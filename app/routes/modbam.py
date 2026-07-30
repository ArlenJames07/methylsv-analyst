import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Annotated

from fastapi import APIRouter, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from starlette.concurrency import run_in_threadpool

from app.services.modbam_validator import validate_modbam


APP_DIR = Path(__file__).resolve().parents[1]
COPY_BUFFER_SIZE = 8 * 1024 * 1024

router = APIRouter(
    prefix="/workflows/modbam",
    tags=["modBAM"],
)

templates = Jinja2Templates(directory=APP_DIR / "templates")

REFERENCE_OPTIONS = {
    "grch38": "GRCh38/hg38",
    "t2t-chm13-v2.0": "T2T-CHM13v2.0",
}


def file_summary(upload: UploadFile) -> dict[str, str]:
    size = upload.size

    return {
        "name": upload.filename or "Unnamed file",
        "size": "Unknown" if size is None else f"{size:,} bytes",
    }


async def read_prefix(upload: UploadFile, length: int) -> bytes:
    prefix = await upload.read(length)
    await upload.seek(0)
    return prefix


def copy_upload_to_path(
    upload: UploadFile,
    destination: Path,
) -> None:
    upload.file.seek(0)

    with destination.open("wb") as output_file:
        shutil.copyfileobj(
            upload.file,
            output_file,
            length=COPY_BUFFER_SIZE,
        )

    upload.file.seek(0)


@router.get(
    "/upload",
    response_class=HTMLResponse,
    name="modbam_upload_form",
)
def upload_form(request: Request):
    return templates.TemplateResponse(
        request=request,
        name="modbam_upload.html",
        context={"reference_options": REFERENCE_OPTIONS},
    )


@router.post(
    "/validate",
    response_class=HTMLResponse,
    name="validate_modbam_upload",
)
async def validate_upload(
    request: Request,
    reference: Annotated[str, Form()],
    bam_files: Annotated[list[UploadFile], File()],
    bai_files: Annotated[list[UploadFile], File()],
):
    errors: list[str] = []
    warnings: list[str] = []
    checks: list[dict[str, str]] = []

    bam_file = bam_files[0] if len(bam_files) == 1 else None
    bai_file = bai_files[0] if len(bai_files) == 1 else None

    try:
        if len(bam_files) != 1:
            errors.append("Exactly one BAM file is required.")

        if len(bai_files) != 1:
            errors.append("Exactly one BAI index file is required.")

        if reference not in REFERENCE_OPTIONS:
            errors.append(
                "The selected reference genome is not supported."
            )

        if bam_file is not None:
            bam_name = bam_file.filename or ""

            if not bam_name.lower().endswith(".bam"):
                errors.append(
                    "The alignment file must end in .bam."
                )

            bam_prefix = await read_prefix(bam_file, 2)

            if bam_prefix != b"\x1f\x8b":
                errors.append(
                    "The BAM does not have the expected compressed "
                    "binary-file signature."
                )

        if bai_file is not None:
            bai_name = bai_file.filename or ""

            if not bai_name.lower().endswith(".bai"):
                errors.append(
                    "The index file must end in .bai."
                )

            bai_prefix = await read_prefix(bai_file, 4)

            if bai_prefix != b"BAI\x01":
                errors.append(
                    "The uploaded index is not a recognizable BAI file."
                )

        if not errors and bam_file is not None and bai_file is not None:
            try:
                with TemporaryDirectory(
                    prefix="methylsv-modbam-"
                ) as temporary_directory:
                    temporary_path = Path(temporary_directory)
                    bam_path = temporary_path / "sample.bam"
                    bai_path = temporary_path / "sample.bam.bai"

                    await run_in_threadpool(
                        copy_upload_to_path,
                        bam_file,
                        bam_path,
                    )

                    await run_in_threadpool(
                        copy_upload_to_path,
                        bai_file,
                        bai_path,
                    )

                    result = await run_in_threadpool(
                        validate_modbam,
                        bam_path,
                        bai_path,
                        reference,
                    )

                    errors.extend(result["errors"])
                    warnings.extend(result["warnings"])
                    checks.extend(result["checks"])

            except OSError:
                errors.append(
                    "The uploaded files could not be staged "
                    "for validation."
                )

        accepted = not errors

        files = []

        if bam_file is not None:
            files.append(file_summary(bam_file))

        if bai_file is not None:
            files.append(file_summary(bai_file))

        return templates.TemplateResponse(
            request=request,
            name="modbam_validation_result.html",
            context={
                "accepted": accepted,
                "errors": errors,
                "warnings": warnings,
                "checks": checks,
                "files": files,
                "reference": REFERENCE_OPTIONS.get(
                    reference,
                    reference,
                ),
            },
            status_code=200 if accepted else 400,
        )
    finally:
        for upload in [*bam_files, *bai_files]:
            await upload.close()