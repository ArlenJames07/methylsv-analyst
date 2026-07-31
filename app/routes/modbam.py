import logging
import shutil
from pathlib import Path
from typing import Annotated

from fastapi import APIRouter, File, Form, Request, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from starlette.concurrency import run_in_threadpool

from app.services.job_store import (
    JOB_TTL_SECONDS,
    create_job,
    delete_job,
    write_job_metadata,
)
from app.services.modbam_validator import validate_modbam

from app.services.job_store import (
    JOB_TTL_SECONDS,
    create_job,
    delete_job,
    load_job,
    write_job_metadata,
)
from app.services.region_reader import (
    RegionValidationError,
    list_contigs,
    read_region,
)

APP_DIR = Path(__file__).resolve().parents[1]
COPY_BUFFER_SIZE = 8 * 1024 * 1024
REGION_TABLE_LIMIT = 200

LOGGER = logging.getLogger(__name__)

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

    destination.chmod(0o600)
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

    job_id: str | None = None

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
            staged_job_id: str | None = None

            try:
                staged_job_id, job_path = await run_in_threadpool(
                    create_job
                )

                bam_path = job_path / "sample.bam"
                bai_path = job_path / "sample.bam.bai"

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

                if errors:
                    await run_in_threadpool(
                        delete_job,
                        staged_job_id,
                    )
                    staged_job_id = None

                else:
                    await run_in_threadpool(
                        write_job_metadata,
                        job_path,
                        {
                            "reference": reference,
                            "reference_label": (
                                REFERENCE_OPTIONS[reference]
                            ),
                            "bam_original_name": (
                                bam_file.filename or "sample.bam"
                            ),
                            "bai_original_name": (
                                bai_file.filename or "sample.bam.bai"
                            ),
                            "validation_checks": checks,
                            "ttl_seconds": JOB_TTL_SECONDS,
                        },
                    )

                    job_id = staged_job_id

            except Exception:
                if staged_job_id is not None:
                    try:
                        await run_in_threadpool(
                            delete_job,
                            staged_job_id,
                        )
                    except OSError:
                        LOGGER.exception(
                            "Temporary job cleanup failed."
                        )

                LOGGER.exception(
                    "The uploaded modBAM could not be staged."
                )

                errors.append(
                    "The uploaded files could not be staged "
                    "for analysis."
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
                "job_id": job_id,
                "job_ttl_minutes": JOB_TTL_SECONDS // 60,
            },
            status_code=200 if accepted else 400,
        )

    finally:
        for upload in [*bam_files, *bai_files]:
            await upload.close()


async def load_available_job(
    job_id: str,
) -> tuple[Path, dict[str, object]] | None:
    try:
        return await run_in_threadpool(load_job, job_id)
    except ValueError:
        # Invalid job-ID format.
        return None


@router.get(
    "/{job_id}/region",
    response_class=HTMLResponse,
    name="modbam_region_form",
)
async def region_form(
    request: Request,
    job_id: str,
):
    loaded_job = await load_available_job(job_id)

    if loaded_job is None:
        return templates.TemplateResponse(
            request=request,
            name="modbam_region.html",
            context={
                "job_available": False,
                "job_id": job_id,
            },
            status_code=404,
        )

    job_path, metadata = loaded_job
    bam_path = job_path / "sample.bam"
    bai_path = job_path / "sample.bam.bai"

    try:
        contigs = await run_in_threadpool(
            list_contigs,
            bam_path,
            bai_path,
        )
    except (FileNotFoundError, OSError):
        LOGGER.exception(
            "The regional input files are unavailable for job %s.",
            job_id,
        )

        return templates.TemplateResponse(
            request=request,
            name="modbam_region.html",
            context={
                "job_available": False,
                "job_id": job_id,
            },
            status_code=404,
        )

    errors: list[str] = []

    if contigs:
        selected_contig = str(contigs[0]["name"])
        start_value: int | str = 1
        end_value: int | str = min(
            int(contigs[0]["length"]),
            1_000,
        )
    else:
        selected_contig = ""
        start_value = ""
        end_value = ""
        errors.append(
            "The BAM header does not declare any reference contigs."
        )

    return templates.TemplateResponse(
        request=request,
        name="modbam_region.html",
        context={
            "job_available": True,
            "job_id": job_id,
            "metadata": metadata,
            "contigs": contigs,
            "selected_contig": selected_contig,
            "start_value": start_value,
            "end_value": end_value,
            "errors": errors,
            "result": None,
            "table_limit": REGION_TABLE_LIMIT,
        },
        status_code=200 if not errors else 500,
    )


@router.post(
    "/{job_id}/region",
    response_class=HTMLResponse,
    name="analyze_modbam_region",
)
async def analyze_region(
    request: Request,
    job_id: str,
    contig: Annotated[str, Form()],
    start_1: Annotated[str, Form()],
    end_1: Annotated[str, Form()],
):
    loaded_job = await load_available_job(job_id)

    if loaded_job is None:
        return templates.TemplateResponse(
            request=request,
            name="modbam_region.html",
            context={
                "job_available": False,
                "job_id": job_id,
            },
            status_code=404,
        )

    job_path, metadata = loaded_job
    bam_path = job_path / "sample.bam"
    bai_path = job_path / "sample.bam.bai"

    try:
        contigs = await run_in_threadpool(
            list_contigs,
            bam_path,
            bai_path,
        )
    except (FileNotFoundError, OSError):
        LOGGER.exception(
            "The regional input files are unavailable for job %s.",
            job_id,
        )

        return templates.TemplateResponse(
            request=request,
            name="modbam_region.html",
            context={
                "job_available": False,
                "job_id": job_id,
            },
            status_code=404,
        )

    errors: list[str] = []
    result: dict[str, object] | None = None

    try:
        parsed_start = int(start_1.strip())
    except (AttributeError, TypeError, ValueError):
        parsed_start = None
        errors.append("The region start must be an integer.")

    try:
        parsed_end = int(end_1.strip())
    except (AttributeError, TypeError, ValueError):
        parsed_end = None
        errors.append("The region end must be an integer.")

    if (
        not errors
        and parsed_start is not None
        and parsed_end is not None
    ):
        try:
            result = await run_in_threadpool(
                read_region,
                bam_path,
                bai_path,
                contig,
                parsed_start,
                parsed_end,
            )
        except RegionValidationError as exc:
            errors.append(str(exc))
        except FileNotFoundError:
            return templates.TemplateResponse(
                request=request,
                name="modbam_region.html",
                context={
                    "job_available": False,
                    "job_id": job_id,
                },
                status_code=404,
            )

    return templates.TemplateResponse(
        request=request,
        name="modbam_region.html",
        context={
            "job_available": True,
            "job_id": job_id,
            "metadata": metadata,
            "contigs": contigs,
            "selected_contig": contig,
            "start_value": start_1,
            "end_value": end_1,
            "errors": errors,
            "result": result,
            "table_limit": REGION_TABLE_LIMIT,
        },
        status_code=400 if errors else 200,
    )            