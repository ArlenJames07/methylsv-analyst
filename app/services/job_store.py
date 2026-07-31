import json
import re
import shutil
import time
from pathlib import Path
from typing import Any
from uuid import uuid4


JOB_ROOT = Path("/tmp/methylsv-jobs")
JOB_TTL_SECONDS = 60 * 60

JOB_ID_PATTERN = re.compile(r"^[0-9a-f]{32}$")


def ensure_job_root() -> Path:
    JOB_ROOT.mkdir(parents=True, exist_ok=True)
    JOB_ROOT.chmod(0o700)
    return JOB_ROOT


def job_path_from_id(job_id: str) -> Path:
    if not JOB_ID_PATTERN.fullmatch(job_id):
        raise ValueError("Invalid analysis job identifier.")

    return ensure_job_root() / job_id


def delete_job(job_id: str) -> None:
    job_path = job_path_from_id(job_id)

    if job_path.is_symlink():
        job_path.unlink(missing_ok=True)
    elif job_path.exists():
        shutil.rmtree(job_path)


def cleanup_expired_jobs() -> None:
    root = ensure_job_root()
    expiration_cutoff = time.time() - JOB_TTL_SECONDS

    for candidate in root.iterdir():
        if not JOB_ID_PATTERN.fullmatch(candidate.name):
            continue

        if candidate.is_symlink() or not candidate.is_dir():
            continue

        try:
            if candidate.stat().st_mtime < expiration_cutoff:
                shutil.rmtree(candidate)
        except FileNotFoundError:
            continue


def create_job() -> tuple[str, Path]:
    cleanup_expired_jobs()

    while True:
        job_id = uuid4().hex
        job_path = ensure_job_root() / job_id

        try:
            job_path.mkdir(mode=0o700)
            return job_id, job_path
        except FileExistsError:
            continue


def write_job_metadata(
    job_path: Path,
    metadata: dict[str, Any],
) -> None:
    payload = dict(metadata)
    payload.setdefault("created_at_epoch", time.time())

    temporary_metadata = job_path / "metadata.json.tmp"
    metadata_path = job_path / "metadata.json"

    temporary_metadata.write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )
    temporary_metadata.chmod(0o600)
    temporary_metadata.replace(metadata_path)


def load_job(
    job_id: str,
) -> tuple[Path, dict[str, Any]] | None:
    cleanup_expired_jobs()

    try:
        job_path = job_path_from_id(job_id)
    except ValueError:
        return None

    if (
        job_path.is_symlink()
        or not job_path.is_dir()
        or not (job_path / "sample.bam").is_file()
        or not (job_path / "sample.bam.bai").is_file()
        or not (job_path / "metadata.json").is_file()
    ):
        return None

    try:
        metadata = json.loads(
            (job_path / "metadata.json").read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError):
        delete_job(job_id)
        return None

    return job_path, metadata