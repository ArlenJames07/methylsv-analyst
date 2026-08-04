from array import array
from pathlib import Path
from typing import Iterable

import pysam


GRCH38_CHR1_LENGTH = 248_956_422
T2T_CHM13_CHR1_LENGTH = 248_387_328

CANONICAL_SEQUENCE = "ACGTCGCGTACGTCGCGT"
CANONICAL_MM = "C+m,0,0,0,0,0,0;"
CANONICAL_ML = (230, 220, 210, 200, 190, 180)

AUTO_MN = object()


def make_read(
    name: str = "synthetic_modbam_read",
    start_0: int = 100,
    sequence: str = CANONICAL_SEQUENCE,
    *,
    flag: int = 0,
    cigar: list[tuple[int, int]] | None = None,
    mm: str | None = CANONICAL_MM,
    ml: Iterable[int] | None = CANONICAL_ML,
    mn: object | int | None = AUTO_MN,
    mapping_quality: int = 60,
    hp: int | None = None,
    ps: int | None = None,
) -> pysam.AlignedSegment:
    """Build one deterministic synthetic alignment record."""

    read = pysam.AlignedSegment()
    read.query_name = name
    read.query_sequence = sequence
    read.flag = flag
    read.reference_id = 0
    read.reference_start = start_0
    read.mapping_quality = mapping_quality
    read.cigartuples = cigar or [(0, len(sequence))]
    read.query_qualities = pysam.qualitystring_to_array(
        "I" * len(sequence)
    )

    if mm is not None:
        read.set_tag("MM", mm, value_type="Z")

    if ml is not None:
        read.set_tag("ML", array("B", ml))

    if mn is AUTO_MN:
        read.set_tag("MN", len(sequence), value_type="i")
    elif mn is not None:
        read.set_tag("MN", int(mn), value_type="i")

    if hp is not None:
        read.set_tag("HP", hp, value_type="i")

    if ps is not None:
        read.set_tag("PS", ps, value_type="i")

    return read


def write_bam(
    directory: Path,
    reads: Iterable[pysam.AlignedSegment],
    *,
    name: str = "sample",
    chromosome_1_length: int = GRCH38_CHR1_LENGTH,
    sort_order: str = "coordinate",
    create_index: bool = True,
) -> tuple[Path, Path]:
    """Write a synthetic BAM and optionally create its BAI."""

    bam_path = directory / f"{name}.bam"
    bai_path = Path(f"{bam_path}.bai")

    header = {
        "HD": {
            "VN": "1.6",
            "SO": sort_order,
        },
        "SQ": [
            {
                "SN": "chr1",
                "LN": chromosome_1_length,
            }
        ],
    }

    sorted_reads = sorted(
        reads,
        key=lambda read: (
            read.reference_id,
            read.reference_start,
            read.query_name,
        ),
    )

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for read in sorted_reads:
            bam.write(read)

    if create_index:
        pysam.index(str(bam_path))

    return bam_path, bai_path
