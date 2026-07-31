from pathlib import Path

import pysam


MAX_REGION_BP = 1_000_000
MAX_READS = 10_000


class RegionValidationError(ValueError):
    """Raised when a requested genomic region is invalid."""


def _check_input_files(
    bam_path: Path,
    bai_path: Path,
) -> None:
    if not bam_path.is_file():
        raise FileNotFoundError(
            f"BAM file not found: {bam_path}"
        )

    if not bai_path.is_file():
        raise FileNotFoundError(
            f"BAI file not found: {bai_path}"
        )


def list_contigs(
    bam_path: Path,
    bai_path: Path,
) -> list[dict[str, str | int]]:
    """Return the contigs declared in the BAM header."""

    _check_input_files(bam_path, bai_path)

    with pysam.AlignmentFile(
        str(bam_path),
        "rb",
        index_filename=str(bai_path),
        require_index=True,
    ) as bam:
        return [
            {
                "name": name,
                "length": length,
            }
            for name, length in zip(
                bam.references,
                bam.lengths,
                strict=True,
            )
        ]


def read_region(
    bam_path: Path,
    bai_path: Path,
    contig: str,
    start_1: int,
    end_1: int,
) -> dict[str, object]:
    """
    Retrieve reads overlapping a genomic interval.

    Public coordinates:
        1-based and inclusive.

    pysam coordinates:
        0-based and half-open.

    Example:
        chr1:101-118 becomes fetch("chr1", 100, 118).
    """

    _check_input_files(bam_path, bai_path)

    if not contig or not contig.strip():
        raise RegionValidationError(
            "A chromosome or contig is required."
        )

    contig = contig.strip()

    if isinstance(start_1, bool) or not isinstance(start_1, int):
        raise RegionValidationError(
            "The region start must be an integer."
        )

    if isinstance(end_1, bool) or not isinstance(end_1, int):
        raise RegionValidationError(
            "The region end must be an integer."
        )

    if start_1 < 1:
        raise RegionValidationError(
            "The region start must be at least 1."
        )

    if end_1 < start_1:
        raise RegionValidationError(
            "The region end must be greater than or equal "
            "to the start."
        )

    region_length = end_1 - start_1 + 1

    if region_length > MAX_REGION_BP:
        raise RegionValidationError(
            f"The requested region contains {region_length:,} bp. "
            f"The current limit is {MAX_REGION_BP:,} bp."
        )

    with pysam.AlignmentFile(
        str(bam_path),
        "rb",
        index_filename=str(bai_path),
        require_index=True,
    ) as bam:
        if contig not in bam.references:
            available_preview = ", ".join(bam.references[:10])

            raise RegionValidationError(
                f"Contig '{contig}' is not present in the BAM. "
                f"Available contigs begin with: {available_preview}"
            )

        contig_length = bam.get_reference_length(contig)

        if end_1 > contig_length:
            raise RegionValidationError(
                f"The requested end ({end_1:,}) exceeds the "
                f"length of {contig} ({contig_length:,} bp)."
            )

        # Convert 1-based inclusive coordinates to the
        # 0-based half-open convention required by pysam.
        start_0 = start_1 - 1
        stop_0 = end_1

        reads: list[dict[str, object]] = []
        truncated = False
        reads_with_modifications = 0
        modified_calls_in_overlapping_reads = 0
        modification_decode_errors = 0

        iterator = bam.fetch(
            contig=contig,
            start=start_0,
            stop=stop_0,
        )

        for read_number, read in enumerate(iterator):
            if read_number >= MAX_READS:
                truncated = True
                break

            has_mm = read.has_tag("MM")
            has_ml = read.has_tag("ML")

            try:
                modified_bases = read.modified_bases
            except (KeyError, ValueError):
                modified_bases = None
                modification_decode_errors += 1

            modification_count = 0

            if modified_bases:
                modification_count = sum(
                    len(calls)
                    for calls in modified_bases.values()
                )

            if modification_count > 0:
                reads_with_modifications += 1
                modified_calls_in_overlapping_reads += (
                    modification_count
                )

            hp = read.get_tag("HP") if read.has_tag("HP") else None
            ps = read.get_tag("PS") if read.has_tag("PS") else None

            read_start_1 = (
                None
                if read.reference_start is None
                else read.reference_start + 1
            )

            read_end_1 = read.reference_end

            reads.append(
                {
                    "name": read.query_name,
                    "start_1": read_start_1,
                    "end_1": read_end_1,
                    "mapping_quality": read.mapping_quality,
                    "strand": "-" if read.is_reverse else "+",
                    "overlap_bp": read.get_overlap(
                        start_0,
                        stop_0,
                    ),
                    "has_mm": has_mm,
                    "has_ml": has_ml,
                    "modified_calls_in_read": modification_count,
                    "hp": hp,
                    "ps": ps,
                    "is_secondary": read.is_secondary,
                    "is_supplementary": read.is_supplementary,
                }
            )

    return {
        "region_label": f"{contig}:{start_1}-{end_1}",
        "contig": contig,
        "contig_length": contig_length,
        "start_1": start_1,
        "end_1": end_1,
        "region_length_bp": region_length,
        "fetch_start_0": start_0,
        "fetch_stop_0": stop_0,
        "read_count": len(reads),
        "reads_with_modifications": reads_with_modifications,
        "modified_calls_in_overlapping_reads": (
            modified_calls_in_overlapping_reads
        ),
        "modification_decode_errors": (
            modification_decode_errors
        ),
        "truncated": truncated,
        "read_limit": MAX_READS,
        "reads": reads,
    }