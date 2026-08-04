from pathlib import Path
from typing import Any

import pysam


MAX_RECORDS_TO_SCAN = 10_000

REFERENCE_FINGERPRINTS = {
    "grch38": {
        "label": "GRCh38/hg38",
        "chromosome_1_length": 248_956_422,
        "chromosome_1_aliases": (
            "chr1",
            "1",
            "NC_000001.11",
            "CM000663.2",
        ),
    },
    "t2t-chm13-v2.0": {
        "label": "T2T-CHM13v2.0",
        "chromosome_1_length": 248_387_328,
        "chromosome_1_aliases": (
            "chr1",
            "1",
            "NC_060925.1",
        ),
    },
}


def validate_modbam(
    bam_path: Path,
    bai_path: Path,
    reference: str,
) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []
    checks: list[dict[str, str]] = []

    fingerprint = REFERENCE_FINGERPRINTS.get(reference)

    if fingerprint is None:
        return {
            "accepted": False,
            "errors": ["The selected reference genome is unsupported."],
            "warnings": [],
            "checks": [],
        }

    try:
        with pysam.AlignmentFile(
            str(bam_path),
            "rb",
            index_filename=str(bai_path),
            require_index=True,
        ) as alignment:
            if not alignment.is_bam:
                errors.append("The uploaded alignment is not a BAM file.")
            else:
                checks.append(
                    {
                        "name": "BAM format",
                        "detail": "The file is a readable BAM.",
                    }
                )

            sort_order = (
                alignment.header
                .to_dict()
                .get("HD", {})
                .get("SO")
            )

            if sort_order != "coordinate":
                errors.append(
                    "The BAM header does not declare coordinate sorting."
                )
            else:
                checks.append(
                    {
                        "name": "Sort order",
                        "detail": "The BAM declares SO:coordinate.",
                    }
                )

            alignment.check_index()

            checks.append(
                {
                    "name": "BAM index",
                    "detail": "The uploaded BAI can be opened by pysam.",
                }
            )

            reference_lengths = dict(
                zip(alignment.references, alignment.lengths)
            )

            chromosome_1_name = next(
                (
                    alias
                    for alias in fingerprint["chromosome_1_aliases"]
                    if alias in reference_lengths
                ),
                None,
            )

            if chromosome_1_name is None:
                errors.append(
                    "Chromosome 1 was not found using a recognized "
                    "contig name, so the selected reference cannot be verified."
                )
            else:
                observed_length = reference_lengths[chromosome_1_name]
                expected_length = fingerprint["chromosome_1_length"]

                if observed_length != expected_length:
                    errors.append(
                        f"The BAM does not match {fingerprint['label']}: "
                        f"{chromosome_1_name} has {observed_length:,} bp, "
                        f"but {expected_length:,} bp were expected."
                    )
                else:
                    checks.append(
                        {
                            "name": "Reference genome",
                            "detail": (
                                f"{fingerprint['label']} matched using "
                                f"{chromosome_1_name} "
                                f"({observed_length:,} bp)."
                            ),
                        }
                    )

            alignment.reset()

            records_scanned = 0
            paired_tag_records = 0
            valid_modified_records = 0
            modified_base_calls = 0
            incomplete_tag_records = 0
            mn_mismatches = 0
            parsing_failures = 0

            for read in alignment.fetch(until_eof=True):
                if records_scanned >= MAX_RECORDS_TO_SCAN:
                    break

                records_scanned += 1

                if read.is_unmapped:
                    continue

                has_mm = read.has_tag("MM")
                has_ml = read.has_tag("ML")

                if not has_mm and not has_ml:
                    continue

                if not has_mm or not has_ml:
                    incomplete_tag_records += 1
                    continue

                paired_tag_records += 1

                if read.has_tag("MN"):
                    sequence_length = read.query_length
                    original_length = read.get_tag("MN")

                    if sequence_length != original_length:
                        mn_mismatches += 1
                        continue

                try:
                    modifications = read.modified_bases
                except (TypeError, ValueError):
                    parsing_failures += 1
                    continue

                if not modifications:
                    parsing_failures += 1
                    continue

                call_count = sum(
                    len(calls)
                    for calls in modifications.values()
                )

                if call_count == 0:
                    parsing_failures += 1
                    continue

                valid_modified_records += 1
                modified_base_calls += call_count

            if records_scanned == 0:
                errors.append("The BAM contains no alignment records.")
            elif valid_modified_records == 0:
                errors.append(
                    "No mapped record with valid paired MM/ML "
                    f"tags was found among the first "
                    f"{records_scanned:,} records."
                )
            else:
                record_label = (
                    "record"
                    if valid_modified_records == 1
                    else "records"
                )

                checks.append(
                    {
                        "name": "Modified-base tags",
                        "detail": (
                            f"{valid_modified_records:,} {record_label} with "
                            f"parsable MM/ML tags and "
                            f"{modified_base_calls:,} modified-base calls "
                            f"were detected."
                        ),
                    }
                )

            if incomplete_tag_records:
                warnings.append(
                    f"{incomplete_tag_records:,} scanned records contained "
                    "MM or ML, but not both."
                )

            if mn_mismatches:
                warnings.append(
                    f"{mn_mismatches:,} scanned records had an MN value "
                    "that did not match their current sequence length."
                )

            if parsing_failures:
                warnings.append(
                    f"{parsing_failures:,} scanned records had MM/ML tags "
                    "that pysam could not parse into modified bases."
                )

    except (OSError, ValueError) as exc:
        errors.append(f"The BAM or BAI could not be validated: {exc}")

    return {
        "accepted": not errors,
        "errors": errors,
        "warnings": warnings,
        "checks": checks,
    }
