import pysam


MODIFICATION_NAMES = {
    ("C", "m"): "5mC",
    ("C", "h"): "5hmC",
    ("C", "f"): "5fC",
    ("C", "c"): "5caC",
    ("A", "a"): "6mA",
    ("G", "o"): "8oxoG",
}


class ModificationCallDecodeError(ValueError):
    """Raised when modified-base calls cannot be mapped safely."""


def _probability_fields(
    ml_value: int,
) -> dict[str, int | float | None]:
    if ml_value == -1:
        return {
            "ml_value": None,
            "probability_lower": None,
            "probability_upper": None,
            "probability_midpoint": None,
        }

    if (
        isinstance(ml_value, bool)
        or not isinstance(ml_value, int)
        or not 0 <= ml_value <= 255
    ):
        raise ModificationCallDecodeError(
            f"Invalid ML value: {ml_value!r}"
        )

    return {
        "ml_value": ml_value,
        "probability_lower": round(ml_value / 256, 6),
        "probability_upper": round((ml_value + 1) / 256, 6),
        "probability_midpoint": round(
            (ml_value + 0.5) / 256,
            6,
        ),
    }


def decode_read_modification_calls(
    read: pysam.AlignedSegment,
    region_start_0: int,
    region_stop_0: int,
) -> dict[str, object]:
    """
    Map MM/ML calls from BAM query positions to reference positions.

    Region coordinates are 0-based and half-open.
    """

    if (
        isinstance(region_start_0, bool)
        or not isinstance(region_start_0, int)
        or region_start_0 < 0
    ):
        raise ModificationCallDecodeError(
            "The internal region start must be a non-negative integer."
        )

    if (
        isinstance(region_stop_0, bool)
        or not isinstance(region_stop_0, int)
        or region_stop_0 <= region_start_0
    ):
        raise ModificationCallDecodeError(
            "The internal region stop must be greater than the start."
        )

    result: dict[str, object] = {
        "read_name": read.query_name,
        "alignment_strand": "-" if read.is_reverse else "+",
        "decoded_call_count": 0,
        "calls_in_region": 0,
        "calls_outside_region": 0,
        "calls_without_reference_position": 0,
        "calls": [],
    }

    query_sequence = read.query_sequence

    if (
        read.is_unmapped
        or query_sequence is None
        or not read.has_tag("MM")
    ):
        return result

    if read.has_tag("MN"):
        mn_length = read.get_tag("MN")

        if mn_length != len(query_sequence):
            raise ModificationCallDecodeError(
                f"Read '{read.query_name}' has MN={mn_length}, "
                f"but its current sequence contains "
                f"{len(query_sequence)} bases. Its MM/ML tags "
                "may be stale after clipping."
            )

    query_to_reference = read.get_reference_positions(
        full_length=True,
    )

    if len(query_to_reference) != len(query_sequence):
        raise ModificationCallDecodeError(
            f"Read '{read.query_name}' has inconsistent query "
            "and reference-position arrays."
        )

    try:
        modified_bases = read.modified_bases
    except (KeyError, TypeError, ValueError) as exc:
        raise ModificationCallDecodeError(
            f"Could not decode MM/ML tags for "
            f"read '{read.query_name}'."
        ) from exc

    if modified_bases is None:
        raise ModificationCallDecodeError(
            f"Could not decode MM/ML tags for "
            f"read '{read.query_name}'."
        )

    if not modified_bases:
        return result

    hp = read.get_tag("HP") if read.has_tag("HP") else None
    ps = read.get_tag("PS") if read.has_tag("PS") else None

    decoded_count = 0
    in_region_count = 0
    outside_region_count = 0
    without_reference_count = 0
    calls: list[dict[str, object]] = []

    for key, encoded_calls in modified_bases.items():
        canonical_base, modification_strand, modification = key

        if modification_strand not in (0, 1):
            raise ModificationCallDecodeError(
                f"Unexpected modification strand "
                f"{modification_strand!r}."
            )

        modification_code = str(modification)
        modification_name = MODIFICATION_NAMES.get(
            (str(canonical_base), modification_code),
            f"{canonical_base}+{modification_code}",
        )

        for query_position_0, ml_value in encoded_calls:
            decoded_count += 1

            if not 0 <= query_position_0 < len(
                query_to_reference
            ):
                raise ModificationCallDecodeError(
                    f"Modified-base position {query_position_0} "
                    f"is outside read '{read.query_name}'."
                )

            reference_position_0 = query_to_reference[
                query_position_0
            ]

            if reference_position_0 is None:
                without_reference_count += 1
                continue

            if not (
                region_start_0
                <= reference_position_0
                < region_stop_0
            ):
                outside_region_count += 1
                continue

            in_region_count += 1

            forward_query_position_0 = query_position_0

            if read.is_reverse:
                forward_query_position_0 = (
                    len(query_sequence)
                    - query_position_0
                    - 1
                )

            call = {
                "read_name": read.query_name,
                "reference_position_0": reference_position_0,
                "reference_position_1": (
                    reference_position_0 + 1
                ),
                "query_position_bam_0": query_position_0,
                "query_position_forward_0": (
                    forward_query_position_0
                ),
                "query_base": query_sequence[
                    query_position_0
                ],
                "canonical_base": str(canonical_base),
                "modification_code": modification_code,
                "modification_name": modification_name,
                "modification_strand": (
                    "+"
                    if modification_strand == 0
                    else "-"
                ),
                "alignment_strand": (
                    "-" if read.is_reverse else "+"
                ),
                "mapping_quality": read.mapping_quality,
                "hp": hp,
                "ps": ps,
            }

            call.update(_probability_fields(ml_value))
            calls.append(call)

    calls.sort(
        key=lambda call: (
            int(call["reference_position_1"]),
            int(call["query_position_bam_0"]),
            str(call["modification_code"]),
        )
    )

    result.update(
        {
            "decoded_call_count": decoded_count,
            "calls_in_region": in_region_count,
            "calls_outside_region": outside_region_count,
            "calls_without_reference_position": (
                without_reference_count
            ),
            "calls": calls,
        }
    )

    return result
