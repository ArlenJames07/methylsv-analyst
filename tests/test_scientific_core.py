import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

from app.services.modbam_validator import validate_modbam
from app.services.modification_calls import (
    ModificationCallDecodeError,
    _probability_fields,
    decode_read_modification_calls,
)
from app.services.region_reader import (
    MAX_REGION_BP,
    RegionValidationError,
    list_contigs,
    read_region,
)
from tests.synthetic_modbam import (
    CANONICAL_ML,
    CANONICAL_MM,
    CANONICAL_SEQUENCE,
    GRCH38_CHR1_LENGTH,
    T2T_CHM13_CHR1_LENGTH,
    make_read,
    write_bam,
)


class TemporaryBamTestCase(unittest.TestCase):
    def setUp(self) -> None:
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        self.temporary_path = Path(temporary_directory.name)


class ModificationCallTests(unittest.TestCase):
    def test_probability_fields_cover_boundaries(self) -> None:
        self.assertEqual(
            _probability_fields(-1),
            {
                "ml_value": None,
                "probability_lower": None,
                "probability_upper": None,
                "probability_midpoint": None,
            },
        )
        self.assertEqual(
            _probability_fields(0),
            {
                "ml_value": 0,
                "probability_lower": 0.0,
                "probability_upper": 0.003906,
                "probability_midpoint": 0.001953,
            },
        )
        self.assertEqual(
            _probability_fields(255),
            {
                "ml_value": 255,
                "probability_lower": 0.996094,
                "probability_upper": 1.0,
                "probability_midpoint": 0.998047,
            },
        )

        for invalid_value in (True, -2, 256, 1.5):
            with self.subTest(invalid_value=invalid_value):
                with self.assertRaises(ModificationCallDecodeError):
                    _probability_fields(invalid_value)  # type: ignore[arg-type]

    def test_forward_read_is_filtered_to_half_open_region(self) -> None:
        read = make_read(hp=1, ps=42)

        result = decode_read_modification_calls(
            read,
            region_start_0=104,
            region_stop_0=111,
        )

        self.assertEqual(result["decoded_call_count"], 6)
        self.assertEqual(result["calls_in_region"], 3)
        self.assertEqual(result["calls_outside_region"], 3)
        self.assertEqual(result["calls_without_reference_position"], 0)
        self.assertEqual(
            result["decoded_call_count"],
            result["calls_in_region"]
            + result["calls_outside_region"]
            + result["calls_without_reference_position"],
        )

        calls = result["calls"]
        self.assertEqual(len(calls), result["calls_in_region"])
        self.assertEqual(
            [call["reference_position_1"] for call in calls],
            [105, 107, 111],
        )
        self.assertEqual(
            [call["query_position_bam_0"] for call in calls],
            [4, 6, 10],
        )
        self.assertEqual(
            [call["ml_value"] for call in calls],
            [220, 210, 200],
        )
        self.assertEqual(calls[0]["probability_midpoint"], 0.861328)
        self.assertTrue(
            all(call["modification_name"] == "5mC" for call in calls)
        )
        self.assertTrue(
            all(call["modification_strand"] == "+" for call in calls)
        )
        self.assertTrue(all(call["hp"] == 1 for call in calls))
        self.assertTrue(all(call["ps"] == 42 for call in calls))

    def test_reverse_read_preserves_orientation_and_ml_order(self) -> None:
        read = make_read(
            name="reverse_read",
            flag=16,
            mapping_quality=42,
            hp=2,
            ps=99,
        )

        result = decode_read_modification_calls(read, 100, 118)
        calls = result["calls"]

        self.assertEqual(
            [call["reference_position_1"] for call in calls],
            [103, 106, 108, 112, 115, 117],
        )
        self.assertEqual(
            [call["query_position_bam_0"] for call in calls],
            [2, 5, 7, 11, 14, 16],
        )
        self.assertEqual(
            [call["query_position_forward_0"] for call in calls],
            [15, 12, 10, 6, 3, 1],
        )
        self.assertEqual(
            [call["ml_value"] for call in calls],
            [180, 190, 200, 210, 220, 230],
        )
        self.assertTrue(all(call["query_base"] == "G" for call in calls))
        self.assertTrue(
            all(call["canonical_base"] == "C" for call in calls)
        )
        self.assertTrue(
            all(call["modification_strand"] == "-" for call in calls)
        )
        self.assertTrue(
            all(call["alignment_strand"] == "-" for call in calls)
        )

    def test_soft_clipped_call_has_no_reference_position(self) -> None:
        read = make_read(
            sequence="CCACCC",
            start_0=300,
            cigar=[(4, 2), (0, 4)],
            mm="C+m,0,1;",
            ml=(64, 192),
        )

        result = decode_read_modification_calls(read, 300, 304)

        self.assertEqual(result["decoded_call_count"], 2)
        self.assertEqual(result["calls_in_region"], 1)
        self.assertEqual(result["calls_outside_region"], 0)
        self.assertEqual(result["calls_without_reference_position"], 1)
        self.assertEqual(
            result["decoded_call_count"],
            result["calls_in_region"]
            + result["calls_outside_region"]
            + result["calls_without_reference_position"],
        )
        self.assertEqual(
            result["calls"],
            [
                {
                    "read_name": "synthetic_modbam_read",
                    "reference_position_0": 301,
                    "reference_position_1": 302,
                    "query_position_bam_0": 3,
                    "query_position_forward_0": 3,
                    "query_base": "C",
                    "canonical_base": "C",
                    "modification_code": "m",
                    "modification_name": "5mC",
                    "modification_strand": "+",
                    "alignment_strand": "+",
                    "mapping_quality": 60,
                    "hp": None,
                    "ps": None,
                    "ml_value": 192,
                    "probability_lower": 0.75,
                    "probability_upper": 0.753906,
                    "probability_midpoint": 0.751953,
                }
            ],
        )

    def test_empty_annotation_is_not_a_decode_error(self) -> None:
        read = make_read(mm="C+m;", ml=())

        result = decode_read_modification_calls(read, 100, 118)

        self.assertEqual(result["decoded_call_count"], 0)
        self.assertEqual(result["calls"], [])

    def test_missing_mm_is_empty_and_missing_ml_has_unknown_probability(self) -> None:
        unmodified = make_read(mm=None, ml=None, mn=None)
        unmodified_result = decode_read_modification_calls(
            unmodified,
            100,
            118,
        )
        self.assertEqual(unmodified_result["decoded_call_count"], 0)
        self.assertEqual(unmodified_result["calls"], [])

        no_ml = make_read(ml=None)
        no_ml_result = decode_read_modification_calls(no_ml, 100, 118)
        self.assertEqual(no_ml_result["decoded_call_count"], 6)
        self.assertTrue(
            all(call["ml_value"] is None for call in no_ml_result["calls"])
        )
        self.assertTrue(
            all(
                call["probability_midpoint"] is None
                for call in no_ml_result["calls"]
            )
        )

    def test_malformed_annotation_is_a_decode_error(self) -> None:
        read = make_read(mm="C+m,99;", ml=(128,))
        previous_verbosity = pysam.set_verbosity(0)

        try:
            with self.assertRaises(ModificationCallDecodeError):
                decode_read_modification_calls(read, 100, 118)
        finally:
            pysam.set_verbosity(previous_verbosity)

    def test_stale_mn_and_invalid_regions_are_rejected(self) -> None:
        stale_read = make_read(
            sequence="ACGCCC",
            cigar=[(5, 2), (0, 6)],
            mm="C+m,0,0,0,0;",
            ml=(220, 210, 200, 190),
            mn=8,
        )

        with self.assertRaisesRegex(
            ModificationCallDecodeError,
            "may be stale after clipping",
        ):
            decode_read_modification_calls(stale_read, 100, 118)

        invalid_regions = (
            (-1, 10),
            (True, 10),
            (10, True),
            (10, 10),
            (11, 10),
        )

        for start_0, stop_0 in invalid_regions:
            with self.subTest(start_0=start_0, stop_0=stop_0):
                with self.assertRaises(ModificationCallDecodeError):
                    decode_read_modification_calls(
                        make_read(),
                        start_0,  # type: ignore[arg-type]
                        stop_0,  # type: ignore[arg-type]
                    )

    def test_modified_bases_exceptions_are_wrapped(self) -> None:
        class BrokenRead:
            query_name = "broken_read"
            query_sequence = "C"
            is_reverse = False
            is_unmapped = False

            def __init__(self, exception: Exception) -> None:
                self.exception = exception

            def has_tag(self, tag: str) -> bool:
                return tag == "MM"

            def get_reference_positions(
                self,
                full_length: bool,
            ) -> list[int]:
                self.full_length = full_length
                return [100]

            @property
            def modified_bases(self) -> object:
                raise self.exception

        for exception in (KeyError("MM"), TypeError("MM"), ValueError("MM")):
            with self.subTest(exception=type(exception).__name__):
                with self.assertRaises(ModificationCallDecodeError):
                    decode_read_modification_calls(
                        BrokenRead(exception),  # type: ignore[arg-type]
                        100,
                        101,
                    )

        class NoneModifiedBasesRead(BrokenRead):
            @property
            def modified_bases(self) -> None:
                return None

        with self.assertRaises(ModificationCallDecodeError):
            decode_read_modification_calls(
                NoneModifiedBasesRead(ValueError()),  # type: ignore[arg-type]
                100,
                101,
            )


class RegionReaderTests(TemporaryBamTestCase):
    def setUp(self) -> None:
        super().setUp()
        self.bam_path, self.bai_path = write_bam(
            self.temporary_path,
            [make_read(hp=1, ps=42)],
        )

    def test_contigs_and_full_region_summary(self) -> None:
        self.assertEqual(
            list_contigs(self.bam_path, self.bai_path),
            [{"name": "chr1", "length": GRCH38_CHR1_LENGTH}],
        )

        result = read_region(
            self.bam_path,
            self.bai_path,
            "chr1",
            101,
            118,
        )

        self.assertEqual(result["region_length_bp"], 18)
        self.assertEqual(result["fetch_start_0"], 100)
        self.assertEqual(result["fetch_stop_0"], 118)
        self.assertEqual(result["read_count"], 1)
        self.assertEqual(result["reads_with_modifications"], 1)
        self.assertEqual(result["reads_with_modifications_in_region"], 1)
        self.assertEqual(result["modified_calls_in_overlapping_reads"], 6)
        self.assertEqual(result["modified_calls_in_region"], 6)
        self.assertEqual(
            result["modified_calls_outside_region_in_overlapping_reads"],
            0,
        )
        self.assertEqual(result["modified_calls_without_reference_position"], 0)
        self.assertEqual(result["modification_decode_errors"], 0)
        self.assertEqual(
            [
                call["reference_position_1"]
                for call in result["modification_calls"]
            ],
            [102, 105, 107, 111, 114, 116],
        )

        read = result["reads"][0]
        self.assertEqual(read["start_1"], 101)
        self.assertEqual(read["end_1"], 118)
        self.assertEqual(read["modified_calls_in_read"], 6)
        self.assertEqual(read["modified_calls_in_region"], 6)
        self.assertEqual(read["hp"], 1)
        self.assertEqual(read["ps"], 42)

    def test_partial_and_one_base_regions_use_inclusive_coordinates(self) -> None:
        partial = read_region(
            self.bam_path,
            self.bai_path,
            "chr1",
            105,
            111,
        )

        self.assertEqual(partial["region_length_bp"], 7)
        self.assertEqual(partial["modified_calls_in_overlapping_reads"], 6)
        self.assertEqual(partial["modified_calls_in_region"], 3)
        self.assertEqual(
            partial["modified_calls_outside_region_in_overlapping_reads"],
            3,
        )
        self.assertEqual(
            [
                call["reference_position_1"]
                for call in partial["modification_calls"]
            ],
            [105, 107, 111],
        )
        self.assertEqual(partial["reads"][0]["modified_calls_in_region"], 3)

        one_base = read_region(
            self.bam_path,
            self.bai_path,
            "chr1",
            102,
            102,
        )
        self.assertEqual(one_base["region_length_bp"], 1)
        self.assertEqual(one_base["modified_calls_in_region"], 1)
        self.assertEqual(
            one_base["modification_calls"][0]["reference_position_1"],
            102,
        )

    def test_reverse_read_survives_bam_round_trip(self) -> None:
        reverse_bam, reverse_bai = write_bam(
            self.temporary_path,
            [
                make_read(
                    name="reverse_read",
                    flag=16,
                    mapping_quality=42,
                    hp=2,
                    ps=99,
                )
            ],
            name="reverse",
        )

        result = read_region(
            reverse_bam,
            reverse_bai,
            "chr1",
            101,
            118,
        )
        calls = result["modification_calls"]

        self.assertEqual(result["read_count"], 1)
        self.assertEqual(result["modified_calls_in_region"], 6)
        self.assertEqual(
            [call["reference_position_1"] for call in calls],
            [103, 106, 108, 112, 115, 117],
        )
        self.assertEqual(
            [call["query_position_bam_0"] for call in calls],
            [2, 5, 7, 11, 14, 16],
        )
        self.assertEqual(
            [call["query_position_forward_0"] for call in calls],
            [15, 12, 10, 6, 3, 1],
        )
        self.assertEqual(
            [call["ml_value"] for call in calls],
            [180, 190, 200, 210, 220, 230],
        )
        self.assertEqual(result["reads"][0]["strand"], "-")
        self.assertEqual(result["reads"][0]["hp"], 2)
        self.assertEqual(result["reads"][0]["ps"], 99)

    def test_decode_errors_are_counted_without_silent_calls(self) -> None:
        error_bam, error_bai = write_bam(
            self.temporary_path,
            [make_read(name="stale_read", mn=len(CANONICAL_SEQUENCE) + 1)],
            name="stale",
        )

        result = read_region(
            error_bam,
            error_bai,
            "chr1",
            101,
            118,
        )

        self.assertEqual(result["read_count"], 1)
        self.assertEqual(result["modification_decode_errors"], 1)
        self.assertEqual(result["modified_calls_in_overlapping_reads"], 0)
        self.assertEqual(result["modified_calls_in_region"], 0)
        self.assertEqual(result["modification_calls"], [])
        self.assertTrue(result["reads"][0]["modification_decode_error"])

    def test_malformed_annotations_are_counted_as_decode_errors(self) -> None:
        malformed_bam, malformed_bai = write_bam(
            self.temporary_path,
            [
                make_read(
                    name="malformed_read",
                    mm="C+m,99;",
                    ml=(128,),
                )
            ],
            name="malformed",
        )
        previous_verbosity = pysam.set_verbosity(0)

        try:
            result = read_region(
                malformed_bam,
                malformed_bai,
                "chr1",
                101,
                118,
            )
        finally:
            pysam.set_verbosity(previous_verbosity)

        self.assertEqual(result["read_count"], 1)
        self.assertEqual(result["modification_decode_errors"], 1)
        self.assertEqual(result["modified_calls_in_overlapping_reads"], 0)
        self.assertEqual(result["modification_calls"], [])
        self.assertTrue(result["reads"][0]["modification_decode_error"])

    def test_truncation_flag_uses_configured_read_limit(self) -> None:
        truncated_bam, truncated_bai = write_bam(
            self.temporary_path,
            [
                make_read(name="read_1", start_0=100),
                make_read(name="read_2", start_0=105),
            ],
            name="truncated",
        )

        with patch("app.services.region_reader.MAX_READS", 1):
            result = read_region(
                truncated_bam,
                truncated_bai,
                "chr1",
                101,
                130,
            )

        self.assertEqual(result["read_count"], 1)
        self.assertTrue(result["truncated"])
        self.assertEqual(result["read_limit"], 1)

    def test_invalid_inputs_and_bounds_are_rejected(self) -> None:
        missing_path = self.temporary_path / "missing.bam"

        with self.assertRaises(FileNotFoundError):
            read_region(missing_path, self.bai_path, "chr1", 1, 1)

        with self.assertRaises(FileNotFoundError):
            read_region(self.bam_path, missing_path, "chr1", 1, 1)

        maximum_region = read_region(
            self.bam_path,
            self.bai_path,
            "chr1",
            1,
            MAX_REGION_BP,
        )
        self.assertEqual(maximum_region["region_length_bp"], MAX_REGION_BP)

        invalid_requests = (
            ("", 1, 1),
            ("   ", 1, 1),
            ("chr1", True, 1),
            ("chr1", 1, True),
            ("chr1", 0, 1),
            ("chr1", 10, 9),
            ("chr1", 1, MAX_REGION_BP + 1),
            ("chr2", 1, 1),
            (
                "chr1",
                GRCH38_CHR1_LENGTH,
                GRCH38_CHR1_LENGTH + 1,
            ),
        )

        for contig, start_1, end_1 in invalid_requests:
            with self.subTest(
                contig=contig,
                start_1=start_1,
                end_1=end_1,
            ):
                with self.assertRaises(RegionValidationError):
                    read_region(
                        self.bam_path,
                        self.bai_path,
                        contig,
                        start_1,  # type: ignore[arg-type]
                        end_1,  # type: ignore[arg-type]
                    )


class ModbamValidatorTests(TemporaryBamTestCase):
    def test_valid_grch38_bam_passes_all_checks(self) -> None:
        bam_path, bai_path = write_bam(
            self.temporary_path,
            [make_read()],
        )

        result = validate_modbam(bam_path, bai_path, "grch38")

        self.assertTrue(result["accepted"])
        self.assertEqual(result["errors"], [])
        self.assertEqual(result["warnings"], [])
        self.assertEqual(
            [check["name"] for check in result["checks"]],
            [
                "BAM format",
                "Sort order",
                "BAM index",
                "Reference genome",
                "Modified-base tags",
            ],
        )
        self.assertIn("6 modified-base calls", result["checks"][-1]["detail"])

    def test_reference_mismatch_and_unsupported_reference_fail(self) -> None:
        bam_path, bai_path = write_bam(
            self.temporary_path,
            [make_read()],
            chromosome_1_length=T2T_CHM13_CHR1_LENGTH,
        )

        mismatch = validate_modbam(bam_path, bai_path, "grch38")
        self.assertFalse(mismatch["accepted"])
        self.assertTrue(
            any("does not match GRCh38/hg38" in error for error in mismatch["errors"])
        )

        unsupported = validate_modbam(
            self.temporary_path / "missing.bam",
            self.temporary_path / "missing.bai",
            "unsupported",
        )
        self.assertEqual(
            unsupported,
            {
                "accepted": False,
                "errors": ["The selected reference genome is unsupported."],
                "warnings": [],
                "checks": [],
            },
        )

    def test_warnings_do_not_hide_a_valid_modified_read(self) -> None:
        reads = [
            make_read(name="valid", start_0=100),
            make_read(
                name="incomplete",
                start_0=200,
                mm=CANONICAL_MM,
                ml=None,
            ),
            make_read(
                name="stale_mn",
                start_0=300,
                mn=len(CANONICAL_SEQUENCE) + 1,
            ),
            make_read(
                name="malformed",
                start_0=400,
                mm="C+m,0,0;",
                ml=(128,),
            ),
        ]
        bam_path, bai_path = write_bam(self.temporary_path, reads)
        previous_verbosity = pysam.set_verbosity(0)

        try:
            result = validate_modbam(bam_path, bai_path, "grch38")
        finally:
            pysam.set_verbosity(previous_verbosity)

        self.assertTrue(result["accepted"])
        self.assertEqual(len(result["warnings"]), 3)
        self.assertTrue(
            any("MM or ML, but not both" in warning for warning in result["warnings"])
        )
        self.assertTrue(
            any("MN value" in warning for warning in result["warnings"])
        )
        self.assertTrue(
            any("could not parse" in warning for warning in result["warnings"])
        )
        modified_base_check = next(
            check
            for check in result["checks"]
            if check["name"] == "Modified-base tags"
        )
        self.assertIn(
            "1 record with parsable MM/ML tags",
            modified_base_check["detail"],
        )
        self.assertIn("6 modified-base calls", modified_base_check["detail"])

    def test_empty_bam_missing_index_and_unsorted_header_fail(self) -> None:
        empty_bam, empty_bai = write_bam(
            self.temporary_path,
            [],
            name="empty",
        )
        empty_result = validate_modbam(empty_bam, empty_bai, "grch38")
        self.assertFalse(empty_result["accepted"])
        self.assertIn(
            "The BAM contains no alignment records.",
            empty_result["errors"],
        )

        no_index_bam, no_index_bai = write_bam(
            self.temporary_path,
            [make_read()],
            name="no_index",
            create_index=False,
        )
        no_index_result = validate_modbam(
            no_index_bam,
            no_index_bai,
            "grch38",
        )
        self.assertFalse(no_index_result["accepted"])
        self.assertTrue(
            any(
                "could not be validated" in error
                for error in no_index_result["errors"]
            )
        )

        unsorted_bam, unsorted_bai = write_bam(
            self.temporary_path,
            [make_read()],
            name="unsorted",
            sort_order="queryname",
        )
        unsorted_result = validate_modbam(
            unsorted_bam,
            unsorted_bai,
            "grch38",
        )
        self.assertFalse(unsorted_result["accepted"])
        self.assertIn(
            "The BAM header does not declare coordinate sorting.",
            unsorted_result["errors"],
        )


if __name__ == "__main__":
    unittest.main()
