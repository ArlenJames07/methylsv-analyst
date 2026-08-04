from pathlib import Path

import pysam

from tests.synthetic_modbam import make_read, write_bam


OUTPUT_DIRECTORY = Path("/tmp/methylsv-test")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    bam_path, bai_path = write_bam(
        OUTPUT_DIRECTORY,
        [make_read()],
        name="valid_grch38_modbam",
    )

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        test_read = next(bam.fetch("chr1"))

        print(f"BAM: {bam_path}")
        print(f"BAI: {bai_path}")
        print(
            "Coordinate sorted: "
            f"{bam.header.to_dict()['HD']['SO']}"
        )
        print(f"MM present: {test_read.has_tag('MM')}")
        print(f"ML present: {test_read.has_tag('ML')}")
        print(
            "Modified bases parsed: "
            f"{bool(test_read.modified_bases)}"
        )


if __name__ == "__main__":
    main()
