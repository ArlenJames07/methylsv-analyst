from array import array
from pathlib import Path

import pysam


OUTPUT_DIR = Path("/tmp/methylsv-test")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

bam_path = OUTPUT_DIR / "valid_grch38_modbam.bam"

header = {
    "HD": {
        "VN": "1.6",
        "SO": "coordinate",
    },
    "SQ": [
        {
            "SN": "chr1",
            "LN": 248_956_422,
        }
    ],
}

sequence = "ACGTCGCGTACGTCGCGT"

read = pysam.AlignedSegment()
read.query_name = "synthetic_modbam_read"
read.query_sequence = sequence
read.flag = 0
read.reference_id = 0
read.reference_start = 100
read.mapping_quality = 60
read.cigartuples = [(0, len(sequence))]
read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))

# Six cytosines are annotated as 5-methylcytosine.
read.set_tag(
    "MM",
    "C+m,0,0,0,0,0,0;",
    value_type="Z",
)

read.set_tag(
    "ML",
    array("B", [230, 220, 210, 200, 190, 180]),
)

# Sequence length when MM/ML were generated.
read.set_tag(
    "MN",
    len(sequence),
    value_type="i",
)

with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
    bam.write(read)

pysam.index(str(bam_path))

with pysam.AlignmentFile(bam_path, "rb") as bam:
    test_read = next(bam.fetch("chr1"))

    print(f"BAM: {bam_path}")
    print(f"BAI: {bam_path}.bai")
    print(f"Coordinate sorted: {bam.header.to_dict()['HD']['SO']}")
    print(f"MM present: {test_read.has_tag('MM')}")
    print(f"ML present: {test_read.has_tag('ML')}")
    print(f"Modified bases parsed: {bool(test_read.modified_bases)}")