# MethylSV Analyst

MethylSV Analyst is a FastAPI research workspace for validating one modBAM/BAI pair and inspecting read-level modified-base calls in a genomic interval.

Version 0.2.0 is a development prototype for research use. The modBAM workflow is functional; the DSS workflow is an informational roadmap and cannot execute analyses yet.

## Requirements

- Python 3.10 or newer
- A POSIX-like environment with writable `/tmp` storage
- A coordinate-sorted modBAM and its matching BAI for real analyses

Create a local environment and install the pinned runtime dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

## Run locally

Start the development server from the repository root:

```bash
uvicorn app.main:app --reload
```

Open <http://127.0.0.1:8000>. The health endpoint is available at <http://127.0.0.1:8000/health>.

`--reload` is intended only for local development. A production server, container image, and deployment configuration are not provided yet.

Keep the server bound to `127.0.0.1`. Do not expose this prototype to an untrusted network: it has no authentication or upload-size limit.

## Implemented workflow

The single-sample modBAM workflow supports:

1. Uploading exactly one `.bam` and one matching `.bai`.
2. Selecting GRCh38/hg38 or T2T-CHM13v2.0.
3. Validating the BAM, index, sort order, reference fingerprint, and MM/ML annotations.
4. Selecting a genomic region with 1-based, inclusive coordinates.
5. Reviewing overlapping reads and decoded modified-base calls.

Results include the genomic call position, modification name, ML-derived probability midpoint, alignment and modification strands, mapping quality, and HP/PS tags when present. The service layer also retains the probability bucket bounds.

Current processing limits:

- Maximum region length: 1,000,000 bp
- Maximum overlapping alignments processed: 10,000
- Maximum validation scan: first 10,000 alignment records
- Maximum reads and calls displayed in each table: 200

The region query includes primary, secondary, supplementary, duplicate, and QC-fail alignments. These categories are not filtered or deduplicated in version 0.2.0.

## Validation behavior

The upload and service layers check:

- `.bam` and `.bai` file extensions and binary signatures
- BAM readability and BAI usability through pysam
- `SO:coordinate` in the BAM header
- The selected reference through a recognized chromosome 1 alias and its exact length
- At least one mapped record with paired, parsable MM/ML tags and modified-base calls among the scanned records

The reference check is a chromosome 1 fingerprint, not full reference-sequence verification. Warnings identify incomplete MM/ML pairs, stale MN lengths, and annotations that pysam cannot decode.

## Automated tests

Run the standard-library regression suite without installing extra test dependencies:

```bash
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover \
  -s tests -p 'test_*.py' -v
```

The suite creates isolated temporary BAMs and covers:

- ML probability boundaries and invalid values
- Forward and reverse modified-base orientation
- Soft clipping and calls without reference positions
- Malformed MM/ML annotations and stale MN lengths
- 0-based/half-open to 1-based/inclusive coordinate conversion
- Region aggregation, decoding errors, and result truncation
- Reference mismatches, warning aggregation, empty BAMs, missing indexes, and unsorted headers

These tests cover the scientific services. HTTP/multipart routes, job lifecycle and permissions, templates, JavaScript, and browser behavior do not yet have automated integration tests.

Generate a reusable manual fixture with:

```bash
python -m tests.create_test_modbam
```

This writes a synthetic GRCh38 modBAM and BAI under `/tmp/methylsv-test`. Select GRCh38 in the UI and query `chr1:101-118`; the expected result is one read with six 5mC calls.

## Data handling and safety

Use only synthetic or properly deidentified data in this development environment. Deidentification must include uploaded filenames because the original BAM and BAI names are stored in job metadata and displayed in the UI.

- Accepted uploads are stored under `/tmp/methylsv-jobs/<random-job-id>`.
- Job directories use mode `0700`; staged files and metadata use mode `0600`.
- Rejected staged uploads are removed.
- Jobs become eligible for cleanup after 60 minutes.
- Cleanup is opportunistic when jobs are created or loaded; it is not a scheduled deletion guarantee.
- The job ID is embedded in the region URL and acts as a bearer secret while the job exists. Anyone with that URL can access the analysis, and normal server access logs can contain the ID.

The application currently has no authentication, TLS termination, encryption-at-rest layer, audit trail, database, configured durable datastore, or production privacy controls. Uploaded files remain on local disk until opportunistic cleanup or operating-system removal. It is not intended for clinical decision-making.

## Known limitations

- DSS differential methylation is planned but not implemented.
- Structural-variant analysis is not implemented.
- There are no plots, downloads, exports, reports, or reproducibility manifests.
- HiPhase/haplotagging and VCF input are not supported.
- Upload size and decoded-call memory limits are not configured.
- Job cleanup runs only in response to application activity.
- Temporary storage uses an absolute path under `/tmp`.
- High-depth regions can accumulate more decoded calls in memory than the 200 rows shown in the UI.
- Container deployment is not currently supported.

## Project layout

```text
app/main.py                       Application and workflow catalog
app/routes/modbam.py              Upload, validation, and region routes
app/services/modbam_validator.py  BAM/index/reference/MM-ML validation
app/services/modification_calls.py Modified-base decoding
app/services/region_reader.py     Regional read and call aggregation
app/templates/                    Server-rendered pages
app/static/                       Styles, JavaScript, and favicon
tests/                            Synthetic fixtures and regression tests
```
