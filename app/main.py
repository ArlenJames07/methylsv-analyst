from pathlib import Path

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from app.routes.modbam import router as modbam_router


BASE_DIR = Path(__file__).resolve().parent

app = FastAPI(
    title="MethylSV Analyst",
    description="Platform for methylation and structural variant analysis",
    version="0.2.0",
)

app.include_router(modbam_router)

app.mount(
    "/static",
    StaticFiles(directory=BASE_DIR / "static"),
    name="static",
)

templates = Jinja2Templates(directory=BASE_DIR / "templates")


WORKFLOWS = {
    "modbam": {
        "id": "modbam",
        "number": "01",
        "title": "Methylation visualization",
        "subtitle": "Single-sample modBAM analysis",
        "description": (
            "Inspect DNA methylation at read level and generate regional "
            "plots using modbamtools."
        ),
        "inputs": [
            "Exactly one coordinate-sorted modBAM file",
            "Matching BAM index: .bai",
            "Reference: GRCh38/hg38 or T2T-CHM13v2.0",
            "Optional indexed VCF for haplotagging",
        ],
        "steps": [
            "Validate the BAM, index, reference and MM/ML tags",
            "Optionally haplotag the BAM with HiPhase",
            "Select a genomic region",
            "Generate read-level methylation visualization",
            "Download plots, logs and the reproducibility manifest",
        ],
    },
    "dss": {
        "id": "dss",
        "number": "02",
        "title": "Differential methylation",
        "subtitle": "Two-group DSS analysis",
        "description": (
            "Identify differentially methylated loci and regions using "
            "per-sample methylated and total-read counts."
        ),
        "inputs": [
            "One metadata CSV file",
            "One DSS-format TSV file per sample",
            "Required DSS columns: chr, pos, N, X",
            "Reference: GRCh38/hg38 or T2T-CHM13v2.0",
        ],
        "steps": [
            "Validate metadata and sample-file correspondence",
            "Validate coverage and methylated-read counts",
            "Perform sample QC and exploratory analysis",
            "Run DSS differential methylation testing",
            "Call, annotate and export DMRs",
        ],
    },
}


@app.get("/", response_class=HTMLResponse, name="home")
def home(request: Request):
    return templates.TemplateResponse(
        request=request,
        name="index.html",
        context={"workflows": list(WORKFLOWS.values())},
    )


@app.get(
    "/workflows/{workflow_id}",
    response_class=HTMLResponse,
    name="workflow_page",
)
def workflow_page(request: Request, workflow_id: str):
    workflow = WORKFLOWS.get(workflow_id)

    if workflow is None:
        raise HTTPException(status_code=404, detail="Unknown workflow")

    return templates.TemplateResponse(
        request=request,
        name="workflow.html",
        context={"workflow": workflow},
    )


@app.get("/health")
def health_check():
    return {
        "status": "healthy",
        "application": "MethylSV Analyst",
        "version": "0.2.0",
    }