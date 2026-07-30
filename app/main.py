from fastapi import FastAPI
from fastapi.responses import HTMLResponse

app = FastAPI(
    title="MethylSV Analyst",
    description="Platform for structural variant and methylation analysis",
    version="0.1.0",
)


@app.get("/", response_class=HTMLResponse)
def home():
    return """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>MethylSV Analyst</title>
    </head>
    <body>
        <h1>MethylSV Analyst</h1>
        <p>Structural variation and DNA methylation analysis platform.</p>
        <p>The FastAPI server is working correctly.</p>
    </body>
    </html>
    """


@app.get("/health")
def health_check():
    return {
        "status": "healthy",
        "application": "MethylSV Analyst",
        "version": "0.1.0",
    }