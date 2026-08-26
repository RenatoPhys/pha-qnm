"""Non-pixel figure contracts for the publication artifact pipeline."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis" / "figure_manifest.json"


@pytest.fixture(scope="module")
def entries():
    payload = json.loads(MANIFEST.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    return payload["figures"]


def test_manifest_names_and_outputs_are_unique(entries):
    names = [entry["figure"] for entry in entries]
    assert len(names) == len(set(names))
    outputs = [output for entry in entries for output in entry["outputs"]]
    assert len(outputs) == len(set(outputs))


@pytest.mark.parametrize("required", ["script", "inputs", "outputs", "axis_labels",
                                      "panel_labels", "semantic_layers", "caption_contract"])
def test_manifest_contract_fields(entries, required):
    assert all(entry.get(required) for entry in entries)


def test_inputs_and_outputs_exist(entries):
    missing = []
    for entry in entries:
        for relative in entry["inputs"] + entry["outputs"]:
            if not (ROOT / relative).exists():
                missing.append(relative)
    assert not missing, f"missing figure products: {missing}"


def test_outputs_are_nonempty_and_current(entries):
    stale = []
    for entry in entries:
        inputs = [ROOT / name for name in entry["inputs"]]
        newest_input = max(path.stat().st_mtime for path in inputs)
        for relative in entry["outputs"]:
            output = ROOT / relative
            assert output.stat().st_size > 1000
            if output.stat().st_mtime + 2.0 < newest_input:
                stale.append(relative)
    assert not stale, f"stale figure products: {stale}"


def test_pdfs_open_and_contain_vector_fonts(entries):
    pdfinfo = shutil.which("pdfinfo")
    assert pdfinfo, "Poppler pdfinfo is required for figure QA"
    for entry in entries:
        pdf = ROOT / entry["outputs"][0]
        result = subprocess.run([pdfinfo, str(pdf)], capture_output=True, text=True)
        assert result.returncode == 0, result.stderr
        assert "Pages:" in result.stdout and "Title:" in result.stdout
        content = pdf.read_bytes()
        assert b"/Font" in content, f"no vector font resources found in {pdf.name}"


def test_declared_pdf_and_png_stems_match(entries):
    for entry in entries:
        pdf, png = map(Path, entry["outputs"])
        assert pdf.suffix == ".pdf" and png.suffix == ".png"
        assert pdf.stem == png.stem == entry["figure"]
