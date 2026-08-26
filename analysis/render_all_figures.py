#!/usr/bin/env python3
"""Single entry point for publication figures, gallery, and contract checks."""

from __future__ import annotations

import argparse
import html
import json
from pathlib import Path
import subprocess
import sys

from PIL import Image, ImageDraw, ImageFont, ImageOps


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "paper" / "figures"
GALLERY = RESULTS / "figure_gallery"
MANIFEST = Path(__file__).with_name("figure_manifest.json")

MAIN_COMMANDS = (
    ("build_phase_diagram.py", "--plot-only"),
    ("plot_background_grid_map.py",),
    ("plot_finite_k_figures.py",),
    ("plot_posterior_uq.py",),
)
AUDIT_COMMAND = (("audit_longitudinal_mode_identity.py", "--hydro-only"),)
APPENDIX_COMMANDS = (
    ("plot_homogeneous_trajectories.py",),
    ("plot_legacy_2018_curves.py",),
    ("plot_cep_reproduction.py",),
    ("plot_qnm_validation.py",),
)


def run_commands(commands) -> None:
    analysis = Path(__file__).parent
    for command in commands:
        subprocess.run([sys.executable, str(analysis / command[0]), *command[1:]],
                       cwd=ROOT, check=True)


def make_gallery() -> None:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    GALLERY.mkdir(parents=True, exist_ok=True)
    tiles = []
    cards = []
    for entry in manifest["figures"]:
        name = entry["figure"]
        source = ROOT / entry["outputs"][1]
        if not source.exists():
            continue
        image = Image.open(source).convert("RGB")
        normal = GALLERY / f"{name}.png"
        image.save(normal)
        thumb = ImageOps.contain(image.copy(), (520, 360))
        thumb_path = GALLERY / f"{name}_thumbnail.png"
        thumb.save(thumb_path)
        gray = ImageOps.grayscale(image)
        gray_path = GALLERY / f"{name}_grayscale.png"
        gray.save(gray_path)
        canvas = Image.new("RGB", (560, 420), "white")
        canvas.paste(thumb, ((560-thumb.width)//2, 36))
        draw = ImageDraw.Draw(canvas)
        draw.text((14, 10), name, fill="#111111", font=ImageFont.load_default())
        tiles.append(canvas)
        cards.append(
            f'<section><h2>{html.escape(name)}</h2>'
            f'<a href="{name}.png"><img src="{name}_thumbnail.png" alt="{name}"></a>'
            f'<img src="{name}_grayscale.png" class="gray" alt="{name} grayscale"></section>'
        )
    if tiles:
        columns = 2; rows = (len(tiles)+1)//2
        sheet = Image.new("RGB", (columns*560, rows*420), "#eef0f2")
        for index, tile in enumerate(tiles):
            sheet.paste(tile, ((index % columns)*560, (index // columns)*420))
        sheet.save(GALLERY / "contact_sheet.png")
    (GALLERY / "index.html").write_text(
        "<!doctype html><meta charset='utf-8'><title>PHA-QNM figure QA</title>"
        "<style>body{font-family:system-ui;margin:2rem;max-width:1200px}"
        "section{border-top:1px solid #ccc;padding:1rem 0}img{max-width:48%;vertical-align:top;margin-right:1%}"
        ".gray{filter:none}</style><h1>PHA-QNM figure QA gallery</h1>" + "".join(cards),
        encoding="utf-8")


def check_contracts() -> None:
    subprocess.run([sys.executable, "-m", "pytest", "-q", "tests/test_figure_contracts.py"],
                   cwd=ROOT, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--main-only", action="store_true")
    group.add_argument("--appendix-only", action="store_true")
    parser.add_argument("--gallery", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.main_only:
        run_commands(MAIN_COMMANDS)
    elif args.appendix_only:
        run_commands(AUDIT_COMMAND + APPENDIX_COMMANDS)
    else:
        run_commands(AUDIT_COMMAND + MAIN_COMMANDS + APPENDIX_COMMANDS)
    if args.gallery or not (args.main_only or args.appendix_only):
        make_gallery()
    if args.check:
        check_contracts()


if __name__ == "__main__":
    main()
