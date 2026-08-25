#!/usr/bin/env python3
"""Reproduce the sample74 critical point with the MUSES line-crossing route."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pha_qnm.thermodynamics import CriticalPointNumerics, locate_critical_point


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lines", type=int, default=800)
    parser.add_argument("--neighbors-only", action="store_true")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "python" / "reference_cep.json")
    args = parser.parse_args()
    result = locate_critical_point(
        options=CriticalPointNumerics(line_count=args.lines),
        progress=lambda message: print(message, flush=True),
        include_non_neighboring=not args.neighbors_only,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
