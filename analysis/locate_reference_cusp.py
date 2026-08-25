#!/usr/bin/env python3
"""Locate the MAP critical endpoint as a cusp of the thermodynamic map."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pha_qnm.thermodynamics import BackgroundNumerics, locate_cusp_critical_point


ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phi-step", type=float, default=1.0e-3)
    parser.add_argument("--charge-step", type=float, default=1.0e-4)
    parser.add_argument("--uv-tolerance", type=float, default=1.0e-4)
    parser.add_argument("--output", type=Path,
                        default=ROOT / "results" / "python" / "reference_cusp.json")
    args = parser.parse_args()
    background_numerics = BackgroundNumerics(
        phiA_tolerance=args.uv_tolerance,
        ricci_tolerance=min(1.2e-3, args.uv_tolerance),
    )
    result = locate_cusp_critical_point(
        phi_step=args.phi_step,
        charge_step=args.charge_step,
        background_numerics=background_numerics,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
