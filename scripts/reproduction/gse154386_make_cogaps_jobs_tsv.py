#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from gse154386_cogaps_sweep_common import parse_int_csv


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Create a headerless jobs.tsv for the GSE154386 CoGAPS Slurm array sweep."
    )
    parser.add_argument(
        "--k-grid",
        required=True,
        help="Comma-separated K grid, for example 6,8,10,12,14",
    )
    parser.add_argument(
        "--seeds",
        required=True,
        help="Comma-separated random seeds, for example 1,2,3,4,5",
    )
    parser.add_argument(
        "--iters",
        required=True,
        help="Comma-separated CoGAPS iteration counts, for example 4000,10000",
    )
    parser.add_argument(
        "--out",
        default=str(Path("GSE154386") / "cogaps_sweep_singleprocess_hpc" / "jobs.tsv"),
        help="Output TSV path. The file has no header and each line is K<TAB>seed<TAB>n_iter.",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()

    k_grid = parse_int_csv(args.k_grid)
    seeds = parse_int_csv(args.seeds)
    iters = parse_int_csv(args.iters)

    out_path = Path(args.out).resolve()
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    for k_value in k_grid:
        for n_iter in iters:
            for seed in seeds:
                lines.append(f"{k_value}\t{seed}\t{n_iter}\n")

    out_path.write_text("".join(lines), encoding="utf-8")
    print(f"[jobs] wrote {out_path} with {len(lines)} rows")


if __name__ == "__main__":
    main()
