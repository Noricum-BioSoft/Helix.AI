#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import random
from pathlib import Path


def _poisson(lam: float, rng: random.Random) -> int:
    # Knuth for small lam; normal approx for large lam
    if lam <= 0:
        return 0
    if lam < 30:
        L = math.exp(-lam)
        k = 0
        p = 1.0
        while p > L:
            k += 1
            p *= rng.random()
        return max(0, k - 1)
    # normal approximation
    x = rng.gauss(lam, math.sqrt(lam))
    return max(0, int(round(x)))


def generate_demo(out_dir: Path, seed: int = 7, n_genes: int = 250) -> tuple[Path, Path]:
    rng = random.Random(seed)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 2×2 factorial, 3 replicates each
    groups = [
        ("uninfected", "11dpi"),
        ("infected", "11dpi"),
        ("uninfected", "33dpi"),
        ("infected", "33dpi"),
    ]
    samples: list[dict] = []
    for infection, time in groups:
        for rep in range(1, 4):
            sample_id = f"{infection[:2]}_{time}_{rep}"
            samples.append(
                {
                    "sample_id": sample_id,
                    "infection_status": infection,
                    "time_point": time,
                }
            )

    # Create effects for a subset of genes
    genes = [f"Gene{idx:04d}" for idx in range(1, n_genes + 1)]
    infection_genes = set(rng.sample(genes, k=max(10, n_genes // 25)))
    time_genes = set(rng.sample([g for g in genes if g not in infection_genes], k=max(10, n_genes // 25)))
    interaction_genes = set(
        rng.sample([g for g in genes if g not in infection_genes and g not in time_genes], k=max(8, n_genes // 30))
    )

    # Baseline expression per gene
    base_means = {g: 20 * math.exp(rng.uniform(0.0, 4.0)) for g in genes}  # ~20..~1100

    # Simple library size factors
    lib_factors = {s["sample_id"]: rng.uniform(0.85, 1.20) for s in samples}

    counts_path = out_dir / "tgondii_counts.csv"
    meta_path = out_dir / "tgondii_metadata.csv"

    # Write metadata
    with meta_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "infection_status", "time_point"])
        w.writeheader()
        for s in samples:
            w.writerow(s)

    # Write counts
    sample_ids = [s["sample_id"] for s in samples]
    with counts_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["gene_id", *sample_ids])
        for g in genes:
            row = [g]
            mu0 = base_means[g]
            for s in samples:
                sid = s["sample_id"]
                mu = mu0 * lib_factors[sid]

                # main effects
                if g in infection_genes and s["infection_status"] == "infected":
                    mu *= rng.uniform(1.4, 2.2)
                if g in time_genes and s["time_point"] == "33dpi":
                    mu *= rng.uniform(1.3, 2.0)

                # interaction: extra boost only in infected@33dpi
                if (
                    g in interaction_genes
                    and s["infection_status"] == "infected"
                    and s["time_point"] == "33dpi"
                ):
                    mu *= rng.uniform(1.6, 3.0)

                # add mild overdispersion by jittering mean
                mu *= rng.uniform(0.85, 1.15)
                row.append(_poisson(mu, rng))
            w.writerow(row)

    return counts_path, meta_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="Output directory for CSVs")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--n-genes", type=int, default=250)
    args = ap.parse_args()

    out_dir = Path(os.path.expanduser(args.out_dir)).resolve()
    counts_path, meta_path = generate_demo(out_dir, seed=args.seed, n_genes=args.n_genes)
    print(str(counts_path))
    print(str(meta_path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

