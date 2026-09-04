#!/usr/bin/env python3
"""Kure boslugu benchmark'ini bir dizi konfigurasyonda kosturur,
results/sphere_scf.csv ve figures/sphere_convergence.png uretir."""

import csv
import pathlib
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN = ROOT / "build" / "cpp" / "fcm3d_sphere"
CSV_PATH = ROOT / "results" / "sphere_scf.csv"
FIG_PATH = ROOT / "figures" / "sphere_convergence.png"

# (p, depth, n_elements, off, R)
CONFIGS = [
    (3, 3, 2, 1.10, 0.15),
    (4, 3, 2, 1.10, 0.15),
    (5, 3, 2, 1.10, 0.15),
    (6, 3, 2, 1.10, 0.15),
    (4, 3, 3, 1.10, 0.15),   # ekvator hucre icinde
    (8, 3, 1, 1.10, 0.15),   # hic arayuz yok
]

FIELDS = ["p", "depth", "n_elements", "off", "R", "dof", "quad_points",
          "scf_avg", "scf_max", "scf_min", "scf_analytic"]


def run(cfg):
    args = [str(BIN)] + [str(v) for v in cfg]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    for line in out.splitlines():
        if line.startswith("CSV,"):
            return dict(zip(FIELDS, line.split(",")[1:]))
    raise RuntimeError(f"CSV satiri yok: {args}")


def main():
    if not BIN.exists():
        sys.exit(f"ikili bulunamadi: {BIN}\nonce derle: cmake --build build")

    rows = []
    for cfg in CONFIGS:
        print("kosuyor", cfg, flush=True)
        rows.append(run(cfg))

    CSV_PATH.parent.mkdir(exist_ok=True)
    with CSV_PATH.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print("yazildi:", CSV_PATH)

    sweep = sorted((r for r in rows if int(r["n_elements"]) == 2),
                   key=lambda r: int(r["p"]))
    dof = [int(r["dof"]) for r in sweep]
    avg = [float(r["scf_avg"]) for r in sweep]
    lo  = [float(r["scf_min"]) for r in sweep]
    hi  = [float(r["scf_max"]) for r in sweep]
    exact = float(sweep[0]["scf_analytic"])

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.axhline(exact, color="k", ls="--", lw=1, label=f"analytic  {exact:.4f}")
    ax.fill_between(dof, lo, hi, alpha=0.25, color="C0",
                    label="angular spread (min-max)")
    ax.plot(dof, avg, "o-", color="C0", label="circumferential mean")
    for r, d, a in zip(sweep, dof, avg):
        ax.annotate(f"p={r['p']}", (d, a), textcoords="offset points",
                    xytext=(6, -12), fontsize=9)
    ax.set_xlabel("degrees of freedom")
    ax.set_ylabel(r"$\sigma_{zz}\,/\,\sigma_\infty$  at  $1.10\,R$")
    ax.set_title("Spherical void in uniaxial tension\n"
                 r"$2^3$ cells, octree depth 3, $\nu=0.3$")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    FIG_PATH.parent.mkdir(exist_ok=True)
    fig.savefig(FIG_PATH, dpi=150)
    print("yazildi:", FIG_PATH)


if __name__ == "__main__":
    main()