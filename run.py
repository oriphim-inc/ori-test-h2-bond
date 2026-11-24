"""Reproduce every result: python run.py [--n N]

Writes results/scan.csv, results/summary.json, results/pes.png, results/ionic_weights.png,
results/reference_point.json (integrals and energies at R = 1.4 bohr).
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "src"))
os.environ.setdefault("MPLCONFIGDIR", str(HERE / ".mplcache"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from h2qc.fci import fci  # noqa: E402
from h2qc.hf import rhf  # noqa: E402
from h2qc.integrals import (BOHR_TO_ANG, HARTREE_TO_EV, HARTREE_TO_KCAL,  # noqa: E402
                            h2_integrals)
from h2qc.scan import run_scan  # noqa: E402


def reference_point(R_bohr=1.4):
    ints = h2_integrals(R_bohr)
    hf = rhf(ints)
    ci = fci(ints, hf.C)
    G = ints.G
    return {
        "R_bohr": R_bohr, "S12": ints.S[0, 1], "T11": ints.T[0, 0], "T12": ints.T[0, 1],
        "V11_total": ints.V[0, 0], "V12_total": ints.V[0, 1],
        "(11|11)": G[0, 0, 0, 0], "(11|22)": G[0, 0, 1, 1],
        "(21|11)": G[1, 0, 0, 0], "(21|21)": G[1, 0, 1, 0],
        "E_RHF": hf.E_total, "E_FCI": ci.E_total, "eps": hf.eps.tolist(),
        "scf_iter": hf.n_iter,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=461, help="number of R points (0.4..5.0 Å)")
    ap.add_argument("--out", default="results")
    args = ap.parse_args()
    out = HERE / args.out
    out.mkdir(exist_ok=True)
    t0 = time.time()

    ref = {k: (float(v) if not isinstance(v, list) else v) for k, v in reference_point().items()}
    rows, summary = run_scan(args.n)
    summary["units"] = {"energy": "hartree", "R": "angstrom (and bohr where noted)",
                        "bohr_to_ang": BOHR_TO_ANG}
    for tag in ("RHF", "FCI"):
        s = summary[tag]
        s["D_e_vs_2H_eV"] = s["D_e_vs_2H_Ha"] * HARTREE_TO_EV
        s["D_e_vs_2H_kcal"] = s["D_e_vs_2H_Ha"] * HARTREE_TO_KCAL
        s["D_e_vs_own_end_eV"] = s["D_e_vs_own_end_Ha"] * HARTREE_TO_EV
    last = rows[-1]
    summary["at_Rmax"] = {k: last[k] for k in ("R_ang", "E_RHF", "E_FCI", "E_triplet",
                                                "FCI_ion_CC", "FCI_ion_L", "RHF_ion_CC",
                                                "RHF_ion_L", "c_gg", "c_uu")}
    summary["all_scf_converged"] = all(r["scf_converged"] for r in rows)
    summary["max_scf_iter"] = max(r["scf_iter"] for r in rows)
    summary["runtime_s"] = time.time() - t0

    with open(out / "scan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    (out / "reference_point.json").write_text(json.dumps(ref, indent=2))

    R = np.array([r["R_ang"] for r in rows])
    Ehf = np.array([r["E_RHF"] for r in rows])
    Eci = np.array([r["E_FCI"] for r in rows])
    E2H = summary["E_2H"]

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(R, Ehf, label="RHF")
    ax.plot(R, Eci, label="FCI (2e, 2 orbitals)")
    ax.axhline(E2H, color="k", ls="--", lw=0.8, label=f"2 E_H(STO-3G) = {E2H:.4f}")
    for tag, c in (("RHF", "C0"), ("FCI", "C1")):
        s = summary[tag]
        ax.plot(s["R_e_ang"], s["E_min"], "o", color=c)
        ax.annotate(f"R_e={s['R_e_ang']:.3f} Å\nE={s['E_min']:.4f}",
                    (s["R_e_ang"], s["E_min"]), textcoords="offset points",
                    xytext=(15, -25 if tag == "FCI" else 10), color=c, fontsize=8)
    ax.set_xlabel("R (Å)")
    ax.set_ylabel("E (hartree)")
    ax.set_ylim(min(Eci.min(), Ehf.min()) - 0.05, max(0.1, Ehf.max()) + 0.02)
    ax.set_title("H$_2$ / STO-3G potential energy curves")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "pes.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.8))
    ax.plot(R, [r["FCI_ion_CC"] for r in rows], "C1-", label="FCI ionic (Chirgwin–Coulson)")
    ax.plot(R, [r["FCI_ion_L"] for r in rows], "C1--", label="FCI ionic (Löwdin)")
    ax.plot(R, [r["RHF_ion_CC"] for r in rows], "C0-", label="RHF ionic (Chirgwin–Coulson)")
    ax.plot(R, [r["RHF_ion_L"] for r in rows], "C0--", label="RHF ionic (Löwdin)")
    ax.plot(R, [r["c_uu"] ** 2 for r in rows], "C2:", label="FCI |c(σu²)|² (MO picture)")
    ax.set_xlabel("R (Å)")
    ax.set_ylabel("weight")
    ax.set_ylim(-0.02, 1.0)
    ax.set_title("Ionic weight of the ground state (covalent = 1 − ionic)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "ionic_weights.png", dpi=150)
    plt.close(fig)

    print(json.dumps({"reference_point": ref, "summary": summary}, indent=2))


if __name__ == "__main__":
    main()
