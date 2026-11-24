"""Bond-length scan, equilibrium search and summary for H2 / STO-3G."""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar

from .fci import fci
from .hf import rhf
from .integrals import ANG_TO_BOHR, h2_integrals, h_atom_energy
from .vb import vb_weights


def point(R_ang: float) -> dict:
    ints = h2_integrals(R_ang * ANG_TO_BOHR)
    hf = rhf(ints)
    ci = fci(ints, hf.C)
    S_ab = float(ints.S[0, 1])
    c_hf = np.zeros((2, 2))
    c_hf[0, 0] = 1.0
    w_ci = vb_weights(hf.C, ci.coeff, S_ab)
    w_hf = vb_weights(hf.C, c_hf, S_ab)
    return {
        "R_ang": R_ang, "R_bohr": R_ang * ANG_TO_BOHR, "S_ab": S_ab,
        "E_RHF": hf.E_total, "E_FCI": ci.E_total, "E_triplet": ci.E_triplet,
        "scf_iter": hf.n_iter, "scf_converged": hf.converged,
        "eps_g": float(hf.eps[0]), "eps_u": float(hf.eps[1]),
        "c_gg": float(ci.coeff[0, 0]), "c_uu": float(ci.coeff[1, 1]),
        "c_gu": float(ci.coeff[0, 1]), "c_ug": float(ci.coeff[1, 0]),
        "FCI_ion_CC": w_ci["ion_CC"], "FCI_cov_CC": w_ci["cov_CC"],
        "FCI_ion_L": w_ci["ion_L"], "FCI_cov_L": w_ci["cov_L"],
        "RHF_ion_CC": w_hf["ion_CC"], "RHF_cov_CC": w_hf["cov_CC"],
        "RHF_ion_L": w_hf["ion_L"], "RHF_cov_L": w_hf["cov_L"],
    }


def find_min(key: str, R_grid: np.ndarray, E_grid: np.ndarray):
    k = int(np.argmin(E_grid))
    if k == 0 or k == len(R_grid) - 1:
        raise RuntimeError(f"{key}: minimum at grid edge")
    lo, hi = R_grid[k - 1], R_grid[k + 1]
    res = minimize_scalar(lambda r: point(r)[key], bounds=(lo, hi), method="bounded",
                          options={"xatol": 1e-7})
    return float(res.x), float(res.fun)


def run_scan(n: int = 461, R_min: float = 0.4, R_max: float = 5.0):
    R = np.linspace(R_min, R_max, n)
    rows = [point(float(r)) for r in R]
    E_H = h_atom_energy()
    summary = {"E_H_atom_computed": E_H, "E_H_atom_stated": -0.4666,
               "E_2H": 2 * E_H, "n_points": n, "R_range_ang": [R_min, R_max]}
    for key, tag in (("E_RHF", "RHF"), ("E_FCI", "FCI")):
        E = np.array([r[key] for r in rows])
        Re, Emin = find_min(key, R, E)
        E_end = float(E[-1])
        summary[tag] = {
            "R_e_ang": Re, "R_e_bohr": Re * ANG_TO_BOHR, "E_min": Emin,
            "D_e_vs_2H_Ha": 2 * E_H - Emin,
            "D_e_vs_own_end_Ha": E_end - Emin,
            "E_at_Rmax": E_end,
            "dissociation_error_Ha": E_end - 2 * E_H,
            "grid_min_R_ang": float(R[np.argmin(E)]),
        }
    return rows, summary
