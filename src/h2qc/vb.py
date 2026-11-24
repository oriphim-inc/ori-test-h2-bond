"""Valence-bond (ionic / covalent) decomposition of a two-electron singlet.

The alpha(1)beta(2) component of sum_ij c_ij |i alpha, j beta| has spatial part
sum_ij c_ij phi_i(1) phi_j(2) = sum_{mu nu} M_{mu nu} chi_mu(1) chi_nu(2),  M = C c C^T.
For H2 with AOs a, b:  Psi = M_ab (ab + ba) + M_aa aa + M_bb bb  (M_aa = M_bb by symmetry).
Structures: Phi_cov = (ab+ba)/N, Phi_ion = (aa+bb)/N, N^2 = 2(1 + S^2);
<Phi_cov|Phi_ion> = 2S/(1+S^2).
"""
from __future__ import annotations

import numpy as np


def vb_weights(C: np.ndarray, c: np.ndarray, S_ab: float):
    M = C @ c @ C.T
    M = 0.5 * (M + M.T)
    N = np.sqrt(2.0 * (1.0 + S_ab**2))
    # coefficients of normalised structures (ionic uses mean of aa, bb)
    x = np.array([M[0, 1] * N, 0.5 * (M[0, 0] + M[1, 1]) * N])  # [cov, ion]
    Sv = np.array([[1.0, 2 * S_ab / (1 + S_ab**2)], [2 * S_ab / (1 + S_ab**2), 1.0]])
    norm = x @ Sv @ x
    x = x / np.sqrt(norm)
    # Chirgwin-Coulson
    cc = x * (Sv @ x)
    # Loewdin (symmetric orthogonalisation of the two structures)
    s, U = np.linalg.eigh(Sv)
    Shalf = U @ np.diag(np.sqrt(s)) @ U.T
    lw = (Shalf @ x) ** 2
    return {
        "cov_CC": float(cc[0]), "ion_CC": float(cc[1]),
        "cov_L": float(lw[0]), "ion_L": float(lw[1]),
        "c_cov": float(x[0]), "c_ion": float(x[1]),
        "asym_aa_bb": float(abs(M[0, 0] - M[1, 1])),
        "norm_before": float(norm),
    }
