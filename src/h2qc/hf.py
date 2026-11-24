"""Closed-shell restricted Hartree-Fock (Roothaan-Hall SCF loop)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

from .integrals import Integrals


@dataclass
class RHFResult:
    E_total: float
    E_elec: float
    C: np.ndarray        # MO coefficients (AO x MO), S-orthonormal
    eps: np.ndarray      # orbital energies
    F: np.ndarray
    P: np.ndarray        # density matrix P = 2 C_occ C_occ^T
    n_iter: int
    converged: bool


def rhf(ints: Integrals, n_elec: int = 2, max_iter: int = 200,
        e_tol: float = 1e-12, d_tol: float = 1e-10) -> RHFResult:
    S, H, G = ints.S, ints.H, ints.G
    nocc = n_elec // 2
    # symmetric orthogonalisation X = S^-1/2
    s, U = np.linalg.eigh(S)
    X = U @ np.diag(s**-0.5) @ U.T

    def diag(F):
        e, Cp = np.linalg.eigh(X.T @ F @ X)
        return e, X @ Cp

    eps, C = diag(H)  # core-Hamiltonian guess
    P = 2.0 * C[:, :nocc] @ C[:, :nocc].T
    E_old = 0.0
    converged = False
    for it in range(1, max_iter + 1):
        J = np.einsum("mnls,ls->mn", G, P)
        K = np.einsum("mlns,ls->mn", G, P)
        F = H + J - 0.5 * K
        E_elec = 0.5 * np.sum(P * (H + F))
        eps, C = diag(F)
        P_new = 2.0 * C[:, :nocc] @ C[:, :nocc].T
        dP = np.sqrt(np.mean((P_new - P) ** 2))
        dE = abs(E_elec - E_old)
        P, E_old = P_new, E_elec
        if dE < e_tol and dP < d_tol:
            converged = True
            break
    # final energy and Fock with converged density
    J = np.einsum("mnls,ls->mn", G, P)
    K = np.einsum("mlns,ls->mn", G, P)
    F = H + J - 0.5 * K
    E_elec = 0.5 * np.sum(P * (H + F))
    eps, C = diag(F)
    # fix MO sign convention: largest-|c| AO coefficient positive on atom A
    for k in range(C.shape[1]):
        if C[0, k] < 0:
            C[:, k] *= -1
    return RHFResult(E_elec + ints.E_nuc, E_elec, C, eps, F, P, it, converged)
