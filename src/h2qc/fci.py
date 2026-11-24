"""Full CI for 2 electrons in the spin orbitals built from the RHF MOs.

Determinants |i alpha, j beta| (M_s = 0 block, n_mo^2 determinants) — exact for the
singlet ground state because H conserves M_s. For two-electron determinants
|ab| = (ab - ba)/sqrt(2) with orthonormal spin orbitals:
  <ab|H|cd> = h_ac d_bd + d_ac h_bd - h_ad d_bc - d_ad h_bc + <ab|cd> - <ab|dc>
with physicists' <pq|rs> = (pr|qs) * spin deltas.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .integrals import Integrals


@dataclass
class FCIResult:
    E_total: float
    coeff: np.ndarray        # (n_mo, n_mo): c[i, j] for |i alpha, j beta|
    energies: np.ndarray     # all M_s=0 eigenvalues (total energies)
    singlet_flags: np.ndarray
    E_triplet: float         # lowest triplet (M_s=0 component) total energy


def mo_integrals(ints: Integrals, C: np.ndarray):
    h = C.T @ ints.H @ C
    g = np.einsum("mnls,mp,nq,lr,st->pqrt", ints.G, C, C, C, C, optimize=True)
    return h, g


def fci(ints: Integrals, C: np.ndarray) -> FCIResult:
    h, g = mo_integrals(ints, C)
    n = h.shape[0]
    # spin orbitals: (spatial index, spin) with spin 0=alpha, 1=beta
    def h1(p, q):
        return h[p[0], q[0]] if p[1] == q[1] else 0.0

    def ov(p, q):
        return 1.0 if p == q else 0.0

    def g2(p, q, r, s):  # <pq|rs> = (pr|qs)
        if p[1] != r[1] or q[1] != s[1]:
            return 0.0
        return g[p[0], r[0], q[0], s[0]]

    dets = [((i, 0), (j, 1)) for i in range(n) for j in range(n)]
    N = len(dets)
    Hm = np.zeros((N, N))
    for I, (a, b) in enumerate(dets):
        for J, (c, d) in enumerate(dets):
            Hm[I, J] = (h1(a, c) * ov(b, d) + ov(a, c) * h1(b, d)
                        - h1(a, d) * ov(b, c) - ov(a, d) * h1(b, c)
                        + g2(a, b, c, d) - g2(a, b, d, c))
    Hm += ints.E_nuc * np.eye(N)
    assert np.allclose(Hm, Hm.T)
    # Spin adaptation: in the M_s=0 block a singlet has symmetric c_ij, a triplet
    # antisymmetric c_ij. Build orthonormal bases of the two subspaces and diagonalise
    # H in each, so near-degenerate singlet/triplet roots (large R) cannot mix.
    sym, anti = [], []
    for i in range(n):
        for j in range(i, n):
            v = np.zeros((n, n))
            if i == j:
                v[i, i] = 1.0
                sym.append(v.ravel())
            else:
                v[i, j] = v[j, i] = 2**-0.5
                sym.append(v.ravel())
                w = np.zeros((n, n))
                w[i, j], w[j, i] = 2**-0.5, -(2**-0.5)
                anti.append(w.ravel())
    Bs, Bt = np.array(sym).T, np.array(anti).T
    Es, Vs = np.linalg.eigh(Bs.T @ Hm @ Bs)
    Et, _ = np.linalg.eigh(Bt.T @ Hm @ Bt)
    E = np.concatenate([Es, Et])
    flags = np.array([True] * len(Es) + [False] * len(Et))
    order = np.argsort(E)
    E, flags = E[order], flags[order]
    c0 = (Bs @ Vs[:, 0]).reshape(n, n)
    if c0[0, 0] < 0:
        c0 = -c0
    return FCIResult(float(Es[0]), c0, E, flags, float(Et.min()))
