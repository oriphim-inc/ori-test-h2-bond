"""Independent reference for the Sceptic's checks (does NOT import h2qc).

Written separately from src/h2qc: vectorised primitive integrals (Szabo & Ostlund App. A)
and the closed-form, symmetry-adapted minimal-basis H2 solution. In a 2-AO basis with
D_inf_h symmetry the MOs are fixed by symmetry,
    g = (a+b)/sqrt(2(S_aa+S)),  u = (a-b)/sqrt(2(S_aa-S)),
(normalised with the actual S_aa, which is 1 - 9.1e-9 for the 8-digit constants, so
the MOs are exactly S-orthonormal), so E_RHF = 2 h_gg + (gg|gg) + 1/R needs no SCF, and the singlet FCI is the 2x2 problem
    [[2h_gg+(gg|gg), (gu|gu)], [(gu|gu), 2h_uu+(uu|uu)]] + 1/R,
the triplet is h_gg + h_uu + (gg|uu) - (gu|gu) + 1/R.
"""
from __future__ import annotations

import numpy as np
from scipy.special import erf

BOHR_ANG = 0.529177210903
# STO-3G H (zeta = 1.24), PLAN.md / Hehre, Stewart & Pople 1969
ALPHA = np.array([3.42525091, 0.62391373, 0.16885540])
D = np.array([0.15432897, 0.53532814, 0.44463454])
CN = D * (2 * ALPHA / np.pi) ** 0.75


def _F0(t):
    t = np.asarray(t, float)
    out = np.ones_like(t)
    m = t > 1e-12
    st = np.sqrt(t[m])
    out[m] = 0.5 * np.sqrt(np.pi) * erf(st) / st
    out[~m] = 1 - t[~m] / 3
    return out


def ao_integrals(R):
    """Return S_ab, h (2x2 core), and the six unique ERIs over AOs a (at 0) and b (at R)."""
    X = np.array([0.0, R])  # 1D positions along bond
    a = ALPHA
    # pair data over (centre i, prim k, centre j, prim l)
    A = X[:, None, None, None]
    B = X[None, None, :, None]
    al = a[None, :, None, None]
    be = a[None, None, None, :]
    c = CN[None, :, None, None] * CN[None, None, None, :]
    p = al + be
    mu = al * be / p
    Rab2 = (A - B) ** 2
    K = np.exp(-mu * Rab2)
    P = (al * A + be * B) / p
    Sp = (np.pi / p) ** 1.5 * K
    Tp = mu * (3 - 2 * mu * Rab2) * Sp
    Vp = np.zeros_like(Sp)
    for C in X:
        Vp += -2 * np.pi / p * K * _F0(p * (P - C) ** 2)
    shp = (2, 3, 2, 3)
    p, P, K, c = (np.broadcast_to(z, shp) for z in (p, P, K, c))
    S = np.einsum("ikjl->ij", c * Sp)
    T = np.einsum("ikjl->ij", c * Tp)
    V = np.einsum("ikjl->ij", c * Vp)
    # ERIs: (ij|kl) over contracted AOs, brute force over primitive quartets
    G = np.zeros((2, 2, 2, 2))
    for i in range(2):
        for j in range(2):
            for k in range(2):
                for l in range(2):
                    p1 = p[i, :, j, :][:, :, None, None]
                    q1 = p[k, :, l, :][None, None, :, :]
                    P1 = P[i, :, j, :][:, :, None, None]
                    Q1 = P[k, :, l, :][None, None, :, :]
                    K1 = (c * K)[i, :, j, :][:, :, None, None]
                    K2 = (c * K)[k, :, l, :][None, None, :, :]
                    pref = 2 * np.pi ** 2.5 / (p1 * q1 * np.sqrt(p1 + q1))
                    G[i, j, k, l] = np.sum(
                        K1 * K2 * pref * _F0(p1 * q1 / (p1 + q1) * (P1 - Q1) ** 2))
    return S, T, V, G


def h_atom():
    # one-centre: T_aa + attraction to own nucleus only
    a = ALPHA
    c = CN[:, None] * CN[None, :]
    p = a[:, None] + a[None, :]
    return float(np.sum(c * 3 * (a[:, None] * a[None, :] / p) * (np.pi / p) ** 1.5)
                 + np.sum(c * (-2 * np.pi / p)))


def solve(R):
    """Closed-form RHF / singlet FCI / triplet energies at R (bohr) plus VB data."""
    S, T, V, G = ao_integrals(R)
    s = S[0, 1]
    saa = S[0, 0]  # = 0.99999999 with the 8-digit STO-3G constants, not exactly 1
    H = T + V
    Cg = np.array([1, 1]) / np.sqrt(2 * (saa + s))
    Cu = np.array([1, -1]) / np.sqrt(2 * (saa - s))
    C = np.column_stack([Cg, Cu])
    h = C.T @ H @ C
    g = np.einsum("mnls,mp,nq,lr,st->pqrt", G, C, C, C, C)
    Enuc = 1 / R
    E_rhf = 2 * h[0, 0] + g[0, 0, 0, 0] + Enuc
    Hci = np.array([[2 * h[0, 0] + g[0, 0, 0, 0], g[0, 1, 0, 1]],
                    [g[0, 1, 0, 1], 2 * h[1, 1] + g[1, 1, 1, 1]]]) + Enuc * np.eye(2)
    w, v = np.linalg.eigh(Hci)
    cg, cu = v[:, 0] * np.sign(v[0, 0])
    E_trip = h[0, 0] + h[1, 1] + g[0, 0, 1, 1] - g[0, 1, 0, 1] + Enuc
    # VB: Psi_spatial = cg g g + cu u u = sum M_mn chi_m chi_n
    M = cg * np.outer(Cg, Cg) + cu * np.outer(Cu, Cu)
    Nst = np.sqrt(2 * (saa * saa + s * s))
    x = np.array([M[0, 1], M[0, 0]]) * Nst  # [cov, ion] coefficients of normalised structures
    o = 2 * s * saa / (saa * saa + s * s)
    Sv = np.array([[1, o], [o, 1]])
    x = x / np.sqrt(x @ Sv @ x)
    cc = x * (Sv @ x)
    ev, U = np.linalg.eigh(Sv)
    lw = (U @ np.diag(np.sqrt(ev)) @ U.T @ x) ** 2
    return dict(S=S, T=T, V=V, G=G, E_RHF=E_rhf, E_FCI=w[0], E_triplet=E_trip,
                c_gg=cg, c_uu=cu, ion_CC=cc[1], ion_L=lw[1], cov_CC=cc[0], cov_L=lw[0],
                eps_g=h[0, 0] + g[0, 0, 0, 0],
                eps_u=h[1, 1] + 2 * g[0, 0, 1, 1] - g[0, 1, 0, 1])
