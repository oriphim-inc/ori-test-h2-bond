"""Analytic integrals over s-type Gaussian primitives and STO-3G contractions.

Formulas: Szabo & Ostlund, Modern Quantum Chemistry (1989), Appendix A.
All quantities in atomic units (bohr, hartree).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import erf

# STO-3G H 1s, zeta = 1.24, as agreed in PLAN.md (Hehre, Stewart & Pople, JCP 51, 2657
# (1969); Szabo & Ostlund 1989 sec. 3.6). Other zeta: exponents scale as (zeta/1.24)**2.
ZETA_H = 1.24
STO3G_ALPHA_H = np.array([3.42525091, 0.62391373, 0.16885540])
STO3G_D = np.array([0.15432897, 0.53532814, 0.44463454])

BOHR_TO_ANG = 0.529177210903  # CODATA 2018
ANG_TO_BOHR = 1.0 / BOHR_TO_ANG
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.5094740631


@dataclass(frozen=True)
class ContractedS:
    """Contracted s function: sum_k d_k * N_k * exp(-alpha_k |r - A|^2)."""

    center: np.ndarray  # (3,) bohr
    alpha: np.ndarray
    d: np.ndarray

    @property
    def norm_coeffs(self) -> np.ndarray:
        """d_k multiplied by the primitive normalisation (2 alpha / pi)^(3/4)."""
        return self.d * (2.0 * self.alpha / math.pi) ** 0.75


def sto3g_h(center, zeta: float = ZETA_H) -> ContractedS:
    return ContractedS(np.asarray(center, float), STO3G_ALPHA_H * (zeta / ZETA_H) ** 2,
                      STO3G_D.copy())


def boys0(t: float) -> float:
    """F0(t) = int_0^1 exp(-t u^2) du."""
    if t < 1e-8:
        return 1.0 - t / 3.0 + t * t / 10.0
    return 0.5 * math.sqrt(math.pi / t) * erf(math.sqrt(t))


# ---- primitive (unnormalised) integrals --------------------------------------
def _s(a, A, b, B):
    p = a + b
    rab2 = float(np.dot(A - B, A - B))
    return (math.pi / p) ** 1.5 * math.exp(-a * b / p * rab2)


def _t(a, A, b, B):
    p = a + b
    rab2 = float(np.dot(A - B, A - B))
    mu = a * b / p
    return mu * (3.0 - 2.0 * mu * rab2) * (math.pi / p) ** 1.5 * math.exp(-mu * rab2)


def _v(a, A, b, B, C, Z):
    p = a + b
    P = (a * A + b * B) / p
    rab2 = float(np.dot(A - B, A - B))
    rpc2 = float(np.dot(P - C, P - C))
    return -2.0 * math.pi / p * Z * math.exp(-a * b / p * rab2) * boys0(p * rpc2)


def _eri(a, A, b, B, c, C, d, D):
    p, q = a + b, c + d
    P = (a * A + b * B) / p
    Q = (c * C + d * D) / q
    rab2 = float(np.dot(A - B, A - B))
    rcd2 = float(np.dot(C - D, C - D))
    rpq2 = float(np.dot(P - Q, P - Q))
    pref = 2.0 * math.pi**2.5 / (p * q * math.sqrt(p + q))
    return pref * math.exp(-a * b / p * rab2 - c * d / q * rcd2) * boys0(p * q / (p + q) * rpq2)


# ---- contracted integrals ----------------------------------------------------
def _contract2(f, g1: ContractedS, g2: ContractedS, *extra):
    c1, c2 = g1.norm_coeffs, g2.norm_coeffs
    tot = 0.0
    for i, a in enumerate(g1.alpha):
        for j, b in enumerate(g2.alpha):
            tot += c1[i] * c2[j] * f(a, g1.center, b, g2.center, *extra)
    return tot


def overlap(basis):
    n = len(basis)
    return np.array([[_contract2(_s, basis[i], basis[j]) for j in range(n)] for i in range(n)])


def kinetic(basis):
    n = len(basis)
    return np.array([[_contract2(_t, basis[i], basis[j]) for j in range(n)] for i in range(n)])


def nuclear(basis, charges, centers):
    n = len(basis)
    V = np.zeros((n, n))
    for Z, C in zip(charges, centers):
        C = np.asarray(C, float)
        V += np.array(
            [[_contract2(_v, basis[i], basis[j], C, Z) for j in range(n)] for i in range(n)]
        )
    return V


def eri(basis):
    """(mu nu | lam sig) in chemists' notation, full 4-index array."""
    n = len(basis)
    G = np.zeros((n, n, n, n))
    for m in range(n):
        for nu in range(n):
            for l in range(n):
                for s in range(n):
                    gm, gn, gl, gs = basis[m], basis[nu], basis[l], basis[s]
                    tot = 0.0
                    for i, a in enumerate(gm.alpha):
                        for j, b in enumerate(gn.alpha):
                            for k, c in enumerate(gl.alpha):
                                for q, d in enumerate(gs.alpha):
                                    tot += (
                                        gm.norm_coeffs[i] * gn.norm_coeffs[j]
                                        * gl.norm_coeffs[k] * gs.norm_coeffs[q]
                                        * _eri(a, gm.center, b, gn.center,
                                               c, gl.center, d, gs.center)
                                    )
                    G[m, nu, l, s] = tot
    return G


@dataclass
class Integrals:
    S: np.ndarray
    T: np.ndarray
    V: np.ndarray
    G: np.ndarray
    E_nuc: float
    basis: list

    @property
    def H(self):
        return self.T + self.V


def h2_integrals(R_bohr: float) -> Integrals:
    A = np.array([0.0, 0.0, 0.0])
    B = np.array([0.0, 0.0, R_bohr])
    basis = [sto3g_h(A), sto3g_h(B)]
    S = overlap(basis)
    T = kinetic(basis)
    V = nuclear(basis, [1.0, 1.0], [A, B])
    G = eri(basis)
    return Integrals(S, T, V, G, 1.0 / R_bohr, basis)


def h_atom_energy(zeta: float = ZETA_H) -> float:
    """One electron in one STO-3G 1s on a bare proton: E = T_11 + V_11."""
    g = [sto3g_h([0.0, 0.0, 0.0], zeta)]
    return float(kinetic(g)[0, 0] + nuclear(g, [1.0], [[0.0, 0.0, 0.0]])[0, 0])
