"""Unit tests for the from-scratch H2 / STO-3G code."""
import math

import numpy as np
import pytest
from scipy.integrate import quad

from h2qc.fci import fci
from h2qc.hf import rhf
from h2qc.integrals import boys0, h2_integrals, h_atom_energy, sto3g_h, overlap
from h2qc.vb import vb_weights


def test_boys_small_and_large():
    for t in [0.0, 1e-9, 1e-3, 0.5, 3.0, 30.0]:
        ref = quad(lambda u: math.exp(-t * u * u), 0, 1)[0]
        assert boys0(t) == pytest.approx(ref, rel=1e-10, abs=1e-14)


def test_basis_normalised():
    S = overlap([sto3g_h([0, 0, 0])])
    assert S[0, 0] == pytest.approx(1.0, abs=1e-6)


def test_h_atom_energy():
    assert h_atom_energy() == pytest.approx(-0.4666, abs=1e-4)


def test_spin_adapted_root_counts():
    # M_s=0 block of 2 orbitals: 3 singlets + 1 triplet component
    ints = h2_integrals(1.4)
    ci = fci(ints, rhf(ints).C)
    assert len(ci.energies) == 4 and ci.singlet_flags.sum() == 3


def test_eri_symmetry():
    G = h2_integrals(1.4).G
    assert np.allclose(G, G.transpose(1, 0, 2, 3))
    assert np.allclose(G, G.transpose(0, 1, 3, 2))
    assert np.allclose(G, G.transpose(2, 3, 0, 1))


def test_large_R_eri_limit():
    # (aa|bb) -> 1/R for well-separated charge clouds
    R = 30.0
    G = h2_integrals(R).G
    assert G[0, 0, 1, 1] == pytest.approx(1 / R, rel=1e-6)


def test_fci_variational_and_singlet():
    for R in [0.8, 1.4, 3.0, 6.0]:
        ints = h2_integrals(R)
        hf = rhf(ints)
        ci = fci(ints, hf.C)
        assert hf.converged
        assert ci.E_total <= hf.E_total + 1e-12
        c = ci.coeff
        assert np.allclose(c, c.T, atol=1e-10)
        assert abs(np.linalg.norm(c) - 1) < 1e-12


def test_vb_weights_sum_to_one():
    ints = h2_integrals(2.0)
    hf = rhf(ints)
    ci = fci(ints, hf.C)
    w = vb_weights(hf.C, ci.coeff, ints.S[0, 1])
    assert w["cov_CC"] + w["ion_CC"] == pytest.approx(1.0, abs=1e-12)
    assert w["cov_L"] + w["ion_L"] == pytest.approx(1.0, abs=1e-12)
