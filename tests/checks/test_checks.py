"""Sceptic's pass/fail checks for H2 / STO-3G RHF + FCI (written before the main run).

Run:  .venv/bin/python -m pytest tests/checks -q
Every test skips until results/ exist. Independent reference: tests/checks/_ref.py
(own integrals + symmetry-fixed closed-form RHF/FCI; does not import h2qc).

Tolerance rationale
- vs Szabo & Ostlund (1989) tabulated values: they print 4 decimals -> 5e-5 rounding;
  allow 1e-4 (2e-3 for orbital energies printed to 3 decimals).
- vs the independent reference: the same physics evaluated two ways; differences should be
  round-off plus SCF convergence (|dE| < 1e-12): energies 1e-9 Ha, coefficients/weights 1e-8.
- R_e: Brent xatol 1e-7 A, plus energy noise ~1e-12 Ha on a curvature of ~1.3 Ha/A^2
  gives ~1e-6 A; allow 2e-5 A.
- literature R_e (HF/STO-3G 0.712 A, FCI/STO-3G 0.735 A, e.g. NIST CCCBDB): quoted to
  3 decimals -> allow 1e-3 A; E_min(RHF) = -1.1175 Ha to 1e-4.
"""
from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from _ref import BOHR_ANG, h_atom, solve

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
pytestmark = pytest.mark.skipif(
    not all((RES / f).exists() for f in ("scan.csv", "summary.json", "reference_point.json")),
    reason="results/ not produced yet (run `python run.py`)")


@pytest.fixture(scope="module")
def rows():
    with open(RES / "scan.csv") as f:
        rd = list(csv.DictReader(f))
    out = {k: np.array([float(r[k]) if r[k] not in ("True", "False") else r[k] == "True"
                        for r in rd]) for k in rd[0]}
    return out


@pytest.fixture(scope="module")
def summ():
    return json.loads((RES / "summary.json").read_text())


@pytest.fixture(scope="module")
def refpt():
    return json.loads((RES / "reference_point.json").read_text())


@pytest.fixture(scope="module")
def ref(rows):
    return [solve(R / BOHR_ANG) for R in rows["R_ang"]]


def col(ref, key):
    return np.array([r[key] for r in ref])


# ---------- 1. anchors: Szabo & Ostlund (1989), R = 1.4 bohr --------------------------
SO = {"S12": 0.6593, "T11": 0.7600, "T12": 0.2365, "(11|11)": 0.7746, "(11|22)": 0.5697,
      "(21|11)": 0.4441, "(21|21)": 0.2970, "E_RHF": -1.1167, "E_FCI": -1.1373,
      # H^core_11 = -1.1204, H^core_12 = -0.9584 (S&O 3.5.2) -> V = H - T
      "V11_total": -1.8804, "V12_total": -1.1949}


@pytest.mark.parametrize("key", list(SO))
def test_szabo_ostlund_anchor(refpt, key):
    assert refpt["R_bohr"] == pytest.approx(1.4, abs=1e-12)
    assert refpt[key] == pytest.approx(SO[key], abs=1e-4), key


def test_szabo_ostlund_orbital_energies(refpt):
    # S&O: eps_1 = -0.578, eps_2 = 0.670 at R = 1.4 bohr
    assert refpt["eps"][0] == pytest.approx(-0.578, abs=2e-3)
    assert refpt["eps"][1] == pytest.approx(0.670, abs=2e-3)


def test_h_atom_energy(summ):
    # STO-3G H atom: -0.46658 Ha (S&O); also equal to independent evaluation
    assert summ["E_H_atom_computed"] == pytest.approx(-0.4666, abs=1e-4)
    assert summ["E_H_atom_computed"] == pytest.approx(h_atom(), abs=1e-10)
    assert summ["E_2H"] == pytest.approx(2 * summ["E_H_atom_computed"], abs=1e-14)


# ---------- 2. grid and units ---------------------------------------------------------
def test_grid_and_units(rows):
    R = rows["R_ang"]
    assert len(R) == 461
    assert R[0] == pytest.approx(0.4, abs=1e-12) and R[-1] == pytest.approx(5.0, abs=1e-12)
    assert np.allclose(np.diff(R), 0.01, atol=1e-12)
    assert np.allclose(rows["R_bohr"] * BOHR_ANG, R, rtol=1e-12)  # bohr, not cgs/other


# ---------- 3. whole curves vs independent reference ----------------------------------
@pytest.mark.parametrize("key,tol", [("S_ab", 1e-10), ("E_RHF", 1e-9), ("E_FCI", 1e-9),
                                     ("E_triplet", 1e-9), ("eps_g", 1e-8), ("eps_u", 1e-8)])
def test_curve_matches_independent_reference(rows, ref, key, tol):
    if key == "S_ab":
        mine = np.array([r["S"][0, 1] for r in ref])
    else:
        mine = col(ref, key)
    err = np.abs(rows[key] - mine)
    assert err.max() < tol, f"{key}: max |diff| {err.max():.2e} at R={rows['R_ang'][err.argmax()]}"


def test_fci_coefficients_vs_reference(rows, ref):
    assert np.abs(rows["c_gg"] - col(ref, "c_gg")).max() < 1e-8
    assert np.abs(rows["c_uu"] - col(ref, "c_uu")).max() < 1e-8


# ---------- 4. SCF ----------------------------------------------------------------------
def test_scf_converged_everywhere(rows, summ):
    assert rows["scf_converged"].all()
    assert summ.get("all_scf_converged", True) is True
    assert rows["scf_iter"].max() < 200  # below max_iter: converged, not stopped


# ---------- 5. variational / symmetry / normalisation ---------------------------------
def test_fci_below_rhf(rows):
    assert np.all(rows["E_FCI"] <= rows["E_RHF"] + 1e-12)


def test_fci_above_exact_h2(summ, rows):
    # exact non-relativistic BO minimum of H2: -1.174475 Ha (Kolos & Wolniewicz);
    # a finite-basis FCI must lie above the exact curve, and above -1 Ha (2 exact H) at large R
    assert summ["FCI"]["E_min"] > -1.174475
    assert np.all(rows["E_FCI"] > -1.174475)
    assert rows["E_FCI"][-1] > -1.0


def test_fci_singlet_structure(rows):
    # g/u mixing forbidden by inversion symmetry; CSF normalised; u^2 enters with opposite sign
    assert np.abs(rows["c_gu"]).max() < 1e-10 and np.abs(rows["c_ug"]).max() < 1e-10
    assert np.allclose(rows["c_gg"] ** 2 + rows["c_uu"] ** 2, 1.0, atol=1e-10)
    assert np.all(rows["c_gg"] > 0) and np.all(rows["c_uu"] < 0)


def test_triplet_above_singlet_and_repulsive(rows, summ):
    assert np.all(rows["E_triplet"] >= rows["E_FCI"] - 1e-12)
    assert np.all(rows["E_triplet"] > summ["E_2H"] - 1e-9)        # unbound in this basis
    assert np.all(np.diff(rows["E_triplet"]) < 0)                 # monotonically repulsive


# ---------- 6. dissociation limits ------------------------------------------------------
def test_fci_dissociates_to_two_atoms(rows, summ):
    assert abs(summ["FCI"]["dissociation_error_Ha"]) < 1e-6
    assert abs(rows["E_FCI"][-1] - summ["E_2H"]) < 1e-6


def test_singlet_triplet_degenerate_at_Rmax(rows):
    assert abs(rows["E_triplet"][-1] - rows["E_FCI"][-1]) < 1e-6
    assert abs(rows["c_gg"][-1] + rows["c_uu"][-1]) < 1e-3   # c_uu -> -c_gg -> -1/sqrt2


def test_rhf_dissociation_error_analytic(rows, summ):
    # RHF -> 2E_H + (aa|aa)/2 - 1/(2R) for non-overlapping atoms (half ionic H+ H- weight);
    # O(S) corrections at 5 A are ~3e-4 Ha, so 1e-3; exact match is checked in section 3.
    R = rows["R_bohr"][-1]
    expected = 0.5 * 0.7746059 - 0.5 / R
    assert summ["RHF"]["dissociation_error_Ha"] == pytest.approx(expected, abs=1e-3)
    assert summ["RHF"]["dissociation_error_Ha"] > 0.3


# ---------- 7. equilibrium and well depth -----------------------------------------------
@pytest.mark.parametrize("tag,lit_Re", [("RHF", 0.712), ("FCI", 0.735)])
def test_equilibrium(summ, rows, tag, lit_Re):
    s = summ[tag]
    key = "E_" + tag
    from scipy.optimize import minimize_scalar
    m = minimize_scalar(lambda x: solve(x / BOHR_ANG)[key], bounds=(0.6, 0.85),
                        method="bounded", options={"xatol": 1e-9})
    assert s["R_e_ang"] == pytest.approx(m.x, abs=2e-5)
    assert s["E_min"] == pytest.approx(m.fun, abs=1e-9)
    assert s["R_e_ang"] == pytest.approx(lit_Re, abs=1e-3)
    assert s["R_e_bohr"] == pytest.approx(s["R_e_ang"] / BOHR_ANG, rel=1e-12)
    assert s["E_min"] <= rows[key].min() + 1e-12                 # Brent beat the grid
    assert abs(s["grid_min_R_ang"] - s["R_e_ang"]) <= 0.0051


def test_rhf_emin_literature(summ):
    assert summ["RHF"]["E_min"] == pytest.approx(-1.1175, abs=1e-4)


@pytest.mark.parametrize("tag", ["RHF", "FCI"])
def test_well_depth_bookkeeping(summ, rows, tag):
    s = summ[tag]
    key = "E_" + tag
    assert s["E_at_Rmax"] == pytest.approx(rows[key][-1], abs=1e-14)
    assert s["D_e_vs_2H_Ha"] == pytest.approx(summ["E_2H"] - s["E_min"], abs=1e-12)
    assert s["D_e_vs_own_end_Ha"] == pytest.approx(s["E_at_Rmax"] - s["E_min"], abs=1e-12)
    assert s["dissociation_error_Ha"] == pytest.approx(s["E_at_Rmax"] - summ["E_2H"], abs=1e-12)
    if "D_e_vs_2H_eV" in s:
        assert s["D_e_vs_2H_eV"] == pytest.approx(s["D_e_vs_2H_Ha"] * 27.211386245988, rel=1e-10)


def test_grid_independence_of_Re(summ):
    # a coarse scan (47 points, 0.1 A) must give the same Brent R_e: result is not a grid artefact
    from h2qc.scan import run_scan
    _, s2 = run_scan(47)
    for tag in ("RHF", "FCI"):
        assert s2[tag]["R_e_ang"] == pytest.approx(summ[tag]["R_e_ang"], abs=2e-5)
        assert s2[tag]["E_min"] == pytest.approx(summ[tag]["E_min"], abs=1e-10)


# ---------- 8. ionic / covalent weights -------------------------------------------------
def test_weights_sum_to_one_and_bounded(rows):
    for m in ("FCI", "RHF"):
        for d in ("CC", "L"):
            ion, cov = rows[f"{m}_ion_{d}"], rows[f"{m}_cov_{d}"]
            assert np.allclose(ion + cov, 1.0, atol=1e-10)
            assert np.all(ion >= -1e-12) and np.all(ion <= 1 + 1e-12)


def test_rhf_weights_exactly_half(rows):
    # sigma_g^2 = [(ab+ba) + (aa+bb)]/norm: both normalised structures have equal coefficient,
    # so CC and Loewdin ionic weights are exactly 1/2 at every R.
    assert np.abs(rows["RHF_ion_CC"] - 0.5).max() < 1e-10
    assert np.abs(rows["RHF_ion_L"] - 0.5).max() < 1e-10


def test_fci_weights_vs_reference(rows, ref):
    assert np.abs(rows["FCI_ion_CC"] - col(ref, "ion_CC")).max() < 1e-8
    assert np.abs(rows["FCI_ion_L"] - col(ref, "ion_L")).max() < 1e-8


def test_fci_ionic_weight_trend(rows):
    R = rows["R_ang"]
    assert np.all(rows["FCI_ion_L"] <= 0.5 + 1e-12) and np.all(rows["FCI_ion_CC"] <= 0.5)
    assert np.all(np.diff(rows["FCI_ion_L"]) <= 1e-12)              # Loewdin: monotone in R
    far = R >= 1.0                                                   # CC: monotone beyond 1 A
    assert np.all(np.diff(rows["FCI_ion_CC"][far]) <= 1e-12)        # (CC peaks near 0.86 A)
    assert rows["FCI_ion_CC"][-1] < 1e-5 and rows["FCI_ion_L"][-1] < 1e-5


# ---------- 9. no hard-coded targets ----------------------------------------------------
def test_no_target_literals_in_code():
    pat = re.compile(r"1\.1167|1\.1373|1\.1175|0\.712\d|0\.73[45]\d|\b0\.2157|\b0\.388")
    for f in list((ROOT / "src").rglob("*.py")) + [ROOT / "run.py"]:
        hits = [ln for ln in f.read_text().splitlines() if pat.search(ln)]
        assert not hits, f"{f}: {hits}"
