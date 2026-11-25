> **The hydrogen molecule bond.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)

# build minimal quantum chemistry

By Oriphim, with Ori.

## What was asked
Build a quantum-chemistry code from scratch for H2 in STO-3G. It should compute the overlap,
kinetic, nuclear-attraction and two-electron integrals analytically over Gaussian primitives,
run RHF with an SCF loop and then full CI in the 2-electron, 2-orbital space, and scan
R = 0.4–5.0 Å. The deliverables are both potential-energy curves, their equilibrium values,
the dissociation-limit error compared with two STO-3G atoms (−0.4666 Ha each), and the ionic
versus covalent weight as a function of R. The agreed details are in [PLAN.md](PLAN.md).

## Answer
| | R_e (Å) | R_e (bohr) | E_min (Ha) | D_e vs 2 E_H (Ha / eV) | E(5.0 Å) − 2 E_H (Ha) |
|---|---|---|---|---|---|
| RHF | **0.7122** | 1.3459 | **−1.117506** | 0.1843 / 5.02 | **+0.3341** |
| FCI | **0.7349** | 1.3887 | **−1.137306** | 0.2041 / 5.55 | **−7×10⁻⁸** (≈ 0) |

- **Separated atoms:** the code's own STO-3G hydrogen atom has E_H = −0.466582 Ha,
  which agrees with the stated −0.4666. So 2 E_H = −0.933164 Ha.
- **FCI** dissociates correctly. At 5.0 Å it is within 1×10⁻⁷ Ha of 2 E_H, and the
  singlet and triplet are degenerate there to 8×10⁻⁸ Ha.
- **RHF** dissociates badly. At 5.0 Å it lies **0.334 Ha (9.1 eV) above two atoms**, and it is
  still rising there, towards about +0.387 Ha as R → ∞. The RHF number is therefore a
  value *at 5.0 Å*, not the limit. Measured from its own value at 5.0 Å, the RHF well would
  be 0.518 Ha deep, which has no physical meaning.
- **Ionic character:** the RHF wavefunction is exactly **50 % ionic at every R**, because a
  doubly occupied σg puts H⁺H⁻ and H·H· in fixed proportion. The FCI ground state removes
  this by mixing in σu² (shown as the dotted curve in the plot). Its σu² coefficient goes
  to −σg², so the ionic weight goes to zero:

| R (Å) | E_RHF | E_FCI | FCI ionic (Chirgwin–Coulson) | FCI ionic (Löwdin) | RHF ionic |
|---|---|---|---|---|---|
| 0.40 | −0.90436 | −0.91415 | 0.059 | 0.441 | 0.500 |
| 0.74 | −1.11676 | −1.13728 | 0.216 | 0.388 | 0.500 |
| 1.00 | −1.06611 | −1.10115 | 0.215 | 0.327 | 0.500 |
| 1.50 | −0.91087 | −0.99815 | 0.121 | 0.168 | 0.500 |
| 2.00 | −0.78379 | −0.94864 | 0.034 | 0.047 | 0.500 |
| 2.50 | −0.70294 | −0.93605 | 0.0066 | 0.0090 | 0.500 |
| 3.00 | −0.65605 | −0.93363 | 0.0011 | 0.0014 | 0.500 |
| 5.00 | −0.59902 | −0.93316 | 1.2×10⁻⁷ | 1.4×10⁻⁷ | 0.500 |

  Beyond about 1.5 Å the two weight definitions agree. At short R the two VB structures
  overlap so strongly that the weight depends on how it is defined: the Chirgwin–Coulson
  curve peaks near 0.86 Å and falls towards R → 0. The trend, not the short-R numbers, is
  the robust result.

Figures: [results/pes.png](results/pes.png) and
[results/ionic_weights.png](results/ionic_weights.png). All the numbers are in
`results/summary.json`, `results/scan.csv` (461 points, 0.01 Å step) and
`results/reference_point.json`.

## Evidence
- **Textbook values:** at R = 1.4 bohr the code reproduces the Szabo & Ostlund values to all
  printed digits: S₁₂ = 0.6593, T₁₁ = 0.7600, T₁₂ = 0.2365, (11|11) = 0.7746, (11|22) = 0.5697,
  (21|11) = 0.4441, (21|21) = 0.2970, E_RHF = −1.1167 Ha, E_FCI = −1.1373 Ha, ε = −0.578 / +0.670.
- **Independent reference:** the Sceptic wrote its own integrals and a closed-form,
  symmetry-adapted RHF/FCI solution (`tests/checks/_ref.py`). Over the whole scan, the
  curves, triplet, orbital energies, CI coefficients and weights agree with it to about
  1×10⁻¹⁴ Ha. All 48 tests and checks pass.
- **Where the minima come from:** R_e is found by Brent minimisation of the energy itself.
  A 47-point scan gives the same R_e and E_min.
- **Spin:** the FCI is diagonalised separately in the singlet and triplet blocks, so the
  ground state stays a pure singlet even where the two states become degenerate.
- **SCF convergence:** the SCF converges at every R. For H₂ in a minimal basis the MOs are
  fixed by symmetry, so it takes 2 iterations.

## Limits
See [LIMITS.md](LIMITS.md). In short:
- These are minimal-basis numbers. The exact H₂ values are E_min = −1.1745 Ha and
  D_e ≈ 4.75 eV, so STO-3G gives a well about 0.8 eV too deep.
- The RHF dissociation error is a value at 5 Å, not the limit.
- The ionic/covalent split depends on how the weight is defined at short R.
- The 8-digit basis constants leave the 1s normalised only to 1 − 9×10⁻⁹. This affects E_H
  and the weights at the 10⁻⁸ level.
- Cosmetic: the FCI label in pes.png overlaps the axis.

## Re-run
```
python3 -m venv .venv && .venv/bin/pip install numpy scipy matplotlib pytest
.venv/bin/python run.py            # 461 points, ~6 s; --n N to change the grid
.venv/bin/python -m pytest tests   # unit tests + the Sceptic's checks
```
The code is in `src/h2qc/`:
- `integrals.py` — analytic integrals over s Gaussians (Boys F₀).
- `hf.py` — RHF SCF.
- `fci.py` — spin-adapted 2×2-orbital FCI.
- `vb.py` — ionic/covalent decomposition.
- `scan.py` — the R scan and equilibrium search.
