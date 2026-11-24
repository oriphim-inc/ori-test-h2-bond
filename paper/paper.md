# Why one determinant cannot break a bond: H$_2$ in a minimal basis, from scratch

*Research question: Oriphim held-out test set. Computation, analysis and drafting: Ori (Oriphim), with an independent review.*

## Abstract

We wrote a small quantum-chemistry code from scratch (numpy/scipy only) for the hydrogen molecule in the minimal STO-3G basis. It computes all one- and two-electron integrals analytically over Gaussian functions, solves restricted Hartree–Fock (RHF) by a self-consistent-field loop, and then solves full configuration interaction (FCI) in the two-orbital space. Both are scanned over 461 bond lengths from 0.4 to 5.0 Å. RHF gives an equilibrium bond length of 0.7122 Å and FCI gives 0.7349 Å. Their minimum energies are -1.11751 and -1.13731 hartree (Ha). At 5.0 Å the FCI curve lies within $10^{-7}$ Ha of two isolated STO-3G hydrogen atoms (computed as -0.466582 Ha each). RHF lies 0.3341 Ha above them and is still rising. The cause shows up in a valence-bond decomposition. The RHF wavefunction is exactly 50 % ionic (H$^+$H$^-$) at every distance. The FCI ionic weight (Chirgwin–Coulson) falls from about 0.2 near equilibrium to 1.2e-07 at 5.0 Å. An independent reference implementation reproduces every energy curve to about $10^{-14}$ Ha, and the ionic weights to about $10^{-8}$.

## 1. Question and why it matters

Most quantum chemistry starts from a single Slater determinant: each pair of electrons shares one orbital (the Hartree–Fock picture). For H$_2$ the two electrons share the bonding orbital $\sigma_g = a + b$, where $a$ and $b$ are the 1s functions on the two atoms. Multiplying this out gives $\sigma_g^2 = (ab + ba) + (aa + bb)$. The first part is covalent (one electron on each atom). The second is ionic (both electrons on one atom). One determinant fixes the two parts in equal proportion. That is acceptable near equilibrium, but when the bond is stretched the molecule should separate into two neutral atoms, not into a 50/50 mix of neutral atoms and H$^+$ + H$^-$. This failure of a single determinant at dissociation is the textbook example of *static correlation*. It is the reason bond-breaking, transition states and many catalytic and magnetic systems need multi-determinant methods. In the smallest basis the problem can be solved exactly, so here we measure the failure directly. We report:

1. both energy curves;
2. their equilibrium bond length and well depth;
3. the error at large separation relative to two STO-3G atoms (stated as $2\times(-0.4666)$ Ha);
4. how the ionic weight of the exact (FCI) ground state changes with bond length.

## 2. Method in brief

- **Basis and integrals.** Each atom carries one STO-3G 1s function: three Gaussians with the standard exponents and coefficients for Slater exponent $\zeta = 1.24$ (Hehre, Stewart & Pople 1969). The overlap, kinetic, nuclear-attraction and electron-repulsion integrals are computed analytically from the closed forms for s-type Gaussians (Szabo & Ostlund 1989, App. A), using the Boys function $F_0$.
- **RHF.** Closed-shell Roothaan equations with symmetric orthogonalisation and a core-Hamiltonian starting guess. The loop stops at $|\Delta E| < 10^{-12}$ Ha.
- **FCI.** All determinants of two electrons in the four spin orbitals built from $\sigma_g$ and $\sigma_u$. The Hamiltonian is diagonalised separately in the singlet and triplet blocks. This keeps the ground state a pure singlet where the two states become degenerate at large separation.
- **Equilibrium.** Bounded Brent minimisation of the energy function itself, starting from the grid minimum. Well depth $D_e = 2E_H - E_\text{min}$, where $E_H$ is the code's own STO-3G atom energy.
- **Ionic weight.** The singlet wavefunction is rewritten in the atomic basis as $C_\text{cov}(ab+ba) + C_\text{ion}(aa+bb)$. Because the two structures are not orthogonal, the weight is not uniquely defined. We report Chirgwin–Coulson weights and give Löwdin-orthogonalised weights alongside.

## 3. Results

| | Bond length (Å) | Minimum energy (Ha) | Well depth vs $2E_H$ (eV) | $E(5\,\text{Å}) - 2E_H$ (Ha) |
|---|---|---|---|---|
| RHF | 0.7122 | -1.117506 | 5.02 | +0.3341 |
| FCI | 0.7349 | -1.137306 | 5.55 | -7.1e-08 |

![Potential energy curves of H$_2$/STO-3G. RHF (blue) and FCI (orange), with the separated-atom energy $2E_H$ (dashed). The two curves agree near equilibrium. Beyond about 1 Å RHF climbs far above the separated-atom energy, while FCI approaches it. (Cosmetic: the FCI annotation is clipped at the bottom and the upper axis limit leaves empty space.)](../results/pes.png)

**Energy curves.** Near equilibrium the two methods are close: -1.1167 vs -1.1373 Ha at 1.4 bohr, and correlation lengthens the bond by about 0.02 Å. At large separation they diverge qualitatively (Fig. 1). FCI reaches the two-atom energy: at 5.0 Å the singlet (-0.9331638 Ha) and the triplet (-0.9331637 Ha) are degenerate to about $10^{-7}$ Ha, as they must be for two independent atoms. RHF ends 0.334 Ha (about 9 eV) too high at 5.0 Å. This is not its limit: the RHF error keeps growing towards half the on-site repulsion, $\tfrac12(aa|aa)$, which the Sceptic evaluated as +0.387 Ha at very large separation.

![Ionic weight of the ground state against bond length. RHF (blue) is exactly 0.5 everywhere. FCI (orange; solid Chirgwin–Coulson, dashed Löwdin) falls to zero as the bond breaks, as the FCI wavefunction mixes in the doubly-excited $\sigma_u^2$ configuration (dotted, its squared coefficient rising to 0.5).](../results/ionic_weights.png)

**Ionic versus covalent.** The RHF ionic weight is exactly 0.5 at all 461 bond lengths (Fig. 2). FCI removes the ionic part by mixing in $\sigma_u^2$ with a coefficient that tends to minus that of $\sigma_g^2$ (at 5.0 Å: 0.7074 and -0.7068). Since $\sigma_g^2 - \sigma_u^2 \propto ab + ba$, the ionic terms cancel:

| Bond length (Å) | FCI ionic (Chirgwin–Coulson) | FCI ionic (Löwdin) | RHF ionic |
|---|---|---|---|
| 0.40 | 0.059 | 0.441 | 0.500 |
| 0.74 | 0.216 | 0.388 | 0.500 |
| 1.50 | 0.121 | 0.168 | 0.500 |
| 2.00 | 0.034 | 0.047 | 0.500 |
| 3.00 | 0.0011 | 0.0014 | 0.500 |
| 5.00 | 1.2e-07 | 1.4e-07 | 0.500 |

The two definitions differ at every distance: between 1.5 and 3 Å the Löwdin weight is 30–40 % larger, and the gap narrows to about 10 % at 5 Å as both fall towards zero. The disagreement is largest at short distance (0.06 vs 0.44 at 0.4 Å), because the two structures become nearly identical when the atoms overlap strongly.

## 4. How it was checked

- **Textbook anchors.** At 1.4 bohr the code reproduces every integral and total energy printed by Szabo & Ostlund to the four decimals printed, e.g. overlap 0.6593, kinetic $T_{11}$ = 0.7600, $E_\text{RHF}$ = -1.1167 Ha, $E_\text{FCI}$ = -1.1373 Ha. The bond lengths agree with tabulated HF/STO-3G (0.712 Å) and FCI/STO-3G (0.735 Å) values to within 0.001 Å.
- **Independent reference.** Before the main run, an independent reviewer (the Sceptic) wrote its own integrals and a closed-form symmetry-adapted RHF/FCI solution (`tests/checks/_ref.py`), together with 40 pass/fail checks fixed in advance (run with 8 unit tests, 48 in all). The checks cover the whole curves, the variational ordering, singlet purity, triplet behaviour, sum rules on the weights, grid independence of the bond lengths, and a search for hard-coded target values in the source.
- **What failed, and what was changed.** The first full run passed 41 of 48 tests. The seven failures were differences of $10^{-8}$–$10^{-7}$ Ha. They were traced to the *reference*: with the 8-digit basis constants the 1s function is normalised only to $1 - 9\times10^{-9}$, and the reference assumed exactly 1. After the reference was corrected to use the true self-overlap, all 48 tests pass and the energy curves agree to about $10^{-14}$ Ha. No tolerance was loosened. Before that run, the code itself was changed in two ways: it now uses the basis constants exactly as listed in the plan, and the FCI is diagonalised in separate spin blocks so that the nearly degenerate singlet and triplet cannot mix at large separation.
- **The Sceptic's remaining flags** (none judged suspect): the code's atom energy and ionic weights also treat the self-overlap as 1, which shifts them by $10^{-8}$ or less; the RHF dissociation error is a value at 5 Å, not a limit; the weights depend on the definition at short range; and the Fig. 1 annotation is clipped.

## 5. What this does not show

- **These are minimal-basis numbers, not real H$_2$.** The well depth from FCI, 5.55 eV, is about 0.8 eV deeper than the exact value of about 4.75 eV, and the FCI minimum lies well above the exact $-1.1745$ Ha. "FCI" here means exact only within the two-orbital STO-3G space.
- **The RHF dissociation error depends on where the scan stops.** RHF has not reached its limit at 5 Å: its excess over $2E_H$ is +0.3341 Ha at 5 Å, +0.3623 at 20 bohr and +0.3868 at 1000 bohr, approaching +0.3873 Ha. The reported RHF error and the RHF depth measured from its own end point are values at 5 Å, not limits, and the RHF well depth in the table is measured against an energy RHF never reaches.
- **The small FCI residual at 5 Å is not an error of the method.** It is real singlet binding left at that distance plus a $\sim 10^{-8}$ Ha basis-normalisation offset.
- **Ionic weights are definition-dependent at short range.** The Chirgwin–Coulson peak near 0.86 Å is a property of the definition, not a physically "most ionic" distance. Only the trends are robust: RHF stays at 0.5, and FCI falls to zero.
- **Model.** Non-relativistic Hamiltonian with clamped nuclei. The valence-bond picture holds only in this minimal basis, and the integral code handles s-type functions only.
- **Not tested:** larger basis sets, unrestricted HF, excited singlets, nuclear motion and zero-point energy, relativistic effects, and bond lengths outside 0.4–5.0 Å (near zero separation the overlap approaches 1 and the orthogonalisation becomes ill-conditioned). The literature anchors come from standard references rather than a cited research note. The full list is in `LIMITS.md`.

## 6. Reproducing it

```
python3 -m venv .venv && .venv/bin/pip install numpy scipy matplotlib pytest
.venv/bin/python run.py            # 461 points, a few seconds
.venv/bin/python -m pytest tests   # unit tests and the Sceptic's checks
```

The code is in `src/h2qc/` (`integrals.py`, `hf.py`, `fci.py`, `vb.py`, `scan.py`). The outputs are `results/summary.json`, `results/scan.csv` and `results/reference_point.json`. The agreed plan is in `PLAN.md`.
