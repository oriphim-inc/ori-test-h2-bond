# PLAN — minimal quantum chemistry for H2 / STO-3G

## Question
Build from scratch (numpy/scipy only, no quantum-chemistry packages) a code that, for H2 in
STO-3G, computes analytic Gaussian integrals, runs RHF (SCF loop) and full CI in the
2-electron / 2-orbital space, and scans R = 0.4 … 5.0 Å. Deliver:
1. both potential energy curves (RHF, FCI), plotted;
2. each curve's equilibrium bond length R_e and well depth D_e;
3. the large-R limit of each curve vs two isolated STO-3G H atoms (2 × −0.4666 Ha);
4. the ionic vs covalent weights of the FCI ground state as R grows.

## Physics / system (agreed)
- Non-relativistic, clamped-nuclei electronic Hamiltonian + nuclear repulsion 1/R (atomic units).
- Basis: STO-3G for H, one contracted 1s per atom, Slater exponent ζ = 1.24.
  Primitive exponents (ζ=1.24) α = 3.42525091, 0.62391373, 0.16885540;
  contraction coefficients d = 0.15432897, 0.53532814, 0.44463454 (for normalized primitives).
  Source: Hehre, Stewart & Pople, J. Chem. Phys. 51, 2657 (1969); tabulated in
  Szabo & Ostlund, *Modern Quantum Chemistry* (1989), §3.6.
- Integrals (s-type primitives, Szabo & Ostlund App. A): overlap, kinetic, nuclear
  attraction (Boys F0), two-electron repulsion (chemists' notation (μν|λσ)), all analytic.
- RHF: closed shell, symmetric orthogonalization, core-Hamiltonian guess, SCF loop to
  |ΔE| < 1e-12 Ha and RMS ΔD < 1e-10 (numerical settings, tunable).
- FCI: all determinants of 2 electrons in the 4 spin orbitals built from the 2 RHF MOs
  (6 determinants; the Hamiltonian is built in the M_s = 0 block, 4 determinants,
  exact since H conserves M_s). Ground state = lowest *singlet* root (singlet ⇔ symmetric
  spatial coefficient matrix); this matters at large R where the triplet becomes degenerate.
- Unit conversion: 1 bohr = 0.529177210903 Å; 1 Ha = 27.211386245988 eV = 627.5095 kcal/mol
  (CODATA 2018).

## Definitions of the deliverables
- R_e, E_min: minimum of each curve, located by bounded Brent minimisation of the
  energy function itself (not interpolation), bracketed by the grid minimum.
- Well depth: D_e = E_ref − E_min with E_ref = 2 E_H(STO-3G), the separated-atom energy
  computed by the code (one electron in one contracted 1s: E_H = T_11 + V_11). Also reported:
  E(5.0 Å) − E_min (the depth relative to the curve's own end point).
- Dissociation-limit error: E(R = 5.0 Å) − 2 E_H for each method; and the code's E_H vs
  the stated −0.4666 Ha.
- Ionic/covalent weights: write the FCI singlet spatial wavefunction in the AO basis,
  Ψ = C_cov (a b + b a) + C_ion (a a + b b) (a, b = the atomic 1s functions), normalise the
  two (non-orthogonal) VB structures, and report **Chirgwin–Coulson** weights
  (w_i = c_i Σ_j c_j S_ij, sum to 1). Löwdin-orthogonalised weights are reported alongside
  because the weight definition for non-orthogonal structures is not unique. The same
  decomposition of the RHF wavefunction is shown for contrast.

## Scan
- R grid: 0.4 … 5.0 Å, step 0.01 Å (461 points) by default; `python run.py --n N` changes
  the number of points. Nothing is scaled down; the whole problem fits on a laptop.

## What would count as the answer (sanity anchors, Szabo & Ostlund 1989)
At R = 1.4 bohr: S_12 ≈ 0.6593, T_11 ≈ 0.7600, T_12 ≈ 0.2365, (11|11) ≈ 0.7746,
(11|22) ≈ 0.5697, (21|11) ≈ 0.4441, (21|21) ≈ 0.2970, E_RHF ≈ −1.1167 Ha,
E_FCI ≈ −1.1373 Ha. Expected qualitative outcome: RHF dissociates far too high
(ionic contamination fixed at ~50%), FCI dissociates to 2 E_H, ionic weight → 0.
