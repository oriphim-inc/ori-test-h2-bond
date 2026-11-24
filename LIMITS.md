# LIMITS: what this project does NOT show

The project computes H2 in the STO-3G basis: RHF and full CI over the two molecular orbitals,
for R = 0.4 to 5.0 Å. All of the numbers are exact *for that model*. In a separate
reference implementation (tests/checks/_ref.py) the same quantities agree to about 1e-14 Ha.
That agreement shows the code solves the model correctly. It says nothing about how close
the model is to real H2.

## Assumptions (taken as given, not tested)
- **Basis constants.** We used the 8-digit STO-3G exponents and contraction coefficients
  listed in PLAN.md (Hehre, Stewart & Pople 1969, as tabulated by Szabo & Ostlund). Because
  they are rounded, the contracted 1s is not exactly normalised: <a|a> = 0.9999999909
  (1 - 9.1e-9). The code takes this into account wherever it solves the Roothaan equations
  with the true S matrix. In three places it treats <a|a> as 1, which shifts results by
  1e-8 or less:
  - E_H = T_aa + V_aa, without dividing by S_aa. This puts 2 E_H about 8.5e-9 Ha too high.
  - The VB structure norms in src/h2qc/vb.py.
  - The definitions in PLAN.md.

  Basis constants with other roundings, or ζ ≠ 1.24, change the results at the 1e-8 level
  or more.
- **Literature anchors.** The pass criteria use these external values:
  - Szabo & Ostlund (1989) integrals and energies at R = 1.4 bohr, printed to 4 decimals.
  - R_e(HF/STO-3G) = 0.712 Å and R_e(FCI/STO-3G) = 0.735 Å (as tabulated in the NIST
    CCCBDB).
  - E_H(STO-3G) = −0.4666 Ha.
  - The exact H2 minimum, −1.174475 Ha (used only as a lower bound).

  The CCCBDB values and the exact H2 minimum come from standard references we know, not from
  a source cited in notes/research.md (that file is empty). They are agreement checks only.
  The code never uses them.
- **Model.** The Hamiltonian is non-relativistic and electronic, with clamped nuclei
  (Born–Oppenheimer). Energies are in hartree. The unit conversions use CODATA 2018 values.

## What the numbers do NOT mean
- **D_e is not the physical well depth of H2.** D_e(FCI) = 0.2041 Ha (5.55 eV) is measured
  from 2 E_H(STO-3G) = −0.93316 Ha. The minimal basis gives an atom energy 0.033 Ha above
  the exact −0.5 Ha, and the molecule gains more than the atoms from that basis error. The
  experimental and exact D_e of H2 is about 4.75 eV, so the STO-3G value is roughly 0.8 eV too
  deep. R_e and E_min are STO-3G values. E_min(FCI) = −1.1373 Ha, compared with the exact
  −1.1745 Ha.
- **"RHF well depth relative to its own end point" (0.518 Ha) has no physical meaning.**
  The RHF curve has not reached its limit at 5.0 Å. It keeps rising, like
  2 E_H + ½(aa|aa) − 1/(2R), towards 2 E_H + 0.3873 Ha. A single targeted check in the
  reference gave the following excess over 2 E_H:

  | R | Excess over 2 E_H |
  |---|---|
  | 5 Å | +0.3341 Ha |
  | 20 bohr | +0.3623 Ha |
  | 1000 bohr | +0.3868 Ha |

  So the reported RHF dissociation error (0.3341 Ha) and the RHF "D_e vs own end" are values
  *at R = 5 Å*, not dissociation limits. A larger R_max would change both of them.
- **The FCI dissociation error at 5 Å (−7.1e-8 Ha) is not an error of the method.** About
  −6.3e-8 of it is real singlet binding that is left at 9.45 bohr: the singlet lies 8e-8 Ha
  below the triplet. The other −8.5e-9 comes from the S_aa normalisation of E_H described
  above. As R → ∞, E_FCI − 2 E_H → −8.5e-9 (the normalisation offset only).
- **Ionic and covalent weights depend on how they are defined, especially at short R.**
  When the two VB structures overlap strongly (S → 1 as R → 0) they become almost linearly
  dependent, and the Chirgwin–Coulson and Löwdin weights differ widely. At 0.4 Å the FCI
  ionic weight is 0.06 by Chirgwin–Coulson and 0.44 by Löwdin. The CC weight peaks near
  0.86 Å. That peak is a property of the definition, not a physical "most ionic" distance.
  Only the qualitative trends are independent of the definition:
  - The FCI ionic weight goes to 0 at large R (below 2e-7 at 5 Å).
  - The RHF weight stays at exactly 0.5 at every R.

  The VB picture holds only in this minimal basis. In a larger basis the AO structures are
  no longer the full wavefunction.
- **"FCI" means the full CI in this 2-orbital basis**, which is exact correlation for STO-3G
  only. It does not cover dynamic correlation from polarisation or diffuse functions, which
  any realistic basis would need.

## Not modelled / not tested
- Basis-set convergence: no larger basis was tried, so there is no basis superposition
  study and no estimate of the basis error beyond the comparison with the exact energies
  above.
- Unrestricted HF (UHF), which would dissociate correctly after the Coulson–Fischer point.
- Excited singlet states (only the lowest singlet and the lowest triplet are reported).
- Nuclear motion: no zero-point energy, so the bond energy D_0 is not computed. No
  vibrational levels.
- Relativistic effects, the non-adiabatic coupling between electronic and nuclear motion
  (Born–Oppenheimer breakdown), and external fields.
- R below 0.4 Å or above 5.0 Å. Near R → 0 the overlap S → 1 and S^{-1/2} becomes
  ill-conditioned. This was not tested.
- Molecules other than H2, and other basis sets. The integral code handles s-type
  primitives only.
