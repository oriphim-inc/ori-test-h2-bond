
## Sceptic checks (2026-11-14T00:47, after 0 runs)

I wrote 40 pytest checks in tests/checks/test_checks.py. They skip until results/ exist. They rest on an independent reference in tests/checks/_ref.py: my own integrals and the closed-form solution that the symmetry of a two-function basis allows. I validated that reference myself: it reproduces Szabo & Ostlund at 1.4 bohr (E_RHF -1.11671, E_FCI -1.13728) and gives RHF R_e 0.7122 A, FCI R_e 0.7349 A, and E_H -0.466582. The literature R_e values (0.712 and 0.735 A) and the exact H2 bound (-1.174475) come from standard references, not from notes/research.md, which is empty.

- check: Szabo & Ostlund anchors at R=1.4 bohr (reference_point.json): S12 0.6593, T11 0.7600, T12 0.2365, V11_total -1.8804, V12_total -1.1949, (11|11) 0.7746, (11|22) 0.5697, (21|11) 0.4441, (21|21) 0.2970, E_RHF -1.1167, E_FCI -1.1373, each within 1e-4 (S&O print 4 dp); eps = -0.578, 0.670 within 2e-3
- check: E_H(STO-3G) = -0.4666 +/- 1e-4 and equal to an independent evaluation to 1e-10; E_2H = 2 E_H
- check: Grid is 461 points, 0.4-5.0 A, step 0.01; R_bohr*0.529177210903 = R_ang (rel 1e-12)
- check: Whole curves vs independent reference (tests/checks/_ref.py: own vectorised integrals + symmetry-fixed closed-form RHF, 2x2 singlet FCI, triplet; no h2qc import): S_ab 1e-10, E_RHF/E_FCI/E_triplet 1e-9 Ha, eps_g/eps_u 1e-8, c_gg/c_uu 1e-8, FCI ionic weights (CC and Loewdin) 1e-8, at every grid point
- check: SCF converged at every R, iterations < max_iter 200
- check: Variational: E_FCI <= E_RHF + 1e-12 everywhere; E_FCI > exact H2 BO minimum -1.174475 everywhere; E_FCI(5 A) > -1.0
- check: Singlet structure: |c_gu|,|c_ug| < 1e-10; c_gg^2+c_uu^2 = 1 (1e-10); c_gg > 0, c_uu < 0
- check: Triplet: E_triplet >= E_FCI everywhere, > 2E_H - 1e-9 (unbound), strictly decreasing with R; |E_triplet - E_FCI| < 1e-6 at 5 A; |c_gg + c_uu| < 1e-3 at 5 A
- check: FCI dissociation error |E_FCI(5 A) - 2E_H| < 1e-6 Ha
- check: RHF dissociation error = (aa|aa)/2 - 1/(2R) within 1e-3 (O(S) corrections ~3e-4) and > 0.3 Ha
- check: R_e: matches an independent Brent minimisation of the reference to 2e-5 A and E_min to 1e-9; literature HF/STO-3G 0.712 A and FCI/STO-3G 0.735 A within 1e-3 A; RHF E_min -1.1175 +/- 1e-4; E_min <= grid minimum; grid minimum within 0.005 A of R_e; R_e_bohr consistent
- check: Bookkeeping: E_at_Rmax = last row; D_e_vs_2H = E_2H - E_min; D_e_vs_own_end = E_at_Rmax - E_min; dissociation_error = E_at_Rmax - E_2H; eV conversion uses CODATA 27.211386245988
- check: Grid independence: run_scan(47) reproduces R_e to 2e-5 A and E_min to 1e-10 for both methods
- check: Weights: ion+cov = 1 (1e-10) for both definitions and methods, all in [0,1]; RHF ionic weight exactly 0.5 at every R (1e-10) by symmetry; FCI ionic <= 0.5; Loewdin FCI ionic weight monotone non-increasing over the full grid; CC FCI ionic weight monotone for R >= 1.0 A (it peaks near 0.86 A, which is physical for CC at large overlap, not a bug); both < 1e-5 at 5 A
- check: No hard-coded targets: none of 1.1167, 1.1373, 1.1175, 0.712x, 0.734x/0.735x, 0.2157, 0.388 appear as literals in src/ or run.py

## Sceptic review (2026-11-13, after the 461-point run)

**Verdict: pass. 48/48 checks pass after correcting the Sceptic's own reference. No tolerance was loosened.**

### The 7 failures were caused by my reference, not by h2qc
I tested the author's diagnosis directly. With the 8-digit constants, `_ref` gives
S_aa = 0.9999999908898. With the old normalisation 1/sqrt(2(1±s)), C^T S C − I has
diagonal elements −4.9e-9 (g) and −7.2e-8 (u). The MOs were therefore not normalised,
which matches the size of the energy differences (7e-9 to 9e-8).

The fix in tests/checks/_ref.py:
- The MOs are now normalised with sqrt(2(S_aa ± s)).
- The VB structure norms and overlap now use S_aa: N² = 2(S_aa² + s²) and
  o = 2 s S_aa / (S_aa² + s²).

After the fix:
- E_RHF, E_FCI, E_triplet and eps_u agree with results/scan.csv to at most 8e-15 Ha.
- All 48 tests pass.

### Code review
- Integrals: they match the reference to 1e-16 (author's claim). The S&O anchors pass.
- RHF: symmetric orthogonalisation with the true S matrix. It converges in 2 iterations,
  because symmetry fixes the MOs.
- FCI: Slater–Condon rules over the M_s = 0 determinants, then separate singlet and triplet
  blocks. This is correct.
- Definitions of D_e and the dissociation error: they match PLAN.md, and the bookkeeping
  tests pass.
- No hard-coded targets.

### Notes and limits (details in LIMITS.md)
- **h_atom_energy and vb.py treat S_aa as 1.** The effect is 1e-8 or less. E.g.
  2 E_H is 8.5e-9 Ha too high: an FCI scan to R = 1000 bohr levels off at 2 E_H − 8.5e-9.
  Note only.
- **RHF at 5 Å has not reached its limit.** It is +0.3341 Ha above 2 E_H, while the limit as
  R → ∞ is +0.3873. The RHF "dissociation error" and "D_e vs own end" therefore depend on
  R_max. This is a limit.
- **The FCI dissociation error of −7.1e-8 is real residual singlet binding at 5 Å plus the
  S_aa offset**, not a failure of the method.
- **The CC and Löwdin ionic weights diverge at short R** (0.06 vs 0.44 at 0.4 Å), because
  the two structures become nearly linearly dependent. Only the trends are robust.
- **Plot cosmetics in pes.png:**
  - The FCI annotation is clipped at the bottom.
  - The upper y-limit of +0.1 leaves a lot of empty space.

## Sceptic review (2026-11-14T00:51, after 2 runs)

Pass: 48 of 48 checks pass. The 7 earlier failures were in the Sceptic's own reference: _ref normalised the MOs assuming <a|a> = 1, but the 8-digit constants give S_aa = 1 − 9.1e-9, so its MOs were off-normalised by 4.9e-9 and 7.2e-8. With the reference corrected to use S_aa, it agrees with h2qc to below 1e-14. No tolerance was loosened. Code and physics are sound. The main limits: RHF has not reached its limit at 5 Å, so its dissociation numbers depend on R_max; D_e is the STO-3G value, not the physical one; and the VB weights depend on the definition at short R.

- note: The 7 failures were caused by the Sceptic's reference, tests/checks/_ref.py. It normalised the MOs with sqrt(2(1±s)), but S_aa = 0.9999999908898, so C^T S C − I had diagonal elements −4.9e-9 and −7.2e-8. Fixed by normalising with sqrt(2(S_aa±s)) and using S_aa in the VB structure norms. Now E_RHF, E_FCI, E_triplet and eps_u agree with scan.csv to at most 8e-15. 48 of 48 pass, with tolerances unchanged.
- note: h_atom_energy (E_H = T_aa + V_aa) and src/h2qc/vb.py (N^2 = 2(1+S^2)) assume S_aa = 1. This makes 2E_H about 8.5e-9 Ha too high (an FCI scan to R = 1000 bohr levels off at 2E_H − 8.5e-9), and the weights shift by 1e-8 or less. Harmless at the precision reported; worth a one-line fix or a comment.
- limit: RHF has not reached its dissociation limit at 5 Å. E_RHF − 2E_H is +0.3341 Ha at 5 Å, +0.3623 at 20 bohr and +0.3868 at 1000 bohr, approaching 2E_H + ½(aa|aa) = +0.3873. So the RHF 'dissociation error' and 'D_e vs own end' (0.518 Ha) are values at R = 5 Å, not asymptotic values.
- limit: D_e is measured against 2E_H(STO-3G). D_e(FCI) = 5.55 eV, compared with the exact H2 value of about 4.75 eV, and E_min(FCI) = −1.1373 compared with the exact −1.1745 Ha. These are minimal-basis numbers, not physical ones.
- limit: The FCI dissociation error at 5 Å (−7.1e-8 Ha) is mostly real singlet binding left at 9.45 bohr (the singlet–triplet gap is 8e-8) plus the 8.5e-9 S_aa offset. It is not an error of the method.
- limit: The ionic/covalent weights depend on the definition at short R, where the structures are nearly linearly dependent: at 0.4 Å the ionic weight is 0.06 by Chirgwin–Coulson and 0.44 by Löwdin, and the CC peak near 0.86 Å is a property of the definition. Only these trends are robust: FCI ionic weight → 0 at large R, RHF = 0.5 at every R.
- limit: Assumptions: the 8-digit STO-3G constants (ζ = 1.24) from PLAN.md, and literature anchors (S&O at 1.4 bohr; R_e of 0.712 and 0.735 Å from CCCBDB; exact −1.174475 Ha) that come from standard references, not from a source in notes/research.md, which is empty. Not tested: other bases, UHF, R outside 0.4–5.0 Å, nuclear motion and zero-point energy.
- note: pes.png: the FCI R_e/E annotation is clipped at the bottom of the axes, and the upper y-limit of +0.1 wastes vertical space. Cosmetic only.

## Sceptic draft (2026-11-24T03:11, after 2 runs)

Draft is mostly accurate: every placeholder points at the right key or row (scan rows 0/34/110/160/260/460 = 0.4/0.74/1.5/2.0/3.0/5.0 Å), and the hand-typed numbers check out (0.8 eV, 4.75 eV, -1.1745, +0.387 = 0.7746/2, 0.86 Å CC peak, 0.06 vs 0.44, 41/48, ~9 eV, ~0.02 Å, 1e-7 singlet-triplet gap). The PDF is 4 pages. Two claims overstate the evidence and must be fixed: (1) the two weight definitions are said to agree beyond 1.5 Å, but they differ by about 30-40 % there; (2) 1e-14 agreement is claimed for 'every curve', but it holds only for the energies (the weights agree to about 1e-8). There are also small wording fixes and a few missing limits.

- suspect: Sec. 3, 'Beyond about 1.5 Å the two weight definitions agree.' is false according to scan.csv: at 1.5 Å CC 0.121 vs Löwdin 0.168, at 2.0 Å 0.034 vs 0.047, at 3.0 Å 0.0011 vs 0.0014 (a 30-40 % relative gap throughout). Replace with: 'The two definitions differ at every distance (0.121 vs 0.168 at 1.5 Å, 0.034 vs 0.047 at 2.0 Å), but both fall towards zero together; the gap is largest at short distance (0.06 vs 0.44 at 0.4 Å), because the two structures become nearly identical when the atoms overlap strongly.'
- suspect: Abstract, 'An independent reference implementation reproduces every curve to about 10^-14 Ha', and Sec. 4, 'the curves agree to about 10^-14 Ha', overstate the agreement. notes/review.md reports at most 8e-15 Ha only for E_RHF, E_FCI, E_triplet and eps_u. The ionic weights are checked at a 1e-8 tolerance, and vb.py assumes S_aa = 1, so they differ at up to about 1e-8. Change to: 'reproduces every energy curve to about 10^-14 Ha (and the ionic weights to about 10^-8)'.
- note: Abstract, 'the FCI curve lies within -7e-08 Ha of two isolated ... atoms', puts a negative number after 'within'. Write: 'lies 7e-08 Ha below two isolated STO-3G hydrogen atoms' (use the magnitude, or reword so the sign is meaningful).
- note: Abstract, 'The FCI ionic weight falls from about 0.2 near equilibrium' holds only for Chirgwin–Coulson (0.216 at 0.74 Å); Löwdin gives 0.39. Write 'The FCI ionic weight (Chirgwin–Coulson) falls from about 0.2 ...'.
- note: Sec. 4, '40 pass/fail checks fixed in advance' followed by '41 of 48 tests' reads as inconsistent. Explain: '40 pass/fail checks fixed in advance, run together with 8 unit tests (48 tests in all)'.
- note: Sec. 4, 'reproduces every integral and energy printed by Szabo & Ostlund to all four decimals': the orbital energies were checked only to 2e-3 (S&O print them to 3 decimals). Write 'every integral and total energy ... to the four decimals printed'.
- note: Sec. 4 'remaining flags' and the Fig. 1 caption say the labels 'overlap the axis'. notes/review.md says the FCI annotation is clipped at the bottom of the axes and the upper y-limit wastes space. Use the review's wording.
- limit: Sec. 5 leaves out some LIMITS.md items. Add one line: 'The Hamiltonian is non-relativistic with clamped nuclei; the valence-bond picture holds only in this minimal basis; near R -> 0 the overlap approaches 1 and the orthogonalisation becomes ill-conditioned (untested below 0.4 Å); the integral code handles s-type functions only.' Also add relativistic effects to the 'Not tested' list.
- limit: Sec. 5 RHF bullet: add the numbers from LIMITS.md so readers can see the drift: the excess over 2E_H is +0.3341 Ha at 5 Å, +0.3623 Ha at 20 bohr and +0.3868 Ha at 1000 bohr, approaching +0.3873 Ha. Optionally note that the RHF well depth in Table 1 (5.02 eV vs 2E_H) is measured against an energy RHF never reaches.

## Sceptic draft (2026-11-24T03:12, after 2 runs)

Revised draft re-checked against results/summary.json, scan.csv, reference_point.json, LIMITS.md and notes/review.md. Most requested fixes are in and correct: 1e-14 is now limited to the energy curves with weights at ~1e-8; the abstract wording, the Chirgwin–Coulson qualifier, '40 checks + 8 unit tests = 48', 'to the four decimals printed', the clipped-annotation wording and all the added LIMITS items (non-relativistic/clamped nuclei, VB only in minimal basis, s-type only, ill-conditioning near R->0, relativistic, RHF drift +0.3341/+0.3623/+0.3868 -> +0.3873) are present and match LIMITS.md and notes/review.md. Placeholders resolve in paper.md. The PDF is 4 pages. One problem remains, and it came from my own earlier review: 'by 30–40 % beyond 1.5 Å' is true only out to about 3 Å. From scan.csv the Löwdin/CC ratio is 1.39 at 1.5 Å, 1.38 at 2.0, 1.31 at 3.0, 1.19 at 4.0 and 1.09 at 5.0 Å. One sentence must be corrected. Everything else passes.

- suspect: Sec. 3, sentence after the ionic-weight table: 'The two definitions differ at every distance (by 30–40 % beyond 1.5 Å) but fall towards zero together.' This overstates the gap at long range. From results/scan.csv the Löwdin/Chirgwin–Coulson ratio is 1.39 (1.5 Å), 1.38 (2.0), 1.37 (2.5), 1.31 (3.0), 1.19 (4.0) and 1.09 (5.0 Å). The table's own 5.0 Å row (1.2e-7 vs 1.4e-7) contradicts '30–40 %'. My previous review introduced this wording. Replace it with: 'The two definitions differ at every distance: between 1.5 and 3 Å the Löwdin weight is 30–40 % larger, and the gap narrows to about 10 % at 5 Å as both fall towards zero.'
- note: Verified correct: abstract numbers (R_e 0.7122/0.7349 Å, E_min, CC ionic weight about 0.2 near R_e (0.216 at 0.74 Å), 1.2e-7 at 5 Å, FCI within 1e-7 Ha of 2E_H (-7.1e-8)); about 9 eV (0.334 Ha); 0.8 eV over 4.75 eV (5.55 eV); CC peak at 0.86 Å; 0.06 vs 0.44 at 0.4 Å; RHF drift numbers match LIMITS.md lines 44-51 and notes/review.md; the PDF has /Count 4, so it is 4 pages.

## Sceptic draft (2026-11-24T03:13, after 2 runs)

Final draft check: pass. The corrected sentence in Sec. 3 matches results/scan.csv. From the CSV, the Löwdin/Chirgwin–Coulson FCI ionic-weight ratio over 1.5–3.0 Å runs from 1.388 down to 1.311, so '30–40 % larger' is correct. At 5.0 Å it is 1.094 ('about 10 %'). Its minimum over all 461 points is 1.094 > 1, so 'differ at every distance' holds. At 0.4 Å the weights are 0.0588 vs 0.4407 ('0.06 vs 0.44'). This fixes the one issue left from the previous draft review. paper.md has no unresolved placeholders, and the rendered sentence and table row (0.40 | 0.059 | 0.441 | 0.500) agree with the source. paper.pdf was rebuilt at the same time as source.md, and its Pages tree has /Count 4, so it is still 4 pages. paper/ is untracked, so I could not diff it against the version I last reviewed. I am relying on Ori's statement that this was the only change, together with a spot check of the rendered output.

- note: Sec. 3 sentence checked against scan.csv: L/CC ratio 1.388 (1.5 Å), 1.383 (2.0), 1.311 (3.0), 1.094 (5.0); minimum over all 461 points 1.094; 0.0588 vs 0.4407 at 0.4 Å. All of the sentence's claims hold.
- note: paper/ is untracked in git, so I could not confirm 'nothing else changed' by diff. I checked that placeholders resolve (0 left in paper.md) and that the PDF page count is 4 (/Pages /Count 4).
- note: The causal clause ('because the two structures become nearly identical when the atoms overlap strongly') explains the result and does not claim a number. It is acceptable as written.
