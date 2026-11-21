# Shifted q-metallic Hankel determinants

This package proves Han–Pedon's Conjecture C (arXiv:2610.36132v1, p. 6).
For every integer n ≥ 1 and shift ℓ ≥ n + 3, each residue subsequence is a
rational polynomial starting at m = 0, integer-valued on all integers.
The valid period is 2n(n+1) for even n and 4n(n+1) for odd n.
Extending Han–Pedon's lower-shift results (2026), the note also proves a degree
bound, rational generating functions, linear recurrences and periodicity from
the first term modulo every positive integer at every shift. At higher shifts,
a valid modular period is rho*M*lcm(1,...,max(1,d)), where d bounds the residue
degrees; rho*M*d! is also valid.

## Contents

- `note.tex`: complete mathematical argument and bibliography.
- `metadata.json`: authors, license, statement, and mathematical coverage.
- `verify.py`: exact-integer finite checks; Python 3.9 or later, standard library only.
- `checks.json`: the resulting finite values and scope.

The readable document is `note.pdf`.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively use a standard LaTeX distribution:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is included directly in the TeX document.

## Finite checks

```sh
python3 verify.py --output reproduced-checks.json
cmp checks.json reproduced-checks.json
```

The checks use fraction-free integer elimination, coefficient recursions,
continued-fraction recursions, complementary minors, and finite differences.
They cover cyclic identities for 1 ≤ n ≤ 40; eight normal determinants and
648 orthogonality equations; 252 complementary determinants for 1 ≤ n ≤ 4,
a = 2, 3, 4, 0 ≤ j ≤ 20; and 120 endpoint coefficient sequences for
1 ≤ n ≤ 6, 0 ≤ h ≤ 4, 1 ≤ K ≤ 7. They also give the first six terms of
all eight residue subsequences at n = 1, ℓ = 4, and 8100 exact modular
Newton-basis identities for 0 ≤ d ≤ 8, 1 ≤ M ≤ 20, -4 ≤ m ≤ 4 and 0 ≤ v ≤ d.

The function `golden_polynomial_identities` gives 220 exact complementary-minor
identities for the golden shift-4 laws. Each residue is evaluated at all
d_r+1 points m = 0,...,d_r, using the proved degree bounds
[50, 37, 26, 17, 10, 5, 2, 65]. Thus the degree theorem and these finite
identities prove the affine formulas for every m ≥ 0; this conclusion does
not follow from the 48 illustrative direct determinants alone. The pattern
is credited to Ovsienko–Pedon (2023; published 2025), equation (1.18).
The same laws evaluate the third-shift RNA determinants and the diagonal
T(j+2,j) of OEIS A123634. The coefficient-sign identity has a written proof;
48 direct RNA determinants give additional exact checks. The earlier number-wall
data are credited to Michael Somos (2006), with the distinct earlier diagonal
formulas recorded in that entry in 2016.

## Mathematical coverage

The exact formal statement is
`D5.S3.Combinatorics.MetallicHankel.MetallicHankelResiduePolynomials.result`
at the pinned source:

https://github.com/the-omega-institute/trureturing/blob/8c7bbde1c8750c7c307bdc6bc62f44c94b7b2c7f/D5/S3/Combinatorics/MetallicHankel/MetallicHankelResiduePolynomials.lean

It proves the residue identities for every allowed n, ℓ, r, and m, and
integer-valuedness at every integer argument. The degree bound and the
explicit composite-modulus period have written proofs in the note.
Appendix A links all 22 pinned Lean modules and states each declaration's
scope, distinguishing definitions from theorems. The rational generating
function, linear recurrence and growth dichotomy also have written proofs;
the golden corollary combines the degree bound and the 220 finite identities.
The source-item coverage table distinguishes the previously proved bounded
shift from the higher-shift theorem, the earlier unboundedness question,
and the bilinear-recurrence searches.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

CC BY 4.0 (Creative Commons Attribution 4.0 International).
