This repository contains the source of the **Q7 Cosmochrony paper**  
*Toward Three Spatial Directions from a Supplied Spin-One Carrier:
Requirements for an Equivariant Bridge*.

This work addresses a central structural question of the Cosmochrony
spectral admissibility programme:

> Are the three admissible directions arising from the neutral sector of the supplied
> spinor carrier (O23, conditional) identical to the three spatial directions of the
> emergent Lorentzian geometry?

The paper states what an identification would require, and proves the requirements that
follow from representation theory on a supplied carrier. It does not identify the two
structures: the required bridge is a **missing identification**, not a proved obstruction.

## Conceptual Overview

1. **Separation of the two "3"s**  
   Two mechanisms produce the integer 3:
    - the geometric 3 from the Carnot structure of $\mathrm{Heis}_3$ (Q5b, whose co-metric
      extraction is conditional on its unestablished hypothesis [H-L], with a degenerate
      central slot at principal-symbol order),
    - the admissible 3 from the neutral sector of a supplied spinor carrier
      (O23 Theorem 3.1, conditional on that carrier) and from the projection space
      $H_{\mathrm{eff}} = \mathbb{C}^3$, which O28 supplies as a parameter consistent with that
      selection rule, and whose covariance rank three O28 itself calls a finite-data value of its
      protocol and not an invariant (it also reports a participation-ratio value $8/3$).

   These structures are *a priori distinct* and cannot be identified directly.

2. **Algebraic obstruction and representation-level bridge**  
   The Lie algebras $\mathfrak{su}(2)$ and $\mathfrak{heis}_3$ are
   non-isomorphic (semisimple vs nilpotent).  
   Any identification must therefore occur at the level of **representations**.

3. **The identification hypothesis, isolated**  
   That the measured space carries the spin-1 module structure,
   $H_{\mathrm{eff}} \simeq \mathrm{Sym}^2(V_\rho)$, is **Hypothesis [ID]** of the paper.
   No source supplies it: O29 states that no calculation there identifies these objects,
   and $d_\rho = 2$ is O26's minimality selection.

4. **Requirements under [ID]**  
   Up to isomorphism, $\mathrm{Sym}^2(V_\rho)$ is the unique irreducible three-dimensional
   complex module of $\mathfrak{su}(2)$, isomorphic to the complexified adjoint. Hence any
   equivariant bridge is **either zero or an isomorphism, unique up to one non-zero complex
   scalar**; reality and positivity of that scalar are a further requirement. The Casimir
   acts as $2 \cdot \mathrm{Id}$, and invariance under the Cartan generator forces isotropy
   *within the horizontal plane* only: it does not force $A_H = A_z$.

5. **Reduction to a symbol criterion**  
   An identification would have to satisfy
   $\sigma_2(L_{\mathrm{eff}})\big|_{\mathrm{Sym}^2(V_\rho)} = A_H(k_X^2 + k_Y^2) + A_z k_Z^2$,
   a test on the continuum operator of Q5b. No such computation exists in the corpus.

## What the O25 Checkpoints Measure

Version 2.0 re-types the numerical section. The three stored basis rows are **pure Fourier
modes** of indices $\{0, -k, +k\}$ (purity 1.000000 on all twelve tested pairs), so every
entry of the compression $\tilde L = B_{\mathrm{eff}} L_{\mathrm{Weil}} B_{\mathrm{eff}}^\dagger$
is a Rayleigh quotient in closed form:

$\langle e_m, L_{\mathrm{Weil}} e_m \rangle = 4 - 2\cos(2\pi m/q).$

- The central value is $4 - 2 - 0 = 2$ exactly, for every $q$ and $c$: the additive $4I$ term
  of $L_{\mathrm{Weil}} = 4I - W_a - W_a^\dagger - W_b - W_b^\dagger$ less the shift eigenvalue
  at the flat mode. It does agree numerically with the $\mathfrak{su}(2)$ Casimir eigenvalue on a
  spin-1 module, but the agreement is not normalisation-invariant — $L/2$ gives 1, $L - 2I$ gives
  0 — and therefore carries no evidential weight.
- The horizontal value is $4 - 2\cos(2\pi k/q)$ for the selected index $k$, so
  $A_H - 2 = 4\sin^2(\pi k/q)$ **identically**. Vanishing of that gap is equivalent to
  $\mathrm{dist}(k/q, \mathbb{Z}) \to 0$, a property of the pipeline's mode-selection rule,
  which no paper of the corpus derives.
- An off-diagonal entry of modulus 1 appears exactly when two stored indices differ by
  $\pm c \bmod q$: this accounts for the two pairs, $(61,3)$ and $(151,5)$, that version 1.3
  reported as finite-size and sector anomalies.

**The closed form describes the stored-basis compression, not directly the published entries.**
The pipeline conjugates that compression by the covariance eigenbasis, and that basis is fixed by
the data only where the covariance spectrum is non-degenerate. Measured: the gap between the second
and third normalised eigenvalues of $C_c$ is at machine precision ($10^{-16}$ to
$2\times10^{-15}$) at nine of the twelve pairs, and resolved at three ($6.2\times10^{-3}$,
$1.2\times10^{-2}$, $1.6\times10^{-2}$). At the nine degenerate pairs the rotation of that plane is arbitrary, but it has
no effect wherever the compression restricted to the plane is **scalar**, which holds at eleven of
the twelve pairs: a rotation of a plane on which an operator acts as a scalar returns the same
scalar, so no solver can split those entries. Exactly one pair has both an arbitrary rotation and a
non-scalar block — $(151,5)$, where the indices $\pm 73$ differ by $-c$ — and there the split
$6.9527 / 5.0257$ and the off-diagonal $0.2677$ are solver-dependent; the ratio of that
off-diagonal to the mean diagonal is the published $\rho_{XY} = 0.057$, in a column captioned as
the isotropy ratio, which evaluates to $0.2772$ under the caption's own formula and to $0.3218$ as
the script computes it. The invariant content is the spectrum of the compression, which the closed
form gives outright at the ten pairs with no index difference equal to $\pm c \bmod q$ and which
is shifted at the two that have one. §6.3 of the paper works this out, and records one published number,
$\rho_{XY} = 0.005$ at $(101,3)$, that no versioned quantity reproduces.

Reproducible diagnostic: `code/q7_mode_diagnostic.py` prints Fourier indices, purities, measured
and analytic values, residuals, the cross-term check, the covariance spectrum with its degeneracy
gap, and the rotated matrix with a warning where the entries are solver-dependent. Largest
residual: $3.6 \times 10^{-15}$.

It runs from the repository alone: `code/data/q7_stage_inputs.npz` (37 kB) ships the inputs it
needs — per pair, the three stored basis rows, the $3\times3$ covariance, and the character pair —
with a SHA-256 digest in `code/data/SHA256SUMS` that the script verifies and prints. This is a
derived extract: reproducing the extraction itself needs the full O25 checkpoints of the Q5a-O5
campaign, which are not redistributed here (`--from-checkpoints` with `--checkpoint-dir` uses
them).

## The One Conditional Obstruction in the Corpus

Under Q5a's depth and weight hypotheses, every subsequential limit of the rescaled admissibility
forms on toric targets with finitely many Fourier modes is the zero form, a bounded multiplication
form with no derivative term, or $+\infty$ (Q5a Theorem 6.1, Remark 6.2). That is a conditional
obstruction to one route, for common scalar normalisations of the published form; it leaves open a
bridge on a different filtration, under an anisotropic or non-scalar renormalisation, or on a
non-toric target, and it proves nothing about the existence of the equivariant map.

## Claims No Longer Asserted in Version 2.0

Q7 withdraws its own assertions and evidential uses; it takes no position on the editorial
status of the papers cited. Notably: the reading of the central value as a Casimir signature;
the fit $|A_H - A_z| \approx 0.364\,q^{-0.52}$ as evidence for asymptotic isotropy; the
transfer of block-diagonality from a discrete compression to cross terms of the continuum
symbol; the attributions to O29 ($H_{\mathrm{eff}} \simeq \mathrm{Sym}^2(V_\rho)$, $d_\rho = 2$),
to O27 (unconditional factorisation through $\mathfrak{su}(2)$), to Q5b Theorem 6.1 (a
co-metric with $A_z > 0$) and Theorem 3.2 (an $O(q^{-1/2})$ rate); and the supports invoked
from Q8, Q10, U1 and Q11 for $A_z = 2$ and $A_H \to 2$.

## What This Paper Does Not Assume

The analysis avoids unnecessary structural assumptions:

- no Lie algebra identification between $\mathfrak{su}(2)$ and $\mathfrak{heis}_3$,
- no arbitrary choice of bridge (uniqueness enforced),
- no dynamical interpretation of the geometry,
- no assumption of isotropy in the central direction (it is not derived here),
- no background spacetime structure.

The results proved here are representation-theoretic, on a supplied carrier, together with
one exact statement about compressions of the discrete Weil Laplacian.

## Keywords

Spatial dimensionality, spectral admissibility, Heisenberg group,
su(2) representation, Casimir operator, sub-Riemannian geometry,
emergent spacetime, isotropy, non-injectivity

## Repository Contents
```
q7/
├── code/    # Reproducible diagnostic (q7_mode_diagnostic.py, q7_symbol_test.py)
├── tex/     # LaTeX sources
└── README.md
```

## Links

- 🌐 Website: https://cosmochrony.org
- 🔗 Zenodo record: https://doi.org/10.5281/zenodo.19802123

## Citation

If you reference this work, please cite:

> J. Beau, *Toward Three Spatial Directions from a Supplied Spin-One Carrier:
> Requirements for an Equivariant Bridge*, 2026.

## Acknowledgements

Portions of the editorial refinement benefited from iterative interactions with
large language models.  
These tools were used as analytical assistants for exploring alternative
formulations, checking internal consistency, and improving clarity.  
All claims, interpretations, and final formulations remain the sole
responsibility of the author.

## Contributions

This repository is intended as a research reference.

Critical feedback, independent analyses, and formal scrutiny are welcome.  
Please open an issue to discuss:

- the mode-selection rule $q \mapsto k(q,c)$, on which the isotropy conjecture now depends,
- Hypothesis [ID] and any route that would supply it,
- potential counterexamples,
- extensions to larger primes or alternative admissible sectors.
