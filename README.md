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
      (O23 Theorem 3.1, conditional on that carrier) and from the measured projection
      space $H_{\mathrm{eff}} = \mathbb{C}^3$ of rank three (O28).

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
  at the flat mode. It is not the $\mathfrak{su}(2)$ Casimir eigenvalue, and it is not
  normalisation-invariant: $L/2$ gives 1, $L - 2I$ gives 0.
- The horizontal value is $4 - 2\cos(2\pi k/q)$ for the selected index $k$, so
  $A_H - 2 = 4\sin^2(\pi k/q)$ **identically**. Vanishing of that gap is equivalent to
  $\mathrm{dist}(k/q, \mathbb{Z}) \to 0$, a property of the pipeline's mode-selection rule,
  which no paper of the corpus derives.
- An off-diagonal entry of modulus 1 appears exactly when two stored indices differ by
  $\pm c \bmod q$: this accounts for the two pairs, $(61,3)$ and $(151,5)$, that version 1.3
  reported as finite-size and sector anomalies.

Reproducible diagnostic: `code/q7_mode_diagnostic.py` prints Fourier indices, purities,
measured and analytic values, residuals and the cross-term check, with the Fourier convention
and the checkpoints documented. Largest residual: $3.6 \times 10^{-15}$.

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
├── pdfs/    # Archived versioned PDFs
├── tex/     # LaTeX sources
└── README.md
```

## Links

- 📄 [Paper PDF](https://github.com/Cosmochrony/.../Q7.pdf)
- 🌐 Website: https://cosmochrony.org

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
