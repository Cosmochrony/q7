This repository contains the source of the **Q7 Cosmochrony paper**  
[*Three Admissible Directions and Three Spatial Dimensions:
A Structural Bridge Candidate between $H_{\mathrm{eff}} \simeq \mathbb{C}^3$
and the Horizontal Geometry of $\mathrm{Heis}_3(\mathbb{R})$*](out/Q7.pdf).

This work addresses a central structural question of the Cosmochrony
spectral admissibility programme:

> Are the three admissible directions arising from quaternionic minimality
> identical to the three spatial directions of the emergent Lorentzian geometry?

The paper shows that this question can be formulated precisely, and that
its resolution reduces to a **single computable spectral criterion**.

## Conceptual Overview

The paper proceeds in four logical steps:

1. **Separation of the two “3”s**  
   Two independent mechanisms produce the integer 3:
    - the geometric 3 from the Carnot structure of $\mathrm{Heis}_3$ (Q5b),
    - the admissible 3 from quaternionic minimality and projection (O23–O29).

   These structures are *a priori distinct* and cannot be identified directly.

2. **Algebraic obstruction and representation-level bridge**  
   The Lie algebras $\mathfrak{su}(2)$ and $\mathfrak{heis}_3$ are
   non-isomorphic (semisimple vs nilpotent).  
   Any identification must therefore occur at the level of **representations**.

3. **Rigidity via symmetric square structure**  
   The admissible space satisfies:
   \[
   H_{\mathrm{eff}} \simeq \mathrm{Sym}^2(V_\rho),
   \]
   the spin-1 irreducible representation of $\mathfrak{su}(2)$.

   By Schur’s lemma:
    - any equivariant bridge is **unique up to scalar**,
    - the identification problem becomes **binary (existence vs no-go)**.

4. **Reduction to a spectral symbol criterion**  
   The identification is equivalent to testing:
   \[
   \sigma_2(L_{\mathrm{eff}})\big|_{\mathrm{Sym}^2(V_\rho)}
   = A_H(k_X^2 + k_Y^2) + A_z k_Z^2.
   \]

   This reduces the problem to:
    - absence of cross terms,
    - equality $A_H = A_z$ (isotropy).

## Core Claims

The paper establishes the following results:

1. **The two “3”s are structurally distinct but comparable**  
   They arise from different algebraic levels (Carnot vs representation),
   and cannot be identified via Lie algebra isomorphism.

2. **The admissible sector is the spin-1 representation of $\mathfrak{su}(2)$**  
   The effective space $H_{\mathrm{eff}}$ is identified with
   $\mathrm{Sym}^2(V_\rho)$.

3. **Any bridge is rigid and unique**  
   By irreducibility, there is at most one equivariant identification
   (up to normalization).

4. **The spatial metric is constrained by the Casimir structure**  
   The $\mathfrak{su}(2)$ Casimir induces a quadratic form that matches
   the spatial block of the effective symbol up to two scalars.

5. **Cross terms vanish structurally**  
   The metaplectic symmetry
   \[
   [F_c, L_{\mathrm{Weil}}] = 0
   \]
   implies block-diagonalization and forbids mixing between sectors.

6. **The problem reduces to a single scalar condition**  
   The identification holds if and only if:
   \[
   A_H = A_z.
   \]

## Numerical Evidence

Using O25 pipeline checkpoints for $q \in \{61, 101, 151\}$:

- Cross terms are numerically zero (consistent with analytic result),
- The central coefficient satisfies:
  \[
  A_z \approx 2 = C_{\mathfrak{su}(2)},
  \]
- The isotropy gap
  \[
  \Delta(q) = |A_H(q) - A_z(q)|
  \]
  is:
    - small ($< 0.05$ in low-energy sector),
    - monotonically decreasing with $q$,
    - consistent with power-law decay $\sim q^{-0.44}$.

This supports the conjecture:

\[
|A_H(q) - A_z(q)| \to 0 \quad \text{as } q \to \infty.
\]

## Core Result

The identification problem is reduced to:

> a single asymptotic scalar equality, for which all current analytical
> and numerical evidence is positive.

If confirmed, this implies:

- the three spatial dimensions correspond to the spin-1 representation,
- the spatial metric is governed by the $\mathfrak{su}(2)$ Casimir:
  \[
  g_{\mathrm{spatial}} \sim C_{\mathfrak{su}(2)}.
  \]

## What This Paper Does Not Assume

The analysis avoids unnecessary structural assumptions:

- no Lie algebra identification between $\mathfrak{su}(2)$ and $\mathfrak{heis}_3$,
- no arbitrary choice of bridge (uniqueness enforced),
- no dynamical interpretation of the geometry,
- no assumption of isotropy (tested, not imposed),
- no background spacetime structure.

The result is entirely **representation-theoretic and spectral**.

## Keywords

Spatial dimensionality, spectral admissibility, Heisenberg group,
su(2) representation, Casimir operator, sub-Riemannian geometry,
emergent spacetime, isotropy, non-injectivity

## Repository Contents
```
paper/
├── pdf/ # Compiled paper PDF
├── tex/ # LaTeX sources
└── README.md
```

## Links

- 📄 [Paper PDF](https://github.com/Cosmochrony/.../Q7.pdf)
- 🌐 Website: https://cosmochrony.org

## Citation

If you reference this work, please cite:

> J. Beau, *Three Admissible Directions and Three Spatial Dimensions:
> A Structural Bridge Candidate between $H_{\mathrm{eff}} \simeq \mathbb{C}^3$
> and the Horizontal Geometry of $\mathrm{Heis}_3(\mathbb{R})$*, 2026.

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

- the asymptotic isotropy conjecture,
- the role of the Casimir in the spatial metric,
- potential counterexamples,
- extensions to larger primes or alternative admissible sectors.
