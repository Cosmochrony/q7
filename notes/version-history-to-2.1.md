# Q7 version history (working note, not part of the paper)

Version 2.1 states the current result only. The version narration removed from the paper and README
(what version 1.3 asserted and why version 2.0 withdrew it) is kept here verbatim, as the removed lines
of the 2.0 -> 2.1 diff. The deposited 1.3 and 2.0 records on Zenodo (concept 19802123) carry the full texts.

```diff
diff --git a/README.md b/README.md
index b025581..043a87f 100644
--- a/README.md
+++ b/README.md
@@ -54,5 +54,5 @@ structures: the required bridge is a **missing identification**, not a proved ob
    Q5b's existing spatial principal symbol has rank two everywhere, whereas the positive Casimir
    form has rank three. No invertible transport can identify them.
-   The revised Conjecture 4.7 therefore asks for a new positive spatial target with a specified
+   Conjecture 4.7 therefore asks for a new positive spatial target with a specified
    effective operator and compatible action. Criterion 5.1 tests its graded positive form;
    it does not construct the target or establish equivariance, nor does it assert equality to the Casimir.
@@ -61,5 +61,5 @@ structures: the required bridge is a **missing identification**, not a proved ob
 ## What the O25 Checkpoints Measure
 
-Version 2.0 re-types the numerical section. The three stored basis rows are **pure Fourier
+The three stored basis rows are **pure Fourier
 modes** of indices $\{0, -k, +k\}$ (purity 1.000000 on all twelve tested pairs), so every
 diagonal entry of the compression $\tilde L = B_{\mathrm{eff}} L_{\mathrm{Weil}} B_{\mathrm{eff}}^\dagger$
@@ -78,8 +78,8 @@ $\langle e_m, L_{\mathrm{Weil}} e_m \rangle = 4 - 2\cos(2\pi m/q).$
   which no paper of the corpus derives.
 - An off-diagonal entry of modulus 1 appears exactly when two stored indices differ by
-  $\pm c \bmod q$: this accounts for the two pairs, $(61,3)$ and $(151,5)$, that version 1.3
-  reported as finite-size and sector anomalies.
+  $\pm c \bmod q$: this accounts for the two pairs, $(61,3)$ and $(151,5)$, whose pipeline output
+  looks anomalous; neither is a finite-size or sector effect.
 
-**The closed form describes the stored-basis compression, not directly the published entries.**
+**The closed form describes the stored-basis compression, not directly the reported entries.**
 The pipeline conjugates that compression by the covariance eigenbasis, and that basis is fixed by
 the data only where the covariance spectrum is non-degenerate. Measured: the gap between the two eigenvalues *inside
@@ -97,8 +97,7 @@ non-scalar block — $(151,5)$, where the indices $\pm 73$ differ by $-c$ — an
 split
 $6.9527 / 5.0257$ and the off-diagonal $0.2677$ are solver-dependent; the ratio of that
-off-diagonal to the mean diagonal is the published $\rho_{XY} = 0.057$, in a column that version 1.3's
-running text called the isotropy ratio and its caption called a cross-term ratio. The running
-text's formula, on the stored-basis compression it names, gives $0$; read on the rotated matrix it
-gives $0.2772$; and the script's own isotropy variable gives $0.3218$. Four readings, four values.
+off-diagonal to the mean diagonal is the pipeline's $\rho_{XY} = 0.057$, a cross-term ratio of the
+rotated matrix. It is not the XY isotropy ratio: on the stored-basis compression that ratio is $0$,
+read on the rotated matrix it gives $0.2772$, and the script's own isotropy variable gives $0.3218$.
 
 The per-pair covariance is the same Hermitian object whose spectrum O29 publishes as
@@ -122,6 +121,5 @@ three of the thirty pairs at $q=61$ and nothing at $q=29$.
 The invariant content is the spectrum of the compression, which the closed
 form gives outright at the ten pairs with no index difference equal to $\pm c \bmod q$ and which
-is shifted at the two that have one. §6.3 of the paper works this out, and records one published number,
-$\rho_{XY} = 0.005$ at $(101,3)$, that no versioned quantity reproduces.
+is shifted at the two that have one. §6.3 of the paper works this out.
 
 Reproducible diagnostic: `code/q7_mode_diagnostic.py` prints Fourier indices, purities, the two
@@ -152,14 +150,17 @@ bridge on a different filtration, under an anisotropic or non-scalar renormalisa
 non-toric target, and it proves nothing about the existence of the equivariant map.
 
-## Claims No Longer Asserted in Version 2.0
+## What the Cited Sources and the Computations Do Not Supply
 
-Q7 withdraws its own assertions and evidential uses; it takes no position on the editorial
-status of the papers cited. Notably: the reading of the central value as a Casimir signature;
-the fit $|A_H - A_z| \approx 0.364\,q^{-0.52}$ as evidence for asymptotic isotropy; the
-transfer of block-diagonality from a discrete compression to cross terms of the continuum
-symbol; the attributions to O29 ($H_{\mathrm{eff}} \simeq \mathrm{Sym}^2(V_\rho)$, $d_\rho = 2$),
-to O27 (unconditional factorisation through $\mathfrak{su}(2)$), to Q5b Theorem 6.1 (a
-co-metric with $A_z > 0$) and Theorem 3.2 (an $O(q^{-1/2})$ rate); and the supports invoked
-from Q8, Q10, U1 and Q11 for $A_z = 2$ and $A_H \to 2$.
+The paper states what it may rely on, without taking a position on the editorial status of the
+papers cited. The central value is not a Casimir signature, and the fit
+$|A_H - A_z| \approx 0.364\,q^{-0.52}$ is no evidence for asymptotic isotropy. Block-diagonality
+of a discrete compression says nothing about cross terms of the continuum symbol. O29 does not
+identify $H_{\mathrm{eff}} \simeq \mathrm{Sym}^2(V_\rho)$ or establish $d_\rho = 2$; O27 gives no
+unconditional factorisation through $\mathfrak{su}(2)$; Q5b supplies no $A_z$ (its co-metric,
+Theorem 5.2, has $0$ in the central slot) and no convergence rate. On established inputs, Q8, Q10,
+U1 and Q11 supply no coefficient value: Q8 leaves the common value $2\lambda$ of an invariant form
+free; Q10 bounds $A_H/A_z$ under [U], [W], [R], [T] and gives an absolute value of $A_H$ only under
+its scale input [S], none of which is established; U1 leaves [U] open; Q11 leaves $A_\tau$
+independent of the spatial coefficients.
 
 ## What This Paper Does Not Assume
@@ -222,5 +223,5 @@ Critical feedback, independent analyses, and formal scrutiny are welcome.
 Please open an issue to discuss:
 
-- the mode-selection rule $q \mapsto k(q,c)$, on which the isotropy conjecture now depends,
+- the mode-selection rule $q \mapsto k(q,c)$, on which the isotropy conjecture depends,
 - Hypotheses [ID], [ACT], [INV] and [ADAPT], and any route that would supply them,
 - potential counterexamples,
diff --git a/tex/Q7.tex b/tex/Q7.tex
index 8aceade..a55e3c0 100644
--- a/tex/Q7.tex
+++ b/tex/Q7.tex
@@ -8,5 +8,5 @@
     Equivariant Bridge},
   pdfauthor={Jérôme Beau},
-  pdfsubject={Cosmochrony spectral geometry programme, paper Q7, version 2.0},
+  pdfsubject={Cosmochrony spectral geometry programme, paper Q7, version 2.1},
 }
 \usepackage{geometry}
@@ -73,12 +73,11 @@ Requirements for an Equivariant Bridge\\[4pt]
     gives horizontal isotropy when compatible with the real rotation structure, but does not force
     $A_H = A_z$.
-    Version~2.0 re-types the numerical section.
-    The quantities previously reported as symbol coefficients are compressions of the discrete
+    The numerical quantities are not symbol coefficients: they are compressions of the discrete
     Weil Laplacian $L_{\mathrm{Weil}} = 4I - W_a - W_a^\dagger - W_b - W_b^\dagger$ onto the three
     stored basis rows, which are pure Fourier modes: the central value is the flat-mode Rayleigh
     quotient $4 - 2 - 0 = 2$, exactly and for every tested pair, and the horizontal value is
     $4 - 2\cos(2\pi k/q)$ for the mode index $k$ selected by the pipeline.
-    Consequently the identification of the central value with the Casimir eigenvalue and the
-    fitted convergence of the isotropy gap are no longer asserted: for these modes
+    Consequently the central value is not identified with the Casimir eigenvalue, and no
+    convergence of the isotropy gap is inferred from a fit: for these modes
     $A_H - 2 = 2(1 - \cos(2\pi k/q))$, so vanishing of the gap is equivalent to
     $\mathrm{dist}(k/q, \mathbb{Z}) \to 0$, a property of the mode-selection rule that remains to
@@ -86,5 +85,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     The broader bridge to an extended target remains a missing construction.
     The existing Q5b spatial principal symbol has rank two and cannot equal the rank-three
-    Casimir form; the revised conjecture requires a new target.
+    Casimir form; the conjecture therefore requires a new target.
     Separately, under Q5a's depth and weight hypotheses, no common scalar normalisation of its
     published admissibility form yields a toric differential limit.
@@ -105,5 +104,5 @@ Requirements for an Equivariant Bridge\\[4pt]
       The Cosmochrony programme derives physical observables from the spectral properties of
       admissible non-injective projections on the Weil representation of $\Heis(\mathbb{Z}/q\mathbb{Z})$.
-      Two sub-programmes have now independently produced the integer~$3$ as a structural output.
+      Two sub-programmes independently produce the integer~$3$ as a structural output.
 
       \paragraph{The geometric ``3'' (Q-series).}
@@ -234,8 +233,11 @@ Requirements for an Equivariant Bridge\\[4pt]
       principal-symbol level); its positivity is stipulated inside [H-L], not derived,
       \item the central slot is \emph{degenerate at this order}: Theorem~5.2
-      of~\cite{Beau2026q5b} puts $0$ there and defers a full-rank spatial block to sub-principal
-      corrections (open problem Q5b-O2). A value $A_z = 2$ is asserted in Q8~\cite{BeauQ8} and
-      relayed by~\cite{Beau2026q5b} outside that theorem; this paper does not use it, for the
-      reasons given in Section~\ref{ssec:withdrawn}.
+      of~\cite{Beau2026q5b} puts $0$ there; no term of lower differential order changes it, and a
+      non-degenerate co-metric requires a new operator (Remark~5.4 of~\cite{Beau2026q5b} and open
+      problem Q5b-O2). No value of $A_z$ is established in the sources cited here: Q5b establishes
+      none, and Q8~\cite{BeauQ8} shows that a positive rank-three form supplied with the
+      compatibility conditions of this paper and invariant under the full $\su(2)$ action has
+      $A_H = A_z$ with the common value free
+      (Section~\ref{ssec:withdrawn}).
     \end{itemize}
     Theorem~6.1 of~\cite{Beau2026q5b} is a separate statement, on the Lorentzian signature, and is
@@ -243,5 +245,6 @@ Requirements for an Equivariant Bridge\\[4pt]
     Throughout this paper, $A_z$ denotes the coefficient of $k_Z^2$ that a full-rank spatial block
     would carry. It is treated here as an unevaluated parameter: no value for it is derived in this
-    paper, and the value asserted downstream is not used (Section~\ref{ssec:withdrawn}).
+    paper, and none of the cited sources supplies one on established inputs
+    (Section~\ref{ssec:withdrawn}).
 
     \begin{remark}[Horizontal vs.\ derived: two sub-strata of the spatial sector]
@@ -394,6 +397,7 @@ Requirements for an Equivariant Bridge\\[4pt]
       part of $\mathrm{End}(\Heff)$, of complex dimension six, whereas
       $\Sym^2_\mathbb{C}(\Vr)$ has complex dimension three; the two spaces must not be conflated.
-      Version~1.3 of this paper inferred [ID] from O29 and reported $d_\rho = 2$ as an O29 result.
-      Neither attribution is retained (Section~\ref{ssec:withdrawn}).
+      Neither [ID] nor $d_\rho = 2$ is an O29 result: O29~\cite{Beau2026a33} identifies no such
+      structure, and $d_\rho = 2$ is O26's minimality selection~\cite{Beau2026a30}
+      (Section~\ref{ssec:withdrawn}).
 
   \subsection{Grading Correspondence}
@@ -468,6 +472,4 @@ Requirements for an Equivariant Bridge\\[4pt]
       A unit scalar is available only when the chosen reference normalisation also fixes its
       modulus to one.
-      Version~1.3 stated uniqueness ``up to positive scalar'' with neither the non-vanishing
-      requirement nor this convention.
     \end{remark}
 
@@ -599,10 +601,10 @@ Requirements for an Equivariant Bridge\\[4pt]
       $|\lambda|^2\mathrm{diag}(A_H,A_H,4A_z)$.
       No common scalar removes the relative factor four when both coefficients are positive.
-      Thus [ADAPT] cannot be replaced by a choice of scalar under the earlier hypotheses.
+      Thus [ADAPT] cannot be replaced by a choice of scalar under [ACT] and [INV] alone.
       Q5b Theorem~5.2 supplies only a spatial principal form of rank two under [H-L] and
       [H-lift], with zero central slot in the adapted frame.
       It does not supply the positive rank-three form assumed here.
       Any such form is an additional construction, not a consequence of this lemma or an
-      established import from Q5b; downstream asserted values are not used here.
+      established import from Q5b; no value of either coefficient is imported from the cited sources.
     \end{remark}
 
@@ -620,9 +622,9 @@ Requirements for an Equivariant Bridge\\[4pt]
       Under [ID], an equivariant identification then relates the measured space to this target;
       its pullback has the two-coefficient form of Lemma~\ref{lem:qform}(3).
-      The map has the same type as in the earlier formulation,
+      The map has the type
       \begin{equation}
         \varphi:\Sym^2(\Vr)\xrightarrow{\sim}W_{\mathrm{sp}}\otimes_\mathbb{R}\mathbb{C},
       \end{equation}
-      but its target form must now come from the new construction specified here.
+      and its target form comes from the construction specified here.
       The construction must specify the operator or effective form producing this target and
       its relation to Q5b's horizontal geometry.
@@ -632,5 +634,5 @@ Requirements for an Equivariant Bridge\\[4pt]
 
     \noindent
-    The conjecture no longer asks for equality with Q5b's existing spatial principal symbol.
+    The conjecture does not ask for equality with Q5b's existing spatial principal symbol.
     Under [H-L] and [H-lift], that symbol is
     \[
@@ -642,14 +644,13 @@ Requirements for an Equivariant Bridge\\[4pt]
     Adding terms of strictly lower ordinary differential order cannot change that principal
     symbol.
-    The previous target contract is therefore ruled out under Q5b's stated hypotheses;
+    A target equal to that symbol is therefore ruled out under Q5b's stated hypotheses;
     the broader construction in Conjecture~\ref{conj:bridge} remains missing, not ruled out.
     Even for a new positive target, equality with the transported Casimir is a stronger condition
     than the graded form: under [ADAPT] it requires $A_H=A_z$ and a scale choice.
-    Neither that equality nor an absolute coefficient is asserted by the revised conjecture.
+    Neither that equality nor an absolute coefficient is asserted by the conjecture.
 
-    Version~1.3 formulated the target as the spatial block of $\sigma_2(L_{\mathrm{eff}})$
-    ``restricted to the admissible sector of $L^2(\mathbb{R})$ in the large-$q$ limit of Q5a''.
-    That object is not supplied by Q5a and the formulation is withdrawn: in the published version
-    of Q5a~\cite{Beau2026q5a}, under its depth hypothesis [H-depth] and its weight hypothesis
+    Nor can the target be taken as the spatial block of $\sigma_2(L_{\mathrm{eff}})$ restricted
+    to the admissible sector of $L^2(\mathbb{R})$ in a large-$q$ limit of Q5a, since Q5a does not
+    supply that object: in Q5a~\cite{Beau2026q5a}, under its depth hypothesis [H-depth] and its weight hypothesis
     [H-$w'$], every subsequential limit of the rescaled admissibility forms on toric targets with
     finitely many Fourier modes is the zero form, a bounded multiplication form containing no
@@ -693,6 +694,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     span a real two-dimensional subspace of the weight-$\pm 1$ sector, and the central direction
     $k_Z$ with $e_0$.
-    (The pairing $k_X \leftrightarrow e_+$, $k_Y \leftrightarrow e_-$ used in version~1.3 matched
-    real cotangent coordinates with complex weight vectors and is corrected here.)
+    (Real cotangent coordinates are not paired directly with the complex weight vectors $e_\pm$.)
     Test whether:
     \begin{equation}
@@ -721,9 +721,9 @@ Requirements for an Equivariant Bridge\\[4pt]
     to a continuum symbol coefficient is exactly the identification that Criterion~\ref{crit:symbol}
     presupposes and that no source supplies.
-    Version~1.3 reported the isotropy condition as resolved, citing Q8~\cite{BeauQ8} for
-    $A_z = C_{\mathfrak{su}(2)} = 2$ and Q10~\cite{Beau2026q10} with U1~\cite{Beau2026u1} for
-    $A_H \to 2$.
-    Those supports are no longer asserted here; see Section~\ref{ssec:withdrawn} for the reasons,
-    which concern what this paper may rely on, not the editorial status of those papers.
+    The cited sources do not resolve the isotropy condition on established inputs:
+    Q8~\cite{BeauQ8} leaves the common value of an invariant form free,
+    Q10~\cite{Beau2026q10} bounds the ratio $A_H/A_z$ only under [U], [W], [R] and [T] and gives an
+    absolute value only under its scale input [S], none of which is established, and
+    U1~\cite{Beau2026u1} leaves [U] open (Section~\ref{ssec:withdrawn}).
   \end{remark}
 
@@ -750,10 +750,9 @@ Requirements for an Equivariant Bridge\\[4pt]
       F_c W_a F_c^\dagger = W_b, \qquad F_c W_b F_c^\dagger = W_a^\dagger .
     \end{equation}
-    (Version~1.3 printed the first relation as $F_c W_a F_c^\dagger = W_b^\dagger$; the sign is
-    corrected here. The diagnostic of Section~\ref{ssec:diagnostic} checks both relations, the
-    commutation $[F_c, L_{\mathrm{Weil}}] = 0$ and the deviation of the $\tilde F$ eigenvalue phases
-    from multiples of $90^\circ$, and prints the residuals. The commutation is in any case
-    unaffected by the sign, since $L_{\mathrm{Weil}}$ contains each generator together with its
-    adjoint.)
+    (The diagnostic of Section~\ref{ssec:diagnostic} checks both relations, the commutation
+    $[F_c, L_{\mathrm{Weil}}] = 0$ and the deviation of the $\tilde F$ eigenvalue phases from
+    multiples of $90^\circ$, and prints the residuals. The commutation does not depend on the sign
+    convention in~\eqref{eq:fc-conj}, since $L_{\mathrm{Weil}}$ contains each generator together
+    with its adjoint.)
 
     \begin{proposition}[Metaplectic symmetry]
@@ -787,7 +786,7 @@ Requirements for an Equivariant Bridge\\[4pt]
 
     \noindent
-    The mathematical content of Corollary~\ref{cor:blockdiag} is unaffected by this revision: it is
-    the standard consequence that an operator commuting with $F_c$ preserves its eigenspaces.
-    Its \emph{application} to the O25 data is withdrawn, on two counts.
+    Corollary~\ref{cor:blockdiag} is the standard consequence that an operator commuting with $F_c$
+    preserves its eigenspaces.
+    It does not \emph{apply} to the O25 data, on two counts.
     First, its hypothesis is not met by those data: the three stored rows are not
     $F_c$-eigenvectors, the off-diagonal entries of
@@ -797,5 +796,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     Second, even with the hypothesis met, the conclusion would concern off-diagonal entries of a
     compression of the discrete operator $L_{\mathrm{Weil}}$ on $\mathbb{C}^q$, not cross terms of
-    $\sigma_2(L_{\mathrm{eff}})$: version~1.3 asserted the latter, which requires the missing
+    $\sigma_2(L_{\mathrm{eff}})$; a statement about the latter would require the missing
     identification of Criterion~\ref{crit:symbol}.
     What does hold for these data, exactly and for the reason given in
@@ -872,5 +871,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     $\rho_{Z|XY} = \max(|(L_{\mathrm{sw}})_{01}|, |(L_{\mathrm{sw}})_{02}|)/\bar\lambda$ and
     $|(L_{\mathrm{sw}})_{12}|/\bar\lambda$, where $\bar\lambda$ is the mean diagonal entry.
-    The quantities historically named $A_z := (L_{\mathrm{sw}})_{00}$ and
+    The quantities labelled $A_z := (L_{\mathrm{sw}})_{00}$ and
     $A_H := \tfrac12\bigl((L_{\mathrm{sw}})_{11} + (L_{\mathrm{sw}})_{22}\bigr)$ are therefore
     diagonal entries of a compression of an operator on $\mathbb{C}^q$, read in the covariance
@@ -920,8 +919,8 @@ Requirements for an Equivariant Bridge\\[4pt]
         exactly when an off-diagonal entry of modulus $1$ appears; at $(61,3)$ this happens for two
         pairs of indices, at $(151,5)$ for one.
-        Historically, $A_z$ and $A_H$ were read off this table after the Stage~A ordering; that
-        ordering puts the flat mode first in eleven of the twelve pairs, and at $(q,c) = (101,3)$
-        puts the mode $k = 5$ first, which is the whole origin of the published row
-        $A_z = 2.096$, $A_H = \tfrac12(2.096 + 2.000) = 2.048$.}
+        The pipeline reads $A_z$ and $A_H$ off this table after the Stage~A ordering, which puts the
+        flat mode first in eleven of the twelve pairs and, at $(q,c) = (101,3)$, puts the mode
+        $k = 5$ first; that ordering alone produces the pipeline row $A_z = 2.096$,
+        $A_H = \tfrac12(2.096 + 2.000) = 2.048$ at that pair.}
       \label{tab:numerical}
       \smallskip
@@ -953,9 +952,9 @@ Requirements for an Equivariant Bridge\\[4pt]
     not in phase, and at $(101,3)$ and $(151,5)$ they differ outright, as
     Section~\ref{ssec:rotation} works out.
-    In particular the two pairs at which version~1.3 reported anomalies, $(61,3)$ and $(151,5)$,
+    In particular the two pairs whose pipeline output looks anomalous, $(61,3)$ and $(151,5)$,
     are exactly the two pairs whose stored indices differ by $\pm c$: at $(61,3)$ the flat mode and
     the mode $-3$ differ by $c = 3$, giving $|\tilde L_{01}| = 1$ and hence the reported
     $\rho_{Z|XY} = 0.485$; at $(151,5)$ the modes $\pm 73$ differ by $146 \equiv -5$, which makes the
-    restriction of $\tilde L$ to their plane non-scalar, and the $0.057$ printed in that row follows
+    restriction of $\tilde L$ to their plane non-scalar, and the $0.057$ the pipeline reports in that row follows
     from that together with the covariance rotation, as Section~\ref{ssec:rotation} works out.
     Neither is a finite-size or finite-coverage effect.
@@ -988,7 +987,6 @@ Requirements for an Equivariant Bridge\\[4pt]
       $j(j+1) = 2$ is untouched, and every ratio of entries is unchanged; under
       $L_{\mathrm{Weil}} \mapsto L_{\mathrm{Weil}} - 2I$ it becomes $0$.
-      Version~1.3 read this coincidence as a spectral signature of the bridge, and described the
-      $(101,3)$ entry as a finite-window correction.
-      Both readings are withdrawn (Section~\ref{ssec:withdrawn}).
+      The coincidence is therefore not a spectral signature of the bridge, and the $(101,3)$ entry
+      is not a finite-window correction (Section~\ref{ssec:withdrawn}).
       What survives is the identity~\eqref{eq:rayleigh} and the observation that any test of
       Criterion~\ref{crit:symbol} must compare normalisation-independent quantities.
@@ -998,5 +996,5 @@ Requirements for an Equivariant Bridge\\[4pt]
       \label{rem:finite}
       The cross-term ratio $\rho_{Z|XY} = 0.485$ at $(q,c)=(61,3)$ is not an outlier and not a
-      finite-size artefact of the short window, as version~1.3 stated.
+      finite-size artefact of the short window.
       The stored indices there are $\{0, -3, +3\}$ and the character is $c = 3$, so the flat mode
       differs from \emph{both} non-flat modes by $\pm c$ ($0 - (-3) = 3$ and $0 - 3 = -3$); by the
@@ -1013,5 +1011,5 @@ Requirements for an Equivariant Bridge\\[4pt]
       For all tested primes, the character $c=5$ gives a mode entry close to $6$, against $2$ at
       the flat mode.
-      The mechanism is now explicit: the stored index at $c=5$ is $k = 29, 49, 73, 102$ for
+      The mechanism is explicit: the stored index at $c=5$ is $k = 29, 49, 73, 102$ for
       $q = 61, 101, 151, 211$, that is $k \approx q/2$, and~\eqref{eq:rayleigh} gives
       $4 - 2\cos(2\pi k/q) \approx 4 + 2 = 6$ at the top of the spectrum of
@@ -1032,6 +1030,6 @@ Requirements for an Equivariant Bridge\\[4pt]
       characters: the metaplectic overlap is $0.72$ ($c=3$), $0.60$ ($c=5$) and $0.38$ ($c=7$),
       all below the threshold $0.75$.
-      Version~1.3 attributed this to low trajectory coverage, by analogy with the $(61,3)$
-      anomaly. That attribution is not retained, and the following is what the quantity is.
+      This is not attributed here to low trajectory coverage, by analogy with the $(61,3)$
+      configuration; the following is what the quantity is.
       Stage~B does not test an entry of $\tilde F = B_{\mathrm{eff}} F_c B_{\mathrm{eff}}^\dagger$;
       it compares two unit vectors of $\mathbb{C}^3$, the dominant eigenvector of the covariance
@@ -1060,6 +1058,5 @@ Requirements for an Equivariant Bridge\\[4pt]
         The values are identities, by~\eqref{eq:gap}, once $k(q)$ is read from the checkpoint; they
         are not measurements of a metric anisotropy.
-        The $c = 3$ column is included because it is not monotone, which the version~1.3 reading
-        did not report: $k(q) = 3, 5, 7, 10$ gives $k/q = 0.049, 0.050, 0.046, 0.047$.}
+        The $c = 3$ column is included because it is not monotone: $k(q) = 3, 5, 7, 10$ gives $k/q = 0.049, 0.050, 0.046, 0.047$.}
       \label{tab:asymp}
       \smallskip
@@ -1074,8 +1071,4 @@ Requirements for an Equivariant Bridge\\[4pt]
         \hline
       \end{tabular}
-
-        {\small The extrapolated row at $q = 307$ printed in version~1.3 is removed: it had no
-      asymptotic justification, and its value was placed in the $A_H$ column while being a value of
-      $\Delta$.}
     \end{table}
 
@@ -1084,5 +1077,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     decreases with it; for $c = 3$ it does not, $k/q$ staying in $[0.046, 0.050]$ across the four
     primes.
-    The four-point fit $\Delta \approx 0.364\,q^{-0.52}$ reported in version~1.3 is a fit to the
+    A four-point fit $\Delta \approx 0.364\,q^{-0.52}$ is a fit to the
     sequence $4\sin^2(\pi k(q)/q)$ for the four selected indices $k = 2,3,4,5$: recomputing it from
     those indices alone returns the same exponent to two decimals.
@@ -1105,6 +1098,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     \noindent
     The equivalence is immediate from~\eqref{eq:gap}, since $\sin^2(\pi x)$ vanishes exactly on the
-    integers, and it separates the conjecture from the numerical support version~1.3 offered for
-    it.
+    integers, and it separates the conjecture from any fit to the tabulated values.
     A power-law fit to four values of $\Delta$ cannot establish the conjecture, because those four
     values are determined by four integers $k$ read off the same checkpoints; the fitted exponent
@@ -1118,5 +1110,5 @@ Requirements for an Equivariant Bridge\\[4pt]
     \begin{corollary}[Status of the bridge identification]
       \label{cor:quasithm}
-      The following is the status of the three claims this corollary carried in version~1.3.
+      The status of three claims about these data is as follows.
       \begin{enumerate}
         \item The vanishing of off-diagonal entries of the compression $\tilde L$ holds, for the
@@ -1127,6 +1119,6 @@ Requirements for an Equivariant Bridge\\[4pt]
         identification presupposed by Criterion~\ref{crit:symbol}.
         \item The flat-mode entry equals $2$ exactly and by construction
-        (Remark~\ref{rem:az2}); its reading as a universal coefficient
-        $A_z = C_{\mathfrak{su}(2)}$ is not retained.
+        (Remark~\ref{rem:az2}); it is not a universal coefficient
+        $A_z = C_{\mathfrak{su}(2)}$.
         \item The claim $|A_H - A_z| < 0.05$ for $c \in \{3,7\}$, decreasing with $q$, is false as
         stated: at $c = 3$ the gap is $0.0947$ at $q = 61$ and $0.0842$ at $q = 151$, and it is not
@@ -1155,15 +1147,9 @@ Requirements for an Equivariant Bridge\\[4pt]
         This does not exclude another constructed target.
       \end{enumerate}
-      Version~1.3 read Table~\ref{tab:numerical} as pointing to case~(1).
-      That reading is withdrawn: the table contains compressions of the discrete operator, whose
+      Table~\ref{tab:numerical} does not point to case~(1): it contains compressions of the
+      discrete operator, whose
       two distinct entries differ by the identity~\eqref{eq:gap}, so it discriminates between the
       three cases only through an identification it does not supply.
       These data exclude no case for a new target; Q5b's existing symbol has its known rank-two obstruction.
-      (Papers Q8~\cite{BeauQ8}, Q9~\cite{Beau2026q9} and Q10~\cite{Beau2026q10} cite this
-      structure as ``Q7 Remark~6.8'': Q9 restates all three cases, Q8 and Q10 name Case~3 and
-      Case~1 respectively. In the deposited version~1.3 and here
-      it is Remark~\ref{rem:threecases}, the number 6.8 belonging to
-      Corollary~\ref{cor:quasithm}. The numbering of
-      this paper is unchanged in version~2.0.)
     \end{remark}
 
@@ -1182,5 +1168,5 @@ Requirements for an Equivariant Bridge\\[4pt]
       At eleven pairs the flat mode carries the largest covariance eigenvalue and is ordered first.
       At $(q,c) = (101,3)$ it carries the \emph{smallest} (normalised eigenvalues
-      $1$, $1$, $0.98758$) and is ordered last, which is the whole origin of the published row
+      $1$, $1$, $0.98758$) and is ordered last, which alone produces the pipeline row
       $A_z = 2.096$, $A_H = \tfrac12(2.096 + 2.000) = 2.048$.
       What is resolved there, by about one per cent, is the flat mode's separation from the plane;
@@ -1231,24 +1217,20 @@ Requirements for an Equivariant Bridge\\[4pt]
       The value $5.989188$ survives only as the mean of the two split entries, by preservation of
       the trace, and the ratio $0.2677/\bar\lambda = 0.0574$ with $\bar\lambda = 4.6595$ is the
-      number printed as $\rho_{XY} = 0.057$ in version~1.3.
+      value $0.057$ of the pipeline column $\rho_{XY}$ at that pair.
       Those two entries and that off-diagonal term are therefore solver-dependent, and this is the
       only pair of which that is true.
 
     \paragraph{The column label.}
-      Version~1.3's running text defined that column as the XY isotropy ratio
-      $|\tilde L_{11} - \tilde L_{22}|/|\tilde L_{11}|$, on the matrix its own Stage~C had just
-      introduced, namely the stored-basis compression
-      $\tilde L = B_{\mathrm{eff}} L_{\mathrm{Weil}} B_{\mathrm{eff}}^\dagger$; its table caption
-      instead called that same column a relative cross-term ratio.
-      Neither description matches the number: on the stored-basis compression the formula gives $0$
-      at every pair, both entries being equal there, while the quantity the script prints is the
-      cross-term ratio $|(L_{\mathrm{sw}})_{12}|/\bar\lambda$ of the covariance-rotated matrix,
-      $0.0574$ at $(151,5)$.
-      Read as an isotropy ratio on that rotated matrix, the numerator is
-      $|(L_{\mathrm{sw}})_{11} - (L_{\mathrm{sw}})_{22}| = 1.9270$; over the
-      $|(L_{\mathrm{sw}})_{11}| = 6.9527$ the formula prescribes it gives $0.2772$, and over the
-      mean $A_H = 5.9892$ the script's own isotropy variable uses it gives $0.3218$.
-      Four readings, four distinct values, and the published column is the one its running text
-      does not describe. The label is corrected here.
+      The pipeline column $\rho_{XY}$ is the cross-term ratio
+      $|(L_{\mathrm{sw}})_{12}|/\bar\lambda$ of the covariance-rotated matrix, $0.0574$ at
+      $(151,5)$.
+      It is not the XY isotropy ratio $|\tilde L_{11} - \tilde L_{22}|/|\tilde L_{11}|$ on the
+      stored-basis compression $\tilde L = B_{\mathrm{eff}} L_{\mathrm{Weil}} B_{\mathrm{eff}}^\dagger$,
+      which is $0$ at every pair, both entries being equal there.
+      Read as an isotropy ratio on the rotated matrix, the numerator is
+      $|(L_{\mathrm{sw}})_{11} - (L_{\mathrm{sw}})_{22}| = 1.9270$; over
+      $|(L_{\mathrm{sw}})_{11}| = 6.9527$ it gives $0.2772$, and over the mean $A_H = 5.9892$ used by
+      the script's isotropy variable it gives $0.3218$.
+      These are distinct quantities, and only the first is the one the column reports.
 
     \paragraph{The invariant content.}
@@ -1266,12 +1248,4 @@ Requirements for an Equivariant Bridge\\[4pt]
       entries depart from the stored-basis entries. They coincide only at $(151,5)$.
 
-    \paragraph{One published number is not accounted for.}
-      The value $\rho_{XY} = 0.005$ printed at $(101,3)$ in version~1.3 is reproduced by no
-      quantity the versioned script computes: at that pair the isotropy ratio is $0.0469$ with
-      $A_H$ in the denominator and $0.0458$ with the caption's, and the cross-term ratios vanish to
-      below $10^{-14}$.
-      We record this rather than reconstruct it: the number's provenance is unknown, and no
-      statement of this paper depends on it.
-
   \subsection{Reproducible diagnostic}
     \label{ssec:diagnostic}
@@ -1399,9 +1373,9 @@ Requirements for an Equivariant Bridge\\[4pt]
     defined effective-form construction, rather than identify its positive form with the
     unchanged rank-two symbol.
-    Version~1.3 reported the problem as resolved at the algebraic level, citing Q8~\cite{BeauQ8}
-    for $A_z = C_{\mathfrak{su}(2)} = 2$, with Q10~\cite{Beau2026q10} and
-    U1~\cite{Beau2026u1} for $A_H \to 2$, and treated the absence of cross terms as a theorem about
-    the continuum symbol.
-    None of these is asserted here; see Section~\ref{ssec:withdrawn}.
+    No cited source resolves it at the algebraic level: Q8~\cite{BeauQ8} leaves the common value
+    of an invariant form free, Q10~\cite{Beau2026q10} bounds $A_H/A_z$ only under named inputs that
+    are not established, U1~\cite{Beau2026u1} leaves [U] open, and the block-diagonality of
+    Corollary~\ref{cor:blockdiag} concerns a discrete compression, not cross terms of the continuum
+    symbol (Section~\ref{ssec:withdrawn}).
     Q5b-O2 is, from the point of view of this paper, open.
 
@@ -1425,72 +1399,64 @@ Requirements for an Equivariant Bridge\\[4pt]
     It would bound the scope of the tested construction.
     It would not establish that the two ``3''s are structurally independent, which would require a
-    proof that no admissible bridge exists; no such proof is available, and version~1.3's wording
-    to the contrary is corrected.
+    proof that no admissible bridge exists; no such proof is available.
 
-  \subsection{Claims No Longer Asserted in Version 2.0}
+  \subsection{What the Cited Sources and the Computations Do Not Supply}
     \label{ssec:withdrawn}
 
-    This paper withdraws its own assertions and its own evidential uses.
-    It takes no position on the editorial status of the papers cited below, whose claims are matters
-    for those papers.
-    Each item states what version~1.3 asserted, and why it is not asserted here.
-
-    \paragraph{False statements about sources.}
-      Version~1.3 stated that O29~\cite{Beau2026a33} identifies
-      $\Heff \simeq \Sym^2(\Vr)$ as an $\su(2)$-module and establishes $d_\rho = 2$, and that the
-      spin-$\tfrac12$ doublet is recovered as the Veronese square-root structure of the data.
-      O29 states that no calculation there identifies these objects; $d_\rho = 2$ is O26's
-      minimality selection~\cite{Beau2026a30}.
-      Version~1.3 also stated that O27~\cite{Beau2026a31} proves unconditionally that any
-      admissible morphism factors through $\su(2)$; O27 proves a real-linear factorisation through
-      the admissible quotient, and the $\su(2)$ reading holds only under its saturation hypothesis.
-      Version~1.3 attributed to Theorem~6.1 of~\cite{Beau2026q5b} a co-metric
-      $\mathrm{diag}(-A_\tau, A_H, A_H, A_z)$ with $A_H, A_z > 0$ established there, and a
-      convergence rate $O(q^{-1/2})$ to its Theorem~3.2.
-      The co-metric is its Theorem~5.2 with $0$ in the central slot, $A_H > 0$ is stipulated inside
-      [H-L], that theorem supplies no $A_z$, and its Theorem~3.2 is a Gromov--Hausdorff statement
+    This section states what the sources and computations a reader may expect this paper to rely on
+    do and do not supply.
+    It concerns what this paper uses, not the editorial status of the cited papers.
+
+    \paragraph{What the cited sources state.}
+      O29~\cite{Beau2026a33} does not identify an $\su(2)$-module structure
+      $\Heff \simeq \Sym^2(\Vr)$,
+      does not establish $d_\rho = 2$, and does not recover the spin-$\tfrac12$ doublet as a
+      Veronese square-root structure of the data; $d_\rho = 2$ is O26's minimality
+      selection~\cite{Beau2026a30}.
+      O27~\cite{Beau2026a31} proves a real-linear factorisation of admissible morphisms through the
+      admissible quotient; the $\su(2)$ reading holds only under its saturation hypothesis.
+      The co-metric of Q5b~\cite{Beau2026q5b} is its Theorem~5.2, with $0$ in the central slot;
+      $A_H > 0$ is stipulated inside [H-L], Q5b supplies no $A_z$ (its Remark~5.4), its Theorem~6.1
+      concerns the Lorentzian signature, and its Theorem~3.2 is a Gromov--Hausdorff statement
       carrying no rate.
-      (Q5b does relay $A_z = 2$ from Q8 elsewhere, in its abstract and in its Remark~5.4; the point
-      here is only that Q7 attributed the value to the theorem it cited.)
-      Version~1.3 located the target of Conjecture~\ref{conj:bridge} in a large-$q$ limit of
-      Q5a~\cite{Beau2026q5a} that Q5a does not supply.
-
-    \paragraph{Invalid arguments and invalid readings of computations.}
-      The proof of Proposition~\ref{prop:sym2} in version~1.3 inferred the module structure from
-      the O28 covariance, conflating the symmetric part of $\mathrm{End}(\Heff)$, of complex
-      dimension six, with $\Sym^2(\Vr)$, of complex dimension three.
-      Corollary~\ref{cor:blockdiag} was applied to data that do not satisfy its hypothesis, and its
-      conclusion was transferred from a compression of $L_{\mathrm{Weil}}$ to cross terms of
-      $\sigma_2(L_{\mathrm{eff}})$.
-      The coincidence $A_z = C_{\mathfrak{su}(2)} = 2$ was read as a spectral signature of the
-      bridge; it is the flat-mode Rayleigh quotient $4 - 2 - 0$ and is not normalisation-invariant
+      Q5a~\cite{Beau2026q5a} supplies no large-$q$ limit carrying the target of
+      Conjecture~\ref{conj:bridge}.
+
+    \paragraph{What the computations support.}
+      The module structure of Proposition~\ref{prop:sym2} is not inferred from the O28 covariance:
+      the symmetric part of $\mathrm{End}(\Heff)$, of complex dimension six, is not $\Sym^2(\Vr)$,
+      of complex dimension three.
+      Corollary~\ref{cor:blockdiag} concerns a compression of $L_{\mathrm{Weil}}$; its hypothesis is
+      not met by the O25 data, and it says nothing about cross terms of $\sigma_2(L_{\mathrm{eff}})$.
+      The coincidence $A_z = C_{\mathfrak{su}(2)} = 2$ is the flat-mode Rayleigh quotient
+      $4 - 2 - 0$, is not normalisation-invariant, and is no spectral signature of the bridge
       (Remark~\ref{rem:az2}).
-      The anomalies at $(61,3)$, $(151,5)$ and $q = 211$ were attributed to a short window, a
-      high-energy sector and low coverage; they follow from the stored mode indices
-      (Remarks~\ref{rem:finite}, \ref{rem:c5}, \ref{rem:q211}).
-      The fit $|A_H - A_z| \approx 0.364\,q^{-0.52}$ was offered as evidence for asymptotic
-      isotropy; it is a fit to $4\sin^2(\pi k(q)/q)$ on four selected indices
+      The configurations at $(61,3)$ and $(151,5)$ follow from the stored mode indices, the $c = 5$
+      values from a selected index $k \approx q/2$, and the Stage~B outcome at $q = 211$ is not
+      attributed to low coverage (Remarks~\ref{rem:finite}, \ref{rem:c5}, \ref{rem:q211}).
+      The fit $|A_H - A_z| \approx 0.364\,q^{-0.52}$ is a fit to $4\sin^2(\pi k(q)/q)$ on four
+      selected indices and is no evidence for asymptotic isotropy
       (Conjecture~\ref{conj:isotropy}).
-      Corollary~\ref{cor:quasithm} listed as proved or numerically supported three claims that are
-      restated there with their actual status, one of them false as stated for $c = 3$.
-      Lemma~\ref{lem:schur} was stated as giving uniqueness up to a \emph{positive} scalar without
-      the requirement that makes the scalar real and positive, and without noting that the zero map
-      is not excluded.
-
-    \paragraph{Evidential uses of downstream results.}
-      Version~1.3 reported the isotropy condition as resolved by Q8~\cite{BeauQ8}
-      ($A_z = C_{\mathfrak{su}(2)} = 2$ by Casimir rigidity) and by Q10~\cite{Beau2026q10} with
-      U1~\cite{Beau2026u1} ($A_H \to 2$ by spectral universality), and treated
-      Criterion~\ref{crit:symbol} as thereby conditionally confirmed.
-      This paper no longer uses those results as support.
-      Two reasons are internal to the present revision and suffice for it.
-      First, those statements concern the same two numbers that Section~\ref{ssec:rayleigh} shows to
-      be members of the family~\eqref{eq:rayleigh}, so they cannot corroborate an identification
-      that the family does not supply.
-      Second, Q8's derivation fixes its scalar by matching $A_H$, and Q10's absolute bound is
-      obtained by inserting $A_z = 2$ from Q8, so using them here would make this paper's support
-      circular through its own Remark~\ref{rem:az2}.
-      Q11~\cite{Beau2026q11} imports the same normalisation for the temporal coefficient; it is not
-      used here either.
+      Corollary~\ref{cor:quasithm} gives the status of three claims about these data, one of them
+      false as stated for $c = 3$.
+      Lemma~\ref{lem:schur} gives uniqueness up to one non-zero complex scalar; a real positive
+      scalar requires the convention of Remark~\ref{rem:rigidity}, and the zero map is not excluded.
+
+    \paragraph{What downstream results supply.}
+      None of Q8, Q10, U1 and Q11 resolves the isotropy condition, and Criterion~\ref{crit:symbol}
+      is not confirmed by them: on established inputs none of them supplies the required values.
+      Q8~\cite{BeauQ8} shows that a positive rank-three form supplied with the compatibility
+      conditions of this paper and invariant under the full $\su(2)$ action has
+      $A_H = A_z = 2\lambda$ with $\lambda > 0$ free; the Casimir eigenvalue
+      fixes the normalisation of one invariant form, not a coefficient.
+      Q10~\cite{Beau2026q10} bounds the ratio $A_H/A_z$ under [U], [W], [R] and [T], and gives an
+      absolute value of $A_H$ under its scale input [S], which supplies $A_z$ independently of $A_H$;
+      none of these inputs is established, and the inputs cited for [U] do not imply it.
+      U1~\cite{Beau2026u1} leaves [U] open.
+      Q11~\cite{Beau2026q11} shows that spatial $\mathrm{SU}(2)$ invariance leaves $A_\tau$ independent of the
+      spatial coefficients, so a Casimir normalisation sets no value of $A_\tau$.
+      Moreover, those statements concern the same two numbers that Section~\ref{ssec:rayleigh} shows
+      to be members of the family~\eqref{eq:rayleigh}, so even an established value would not
+      corroborate an identification that the family does not supply.
 
   \section{Conclusion}
@@ -1513,6 +1479,5 @@ Requirements for an Equivariant Bridge\\[4pt]
   The open question is thereby given a precise contract rather than a resolution.
 
-  The numerical section has been re-typed in version~2.0.
-  The quantities it computes are compressions of the discrete Weil Laplacian onto three stored
+  The quantities the numerical section computes are compressions of the discrete Weil Laplacian onto three stored
   basis rows, and those rows are pure Fourier modes, so every diagonal entry is the Rayleigh quotient
   $4 - 2\cos(2\pi m/q)$: exactly $2$ at the flat mode, and $2 + 4\sin^2(\pi k/q)$ at the selected
@@ -1523,13 +1488,13 @@ Requirements for an Equivariant Bridge\\[4pt]
   $\mathrm{dist}(k/q, \mathbb{Z}) \to 0$, a property of the pipeline's mode-selection rule that no
   paper of the corpus derives.
-  Accordingly the identification of the central value with the $\su(2)$ Casimir eigenvalue, the
-  fitted convergence of the isotropy gap, and the transfer of the block-diagonality of a discrete
-  compression to the cross terms of a continuum symbol are no longer asserted, and neither are the
-  downstream supports version~1.3 invoked for them.
+  Accordingly the central value is not identified with the $\su(2)$ Casimir eigenvalue, no
+  convergence of the isotropy gap follows from a fit, the block-diagonality of a discrete
+  compression says nothing about the cross terms of a continuum symbol, and no cited source
+  supplies the coefficient values on established inputs.
 
-  What remains open is what was open before those supports were claimed: the bridge is a missing
-  identification, in the sense of verdict~3 of the programme's bridge discipline.
+  What remains open is the bridge itself: it is a missing identification, in the sense of
+  verdict~3 of the programme's bridge discipline.
   The exact rank-two principal symbol of Q5b is incompatible with a rank-three positive form;
-  the revised conjecture leaves open a different, explicitly constructed target.
+  the conjecture leaves open a different, explicitly constructed target.
   A separate conditional obstruction applies under the depth and
   weight hypotheses of Q5a~\cite{Beau2026q5a}: no common scalar normalisation of the published
```
