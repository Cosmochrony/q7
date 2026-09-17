#!/usr/bin/env python3
"""
Q7 mode diagnostic: what the O25 checkpoints supply, and what the Section 6
quantities are as a closed form.

Purpose
-------
Section 6 of Q7 reports two numbers per admissible pair, historically called
A_z and A_H, obtained as diagonal entries of

    L_tilde = B_eff L_Weil B_eff^dag,
    L_Weil  = 4 I - W_a - W_a^dag - W_b - W_b^dag,

where B_eff holds the first three rows of `basis_c` stored by the O25 pipeline,
W_a is the cyclic shift (W_a f)(k) = f(k-1 mod q) and W_b is the character-c
phase diagonal (W_b f)(k) = exp(2 pi i c k / q) f(k).

This script establishes, for every stored pair, that those three rows are pure
Fourier modes and that each diagonal entry equals a Rayleigh quotient in closed
form.  It prints, per (q, c): the Fourier index of each row, its spectral purity,
the measured diagonal entries, the analytic values, the residuals, and an
explicit check of the cross-term rule.

Derivation of the closed form
-----------------------------
Fourier convention used throughout (and matching numpy.fft):

    e_m(k) = q^{-1/2} exp(2 pi i m k / q),   m, k in Z/qZ,

so that numpy.fft.fft(v)_m = sum_k v(k) exp(-2 pi i m k / q).  The purity of a
unit vector v is max_m |fft(v)_m|^2 / q, equal to 1 exactly when v is a single
mode up to a phase.

The two generators act on modes as

    W_a e_m = exp(-2 pi i m / q) e_m        (diagonal: e_m is an eigenvector),
    W_b e_m = e_{m + c}                     (a shift of the Fourier index by c).

Hence for a single mode

    <e_m, W_a e_m>       = exp(-2 pi i m / q),
    <e_m, W_b e_m>       = <e_m, e_{m+c}> = 0        (c != 0 mod q),

and therefore

    <e_m, L_Weil e_m> = 4 - 2 cos(2 pi m / q) - 0 = 4 - 2 cos(2 pi m / q).      (*)

Two consequences, both tested below.

1.  The m = 0 (flat) mode gives exactly 4 - 2 - 0 = 2, for every q and every c.
    The 2 is the additive `4 I` term of L_Weil minus the shift eigenvalue at
    m = 0; it is independent of q and of c by construction, and it is not a
    limit.  Rescaling L_Weil moves it (L/2 gives 1, L - 2I gives 0) while
    leaving every ratio unchanged.

2.  For a mode of index m = k, (*) gives 4 - 2 cos(2 pi k / q), so

        value - 2 = 2 (1 - cos(2 pi k / q)) = 4 sin^2(pi k / q),

    which tends to 0 exactly when dist(k/q, Z) tends to 0.  Whether the stored
    mode indices satisfy that is a property of the pipeline's selection rule,
    not of the spectral geometry.

Off-diagonal entries: <e_m, L_Weil e_n> = -<e_m, W_b e_n> - <e_m, W_b^dag e_n>
is non-zero exactly when m - n = +/- c mod q, in which case its modulus is 1.
This is the cross-term rule checked below.

Inputs
------
By default the script reads the compact bundle shipped with this repository,
`code/data/q7_stage_inputs.npz`, after verifying its SHA-256 digest against
`code/data/SHA256SUMS`.  The bundle holds, per (q, c): the three stored basis
rows, the 3x3 covariance C_c of the per-shell projections, and the conjugate
character pair.  It is a DERIVED extract: reproducing the extraction itself
needs the full O25 checkpoints, which are not redistributed here.

With --from-checkpoints the script recomputes everything from those checkpoints.
Default directory: the O25 outputs of the Q5a-O5 campaign,
`simulation/gravity/q5a-o5/o25_outputs/q{q}_o25.npz`, relative to the workspace
root.  Each file stores `pairs` (conjugate character pairs), `basis_c` (the
accumulated Gram-Schmidt basis per pair; its first three rows are the
"admissible directions" used by Section 6) and `pi_c` (the per-shell
H_eff-projections).  Pass --checkpoint-dir to point elsewhere.

Usage
-----
    python3 q7_mode_diagnostic.py [--primes 61 101 151 211] [--chars 3 5 7]
    python3 q7_mode_diagnostic.py --from-checkpoints --checkpoint-dir <dir>
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys

import numpy as np

HEFF_DIM = 3
DEFAULT_PRIMES = [61, 101, 151, 211]
DEFAULT_CHARS = [3, 5, 7]
# Resolved against this file, not the current working directory, so the script
# runs from anywhere inside a full workspace checkout.  The checkpoints
# themselves are produced by the O25 pipeline and are NOT distributed with this
# repository: pass --checkpoint-dir to point at your own copy.
BUNDLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data",
                      "q7_stage_inputs.npz")
DEFAULT_DIR = os.path.normpath(
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "simulation", "gravity", "q5a-o5", "o25_outputs",
    )
)

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from q7_symbol_test import (  # noqa: E402
    extract_Beff_from_checkpoint,
    extract_pi_c_from_checkpoint,
    load_checkpoint,
    weil_laplacian,
    weil_phase,
    weil_shift,
)


def fourier_index_and_purity(v: np.ndarray) -> tuple[int, float]:
    """Return (signed Fourier index, purity) of a unit vector v in C^q.

    With e_m(k) = q^{-1/2} exp(2 pi i m k / q), the coefficient of a unit vector v
    on e_m has modulus q^{-1/2} |numpy.fft.fft(v)_{-m}|, so |fft(v)|^2 / q is the
    modulus squared of a coefficient and its maximum over the index is the purity.
    Purity 1 means v is a single mode up to a phase.  (The index printed below is
    read off the same spectrum, so it is consistent with this convention.)
    """
    q = v.size
    spec = np.abs(np.fft.fft(v)) ** 2 / q
    m = int(np.argmax(spec))
    signed = m if m <= q // 2 else m - q
    return signed, float(spec[m])


def analytic_rayleigh(q: int, m: int) -> float:
    """Closed form (*) of the Rayleigh quotient of L_Weil on the mode e_m."""
    return 4.0 - 2.0 * np.cos(2.0 * np.pi * m / q)


def check_mode_actions(q: int, c: int, m: int) -> tuple[float, float]:
    """Return the two residuals of the mode action identities at index m.

    Checks  W_a e_m = exp(-2 pi i m / q) e_m  and  W_b e_m = e_{m + c}.
    """
    k = np.arange(q)
    e = np.exp(2j * np.pi * m * k / q) / np.sqrt(q)
    e_shift = np.exp(2j * np.pi * ((m + c) % q) * k / q) / np.sqrt(q)
    res_a = float(
        np.abs(weil_shift(q) @ e - np.exp(-2j * np.pi * m / q) * e).max()
    )
    res_b = float(np.abs(weil_phase(q, c) @ e - e_shift).max())
    return res_a, res_b


def load_bundle(path: str = BUNDLE):
    """Load the compact inputs shipped with this repository, checksum first.

    The bundle holds, per (q, c): the three stored basis rows, the 3x3 covariance
    C_c of the per-shell projections, and the conjugate character pair.  It is
    what the statements of Section 6 need; the full O25 checkpoints are not
    redistributed here.  SHA256SUMS sits next to it.
    """
    sums = os.path.join(os.path.dirname(path), "SHA256SUMS")
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    if os.path.exists(sums):
        expected = open(sums).read().split()[0]
        if digest != expected:
            raise SystemExit(
                f"bundle checksum mismatch\n  expected {expected}\n  got      {digest}"
            )
        print(f"Bundle: {path}\n  sha256 {digest} (verified against SHA256SUMS)")
    else:
        print(f"Bundle: {path}\n  sha256 {digest}\n"
              f"  WARNING: {os.path.basename(sums)} absent, digest NOT verified")
    return np.load(path)


def report_pair(q: int, c: int, data, bundle=None) -> dict | None:
    if bundle is not None:
        key = f"basis_{q}_{c}"
        if key not in bundle.files:
            print(f"  q={q} c={c}: not in bundle")
            return None
        B = bundle[key]
        cov = bundle[f"cov_{q}_{c}"]
        pair = tuple(int(x) for x in bundle[f"pair_{q}_{c}"])
    else:
        try:
            B = extract_Beff_from_checkpoint(data, c)
        except (KeyError, ValueError) as exc:
            print(f"  q={q} c={c}: {exc}")
            return None
        w = np.asarray(extract_pi_c_from_checkpoint(data, c), dtype=complex)
        cov = w.conj().T @ w / len(w)
        pair = tuple(int(x) for x in data["pairs"][_pair_row(data, c)])

    L = weil_laplacian(q, c)
    L_tilde = B @ L @ B.conj().T

    rows = []
    for i, v in enumerate(B):
        m, purity = fourier_index_and_purity(v)
        measured = float(L_tilde[i, i].real)
        analytic = analytic_rayleigh(q, m)
        res_a, res_b = check_mode_actions(q, c, m)
        rows.append(
            {
                "row": i,
                "index": m,
                "purity": purity,
                "measured": measured,
                "analytic": analytic,
                "residual": measured - analytic,
                "res_Wa": res_a,
                "res_Wb": res_b,
            }
        )

    print(f"\n  q={q}  c={c}   (conjugate pair stored as {pair})")
    print("    row  index   purity      measured      analytic      residual"
          "    |W_a id|   |W_b id|")
    for r in rows:
        print(f"    {r['row']:3d}  {r['index']:5d}   {r['purity']:.6f}   "
              f"{r['measured']:11.6f}   {r['analytic']:11.6f}   "
              f"{r['residual']:+.2e}   {r['res_Wa']:.1e}   {r['res_Wb']:.1e}")

    # Cross-term rule: |L_tilde[i,j]| is 1 iff m_i - m_j = +/- c mod q, else 0.
    print("    cross terms (i,j): |measured|  predicted  index difference")
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            diff = (rows[i]["index"] - rows[j]["index"]) % q
            predicted = 1.0 if diff in (c % q, (-c) % q) else 0.0
            meas = float(abs(L_tilde[i, j]))
            flag = "ok" if abs(meas - predicted) < 1e-8 else "MISMATCH"
            print(f"      ({i},{j}):      {meas:.6f}    {predicted:.1f}      "
                  f"{diff:4d}   {flag}")

    # Covariance spectrum: this, and not the equality of the Rayleigh quotients,
    # is what decides whether the Stage-A eigenbasis is free to rotate.
    ev = np.linalg.eigvalsh(cov)[::-1]
    ev_n = ev / ev[0]
    # The covariance is block-diagonal with respect to {flat} + {e_{-k}, e_{+k}}.
    # What decides whether the eigenbasis of the horizontal plane is determined is
    # the gap INSIDE that plane, not the gap between the second and third
    # eigenvalues of the whole matrix: at (101,3) the latter is the flat mode's.
    plane = np.linalg.eigvalsh(cov[1:, 1:])[::-1] / ev[0]
    gap = float(abs(plane[0] - plane[1]))
    gap23 = float(abs(ev_n[1] - ev_n[2]))
    U = np.linalg.eigh(cov)[1][:, ::-1]
    L_sw = U.conj().T @ L_tilde @ U
    off_sw = float(
        np.abs(L_sw - np.diag(np.diag(L_sw))).max()
    )
    print("    covariance spectrum (normalised): "
          f"{np.round(ev_n, 6)}   gap(2,3) = {gap23:.2e}")
    print(f"    gap inside the horizontal plane: {gap:.2e}   "
          f"{'degenerate to numerical precision' if gap < 1e-12 else 'resolved'}")
    print("    rotated by the covariance eigenbasis: diag "
          f"{np.round(np.diag(L_sw).real, 6)}   max |off-diagonal| = {off_sw:.6f}")
    blk = L_tilde[1:, 1:]
    scalar_block = bool(
        np.allclose(blk, blk[0, 0] * np.eye(blk.shape[0]), atol=1e-9)
    )
    if gap < 1e-12 and not scalar_block:
        print("      (plane degenerate to numerical precision AND block non-scalar:"
              " these two numbers are solver-dependent; only the spectrum is"
              " invariant)")
    elif gap < 1e-12:
        print("      (block is scalar on that plane, so no rotation can split it,"
              " however the degeneracy is classified: these entries are invariant)")
    print("    spectrum of the compression (invariant): "
          f"{np.round(np.linalg.eigvalsh(L_tilde).real, 6)}")

    return {"q": q, "c": c, "rows": rows, "L_tilde": L_tilde, "cov_gap": gap}


def _pair_row(data, c: int) -> int:
    pairs = data["pairs"]
    idx = np.where(pairs[:, 0] == c)[0]
    if len(idx) == 0:
        idx = np.where(pairs[:, 1] == c)[0]
    return int(idx[0])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--primes", nargs="+", type=int, default=DEFAULT_PRIMES)
    ap.add_argument("--chars", nargs="+", type=int, default=DEFAULT_CHARS)
    ap.add_argument("--checkpoint-dir", type=str, default=DEFAULT_DIR)
    ap.add_argument("--from-checkpoints", action="store_true",
                    help="use the full O25 checkpoints instead of the bundled inputs")
    args = ap.parse_args()

    print(__doc__)
    bundle = None
    if not args.from_checkpoints and os.path.exists(BUNDLE):
        bundle = load_bundle()
    else:
        print(f"Checkpoint directory: {os.path.abspath(args.checkpoint_dir)}")

    worst_purity = 1.0
    worst_residual = 0.0
    pairs_done = 0
    for q in args.primes:
        data = None
        if bundle is None:
            try:
                data, path = load_checkpoint(q, args.checkpoint_dir)
            except FileNotFoundError as exc:
                print(f"\n  q={q}: {exc}")
                continue
            print(f"\n== q={q}  checkpoint {os.path.basename(path)}")
        else:
            print(f"\n== q={q}  (bundled inputs)")
        for c in args.chars:
            out = report_pair(q, c, data, bundle)
            if out is None:
                continue
            pairs_done += 1
            for r in out["rows"]:
                worst_purity = min(worst_purity, r["purity"])
                worst_residual = max(worst_residual, abs(r["residual"]))

    print("\nSummary")
    print(f"  pairs processed                              : {pairs_done}")
    if pairs_done == 0:
        print("  NOTHING WAS PROCESSED: no bundle and no readable checkpoint.")
        print("  The statements below are vacuous; re-run with --checkpoint-dir "
              "or restore code/data/q7_stage_inputs.npz.")
        return 2
    print(f"  lowest Fourier purity over all reported rows : {worst_purity:.6f}")
    print(f"  largest |measured - analytic| residual       : {worst_residual:.2e}")
    print("  Every reported diagonal entry is the Rayleigh quotient "
          "4 - 2 cos(2 pi m / q)")
    print("  of L_Weil on a single Fourier mode: 2 exactly at m = 0, and "
          "2 + 4 sin^2(pi k / q) at m = k.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
