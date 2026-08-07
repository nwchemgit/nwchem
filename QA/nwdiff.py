#!/usr/bin/env python3
# ---------------------------------------------------------------------
# NOTE: This script was generated with the assistance of
#       Claude Sonnet 4.5 by Anthropic (https://www.anthropic.com).
# Companion script to nwparse.py, itself a Python conversion of
# nwparse.pl [1] by Ricky A. Kendall, Pacific Northwest National Laboratory
# ---------------------------------------------------------------------
"""
nwdiff.py  --  numerically-aware diff for nwparse output files.

nwparse.py writes all numerical values with ONE extra decimal digit
beyond the originally requested accuracy. nwdiff.py therefore applies
tolerances at the ORIGINAL accuracy level (one decimal place coarser
than what is stored), so that values agreeing within the original
requested accuracy always compare as equal.

Stored digits -> original accuracy -> tolerance applied here:
    stored   original   tolerance
    ------   --------   ---------
      11       10        1e-9     (MR-BWCCSD / MkCCSD energies)
      10        9        1e-8     (1-e int)
       8        7        1e-6     (MBPT/CCSD/... total energies)
       6        5        1e-4     (SCF/DFT/... total energies)
       5        4        1e-3     (PSPW)
       4        3        1e-2     (3-digit values)
       3        2        1e-1     (2-digit values)
       1        0        1.0      (frequencies)

Usage:
    python nwdiff.py [--diff-output] [--verbose] reference.nwparse new.nwparse

Options:
    --diff-output   Show diff-style output for differing lines
                    (off by default: only a summary is printed)
    --verbose       Print detailed token-level difference information
                    (implies --diff-output)

Exit codes (mirrors Unix diff behavior):
    0   No differences found
    1   Differences found
    2   Error (e.g. file not found, bad arguments)
"""

import sys
import re


# Tolerances are set at the ORIGINAL requested accuracy level,
# keyed by the number of decimal digits STORED in the .nwparse file
# (which is original_digits + 1).
_TOLS = {
    11: 1e-9,    # MR-BWCCSD  (originally 10 digits)
    10: 1e-8,    # 1-e int    (originally  9 digits)
     8: 1e-6,    # MBPT/CCSD  (originally  7 digits)
     6: 1e-4,    # DFT/SCF    (originally  5 digits)
     5: 1e-3,    # PSPW       (originally  4 digits)
     4: 1e-2,    # 3-digit    (originally  3 digits)
     3: 1e-1,    # 2-digit    (originally  2 digits)
     1: 1.0,     # frequencies(originally  0 digits)
}
_DEFAULT_TOL = 1e-4   # fallback


def _tolerance(token: str) -> float:
    """
    Return the numerical tolerance for *token* based on its printed width.

    If the token contains a decimal point, the number of decimal places
    determines the tolerance from _TOLS.

    If the token has NO decimal point but is a valid number
    (integer-formatted), use _TOLS[0] = 1.0 as a safety fallback.
    """
    m = re.search(r'\.(\d+)$', token)
    if m:
        ndec = len(m.group(1))
        return _TOLS.get(ndec, _DEFAULT_TOL)
    elif _is_float(token):
        # Integer-formatted number: no decimal point present.
        # Should not normally occur with the extra-digit approach
        # but kept as a safety fallback.
        return _TOLS.get(0, 1.0)
    return _DEFAULT_TOL


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _normalize_zero(s: str) -> str:
    """
    Normalize signed zero strings so that '-0', '-0.0', '-0.00',
    '-0.000' etc. compare equal to their positive counterparts.
    Returns the input string unchanged if it is not a signed zero.
    """
    if _is_float(s) and float(s) == 0.0:
        m = re.search(r'\.(\d+)$', s)
        if m:
            ndec = len(m.group(1))
            return f"{0.0:.{ndec}f}"
        else:
            return '0'
    return s


def compare_files(ref_path: str, new_path: str,
                  diff_output: bool = False,
                  verbose: bool = False) -> int:
    """
    Compare *ref_path* and *new_path* token by token.

    Args:
        ref_path:    Path to the reference .nwparse file
        new_path:    Path to the new .nwparse file
        diff_output: If True, print diff-style < / > lines for
                     differing lines. If False (default), only a
                     summary count is printed.
        verbose:     If True, print token-level details for every
                     differing line (implies diff_output).

    Returns:
        0  no differences found
        1  differences found
        2  error opening files
    """
    # verbose implies diff_output
    if verbose:
        diff_output = True

    try:
        with open(ref_path) as fref:
            ref_lines = fref.readlines()
    except IOError:
        print(f"nwdiff: error: cannot open reference file: {ref_path}",
              file=sys.stderr)
        return 2

    try:
        with open(new_path) as fnew:
            new_lines = fnew.readlines()
    except IOError:
        print(f"nwdiff: error: cannot open new file: {new_path}",
              file=sys.stderr)
        return 2

    ndiff = 0
    max_lines = max(len(ref_lines), len(new_lines))

    for lineno in range(max_lines):

        # ---- handle files of different lengths ---------------------------
        if lineno >= len(ref_lines):
            ndiff += 1
            if diff_output:
                print(f"> {new_lines[lineno]}", end='')
                if verbose:
                    print(f"  [line {lineno+1}: present only in NEW file]")
            continue

        if lineno >= len(new_lines):
            ndiff += 1
            if diff_output:
                print(f"< {ref_lines[lineno]}", end='')
                if verbose:
                    print(f"  [line {lineno+1}: present only in REF file]")
            continue

        ref_line = ref_lines[lineno].rstrip('\n')
        new_line = new_lines[lineno].rstrip('\n')

        if ref_line == new_line:
            continue  # fast path – identical strings

        ref_toks = ref_line.split()
        new_toks = new_line.split()

        # ---- different number of tokens ----------------------------------
        if len(ref_toks) != len(new_toks):
            ndiff += 1
            if diff_output:
                print(f"< {ref_line}")
                print(f"> {new_line}")
                if verbose:
                    print(f"  [line {lineno+1}: token count differs:"
                          f" ref={len(ref_toks)} new={len(new_toks)}]")
            continue

        # ---- compare token by token --------------------------------------
        line_differs = False
        diff_details = []

        for col, (rtok, ntok) in enumerate(zip(ref_toks, new_toks), start=1):

            # Normalize signed zeros before any comparison so that
            # '-0.0000' and '0.0000' are treated as identical
            rtok_norm = _normalize_zero(rtok)
            ntok_norm = _normalize_zero(ntok)

            if rtok_norm == ntok_norm:
                if verbose and rtok != ntok:
                    diff_details.append(
                        f"    token {col}: ref='{rtok}' new='{ntok}'"
                        f" -> SIGNED ZERO (normalized to equal)")
                continue

            if _is_float(rtok_norm) and _is_float(ntok_norm):
                tol = _tolerance(rtok_norm)
                delta = abs(float(rtok_norm) - float(ntok_norm))
                if delta <= tol:
                    if verbose:
                        diff_details.append(
                            f"    token {col}: ref={rtok} new={ntok}"
                            f" |delta|={delta:.2e} tol={tol:.2e}"
                            f" -> WITHIN TOLERANCE (ok)")
                    continue

                # Exceeds tolerance -> real difference
                line_differs = True
                if verbose:
                    diff_details.append(
                        f"    token {col}: ref={rtok} new={ntok}"
                        f" |delta|={delta:.2e} tol={tol:.2e}"
                        f" -> EXCEEDS TOLERANCE ***")
            else:
                # Non-numeric token mismatch
                line_differs = True
                if verbose:
                    diff_details.append(
                        f"    token {col}: ref='{rtok}' new='{ntok}'"
                        f" -> STRING MISMATCH ***")

        if line_differs:
            ndiff += 1
            if diff_output:
                print(f"< {ref_line}")
                print(f"> {new_line}")
                if verbose and diff_details:
                    print(f"  [line {lineno+1} details:]")
                    for detail in diff_details:
                        print(detail)
        elif verbose and diff_details:
            # Line had token differences but all were within tolerance
            # or were signed-zero normalizations.
            print(f"  [line {lineno+1}: all differences within tolerance"
                  f" or signed-zero normalization]")
            for detail in diff_details:
                print(detail)

    return 1 if ndiff > 0 else 0

def main():
    # ---- argument parsing ------------------------------------------------
    args = sys.argv[1:]
    verbose = False
    diff_output = False

    if '--verbose' in args:
        verbose = True
        args = [a for a in args if a != '--verbose']

    if '--diff-output' in args:
        diff_output = True
        args = [a for a in args if a != '--diff-output']

    if len(args) != 2:
        print("Usage: python nwdiff.py [--diff-output] [--verbose] "
              "reference.nwparse new.nwparse",
              file=sys.stderr)
        sys.exit(2)

    ref_path = args[0]
    new_path = args[1]

    # ---- run comparison --------------------------------------------------
    exit_code = compare_files(ref_path, new_path,
                              diff_output=diff_output,
                              verbose=verbose)

    if exit_code == 0:
        # Only print the success message when --diff-output is requested
        if diff_output:
            print("nwdiff: no numerical differences found.")
    elif exit_code == 1:
        # Always print the failure message so the caller knows why
        # the exit code is non-zero
        print("nwdiff: differences found.")

    sys.exit(exit_code)

if __name__ == '__main__':
    main()
