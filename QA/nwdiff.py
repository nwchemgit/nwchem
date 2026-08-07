#!/usr/bin/env python3
# ---------------------------------------------------------------------
# NOTE: This script was generated with the assistance of
#       Claude Sonnet 4.5 by Anthropic (https://www.anthropic.com).
# Companion script to nwparse.py, itself a Python conversion of
# nwparse.pl [1] by Ricky A. Kendall, Pacific Northwest National Laboratory
# ---------------------------------------------------------------------
"""
nwdiff.py  --  numerically-aware diff for nwparse output files.

For each line pair the script tokenises both lines and compares
tokens. Tokens that look like floating-point numbers are compared
with a per-field tolerance derived from the number of decimal places
present in the *reference* (ok) token:

    digits  tolerance
    ------  ---------
       10   1e-9          (MR-BWCCSD / MkCCSD energies)
        9   1e-8          (1-e int)
        8   1e-7          (1-e int stdout)
        7   1e-6          (MBPT/CCSD/... total energies)
        5   1e-4          (SCF/DFT/... total energies)
        4   1e-3
        3   1e-2
        2   1e-1
        0   1             (frequencies)

Usage:
    python nwdiff.py [--verbose] reference.nwparse new.nwparse

Options:
    --verbose   Print detailed token-level difference information

Exit codes (mirrors Unix diff behavior):
    0   No differences found
    1   Differences found
    2   Error (e.g. file not found, bad arguments)
"""

import sys
import re


# Tolerance = 1 ULP at the printed precision
_TOLS = {
    10: 1e-9,
    9:  1e-8,
    8:  1e-7,
    7:  1e-6,
    5:  1e-4,
    4:  1e-3,
    3:  1e-2,
    2:  1e-1,
    0:  1.0,
}
_DEFAULT_TOL = 1e-4   # fallback – covers the common 5-digit DFT case


def _tolerance(token: str) -> float:
    """Return the numerical tolerance for *token* based on its printed width."""
    m = re.search(r'\.(\d+)$', token)
    if m:
        ndec = len(m.group(1))
        return _TOLS.get(ndec, _DEFAULT_TOL)
    return _DEFAULT_TOL


def _is_float(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def compare_files(ref_path: str, new_path: str, verbose: bool = False) -> int:
    """
    Compare *ref_path* and *new_path* token by token.

    Returns:
        0  no differences found
        1  differences found
        2  error opening files
    """
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
            print(f"> {new_lines[lineno]}", end='')
            if verbose:
                print(f"  [line {lineno+1}: present only in NEW file]")
            continue

        if lineno >= len(new_lines):
            ndiff += 1
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
            print(f"< {ref_line}")
            print(f"> {new_line}")
            if verbose:
                print(f"  [line {lineno+1}: token count differs:"
                      f" ref={len(ref_toks)} new={len(new_toks)}]")
            continue

        # ---- compare token by token --------------------------------------
        line_differs = False
        diff_details = []   # populated only when verbose=True

        for col, (rtok, ntok) in enumerate(zip(ref_toks, new_toks), start=1):
            if rtok == ntok:
                continue

            if _is_float(rtok) and _is_float(ntok):
                tol = _tolerance(rtok)
                delta = abs(float(rtok) - float(ntok))
                if delta <= tol:
                    # Within tolerance – not a real difference.
                    # Record in verbose mode so the user can see it was checked.
                    if verbose:
                        diff_details.append(
                            f"    token {col}: ref={rtok} new={ntok}"
                            f" |delta|={delta:.2e} tol={tol:.2e}"
                            f" -> WITHIN TOLERANCE (ok)")
                    continue

                # Exceeds tolerance → real difference
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
            print(f"< {ref_line}")
            print(f"> {new_line}")
            if verbose and diff_details:
                print(f"  [line {lineno+1} details:]")
                for detail in diff_details:
                    print(detail)
        elif verbose and diff_details:
            # Line had token differences but all were within tolerance.
            # Only shown in verbose mode so the user knows they were checked.
            print(f"  [line {lineno+1}: all numerical differences within"
                  f" tolerance]")
            for detail in diff_details:
                print(detail)

    return 1 if ndiff > 0 else 0


def main():
    # ---- argument parsing ------------------------------------------------
    args = sys.argv[1:]
    verbose = False

    if '--verbose' in args:
        verbose = True
        args = [a for a in args if a != '--verbose']

    if len(args) != 2:
        print("Usage: python nwdiff.py [--verbose] "
              "reference.nwparse new.nwparse",
              file=sys.stderr)
        sys.exit(2)

    ref_path = args[0]
    new_path = args[1]

    # ---- run comparison --------------------------------------------------
    exit_code = compare_files(ref_path, new_path, verbose=verbose)

    if exit_code == 0:
        if verbose:
            print("nwdiff: no numerical differences found.")
    elif exit_code == 1:
        if verbose:
            print(f"nwdiff: differences found.")

    sys.exit(exit_code)


if __name__ == '__main__':
    main()
