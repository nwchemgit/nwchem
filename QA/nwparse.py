#!/usr/bin/env python3
# ---------------------------------------------------------------------
# NOTE: This script was converted from Perl to Python with the
#       assistance of Claude Sonnet 4.5 by Anthropic
#       (https://www.anthropic.com).
# Original Perl script: nwparse.pl [1]
# Original author: Ricky A. Kendall, Pacific Northwest National Laboratory
# ---------------------------------------------------------------------
"""
Python script to parse nwchem output files.

Usage: python nwparse.py [-h||-H||-help] [-d] [-q] [-s suffix]
           nwchem_output_file_1 [nwchem_output_file_2 ...]

  -d := debug mode
  -q := quiet mode (nothing to stdout) **
  -s := override default suffix of .nwparse to user supplied 'suffix'
  -h := prints this help message (equivalent to -help or -H)

 **:Note: if -d is set -q is ignored

All numerical values are written with ONE extra decimal digit beyond
the requested accuracy so that boundary rounding cases (e.g. 3680.50
vs 3680.51 at 0-digit accuracy) are stored as 3680.5 vs 3680.5 and
compare equal. nwdiff.py applies tolerances at the ORIGINAL accuracy
level when comparing .nwparse files.
"""

import sys
import re


def usage():
    print("\n\nUsage: python nwparse.py [-h||-H||-help] [-d] [-q] [-s suffix]"
          "  nwchem_output_file_1 [nwchem_output_file_2 ...]\n")
    print(" -d := debug mode")
    print(" -q := quiet mode (nothing to stdout) **")
    print(" -s := override default suffix of .nwparse to user supplied"
          " 'suffix'")
    print(" -h := prints this help message (equivalent to -help or -H)")
    print("\n **:Note: if -d is set -q is ignored")


def set_to_digits(value, digits):
    """
    Rounds a float to a specified number of decimal digits using
    Python's built-in round(), which uses IEEE 754 banker's rounding
    and avoids floating-point accumulation errors that the manual
    multiply/truncate/divide approach in the original Perl script [1]
    suffers from.
    Signed zero (-0.0) is normalized to 0.0 to avoid -0.000 vs 0.000
    spurious differences in output.
    Values are rounded to (digits+1) decimal places so that boundary
    rounding cases are preserved with one extra digit of precision.
    nwdiff.py applies tolerances at the original 'digits' accuracy level.
    """
    result = round(float(value), digits + 1)
    # Normalize -0.0 to 0.0 to avoid "-0.000" vs "0.000" diffs
    if result == 0.0:
        result = 0.0
    return result


def _is_float(s: str) -> bool:
    """Return True if string s can be converted to a float."""
    try:
        float(s)
        return True
    except ValueError:
        return False


def main():
    quiet = 0
    debug = 0
    suffix = '.nwparse'
    files_to_parse = []
    args = sys.argv[1:]
    num_argv = len(args)

    if num_argv == 0:
        usage()
        sys.exit("fatal error: no file to parse")

    # Parse command-line arguments
    get_suffix = False
    for argument in args:
        if get_suffix:
            suffix = argument
            if not suffix.startswith('.'):
                suffix = '.' + suffix
            get_suffix = False
        elif argument in ('-h', '-help', '-H'):
            usage()
            sys.exit(0)
        elif argument == '-d':
            print("debug: debug turned on at command line")
            debug = 1
        elif argument == '-s':
            get_suffix = True
        elif argument == '-q':
            quiet = 1
        elif argument.startswith('-'):
            print(f"\n\nUnrecognized argument: {argument}")
            sys.exit("fatal error")
        else:
            files_to_parse.append(argument)

    if debug:
        quiet = 0

    if debug:
        print(f"\ndebug:number of arguments: {num_argv}\n")
        print(f"debug: arguments {' '.join(args)}")
        print(f"\ndebug: suffix is {suffix}")
        print(f"\ndebug: files to parse {' '.join(files_to_parse)}")

    # Process each file
    for filename in files_to_parse:
        atoms = []
        coords = []
        grads = []

        if debug:
            print(f"\ndebug: file to open is {filename}")

        try:
            file_to_parse = open(filename, 'r')
        except IOError:
            sys.exit(f"fatal error: Could not open file:{filename}")

        fileout = filename + suffix
        if debug:
            print(f"\ndebug: file for parsed output is: {fileout}")

        try:
            file_output = open(fileout, 'w')
        except IOError:
            sys.exit(f"fatal error: Could not open file:{fileout}")

        sgroup = 0
        selcipt_block = 0
        gradient_block = 0
        dirdyv_block = 0
        lines = 0
        ci_energy = []
        pt_correc = []
        cipt_ene = []
        pt_norm = []

        for line in file_to_parse:
            lines += 1

            # ----------------------------------------------------------------
            # Handle blank lines: close selcipt_block or gradient_block
            # ----------------------------------------------------------------
            if re.search(r'^\s*$', line):

                if selcipt_block:
                    selcipt_block = 0
                    num_energies = len(ci_energy)

                    if len(pt_correc) != num_energies:
                        sys.exit("number of ci+pt energies different than "
                                 "number of corrections")
                    if len(cipt_ene) != num_energies:
                        sys.exit("number of ci+pt energies different than "
                                 "number of summed ci+pt energies")
                    if len(pt_norm) != num_energies:
                        sys.exit("number of ci+pt energies different than "
                                 "number of pt norms")

                    if not quiet:
                        print(" ci energy   pt correction ci+pt energy PT norm")
                        print(" ----------  ------------- ------------ -------")
                    file_output.write("ci energy   pt correction"
                                      " ci+pt energy PT norm\n")
                    file_output.write("----------  -------------"
                                      " ------------ -------\n")

                    for itok in range(num_energies):
                        if not quiet:
                            print(
                                f"{set_to_digits(ci_energy[itok], 5):12.6f}"
                                f" {set_to_digits(pt_correc[itok], 5):14.6f}"
                                f" {set_to_digits(cipt_ene[itok], 5):13.6f}"
                                f" {set_to_digits(pt_norm[itok], 3):8.4f}")
                        file_output.write(
                            f"{set_to_digits(ci_energy[itok], 5):12.6f}"
                            f" {set_to_digits(pt_correc[itok], 5):14.6f}"
                            f" {set_to_digits(cipt_ene[itok], 5):13.6f}"
                            f" {set_to_digits(pt_norm[itok], 3):8.4f}\n")

                if gradient_block:
                    gradient_block = 0
                    num_atoms = len(atoms)
                    num_grads = len(grads)
                    num_coords = len(coords)

                    if (num_grads / 3) != num_atoms:
                        print(f" num_grads = {num_grads}")
                        print(f" num_atoms = {num_atoms}")
                        sys.exit(" fatal error ")
                    if (num_coords / 3) != num_atoms:
                        print(f" num_coords = {num_coords}")
                        print(f" num_atoms  = {num_atoms}")
                        sys.exit(" fatal error ")

                    if debug:
                        print(f"debug: number of atoms:"
                              f" {num_atoms} {atoms}")
                        print(f"debug: number of grads:"
                              f" {num_grads} {grads}")
                        print(f"debug: number of coords:"
                              f" {num_coords} {coords}")

                    if not quiet:
                        print("   Atoms             Coordinates:")
                    file_output.write("   Atoms             Coordinates:\n")

                    for iatom in range(num_atoms):
                        indx1 = iatom * 3
                        indx2 = indx1 + 1
                        indx3 = indx1 + 2
                        if not quiet:
                            print(
                                f" {atoms[iatom]:>10s}"
                                f" {set_to_digits(coords[indx1], 4):10.4f}"
                                f" {set_to_digits(coords[indx2], 4):10.4f}"
                                f" {set_to_digits(coords[indx3], 4):10.4f}")
                        file_output.write(
                            f" {atoms[iatom]:>10s}"
                            f" {set_to_digits(coords[indx1], 4):10.4f}"
                            f" {set_to_digits(coords[indx2], 4):10.4f}"
                            f" {set_to_digits(coords[indx3], 4):10.4f}\n")

                    if not quiet:
                        print("   Atoms              Gradients:")
                    file_output.write("   Atoms              Gradients:\n")

                    for iatom in range(num_atoms):
                        indx1 = iatom * 3
                        indx2 = indx1 + 1
                        indx3 = indx1 + 2
                        if not quiet:
                            print(
                                f" {atoms[iatom]:>10s}"
                                f" {set_to_digits(grads[indx1], 3):10.4f}"
                                f" {set_to_digits(grads[indx2], 3):10.4f}"
                                f" {set_to_digits(grads[indx3], 3):10.4f}")
                        file_output.write(
                            f" {atoms[iatom]:>10s}"
                            f" {set_to_digits(grads[indx1], 3):10.4f}"
                            f" {set_to_digits(grads[indx2], 3):10.4f}"
                            f" {set_to_digits(grads[indx3], 3):10.4f}\n")

                    atoms = []
                    coords = []
                    grads = []

                # Skip blank lines for all further processing
                continue

            # ----------------------------------------------------------------
            # Warn on failed/warning lines
            # ----------------------------------------------------------------
            if re.search(r'failed', line, re.IGNORECASE) or \
               re.search(r'warning', line, re.IGNORECASE):
                print(line, end='')

            # ----------------------------------------------------------------
            # Detect GA subgroups
            # ----------------------------------------------------------------
            if re.search(r'^ Creating groups', line):
                sgroup = 1

            # ----------------------------------------------------------------
            # P.Frequency block
            # Original Perl [1]: printf "%10.0f", set_to_digits($tok, 0)
            # Now writes one extra digit: %.1f with set_to_digits(tok, 0)
            # ----------------------------------------------------------------
            if re.search(r'^ P\.Frequency', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if not quiet:
                    print(line_tokens[0], end='')
                file_output.write(line_tokens[0])
                for itok in range(1, num_line_tokens):
                    if not quiet:
                        print(
                            f"{set_to_digits(line_tokens[itok], 0):11.1f} ",
                            end='')
                    file_output.write(
                        f"{set_to_digits(line_tokens[itok], 0):11.1f} ")
                if not quiet:
                    print()
                file_output.write("\n")

            # ----------------------------------------------------------------
            # 1-e int block
            # Original Perl [1]: printf "%17.9f", set_to_digits($tok, 9)
            # Now writes one extra digit: %.10f with set_to_digits(tok, 9)
            # ----------------------------------------------------------------
            if re.search(r'^1-e int', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if not quiet:
                    print(f"{line_tokens[0]} {line_tokens[1]} ", end='')
                    print(
                        f"{int(line_tokens[2]):5d}{int(line_tokens[3]):5d} ",
                        end='')
                    print(f"  {set_to_digits(line_tokens[4], 8):17.9f}")
                file_output.write(
                    f"{line_tokens[0]} {line_tokens[1]} "
                    f"{int(line_tokens[2]):5d}{int(line_tokens[3]):5d} "
                    f"  {set_to_digits(line_tokens[4], 9):18.10f}\n")

            # ----------------------------------------------------------------
            # Root singlet/triplet block
            # ----------------------------------------------------------------
            if re.search(r'^  Root ', line) and \
               (re.search(r' singlet ', line) or
                re.search(r' triplet ', line)):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if not quiet:
                    print(f"{line_tokens[0]} {int(line_tokens[1]):d}"
                          f" {line_tokens[2]}", end='')
                    print(f" {set_to_digits(line_tokens[4], 3):0.4f}"
                          f" {line_tokens[5]}", end='')
                    print(f" {set_to_digits(line_tokens[6], 2):0.3f}"
                          f" {line_tokens[7]}")
                file_output.write(
                    f"{line_tokens[0]} {int(line_tokens[1]):d}"
                    f" {line_tokens[2]}"
                    f" {set_to_digits(line_tokens[4], 3):0.4f}"
                    f" {line_tokens[5]}"
                    f" {set_to_digits(line_tokens[6], 2):0.3f}"
                    f" {line_tokens[7]}\n")

            # ----------------------------------------------------------------
            # Root (no singlet/triplet) block
            # ----------------------------------------------------------------
            elif re.search(r'^  Root ', line) and \
                 not (re.search(r' singlet ', line) or
                      re.search(r' triplet ', line)):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if not quiet:
                    print(f"{line_tokens[0]} {int(line_tokens[1]):d}", end='')
                    print(f" {set_to_digits(line_tokens[3], 3):0.4f}"
                          f" {line_tokens[4]}", end='')
                    print(f" {set_to_digits(line_tokens[5], 2):0.3f}"
                          f" {line_tokens[6]}")
                file_output.write(
                    f"{line_tokens[0]} {int(line_tokens[1]):d}"
                    f" {set_to_digits(line_tokens[3], 3):0.4f}"
                    f" {line_tokens[4]}"
                    f" {set_to_digits(line_tokens[5], 2):0.3f}"
                    f" {line_tokens[6]}\n")

            # ----------------------------------------------------------------
            # Zero-Point correction to Energy
            # Original Perl [1]: printf "%.3f", set_to_digits($tok, 3)
            # Now writes one extra digit: %.4f with set_to_digits(tok, 3)
            # ----------------------------------------------------------------
            if re.search(r'Zero-Point correction to Energy', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                for itok in range(num_line_tokens - 5):
                    if not quiet:
                        print(f"{line_tokens[itok]} ", end='')
                    file_output.write(f"{line_tokens[itok]} ")
                itok = num_line_tokens - 5
                if not quiet:
                    print(f"{set_to_digits(line_tokens[itok], 3):.4f}")
                file_output.write(
                    f"{set_to_digits(line_tokens[itok], 3):.4f}\n")

            # ----------------------------------------------------------------
            # Nuclear repulsion energy
            # Original Perl [1]: printf "%.2f", set_to_digits($tok, 2)
            # Now writes one extra digit: %.3f with set_to_digits(tok, 2)
            # ----------------------------------------------------------------
            if re.search(r'nuclear', line) and \
               re.search(r'repulsion', line) and \
               re.search(r'energy', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                for itok in range(num_line_tokens - 1):
                    if not quiet:
                        print(f"{line_tokens[itok]} ", end='')
                    file_output.write(f"{line_tokens[itok]} ")
                itok = num_line_tokens - 1
                if not quiet:
                    print(f"{set_to_digits(line_tokens[itok], 3):.4f}")
                file_output.write(
                    f"{set_to_digits(line_tokens[itok], 2):.3f}\n")

            # ----------------------------------------------------------------
            # Total energy: SCF/DFT/CCSD/MP2/MCSCF/RIMP2/RISCF/BAND/PAW/
            #               WFN1/xTB
            # Original Perl [1]: printf "%.5f", set_to_digits($tok, 5)
            # Now writes one extra digit: %.6f with set_to_digits(tok, 5)
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(r'Total', line) and \
                   re.search(r'energy', line):
                    if (re.search(r'SCF', line) or
                            re.search(r'DFT', line) or
                            re.search(r'CCSD', line) or
                            re.search(r'MP2', line) or
                            re.search(r'MCSCF', line) or
                            re.search(r'RIMP2', line) or
                            re.search(r'RISCF', line) or
                            re.search(r'BAND', line) or
                            re.search(r'PAW', line) or
                            re.search(r'WFN1', line) or
                            re.search(r'xTB', line)):
                        if debug:
                            print(f"\ndebug: {line}", end='')
                        line_tokens = line.split()
                        num_line_tokens = len(line_tokens)
                        if debug:
                            print(f"debug:line_tokens: {line_tokens}")
                            print(f"debug:number     : {num_line_tokens}")
                        for itok in range(num_line_tokens - 1):
                            if not quiet:
                                print(f"{line_tokens[itok]} ", end='')
                            file_output.write(f"{line_tokens[itok]} ")
                        itok = num_line_tokens - 1
                        if not quiet:
                            print(f"{set_to_digits(line_tokens[itok], 5):.6f}")
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 5):.6f}\n")

            # ----------------------------------------------------------------
            # Total energy: PSPW
            # Original Perl [1]: printf "%.5f", set_to_digits($tok, 4)
            # Now writes one extra digit: %.5f with set_to_digits(tok, 4)
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(r'Total', line) and \
                   re.search(r'energy', line):
                    if re.search(r'PSPW', line):
                        if debug:
                            print(f"\ndebug: {line}", end='')
                        line_tokens = line.split()
                        num_line_tokens = len(line_tokens)
                        if debug:
                            print(f"debug:line_tokens: {line_tokens}")
                            print(f"debug:number     : {num_line_tokens}")
                        for itok in range(num_line_tokens - 1):
                            if not quiet:
                                print(f"{line_tokens[itok]} ", end='')
                            file_output.write(f"{line_tokens[itok]} ")
                        itok = num_line_tokens - 1
                        if not quiet:
                            print(f"{set_to_digits(line_tokens[itok], 4):.5f}")
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 4):.5f}\n")

            # ----------------------------------------------------------------
            # total energy: MBPT/LCCD/CCD/LCCSD/CCSD/CCSDT/CCSDTQ/
            #               QCISD/CISD/CISDT/CISDTQ
            # Original Perl [1]: printf "%.7f", set_to_digits($tok, 7)
            # Now writes one extra digit: %.8f with set_to_digits(tok, 7)
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(r'total', line) and \
                   re.search(r'energy', line):
                    if (re.search(r'MBPT', line) or
                            re.search(r'LCCD', line) or
                            re.search(r'CCD', line) or
                            re.search(r'LCCSD', line) or
                            re.search(r'CCSD', line) or
                            re.search(r'CCSDT', line) or
                            re.search(r'CCSDTQ', line) or
                            re.search(r'QCISD', line) or
                            re.search(r'CISD', line) or
                            re.search(r'CISDT', line) or
                            re.search(r'CISDTQ', line)):
                        if debug:
                            print(f"\ndebug: {line}", end='')
                        line_tokens = line.split()
                        num_line_tokens = len(line_tokens)
                        if debug:
                            print(f"debug:line_tokens: {line_tokens}")
                            print(f"debug:number     : {num_line_tokens}")
                        for itok in range(num_line_tokens - 1):
                            if not quiet:
                                print(f"{line_tokens[itok]} ", end='')
                            file_output.write(f"{line_tokens[itok]} ")
                        itok = num_line_tokens - 1
                        if not quiet:
                            print(f"{set_to_digits(line_tokens[itok], 7):.8f}")
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 7):.8f}\n")

            # ----------------------------------------------------------------
            # @GW block
            # Original Perl [1]: printf "%.2f", set_to_digits($tok, 2)
            # Now writes one extra digit: %.3f with set_to_digits(tok, 2)
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(r'@GW', line):
                    if debug:
                        print(f"\ndebug: {line}", end='')
                    line_tokens = line.split()
                    num_line_tokens = len(line_tokens)
                    if debug:
                        print(f"debug:line_tokens: {line_tokens}")
                        print(f"debug:number     : {num_line_tokens}")
                    for itok in range(num_line_tokens - 2):
                        if not quiet:
                            print(f"{line_tokens[itok]} ", end='')
                        file_output.write(f"{line_tokens[itok]} ")
                    itok = num_line_tokens - 2
                    if not quiet:
                        print(f"{set_to_digits(line_tokens[itok], 2):.3f}")
                    file_output.write(
                        f"{set_to_digits(line_tokens[itok], 2):.3f}\n")

            # ----------------------------------------------------------------
            # rt_tddft Etot block
            # Original Perl [1]: printf "%.4f", set_to_digits($tok, 4)
            # Now writes one extra digit: %.5f with set_to_digits(tok, 4)
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(r'rt_tddft', line) and \
                   re.search(r'Etot', line):
                    if debug:
                        print(f"\ndebug: {line}", end='')
                    line_tokens = line.split()
                    num_line_tokens = len(line_tokens)
                    if debug:
                        print(f"debug:line_tokens: {line_tokens}")
                        print(f"debug:number     : {num_line_tokens}")
                    for itok in range(num_line_tokens - 3):
                        if not quiet:
                            print(f"{line_tokens[itok]} ", end='')
                        file_output.write(f"{line_tokens[itok]} ")
                    itok = num_line_tokens - 3
                    if not quiet:
                        print(f"{set_to_digits(line_tokens[itok], 4):.5f}")
                    file_output.write(
                        f"{set_to_digits(line_tokens[itok], 4):.5f}\n")

            # ----------------------------------------------------------------
            # Excitation energy / Rotatory / IBO loc block
            # Original Perl [1]: printf "%.3f", set_to_digits($tok, 3)
            # Now writes one extra digit: %.4f with set_to_digits(tok, 3)
            # ----------------------------------------------------------------
            if (re.search(r'Excitation energy', line) or
                    re.search(r'Rotatory ', line) or
                    re.search(r'IBO loc: largest element in ', line)):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                for itok in range(num_line_tokens - 1):
                    if not quiet:
                        print(f"{line_tokens[itok]} ", end='')
                    file_output.write(f"{line_tokens[itok]} ")
                itok = num_line_tokens - 1
                if not quiet:
                    print(f"{set_to_digits(line_tokens[itok], 3):.4f}")
                file_output.write(
                    f"{set_to_digits(line_tokens[itok], 3):.4f}\n")

            # ----------------------------------------------------------------
            # MR-BWCCSD / BW-MRCCSD / MR-MkCCSD energy block
            # Original Perl [1]: printf "%.10f", set_to_digits($tok, 10)
            # Now writes one extra digit: %.11f with set_to_digits(tok, 10)
            # ----------------------------------------------------------------
            if (re.search(r'MR-BWCCSD energy', line) or
                    (re.search(r'BW-MRCCSD', line) and
                     re.search(r'a posteriori', line)) or
                    re.search(r'MR-MkCCSD energy', line)):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                for itok in range(num_line_tokens - 1):
                    if not quiet:
                        print(f"{line_tokens[itok]} ", end='')
                    file_output.write(f"{line_tokens[itok]} ")
                itok = num_line_tokens - 1
                if not quiet:
                    print(f"{set_to_digits(line_tokens[itok], 10):.11f}")
                file_output.write(
                    f"{set_to_digits(line_tokens[itok], 10):.11f}\n")

            # ----------------------------------------------------------------
            # Isotropic / anisotropy block
            # Original Perl [1]: printf "%.3f", set_to_digits($tok, 3)
            # Now writes one extra digit: %.4f with set_to_digits(tok, 3)
            # Searches backwards for last numeric token to handle lines
            # ending with units labels like 'au'
            # ----------------------------------------------------------------
            if re.search(r'sotropic =', line) or \
               re.search(r'anisotropy =', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")

                # Find the last numeric token, skipping trailing
                # non-numeric tokens such as units labels (e.g. 'au')
                last_numeric_idx = None
                for i in range(num_line_tokens - 1, -1, -1):
                    if _is_float(line_tokens[i]):
                        last_numeric_idx = i
                        break

                if last_numeric_idx is None:
                    if debug:
                        print(f"debug: no numeric token found in "
                              f"isotropic/anisotropy line, skipping:"
                              f" {line}", end='')
                else:
                    for itok in range(last_numeric_idx):
                        if not quiet:
                            print(f"{line_tokens[itok]} ", end='')
                        file_output.write(f"{line_tokens[itok]} ")
                    if not quiet:
                        print(f"{set_to_digits(line_tokens[last_numeric_idx], 3):.4f}")
                    file_output.write(
                        f"{set_to_digits(line_tokens[last_numeric_idx], 3):.4f}\n")

            # ----------------------------------------------------------------
            # average: block
            # Original Perl [1]: printf "%.3f", with 5-token format
            # Now writes one extra digit: %.4f
            # ----------------------------------------------------------------
            if re.search(r' average: ', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if num_line_tokens == 5:
                    # Standard case: "average: <val> + I <val>"
                    # zeros added to avoid -0.000 vs 0.000 diffs
                    if not quiet:
                        print(f"{line_tokens[0]}"
                              f" {float(line_tokens[1]) + 0:.4f}"
                              f" {line_tokens[2]}"
                              f" {line_tokens[3]}"
                              f" {float(line_tokens[4]) + 0:.4f}")
                    file_output.write(
                        f"{line_tokens[0]}"
                        f" {float(line_tokens[1]) + 0:.4f}"
                        f" {line_tokens[2]}"
                        f" {line_tokens[3]}"
                        f" {float(line_tokens[4]) + 0:.4f}\n")
                elif num_line_tokens == 2 and _is_float(line_tokens[1]):
                    # Shorter case: "average: <val>"
                    if not quiet:
                        print(f"{line_tokens[0]}"
                              f" {float(line_tokens[1]) + 0:.4f}")
                    file_output.write(
                        f"{line_tokens[0]}"
                        f" {float(line_tokens[1]) + 0:.4f}\n")
                else:
                    # Fallback: write line as-is to avoid data loss
                    if not quiet:
                        print(line, end='')
                    file_output.write(line)

            # ----------------------------------------------------------------
            # DMX / DMY / DMZ block
            # Original Perl [1]: printf "%.3f", set_to_digits($tok, 3)
            # Now writes one extra digit: %.4f with set_to_digits(tok, 3)
            # ----------------------------------------------------------------
            if (re.search(r'DMX', line) or
                    re.search(r'DMY', line) or
                    re.search(r'DMZ', line)):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if num_line_tokens == 4:
                    if not quiet:
                        print(f"{line_tokens[0]} ", end='')
                        print(f"{set_to_digits(line_tokens[1], 3):.4f}")
                        print(f"{line_tokens[2]} ", end='')
                        print(f"{set_to_digits(line_tokens[3], 3):.4f}")
                    file_output.write(f"{line_tokens[0]} ")
                    file_output.write(
                        f"{set_to_digits(line_tokens[1], 3):.4f}\n")
                    file_output.write(f"{line_tokens[2]} ")
                    file_output.write(
                        f"{set_to_digits(line_tokens[3], 3):.4f}\n")

            # ----------------------------------------------------------------
            # XX/YY/ZZ/XY/XZ/YZ tensor block (not Transition)
            # Uses \b word boundaries to avoid matching substrings like XYZ
            # Original Perl [1]: printf "%.3f", set_to_digits($tok, 3)
            # Now writes one extra digit: %.4f with set_to_digits(tok, 3)
            # ----------------------------------------------------------------
            if re.search(r'\b(XX|YY|ZZ|XY|XZ|YZ)\b', line) and \
               not re.search(r'Transition', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                if num_line_tokens == 4:
                    if not quiet:
                        print(f"{line_tokens[0]} ", end='')
                        for itok in range(1, num_line_tokens):
                            print(
                                f"{set_to_digits(line_tokens[itok], 3):.4f} ",
                                end='')
                        print()
                    file_output.write(f"{line_tokens[0]} ")
                    for itok in range(1, num_line_tokens):
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 3):.4f} ")
                    file_output.write("\n")

            # ----------------------------------------------------------------
            # Gradient block: reading atom/coord/grad data
            # ----------------------------------------------------------------
            if gradient_block == 2:
                if debug:
                    print(f"debug:g3: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:num tok: {num_line_tokens}")
                if num_line_tokens == 8:
                    atoms.append(line_tokens[1])
                    coords.extend(line_tokens[2:5])
                    grads.extend(line_tokens[5:8])
                    if debug:
                        print(f" number of atoms: {len(atoms)} {atoms}")
                        print(f" number of grads: {len(grads)} {grads}")
                        print(f" number of coords: {len(coords)} {coords}")
                else:
                    print("possible bad gradient block")

            # ----------------------------------------------------------------
            # Gradient block: detect header line
            # ----------------------------------------------------------------
            if not sgroup:
                if re.search(
                        r'atom               coordinates'
                        r'                        gradient', line):
                    atoms = []
                    coords = []
                    grads = []
                    gradient_block = 1
                    if debug:
                        print(f"debug:g1: {line}", end='')

            # ----------------------------------------------------------------
            # Gradient block: detect column header line
            # ----------------------------------------------------------------
            if re.search(
                    r'x          y          z'
                    r'           x          y          z', line):
                if debug:
                    print(f"debug:g2: gradient_block is {gradient_block}")
                if gradient_block == 1:
                    gradient_block += 1
                    if debug:
                        print(f"debug:g2: {line}", end='')

            # ----------------------------------------------------------------
            # selci/selci+pt block: reading energy data
            # ----------------------------------------------------------------
            if selcipt_block == 2:
                if debug:
                    print(f"debug:selci get info block: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:num tok: {num_line_tokens}")
                if num_line_tokens == 5:
                    ci_energy.append(float(line_tokens[1]))
                    pt_correc.append(float(line_tokens[2]))
                    cipt_ene.append(float(line_tokens[3]))
                    pt_norm.append(float(line_tokens[4]))
                else:
                    print("possible bad selci or selci+pt energy block")

            # ----------------------------------------------------------------
            # selci+pt block: detect EN| or MP| lines
            # ----------------------------------------------------------------
            if re.search(r'^ EN\|', line) or re.search(r'^ MP\|', line):
                if selcipt_block == 1:
                    ci_energy = []
                    pt_correc = []
                    cipt_ene = []
                    pt_norm = []
                selcipt_block += 1
                if debug:
                    print(f"debug:selcipt inc:{selcipt_block}:"
                          f" line: {line}", end='')

            # ----------------------------------------------------------------
            # Root final energy block
            # Original Perl [1]: printf "%.5f", set_to_digits($tok, 5)
            # Now writes one extra digit: %.6f with set_to_digits(tok, 5)
            # ----------------------------------------------------------------
            if re.search(r'^ Root', line) and \
               re.search(r'final energy', line):
                if debug:
                    print(f"\ndebug: {line}", end='')
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if debug:
                    print(f"debug:line_tokens: {line_tokens}")
                    print(f"debug:number     : {num_line_tokens}")
                for itok in range(num_line_tokens - 1):
                    if itok == 1:
                        if not quiet:
                            print(f"{int(line_tokens[itok]):4d} ", end='')
                        file_output.write(f"{int(line_tokens[itok]):4d} ")
                    else:
                        if not quiet:
                            print(f"{line_tokens[itok]} ", end='')
                        file_output.write(f"{line_tokens[itok]} ")
                itok = num_line_tokens - 1
                if not quiet:
                    print(f"{set_to_digits(line_tokens[itok], 5):.6f}")
                file_output.write(
                    f"{set_to_digits(line_tokens[itok], 5):.6f}\n")

            # ----------------------------------------------------------------
            # DIRDYVTST block: end
            # ----------------------------------------------------------------
            if dirdyv_block and \
               re.search(r'drdy_NWChem has finished', line):
                dirdyv_block = 0

            # ----------------------------------------------------------------
            # DIRDYVTST block: reading data
            # Original Perl [1]: printf "%.5f", set_to_digits($tok, 5)
            # Now writes one extra digit: %.6f with set_to_digits(tok, 5)
            # ----------------------------------------------------------------
            if dirdyv_block:
                line_tokens = line.split()
                num_line_tokens = len(line_tokens)
                if num_line_tokens != 4:
                    file_output.write(f"{line_tokens[0]} ")
                    for itok in range(1, num_line_tokens - 1):
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 5):.6f} ")
                    itok = num_line_tokens - 1
                    file_output.write(
                        f"{set_to_digits(line_tokens[itok], 5):.6f}\n")
                else:
                    for itok in range(num_line_tokens - 1):
                        file_output.write(
                            f"{set_to_digits(line_tokens[itok], 5):.6f} ")
                    itok = num_line_tokens - 1
                    file_output.write(
                        f"{set_to_digits(line_tokens[itok], 5):.6f}\n")

            # ----------------------------------------------------------------
            # DIRDYVTST block: start
            # ----------------------------------------------------------------
            if re.search(
                    r's \(au\)                      '
                    r'frequencies \(cm\^-1\)', line):
                dirdyv_block = 1

        # End of file loop
        if not quiet:
            print(f"nwparse.py: parsed {lines} in file {filename}"
                  f" sent output to {fileout} ")

        file_to_parse.close()
        file_output.close()


if __name__ == '__main__':
    main()
