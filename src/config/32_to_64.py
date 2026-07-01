#!/usr/bin/env python3
"""
Python script to do transliteration from "single" values to "double" values
Usage: python 32_to_64.py file1.f [file2.f ...]
Original Perl script:
    Written:  3/14/97
    By:       Ricky A. Kendall
    High Performance Computational Chemistry Group
    Theory Modeling and Simulation Program [2]

Converted from Perl to Python using Claude (claude-4-5-opus-4-5-20251101-v1)
Optimizations applied:
    - Early exits to skip unnecessary processing
    - Reduced redundant regex passes
    - Combined line-type checks
    - Precompiled regex patterns
    - Lookahead for correct underscore boundary matching
"""

import os
import sys
import re

debug = False


def copy_case(from_str, to_str):
    """Take case from 'from_str' and apply it to 'to_str' and return that new string"""
    result = []
    for i, char in enumerate(to_str):
        if i < len(from_str):
            if from_str[i].isupper():
                result.append(char.upper())
            elif from_str[i].islower():
                result.append(char.lower())
            else:
                result.append(char)
        else:
            result.append(char)
    return ''.join(result)


class ConversionPattern:
    """Holds precompiled regex patterns for a single from/to conversion pair"""
    
    def __init__(self, from_str, to_str):
        self.from_str = from_str
        self.to_str = to_str
        self.from_str_upper = from_str.upper()
        # Precompile all regex patterns
        escaped = re.escape(from_str)
        
        # Pattern for finding the from_str (case insensitive)
        self.search_pattern = re.compile(escaped, re.IGNORECASE)
        
        # Pattern for Fortran continuation lines (5 spaces + non-space)
        self.fortran_pattern = re.compile(
            r'([ ]{5}.)' + escaped + r'(?=_|\W)',
            re.IGNORECASE
        )

        # Pattern for general substitution
        # Use capturing group for preceding \W (consumed, restored via group(1))
        # Use lookahead (?=_|\W) so underscore after pattern is not consumed [2]
        # This correctly handles: void ycopy_() -> void dcopy_()
        self.general_pattern = re.compile(
            r'(\W)' + escaped + r'(?=_|\W)',
            re.IGNORECASE
        )


def load_data_file():
    """Load the from/to conversion pairs from the data file and precompile patterns"""
    patterns = []
    
    data_path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(data_path, "data.64_to_32")
    
    if debug:
        print(f"Data path: {data_path}")
    
    try:
        with open(data_path, 'r') as data_file:
            for line in data_file:
                if line and not line.startswith('#'):
                    tokens = line.split()
                    if len(tokens) >= 2:
                        # For 32_to_64: tokens[1] is "from", tokens[0] is "to" [2]
                        patterns.append(ConversionPattern(tokens[1], tokens[0]))
    except IOError:
        sys.exit(f"Unable to open: {data_path}")
    
    if debug:
        print(f"Loaded {len(patterns)} conversion patterns")
    
    return patterns


def process_line(line, patterns, line_upper):
    """
    Process a single line with all conversion patterns.
    Uses early exits and optimized regex passes [2]
    """
    # Early exit: Skip comment lines (Fortran style) or empty lines [2]
    if not line or line[0] in ('c', 'C', '*') or line.strip() == '':
        return line

    # Early exit: Quick check if any pattern might match using uppercase comparison
    has_potential_match = any(conv.from_str_upper in line_upper for conv in patterns)
    if not has_potential_match:
        return line

    # Determine line type once (not for each pattern) [2]
    is_fortran_continuation = (
        line.startswith('     ') and len(line) > 5 and not line[5].isspace()
    )
    # Matches lines starting with space or tab [2]
    is_general_line = len(line) > 0 and line[0] in ' \t'
    # Matches any line starting with space or non-whitespace char [2]
    # ^[ \S] in Perl matches virtually every non-empty line
    # including C declarations like "void ycopy_()"
    is_declaration_line = len(line) > 0 and (line[0] == ' ' or not line[0].isspace())

    # Process each conversion pattern
    for conv in patterns:
        # Early exit: Skip if pattern not in line (case-insensitive quick check)
        if conv.from_str_upper not in line_upper:
            continue

        match = conv.search_pattern.search(line)
        if not match:
            continue

        # Found a match - compute replacement
        froom = match.start()
        toot = copy_case(
            line[froom:froom + len(conv.from_str)],
            conv.to_str
        )

        # Apply Fortran continuation substitution [2]
        if is_fortran_continuation:
            line = conv.fortran_pattern.sub(
                lambda m: m.group(1) + toot,
                line
            )

        # Apply general substitution for tab/space lines [2]
        if is_general_line:
            line = conv.general_pattern.sub(
                lambda m: m.group(1) + toot,
                line
            )

        # Apply declaration substitution for C-style lines [2]
        # Handles cases like "void ycopy_()" where line starts with non-space
        if is_declaration_line:
            line = conv.general_pattern.sub(
                lambda m: m.group(1) + toot,
                line
            )

        # Update line_upper after substitution for subsequent patterns
        line_upper = line.upper()

    return line


def process_file(file_path, patterns):
    """Process a single file and perform transliteration"""
    pid = os.getpid()
    orgfile = file_path
    filebak = f"{file_path}.{pid}"
    
    if debug:
        print(f"Processing: {file_path}")
        print(f"Backup file: {filebak}")
    
    # Rename original to backup [2]
    os.rename(file_path, filebak)
    
    try:
        with open(filebak, 'r') as fh_in:
            with open(orgfile, 'w') as fh_out:
                for line in fh_in:
                    line_upper = line.upper()
                    processed_line = process_line(line, patterns, line_upper)
                    fh_out.write(processed_line)
    except IOError as e:
        sys.exit(f"Can't open file: {e}")
    
    # Remove backup [2]
    os.unlink(filebak)


def main():
    if debug:
        print(f"Arguments: {sys.argv[1:]}")
    
    # Load patterns once (precompiled)
    patterns = load_data_file()
    
    if len(patterns) == 0:
        sys.exit("Fatal sngl2dbl error: No conversion patterns loaded")
    
    files = sys.argv[1:]

    if len(files) == 0:
        print("Usage: python 32_to_64.py file1.f [file2.f ...]")
        sys.exit(1)

    # Process each file sequentially
    for file_path in files:
        process_file(file_path, patterns)


if __name__ == "__main__":
    main()
