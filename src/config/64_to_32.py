#!/usr/bin/env python3
"""
Python script to do transliteration from "double" values to "8wrap" values
Usage: python 64_to_32.py file1.f [file2.f ...]

Original Perl script:
    Written:  3/14/97
    By:       Ricky A. Kendall

Converted from Perl to Python using Claude (claude-4-5-opus-4-5-20251101-v1)
Optimizations applied:
    - Early exits to skip unnecessary processing
    - Reduced redundant regex passes
    - Combined line-type checks
    - Precompiled regex patterns
"""

import os
import sys
import re
import subprocess # nosec B404

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
            r'([ ]{5}.)' + escaped + r'(\W)',
            re.IGNORECASE
        )
        
        # Pattern for general substitution (word boundaries)
        self.general_pattern = re.compile(
            r'(\W)' + escaped + r'(\W)',
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
                        # For 64_to_32: tokens[0] is "from", tokens[1] is "to" [1]
                        patterns.append(ConversionPattern(tokens[0], tokens[1]))
    except IOError:
        sys.exit(f"Unable to open: {data_path}")
    
    if debug:
        print(f"Loaded {len(patterns)} conversion patterns")
    
    return patterns


# Precompile line-type detection patterns (used for every line)
FORTRAN_CONTINUATION = re.compile(r'^[ ]{5}[^\s]')


def process_line(line, patterns, line_upper):
    """
    Process a single line with all conversion patterns.
    Uses early exits and optimized regex passes.
    """
    # Early exit: Skip comment lines (Fortran style) or empty lines [1]
    if not line or line[0] in ('c', 'C', '*') or line.strip() == '':
        return line
    
    # Early exit: Quick check if any pattern might match using uppercase comparison
    # This avoids regex overhead for lines with no potential matches
    has_potential_match = False
    for conv in patterns:
        if conv.from_str_upper in line_upper:
            has_potential_match = True
            break
    
    if not has_potential_match:
        return line
    
    # Determine line type once (not for each pattern) [1]
    is_fortran_continuation = line.startswith('     ') and len(line) > 5 and not line[5].isspace()
    is_general_line = len(line) > 0 and line[0] in ' \t'
    
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
        
        # Apply appropriate substitution based on line type (single pass) [1]
        if is_fortran_continuation:
            line = conv.fortran_pattern.sub(
                lambda m: m.group(1) + toot + m.group(2),
                line
            )
        
        if is_general_line:
            line = conv.general_pattern.sub(
                lambda m: m.group(1) + toot + m.group(2),
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
    
    # Runs the Linux 'file' command
    
    result = subprocess.run(['file', file_path], shell=False, capture_output=True, text=True) # nosec

    # DOS files output contains "with CRLF line terminators"
    if "CRLF" in result.stdout:
        print(f"ERROR: DOS file with CRLF {file_path}", file=sys.stderr)
        sys.exit(123)  
    
    # Rename original to backup [1]
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
    
    # Remove backup [1]
    os.unlink(filebak)


def main():
    if debug:
        print(f"Arguments: {sys.argv[1:]}")
    
    # Load patterns once (precompiled)
    patterns = load_data_file()
    
    if len(patterns) == 0:
        sys.exit("Fatal dbl2sngl error: No conversion patterns loaded")
    
    files = sys.argv[1:]
    
    if len(files) == 0:
        print("Usage: python 64_to_32.py file1.f [file2.f ...]")
        sys.exit(1)
    
    # Process each file sequentially
    for file_path in files:
        process_file(file_path, patterns)


if __name__ == "__main__":
    main()
