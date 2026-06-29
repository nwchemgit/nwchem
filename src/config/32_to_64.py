#!/usr/bin/env python3
"""
Python script to do transliteration from "single" values to "double" values
Usage: python 32_to_64.py file1.f [file2.f ...]
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
                        # Note: For 32_to_64, tokens are reversed compared to 64_to_32
                        # tokens[1] is "from" and tokens[0] is "to" [2]
                        patterns.append(ConversionPattern(tokens[1], tokens[0]))
    except IOError:
        sys.exit(f"Unable to open: {data_path}")
    
    if debug:
        print(f"Loaded {len(patterns)} conversion patterns")
    
    return patterns


# Precompile line-type detection patterns (used for every line)
COMMENT_PATTERN = re.compile(r'^[cC*]')
FORTRAN_CONTINUATION = re.compile(r'^[ ]{5}[^\s]')
TAB_START = re.compile(r'^[ \t]')
DECLARATION_START = re.compile(r'^[ \S]')
DIGIT_START = re.compile(r'^[ \d]')


def process_file(file_path, patterns):
    """Process a single file and perform transliteration"""
    pid = os.getpid()
    orgfile = file_path
    filebak = f"{file_path}.{pid}"
    
    if debug:
        print(f"Processing: {file_path}")
        print(f"Backup file: {filebak}")
    
    # Rename original to backup
    os.rename(file_path, filebak)
    
    try:
        with open(filebak, 'r') as fh_in:
            with open(orgfile, 'w') as fh_out:
                for line in fh_in:
                    # Skip comment lines (Fortran style) or empty lines
                    if COMMENT_PATTERN.match(line) or line.strip() == '':
                        fh_out.write(line)
                        continue
                    
                    # Process each conversion pattern
                    for conv in patterns:
                        match = conv.search_pattern.search(line)
                        if match:
                            # Find where the string starts
                            froom = match.start()
                            toot = copy_case(
                                line[froom:froom + len(conv.from_str)],
                                conv.to_str
                            )
                            
                            # Fortran continuation lines (5 spaces + non-space)
                            if FORTRAN_CONTINUATION.match(line):
                                line = conv.fortran_pattern.sub(
                                    lambda m: m.group(1) + toot + m.group(2),
                                    line
                                )
                            
                            # Handle tab chars in C [2]
                            if TAB_START.match(line):
                                line = conv.general_pattern.sub(
                                    lambda m: m.group(1) + toot + m.group(2),
                                    line
                                )
                            
                            # Handle declarations in C [2]
                            if DECLARATION_START.match(line):
                                line = conv.general_pattern.sub(
                                    lambda m: m.group(1) + toot + m.group(2),
                                    line
                                )
                            
                            # Handle digit starts
                            if DIGIT_START.match(line):
                                line = conv.general_pattern.sub(
                                    lambda m: m.group(1) + toot + m.group(2),
                                    line
                                )
                    
                    fh_out.write(line)
                    
    except IOError as e:
        sys.exit(f"Can't open file: {e}")
    
    # Remove backup
    os.unlink(filebak)


def main():
    if debug:
        print(f"Arguments: {sys.argv[1:]}")
    
    # Load patterns once (precompiled)
    patterns = load_data_file()
    
    if len(patterns) == 0:
        sys.exit("Fatal sngl2dbl error: No conversion patterns loaded")
    
    # Process each file
    for file_path in sys.argv[1:]:
        process_file(file_path, patterns)


if __name__ == "__main__":
    main()
