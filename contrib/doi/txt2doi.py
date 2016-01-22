#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vi:ts=4:et

# Collect DOI-s from text files and write the list of unique DOI-s.
# DOI-s are indicated by the preceeding "DOI:" or "doi:" string.

from argparse import ArgumentParser
import sys
import io
import re
try:
    from io import BytesIO
except ImportError:
    from StringIO import StringIO as BytesIO

def parse_args():
    """Parse the command line arguments"""
    parser = ArgumentParser(description="Read text files and extract DOI-s")
    parser.add_argument("filename",help="A text file with embedded DOI-s", nargs="+")
    args = parser.parse_args()
    return vars(args)

def parse_files(filenames):
    """Go through every file and add DOI-s to the doi-table."""
    doi_table = []
    #nfiles = len(filenames)
    #if nfiles < 1:
    #    raise Exception("Need at least 1 file name to read DOI-s from.")
    expr = re.compile("[dD][oO][iI]:")
    for filename in filenames:
        try:
            file = io.open(filename,'r',encoding='utf-8')
            file_content = file.read()
            file.close()
            # message = filename + ": utf-8"
            # print(message)
        except UnicodeDecodeError:
            file = io.open(filename,'r',encoding='utf-16')
            file_content = file.read()
            file.close()
            file_content = file_content.replace(u"\ufeff"," ")
            # message = filename + ": utf-16"
            # print(message)
        file_tokens = file_content.split()
        length = len(file_tokens)
        ii = 0
        while ii < length:
           token = file_tokens[ii]
           if expr.search(token):
               # print(token)
               # print(file_tokens[ii+1])
               if len(token) > 4:
                   # Token is e.g. doi:10.1016/S0883-2927(03)00115-X
                   # i.e. not space between "doi:" and the actual DOI.
                   token = token[4:]
               else:
                   # Token is "doi:" so the actual DOI is in the next token
                   ii += 1
                   token = file_tokens[ii]
               # message = filename + ": " + token
               # print(message)
               doi_table.append(token)
           ii += 1
        # print(doi_table)
    return doi_table

def remove_duplicates(doi_table):
    """Remove all duplicate DOI entries. We exploit that each element in a set
    must be unique"""
    doi_set = set(doi_table)
    doi_table = list(doi_set)
    doi_table.sort()
    return doi_table

def print_dois(doi_table):
    """Print the list of DOI-s."""
    for doi in doi_table:
        print( doi )

def main ():
    """Run the whole thing (limited capability at the moment)"""
    opts = []
    args = []
    doi_table = []
    args = parse_args()
    doi_table = parse_files(args["filename"])
    doi_table = remove_duplicates(doi_table)
    print_dois(doi_table)

if __name__ == "__main__":
    main()

