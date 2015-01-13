README
======

The script

   doi.py

is a little Python tool that takes one or more files. Each of the files is 
assumed to be a list of DOI-s with one DOI per line. The script compiles a
list of all unique DOI-s, queries the CrossRef infrastructure to obtain 
the full bibliographic information for each, and finally prints the BibTeX 
bibliography file on the standard output. 

This tool requires:
	PycURL   - For the curl functionality
	ArgParse - For command line parsing

As an example try running

   ./doi.py test.doi

This should print a BibTeX list with two references that is identical to the
contents of test.bib.

Note that this tool is extremely limited at the moment. In particular it 
does not do anything clever if there are network problems, for example.
