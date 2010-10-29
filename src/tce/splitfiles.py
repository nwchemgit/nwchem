# Usage: python splitfiles.py < inputfile.F
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)
# $Id$
 
import string
import copy
import sys

source = sys.stdin.readlines()

if (not source):
   print "Usage: python splitfiles.py < inputfile.F"

nfiles = 0
for line in source:
   if (string.find(line,"SUBROUTINE") != -1):
      if (nfiles):
         file = open(filename+".F","w")
         for newline in filecontent:
            file.write(newline)
      filename = string.split(line[string.find(line,"SUBROUTINE")+11:],"(")[0]
      print filename+".o\\"
      nfiles = nfiles + 1
      filecontent = [line]
   else:
      filecontent.append(line)
# don't forget to dump the last subroutine
file = open(filename+".F","w")
for newline in filecontent:
   file.write(newline)
print "Number of files generated:",nfiles
