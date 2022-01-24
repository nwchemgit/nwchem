# Usage: python splitfiles.py < inputfile.F
# (c) All rights reserved by Battelle & Pacific Northwest Nat'l Lab (2002)
# $Id$

import sys

source = sys.stdin.readlines()

if (not source):
    print("Usage: python splitfiles.py < inputfile.F")

nfiles = 0
filename = " "
filecontent = " "
for line in source:
   if (line.find("SUBROUTINE") != -1):
      if (nfiles):
         file = open(filename+".F", "w")
         for newline in filecontent:
            file.write(newline)
      filename = line[line.find("SUBROUTINE")+11:].split("(", 99999)[0]
      print(filename+".o\\")
      nfiles = nfiles + 1
      filecontent = [line]
   else:
      filecontent.append(line)
# don't forget to dump the last subroutine
file = open(filename+".F", "w")
for newline in filecontent:
   file.write(newline)
print("Number of files generated:", nfiles)
