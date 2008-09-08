#!/usr/bin/env python
from math import *
import sys
import os

name=""
for arg in sys.argv: 
    name=arg

outname=name+'.low'

inf=open(name,'r').readlines()
out=open(outname,'w')

for i in inf:
  x=i.lower()
  out.write(x)

