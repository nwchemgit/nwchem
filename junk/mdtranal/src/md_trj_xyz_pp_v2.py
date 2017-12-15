import sys
import numpy as np
import scipy.io

ne = 168 # number of elements, TODO: parameterize later
fh = "  {}".format(ne) # frame header

def skip_n(k):
  for i in range(k-1):
    sys.stdin.readline()
  return sys.stdin.readline()

line = sys.stdin.readline()
# initial skip
while line:
  if line.startswith(fh): break
  line = sys.stdin.readline()

frames = []
while line:
  f = []
  line = skip_n(2) # d onot understand this line
  while line:
    if line.startswith(fh): break
    f.append(map(float, line.split()[1:4]))
    line = sys.stdin.readline()
  print f
  frames.append(f)
trace = np.array(frames)
scipy.io.savemat('dump.mat', mdict={'trace' : trace})
