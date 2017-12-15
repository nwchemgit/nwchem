import sys
import numpy as np
import scipy.io

def skip_n(k):
  for i in range(k-1):
    sys.stdin.readline()
  return sys.stdin.readline()

line = sys.stdin.readline()
# initial skip
while line:
  if line.startswith('frame'): break
  line = sys.stdin.readline()

frames = []
while line:
  f = []
  line = skip_n(6)
  while line:
    if line.startswith('frame'): break
    f.append(map(float, line.split()))
    line = sys.stdin.readline()
  frames.append(f)
trace = np.array(frames)
scipy.io.savemat('dump.mat', mdict={'trace' : trace})
