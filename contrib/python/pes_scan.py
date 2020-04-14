from nwchem import *
from math import *

def pes_scan(input,start,end,nstep,theory,task):
  # 
  # Does a true multidimensional potential energy surface
  # scan, with user input specifying the minimum and maximum
  # values of each variable, and number of times to step each.
  # The number of calculations done SCALES QUICKLY; according
  # to (nstep+1)^(# of variables to scan).  This code can handle
  # any computationally feasible number of variables, including
  # just one.
  #
  # Theory selects the electronic wavefunction and task can be
  # any defined task but most likely is one of task_energy or
  # task_optimize.  If task optimize is selected, then it makes
  # most sense for one or more of the parameters to be frozen
  # in a geometry. 
  #
  # Calculations ARE performed at the end points so specifying
  # nstep=1 does one calculation at the start point, and one at
  # the end point.
  #
  # Returned is the tuple
  #
  # [(param-1,results-1), (param-2,results-2), ...]
  #
  # where param-n are the parameters of the n-th step and
  # results-n are the corresponding results returned by task()
  #
  # Example.  Scan a bond & angle computing the SCF energy
  #
  #  geom = '''
  #    geometry noprint adjust
  #      zcoord
  #        bond  3 2   %f oh
  #        angle 3 2 1 %f hon constant
  #      end
  #    end
  #  '''
  #  results = pes_scan(geom, \
  #                      [0.967, 103.3], 
  #                      [2.109,  26.96],
  #                      10, 'scf', task_energy)
  #
  # in this example using pes_scan, over 120 single point energy
  # calculations would be done.  using scan_input, only 10 would
  # be done.

  results = []
 
  if (len(start) != len(end)):
    raise NWChemError('pes_scan: inconsistent #parameters')

  npoint = (nstep+1)**len(start)

  if (ga_nodeid() == 0): 
    print(' ')
    print(' Doing a PES Scan on input ')
    print(' -------------------------')
    print(' ')
    print(input)
    print(' ')
    print(' Number of points ', npoint)
    print(' Minimum values ', start)
    print(' Maximum values ', end)
    print(' ')

  step = []
  for i in range (0, len(start)):
    step.append((end[i]-start[i])/nstep)


  ylist = []
  for i in range (0, len(start)):
    xlist=[]
    for j in range (0, nstep+1):
      xlist.append(start[i]+j*step[i])
    ylist.append(xlist)

  
  zlist = []
  indexes = len(ylist) * [0]
  base = len(ylist[0])

  while 1:
    elt = []
    for i in range(0, len(indexes)):
        elt.append(ylist[i][indexes[i]])
    zlist.append(elt)
    for i in range(len(indexes)-1, -1, -1):
        if indexes[i] < base-1:
            indexes[i] = indexes[i]+1
            for j in range(i+1, len(indexes)):
                indexes[j] = 0
            break
    else: break


  for i in range(0,len(zlist)):
     new = zlist[i]
    
       
     if (ga_nodeid() == 0):
        print(' ')
        print(' Scanning NWChem input - point %d of %d ' % (i+1,npoint))
        print(' ')
        print(input % tuple(new))
        print(' ')

     input_parse(input % tuple(new))
     result = task(theory)
     if (ga_nodeid() == 0):
        print(' ')
        print(' Scanning NWChem input - results from point ', i+1)
        print(' ')
        print(result)
        print(' ')
     results.append((new,result));

  if (ga_nodeid() == 0):
     print(' ')
     print(' Python Scan Output ')
     print(' ')
     for i in range(0,len(results)):
       print(results[i][0], results[i][1])
     print(' ')
     print(' Python Scan Output Finished ')
  return tuple(results)
