from nwchem import *
from math import *

def geom_get_coords(name):
  #
  # This routine returns a list with the cartesian
  # coordinates in user input units for the geometry
  # of given name
  #
  try:
    actualname = rtdb_get(name)
    if (actualname == None):
      actualname = name
  except NWChemError:
    actualname = name
  if actualname is None:
    actualname = name
  coords = rtdb_get('geometry:' + actualname + ':coords')
  units = rtdb_get('geometry:'+actualname+':user units')
  if (units == 'a.u.'):
    factor = 1.0
  elif (units == 'angstroms'):
    factor = rtdb_get('geometry:' + str(actualname) + ':angstrom_to_au')
  else:
    raise NWChemError('unknown units')
  i = 0
  while (i < len(coords)):
    coords[i] = coords[i] / factor
    i = i + 1
  return coords

def geom_set_coords(name,coords):
  #
  # This routine, given a list with the cartesian
  # coordinates in user input units set them in 
  # the geometry of given name.
  #
  try:
    actualname = rtdb_get(name)
  except NWChemError:
    actualname = name
  if actualname is None:
    actualname = name
  units = rtdb_get('geometry:'+actualname+':user units')
  if (units == 'a.u.'):
    factor = 1.0
  elif (units == 'angstroms'):
    factor = rtdb_get('geometry:' + str(actualname) + ':angstrom_to_au')
  else:
    raise NWChemError('unknown units')
  coords = list(coords)
  i = 0
  while (i < len(coords)):
    coords[i] = coords[i] * factor
    i = i + 1
  rtdb_put('geometry:' + str(actualname) + ':coords',coords)

def bond_length(i,j):  # atoms numbered 1,2,...
  #
  # Return the distance betwen atoms i and j in user
  # units in the default geometry
  #
  coords = geom_get_coords('geometry')
  x = coords[(i-1)*3  ]-coords[(j-1)*3  ]
  y = coords[(i-1)*3+1]-coords[(j-1)*3+1]
  z = coords[(i-1)*3+2]-coords[(j-1)*3+2]
  return sqrt(x*x + y*y + z*z)

def minimize1d(f, xlo, xhi, xtol, maxeval):
  #
  # Find the minimum value of function(x) in [xlo,xhi]
  # to a precision in x of xtol.  Maxeval specifies
  # the maximum no. of function evaluations.
  #
  # If you want to maximize f() then minimize -f()
  #
  # It returns (xmin, fmin) where 
  #    xmin = position of lowest value computed
  #    fmin = f(xmin)
  # Also it guarantees that the last point evaluated
  # was xmin so any external state is consistent.

  gold = 0.38197
  neval = 0

  if (xhi < xlo):
    tmp = xhi
    xhi = xlo
    xlo = tmp      

  if (ga_nodeid() == 0):
     print('   Mode     xlo       xmid      xhi        flo          fmid          fhi')
     print(' ------- --------- --------- --------- ------------ ------------ ------------')

  xmid = xlo + (xhi - xlo)*gold
  flo  = f(xlo)
  if (ga_nodeid() == 0):
     print(' startup%10.4f                    %13.6f'   % (xlo, flo))

  fmid = f(xmid)
  if (ga_nodeid() == 0):
     print(' startup%10.4f%10.4f          %13.6f%13.6f' % (xlo, xmid, flo, fmid))

  fhi  = f(xhi)
  neval = neval + 3

  xlast = xhi    # Tracks last point of function evaluation

  # First bracket the minimum
  while (not ((fmid<flo) and (fmid<fhi))):
    if (ga_nodeid() == 0):
       print(' bracket%10.4f%10.4f%10.4f%13.6f%13.6f%13.6f' % (xlo, xmid, xhi, flo, fmid, fhi))
    if (neval >= maxeval):
       raise NWChemError('min1d: too many evaluations')
    if ((flo>fmid) and (fmid>fhi)):
       xlo  = xmid
       flo  = fmid
    elif ((fhi>fmid) and (fmid>flo)):
       xhi  = xmid
       fhi  = fmid
    elif ((fmid>flo) and (fmid>fhi)):
       if (flo < fhi):
          xhi  = xmid
          fhi  = fmid
       else:
          xlo  = xmid
          flo  = fmid
    else:
       raise NWChemError('unanticipated')
    xmid = xlo + (xhi - xlo)*gold
    fmid = f(xmid)
    neval = neval + 1

  # The minimum is now bracketed.  
  # Parabolic fit or Golden section search.
  mode = 'search'
  dx = xhi-xlo
  while (dx > xtol):
     if (ga_nodeid() == 0):
         print(' %s %10.4f%10.4f%10.4f%13.6f%13.6f%13.6f' % (mode, xlo, xmid, xhi, flo, fmid, fhi))
     if (neval >= maxeval):
         raise NWChemError('min1d: too many evaluations')

     # Try a parabolic fit
     d1 = (fmid - flo) / (xmid-xlo)
     d2 = (fhi  - flo) / (xhi -xlo)
     a  = (d2 - d1) / (xhi - xmid)
     b  = d1 - a*(xmid-xlo)
     c  = flo
     if (a == 0.0):
       a = -1.0
     xtest  = xlo + -b / (2.0*a)
     ftestp = a*(xtest-xlo)*(xtest-xlo) +  b*(xtest-xlo) + c
     ftest = 0.0
     if ((a > 0) and (b < 0) and (xtest>xlo) and (xtest<xhi)):
        mode = 'newton'
        ftest = f(xtest)
        neval = neval + 1
        if ((ftest > fmid)):
           print(' Rejecting Newton step due to uphill motion %13.6f' % ftest)
           mode = 'search'
        elif (xtest > xmid):
           xlo = xmid
           flo = fmid
           xmid= xtest
           fmid= ftest
        else:
           xhi = xmid
           fhi = fmid
           xmid= xtest
           fmid= ftest

     if (mode == 'search'):
        mode = 'search'
        if ((xhi - xmid) > (xmid - xlo)):
            xtest = xmid + (xhi-xmid)*gold
            ftest = f(xtest)
            neval = neval + 1
            if (ftest < fmid):
                xlo  = xmid
                flo  = fmid
                xmid = xtest
                fmid = ftest
            else:
                xhi = xtest
                fhi = ftest
        else:
            xtest = xmid + (xlo-xmid)*gold
            ftest = f(xtest)
            neval = neval + 1
            if (ftest < fmid):
                xhi  = xmid
                fhi  = fmid
                xmid = xtest
                fmid = ftest
            else:
                xlo = xtest
                flo = ftest
     dx = fabs(xlast-xtest) 
     xlast = xtest
    
  print(' %s %10.4f%10.4f%10.4f%13.6f%13.6f%13.6f' % (mode, xlo, xmid, xhi, flo, fmid, fhi))
  if (xlast != xmid):
     print(' Re-evaluating at final coordinate')
     fmid = f(xmid)

  return (xmid, fmid)


def scan_input(input,start,end,nstep,theory,task):
  # 
  # Scan some NWChem input in nsteps from the parameters in
  # start[] to those in end[].  The parameters are substituted
  # in order into input.  Theory selects the electronic
  # wavefunction and task can be any defined task but most
  # likely is one of task_energy or task_optimize.  If task
  # optimize is selected, then it makes most sense for one
  # or more of the parameters to be frozen in a geometry. 
  #
  # Calculations are NOT performed at the end points
  # so specifying nstep=1 does one calculation at the mid point.
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
  #  results = scan_input(geom, \
  #                      [0.967, 103.3], 
  #                      [2.109,  26.96],
  #                      10, 'scf', task_energy)

  results = []
 
  if (len(start) != len(end)):
    raise NWChemError('scan_input: inconsistent #parameters')

  if (ga_nodeid() == 0): 
    print(' ')
    print(' Scanning NWChem input ')
    print(' ---------------------')
    print(' ')
    print(input)
    print(' ')
    print(' Nstep ', nstep)
    print(' Start ', start)
    print(' End   ', end)
    print(' ')

  for i in range(1,nstep+1):
     alpha = (1.0*i)/(nstep+1)
     new = []
     for j in range(0,len(start)):
        new.append((1.0-alpha)*start[j] + alpha*end[j])
       
     if (ga_nodeid() == 0):
        print(' ')
        print(' Scanning NWChem input - step %d of %d ' % (i,nstep))
        print(' ')
        print(input % tuple(new))
        print(' ')

     input_parse(input % tuple(new))
     result = task(theory)
     if (ga_nodeid() == 0):
        print(' ')
        print(' Scanning NWChem input - results from step ', i)
        print(' ')
        print(result)
        print(' ')
     results.append((new,result));

  return tuple(results)
