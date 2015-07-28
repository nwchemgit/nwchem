#!/usr/bin/python -E
#
# -E  Suppresses all PYTHON* environment variables. In particular
#     PYTHONHOME is a source of disaster.
#
#
# === README ===
#
# This program replaces function calls of density functionals with the
# corresponding results of a density functional subroutine call. In addition
# the appropriate subroutine call is also spliced in at the right place.
#
#
# BACKGROUND
# ----------
#
# In the normal operation of Maxima generating source code one would specify
# the expression for the density functional in full. The resulting (typically
# very large) expression is then differentiated in every required way, the
# catalogue of expressions is optimized and everything is transformed into 
# Fortran. This leads to very large, illegible, subroutines that may take 
# hours to generate, and tens of minutes to compile (e.g. nwxcm_c_tpss03_d3
# is one subroutine of 20857 lines and 1.25 MB in size). Clearly this is 
# undesirable.
#
# In practice the functionals implemented in these routines can be constructed
# from other functionals. For example, the TPSS correlation functional
# mentioned above is expressed in terms of a modified PBE functional, which
# in turn is specified in terms of the PW91 LDA functional. When hand coding
# functionals this structure is maintained to make functionals manageable. 
# Ideally we would like to preserve this structure while using symbolic algebra
# to generate the source code.
#
# Using symbolic algebra it is in principle possible to do this as one can 
# specify a function in terms of other functions that Maxima does not know.
# For example, one can specify a function as
#
#    f(x,y) := (x+y)^2*h(y);
#
# where Maxima does not know what h(y) is. Differentiating f(x,y) wrt. y
# Maxima will generate
#
#    diff(f(x,y),y) := 2*(x+y)*h(y)+(x+y)^2*'diff(h(y),y,1)
#
# where 'diff indicates that a derivative is required. The Fortran generated
# by Maxima will reference h(y) and 'diff(h(y),y,1) which is not necessarily
# valid Fortran. If we know what h(y) and 'diff(h(y),y,1) stand for this problem
# can resolved by replacing these entities with valid Fortran functions or
# variables. In that fashion Maxima can be used to generate source code while
# preserving the structure of the original functional implementation. 
#
# This script implements these source code transformations.
#
#
# APPROACH
# --------
#
# The main complication in this program is that the expression optimization 
# that Maxima performs messes things up as in the example above it would 
# break out the references to h(y) as in
#
#    t1 = h(y)
#    diff(f(x,y),y) := 2*(x+y)*t1+(x+y)^2*'diff(t1,y,1)
#
# Also note that h may be invoked multiple times as in, for example,
#
#    f(x,y) := (x+y)^2*h(x,y) - x^2*h(x,0) - y^2*h(0,y)
#
# so we need to distinguish every invokation based on the arguments passed in.
# In addition, in our case, h is actually implemented as a subroutine that
# returns the function value and the values of the various derivatives. The
# inputs and outputs are stored in arrays. The steps required are:
#
# 1. All Fortran lines need to be unwrapped so that simple text analysis 
#    can reliably find references to h(x,y), 'diff(h(x,y)), etc.
# 2. Where h(x,y) is invoked we need to generate arrays x and y, as well as
#    arrays for the required derivatives and initialize them.
# 3. We need to insert the array declarations in the type declarations.
# 4. We need to insert the subroutine calls at the appropriate places (i.e.
#    right before the first place where results are used as we need to make 
#    sure that the input data is defined).
# 5. Where h(x,y) or its derivatives are referenced we need to substitute 
#    the appropriate variables.
# 6. The Fortran lines need to be rewrapped such that we use at most 72 columns.
#
import re
import sys
import string

type_autoxc   = 1
type_autoxcDs = 2

def usage(code):
   """
   Print usage information
   """
   sys.stderr.write("Insert subroutine calls to functionals in Maxima generated code.")
   sys.stderr.write("")
   sys.stderr.write("  %s [-h] < <filein> > <fileout>"%sys.argv[0])
   sys.stderr.write("")
   sys.stderr.write("-h         Print this information")
   sys.stderr.write("<filein>   Raw Fortran from Maxima")
   sys.stderr.write("<fileout>  Fortran from Maxima with subroutine calls")
   sys.stderr.write("           for functionals")
   exit(code)

def unwrap_lines(lines_in):
   """
   Take a list of lines of Fortran and unwrap all continuation lines. 
   The result is stored in a new list which is returned.
   All the input lines end with a '\n' character, the output lines have no
   newline characters.
   """
   pattern = re.compile("     [0-9:;<=>?@+]")
   lines_out = []
   longline = ""
   for line in lines_in:
      length = len(line)
      if pattern.match(line):
         shortline = line[6:length-1]
         longline += shortline.lstrip()
      elif length > 1:
         if len(longline) > 0:
            lines_out.append(longline)
         longline = line[:length-1]
   if len(longline) > 0:
      lines_out.append(longline)
   return lines_out

def find_subroutine(lines,lineno):
   """
   Given the starting line number (lineno) and a list of lines find the 
   line numbers where the next subroutine starts and ends. Here we use
   Fortran90 coding styles so the starting and end points can be identified
   by looking for
   - subroutine
   - end subroutine
   A tuple with the corresponding starting and end line numbers is returned.
   """
   length = len(lines)
   line = lineno
   lineno_start = -1
   lineno_end   = -1
   pattern = re.compile("      subroutine")
   while !pattern.match(lines[line]):
      line += 1
      if line > length:
         break
   if line <= length:
      lineno_start = line
   pattern = re.compile("      end subroutine")
   while !pattern.match(lines[line]):
      line += 1
      if line > length:
         break
   if line <= length:
      lineno_end = line
   return (lineno_start,lineno_end)

def find_code_skeleton_type(lines,subr_lines):
   """
   Scan the lines of a subroutine to find out whether it is an autoxc or and
   autoxc-Ds generated subroutine.
   """
   (lineno_start,lineno_end) = subr_lines
   pattern = re.compile("if (taua.gt.tol_rho) then")
   subr_type = type_autoxc
   line = lineno_start
   while line <= lineno_end:
      if pattern.match(lines[line]):
         subr_type = type_autoxcDs
         break
      line += 1
   return subr_type

def find_autoxc_code_skeleton(lines,subr_lines):
   """
   Assuming that the subroutine is of the autoxc code skeleton type find all
   the if-branches and return a list of all the (start,end) tuples.
   """
   (lineno_start,lineno_end) = subr_lines
   ifbranches = []
   line = lineno_start
   ifstart = -1
   ifend   = -1
   #
   pattern = re.compile("if (rhoa.gt.tol_rho) then")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifstart = line+1
         break
      line += 1
   #
   pattern = re.compile("endif ! rhoa.gt.tol_rho")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifend   = line-1
         break
      line += 1
   #
   ifbranches.append((ifstart,ifend))
   #
   ifstarta = -1
   ifenda   = -1
   ifstartb = -1
   ifendb   = -1
   ifstartc = -1
   ifendc   = -1
   pattern = re.compile("if (rhoa.gt.tol_rho.and.rhob.gt.tol_rho) then")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifstarta = line+1
         break
      line += 1
   #
   pattern = re.compile("elseif (rhoa.gt.tol_rho.and.rhob.le.tol_rho) then")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifenda   = line-1
         ifstartb = line+1
         break
      line += 1
   #
   pattern = re.compile("elseif (rhoa.le.tol_rho.and.rhob.gt.tol_rho) then")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifendb   = line-1
         ifstartc = line+1
         break
      line += 1
   #
   pattern = re.compile("endif ! rhoa.gt.tol_rho.and.rhob.gt.tol_rho")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifendc   = line-1
         break
      line += 1
   #
   ifbranches.append((ifstarta,ifenda))
   ifbranches.append((ifstartb,ifendb))
   ifbranches.append((ifstartc,ifendc))
   #
   return ifbranches

def find_autoxcDs_code_skeleton(lines,subr_lines):
   """
   Assuming that the subroutine is of the autoxc-Ds code skeleton type find all
   the if-branches and return a list of all the (start,end) tuples.
   """
   (lineno_start,lineno_end) = subr_lines
   ifbranches = []
   line = lineno_start
   ifstarta = -1
   ifenda   = -1
   ifstartb = -1
   ifendb   = -1
   #
   pattern = re.compile("if (taua.gt.tol_rho) then")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifstarta = line+1
         break
      line += 1
   #
   pattern = re.compile("else")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifenda   = line-1
         ifstartb = line+1
         break
      line += 1
   #
   pattern = re.compile("endif")
   while line <= lineno_end:
      if pattern.match(lines[line])
         ifendb  = line-1
         break
      line += 1
   #
   ifbranches.append((ifstarta,ifenda))
   ifbranches.append((ifstartb,ifendb))
   #
   for x in range(0,3):
      ifstarta = -1
      ifenda   = -1
      ifstartb = -1
      ifendb   = -1
      ifstartc = -1
      ifendc   = -1
      ifstartd = -1
      ifendd   = -1
      pattern = re.compile("if (taua.gt.tol_rho.and.taub.gt.tol_rho) then")
      while line <= lineno_end:
         if pattern.match(lines[line])
            ifstarta = line+1
            break
         line += 1
      #
      pattern = re.compile("elseif (taua.gt.tol_rho.and.taub.le.tol_rho) then")
      while line <= lineno_end:
         if pattern.match(lines[line])
            ifenda   = line-1
            ifstartb = line+1
            break
         line += 1
      #
      pattern = re.compile("elseif (taua.le.tol_rho.and.taub.gt.tol_rho) then")
      while line <= lineno_end:
         if pattern.match(lines[line])
            ifendb   = line-1
            ifstartc = line+1
            break
         line += 1
      #
      pattern = re.compile("else")
      while line <= lineno_end:
         if pattern.match(lines[line])
            ifendc   = line-1
            ifstartd = line+1
            break
         line += 1
      #
      pattern = re.compile("endif")
      while line <= lineno_end:
         if pattern.match(lines[line])
            ifendd   = line-1
            break
         line += 1
      #
      ifbranches.append((ifstarta,ifenda))
      ifbranches.append((ifstartb,ifendb))
      ifbranches.append((ifstartc,ifendc))
      ifbranches.append((ifstartd,ifendd))
      #
   return ifbranches

def find_type_declaration_insertion_point(lines,subr_lines):
   """
   Find the point in the subroutine where the array declarations for the
   subroutine calls can be inserted.
   """
   (lineno_start,lineno_end) = subr_lines
   pattern = re.compile("#include \"nwxc_param.fh\"")
   line_insert = -1
   line = lineno_start
   while line <= lineno_end:
      if pattern.match(lines[line])
         line_insert = line
         break
      line += 1
   return line_insert

def collect_subroutine_calls(lines,ifbranch_lines):
   """
   Collect a dictionary of all subroutine call instances. In the source code
   on entry the "subroutine call" is given in the form of a function call, e.g.
   t5 = nwxc_c_Mpbe(rhoa,0.0d+0,gammaaa,0.0d+0,0.0d+0)
   we need to know the "nwxc_c_Mpbe(rhoa,0.0d+0,gammaaa,0.0d+0,0.0d+0)" part.
   The key in the dictionary is going to be the variable name (i.e. "t5") as we
   will have to replace those variable instances with an array element 
   reference. The resulting dictionary is returned.
   """
   (lineno_start,lineno_end) = ifbranch_lines
   dict = {}
   pattern = re.compile("nwxc")
   line = lineno_start
   while line <= lineno_end:
      if pattern.match(lines[line]):
         aline = lines[line]
         aline.split(" = ")
         dict[aline[0]] = aline[1]
      line += 1
   return dict

def delete_lines(lines,ifbranches):
   """
   Collect a list of lines that should be deleted (i.e. skipped) when writing
   the output subroutine. The lines in question are the ones generated from the
   "at" command, e.g.
   t10(1) = (gammaaa = gammaaa)
   t10(2) = (gammaab = gammaaa)
   t10(3) = (gammabb = gammaaa)
   t10(4) = (rhoa = rhoa)
   t10(5) = (rhob = rhoa)
   t10(6) = (taua = taua)
   t10(7) = (taub = taua)
   The list of line numbers is returned.
   """
   dlist = []
   pattern = re.compile("nwxc")
   for ifbranch in ifbranches:
      (lineno_start,lineno_end) = ifbranch
      line = lineno_start
      while line <= lineno_end:
         aline = lines[line]
         aline.split(" = ")
         if len(aline) == 3:
            # line: t10(2) = (gammaab = gammaaa)
            dlist.append(line)
         line += 1
   return dlist

def find_maxno_calls(lines,ifbranches):
   """
   In every if-branch some subroutine calls will be inserted. For the inputs
   we need one set of arrays for rho, gamma and tau. However, for the outputs
   we need a separate set of arrays for each separate call as the results may
   appear in multiple place throughout the remained of the if-branch. In 
   addition there is no guarantee that the results are not used in overlapping
   code segments. Therefore we need to count how many different sets of 
   output variables we need to declare. This routine returns the resulting
   number.
   """
   varsets = 0
   for ifbranch in ifbranches:
      dict = collect_subroutine_calls(lines,ifbranch)
      numcalls = len(dict)
      if numcalls > varsets:
         varsets = numcalls
   return varsets

def append_declarations(olines,varsets,orderdiff):
   """
   Given the list of output lines so far, the number of output variable sets
   and the maximum order of differentiation append the required array
   declarations. The new list of lines is returned.
   """
   olines.append("      double precision sr(NCOL_RHO)")
   olines.append("      double precision sg(NCOL_GAMMA)")
   olines.append("      double precision st(NCOL_TAU)")
   for ii in range(1,varsets+1):
      line = "      double precision s"+str(ii)+"f"
      olines.append(line)
      line = "      double precision s"+str(ii)+"a(NCOL_AMAT)"
      olines.append(line)
      line = "      double precision s"+str(ii)+"c(NCOL_CMAT)"
      olines.append(line)
      line = "      double precision s"+str(ii)+"m(NCOL_MMAT)"
      olines.append(line)
      if orderdiff > 1:
         line = "      double precision s"+str(ii)+"a2(NCOL_AMAT2)"
         olines.append(line)
         line = "      double precision s"+str(ii)+"c2(NCOL_CMAT2)"
         olines.append(line)
         line = "      double precision s"+str(ii)+"m2(NCOL_MMAT2)"
         olines.append(line)
      if orderdiff > 2:
         line = "      double precision s"+str(ii)+"a3(NCOL_AMAT3)"
         olines.append(line)
         line = "      double precision s"+str(ii)+"c3(NCOL_CMAT3)"
         olines.append(line)
         line = "      double precision s"+str(ii)+"m3(NCOL_MMAT3)"
         olines.append(line)
   return olines

def find_max_order_diff(lines,subr_lines):
   """
   Work out the maximum order of differentiation for a given subroutine.
   The maximum order is returned.
   """
   (lineno_start,lineno_end) = subr_lines
   aline = lines[lineno_start]
   orderdiff = 0
   pattern_d2 = re.compile("_d2(")
   pattern_d3 = re.compile("_d3(")
   if pattern_d2.match(aline):
      orderdiff = 2
   elif pattern_d3.match(aline):
      orderdiff = 3
   else:
      orderdiff = 1
   return orderdiff


   

if len(sys.argv) == 2:
   if sys.argv[1] == "-h":
      usage(0)
   else:
      usage(1)
elif len(sys.argv) > 2:
   usage(1)

lines = sys.stdin.readlines()
lines = unwrap_lines(lines)
for line in lines:
   sys.stdout.write("%s\n"%line)
