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

ifbranch_closedshell = 1
ifbranch_openshell   = 2

func_invalid  = -1
func_lda      = 1
func_gga      = 2
func_mgga     = 3

version = "$Id$"
version = version.split()
revision = version[1]+" revision "+version[2]+" "+version[3]

def usage(code):
   """
   Print usage information
   """
   sys.stderr.write("Insert subroutine calls to functionals in Maxima generated code.")
   sys.stderr.write("")
   sys.stderr.write("  %s [-h] [-v|--version] < <filein> > <fileout>"%sys.argv[0])
   sys.stderr.write("")
   sys.stderr.write("-h            Print this information")
   sys.stderr.write("-v|--version  Print the version data")
   sys.stderr.write("<filein>      Raw Fortran from Maxima")
   sys.stderr.write("<fileout>     Fortran from Maxima with subroutine calls")
   sys.stderr.write("              for functionals")
   sys.stderr.write("")
   sys.stderr.write("$Id$")
   sys.exit(code)

def var_to_int(var):
   """
   Convert a variable name (such as "t5") to an integer (in this case 5).
   Variables like "Amat(iq,D1_RA)" should never show up here.
   """
   return int(var[1:])

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
   while not pattern.match(lines[line]):
      line += 1
      if line >= length:
         break
   if line < length:
      lineno_start = line
   # next statement needed to guarantee that line is in the valid range
   if line >= length:
      line = length-1
   pattern = re.compile("      end subroutine")
   while not pattern.match(lines[line]):
      line += 1
      if line >= length:
         break
   if line < length:
      lineno_end = line
   return (lineno_start,lineno_end)

def find_code_skeleton_type(lines,subr_lines):
   """
   Scan the lines of a subroutine to find out whether it is an autoxc or and
   autoxc-Ds generated subroutine.
   """
   (lineno_start,lineno_end) = subr_lines
   pattern = re.compile("if \(taua\.gt\.tol_rho\) then")
   subr_type = type_autoxc
   line = lineno_start
   while line <= lineno_end:
      if pattern.search(lines[line]):
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
   pattern = re.compile("if \(rhoa\.gt\.tol_rho\) then")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifstart = line+1
         break
      line += 1
   #
   pattern = re.compile("endif ! rhoa\.gt\.tol_rho")
   while line <= lineno_end:
      if pattern.search(lines[line]):
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
   pattern = re.compile("if \(rhoa\.gt\.tol_rho\.and\.rhob\.gt\.tol_rho\) then")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifstarta = line+1
         break
      line += 1
   #
   pattern = re.compile("elseif \(rhoa\.gt\.tol_rho\.and\.rhob\.le\.tol_rho\) then")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifenda   = line-1
         ifstartb = line+1
         break
      line += 1
   #
   pattern = re.compile("elseif \(rhoa\.le\.tol_rho\.and\.rhob\.gt\.tol_rho\) then")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifendb   = line-1
         ifstartc = line+1
         break
      line += 1
   #
   pattern = re.compile("endif ! rhoa\.gt\.tol_rho\.and\.rhob\.gt\.tol_rho")
   while line <= lineno_end:
      if pattern.search(lines[line]):
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
   pattern = re.compile("if \(taua\.gt\.tol_rho\) then")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifstarta = line+1
         break
      line += 1
   #
   pattern = re.compile("else")
   while line <= lineno_end:
      if pattern.search(lines[line]):
         ifenda   = line-1
         ifstartb = line+1
         break
      line += 1
   #
   pattern = re.compile("endif")
   while line <= lineno_end:
      if pattern.search(lines[line]):
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
      pattern = re.compile("if \(taua\.gt\.tol_rho\.and\.taub\.gt\.tol_rho\) then")
      while line <= lineno_end:
         if pattern.search(lines[line]):
            ifstarta = line+1
            break
         line += 1
      #
      pattern = re.compile("elseif \(taua\.gt\.tol_rho\.and\.taub\.le\.tol_rho\) then")
      while line <= lineno_end:
         if pattern.search(lines[line]):
            ifenda   = line-1
            ifstartb = line+1
            break
         line += 1
      #
      pattern = re.compile("elseif \(taua\.le\.tol_rho\.and\.taub\.gt\.tol_rho\) then")
      while line <= lineno_end:
         if pattern.search(lines[line]):
            ifendb   = line-1
            ifstartc = line+1
            break
         line += 1
      #
      # needed because "else" is a substring of "elseif ..."
      line = ifstartc
      pattern = re.compile("else")
      while line <= lineno_end:
         if pattern.search(lines[line]):
            ifendc   = line-1
            ifstartd = line+1
            break
         line += 1
      #
      pattern = re.compile("endif")
      while line <= lineno_end:
         if pattern.search(lines[line]):
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
      if pattern.match(lines[line]):
         line_insert = line
         break
      line += 1
   return line_insert

def find_subroutine_call_insertion_point(lines,ifbranch,varname):
   """
   The assumption is that the symbolic algebra optimization will always break
   out the function evaluation. Therefore there always is a line where the
   functional value is assigned to a variable. The actual functional 
   subroutine call, at the latest, has to happen on the line before.
   This routine returns the line where the variable first appears.
   """
   (lineno_start,lineno_end) = ifbranch
   line = lineno_start
   insert_point = -1
   while line <= lineno_end:
      aline = lines[line]
      aline = aline.split(" = ")
      key = aline[0].lstrip()
      if key == varname:
         insert_point = line
         break
      line += 1
   return insert_point

def collect_subroutine_calls(lines,ifbranch_lines):
   """
   Collect a dictionary of all subroutine call instances. In the source code
   on entry the "subroutine call" is given in the form of a function call, e.g.
   t5 = nwxc_c_Mpbe(rhoa,0.0d+0,gammaaa,0.0d+0,0.0d+0)
   we need to know the "nwxc_c_Mpbe(rhoa,0.0d+0,gammaaa,0.0d+0,0.0d+0)" part.
   The key in the dictionary is going to be the variable name (i.e. "t5") as we
   will have to replace those variable instances with an array element 
   reference.
   The resulting dictionary is returned.
   """
   (lineno_start,lineno_end) = ifbranch_lines
   dict = {}
   pattern = re.compile("nwxc")
   line = lineno_start
   while line <= lineno_end:
      if pattern.search(lines[line]):
         aline = lines[line]
         aline = aline.split(" = ")
         key   = aline[0].lstrip()
         dict[key] = aline[1]
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
   for ifbranch in ifbranches:
      (lineno_start,lineno_end) = ifbranch
      line = lineno_start
      while line <= lineno_end:
         aline = lines[line]
         aline = aline.split(" = ")
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

def append_subroutine_call(olines,subrname,arglist,orderdiff,ifbranch_kind,funckind,num,indent):
   """
   Given the list of output lines so far, the subroutine name, the input
   argument list, the order of differentiation, the functional kind as well as
   the number of the output variable set, generate the subroutine call code.
   The subroutine call code involves initializing the input and output 
   arguments, and constructing the actual call itself.
   The list of output lines is returned.
   """
   #DEBUG
   #print "append_subroutine_call: arglist:",arglist
   #DEBUG
   if ifbranch_kind == ifbranch_closedshell:
      if arglist[1] == '0.0d+0' or arglist[2] == '0.0d+0':
         # We are dealing with the open shell term of the Stoll partitioning
         # of the correlation energy
         line = indent+"sr(R_A) = "+arglist[1]
         olines.append(line)
         line = indent+"sr(R_B) = "+arglist[2]
         olines.append(line)
         if funckind >= func_gga:
            line = indent+"sg(G_AA) = "+arglist[3]
            olines.append(line)
            line = indent+"sg(G_AB) = "+arglist[4]
            olines.append(line)
            line = indent+"sg(G_BB) = "+arglist[5]
            olines.append(line)
         if funckind >= func_mgga:
            line = indent+"st(T_A) = "+arglist[6]
            olines.append(line)
            line = indent+"st(T_B) = "+arglist[7]
            olines.append(line)
      else:
         # We are dealing with a regular closed shell call
         line = indent+"sr(R_T) = 2.0d0*"+arglist[1]
         olines.append(line)
         if funckind >= func_gga:
            line = indent+"sg(G_TT) = 4.0d0*"+arglist[3]
            olines.append(line)
         if funckind >= func_mgga:
            line = indent+"st(T_T) = 2.0d0*"+arglist[6]
            olines.append(line)
   elif ifbranch_kind == ifbranch_openshell:
      # We are dealing with a regular open shell call
      line = indent+"sr(R_A) = "+arglist[1]
      olines.append(line)
      line = indent+"sr(R_B) = "+arglist[2]
      olines.append(line)
      if funckind >= func_gga:
         line = indent+"sg(G_AA) = "+arglist[3]
         olines.append(line)
         line = indent+"sg(G_AB) = "+arglist[4]
         olines.append(line)
         line = indent+"sg(G_BB) = "+arglist[5]
         olines.append(line)
      if funckind >= func_mgga:
         line = indent+"st(T_A) = "+arglist[6]
         olines.append(line)
         line = indent+"st(T_B) = "+arglist[7]
         olines.append(line)
   else:
      sys.stderr.write("append_subroutine_call: invalid ifbranch_kind: %d\n"%ifbranch_kind)
      sys.exit(20)
   line = indent+"s"+str(num)+"f = 0.0d0"
   olines.append(line)
   line = indent+"call dcopy(NCOL_AMAT,0.0d0,0,s"+str(num)+"a,1)"
   olines.append(line)
   if orderdiff >= 2:
      line = indent+"call dcopy(NCOL_AMAT2,0.0d0,0,s"+str(num)+"a2,1)"
      olines.append(line)
   if orderdiff >= 3:
      line = indent+"call dcopy(NCOL_AMAT3,0.0d0,0,s"+str(num)+"a3,1)"
      olines.append(line)
   if funckind >= func_gga:
      line = indent+"call dcopy(NCOL_CMAT,0.0d0,0,s"+str(num)+"c,1)"
      olines.append(line)
      if orderdiff >= 2:
         line = indent+"call dcopy(NCOL_CMAT2,0.0d0,0,s"+str(num)+"c2,1)"
         olines.append(line)
      if orderdiff >= 3:
         line = indent+"call dcopy(NCOL_CMAT3,0.0d0,0,s"+str(num)+"c3,1)"
         olines.append(line)
   if funckind >= func_mgga:
      line = indent+"call dcopy(NCOL_MMAT,0.0d0,0,s"+str(num)+"m,1)"
      olines.append(line)
      if orderdiff >= 2:
         line = indent+"call dcopy(NCOL_MMAT2,0.0d0,0,s"+str(num)+"m2,1)"
         olines.append(line)
      if orderdiff >= 3:
         line = indent+"call dcopy(NCOL_MMAT3,0.0d0,0,s"+str(num)+"m3,1)"
         olines.append(line)
   line = indent+"call "+subrname
   if orderdiff == 2:
      line = line+"_d2"
   elif orderdiff == 3:
      line = line+"_d3"
   line = line+"("+arglist[0]+",tol_rho,"
   if arglist[1] == '0.0d+0' or arglist[2] == '0.0d+0':
      # We are dealing with the Stoll partitioning of the correlation energy
      line = line+"2"
   else:
      line = line+"ipol"
   line = line+",1,1.0d0,sr"
   if funckind >= func_gga:
      line = line+",sg"
   if funckind >= func_mgga:
      line = line+",st"
   line = line+",s"+str(num)+"f,s"+str(num)+"a"
   if orderdiff >= 2:
      line = line+",s"+str(num)+"a2"
   if orderdiff >= 3:
      line = line+",s"+str(num)+"a3"
   if funckind >= func_gga:
      line = line+",s"+str(num)+"c"
      if orderdiff >= 2:
         line = line+",s"+str(num)+"c2"
      if orderdiff >= 3:
         line = line+",s"+str(num)+"c3"
   if funckind >= func_mgga:
      line = line+",s"+str(num)+"m"
      if orderdiff >= 2:
         line = line+",s"+str(num)+"m2"
      if orderdiff >= 3:
         line = line+",s"+str(num)+"m3"
   line = line+")"
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
   pattern_d2 = re.compile("_d2\(")
   pattern_d3 = re.compile("_d3\(")
   if pattern_d2.search(aline):
      orderdiff = 2
   elif pattern_d3.search(aline):
      orderdiff = 3
   else:
      orderdiff = 1
   return orderdiff

def find_functional_kind(funccall):
   """
   Given the way the functional was invoked as a function work out whether
   the functional is an LDA, GGA or Meta-GGA functional. We establish this
   by counting the number of arguments which is 3 for LDA, 6 for GGA and 
   8 for a Meta-GGA functional. The numbers include the "param" argument
   as well as the density dependent quantities, i.e. for LDA we have 3
   arguments: param + rhoa + rhob.
   The functional kind is returned.
   """
   data = funccall
   data = data.split(",")
   length = len(data)
   funckind = func_invalid
   if   length == 3:
      funckind = func_lda
   elif length == 6:
      funckind = func_gga
   elif length == 8:
      funckind = func_mgga
   return funckind

def make_functional_name(funccall):
   """
   Given the way the functional was invoked as a function work out what the
   name of the subroutine to call is.
   The subroutine name is returned.
   """
   data = funccall
   data = data.split("(")
   data = data[0]
   data = data.replace("nwxc_","nwxcm_")
   return data

def make_input_args_list(funccall):
   """
   Given the way the functional was invoked as a function work out what the
   argument list is. This list will be used to initialize the input data to
   the actual subroutine call.
   The list of arguments is returned.
   """
   data = funccall.split("(")
   data = data[1].split(")")
   data = data[0].split(",")
   return data

def find_varname(dict,diffstr):
   """
   Given the subroutine call dictionary and the diffstr representing the
   "'diff(...)" command work out which variable is indicated. If
   diffstr is e.g. 'diff(t5,gammaa,2,rhob,1) then we return e.g. 
   s3c3(D3_RB_GAA_GAA).
   If diffstr is just a variable name, e.g. t5. then it represents the
   functional value rather than a derivative and we return e.g. s3f.
   The indicated variable is returned as a string.
   """
   data = diffstr
   #DEBUG
   #print "find_varname: data:",data
   #DEBUG
   if "%at(" == data[:4]:
      data = data[4:-1]
      iend = data.rfind(",")
      data = data[:iend]
      #DEBUG
      #print "find_varname: at:",data
      #DEBUG
   if "'diff(" == data[:6]:
      data = data[6:-1]
      #DEBUG
      #print "find_varname: diff:",data
      #DEBUG
   data = data.split(",")
   callref = data[0]
   list = dict.keys()
   list = sorted(list,key=var_to_int)
   num = -1
   if len(list) > 0:
      num = 0
      #DEBUG
      #print "find_varname:",callref,list
      #DEBUG
      while callref != list[num]:
         num += 1
      #DEBUG
      #print "find_varname: num:",num
      #DEBUG
   lengthl = len(list)
   if num == lengthl:
      sys.stdout.write("entity %s not found\n"%callref)
      for jj in range(0,length):
         sys.stdout.write("list %d: %s\n"%(jj,list[jj]))
      sys.exit(10)
   # num is now the variable set number
   lengthd = len(data)
   #DEBUG
   #print "find_varname: lengthd:",lengthd
   #DEBUG
   if lengthd == 1:
      # this is the energy functional value
      var_name = "s"+str(num+1)+"f"
      return var_name
   orderdiff = 0
   if lengthd >= 3:
      orderdiff += int(data[2])
   if lengthd >= 5:
      orderdiff += int(data[4])
   if lengthd >= 7:
      orderdiff += int(data[6])
   #DEBUG
   #print "find_varname: orderdiff:",orderdiff
   #DEBUG
   # orderdiff is now the order of differentiation
   patternc = re.compile("gamma")
   patternt = re.compile("tau")
   var_func = func_lda
   if lengthd >= 3:
      if patternc.match(data[1]):
         var_func = func_gga
      if patternt.match(data[1]):
         var_func = func_mgga
   if lengthd >= 5:
      if var_func < func_gga and patternc.match(data[3]):
         var_func = func_gga
      if var_func < func_mgga and patternt.match(data[3]):
         var_func = func_mgga
   if lengthd >= 7:
      if var_func < func_gga and patternc.match(data[5]):
         var_func = func_gga
      if var_func < func_mgga and patternt.match(data[5]):
         var_func = func_mgga
   #DEBUG
   #print "find_varname: var_func:",var_func
   #DEBUG
   var_char = "invalid"
   if var_func == func_lda:
      var_char = "a"
   elif var_func == func_gga:
      var_char = "c"
   elif var_func == func_mgga:
      var_char = "m"
   # var_char is now the character representing the output matrix, 
   # a for Amat, c for Cmat, and m for Mmat
   var_field = "D"+str(orderdiff)
   if lengthd >= 3:
      if "rhoa" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_RA"
   if lengthd >= 5:
      if "rhoa" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_RA"
   if lengthd >= 7:
      if "rhoa" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_RA"
   #
   if lengthd >= 3:
      if "rhob" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_RB"
   if lengthd >= 5:
      if "rhob" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_RB"
   if lengthd >= 7:
      if "rhob" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_RB"
   #
   if lengthd >= 3:
      if "gammaaa" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_GAA"
   if lengthd >= 5:
      if "gammaaa" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_GAA"
   if lengthd >= 7:
      if "gammaaa" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_GAA"
   #
   if lengthd >= 3:
      if "gammaab" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_GAB"
   if lengthd >= 5:
      if "gammaab" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_GAB"
   if lengthd >= 7:
      if "gammaab" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_GAB"
   #
   if lengthd >= 3:
      if "gammabb" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_GBB"
   if lengthd >= 5:
      if "gammabb" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_GBB"
   if lengthd >= 7:
      if "gammabb" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_GBB"
   #
   if lengthd >= 3:
      if "taua" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_TA"
   if lengthd >= 5:
      if "taua" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_TA"
   if lengthd >= 7:
      if "taua" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_TA"
   #
   if lengthd >= 3:
      if "taub" == data[1]:
         for ii in range(0,int(data[2])):
            var_field = var_field + "_TB"
   if lengthd >= 5:
      if "taub" == data[3]:
         for ii in range(0,int(data[4])):
            var_field = var_field + "_TB"
   if lengthd >= 7:
      if "taub" == data[5]:
         for ii in range(0,int(data[6])):
            var_field = var_field + "_TB"
   # var_field now contains the array field, e.g. D1_RA, D3_GAA_TB_TB, etc.
   var_name = "s"+str(num+1)+var_char
   if orderdiff >= 2:
      var_name = var_name+str(orderdiff)
   var_name = var_name+"("+var_field+")"
   return var_name

def line_contains_var(line,var):
   """
   Check whether a variable specified by "var" is contained within the line
   specified by "line". E.g. var="t2" and line="a=b+t2*t4" should return
   True whereas line="a=b+t21" should return False.
   """
   regular  = var+"[^0-9]"
   pattern  = re.compile(regular)
   longline = line+" "
   result   = pattern.match(longline)
   return (result != None)

def find_indent(line):
   """
   Given a line of source code work the indentation out and return a string
   containing as many spaces as the indentation.
   """
   pattern = re.compile("\s*")
   obj     = pattern.match(line)
   indent  = obj.group()
   return indent

def find_var_in_line(line,var_name):
   """
   Find all locations where the variable given by "var_name" is used.
   The locations are stored in a list which is returned.
   """
   patterne = re.compile(" = ")
   # If var_name = t2 we need to make sure we do not replace
   # cmat2 as well. Patternw is needed to achieve that. So we look for
   # t2 as well at2, if the end-points are equal the string found does
   # not match the variable t2.
   patternv = re.compile(var_name+"[^0-9]")
   patternw = re.compile("a"+var_name+"[^0-9]")
   aline = line+" "
   found = patterne.search(aline)
   list = []
   if found:
      pos = found.end(0)
      obj = patternv.search(aline,pos)
      while obj:
         objw = patternw.search(aline,pos)
         if objw:
            if objw.end(0) != obj.end(0):
               list.append((obj.start(0),obj.end(0)-1))
         else:
            list.append((obj.start(0),obj.end(0)-1))
         pos = obj.end(0)
         obj = patternv.search(aline,pos)
   return list

def expand_var_in_line(line,list):
   """
   The variable references are given in list in the form of tuples so that
   line(begin:end) matches the variable exactly. However, the string we need
   to replace might be larger, e.g. in the case we have an entity such as
   'diff(t5,rhoa,1,gammaaa,2), or worse %at('diff(t5,rhoa,1),t10). 
   Here we expand the tuples in the list to cover these strings in full.
   The list with updated tuples is returned.
   """
   pattern = re.compile("\)")
   item = 0
   length = len(list)
   while item < length:
      (ibegin,iend) = list[item]
      iibegin = ibegin-6
      if line[iibegin:ibegin] == "'diff(":
         iiend = pattern.search(line,iend).end(0)
         list[item] = (iibegin,iiend)
      (ibegin,iend) = list[item]
      iibegin = ibegin-4
      if line[iibegin:ibegin] == "%at(":
         iiend = pattern.search(line,iend).end(0)
         list[item] = (iibegin,iiend)
      item += 1
   return list

def find_replace_var_in_range(lines,iline_begin,iline_end,dict,var):
   """
   Find all locations in the range of lines (ibegin,iend) where the variable
   "var" is referenced and replace those references with the proper functional
   output variable.
   The updated list of lines "lines" is returned.
   """
   iline = iline_begin
   while iline <= iline_end:
      list = find_var_in_line(lines[iline],var)
      list = expand_var_in_line(lines[iline],list)
      aline = lines[iline]
      bline = ""
      oend = 0
      for tuple in list:
         (ibegin,iend) = tuple
         diffstr = aline[ibegin:iend]
         #DEBUG
         #print "diffstr:",diffstr
         #DEBUG
         var_name = find_varname(dict,diffstr)
         #DEBUG
         #print "var_name:",var_name
         #DEBUG
         bline = bline+aline[oend:ibegin]+var_name
         oend = iend
      bline = bline+aline[oend:]
      lines[iline] = bline
      iline += 1
   return lines

def rewrap_line(longline):
  """
  Break a given long line "longline" up into 72 character long chunks that
  conform the Fortran77 standard. The chunks are written to standard output.
  In addition we do not want to break the line in the middle of numbers.
  """
  pattern = re.compile("\S")
  i = (pattern.search(longline)).start()
  indent = longline[:i]
  indent = indent[:5]+"+"+indent[7:]+"   "
  while len(longline) > 72:
    i = -1
    # wrap before * / ( ) + or -
    i = max(i,string.rfind(longline,",",0,70)+1)
    i = max(i,string.rfind(longline,"*",0,71))
    i = max(i,string.rfind(longline,"/",0,71))
    i = max(i,string.rfind(longline,"(",0,71))
    i = max(i,string.rfind(longline,")",0,71))
    # wrap before + but not in the middle of a numerical constant...
    j = string.rfind(longline,"+",0,71)
    k = string.rfind(longline,"d+",0,71)
    if j-1 == k:
      j = string.rfind(longline,"+",0,k)
    i = max(i,j)
    # wrap before - but not in the middle of a numerical constant...
    j = string.rfind(longline,"-",0,71)
    k = string.rfind(longline,"d-",0,71)
    if j-1 == k:
      j = string.rfind(longline,"-",0,k)
    i = max(i,j)
    if i == -1:
      sys.stderr.write("No sensible break point found in:\n")
      sys.stderr.write(longline)
      sys.exit(1)
    elif i == 6:
      sys.stderr.write("Same break point found repeatedly in:\n")
      sys.stderr.write(longline)
      sys.exit(1)
    sys.stdout.write(longline[:i]+"\n")
    longline = indent + longline[i:]
  sys.stdout.write(longline+"\n")


if len(sys.argv) == 2:
   if sys.argv[1] == "-h":
      usage(0)
   elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
      sys.stdout.write("%s\n"%revision)
      sys.exit(0)
   else:
      usage(1)
elif len(sys.argv) > 2:
   usage(1)

ilines = sys.stdin.readlines()
ilines = unwrap_lines(ilines)
nlines = len(ilines)
#DEBUG
#file = open("junkjunk",'w')
#for line in ilines:
#   file.write("%s\n"%line)
#file.close()
#DEBUG
olines = []
line_start = 0
subr_lines = find_subroutine(ilines,line_start)
(subr_lines_start,subr_lines_end)=subr_lines
#DEBUG
#print "subrs: start,end:",subr_lines_start,subr_lines_end
#DEBUG
while subr_lines_start != -1 and subr_lines_end != -1:
   #
   # Roll forward to the beginning of the subroutine
   #
   while line_start < subr_lines_start:
      olines.append(ilines[line_start])
      line_start += 1
   #
   # Work the subroutine type out and then work the if-branches out
   #
   code_skel_type = find_code_skeleton_type(ilines,subr_lines)
   if   code_skel_type == type_autoxc:
      ifbranches = find_autoxc_code_skeleton(ilines,subr_lines)
   elif code_skel_type == type_autoxcDs:
      ifbranches = find_autoxcDs_code_skeleton(ilines,subr_lines)
   else:
      sys.stderr.write("Unexpected code_type\n")
      exit(20)
   #
   # Work the order of differentiation out
   #
   orderdiff = find_max_order_diff(ilines,subr_lines)
   #DEBUG
   #print "orderdiff:",orderdiff
   #DEBUG
   #
   # Work the declaration insertion point out
   #
   type_decl_insertion_point = find_type_declaration_insertion_point(ilines,subr_lines)
   #
   # How many additional variables do we need?
   #
   max_calls = find_maxno_calls(ilines,ifbranches)
   #DEBUG
   #print "max_calls:",max_calls
   #DEBUG
   #
   # Which lines do we need to drop?
   #
   delete_lines_list = delete_lines(ilines,ifbranches)
   #DEBUG
   #print "delete_list: ",delete_lines_list
   #DEBUG
   #
   # We know what additional variables we need to declare. 
   # So copy more lines from the input file to the output until we reach
   # the declaration insertion point. Then call the function to the inject
   # the additional declarations in the output routine.
   #
   while line_start <= type_decl_insertion_point:
      olines.append(ilines[line_start])
      line_start += 1
   olines = append_declarations(olines,max_calls,orderdiff)
   #
   # Go through branches and mess with subroutine calls
   #
   i_ifbranch = 0
   for ifbranch in ifbranches:
      i_ifbranch += 1
      if   code_skel_type == type_autoxc and i_ifbranch <= 1:
         ifbranch_kind = ifbranch_closedshell
      elif code_skel_type == type_autoxcDs and i_ifbranch <= 2:
         ifbranch_kind = ifbranch_closedshell
      else:
         ifbranch_kind = ifbranch_openshell
      (line_if_start,line_if_end) = ifbranch
      call_lines = collect_subroutine_calls(ilines,ifbranch)
      #DEBUG
      #print "calls: start,end:",line_if_start,line_if_end
      #print "calls: call_lines:",call_lines
      #DEBUG
      call_vars = call_lines.keys()
      #DEBUG
      #print "calls: call_vars :",call_vars
      #DEBUG
      call_vars = sorted(call_vars,key=var_to_int)
      num = 0
      for var in call_vars:
         num += 1
         #DEBUG
         #print "looping over vars:",num,var
         #DEBUG
         call_insert = find_subroutine_call_insertion_point(ilines,ifbranch,var)
         indent = find_indent(ilines[call_insert])
         ilines = find_replace_var_in_range(ilines,call_insert,line_if_end,call_lines,var)
         # var_name = find_varname()
         #DEBUG
         #print "ilines:",call_insert,":",ilines[call_insert]
         #DEBUG
         ilines[call_insert] = indent+var+" = s"+str(num)+"f"
         # Roll forward to just before the insertion point
         while line_start < call_insert:
            if not (line_start in delete_lines_list):
               olines.append(ilines[line_start])
            line_start += 1
         funckind = find_functional_kind(call_lines[var])
         funcname = make_functional_name(call_lines[var])
         arglist  = make_input_args_list(call_lines[var])
         #DEBUG
         #print "funckind: ",funckind
         #print "funcname: ",funcname
         #print "arglist : ",arglist
         #DEBUG
         olines = append_subroutine_call(olines,funcname,arglist,orderdiff,ifbranch_kind,funckind,num,indent)
         #DEBUG
         #print "var: ",var,call_insert
         #DEBUG
      #break
      while line_start <= line_if_end:
         if not (line_start in delete_lines_list):
            olines.append(ilines[line_start])
         line_start += 1
   #break
   while line_start <= subr_lines_end:
      olines.append(ilines[line_start])
      line_start += 1
   subr_lines = find_subroutine(ilines,line_start)
   (subr_lines_start,subr_lines_end)=subr_lines
   #DEBUG
   #print "subrt: start,end:",subr_lines_start,subr_lines_end
   #DEBUG
while line_start < nlines:
    olines.append(ilines[line_start])
    line_start += 1
for line in olines:
   rewrap_line(line)
   #sys.stdout.write("%s\n"%line)
