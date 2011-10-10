#FC=ifc -w -132
.SUFFIXES:

.SUFFIXES: .o .f90

.f90.o :
	$(FC) $(F90FLAGS) $(LIB_INCLUDES)  -I$(TOPDIR)/src/include -c $<

# $Id$
