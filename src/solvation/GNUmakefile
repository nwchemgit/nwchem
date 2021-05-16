#
# $Id$
#
   OBJ = mnsol_mem.o\
         mnsol_interface.o\
         mnsol.o\
         solv_data.o\
         cosmo_smd_output.o\
         cosmo_vem_output.o hnd_cosmo_lib.o

   OBJ_OPTIMIZE = cosmo_input.o \
         cosmo_initialize.o\
         cosmo_tidy.o\
         cosmo_def_radii.o\
         cosmo_charges.o\
         grad_hnd_cos.o

   LIBRARY = libsolvation.a

# If you add -DJEFF_OUTPUT to LIB_DEFINES, the dipole polarizability will be printed in the same format as the TCE.
# This is very useful (1) if you wrote scripts to parse that particular formatting and (2) if you want the anisotrophy.
#   THE LAPACK CALL ENABLED BY THIS IS BROKEN ON SOME ARCHITECTURES
#   LIB_DEFINES = -DJEFF_OUTPUT

   LIB_INCLUDES = -I../nwdft/include -I../ddscf

USES_BLAS = hnd_cosmo_lib.F

HEADERS = ../property/prop.fh ../property/atomdata.fh cosmoP.fh cosmo_params.fh

include ../config/makefile.h
include ../config/makelib.h

#FDEBUG += -O0
