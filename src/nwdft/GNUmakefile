# $Id: GNUmakefile,v 1.3 1994-10-25 23:44:18 og845 Exp $

#   OBJ = nwdft.o

   LIBRARY = libnwdft.a
#   SUBDIRS = input_dft symmetry scf_dft coulomb xc util
   SUBDIRS = input_dft scf_dft coulomb xc util

   HEADERS = \
      include/atoms.fh  \
      include/handles.fh \
      include/timer.fh   \
      include/CD_basis_set.fh \
      include/bas_dft.fh  \
      include/io.fh    \
      include/types.fh \
      include/quad.fh  \
      include/XC_basis_set.fh \
      include/dimensions.fh  \
      include/stdio.fh    \
      include/bas_new.fh


include ../config/makefile.h
include ../config/makelib.h

#
# add any include dependencies here...
# e.g. ...
# rhf.o:		schwarz.fh
#

nwdft.o:  \
      include/atoms.fh  \
      include/handles.fh \
      include/timer.fh   \
      include/CD_basis_set.fh \
      include/bas_dft.fh  \
      include/io.fh    \
      include/types.fh \
      include/quad.fh  \
      include/XC_basis_set.fh \
      include/dimensions.fh  \
      include/stdio.fh
