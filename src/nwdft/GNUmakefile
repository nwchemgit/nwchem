# $Id: GNUmakefile,v 1.4 1994-10-26 08:48:32 og845 Exp $

   OBJ = nwdft.o

   LIBRARY = libnwdft.a
#   SUBDIRS = input_dft symmetry scf_dft coulomb xc util
   SUBDIRS = input_dft scf_dft coulomb xc util

   HEADERS = \
      include/atoms.fh  \
      include/handles.fh \
      include/timer.fh   \
      include/io.fh    \
      include/types.fh \
      include/quad.fh  \
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
      include/io.fh    \
      include/types.fh \
      include/quad.fh  \
      include/dimensions.fh  \
      include/stdio.fh
