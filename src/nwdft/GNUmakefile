#$Id: GNUmakefile,v 1.11 1995-04-26 18:09:25 og845 Exp $
# makefile,v 1.8 1994/12/05 20:37:08 og845 Exp

   OBJ = nwdft.o

   LIBRARY = libnwdft.a
   LIB_DEFINES = -DDFT_TIMER
   SUBDIRS = input_dft scf_dft coulomb xc util 3d

   HEADERS = \
      include/timer.fh   \
      include/quad.fh  \
      include/dimensions.fh  


include ../config/makefile.h
include ../config/makelib.h

#
# add any include dependencies here...
# e.g. ...
# rhf.o:		schwarz.fh
#

nwdft.o:  \
      include/timer.fh   \
      include/quad.fh  \
      include/dimensions.fh  
