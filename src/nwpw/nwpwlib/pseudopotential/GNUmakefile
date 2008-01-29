#
# $Id: GNUmakefile,v 1.14 2008-01-29 19:55:56 d3p708 Exp $
#
      SUBDIRS = paw_atom

      OBJ = \
                psp_generator_input.o \
                pspsolve.o \
                hgh_parse.o \
                cpi_parse.o \
                teter_parse.o \
                qmmm_parse.o \
                get_word.o name.o\
                loggrid.o pred_cor.o spline.o grids.o atom.o  \
                schrodin.o \
                pauli.o \
                dirac.o \
                dft.o \
                hartree.o \
                dirac_exchange.o \
                perdew_wang.o perdew_zunger.o vosko.o \
                pbe_exchange.o pbe_correlation.o \
                revpbe_exchange.o revpbe_correlation.o \
                becke_exchange.o lyp_correlation2.o \
                hamann.o \
                troullier.o xpansion.o gaussj.o \
                vanderbilt.o xpansion2.o \
                psp1d.o generate_rho_semicore.o \
		debug.o rhamann.o rpsp.o        

      LIBRARY = libnwpwlib.a

    USES_BLAS = 

  LIB_DEFINES = -g


include ../../../config/makefile.h
include ../../../config/makelib.h

