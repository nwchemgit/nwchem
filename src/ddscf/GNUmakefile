# $Id: GNUmakefile,v 1.1 1994-04-22 02:15:11 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_extract.o           rhf_energy.o \
        comp4_insert.o            rhf_fock_2e_a.o \
                                  rhf_fock_2e_b.o \
        dabsmax.o                 rhf_fock_2e_c.o \
        rhf_do_canon.o            rhf_get_param.o \
        rhf_fock_2e.o             rhf_nr_solve.o \
        ga_conform.o              rhf_precond.o \
        ga_create_atom_blocked.o  scf.o \
        ga_maxelt.o               schwarz_atom.o \
        get_atom_block.o          schwarz_init.o \
        rhf_hessv.o               schwarz_max.o \
        int_1e_ga.o               schwarz_print.o \
        k2cf.o                    schwarz_shell.o \
                                  schwarz_tidy.o \
        rhf.o                     upd_atom_block.o \
        nxtask.o

    LIBRARY = libddscf.a


include ../config/makefile.h
include ../config/makelib.h


rhf.o:		schwarz.fh
schwarz.o:	schwarz.fh cschwarz.fh
context.o:      context.fh















