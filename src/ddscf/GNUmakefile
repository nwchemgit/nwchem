# $Id: GNUmakefile,v 1.9 1994-08-15 17:28:25 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_extract.o           rhf_energy.o \
        comp4_insert.o            rhf_fock_2e_a.o \
                                  rhf_fock_2e_b.o \
        dabsmax.o                 rhf_fock_2e_c.o \
        rhf_do_canon.o            rhf_get_param.o \
        rhf_fock_2e.o             \
        ga_create_atom_blocked.o  scf.o \
                                  schwarz_atom.o \
        get_atom_block.o          schwarz_init.o \
        rhf_hessv.o               schwarz_max.o \
        int_1e_ga.o               schwarz_print.o \
        k2cf.o                    schwarz_shell.o \
                                  schwarz_tidy.o \
        rhf.o                     upd_atom_block.o \
        nxtask.o                  bad.o rhf_nr_solve.o \
        vectors.o \
	ri_invert_v.o             ga_create_cont_blocked.o \
        ri_v_fit_cd.o             ri_v_fit_j.o \
        ri_drive.o                dabssum.o \
        matrix_exp.o
      

    LIBRARY = libddscf.a

  USES_BLAS = *.F

HEADERS = schwarz.fh

include ../config/makefile.h
include ../config/makelib.h


rhf.o:		schwarz.fh
schwarz.o:	schwarz.fh cschwarz.fh
rhf_hessv.o rhf_nr_solve.o:	crhf.fh
