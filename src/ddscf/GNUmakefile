# $Id: GNUmakefile,v 1.5 1994-05-27 02:01:23 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_extract.o           rhf_energy.o \
        comp4_insert.o            rhf_fock_2e_a.o \
                                  rhf_fock_2e_b.o \
        dabsmax.o                 rhf_fock_2e_c.o \
        rhf_do_canon.o            rhf_get_param.o \
        rhf_fock_2e.o             \
        ga_conform.o              rhf_precond.o \
        ga_create_atom_blocked.o  scf.o \
        ga_maxelt.o               schwarz_atom.o \
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
        ri_drive.o                dabssum.o
      

    LIBRARY = libddscf.a

  USES_BLAS = rhf_fock_2e.F      rhf_fock_2e_a.F \
              rhf_fock_2e_b.F    rhf_fock_2e_c.F \
              rhf_get_param.F    ga_conform.F \
              rhf_hessv.F        rhf_nr_solve.F \
              ga_maxelt.F        rhf_precond.F \
              scf.F              schwarz_atom.F \
              schwarz_init.F     k2cf.F \
              schwarz_max.F \
              schwarz_print.F    rhf.F \
              rhf_dens_to_mo.F   schwarz_tidy.F \
              rhf_do_canon.F     rhf_energy.F



include ../config/makefile.h
include ../config/makelib.h


rhf.o:		schwarz.fh
schwarz.o:	schwarz.fh cschwarz.fh
context.o:      context.fh















