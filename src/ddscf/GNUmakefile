# $Id: GNUmakefile,v 1.19 1994-12-23 03:28:16 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_insert.o            rhf_energy.o \
        rhf_canon.o               rhf_fock_2e_a.o \
        rhf_do_canon.o            rhf_get_param.o \
        rhf_fock_2e.o             rhf_fock_2e_b.o \
        ga_create_atom_blocked.o  scf.o \
        get_atom_block.o          schwarz_init.o \
        rhf_hessv.o               schwarz_max.o \
        int_1e_ga.o               schwarz_print.o \
        k2cf.o                    \
        ga_orthog_mos.o           schwarz_tidy.o \
        rhf.o                     upd_atom_block.o \
        nxtask.o                  bad.o rhf_nr_solve.o \
        vectors.o \
	ri_invert_v.o             ga_create_cont_blocked.o \
        ri_v_fit_cd.o             ri_v_fit_j.o \
        ri_drive.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_project.o


# rhf_fock_2e_b compiles incorrectly with optimize on the SGITFP

  OBJ_OPTIMIZE =  rhf_fock_2e_c.o dabsmax.o dabssum.o  \
                  schwarz_atom.o comp4_extract.o schwarz_shell.o

    LIBRARY = libddscf.a

  USES_BLAS = ga_conform.F ga_mxv.F k2cf.F matrix_exp.F rhf_dens_to_mo.F \
              rhf_do_canon.F \
              rhf_energy.F rhf_fock_2e.F rhf_hessv.F rhf_nr_solve.F \
              ri_invert_v.F \
              ri_v_fit_cd.F ri_v_fit_j.F

HEADERS = schwarz.fh

include ../config/makefile.h
include ../config/makelib.h

rhf.o:		schwarz.fh
schwarz.o:	schwarz.fh cschwarz.fh
rhf_hessv.o rhf_nr_solve.o:	crhf.fh

int_1e_ga.o matrix_exp.o rhf_fock_2e.o rhf_fock_2e_b.o \
rhf_nr_solve.o scf_pstat.o rhf_dens_to_mo.o \
rhf_do_canon.o vectors.o:      cscfps.fh

rhf_fock_2e_a.o rhf_get_param.o:        cfock.fh
