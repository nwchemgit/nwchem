# $Id: GNUmakefile,v 1.29 1995-04-06 21:59:35 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_insert.o            \
        rhf_fock_2e.o \
        ga_create_atom_blocked.o  scf.o \
        schwarz_init.o \
        schwarz_max.o \
        int_1e_ga.o               schwarz_print.o \
        ga_orthog_mos.o           schwarz_tidy.o \
        nxtask.o  \
        vectors.o \
	ga_create_cont_blocked.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_project.o \
	fock_2e.o \
	fock_2e_a.o \
	fock_2e_b.o \
	fock_upd_atom_blk.o  \
	fock_get_atom_blk.o      scf_vectors_guess.o \
        rohf.o                   rohf_nr_solve.o     \
        rohf_wrap.o             rohf_canon.o \
        rohf_k2cf.o             rohf_search_precond.o \
        rohf_make_shift.o       rohf_step_energy.o \
        rohf_rotate_movecs.o    scf_get_info.o \
        rohf_hdiag_scale.o \
        rohf_hessv_xx.o         rohf_fock.o \
        rohf_fock2grad.o        scf_movecs.o \
        rohf_hessv_1e.o         rohf_hessv_2e.o \
        rohf_dens.o             movecs_fix_phase.o \
        dabsmax.o               uhf.o




# rhf_fock_2e_b compiles incorrectly with optimize on the SGITFP

 OBJ_OPTIMIZE =	\
		fock_2e_c.o \
		fock_eri_scale.o \
                fock_density_screen.o \
		schwarz_atom.o \
		comp4_extract.o \
		schwarz_shell.o

      LIBRARY =	libddscf.a

    USES_BLAS =	movecs_fix_phase.F rohf_k2cf.F fock_2e.F k2cf.F \
                uhf.F scf_movecs.F scf_vectors_guess.F


#
 LIB_DEFINES = 
#


HEADERS = schwarz.fh cscfps.fh

include ../config/makefile.h
include ../config/makelib.h

fock_2e_a.o rhf_fock_2e_a.o rhf_get_param.o rhf_hessv.o scf_get_info.o:	 cfock.fh
rhf_hessv.o rhf_nr_solve.o:	 crhf.fh
rohf.o rohf_make_shift.o rohf_nr_solve.o rohf_search_precond.o rohf_step_energy.o rohf_wrap.o scf_movecs.o:	 crohf.fh
rohf_make_shift.o rohf_nr_solve.o rohf_rotate_movecs.o rohf_search_precond.o rohf_step_energy.o rohf_wrap.o scf_get_info.o scf_movecs.o scf_vectors_guess.o:	 cscf.fh
fock_2e.o fock_2e_b.o ga_orthog_mos.o int_1e_ga.o matrix_exp.o rhf_canon.o rhf_canon_save.o rhf_dens_to_mo.o rhf_do_canon.o rhf_fock_2e.o rhf_fock_2e_b.o rhf_nr_solve.o rohf_nr_solve.o rohf_search_precond.o rohf_step_energy.o scf.o scf_pstat.o vectors.o:	 cscfps.fh
schwarz_atom.o schwarz_init.o schwarz_max.o schwarz_shell.o schwarz_tidy.o:	 cschwarz.fh
rhf_dens_to_mo.o rhf_energy.o rhf_hessv.o:	 cskeleton.fh
fock_2e.o fock_2e_a.o fock_2e_b.o rhf_fock_2e.o rhf_fock_2e_a.o rhf_fock_2e_b.o schwarz_atom.o schwarz_init.o schwarz_max.o schwarz_print.o schwarz_shell.o schwarz_tidy.o:	 schwarz.fh
