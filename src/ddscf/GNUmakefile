# $Id: GNUmakefile,v 1.74 1996-07-08 23:39:16 rg240 Exp $

#   LIB_DEFINES = -DDDSCF_TEXAS


  OBJ = comp4_bins.o              rhf_dens_mo.o \
        comp4_insert.o            scf_sym_adap.o \
        rhf_fock_2e.o             scf_lindep.o \
        ga_atom_blk.o  scf.o \
        schwarz_init.o \
        int_1e_ga.o               schwarz_pr.o \
        ga_orth_mos.o           schwarz_tidy.o \
        nxtask.o  \
        vectors.o \
	ga_cont_blk.o           mull_pop.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_proj.o \
        ao_fock_2e.o              movecs_lock.o \
	fock_2e.o                 fock_2e_et.o \
	fock_2e_a.o               \
	fock_2e_b.o               \
	fock_upd_blk.o       fock_2e_cdf.o \
	fock_get_blk.o      scf_vec_guess.o \
        rohf.o                   rohf_nr_sol.o     \
        rohf_wrap.o             rohf_canon.o \
        rohf_k2cf.o             rohf_precon.o \
        rohf_shift.o       rohf_step.o \
        rohf_rotate.o    scf_get_info.o \
        rohf_hscale.o      fock_2e_file.o \
        rohf_hxx.o         rohf_fock.o \
        rohf_fock2g.o        scf_movecs.o \
        rohf_h1e.o         rohf_h2e.o \
        rohf_dens.o             movecs_phase.o \
        riscf_fock.o         riscf_init.o \
        riscf_list.o            riscf_tidy.o \
        riscf_trans.o 		do_riscf.o \
	riscf_dovecs.o \
        uhf.o                   uhf_hessv.o fock_2e_rep.o \
        print_integ.o \
	riscf_square.o 		riscf_io.o \
        rohf_lagr.o 		rohf_diis.o \
	rhf_fock.o


 OBJ_OPTIMIZE = \
                fock_2e_task.o \
                fock_2e_slab.o \
                fock_2e_lab.o fock_labscal.o \
		fock_2e_c.o \
		fock_scale.o \
                fock_densc.o \
		schwarz_atom.o \
		comp4_ext.o \
		schwarz_sh.o \
                int2e_file.o \
                ao_replicated.o

      LIBRARY =	libddscf.a

    USES_BLAS =	movecs_phase.F rohf_k2cf.F fock_2e.F \
                uhf.F scf_movecs.F scf_vec_guess.F riscf_trans.F \
		int2e_file.F fock_2e_slab.F fock_2e_a.F \
		rhf_dens_mo.F scf.F schwarz_init.F vectors.F \
		movecs_proj.F fock_2e_b.F rohf.F rohf_nr_sol.F \
		rohf_precon.F \
		fock_2e_file.F rohf_fock.F rohf_h2e.F \
		scf_get_info.F ao_fock_2e.F rohf_diis.F movecs_lock.F \
                ao_replicated.F

HEADERS = schwarz.fh cscfps.fh cfock.fh

include ../config/makefile.h
include ../config/makelib.h

#activate disk array use of RI-SCF
#  DEFINES += -DRISCF_USE_DA

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
