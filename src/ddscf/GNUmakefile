# $Id$

#   LIB_DEFINES = -DDDSCF_TEXAS


  OBJ = comp4_bins.o              rhf_dens_mo.o \
        comp4_insert.o            scf_sym_adap.o \
        rhf_fock_2e.o             scf_lindep.o \
        ga_atom_blk.o  scf.o      scf_input.o \
        schwarz_init.o            vectors_inp.o \
        int_1e_ga.o               schwarz_pr.o \
        ga_orth_mos.o           schwarz_tidy.o \
        int2e_balance.o \
        vectors.o               scf_sym_scr.o \
	ga_cont_blk.o           mull_pop.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_proj.o \
        movecs_rotate.o movecs_molden.o \
        ao_fock_2e.o              movecs_lock.o \
	fock_2e.o                 fock_2e_et.o \
	fock_2e_a.o               print_aoints.o \
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
        rohf_hessv2.o  \
        rohf_dens.o             movecs_phase.o \
        riscf_fock.o         riscf_init.o \
        riscf_list.o            riscf_tidy.o \
        riscf_trans.o 		do_riscf.o \
	riscf_dovecs.o          rohf_scat.o \
        uhf.o                   uhf_hessv.o fock_2e_rep.o \
        uhf_hessv2.o \
        print_integ_list.o           print_1eints.o \
	riscf_square.o 		riscf_io.o movecs_pr_anal.o \
        rohf_lagr.o 		rohf_diis.o \
	rhf_fock.o              multipole.o ga_mat_to_vec.o \
        movecs_frag.o  localize.o int_dip_ga.o \
        fock_xc.o fock_j_fit.o rohf_jkfac.o uhf_jkfac.o \
        rohf_hessv3.o \
	fock_2e_cam.o ga_get2eri.o scf_frozemb.o print_integrals.o \
	rohf_hessv2_ext.o rohf_hessv3_ext.o uhf_hessv2_ext.o uhf_precond_ext.o \
	scf_dump_evp_input.o


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

    USES_BLAS = ao_fock_2e.F ao_replicated.F fock_2e_file.F fock_2e_slab.F int2e_file.F \
                movecs_lock.F movecs_phase.F riscf_fock.F riscf_trans.F rohf_diis.F \
                rohf_k2cf.F scf_movecs.F scf_vec_guess.F uhf.F movecs_frag.F localize.F fock_2e.F \
		fock_2e_cam.F int_1e_ga.F movecs_rotate.F vectors.F ga_get2eri.F

HEADERS = schwarz.fh cscfps.fh cfock.fh cscf.fh cint2efile.fh cfockmul.fh crohf.fh cosmo.fh cuhf.fh frozemb.fh

#activate disk array use of RI-SCF
  DEFINES += -DRISCF_USE_DRA


SUBDIRS = fast

include ../config/makefile.h
include ../config/makelib.h
