# $Id: GNUmakefile,v 1.51 1995-12-01 16:58:10 d3g681 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_insert.o            scf_movecs_sym_adapt.o \
        rhf_fock_2e.o \
        ga_create_atom_blocked.o  scf.o \
        schwarz_init.o \
        int_1e_ga.o               schwarz_print.o \
        ga_orthog_mos.o           schwarz_tidy.o \
        nxtask.o  \
        vectors.o \
	ga_create_cont_blocked.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_project.o \
        ao_fock_2e.o \
	fock_2e.o                 fock_2e_est_task.o \
	fock_2e_a.o               \
	fock_2e_b.o               \
	fock_upd_atom_blk.o       fock_2e_cache_dens_fock.o \
	fock_get_atom_blk.o      scf_vectors_guess.o \
        rohf.o                   rohf_nr_solve.o     \
        rohf_wrap.o             rohf_canon.o \
        rohf_k2cf.o             rohf_search_precond.o \
        rohf_make_shift.o       rohf_step_energy.o \
        rohf_rotate_movecs.o    scf_get_info.o \
        rohf_hdiag_scale.o      fock_2e_from_file.o \
        rohf_hessv_xx.o         rohf_fock.o \
        rohf_fock2grad.o        scf_movecs.o \
        rohf_hessv_1e.o         rohf_hessv_2e.o \
        rohf_dens.o             movecs_fix_phase.o \
        riscf_fock_2e.o         riscf_init.o \
        riscf_list.o            riscf_tidy.o \
        riscf_trans_int.o \
        uhf.o                   fock_2e_report.o fock_2e_task.o \
        print_integ_list.o \
        do_ri.o rohf_lagr.o 


# rhf_fock_2e_b compiles incorrectly with optimize on the SGITFP

 OBJ_OPTIMIZE =	\
                fock_2e_shell_to_labels.o \
                fock_2e_label.o fock_eri_label_scale.o \
		fock_2e_c.o \
		fock_eri_scale.o \
                fock_density_screen.o \
		schwarz_atom.o \
		comp4_extract.o \
		schwarz_shell.o \
                riscf_square.o \
                int2e_file.o 

      LIBRARY =	libddscf.a

    USES_BLAS =	movecs_fix_phase.F rohf_k2cf.F fock_2e.F \
                uhf.F scf_movecs.F scf_vectors_guess.F riscf_trans_int.F \
		int2e_file.F fock_2e_shell_to_labels.F fock_2e_a.F \
		rhf_dens_to_mo.F scf.F schwarz_init.F vectors.F \
		movecs_project.F fock_2e_b.F rohf.F rohf_nr_solve.F \
		rohf_search_precond.F \
		fock_2e_from_file.F rohf_fock.F rohf_hessv_2e.F \
		scf_get_info.F ao_fock_2e.F

#
 LIB_DEFINES = 
#


HEADERS = schwarz.fh cscfps.fh cfock.fh

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(ao_fock_2e.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(comp4_bins.o):	 $(INCDIR)/bitops_funcs.fh $(INCDIR)/bitops_decls.fh 
$(LIBRARY_PATH)(comp4_extract.o):	 $(INCDIR)/bitops_funcs.fh $(INCDIR)/bitops_decls.fh 
$(LIBRARY_PATH)(comp4_insert.o):	 $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(do_ri.o):	 ./cfock.fh 
$(LIBRARY_PATH)(fock_2e.o):	 ./cscf.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_a.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_b.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/mafdecls.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_2e_cache_dens_fock.o):	 ./cscfps.fh 
$(LIBRARY_PATH)(fock_2e_est_task.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_2e_from_file.o):	 $(INCDIR)/util.fh ./cscfps.fh 
$(LIBRARY_PATH)(fock_2e_report.o):	 $(INCDIR)/global.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_task.o):	 $(INCDIR)/sym.fh ./schwarz.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_get_atom_blk.o):	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_upd_atom_blk.o):	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_create_atom_blocked.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_create_cont_blocked.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_orthog_mos.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int2e_file.o):	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cint2efile.fh $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(int_1e_ga.o):	 $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(matrix_exp.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh 
$(LIBRARY_PATH)(movecs_fix_phase.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(movecs_project.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(nxtask.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rhf_dens_to_mo.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_fock_2e.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_init.o):	 $(INCDIR)/tcgmsg.fh ./schwarz.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_list.o):	 ./riscf.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(riscf_square.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_trans_int.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf_canon.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_dens.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_fock.o):	 $(INCDIR)/util.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf_fock2grad.o):	 $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_hdiag_scale.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_hessv_1e.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_hessv_2e.o):	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh 
$(LIBRARY_PATH)(rohf_hessv_xx.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_k2cf.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_lagr.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_make_shift.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_nr_solve.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_rotate_movecs.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh 
$(LIBRARY_PATH)(rohf_search_precond.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_step_energy.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_wrap.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cscf.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_get_info.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_movecs.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(scf_movecs_sym_adapt.o):	 $(INCDIR)/sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_pstat.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh ./cscfps.fh ./cscf.fh 
$(LIBRARY_PATH)(scf_vectors_guess.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(schwarz_atom.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_init.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cschwarz.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(schwarz_max.o):	 ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_print.o):	 $(INCDIR)/tcgmsg.fh ./schwarz.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(schwarz_shell.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_tidy.o):	 $(INCDIR)/mafdecls.fh ./cschwarz.fh 
$(LIBRARY_PATH)(uhf.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cuhf.fh ./cuhf.fh ./cscfps.fh ./cscf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(vectors.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh 
