# $Id: GNUmakefile,v 1.43 1995-11-02 01:52:14 rg240 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_insert.o            scf_movecs_sym_adapt.o \
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
        ao_fock_2e.o \
	fock_2e.o                 \
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
		rohf_search_precond.F scf_vectors_guess.F \
		fock_2e_from_file.F rohf_fock.F rohf_hessv_2e.F \
		scf_get_info.F ao_fock_2e.F

#
 LIB_DEFINES = 
#


HEADERS = schwarz.fh cscfps.fh cfock.fh

include ../config/makefile.h
include ../config/makelib.h

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
ao_fock_2e.o :	 schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscfps.fh cfock.fh $(INCDIR)/bas.fh 
bad.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
do_ri.o :	 cfock.fh 
fock_2e.o :	 cscf.fh cfock.fh 
fock_2e_a.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh schwarz.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cfock.fh 
fock_2e_a.semi.save.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh schwarz.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cfock.fh 
fock_2e_b.o :	 $(INCDIR)/util.fh schwarz.fh $(INCDIR)/mafdecls.fh cscfps.fh cfock.fh $(INCDIR)/bas.fh 
fock_2e_from_file.o :	 $(INCDIR)/util.fh cscfps.fh 
fock_2e_report.o :	 $(INCDIR)/global.fh cfock.fh 
fock_2e_task.o :	 cfock.fh $(INCDIR)/bas.fh 
fock_get_atom_blk.o :	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
fock_upd_atom_blk.o :	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
ga_create_atom_blocked.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
ga_create_cont_blocked.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
ga_mxv.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
ga_orthog_mos.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh $(INCDIR)/bas.fh 
int2e_file.o :	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh cint2efile.fh 
int_1e_ga.o :	 $(INCDIR)/sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh cscfps.fh $(INCDIR)/bas.fh 
k2cf.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
matrix_exp.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh 
movecs_fix_phase.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
movecs_project.o :	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
nxtask.o :	 $(INCDIR)/global.fh 
rhf_dens_to_mo.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscfps.fh $(INCDIR)/bas.fh 
riscf_fock_2e.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh $(INCDIR)/bas.fh 
riscf_init.o :	 $(INCDIR)/tcgmsg.fh schwarz.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh $(INCDIR)/bas.fh 
riscf_list.o :	 riscf.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
riscf_square.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/bas.fh 
riscf_trans_int.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
rohf.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh cscf.fh crohf.fh 
rohf_canon.o :	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh cscf.fh crohf.fh 
rohf_dens.o :	 $(INCDIR)/global.fh 
rohf_fock.o :	 $(INCDIR)/util.fh $(INCDIR)/global.fh cscf.fh crohf.fh $(INCDIR)/bas.fh 
rohf_fock2grad.o :	 $(INCDIR)/global.fh cscf.fh crohf.fh 
rohf_hdiag_scale.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rohf_hessv_1e.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rohf_hessv_2e.o :	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rohf_hessv_xx.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rohf_k2cf.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
rohf_make_shift.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscf.fh crohf.fh 
rohf_nr_solve.o :	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscfps.fh cscf.fh crohf.fh 
rohf_rotate_movecs.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscf.fh 
rohf_search_precond.o :	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh cscf.fh crohf.fh 
rohf_step_energy.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscfps.fh cscf.fh crohf.fh 
rohf_wrap.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscf.fh crohf.fh $(INCDIR)/bas.fh 
scf.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscfps.fh cscf.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
scf_get_info.o :	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cscf.fh cfock.fh $(INCDIR)/bas.fh 
scf_movecs.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh cscf.fh crohf.fh 
scf_movecs_sym_adapt.o :	 $(INCDIR)/sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
scf_pstat.o :	 $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh cscfps.fh cscf.fh 
scf_vectors_guess.o :	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh cscf.fh $(INCDIR)/bas.fh 
schwarz_atom.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh cschwarz.fh 
schwarz_init.o :	 $(INCDIR)/util.fh schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cschwarz.fh $(INCDIR)/bas.fh 
schwarz_max.o :	 cschwarz.fh 
schwarz_print.o :	 $(INCDIR)/tcgmsg.fh schwarz.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
schwarz_shell.o :	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh cschwarz.fh 
schwarz_tidy.o :	 $(INCDIR)/mafdecls.fh cschwarz.fh 
uhf.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh cuhf.fh cuhf.fh cscfps.fh cscf.fh 
vectors.o :	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh cscfps.fh 
