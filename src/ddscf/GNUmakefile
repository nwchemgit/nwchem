# $Id: GNUmakefile,v 1.66 1996-02-27 00:29:14 d3e129 Exp $

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
	ga_cont_blk.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_proj.o \
        ao_fock_2e.o \
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
        uhf.o                   fock_2e_rep.o \
        print_integ.o \
	riscf_square.o 		riscf_io.o \
        rohf_lagr.o 		rohf_diis.o 


 OBJ_OPTIMIZE =	\
                fock_2e_task.o \
                fock_2e_slab.o \
                fock_2e_lab.o fock_labscal.o \
		fock_2e_c.o \
		fock_scale.o \
                fock_densc.o \
		schwarz_atom.o \
		comp4_ext.o \
		schwarz_sh.o \
                int2e_file.o 

      LIBRARY =	libddscf.a

    USES_BLAS =	movecs_phase.F rohf_k2cf.F fock_2e.F \
                uhf.F scf_movecs.F scf_vec_guess.F riscf_trans.F \
		int2e_file.F fock_2e_slab.F fock_2e_a.F \
		rhf_dens_mo.F scf.F schwarz_init.F vectors.F \
		movecs_proj.F fock_2e_b.F rohf.F rohf_nr_sol.F \
		rohf_precon.F \
		fock_2e_file.F rohf_fock.F rohf_h2e.F \
		scf_get_info.F ao_fock_2e.F rohf_diis.F

HEADERS = schwarz.fh cscfps.fh cfock.fh

include ../config/makefile.h
include ../config/makelib.h

#activate disk array use of RI-SCF
#  DEFINES += -DRISCF_USE_DA

# DO NOT EDIT BENEATH THIS LINE ... GENERATED AUTOMATICALLY
$(LIBRARY_PATH)(ao_fock_2e.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(comp4_bins.o):	 $(INCDIR)/bitops_funcs.fh $(INCDIR)/bitops_decls.fh 
$(LIBRARY_PATH)(comp4_ext.o):	 $(INCDIR)/bitops_funcs.fh $(INCDIR)/bitops_decls.fh 
$(LIBRARY_PATH)(comp4_insert.o):	 $(INCDIR)/bitops.fh 
$(LIBRARY_PATH)(do_riscf.o):	 ./cfock.fh 
$(LIBRARY_PATH)(fock_2e.o):	 ./cscf.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_a.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_b.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/mafdecls.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_2e_cdf.o):	 ./cscfps.fh 
$(LIBRARY_PATH)(fock_2e_et.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_2e_file.o):	 $(INCDIR)/util.fh ./cscfps.fh 
$(LIBRARY_PATH)(fock_2e_lab.o):	 ./cfockmul.fh 
$(LIBRARY_PATH)(fock_2e_rep.o):	 $(INCDIR)/global.fh ./cfock.fh 
$(LIBRARY_PATH)(fock_2e_task.o):	 $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/geom.fh ./cscfps.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_get_blk.o):	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(fock_upd_blk.o):	 $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_atom_blk.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_cont_blk.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(ga_orth_mos.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int2e_file.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cint2efile.fh ./cfock.fh $(INCDIR)/bitops.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(int_1e_ga.o):	 $(INCDIR)/sym.fh ./schwarz.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh ./cint1cache.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(matrix_exp.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh 
$(LIBRARY_PATH)(movecs_phase.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(movecs_proj.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(nxtask.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rhf_dens_mo.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_dovecs.o):	 ./cfock.fh 
$(LIBRARY_PATH)(riscf_fock.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh ./cscfps.fh ./cscf.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_init.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/stdio.fh ./schwarz.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_io.o):	 $(INCDIR)/msgids.fh $(INCDIR)/global.fh $(INCDIR)/da.fh 
$(LIBRARY_PATH)(riscf_list.o):	 ./riscf.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh 
$(LIBRARY_PATH)(riscf_square.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_tidy.o):	 $(INCDIR)/geom.fh $(INCDIR)/da.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(riscf_trans.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/da.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf_canon.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_dens.o):	 $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_diis.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_fock.o):	 $(INCDIR)/util.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(rohf_fock2g.o):	 $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_h1e.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_h2e.o):	 $(INCDIR)/util.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh 
$(LIBRARY_PATH)(rohf_hscale.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_hxx.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_k2cf.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(rohf_lagr.o):	 $(INCDIR)/tcgmsg.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_nr_sol.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_precon.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_rotate.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh 
$(LIBRARY_PATH)(rohf_shift.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_step.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscfps.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(rohf_wrap.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./crohf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscfps.fh ./cscf.fh $(INCDIR)/context.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_get_info.o):	 $(INCDIR)/util.fh $(INCDIR)/sym.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cscf.fh ./cfock.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_lindep.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_movecs.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/global.fh ./cscf.fh ./crohf.fh 
$(LIBRARY_PATH)(scf_pstat.o):	 $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/mafdecls.fh ./cscfps.fh ./cscf.fh 
$(LIBRARY_PATH)(scf_sym_adap.o):	 $(INCDIR)/sym.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(scf_vec_guess.o):	 $(INCDIR)/util.fh $(INCDIR)/rtdb.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(schwarz_atom.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_init.o):	 $(INCDIR)/util.fh ./schwarz.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cschwarz.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(schwarz_max.o):	 ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_pr.o):	 $(INCDIR)/tcgmsg.fh ./schwarz.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/global.fh 
$(LIBRARY_PATH)(schwarz_sh.o):	 $(INCDIR)/mafdecls.fh $(INCDIR)/itri.fh ./cschwarz.fh 
$(LIBRARY_PATH)(schwarz_tidy.o):	 $(INCDIR)/mafdecls.fh ./cschwarz.fh 
$(LIBRARY_PATH)(uhf.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/pstat.fh $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh $(INCDIR)/geom.fh ./cuhf.fh ./cuhf.fh ./cscfps.fh ./cscf.fh $(INCDIR)/bas.fh 
$(LIBRARY_PATH)(vectors.o):	 $(INCDIR)/util.fh $(INCDIR)/tcgmsg.fh $(INCDIR)/rtdb.fh $(INCDIR)/msgtypesf.h $(INCDIR)/msgids.fh $(INCDIR)/mafdecls.fh $(INCDIR)/inp.fh $(INCDIR)/global.fh ./cscfps.fh 
