# $Id: GNUmakefile,v 1.25 1995-03-14 02:00:02 rg240 Exp $


  OBJ = comp4_bins.o              rhf_dens_to_mo.o \
        comp4_insert.o            rhf_energy.o \
        rhf_canon.o \
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
	ga_create_cont_blocked.o \
        matrix_exp.o              scf_pstat.o \
        schwarz_max.o             movecs_project.o \
	fock_2e.o \
	fock_2e_a.o \
	fock_2e_b.o \
	fock_upd_atom_blk.o  \
	fock_get_atom_blk.o 


# rhf_fock_2e_b compiles incorrectly with optimize on the SGITFP

 OBJ_OPTIMIZE =	\
		fock_2e_c.o \
		fock_eri_scale.o \
		rhf_fock_2e_a.o \
		rhf_fock_2e_c.o \
		dabsmax.o \
		dabssum.o  \
		schwarz_atom.o \
		comp4_extract.o \
		schwarz_shell.o

      LIBRARY =	libddscf.a

    USES_BLAS =	\
		fock_2e.F \
		ga_conform.F \
		ga_mxv.F \
		k2cf.F \
		matrix_exp.F \
		rhf_dens_to_mo.F \
		rhf_do_canon.F \
		rhf_energy.F \
		rhf_fock_2e.F \
		rhf_hessv.F \
		rhf_nr_solve.F \
		int_1e_ga.F

#
# Possible defines
#        GENERIC_FOCK        : Use wrapper rhf_fock_2e to call generic fock routines
#
#
 LIB_DEFINES = 

#
#


HEADERS = schwarz.fh cscfps.fh

include ../config/makefile.h
include ../config/makelib.h

rhf.o:		schwarz.fh
schwarz.o:	schwarz.fh cschwarz.fh
rhf_hessv.o rhf_nr_solve.o:	crhf.fh

int_1e_ga.o matrix_exp.o rhf_fock_2e.o rhf_fock_2e_b.o \
rhf_nr_solve.o scf_pstat.o rhf_dens_to_mo.o \
rhf_do_canon.o vectors.o:      cscfps.fh

rhf_fock_2e_a.o rhf_get_param.o:        cfock.fh


fock_2e.o: cscfps.fh

fock_2e_a.o: cfock.fh

fock_2e_b.o: cscfps.fh


