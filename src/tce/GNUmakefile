#$Id: GNUmakefile,v 1.54 2008-03-21 21:46:58 jhammond Exp $
#

SUBDIRS = include ccsd ccsd_lambda eomccsd \
          ccsd_t cr-eomccsd_t lrh emb diis2 \
          response  ccsd_yr ccsd_alpha ccsd_beta \
          ccd lccd lccsd qcisd cc3 \
          mbpt1 mbpt2 mbpt3 mbpt4 cisd cisdt cisdtq \
          ccsdt ccsdt_lambda eomccsdt ccsdt_alpha_offdiag \
          ccsdtq ccsdtq_lambda eomccsdtq ccsdtq_alpha_offdiag

OBJ_OPTIMIZE = tce_input.o tce_init.o tce_tile.o tce_ao1e.o tce_mo1e.o tce_ao2e.o \
      tce_filename.o createfile.o deletefile.o reconcilefile.o dratoga.o gatodra.o gatoeaf.o \
      get_block.o put_block.o add_block.o tce_tidy.o sf_test.o hbar.o \
      tce_mo2e.o tce_mo1e_offset.o tce_mo2e_offset.o \
      tce_e_offset.o tce_t1_offset.o tce_t2_offset.o tce_t3_offset.o tce_t4_offset.o \
      tce_guess_t1.o tce_guess_t2.o tce_zero.o tce_print_x1.o tce_print_x2.o \
      tce_residual_t1.o tce_residual_t2.o tce_residual_t3.o tce_residual_t4.o \
      tce_jacobi_t1.o tce_jacobi_t2.o tce_jacobi_t3.o tce_jacobi_t4.o \
      tce_jacobi_l1.o tce_jacobi_l2.o tce_jacobi_l3.o tce_jacobi_l4.o \
      ma_pack.o sf_print.o ma_print_compact.o ma_transpose.o ma_scale.o \
      tce_diis_init.o tce_diis_tidy.o tce_diis_switch.o tce_diis.o \
      copyfile.o ddotfile.o daxfile.o daxpyfile.o dscalfile.o \
      tce_aod1.o tce_nud1.o tce_mod1.o tce_mod1_offset.o sym_abelian_axis.o \
      tce_eom_init.o tce_eom_tidy.o tce_eom_xguess.o tce_eom_xdiagon.o tce_eom_yguess.o tce_eom_ydiagon.o \
      tce_x1_offset.o tce_x2_offset.o tce_x3_offset.o tce_x4_offset.o \
      tce_y1_offset.o tce_y2_offset.o tce_y3_offset.o tce_y4_offset.o \
      tce_jacobi_x1.o tce_jacobi_x2.o tce_jacobi_x3.o tce_jacobi_x4.o \
      tce_jacobi_y1.o tce_jacobi_y2.o tce_jacobi_y3.o tce_jacobi_y4.o \
      tce_mbpt2.o tce_energy.o tce_energy_fragment.o tce_sort.o tce_sortacc.o tce_restricted.o \
      tce_hash.o get_hash_block.o put_hash_block.o add_hash_block.o tce_sort_4kg.o \
      tce_residual_t3a.o tce_t3a_offset.o tce_jacobi_t3a.o \
      tce_x3a_offset.o tce_jacobi_x3a.o \
      e2_offset_2_offset.o get_block_ind.o tce_hash_ind.o \
      tce_mo2e_offset_intorb.o tce_mo2e_offset_size.o tce_mo2e_trans.o \
      atpart_nalength.o tce_mo2e_zones_4a_disk_ga.o tce_mo2e_zones_4a_disk_ga_chop.o \
      tce_mo2e_zones_4a_disk_2s_new.o\
      tce_dens_hh_offset.o tce_dens_hp_offset.o tce_dens_pp_offset.o tce_dens_ph_offset.o \
      get_mo_rdm_hh.o get_mo_rdm_hp.o get_mo_rdm_ph.o get_mo_rdm_pp.o \
      ao_rdm_write.o tce_prop_start.o tce_property.o tce_diagnose_t1.o \
      tce_ao_mp1.o tce_mo_mp1.o tce_mp_filenames.o tce_mp_irreps.o tce_mo_mp1_offsets.o \
      tce_fourindex_driver.o tce_fock_create.o tce_multipole_create.o tce_reference_check.o \
      tce_fourindex_cleanup.o tce_fock_destroy.o tce_multipole_destroy.o \
      tce_create_t1.o tce_create_t2.o tce_create_t3.o tce_create_t4.o \
      tce_create_x1.o tce_create_x2.o tce_create_x3.o tce_create_x4.o \
      tce_create_y1.o tce_create_y2.o tce_create_y3.o tce_create_y4.o \
      tce_clone_x1.o tce_clone_x2.o tce_clone_x3.o tce_clone_x4.o \
      tce_clone_y1.o tce_clone_y2.o tce_clone_y3.o tce_clone_y4.o \
      tce_t1_offset_new.o tce_t2_offset_new.o tce_t3_offset_new.o tce_t4_offset_new.o \
      tce_x1_offset_new.o tce_x2_offset_new.o tce_x3_offset_new.o tce_x4_offset_new.o \
      tce_y1_offset_new.o tce_y2_offset_new.o tce_y3_offset_new.o tce_y4_offset_new.o \
      tce_ccsd_driver.o tce_ccsdt_driver.o tce_ccsdtq_driver.o \
      tce_ccsd_perturbative.o \
      tce_ccsd_lambda_driver.o tce_ccsdt_lambda_driver.o tce_ccsdtq_lambda_driver.o \
      tce_ccsd_response_driver.o tce_ccsdt_response_driver.o tce_ccsdtq_response_driver.o \
      tce_ccsd_symm_polar.o tce_ccsdt_symm_polar.o tce_ccsdtq_symm_polar.o \
      tce_ccsd_lambda_response_driver.o tce_ccsd_asym_polar.o tce_ccsd_hyperpolar.o \
      tce_ao2e_new.o tce_tile_new.o tce_init_new.o #tce_mo2e_new.o
      
ifdef NEW_PARIO
      OBJ_OPTIMIZE += tce_mo2e_zones_4a_disk_4s_new_patch.o\
                      tce_mo2e_zones_4a_disk_4s_new_patch_m2.o 
endif

LIB_INCLUDES = -I./include

LIBRARY = libtce.a

USES_BLAS = tce_mbpt2.F tce_residual_t1.F tce_residual_t2.F \
            tce_residual_t3.F tce_residual_t4.F tce_residual_t3a.F \
            tce_ao2e.F tce_mo2e.F tce_guess_t1.F tce_diis.F tce_mod1.F \
            tce_eom_xdiagon.F dscalfile.F tce_eom_ydiagon.F\
            tce_ao1e.F tce_aod1.F tce_mo1e.F tce_ao_mp1.F tce_mo_mp1.F\
            tce_mo2e_zones_4a_disk_ga.F tce_mo2e_zones_4a_disk_ga_chop.F\
            tce_mo2e_zones_4a_disk_4s_new_patch.F\
            tce_mo2e_zones_4a_disk_2s_new.F tce_mo2e_zones_4a_disk_4s_new_patch_m2.F\
            tce_diagnose_t1.F ddotfile.F

LIB_DEFINES = -DDEBUG_PRINT
#
# Possible #defines
#
#              -DDEBUG_PRINT          /* Debug printing routines */
#	       -DNOCOMMS              /* Turn off ALL communication calls */
#              -DBLOCK_TRANSF         /* Use block-transfer for comms calls */
#              -DOLD_AOBLOCK          /* Recover original AO blocking scheme */
#

HEADERS = 

include ../config/makefile.h
include ../config/makelib.h

