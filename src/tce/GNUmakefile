#$Id$
#

OBJ = tce_input.o tce_init.o tce_energy.o tce_energy_fragment.o

OBJ_OPTIMIZE = tce_tile.o tce_ao1e.o tce_mo1e.o tce_ao2e.o \
      tce_filename.o createfile.o deletefile.o reconcilefile.o dratoga.o gatodra.o gatoeaf.o \
      get_block.o put_block.o add_block.o tce_tidy.o sf_test.o hbar.o \
      tce_mo2e.o tce_mo1e_offset.o tce_mo2e_offset.o tce_e_offset.o \
      tce_guess_t1.o tce_guess_t2.o tce_zero.o tce_print_x1.o tce_print_x2.o \
      tce_residual_t1.o tce_residual_t2.o tce_residual_t3.o tce_residual_t4.o \
      ma_pack.o sf_print.o ma_print_compact.o ma_scale.o \
      copyfile.o ddotfile.o daxfile.o daxpyfile.o dscalfile.o \
      tce_aod1.o tce_nud1.o tce_mod1.o tce_mod1_offset.o \
      sym_abelian_axis.o tce_mod1_new.o tce_mod1_offset_new.o \
      tce_eom_init.o tce_eom_tidy.o tce_eom_xdiagon.o tce_eom_ydiagon.o \
      tce_t1_offset.o tce_t2_offset.o tce_t3_offset.o tce_t4_offset.o \
      tce_x1_offset.o tce_x2_offset.o tce_x3_offset.o tce_x4_offset.o \
      tce_y1_offset.o tce_y2_offset.o tce_y3_offset.o tce_y4_offset.o \
      tce_mbpt2.o tce_restricted.o \
      tce_hash.o get_hash_block.o put_hash_block.o add_hash_block.o \
      tce_residual_t3a.o tce_t3a_offset.o tce_x3a_offset.o \
      tce_cholesky.o tce_ao2e_cholesky.o tce_ao2e_cholesky_2d.o tce_ao2e_cholesky_4d.o \
      e2_offset_2_offset.o get_block_ind.o get_block_ind_sf.o tce_hash_ind.o \
      tce_mo2e_offset_intorb.o tce_mo2e_offset_size.o tce_mo2e_trans.o \
      atpart_nalength.o tce_mo2e_zones_4a_disk_ga.o tce_mo2e_zones_4a_disk_ga_chop.o \
      tce_mo2e_zones_4a_disk_2s_new.o \
      tce_restart.o tce_restart_triples.o tce_restart_quadruples.o \
      tce_dens_hh_offset.o tce_dens_hp_offset.o tce_dens_pp_offset.o tce_dens_ph_offset.o \
      get_mo_rdm_hh.o get_mo_rdm_hp.o get_mo_rdm_ph.o get_mo_rdm_pp.o \
      ao_rdm_write.o tce_prop_start.o tce_property.o tce_diagnose_t1.o \
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
      tce_ccsd_perturbative.o tce_ccsd_lambda_perturbative.o \
      tce_ccsd_lambda_driver.o tce_ccsdt_lambda_driver.o tce_ccsdtq_lambda_driver.o \
      tce_ccsd_response_driver.o tce_ccsdt_response_driver.o tce_ccsdtq_response_driver.o \
      tce_ccsd_symm_polar.o tce_ccsdt_symm_polar.o tce_ccsdtq_symm_polar.o \
      tce_ccsd_lambda_response_driver.o tce_ccsd_asym_polar.o tce_ccsd_hyperpolar.o \
      tce_ao2e_new.o tce_mo2e_new.o tce_tile_new.o tce_init_new.o \
      tce_mo2e_zones_4a_disk_4s_new_patch.o tce_mo2e_zones_4a_disk_4s_new_patch_m2.o \
      tce_mo2e_zones_4a_disk_4s_new_patch_m3.o tce_ao2e_disk.o tce_mo2e_disk.o \
      tce_mo2e_disk_2eorb.o tce_mo2e_incore_2eorb.o tce_mo2e_incore_2eorb_split.o \
      tce_mo2e_hybrid_2eorb_split.o tce_mo2e_offset_intorb_hybrid.o \
      tce_mo2e_2emet3.o tce_mo2e_2emet3_n5.o tce_mo2e_2emet4.o tce_mo2e_2emet4_n5.o \
      get_hash_block_ma.o \
      ccsd_energy_loc.o cc2_energy.o ccsd_lambda.o ccsd_lr_alpha.o \
      ccsd_qr_beta.o tce_dgeev.o tce_schmidt.o tensor_read_write.o \
      tce_diis_act.o tce_mo2e_zones_4a_disk_ga_act.o tce_mo2e_zones_4a_disk_ga_chop_act.o \
      tce_t2a_offset.o tce_x2a_offset.o new_ga4ind_N5.o new4ind_N5_osu.o

OBJ_OPTIMIZE += tce_mo2e_zones_4a_disk_ga_N5.o tce_mo2e_zones_4a_disk_ga_chop_N5.o

ifdef USE_OPENMP_TASKS
    OBJ_OPTIMIZE += threadsafe_ga_wrappers.o
    #LIB_DEFINES += -DTCE_TIME_COMMS
endif

ifeq ($(NWCHEM_TARGET),$(findstring $(NWCHEM_TARGET),LAPI64 IBM64))
      #OBJ += tce_eom_xguess_karol.o tce_eom_yguess_karol.o
      OBJ += tce_eom_xguess.o tce_eom_yguess.o
else
      #OBJ_OPTIMIZE += tce_eom_xguess_karol.o tce_eom_yguess_karol.o
      OBJ_OPTIMIZE += tce_eom_xguess.o tce_eom_yguess.o
endif


ifdef MA_TRANS_BLOCKED
      LIB_DEFINES += -DMA_TRANS_BLOCKED
else
      OBJ_OPTIMIZE += ma_transpose.o
endif

LIB_INCLUDES = -I./include

LIBRARY = libtce.a

USES_BLAS = tce_mbpt2.F tce_residual_t1.F tce_residual_t2.F \
            tce_residual_t3.F tce_residual_t4.F tce_residual_t3a.F \
            tce_ao2e.F tce_mo2e.F tce_guess_t1.F tce_mod1.F  \
            tce_eom_xdiagon.F dscalfile.F tce_eom_ydiagon.F \
            tce_ao1e.F tce_aod1.F tce_mo1e.F \
            tce_mo2e_zones_4a_disk_ga.F tce_mo2e_zones_4a_disk_ga_chop.F \
            tce_mo2e_zones_4a_disk_4s_new_patch.F \
            tce_mo2e_zones_4a_disk_2s_new.F tce_mo2e_zones_4a_disk_4s_new_patch_m2.F \
            tce_mo2e_2emet3.F tce_mo2e_2emet3_n5.F tce_mo2e_2emet4.F tce_mo2e_2emet4_n5.F \
            tce_mo2e_disk_2eorb.F tce_mo2e_incore_2eorb.F tce_mo2e_incore_2eorb_split.F \
            tce_mo2e_hybrid_2eorb_split.F \
            tce_diagnose_t1.F ddotfile.F tce_ao2e_cholesky.F tce_energy.F \
            ccsd_lr_alpha.F tensor_read_write.F \
            tce_diis_act.F tce_mo2e_zones_4a_disk_ga_act.F tce_mo2e_zones_4a_disk_ga_chop_act.F \
            ccsd_qr_beta.F daxpyfile.F tce_ao2e_cholesky_orb.F tce_ao2e_disk.F \
            tce_ao2e_new.F tce_dgeev.F tce_mo2e_disk.F tce_mo2e_new.F \
            tce_mo2e_zones_4a_disk_4s_new_patch_m3.F tce_mo2e_zones_4a_disk_ga_N5.F \
            tce_mo2e_zones_4a_disk_ga_chop_N5.F tce_mod1_new.F new_ga4ind_N5.F new4ind_N5_osu.F \
            tce_read_t.F

SUBDIRS = include diis diis2 gradients \
          ccsd ccsd_lambda eomccsd \
          ccsd_t cr-eomccsd_t lrh emb \
          ccd lccd lccsd qcisd cc3 \
          mbpt1 mbpt2 mbpt3 mbpt4 cis cisd cisdt cisdtq \
          ccsd_act eomccsd_act cr_eomccsd_t_act \
	      sort ducc ucc4
          #response ccsd_yr ccsd_alpha ccsd_beta \

ifdef EACCSD
      SUBDIRS +=eaccsd
      LIB_DEFINES += -DEACCSD
endif
ifdef IPCCSD
      SUBDIRS +=ipccsd
      LIB_DEFINES += -DIPCCSD
endif

# to be backwards compatible, these must be on by default
CCSDLR=1
CCSDT=1

ifdef CCSDTLR
    CCSDT=1
    CCSDLR=1
    SUBDIRS += ccsdt_alpha_offdiag ccsdt_alpha_new ccsdt_yr ccsdt_beta
    LIB_DEFINES += -DCCSDTLR
endif

ifdef CCSDT
    SUBDIRS += ccsdt ccsdt_lambda eomccsdt
    LIB_DEFINES += -DCCSDT
endif

ifdef CCSDLR
    SUBDIRS += ccsd_alpha ccsd_beta ccsd_yr #ccsd_qr
    SUBDIRS += response
    LIB_DEFINES += -DCCSDLR
endif

ifdef CCSDTQ
    SUBDIRS += ccsdtq ccsdtq_lambda eomccsdtq ccsdtq_alpha_offdiag
    SUBDIRS += response
    LIB_DEFINES += -DCCSDTQ
endif

ifdef MRCC_METHODS
#ckbn mrcc-r-1 -3
      SUBDIRS += mrcc
      LIB_DEFINES += -DMRCC_METHODS
endif

ifdef TCE_CUDA
      LIB_DEFINES += -DTCE_CUDA
      ifdef USE_TTLG
        SUBDIRS += ttlg
      endif
endif

ifdef TCE_HIP
      LIB_DEFINES += -DTCE_HIP
endif

ifdef USE_PSTAT
      LIB_DEFINES += -DUSE_PSTAT
endif



#LIB_DEFINES += -DDEBUG_PRINT

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

ga_test: ga_test.o
	$(LINK.f) $< $(LIBS) -o $@
