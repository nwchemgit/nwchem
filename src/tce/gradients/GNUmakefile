#$Id: GNUmakefile,v 1.3 2009-02-28 01:38:32 d3p975 Exp $

OBJ_OPTIMIZE = tce_gradient.o ccsd_energy.o ccsd_left.o ccsd_gradient.o \
               ccsd_1pdm_hh_mo.o ccsd_1pdm_pp_mo.o ccsd_1pdm_hp_mo.o ccsd_1pdm_ph_mo.o \
               ccsd_2pdm_hhhh_mo.o ccsd_2pdm_hhhp_mo.o ccsd_2pdm_hhpp_mo.o \
               ccsd_2pdm_hphh_mo.o ccsd_2pdm_hphp_mo.o ccsd_2pdm_hppp_mo.o \
               ccsd_2pdm_pphh_mo.o ccsd_2pdm_pphp_mo.o ccsd_2pdm_pppp_mo.o \
               zvec_x_a_a.o   zvec_x_a_b_a.o zvec_x_a_b_b.o zvec_x_a_c_a.o zvec_x_a_c_b.o \
               zvec_x_a_d.o   zvec_x_a_e_a.o zvec_x_a_e_b.o zvec_x_a_f_a.o zvec_x_a_f_b.o \
               zvec_x_b_a_a.o zvec_x_b_a_b.o zvec_x_b_b_a.o zvec_x_b_b_b.o zvec_x_b_c.o \
               zvec_x_b_d_a.o zvec_x_b_d_b.o zvec_x_b_e_a.o zvec_x_b_e_b.o zvec_x_b_f.o \
               zvec_x_c_a.o   zvec_x_c_b.o   ccsd_zvec.o \
               wdm_hh_mo_a.o   wdm_hh_mo_b_a_a.o wdm_hh_mo_b_a_b.o wdm_hh_mo_b_b_a.o wdm_hh_mo_b_b_b.o \
               wdm_hh_mo_c_a.o wdm_hh_mo_c_b_a.o wdm_hh_mo_c_b_b.o wdm_hh_mo_c_c_a.o wdm_hh_mo_c_c_b.o \
               wdm_hh_mo_c_d.o wdm_hh_mo_c_e_a.o wdm_hh_mo_c_e_b.o wdm_hh_mo_c_f_a.o wdm_hh_mo_c_f_b.o \
               wdm_hh_mo_d.o \
               wdm_hp_mo_a.o   wdm_hp_mo_b_a_a.o wdm_hp_mo_b_a_b.o wdm_hp_mo_b_b_a.o wdm_hp_mo_b_b_b.o \
               wdm_hp_mo_b_c.o wdm_hp_mo_b_d_a.o wdm_hp_mo_b_d_b.o wdm_hp_mo_b_e_a.o wdm_hp_mo_b_e_b.o \
               wdm_hp_mo_b_f.o wdm_hp_mo_c.o \
               wdm_pp_mo_a.o   wdm_pp_mo_b_a_a.o wdm_pp_mo_b_a_b.o wdm_pp_mo_b_b_a.o wdm_pp_mo_b_b_b.o \
               wdm_pp_mo_b_c.o wdm_pp_mo_b_d_a.o wdm_pp_mo_b_d_b.o wdm_pp_mo_b_e_a.o wdm_pp_mo_b_e_b.o \
               wdm_pp_mo_b_f.o \
               tce_grad_tidy.o tce_ccsd_1prdm.o \
               tce_dens_hhhh_offset.o tce_dens_hhhp_offset.o tce_dens_hhpp_offset.o \
               tce_dens_hphh_offset.o tce_dens_hphp_offset.o tce_dens_hppp_offset.o \
               tce_dens_pphh_offset.o tce_dens_pphp_offset.o tce_dens_pppp_offset.o \
               tce_zvec_x1_offset.o tce_2pdm_offset.o \
               put_2pdm_hhhh.o put_2pdm_hhhp.o put_2pdm_hhpp.o \
               put_2pdm_hphh.o put_2pdm_hphp.o put_2pdm_hppp.o \
               put_2pdm_pphh.o put_2pdm_pphp.o put_2pdm_pppp.o \
               tce_1pdm_offset.o put_1pdm_hh.o put_1pdm_hp.o put_1pdm_ph.o \
               put_1pdm_ph_2.o put_1pdm_pp.o \
               put_1pwdm_hp.o ao_tiles.o btrans1.o ao_1prdm_write.o \
               tce_gradients_doit.o tce_grad_force.o tce_grad1.o tce_grad2.o \
               sym_2pdm_ao.o btrans2.o \
               tce_dplot.o
                 

LIB_INCLUDES = -I../include

LIBRARY = libtce.a

USES_BLAS = ccsd_1pdm_hh_mo.F   ccsd_1pdm_pp_mo.F ccsd_1pdm_hp_mo.F ccsd_1pdm_ph_mo.F \
            ccsd_2pdm_hhhh_mo.F ccsd_2pdm_hhhp_mo.F ccsd_2pdm_hhpp_mo.F \
            ccsd_2pdm_hphh_mo.F ccsd_2pdm_hphp_mo.F ccsd_2pdm_hppp_mo.F \
            ccsd_2pdm_pphh_mo.F ccsd_2pdm_pphp_mo.F ccsd_2pdm_pppp_mo.F \
            zvec_x_a_a.F   zvec_x_a_b_a.F zvec_x_a_b_b.F zvec_x_a_c_a.F zvec_x_a_c_b.F \
            zvec_x_a_d.F   zvec_x_a_e_a.F zvec_x_a_e_b.F zvec_x_a_f_a.F zvec_x_a_f_b.F \
            zvec_x_b_a_a.F zvec_x_b_a_b.F zvec_x_b_b_a.F zvec_x_b_b_b.F zvec_x_b_c.F \
            zvec_x_b_d_a.F zvec_x_b_d_b.F zvec_x_b_e_a.F zvec_x_b_e_b.F zvec_x_b_f.F \
            zvec_x_c_a.F   zvec_x_c_b.F   ccsd_zvec.F \
            wdm_hh_mo_a.F   wdm_hh_mo_b_a_a.F wdm_hh_mo_b_a_b.F wdm_hh_mo_b_b_a.F wdm_hh_mo_b_b_b.F \
            wdm_hh_mo_c_a.F wdm_hh_mo_c_b_a.F wdm_hh_mo_c_b_b.F wdm_hh_mo_c_c_a.F wdm_hh_mo_c_c_b.F \
            wdm_hh_mo_c_d.F wdm_hh_mo_c_e_a.F wdm_hh_mo_c_e_b.F wdm_hh_mo_c_f_a.F wdm_hh_mo_c_f_b.F \
            wdm_hh_mo_d.F \
            wdm_hp_mo_a.F   wdm_hp_mo_b_a_a.F wdm_hp_mo_b_a_b.F wdm_hp_mo_b_b_a.F wdm_hp_mo_b_b_b.F \
            wdm_hp_mo_b_c.F wdm_hp_b_d_a.F wdm_hp_b_d_b.F wdm_hp_b_e_a.F wdm_hp_b_e_b.F \
            wdm_hp_mo_b_f.F wdm_hp_mo_c.F \
            wdm_pp_mo_a.F   wdm_pp_mo_b_a_a.F wdm_pp_mo_b_a_b.F wdm_pp_mo_b_b_a.F wdm_pp_mo_b_b_b.F \
            wdm_pp_mo_b_c.F wdm_pp_mo_b_d_a.F wdm_pp_mo_b_d_b.F wdm_pp_mo_b_e_a.F wdm_pp_mo_b_e_b.F \
            wdm_pp_mo_b_f.F 

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

include ../../config/makefile.h
include ../../config/makelib.h

