#$Id$

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
               eomccsd_energy.o eomccsd_left.o eomccsd_gradients.o \
               eomccsd_1pdm_hh_mo.o eomccsd_1pdm_pp_mo.o tce_eomccsd_1prdm.o \
               eomccsd_2pdm_hhhh_mo.o eomccsd_2pdm_hhhp_mo.o eomccsd_2pdm_hhpp_mo.o \
               eomccsd_2pdm_hphh_mo.o eomccsd_2pdm_hphp_mo.o eomccsd_2pdm_hppp_mo.o \
               eomccsd_2pdm_pphh_mo.o eomccsd_2pdm_pphp_mo.o eomccsd_2pdm_pppp_mo.o \
               eomccsd_z1.o eomccsd_z2.o eomccsd_2pdm_pphh_zvec_mo.o eomccsd_1pdm_hp_mo.o \
               eomccsd_1pdm_ph_mo.o eomccsd_x1_grad.o eomccsd_x2_grad.o tce_eom_xdiagon_grad.o
                 

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
            wdm_hp_mo_b_c.F wdm_hp_mo_b_d_a.F wdm_hp_mo_b_d_b.F wdm_hp_mo_b_e_a.F wdm_hp_mo_b_e_b.F \
            wdm_hp_mo_b_f.F wdm_hp_mo_c.F \
            wdm_pp_mo_a.F   wdm_pp_mo_b_a_a.F wdm_pp_mo_b_a_b.F wdm_pp_mo_b_b_a.F wdm_pp_mo_b_b_b.F \
            wdm_pp_mo_b_c.F wdm_pp_mo_b_d_a.F wdm_pp_mo_b_d_b.F wdm_pp_mo_b_e_a.F wdm_pp_mo_b_e_b.F \
            wdm_pp_mo_b_f.F \
            eomccsd_1pdm_hh_mo.F eomccsd_1pdm_pp_mo.F tce_eomccsd_1prdm.F \
            eomccsd_2pdm_hhhh_mo.F eomccsd_2pdm_hhhp_mo.F eomccsd_2pdm_hhpp_mo.F \
            eomccsd_2pdm_hphh_mo.F eomccsd_2pdm_hphp_mo.F eomccsd_2pdm_hppp_mo.F \
            eomccsd_2pdm_pphh_mo.F eomccsd_2pdm_pphp_mo.F eomccsd_2pdm_pppp_mo.F \
            eomccsd_z1.F eomccsd_z2.F eomccsd_2pdm_pphh_zvec_mo.F eomccsd_1pdm_hp_mo.F \
            eomccsd_1pdm_ph_mo.F eomccsd_x1_grad.F eomccsd_x2_grad.F tce_eom_xdiagon_grad.F \
            tce_cisdt_denominator.F cisdtq_2pdm_hhhh_mo.F cisdtq_2pdm_hphp_mo.F cisdtq_2pdm_hppp_mo.F \
            cisdtq_2pdm_pphp_mo.F cisdtq_2pdm_pppp_mo.F mp4_1pdm_pp_mo_III_b.F tce_cis_denominator.F \
            cisd_1pdm_hh_mo.F ccsdt_2pdm_hphh_mo.F ccsdt_2pdm_hhhp_mo.F mp3_1pdm_hh_mo.F \
            ccsdt_2pdm_pphh_mo.F tce_cisdtq_denominator.F tce_mp4_t24.F cis_x1.F \
            cisd_2pdm_hhhh_mo.F mp4_1pdm_hh_mo_II.F cisd_2pdm_hphp_mo.F cisdtq_1pdm_hh_mo.F \
            cisd_2pdm_hppp_mo.F cisd_2pdm_pphp_mo.F cis_1pdm_pp_mo.F cisd_2pdm_pppp_mo.F \
            ccsdt_1pdm_pp_mo.F cisdt_1pdm_pp_mo.F mp4_2pdm_hphh_mo.F mp2_1pdm_hh_mo.F \
            mp4_2pdm_pphh_mo.F ccsdtq_2pdm_hphh_mo.F ccsdtq_2pdm_hhhp_mo.F cis_2pdm_hphp_mo.F \
            tce_cisd_denominator.F ccsdtq_2pdm_pphh_mo.F cisdt_2pdm_hhhh_mo.F mp4_1pdm_hh_mo_III_a.F \
            ccsdtq_1pdm_pp_mo.F cisdt_2pdm_hphp_mo.F btrans1.F cisdt_2pdm_hppp_mo.F \
            cisdt_2pdm_pphp_mo.F cisdt_2pdm_pppp_mo.F cisdtq_2pdm_hphh_mo.F cisdtq_2pdm_hhhp_mo.F \
            wdm_hp_mo_b_d_a.F cisdtq_2pdm_pphh_mo.F cisdtq_2pdm_hhpp_mo.F mp4_1pdm_hh_mo_III_b.F \
            btrans2.F ccsdt_2pdm_hhhh_mo.F wdm_hp_mo_b_d_b.F ccsdt_2pdm_hphp_mo.F \
            cisd_1pdm_pp_mo.F ccsdt_2pdm_hppp_mo.F ccsdt_2pdm_pphp_mo.F mp3_1pdm_pp_mo.F \
            ccsdt_2pdm_pppp_mo.F tce_mp2_norm.F cisd_2pdm_hphh_mo.F cisd_2pdm_hhhp_mo.F \
            cis_1pdm_hh_mo.F cisd_2pdm_pphh_mo.F cisd_2pdm_hhpp_mo.F wdm_hp_mo_b_e_a.F \
            ccsdt_1pdm_hh_mo.F mp4_1pdm_pp_mo_II.F cisdt_1pdm_hh_mo.F cisdtq_1pdm_pp_mo.F \
            ccsdt_1pdm_ph_mo.F mp4_2pdm_hhhh_mo.F cis_xdiagon.F mp4_2pdm_hphp_mo.F \
            ccsdtq_2pdm_hhhh_mo.F mp4_2pdm_pphp_mo.F mp2_1pdm_pp_mo.F mp3_2pdm_hhhh_mo.F \
            ccsdtq_1pdm_hh_mo.F mp4_2pdm_pppp_mo.F ccsdtq_2pdm_hphp_mo.F wdm_hp_mo_b_e_b.F \
            ccsdtq_1pdm_ph_mo.F ccsdtq_2pdm_hppp_mo.F ccsdtq_2pdm_pphp_mo.F cisdt_2pdm_hphh_mo.F \
            cisdt_2pdm_hhhp_mo.F mp3_2pdm_hphp_mo.F ccsdtq_2pdm_pppp_mo.F cisdt_2pdm_pphh_mo.F \
            cisdt_2pdm_hhpp_mo.F mp4_1pdm_pp_mo_III_a.F mp3_2pdm_pppp_mo.F 


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

