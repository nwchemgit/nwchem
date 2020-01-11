bash runtests.mpi.unix auh2o autosym dft_he2+ h2mp2 h2o hess_h2o prop_h2o pyqa
bash runtests.mpi.unix rimp2_ne sadsmall scf_feco5 small_intchk tagcheck testtab
bash runtests.mpi.unix h2o_dk u_sodft cosmo_h2o ch5n_nbo h2s_finite startag
bash runtests.mpi.unix cosmo_trichloroethene esp dplot dft_bsse
bash runtests.mpi.unix et_zn_dimer prop_uhf_h2o
bash runtests.mpi.unix vectors_rotate

REM ---   small tests that should fail!
echo " "
echo "The oh2 test is testing the perl parsing script and SHOULD fail"
bash runtests.mpi.unix oh2
REM ---   medium tests
bash runtests.mpi.unix dft_feco5 
bash runtests.mpi.unix dft_ozone 
bash runtests.mpi.unix dft_siosi3
bash runtests.mpi.unix dft_x
bash runtests.mpi.unix dft_mpw1k
bash runtests.mpi.unix dft_li2freq
bash runtests.mpi.unix dielsalder
bash runtests.mpi.unix grad_ozone
bash runtests.mpi.unix hess_c6h6
bash runtests.mpi.unix intchk
bash runtests.mpi.unix sadbig
bash runtests.mpi.unix br2_dk
bash runtests.mpi.unix uo2_sodft
bash runtests.mpi.unix uo2_sodft_grad
bash runtests.mpi.unix si2cl6_gc
bash runtests.mpi.unix pspw
bash runtests.mpi.unix pspw_SiC
bash runtests.mpi.unix pspw_md
bash runtests.mpi.unix pspw_polarizability
bash runtests.mpi.unix pspw_blyp_h2o
bash runtests.mpi.unix pspw_vs98_h2o
bash runtests.mpi.unix pspw_revpbe_h2o
bash runtests.mpi.unix pspw_pbesol_h2o
bash runtests.mpi.unix pspw_pbe0_h2o
bash runtests.mpi.unix pspw_hse_h2o
bash runtests.mpi.unix pspw_tpss03_h2o
bash runtests.mpi.unix pspw_scan_h2o
bash runtests.mpi.unix pspw_acm_h2o
bash runtests.mpi.unix pspw_becke97_h2o
bash runtests.mpi.unix pspw_becke98_h2o
bash runtests.mpi.unix pspw_hcth120_h2o
bash runtests.mpi.unix pspw_hcth147_h2o
bash runtests.mpi.unix pspw_hcth407_h2o
bash runtests.mpi.unix pspw_hcth_h2o  
bash runtests.mpi.unix pspw_mpw1k_h2o 
bash runtests.mpi.unix pspw_sogga_h2o
bash runtests.mpi.unix pspw_sogga11-x_h2o
bash runtests.mpi.unix pspw_b3lyp_h2o
bash runtests.mpi.unix pspw_beef_h2o
bash runtests.mpi.unix tddft_h2o
bash runtests.mpi.unix tddft_n2+
bash runtests.mpi.unix tce_h2o
bash runtests.mpi.unix oniom1
bash runtests.mpi.unix h2o_vscf
REM ---   long  may not run on workstations
bash runtests.mpi.unix aump2 
bash runtests.mpi.unix uoverlap
bash runtests.mpi.unix grad_nh3_trimer 
bash runtests.mpi.unix hess_nh3
bash runtests.mpi.unix hess_nh3_dimer
bash runtests.mpi.unix mp2_si2h6
bash runtests.mpi.unix pbo_nesc1e
bash runtests.mpi.unix oniom3
