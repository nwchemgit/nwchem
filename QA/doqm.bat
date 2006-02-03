tcsh runtests.unix auh2o autosym dft_he2+ h2mp2 h2o hess_h2o prop_h2o pyqa
tcsh runtests.unix rimp2_ne sadsmall scf_feco5 small_intchk tagcheck testtab
tcsh runtests.unix h2o_dk u_sodft cosmo_h2o ch5n_nbo h2s_finite startag
tcsh runtests.unix cosmo_trichloroethene esp dplot dft_bsse
tcsh runtests.unix et_zn_dimer prop_uhf_h2o
tcsh runtests.unix vectors_rotate

REM ---   small tests that should fail!
echo " "
echo "The oh2 test is testing the perl parsing script and SHOULD fail"
tcsh runtests.unix oh2
REM ---   medium tests
tcsh runtests.unix dft_feco5 
tcsh runtests.unix dft_ozone 
tcsh runtests.unix dft_siosi3
tcsh runtests.unix dft_x
tcsh runtests.unix dft_mpw1k
tcsh runtests.unix dft_li2freq
tcsh runtests.unix dielsalder
tcsh runtests.unix grad_ozone
tcsh runtests.unix hess_c6h6
tcsh runtests.unix intchk
tcsh runtests.unix sadbig
tcsh runtests.unix br2_dk
tcsh runtests.unix uo2_sodft
tcsh runtests.unix uo2_sodft_grad
tcsh runtests.unix si2cl6_gc
tcsh runtests.unix pspw
tcsh runtests.unix pspw_SiC
tcsh runtests.unix pspw_md
tcsh runtests.unix pspw_polarizability
tcsh runtests.unix tddft_h2o
tcsh runtests.unix tddft_n2+
tcsh runtests.unix tce_h2o
tcsh runtests.unix oniom1
tcsh runtests.unix h2o_vscf
REM ---   long  may not run on workstations
tcsh runtests.unix aump2 
tcsh runtests.unix uoverlap
tcsh runtests.unix grad_nh3_trimer 
tcsh runtests.unix hess_nh3
tcsh runtests.unix hess_nh3_dimer
tcsh runtests.unix mp2_si2h6
tcsh runtests.unix pbo_nesc1e
tcsh runtests.unix oniom3
