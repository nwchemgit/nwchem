csh runtests.unix auh2o autosym dft_he2+ h2mp2 h2o hess_h2o prop_h2o pyqa
csh runtests.unix rimp2_ne sadsmall scf_feco5 small_intchk tagcheck testtab
csh runtests.unix h2o_dk u_sodft cosmo_h2o ch5n_nbo h2s_finite startag
csh runtests.unix neda_nbo

REM ---   small tests that should fail!
echo " "
echo "The oh2 test is testing the perl parsing script and SHOULD fail"
csh runtests.unix oh2
REM ---   medium tests
csh runtests.unix dft_feco5 
csh runtests.unix dft_ozone 
csh runtests.unix dft_siosi3
csh runtests.unix dielsalder
csh runtests.unix grad_ozone
csh runtests.unix hess_c6h6
csh runtests.unix intchk
csh runtests.unix sadbig
csh runtests.unix br2_dk
csh runtests.unix uo2_sodft
csh runtests.unix uo2_sodft_grad
csh runtests.unix si2cl6_gc
REM ---   long  may not run on workstations
csh runtests.unix aump2 
csh runtests.unix uoverlap
csh runtests.unix grad_nh3_trimer 
csh runtests.unix hess_nh3
csh runtests.unix hess_nh3_dimer
csh runtests.unix mp2_si2h6
csh runtests.unix pbo_nesc1e
