# This configuration generated automatically on cowboys at Wed Jul 12 14:26:50 PDT 2000
# Request modules from user: all

NW_CORE_SUBDIRS = basis geom inp input \
	pstat rtdb task symmetry util peigs lapack

NW_MODULE_SUBDIRS = NWints atomscf ddscf gradients moints nwdft rimp2 stepper driver dftgrad cphf ccsd vib mcscf prepar nwargos esp hessian selci dplot mp2_grad property fft analyz nwmd cafe space 

UTIL_LIBS = util.lib pario.lib global.lib ma.lib peigs.lib linalg.lib \
tcgmsg-mpi.lib armci.lib cvwmpi.lib wsock32.lib lapack.lib $(LINK_F90)

LIBS = nwctask.lib ccsd.lib mcscf.lib selci.lib mp2.lib moints.lib stepper.lib driver.lib dftgrad.lib nwdft.lib gradients.lib cphf.lib esp.lib ddscf.lib guess.lib hessian.lib vib.lib util.lib rimp2.lib property.lib nwints.lib prepar.lib nwargos.lib nwmd.lib cafe.lib space.lib analyze.lib pfft.lib dplot.lib $(UTIL_LIBS)

#EXCLUDED_SUBDIRS = develop ideaz scfaux hessian plane_wave oimp2 gapss pspw rimp2_grad nbo python vscf uccsdt drdy
#CONFIG_LIBS = 
