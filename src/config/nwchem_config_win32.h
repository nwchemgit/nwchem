#
#	$Id: nwchem_config_win32.h,v 1.1 1999-11-17 18:36:49 bjohnson Exp $
#

# comment out tools temporarily to get build working
#NW_CORE_SUBDIRS = tools basis geom inp input \
NW_CORE_SUBDIRS = basis geom inp input \
	pstat rtdb task symmetry util peigs lapack

NW_MODULE_SUBDIRS = NWints atomscf ddscf gradients moints nwdft rimp2 \
	stepper driver dftgrad cphf ccsd vib mcscf esp hessian selci dplot \
	mp2_grad property

UTIL_LIBS = util.lib chemio.lib global.lib ma.lib peigs.lib linalg.lib \
	tcgmsg-mpi.lib armci.lib cvwmpi.lib wsock32.lib lapack.lib smaths.lib

LIBS = nwctask.lib ccsd.lib mcscf.lib selci.lib mp2.lib moints.lib \
	stepper.lib driver.lib dftgrad.lib nwdft.lib gradients.lib cphf.lib \
	esp.lib ddscf.lib guess.lib hessian.lib vib.lib util.lib rimp2.lib \
	property.lib nwints.lib dplot.lib $(UTIL_LIBS)
