# $Id: nwchem_config_win32.h,v 1.9 2000-10-28 15:39:54 bjohnson Exp $

!IFDEF PYTHONHOME
PYTHON_NWLIB = python.lib
PYTHON_SYSLIB = /libpath:$(PYTHONHOME)\libs
PYTHON_SUBDIR = python
!ENDIF

NW_CORE_SUBDIRS = basis geom inp input \
	pstat rtdb task symmetry util peigs lapack blas

NW_MODULE_SUBDIRS = NWints atomscf ddscf gradients moints nwdft rimp2 stepper driver dftgrad cphf ccsd vib mcscf prepar nwargos esp hessian selci dplot mp2_grad property fft analyz nwmd cafe space drdy vscf qmmm qmd $(PYTHON_SUBDIR)

!IFNDEF LINK_F90
LINK_F90 = smathd.lib smaths.lib
!ENDIF

!IFNDEF NWCHEM_EXTRA_LIBS
NWCHEM_EXTRA_LIBS = cvwmpi.lib
!ENDIF

!IFDEF USECXML
CXMLLIB = cxml.lib
!ENDIF

UTIL_LIBS = util.lib $(PYTHON_NWLIB) pario.lib global.lib ma.lib peigs.lib \
tcgmsg-mpi.lib armci.lib $(CXMLLIB) lapack.lib blas.lib $(NWCHEM_EXTRA_LIBS) wsock32.lib \
$(LINK_F90)

LIBS = nwctask.lib ccsd.lib mcscf.lib selci.lib mp2.lib moints.lib stepper.lib driver.lib dftgrad.lib nwdft.lib gradients.lib cphf.lib esp.lib ddscf.lib guess.lib hessian.lib vib.lib util.lib rimp2.lib property.lib nwints.lib prepar.lib nwargos.lib nwmd.lib cafe.lib space.lib analyze.lib pfft.lib dplot.lib drdy.lib vscf.lib qmmm.lib qmd.lib $(UTIL_LIBS)

#EXCLUDED_SUBDIRS = develop ideaz scfaux hessian plane_wave oimp2 gapss pspw rimp2_grad nbo uccsdt
#CONFIG_LIBS = 
