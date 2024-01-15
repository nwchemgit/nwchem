      subroutine xtb_restart_read(mol, calc, wfn)
      use mctc_io, only : structure_type
      use tblite_wavefunction_type, only : wavefunction_type
      use tblite_xtb_calculator, only: xtb_calculator
      implicit none
#include "xtb.fh"
      type(structure_type) :: mol
      type(wavefunction_type) :: wfn
      type(xtb_calculator) :: calc
      integer :: ispin,unitno

      parameter(unitno=66)

      open(unit=unitno,status='old',form='unformatted',file=restart_out)
      do ispin=1,nspin
        call sread(unitno, wfn%qsh(1,ispin), int(calc%bas%nsh))
        call sread(unitno, wfn%qat(1,ispin), int(mol%nat))
        call sread(unitno, wfn%dpat(1,1,ispin), int(3*mol%nat))
        call sread(unitno, wfn%qpat(1,1,ispin), int(6*mol%nat))
      enddo
      close(unit=unitno)
      return

      end
