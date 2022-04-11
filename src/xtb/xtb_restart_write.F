      subroutine xtb_restart_write(mol, calc, wfn)
      use mctc_io, only : structure_type
      use tblite_xtb_calculator, only: xtb_calculator
      use tblite_wavefunction_type, only : wavefunction_type
      implicit none
#include "xtb.fh"

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(wavefunction_type) :: wfn

      integer i, unitno, ispin
      parameter(unitno=66)

      open(unit=unitno,status='unknown',form='unformatted',
     $     file=restart_out)
      do ispin=1,nspin
c        do i=1,nao
c          call swrite(unitno, dmat(1,i,ispin), int(calc%bas%nao))
c        enddo
        call swrite(unitno, wfn%qsh(1,ispin), int(calc%bas%nsh))
        call swrite(unitno, wfn%qat(1,ispin), int(mol%nat))
        call swrite(unitno, wfn%dpat(1,1,ispin), int(3*mol%nat))
        call swrite(unitno, wfn%qpat(1,1,ispin), int(6*mol%nat))
      enddo
      close(unit=unitno)
      end
