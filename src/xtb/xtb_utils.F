      subroutine xtb_deallocate_wavefunction(wfn)
      use tblite_wavefunction_type, only : wavefunction_type
      implicit none
      type(wavefunction_type) :: wfn

      deallocate(wfn%homo)
      deallocate(wfn%nel)
      deallocate(wfn%n0at)
      deallocate(wfn%n0sh)
      deallocate(wfn%density)
      deallocate(wfn%coeff)
      deallocate(wfn%emo)
      deallocate(wfn%focc)
      deallocate(wfn%qat)
      deallocate(wfn%qsh)
      deallocate(wfn%dpat)
      deallocate(wfn%qpat)

      end subroutine

      subroutine xtb_deallocate_structure(struct)
      use mctc_io, only: structure_type
      implicit none
      type(structure_type) :: struct

      deallocate(struct%id)
      deallocate(struct%xyz)
      deallocate(struct%lattice)
      deallocate(struct%periodic)
      deallocate(struct%num)
      deallocate(struct%sym)

      end
