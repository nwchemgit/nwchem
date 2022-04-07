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

      subroutine xtb_deallocate_params(param)
      use tblite_param, only : param_record
      implicit none
      type(param_record) :: param
      deallocate(param%name)
      deallocate(param%reference)
      deallocate(param%dispersion)
      deallocate(param%repulsion)
      deallocate(param%charge)
      deallocate(param%multipole)
      deallocate(param%halogen)
      deallocate(param%thirdorder)
      deallocate(param%record)
      end
