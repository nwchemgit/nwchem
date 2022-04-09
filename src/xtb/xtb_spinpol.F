      subroutine xtb_spinpol(mol,calc,bas)
      use mctc_io, only: structure_type
      use tblite_container, only : container_type
      use tblite_spin, only : spin_polarization, new_spin_polarization
      use tblite_data_spin, only : get_spin_constant
      use tblite_basis_type, only : basis_type
      use tblite_xtb_calculator, only: xtb_calculator

      implicit none

      type(structure_type) :: mol
      type(xtb_calculator) :: calc
      type(basis_type)     :: bas
      type(spin_polarization), allocatable :: spin
      class(container_type),   allocatable :: cont

      integer i, ish, jsh
      double precision, allocatable :: wll(:,:,:)

      allocate(spin) 
      allocate(wll(bas%nsh,bas%nsh,mol%nid),source=0.0d0)
      do i=1,mol%nid
        do ish=1,bas%nsh_id(i)
          do jsh=1,bas%nsh_id(i)
            wll(ish,jsh,i) = get_spin_constant(
     $                           bas%cgto(ish,i)%ang,
     $                           bas%cgto(jsh,i)%ang,
     $                           mol%num(i))
          enddo
        enddo
      enddo
      call new_spin_polarization(spin, mol, wll, bas%nsh_id)
      call move_alloc(spin, cont)
      call calc%push_back(cont)

      !
      ! Both xTB objects should be deallocated by  now
      !
      if (allocated(spin)) deallocate(spin)
      if (allocated(cont)) deallocate(cont)

      end subroutine
