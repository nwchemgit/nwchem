      subroutine xtb_guess(mol,qat,qshl,param,nat,nshl,nshl_at,atnum)
      use mctc_io, only: structure_type
      use tblite_param, only : param_record
      use multicharge, only : mchrg_model_type, new_eeq2019_model,
     $                        get_covalent_rad, get_coordination_number,
     $                        get_lattice_points 
      implicit none
      type(structure_type) :: mol
      type(param_record) :: param
      integer*4 :: nat, nshl
      integer*4 :: nshl_at(nat)
      integer*4 :: atnum(nat)
      double precision :: qat(nat), qshl(nshl)

      type(mchrg_model_type) :: model
      integer :: i, j, ksh
      double precision :: tot
      double precision, allocatable :: rcov(:), cn(:), trans(:,:)


      call new_eeq2019_model(mol, model)
      call get_lattice_points(mol%periodic,mol%lattice,40.0d0,trans)

      allocate(cn(nat))
      rcov = get_covalent_rad(mol%num)

      call get_coordination_number(mol, trans, 40.0d0, rcov, cn)

      ! get atomic charges guess
      call model%solve(mol, cn, qvec=qat)

      ! transform atomic charges into shell charges
      ksh = 0
      do i=1,nat
        tot = sum(param%record(atnum(i))%refocc(:))
        do j=1,nshl_at(i)
          ksh = ksh + 1
          qshl(ksh) = param%record(atnum(i))%refocc(j)/tot*qat(i)
        enddo
      enddo

      deallocate(rcov,cn,trans)

      end subroutine


