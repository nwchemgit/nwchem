      subroutine xtb_guess(mol,qat,qshl,dpat,qpat,param,nat,nshl,
     $                     nshl_at,atnum)
      use mctc_io, only: structure_type
      use tblite_param, only : param_record
      use multicharge, only : mchrg_model_type, new_eeq2019_model,
     $                        get_covalent_rad, get_coordination_number,
     $                        get_lattice_points 
      implicit none
#include "util_intsize.fh"
      type(structure_type) :: mol
      type(param_record) :: param
      INTGR4 :: nat, nshl
      INTGR4 :: nshl_at(nat)
      INTGR4 :: atnum(nat)
      double precision :: qat(nat), qshl(nshl), dpat(3,nat), qpat(6,nat)

      type(mchrg_model_type) :: model
      integer :: i, j, ksh
      double precision :: tot, vec(3), tr
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


      do i=1,nat
        dpat(:,i) = qat(i)*mol%xyz(:,i)
        qpat(1,i) = qat(i)*mol%xyz(1,i)**2
        qpat(2,i) = qat(i)*mol%xyz(1,i)*mol%xyz(2,i)
        qpat(3,i) = qat(i)*mol%xyz(2,i)**2
        qpat(4,i) = qat(i)*mol%xyz(3,i)*mol%xyz(1,i)
        qpat(5,i) = qat(i)*mol%xyz(3,i)*mol%xyz(2,i)
        qpat(6,i) = qat(i)*mol%xyz(3,i)**2
      enddo

      dpat = 0d0
      qpat = 0d0

      deallocate(rcov,cn,trans)

      end subroutine


