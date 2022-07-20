      program cif2xyz
      implicit none
      real*8, allocatable :: xyz(:,:)
      character*6, allocatable :: label(:)
      real*8 :: lattice_vec(3)
      integer natom, iatm, newatom
      read(*,*)natom
      read(*,*)lattice_vec
      allocate(xyz(3,4*natom),label(4*natom))
!     Read coordinates from "cif" file as fractional coordinates
      do iatm = 1, natom
        read(*,*)label(iatm),xyz(1:3,iatm)
      enddo
!     Convert coordinates to real space coordinates
      do iatm = 1, natom
        xyz(1,iatm) = xyz(1,iatm)*lattice_vec(1)
        xyz(2,iatm) = xyz(2,iatm)*lattice_vec(2)
        xyz(3,iatm) = xyz(3,iatm)*lattice_vec(3)
      enddo
!     Make a supercell
!     Add in the x-dimension
      newatom = natom
      do iatm = 1, natom
        newatom = newatom + 1
        label(newatom) = label(iatm)
        xyz(1,newatom) = xyz(1,iatm) + lattice_vec(1)
        xyz(2,newatom) = xyz(2,iatm)
        xyz(3,newatom) = xyz(3,iatm)
      enddo
      natom = newatom
!     Add in the z-dimension
      newatom = natom
      do iatm = 1, natom
        newatom = newatom + 1
        label(newatom) = label(iatm)
        xyz(1,newatom) = xyz(1,iatm)
        xyz(2,newatom) = xyz(2,iatm)
        xyz(3,newatom) = xyz(3,iatm) + lattice_vec(3)
      enddo
      natom = newatom
!     Write xyz file
      write(*,*)natom
      write(*,*)
      do iatm = 1, natom
        write(*,'(a6,3f24.10)')label(iatm),xyz(1:3,iatm)
      enddo
      end
