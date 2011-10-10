c     ==========================================================
c     == Delocalized internals utility functions and routines ==
c     ==========================================================
c
c     == get the covalent radii of a collection of atoms ==
      logical function getcovalentradii(igeom, natom, covrad)
c
      implicit none
c
c     == general includes ==
#include "geom.fh"
#include "stdio.fh"
c
      integer igeom  ! geometry handle
      integer iatom, natom
c
      double precision covrad(*)
      double precision xtmp(3), qtmp, rtmp  ! local temp variables
c
      character*16 symb
c
c     == preliminaries ==
      getcovalentradii = .true.
c
c     == initial checks ==
      if (.not. geom_check_handle(igeom, "Check geometry handle")) then
       getcovalentradii = .false.
       return
      end if
      if (natom.le.0) getcovalentradii = .false.
c
c      write(luout,*) "natoms:",natom
c     == loop over the atoms and get the covalent radii ==
      do iatom = 1,natom
        covrad(iatom) = 0.d0  ! default radius
        if (geom_cent_get(igeom, iatom, symb, xtmp, qtmp)) then
          if (geom_tag_to_covalent_radius(symb, rtmp)) then
             covrad(iatom) = rtmp
c             write(luout,*) iatom, covrad(iatom)
          end if
        end if
      end do
c
      return
      end 
c
c     == get the atomic numbers of a collection of atoms ==
      logical function getatomicnumbers(igeom, natom, atnum)
c
      implicit none
c
c     == general includes ==
#include "geom.fh"
#include "stdio.fh"
c
      integer igeom  ! geometry handle
      integer iatom, natom
      integer atnum(*),atn
c
      double precision covrad(natom)
      double precision xtmp(3), qtmp, rtmp  ! local temp variables
c
      character*2 tag_symb
      character*16 tag_elem, tag
c
c     == preliminaries ==
      getatomicnumbers = .true.
c
c     == initial checks ==
      if (.not. geom_check_handle(igeom, "Check geometry handle")) then
       getatomicnumbers = .false.
       return
      end if
      if (natom.le.0) getatomicnumbers = .false.
c
c      write(luout,*) "natoms:",natom
c     == loop over the atoms and get the atomic numbers ==
      do iatom = 1,natom
        atnum(iatom) = 0    ! default atomic number
        if (geom_cent_get(igeom, iatom, tag, xtmp, qtmp)) then
          if (geom_tag_to_element(tag, tag_symb, tag_elem, atn)) then
             atnum(iatom) = atn
c             write(luout,*) iatom, atnum(iatom)
          end if
        end if
      end do
c
      return
      end
c
c     == convergence check ==
      logical function hasconverged(ndim,grad)
c
      implicit none
c
#include "deloc_data.fh"
c
      integer ndim
      double precision grad(*)
      double precision rmsvalue
c
      double precision scalarprod
      external scalarprod
c
c     == preliminaries ==
      hasconverged = .false.
      if (ndim.le.0) then
        hasconverged = .false.
        return
      end if
c
c     == calculate the rms value ==
      rmsvalue = scalarprod(ndim,grad,grad)
      rmsvalue = dsqrt(rmsvalue/ndim)
c
      if (rmsvalue.le.smallvalue) hasconverged = .true.
c
      return
      end
c
c     == finds the scalar product of 2 vectors ==
      double precision function scalarprod(n,v1,v2)
c
      implicit none
c
      integer i,n
      double precision v1(*),v2(*)
c
      double precision zero
      parameter (zero=0.0d0)
c
      scalarprod = zero
      do i=1,n
        scalarprod = scalarprod + v1(i)*v2(i)
      end do
c
      return
      end
c $Id$
