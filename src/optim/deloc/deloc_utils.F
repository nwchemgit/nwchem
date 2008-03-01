c     ==========================================================
c     == Delocalized internals utility functions and routines ==
c     ==========================================================
c
c     == get the covalent radii of the atoms in the system ==
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
      double precision covrad(natom)
      double precision xtmp(3), qtmp, rtmp  ! local temp variables
c
      character*16 symb
c
c     == preliminaries ==
      if (natom.le.0) getcovalentradii = .false.
c
      write(luout,*) "natoms:",natom
c
c     == loop over the atoms and get the covalent radii ==
      do iatom = 1,natom
        covrad(iatom) = 0.d0  ! default radius
        if (geom_cent_get(igeom, iatom, symb, xtmp, qtmp)) then
          if (geom_tag_to_covalent_radius(symb, rtmp)) then
             covrad(iatom) = rtmp
             write(luout,*) iatom, covrad(iatom)
          end if
        end if
      end do
c
      getcovalentradii = .true.
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
      double precision grad(ndim)
      double precision rmsvalue
c
      double precision scalarprod
      external scalarprod
c
      rmsvalue = scalarprod(ndim,grad,grad)
      rmsvalue = dsqrt(rmsvalue/ndim)
c
      hasconverged = .false.
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
      double precision v1(n),v2(n)
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
