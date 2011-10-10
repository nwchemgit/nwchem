c     ==================================================
c     == Delocalized internals optimizer main routine ==
c     ==================================================
c
      subroutine deloc(periodic, natom, ndim,
     &                 coord_in, grad_in, covrad,
     &                 converged, ierror)
c
      implicit none
c
c     == general includes ==
#include "stdio.fh"
c
c     == optimizer specific includes ==
c
      logical periodic           ! periodic/molecular flag
      integer natom, ndim        ! number of atoms and dimensions
      double precision coord_in(*)  ! incoming coordinates
      double precision grad_in(*)   ! incoming gradients
      double precision covrad(*)    ! covalent radii
      logical converged             ! convergence flag
      integer ierror             ! return status (0: success, 1: failure)
c
      logical hasconverged
      external hasconverged
c
      integer i
c
c     == check convergence on entry ==
      converged = .false.
      ierror = 0
      if (hasconverged(ndim,grad_in)) then
         converged = .true.
         return
      end if
c
c     == test print out ==
      write(luout,*) "== natom ==",natom
      write(luout,*) "== coordinates =="
      do i = 1,natom
       write(luout,*) coord_in(3*i-2),coord_in(3*i-1),coord_in(3*i)
      end do
      write(luout,*) "== gradients =="
      do i = 1,natom
       write(luout,*) grad_in(3*i-2),grad_in(3*i-1),grad_in(3*i)
      end do
      write(luout,*) "== covalent radii =="
      do i = 1,natom
       write(luout,*) i, covrad(i)
      end do
c
      return
      end
c $Id$
