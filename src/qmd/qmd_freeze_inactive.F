c
c     qmd_freeze_inactive: zero out velocities on inactive atoms
c
      subroutine qmd_freeze_inactive(n,m,v,act_at)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
#include "qmd_common.fh"
c
      integer n                ! number of atoms
      logical act_at(n)        ! active atoms
      double precision m(n)    ! atomic masses
      double precision v(3,n)  ! atomic velocities
      double precision KE_in, KE_after
c
      integer i,j
c
      character*32 pname
      pname='qmd_randomize: '
c
      call qmd_kinetic(n,m,v,KE_in)
      do i = 1,n
       if (.not.act_at(i)) then
        do j = 1,3
          v(j,i) = 0.d0
        end do
       end if
      end do
      call qmd_kinetic(n,m,v,KE_after)
c     rescale velocities to match kinetic energy on entry
      call qmd_rescale(n,m,v,KE_after,KE_in)
c
      return
      end
