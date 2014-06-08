c
c     Create the global array to hold R (rhs of linear system).
c
      subroutine tddft_grad_create_r(ipol,naoc,nav,nroot,g_r)
c
c     This array is stored as a 3 dimensional object with indeces
c     (ivec,iocc,ivirt) and the first dimension has a chunk size equal
c     to the number of vectors so that each processor holds the same
c     part of all vectors.
c
c     Written by Huub van Dam, May 2006.
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nroot     ! the number of roots to treat
      integer g_r(2)    ! global array handle for (X+Y)
c
c     Local
c
      integer idim(3)   ! dimension array for the nga calls
      integer ichunk(3) ! chunk size array for the nga calls
      integer i         ! just a counter
c
      do i = 1, ipol
        idim(1)   = nroot
        idim(2)   = naoc(i)
        idim(3)   = nav(i)
        ichunk(1) = nroot
        ichunk(2) = -1
        ichunk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vectors R',ichunk,
     +                      g_r(i)))
     +    call errquit('tddft_grad_create_r: failed to create g_r',
     +                 0, GA_ERR)
      enddo
c
      end
c
c     Destroy the global arrays holding R
c
      subroutine tddft_grad_destroy_r(ipol,g_r)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer ipol    ! =1 (restricted), =2 (unrestricted)
      integer g_r(2)  ! global array handle for (X+Y)
c
c     Local
c
      integer i ! just a counter
c
      do i = 1, ipol
        if (.not.ga_destroy(g_r(i))) call errquit(
     +         'tddft_grad_destroy_r: failed to destroy g_r',
     +         0, GA_ERR)
      enddo
c
      end
c $Id$
