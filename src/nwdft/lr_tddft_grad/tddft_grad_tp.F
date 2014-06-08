c
c     Create and destroy global arrays to hold T, T+(X+Y), or T+V(=P).
c     HvD 5/2006, NG 11/2012
c
      subroutine tddft_grad_create_tp(ipol,naoc,nav,nroot,g_tp)
c
c     These arrays are stored as 3 dimensional objects with indeces
c     (ivec,imo,imo) and the first dimension has a chunk size equal
c     to the number of vectors so that each processor holds the same
c     part of all vectors.
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nroot     ! the number of roots to treat
      integer g_tp(2)   ! global array handle for T
c
c     Local
c
      integer idim(3)   ! dimension array for the nga calls
      integer ichunk(3) ! chunk size array for the nga calls
      integer i         ! just a counter
      character*32 pname
c
      pname ="tddft_grad_create_tp: "      
c
      do i = 1, ipol
        idim(1)   = nroot
        idim(2)   = naoc(i)+nav(i)
        idim(3)   = naoc(i)+nav(i)
        ichunk(1) = nroot
        ichunk(2) = -1
        ichunk(3) = -1
        if (.not.nga_create(mt_dbl,3,idim,'vectors T',ichunk,g_tp(i)))
     +    call errquit(pname//'failed to create g_tp',0,GA_ERR)
      enddo
c
      end
c
      subroutine tddft_grad_destroy_tp(ipol,g_tp)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer g_tp(2)   ! global array handle for (X+Y)
c
c     Local
c
      integer i ! just a counter
      character*32 pname
c
      pname="tddft_grad_destroy_tp: "
c
      do i = 1,ipol
        if (.not.ga_destroy(g_tp(i))) 
     +    call errquit(pname//'failed to destroy g_tp',0,GA_ERR)
      enddo
c
      end
c $Id$
