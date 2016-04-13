c
c     Create, Destroy global arrays (X+Y) and (X-Y).
c     HVD 5/2006, NG 11/2012
c
      subroutine tddft_grad_create_xy(tda,ipol,naoc,nav,nroot,
     +                                g_xpy,g_xmy)
c
c     These arrays are stored as 3 dimensional objects with indices
c     (ivec,iocc,ivirt) and the first dimension has a chunk size equal
c     to the number of vectors so that each processor holds the same
c     part of all vectors.
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
c
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer naoc(2)   ! the number of active occupied orbitals
      integer nav(2)    ! the number of active virtual orbitals
      integer nroot     ! the number of roots to treat
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
c     Local
c
      integer idim(3)   ! dimension array for the nga calls
      integer ichunk(3) ! chunk size array for the nga calls
      integer i         ! just a counter
      character*32 pname
c
      pname ="tddft_grad_create_xy: "

      if (ga_nodeid().eq.0) write(luout,*) "nroot: ", nroot
c
c     Check if full TDDFT or TDA
      if (tda) then    ! Tamm-Dancoff: B=0 and hence Y=0
c
c       Only need g_xpy (which equals g_x and therefore g_xmy)
        do i = 1, ipol
         idim(1)   = nroot
         idim(2)   = naoc(i)
         idim(3)   = nav(i)
         ichunk(1) = nroot
         ichunk(2) = -1
         ichunk(3) = -1
         if (.not.nga_create(mt_dbl,3,idim,'vectors X',ichunk,g_xpy(i)))
     +      call errquit(pname//'failed to create g_xpy',0,GA_ERR)
         call ga_zero(g_xpy(i))
        enddo
c
      else
c
c       Only need both g_xpy and g_xmy: Full non-Hermitian matrix
        do i = 1, ipol
         idim(1)   = nroot
         idim(2)   = naoc(i)
         idim(3)   = nav(i)
         ichunk(1) = nroot
         ichunk(2) = -1
         ichunk(3) = -1
         if (.not.nga_create(mt_dbl,3,idim,'vectors X+Y',ichunk,
     +             g_xpy(i)))
     +    call errquit(pname//'failed to create g_xpy',0,GA_ERR)
         call ga_zero(g_xpy(i))
c
         if (.not.nga_create(mt_dbl,3,idim,'vectors X-Y',ichunk,
     +             g_xmy(i)))
     +    call errquit(pname//'failed to create g_xmy',0,GA_ERR)
         call ga_zero(g_xmy(i))
        enddo
c
      endif
      end
c
      subroutine tddft_grad_destroy_xy(tda,ipol,g_xpy,g_xmy)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      logical tda       ! .true. if Tamm-Dancoff approximation is used
      integer ipol      ! =1 (restricted), =2 (unrestricted)
      integer g_xpy(2)  ! global array handle for (X+Y)
      integer g_xmy(2)  ! global array handle for (X-Y)
c
c     Local
c
      integer i ! just a counter
      character*32 pname
c
      pname = "tddft_grad_destroy_xy: "
c
      if (tda) then 
c
c       We only have g_xpy (which equals g_x and therefore g_xmy)
c
        do i = 1, ipol
         if (.not.ga_destroy(g_xpy(i))) 
     +     call errquit(pname//'failed to destroy g_xpy',0,GA_ERR)
        enddo
c
      else
c
c       We only need both g_xpy and g_xmy
c
        do i = 1, ipol
         if (.not.ga_destroy(g_xpy(i))) 
     +     call errquit(pname//'failed to destroy g_xpy',0,GA_ERR)
         if (.not.ga_destroy(g_xmy(i))) 
     +     call errquit(pname//'failed to destroy g_xmy',0,GA_ERR)
        enddo
c
      endif ! tda
      end
c $Id$
