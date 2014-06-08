      subroutine tddft_grad_store_gradient(rtdb,nat,nroots,iroots,g_g)
      implicit none
c
c     Write the final TDDFT gradient to the runtime database for use
c     elsewhere in the code.
c
c     Written by Huub van Dam, June 2007.
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
c
c     Input:
c
      integer rtdb           ! the runtime database handle
      integer nat            ! the number of atoms
      integer nroots         ! the number of excited states
      integer iroots(nroots) ! the excited state labels
      integer g_g            ! the global array holding the gradient
c
c     Local:
c
      integer lforce           ! the length of a gradient array
      integer l_force, k_force ! memory for the gradient array
      integer ilo(3), ihi(3)   ! the dimensions of the global array
      integer ld(3)            ! the leading dimensions of the array
      integer ir               ! counter over roots
      logical status           ! return value of RTDB calls
      integer iproc            ! the rank of the current processor
c
      character*26 rtdb_string ! name for the data on the RTDB
c
c     Code:
c
      call ga_sync
c
      status = rtdb_parallel(.false.)
c
      lforce = nat * 3 
      iproc = ga_nodeid()
      if (iproc.eq.0) then
c
        if (.not.ma_push_get(mt_dbl,lforce,'forces',l_force,k_force))
     &    call errquit(
     &    'tddft_grad_store_gradient:could not allocate l_force',lforce,
     &    MA_ERR)
c
        do ir = 1, nroots
c
          ilo(1) = ir
          ihi(1) = ir
          ilo(2) = 1
          ihi(2) = 3
          ilo(3) = 1
          ihi(3) = nat
          ld(1)  = 1
          ld(2)  = 3
          call nga_get(g_g,ilo,ihi,dbl_mb(k_force),ld)
c
c         write(rtdb_string,'("tddft:gradient:",i2)')iroots(ir)
          write(rtdb_string,'("tddft:gradient")')
          if (.not. rtdb_put(rtdb, rtdb_string, mt_dbl, lforce,
     &      dbl_mb(k_force)))call errquit
     &      ('tddft_grad_store_gradients: could not store gradients',
     &       1, RTDB_ERR)
c
        enddo
c
        if (.not.ma_pop_stack(l_force)) call errquit(
     &    'tddft_grad_store_gradient:could not deallocate l_force',
     &    0, MA_ERR)
c
      endif
c
      call ga_sync
c
      status = rtdb_parallel(.true.)
c
      end

c $Id$
