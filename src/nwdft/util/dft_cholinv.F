      subroutine dft_cholinv(g_a,g_b)
c     cholesky followed by inversion
c     assumes square matrices
      implicit none
      integer g_a
      integer g_b
c
      integer info,itype
      integer n
      character*13 pname
#include "global.fh"
#include "dftps.fh"
#include "stdio.fh"
#ifdef SCALAPACK
#include "dft_scaleig.fh"
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
#endif
c     
      pname = "dft_cholinv:"
      if (odftps) call pstat_on(ps_chol)
      call ga_inquire(g_b, itype, n, n)
      call ga_zero(g_b)
      info = 0
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
      call ga_copy(g_a, g_b)
      call ga_sync()
      info= ga_cholesky('U',g_b)
#else
      call ga_chol(g_a, g_b, info)
#endif
#else 
      call ga_chol_seq(g_a, g_b, info)
#endif
      if (odftps) call pstat_off(ps_chol)
      if (info.ne.0)then
         if (ga_nodeid().eq.0)then
            write(LuOut,*)' Problem in performing a Cholesky '
            write(LuOut,*)' decomposition of the 2-ctr ERI '
            write(LuOut,*)' matrix using CD fitting basis. '
            write(LuOut,*)' Attempting a diag/inverse. '
         endif
      endif
      if (info.eq.0) then
      if (odftps) call pstat_on(ps_inv)
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
         info = ga_llt_i('U',g_b,-1)
         if (info.ne.0)then
            if (ga_nodeid().eq.0)then
               write(LuOut,*)' Problem in performing a Invers. '
               write(LuOut,*)' of the 2-ctr ERI '
            endif
            call ga_sync
            call errquit(pname//'Inverse failed ',0,0)
         endif
#else
         call ga_inverse(g_b, g_a)
#endif
         if (odftps) call pstat_on(ps_bdgemm)
#ifndef SCALAPACK
         call ga_dgemm('T', 'N', n, n, n, 1.d0, 
     &        g_a, g_a, 0.d0, g_b)
#endif
         if (odftps) call pstat_off(ps_bdgemm)
#else 
         call ga_copy(g_b, g_a)
         call ga_inv_seq(g_a, g_b)
#endif
      else
         call dft_invdiag(g_a, g_b, n)
      endif
      if (odftps) call pstat_off(ps_inv)
      return
      end
