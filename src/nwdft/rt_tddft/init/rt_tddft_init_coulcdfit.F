C
C     Initialization for charge density Coulomb fitting.
C
      subroutine rt_tddft_init_coulcdfit (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "rt_tddft.fh"
#ifdef SCALAPACK
            integer ga_cholesky,ga_llt_i
            external ga_cholesky,ga_llt_i
#endif


C     == In/out ==
      type(rt_params_t) params  !cdfit params stored in here

      
C     == Parameters ==
      integer, parameter   :: npol = 1
      character(*), parameter :: pname = "rt_tddft_init_coulcdfit: "


C     == External ==
      logical dft_mem3c
      external dft_mem3c

      
C     == Variables ==
      logical new_incore
      integer n_batch
      integer n3c_int, n3c_dbl, n_semi_bufs
      integer iwhat_max
      integer l_3ceri
      integer k_3ceri
      integer l_3cwhat
      integer k_3cwhat
      integer fd

C     (for part from dft_main0d.F)
      integer g_tmpb
      integer info
      logical IOLGC
      integer lmiss
      integer me


      me = ga_nodeid ()


C
C     (ripped from dft_main0d.F)
C
      IOLGC = .TRUE.

C      if (noio.eq.1) IOLGC = .FALSE.

      if (params%nodisk) then 
         call errquit (pname//"not working with nodisk yet",0,0)
         IOLGC = .false.
      endif


C      call errquit (pname//"XXX check npol, always 1?",0,0)
C      call rt_tddft_print_warning ("CHECK THAT CD FITTING WORKING")


      if (CDFIT)then
c     
c        Determine the characteristics of the CD Gaussian basis set.
c     
c     
c        Compute the matrix inverse of the CD 2-ctr ERIs.
c     
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERI', 
     &                       0, nbf_cd, g_2ceri)) 
     &      call errquit(pname//'Error creating g_2ceri',0,
     &       GA_ERR)
         call ga_zero(g_2ceri)
         call dft_get2eri(CD_bas_han, g_2ceri,oskel)
         if (oskel)call 
     .        sym_symmetrize(geom,cd_bas_han,.false.,g_2ceri)
         call ga_sync()
         if (.not. ga_duplicate(g_2ceri, g_cdinv, 'CD 2cERInv'))
     &    call errquit(pname//'Error creating g_cdinv',0, GA_ERR)
         call ga_zero(g_cdinv)
         lmiss = 1
c         if (lmiss.eq.1) then
            call ga_zero(g_cdinv)
            info = 0
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
            call ga_copy(g_2ceri, g_cdinv)
            call ga_sync()
            info= ga_cholesky('U',g_cdinv)
#else
            call ga_chol(g_2ceri, g_cdinv, info)
#endif
#else 
            call ga_chol_seq(g_2ceri, g_cdinv, info)
#endif
            if (info.ne.0)then
               if (me.eq.0)then
                  write(LuOut,*)' Problem in performing a Cholesky '
                  write(LuOut,*)' decomposition of the 2-ctr ERI '
                  write(LuOut,*)' matrix using CD fitting basis. '
                  write(LuOut,*)' Attempting a diag/inverse. '
               endif
            endif
            if (info.eq.0) then
               g_tmpb = g_2ceri
#if defined(PARALLEL_DIAG)
#ifdef SCALAPACK
               info = ga_llt_i('U',g_cdinv,-1)
            if (info.ne.0)then
               if (me.eq.0)then
                  write(LuOut,*)' Problem in performing a Invers. '
                  write(LuOut,*)' of the 2-ctr ERI '
               endif
               call ga_sync
               call errquit(pname//'Inverse failed ',0,0)
            endif
               
#else
               call ga_inverse(g_cdinv, g_tmpb)
#endif
#ifndef SCALAPACK
               call ga_dgemm('T', 'N', nbf_cd, nbf_cd, nbf_cd, 1.d0, 
     &              g_tmpb, g_tmpb, 0.d0, g_cdinv)
#endif
#else 
               call ga_copy(g_cdinv, g_tmpb)
               call ga_inv_seq(g_tmpb, g_cdinv)
#endif
            else
               call dft_invdiag(g_2ceri, g_cdinv, 
     &                          nbf_cd)
            endif
#ifndef SCALAPACK            
c
c     second build of g_2ceri needed becuase previous calls destroyed it
c     
            call ga_zero(g_2ceri)
            call dft_get2eri(CD_bas_han, g_2ceri,oskel)
            if (oskel)call 
     .           sym_symmetrize(geom,cd_bas_han,.false.,g_2ceri)
#endif
            if (IOLGC.and.(me.eq.0)) then
               lmiss = 0
               call dft_invio('CDI', g_cdinv, nbf_cd, 'WRITE', lmiss)
               if (lmiss.ne.0)call errquit
     &         (pname//'dft_invio - abnormal write of CDI ', 0,
     &       DISK_ERR)
               lmiss = 0
               call dft_invio('CD', g_2ceri, nbf_cd, 'WRITE', lmiss)
               if (lmiss.ne.0)call errquit
     &         (pname//'dft_invio - abnormal write of CD ', 0,
     &       DISK_ERR)
            endif
c         endif
         if (IOLGC) then
            if (.not. ga_destroy(g_cdinv)) call errquit
     &         (pname//'Could not destroy g_xcinv', 0, GA_ERR)
            if (.not. ga_destroy(g_2ceri)) call errquit
     &         (pname//'Could not destroy g_xcinv', 0, GA_ERR)
         endif
      endif
C
C     (end from dft_main0d.F)
C


C
C     Set up three center integrals.  This was already done in SCF, but
C     we have to redo it here.
C     
      if (dft_mem3c(params%rtdb,
     $     params%natoms, npol, .false., .false.,
     $     n3c_int, n3c_dbl, !! n_semi_bufs,
     $     l_3ceri, k_3ceri, l_3cwhat, k_3cwhat)) then
         
         incore=.false.
         call dft_3cincor (n_batch, n3c_int, int_mb(k_3cwhat),
     $        dbl_mb(k_3cERI), n3c_dbl)
         
         new_incore = .true.
      else
         new_incore = .false.
      endif

      
C
C     Replace the "incore" flag in cdft.fh
C
      incore = new_incore

C      if (new_incore .neqv. incore)
C     $     call errquit (pname//"incore different from SCF", 0, 0)


C
C     Load params necessary for CD fitting into "params" struct.
C
      params%n_batch = n_batch
      params%k_3ceri = k_3ceri
      params%l_3ceri = l_3ceri
      params%n3c_int = n3c_int
      params%k_3cwhat = k_3cwhat
      params%l_3cwhat = l_3cwhat
      params%n3c_dbl = n3c_dbl
      params%iwhat_max = iwhat_max
      params%n_semi_bufs = n_semi_bufs
      params%fd = fd

      
      end subroutine rt_tddft_init_coulcdfit



c $Id$
