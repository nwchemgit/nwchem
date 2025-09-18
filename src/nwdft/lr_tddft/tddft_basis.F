      subroutine tddft_basis(rtdb,geom_in,aobasis,lhascd)
      implicit none
#include "rtdb.fh"
#include "cdft.fh"
#include "errquit.fh"
#include "bas.fh"
#include "geom.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      logical lhascd
      integer rtdb, geom_in
      integer aobasis
      character(len=13),parameter :: pname = 'tddft_basis: '

      logical ok, iolgc
      logical, external :: int_normalize, int_norm_2c

      integer g_tmpb,lmiss,me,info

      lhascd = .false.
      me = ga_nodeid()
      iolgc = noio.ne.1

      geom = geom_in
c
c     Define the contracted Gaussian functions for the AO basis set.
c
      ok = bas_destroy(AO_bas_han)
      ok = bas_destroy(CD_bas_han)
      if (.not. bas_create(AO_bas_han, 'ao basis'))
     &   call errquit(pname//'bas_create?', 86, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, AO_bas_han, 'ao basis'))
     &   call errquit(pname//'no ao basis ', 86, BASIS_ERR)
      if (.not. int_normalize(rtdb,AO_bas_han))
     &   call errquit(pname//'int_normalize failed', 86, INT_ERR)
      if (.not. bas_numbf(AO_bas_han, nbf))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      nbf_ao = nbf
      if (.not. bas_nprim_cn_max(AO_bas_han,nbf_ao_mxprim))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_high_angular(AO_bas_han,nbf_ao_mxang))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_ncontr_cn_max(AO_bas_han,nbf_ao_mxcont))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_cn_max(AO_bas_han,nbf_ao_mxnbf_cn))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)
      if (.not. bas_numcont(AO_bas_han,nshells_ao))
     &   call errquit(pname//'basis set error:', 86, BASIS_ERR)

      aobasis = ao_bas_han


      if (.not. bas_create(CD_bas_han, 'cd basis'))
     &   call errquit(pname//'bas_create?',102, BASIS_ERR)
      CDFIT=bas_rtdb_load(rtdb, geom, CD_bas_han, 'cd basis')
      if (CDFIT)then
         call int_init(rtdb, 1, cd_bas_han)         
         if (.not. int_norm_2c(rtdb,CD_bas_han))
     &      call errquit('rdinput: int_norm_2c failed', 208, BASIS_ERR)
         call int_terminate()
         if (.not. bas_numbf(CD_bas_han, nbf_cd))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nprim_cn_max(CD_bas_han,nbf_cd_mxprim))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_high_angular(CD_bas_han,nbf_cd_mxang))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_ncontr_cn_max(CD_bas_han,nbf_cd_mxcont))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_cn_max(CD_bas_han,nbf_cd_mxnbf_cn))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_nbf_ce_max(CD_bas_han,nbf_cd_mxnbf_ce))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. bas_numcont(CD_bas_han,nshells_cd))
     &      call errquit('rdinput: basis set error:', 86, BASIS_ERR)
         if (.not. rtdb_put(rtdb, 'dft:cdfit', mt_log, 1, .true.))
     &      call errquit('rdinput: failed put_rtdb:', 99, RTDB_ERR)
c
c     texas breaks when AO_spherical and CD_cartesian
c
         if(bas_is_spherical(AO_bas_han).and.
     N        (.not.bas_is_spherical(CD_bas_han))) 
     T        call int_app_set_no_texas(rtdb)

      else
         if (.not. bas_destroy(CD_bas_han))
     &      call errquit('rdinput:not able to destroy CD_bas:',86,
     &       BASIS_ERR)
         CD_bas_han=-99999
         nbf_cd = 0
         nbf_cd_mxprim = 0
         nbf_cd_mxang = 0
         nbf_cd_mxcont = 0
         nbf_cd_mxnbf_cn = 0
         nbf_cd_mxnbf_ce = 0
         nshells_cd = 0
      endif

      if (cd_bas_han.eq.-99999) then
        call int_init(rtdb,1,ao_bas_han)
      else
        call int_init(rtdb,2,(/ao_bas_han,cd_bas_han/))
      endif

      if (CDFIT)then
         if (.not. ga_create(mt_dbl, nbf_cd, nbf_cd, 'CD 2cERI', 
     &                       0, nbf_cd, g_2ceri)) 
     &      call errquit(pname//'Error creating g_2ceri',0,
     &       GA_ERR)
         call ga_zero(g_2ceri)
         call dft_get2eri(CD_bas_han, g_2ceri,.false.)
         call ga_sync()
         if (.not. ga_duplicate(g_2ceri, g_cdinv, 'CD 2cERInv'))
     &        call errquit(pname//'Error creating g_cdinv',0, GA_ERR)
         call dft_cholinv(g_2ceri,g_cdinv)
#ifndef SCALAPACK            
c
c     second build of g_2ceri needed because previous calls destroyed it
c     
            call ga_zero(g_2ceri)
            call dft_get2eri(CD_bas_han, g_2ceri,.false.)
#endif
         if (IOLGC.and.(me.eq.0)) then
           lmiss = 0
           call dft_invio('CDI', g_cdinv, nbf_cd, 'WRITE', lmiss)
           if (lmiss.ne.0)call errquit
     &       (pname//'dft_invio - abnormal write of CDI ', 0,DISK_ERR)
           lmiss = 0
           call dft_invio('CD', g_2ceri, nbf_cd, 'WRITE', lmiss)
           if (lmiss.ne.0)call errquit
     &       (pname//'dft_invio - abnormal write of CD ', 0,DISK_ERR)
         endif
         if (IOLGC) then
            if (.not. ga_destroy(g_cdinv)) call errquit
     &         (pname//'Could not destroy g_xcinv', 0, GA_ERR)
            if (.not. ga_destroy(g_2ceri)) call errquit
     &         (pname//'Could not destroy g_xcinv', 0, GA_ERR)
         endif

         lhascd = .true.
      endif

      return
      end
