      subroutine bse_ri_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "bas.fh"
#include "bse.fh"
#include "stdio.fh"
#include "util.fh"      
#include "global.fh"
#ifdef SCALAPACK
#include "scaleig.fh"
#endif
c
      type(bse_params_t) :: pars
c
      character(len=*), parameter :: pname = "bse_ri_init: "
c
      logical iolgc
      integer g_2ceris, g_2cinv, info, lmiss
c
      integer ilo,ihi,jlo,jhi
      integer lbuf,lscr,ltmp,leri,lvecs
      integer l_buf,k_buf,l_scr,k_scr,l_tmp,k_tmp,k_eri,k_mo,l_mo
      integer l_vecs,k_vecs,l_idx,k_idx
      integer l_shlbuf, k_shlbuf, l_mobuf, k_mobuf
      integer me, isp, addr
      integer g_3ceris,shlo,shhi,offset,ld
      integer nmo,nvir,rihi,rilo,nri,nri_me,nocc,offk
      integer ish,jsh,ksh,ni,nj,nk,nij
      integer ifirst,jfirst,kfirst
      integer ilast,jlast,klast,imaxthreads
      integer aopairs,mopairs(2),idum,totmopairs
      logical IeqJ,ok,oactive
      double precision eri_est,timer,timing
c
      double precision schwarz_shell
      external schwarz_shell
c
#ifdef USE_OPENMP
      integer,external :: omp_get_max_threads
#endif

      timer = util_wallsec()
c
c     Create Global Arrays
c
      if (.not.ga_create(mt_dbl,nbf_cd,nbf_cd,'2c ERI',
     $     0, nbf_cd, g_2ceris))
     $  call errquit(pname//'Error creating 2c ERI',0,GA_ERR)
      if (.not.ga_create(mt_dbl,nbf_cd,nbf_cd,'2c INV',
     $     0, nbf_cd, g_2cinv))
     $  call errquit(pname//'Error creating 2c ERI',0,GA_ERR)
      call ga_zero(g_2ceris)
c
c     Compute 2-center ERIs
c
      if (pars%me.eq.0) write(luout,9000)
 9000 format(/,10x,'Computing 2-center integrals')
c
      call dft_get2eri(cd_bas_han, g_2ceris, .false.)
c
c     Obtain inverse Cholesky factor
c
      if (pars%me.eq.0) write(luout,9010)
 9010 format(10x,'Computing Inverse Cholesky factor')
c
      call ga_sync()
      call ga_chol_inv_seq(g_2ceris, g_2cinv, info)
      if (info.ne.0) then
        if (pars%me.eq.0) then
          write(luout,*)
          write(luout,*)' Problem in performing Cholesky'
          write(luout,*)' Obtaining square root via truncated '
          write(luout,*)' eigenvalue decomposition '
          write(luout,*)
        endif
        call dft_invdiag(g_2ceris, g_2cinv, nbf_cd)
      endif
c      
      if (.not. ga_destroy(g_2ceris))
     $  call errquit(pname//'failed to destrou g_2ceri',0,GA_ERR)
c
c     Allocate storage
c
      nri = nbf_cd
      nmo = pars%nmo
      aopairs = (nbf_ao*(nbf_ao+1))/2

      do isp=1,pars%ipol
        mopairs(isp) = pars%nmo*pars%nmo
      enddo
      totmopairs = sum(mopairs(1:pars%ipol))

      if (.not.ga_create(mt_dbl,1,nri,'3c ERI',aopairs,0,
     $     g_3ceris))
     $  call errquit(pname//'can''t create 3c ERIs',0,GA_ERR)
      call ga_distribution(g_3ceris,pars%me,ilo,ihi,rilo,rihi)
      nri_me = rihi - rilo + 1
      ok = ga_destroy(g_3ceris)

c      call ga_access(g_3ceris,ilo,ihi,rilo,rihi,k_eri,ld)
c      call dfill(aopairs*nri_me,0d0,dbl_mb(k_eri),1)

      shlo = 0
      do ksh=1,nshells_cd
        ok = bas_cn2bfr(cd_bas_han,ksh,kfirst,klast)
        if (klast.lt.rilo) cycle
        if (shlo.eq.0) then
          offset = (rilo-kfirst)
          shlo = ksh
        endif
        shhi = ksh
        if (klast.ge.rihi) exit
      enddo

      oactive = nri_me.gt.0

      call int_mem_2e3c(lbuf, lscr)
      if (.not.ma_alloc_get(mt_dbl,lbuf, 'RI buffer',l_buf, k_buf))
     $  call errquit(pname//'can''t get buffer space',0,MA_ERR)
      lscr = max(lscr,nbf_cd*nbf_ao_mxnbf_cn**2)
      if (.not.ma_alloc_get(mt_dbl,lscr, 'RI scratch', l_scr, k_scr))
     $  call errquit(pname//'can''t get scratch space',0,MA_ERR)

      if (.not.ma_alloc_get(mt_dbl,nri_me*nbf_ao_mxnbf_cn**2,
     $                      'Shell buffer',l_shlbuf, k_shlbuf))
     $  call errquit(pname//'can''t allocate shell buffer',0,MA_ERR)
      if (.not.ma_alloc_get(mt_dbl,
     $          nri_me*nbf_ao_mxnbf_cn*nmo*pars%ipol,
     $         'MO buffer',l_mobuf, k_mobuf))
     $  call errquit(pname//'can''t allocate shell buffer',0,MA_ERR)
c
c     Compute three-center ERIs in MO representation
c
      call ga_sync()
c
      if (pars%me.eq.0) write(luout,9020)
 9020 format(10x,'Computing 3-center integrals')

      lvecs = nbf_ao*nmo*pars%ipol
      if (.not.ma_alloc_get(mt_dbl,lvecs,'MOVECS',l_vecs,k_vecs))
     $  call errquit(pname//'can''t get movecs space',0,MA_ERR)
      if (.not.ma_alloc_get(mt_dbl,totmopairs*nri_me,'MO Eris',
     $                      l_mo,k_mo))
     $  call errquit(pname//'can''t get MO eris space',0,MA_ERR)
      call dfill(totmopairs*nri_me,0.0d0,dbl_mb(k_mo),1)

      do isp=1,pars%ipol
        call ga_get(pars%g_movecs(isp),1,nbf_ao,1,nmo,
     $              dbl_mb(k_vecs+(isp-1)*nbf_ao*nmo),nbf_ao)
      enddo
c
#ifdef USE_OPENMP
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif

      do ish=1,nshells_ao
        ok = bas_cn2bfr(ao_bas_han,ish,ifirst,ilast)
        ni = ilast - ifirst + 1
        call dfill(ni*nmo*pars%ipol*nri_me,0d0,dbl_mb(k_mobuf),1)
        do jsh=1,nshells_ao
          eri_est = schwarz_shell(ish,jsh)
          if (eri_est.lt.pars%tol2e) cycle
          IeqJ = ish.eq.jsh
          ok = bas_cn2bfr(ao_bas_han,jsh,jfirst,jlast)
          nj = jlast - jfirst + 1
          nij = ni*nj
          do ksh=shlo,shhi
            ok = bas_cn2bfr(cd_bas_han, ksh, kfirst, klast)
            call int_2e3c(cd_bas_han,ksh,ao_bas_han,ish,jsh,lscr,
     $                    dbl_mb(k_scr),lbuf,dbl_mb(k_buf))

            if (ksh.eq.shlo) then
              kfirst = rilo
              offk = k_buf + offset*nij
            else
              offk = k_buf
            endif

            if (ksh.eq.shhi) klast = rihi

            nk = klast - kfirst + 1
            kfirst = kfirst - rilo + 1
            call dcopy(ni*nj*nk,dbl_mb(offk),1,
     $                 dbl_mb(k_shlbuf+(kfirst-1)*ni*nj),1)
          enddo
          do isp=1,pars%ipol
            call gw_mo2(dbl_mb(k_shlbuf),
     $               dbl_mb(k_mobuf+(isp-1)*nmo*ni*nri_me),
     $               dbl_mb(k_vecs+(isp-1)*nbf_ao*nmo),nri_me,
     $               nbf_ao,ni,nj,nmo,nmo,ifirst,
     $               jfirst,ieqj)
          enddo
        enddo
        do isp=1,pars%ipol
          call gw_mo(dbl_mb(k_mobuf+(isp-1)*nmo*ni*nri_me),
     $                dbl_mb(k_mo+(isp-1)*mopairs(1)*nri_me),
     $                dbl_mb(k_vecs+(isp-1)*nbf_ao*nmo),nri_me,
     $                nbf_ao,ni,nmo,nmo,ifirst)
        enddo
      enddo
c
      ok = ma_free_heap(l_buf) .and. ma_free_heap(l_scr) .and.
     $     ma_free_heap(l_shlbuf) .and.
     $     ma_free_heap(l_mobuf)
      if (.not.ok) 
     $  call errquit(pname//'can''t free heap',165,MA_ERR)
c
c
#ifdef GWDEBUG
      if (pars%me.eq.0) write(*,*) '  Transform ERIs'
#endif

      call ga_sync()

      ok = ma_free_heap(l_vecs)
      if (.not.ok)
     $  call errquit(pname//'can''t free heap',176,MA_ERR)
c
      call ga_sync()
c
c     Distribute and orthonormalize ERIs
c
      timer = util_wallsec()
#ifdef GWDEBUG
      if (pars%me.eq.0) write(*,*) '  Distribute ERIs'
#endif

      addr = k_mo
      do isp=1,pars%ipol
        if (nmo.eq.0) cycle
        if (isp.eq.2) addr = k_mo + mopairs(1)*nri_me

        ! OO block
        call bse_puteris(dbl_mb(addr),nmo,nmo,nri,nri_me,
     $                  rilo,1,pars%nocc(isp),1,pars%nocc(isp),
     $                  pars%g_erioo(isp),'oo',pars%oolo(isp),
     $                  pars%oohi(isp),idum)

        ! OV block
        call bse_puteris(dbl_mb(addr),nmo,nmo,nri,nri_me,
     $                rilo,pars%nocc(isp)+1,nmo,1,pars%nocc(isp),
     $                pars%g_eriov(isp),'ov',pars%ovlo(isp),
     $                pars%ovhi(isp),pars%mynpoles(isp))

        ! VV block
        call bse_puteris(dbl_mb(addr),nmo,nmo,nri,nri_me,
     $                  rilo,pars%nocc(isp)+1,nmo,
     $                  pars%nocc(isp)+1,nmo,pars%g_erivv(isp),
     $                  'vv',pars%vvlo(isp),pars%vvhi(isp),idum)

      enddo

      ok = ma_free_heap(l_mo)
#ifdef GWDEBUG
      if (pars%me.eq.0) write(*,*) '  Read inverse'
#endif
c
      if (.not.ma_alloc_get(mt_dbl,nri**2,'RI temporary',l_tmp,k_tmp))
     $  call errquit(pname//'can''t get temp space',0,MA_ERR)
      if (pars%me.eq.0) then
        call ga_get(g_2cinv,1,nbf_cd,1,nbf_cd,dbl_mb(k_tmp),nbf_cd)
      endif
      call ga_sync()
      if (.not. ga_destroy(g_2cinv))
     $  call errquit(pname//'failed to destroy g_2cinv',0,GA_ERR)
      call ga_brdcst(1038,dbl_mb(k_tmp),ma_sizeof(mt_dbl,nri**2,
     $               mt_byte),0)

#ifdef GWDEBUG
      if (pars%me.eq.0) write(*,*) '  Orthogonalize ERIs'
#endif
      do isp=1,pars%ipol
        if (nmo.eq.0) cycle
c
c     Transform OO block
c
        call ga_distribution(pars%g_erioo(isp),pars%me,ilo,ihi,jlo,jhi)
        if (jlo.gt.jhi) goto 101
        call ga_access(pars%g_erioo(isp),ilo,ihi,jlo,jhi,k_eri,ld)
        nij = jhi-jlo+1
        call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $              nbf_cd,dbl_mb(k_eri),nbf_cd)
        call ga_release_update(pars%g_erioo(isp),ilo,ihi,jlo,jhi)
c
c     Transform OV block
c
 101    continue
        call ga_distribution(pars%g_eriov(isp),pars%me,ilo,ihi,jlo,jhi)
        if (jlo.gt.jhi) goto 102
        call ga_access(pars%g_eriov(isp),ilo,ihi,jlo,jhi,k_eri,ld)
        nij = jhi-jlo+1
        call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $              nbf_cd,dbl_mb(k_eri),nbf_cd)
        call ga_release_update(pars%g_eriov(isp),ilo,ihi,jlo,jhi)
c
c     Transform VV block
c
 102    continue
        if (nmo.gt.pars%nocc(isp)) then
          call ga_distribution(pars%g_erivv(isp),pars%me,ilo,ihi,
     $                         jlo,jhi)
          if (jlo.gt.jhi) goto 103
          call ga_access(pars%g_erivv(isp),ilo,ihi,jlo,jhi,k_eri,ld)
          nij = jhi-jlo+1
          call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $                nbf_cd,dbl_mb(k_eri),nbf_cd)
          call ga_release_update(pars%g_erivv(isp),ilo,ihi,jlo,jhi)
        endif
 103    continue
      enddo
c
c     Free heap
c
      ok = ma_free_heap(l_tmp)
      if (.not.ok)
     $  call errquit(pname//'failed to free heap',224,MA_ERR)
c
#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif
      return
      end subroutine
