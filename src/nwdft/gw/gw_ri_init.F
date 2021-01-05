      subroutine gw_ri_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "bas.fh"
#include "gw.fh"
#include "stdio.fh"
#include "util.fh"      
#ifdef SCALAPACK
#include "dft_scaleig.fh"
      integer ga_cholesky, ga_llt_i
      external ga_cholesky, ga_llt_i
#endif
c
      type(gw_params_t) :: pars
c
      character(len=12), parameter :: pname = "gw_ri_init: "
c
      logical iolgc
      integer g_2ceris, g_2cinv, info, lmiss
c
      integer ilo,ihi,jlo,jhi
      integer lbuf,lscr,ltmp,leri,lvecs
      integer l_buf,k_buf,l_scr,k_scr,l_tmp,k_tmp,k_eri,k_mo,l_mo
      integer l_vecs,k_vecs,l_idx,k_idx
      integer me
      integer np,nstates
      integer g_3ceris,shlo,shhi,offset,ld
      integer nmo,nvir,rihi,rilo,nri,nri_me,nocc,offk
      integer ish,jsh,ksh,ni,nj,nk,nij
      integer ifirst,jfirst,kfirst
      integer ilast,jlast,klast
      integer aopairs,mopairs
      logical IeqJ,ok,oactive
      double precision eri_est,timer,timing
c
      logical ga_create, ga_duplicate, ga_destroy, ga_create_irreg
      external ga_create, ga_duplicate, ga_destroy, ga_create_irreg
c
      double precision schwarz_shell
      external schwarz_shell
c
      integer ga_nnodes
      external ga_nnodes

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
      if (pars%me.eq.0) write(luout,*)
c
      if (pars%me.eq.0) write(luout,9000)
 9000 format(10x,'Computing 2-center integrals')
c
      call dft_get2eri(cd_bas_han, g_2ceris, oskel)
      if (oskel) then
        call sym_symmetrize(geom, cd_bas_han,.false., g_2ceris)
      endif
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
      np = ga_nnodes()
      nri = nbf_cd
      nmo = pars%nmo(1)
      nocc = pars%nocc(1)
      nstates = pars%nstates
      mopairs = nstates*nmo
      aopairs = (nbf_ao*(nbf_ao+1))/2

      if (.not.ga_create(mt_dbl,aopairs,nri,'3c ERI',aopairs,0,
     $     g_3ceris))
     $  call errquit(pname//'can''t create 3c ERIs',0,GA_ERR)
      call ga_distribution(g_3ceris,pars%me,ilo,ihi,rilo,rihi)
      nri_me = rihi - rilo + 1

      call ga_access(g_3ceris,1,aopairs,rilo,rihi,k_eri,ld)
      call dfill(aopairs*nri_me,0d0,dbl_mb(k_eri),1)

      shlo = 0
      do ksh=1,nshells_cd
        ok = bas_cn2bfr(cd_bas_han,ksh,kfirst,klast)
        if (klast.lt.rilo) cycle
        if (shlo.eq.0) then
          offset = (kfirst-rilo)
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
      if (.not.ma_alloc_get(mt_int,lscr/nri,'Index',l_idx,k_idx))
     $  call errquit(pname/'can''t allocate index space',0,MA_ERR)
c
c     Compute three-center ERIs in MO representation
c
      call ga_sync()
c
      if (pars%me.eq.0) write(luout,9020)
 9020 format(10x,'Computing 3-center integrals')
c
      do ish=1,nshells_ao
        ok = bas_cn2bfr(ao_bas_han,ish,ifirst,ilast)
        ni = ilast - ifirst + 1
        do jsh=1,ish
          eri_est = schwarz_shell(ish,jsh)
          if (eri_est.lt.pars%tol2e*0.1d0) cycle
          IeqJ = ish.eq.jsh
          ok = bas_cn2bfr(ao_bas_han,jsh,jfirst,jlast)
          nj = jlast - jfirst + 1
          nij = ni*nj
          offk = k_buf+offset*nij
          call gw_ri_idx(int_mb(k_idx),ifirst,ilast,jfirst,jlast,nij,
     $                   IeqJ)
          do ksh=shlo,shhi
            ok = bas_cn2bfr(cd_bas_han, ksh, kfirst, klast)
            call int_2e3c(cd_bas_han,ksh,ao_bas_han,ish,jsh,lscr,
     $                    dbl_mb(k_scr),lbuf,dbl_mb(k_buf))
            if (ksh.eq.shlo) kfirst = rilo
            if (ksh.eq.shhi) klast = rihi
            nk = plast - kfirst + 1
            call gw_put_int(dbl_mb(offk),dbl_mb(k_eri),int_mb(k_idx),
     $                      nk,ni,nj,nij,aopairs,nri_me,kfirst,IeqJ)
          enddo
        enddo
      enddo
c
      ok = ma_free_heap(l_buf) .and. ma_free_heap(l_scr)
      if (.not.ok) 
     $  call errquit(pname//'can''t free heap',165,MA_ERR)
c
      lvecs = nbf_ao*nmo
      if (.not.ma_alloc_get(mt_dbl,lvecs,'MOVECS',l_vecs,k_vecs))
     $  call errquit(pname//'can''t get movecs space',0,MA_ERR)
      if (.not.ma_alloc_get(mt_dbl,mopairs*nri_me,'MO Eris',l_mo,k_mo))
     $  call errquit(pname//'can''t get MO eris space',0,MA_ERR)
c
      call ga_get(pars%g_movecs(1),1,nbf_ao,1,nmo,dbl_mb(k_vecs),nbf_ao)
      call gw_mo(dbl_mb(k_eri),dbl_mb(k_vecs),dbl_mb(k_mo),nbf_ao,nmo,
     $           nstates,aopairs,nri_me)

      ok = ma_free_heap(l_vecs) .and. ga_destroy(g_3ceris)
      if (.not.ok)
     $  call errquit(pname//'can''t free heap',176,MA_ERR)
c
      pars%timer%eris = util_wallsec() - timer
c
c     Distribute and orthonormalize ERIs
c
      timer = util_wallsec()
      call gw_puteris(pars,dbl_mb(k_mo),nmo,nstates,nri,nri_me,rilo)
      ok = ma_free_heap(l_mo)
c
      if (.not.ma_alloc_get(mt_dbl,nri**2,'RI temporary',l_tmp,k_tmp))
     $  call errquit(pname//'can''t get temp space',0,MA_ERR)
      call ga_get(g_2cinv,1,nbf_cd,1,nbf_cd,dbl_mb(k_tmp),nbf_cd)
      if (.not. ga_destroy(g_2cinv))
     $  call errquit(pname//'failed to destroy g_2cinv',0,GA_ERR)
c
c     Transform OO block
c
      call ga_distribution(pars%g_erioo,pars%me,ilo,ihi,jlo,jhi)
      call ga_access(pars%g_erioo,ilo,ihi,jlo,jhi,k_eri,ld)
      nij = jhi-jlo+1
      call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $            nbf_cd,dbl_mb(k_eri),nbf_cd)
c
c     Transform OV block
c
      call ga_distribution(pars%g_eriov,pars%me,ilo,ihi,jlo,jhi)
      call ga_access(pars%g_eriov,ilo,ihi,jlo,jhi,k_eri,ld)
      nij = jhi-jlo+1
      call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $            nbf_cd,dbl_mb(k_eri),nbf_cd)
c
c     Transform VV block
c
      if (nstates.gt.nocc) then
        call ga_distribution(pars%g_erivv,pars%me,ilo,ihi,jlo,jhi)
        call ga_access(pars%g_erivv,ilo,ihi,jlo,jhi,k_eri,ld)
        nij = jhi-jlo+1
        call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $              nbf_cd,dbl_mb(k_eri),nbf_cd)
      endif
c
      pars%timer%distr = util_wallsec() - timer
c
c     Free heap
c
      ok = ma_free_heap(l_tmp)
      if (.not.ok)
     $  call errquit(pname//'failed to free heap',224,MA_ERR)
c
      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_ri_idx(idx,ifirst,ilast,jfirst,jlast,nij,IeqJ)
      implicit none
      integer,intent(in) :: ifirst,jfirst,ilast,jlast,nij
      logical,intent(in) :: IeqJ
      integer,intent(out),dimension(nij) :: idx

      integer :: i, j, ij, jupper

      ij = 0
      jupper = jlast
      do i=ifirst,ilast
        idx(1:nij) = (i*(i-1))/2
        if (IeqJ) jupper = i
        do j=jfirst,jupper
          ij = ij + 1
          idx(ij) = idx(ij) + j
        enddo
      enddo
      return
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_put_int(buf,tmp,idx,nk,ni,nj,nij,aopairs,nri,kfirst,
     $                      IeqJ)
      implicit none
      integer,intent(in) :: nk, ni, nj, nij, aopairs, nri, kfirst
      integer,intent(in),dimension(nij) :: idx
      logical,intent(in) :: ieqj
      double precision,intent(in),dimension(nj,ni,nk) :: buf
      double precision,intent(out),dimension(aopairs,nri) :: tmp
c
      integer i, j, k, kri, ij

      if (IeqJ) then
        do k=1,nk
          kri = k + kfirst - 1
          ij = 0
          do i=1,ni
            tmp(idx(ij+1:ij+i),kri) = buf(1:i,i,k)
          enddo
        enddo
      else
        do k=1,nk
          kri = k + kfirst - 1
          ij = 0
          do i=1,ni
            tmp(idx(ij+1:ij+nj),kri) = buf(1:nj,i,k)
            ij = ij + nj
          enddo
        enddo
      endif
c
      end subroutine
c
c     *****************************************************************
c
      subroutine gw_mo(eriao,mos,erimo,nbf,nmo,nstates,npairs,nri)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
      integer,intent(in) :: nbf, nmo, nstates, npairs, nri
      double precision,intent(in),dimension(npairs,nri) :: eriao
      double precision,intent(in),dimension(nbf,nmo) :: mos
      double precision,intent(in),dimension(nstates,nmo,nri) :: erimo

      integer k_tmp, l_tmp, k_sym, l_sym
      integer istate, iri, ibf, istart, jstart

      if (.not.ma_push_get(mt_dbl,nstates*nbf,'temp',l_tmp,k_tmp))
     $  call errquit('gw_mo: cannot allocate temporart space',0,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf*nbf,'symm',l_sym,k_sym))
     $  call errquit('gw_mo: cannot allocate space',0,MA_ERR)
      call dfill(nmo*nmo,0.0d0,dbl_mb(k_sym),1)

      do iri=1,nri
c        
        do ibf=1,nbf
          istart = (ibf*(ibf-1))/2 + 1
          jstart = k_sym + nbf*(ibf-1)
          call dcopy(ibf,eriao(istart,iri),1,dbl_mb(jstart),1)
        enddo
c        
        call dsymm('l','u',nbf,nstates,1d0,dbl_mb(k_sym),nbf,mos,nbf,
     $              0d0,dbl_mb(k_tmp),nbf)
        call dgemm('t','n',nstates,nmo,nbf,1d0,dbl_mb(k_tmp),nbf,mos,
     $              nbf,0d0,erimo(1,1,iri),nstates)
c
      enddo
c
      if (.not.ma_chop_stack(l_tmp))
     $  call errquit('gw_mo: cannot chop stack',0,MA_ERR)

      return
      end subroutine
c
c     *****************************************************************
      subroutine gw_puteris(pars,erimo,nmo,nstates,nri,nri_me,rilo)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "gw.fh"

      type(gw_params_t) :: pars
      integer,intent(in) :: nmo, nstates, nri, nri_me, rilo
      double precision,intent(in),dimension(nstates,nmo,nri_me) :: erimo

      character(len=12),parameter :: pname = 'gw_puteris: '

      integer ilo,ihi,jlo,jhi,imo,jmo,ipair,npairs,nocc,nvir,rihi

      logical ok
      logical,external :: ga_create


      nocc = pars%nocc(1)
      nvir = nmo - nocc
      rihi = rilo + nri_me - 1

      npairs = nocc**2
      if (.not.ga_create(mt_dbl,nri,npairs,'eris oo',nri,0,
     $  pars%g_erioo))
     $  call errquit(pname//'can''t create OO block',0,GA_ERR)

      npairs = nocc*nvir
      if (.not.ga_create(mt_dbl,nri,npairs,'eris ov',nri,0,
     $  pars%g_eriov))
     $  call errquit(pname//'can''t create OV block',0,GA_ERR)

      if (nstates.gt.nocc) then
        npairs = nvir*(nstates-nocc)
        if (.not.ga_create(mt_dbl,nri,npairs,'eris vv',nri,0,
     $    pars%g_erivv))
     $    call errquit(pname//'can''t create VV block',0,GA_ERR)
      endif

      ! Fill OO block
      ipair = 0
      do imo=1,nocc
        do jmo=1,nocc
          ipair = ipair + 1
          call ga_put(pars%g_erioo,rilo,rihi,ipair,ipair,
     $                erimo(jmo,imo,:),nri)
        enddo
      enddo

      ! Fill OV block
      ipair = 0
      do imo=1,nocc
        do jmo=nocc+1,nmo
          ipair = ipair + 1
          call ga_put(pars%g_eriov,rilo,rihi,ipair,ipair,
     $                erimo(imo,jmo,:),nri)
        enddo
      enddo

      ! Fill VV block
      ipair = 0
      do imo=1,nocc
        do jmo=nocc+1,nmo
          ipair = ipair + 1
          call ga_put(pars%g_eriov,rilo,rihi,ipair,ipair,
     $                erimo(imo,jmo,:),nri)
        enddo
      enddo

      ! Fill VV block
      ipair = 0
      do imo=nocc+1,nstates
        do jmo=nocc+1,nmo
          ipair = ipair + 1
          call ga_put(pars%g_erivv,rilo,rihi,ipair,ipair,
     $                erimo(imo,jmo,:),nri)
        enddo
      enddo

      end subroutine

