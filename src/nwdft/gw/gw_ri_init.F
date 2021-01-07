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
      integer aopairs,mopairs,idum
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

      call ga_access(g_3ceris,ilo,ihi,rilo,rihi,k_eri,ld)
      call dfill(aopairs*nri_me,0d0,dbl_mb(k_eri),1)

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
      if (.not.ma_alloc_get(mt_int,lscr/nri,'Index',l_idx,k_idx))
     $  call errquit(pname//'can''t allocate index space',0,MA_ERR)
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
          call gw_ri_idx(int_mb(k_idx),ifirst,ilast,jfirst,jlast,
     $                   IeqJ)
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
            call gw_put_int(dbl_mb(offk),dbl_mb(k_eri),int_mb(k_idx),
     $                      nk,ni,nj,nij,aopairs,nri_me,kfirst,IeqJ)
          enddo
        enddo
      enddo
c
      ok = ma_free_heap(l_buf) .and. ma_free_heap(l_scr) .and.
     $     ma_free_heap(l_idx) 
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

      ! OO block
      call gw_puteris(dbl_mb(k_mo),nmo,nstates,nri,nri_me,rilo,1,nocc,
     $                1,nocc,pars%g_erioo,pars%idx_oo,'oo',idum,idum,
     $                idum)

      ! OV block
      call gw_puteris(dbl_mb(k_mo),nmo,nstates,nri,nri_me,rilo,1,nocc,
     $                nocc+1,nmo,pars%g_eriov,pars%idx_ov,'ov',
     $                pars%polelo,pars%polehi,pars%mynpoles)

      ! VV block
      call gw_puteris(dbl_mb(k_mo),nmo,nstates,nri,nri_me,rilo,nocc+1,
     $                nstates,nocc+1,nmo,pars%g_erivv,pars%idx_vv,'vv',
     $                idum,idum,idum)


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
      subroutine gw_ri_idx(idx,ifirst,ilast,jfirst,jlast,IeqJ)
      implicit none
      integer,intent(in) :: ifirst,jfirst,ilast,jlast
      logical,intent(in) :: IeqJ
      integer,intent(out),dimension(*) :: idx

      integer :: i, j, ij, jupper, itmp

      ij = 0
      jupper = jlast
      do i=ifirst,ilast
        itmp = (i*(i-1))/2
        if (IeqJ) jupper = i
        do j=jfirst,jupper
          ij = ij + 1
          idx(ij) = itmp + j
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
            ij = ij + i
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
      subroutine gw_puteris(erimo,nmo,nstates,nri,nri_me,rilo,ilower,
     $                      iupper,jlower,jupper,eri,idx,blk,polelo,
     $                      polehi,mynpoles)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"

      character(len=2),intent(in) :: blk
      integer,intent(in) :: nmo, nstates, nri, nri_me, rilo
      integer,intent(in) :: ilower, iupper, jlower, jupper
      integer,intent(out) :: eri,idx,polelo,polehi,mynpoles
      double precision,intent(in),dimension(nstates,nmo,nri_me) :: erimo

      character(len=12),parameter :: pname = 'gw_puteris: '
      character(len=6) :: erilabel,idxlabel

      integer ilo,ihi,jlo,jhi,imo,jmo,ipair,npairs,me,rihi
      integer ni,nj,l_tmp,k_tmp,k_loc,ld

      logical ok
      logical,external :: ga_create
      integer,external :: ga_nodeid, gw_get_pair

      integer l_oo,k_oo,l_ov,k_ov,l_vv,k_vv

      ni = iupper - ilower + 1
      if (ni.lt.1) return

      nj = jupper - jlower + 1
      rihi = rilo + nri_me - 1
      me = ga_nodeid()

      erilabel = 'eri '//blk
      idxlabel = 'idx '//blk

      npairs = ni*nj
      if (.not.ga_create(mt_dbl,nri,npairs,erilabel,nri,0,eri))
     $  call errquit(pname//'can''t create '//blk//' block',0,GA_ERR)
      if (.not.ma_alloc_get(mt_int,2*npairs,idxlabel,l_tmp,k_tmp))
     $  call errquit(pname//'can''t create '//blk//' index',0,MA_ERR)
      if (.not.ga_create(mt_int,npairs,2,idxlabel,0,2,idx))
     $  call errquit(pname//'can''t create '//blk//' index',0,GA_ERR)

      call gw_get_indeces(int_mb(k_tmp),npairs,ilower,iupper,jlower,
     $                    jupper)

      do ipair=1,npairs
        imo = gw_get_pair(int_mb(k_tmp),ipair,1,npairs)
        jmo = gw_get_pair(int_mb(k_tmp),ipair,2,npairs)
        call ga_put(eri,rilo,rihi,ipair,ipair,erimo(imo,jmo,1:nri_me),
     $              nri_me)
      enddo

      call ga_distribution(idx,me,ilo,ihi,jlo,jhi)
      call ga_access(idx,ilo,ihi,1,2,k_loc,ld)
      call dcopy(ld,int_mb(k_tmp+ilo-1),1,int_mb(k_loc),1)
      call dcopy(ld,int_mb(k_tmp+npairs+ilo-1),1,int_mb(k_loc+ld),1)
      ok = ma_free_heap(l_tmp)

      if (blk .eq. 'ov') then
        polelo = ilo
        polehi = ihi
        mynpoles = ihi - ilo + 1
      endif

      end subroutine

c     *****************************************************************      

      subroutine gw_get_indeces(idx,npairs,ilower,iupper,jlower,jupper)
      implicit none
      integer,intent(in) :: npairs,ilower,iupper,jlower,jupper

      integer,intent(out),dimension(npairs,2) :: idx

      integer iproc, ipair, jproc, nprocs
      integer imo, jmo
      integer,external :: ga_nnodes

      nprocs = ga_nnodes()

      ipair = 0
      do iproc=0,nprocs-1
        jproc = 0
        do imo=ilower,iupper
          do jmo=jlower,jupper
            if (jproc.eq.iproc) then
              ipair = ipair + 1
              idx(ipair,1) = imo
              idx(ipair,2) = jmo
            endif
            jproc = jproc + 1
            if (jproc.gt.nprocs-1) jproc = 0
          enddo
        end do
      enddo
      return
      end subroutine


      integer function gw_get_pair(idx,ipair,col,npairs)
      implicit none
      integer,intent(in) :: ipair,col,npairs
      integer,intent(in),dimension(npairs,2) :: idx
      gw_get_pair = idx(ipair,col)
      return
      end function
