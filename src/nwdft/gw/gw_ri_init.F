      subroutine gw_ri_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "bas.fh"
#include "gw.fh"
#include "stdio.fh"
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
      integer ilo,ihi,jlo,jhi,g_a
      integer lbuf,lscr,ltmp,leri,lvecs,lvtmp,laomo
      integer l_buf,k_buf,l_scr,k_scr,l_tmp,k_tmp,l_eri,k_eri
      integer l_vecs,k_vecs,l_vtmp,k_vtmp,l_aomo,k_aomo
      integer cfirst,clast,dfirst,dlast,pfirst,plast,ip
      integer nfc,nfd,nfp,nij,eri_offset,ishc,ishd,ishp,nactive,tmp,me
      integer imo,istart,iend,l_map,k_map,np,nstates,npairs,nbf_cd_me
      logical CeqD,ok,oactive
      double precision eri_est
c
      logical ga_create, ga_duplicate, ga_destroy, ga_create_irreg
      external ga_create, ga_duplicate, ga_destroy, ga_create_irreg
c
      double precision schwarz_shell
      external schwarz_shell
c
      integer ga_nnodes
      external ga_nnodes
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
      write(luout,*)
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
      nstates = pars%nstates
      npairs = nstates*pars%nmo(1)

      if (.not.ma_alloc_get(mt_int,np,'ga map',l_map,k_map))
     &  call errquit(pname//'can''t get map space',0,MA_ERR)

      call ifill(np,pars%nmo(1)/np,int_mb(k_map),1)
      do ip=0,mod(pars%nmo(1),np)-1
        int_mb(k_map+ip) = int_mb(k_map+ip) + 1
      enddo

      oactive = int_mb(k_map+pars%me) > 0
      if (oactive) then
        if (pars%me.eq.0) then
          istart = 1
        else
          istart = sum(int_mb(k_map:k_map+pars%me-1)) + 1
        endif
        if (pars%me.eq.np-1) then
          iend = pars%nmo(1)
        else
          iend = sum(int_mb(k_map:k_map+pars%me))
        endif
      endif

      if (.not.oactive) goto 200

      leri = nbf_cd*(iend-istart+1)*nstates
      if (.not.ma_alloc_get(mt_dbl,leri, 'ERI buffer',l_eri,k_eri))
     $  call errquit(pname//'can''t allocate eri buffer',0,MA_ERR) 
      call int_mem_2e3c(lbuf, lscr)
      if (.not.ma_alloc_get(mt_dbl,lbuf, 'RI buffer',l_buf, k_buf))
     $  call errquit(pname//'can''t get buffer space',0,MA_ERR)
      lscr = max(lscr,nbf_cd*nbf_ao_mxnbf_cn**2)
      if (.not.ma_alloc_get(mt_dbl,lscr, 'RI scratch', l_scr, k_scr))
     $  call errquit(pname//'can''t get scratch space',0,MA_ERR)
      ltmp = max(nbf_cd*nbf_ao_mxnbf_cn**2,nbf_cd**2)
      if (.not.ma_alloc_get(mt_dbl,ltmp,'RI temporary',l_tmp,k_tmp))
     $  call errquit(pname//'can''t get temp space',0,MA_ERR)
      lvecs = nbf_ao*pars%nmo(1)
      if (.not.ma_alloc_get(mt_dbl,lvecs,'MOVECS',l_vecs,k_vecs))
     $  call errquit(pname//'can''t get movecs space',0,MA_ERR)
      lvtmp = nbf_ao_mxnbf_cn*pars%nmo(1)
      if (.not.ma_alloc_get(mt_dbl,lvtmp,'MOVECS tmp',l_vtmp,k_vtmp))
     $  call errquit(pname//'can''t get temp movecs space',0,MA_ERR)
      laomo = nbf_cd*nbf_ao_mxnbf_cn*nstates
      if (.not.ma_alloc_get(mt_dbl,laomo,'AOMO tmp',l_aomo,k_aomo))
     $  call errquit(pname//'can''t get temp aomo space',0,MA_ERR)

c
c     Compute three-center ERIs in MO representation
c
      if (pars%me.eq.0) write(luout,9020)
 9020 format(10x,'Computing 3-center integrals')
c
      call ga_get(pars%g_movecs(1),1,nbf_ao,1,pars%nmo(1),
     $            dbl_mb(k_vecs),nbf_ao)
c
      call dfill(leri,0.0d0,dbl_mb(k_eri),1)
      do ishc=1,nshells_ao
        if (.not.bas_cn2bfr(ao_bas_han, ishc, cfirst, clast))
     $    call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
        nfc = clast - cfirst + 1
c
        call dfill(laomo,0.0d0,dbl_mb(k_aomo),1)
        do ishd=1,nshells_ao
c
          eri_est = schwarz_shell(ishc,ishd)
          if (eri_est.lt.pars%tol2e*0.001d0) cycle
c
          if (.not.bas_cn2bfr(ao_bas_han, ishd, dfirst, dlast))
     $      call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
          nfd = dlast - dfirst + 1
c
          do ishp=1,nshells_cd
            if (.not.bas_cn2bfr(cd_bas_han, ishp, pfirst, plast))
     $        call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
            nfp = plast - pfirst + 1
c
            call int_2e3c(cd_bas_han, ishp, ao_bas_han,
     $                    ishc, ishd, lscr, dbl_mb(k_scr), lbuf,
     $                    dbl_mb(k_buf))
            call put_int(dbl_mb(k_buf), pfirst, nfp, nfc, nfd,
     $                   dbl_mb(k_tmp), nbf_cd, .false.)
          enddo
c          
          do imo=1,nstates
            call dcopy(nfd,dbl_mb(k_vecs+(imo-1)*nbf_ao+dfirst-1),
     $                 1,dbl_mb(k_vtmp+(imo-1)*nfd),1)
          enddo
c
          call gw_aomo(dbl_mb(k_tmp),nbf_cd,nfc,nfd,
     $                 dbl_mb(k_aomo),dbl_mb(k_vtmp),nstates)
c
        enddo
        do imo=1,pars%nmo(1)
          call dcopy(nfc,dbl_mb(k_vecs+(imo-1)*nbf_ao+cfirst-1),1,
     $               dbl_mb(k_vtmp+(imo-1)*nfc),1)
        enddo
        call gw_mo(dbl_mb(k_aomo),nbf_cd,nfc,nstates,pars%nmo(1),
     $             dbl_mb(k_eri),dbl_mb(k_vtmp),istart,iend)
      enddo
  200 continue
c
c     Orthonormalize 3-center ERIs
c
      if (pars%me.eq.0) write(luout,9030)
 9030 format(10x,'Orthonormalization of 3-center integrals')
c
      call ga_get(g_2cinv,1,nbf_cd,1,nbf_cd,dbl_mb(k_tmp),
     $            nbf_cd)
c
      nij = leri/nbf_cd
      call dtrmm('l','l','n','n',nbf_cd,nij,1.0d0,dbl_mb(k_tmp),
     $            nbf_cd,dbl_mb(k_eri),nbf_cd)
c
c     Destroy 2-ctr global array
c
      if (.not. ga_destroy(g_2cinv))
     $  call errquit(pname//'failed to destroy g_2cinv',0,GA_ERR)

      npairs = pars%nocc(1)**2
      if (.not.ga_create(mt_dbl,nbf_cd,npairs,'eris oo',nbf_cd,-1,
     $  pars%g_erioo))
     $  call errquit(pname//'can''t create O-O block',0,GA_ERR)

      npairs = pars%nocc(1)*pars%nvir(1)
      if (.not.ga_create(mt_dbl,nbf_cd,npairs,'eris ov',nbf_cd,-1,
     $  pars%g_eriov))
     $  call errquit(pname//'can''t create O-V block',0,GA_ERR)

      if (nstates.gt.pars%nocc(1)) then
        npairs = pars%nvir(1)*(nstates-pars%nocc(1))
        if (.not.ga_create(mt_dbl,nbf_cd,npairs,'eris vv',nbf_cd,-1,
     $    pars%g_erivv))
     $    call errquit(pname//'can''t create V-V block',0,GA_ERR)
      endif
c
c     Put MO ERIs in the global array
c
      if (.not.oactive) return

      call putmoeris_ga(pars%g_erioo,pars%g_eriov,pars%g_erivv,
     &                  dbl_mb(k_eri),pars%nocc(1),pars%nvir(1),nstates,
     &                  nbf_cd,pars%me,istart,iend)
c
c     Free heap
c
      ok = ma_free_heap(l_buf) .and. ma_free_heap(l_tmp) .and.
     $     ma_free_heap(l_scr) .and. ma_free_heap(l_vtmp) .and.
     $     ma_free_heap(l_vecs) .and. ma_free_heap(l_aomo) .and.
     $     ma_free_heap(l_eri)  .and. ma_free_heap(l_map)
      if (.not.ok)
     $  call errquit(pname//'failed to free heap',0,MA_ERR)
c
      return
      end subroutine
c
      subroutine gw_aomo(buf,nri,ni,nj,res,movecs,nstates)
      implicit none
      integer nri,ni,nj,nstates
      double precision buf(nri,nj,ni),res(nri,ni,nstates)
      double precision movecs(nj,nstates)
c
      integer i,j,s,imo,ib,ir
c
      do imo=1,nstates
        do i=1,ni
c          call dgemv('n',nri,nj,1.0d0,buf(1,ib+1),nri,movecs(1,imo),1,
c     $                1.0d0,res(1,ir),1)
          do j=1,nj
            res(:,i,imo) = res(:,i,imo) + buf(:,j,i)*movecs(j,imo)
          enddo
        enddo
      enddo
c
      return
      end subroutine
c
      subroutine gw_mo(buf,nri,ni,nstates,nmo,eri,movecs,istart,iend)
      implicit none
      integer nri,ni,nstates,nmo,istart,iend
      double precision buf(nri,ni,nstates),eri(nri,*)
      double precision movecs(ni,nmo)
      integer imo,jmo,i,ieri,ib

      ieri = 0
      do jmo=1,nstates
        do imo=istart,iend
          ieri = ieri + 1
          do i=1,ni
            eri(:,ieri) = eri(:,ieri) + buf(:,i,jmo)*movecs(i,imo)
          enddo
        enddo
      enddo
      return
      end subroutine

      subroutine putmoeris_ga(oo,ov,vv,eris,nocc,nvir,nstates,nri,me,
     $                        istart,iend)
      implicit none
      integer,intent(in) :: oo,ov,vv
      integer,intent(in) :: nocc,nvir,nstates,nri,istart,iend,me
      double precision,intent(in) :: eris(nri,*)

      integer myend,mystart,mynmo
      integer imo,jmo,ieri,ipair

      mynmo = iend - istart + 1

      ! Save occupied-occupied block
      if (istart.le.nocc) then
        myend = min(iend,nocc)
        do imo=1,nocc
          ieri = (imo-1)*mynmo
          ipair = (imo-1)*nocc + istart - 1
          do jmo=istart,myend
            ieri = ieri + 1
            ipair = ipair + 1
            call ga_put(oo,1,nri,ipair,ipair,eris(1,ieri),1)
          enddo
        enddo
      endif

      ! Save occupied-virtual block
      if (iend.ge.nocc+1) then
        mystart = max(nocc+1,istart)
        do imo=1,nocc
          ieri = (imo-1)*mynmo + mystart - istart
          ipair = (imo-1)*nvir + mystart - nocc - 1
          do jmo=mystart,iend
            ieri = ieri + 1
            ipair = ipair + 1
            call ga_put(ov,1,nri,ipair,ipair,eris(1,ieri),1)
          enddo
        enddo
      endif

      ! Save virtual-virtual block, if needed
      if (nstates.gt.nocc) then
        if (iend.ge.nocc+1) then
          mystart = max(nocc+1,istart)
          do imo=nocc+1,nstates
            ieri = (imo-1)*mynmo + mystart - istart
            ipair = (imo-nocc-1)*nvir + mystart - nocc - 1
            do jmo=mystart,iend
              ieri = ieri + 1
              ipair = ipair + 1
              call ga_put(vv,1,nri,ipair,ipair,eris(1,ieri),1)
            enddo
          enddo
        endif
      endif

      end subroutine
