      subroutine gw_coul_init(pars)
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
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
      character(len=14), parameter :: pname = "gw_coul_init: "
c
      logical oskel, iolgc
      integer g_2ceri, g_2cinv, info, lmiss
c
      integer lbuf,lscr,ltmp,leri,lvecs,lvtmp,laomo
      integer l_buf,k_buf,l_scr,k_scr,l_tmp,k_tmp,l_eri,k_eri
      integer l_vecs,k_vecs,l_vtmp,k_vtmp,l_aomo,k_aomo
      integer cfirst,clast,dfirst,dlast,pfirst,plast
      integer nfc,nfd,nfp,nij,eri_offset,ishc,ishd,ishp
      integer imo
      logical CeqD,ok
      double precision eri_est
c
      logical ga_create, ga_duplicate, ga_destroy
      external ga_create, ga_duplicate, ga_destroy
c
      double precision schwarz_shell
      external schwarz_shell
c
c     Create Global Arrays
c
      if (.not.ga_create(mt_dbl,pars%nbf_ri,pars%nbf_ri,'2c ERI',
     $     0, pars%nbf_ri, g_2ceri))
     $  call errquit(pname//'Error creating 2c ERI',0,GA_ERR)
      if (.not.ga_create(mt_dbl,pars%nbf_ri,pars%nbf_ri,'2c INV',
     $     0, pars%nbf_ri, g_2cinv))
     $  call errquit(pname//'Error creating 2c ERI',0,GA_ERR)
      call ga_zero(g_2ceri)
c
c     Compute 2-center ERIs
c
      write(luout,*)
c
      if (pars%me.eq.0) write(luout,9000)
 9000 format(10x,'Computing 2-center integrals')
c
      call dft_get2eri(pars%ri_bas_han, g_2ceri, oskel)
      if (oskel) then
        call sym_symmetrize(pars%geom, pars%ri_bas_han,.false., g_2ceri)
      endif
c
c     Obtain inverse Cholesky factor
c
      if (pars%me.eq.0) write(luout,9010)
 9010 format(10x,'Computing Inverse Cholesky factor')
c
      call ga_sync()
      call ga_chol_inv_seq(g_2ceri, g_2cinv, info)
      if (info.ne.0) then
        if (pars%me.eq.0) then
          write(luout,*)
          write(luout,*)' Problem in performing Cholesky'
          write(luout,*)' Obtaining square root via truncated '
          write(luout,*)' eigenvalue decomposition '
          write(luout,*)
        endif
        call dft_invdiag(g_2ceri, g_2cinv, pars%nbf_ri)
      endif
c
c     Allocate storage
c
      leri = pars%nbf_ri*(pars%nmo(1)*(pars%nmo(1)+1)/2)
      if (.not.ma_alloc_get(mt_dbl,leri,'3 centers',l_eri,k_eri))
     $  call errquit(pname//'can''t get eri space',0,MA_ERR)
      call int_mem_2e3c(lbuf, lscr)
      if (.not.ma_alloc_get(mt_dbl,lbuf, 'RI buffer',l_buf, k_buf))
     $  call errquit(pname//'can''t get buffer space',0,MA_ERR)
      lscr = max(lscr,pars%nbf_ri*pars%nbf_ao_mxnbf_cn**2)
      if (.not.ma_alloc_get(mt_dbl,lscr, 'RI scratch', l_scr, k_scr))
     $  call errquit(pname//'can''t get scratch space',0,MA_ERR)
      ltmp = max(pars%nbf_ri*pars%nbf_ao_mxnbf_cn**2,pars%nbf_ri**2)
      if (.not.ma_alloc_get(mt_dbl,ltmp,'RI temporary',l_tmp,k_tmp))
     $  call errquit(pname//'can''t get temp space',0,MA_ERR)
      lvecs = pars%nbf_ao*pars%nmo(1)
      if (.not.ma_alloc_get(mt_dbl,lvecs,'MOVECS',l_vecs,k_vecs))
     $  call errquit(pname//'can''t get movecs space',0,MA_ERR)
      lvtmp = pars%nbf_ao_mxnbf_cn*pars%nmo(1)
      if (.not.ma_alloc_get(mt_dbl,lvtmp,'MOVECS tmp',l_vtmp,k_vtmp))
     $  call errquit(pname//'can''t get temp movecs space',0,MA_ERR)
      laomo = pars%nbf_ri*lvtmp
      if (.not.ma_alloc_get(mt_dbl,laomo,'AOMO tmp',l_aomo,k_aomo))
     $  call errquit(pname//'can''t get temp aomo space',0,MA_ERR)

c
c     Compute three-center ERIs in MO representation
c
      if (pars%me.eq.0) write(luout,9020)
 9020 format(10x,'Computing 3-center integrals')
c
      call ga_get(pars%g_movecs(1),1,pars%nbf_ao,1,pars%nmo(1),
     $            dbl_mb(k_vecs),pars%nbf_ao)
c
      call dfill(leri,0.0d0,dbl_mb(k_eri),1)
      do ishc=1,pars%nshells_ao
        if (.not.bas_cn2bfr(pars%ao_bas_han, ishc, cfirst, clast))
     $    call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
        nfc = clast - cfirst + 1
c
        call dfill(laomo,0.0d0,dbl_mb(k_aomo),1)
        do ishd=1,pars%nshells_ao
c
          eri_est = schwarz_shell(ishc,ishd)
          if (eri_est.lt.pars%tol2e*0.001d0) cycle
c
          if (.not.bas_cn2bfr(pars%ao_bas_han, ishd, dfirst, dlast))
     $      call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
          nfd = dlast - dfirst + 1
c
          do imo=1,pars%nmo(1)
            call dcopy(nfd,dbl_mb(k_vecs+(imo-1)*pars%nbf_ao+dfirst-1),
     $                 1,dbl_mb(k_vtmp+(imo-1)*nfd),1)
          enddo
c
          do ishp=1,pars%nshells_ri
            if (.not.bas_cn2bfr(pars%ri_bas_han, ishp, pfirst, plast))
     $        call errquit(pname//'bas_cn2bfr failed',0,BASIS_ERR)
            nfp = plast - pfirst + 1
c
            call int_2e3c(pars%ri_bas_han, ishp, pars%ao_bas_han,
     $                    ishc, ishd, lscr, dbl_mb(k_scr), lbuf,
     $                    dbl_mb(k_buf))
            call put_int(dbl_mb(k_buf), pfirst, nfp, nfc, nfd,
     $                   dbl_mb(k_tmp), pars%nbf_ri, .false.)
          enddo
c
          call gw_aomo(dbl_mb(k_tmp),pars%nbf_ri,nfc,nfd,
     $                 dbl_mb(k_aomo),dbl_mb(k_vtmp),pars%nmo(1))
c
        enddo
        do imo=1,pars%nmo(1)
          call dcopy(nfc,dbl_mb(k_vecs+(imo-1)*pars%nbf_ao+cfirst-1),1,
     $               dbl_mb(k_vtmp+(imo-1)*nfc),1)
        enddo
        call gw_mo(dbl_mb(k_aomo),pars%nbf_ri,nfc,pars%nmo(1),
     $             dbl_mb(k_eri),dbl_mb(k_vtmp))
      enddo
c
c     Orthonormalize 3-center ERIs
c
      if (pars%me.eq.0) write(luout,9030)
 9030 format(10x,'Orthonormalization of 3-center integrals')
c
      call ga_get(g_2cinv,1,pars%nbf_ri,1,pars%nbf_ri,dbl_mb(k_tmp),
     $            pars%nbf_ri)
c
      nij = leri/pars%nbf_ri
      call dtrmm('l','l','n','n',pars%nbf_ri,nij,1.0d0,dbl_mb(k_tmp),
     $            pars%nbf_ri,dbl_mb(k_eri),pars%nbf_ri)
c
c     Destroy global arrays
c
      if (.not. ga_destroy(g_2cinv))
     $  call errquit(pname//'failed to destroy g_2cinv',0,GA_ERR)
      if (.not. ga_destroy(g_2ceri))
     $  call errquit(pname//'failed to destrou g_2ceri',0,GA_ERR)
c
c     Free heap
c
      ok = ma_free_heap(l_buf) .and. ma_free_heap(l_tmp) .and.
     $     ma_free_heap(l_scr) .and. ma_free_heap(l_vtmp) .and.
     $     ma_free_heap(l_vecs) .and. ma_free_heap(l_aomo)
      if (.not.ok)
     $  call errquit(pname//'failed to free heap',0,MA_ERR)
c
c      Save 3-center ERIs information
c
       pars%l_3ceri = l_eri
       pars%k_3ceri = k_eri
c
      return
      end subroutine
c
      subroutine gw_aomo(buf,nri,ni,nj,res,movecs,nmo)
      implicit none
      integer nri,ni,nj,nmo
      double precision buf(nri,nj,ni),res(nri,ni,nmo),movecs(nj,nmo)
c
      integer i,j,s,imo,ib,ir
c
      do imo=1,nmo
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
      subroutine gw_mo(buf,nri,ni,nmo,eri,movecs)
      implicit none
      integer nri,ni,nmo
      double precision buf(nri,ni,nmo),eri(nri,*),movecs(ni,nmo)
      integer imo,jmo,i,ieri,ib

      ieri = 0
      do imo=1,nmo
        do jmo=imo,nmo
          ieri = ieri + 1
          do i=1,ni
            eri(:,ieri) = eri(:,ieri) + buf(:,i,jmo)*movecs(i,imo)
          enddo
        enddo
      enddo
      return
      end subroutine
