      subroutine bse_davidson_kvec(pars,vec,wia,kv,npoles,ntrials,
     $                             dohartree,sgn)
      implicit none 
#include "bse.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"

      type(bse_params_t) :: pars

      integer npoles,ntrials
      integer vec,wia,kv,mv
      double precision sgn

      logical dohartree

      integer ipol,npoles1,npoles2,uli,lli,ipole,off1
      integer nocc(2),nvir(2),amo,imo,itrial,ula,lla,nri
      integer ilo,ihi,jlo,jhi,ld,ld2, myld

      integer l_tmp1,k_tmp1,k_eriov,isp,k_vec,k_erivv
      integer g_tmp1, g_tmp2, g_tmp3, g_tmp4, k_tmp
      integer subscript(2), k_tmp3
      double precision temp(npoles),factor

      ipol = pars%ipol
      nri = pars%nri
      nocc = pars%nocc
      nvir = pars%nvir
      npoles1 = nocc(1)*nvir(1)
      npoles2 = nocc(2)*nvir(2)

      call ga_copy(vec,kv)

      if (pars%tda .and. (sgn.eq.-1d0)) return

      call ga_scale_rows(kv,wia)

      ! Screened Coulomb A
      do isp=1,pars%ipol
        off1 = (isp-1)*pars%npoles(1)
        subscript = (/pars%me+(isp-1)*ga_nnodes(),0/)

        if(.not.ga_create(mt_dbl,pars%nvir(isp),pars%nocc(isp),'g_tmp1',
     $                    0,pars%nocc(isp),g_tmp1))
     $    call errquit('bse_kvec: ga_create failed',42,GA_ERR)
        if(.not.ga_create(mt_dbl,nri,pars%nvir(isp),'g_tmp2',
     $                    nri,0,g_tmp2))
     $    call errquit('bse_kvec: ga_create failed',45,GA_ERR)
        if(.not.ga_create(mt_dbl,nri*pars%nvir(isp),1,'g_tmp3',
     $                    0,1,g_tmp3))
     $    call errquit('bse_kvec: ga_create failed',48,GA_ERR)
        if(.not.ga_create(mt_dbl,pars%npoles(isp),1,'g_tmp4',
     $                    0,1,g_tmp4))   
     $    call errquit('bse_kvec: ga_create failed',48,GA_ERR)

        call nga_access_block_grid(vec,subscript,k_vec,ld2)
        call ga_distribution(pars%g_erivv(isp),pars%me,ilo,ihi,jlo,jhi)
        call ga_access(pars%g_erivv(isp),ilo,ihi,jlo,jhi,k_erivv,ld)
        call ga_access(g_tmp3,ilo,ihi,1,1,k_tmp3,ld)

        do itrial=1,ntrials
          call nga_copy_patch('n',vec,(/off1+1,itrial/),
     $            (/off1+pars%npoles(isp),itrial/),g_tmp4,
     $            (/1,1/),(/pars%npoles(isp),1/))
          call ga_copy_patch('n',g_tmp4,1,pars%npoles(isp),1,1,
     $            g_tmp1,1,pars%nvir(isp),1,pars%nocc(isp))
          do imo=1,pars%nocc(isp)
            lli = (imo-1)*pars%nocc(isp) + 1
            uli = imo*pars%nocc(isp)
            call ga_matmul_patch('n','t',1d0,0d0,
     $           pars%g_erioo(isp),1,nri,lli,uli,
     $           g_tmp1,1,pars%nocc(isp),1,pars%nvir(isp),
     $           g_tmp2,1,nri,1,pars%nvir(isp))
            call ga_copy_patch('n',g_tmp2,1,nri,1,pars%nvir(isp),
     $                             g_tmp3,1,nri*pars%nvir(isp),1,1)
            call dgemv('t',ld,pars%nvir(isp),-1d0,dbl_mb(k_erivv),ld,
     $       dbl_mb(k_tmp3),1,0d0,temp((imo-1)*pars%nvir(isp)+1),1)
          enddo
          call nga_acc(kv,(/1+off1,itrial/),
     $       (/pars%npoles(isp)+off1,itrial/),temp,npoles,1d0)
        enddo
        if(.not. (ga_destroy(g_tmp1).and.ga_destroy(g_tmp2).and.
     $            ga_destroy(g_tmp3).and.ga_destroy(g_tmp4)))
     $    call errquit('bse_davidson_kvec: ga_destroy_failed',68,GA_ERR)
        
        call ga_release(pars%g_erivv(isp),ilo,ihi,jlo,jhi)
        call ga_release_update(g_tmp3,ilo,ihi,1,1)
        call nga_release_block_grid(vec,subscript)
      enddo


      if (.not.ma_push_get(mt_dbl,nri*ntrials,'tmp1',l_tmp1,k_tmp1))
     $  call errquit('bse_davidson_kvec: MA push failed',73,MA_ERR)

      ! Hartree
      if (dohartree) then
        factor = 4d0/pars%ipol
        if (pars%tda) factor = factor/2d0
        call dcopy(nri*ntrials,0d0,0,dbl_mb(k_tmp1),1)
        do isp=1,pars%ipol
         off1 = (isp-1)*pars%npoles(1)
         subscript = (/pars%me+(isp-1)*ga_nnodes(),0/)
         call ga_distribution(pars%g_eriov(isp),pars%me,ilo,ihi,jlo,jhi)
         call ga_access(pars%g_eriov(isp),ilo,ihi,jlo,jhi,k_eriov,ld)
         call nga_access_block_grid(vec,subscript,k_vec,ld)
         call dgemm('n','n',nri,ntrials,ld,1d0,dbl_mb(k_eriov),nri,
     $                dbl_mb(k_vec),ld,1d0,dbl_mb(k_tmp1),nri)
         call ga_release(pars%g_eriov(isp),ilo,ihi,jlo,jhi)
         call nga_release_block_grid(vec,subscript)
        enddo
        call ga_dgop((/93/),dbl_mb(k_tmp1),nri*ntrials,'+')

        do isp=1,pars%ipol
         off1 = (isp-1)*pars%npoles(1)
         subscript = (/pars%me+(isp-1)*ga_nnodes(),0/)
         call ga_distribution(pars%g_eriov(isp),pars%me,ilo,ihi,jlo,jhi)
         call ga_access(pars%g_eriov(isp),ilo,ihi,jlo,jhi,k_eriov,ld)
         call nga_access_block_grid(kv,subscript,k_vec,ld)
         call dgemm('t','n',ld,ntrials,nri,factor,dbl_mb(k_eriov),nri,
     $               dbl_mb(k_tmp1),nri,1d0,dbl_mb(k_vec),ld)
         call ga_release(pars%g_eriov(isp),ilo,ihi,jlo,jhi)
         call nga_release_update_block_grid(kv,subscript)
        enddo
      endif

      ! Screened Coulomb B
      if(.not.pars%tda) then

      do isp=1,pars%ipol
        off1 = (isp-1)*pars%npoles(1)
        subscript = (/pars%me+(isp-1)*ga_nnodes(),0/)
        call ga_distribution(pars%g_wov(isp),pars%me,ilo,ihi,jlo,jhi)
        call ga_access(pars%g_wov(isp),ilo,ihi,jlo,jhi,k_eriov,ld)
        call nga_access_block_grid(vec,subscript,k_vec,ld)
        call dgemm('n','n',nri,ntrials,ld,-1d0,dbl_mb(k_eriov),nri,
     $              dbl_mb(k_vec),ld,0d0,dbl_mb(k_tmp1),nri)
        call nga_release_block_grid(vec,subscript)

        call ga_dgop((/115/),dbl_mb(k_tmp1),nri*ntrials,'+')

        call nga_access_block_grid(kv,subscript,k_vec,ld)
        call dgemm('t','n',ld,ntrials,nri,sgn,dbl_mb(k_eriov),nri,
     $              dbl_mb(k_tmp1),nri,1d0,dbl_mb(k_vec),ld)
        call ga_release(pars%g_wov(isp),ilo,ihi,jlo,jhi)
        call nga_release_update_block_grid(kv,subscript)
      enddo

      endif

      if (.not.ma_chop_stack(l_tmp1))
     $  call errquit('bse_kvec: failed chop stack',122,MA_ERR)

      end
