      subroutine bse_davidson_diagonal(pars,wia,maxpoles)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bse.fh"
#include "util.fh"

      type(bse_params_t) :: pars

      integer maxpoles, npoles
      double precision wia(maxpoles,pars%ipol)

      character(*),parameter :: pname = 'bse_casida: '

      integer isp,jsp,nri,off1,off2
      integer ilo,ihi,jlo,jhi
      integer nocc2(2),nvir2(2)
      integer ipole,imo,jmo,amo,bmo
      integer k,l,kproc,lproc,klocal,llocal,kglobal,lglobal
      integer g_diag, g_tmp, g_amb, ktmp, ld, info, i, j
      integer ij,jb,llia,ulia,ltmp2,ltmp3,ktmp2,ktmp3
      double precision temp(maxpoles),factor

      double precision, external :: ddot

#ifdef USE_OPENMP
      call util_blas_set_num_threads(pars%iMaxthreads)
#endif

      nri = pars%nri
      npoles = pars%npoles(1)
      if (pars%ipol.gt.1) npoles = npoles + pars%npoles(2)
      if (.not.nga_create(mt_dbl,1,npoles,'apb',0,pars%g_apb))
     $  call errquit(pname//'could not create wia',0,GA_ERR)
      if (.not.nga_create(mt_dbl,1,npoles,'amb',0,pars%g_amb))
     $  call errquit(pname//'could not create wia',0,GA_ERR)
      if(.not.ma_push_get(mt_dbl,nri,'dia',ltmp2,ktmp2))
     &  call errquit(pname//'failed to allocate Delta_ia',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nri,'dia',ltmp3,ktmp3))
     &  call errquit(pname//'failed to allocate Delta_ia',0,MA_ERR)

      factor = 2d0
      if ((pars%singlet) .and. (pars%ipol.eq.1)) factor = 4d0
      if ((pars%triplet) .and. (pars%ipol.eq.1)) factor = 0d0
      if (pars%tda) factor = factor/2d0

      do isp=1,pars%ipol
        off1 = (isp-1)*pars%npoles(1)

        ! Eigenvalue differences
        call ga_distribution(pars%g_eriov(isp),pars%me,ilo,ihi,jlo,jhi)
        call dcopy(jhi-jlo+1,wia(1,isp),1,temp,1)

        ! Screened Coulomb from A
        do ipole=jlo,jhi
          imo = (ipole-1)/pars%nvir(isp) + 1
          amo = ipole - (imo-1)*pars%nvir(isp)
          imo = imo + (imo-1)*pars%nocc(isp)
          amo = amo + (amo-1)*pars%nvir(isp)
          call ga_get(pars%g_erioo(isp),1,nri,imo,imo,dbl_mb(ktmp2),nri)
          call ga_get(pars%g_erivv(isp),1,nri,amo,amo,dbl_mb(ktmp3),nri)
          temp(ipole-jlo+1) = temp(ipole-jlo+1) - 
     $        ddot(nri,dbl_mb(ktmp2),1,dbl_mb(ktmp3),1)
        enddo
        call nga_put(pars%g_apb,jlo+off1,jhi+off1,temp,jhi-jlo+1)
        call nga_put(pars%g_amb,jlo+off1,jhi+off1,temp,jhi-jlo+1)

        ! Screened Coulomb from B
        if (.not.pars%tda) then
          call ga_access(pars%g_wov(isp),ilo,ihi,jlo,jhi,ktmp,ld)
          do ipole=jlo,jhi
            temp(ipole-jlo+1) = ddot(nri,dbl_mb(ktmp+(ipole-jlo)*nri),1,
     $                          dbl_mb(ktmp+(ipole-jlo)*nri),1)
          enddo
          call nga_acc(pars%g_apb,jlo+off1,jhi+off1,temp,jhi-jlo+1,-1d0)
          call nga_acc(pars%g_amb,jlo+off1,jhi+off1,temp,jhi-jlo+1,1d0)
          call ga_release(pars%g_wov(isp),ilo,ihi,jlo,jhi)
        endif

        ! Hartree term
        call ga_access(pars%g_eriov(isp),ilo,ihi,jlo,jhi,ktmp,ld)
        do ipole=jlo,jhi
          temp(ipole-jlo+1) =
     $       factor*ddot(nri,dbl_mb(ktmp+(ipole-jlo)*nri),1,
     $                       dbl_mb(ktmp+(ipole-jlo)*nri),1)
        enddo
        call nga_acc(pars%g_apb,jlo+off1,jhi+off1,temp,jhi-jlo+1,1d0)
        call ga_release(pars%g_eriov(isp),ilo,ihi,jlo,jhi)

      enddo

      if(.not.ma_chop_stack(ltmp2))
     $  call errquit('bse_davidson_diagonal: chop stack',0,MA_ERR)

      call ga_sync()

      end subroutine
