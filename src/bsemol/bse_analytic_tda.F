      subroutine bse_analytic_tda(pars,wia,omega,maxpoles,totpoles,
     $                            singlet)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "bse.fh"
#ifdef SCALAPACK
#include "scaleig.fh"
      integer ga_cholesky,ga_llt_i
      external ga_cholesky,ga_llt_i
#endif

      type(bse_params_t) :: pars

      integer maxpoles, totpoles
      double precision wia(maxpoles,pars%ipol)
      double precision omega(totpoles)
      logical singlet

      character(*),parameter :: pname = 'bse_analytic_tda: '

      integer nri,isp,jsp,npoles
      integer ilo,ihi,jlo,jhi
      integer nocc2(2),nvir2(2)
      integer ipole,imo,jmo,amo,bmo
      integer k,l,kproc,lproc,klocal,llocal,kglobal,lglobal
      integer g_diag, g_tmp, g_amb, ktmp, ld, info, i, j
      integer modtn,tovern,off1,off2
      double precision alpha
      integer ij,jb,llia,ulia

#ifdef USE_OPENMP
      call util_blas_set_num_threads(pars%iMaxthreads)
#endif

      npoles = pars%npoles(1)
      if (pars%ipol.gt.1) npoles = npoles + pars%npoles(2)
      nri = pars%nri

      if(.not.ga_create(mt_dbl,totpoles,totpoles,'apb',0,npoles,
     &   pars%g_apb))
     &  call errquit(pname//'could not create A+B GA',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,npoles,'wia',0,g_diag))
     $  call errquit(pname//'could not create wia',0,GA_ERR)

      if (pars%ipol.eq.1) then
        nocc2(1) = pars%nocc(1)**2
        nvir2(1) = pars%nvir(1)**2

        call ga_distribution(pars%g_apb,pars%me,ilo,ihi,jlo,jhi)
        if (jhi.ge.jlo) call nga_put(g_diag,ilo,ihi,wia,ihi-ilo+1)


        ! Screened Coulomb contribution
        if(.not.ga_create(mt_dbl,nvir2(1),nocc2(1),'tmp',nvir2(1),0,
     &                    g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('t','n',nvir2(1),nocc2(1),pars%nri,-1d0,
     %                 pars%g_erivv(1),pars%g_erioo(1),0d0,g_tmp)
        call ga_distribution(g_tmp,pars%me,ilo,ihi,jlo,jhi)
        call ga_access(g_tmp,ilo,ihi,jlo,jhi,ktmp,ld)
        ! Need to take care about indexing
        do ij=jlo,jhi
          jmo = (ij-1)/pars%nocc(1) + 1
          imo = ij - (jmo-1)*pars%nocc(1)
          do bmo = 1,pars%nvir(1)
            llia = (imo-1)*pars%nvir(1) + 1
            ulia = imo*pars%nvir(1)
            jb = (jmo-1)*pars%nvir(1) + bmo
            call ga_put(pars%g_apb,llia,ulia,jb,jb,
     &           dbl_mb(ktmp+(ij-jlo)*nvir2(1)+(bmo-1)*pars%nvir(1)),
     &           pars%nvir(1))
          enddo
        enddo
        call ga_release(g_tmp,ilo,ihi,jlo,jhi)
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

        ! Hartree contribution
        if (singlet) then
          call ga_dgemm('t','n',npoles,npoles,pars%nri,2d0,
     &                   pars%g_eriov(1),pars%g_eriov(1),1d0,pars%g_apb)
        endif

        ! Diagonal term
        call ga_add_diagonal(pars%g_apb,g_diag)

        ! Destroy ERIs
        if (.not.ga_destroy(pars%g_erioo(1)))
     &    call errquit(pname//'could not destroy erioo GA',1,GA_ERR) 
        if (.not.ga_destroy(pars%g_erivv(1)))
     &    call errquit(pname//'could not destroy erivv GA',1,GA_ERR) 
        if (.not.ga_destroy(pars%g_eriov(1)))
     &    call errquit(pname//'could not destroy eriov GA',1,GA_ERR) 


        ! Diagonalize (A+B) matrix
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
        call dft_scaleig(pars%g_apb,pars%g_apb,omega,0)
#else
        call ga_diag_std(pars%g_apb,pars%g_apb,omega)
#endif
#else
        call ga_diag_std_seq(pars%g_apb,pars%g_apb,omega)
#endif

        ! Destroy unused GAs
        if (.not.ga_destroy(g_diag))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 
        
      else
        nocc2(1) = pars%nocc(1)**2
        nvir2(1) = pars%nvir(1)**2
        nocc2(2) = pars%nocc(2)**2
        nvir2(2) = pars%nvir(2)**2


        ! Gather eigenvalue differences
        call ga_distribution(pars%g_eriov(1),pars%me,ilo,ihi,jlo,jhi)
        call nga_put(g_diag,jlo,jhi,wia(1,1),jhi-jlo+1)
        call ga_distribution(pars%g_eriov(2),pars%me,ilo,ihi,jlo,jhi)
        call nga_put(g_diag,jlo+pars%npoles(1),jhi+pars%npoles(1),
     $               wia(1,2),jhi-jlo+1)

        ! Form RPA matrix
        if (.not.singlet) then
          call ga_zero(pars%g_apb)
        else
          do isp=1,pars%ipol
            off1 = (isp-1)*pars%npoles(1)
            do jsp=1,pars%ipol
              off2 = (jsp-1)*pars%npoles(1)
              call ga_matmul_patch('t','n',1d0,0d0,
     $         pars%g_eriov(isp),1,pars%npoles(isp),1,nri,
     $         pars%g_eriov(jsp),1,nri,1,pars%npoles(jsp),
     $         pars%g_apb,1+off1,pars%npoles(isp)+off1,
     $                    1+off2,pars%npoles(jsp)+off2)
            enddo
          enddo
        endif
        call ga_add_diagonal(pars%g_apb,g_diag)

        !Add Screened Coulomb contribution
        do isp=1,pars%ipol
          off1 = (isp-1)*pars%npoles(1)
          if (.not.ga_destroy(pars%g_eriov(isp)))
     &      call errquit(pname//'could not destroy eriov GA',147,GA_ERR)

          if (.not.ga_create(mt_dbl,nvir2(isp),nocc2(isp),'tmp',
     $                       nvir2(isp),0,g_tmp))
     $    call errquit('bse_analytic_tda: could not create tmp',145,
     $                  GA_ERR)

          call ga_dgemm('t','n',nvir2(isp),nocc2(isp),nri,-1d0,
     $      pars%g_erivv(isp),pars%g_erioo(isp),0d0,g_tmp)
          call ga_distribution(g_tmp,pars%me,ilo,ihi,jlo,jhi)
          call ga_access(g_tmp,ilo,ihi,jlo,jhi,ktmp,ihi-ilo+1)
          do ij=jlo,jhi
            jmo = (ij-1)/pars%nocc(isp) + 1
            imo = ij - (jmo-1)*pars%nocc(isp)
            do bmo = 1,pars%nvir(isp)
              llia = (imo-1)*pars%nvir(isp) + 1 + off1
              ulia = imo*pars%nvir(isp) + off1
              jb = (jmo-1)*pars%nvir(isp) + bmo + off1
              call ga_acc(pars%g_apb,llia,ulia,jb,jb,
     &             dbl_mb(ktmp+(ij-jlo)*nvir2(isp)+
     &                    (bmo-1)*pars%nvir(isp)),pars%nvir(isp),1d0)
            enddo
          enddo
          call ga_release(g_tmp,ilo,ihi,jlo,jhi)
          if (.not.ga_destroy(g_tmp))
     &      call errquit(pname//'could not destroy temp GA',172,GA_ERR) 
          ! Destroy ERIs
          if (.not.ga_destroy(pars%g_erioo(isp)))
     &      call errquit(pname//'could not destroy erioo GA',175,GA_ERR)
          if (.not.ga_destroy(pars%g_erivv(isp)))
     &      call errquit(pname//'could not destroy erivv GA',177,GA_ERR)
        enddo

        ! Diagonalize (A+B) matrix
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK
        call dft_scaleig(pars%g_apb,pars%g_apb,omega,0)
#else   
        call ga_diag_std(pars%g_apb,pars%g_apb,omega)
#endif
#else
        call ga_diag_std_seq(pars%g_apb,pars%g_apb,omega)
#endif

        ! Destroy unused GAs
        if (.not.ga_destroy(g_diag))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR)


      endif

      end subroutine
