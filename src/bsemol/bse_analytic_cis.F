      subroutine bse_analytic_cis(pars,wia,omega,maxpoles,totpoles,
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

      logical singlet
      integer maxpoles, totpoles, npoles
      double precision wia(maxpoles,pars%ipol)
      double precision omega(totpoles)

      character(*),parameter :: pname = 'bse_casida: '

      integer isp,jsp,nri,off1,off2
      integer ilo,ihi,jlo,jhi
      integer nocc2(2),nvir2(2)
      integer ipole,imo,jmo,amo,bmo
      integer k,l,kproc,lproc,klocal,llocal,kglobal,lglobal
      integer g_diag, g_tmp, g_amb, ktmp, ld, info, i, j
      integer ij,jb,llia,ulia

#ifdef USE_OPENMP
      call util_blas_set_num_threads(pars%iMaxthreads)
#endif

      nri = pars%nri
      npoles = pars%npoles(1)
      if (pars%ipol.gt.1) npoles = npoles + pars%npoles(2)

      if(.not.ga_create(mt_dbl,npoles,npoles,'apb',0,npoles,
     &   pars%g_apb))
     &  call errquit(pname//'could not create A+B GA',0,GA_ERR) 
      if(.not.ga_create(mt_dbl,npoles,npoles,'amb',0,npoles,
     &   g_amb))
     &  call errquit(pname//'could not create A-B GA',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,npoles,'wia',0,g_diag))
     $  call errquit(pname//'could not create wia',0,GA_ERR)

      if (pars%ipol.eq.1) then
        nocc2(1) = pars%nocc(1)**2
        nvir2(1) = pars%nvir(1)**2

        call ga_distribution(pars%g_apb,pars%me,ilo,ihi,jlo,jhi)
        if (jhi.ge.jlo) call nga_put(g_diag,ilo,ihi,wia,ihi-ilo+1)

        ! Screened Coulomb contribution from A
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
          llia = (imo-1)*pars%nvir(1) + 1
          ulia = imo*pars%nvir(1)
          do bmo = 1,pars%nvir(1)
            jb = (jmo-1)*pars%nvir(1) + bmo
            call ga_put(pars%g_apb,llia,ulia,jb,jb,
     &           dbl_mb(ktmp+(ij-jlo)*nvir2(1)+(bmo-1)*pars%nvir(1)),
     &           pars%nvir(1))
          enddo
        enddo
        call ga_sync()
        call ga_copy(pars%g_apb,g_amb)
        call ga_release(g_tmp,ilo,ihi,jlo,jhi)
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

        ! Screened-Coulomb contribution from B
        if(.not.ga_create(mt_dbl,npoles,npoles,'tmp',0,npoles,
     &     g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('t','n',npoles,npoles,pars%nri,-1d0,
     &                 pars%g_wov(1),pars%g_wov(1),0d0,g_tmp)
        call ga_add(1d0,g_tmp,1d0,pars%g_apb,pars%g_apb)
        call ga_add(-1d0,g_tmp,1d0,g_amb,g_amb)
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

        if (singlet) then
        ! Hartree contribution
          call ga_dgemm('t','n',npoles,npoles,pars%nri,4d0,
     &                   pars%g_eriov(1),pars%g_eriov(1),1d0,pars%g_apb)
        endif

        ! Diagonal term
        call ga_add_diagonal(pars%g_apb,g_diag)
        call ga_add_diagonal(g_amb,g_diag)

        ! Destroy ERIs
        if (.not.ga_destroy(pars%g_erioo(1)))
     &    call errquit(pname//'could not destroy erioo GA',1,GA_ERR) 
        if (.not.ga_destroy(pars%g_erivv(1)))
     &    call errquit(pname//'could not destroy erivv GA',1,GA_ERR) 
        if (.not.ga_destroy(pars%g_eriov(1)))
     &    call errquit(pname//'could not destroy eriov GA',1,GA_ERR) 
        if (.not.ga_destroy(pars%g_wov(1)))
     &    call errquit(pname//'could not destroy wov GA',1,GA_ERR) 

        ! First step is to get (A-B)^(1/2)
        call ga_sync()
#ifdef SCALAPACK
        info= ga_cholesky('L',g_amb)
#elif defined(PEIGS)
        call ga_chol(g_amb, g_amb, info)
#else 
        call ga_chol_seq(g_amb, g_amb, info)
#endif

        ! Transform (A+B) matrix with Cholesky factors
        if(.not.ga_create(mt_dbl,npoles,npoles,'tmp',0,npoles,
     &     g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('n','n',npoles,npoles,npoles,1d0,
     &                 pars%g_apb,g_amb,0d0,g_tmp)
        call ga_dgemm('t','n',npoles,npoles,npoles,1d0,
     &                 g_amb,g_tmp,0d0,pars%g_apb)   
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

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

        !
        ! Since we transformed the eigenvalue equation to an Hermitian
        ! form, we obtain the square of the neutral excitations
        ! (Equation (32)).
        !
!$omp   parallel do
        do ipole=1,npoles
          omega(ipole) = dsqrt(omega(ipole))
        enddo
!$omp   end parallel do

        ! Transform eigenvectors with (A-B)^(1/2)
        if(.not.ga_create(mt_dbl,npoles,npoles,'tmp',0,npoles,
     &     g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('n','n',npoles,npoles,npoles,1d0,
     &                 g_amb,pars%g_apb,0d0,g_tmp)

        ! Gather neutral excitations vectors
        call nga_distribution(g_diag,pars%me,ilo,ihi)
        if(ihi.ge.ilo) then
          call nga_put(g_diag,ilo,ihi,1d0/dsqrt(omega(ilo:ihi)),
     &                 ihi-ilo+1)
        endif
        call ga_sync()
        call ga_scale_cols(g_tmp,g_diag)

        ! Destroy unused GAs
        call ga_copy(g_tmp,pars%g_apb)
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 
        if (.not.ga_destroy(g_amb))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 
        if (.not.ga_destroy(g_diag))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 
        
      else
        nocc2(1) = pars%nocc(1)**2
        nvir2(1) = pars%nvir(1)**2
        nocc2(2) = pars%nocc(2)**2
        nvir2(2) = pars%nvir(2)**2

        call ga_zero(pars%g_apb)

        ! Gather eigenvalue differences
        call ga_distribution(pars%g_eriov(1),pars%me,ilo,ihi,jlo,jhi)
        call nga_put(g_diag,jlo,jhi,wia(1,1),jhi-jlo+1)
        call ga_distribution(pars%g_eriov(2),pars%me,ilo,ihi,jlo,jhi)
        call nga_put(g_diag,jlo+pars%npoles(1),jhi+pars%npoles(1),
     $               wia(1,2),jhi-jlo+1)

        do isp=1,pars%ipol
          off1 = (isp-1)*pars%npoles(1)

          ! Add Screened Coulomb contribution from A
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
              call ga_put(pars%g_apb,llia,ulia,jb,jb,
     &             dbl_mb(ktmp+(ij-jlo)*nvir2(isp)+
     &                    (bmo-1)*pars%nvir(isp)),pars%nvir(isp))
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

          call ga_copy_patch('n',pars%g_apb,1+off1,
     &           pars%npoles(isp)+off1,1+off1,pars%npoles(isp)+off1,
     &           g_amb,1+off1,pars%npoles(isp)+off1,1+off1,
     &           pars%npoles(isp)+off1)


          ! Add Screened Coulomb contribution from B
          if(.not.ga_create(mt_dbl,pars%npoles(isp),pars%npoles(isp),
     &                     'tmp',0,pars%npoles(isp),g_tmp))
     &        call errquit(pname//'could not create temp GA',234,GA_ERR)
          call ga_dgemm('t','n',pars%npoles(isp),pars%npoles(isp),
     &                   pars%nri,-1d0,pars%g_wov(isp),pars%g_wov(isp),
     &                   0d0,g_tmp)
          call ga_add_patch(1d0,g_tmp,1,pars%npoles(isp),
     &                                1,pars%npoles(isp),
     &                 1d0,pars%g_apb,1+off1,pars%npoles(isp)+off1,
     &                                1+off1,pars%npoles(isp)+off1,
     &                     pars%g_apb,1+off1,pars%npoles(isp)+off1,
     &                                1+off1,pars%npoles(isp)+off1)
          call ga_add_patch(-1d0,g_tmp,1,pars%npoles(isp),
     &                                 1,pars%npoles(isp),
     &                       1d0,g_amb,1+off1,pars%npoles(isp)+off1,
     &                                 1+off1,pars%npoles(isp)+off1,
     &                           g_amb,1+off1,pars%npoles(isp)+off1,
     &                                 1+off1,pars%npoles(isp)+off1)

          if (.not.ga_destroy(g_tmp))
     &      call errquit(pname//'could not destroy temp GA',1,GA_ERR) 
          if (.not.ga_destroy(pars%g_wov(isp)))
     &      call errquit(pname//'could not destroy wov GA',1,GA_ERR) 
        enddo

        ! Form RPA matrix
        if (singlet) then
          do isp=1,pars%ipol
            off1 = (isp-1)*pars%npoles(1)
            do jsp=1,pars%ipol
             off2 = (jsp-1)*pars%npoles(1)
             call ga_matmul_patch('t','n',2d0,1d0,
     $         pars%g_eriov(isp),1,pars%npoles(isp),1,nri,
     $         pars%g_eriov(jsp),1,nri,1,pars%npoles(jsp),
     $         pars%g_apb,1+off1,pars%npoles(isp)+off1,
     $                    1+off2,pars%npoles(jsp)+off2)
            enddo
          enddo
        endif
        call ga_add_diagonal(pars%g_apb,g_diag)
        call ga_add_diagonal(g_amb,g_diag)


        ! Destroy ERIs
        if (.not.ga_destroy(pars%g_eriov(1)))
     &    call errquit(pname//'could not destroy eriov GA',275,GA_ERR) 
        if (.not.ga_destroy(pars%g_eriov(2)))
     &    call errquit(pname//'could not destroy eriov GA',277,GA_ERR) 


        ! First step is to get (A-B)^(1/2)
        call ga_sync()
#ifdef SCALAPACK
        info= ga_cholesky('L',g_amb)
#elif defined(PEIGS)
        call ga_chol(g_amb, g_amb, info)
#else 
        call ga_chol_seq(g_amb, g_amb, info)
#endif
        if (info.ne.0) call errquit("Cholesky Failed",298,GA_ERR)

        ! Transform (A+B) matrix with Cholesky factors
        if(.not.ga_create(mt_dbl,npoles,npoles,'tmp',0,npoles,
     &     g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('n','n',npoles,npoles,npoles,1d0,
     &                 pars%g_apb,g_amb,0d0,g_tmp)
        call ga_dgemm('t','n',npoles,npoles,npoles,1d0,
     &                 g_amb,g_tmp,0d0,pars%g_apb)   
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

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

        !
        ! Since we transformed the eigenvalue equation to an Hermitian
        ! form, we obtain the square of the neutral excitations
        ! (Equation (32)).
        !
!$omp   parallel do
        do ipole=1,npoles
          omega(ipole) = dsqrt(omega(ipole))
        enddo
!$omp   end parallel do

        ! Transform eigenvectors with (A-B)^(1/2)
        if(.not.ga_create(mt_dbl,npoles,npoles,'tmp',0,npoles,
     &     g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('n','n',npoles,npoles,npoles,1d0,
     &                 g_amb,pars%g_apb,0d0,g_tmp)

        ! Gather neutral excitations vectors
        call nga_distribution(g_diag,pars%me,ilo,ihi)
        if(ihi.ge.ilo) then
          call nga_put(g_diag,ilo,ihi,1d0/dsqrt(omega(ilo:ihi)),
     &                 ihi-ilo+1)
        endif
        call ga_sync()
        call ga_scale_cols(g_tmp,g_diag)
        call ga_copy(g_tmp,pars%g_apb)

        ! Destroy unused GAs
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',344,GA_ERR) 
        if (.not.ga_destroy(g_amb))
     &    call errquit(pname//'could not destroy temp GA',346,GA_ERR) 
        if (.not.ga_destroy(g_diag))
     &    call errquit(pname//'could not destroy temp GA',348,GA_ERR) 

      endif

      end subroutine
