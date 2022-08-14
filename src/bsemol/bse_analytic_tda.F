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

      integer ilo,ihi,jlo,jhi
      integer nocc2(3),nvir2(3)
      integer ipole,imo,jmo,amo,bmo
      integer k,l,kproc,lproc,klocal,llocal,kglobal,lglobal
      integer g_diag, g_tmp, g_amb, ktmp, ld, info, i, j
      integer modtn,tovern
      double precision alpha

#ifdef USE_OPENMP
      call util_blas_set_num_threads(pars%iMaxthreads)
#endif

      if(.not.ga_create(mt_dbl,totpoles,totpoles,'apb',0,totpoles,
     &   pars%g_apb))
     &  call errquit(pname//'could not create A+B GA',0,GA_ERR) 
      if (.not.nga_create(mt_dbl,1,totpoles,'wia',0,g_diag))
     $  call errquit(pname//'could not create wia',0,GA_ERR)

      if (pars%ipol.eq.1) then
        nocc2(1) = pars%nocc(1)**2
        nvir2(1) = pars%nvir(1)**2

        call ga_distribution(pars%g_apb,pars%me,ilo,ihi,jlo,jhi)
        if (jhi.ge.jlo) call nga_put(g_diag,jlo,jhi,wia,jhi-jlo+1)

        ! Screened Coulomb contribution
        if(.not.ga_create(mt_dbl,nocc2(1),nvir2(1),'tmp',nocc2(1),0,
     &                    g_tmp))
     &    call errquit(pname//'could not create temp GA',1,GA_ERR) 
        call ga_dgemm('t','n',nocc2(1),nvir2(1),pars%nri,-1d0,
     %                 pars%g_erioo(1),pars%g_erivv(1),0d0,g_tmp)
        call ga_distribution(g_tmp,pars%me,ilo,ihi,jlo,jhi)
        call ga_access(g_tmp,ilo,ihi,jlo,jhi,ktmp,ld)

        ! Need to take care about some crazy indexing
        tovern = totpoles/pars%nprocs
        modtn  = mod(totpoles,pars%nprocs)
        do i=1,nocc2(1)
          jmo = (i-1)/pars%nocc(1) + 1
          imo = i - (jmo-1)*pars%nocc(1)
          do j=jlo,jhi
            bmo = (j-1)/pars%nvir(1) + 1
            amo = j - (bmo-1)*pars%nvir(1)
            k = amo + (imo-1)*pars%nvir(1)
            l = bmo + (jmo-1)*pars%nvir(1)
            kproc = mod(k-1,pars%nprocs)
            lproc = mod(l-1,pars%nprocs)
            klocal = (k-1)/pars%nprocs + 1
            llocal = (l-1)/pars%nprocs + 1
            kglobal = klocal + kproc*tovern + min(kproc, modtn)
            lglobal = llocal + lproc*tovern + min(lproc, modtn)
            call ga_put(pars%g_apb,kglobal,kglobal,lglobal,lglobal,
     &                  dbl_mb(ktmp+(j-jlo)*nocc2(1)+i-1),1)
          enddo
        enddo
        call ga_release(g_tmp,ilo,ihi,jlo,jhi)
        if (.not.ga_destroy(g_tmp))
     &    call errquit(pname//'could not destroy temp GA',1,GA_ERR) 

        ! Hartree contribution
        if (singlet) then
          call ga_dgemm('t','n',totpoles,totpoles,pars%nri,2d0,
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

      endif

      end subroutine
