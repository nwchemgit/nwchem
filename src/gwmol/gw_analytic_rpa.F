      subroutine gw_analytic_rpa(pars,eia,omega,maxpoles,maxnpoles,nri)
C
C     Build and diagonalize the RPA polarizability Casida matrix.
C     Transform the ERIs with the eigenvectors obtained.
C
C     We convert the full Casida matrix (Equation (12) and (13)), into
C     the Hermitian matrix `A+B` appearing on Equations (32) and (33).
C     The matrix `A-B` and, in principle, should be positive definite.
C
C     See Equations (12), (13), (34), (36) of
C     JCTC 17, 7504 (2021).
C
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"      
#include "gw.fh"
#ifdef SCALAPACK
#include "scaleig.fh"
#endif
      type(gw_params_t)  :: pars
      integer,intent(in) :: maxpoles, maxnpoles, nri

      double precision, intent(in)  :: eia(maxpoles,*)
      double precision, intent(out) :: omega(2*maxnpoles)

      integer g_rpa, g_xpy, g_diag, g_diag1, g_diag2
      integer npoles,mynpoles,ipole,ilo,ihi,jlo,jhi
      integer kRPA,kOV1,kOV2,kXPY,ld
      integer iMaxThreads
      integer lo,hi
      logical ok

      character(len=17),parameter :: pname = 'gw_analytic_rpa: '

#ifdef USE_OPENMP
      integer,external :: omp_get_max_threads
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif      

      !
      ! Closed-shell case
      !
      if (pars%ipol.eq.1) then

        npoles = pars%npoles(1)
        mynpoles = pars%mynpoles(1)

        if (.not.ga_create(mt_dbl,npoles,npoles,'rpa',0,npoles,g_rpa))
     $    call errquit(pname//'can''t create RPA GA',0,GA_ERR)

        call ga_distribution(g_rpa,pars%me,ilo,ihi,jlo,jhi)
        if (pars%ovlo(1).ne.ilo .or. pars%ovhi(1).ne.ihi)
     $    call errquit('gw_rpa: wrong distribution',0,GA_ERR)

        call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $                 kOV1,ld)                

        !
        ! In order to obtain the matrix (A-B)^{1/2}(A+B)(A-B)^{1/2} is
        ! convenient to scale the ERIs with (A-B)^{1/2} and then use
        ! DYSRK.
        !
        call gw_analytic_ovscale(1,dbl_mb(kOV1),eia,nri,mynpoles)
        call ga_release_update(pars%g_eriov(1),1,nri,pars%ovlo(1),
     $                         pars%ovhi(1))   

        !
        ! We use GA_DGEMM instead of the missing GA_DSYRK
        !
        call ga_dgemm('t','n',pars%npoles(1),pars%npoles(1),nri,4d0,
     $                 pars%g_eriov(1),pars%g_eriov(1),0d0,g_rpa)

        !
        ! We still need to add the square of the eigenvalue differences
        ! (Equation (34) scaled by (A-B)^{1/2} on both sides)
        !
        call ga_access(g_rpa,pars%ovlo(1),pars%ovhi(1),1,npoles,
     $                 kRPA,ld)   
        call gw_analytic_addrpa(dbl_mb(kRPA),eia,mynpoles,pars%ovlo(1))
        call ga_release_update(g_rpa,ilo,ihi,jlo,jhi)


        !
        ! Obtain the RPA correlation energy from the diagonal of the
        ! Casida matrix.
        !
        call ga_access(g_rpa,pars%ovlo(1),pars%ovhi(1),1,npoles,kRPA,ld)
        call gw_analytic_ecrpa(dbl_mb(kRPA),eia,mynpoles,pars%ovlo(1),
     $                         pars%ecrpa)
        call ga_release(g_rpa,ilo,ihi,jlo,jhi)

        !
        ! Perform the Eigendecomposition of the Casida matrix
        !
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK        
        call dft_scaleig(g_rpa,g_rpa,omega,0)
#else
        call ga_diag_std(g_rpa,g_rpa,omega)
#endif
#else
        call ga_diag_std_seq(g_rpa,g_rpa,omega)
#endif

        !
        ! Since we transform the eigenvalue equation to an Hermitian
        ! form, we obtain the square of the neutral excitations
        ! (Equation (32)).
        !
!$omp   parallel do
        do ipole=1,npoles
          omega(ipole) = dsqrt(omega(ipole))
        enddo
!$omp   end parallel do

        pars%ecrpa = pars%ecrpa + 0.5d0*sum(omega(1:npoles))


        !
        ! Include 1/\Omega^{1/2} factor in the eigenvectors.
        !
        ! The manuscript has the 1/\Omega^{1/2} missing!
        ! (see Equations (42a) and (42b) from Comp. Phys. Comm. 208, 149
        ! (2016)).
        !
        if (.not.ga_duplicate(pars%g_eriov(1),pars%g_erim,
     $                        'transformed eri alpha'))
     $    call errquit(pname//'failed to duplicate array',0,GA_ERR)

        call ga_access(g_rpa,pars%ovlo(1),pars%ovhi(1),1,npoles,
     $                 kRPA,ld)   
        call gw_analytic_scalxpy(dbl_mb(kRPA),omega,mynpoles,npoles)
        call ga_release_update(g_rpa,pars%ovlo(1),pars%ovhi(1),1,
     $                         npoles)   


        !
        ! We do not need to multiply Equation (33) by (A-B)^{1/2} since
        ! this factor is already present in the ERIs. Therefore, we just
        ! transform the scaled ERIs by the scaled T, to obtain the
        ! first summation of Equation (36).
        !
        call ga_dgemm('n','n',nri,npoles,npoles,1d0,pars%g_eriov(1),
     $                 g_rpa,0d0,pars%g_erim)


        !
        ! We will need the original ERIs, so we remove the (A-B)^{1/2}
        ! factor.
        !
        call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $                 kOV1,ld)
        call gw_analytic_ovscale(2,dbl_mb(kOV1),eia,nri,mynpoles)
        call ga_release_update(pars%g_eriov(1),1,nri,pars%ovlo(1),
     $                         pars%ovhi(1))   

        if (.not.(ga_destroy(g_rpa)))
     $    call errquit(pname//'can''t destroy RPA GA',0,GA_ERR)


      !
      ! Open-shell case
      !
      else

        npoles = pars%npoles(1) + pars%npoles(2)
        omega(1:npoles) = 0d0

        if (.not.ga_create(mt_dbl,npoles,npoles,'rpa',0,npoles,g_rpa))
     $    call errquit(pname//'can''t create RPA GA',0,GA_ERR)
        if (.not.nga_create(mt_dbl,1,pars%npoles(1),'RPA diag 1',0,
     $       g_diag1))
     $    call errquit(pname//'can''t create RPA diag 1',0,GA_ERR)
        if (.not.nga_create(mt_dbl,1,pars%npoles(2),'RPA diag 2',0,
     $       g_diag2))
     $    call errquit(pname//'can''t create RPA diag 2',0,GA_ERR)
        if (.not.nga_create(mt_dbl,1,npoles,'RPA diag',0,g_diag))
     $    call errquit(pname//'can''t create RPA diag',0,GA_ERR)

c       Gather eigenvalue differences in a GA vector
        call nga_distribution(g_diag1,pars%me,lo,hi)
        if (pars%ovlo(1).ne.lo .or. pars%ovhi(1).ne.hi)
     $    call errquit('gw_rpa: wrong distribution diag1',0,GA_ERR)
        if (hi.ge.lo) then
          call nga_put(g_diag1,lo,hi,eia(1,1),ld)
        endif

        call nga_distribution(g_diag2,pars%me,lo,hi)
        if (pars%ovlo(2).ne.lo .or. pars%ovhi(2).ne.hi)
     $    call errquit('gw_rpa: wrong distribution diag2',0,GA_ERR)
        if (hi.ge.lo) then
          call nga_put(g_diag2,lo,hi,eia(1,2),ld)
        endif

        call ga_sync()
        call nga_copy_patch('n',g_diag1,1,pars%npoles(1),
     $                      g_diag,1,pars%npoles(1))
        call nga_copy_patch('n',g_diag2,1,pars%npoles(2),
     $                      g_diag,pars%npoles(1)+1,npoles)   

c       Square eigenvalue differences        
        call ga_elem_multiply(g_diag,g_diag,g_diag)

c       Scale integrals with eigenvalue differences
        call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $                 kOV1,ld)   
        call gw_analytic_ovscale(1,dbl_mb(kOV1),eia,nri,
     $                           pars%mynpoles(1))
        call ga_release_update(pars%g_eriov(1),1,nri,pars%ovlo(1),
     $                         pars%ovhi(1))   


        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV2,ld)   
        call gw_analytic_ovscale(1,dbl_mb(kOV2),eia(1,2),nri,
     $                           pars%mynpoles(2))
        call ga_release_update(pars%g_eriov(2),1,nri,pars%ovlo(2),
     $                         pars%ovhi(2))   

c       Form open-shell RPA matrix        
        call ga_zero(g_rpa)
        call ga_matmul_patch('t','n',2d0,0d0,
     $                        pars%g_eriov(1),1,pars%npoles(1),1,nri,
     $                        pars%g_eriov(1),1,nri,1,pars%npoles(1),
     $                        g_rpa,1,pars%npoles(1),1,pars%npoles(1))
        call ga_matmul_patch('t','n',2d0,0d0,
     $                        pars%g_eriov(1),1,pars%npoles(1),1,nri,
     $                        pars%g_eriov(2),1,nri,1,pars%npoles(2),
     $                        g_rpa,1,pars%npoles(1),
     $                              pars%npoles(1)+1,npoles)
        call ga_matmul_patch('t','n',2d0,0d0,
     $                        pars%g_eriov(2),1,pars%npoles(2),1,nri,
     $                        pars%g_eriov(1),1,nri,1,pars%npoles(1),
     $                        g_rpa,pars%npoles(1)+1,npoles,
     $                              1,pars%npoles(1))
        call ga_matmul_patch('t','n',2d0,0d0,
     $                        pars%g_eriov(2),1,pars%npoles(2),1,nri,
     $                        pars%g_eriov(2),1,nri,1,pars%npoles(2),
     $                        g_rpa,pars%npoles(1)+1,npoles,
     $                              pars%npoles(1)+1,npoles)
        call ga_add_diagonal(g_rpa,g_diag)
        call ga_sync()

c       ECRPA pending
c       pars%ecrpa =      

c       Diagonalize RPA matrix
#ifdef PARALLEL_DIAG
#ifdef SCALAPACK        
        call dft_scaleig(g_rpa,g_rpa,omega,0)
#else
        call ga_diag_std(g_rpa,g_rpa,omega)
#endif
#else
        call ga_diag_std_seq(g_rpa,g_rpa,omega)
#endif

        ok = ma_verify_allocator_stuff()

!$omp   parallel do
        do ipole=1,npoles
          omega(ipole) = dsqrt(omega(ipole))
        enddo
!$omp   end parallel do

c       Create transformed ERI array        
        if (.not.ga_create(mt_dbl,nri,npoles,'transformed eris',nri,0,
     $       pars%g_erim))
     $    call errquit(pname//'failed to create array',0,GA_ERR)
          
c       Gather neutral excitations vector
        call nga_distribution(g_diag,pars%me,lo,hi)
        if (hi.ge.lo) then
          call nga_put(g_diag,lo,hi,1d0/dsqrt(omega(lo:hi)),hi-lo+1)
        endif
        call ga_sync()
        call ga_scale_cols(g_rpa,g_diag)

c       Transform ERIs with excitation vectors
        call ga_matmul_patch('n','n',1d0,0d0,
     $          pars%g_eriov(1),1,nri,1,pars%npoles(1),
     $          g_rpa,1,pars%npoles(1),1,npoles,
     $          pars%g_erim,1,nri,1,npoles)
        call ga_matmul_patch('n','n',1d0,1d0,
     $          pars%g_eriov(2),1,nri,1,pars%npoles(2),
     $          g_rpa,pars%npoles(1)+1,npoles,1,npoles,
     $          pars%g_erim,1,nri,1,npoles)   

c       Scale back ERIs        
        call ga_access(pars%g_eriov(1),1,nri,pars%ovlo(1),pars%ovhi(1),
     $                 kOV1,ld)   
        call gw_analytic_ovscale(2,dbl_mb(kOV1),eia,nri,
     $                           pars%mynpoles(1))
        call ga_release_update(pars%g_eriov(1),1,nri,pars%ovlo(1),
     $                         pars%ovhi(1))   

        call ga_access(pars%g_eriov(2),1,nri,pars%ovlo(2),pars%ovhi(2),
     $                 kOV2,ld)   
        call gw_analytic_ovscale(2,dbl_mb(kOV2),eia(1,2),nri,
     $                           pars%mynpoles(2))
        call ga_release_update(pars%g_eriov(2),1,nri,pars%ovlo(2),
     $                         pars%ovhi(2))   

c       Destroy GAs        
        if (.not.(ga_destroy(g_rpa)))
     $    call errquit(pname//'can''t destroy RPA GA',0,GA_ERR)
        if (.not.(ga_destroy(g_diag)))
     $    call errquit(pname//'can''t destroy RPA diag',0,GA_ERR)
        if (.not.(ga_destroy(g_diag1)))
     $    call errquit(pname//'can''t destroy RPA diag',0,GA_ERR)
        if (.not.(ga_destroy(g_diag2)))
     $    call errquit(pname//'can''t destroy RPA diag',0,GA_ERR)

      endif

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif

      end subroutine
