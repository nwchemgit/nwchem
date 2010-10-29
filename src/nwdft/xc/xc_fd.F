c
c     Routines for calculating XC functional second derivatives
c     by finite difference
c
c     BGJ - 1/02
c
c     $Id$
c
c     Finite difference step size
c
#define STEP_SIZE 0.0001
c
      subroutine xc_setup_fd(tol_rho, rho, delrho, qwght, nq, ipol, GC,
     &     l_storage, i_prho, i_pdelrho, i_pAmat, i_pCmat, i_pfunc,
     &     i_qwght_copy)
c
      implicit none
#include "errquit.fh"
c
      integer nq, ipol, l_storage, i_prho, i_pdelrho, i_pAmat, i_pCmat,
     &     i_pfunc, i_qwght_copy
      double precision tol_rho, rho(nq,*), delrho(nq,3,*), qwght(nq)
      logical GC
c
c     Local variables
c
      integer npert, len_prho, len_pAmat, len_prho_tot, len_pAmat_tot
      integer i
c
#include "mafdecls.fh"
c
      npert = 2
      len_prho = 2*3*nq
      len_pAmat = npert*npert*nq
      if (GC) then
         npert = 5
         len_prho_tot = len_prho + 2*6*nq
         len_pAmat_tot = npert*npert*nq
      else
         len_prho_tot = len_prho
         len_pAmat_tot = len_pAmat
      endif
      if (.not.
     &     MA_Alloc_Get(MT_DBL,len_prho_tot+len_pAmat_tot+4*nq*npert,
     &     'xc_fd',l_storage,i_prho))
     &     call errquit('xc_setup_fd: cannot allocate storage',0,
     &       MA_ERR)
      call dfill(len_prho_tot+len_pAmat_tot+4*nq*npert, 0d0,
     &     dbl_mb(i_prho), 1)
      i_pAmat = i_prho + len_prho_tot
      i_pdelrho = i_prho + len_prho
      i_pCmat = i_pAmat + len_pAmat
      i_pfunc = i_pAmat + len_pAmat_tot
      i_qwght_copy = i_pfunc + 2*nq*npert
c
c     Copy weights
c
      do i = 1, 2*npert
         call dcopy(nq, qwght, 1, dbl_mb(i_qwght_copy+(i-1)*nq), 1)
      enddo
c
c     Set up perturbed density values
c
      call xc_pert_rho(tol_rho, rho, delrho, nq, ipol, GC,
     &     dbl_mb(i_prho), dbl_mb(i_pdelrho))
c
      return
      end
c
      subroutine xc_pert_rho(tol_rho, rho, delrho, nq, ipol, npert,
     &     prho, pdelrho)
c
      implicit none
c
      integer nq, ipol, npert
      double precision tol_rho, rho(nq,*), delrho(nq,3,*),
     &     prho(nq,3,2,npert), pdelrho(nq,3,2,2,npert)
c
c     Local variables
c
      double precision h
      parameter (h = STEP_SIZE)
      logical GC, DoA, DoB
      integer i, ipert
      double precision ra, drax, dray, draz, Gaa, Gab, Gbb,
     &     mod_drhoa, mod_drhob, cos_theta, x, mod_drhob_p,
     &     cos_theta_p, sin_theta_p
c
      GC = npert.eq.5
c
c     Initial copy of density values - copy as unrestricted
c
      call dfill(nq*3, 0d0, prho, 1)
      call dfill(nq*6, 0d0, pdelrho, 1)
      if (ipol.eq.1) then
         do i = 1, nq
            if (rho(i,1).gt.tol_rho) then
               ra = rho(i,1) * 0.5d0
               prho(i,1,1,1) = rho(i,1)
               prho(i,2,1,1) = ra
               prho(i,3,1,1) = ra
            endif
         enddo
         if (GC) then
            do i = 1, nq
               drax = delrho(i,1,1) * 0.5d0
               dray = delrho(i,2,1) * 0.5d0
               draz = delrho(i,3,1) * 0.5d0
               pdelrho(i,1,1,1,1) = drax
               pdelrho(i,2,1,1,1) = dray
               pdelrho(i,3,1,1,1) = draz
               pdelrho(i,1,2,1,1) = drax
               pdelrho(i,2,2,1,1) = dray
               pdelrho(i,3,2,1,1) = draz
            enddo
         endif
      else
         do i = 1, nq
            if (rho(i,1).gt.tol_rho) then
               prho(i,1,1,1) = rho(i,1)
               prho(i,2,1,1) = rho(i,2)
               prho(i,3,1,1) = rho(i,3)
            endif
         enddo
         if (GC) then
            call dcopy(nq*6, delrho, 1, pdelrho, 1)
         endif
      endif
c
c     Copy to all other perturbation locations
c
      call dcopy(nq*3, prho, 1, prho(1,1,2,1), 1)
      if (GC) then
         call dcopy(nq*6, pdelrho, 1, pdelrho(1,1,1,2,1), 1)
      endif
      do ipert = 2, npert
         do i = 1, 2
            call dcopy(nq*3, prho, 1, prho(1,1,i,ipert), 1)
            if (GC) then
               call dcopy(nq*6, pdelrho, 1, pdelrho(1,1,1,i,ipert), 1)
            endif
         enddo
      enddo
c
c     Perturb the density parameter values - we aren't concerned if
c     some density values go negative as a result of this since the
c     functional implementation routines map negative densities to
c     a zero value for the functional and its derivatives (e.g. like
c     analytic continuation)
c
      do ipert = 1, 2
         do i = 1, nq
            if (prho(i,1,1,ipert).gt.tol_rho) then
               prho(i,1,1,ipert) = prho(i,1,1,ipert) + h
               prho(i,ipert+1,1,ipert) = prho(i,ipert+1,1,ipert) + h
               prho(i,1,2,ipert) = prho(i,1,2,ipert) - h
               prho(i,ipert+1,2,ipert) = prho(i,ipert+1,2,ipert) - h
            endif
         enddo
      enddo
      if (GC) then
c
c     Perturb gamma values
c
         do i = 1, nq
            Gaa = pdelrho(i,1,1,1,3)*pdelrho(i,1,1,1,3)
     &          + pdelrho(i,2,1,1,3)*pdelrho(i,2,1,1,3)
     &          + pdelrho(i,3,1,1,3)*pdelrho(i,3,1,1,3)
            Gab = pdelrho(i,1,1,1,3)*pdelrho(i,1,2,1,3)
     &          + pdelrho(i,2,1,1,3)*pdelrho(i,2,2,1,3)
     &          + pdelrho(i,3,1,1,3)*pdelrho(i,3,2,1,3)
            Gbb = pdelrho(i,1,2,1,3)*pdelrho(i,1,2,1,3)
     &          + pdelrho(i,2,2,1,3)*pdelrho(i,2,2,1,3)
     &          + pdelrho(i,3,2,1,3)*pdelrho(i,3,2,1,3)
            mod_drhoa = sqrt(Gaa)
            mod_drhob = sqrt(Gbb)
c
c     Ensure perturbed gammas are always > 0, otherwise take derivative as 0
c
            DoA = Gaa.gt.h
            DoB = DoA.and.Gbb.gt.h
            if (DoB) then
               cos_theta = Gab/(mod_drhoa*mod_drhob)
            else
               cos_theta = 1d0
            endif
c
c     Now that we have the vital parameters clear out locations
c     in the perturbed density gradient array - we only need three
c     locations to construct the proper perturbed gamma values
c
            pdelrho(i,1,1,1,3) = 0d0
            pdelrho(i,2,1,1,3) = 0d0
            pdelrho(i,3,1,1,3) = 0d0
            pdelrho(i,1,2,1,3) = 0d0
            pdelrho(i,2,2,1,3) = 0d0
            pdelrho(i,3,2,1,3) = 0d0
            pdelrho(i,1,1,1,4) = 0d0
            pdelrho(i,2,1,1,4) = 0d0
            pdelrho(i,3,1,1,4) = 0d0
            pdelrho(i,1,2,1,4) = 0d0
            pdelrho(i,2,2,1,4) = 0d0
            pdelrho(i,3,2,1,4) = 0d0
            pdelrho(i,1,1,1,5) = 0d0
            pdelrho(i,2,1,1,5) = 0d0
            pdelrho(i,3,1,1,5) = 0d0
            pdelrho(i,1,2,1,5) = 0d0
            pdelrho(i,2,2,1,5) = 0d0
            pdelrho(i,3,2,1,5) = 0d0
            if (DoA) then
c
c     Perturb Gaa up
c
               x = mod_drhoa*cos_theta
               pdelrho(i,1,1,1,3) = x
               pdelrho(i,2,1,1,3) = sqrt(Gaa-x*x+h)
               pdelrho(i,1,2,1,3) = mod_drhob
            endif
            if (DoB) then
c
c     Perturb Gab up
c
               x = (Gab+h)/mod_drhob
               pdelrho(i,1,1,1,4) = x
               pdelrho(i,2,1,1,4) = sqrt(Gaa-x*x)
               pdelrho(i,1,2,1,4) = mod_drhob
c
c     Perturb Gbb up
c
               mod_drhob_p = sqrt(Gbb+h)
               cos_theta_p = (mod_drhob/mod_drhob_p)*cos_theta
               sin_theta_p = sqrt(1d0-cos_theta_p*cos_theta_p)
               pdelrho(i,1,1,1,5) = mod_drhoa*cos_theta_p
               pdelrho(i,2,1,1,5) = mod_drhoa*sin_theta_p
               pdelrho(i,1,2,1,5) = mod_drhob_p
            endif
            if (DoA) then
c
c     Perturb Gaa back
c
               x = mod_drhoa*cos_theta
               pdelrho(i,1,1,2,3) = x
               pdelrho(i,2,1,2,3) = sqrt(Gaa-x*x-h)
               pdelrho(i,1,2,2,3) = mod_drhob
            endif
            if (DoB) then
c
c     Perturb Gab back
c
               x = (Gab-h)/mod_drhob
               pdelrho(i,1,1,2,4) = x
               pdelrho(i,2,1,2,4) = sqrt(Gaa-x*x)
               pdelrho(i,1,2,2,4) = mod_drhob
c
c     Perturb Gbb back
c
               mod_drhob_p = sqrt(Gbb-h)
               cos_theta_p = (mod_drhob/mod_drhob_p)*cos_theta
               sin_theta_p = sqrt(1d0-cos_theta_p*cos_theta_p)
               pdelrho(i,1,1,2,5) = mod_drhoa*cos_theta_p
               pdelrho(i,2,1,2,5) = mod_drhoa*sin_theta_p
               pdelrho(i,1,2,2,5) = mod_drhob_p
            endif
         enddo
      endif
c
      return
      end
c
      subroutine xc_make_fd(Amat2, Cmat2, nq, GC, pAmat, pCmat)
c
      implicit none
c
      integer nq
      double precision Amat2(nq,*), Cmat2(nq,*),
     &     pAmat(nq,2,2,*), pCmat(nq,3,2,*)
      logical GC
c
c     Local variables
c
      double precision h, r2h
      parameter (h = STEP_SIZE)
      integer i, npert, ipert
c
#include "dft2drv.fh"
c
      if (GC) then
         npert = 5
      else
         npert = 2
      endif
      r2h = 0.5d0*h
c
c     Construct finite differences in the temporary arrays
c
      do ipert = 1, npert
         do i = 1, nq
            pAmat(i,1,1,ipert) = (pAmat(i,1,1,ipert)
     &                           -pAmat(i,1,2,ipert))*r2h
            pAmat(i,2,1,ipert) = (pAmat(i,2,1,ipert)
     &                           -pAmat(i,2,2,ipert))*r2h
         enddo
         if (GC) then
            do i = 1, nq
               pCmat(i,1,1,ipert) = (pCmat(i,1,1,ipert)
     &                              -pCmat(i,1,2,ipert))*r2h
               pCmat(i,2,1,ipert) = (pCmat(i,2,1,ipert)
     &                              -pCmat(i,2,2,ipert))*r2h
               pCmat(i,3,1,ipert) = (pCmat(i,3,1,ipert)
     &                              -pCmat(i,3,2,ipert))*r2h
            enddo
         endif
      enddo
c
c     Now scatter to Amat2 and Cmat2
c
      call dcopy(nq, pAmat(1,1,1,1), 1, Amat2(1,D2_RA_RA), 1)
      call dcopy(nq, pAmat(1,2,1,1), 1, Amat2(1,D2_RA_RB), 1)
c     Or could be as below, etc., etc.
c      call dcopy(nq, pAmat(1,1,1,2), 1, Amat2(1,D2_RA_RB), 1)
      call dcopy(nq, pAmat(1,2,1,2), 1, Amat2(1,D2_RB_RB), 1)
      if (GC) then
         call dcopy(nq, pCmat(1,D1_GAA,1,1), 1, Cmat2(1,D2_RA_GAA), 1)
c     Or could be as below, etc., etc.
c         call dcopy(nq, pAmat(1,1,1,3), 1, Cmat2(1,D2_RA_GAA), 1)
         call dcopy(nq, pCmat(1,D1_GAA,1,2), 1, Cmat2(1,D2_RB_GAA), 1)
         call dcopy(nq, pCmat(1,D1_GAA,1,3), 1, Cmat2(1,D2_GAA_GAA), 1)
c
         call dcopy(nq, pCmat(1,D1_GAB,1,1), 1, Cmat2(1,D2_RA_GAB), 1)
         call dcopy(nq, pCmat(1,D1_GAB,1,2), 1, Cmat2(1,D2_RB_GAB), 1)
         call dcopy(nq, pCmat(1,D1_GAB,1,3), 1, Cmat2(1,D2_GAA_GAB), 1)
         call dcopy(nq, pCmat(1,D1_GAB,1,4), 1, Cmat2(1,D2_GAB_GAB), 1)
c
         call dcopy(nq, pCmat(1,D1_GBB,1,1), 1, Cmat2(1,D2_RA_GBB), 1)
         call dcopy(nq, pCmat(1,D1_GBB,1,2), 1, Cmat2(1,D2_RB_GBB), 1)
         call dcopy(nq, pCmat(1,D1_GBB,1,3), 1, Cmat2(1,D2_GAA_GBB), 1)
         call dcopy(nq, pCmat(1,D1_GBB,1,4), 1, Cmat2(1,D2_GAB_GBB), 1)
         call dcopy(nq, pCmat(1,D1_GBB,1,5), 1, Cmat2(1,D2_GBB_GBB), 1)
      endif
c
      return
      end
