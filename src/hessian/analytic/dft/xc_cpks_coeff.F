c
c     Precomputes coefficients of basis functions pairs for use in
c     making multiple XC matrices for CPKS
c
c     BGJ - 8/98
c
c     $Id: xc_cpks_coeff.F,v 1.7 2002-03-22 02:29:47 sohirata Exp $
c
c     2/02 So Hirata --- added a code for triplet excitation or
c     instability calculation and the corresponding logical argu-
c     ment TRIPLET, which must be false for hessian calculation.
c
      Subroutine xc_cpks_coeff(prho, pdelrho, Amat2, Cmat2, Cmat,
     &                         delrho, npert, ipol, nq, GRAD, TRIPLET)
      implicit none
#include "rtdb.fh"
#include "dft2drv.fh"
c !!! BGJ test
#include "bgj.fh"
c !!! BGJ test
c
      integer npert    ! Number CPKS perturbations [input]
      integer ipol     ! [input]
      integer nq       ! [input]
      logical GRAD     ! Whether gradient-corrected [input]
      logical TRIPLET  ! Whether triplet [input]
c
c     Current approximate perturbed spin densities and density gradients
c     These are overwritten with the XC matrix coefficients to save space
c
      double precision prho(nq,ipol,npert), pdelrho(nq,3,ipol,npert)
c
c     Second derivatives of XC functional [input]
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
c
c     GC first derivatives of XC functional [input]
c
      double precision Cmat(nq,*)
c
c     Gradients of spin densities [input]
c
      double precision delrho(nq,3,ipol)
c
      integer ipert, n ! Loop indices
      double precision pz(5), pdra(3), pdrb(3), t(2)
c
      if (ipol.eq.1) then
c
c     Since the total densities are evaluated in the restricted case,
c     scale them by a factor of 0.5 so that the correct CPKS matrices
c     will be produced.
c
         call dscal(nq*ipol*npert,0.5d0,prho,1)
         if (GRAD) then
            call dscal(nq*3*ipol*npert,0.5d0,pdelrho,1)
            call dscal(nq*3*ipol,0.5d0,delrho,1)
         endif
      endif
c
      do ipert = 1, npert
c
c        !!! Put in cutoffs here similar to xc_tabcd? !!!
c
         if (ipol.eq.2) then
c
c           Unrestricted case
c           (DFT analytic hessian)
c           (TDDFT excitation energies)
c           (DFT instability check ... NYI)
c
            if (.not. GRAD) then
c
               do n = 1, nq
                  pz(1) = prho(n,1,ipert)
                  pz(2) = prho(n,2,ipert)
                  prho(n,1,ipert) = Amat2(n,D2_RA_RA)*pz(1)
     &                            + Amat2(n,D2_RA_RB)*pz(2)
                  prho(n,2,ipert) = Amat2(n,D2_RA_RB)*pz(1)
     &                            + Amat2(n,D2_RB_RB)*pz(2)
               enddo
c
            else
c
               do n = 1, nq
c
c                 Perturbed functional parameters
c
                  pz(1) = prho(n,1,ipert)
                  pz(2) = prho(n,2,ipert)
                  pz(3) =(delrho(n,1,1)*pdelrho(n,1,1,ipert)
     &                  + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     &                  + delrho(n,3,1)*pdelrho(n,3,1,ipert))*2d0
                  pz(4) = delrho(n,1,1)*pdelrho(n,1,2,ipert)
     &                  + delrho(n,2,1)*pdelrho(n,2,2,ipert)
     &                  + delrho(n,3,1)*pdelrho(n,3,2,ipert)
     &                  + delrho(n,1,2)*pdelrho(n,1,1,ipert)
     &                  + delrho(n,2,2)*pdelrho(n,2,1,ipert)
     &                  + delrho(n,3,2)*pdelrho(n,3,1,ipert)
                  pz(5) =(delrho(n,1,2)*pdelrho(n,1,2,ipert)
     &                  + delrho(n,2,2)*pdelrho(n,2,2,ipert)
     &                  + delrho(n,3,2)*pdelrho(n,3,2,ipert))*2d0
c
c                 Local CPKS coefficients
c
                  prho(n,1,ipert) = Amat2(n,D2_RA_RA) *pz(1)
     &                            + Amat2(n,D2_RA_RB) *pz(2)
     &                            + Cmat2(n,D2_RA_GAA)*pz(3)
     &                            + Cmat2(n,D2_RA_GAB)*pz(4)
     &                            + Cmat2(n,D2_RA_GBB)*pz(5)
                  prho(n,2,ipert) = Amat2(n,D2_RA_RB) *pz(1)
     &                            + Amat2(n,D2_RB_RB) *pz(2)
     &                            + Cmat2(n,D2_RB_GAA)*pz(3)
     &                            + Cmat2(n,D2_RB_GAB)*pz(4)
     &                            + Cmat2(n,D2_RB_GBB)*pz(5)
c
c                 GC CPKS coefficients
c
                  pdra(1) = pdelrho(n,1,1,ipert)
                  pdra(2) = pdelrho(n,2,1,ipert)
                  pdra(3) = pdelrho(n,3,1,ipert)
                  pdrb(1) = pdelrho(n,1,2,ipert)
                  pdrb(2) = pdelrho(n,2,2,ipert)
                  pdrb(3) = pdelrho(n,3,2,ipert)
c
c                 Perturbed density gradient contribution
c
                  pdelrho(n,1,1,ipert) = Cmat(n,D1_GAA)*pdra(1)*2d0
     &                                 + Cmat(n,D1_GAB)*pdrb(1)
                  pdelrho(n,2,1,ipert) = Cmat(n,D1_GAA)*pdra(2)*2d0
     &                                 + Cmat(n,D1_GAB)*pdrb(2)
                  pdelrho(n,3,1,ipert) = Cmat(n,D1_GAA)*pdra(3)*2d0
     &                                 + Cmat(n,D1_GAB)*pdrb(3)
                  pdelrho(n,1,2,ipert) = Cmat(n,D1_GBB)*pdrb(1)*2d0
     &                                 + Cmat(n,D1_GAB)*pdra(1)
                  pdelrho(n,2,2,ipert) = Cmat(n,D1_GBB)*pdrb(2)*2d0
     &                                 + Cmat(n,D1_GAB)*pdra(2)
                  pdelrho(n,3,2,ipert) = Cmat(n,D1_GBB)*pdrb(3)*2d0
     &                                 + Cmat(n,D1_GAB)*pdra(3)
c
c                 Perturbed functional argument contribution
c
                  t(1) =(Cmat2(n,D2_RA_GAA) *pz(1)
     &                 + Cmat2(n,D2_RB_GAA) *pz(2)
     &                 + Cmat2(n,D2_GAA_GAA)*pz(3)
     &                 + Cmat2(n,D2_GAA_GAB)*pz(4)
     &                 + Cmat2(n,D2_GAA_GBB)*pz(5))*2d0
                  t(2) = Cmat2(n,D2_RA_GAB) *pz(1)
     &                 + Cmat2(n,D2_RB_GAB) *pz(2)
     &                 + Cmat2(n,D2_GAA_GAB)*pz(3)
     &                 + Cmat2(n,D2_GAB_GAB)*pz(4)
     &                 + Cmat2(n,D2_GAB_GBB)*pz(5)
                  pdelrho(n,1,1,ipert) = pdelrho(n,1,1,ipert)
     &                 + delrho(n,1,1)*t(1) + delrho(n,1,2)*t(2)
                  pdelrho(n,2,1,ipert) = pdelrho(n,2,1,ipert)
     &                 + delrho(n,2,1)*t(1) + delrho(n,2,2)*t(2)
                  pdelrho(n,3,1,ipert) = pdelrho(n,3,1,ipert)
     &                 + delrho(n,3,1)*t(1) + delrho(n,3,2)*t(2)
                  t(1) = t(2)
                  t(2) =(Cmat2(n,D2_RA_GBB) *pz(1)
     &                 + Cmat2(n,D2_RB_GBB) *pz(2)
     &                 + Cmat2(n,D2_GAA_GBB)*pz(3)
     &                 + Cmat2(n,D2_GAB_GBB)*pz(4)
     &                 + Cmat2(n,D2_GBB_GBB)*pz(5))*2d0
                  pdelrho(n,1,2,ipert) = pdelrho(n,1,2,ipert)
     &                 + delrho(n,1,1)*t(1) + delrho(n,1,2)*t(2)
                  pdelrho(n,2,2,ipert) = pdelrho(n,2,2,ipert)
     &                 + delrho(n,2,1)*t(1) + delrho(n,2,2)*t(2)
                  pdelrho(n,3,2,ipert) = pdelrho(n,3,2,ipert)
     &                 + delrho(n,3,1)*t(1) + delrho(n,3,2)*t(2)
               enddo
c
            endif
c
         else if (TRIPLET) then
c
c           Restricted triplet case
c           (TDDFT triplet excitation energies)
c           (DFT triplet instability check ... NYI)
c
            if (.not. GRAD) then
c
               do n = 1, nq
                  prho(n,1,ipert) =(Amat2(n,D2_RA_RA)-Amat2(n,D2_RA_RB))
     &                            * prho(n,1,ipert)
               enddo
c
            else
c
               do n = 1, nq
c
c                 Perturbed functional parameters
c
                  pz(1) = prho(n,1,ipert)
                  pz(3) =(delrho(n,1,1)*pdelrho(n,1,1,ipert)
     &                  + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     &                  + delrho(n,3,1)*pdelrho(n,3,1,ipert))*2d0
c
c                 Local CPKS coefficient
c
                  prho(n,1,ipert) =(Amat2(n,D2_RA_RA)
     &                            - Amat2(n,D2_RA_RB)) *pz(1)
     &                            +(Cmat2(n,D2_RA_GAA)
     &                            - Cmat2(n,D2_RA_GBB))*pz(3)
c
c                 GC CPKS coefficients
c
                  pdra(1) = pdelrho(n,1,1,ipert)
                  pdra(2) = pdelrho(n,2,1,ipert)
                  pdra(3) = pdelrho(n,3,1,ipert)
c
c                 Perturbed density gradient contribution
c
                  pdelrho(n,1,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 - Cmat(n,D1_GAB))*pdra(1)
                  pdelrho(n,2,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 - Cmat(n,D1_GAB))*pdra(2)
                  pdelrho(n,3,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 - Cmat(n,D1_GAB))*pdra(3)
c
c                 Perturbed functional argument contribution
c
                  t(1) =((Cmat2(n,D2_RA_GAA)
     &                 -  Cmat2(n,D2_RA_GBB)) *pz(1)
     &                 + (Cmat2(n,D2_GAA_GAA)
     &                 -  Cmat2(n,D2_GAA_GBB))*pz(3))*2d0
                  pdelrho(n,1,1,ipert) = pdelrho(n,1,1,ipert)
     &                 + delrho(n,1,1)*t(1)
                  pdelrho(n,2,1,ipert) = pdelrho(n,2,1,ipert)
     &                 + delrho(n,2,1)*t(1)
                  pdelrho(n,3,1,ipert) = pdelrho(n,3,1,ipert)
     &                 + delrho(n,3,1)*t(1)
               enddo
            endif
c
         else
c
c           Restricted case
c           (DFT analytic hessian)
c           (TDDFT singlet excitation energies)
c           (DFT singlet instability check ... NYI)
c
            if (.not. GRAD) then
c
               do n = 1, nq
                  prho(n,1,ipert) =(Amat2(n,D2_RA_RA)+Amat2(n,D2_RA_RB))
     &                            * prho(n,1,ipert)
               enddo
c
            else
c
               do n = 1, nq
c
c                 Perturbed functional parameters
c
                  pz(1) = prho(n,1,ipert)
                  pz(3) =(delrho(n,1,1)*pdelrho(n,1,1,ipert)
     &                  + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     &                  + delrho(n,3,1)*pdelrho(n,3,1,ipert))*2d0
c
c                 Local CPKS coefficient
c
                  prho(n,1,ipert) =(Amat2(n,D2_RA_RA)
     &                            + Amat2(n,D2_RA_RB)) *pz(1)
     &                            +(Cmat2(n,D2_RA_GAA)
     &                            + Cmat2(n,D2_RA_GAB)
     &                            + Cmat2(n,D2_RA_GBB))*pz(3)
c
c                 GC CPKS coefficients
c
                  pdra(1) = pdelrho(n,1,1,ipert)
                  pdra(2) = pdelrho(n,2,1,ipert)
                  pdra(3) = pdelrho(n,3,1,ipert)
c
c                 Perturbed density gradient contribution
c
                  pdelrho(n,1,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 + Cmat(n,D1_GAB))*pdra(1)
                  pdelrho(n,2,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 + Cmat(n,D1_GAB))*pdra(2)
                  pdelrho(n,3,1,ipert) =(Cmat(n,D1_GAA)*2d0
     &                                 + Cmat(n,D1_GAB))*pdra(3)
c
c                 Perturbed functional argument contribution
c
                  t(1) =((Cmat2(n,D2_RA_GAA)
     &                 +  Cmat2(n,D2_RA_GBB)) *pz(1)
     &                 + (Cmat2(n,D2_GAA_GAA)
     &                 +  Cmat2(n,D2_GAA_GAB)
     &                 +  Cmat2(n,D2_GAA_GBB))*pz(3))*2d0
                  t(2) = Cmat2(n,D2_RA_GAB)*2d0*pz(1)
     &                 +(Cmat2(n,D2_GAA_GAB)*2d0
     &                 + Cmat2(n,D2_GAB_GAB))*pz(3)
                  pdelrho(n,1,1,ipert) = pdelrho(n,1,1,ipert)
     &                 + delrho(n,1,1)*(t(1) + t(2))
                  pdelrho(n,2,1,ipert) = pdelrho(n,2,1,ipert)
     &                 + delrho(n,2,1)*(t(1) + t(2))
                  pdelrho(n,3,1,ipert) = pdelrho(n,3,1,ipert)
     &                 + delrho(n,3,1)*(t(1) + t(2))
               enddo
            endif
c
         endif
c
      enddo
c
      if (ipol.eq.1 .and. GRAD) then
c
c     Put delrho back the way it was since it may be used later on
c
         call dscal(nq*3*ipol,2d0,delrho,1)
      endif
c
#if 0
      if (bgj_print() .gt. 1) then
      write(6,*) 'xc_cpks_coeff: prho out'
      call output(prho, 1, nq, 1, ipol*npert, nq, ipol*npert, 1)
      if (GRAD) then
         write(6,*) 'xc_cpks_coeff: pdelrho out'
         call output(pdelrho, 1, nq, 1, 3*ipol*npert, nq, 3*ipol*npert,
     &               1)
      endif
      endif
#endif
      return
      end
