c
c     Precomputes coefficients of basis functions pairs for use in
c     making multiple XC matrices for CPKS
c
c     BGJ - 8/98
c
c     $Id$
c
c     2/02 So Hirata --- added a code for triplet excitation or
c     instability calculation and the corresponding logical argu-
c     ment TRIPLET, which must be false for hessian calculation.
c
      subroutine xc_cpks_coeff(prho, pdelrho, pttau, 
     &        Amat2, Cmat2, Cmat, Mmat2, 
     &        delrho, npert, ipol, nq, grad, kske, 
     &        triplet)
c
      implicit none
#include "rtdb.fh"
#include "dft2drv.fh"
#include "stdio.fh"
c !!! BGJ test
#include "bgj.fh"
c !!! BGJ test
c
      integer npert    ! Number CPKS perturbations [input]
      integer ipol     ! [input]
      integer nq       ! [input]
      logical grad     ! Whether gradient-corrected [input]
      logical triplet  ! Whether triplet [input]
c
c     Current approximate perturbed spin densities and density gradients
c     These are overwritten with the XC matrix coefficients to save space
c
      double precision prho(nq,ipol,npert)
      double precision pdelrho(nq,3,ipol,npert)
      double precision pttau(nq,ipol,npert)
c
c     Second derivatives of XC functional [input]
c
      double precision Amat2(nq,NCOL_AMAT2)
      double precision Cmat2(nq,NCOL_CMAT2)
      double precision Mmat2(nq,NCOL_MMAT2)
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
c
      logical kske
c
      double precision ptmp(5),t(2)
      double precision pdra(3), pdrb(3)
      double precision term_g,term_rr,term_rg,term_gg
      double precision term_rt,term_gt,term_tt
      double precision term_prho, term_pdelrho, term_pttau
      double precision gaa,gbb,gab
      double precision rara,rarb,rbra,rbrb
      double precision ragaa, ragab, ragbb, rbgaa, rbgab, rbgbb
      double precision gaagaa,gaagab,gaagbb,gabgab,gabgbb,gbbgbb
c
c     preliminaries
c
      do n = 1,5
         ptmp(n) = 0.d0
      end do
c
      if (ipol.eq.1) then
c
c     Since the total densities are evaluated in the restricted case,
c     scale them by a factor of 0.5 so that the correct CPKS matrices
c     will be produced.
c
         call dscal(nq*ipol*npert,0.5d0,prho,1)
         if (grad) then
            call dscal(nq*3*ipol*npert,0.5d0,pdelrho,1)
            call dscal(nq*3*ipol,0.5d0,delrho,1)
         endif
      endif
c
      do ipert = 1, npert
c
c        !!! Put in cutoffs here similar to xc_tabcd? !!!
c
         if (ipol.eq.2) then  ! unrestricted
c
c           Unrestricted case
c           (DFT analytic hessian)
c           (TDDFT excitation energies)
c           (DFT instability check ... NYI)
c
            if (.not. grad) then
               do n = 1, nq
                  ptmp(1) = prho(n,1,ipert)
                  ptmp(2) = prho(n,2,ipert)
                  prho(n,1,ipert) = Amat2(n,D2_RA_RA)*ptmp(1)
     &                            + Amat2(n,D2_RA_RB)*ptmp(2)
                  prho(n,2,ipert) = Amat2(n,D2_RA_RB)*ptmp(1)
     &                            + Amat2(n,D2_RB_RB)*ptmp(2)
               enddo
            else ! gradient dependent functionals
c
CDIR$ NEXTSCALAR
               do n = 1, nq
c
c                 Perturbed functional parameters
c
                  ptmp(1) = prho(n,1,ipert)
                  ptmp(2) = prho(n,2,ipert)
c
                  ptmp(3) =(delrho(n,1,1)*pdelrho(n,1,1,ipert) +
     &                      delrho(n,2,1)*pdelrho(n,2,1,ipert) +
     &                      delrho(n,3,1)*pdelrho(n,3,1,ipert))*2d0
c
                  ptmp(4) = delrho(n,1,1)*pdelrho(n,1,2,ipert) +
     &                      delrho(n,2,1)*pdelrho(n,2,2,ipert) +
     &                      delrho(n,3,1)*pdelrho(n,3,2,ipert) + 
     &                      delrho(n,1,2)*pdelrho(n,1,1,ipert) +
     &                      delrho(n,2,2)*pdelrho(n,2,1,ipert) +
     &                      delrho(n,3,2)*pdelrho(n,3,1,ipert)
c
                  ptmp(5) =(delrho(n,1,2)*pdelrho(n,1,2,ipert) +
     &                      delrho(n,2,2)*pdelrho(n,2,2,ipert) +
     &                      delrho(n,3,2)*pdelrho(n,3,2,ipert))*2d0
c
c                 Local CPKS coefficients
c
                  prho(n,1,ipert) = Amat2(n,D2_RA_RA) *ptmp(1) +
     &                              Amat2(n,D2_RA_RB) *ptmp(2) +
     &                              Cmat2(n,D2_RA_GAA)*ptmp(3) +
     &                              Cmat2(n,D2_RA_GAB)*ptmp(4) +
     &                              Cmat2(n,D2_RA_GBB)*ptmp(5)
c
                  prho(n,2,ipert) = Amat2(n,D2_RA_RB) *ptmp(1) +
     &                              Amat2(n,D2_RB_RB) *ptmp(2) +
     &                              Cmat2(n,D2_RB_GAA)*ptmp(3) +
     &                              Cmat2(n,D2_RB_GAB)*ptmp(4) +
     &                              Cmat2(n,D2_RB_GBB)*ptmp(5)
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
                  pdelrho(n,1,1,ipert) = Cmat(n,D1_GAA)*pdra(1)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdrb(1)
                  pdelrho(n,2,1,ipert) = Cmat(n,D1_GAA)*pdra(2)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdrb(2)
                  pdelrho(n,3,1,ipert) = Cmat(n,D1_GAA)*pdra(3)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdrb(3)
                  pdelrho(n,1,2,ipert) = Cmat(n,D1_GBB)*pdrb(1)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdra(1)
                  pdelrho(n,2,2,ipert) = Cmat(n,D1_GBB)*pdrb(2)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdra(2)
                  pdelrho(n,3,2,ipert) = Cmat(n,D1_GBB)*pdrb(3)*2d0 +
     &                                   Cmat(n,D1_GAB)*pdra(3)
c
c                 Perturbed functional argument contribution
c
                  t(1) =(Cmat2(n,D2_RA_GAA) *ptmp(1) +
     &                   Cmat2(n,D2_RB_GAA) *ptmp(2) +
     &                   Cmat2(n,D2_GAA_GAA)*ptmp(3) +
     &                   Cmat2(n,D2_GAA_GAB)*ptmp(4) +
     &                   Cmat2(n,D2_GAA_GBB)*ptmp(5))*2d0
c
                  t(2) = Cmat2(n,D2_RA_GAB) *ptmp(1) +
     &                   Cmat2(n,D2_RB_GAB) *ptmp(2) +
     &                   Cmat2(n,D2_GAA_GAB)*ptmp(3) +
     &                   Cmat2(n,D2_GAB_GAB)*ptmp(4) +
     &                   Cmat2(n,D2_GAB_GBB)*ptmp(5)
c
                  pdelrho(n,1,1,ipert) = pdelrho(n,1,1,ipert) +
     &                                   delrho(n,1,1)*t(1) + 
     &                                   delrho(n,1,2)*t(2)
c
                  pdelrho(n,2,1,ipert) = pdelrho(n,2,1,ipert) +
     &                                   delrho(n,2,1)*t(1) + 
     &                                   delrho(n,2,2)*t(2)
c
                  pdelrho(n,3,1,ipert) = pdelrho(n,3,1,ipert) +
     &                                   delrho(n,3,1)*t(1) + 
     &                                   delrho(n,3,2)*t(2)
c
                  t(1) = t(2)
                  t(2) =(Cmat2(n,D2_RA_GBB) *ptmp(1) +
     &                   Cmat2(n,D2_RB_GBB) *ptmp(2) +
     &                   Cmat2(n,D2_GAA_GBB)*ptmp(3) +
     &                   Cmat2(n,D2_GAB_GBB)*ptmp(4) +
     &                   Cmat2(n,D2_GBB_GBB)*ptmp(5))*2d0
c
                  pdelrho(n,1,2,ipert) = pdelrho(n,1,2,ipert) +
     &                                   delrho(n,1,1)*t(1) + 
     &                                   delrho(n,1,2)*t(2)
c
                  pdelrho(n,2,2,ipert) = pdelrho(n,2,2,ipert) +
     &                                   delrho(n,2,1)*t(1) + 
     &                                   delrho(n,2,2)*t(2)
c
                  pdelrho(n,3,2,ipert) = pdelrho(n,3,2,ipert) +
     &                                   delrho(n,3,1)*t(1) + 
     &                                   delrho(n,3,2)*t(2)
               enddo
            endif
c
         else if (TRIPLET) then ! restricted triplet case
c
c          --------------------------------------------------------------
c          Restricted triplet case
c          (TDDFT triplet excitation energies)
c          (DFT triplet instability check ... NYI)
c          --------------------------------------------------------------
c
            if (.not. grad) then                                     ! local functionals
               do n = 1, nq
                 term_rr = Amat2(n,D2_RA_RA) - Amat2(n,D2_RA_RB)     ! rara - rarb
                 prho(n,1,ipert) =term_rr*prho(n,1,ipert)
               enddo
            else      ! gradient dependent functionals
               do n = 1, nq
                ptmp(1) = prho(n,1,ipert)                              ! perturbed density
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert) +         ! delrho*perturbed density 
     &                    delrho(n,2,1)*pdelrho(n,2,1,ipert) +
     &                    delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                term_g  = Cmat(n,D1_GAA)*2.0d0 - Cmat(n,D1_GAB)      ! gaa - gab
                term_rr = Amat2(n,D2_RA_RA)    - Amat2(n,D2_RA_RB)   ! rara - rarb
                term_rg = Cmat2(n,D2_RA_GAA)   - Cmat2(n,D2_RA_GBB)  ! ragaa - ragbb
                term_gg = Cmat2(n,D2_GAA_GAA)  - Cmat2(n,D2_GAA_GBB) ! gaagaa - gaagbb
c
                if (kske) then  ! tau dependent functionals
                 term_rt = Mmat2(n,D2_RA_TA)  - Mmat2(n,D2_RA_TB)     ! rata - ratb
                 term_gt = Mmat2(n,D2_TA_GAA) - Mmat2(n,D2_TA_GBB)    ! gaata - gbbta
                 term_tt = Mmat2(n,D2_TA_TA)  - Mmat2(n,D2_TA_TB)     ! tata - tatb
                 term_prho = term_rr*ptmp(1) +  
     &                       term_rg*ptmp(2)*2.d0 + 
     &                       term_rt*pttau(n,1,ipert)
                 term_pdelrho = term_rg*ptmp(1) + 
     &                          term_gg*ptmp(2)*2.d0 +
     &                          term_gt*pttau(n,1,ipert)
                 term_pttau = term_tt*pttau(n,1,ipert) + 
     &                        term_gt*ptmp(2) +
     &                        term_rt*ptmp(1)
                 pttau(n,1,ipert) = term_pttau
                else 
                 term_prho = term_rr*ptmp(1) + term_rg*ptmp(2)*2.d0
                 term_pdelrho = term_rg*ptmp(1) + term_gg*ptmp(2)*2.d0
                end if
c
                prho(n,1,ipert) = term_prho
                pdelrho(n,1,1,ipert) = term_g*pdelrho(n,1,1,ipert) +
     &                                 term_pdelrho*delrho(n,1,1)*2.d0
                pdelrho(n,2,1,ipert) = term_g*pdelrho(n,2,1,ipert) +
     &                                 term_pdelrho*delrho(n,2,1)*2.d0
                pdelrho(n,3,1,ipert) = term_g*pdelrho(n,3,1,ipert) +
     &                                 term_pdelrho*delrho(n,3,1)*2.d0
               enddo
            endif
c
         else   ! singlet case
c
c          --------------------------------------------------------------
c          Restricted case
c          (DFT analytic hessian)
c          (TDDFT singlet excitation energies)
c          (DFT singlet instability check ... NYI)
c          --------------------------------------------------------------
c
            if (.not. grad) then                                   ! local functionals
               do n = 1, nq
                 term_rr = Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB)   ! rara + rarb
                 prho(n,1,ipert) =term_rr*prho(n,1,ipert)
               enddo
            else                                                   ! gradient dependent functionals
               do n = 1, nq
                ptmp(1) = prho(n,1,ipert)                          ! perturbed density
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert) +     ! delrho*perturbed density 
     &                    delrho(n,2,1)*pdelrho(n,2,1,ipert) +
     &                    delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                term_g  = Cmat(n,D1_GAA)*2.0d0 + Cmat(n,D1_GAB)    ! gaa + gab
                term_rr = Amat2(n,D2_RA_RA)    + Amat2(n,D2_RA_RB) ! rara + rarb
                term_rg = Cmat2(n,D2_RA_GAA) +                     ! ragaa + ragbb + ragab
     &                    Cmat2(n,D2_RA_GBB) +                     
     &                    Cmat2(n,D2_RA_GAB)                      
                term_gg = Cmat2(n,D2_GAA_GAA) +                    ! gaagaa + gaagbb + gaagab + gabgab
     &                    Cmat2(n,D2_GAA_GBB) + 
     &                    Cmat2(n,D2_GAA_GAB)*2.0d0 +
     &                    Cmat2(n,D2_GAB_GAB)*0.5d0
c
                if (kske) then
                 term_rt = Mmat2(n,D2_RA_TA)  + Mmat2(n,D2_RA_TB)   ! rata + ratb
                 term_gt = Mmat2(n,D2_TA_GAA) +                     ! gaata + gabta + gbbta
     &                     Mmat2(n,D2_TA_GAB) + 
     &                     Mmat2(n,D2_TA_GBB)
                 term_tt = Mmat2(n,D2_TA_TA)  + Mmat2(n,D2_TA_TB)   ! tata + tatb
                 term_prho = term_rr*ptmp(1) + 
     &                       term_rg*ptmp(2)*2.d0 + 
     &                       term_rt*pttau(n,1,ipert)
                 term_pdelrho = term_rg*ptmp(1) + 
     &                          term_gg*ptmp(2)*2.d0 +
     &                          term_gt*pttau(n,1,ipert)
                 term_pttau = term_tt*pttau(n,1,ipert) + 
     &                        term_gt*ptmp(2) +
     &                        term_rt*ptmp(1)
                 pttau(n,1,ipert) = term_pttau
                else
                 term_prho = term_rr*ptmp(1) + term_rg*ptmp(2)*2.d0
                 term_pdelrho = term_rg*ptmp(1) + term_gg*ptmp(2)*2.d0
                end if
c
                prho(n,1,ipert) = term_prho
                pdelrho(n,1,1,ipert) = term_g*pdelrho(n,1,1,ipert) +
     &                                 term_pdelrho*delrho(n,1,1)*2.d0
                pdelrho(n,2,1,ipert) = term_g*pdelrho(n,2,1,ipert) +
     &                                 term_pdelrho*delrho(n,2,1)*2.d0
                pdelrho(n,3,1,ipert) = term_g*pdelrho(n,3,1,ipert) +
     &                                 term_pdelrho*delrho(n,3,1)*2.d0
               enddo
            endif
         endif
      enddo
c
      if (ipol.eq.1 .and. grad) then
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
      if (grad) then
       write(6,*) 'xc_cpks_coeff: pdelrho out'
       call output(pdelrho, 1, nq, 1, 3*ipol*npert, nq, 3*ipol*npert, 1)
      endif
      endif
#endif
      return
      end
