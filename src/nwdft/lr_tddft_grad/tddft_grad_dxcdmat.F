      subroutine tddft_grad_dxcdmat(func, prho,
     1  pdelrho, nqpts, ipol, ngxc, grad, rho, delrho, 
     2  Amat, Cmat, Amat2, Cmat2, qwght, triplet, calc_type)
c      
c     $Id$
c
c
c Contracts the perturbed electron density with the XC-potential or
c XC-kernel, for  building the derivatives of the quadrature weights 
c in analytical  TDDFT gradients.  This routine mimics the code in 
c xc_cpks_coef.F (in src/hessian/analytic/dft/) for the case of the 
c 2nd derivative.
c
c Called from: dftg_gridv0b
c
      implicit none
c
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
#include "dft2drv.fh"
c
      integer ipol      ! 1=restricted, 2=unrestricted
      integer nqpts     ! Number of quadrature points
      integer ngxc      ! Number of perturbations
      integer ipert, n  ! Looping variables
      integer calc_type ! 1 = Construct XC-potential
                        ! 2 = Construct XC-kernel
c
      logical grad     ! Determines if the functional is a GGA
      logical triplet  ! Whether a triplet calculation is being done
c
      double precision func(nqpts)
      double precision prho(nqpts,ipol,ngxc)
      double precision pdelrho(nqpts,3,ipol,ngxc)
c
c Spin densities [input]
c
      double precision rho(nqpts,ipol)
c
c Gradients of spin densities [input]
c
      double precision delrho(nqpts,3,ipol)
c
c Quadrature Weights
c
      double precision qwght(nqpts)
c
c Second derivatives of XC functional [input]
c
      double precision Amat2(nqpts,NCOL_AMAT2)
      double precision Cmat2(nqpts,NCOL_CMAT2)
      double precision Mmat2(nqpts,NCOL_MMAT2)
c
c GC first derivatives of XC functional [input]
c
      double precision Amat(nqpts,ipol)
      double precision Cmat(nqpts,*)
c
      double precision term_r, term_g
      double precision term_rr, term_rg, term_gg
      double precision ptmp(9)    
      double precision pdra(3), pdrb(3)
      double precision t(2)
      double precision tm1a, tm1b
      double precision tm2a(3), tm2b(3)
c
      character*32 pname
c
      pname = "tddft_grad_dxcdmat: "
c
c Initialize ptmp
      do n = 1, 9
         ptmp(n) = 0.0d0
      end do
c
      if (ipol.eq.1) then
c
c Since the total densities are evaluated in the restricted case,
c scale them by a factor of 0.5 so that the correct CPKS matrices
c will be produced.
c
        call dscal(nqpts*ipol*ngxc,0.5d0,prho,1)
        if (grad) then
          call dscal(nqpts*3*ipol*ngxc,0.5d0,pdelrho,1)
          call dscal(nqpts*3*ipol,0.5d0,delrho,1)
        endif
      endif
c
c NOTE!!! For unrestricted jobs, the total perturbed density is given as:
c
c prho_tot = -1/sqrt(2)*(prho_alpha + prho_beta)
c
      do ipert = 1, ngxc
        if (ipol.eq.2) then
c -----------------
c Unrestricted case
c -----------------
          if (.not.grad) then ! local functionals
            do n = 1, nqpts
              if (calc_type.eq.1) then ! do nothing
c This part will be removed eventually
                func(n) = 0.0d0
              else if (calc_type.eq.2) then
                ptmp(1) = prho(n,1,ipert)
                ptmp(2) = prho(n,2,ipert)
c
                func(n) = Amat2(n,D2_RA_RA)*ptmp(1)*ptmp(1)
     1                  + 2.0d0*Amat2(n,D2_RA_RB)*ptmp(1)*ptmp(2)
     2                  + Amat2(n,D2_RB_RB)*ptmp(2)*ptmp(2)
              endif
            enddo
          else ! gradient dependent functionals
            do n = 1, nqpts
              if (calc_type.eq.1) then ! do nothing
c This part will be removed eventually
                func(n) = 0.0d0
              else if (calc_type.eq.2) then
c
c Perturbed densities
c
                ptmp(1)  = prho(n,1,ipert)                           ! prho_alpha
c
                ptmp(2)  = prho(n,2,ipert)                           ! prho_beta
c
c Products of the ground state density gradient and perturbed density
c gradient
c
                ptmp(3)  = delrho(n,1,1)*pdelrho(n,1,1,ipert)        ! delrho_alpha * pdelrho_alpha
     1                   + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     2                   + delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                ptmp(4)  = delrho(n,1,1)*pdelrho(n,1,2,ipert)        ! delrho_alpha * pdelrho_beta
     1                   + delrho(n,2,1)*pdelrho(n,2,2,ipert)
     2                   + delrho(n,3,1)*pdelrho(n,3,2,ipert)
c
                ptmp(5)  = delrho(n,1,2)*pdelrho(n,1,1,ipert)        ! delrho_beta * pdelrho_alpha
     1                   + delrho(n,2,2)*pdelrho(n,2,1,ipert)
     2                   + delrho(n,3,2)*pdelrho(n,3,1,ipert)
c
                ptmp(6)  = delrho(n,1,2)*pdelrho(n,1,2,ipert)        ! delrho_beta * pdelrho_beta
     1                   + delrho(n,2,2)*pdelrho(n,2,2,ipert)
     2                   + delrho(n,3,2)*pdelrho(n,3,2,ipert)
c
c Products of the perturbed density gradients
c
                ptmp(7) = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! pdelrho_alpha * pdelrho_alpha
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c
                ptmp(8) = pdelrho(n,1,1,ipert)*pdelrho(n,1,2,ipert) ! pdelrho_alpha * pdelrho_beta
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,2,ipert) ! OR
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,2,ipert) ! pdelrho_beta * pdelrho_alpha
c
                ptmp(9) = pdelrho(n,1,2,ipert)*pdelrho(n,1,2,ipert) ! pdelrho_beta * pdelrho_beta
     1                  + pdelrho(n,2,2,ipert)*pdelrho(n,2,2,ipert)
     2                  + pdelrho(n,3,2,ipert)*pdelrho(n,3,2,ipert)
c -------------------------------------------------------------------
c Second derivative contributions: drdr
c -------------------------------------------------------------------
              term_rr = Amat2(n,D2_RA_RA)*ptmp(1)*ptmp(1)
     1                + 2.0d0*Amat2(n,D2_RA_RB)*ptmp(1)*ptmp(2)
     2                + Amat2(n,D2_RB_RB)*ptmp(2)*ptmp(2)
c -------------------------------------------------------------------
c Second derivative contributions: drdg
c -------------------------------------------------------------------
              term_rg = 2.0d0*( 2.0d0*Cmat2(n,D2_RA_GAA)*ptmp(3)
     1                        + Cmat2(n,D2_RA_GAB)*ptmp(5) )*ptmp(1)     ! aa
     2                + 2.0d0*( 2.0d0*Cmat2(n,D2_RA_GBB)*ptmp(5)*ptmp(2)
     3                        + Cmat2(n,D2_RA_GAB)*ptmp(3)*ptmp(2)
     4                        + 2.0d0*Cmat2(n,D2_RB_GAA)*ptmp(4)*ptmp(1)
     5                        + Cmat2(n,D2_RB_GAB)*ptmp(6)*ptmp(1) )     ! ab and ba
     6                + 2.0d0*( 2.0d0*Cmat2(n,D2_RB_GBB)*ptmp(6)
     7                        + Cmat2(n,D2_RB_GAB)*ptmp(4) )*ptmp(2)     ! bb 
c -------------------------------------------------------------------
c Second derivative contributions: dgdg
c -------------------------------------------------------------------
              term_gg = ( 4.0d0*Cmat2(n,D2_GAA_GAA)*ptmp(3)*ptmp(3)
     1                  + 4.0d0*Cmat2(n,D2_GAA_GAB)*ptmp(5)*ptmp(3)
     2                  + Cmat2(n,D2_GAB_GAB)*ptmp(5)*ptmp(5) )          ! aa 
     3                + 2.0d0*( 4.0d0*Cmat2(n,D2_GAA_GBB)
     4                          *ptmp(5)*ptmp(4)
     5                        + 2.0d0*Cmat2(n,D2_GAA_GAB)
     6                          *ptmp(3)*ptmp(4)
     7                        + 2.0d0*Cmat2(n,D2_GAB_GBB)
     8                          *ptmp(5)*ptmp(6)
     9                        + Cmat2(n,D2_GAB_GAB)
     A                          *ptmp(3)*ptmp(6) )                       ! ab and ba
     B                + ( 4.0d0*Cmat2(n,D2_GBB_GBB)*ptmp(6)*ptmp(6)
     C                  + 4.0d0*Cmat2(n,D2_GAB_GBB)*ptmp(4)*ptmp(6)
     D                  + Cmat2(n,D2_GAB_GAB)*ptmp(4)*ptmp(4) )          ! bb
c -------------------------------------------------------------------
c First derivative contributions: dg
c -------------------------------------------------------------------
             term_g = 2.0d0*( Cmat(n,D1_GAA)*ptmp(7)
     1                      + Cmat(n,D1_GAB)*ptmp(8)
     2                      + Cmat(n,D1_GBB)*ptmp(9) )
c Build func (all alpha terms + all beta terms).  All of the work
c was done above so this should function as expected.
                func(n) = term_rr + term_rg + term_gg + term_g 
              endif
            enddo
          endif
        else if (triplet) then
c -----------------------
c Restricted Triplet case
c -----------------------
          if (.not.grad) then ! local functionals
            do n = 1, nqpts
              if (calc_type.eq.1) then ! do nothing
c This part will be removed eventually
                func(n) = 0.0d0
              else if (calc_type.eq.2) then
                term_rr = Amat2(n,D2_RA_RA) - Amat2(n,D2_RA_RB)     ! rara - rarb
                func(n) = term_rr*prho(n,1,ipert)*prho(n,1,ipert)
              endif
            enddo
          else ! gradient dependent functionals
            do n = 1, nqpts
              if (calc_type.eq.1) then ! do nothing
c This part will be removed eventually
                func(n) = 0.0d0
              else if (calc_type.eq.2) then
                ptmp(1) = prho(n,1,ipert)                           ! perturbed density
c
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert) +      ! delrho*perturbed density 
     1                    delrho(n,2,1)*pdelrho(n,2,1,ipert) +
     2                    delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                ptmp(3) = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! del-perturbed density*del-perturbed density
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c
                term_g  = Cmat(n,D1_GAA)*2.0d0 - Cmat(n,D1_GAB)      ! 2*gaa - gab
                term_rr = Amat2(n,D2_RA_RA)    - Amat2(n,D2_RA_RB)   ! rara - rarb
                term_rg = Cmat2(n,D2_RA_GAA)   - Cmat2(n,D2_RA_GBB)  ! ragaa - ragbb
                term_gg = Cmat2(n,D2_GAA_GAA)  - Cmat2(n,D2_GAA_GBB) ! gaagaa - gaagbb
c
                func(n) = term_rr*ptmp(1)*ptmp(1)       ! term_rr*prho*prho
     1                  + term_rg*ptmp(2)*ptmp(1)*2.0d0 ! term_rg*(delprho*delrho)*prho
     2                  + term_rg*ptmp(1)*ptmp(2)*2.0d0 ! term_rg*(delprho*delrho)*prho
     3                  + term_gg*ptmp(2)*ptmp(2)*4.0d0 ! term_gg*(delprho*delrho)*(delprho*delrho)
     4                  + term_g*ptmp(3)                ! term_g*(delprho*delprho)
              endif
            enddo
          endif
        else
c ---------------
c Restricted case
c ---------------
          if (.not.grad) then ! local functionals 
            do n = 1, nqpts
c ---------------------------------------------------------------------
c XC-potential contracted with the relaxed one-particle difference  
c density (P)
c ---------------------------------------------------------------------
              if (calc_type.eq.1) then
c This part will be removed eventually
c                term_r = Amat(n,1)
c                func(n) = term_r*prho(n,1,ipert)                
                func(n) = 0.0d0
c ---------------------------------------------------------------------
c XC-kernel contracted with the square of the excitation vector, (X+Y), 
c density (also called the transition density).
c ---------------------------------------------------------------------
              else if (calc_type.eq.2) then
                term_rr = 
     1            Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB)
                func(n) = term_rr*prho(n,1,ipert)*prho(n,1,ipert)
              endif
            enddo
          else ! gradient dependent functionals 
            do n = 1, nqpts
c ---------------------------------------------------------------------
c XC-potential contracted with the relaxed one-particle difference  
c density (P)
c ---------------------------------------------------------------------
              if (calc_type.eq.1) then
c This part will be removed eventually
c
c                ptmp(1) = prho(n,1,ipert)                       ! perturbed density
cc
c                ptmp(2) = pdelrho(n,1,1,ipert)*delrho(n,1,1)    ! del-perturbed density*delrho
c     1                  + pdelrho(n,2,1,ipert)*delrho(n,2,1)
c     2                  + pdelrho(n,3,1,ipert)*delrho(n,3,1)
cc
c                term_r = Amat(n,1)                              ! ra
cc
c                term_g = 2.0d0*Cmat(n,D1_GAA) + Cmat(n,D1_GAB)  ! 2*gaa + gab
cc
cc I'm not sure why a factor of 0.5 is needed for the gradient term, but
cc I must have lost track of a factor of 2 somewhere in the code.  This
cc gives the correct answer for everything I've tested so far.
c                func(n) = term_r*ptmp(1) + 0.50d0*term_g*ptmp(2)
                func(n) = 0.0d0
c ---------------------------------------------------------------------
c XC-kernel contracted with the excitation vector density matrix,
c (X+Y)
c ---------------------------------------------------------------------
              elseif (calc_type.eq.2) then
                ptmp(1) = prho(n,1,ipert)                    ! perturbed density
c
                ptmp(2) = pdelrho(n,1,1,ipert)*delrho(n,1,1) ! del-perturbed density*delrho
     1                  + pdelrho(n,2,1,ipert)*delrho(n,2,1)
     2                  + pdelrho(n,3,1,ipert)*delrho(n,3,1)
c
                ptmp(3) = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! del-perturbed density*del-perturbed density
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c
                term_rr = Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB) ! rara + rarb
c
                term_rg = Cmat2(n,D2_RA_GAA)                    ! ragaa + ragab + ragbb
     1                  + Cmat2(n,D2_RA_GAB)
     2                  + Cmat2(n,D2_RA_GBB)
c
                term_gg = Cmat2(n,D2_GAA_GAA)                   ! gaagaa + gaagbb + 2*gaagab + 0.5*gabgab
     1                  + Cmat2(n,D2_GAA_GBB)
     2                  + 2.0d0*Cmat2(n,D2_GAA_GAB)
     3                  + 0.5d0*Cmat2(n,D2_GAB_GAB)
c
                term_g = 2.0d0*Cmat(n,D1_GAA) + Cmat(n,D1_GAB)  ! 2*gaa + gab
c
                func(n) = term_rr*ptmp(1)*ptmp(1)       ! term_rr*prho*prho
     1                  + term_rg*ptmp(2)*ptmp(1)*2.0d0 ! term_rg*(delprho*delrho)*prho
     2                  + term_rg*ptmp(1)*ptmp(2)*2.0d0 ! term_rg*(delprho*delrho)*prho
     3                  + term_gg*ptmp(2)*ptmp(2)*4.0d0 ! term_gg*(delprho*delrho)*(delprho*delrho)
     4                  + term_g*ptmp(3)                ! term_g*(delprho*delprho)
              endif
            enddo
          endif
        endif
      enddo ! ipert

      if (ipol.eq.1) then
c
c Put delrho back the way it was since it may be used later on
c
        call dscal(nqpts*ipol*ngxc,2.0d0,prho,1)
        if (grad) then
          call dscal(nqpts*3*ipol*ngxc,2.0d0,pdelrho,1)
          call dscal(nqpts*3*ipol,2.0d0,delrho,1)
        endif
      endif

c ------
c Return 
c ------
      return
      end
