c
c Calculates the terms related to the gradient of the XC-kernel.  
c Called from: dftg_gridv0b
c
      subroutine tddft_grad_dfxc(prho, pdelrho,
     1  Amat3, Cmat3, Amat2, Cmat2, Amat, Cmat, 
     2  delrho, npert, ipol, nq, 
     3  grad, crho, cdelrho, calctype, triplet)
c
      implicit none
#include "rtdb.fh"
#include "dft2drv.fh"
#include "dft3drv.fh"
#include "stdio.fh"
#include "errquit.fh"
c
      integer npert    ! Number CPKS perturbations [input]
      integer ipol     ! Restricted or unrestricted [input]
      integer nq       ! Number of grid points [input]
      integer calctype ! 1 = XC-coefficient to be multiplied with the GS density gradient
                       ! 2 = XC-coefficient to be multiplied with the perturbed density gradient
c
      logical grad     ! True if a GGA is being used [input]
      logical triplet  ! True if a triplet calculation is being done.
c
c Current approximate perturbed spin densities and density gradients
c
      double precision prho(nq,ipol,npert)
      double precision pdelrho(nq,3,ipol,npert)
c
c Matrices for overwriting with XC matrix coefficients
c
      double precision crho(nq,ipol,npert)
      double precision cdelrho(nq,3,ipol,npert)
c
c Third derivatives of XC functional [input]
c
      double precision Amat3(nq,NCOL_AMAT3)
      double precision Cmat3(nq,NCOL_CMAT3)
c
c Second derivatives of XC functional [input]
c
      double precision Amat2(nq,NCOL_AMAT2)
      double precision Cmat2(nq,NCOL_CMAT2)
c
c First derivatives of XC functional [input]
c
      double precision Amat(nq,ipol)
c      double precision Cmat(nq,*)
      double precision Cmat(nq,3)
c
c Gradients of spin densities [input]
c
      double precision delrho(nq,3,ipol)
c
      integer ipert, n ! Loop indices
c
      double precision ptmp(10)
      double precision t(2)
      double precision pdra(3), pdrb(3)
      double precision term_rrr, term_rrg, term_rgg, term_ggg
      double precision term_rg, term_gg, term_rr, term_r, term_g
      double precision term_prho, term_pdelrho
      double precision term_rrg1, term_rrg2, term_rrg3
      double precision term_rgg1, term_rgg2, term_rgg3
      double precision term_rg1, term_rg2, term_rg3
      double precision term_gg1, term_gg2, term_gg3
c
c Unrestricted specific variables
c
      double precision term_rrr1a, term_rrr1b
      double precision term_rrg1a, term_rrg2a, term_rrg3a, term_rrg4a
      double precision term_rrg1b, term_rrg2b, term_rrg3b, term_rrg4b
      double precision term_rgg1a, term_rgg2a, term_rgg3a, term_rgg4a,
     1                 term_rgg5a, term_rgg6a
      double precision term_rgg1b, term_rgg2b, term_rgg3b, term_rgg4b,
     1                 term_rgg5b, term_rgg6b
      double precision term_ggg1a, term_ggg2a, term_ggg3a, term_ggg4a,
     1                 term_ggg5a, term_ggg6a, term_ggg7a, term_ggg8a
      double precision term_ggg1b, term_ggg2b, term_ggg3b, term_ggg4b,
     1                 term_ggg5b, term_ggg6b, term_ggg7b, term_ggg8b
      double precision term_rg1a, term_rg2a, term_rg3a
      double precision term_rg1b, term_rg2b, term_rg3b
      double precision term_gg1a, term_gg2a, term_gg3a, term_gg4a
      double precision term_gg1b, term_gg2b, term_gg3b, term_gg4b
      double precision term_pdelrho1a, term_pdelrho1b
      double precision term_pdelrho2a, term_pdelrho2b
c
      character*32 pname
c
      pname = "tddft_grad_dfxc: "
c
c Initialize ptmp
c
      do n = 1,10
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
c Check calctype
c
      if (.not.(calctype.eq.1.or.calctype.eq.2))
     1  call errquit(pname//'calctype should be 1 or 2', calctype,
     2    CAPMIS_ERR)
c
c NOTE!!! For unrestricted jobs, the total perturbed density is given as:
c
c prho_tot = -1/sqrt(2)*(prho_alpha + prho_beta)
c
      do ipert = 1, npert
c
c !!! Put in cutoffs here similar to xc_tabcd? !!!
c
        if (ipol.eq.2) then ! unrestricted case
c --------------------------------------------------------------
c Unrestricted case
c (XC-kernel gradients)
c --------------------------------------------------------------
          if (.not. grad) then ! local functionals
            do n = 1, nq
              if (calctype.eq.1) then
                ptmp(1) = prho(n,1,ipert)*prho(n,1,ipert)
                ptmp(2) = prho(n,1,ipert)*prho(n,2,ipert)
                ptmp(3) = prho(n,2,ipert)*prho(n,2,ipert)
c Here, we multiply by the perturbed density twice, mimicing what was
c done for the XC-kernel where the perturbed density is multiplied in
c once. 
                crho(n,1,ipert) = Amat3(n,D3_RA_RA_RA)*ptmp(1)
     1                          + 2.0d0*Amat3(n,D3_RA_RA_RB)*ptmp(2)
     2                          + Amat3(n,D3_RA_RB_RB)*ptmp(3)
c
                crho(n,2,ipert) = Amat3(n,D3_RA_RA_RB)*ptmp(1)
     1                          + 2.0d0*Amat3(n,D3_RA_RB_RB)*ptmp(2)
     2                          + Amat3(n,D3_RB_RB_RB)*ptmp(3)
              else if (calctype.eq.2) then
                ptmp(1) = prho(n,1,ipert)
                ptmp(2) = prho(n,2,ipert)
                crho(n,1,ipert) = Amat2(n,D2_RA_RA)*ptmp(1)
     &                          + Amat2(n,D2_RA_RB)*ptmp(2)
                crho(n,2,ipert) = Amat2(n,D2_RA_RB)*ptmp(1)
     &                          + Amat2(n,D2_RB_RB)*ptmp(2)
              endif
            enddo
          else ! gradient dependent functionals
c            call errquit(
c     1        pname//'Unrestricted GGA calculations not implemented',
c     2        0, MA_ERR)
            do n = 1, nq
              if (calctype.eq.1) then
c
c Perturbed densities
c
                ptmp(1)  = prho(n,1,ipert)                           ! prho_alpha
c
                ptmp(2)  = prho(n,2,ipert)                           ! prho_beta
c
                ptmp(3)  = prho(n,1,ipert) + prho(n,2,ipert)         ! prho_alpha + prho_beta
c
c Products of the ground state density gradient and perturbed density
c gradient
c
                ptmp(4)  = delrho(n,1,1)*pdelrho(n,1,1,ipert)        ! delrho_alpha * pdelrho_alpha
     1                   + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     2                   + delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                ptmp(5)  = delrho(n,1,1)*pdelrho(n,1,2,ipert)        ! delrho_alpha * pdelrho_beta
     1                   + delrho(n,2,1)*pdelrho(n,2,2,ipert)
     2                   + delrho(n,3,1)*pdelrho(n,3,2,ipert)
c
                ptmp(6)  = delrho(n,1,2)*pdelrho(n,1,1,ipert)        ! delrho_beta * pdelrho_alpha
     1                   + delrho(n,2,2)*pdelrho(n,2,1,ipert)
     2                   + delrho(n,3,2)*pdelrho(n,3,1,ipert)
c
                ptmp(7)  = delrho(n,1,2)*pdelrho(n,1,2,ipert)        ! delrho_beta * pdelrho_beta
     1                   + delrho(n,2,2)*pdelrho(n,2,2,ipert)
     2                   + delrho(n,3,2)*pdelrho(n,3,2,ipert)
c
c Products of the perturbed density gradients
c
                ptmp(8)  = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! pdelrho_alpha * pdelrho_alpha
     1                   + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                   + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c
                ptmp(9)  = pdelrho(n,1,1,ipert)*pdelrho(n,1,2,ipert) ! pdelrho_alpha * pdelrho_beta
     1                   + pdelrho(n,2,1,ipert)*pdelrho(n,2,2,ipert) ! OR
     2                   + pdelrho(n,3,1,ipert)*pdelrho(n,3,2,ipert) ! pdelrho_beta * pdelrho_alpha
c
                ptmp(10) = pdelrho(n,1,2,ipert)*pdelrho(n,1,2,ipert) ! pdelrho_beta * pdelrho_beta
     1                   + pdelrho(n,2,2,ipert)*pdelrho(n,2,2,ipert)
     2                   + pdelrho(n,3,2,ipert)*pdelrho(n,3,2,ipert)
c -------------------------------------------------------------------
c Third derivative contributions: drdrdr
c -------------------------------------------------------------------
                term_rrr1a = Amat3(n,D3_RA_RA_RA)*ptmp(1)*ptmp(1)
     1                     + 2.0d0*Amat3(n,D3_RA_RA_RB)*ptmp(1)*ptmp(2)
     2                     + Amat3(n,D3_RA_RB_RB)*ptmp(2)*ptmp(2)
c
                term_rrr1b = Amat3(n,D3_RA_RA_RB)*ptmp(1)*ptmp(1)
     1                     + 2.0d0*Amat3(n,D3_RA_RB_RB)*ptmp(1)*ptmp(2)
     2                     + Amat3(n,D3_RB_RB_RB)*ptmp(2)*ptmp(2)
c -------------------------------------------------------------------
c Third derivative contributions: drdrdg
c -------------------------------------------------------------------
c Alpha spin terms
                term_rrg1a = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_RA_GAA)*ptmp(4)
     1                             + Cmat3(n,D3_RA_RA_GAB)
     2                               *( ptmp(6) + ptmp(5) )
     3                             + 2.0d0*Cmat3(n,D3_RA_RA_GBB)
     4                               *ptmp(7) )
c
                term_rrg2a = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_RB_GAA)*ptmp(4)
     1                             + Cmat3(n,D3_RA_RB_GAB)
     2                               *( ptmp(6) + ptmp(5) )
     3                             + 2.0d0*Cmat3(n,D3_RA_RB_GBB)
     4                               *ptmp(7) )
c
                term_rrg3a = 2.0d0*( Cmat3(n,D3_RA_RA_GAA)
     1                               *ptmp(1)*ptmp(1)
     2                             + 2.0d0*Cmat3(n,D3_RA_RB_GAA)
     3                               *ptmp(2)*ptmp(1)
     4                             + Cmat3(n,D3_RB_RB_GAA)
     5                               *ptmp(2)*ptmp(2) )
c
                term_rrg4a = Cmat3(n,D3_RA_RA_GAB)*ptmp(1)*ptmp(1)
     1                     + 2.0d0*Cmat3(n,D3_RA_RB_GAB)*ptmp(2)*ptmp(1)
     2                     + Cmat3(n,D3_RB_RB_GAB)*ptmp(2)*ptmp(2)
c
c Beta spin terms
                term_rrg1b = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_RB_GAA)*ptmp(4)
     1                             + Cmat3(n,D3_RA_RB_GAB)
     2                               *( ptmp(6) + ptmp(5) )
     3                             + 2.0d0*Cmat3(n,D3_RA_RB_GBB)
     4                               *ptmp(7) )
c
                term_rrg2b = 2.0d0*( 2.0d0*Cmat3(n,D3_RB_RB_GAA)*ptmp(4)
     1                             + Cmat3(n,D3_RB_RB_GAB)
     2                               *( ptmp(6) + ptmp(5) )
     3                             + 2.0d0*Cmat3(n,D3_RB_RB_GBB)
     4                               *ptmp(7) )
c
                term_rrg3b = Cmat3(n,D3_RA_RA_GAB)*ptmp(1)*ptmp(1)
     1                     + 2.0d0*Cmat3(n,D3_RA_RB_GAB)*ptmp(2)*ptmp(1)
     2                     + Cmat3(n,D3_RB_RB_GAB)*ptmp(2)*ptmp(2)
c 
                term_rrg4b = 2.0d0*( Cmat3(n,D3_RA_RA_GBB)
     1                               *ptmp(1)*ptmp(1)
     2                             + 2.0d0*Cmat3(n,D3_RA_RB_GBB)
     3                               *ptmp(2)*ptmp(1)
     4                             + Cmat3(n,D3_RB_RB_GBB)
     5                               *ptmp(2)*ptmp(2) )
c -------------------------------------------------------------------
c Third derivative contributions: drdgdg
c -------------------------------------------------------------------
c Alpha spin terms
                term_rgg1a = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_GAA_GAA)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RA_GAA_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_RA_GAA_GBB)
     5                               *ptmp(7) )
c
                term_rgg2a = 2.0d0*Cmat3(n,D3_RA_GAA_GAB)*ptmp(4)
     1                     + Cmat3(n,D3_RA_GAB_GAB)
     2                       *( ptmp(6) + 2.0d0*ptmp(5) )
     3                     + 4.0d0*Cmat3(n,D3_RA_GAB_GBB)*ptmp(7)
c
                term_rgg3a = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_RA_GAB_GBB)*ptmp(5) )
c
                term_rgg4a = 2.0d0*Cmat3(n,D3_RA_GAB_GBB)*ptmp(7)
     1                     + Cmat3(n,D3_RA_GAB_GAB)*ptmp(5)
c
                term_rgg5a = 4.0d0*( 2.0d0*Cmat3(n,D3_RA_GAA_GAA)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RA_GAA_GAB)
     3                               *( ptmp(6) + ptmp(5) )
     4                             + 2.0d0*Cmat3(n,D3_RA_GAA_GBB)
     5                               *ptmp(7) )*ptmp(1)
     6                     + 4.0d0*( 2.0d0*Cmat3(n,D3_RB_GAA_GAA)
     7                               *ptmp(4)
     8                             + Cmat3(n,D3_RB_GAA_GAB)
     9                               *( ptmp(6) + ptmp(5) )
     A                             + 2.0d0*Cmat3(n,D3_RB_GAA_GBB)
     B                               *ptmp(7) )*ptmp(2)
c
                term_rgg6a = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_GAA_GAB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RA_GAB_GAB)
     3                               *( ptmp(6) + ptmp(5) )
     4                             + 2.0d0*Cmat3(n,D3_RA_GAB_GBB)
     5                               *ptmp(7) )*ptmp(1)
     6                     + 2.0d0*( 2.0d0*Cmat3(n,D3_RB_GAA_GAB)
     7                               *ptmp(4)
     8                             + Cmat3(n,D3_RB_GAB_GAB)
     9                               *( ptmp(6) + ptmp(5) )
     A                             + 2.0d0*Cmat3(n,D3_RB_GAB_GBB)
     B                               *ptmp(7) )*ptmp(2)
c
c Beta spin terms
                term_rgg1b = 2.0d0*( 2.0d0*Cmat3(n,D3_RB_GAA_GAA)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RB_GAA_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_RB_GAA_GBB)
     5                               *ptmp(7) )
c
                term_rgg2b = 2.0d0*Cmat3(n,D3_RB_GAA_GAB)*ptmp(4)
     1                     + Cmat3(n,D3_RB_GAB_GAB)
     2                       *( ptmp(6) + 2.0d0*ptmp(5) )
     3                     + 4.0d0*Cmat3(n,D3_RB_GAB_GBB)*ptmp(7)
c
                term_rgg3b = 2.0d0*( 2.0d0*Cmat3(n,D3_RB_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_RB_GAB_GBB)*ptmp(5) )
c
                term_rgg4b = 2.0d0*Cmat3(n,D3_RB_GAB_GBB)*ptmp(7)
     1                     + Cmat3(n,D3_RB_GAB_GAB)*ptmp(5)
c
                term_rgg5b = 4.0d0*( 2.0d0*Cmat3(n,D3_RA_GAA_GBB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RA_GAB_GBB)
     3                               *( ptmp(6) + ptmp(5) )
     4                             + 2.0d0*Cmat3(n,D3_RA_GBB_GBB)
     5                               *ptmp(7) )*ptmp(1)
     6                     + 4.0d0*( 2.0d0*Cmat3(n,D3_RB_GAA_GBB)
     7                               *ptmp(4)
     8                             + Cmat3(n,D3_RB_GAB_GBB)
     9                               *( ptmp(6) + ptmp(5) )
     A                             + 2.0d0*Cmat3(n,D3_RB_GBB_GBB)
     B                               *ptmp(7) )*ptmp(2)
c
                term_rgg6b = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_GAA_GAB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_RA_GAB_GAB)
     3                               *( ptmp(6) + ptmp(5) )
     4                             + 2.0d0*Cmat3(n,D3_RA_GAB_GBB)
     5                               *ptmp(7) )*ptmp(1)
     6                     + 2.0d0*( 2.0d0*Cmat3(n,D3_RB_GAA_GAB)
     7                               *ptmp(4)
     8                             + Cmat3(n,D3_RB_GAB_GAB)
     9                               *( ptmp(6) + ptmp(5) )
     A                             + 2.0d0*Cmat3(n,D3_RB_GAB_GBB)
     B                               *ptmp(7) )*ptmp(2)
c              term_rgg6b = term_rgg6a
c -------------------------------------------------------------------
c Third derivative contributions: dgdgdg
c -------------------------------------------------------------------
c Alpha spin terms
                term_ggg1a = 4.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GAA)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAA_GAA_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAA_GAA_GBB)
     5                               *ptmp(7) )
c
                term_ggg2a = 2.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GAB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAA_GAB_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     5                               *ptmp(7) )
c
                term_ggg3a = 4.0d0*( 2.0d0*Cmat3(n,D3_GAA_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAA_GAB_GBB)*ptmp(5) )
c
                term_ggg4a = 2.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAA_GAB_GAB)*ptmp(5) )
c
                term_ggg5a = 2.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GAB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAA_GAB_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     5                               *ptmp(7) )
c
                term_ggg6a = 2.0d0*Cmat3(n,D3_GAA_GAB_GAB)*ptmp(4)
     1                     + Cmat3(n,D3_GAB_GAB_GAB)
     2                       *( ptmp(6) + 2.0d0*ptmp(5) )
     3                     + 4.0d0*Cmat3(n,D3_GAB_GAB_GBB)*ptmp(7)
c
                term_ggg7a = 2.0d0*( 2.0d0*Cmat3(n,D3_GAB_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAB_GAB_GBB)
     3                               *ptmp(5) )
c
                term_ggg8a = 2.0d0*Cmat3(n,D3_GAB_GAB_GBB)*ptmp(7)
     1                     + Cmat3(n,D3_GAB_GAB_GAB)*ptmp(5)
c
c Beta spin terms
                term_ggg1b = 2.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GAB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAA_GAB_GAB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     5                               *ptmp(7) )
c
                term_ggg2b = 2.0d0*Cmat3(n,D3_GAA_GAB_GAB)*ptmp(4)
     1                     + Cmat3(n,D3_GAB_GAB_GAB)
     2                       *( ptmp(6) + 2.0d0*ptmp(5) )
     3                     + 4.0d0*Cmat3(n,D3_GAB_GAB_GBB)*ptmp(7)
c
                term_ggg3b = 2.0d0*( 2.0d0*Cmat3(n,D3_GAB_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAB_GAB_GBB)*ptmp(5) )
c
                term_ggg4b = 2.0d0*Cmat3(n,D3_GAB_GAB_GBB)*ptmp(7)
     1                     + Cmat3(n,D3_GAB_GAB_GAB)*ptmp(5)
c
                term_ggg5b = 4.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GBB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAA_GAB_GBB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAA_GBB_GBB)
     5                               *ptmp(7) )
c
                term_ggg6b = 2.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     1                               *ptmp(4)
     2                             + Cmat3(n,D3_GAB_GAB_GBB)
     3                               *( ptmp(6) + 2.0d0*ptmp(5) )
     4                             + 4.0d0*Cmat3(n,D3_GAB_GBB_GBB)
     5                               *ptmp(7) )
c
                term_ggg7b = 4.0d0*( 2.0d0*Cmat3(n,D3_GBB_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAB_GBB_GBB)*ptmp(5) )
c
                term_ggg8b = 2.0d0*( 2.0d0*Cmat3(n,D3_GAB_GBB_GBB)
     1                               *ptmp(7)
     2                             + Cmat3(n,D3_GAB_GAB_GBB)*ptmp(5) )
c -------------------------------------------------------------------
c Third derivative contributions: drdg
c -------------------------------------------------------------------
c Alpha spin terms
                term_rg1a = 2.0d0*( Cmat2(n,D2_RA_GAA)*ptmp(8)
     1                            + Cmat2(n,D2_RA_GAB)*ptmp(9)
     2                            + Cmat2(n,D2_RA_GBB)*ptmp(10) )
c
                term_rg2a = 4.0d0*( Cmat2(n,D2_RA_GAA)*ptmp(1)
     1                            + Cmat2(n,D2_RB_GAA)*ptmp(2) )
c
                term_rg3a = 2.0d0*( Cmat2(n,D2_RA_GAB)*ptmp(1)
     1                            + Cmat2(n,D2_RB_GAB)*ptmp(2) )
c
c Beta spin terms
                term_rg1b = 2.0d0*( Cmat2(n,D2_RB_GAA)*ptmp(8)
     1                            + Cmat2(n,D2_RB_GAB)*ptmp(9)
     2                            + Cmat2(n,D2_RB_GBB)*ptmp(10) )
c
                term_rg2b = 2.0d0*( Cmat2(n,D2_RA_GAB)*ptmp(1)
     1                            + Cmat2(n,D2_RB_GAB)*ptmp(2) )
c
                term_rg3b = 4.0d0*( Cmat2(n,D2_RA_GBB)*ptmp(1)
     1                            + Cmat2(n,D2_RB_GBB)*ptmp(2) )
c -------------------------------------------------------------------
c Third derivative contributions: dgdg
c -------------------------------------------------------------------
c Alpha spin terms
                term_gg1a = 4.0d0*( 2.0d0*Cmat2(n,D2_GAA_GAA)*ptmp(4)
     1                            + Cmat2(n,D2_GAA_GAB)
     2                              *( ptmp(6) + ptmp(5) )
     3                            + 2.0d0*Cmat2(n,D2_GAA_GBB)*ptmp(7) )
c
                term_gg2a = 4.0d0*( Cmat2(n,D2_GAA_GAA)*ptmp(8)
     1                            + Cmat2(n,D2_GAA_GAB)*ptmp(9)
     2                            + Cmat2(n,D2_GAA_GBB)*ptmp(10) )
c
                term_gg3a = 2.0d0*( 2.0d0*Cmat2(n,D2_GAA_GAB)*ptmp(4)
     1                            + Cmat2(n,D2_GAB_GAB)
     2                              *( ptmp(6) + ptmp(5) )
     3                            + 2.0d0*Cmat2(n,D2_GAB_GBB)*ptmp(7) )
c
                term_gg4a = 2.0d0*( Cmat2(n,D2_GAA_GAB)*ptmp(8)
     1                            + Cmat2(n,D2_GAB_GAB)*ptmp(9)
     2                            + Cmat2(n,D2_GAB_GBB)*ptmp(10) )
c
c Beta spin terms
                term_gg1b = 2.0d0*( 2.0d0*Cmat2(n,D2_GAA_GAB)*ptmp(4)
     1                            + Cmat2(n,D2_GAB_GAB)
     2                              *( ptmp(6) + ptmp(5) )
     3                            + 2.0d0*Cmat2(n,D2_GAB_GBB)*ptmp(7) )
c
                term_gg2b = 2.0d0*( Cmat2(n,D2_GAA_GAB)*ptmp(8)
     1                            + Cmat2(n,D2_GAB_GAB)*ptmp(9)
     2                            + Cmat2(n,D2_GAB_GBB)*ptmp(10) )
c
                term_gg3b = 4.0d0*( 2.0d0*Cmat2(n,D2_GAA_GBB)*ptmp(4)
     1                            + Cmat2(n,D2_GAB_GBB)
     2                              *( ptmp(6) + ptmp(5) )
     3                            + 2.0d0*Cmat2(n,D2_GBB_GBB)*ptmp(7) )
c
                term_gg4b = 4.0d0*( Cmat2(n,D2_GAA_GBB)*ptmp(8)
     1                            + Cmat2(n,D2_GAB_GBB)*ptmp(9)
     2                            + Cmat2(n,D2_GBB_GBB)*ptmp(10) )
c -------------------------------------------------------------------
c Terms that do not involve the gradient with respect to the basis
c functions.  
c -------------------------------------------------------------------
c
c Recall that:
c  ptmp(1)  = prho_alpha
c  ptmp(2)  = prho_beta
c  ptmp(3)  = prho_alpha + prho_beta
c  ptmp(4)  = delrho_alpha * pdelrho_alpha
c  ptmp(5)  = delrho_alpha * pdelrho_beta
c  ptmp(6)  = delrho_beta  * pdelrho_alpha
c  ptmp(7)  = delrho_beta  * pdelrho_beta 
c  ptmp(8)  = pdelrho_alpha * pdelrho_alpha
c  ptmp(9)  = pdelrho_alpha * pdelrho_beta = pdelrho_beta * pdelrho_alpha
c  ptmp(10) = pdelrho_beta * pdelrho_beta
c
              crho(n,1,ipert) = term_rrr1a
     1                        + term_rrg1a*ptmp(1) + term_rrg2a*ptmp(2)
     2                        + term_rgg1a*ptmp(4) + term_rgg2a*ptmp(6)
     3                        + term_rgg3a*ptmp(7) + term_rgg4a*ptmp(5)
     4                        + term_rg1a
c
              crho(n,2,ipert) = term_rrr1b
     1                        + term_rrg1b*ptmp(1) + term_rrg2b*ptmp(2)
     2                        + term_rgg1b*ptmp(4) + term_rgg2b*ptmp(6)
     3                        + term_rgg3b*ptmp(7) + term_rgg4b*ptmp(5)
     4                        + term_rg1b
c -------------------------------------------------------------------
c Terms that involve the gradient with respect to the basis 
c functions.
c -------------------------------------------------------------------
c These terms are multiplied by the ground state density gradient
c Multiplied by delrho_alpha (alpha coefficient)
                term_pdelrho1a = term_rrg3a
     1                         + term_rgg5a
     2                         + term_ggg1a*ptmp(4) + term_ggg2a*ptmp(6)
     3                         + term_ggg3a*ptmp(7) + term_ggg4a*ptmp(5)
     4                         + term_gg2a
c Multiplied by delrho_beta (alpha coefficient)
                term_pdelrho1b = term_rrg4a
     1                         + term_rgg6a
     2                         + term_ggg5a*ptmp(4) + term_ggg6a*ptmp(6)
     3                         + term_ggg7a*ptmp(7) + term_ggg8a*ptmp(5)
     4                         + term_gg4a
c Multiplied by delrho_alpha (beta coefficient)
                term_pdelrho2a = term_rrg3b
     1                         + term_rgg5b
     2                         + term_ggg1b*ptmp(4) + term_ggg2b*ptmp(6)
     3                         + term_ggg3b*ptmp(7) + term_ggg4b*ptmp(5)
     4                         + term_gg2b
c Multiplied by delrho_beta (beta coefficient)
                term_pdelrho2b = term_rrg4b
     1                         + term_rgg6b
     2                         + term_ggg5b*ptmp(4) + term_ggg6b*ptmp(6)
     3                         + term_ggg7b*ptmp(7) + term_ggg8b*ptmp(5)
     4                         + term_gg4b
c -------------------------------------------------------------------
c Daniel (3-19-13): These can be written more efficiently
                pdra(1) = pdelrho(n,1,1,ipert)
                pdra(2) = pdelrho(n,2,1,ipert)
                pdra(3) = pdelrho(n,3,1,ipert)
                pdrb(1) = pdelrho(n,1,2,ipert)
                pdrb(2) = pdelrho(n,2,2,ipert)
                pdrb(3) = pdelrho(n,3,2,ipert)
c
                cdelrho(n,1,1,ipert) = term_rg2a*pdra(1)
     1                               + term_rg3a*pdrb(1)
     2                               + term_gg1a*pdra(1)
     3                               + term_gg3a*pdrb(1)
     4                               + term_pdelrho1a*delrho(n,1,1)
     5                               + term_pdelrho1b*delrho(n,1,2)
c
                cdelrho(n,2,1,ipert) = term_rg2a*pdra(2)
     1                               + term_rg3a*pdrb(2)
     2                               + term_gg1a*pdra(2)
     3                               + term_gg3a*pdrb(2)
     4                               + term_pdelrho1a*delrho(n,2,1)
     5                               + term_pdelrho1b*delrho(n,2,2)
c
                cdelrho(n,3,1,ipert) = term_rg2a*pdra(3)
     1                               + term_rg3a*pdrb(3)
     2                               + term_gg1a*pdra(3)
     3                               + term_gg3a*pdrb(3)
     4                               + term_pdelrho1a*delrho(n,3,1)
     5                               + term_pdelrho1b*delrho(n,3,2)
c
                cdelrho(n,1,2,ipert) = term_rg2b*pdra(1)
     1                               + term_rg3b*pdrb(1)
     2                               + term_gg1b*pdra(1)
     3                               + term_gg3b*pdrb(1)
     4                               + term_pdelrho2a*delrho(n,1,1)
     5                               + term_pdelrho2b*delrho(n,1,2)
c
                cdelrho(n,2,2,ipert) = term_rg2b*pdra(2)
     1                               + term_rg3b*pdrb(2)
     2                               + term_gg1b*pdra(2)
     3                               + term_gg3b*pdrb(2)
     4                               + term_pdelrho2a*delrho(n,2,1)
     5                               + term_pdelrho2b*delrho(n,2,2)
c
                cdelrho(n,3,2,ipert) = term_rg2b*pdra(3)
     1                               + term_rg3b*pdrb(3)
     2                               + term_gg1b*pdra(3)
     3                               + term_gg3b*pdrb(3)
     4                               + term_pdelrho2a*delrho(n,3,1)
     5                               + term_pdelrho2b*delrho(n,3,2)
c
              else if (calctype.eq.2) then
c
c Perturbed density and density gradients
c
                  ptmp(1) = prho(n,1,ipert)
                  ptmp(2) = prho(n,2,ipert)
c
                  ptmp(3) = ( delrho(n,1,1)*pdelrho(n,1,1,ipert) 
     1                      + delrho(n,2,1)*pdelrho(n,2,1,ipert) 
     2                      + delrho(n,3,1)*pdelrho(n,3,1,ipert) )*2.0d0
c
                  ptmp(4) = delrho(n,1,1)*pdelrho(n,1,2,ipert) 
     1                    + delrho(n,2,1)*pdelrho(n,2,2,ipert) 
     2                    + delrho(n,3,1)*pdelrho(n,3,2,ipert) 
     3                    + delrho(n,1,2)*pdelrho(n,1,1,ipert) 
     4                    + delrho(n,2,2)*pdelrho(n,2,1,ipert) 
     5                    + delrho(n,3,2)*pdelrho(n,3,1,ipert)
c
                  ptmp(5) = ( delrho(n,1,2)*pdelrho(n,1,2,ipert) 
     1                      + delrho(n,2,2)*pdelrho(n,2,2,ipert) 
     2                      + delrho(n,3,2)*pdelrho(n,3,2,ipert) )*2.0d0
c
c Local perturbed density coefficients
c
                  crho(n,1,ipert) = Amat2(n,D2_RA_RA) *ptmp(1) 
     1                            + Amat2(n,D2_RA_RB) *ptmp(2) 
     2                            + Cmat2(n,D2_RA_GAA)*ptmp(3) 
     3                            + Cmat2(n,D2_RA_GAB)*ptmp(4) 
     4                            + Cmat2(n,D2_RA_GBB)*ptmp(5)
c
                  crho(n,2,ipert) = Amat2(n,D2_RA_RB) *ptmp(1) 
     1                            + Amat2(n,D2_RB_RB) *ptmp(2) 
     2                            + Cmat2(n,D2_RB_GAA)*ptmp(3) 
     3                            + Cmat2(n,D2_RB_GAB)*ptmp(4) 
     4                            + Cmat2(n,D2_RB_GBB)*ptmp(5)
c
c Gradient corrected CPKS coefficients
c
                  pdra(1) = pdelrho(n,1,1,ipert)
                  pdra(2) = pdelrho(n,2,1,ipert)
                  pdra(3) = pdelrho(n,3,1,ipert)
                  pdrb(1) = pdelrho(n,1,2,ipert)
                  pdrb(2) = pdelrho(n,2,2,ipert)
                  pdrb(3) = pdelrho(n,3,2,ipert)
c
c Perturbed density gradient contribution
c
                  cdelrho(n,1,1,ipert) = Cmat(n,D1_GAA)*pdra(1)*2.0d0 
     1                                 + Cmat(n,D1_GAB)*pdrb(1)
                  cdelrho(n,2,1,ipert) = Cmat(n,D1_GAA)*pdra(2)*2.0d0
     1                                 + Cmat(n,D1_GAB)*pdrb(2)
                  cdelrho(n,3,1,ipert) = Cmat(n,D1_GAA)*pdra(3)*2.0d0
     1                                 + Cmat(n,D1_GAB)*pdrb(3)
                  cdelrho(n,1,2,ipert) = Cmat(n,D1_GBB)*pdrb(1)*2.0d0
     1                                 + Cmat(n,D1_GAB)*pdra(1)
                  cdelrho(n,2,2,ipert) = Cmat(n,D1_GBB)*pdrb(2)*2.0d0
     1                                 + Cmat(n,D1_GAB)*pdra(2)
                  cdelrho(n,3,2,ipert) = Cmat(n,D1_GBB)*pdrb(3)*2.0d0
     1                                 + Cmat(n,D1_GAB)*pdra(3)
c
c Perturbed functional argument contribution
c
                  t(1) = ( Cmat2(n,D2_RA_GAA) *ptmp(1) 
     1                   + Cmat2(n,D2_RB_GAA) *ptmp(2) 
     2                   + Cmat2(n,D2_GAA_GAA)*ptmp(3) 
     3                   + Cmat2(n,D2_GAA_GAB)*ptmp(4) 
     4                   + Cmat2(n,D2_GAA_GBB)*ptmp(5) )*2.0d0
c
                  t(2) = Cmat2(n,D2_RA_GAB) *ptmp(1)
     1                 + Cmat2(n,D2_RB_GAB) *ptmp(2)
     2                 + Cmat2(n,D2_GAA_GAB)*ptmp(3) 
     3                 + Cmat2(n,D2_GAB_GAB)*ptmp(4) 
     4                 + Cmat2(n,D2_GAB_GBB)*ptmp(5)
c
                  cdelrho(n,1,1,ipert) = cdelrho(n,1,1,ipert) 
     1                                 + delrho(n,1,1)*t(1) 
     2                                 + delrho(n,1,2)*t(2)
c
                  cdelrho(n,2,1,ipert) = cdelrho(n,2,1,ipert) 
     1                                 + delrho(n,2,1)*t(1) 
     2                                 + delrho(n,2,2)*t(2)
c
                  cdelrho(n,3,1,ipert) = cdelrho(n,3,1,ipert) 
     1                                 + delrho(n,3,1)*t(1) 
     2                                 + delrho(n,3,2)*t(2)
c
                  t(1) = t(2)
                  t(2) = ( Cmat2(n,D2_RA_GBB) *ptmp(1) 
     1                   + Cmat2(n,D2_RB_GBB) *ptmp(2) 
     2                   + Cmat2(n,D2_GAA_GBB)*ptmp(3) 
     3                   + Cmat2(n,D2_GAB_GBB)*ptmp(4) 
     4                   + Cmat2(n,D2_GBB_GBB)*ptmp(5) )*2.0d0
c
                  cdelrho(n,1,2,ipert) = cdelrho(n,1,2,ipert) 
     1                                 + delrho(n,1,1)*t(1) 
     2                                 + delrho(n,1,2)*t(2)
c
                  cdelrho(n,2,2,ipert) = cdelrho(n,2,2,ipert) 
     1                                 + delrho(n,2,1)*t(1) 
     2                                 + delrho(n,2,2)*t(2)
c
                  cdelrho(n,3,2,ipert) = cdelrho(n,3,2,ipert) 
     1                                 + delrho(n,3,1)*t(1) 
     2                                 + delrho(n,3,2)*t(2)
              endif
            enddo
          endif 
        else if (triplet) then ! triplet case
c --------------------------------------------------------------
c Restricted Triplet case
c (XC-kernel gradients)
c --------------------------------------------------------------
          if (.not. grad) then ! local functionals
            do n = 1, nq
              if (calctype.eq.1) then
c Note that there is a cancellation of terms here (see comment at
c right).
                term_rrr = Amat3(n,D3_RA_RA_RA)                     !   rarara + rararb 
     1                   - Amat3(n,D3_RA_RB_RB)                     ! - rararb - rarbrb
c Here, we multiply by the perturbed density twice, mimicing what was
c done for the XC-kernel where the perturbed density is multiplied in
c once. 
                crho(n,1,ipert) =
     1            0.50d0*term_rrr*prho(n,1,ipert)*prho(n,1,ipert)
              else if (calctype.eq.2) then
                term_rr = Amat2(n,D2_RA_RA) - Amat2(n,D2_RA_RB)
                crho(n,1,ipert) = term_rr*prho(n,1,ipert)
              endif
            enddo
          else ! gradient dependent functionals
c            call errquit(
c     1        pname//'Triplet GGA calculations not implemented',
c     2        0, MA_ERR)
            do n = 1, nq
              if (calctype.eq.1) then
                ptmp(1) = prho(n,1,ipert)                           ! perturbed density
c
c Here we collect the various dot products needed to construct the 
c total functional derivative for a GGA.
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert)        ! delrho*del-perturbed density 
     1                  + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     2                  + delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                ptmp(3) = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! del-perturbed density*del-perturbed density
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c -------------------------------------------------------------------
c Third derivative contributions: drdrdr
c -------------------------------------------------------------------
                term_rrr = Amat3(n,D3_RA_RA_RA) - Amat3(n,D3_RA_RB_RB)
c -------------------------------------------------------------------
c Third derivative contributions: drdrdg
c -------------------------------------------------------------------
                term_rrg1 = 2.0d0*( 2.0d0*Cmat3(n,D3_RA_RA_GAA)
     1                            + Cmat3(n,D3_RA_RA_GAB)
     2                            - 2.0d0*Cmat3(n,D3_RA_RB_GBB)
     3                            - Cmat3(n,D3_RA_RB_GAB) )
c
                term_rrg2 = 2.0d0*Cmat3(n,D3_RA_RA_GAA)
     1                    + Cmat3(n,D3_RA_RA_GAB)
     2                    - 2.0d0*Cmat3(n,D3_RA_RB_GBB)
     3                    - Cmat3(n,D3_RA_RB_GAB)
c -------------------------------------------------------------------
c Third derivative contributions: drdgdg
c -------------------------------------------------------------------
                term_rgg1 = 4.0d0*( Cmat3(n,D3_RA_GAA_GAA)
     1                            + Cmat3(n,D3_RA_GAA_GAB)
     2                            - Cmat3(n,D3_RA_GBB_GBB)
     3                            - Cmat3(n,D3_RA_GAB_GBB) )
c
                term_rgg2 = 8.0d0*( Cmat3(n,D3_RA_GAA_GAA)
     1                            + Cmat3(n,D3_RA_GAA_GAB)
     2                            - Cmat3(n,D3_RA_GAA_GBB)
     3                            - Cmat3(n,D3_RA_GAB_GBB) )
c -------------------------------------------------------------------
c Third derivative contributions: dgdgdg
c -------------------------------------------------------------------
                term_ggg = 8.0d0*( Cmat3(n,D3_GAA_GAA_GAA)
     1                           + Cmat3(n,D3_GAA_GAA_GAB)
     3                           - Cmat3(n,D3_GAA_GAA_GBB)
     4                           - Cmat3(n,D3_GAA_GAB_GBB) )
c -------------------------------------------------------------------
c Second derivative contributions: drdg
c -------------------------------------------------------------------
                term_rg1 = 2.0d0*( Cmat2(n,D2_RA_GAA)
     1                           - Cmat2(n,D2_RA_GBB) )
c
                term_rg2 = 4.0d0*( Cmat2(n,D2_RA_GAA)
     1                           - Cmat2(n,D2_RA_GBB) )
c -------------------------------------------------------------------
c Second derivative contributions: dgdg
c -------------------------------------------------------------------
                term_gg1 = 8.0d0*( Cmat2(n,D2_GAA_GAA)
     1                           - Cmat2(n,D2_GAA_GBB) )
c
                term_gg2 = 4.0d0*( Cmat2(n,D2_GAA_GAA)
     1                           - Cmat2(n,D2_GAA_GBB) )
c -------------------------------------------------------------------
c Terms that do not involve the gradient with respect to the basis
c functions.  
c -------------------------------------------------------------------
c
c Recall that:
c  ptmp(1) = prho
c  ptmp(2) = pdelrho*delrho
c  ptmp(3) = pdelrho*pdelrho
c
                term_prho = term_rrr*ptmp(1)*ptmp(1)
     1                    + term_rrg1*ptmp(2)*ptmp(1)
     2                    + term_rgg1*ptmp(2)*ptmp(2)
     3                    + term_rg1*ptmp(3)
c -------------------------------------------------------------------
c Terms that involve the gradient with respect to the basis 
c functions.
c -------------------------------------------------------------------
                term_pdelrho = term_rrg2*ptmp(1)*ptmp(1)
     1                       + term_rgg2*ptmp(2)*ptmp(1)
     2                       + term_ggg*ptmp(2)*ptmp(2)
     3                       + term_gg2*ptmp(3)
c -------------------------------------------------------------------
                crho(n,1,ipert) = term_prho * 0.50d0
c
                cdelrho(n,1,1,ipert) =
     1              ( term_rg2*pdelrho(n,1,1,ipert)*ptmp(1)
     2              + term_gg1*pdelrho(n,1,1,ipert)*ptmp(2)
     3              + term_pdelrho*delrho(n,1,1) )*0.50d0
                cdelrho(n,2,1,ipert) =
     1              ( term_rg2*pdelrho(n,2,1,ipert)*ptmp(1)
     2              + term_gg1*pdelrho(n,2,1,ipert)*ptmp(2)
     2              + term_pdelrho*delrho(n,2,1) )*0.50d0
                cdelrho(n,3,1,ipert) =
     1              ( term_rg2*pdelrho(n,3,1,ipert)*ptmp(1)
     2              + term_gg1*pdelrho(n,3,1,ipert)*ptmp(2)
     3              + term_pdelrho*delrho(n,3,1) )*0.50d0
              else if (calctype.eq.2) then
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
                term_prho = term_rr*ptmp(1) + term_rg*ptmp(2)*2.0d0
                term_pdelrho = term_rg*ptmp(1) + term_gg*ptmp(2)*2.0d0
c
                crho(n,1,ipert) = term_prho
                cdelrho(n,1,1,ipert) = term_g*pdelrho(n,1,1,ipert) +
     &                                 term_pdelrho*delrho(n,1,1)*2.0d0
                cdelrho(n,2,1,ipert) = term_g*pdelrho(n,2,1,ipert) +
     &                                 term_pdelrho*delrho(n,2,1)*2.0d0
                cdelrho(n,3,1,ipert) = term_g*pdelrho(n,3,1,ipert) +
     &                                 term_pdelrho*delrho(n,3,1)*2.0d0
              endif
            enddo
          endif 
        else   ! singlet case
c --------------------------------------------------------------
c Restricted case
c (XC-kernel gradients)
c --------------------------------------------------------------
          if (.not. grad) then ! local functionals
            do n = 1, nq
              if (calctype.eq.1) then
                term_rrr = Amat3(n,D3_RA_RA_RA)                  ! rarara + 2*rararb + rarbrb
     1                   + 2.0d0*Amat3(n,D3_RA_RA_RB)
     2                   + Amat3(n,D3_RA_RB_RB)
c Factor of 0.50 seems to imply there's an extra factor of 2 somewhere,
c or an inconsistency in the third derivatives with the remaining code.
c Oddly the third derivatives work correctly for generating the
c Z-vector, so I think the 0.50 is needed to make everything consistent.
                crho(n,1,ipert) = 
     1            0.50d0*term_rrr*prho(n,1,ipert)*prho(n,1,ipert)
              else if (calctype.eq.2) then
                term_rr = Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB)
                crho(n,1,ipert) =term_rr*prho(n,1,ipert)
              endif
            enddo
          else ! gradient dependent functionals
            do n = 1, nq
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Terms to be contracted with the nuclear gradient of the GS density.
c
              if (calctype.eq.1) then
c Here, we construct the two contributions to the gradient of fxc:
c
c 1. crho => part that only involves terms multiplied by up to the 
c   gradient of the basis functions.
c 2. cdelrho => part that involves terms multiplied by the Hessian 
c   or gradient dot product of the basis functions.
c
                ptmp(1) = prho(n,1,ipert)                           ! perturbed density
c
c Here we collect the various dot products needed to construct the 
c total functional derivative for a GGA.
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert)        ! delrho*del-perturbed density 
     1                  + delrho(n,2,1)*pdelrho(n,2,1,ipert)
     2                  + delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                ptmp(3) = pdelrho(n,1,1,ipert)*pdelrho(n,1,1,ipert) ! del-perturbed density*del-perturbed density
     1                  + pdelrho(n,2,1,ipert)*pdelrho(n,2,1,ipert)
     2                  + pdelrho(n,3,1,ipert)*pdelrho(n,3,1,ipert)
c
c ----------------------------------------------------------------------
c Third derivative contributions: drdrdr
c ----------------------------------------------------------------------
c Correct terms: term_rrr
                term_rrr = Amat3(n,D3_RA_RA_RA)
     1                   + 2.0d0*Amat3(n,D3_RA_RA_RB)
     2                   + Amat3(n,D3_RA_RB_RB)
c ----------------------------------------------------------------------
c Third derivative contributions: drdrdg
c ----------------------------------------------------------------------
c Correct terms: term_rrg1, term_rrg2, term_rrg3
                term_rrg1 = 2.0d0*( Cmat3(n,D3_RA_RA_GAA)
     1                            + 1.5d0*Cmat3(n,D3_RA_RA_GAB)
     2                            + 2.0d0*Cmat3(n,D3_RA_RA_GBB)
     3                            + Cmat3(n,D3_RA_RB_GBB)
     4                            + 0.5d0*Cmat3(n,D3_RA_RB_GAB))
c
                term_rrg2 = 2.0d0*( Cmat3(n,D3_RA_RA_GAA)
     1                            + 0.5d0*Cmat3(n,D3_RA_RA_GAB)
     2                            + 2.0d0*Cmat3(n,D3_RA_RB_GAA)
     3                            + 1.5d0*Cmat3(n,D3_RA_RB_GAB)
     4                            + Cmat3(n,D3_RA_RB_GBB))
c
                term_rrg3 = 2.0d0*( Cmat3(n,D3_RA_RA_GAA)
     1                            + Cmat3(n,D3_RA_RA_GAB)
     2                            + 2.0d0*Cmat3(n,D3_RA_RB_GAA)
     3                            + Cmat3(n,D3_RA_RB_GAB)
     4                            + Cmat3(n,D3_RA_RA_GBB))
c ----------------------------------------------------------------------
c Third derivative contributions: drdgdg
c ----------------------------------------------------------------------
c Correct terms: term_rgg1, term_rgg2, term_rgg3
                term_rgg1 = 4.0d0*( Cmat3(n,D3_RA_GAA_GAA)
     1                            + 2.0d0*Cmat3(n,D3_RA_GAA_GAB)
     2                            + Cmat3(n,D3_RA_GAB_GAB)
     3                            + 2.0d0*Cmat3(n,D3_RA_GAA_GBB)
     4                            + 2.0d0*Cmat3(n,D3_RA_GAB_GBB)
     5                            + Cmat3(n,D3_RA_GBB_GBB))
c
                term_rgg2 = 4.0d0*( Cmat3(n,D3_RA_GAA_GAA)
     1                            + 2.0d0*Cmat3(n,D3_RA_GAA_GAB)
     2                            + Cmat3(n,D3_RA_GAB_GAB)
     3                            + 3.0d0*Cmat3(n,D3_RA_GAA_GBB)
     4                            + 2.0d0*Cmat3(n,D3_RA_GAB_GBB))
c
                term_rgg3 = 4.0d0*( Cmat3(n,D3_RA_GAA_GAA)
     1                            + Cmat3(n,D3_RA_GAA_GAB)
     2                            + Cmat3(n,D3_RA_GAB_GAB)
     3                            + 2.0d0*Cmat3(n,D3_RA_GBB_GBB)
     4                            + 3.0d0*Cmat3(n,D3_RA_GAB_GBB)
     5                            + Cmat3(n,D3_RA_GAA_GBB))
c ----------------------------------------------------------------------
c Third derivative contributions: dgdgdg
c ----------------------------------------------------------------------
c Correct terms: term_ggg
                term_ggg = 4.0d0*( 2.0d0*Cmat3(n,D3_GAA_GAA_GAA)
     1                           + 6.0d0*Cmat3(n,D3_GAA_GAA_GAB)
     2                           + 6.0d0*Cmat3(n,D3_GAA_GAB_GAB)
     3                           + Cmat3(n,D3_GAB_GAB_GAB)
     4                           + 4.0d0*Cmat3(n,D3_GAA_GAA_GBB)
     5                           + 6.0d0*Cmat3(n,D3_GAA_GAB_GBB)
     6                           + 2.0d0*Cmat3(n,D3_GAA_GBB_GBB))
c ----------------------------------------------------------------------
c Second derivative contributions: drdg
c ----------------------------------------------------------------------
c Correct terms: term_rg1, term_rg2, term_rg3
                term_rg1 = 2.0d0*( Cmat2(n,D2_RA_GAA)
     1                           + Cmat2(n,D2_RA_GAB)
     2                           + Cmat2(n,D2_RA_GBB))
c
                term_rg2 = 2.0d0*( Cmat2(n,D2_RA_GAA)
     1                           + 1.5d0*Cmat2(n,D2_RA_GAB))
c
                term_rg3 = 2.0d0*( Cmat2(n,D2_RA_GAA)
     1                           + 0.5d0*Cmat2(n,D2_RA_GAB)
     2                           + 2.0d0*Cmat2(n,D2_RA_GBB))
c ----------------------------------------------------------------------
c Second derivative contributions: dgdg
c ----------------------------------------------------------------------
c Correct terms: term_gg1, term_gg2, term_gg3
                term_gg1 = 4.0d0*Cmat2(n,D2_GAA_GAA)
     1                   + 4.0d0*Cmat2(n,D2_GAA_GBB)
     2                   + 8.0d0*Cmat2(n,D2_GAA_GAB)
     3                   + 2.0d0*Cmat2(n,D2_GAB_GAB)
c
                term_gg2 = 4.0d0*Cmat2(n,D2_GAA_GAA)
     1                   + 8.0d0*Cmat2(n,D2_GAA_GAB)
     2                   + 3.0d0*Cmat2(n,D2_GAB_GAB)
c
                term_gg3 = 4.0d0*Cmat2(n,D2_GAA_GAA)
     1                   + 8.0d0*Cmat2(n,D2_GAA_GBB)
     2                   + 8.0d0*Cmat2(n,D2_GAA_GAB)
     3                   + Cmat2(n,D2_GAB_GAB)
c ----------------------------------------------------------------------
c Terms that do not involve the gradient with respect to the basis
c functions.  Recall that:
c  ptmp(1) = prho
c  ptmp(2) = pdelrho*delrho
c  ptmp(3) = pdelrho*pdelrho
c
                term_prho = term_rrr*ptmp(1)*ptmp(1)
     1                    + (term_rrg1 + term_rrg2)*ptmp(2)*ptmp(1)
     2                    + term_rgg1*ptmp(2)*ptmp(2)
     3                    + term_rg1*ptmp(3)
c ----------------------------------------------------------------------
c Terms that involve the gradient with respect to the basis 
c functions.
c ----------------------------------------------------------------------
                term_pdelrho = term_rrg3*ptmp(1)*ptmp(1)
     1                       + (term_rgg2 + term_rgg3)*ptmp(2)*ptmp(1)
     2                       + term_ggg*ptmp(2)*ptmp(2)
     3                       + term_gg1*ptmp(3)
c ----------------------------------------------------------------------
c Above we needed a factor of 0.5 to get the correct contribution for
c LDA, so I assume we need it here also...
                crho(n,1,ipert) = term_prho*0.50d0
                cdelrho(n,1,1,ipert) =
     1              ( (term_rg2 + term_rg3)*pdelrho(n,1,1,ipert)*ptmp(1)
     2              + (term_gg2 + term_gg3)*pdelrho(n,1,1,ipert)*ptmp(2)
     3              + term_pdelrho*delrho(n,1,1) )*0.50d0
                cdelrho(n,2,1,ipert) =
     1              ( (term_rg2 + term_rg3)*pdelrho(n,2,1,ipert)*ptmp(1)
     2              + (term_gg2 + term_gg3)*pdelrho(n,2,1,ipert)*ptmp(2)
     3              + term_pdelrho*delrho(n,2,1) )*0.50d0
                cdelrho(n,3,1,ipert) =
     1              ( (term_rg2 + term_rg3)*pdelrho(n,3,1,ipert)*ptmp(1)
     2              + (term_gg2 + term_gg3)*pdelrho(n,3,1,ipert)*ptmp(2)
     3              + term_pdelrho*delrho(n,3,1) )*0.50d0
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Terms to be contracted with the nuclear gradient of the perturbed
c density.
c
              else if (calctype.eq.2) then
                ptmp(1) = prho(n,1,ipert)                          ! perturbed density
                ptmp(2) = delrho(n,1,1)*pdelrho(n,1,1,ipert) +     ! delrho*perturbed density 
     1                    delrho(n,2,1)*pdelrho(n,2,1,ipert) +
     2                    delrho(n,3,1)*pdelrho(n,3,1,ipert)
c
                term_g  = Cmat(n,D1_GAA)*2.0d0 + Cmat(n,D1_GAB)    ! 2*gaa + gab
                term_rr = Amat2(n,D2_RA_RA)    + Amat2(n,D2_RA_RB) ! rara + rarb
                term_rg = Cmat2(n,D2_RA_GAA) +                     ! ragaa + ragbb + ragab
     1                    Cmat2(n,D2_RA_GBB) +
     2                    Cmat2(n,D2_RA_GAB)
                term_gg = Cmat2(n,D2_GAA_GAA) +                    ! gaagaa + gaagbb + 2*gaagab + 0.5*gabgab
     1                    Cmat2(n,D2_GAA_GBB) +
     2                    Cmat2(n,D2_GAA_GAB)*2.0d0 +
     3                    Cmat2(n,D2_GAB_GAB)*0.5d0
c
                term_prho = term_rr*ptmp(1) + term_rg*ptmp(2)*2.0d0
                term_pdelrho = term_rg*ptmp(1) + term_gg*ptmp(2)*2.0d0
c
                crho(n,1,ipert) = term_prho
c
                cdelrho(n,1,1,ipert) = term_g*pdelrho(n,1,1,ipert)
     1                               + term_pdelrho*delrho(n,1,1)*2.0d0
                cdelrho(n,2,1,ipert) = term_g*pdelrho(n,2,1,ipert)
     1                               + term_pdelrho*delrho(n,2,1)*2.0d0
                cdelrho(n,3,1,ipert) = term_g*pdelrho(n,3,1,ipert)
     1                               + term_pdelrho*delrho(n,3,1)*2.0d0
              endif
            enddo
          endif
        endif
      enddo
c
c     Put delrho back the way it was since it may be used later on
      if (ipol.eq.1) then
        call dscal(nq*ipol*npert,2.0d0,prho,1)
        if (grad) then
          call dscal(nq*3*ipol*npert,2.0d0,pdelrho,1)
          call dscal(nq*3*ipol,2.0d0,delrho,1)
        endif
      endif
c ------
c Return
c ------
      return
      end
c $Id$
