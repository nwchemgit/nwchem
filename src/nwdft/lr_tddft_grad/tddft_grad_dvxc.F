      subroutine tddft_grad_dvxc(prho, pdelrho,
     1  Amat2, Cmat2, Amat, Cmat,  
     2  delrho, npert, ipol, nq, grad,
     3  crho, cdelrho, calctype)
c
c     $Id$
c
c Routine for building the gradient of the XC-potential, based on the
c recipe outlined in J. Chem. Phys., 124, 094107 (2006).
c
c Called from: dftg_gridv0b
c
      implicit none
#include "rtdb.fh"
#include "dft2drv.fh"
#include "stdio.fh"
#include "errquit.fh"
c
      integer npert    ! Number of density matrices
      integer ipol     ! 1 = restricted, 2 = unrestricted
      integer nq       ! Number of quadrature points
      integer calctype ! 1 = XC-coefficient to be multiplied with GS density gradient
                       ! 2 = XC-coefficient to be multiplied with perturbed density gradient
c
      logical grad     ! True if a GGA is used 
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
c Second derivatives of XC functional
c
      double precision Amat2(nq,NCOL_AMAT2)
      double precision Cmat2(nq,NCOL_CMAT2)
c
c First derivatives of XC functional
c
c      double precision Amat(nq,*)
c      double precision Cmat(nq,*)
      double precision Amat(nq,ipol)
      double precision Cmat(nq,3)
c
c Gradients of the GS spin densities 
c
      double precision delrho(nq,3,ipol)
c
      integer ipert, n ! Loop indices
c
      double precision ptmp(5)
      double precision term_r,term_g,term_rr,term_rg,term_gg
      double precision term_prho, term_pdelrho
c
c Initialize ptmp
c
      do n = 1, 5
        ptmp(n) = 0.d0
      end do
c
      if (ipol.eq.1) then
c
c Since the total densities are evaluated in the restricted case,
c scale them so that the correct CPKS matrices will be produced.
c
        call dscal(nq*ipol*npert,0.25d0,prho,1)
        if (grad) then
          call dscal(nq*3*ipol*npert,0.25d0,pdelrho,1)
          call dscal(nq*3*ipol,0.25d0,delrho,1)
        endif
      endif

      do ipert = 1, npert
        if (ipol.eq.2) then  ! unrestricted
          call errquit(
     1      'tddftg_dvxc: Cannot do unrestricted TDDFT gradients',
     2      0, CAPMIS_ERR)
c --------------------------------------------------------------
c Unrestricted case
c (XC-potential gradients)
c --------------------------------------------------------------
        else ! restricted 
c --------------------------------------------------------------
c Restricted case
c (TDDFT XC-potential gradients)
c --------------------------------------------------------------
c Here, we assemble/reassign the matrices without overwriting the 
c perturbed density.  This is a consequence of the gradient of the 
c XC-kernel requiring the perturbed density matrix, so this routine
c is written to parallel the XC-kernel gradient construction.
          if (.not. grad) then                                   ! local functionals
            do n = 1, nq
              if (calctype.eq.1) then 
                term_rr = Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB)    ! rara + rarb
                crho(n,1,ipert) = term_rr*prho(n,1,ipert)
c                crho(n,1,ipert) = 0.0d0
              else if (calctype.eq.2) then
                term_r = 0.50d0*Amat(n,1)
                crho(n,1,ipert) = term_r
c                crho(n,1,ipert) = 0.0d0
              endif
            enddo
          else                                                   ! gradient dependent functionals
            do n = 1, nq
              if (calctype.eq.1) then
c Here, we construct the two contributions to the gradient of Vxc:
c
c 1. crho => part that only involves terms multiplied by up to the 
c   gradient of the basis functions.
c 2. cdelrho => part that involves terms multiplied by the Hessian 
c   or gradient dot product of the basis functions.
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
              else if (calctype.eq.2) then
                term_r = Amat(n,1)                                 ! ra
                term_g = 2.0d0*Cmat(n,D1_GAA) + Cmat(n,D1_GAB)     ! 2*gaa + gab
c
                crho(n,1,ipert) = term_r
c
                cdelrho(n,1,1,ipert) = term_g*delrho(n,1,1)
                cdelrho(n,2,1,ipert) = term_g*delrho(n,2,1)
                cdelrho(n,3,1,ipert) = term_g*delrho(n,3,1)
c
c                write(6,*) 'Cmat_GAA', Cmat(n,D1_GAA)
c                write(6,*) 'term_r', term_r                  ! correct
c                write(6,*) 'term_g', term_g                  ! incorrect ?
c                write(6,*) 'delrho_x', delrho(n,1,1)         ! correct
c                write(6,*) 'delrho_y', delrho(n,2,1)         ! correct
c                write(6,*) 'delrho_z', delrho(n,3,1)         ! correct
c                write(6,*) 'pdelrho_x', pdelrho(n,1,1,ipert)
c                write(6,*) 'pdelrho_y', pdelrho(n,2,1,ipert)
c                write(6,*) 'pdelrho_z', pdelrho(n,3,1,ipert)
c                write(6,*) 'crho', crho(n,1,ipert)           ! correct
c                write(6,*) 'cdelrho_x', cdelrho(n,1,1,ipert) ! incorrect ?
c                write(6,*) 'cdelrho_y', cdelrho(n,2,1,ipert) ! incorrect ?
c                write(6,*) 'cdelrho_z', cdelrho(n,3,1,ipert) ! incorrect ?
              endif
            enddo
          endif
        endif
      enddo
c
c Put densities back the way they were
c
      if (ipol.eq.1) then
        call dscal(nq*ipol*npert,4.0d0,prho,1)
        if (grad) then
          call dscal(nq*3*ipol*npert,4.0d0,pdelrho,1)
          call dscal(nq*3*ipol,4.0d0,delrho,1)
        endif
      endif
c ------
c Return
c ------
      return
      end
