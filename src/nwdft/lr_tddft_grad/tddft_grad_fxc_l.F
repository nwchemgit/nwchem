c
c     Compute the TDDFT contributions to the gradient. This routine
c     considers one excited state at a time as the final gradient 
c     evaluator tddft_grad_1e loops over excited states. This could
c     be changed for efficiency reasons but that requires quite some
c     code adjustments in the DFT drivers.
c
      subroutine tddft_grad_fxc_l(ldew,ipol,nq,natoms,oactive,
     &                            nactive,ictr,
     &                            rho,drho,amat2,amat3,force,exc)
      implicit none
c
c     $Id$
c
c     The general expression to evaluate 
c
c       dE   dE   d d2f*(X+Y)*(X+Y)
c       -- = -- + -----------------
c       dl   dl        dl
c
c     The factors in this expression can be classified into two main
c     groups:
c
c       1. density like factors
c       2. derivatives of density like factors wrt. nuclear coordinates
c
c     There is of course a tricky bit involved with the derivatives of
c     the weights. I.e. the "energy" factor that gets multiplied with
c     the derivatives of the weights needs to be computed as well.
c     
c     Written by Huub van Dam, November 2006.
c
#include "dft2drv.fh"
#include "dft3drv.fh"
c
c     Input:
c
      logical ldew   ! .true. if derivatives of weights required
      integer ipol   ! =1 (restricted), =2 (unrestricted)
      integer nq     ! the number of quadrature points in this batch
      integer natoms ! the number of atoms
      integer ictr   ! the centre on which the current batch of grid
                     ! points reside
c
      integer nactive          ! the number of "active" atoms
      logical oactive(nactive) ! .true. if atom active in current batch
c
c     density factors rho( *, *, density/(X+Y)_ao )
      double precision rho(nq,ipol,2)
c
c     derivatives of the density wrt. nuclear coordinates
c     drho( *, *, x/y/z, density/(X+Y)_ao )
      double precision drho(nq,ipol,3*nactive,2) 
c
c     second order functional derivatives wrt. rho
      double precision amat2(nq,NCOL_AMAT2) 
c
c     third order functional derivatives wrt. rho
      double precision amat3(nq,NCOL_AMAT3) 
c
c     Input/Output:
c
      double precision exc(nq) ! "energy" terms needed for derivatives
                               ! of weights contribution
c
c     Output:
c
      double precision force(3,natoms) ! the TDDFT force contribution
c
c     Local:
c
      integer RA, RB, XPY, D0
      parameter (RA  = 1)  ! Alpha density component
      parameter (RB  = 2)  ! Beta  density component
      parameter (D0  = 1)  ! Ground state density component
      parameter (XPY = 2)  ! (X+Y) density component
c
      logical ldo  ! do we need to do this atom
      integer iq   ! counter over grid points
      integer ix   ! counter over atomic coordinates x,y,z
      integer ic   ! counter over atomic coordinates A(x,y,z),B(x,y,z)..
      integer iat  ! counter over atoms
      integer icur ! the current atom in the compressed list of atoms
c
      double precision out(3) ! buffer for contributions to atom
c
c     Code:
c
c     First do the gradient part
c
cDEBUG
#include "global.fh"
      logical oroot
      oroot = ga_nodeid().eq.0
      if (oroot) then
c     write(*,*)'*** fxc_l loc(rho(XPY))  = ',loc(rho(1,RA,XPY))
c     write(*,*)'*** fxc_l loc(drho(D0))  = ',loc(drho(1,1,1,D0))
c     write(*,*)'*** fxc_l loc(drho(XPY)) = ',loc(drho(1,1,1,XPY))
      endif
cDEBUG
      icur = 0
      if (ipol.eq.1) then
c
c       Note: in the closed shell case RA refers to the total density
c
        do iat = 1, natoms
          if (oactive(iat)) then
            icur = icur + 1
            ldo  = (.not.ldew).or.(iat.ne.ictr)
            if (ldo) then
              do ix = 1, 3
                out(ix) = 0.0d0
              enddo
cDEBUG
c             do iq = 1, nq
c               if (oroot) write(*,150)iq,rho(iq,RA,D0)
c             enddo
c150          format(i4,' rho(D0)   ',f30.10)
c
c             do iq = 1, nq
c               if (oroot) write(*,160)iq,rho(iq,RA,XPY)
c             enddo
c160          format(i4,' rho(XPY)  ',f30.10)
c
c             do iq = 1, nq
c               do ix = 1, 3
c                 ic = 3*(icur-1)+ix
c                 if (oroot) write(*,100)iq,ic,drho(iq,RA,ic,D0)
c               enddo
c             enddo
c100          format(2i4,' drho(D0)',f30.10)
c
c             do iq = 1, nq
c               do ix = 1, 3
c                 ic = 3*(icur-1)+ix
c                 if (oroot) write(*,110)iq,ic,drho(iq,RA,ic,XPY)
c               enddo
c             enddo
c110          format(2i4,' drho(XPY)',f30.10)
c
c             do iq = 1, nq
c               if (oroot) write(*,120)iq,amat2(iq,D2_RA_RA)
c             enddo
c120          format(i4,' amat2(A2)',f30.10)
c
c             do iq = 1, nq
c               if (oroot) write(*,130)iq,amat3(iq,D3_RA_RA_RA)
c             enddo
c130          format(i4,' amat3(A3)',f30.10)
c
c             do iq = 1, nq
c               if (oroot) write(*,140)iq,amat3(iq,D3_RA_RA_RB)
c             enddo
c140          format(i4,' amat3(A2B)',f30.10)
cDEBUG
              do iq = 1, nq
                do ix = 1, 3
c
c                 Do Fijs,kls first
c
                  ic = 3*(icur-1)+ix
                  out(ix) = out(ix) 
     &            + 0.5d0*amat2(iq,D2_RA_RA)*rho(iq,RA,XPY)
     &              *drho(iq,RA,ic,XPY)
chvd &            + 0.125d0*amat3(iq,D3_RA_RA_RA)*rho(iq,RA,XPY)**2
     &            + 0.0625d0*amat3(iq,D3_RA_RA_RA)*rho(iq,RA,XPY)**2
     &              *drho(iq,RA,ic,D0)
     &            + 0.125d0*amat3(iq,D3_RA_RA_RB)*rho(iq,RA,XPY)**2
chvd &            + 0.0625d0*amat3(iq,D3_RA_RA_RB)*rho(iq,RA,XPY)**2
     &              *drho(iq,RA,ic,D0)
c
c                 Do Fijs,klt next
c
                  out(ix) = out(ix) 
     &            + 0.25d0*amat3(iq,D3_RA_RA_RB)*rho(iq,RA,XPY)**2
     &              *drho(iq,RA,ic,D0)
     &            + 0.5d0*amat2(iq,D2_RA_RB)*rho(iq,RA,XPY)
     &              *drho(iq,RA,ic,XPY)
                enddo
              enddo
              do ix = 1, 3
                force(ix,iat) = force(ix,iat) + out(ix)
              enddo
              if (ldew.and.oactive(ictr)) then
                do ix = 1, 3
                  force(ix,ictr) = force(ix,ictr) - out(ix)
                enddo
              endif
            endif
          endif
        enddo
      else
        do iat = 1, natoms
          if (oactive(iat)) then
            icur = icur + 1
            ldo  = (.not.ldew).or.(iat.ne.ictr)
            if (ldo) then
              do ix = 1, 3
                out(ix) = 0.0d0
              enddo
              do iq = 1, nq
                do ix = 1, 3
c
c                 First the Fija,kla terms
c
                  ic = 3*(icur-1)+ix
                  out(ix) = out(ix) 
     &            + 2*amat2(iq,D2_RA_RA)*rho(iq,RA,XPY)
     &              *drho(iq,RA,ic,XPY)
     &            + amat3(iq,D3_RA_RA_RA)*rho(iq,RA,XPY)**2
     &              *drho(iq,RA,ic,D0)
     &            + amat3(iq,D3_RA_RA_RB)*rho(iq,RA,XPY)**2
     &              *drho(iq,RB,ic,D0)
c
c                 Then the Fijb,klb terms
c
                  out(ix) = out(ix) 
     &            + 2*amat2(iq,D2_RB_RB)*rho(iq,RB,XPY)
     &              *drho(iq,RB,ic,XPY)
     &            + amat3(iq,D3_RA_RB_RB)*rho(iq,RB,XPY)**2
     &              *drho(iq,RA,ic,D0)
     &            + amat3(iq,D3_RB_RB_RB)*rho(iq,RB,XPY)**2
     &              *drho(iq,RB,ic,D0)
c
c                 Then the Fija,klb terms
c
                  out(ix) = out(ix) 
     &            + amat3(iq,D3_RA_RA_RB)*rho(iq,RA,XPY)*rho(iq,RB,XPY)
     &              *drho(iq,RA,ic,D0)
     &            + amat3(iq,D3_RA_RB_RB)*rho(iq,RA,XPY)*rho(iq,RB,XPY)
     &              *drho(iq,RB,ic,D0)
     &            + amat2(iq,D2_RA_RB)*rho(iq,RB,XPY)*drho(iq,RA,ic,XPY)
     &            + amat2(iq,D2_RA_RB)*rho(iq,RA,XPY)*drho(iq,RB,ic,XPY)
c
c                 Last the Fijb,kla terms
c
                  out(ix) = out(ix) 
     &            + amat3(iq,D3_RA_RB_RB)*rho(iq,RB,XPY)*rho(iq,RA,XPY)
     &              *drho(iq,RB,ic,D0)
     &            + amat3(iq,D3_RA_RA_RB)*rho(iq,RB,XPY)*rho(iq,RA,XPY)
     &              *drho(iq,RA,ic,D0)
     &            + amat2(iq,D2_RA_RB)*rho(iq,RA,XPY)*drho(iq,RB,ic,XPY)
     &            + amat2(iq,D2_RA_RB)*rho(iq,RB,XPY)*drho(iq,RA,ic,XPY)
                enddo
              enddo
              do ix = 1, 3
                force(ix,iat) = force(ix,iat) + out(ix)
              enddo
              if (ldew.and.oactive(ictr)) then
                do ix = 1, 3
                  force(ix,ictr) = force(ix,ictr) - out(ix)
                enddo
              endif
            endif
          endif
        enddo
      endif
c
c     Second do the "energy" part that is to be multiplied with the
c     derivatives of the weights
c
      if (ldew) then
        if (ipol.eq.1) then
          do iq = 1, nq
            exc(iq) = 
     &        0.25d0*amat2(iq,D2_RA_RA)*rho(iq,RA,XPY)**2
     &      + 0.25d0*amat2(iq,D2_RA_RB)*rho(iq,RA,XPY)**2
c           exc(iq) = 2*exc(iq)
cDEBUG
c           exc(iq) = 0.d0
cDEBUG
          enddo
        else
          do iq = 1, nq
            exc(iq) = 
     &          amat2(iq,D2_RA_RA)*rho(iq,RA,XPY)**2
     &      + 2*amat2(iq,D2_RA_RB)*rho(iq,RA,XPY)*rho(iq,RB,XPY)
     &      +   amat2(iq,D2_RB_RB)*rho(iq,RB,XPY)**2
          enddo
        endif
      endif
      end
