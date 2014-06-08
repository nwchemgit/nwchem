      subroutine tddft_grad_gxc(ipol,nroot,nbf,ibf,nq,chi,dchi,xpy_rho,
     &                          amat3,scr,scrmat)
      implicit none
c
c     $Id$
c
c     Compute the DFT contributions to R and W in AO basis. This routine
c     will be called from within the DFT drivers, so it has to be
c     adapted to the idiosyncracies thereof. 
c
c     The general expression to evaluate is
c
c        Rijs = Rijs + Sum_klt Sum_mnu g_ijskltmnu*(X+Y)_klt*(X+Y)_mnu
c
c     where i,j,k,l,m and n are AO labels and s,t and u are spin labels.
c     The quantity g is the third derivative of the functional with
c     respect to the density matrix, i.e.:
c
c                      d3 f(ra,rb,gaa,gab,gbb,ta,tb)
c        g_ijskltmnu = -----------------------------
c                      d D_ijs * d D_klt * d D_mnu
c
c     Grouping factors in the expression at the top it can be reduced
c     to an equation involving 3 different kinds of factors:
c
c       1. partial derivatives of the functional
c       2. "basis function" factors
c       3. "density" (like) factors
c
c     The first kind of factors are obvious, the basis function factors
c     are:
c
c       X_i(r)*X_j(r)                                          (101)
c       ... more for (meta) GGAs ...
c
c     The density factors are:
c
c       Sum_ij X_i(r)*X_j(r)*(X+Y)_ijs                         (201)
c       ... more for (meta) GGAs ...
c
c     In particular computing the density factors up front helps 
c     eliminating indeces early and is therefore recommended to 
c     enhance the efficiency. As these expressions are identical to
c     the normal density expressions the same machinery can be used
c     to evaluate them. Note that typical density evaluation routines
c     will assume that the density presented to them is symmetric. 
c     Due to the form of the expressions this can be dealt with by
c     symmetrizing (X+Y) before passing it into the density evaluators.
c
c     Written by Huub van Dam, November 2006.
c
#include "dft3drv.fh"
#include "global.fh"
c
c     Input:
c
      integer ipol     ! =1 (restricted), =2 (unrestricted)
      integer nroot    ! the number of roots
      integer nbf      ! the number of active basis functions
      integer ibf(nbf) ! the table giving the true basis function number
                       ! for every element of compressed list
      integer nq       ! the number of grid points
c
      double precision chi(nq,nbf)            ! basis function values
      double precision dchi(nq,3,nbf)         ! basis function gradients
      double precision xpy_rho(nq,ipol,nroot) ! density factor eq. (201)
c
c     3rd order functional derivatives wrt rho
      double precision amat3(nq,NCOL_AMAT3)
c                                       
c     Input/Output:
c
c     buffer to accumulate R matrix contributions into:
c
      double precision scrmat(ipol*nroot*nbf*nbf) 
c
c     integer g_r ! the global array holding the TDDFT gradient
c                 ! right-hand-side. The contribitions are added onto
c                 ! this.
c
c     Workspace:
c
      double precision scr(nq) ! buffer to minimize multiplications
c
c     Local:
c
      integer kbf, lbf ! counters over basis functions
      integer iq       ! counter over quadrature points
      integer ir       ! counter over roots
c
      double precision t ! a term
c
c     Code:
c
cDEBUG
      logical oroot
      oroot = ga_nodeid().eq.0
c     if (oroot) write(*,*)'*** tddft_grad_gxc begin'
c     if (oroot) write(*,*)'*** basis functions'
c     do kbf=1,nbf
c       do iq=1,nq
c         if (oroot) write(*,*)'*** chi=',kbf,iq,chi(iq,kbf)
c       enddo
c     enddo
c     if (oroot) write(*,*)'*** amat3'
c     do iq=1,nq
c       if (oroot) write(*,*)'*** amat3=',amat3(iq,D3_RA_RA_RA)
c     enddo
c     if (oroot) write(*,*)'*** xpy_rho'
c     do ir=1,nroot
c       do iq=1,nq
c         if (oroot) write(*,*)'*** xpy_rho=',ir,iq,xpy_rho(iq,1,ir)
c       enddo
c     enddo
c     if (oroot) write(*,*)'*** tddft_grad_gxc end'
cDEBUG
      if (ipol.eq.1) then
c
c       Only one spin component to do
c
        do ir = 1, nroot
          do kbf = 1, nbf
            do iq = 1, nq
              scr(iq) = 
     &            chi(iq,kbf)*amat3(iq,D3_RA_RA_RA)* ! gaaa
     &            xpy_rho(iq,1,ir)*xpy_rho(iq,1,ir)
     &        + 2*chi(iq,kbf)*amat3(iq,D3_RA_RA_RB)* ! gaab+gaba
     &            xpy_rho(iq,1,ir)*xpy_rho(iq,1,ir)
     &        +   chi(iq,kbf)*amat3(iq,D3_RA_RA_RB)* ! gabb
     &            xpy_rho(iq,1,ir)*xpy_rho(iq,1,ir)
              scr(iq) = 0.25d0*scr(iq) ! xpy_rho: total rho -> alpha rho
              scr(iq) = 2.0d0*scr(iq)  ! 2*gxc
cDEBUG
c     scr(iq) = 0.d0
cDEBUG
            enddo
            do lbf = 1, kbf-1
              t = 0.0d0
              do iq = 1, nq
                t = t + chi(iq,lbf)*scr(iq)
              enddo
              scrmat(ir+(kbf-1)*nroot+(lbf-1)*nroot*nbf) = t
              scrmat(ir+(lbf-1)*nroot+(kbf-1)*nroot*nbf) = t
            enddo
            t = 0.0d0
            do iq = 1, nq
              t = t + chi(iq,kbf)*scr(iq)
            enddo
            scrmat(ir+(kbf-1)*nroot+(kbf-1)*nroot*nbf) = t
          enddo
        enddo
c
      else if (ipol.eq.2) then
c
c       Alpa density component first
c
        do ir = 1, nroot
          do kbf = 1, nbf
            do iq = 1, nq
              scr(iq) = chi(iq,kbf)*(
     &                    amat3(iq,D3_RA_RA_RA)*
     &                    xpy_rho(iq,1,ir)*xpy_rho(iq,1,ir)
     &                 +2*amat3(iq,D3_RA_RA_RB)*
     &                    xpy_rho(iq,1,ir)*xpy_rho(iq,2,ir)
     &                 +  amat3(iq,D3_RA_RB_RB)*
     &                    xpy_rho(iq,2,ir)*xpy_rho(iq,2,ir)
     &                  )
            enddo
            do lbf = 1, kbf-1
              t = 0.0d0
              do iq = 1, nq
                t = t + chi(iq,lbf)*scr(iq)
              enddo
              scrmat(ir+(lbf-1)*nroot*ipol+(kbf-1)*nroot*ipol*nbf) = t
              scrmat(ir+(kbf-1)*nroot*ipol+(lbf-1)*nroot*ipol*nbf) = t
            enddo
            t = 0.0d0
            do iq = 1, nq
              t = t + chi(iq,kbf)*scr(iq)
            enddo
            scrmat(ir+(kbf-1)*nroot*ipol+(kbf-1)*nroot*ipol*nbf) = t
          enddo
        enddo
c
c       Beta density component next
c
        do ir = 1, nroot
          do kbf = 1, nbf
            do iq = 1, nq
              scr(iq) = chi(iq,kbf)*(
     &                    amat3(iq,D3_RA_RA_RB)*
     &                    xpy_rho(iq,1,ir)*xpy_rho(iq,1,ir)
     &                 +2*amat3(iq,D3_RA_RB_RB)*
     &                    xpy_rho(iq,1,ir)*xpy_rho(iq,2,ir)
     &                 +  amat3(iq,D3_RB_RB_RB)*
     &                    xpy_rho(iq,2,ir)*xpy_rho(iq,2,ir)
     &                  )
            enddo
            do lbf = 1, kbf-1
              t = 0.0d0
              do iq = 1, nq
                t = t + chi(iq,lbf)*scr(iq)
              enddo
              scrmat(ir+(ipol-1)*nroot+(lbf-1)*nroot*ipol
     &               +(kbf-1)*nroot*ipol*nbf) = t
              scrmat(ir+(ipol-1)*nroot+(kbf-1)*nroot*ipol
     &               +(lbf-1)*nroot*ipol*nbf) = t
            enddo
            t = 0.0d0
            do iq = 1, nq
              t = t + chi(iq,kbf)*scr(iq)
            enddo
            scrmat(ir+(ipol-1)*nroot+(kbf-1)*nroot*ipol
     &             +(kbf-1)*nroot*ipol*nbf) = t
          enddo
        enddo
      endif
c
      end
