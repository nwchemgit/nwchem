c
c     Precomputes coefficients of basis functions pairs for use in
c     making multiple XC matrices for CPKS
c
c     BGJ - 8/97
c
c     $Id: xc_cpks_coeff.F,v 1.1 1998-09-09 13:21:12 bjohnson Exp $
c
      Subroutine xc_cpks_coeff(rho, delrho, Amat2, Cmat2, npert, ipol,
     &                         nq, GRAD)
      implicit none
#include "dft2drv.fh"
c
      integer npert ! Number CPKS perturbations [input]
      integer ipol  ! [input]
      integer nq    ! [input]
      logical GRAD  ! Whether gradient-corrected [input]
c
c     Current approximate perturbed spin densities and density gradients
c     These are overwritten with the XC matrix coefficients to save space
c
      double precision rho(nq,ipol,npert), delrho(nq,3,ipol,npert)
c
c     Second derivatives of XC functional [input, destroyed]
c
      double precision Amat2(nq,NCOL_AMAT2), Cmat2(nq,NCOL_CMAT2)
c
      integer ipert, n ! Loop indices
      double precision rhoa, rhob
c
      if (GRAD) then
         call errquit('xc_cpks_coeff: Not yet implemented for GC',0)
      endif
c
      if (ipol.eq.1) then
c
c     Since the total densities are evaluated in the restricted case,
c     scale the functional derivatives by a factor of 0.5 so that the
c     correct CPKS matrices will be produced.
c
         call dscal(nq*NCOL_AMAT2,0.5d0,Amat2,1)
         if (GRAD) then
            call errquit('xc_cpks_coeff: Scale for restricted GC?',0)
c     !!! Is 1/2 the right factor?  1/4?
            call dscal(nq*NCOL_CMAT2,0.5d0,Cmat2,1)
         endif
      endif
c
      do ipert = 1,npert
c
c        !!! Put in cutoffs here similar to xc_tabcd? !!!
c
         if (ipol.eq.2) then
c
c           Unrestricted case
c
            do n = 1, nq
               rhoa = rho(n,1,ipert)
               rhob = rho(n,2,ipert)
               rho(n,1,ipert) = Amat2(n,D2_RA_RA)*rhoa
     &                        + Amat2(n,D2_RA_RB)*rhob
               rho(n,2,ipert) = Amat2(n,D2_RA_RB)*rhoa
     &                        + Amat2(n,D2_RB_RB)*rhob
            enddo
c
            if (GRAD) then
               do n = 1, nq
c                 !!! GC part goes here !!!
               enddo
            endif
c
         else
c
c           Restricted case
c
            do n = 1, nq
               rho(n,1,ipert) = (Amat2(n,D2_RA_RA) + Amat2(n,D2_RA_RB))
     &                        * rho(n,1,ipert)
            enddo
c
            if (GRAD) then
               do n = 1, nq
c                 !!! GC part goes here !!!
               enddo
            endif
c
         endif
c
      enddo
c
c      write(6,*) 'xc_cpks_coeff: rho out '
c      call output(rho, 1, nq, 1, ipol*npert, nq, ipol*npert, 1)
      return
      end
