c
c $Id: paw_spher_func.F,v 1.2 2003-02-10 02:33:32 marat Exp $
c
c $Log: not supported by cvs2svn $
c Revision 1.1  2003/02/06 06:07:12  marat
c initial version
c

!**************************************************
!
!	Name: spher_harmonics_generate
!
!	Purpose: 
!
!	Created:	8/05/2002
!**************************************************
      subroutine spher_harmonics_generate(l,m,nfft3d,
     >                                    Gx,Gy,Gz,
     >                                    Ylm)
      implicit none
      integer    l,m,nfft3d
      real*8     Gx(*),Gy(*),Gz(*)
      complex*16 Ylm(*)

*     **** local variables ****
      real*8 small
      parameter (small=1.0d-9)

      integer i,k,mabs
      real*8     g,gxy,cs_theta,sgn,coeff,fourpi
      complex*16 exp_phi

*     **** external functions ****
      real*8   legendre_lm
      external legendre_lm


      fourpi = 16.0d0*datan(1.0d0)
      sgn    = 1.0d0
      if (m.ge.0) then
        mabs = m
      else
        mabs = -m
        if (mod(mabs,2).eq.1) sgn = -1.0d0
      end if
     

      !*** find coefficient ***
      if(mabs.eq.0) then
         coeff = 1.0d0
      else
         coeff = 1.0d0
         do i=1,2*mabs
           coeff = coeff/(l-mabs+i)
         end do
      end if
      coeff = coeff*(2*l+1.0d0)/fourpi
      coeff = sgn*dsqrt(coeff)

      if (m.ge.0) then
        do k=1,nfft3d
          gxy = Gx(k)**2 + Gy(k)**2
          g   = gxy + Gz(k)**2
          gxy = dsqrt(gxy)
          g   = dsqrt(g) 

          if (gxy.gt.small)then
            exp_phi = (dcmplx(Gx(k),Gy(k))/gxy)**mabs 
          else
            exp_phi = dcmplx(1.0d0,0.0d0)
          end if 

          if (g.gt.small) then
            cs_theta = Gz(k)/g
          else
            cs_theta = 0.0d0
          end if 

          Ylm(k) = coeff*legendre_lm(l,mabs,cs_theta)
     >                  *exp_phi
        end do
      else
        do k=1,nfft3d
          gxy = Gx(k)**2 + Gy(k)**2
          g   = gxy + Gz(k)**2
          gxy = dsqrt(gxy)
          g   = dsqrt(g) 

          if (gxy.gt.small)then
            exp_phi = (dcmplx(Gx(k),Gy(k))/gxy)**mabs 
          else
            exp_phi = dcmplx(1.0d0,0.0d0)
          end if 

          if (g.gt.small) then
            cs_theta = Gz(k)/g
          else
            cs_theta = 0.0d0
          end if 
       
          Ylm(k) = coeff*legendre_lm(l,mabs,cs_theta)
     >                  *dconjg(exp_phi)
        end do
      end if
     
      return
      end

!**************************************************
!
!	Name: legendre_lm
!
!	Purpose: 
!
!	Created:	8/05/2002
!**************************************************
    
      function legendre_lm(l,m,x)
      implicit none
      integer l,m
      real*8 x
      double precision legendre_lm

*     *** local variables ****
      integer i,ll
      real*8  fact,pll,pmm,pmmp1,somx2,Plm

      if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.0d0))
     >  call errquit('legendre_lm: bad argments',0)

      pmm = 1.0d0
      if (m.gt.0) then
        somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
        fact  = 1.0d0
        do i=1,m
          pmm  = -pmm*fact*somx2
          fact = fact + 2.0d0
        end do
      end if

      if (l.eq.m) then
        Plm = pmm
      else
         pmmp1 = x*(2*m+1.0d0)*pmm
         if (l.eq.(m+1)) then
           Plm = pmmp1
         else
            do ll=m+2,l
              pll   = (x*(2*ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/dble(ll-m)
              pmm   = pmmp1
              pmmp1 = pll
            end do
            Plm = pll
         end if
      end if

      legendre_lm = Plm
      return
      end


!  *************************************************
!  
!     Name    : theta_lm
!  
!  
!     Purpose : calculates theta_lm for a scalar cos_theta
!               such that
!  
!         Y_lm(cos_theta,phi)=theta_lm(cos_theta)*exp(i*m*phi)
!  
!  
!     Created : 1/99
!  
!  *************************************************
      function theta_lm (l, m, cos_theta)
      implicit none


      integer l
      integer m
      double precision cos_theta
      double precision theta_lm
!    *** LOCAL VARIABLES ***
      integer i
      integer mod_m
      double precision coeff
      
#include "paw_params.fh"
      double precision legendre_lm
      external legendre_lm

      if (m .gt. l) call errquit(
     1   ' PARAMETER OUT OF ORDER IN FUNCTION THETA_LM ',1)
      mod_m = abs(m)
!    *** FIND COEFFICIENT ***
      if (mod_m .eq. 0) then
         coeff = 1.0
      else if (mod_m .gt. 0) then
         coeff = 1.0
         do i = 1, 2*mod_m
            coeff = coeff/(l - mod_m + i)
         end do
      endif
      coeff = coeff*(2*l + 1)/fourpi
      coeff = sqrt(coeff)
      theta_lm = coeff*legendre_lm(l,mod_m,cos_theta)
      if (m .lt. 0) then
         if (mod(mod_m,2) .eq. 1) theta_lm = -theta_lm
      endif
      end 
 
