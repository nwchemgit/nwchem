c
c $Id: paw_spher_func.F,v 1.3 2003-02-21 03:02:19 bylaska Exp $
c
c $Log: not supported by cvs2svn $
c Revision 1.2  2003/02/10 02:33:32  marat
c fixed bug with undeclared variables in theta_lm
c
c
c MV
c
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
 
!**************************************************
!
!	Name: find_Y_lm
!
!	Purpose: 
!
!	Created:	8/05/2002
!**************************************************
      function find_Y_lm(l,m,cos_theta,phi)
      implicit none
      integer l,m
      real*8 cos_theta,phi
      complex*16 find_Y_lm !* RESULT*

      !*** external function ****
      real*8   theta_lm
      external theta_lm

      find_Y_lm = theta_lm(l,m,cos_theta)
     >            *dcmplx(dcos(m*phi),dsin(m*phi))
      return
      end 



c  !*************************************************
c  !
c  !   Name    : find_Y_lm_deriv
c  !
c  !
c  !   Purpose : to calculate derivative of the spherical harmonic
c  !             for scalar argument,
c  !             namely g_r*d(Ylm)/dx
c  !                    g_r*d(Ylm)/dy
c  !                    g_r*d(Ylm)/dz
c  !
c  !
c  !
c  !   Created : 1999
c  !
c  !   Notes:
c  !
c  !
c  !*************************************************
c      subroutine find_Y_lm_deriv(l,m,cos_theta,angle_phi,
c     >                           Y_lm_x,Y_lm_y,Y_lm_z)
c      implicit none
c      integer     l,m
c      real*8      cos_theta
c      real*8      angle_phi
c      complex*16  Y_lm_x
c      complex*16  Y_lm_y
c      complex*16  Y_lm_z
c
c      !*** local variables ***
c      complex*16  f1, f2
c
c      !**** external functions ****
c      real*8   dtheta_lm,theta_lm_div
c      external dtheta_lm,theta_lm_div
c
c
c      !*** derivative with respect to theta ***
c      f1 = dtheta_lm(l,m,cos_theta)
c     >     *dcmplx(cos(m*angle_phi),sin(m*angle_phi))
c
c      !*** derivative with respect to phi ***
c      IF(m.eq.) THEN
c        f2 = dcmplx(0.0d0,0.0d0)
c      ELSE
c        f2 = theta_lm_div(l,m,cos_theta)
c     >       *dcmplx(cos(m*angle_phi),sin(m*angle_phi))
c     >       m*dcmplx(0.0d0,1.0d0)
c      END IF
c
c      !*** cartesian derivatives (see Ashkroft p.108, 3rd. ed)
c      Y_lm_x =  f1*cos_theta*COS(angle_phi) - SIN(angle_phi)*f2
c      Y_lm_y =  f1*cos_theta*SIN(angle_phi) + COS(angle_phi)*f2
c      Y_lm_z = -f1*SQRT(1.0_DP - cos_theta**2)
c
c      return
c      end
c
c
c!**************************************************
c!
c!	Name: dtheta_lm
c!
c!	Purpose: 
c!
c!	Created:	8/05/2002
c!**************************************************
c      real*8 FUNCTION dtheta_lm(l,m,cos_theta)
c      implicit none
c      INTEGER l
c      INTEGER m
c      real*8  cos_theta
c
c      !*** local variables ***
c      integer i
c      integer mod_m
c      real*8  coeff,tmp,fourpi
c
c      fourpi = 16.0d0*datan(1.0d0)
c
c      IF( m .gt.l ) THEN
c      CALL errquit("parameter out of order in function theta_lm",1)
c      END IF
c
c      mod_m = ABS(m)
c
c      !*** find coefficient ***
c      if (mod_m.eq.0) then
c        coeff = 1.0_DP
c      else if(mod_m.gt.0) then
c        coeff = 1.0_DP
c        do i=1,2*mod_m
c          coeff = coeff/(l-mod_m+i)
c        end do
c      end if
c
c      coeff = coeff*(2*l+1)/(fourpi)
c      coeff = dsqrt(coeff)
c
c      tmp = coeff*dlegendre_lm(l,mod_m,cos_theta)
c
c      if (m.lt.0) then
c        if( MOD(mod_m,2).eq.1) then
c           tmp = -tmp
c        end if
c      end if
c
c      dtheta_lm = tmp
c      return
c      end



