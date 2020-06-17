c
c $Id$
c
c $Log: not supported by cvs2svn $
c Revision 1.7  2003/10/21 02:05:17  marat
c switched to new errquit by running global replace operation
c see the script below (note it will not work on multiline errquit calls)
c *********************************************************
c #!/bin/sh
c
c e=`find . -name "*F" -print`
c
c for f in $e
c do
c cp $f $f.bak
c sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
c #rm $f.bak
c done
c **********************************************************
c
c Revision 1.6  2003/03/25 22:22:20  edo
c bug fix in paw_spher_func.F (uninitialized variable)
c
c Revision 1.2.2.2  2003/03/25 22:21:03  edo
c bug fix in paw_spher_func.F (uninitialized variable)
c
c files copied over from main branch
c
c Revision 1.2.2.1  2003/03/24 02:15:29  bylaska
c PAW added to release 4.5 tree....EJB
c
c Revision 1.5  2003/03/18 01:00:17  marat
c fixed error in evaluating legendre function
c for |m|>l
c MV
c
c Revision 1.4  2003/02/21 19:27:11  bylaska
c find_Y_lm_deriv subourine added...EJB
c
c Revision 1.3  2003/02/21 03:02:19  bylaska
c find_Y_lm added...EJB
c
c Revision 1.2  2003/02/10 02:33:32  marat
c fixed bug with undeclared variables in theta_lm
c
c
c MV
c
c Revision 1.1  2003/02/06 06:07:12  marat
c initial version
c

c!**************************************************
c!
c!      Name: spher_harmonics_generate
c!
c!      Purpose: 
c!
c!      Created:        8/05/2002
c!**************************************************
      subroutine spher_harmonics_generate(l,m,nfft3d,
     >                                    Gx,Gy,Gz,
     >                                    Ylm)
      implicit none
      integer          l,m,nfft3d
      double precision Gx(*),Gy(*),Gz(*)
      double complex   Ylm(*)

*     **** local variables ****
      double precision small
      parameter (small=1.0d-9)

      integer          i,k,mabs
      double precision g,gxy,cs_theta,sgn,coeff,fourpi
      double complex   exp_phi

*     **** external functions ****
      double precision legendre_lm
      external         legendre_lm


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


C      DOUBLE PRECISION FUNCTION LEGENDRE_LM (L, M, X)
C 
C      INTEGER L
C      INTEGER M
C      DOUBLE PRECISION X
C 
C 
C!    *** local variables ***
C      INTEGER I
C      INTEGER MOD_M
C      DOUBLE PRECISION P_MM
C      DOUBLE PRECISION P_MP1M
C      DOUBLE PRECISION TMP
C      DOUBLE PRECISION FACT
C      DOUBLE PRECISION COEFF
C 
C 
C      MOD_M = ABS(M)
C      IF (ABS(X) .GT. 1) CALL errquit(
C     1   'parameter ot of range in legendre_lm',0,0)
C 
C      IF (MOD_M .GT. L) THEN
C         LEGENDRE_LM = 0d0
C         RETURN 
C      ENDIF
C 
C!    *** find P(mod_m,mod_m) for mod_m=0 ***
C      P_MM = 1.0D0
C 
C!    *** find P(mod_m,mod_m) for mod_m > 0 ***
C      IF (MOD_M .GT. 0) THEN
C 
C         TMP = SQRT((1.0D0 - X)*(1.0D0 + X))
C 
C         FACT = 1.0D0
C 
C         DO I = 1, MOD_M
C 
C            P_MM = -P_MM*FACT*TMP
C 
C            FACT = FACT + 2.0D0
C 
C         END DO
C 
C      ENDIF
C 
C!    *** find P(l,mod_m) ***
C      IF (MOD_M .EQ. L) THEN
C 
C         LEGENDRE_LM = P_MM
C 
C      ELSE
C 
C      !*** find P(mod_m+1,mod_m) ***
C         P_MP1M = X*(2*MOD_M + 1)*P_MM
C 
C         IF (L .EQ. MOD_M+1) THEN
C 
C            LEGENDRE_LM = P_MP1M
C 
C         ELSE
C 
C            DO I = MOD_M + 2, L
C 
C               TMP=(X*(2*I-1)*P_MP1M-(I+MOD_M-1)*P_MM)/(I-MOD_M)
C 
C               P_MM = P_MP1M
C               P_MP1M = TMP
C 
C            END DO
C 
C            LEGENDRE_LM = TMP
C 
C         ENDIF
C 
C      ENDIF
C 
Cc!    *** negative m - this routine is only call with negative m from dtheta_lm and ddtheta_lm ***
C      IF (M .LT. 0) THEN
C 
C         COEFF = 1.0D0
C 
C         DO I = 1, 2*MOD_M
C 
C            COEFF = COEFF/(L - MOD_M + I)
C 
C         END DO
C 
C         COEFF = COEFF*(-1)**MOD_M
C 
C         LEGENDRE_LM = LEGENDRE_LM*COEFF
C 
C      ENDIF
C 
C 
C      END 
 
cc!**************************************************
cc!
cc!     Name: legendre_lm
cc!
cc!     Purpose: 
cc!
cc!     Created:        8/05/2002
cc!**************************************************
c    
c      function legendre_lm1(l,m,x)
c      implicit none
c      integer l,m
c      real*8 x
c      double precision legendre_lm1
c
c*     *** local variables ****
c      integer i,ll
c      real*8  fact,pll,pmm,pmmp1,somx2,Plm
c
c      if ((m.lt.0).or.(m.gt.l).or.(dabs(x).gt.1.0d0))
c     >  call errquit('legendre_lm: bad argments',0,0)
c
c      pmm = 1.0d0
c      if (m.gt.0) then
c        somx2 = dsqrt((1.0d0-x)*(1.0d0+x))
c        fact  = 1.0d0
c        do i=1,m
c          pmm  = -pmm*fact*somx2
c          fact = fact + 2.0d0
c        end do
c      end if
c
c      if (l.eq.m) then
c        Plm = pmm
c      else
c         pmmp1 = x*(2*m+1.0d0)*pmm
c         if (l.eq.(m+1)) then
c           Plm = pmmp1
c         else
c            do ll=m+2,l
c              pll   = (x*(2*ll-1.0d0)*pmmp1-(ll+m-1.0d0)*pmm)/dble(ll-m)
c              pmm   = pmmp1
c              pmmp1 = pll
c            end do
c            Plm = pll
c         end if
c      end if
c
c      legendre_lm1 = Plm
c      return
c      end


Cc!  *************************************************
Cc!  
Cc!     Name    : theta_lm
Cc!  
Cc!  
Cc!     Purpose : calculates theta_lm for a scalar cos_theta
Cc!               such that
Cc!  
Cc!         Y_lm(cos_theta,phi)=theta_lm(cos_theta)*exp(i*m*phi)
Cc!  
Cc!  
Cc!     Created : 1/99
Cc!  
Cc!  *************************************************
C      double precision function theta_lm(l,m,cos_theta)
C      implicit none
C      integer l,m
C      double precision cos_theta
Cc      double precision theta_lm !*RESULT*
C
C      !*** local variables ***
C      integer i,mod_m
C      double precision coeff,fourpi,f
C
C      !*** external functions *** 
C      double precision legendre_lm
C      external         legendre_lm
C
C      fourpi = 16.0d0*datan(1.0d0)
C
C      if (m .gt. l) 
C     > call errquit(
C     >      'parameter out of order in function theta_lm',0,0)
C
C      mod_m = abs(m)
C
C!     *** find coefficient ***
C      if (mod_m.eq.0) then
C         coeff= 1.0d0
C      else if (mod_m.gt.0) then
C         coeff= 1.0d0
C         do i=1,2*mod_m
C            coeff = coeff/dble(l-mod_m+i)
C         end do
C      endif
C      coeff = coeff*(2*l+1)/fourpi
C      coeff = sqrt(coeff)
C      f     = coeff*legendre_lm(l,mod_m,cos_theta)
C      if (m.lt.0) then
C         if (mod(mod_m,2).eq.1) f=-f
C      endif
C
C      theta_lm = f
C      return
C      end 
 
!**************************************************
!
!       Name: find_Y_lm
!
!       Purpose: 
!
!       Created:        8/05/2002
!**************************************************
      double complex function find_Y_lm(l,m,cos_theta,phi)
      implicit none
      integer l,m
      double precision cos_theta,phi
c      complex*16 find_Y_lm !* RESULT*

      !*** external function ****
      double precision theta_lm
      external         theta_lm

      find_Y_lm = theta_lm(l,m,cos_theta)
     >            *dcmplx(dcos(m*phi),dsin(m*phi))
      return
      end 



c !*************************************************
c !
c !   Name    : find_Y_lm_deriv
c !
c !
c !   Purpose : to calculate derivative of the spherical harmonic
c !             for scalar argument,
c !             namely g_r*d(Ylm)/dx
c !                    g_r*d(Ylm)/dy
c !                    g_r*d(Ylm)/dz
c !
c !
c !
c !   Created : 1999
c !
c !   Notes:
c !
c !
c !*************************************************
      subroutine find_Y_lm_deriv(l,m,cos_theta,angle_phi,
     >                           Y_lm_x,Y_lm_y,Y_lm_z)
      implicit none
      integer     l,m
      real*8      cos_theta
      real*8      angle_phi
      complex*16  Y_lm_x
      complex*16  Y_lm_y
      complex*16  Y_lm_z

      !*** local variables ***
      complex*16  f1, f2

      !**** external functions ****
      real*8   dtheta_lm,theta_lm_div
      external dtheta_lm,theta_lm_div


      !*** derivative with respect to theta ***
      f1 = dtheta_lm(l,m,cos_theta)
     >     *dcmplx(cos(m*angle_phi),sin(m*angle_phi))

      !*** derivative with respect to phi ***
      if(m.eq.0) THEN
        f2 = dcmplx(0.0d0,0.0d0)
      else
        f2 = theta_lm_div(l,m,cos_theta)
     >       *dcmplx(cos(m*angle_phi),sin(m*angle_phi))
     >       *m*dcmplx(0.0d0,1.0d0)
      end if

      !*** cartesian derivatives (see Arfken p.108, 3rd. ed) ***
      Y_lm_x =  f1*cos_theta*COS(angle_phi)-sin(angle_phi)*f2
      Y_lm_y =  f1*cos_theta*SIN(angle_phi)+cos(angle_phi)*f2
      Y_lm_z = -f1*dsqrt(1.0d0-cos_theta**2)

      return
      end


!**************************************************
!
!       Name: dtheta_lm
!
!       Purpose: 
!
!       Created:        8/05/2002
!**************************************************
      double precision function dtheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision  cos_theta

      !*** local variables ***
      integer i
      integer mod_m
      double precision  coeff,tmp,fourpi

      !*** external functions ***
      double precision   dlegendre_lm
      external           dlegendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      IF( m .gt.l ) THEN
      CALL errquit("parameter out of order in function dtheta_lm",0,1)
      END IF

      mod_m = ABS(m)

      !*** find coefficient ***
      if (mod_m.eq.0) then
        coeff = 1.0d0
      else if(mod_m.gt.0) then
        coeff = 1.0d0
        do i=1,2*mod_m
          coeff = coeff/dble(l-mod_m+i)
        end do
      end if

      coeff = coeff*(2*l+1)/(fourpi)
      coeff = dsqrt(coeff)

      

      tmp = coeff*dlegendre_lm(l,mod_m,cos_theta)

      if (m.lt.0) then
        if( MOD(mod_m,2).eq.1) then
           tmp = -tmp
        end if
      end if

      dtheta_lm = tmp
      return
      end



!**************************************************
!
!       Name: ddtheta_lm
!
!       Purpose:
!
!       Created:        8/05/2002
!**************************************************
      double precision function ddtheta_lm(l,m,cos_theta)
      implicit none
      integer l,m
      double precision  cos_theta

      !*** local variables ***
      integer i
      integer mod_m
      double precision coeff,tmp,fourpi

      !*** external functions ***
      double precision ddlegendre_lm
      external         ddlegendre_lm

      fourpi = 16.0d0*datan(1.0d0)

      IF( m .gt.l ) THEN
      CALL errquit("parameter out of order in function ddtheta_lm",0,1)
      END IF

      mod_m = ABS(m)

      !*** find coefficient ***
      if (mod_m.eq.0) then
        coeff = 1.0d0
      else if(mod_m.gt.0) then
        coeff = 1.0d0
        do i=1,2*mod_m
          coeff = coeff/dble(l-mod_m+i)
        end do
      end if

      coeff = coeff*(2*l+1)/(fourpi)
      coeff = dsqrt(coeff)


      tmp = coeff*ddlegendre_lm(l,mod_m,cos_theta)

      if (m.lt.0) then
        if( MOD(mod_m,2).eq.1) then
           tmp = -tmp
        end if
      end if

      ddtheta_lm = tmp
      return
      end




Cc !*************************************************
Cc !
Cc !   Name    : dlegendre_lm
Cc !
Cc !
Cc !   Purpose : calculates the derivative of the
Cc !             associated Legendre with respect to angle theta
Cc !             (x=cos(theta))
Cc !
Cc !
Cc !   Created : 10/99
Cc !
Cc !
Cc !   Notes: Arfken definition of P(l,m,x) is different from mine by a
Cc !          factor (-1)**m. Thus all recursive relationships should be modified
Cc !          accordingly
Cc !*************************************************
C      double precision function dlegendre_lm(l,m,x) 
C      implicit none
C      integer l,m
C      double precision  x
Cc      real*8  dlegendre_lm
C
C      !*** local variables ***
C      double precision f1,f2
C
C      !*** external functions ****
C      double precision legendre_lm
C      external         legendre_lm
C
C      !*** check the arguments ***
C      if ((m.gt.l) .or. (dabs(x).gt.1.0d0)) then
C        call errquit("parameter ot of range in dlegendre_lm",0,0)
C      end if
C
C      if (m.eq.l) then
C        f1 = 0.0d0
C      else
C        f1 = legendre_lm(l,m+1,x)
C      END IF
C
C      f2 = legendre_lm(l,m-1,x)
C
C      !*** see Arfken p.669 3rd ed, p. 725 4th ed.                             ***
C      !*** sqrt(1-x*x)*dP_l^m/dx = 0.5*P_l^(m+1) - 0.5*(l+m)*(l-m+1)*P_l^(m-1) ***
C      !*** Note this recursion relationship doesn't have the problems at       ***
C      !*** theta=0,pi that the more traditional recursion relationship in term ***
C      !*** of l (i.e. p. 115 of Jahnke and Emde) have                          ***
C      !*** (1-x*x)*dP_l^m/dx = (l+1)*x*P_l^m - (l-m+1)*P_(l+1)^m               ***
C
C      dlegendre_lm = 0.5d0*f1 - 0.5d0*(l+m)*(l-m+1)*f2
C
C      return
C      end 


c !*************************************************
c !
c !   Name    : ddlegendre_lm
c !
c !
c !   Purpose : calculates the derivative of the
c !             associated Legendre with respect to angle theta
c !             (x=cos(theta))
c !
c !
c !   Created : 10/99
c !
c !
c !   Notes: Arfken definition of P(l,m,x) is different from mine by a
c !          factor (-1)**m. Thus all recursive relationships should be modified
c !          accordingly
c !*************************************************
      double precision function ddlegendre_lm(l,m,x)
      implicit none
      integer l,m
      double precision x

      !*** local variables ***
      double precision f1,f2

      !*** external functions ****
      double precision dlegendre_lm
      external         dlegendre_lm

      !*** check the arguments ***
      if ((m.gt.l) .or. (dabs(x).gt.1.0d0)) then
        write(*,*) "parameter ot of range in ddlegendre_lm"
        stop
      end if

      if (m.eq.l) then
        f1 = 0.0d0
      else
        f1 = dlegendre_lm(l,m+1,x)
      END IF

      f2 = dlegendre_lm(l,m-1,x)

      !*** derived from by differentiation of recursion formula used for dlegendre_lm ***
      !*** see Arfken p.669 3rd ed, p. 725 4th ed                                     ***
      ddlegendre_lm = 0.5d0*f1 - 0.5d0*(l+m)*(l-m+1)*f2

      return
      end




c  !*************************************************
c  !
c  !   Name    : theta_lm_div
c  !
c  !
c  !   Purpose : calculates theta_lm/sin(theta)
c  !
c  !
c  !   Created : 10/99
c  !
c  !*************************************************
      double precision function theta_lm_div(l,m,cos_theta)
      implicit none
      integer l,m
      double precision cos_theta
c      real*8 theta_lm_div !*RESULT*

      !*** local variables ***
      integer i,mod_m
      double precision coeff,fourpi,f

      !*** external functions ***
      double precision legendre_lm_div
      external         legendre_lm_div

      fourpi = 16.0d0*datan(1.0d0)

      if ( m.gt.l) then
        call errquit(
     >   'parameter out of order in function theta_lm_div',0,0)
      end if

      mod_m = abs(m)

      !*** find coefficient ***
      if (mod_m.eq.0) then
        coeff = 1.0d0
      else if (mod_m.gt.0) then
        coeff = 1.0d0
        do i=1,2*mod_m
          coeff=coeff/dble(l-mod_m+i)
        end do
      end if

      coeff = coeff*(2*l+1)/fourpi
      coeff = dsqrt(coeff)

      f = coeff*legendre_lm_div(l,mod_m,cos_theta)

      if (m.lt.0) then
        if (mod(mod_m,2).eq.1) then
           f = -f
        END IF
      END IF

      theta_lm_div = f
      return
      end 


Cc  !*************************************************
Cc  !
Cc  !   Name    : legendre_lm_div
Cc  !
Cc  !   Purpose : calculates associated Legendre
Cc  !             polynomial divided by sin(theta)
Cc  !             for a scalar argument
Cc  !
Cc  !   Created : 10/99
Cc  !
Cc  !*************************************************
C      double precision function legendre_lm_div(l,m,x) 
C      implicit none
C      integer l,m
C      double precision x
Cc      real*8 legendre_lm_div !*RESULT*
C
C      !*** local variables ***
C      integer i
C      double precision p_mm,p_mp1m,tmp,fact,f
C
C      !*** check the arguments ***
C      if ((m.lt.0)           .or. 
C     >    (m.gt.l)           .or. 
C     >    (dabs(x).gt.1.0d0) .or.
C     >    (m.eq.0) ) then
C        call errquit(
C     >        'parameter ot of range in legendre_lm_div',0,0)
C      end if
C
C      !*** P(1,1,x)/sqrt(1-x**2) ***
C      p_mm = -1.0d0
C
C      fact = 3.0d0
C      tmp = dsqrt((1.0d0-x)*(1.0d0+x))
C
C      do i=2,m
C        p_mm = -p_mm*fact*tmp
C        fact = fact + 2.0d0
C      end do
C
C      !*** find P(l,m) ***
C      if (m.eq.l) then
C        f = p_mm
C
C      else
C          !*** find P(m+1,m) ***
C          p_mp1m = x*(2*m+1)*p_mm
C          if (l.eq.(m+1)) then
C            f = p_mp1m
C          else
C            do i=m+2,l
C              tmp = (x*(2*i-1)*p_mp1m - (i+m-1)*p_mm)/(i-m)
C              p_mm = p_mp1m
C              p_mp1m = tmp
C            end do
C            f = tmp
C          end if
C      end if
C
C      legendre_lm_div = f
C      return
C      end




