
!*************************************************
!
!   Name    :  gaussian_mult_interaction
!
!
!   Purpose : calculates the Coulomb interaction energy
!             between two multipole gaussian charge
!             distributions. First charge distribution
!             is a multipole of order l1,m1, it is
!             located at R1 and the width of the
!             gaussian is sigma1. The same goes to
!             for the second charge.
!
!
!   Created : Spring, 1998
!
!*************************************************

        FUNCTION paw_mult_interaction(l1,m1,sigma1,R1,l2,m2,sigma2,R2)

        implicit none


        integer     l1
        integer     m1
        integer     l2
        integer     m2
        double precision    sigma1
        double precision    sigma2
        double precision    R1(3)
        double precision    R2(3)
        complex*16  paw_mult_interaction
        complex*16  f

#include "paw_gaunt.fh"
#include "paw_hypergeom.fh"
#include "paw_double_factorial.fh"
#include "paw_params.fh"
#include "paw_spher_func.fh"

        double precision paw_get_mult_coeff
        external paw_get_mult_coeff

        !*** local variables ***
        INTEGER  l,m
        double precision  cos_theta,phi,R,sigma
        double precision  gaunt_coeff
        double precision  r_over_s
        double precision  phase_factor

        sigma = SQRT((sigma1**2+sigma2**2)/2)
c        R     = SQRT(DOT_PRODUCT(R1-R2,R1-R2))
        R = SQRT( (R1(1)-R2(1))**2 +
     >            (R1(2)-R2(2))**2 +
     >            (R1(3)-R2(3))**2 )

        IF(R .eq. 0.0d0) THEN

        IF(l1 .ne. l2 .OR. m1 .ne. -m2 ) THEN

        f = dcmplx(0.0,0.0)

        ELSE

        f =  (-1)**m1*4.0*SQRT(2.0d0*PI)/(
     >              (2*l1+1)*                  
     >              paw_double_factorial(2*l1+1)*  
     >              sigma**(2*l1+1)            
     >              )

        END IF

        ELSE IF(R .ne. 0.0d0) THEN

        cos_theta = (R1(3)-R2(3))/R

        IF(R1(2) .EQ. R2(2) .AND. R1(1).EQ.R2(1)) THEN
        phi = 0.0d0
        ELSE
        phi = ATAN2(R1(2)-R2(2),R1(1)-R2(1))
        END IF

        f = dcmplx(0.0,0.0)

        m = m1 + m2

        phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
        r_over_s = R/sigma

        DO l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2

          
          gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

            f = f + paw_get_mult_coeff(l1,l2,l)*phase_factor*
     >           r_over_s**l*gaunt_coeff*
     >           theta_lm(l,m,cos_theta)*
     >           exp(dcmplx(0.0,m*phi))*
     >           paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s)**2/2)

            phase_factor = -phase_factor

        END DO

        END IF

        paw_mult_interaction = f

        END ! gaussian_mult_interaction


        SUBROUTINE paw_mult_force(l1,m1,sigma1,R1,
     >                                 l2,m2,sigma2,R2,
     >                                 fx,fy,fz)


        integer          l1
        integer          m1
        integer          l2
        integer          m2
        double precision         sigma1
        double precision         sigma2
        double precision         R1(3)
        double precision         R2(3)
        complex*16    fx
        complex*16    fy
        complex*16    fz

        !*** local variables ***
        double precision  phase_factor
        INTEGER  l,m
        double precision  sin_phi, cos_phi
        double precision  cos_theta,sin_theta,phi,R,sigma
        double precision  gaunt_coeff
        double precision  r_over_s,tmp_real
        double precision  kummer, dkummer
        complex*16  prefactor, dprefactor,tmp_cmplx
        complex*16  Y_lm_x,Y_lm_y,Y_lm_z,Y_lm

        sigma = SQRT((sigma1**2+sigma2**2)/2)
c        R     = SQRT(DOT_PRODUCT(R1-R2,R1-R2))
        R = SQRT( (R1(1)-R2(1))**2 +
     >            (R1(2)-R2(2))**2 +
     >            (R1(3)-R2(3))**2 )


        IF(R .eq. 0.0d0) THEN

        fx = dcmplx(0.0,0.0)
        fy = dcmplx(0.0,0.0)
        fz = dcmplx(0.0,0.0)
        RETURN
        END IF

        cos_theta = (R1(3)-R2(3))/R
        !sin_theta = SQRT(1.0d0-cos_theta**2)
        sin_theta = ((R1(1)-R2(1))**2 + (R1(2)-R2(2))**2)/R**2
        sin_theta = SQRT(sin_theta)

        IF(R1(2) .EQ. R2(2) .AND. R1(1).EQ.R2(1)) THEN
        phi = 0.0d0
        ELSE
        phi = ATAN2(R1(2)-R2(2),R1(1)-R2(1))
        END IF

        sin_phi = SIN(phi)
        cos_phi = COS(phi)

        fx = dcmplx(0.0,0.0)
        fy = dcmplx(0.0,0.0)
        fz = dcmplx(0.0,0.0)


        m = m1 + m2

        phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
        r_over_s = R/sigma

        DO l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2


        gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

            
        CALL find_Y_lm_deriv(l,m,cos_theta,phi,Y_lm_x,Y_lm_y,Y_lm_z)
        Y_lm =  find_Y_lm(l,m,cos_theta,phi)

        kummer = paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s)**2/2)

        dkummer = (-R/sigma**2)*((l1+l2+l+1.0d0)/2)/(l+1.5d0)*
     >                 paw_hypergeom1F1((l1+l2+l)
     >                 /2+1,l+1,-(r_over_s)**2/2)

           prefactor = phase_factor*paw_get_mult_coeff(l1,l2,l)*
     >                 (r_over_s)**l

           dprefactor = l*prefactor/R

           tmp_real = prefactor*gaunt_coeff*kummer/R
           fx = fx + tmp_real*Y_lm_x
           fy = fy + tmp_real*Y_lm_y
           fz = fz + tmp_real*Y_lm_z


           tmp_cmplx = dprefactor*gaunt_coeff*Y_lm*kummer 
           fx = fx + tmp_cmplx*sin_theta*cos_phi
           fy = fy + tmp_cmplx*sin_theta*sin_phi
           fz = fz + tmp_cmplx*cos_theta
        

           tmp_cmplx = prefactor*gaunt_coeff*Y_lm*dkummer
           fx = fx + tmp_cmplx*sin_theta*cos_phi
           fy = fy + tmp_cmplx*sin_theta*sin_phi
           fz = fz + tmp_cmplx*dkummer*cos_theta

           phase_factor = -phase_factor

        END DO

        END 



!*************************************************
!
!   Name    :  paw_get_mult_coeff
!
!
!   Purpose :  calculates the
!              coefficient required for 
!              FUNCTION gaussian_mult_interaction
!
!
!   Created :  Spring,1998
!
!*************************************************   
      FUNCTION paw_get_mult_coeff(l1,l2,l) 
      implicit none

        INTEGER  l1,l2,l
        double precision  paw_get_mult_coeff
#include "paw_gamma_function.fh"

        double precision  c,pi
        pi=acos(-1.d0)

        c = sqrt((2.0d0*pi)**5)/sqrt(2.0d0**(2+l1+l2+l))*
     >      exp(paw_ln_gamma((l1+l2+l+1)/2.0d0)- paw_ln_gamma(l1+1.5d0)-
     >      paw_ln_gamma(l2+1.5d0) - paw_ln_gamma(l+1.5d0) )

        paw_get_mult_coeff = c

              
      END 


