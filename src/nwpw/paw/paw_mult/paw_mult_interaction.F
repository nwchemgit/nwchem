*
* $Id$
*
* $Log: not supported by cvs2svn $
* Revision 1.12  2004/07/29 15:55:37  bylaska
* Temporary array for Gaunt coefficients added.  Speeds up the program considerably, but it is extrememely memory intensive.  Also added timing routines to multipole calculations and fixed the initial total charge calculation.
*
*  ...EJB
*
* Revision 1.11  2003/03/24 20:31:08  bylaska
* ....EJB
*
* Revision 1.10  2003/03/22 01:03:46  bylaska
* bug in z-component of multipole force fixed....EJB
*
* Revision 1.9  2003/02/25 02:36:49  bylaska
* debug write statements removed...EJB
*
* Revision 1.8  2003/02/22 03:10:44  bylaska
* debugging multipole coefficients...There is currently a bug in
* find_comp_coeff...EJB
*
* Revision 1.7  2003/02/21 19:44:22  bylaska
* Routines for computing the mult_energy_coeff have been added to paw_mult
* ...EJB
*


!*************************************************
!
!   Name    :  paw_gaussian_mult_interaction
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
!   Created : Spring, 1998
!*************************************************

      complex*16 FUNCTION paw_mult_interaction(
     >     l1,m1,sigma1,R1,l2,m2,sigma2,R2)
      implicit none
      integer     l1
      integer     m1
      integer     l2
      integer     m2
      double precision    sigma1
      double precision    sigma2
      double precision    R1(3)
      double precision    R2(3)
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

        call nwpw_timing_start(14)
        sigma = dsqrt((sigma1**2+sigma2**2)/2)
c        R     = SQRT(DOT_PRODUCT(R1-R2,R1-R2))
        R = dsqrt( (R1(1)-R2(1))**2 +
     >            (R1(2)-R2(2))**2 +
     >            (R1(3)-R2(3))**2 )

        IF(R .eq. 0.0d0) THEN

        IF(l1 .ne. l2 .OR. m1 .ne. -m2 ) THEN

        f = dcmplx(0.0d0,0.0d0)

        ELSE

        f =  (-1)**m1*4.0d0*dsqrt(2.0d0*PI)/(
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
        phi = DATAN2(R1(2)-R2(2),R1(1)-R2(1))
        END IF

        f = dcmplx(0.0,0.0)

        m = m1 + m2

        phase_factor = (-1)**(m1+l1)/sigma**(l1+l2+1)
        r_over_s = R/sigma

        DO l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2

          gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)
          !write(*,*) "LL,MM,coeff=",l,m,l2,m2,l1,m1,gaunt_coeff

            f = f + paw_get_mult_coeff(l1,l2,l)*phase_factor*
     >           r_over_s**l*gaunt_coeff*
     >           theta_lm(l,m,cos_theta)*
     >           exp(dcmplx(0.0,m*phi))*
     >           paw_hypergeom1F1((l1+l2+l)/2,l,-(r_over_s)**2/2)

            phase_factor = -phase_factor

        END DO

        END IF
        call nwpw_timing_end(14)

        paw_mult_interaction = f

        return
        end ! gaussian_mult_interaction



!*************************************************
!
!   Name    :  paw_mult_force
!
!   Purpose : 
!
!   Created : Spring, 1998
!*************************************************

        SUBROUTINE paw_mult_force(l1,m1,sigma1,R1,
     >                            l2,m2,sigma2,R2,
     >                            fx,fy,fz)
        implicit none
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

        !*** external functions ***
        double precision paw_get_mult_coeff,gen_gaunt_coeff
        double precision paw_hypergeom1F1
        complex*16       find_Y_lm
        external         paw_get_mult_coeff,gen_gaunt_coeff
        external         paw_hypergeom1F1
        external         find_Y_lm

        call nwpw_timing_start(14)

        sigma = dsqrt((sigma1**2+sigma2**2)/2)
c        R     = SQRT(DOT_PRODUCT(R1-R2,R1-R2))
        R = dsqrt((R1(1)-R2(1))**2 +
     >            (R1(2)-R2(2))**2 +
     >            (R1(3)-R2(3))**2 )


        IF(R .eq. 0.0d0) THEN

        fx = dcmplx(0.0,0.0)
        fy = dcmplx(0.0,0.0)
        fz = dcmplx(0.0,0.0)
        RETURN
        END IF

        cos_theta = (R1(3)-R2(3))/R
        !sin_theta = dsqrt(1.0d0-cos_theta**2)
        sin_theta = ((R1(1)-R2(1))**2 + (R1(2)-R2(2))**2)/R**2
        sin_theta = dsqrt(sin_theta)

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


        !*** gaunt coefficent is non-zero ***
        if ((m.eq.(m2+m1))       .and.
     >      (mod(l+l2+l1,2).eq.0).and.
     >      (l1.le.(l+l2))       .and.
     >      (l1.ge.abs(l-l2)))    then

        gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

            
        call find_Y_lm_deriv(l,m,cos_theta,phi,Y_lm_x,Y_lm_y,Y_lm_z)
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
           fz = fz + tmp_cmplx*cos_theta

           end if

           phase_factor = -phase_factor

        END DO
        call nwpw_timing_end(14)

        return
        end 



!*************************************************
!
!   Name    :  paw_get_mult_coeff
!
!   Purpose :  calculates the
!              coefficient required for 
!              FUNCTION gaussian_mult_interaction
!
!   Created :  Spring,1998
!*************************************************   
      double precision FUNCTION paw_get_mult_coeff(l1,l2,l) 
      implicit none

        INTEGER  l1,l2,l
#include "paw_gamma_function.fh"

        double precision  c,pi
        pi=acos(-1.d0)

        c = sqrt((2.0d0*pi)**5)/sqrt(2.0d0**(2+l1+l2+l))*
     >      exp(paw_ln_gamma((l1+l2+l+1)/2.0d0)- paw_ln_gamma(l1+1.5d0)-
     >      paw_ln_gamma(l2+1.5d0) - paw_ln_gamma(l+1.5d0) )

        paw_get_mult_coeff = c

      return              
      END 

!**************************************************
!
!       Name: paw_set_mult_energy_coeff
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_set_mult_energy_coeff()
      implicit none


#include "bafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer nion,lmax2

      lmax2 = (paw_basis_max_mult_l()+1)**2
      nion  = ion_nion()
      call find_mult_energy_coeff(sigma_smooth,
     >                            nshl3d,dbl_mb(rcell(1)),
     >                            nion,lmax2,
     >                            dcpl_mb(mult_energy_coeff(1)))
      return
      end



!**************************************************
!
!       Name: find_mult_energy_coeff
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine find_mult_energy_coeff(sigma_smooth,
     >                            nshl3d,rcell,
     >                            nion,mult_l_max,mult_energy_coeff)
      implicit none
      real*8 sigma_smooth
      integer nshl3d
      real*8 rcell(3,nshl3d)
      integer nion,mult_l_max
      complex*16 mult_energy_coeff(nion,mult_l_max,nion,mult_l_max)


#include "paw_gaunt.fh"
#include "paw_hypergeom.fh"
#include "paw_spher_func.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer taskid,np,tcount
      integer il,ilp
      integer lshell
      integer ii,ia,jj,ja
      integer l,l1,l2,m,m1,m2
      integer mult_l,mult_l_prime,p1

      real*8 sigma_a, sigma_a_prime
      real*8 ra(3), ra_prime(3)
      real*8 R

      real*8 sigma4,sigma3,sigma2,sigma1
      real*8 cos_theta,phi
      real*8 gaunt_coeff
      real*8 r_over_s4, r_over_s3, r_over_s2, r_over_s1
      real*8 phase_factor4,phase_factor3
      real*8 phase_factor2,phase_factor1, phase_factor
      real*8 const4_l,const2_l,const1_l,const,twopi

      complex*16 f,f2,f3,f4,f5,f6

      !**** external functions ****
      integer    paw_double_factorial
      real*8     paw_get_mult_coeff
      complex*16 find_Y_lm
      external   paw_double_factorial
      external   paw_get_mult_coeff
      external   find_Y_lm


      call nwpw_timing_start(14)
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      twopi = 8.0d0*datan(1.0d0)

      call dcopy(2*(nion*mult_l_max)**2,0.0d0,0,mult_energy_coeff,1)

      tcount = 0
      do ii=1,nion
        ia   = ion_katm(ii)
        ra(1)   = ion_rion(1,ii)
        ra(2)   = ion_rion(2,ii)
        ra(3)   = ion_rion(3,ii)
        sigma_a = paw_basis_sigma(ia)
        mult_l  = paw_basis_mult_l(ia)

        do jj=1,ii-1
        do lshell=1,nshl3d
          ja = ion_katm(jj)
          ra_prime(1)   = ion_rion(1,jj) + rcell(1,lshell)
          ra_prime(2)   = ion_rion(2,jj) + rcell(2,lshell)
          ra_prime(3)   = ion_rion(3,jj) + rcell(3,lshell)
          sigma_a_prime = paw_basis_sigma(ja)
          mult_l_prime  = paw_basis_mult_l(ja)
        
          R = dsqrt( (ra(1)-ra_prime(1))**2
     >              +(ra(2)-ra_prime(2))**2
     >              +(ra(3)-ra_prime(3))**2)

          if (R.lt.(4*sigma_smooth)) then
          if (tcount.eq.taskid) then
            sigma4 = dsqrt((sigma_smooth**2+sigma_a_prime**2)/2.0d0)
            sigma3 = dsqrt((sigma_a**2+sigma_smooth**2)/2.0d0)
            sigma2 = dsqrt((sigma_smooth**2+sigma_smooth**2)/2.0d0)
            sigma1 = dsqrt((sigma_a**2+sigma_a_prime**2)/2.0d0)

            r_over_s4 = R/sigma4
            r_over_s3 = R/sigma3
            r_over_s2 = R/sigma2
            r_over_s1 = R/sigma1

            cos_theta = (ra(3)-ra_prime(3))/R

            if (      (ra(2).eq.ra_prime(2)) 
     >          .and. (ra(1).eq.ra_prime(1))) then
              phi = 0.0d0
            else
              phi = datan2(ra(2)-ra_prime(2),ra(1)-ra_prime(1))
            end if

            il = 1
            do l1=0,mult_l
            do m1=-l1,l1

              ilp = 1
              do l2=0,mult_l_prime
              do m2=-l2,l2
                      
                f = dcmplx(0.0d0,0.0d0)
                m = m1 + m2
                phase_factor4 = (-1)**(m1+l1)/sigma4**(l1+l2+1)
                phase_factor3 = (-1)**(m1+l1)/sigma3**(l1+l2+1)
                phase_factor2 = (-1)**(m1+l1)/sigma2**(l1+l2+1)
                phase_factor1 = (-1)**(m1+l1)/sigma1**(l1+l2+1)
                phase_factor  = 1.0d0
                do l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2

                !*** gaunt coefficent is non-zero ***
                if ((m.eq.(m2+m1))       .and.
     >              (mod(l+l2+l1,2).eq.0).and.
     >              (l1.le.(l+l2))       .and.
     >              (l1.ge.abs(l-l2)))    then
                gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

                p1 = (l1+l2+l)/2
                f = f + paw_get_mult_coeff(l1,l2,l)*phase_factor
     >                  *gaunt_coeff*find_Y_lm(l,m,cos_theta,phi)
     >                  *(
     >  -(r_over_s4**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s4)**2)*phase_factor4
     > - (r_over_s3**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s3)**2)*phase_factor3
     > + (r_over_s2**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s2)**2)*phase_factor2
     > + (r_over_s1**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s1)**2)*phase_factor1
     >                   ) 

                end if
                phase_factor = -phase_factor

                end do !*l*
                
                mult_energy_coeff(ii,il,jj,ilp) 
     >           = mult_energy_coeff(ii,il,jj,ilp) + f 


                ilp = ilp + 1
              end do !*m2*
              end do !*l2*

              il = il + 1
            end do !*m1*
            end do !*l1*

          end if
          tcount = mod(tcount+1,np)
          end if

        end do !*lshell*
        end do !*jj*
      end do !*ii*
      if (np.gt.1) 
     >call D3dB_Vector_SumAll(2*(nion*mult_l_max)**2,mult_energy_coeff)


      do ii=1,nion
       ia      = ion_katm(ii)
       mult_l  = paw_basis_mult_l(ia)
       sigma_a = paw_basis_sigma(ia)

       sigma4 = dsqrt((sigma_smooth**2+sigma_a**2)/2)
       sigma3 = dsqrt((sigma_a**2+sigma_smooth**2)/2)
       sigma2 = dsqrt((sigma_smooth**2+sigma_smooth**2)/2)
       sigma1 = dsqrt((sigma_a**2+sigma_a**2)/2)

        il = 1
        do l1=0,mult_l

          const = 4.0d0*dsqrt(twopi)
     >           / dble((2*l1+1)*paw_double_factorial(2*l1+1))

          const4_l = const/(sigma4**(2*l1+1))
          const2_l = const/(sigma2**(2*l1+1))
          const1_l = const/(sigma1**(2*l1+1))
          do m1=-l1,l1
            f = (-1)**m1*(-2*const4_l+const2_l+const1_l)

            ilp = 1
            do l2=0,mult_l
            do m2=-l2,l2

              if ((l1.eq.l2).and.(m1.eq.(-m2)) ) then
                mult_energy_coeff(ii,il,ii,ilp) = 
     >          mult_energy_coeff(ii,il,ii,ilp) + f
              end if

              f2 = dcmplx(0.0d0,0.0d0)
              do lshell=1,nshl3d
                 R = dsqrt( (rcell(1,lshell))**2
     >                    + (rcell(2,lshell))**2
     >                    + (rcell(3,lshell))**2)


                 if ((R.gt.1.0d-9).and.(R.lt.(4*sigma_smooth))) then
                 cos_theta = (rcell(3,lshell))/R
                 if (  (dabs(rcell(2,lshell)).lt.1.0d-9) 
     >           .and. (dabs(rcell(1,lshell)).lt.1.0d-9)) then
                    phi = 0.0d0
                 else
                    phi = datan2(rcell(2,lshell),rcell(1,lshell))
                 end if

                 r_over_s4 = R/sigma4
                 r_over_s3 = R/sigma3
                 r_over_s2 = R/sigma2
                 r_over_s1 = R/sigma1

                 m = m1 + m2
                 phase_factor4 = (-1)**(m1+l1)/sigma4**(l1+l2+1)
                 phase_factor3 = (-1)**(m1+l1)/sigma3**(l1+l2+1)
                 phase_factor2 = (-1)**(m1+l1)/sigma2**(l1+l2+1)
                 phase_factor1 = (-1)**(m1+l1)/sigma1**(l1+l2+1)
                 phase_factor  = 1.0d0

                 do l=l1+l2,MAX(ABS(l1-l2),ABS(m1+m2)),-2

                 !*** gaunt coefficent is non-zero ***
                if ((m.eq.(m2+m1))       .and.
     >              (mod(l+l2+l1,2).eq.0).and.
     >              (l1.le.(l+l2))       .and.
     >              (l1.ge.abs(l-l2)))    then
                gaunt_coeff = gen_gaunt_coeff(l,m,l2,m2,l1,-m1)

                p1 = (l1+l2+l)/2
                f2=f2+  paw_get_mult_coeff(l1,l2,l)*phase_factor
     >                  *gaunt_coeff*find_Y_lm(l,m,cos_theta,phi)
     >                  *(
     >  -(r_over_s4**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s4)**2)*phase_factor4
     > - (r_over_s3**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s3)**2)*phase_factor3
     > + (r_over_s2**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s2)**2)*phase_factor2
     > + (r_over_s1**l)
     >   *paw_hypergeom1F1(p1,l,-0.5d0*(r_over_s1)**2)*phase_factor1
     >                   ) 

                end if

                 phase_factor = -phase_factor
                 end do !*l*


               end if !*R!=0*
              end do !*lshell*


              mult_energy_coeff(ii,il,ii,ilp) = 
     >        mult_energy_coeff(ii,il,ii,ilp) + f2*2.0d0


              ilp = ilp+1
            end do !*m2*
            end do !*l2*

            il = il+1
          end do !*m1*

        end do !*l1*

      end do !*ii*


      do ii=1,nion
      do jj=ii+1,nion
        do il=1,mult_l_max
        do ilp=1,mult_l_max
              mult_energy_coeff(ii,il,jj,ilp) 
     >         = mult_energy_coeff(jj,ilp,ii,il)
        end do
        end do
      end do
      end do

      call nwpw_timing_end(14)
      return
      end

