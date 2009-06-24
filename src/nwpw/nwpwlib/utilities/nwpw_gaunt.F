*     ***************************
*     *				*
*     *   nwpw_factorial_init   *
*     *				*
*     ***************************

      subroutine nwpw_factorial_init()
      implicit none


*     **** Factorial common block ****
      integer FMAX
      parameter (FMAX=30)
      real*8 factorial(0:FMAX)
      common / FACTORIAL_Block / factorial

*     **** local variables ****
      integer n
      real*8  x


      factorial(0) = 1.0d0
      do n=1,FMAX
         x = dble(n)
         factorial(n) = x*factorial(n-1)
      end do

      return
      end

*     ***************************
*     *				*
*     *      nwpw_clebsch	*
*     *				*
*     ***************************
*
* Computes the Clebsch-Gordon coefficients using Racah's formula.
*

      real*8 function nwpw_clebsch(l1,l2,l,m1,m2,m)
      implicit none
      integer l1,l2,l
      integer m1,m2,m


*     **** Factorial common block ****
      integer FMAX
      parameter (FMAX=30)
      real*8 factorial(0:FMAX)
      common / FACTORIAL_Block / factorial

*     **** local variables ****
      integer z,zmin,zmax
      real*8 tmp1,tmp2,sign

*     **************************************
*     **** check the triangle condition ****
*     **************************************

*     *** a zero Clebsch-Gordon coefficient ***
      if ( ((m1+m2).ne.m)      .OR.
     >     (abs(m1).gt.abs(l1)).OR.
     >     (abs(m2).gt.abs(l2)).OR.
     >     (abs(m) .gt.abs(l) ).OR.
     >     (l.gt.(l1+l2))      .OR.
     >     (l.lt.abs(l1-l2))   .OR.
     >     (mod((l1+l2+l),2).ne. 0) ) then
         tmp2 = 0.0d0

*     *** a non-zero Clebsch-Gordon coeffient ****
      else
         zmin = 0
         if ((l-l2+m1)     .lt.0) zmin = -l + l2 - m1
         if ((l-l1-m2+zmin).lt.0) zmin = -l + l1 + m2
         zmax = l1+l2-l
         if ((l2+m2-zmax).lt.0) zmax = l2 + m2
         if ((l1-m1-zmax).lt.0) zmax = l1 - m1

         tmp1 = 0.0d0
         do z=zmin,zmax
            sign = 1.0d0
            if (mod(z,2).ne.0) sign = -1.0d0
            tmp1 = tmp1 
     >           + sign
     >             / ( factorial(  (z)        )
     >               * factorial( (l1+l2-l-z) )
     >               * factorial( (l1-m1-z)   )
     >               * factorial( (l2+m2-z)   )
     >               * factorial( (l-l2+m1+z) )
     >               * factorial( (l-l1-m2+z) ))
         end do
         tmp2 = dsqrt( dble(2*l+1) * factorial( (l1+l2-l)   )
     >                           * factorial( (l1-l2+l)   )
     >                           * factorial( (-l1+l2+l)  )
     >                           * factorial( (l1+m1)     )
     >                           * factorial( (l1-m1)     )
     >                           * factorial( (l2+m2)     )
     >                           * factorial( (l2-m2)     )
     >                           * factorial( (l+m)       )
     >                           * factorial( (l-m)       )
     >                           / factorial( (l1+l2+l+1) )
     >               )*tmp1
      end if

      nwpw_clebsch = tmp2
      return
      end



*     ***************************
*     *				*
*     *      nwpw_Gaunt		*
*     *				*
*     ***************************
*
* Computes the Gaunt coefficient
*
*  Gaunt(l,m,l1,m1,l2,m2) =

*      /2pi  / pi
*     |     |
*   = |     | dconjg(Y_lm(theta,phi)) * Y_l1m1(theta,phi) * dconjg(Y_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*      /2pi  / pi
*     |     |
*   = |     | Y_lm(theta,phi) * dconjg(Y_l1m1(theta,phi)) * Y_l2m2(theta,phi)  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*   = sqrt( (2*l+1)*(2*l2+1)/(4*pi*(2*l1+1)) ) * Clebsch(l l2 l1; m m2 m1) * Clebsch(l l2 l1; 0 0 0)
*
*

      real*8 function nwpw_Gaunt(l,m,l1,m1,l2,m2)
      implicit none
      integer l,m
      integer l1,m1
      integer l2,m2

*     **** local variables ****
      real*8 coeff,fourpi
      real*8 tmp
   
*     **** external functions ***
      real*8   Clebsch
      external Clebsch


*     *** a zero Gaunt coefficient ***
      if ( ((m+m2).ne.m1)      .OR.
     >     (abs(m1).gt.abs(l1)).OR.
     >     (abs(m2).gt.abs(l2)).OR.
     >     (abs(m) .gt.abs(l) ).OR.
     >     (l1.gt.(l+l2))      .OR.
     >     (l1.lt.abs(l-l2))   .OR.
     >     (mod((l1+l2+l),2).ne. 0 )) then
         tmp = 0.0d0

*     *** a non-zero Gaunt Coefficient ***
      else

         fourpi = 16.0d0*datan(1.0d0)
         coeff  = dsqrt( dble((2*l+1)*(2*l2+1))/(fourpi*dble(2*l1+1)))
         tmp    = coeff*Clebsch(l,l2,l1, m,m2,m1)
     >                 *Clebsch(l,l2,l1, 0,0, 0)
      end if

      nwpw_Gaunt = tmp
      return
      end


*     ***************************
*     *				*
*     *      nwpw_Taunt		*
*     *				*
*     ***************************
*
* Computes the Taunt coefficient
*
*  Taunt(l,m,l1,m1,l2,m2) =

*      /2pi  / pi
*     |     |
*   = |     | T_lm(theta,phi)) * T_l1m1(theta,phi) * T_l2m2(theta,phi))  sin(theta) dtheta dphi
*     |     |
*    /0    / 0
*
*
*   = sqrt( (2*l+1)*(2*l2+1)/(2*pi*(2*l1+1)) ) * Clebsch(l l2 l1; m m2 m1) * Clebsch(l l2 l1; 0 0 0)
*
*

      real*8 function nwpw_Taunt(l,m,l1,m1,l2,m2)
      implicit none
      integer l,m
      integer l1,m1
      integer l2,m2

*     **** local variables ****
      real*8 coeff,twopi
      real*8 tmp
   
*     **** external functions ***
      real*8   Clebsch
      external Clebsch


*     *** a zero Taunt coefficient ***
      if ( ((m+m2).ne.m1)      .OR.
     >     (abs(m1).gt.abs(l1)).OR.
     >     (abs(m2).gt.abs(l2)).OR.
     >     (abs(m) .gt.abs(l) ).OR.
     >     (l1.gt.(l+l2))      .OR.
     >     (l1.lt.abs(l-l2))   .OR.
     >     (mod((l1+l2+l),2).ne. 0 )) then
         tmp = 0.0d0

*     *** a non-zero Taunt Coefficient ***
      else

         twopi  = 8.0d0*datan(1.0d0)
         coeff  = dsqrt( dble((2*l+1)*(2*l2+1))/(twopi*dble(2*l1+1)))
         tmp    = coeff*Clebsch(l,l2,l1, m,m2,m1)
     >                 *Clebsch(l,l2,l1, 0,0, 0)
      end if

      nwpw_Taunt = tmp
      return
      end

