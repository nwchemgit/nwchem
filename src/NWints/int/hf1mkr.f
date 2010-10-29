      Subroutine hf1mkr(Axyz,Bxyz,Cxyz,zan,exinv,ncenters,
     &                  alpha,Pxyz,RS,PC,ff,R,
     &                  R0,R0C,IJK,NPP,Lp,Lp3,CENTER)
c $Id$

      Implicit none
c::passed
      Integer ncenters                  ! number of nuclei
      Integer NPP                       ! number of primitive products
      Integer Lp                        ! angular momentum sum La+Lb+MXD
      Integer Lp3                       ! number of angular momentum functions
      Double precision Axyz(3),Bxyz(3)  ! Cartesian coordinates of centers 
      Double precision Cxyz(3,ncenters) ! Nuclear cartesian coordinates
      Double precision zan(ncenters)    ! Nuclear charges
      Double precision exinv(ncenters)  ! Nuclear inverse exponents
      Double precision alpha(2,NPP)     ! basis function exponents copy
      Double precision Pxyz(3,NPP)      ! coordinates of product center
      Double precision RS(NPP)          ! scale factor
      Double precision PC(NPP,3)        ! distance from product center to nucleus
      Double precision ff(2,NPP)        ! recursion factors
      Double precision R(NPP,0:Lp,Lp3)  ! R/R000j/RIKJj functions
      Double precision R0(NPP,Lp3)      ! R integrals
      Double precision R0C(ncenters,NPP,Lp3) ! R integrals per center
      integer IJK(0:Lp,0:Lp,0:Lp)       ! index array 
      Logical CENTER
c::local
      integer ic,j,m,n
      Double precision a,b,f1,f2,alpha_t,beta,PCx,PCy,PCz
      Double precision PI, PI4          ! pi, 4/pi
      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)
c
c Define the auxiliary function integrals necessary to compute the nuclear
c attraction integrals (NAIs). These integrals are scaled by an appropriate
c factor, RS, defined as
c
c         / a + b \ 1/2
c   RS = | ------- |
c         \  PI/4 /
c
c The scale factor for the Hermite expansion coefficients is assumed to be
c
c         /   PI  \ 3/2     /   a b   __2 \
c   ES = | ------- |    EXP| - -----  AB   |
c         \ a + b /         \  a + b      /
c
c Therefore,
c
c            2 PI        /   a b   __2 \
c   ES RS = -------  EXP| - -----  AB   |
c            a + b       \  a + b      /
c
c******************************************************************************

      do 101 m = 1,NPP
       do 100 j = 1,Lp3
        R0(m,j) = 0.D0
  100  continue
  101 continue

      do 110 m = 1,NPP

c Define the center "P".

       a = alpha(1,m)
       b = alpha(2,m)

       f1 = a/(a+b)
       f2 = b/(a+b)

       Pxyz(1,m) = f1*Axyz(1) + f2*Bxyz(1)
       Pxyz(2,m) = f1*Axyz(2) + f2*Bxyz(2)
       Pxyz(3,m) = f1*Axyz(3) + f2*Bxyz(3)

c Define the scaling factor.

       RS(m) = sqrt((a+b)*PI4)

  110 continue

c Sum over all centers.

      do 150 ic = 1,ncenters

       beta = exinv(ic)

c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.

       do 120 m = 1,NPP

        alpha_t = alpha(1,m) + alpha(2,m)
        alpha_t = alpha_t/(1.0d0+alpha_t*beta)

        ff(1,m) = RS(m)
        ff(2,m) = -2.D0*alpha_t

        PCx = Pxyz(1,m) - Cxyz(1,ic)
        PCy = Pxyz(2,m) - Cxyz(2,ic)
        PCz = Pxyz(3,m) - Cxyz(3,ic)

        R(m,0,1) = alpha_t*(PCx**2 + PCy**2 + PCz**2)

        PC(m,1) = PCx
        PC(m,2) = PCy
        PC(m,3) = PCz

  120  continue

c Evaluate the incomplete gamma function.

       call igamma(R,NPP,Lp)

c Scale for finite nuclear size if necessary

       if (beta .ne. 0.0D0) then
        do 125 m = 1,NPP
         alpha_t = sqrt(1.0d0+(alpha(1,m)+alpha(2,m))*beta)
         do 126 j=0,Lp
          R(m,j,1) = R(m,j,1)/alpha_t
  126    continue
  125   continue
       end if

c Define the initial auxiliary functions (i.e., R000j, j=1,Lr).

       do 135 j = 0,Lp
        do 130 m = 1,NPP
         R(m,j,1) = ff(1,m)*R(m,j,1)
         ff(1,m) = ff(1,m)*ff(2,m)
  130   continue
  135  continue

c Recursively build the remaining auxiliary functions (i.e., RIJKj, j=0).

       call hfmkr(R,IJK,PC,NPP,Lp,Lp3)

c Transfer to R0 array.

       if( CENTER )then
        do 141 n = 1,Lp3
         do 140 m = 1,NPP
          R0C(ic,m,n) = -zan(ic)*R(m,0,n)
          R0(m,n) = R0(m,n) + R0C(ic,m,n)
  140    continue
  141   continue
       else
        do 146 n = 1,Lp3
         do 145 m = 1,NPP
          R0(m,n) = R0(m,n) - zan(ic)*R(m,0,n)
  145    continue
  146   continue
       end if

  150 continue

      end
