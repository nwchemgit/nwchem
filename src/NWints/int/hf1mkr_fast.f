      Subroutine hf1mkr_fast(Axyz,Bxyz,Cxyz,zan,ncenters,
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
      Double precision alpha(2,NPP)     ! basis function exponents copy
      Double precision Pxyz(3,NPP)      ! coordinates of product center
      Double precision RS(NPP)          ! scale factor
      Double precision PC(NPP,3)        ! distance from product center to nucleus
      Double precision ff(NPP,3)        ! recursion factors
      Double precision R(NPP,0:Lp,Lp3)  ! R/R000j/RIKJj functions
      Double precision R0(NPP,Lp3)      ! R integrals
      Double precision R0C(ncenters,NPP,Lp3) ! R integrals per center
      integer IJK(0:Lp,0:Lp,0:Lp)       ! index array 
      Logical CENTER
c::local
      integer ic,j,m,n
      Double precision a,b,f1,f2,alpha_t,PCx,PCy,PCz
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
       alpha_t = a + b
       ff(m,2) = -2.D0*alpha_t

       f1 = a/(a+b)
       f2 = b/(a+b)

       Pxyz(1,m) = f1*Axyz(1) + f2*Bxyz(1)
       Pxyz(2,m) = f1*Axyz(2) + f2*Bxyz(2)
       Pxyz(3,m) = f1*Axyz(3) + f2*Bxyz(3)

c Define the scaling factor.

       RS(m) = sqrt((a+b)*PI4)


 110  continue

c Sum over all centers.


      do 150 ic = 1,ncenters


c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.

       do 120 m = 1,NPP


        ff(m,1) = RS(m)

        PCx = Pxyz(1,m) - Cxyz(1,ic)
        PCy = Pxyz(2,m) - Cxyz(2,ic)
        PCz = Pxyz(3,m) - Cxyz(3,ic)

        R(m,0,1) = -0.5*ff(m,2)*(PCx**2 + PCy**2 + PCz**2)

        PC(m,1) = PCx
        PC(m,2) = PCy
        PC(m,3) = PCz

  120  continue

c Evaluate the incomplete gamma function.

       call igamma(R,NPP,Lp)

c Scale for finite nuclear size if necessary


c Define the initial auxiliary functions (i.e., R000j, j=1,Lr).

       if(lp.eq.0) then
        do  m = 1,NPP
         R(m,0,1) = ff(m,1)*R(m,0,1)
         ff(m,1) = ff(m,1)*ff(m,2)
      enddo
       else
       do 135 j = 0,Lp
        do 130 m = 1,NPP
         R(m,j,1) = ff(m,1)*R(m,j,1)
         ff(m,1) = ff(m,1)*ff(m,2)
  130   continue
  135  continue
       endif
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
         if(lp3.eq.1) then
            if (NPP.eq.1) then
                  R0(1,1) = R0(1,1) - zan(ic)*R(1,0,1)
               else
                  do m = 1,NPP
                     R0(m,1) = R0(m,1) - zan(ic)*R(m,0,1)
                  enddo
               endif
         else
            if(NPP.eq.1) then
               do  n = 1,Lp3
                  R0(1,n) = R0(1,n) - zan(ic)*R(1,0,n)
               enddo
            else
               do n = 1,Lp3
                  do m = 1,NPP
                     R0(m,n) = R0(m,n) - zan(ic)*R(m,0,n)
                  enddo
               enddo
            endif
         endif
       end if

  150 continue

      end
