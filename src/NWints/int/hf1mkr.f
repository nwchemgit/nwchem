      Subroutine hf1mkr(Axyz,Bxyz,Cxyz,zan,ncenters,
     &                  alpha,Pxyz,RS,PC,ff,R,
     &                  R0,R0C,IJK,NPP,Lp,Lp3,CENTER)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical CENTER

      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)

c--> Cartesian Coordinates

      Dimension Axyz(3),Bxyz(3)

c--> Nuclear Cartesian Coordinates & Charges

      Dimension Cxyz(3,ncenters),zan(ncenters)

c--> Exponents

      Dimension alpha(2,NPP)

c--> Auxiliary Function Integrals & Index

      Dimension R0(NPP,Lp3),R0C(ncenters,NPP,Lp3),IJK(0:Lp,0:Lp,0:Lp)

c--> Scratch Space

      Dimension Pxyz(3,NPP),RS(NPP),PC(NPP,3),ff(2,NPP),R(NPP,0:Lp,Lp3)

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

      do 100 mp = 1,NPP*Lp3
       R0(mp,1) = 0.D0
  100 continue

      do 110 mp = 1,NPP

c Define the center "P".

       a = alpha(1,mp)
       b = alpha(2,mp)

       f1 = a/(a+b)
       f2 = b/(a+b)

       Pxyz(1,mp) = f1*Axyz(1) + f2*Bxyz(1)
       Pxyz(2,mp) = f1*Axyz(2) + f2*Bxyz(2)
       Pxyz(3,mp) = f1*Axyz(3) + f2*Bxyz(3)

c Define the scaling factor.

       RS(mp) = sqrt((a+b)*PI4)

  110 continue

c Sum over all centers.

      do 150 ic = 1,ncenters

c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.

       do 120 m = 1,NPP

        alpha_t = alpha(1,m) + alpha(2,m)

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
