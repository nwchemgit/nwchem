      Subroutine hf1PRtoAB(dP,dR,dA,dB,alpha,ipair,ff,NPP,Nint)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

c--> Derivative integrals wrt to (P,R).

      Dimension dP(NPP,Nint),dR(NPP,Nint)

c--> Derivative integrals wrt to (A,B).

      Dimension dA(Nint),dB(Nint)

c--> Exponents & Pair Index

      Dimension alpha(2,NPP),ipair(2,NPP)

c--> Scratch space.

      Dimension ff(2,NPP)
c
c Transform derivative integrals wrt (P,R) to (A,B).
c
c N.B. It is assumed that the product of contraction coefficients has been
c      factored into each primitive (P,R) integral derivative. Thus, this 
c      routine currently transforms primitive (P,R) integral derivatives to 
c      contracted (A,B) integral derivatives.
c
c*******************************************************************************

c Initialize derivative integrals wrt to (A,B).

      do 10 nn = 1,Nint
       dA(nn) = 0.D0
       dB(nn) = 0.D0
   10 continue

c Compute exponent ratios.

      do 20 mp = 1,NPP
       ff(1,mp) = alpha(1,mp)/(alpha(1,mp) + alpha(2,mp))
       ff(2,mp) = alpha(2,mp)/(alpha(1,mp) + alpha(2,mp))
   20 continue

c Transform.

      do 40 nn = 1,Nint

       do 30 mp = 1,NPP
        dA(nn) = dA(nn) + (ff(1,mp)*dP(mp,nn) + dR(mp,nn))
        dB(nn) = dB(nn) + (ff(2,mp)*dP(mp,nn) - dR(mp,nn))
   30  continue

   40 continue

      end
