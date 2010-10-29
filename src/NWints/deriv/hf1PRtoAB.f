      Subroutine hf1PRtoAB(dP,dR,dA,dB,alpha,ipair,ff,NPP,Nint,
     &       ictrA,ictrB)
c $Id$

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
c--> locals
      double precision sumA, sumB
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

      if(ictra.eq.ictrb) then
        do 00100 nn = 1,Nint
          sumA = 0.0d00
          do 00200 mp = 1,NPP
            sumA = sumA + dP(mp,nn)
00200     continue
          dA(nn) = sumA
00100   continue
        call dlaset(' ',nint,1,0.0d00,0.0d00,dB,nint)
      else
c Compute exponent ratios.

      do 20 mp = 1,NPP
       ff(1,mp) = alpha(1,mp)/(alpha(1,mp) + alpha(2,mp))
       ff(2,mp) = alpha(2,mp)/(alpha(1,mp) + alpha(2,mp))
   20 continue

c Transform.

      do 40 nn = 1,Nint

        sumA = 0.0d00
        sumB = 0.0d00
        do 30 mp = 1,NPP
          sumA = sumA + (ff(1,mp)*dP(mp,nn) + dR(mp,nn))
          sumB = sumB + (ff(2,mp)*dP(mp,nn) - dR(mp,nn))
 30     continue
        dA(nn) = sumA
        dB(nn) = sumB
 40   continue
      endif

      end
