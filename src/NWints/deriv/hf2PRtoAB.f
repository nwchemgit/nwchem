      Subroutine hf2PRtoAB(dP,dR,dA,dB,alpha,ipair,ff,NPP,NPQ,Nint3,
     &       ictra,ictrb)
c $Id: hf2PRtoAB.f,v 1.4 1994-06-07 00:31:15 d3e129 Exp $

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      integer ictra,ictrb

c--> Derivative Integrals WRT (P,R)

      Dimension dP((NPQ*NPP),Nint3),dR((NPQ*NPP),Nint3)

c--> Derivative Integrals WRT (A,B)

      Dimension dA(Nint3),dB(Nint3)

c--> Exponents & Pair Index

      Dimension alpha(2,NPP),ipair(2,NPP)

c--> Scratch space

      Dimension ff(2,(NPQ*NPP))

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

      do 10 nn = 1,Nint3
       dA(nn) = 0.D0
       dB(nn) = 0.D0
   10 continue

      if (ictra.eq.ictrb) then
        do 00100 nn=1,Nint3
          do 00200 mr = 1,(NPQ*NPP)
            dA(nn) = dA(nn) + dP(mr,nn)
00200     continue
00100   continue
      else
c Compute exponent ratios.

        mr = 0
        
        do 25 mp = 1,NPP
          
          ap = alpha(1,mp)/(alpha(1,mp) + alpha(2,mp))
          bp = alpha(2,mp)/(alpha(1,mp) + alpha(2,mp))
          
          do 20 mq = 1,NPQ
            
            mr = mr + 1
            
            ff(1,mr) = ap
            ff(2,mr) = bp
            
20        continue
          
25      continue
        
c Transform.
        
        do 40 nn = 1,Nint3
          
          do 30 mr = 1,(NPQ*NPP)
            dA(nn) = dA(nn) + (ff(1,mr)*dP(mr,nn) + dR(mr,nn))
            dB(nn) = dB(nn) + (ff(2,mr)*dP(mr,nn) - dR(mr,nn))
30        continue
          
40      continue
      endif
      end
