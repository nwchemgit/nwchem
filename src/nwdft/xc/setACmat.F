      Subroutine setACmat(Amat,Cmat,qwght,Vx,Vc,delrho,
     &                    ipol,nExc,nVxc,nq,GRAD)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical GRAD

c--> Sampling Matrices for the XC Potential & Energy

      Dimension Amat(ipol+nExc,nq),Cmat(ipol,3,nq)

c--> Quadrature Weights

      Dimension qwght(nq)

c--> XC Potential & Energy Partial Derivatives

      Dimension Vx(nVxc,nq),Vc(nVxc,nq)

c--> Charge Density Gradient

      Dimension delrho(ipol,3,nq)
c
c Define the sampling matrices for the XC potential and energy.
c
c*******************************************************************************

      if( GRAD )then

       if( nExc.eq.1 )then

        if( ipol.eq.1 )then

         do 110 n = 1,nq

          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)

          g = 0.5D0*(Vx(3,n) + Vc(3,n))*qwght(n)
        
          Cmat(1,1,n) = g*delrho(1,1,n)
          Cmat(1,2,n) = g*delrho(1,2,n)
          Cmat(1,3,n) = g*delrho(1,3,n)

  110    continue

        else
 
         do 120 n = 1,nq

          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)
          Amat(3,n) = (Vx(3,n) + Vc(3,n))*qwght(n)

          gaa = 2.D0*(Vx(4,n) + Vc(4,n))*qwght(n)
          gbb = 2.D0*(Vx(5,n) + Vc(5,n))*qwght(n)
          gab =      (Vx(6,n) + Vc(6,n))*qwght(n)
        
          Cmat(1,1,n) = gaa*delrho(1,1,n) + gab*delrho(2,1,n)
          Cmat(2,1,n) = gbb*delrho(2,1,n) + gab*delrho(1,1,n)
          Cmat(1,2,n) = gaa*delrho(1,2,n) + gab*delrho(2,2,n)
          Cmat(2,2,n) = gbb*delrho(2,2,n) + gab*delrho(1,2,n)
          Cmat(1,3,n) = gaa*delrho(1,3,n) + gab*delrho(2,3,n)
          Cmat(2,3,n) = gbb*delrho(2,3,n) + gab*delrho(1,3,n)

  120    continue

        end if

       else

        if( ipol.eq.1 )then

         do 130 n = 1,nq

          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n)          )*qwght(n)
          Amat(3,n) = (          Vc(2,n))*qwght(n)

          g = 0.5D0*(Vx(3,n) + Vc(3,n))*qwght(n)
        
          Cmat(1,1,n) = g*delrho(1,1,n)
          Cmat(1,2,n) = g*delrho(1,2,n)
          Cmat(1,3,n) = g*delrho(1,3,n)

  130    continue

        else
 
         do 140 n = 1,nq

          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)
          Amat(3,n) = (Vx(3,n)          )*qwght(n)
          Amat(4,n) = (          Vc(3,n))*qwght(n)

          gaa = 2.D0*(Vx(4,n) + Vc(4,n))*qwght(n)
          gbb = 2.D0*(Vx(5,n) + Vc(5,n))*qwght(n)
          gab =      (Vx(6,n) + Vc(6,n))*qwght(n)
        
          Cmat(1,1,n) = gaa*delrho(1,1,n) + gab*delrho(2,1,n)
          Cmat(2,1,n) = gbb*delrho(2,1,n) + gab*delrho(1,1,n)
          Cmat(1,2,n) = gaa*delrho(1,2,n) + gab*delrho(2,2,n)
          Cmat(2,2,n) = gbb*delrho(2,2,n) + gab*delrho(1,2,n)
          Cmat(1,3,n) = gaa*delrho(1,3,n) + gab*delrho(2,3,n)
          Cmat(2,3,n) = gbb*delrho(2,3,n) + gab*delrho(1,3,n)

  140    continue

        end if

       end if

      else

       if( nExc.eq.1 )then

        if( ipol.eq.1 )then

         do 210 n = 1,nq
          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)
  210    continue

        else
 
         do 220 n = 1,nq
          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)
          Amat(3,n) = (Vx(3,n) + Vc(3,n))*qwght(n)
  220    continue

        end if

       else

        if( ipol.eq.1 )then

         do 230 n = 1,nq
          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n)          )*qwght(n)
          Amat(3,n) = (          Vc(2,n))*qwght(n)
  230    continue

        else
 
         do 240 n = 1,nq
          Amat(1,n) = (Vx(1,n) + Vc(1,n))*qwght(n)
          Amat(2,n) = (Vx(2,n) + Vc(2,n))*qwght(n)
          Amat(3,n) = (Vx(3,n)          )*qwght(n)
          Amat(4,n) = (          Vc(3,n))*qwght(n)
  240    continue

        end if

       end if

      end if

      return
      end
