      Subroutine hf2d(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &                Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &                Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,
     &                Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,
     &                bERI,Nint,canAB,canCD,canPQ,W0,maxW0)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB,canCD,canPQ

      Logical GENCON

#if defined(MSS_TIMER)
#include <timer.fh>
#endif

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
      Dimension Cxyz(3),Cprims(NPC),Ccoefs(NPC,NCC)
      Dimension Dxyz(3),Dprims(NPD),Dcoefs(NPD,NCD)

c--> Block of Electron Repulsion Integrals

      Dimension bERI(Nint,*)

c--> Scratch Space

      Dimension W0(maxW0)

c--> Derivative Indices

      Parameter (mxder=3)
      Dimension inder1(6,mxder),inder2(3,mxder),inder3(3,mxder)
c
c Compute gradients of 4-ctr electron repulsion integrals (ERI) for 4 shells 
c of contracted Gaussian functions.
c
c******************************************************************************

      MXD = 1

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB
      NCQ = NCC*NCD

      GENCON = (NCP.ne.1) .or. (NCQ.ne.1)

      if( GENCON )then
       write(*,*) 'HF2D:  Not prepared to handle general contractions.'
       stop
      end if

c Define the angular momentum of the overlap distributions.

      Lp = La + Lb
      Lq = Lc + Ld
      Lr = Lp + Lq

c Increment "Lr" to account for the order of differentiation.

      Lr = Lr + MXD

c Define the accumulated number of angular momentum functions <= Lr.

      Lr3 = ((Lr+1)*(Lr+2)*(Lr+3))/6

c Define the prefactor of the overlap distribution "P".
 
#if defined(MSS_TIMER)
                     call second(t1)
#endif

c Assign pointers to scratch space.

      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1

      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp

      if( i_left.ge.i_right )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '      allocated ',maxW0

       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp 

       stop

      end if

      call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &           Bxyz,Bprims,Bcoefs,NPB,NCB,
     &           GENCON,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)
 
c Define the prefactor of the overlap distribution "Q".

c Assign pointers to scratch space.

      i_ALPHAq = i_IPAIRp + 2*(NPA*NPB)
      i_IPAIRq = i_ALPHAq + 2*(NPC*NPD)
      i_left   = i_IPAIRq + 2*(NPC*NPD) - 1

      i_ESq   = i_right - 3*(NPC*NPD)
      i_right = i_ESq

      if( i_left.ge.i_right )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '      allocated ',maxW0

       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp 
       write(*,*) 'ESq   :  ',i_ESq 

       stop

      end if

      call hfset(Cxyz,Cprims,Ccoefs,NPC,NCC,
     &           Dxyz,Dprims,Dcoefs,NPD,NCD,
     &           GENCON,W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ESq),NPQ)

#if defined(MSS_TIMER)
                     call second(t2)
                     if(itime.eq.1) time(1) = time(1) + (t2-t1)
#endif

c Define the Hermite linear expansion coefficients.

#if defined(MSS_TIMER)
                     call second(t1)
#endif

c Assign pointers to scratch space.

      NPR = NPP*NPQ

      i_ERI  = i_IPAIRq + 2*(NPC*NPD)
      i_Ep   = i_ERI    + NPR*(Nint*12)
      i_Eq   = i_Ep     + 3*NPP*(MXD+1)*((Lp+1)*(La+1)*(Lb+1))
      i_pf   = i_Eq     + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_left = i_pf     + 2*max(NPP,NPQ) - 1

      if( i_left.ge.i_right )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '      allocated ',maxW0

       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'ERI   :  ',i_ERI
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'Eq    :  ',i_Eq
       write(*,*) 'pf    :  ',i_pf
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp 
       write(*,*) 'ESq   :  ',i_ESq 

       stop

      end if

      do 100 nd = 0,MXD

       call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &            nd,NPP,MXD,La,Lb)

       call hfmke(Cxyz,Dxyz,W0(i_ALPHAq),W0(i_ESq),W0(i_Eq),W0(i_pf),
     &            nd,NPQ,MXD,Lc,Ld)

  100 continue

#if defined(MSS_TIMER)
                     call second(t2)
                     if(itime.eq.1) time(2) = time(2) + (t2-t1)
#endif

c Evaluate the auxiliary function integrals.

#if defined(MSS_TIMER)
                     call second(t1)
#endif

c Assign pointers to scratch space.

      i_R0  = i_Eq  + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_IJK = i_R0  + NPR*Lr3
      i_P   = i_IJK + (Lr+1)**3
      i_Q   = i_P   + 3*NPP
      i_PQ  = i_Q   + 3*NPQ
      i_ff  = i_PQ  + NPR*3
      i_Rj  = i_ff  + 2*NPR
      i_top = i_Rj  + NPR*(Lr+1)*Lr3 - 1

      if( i_top.gt.maxW0 )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_top
       write(*,*) '      allocated ',maxW0

       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'ERI   :  ',i_ERI
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'Eq    :  ',i_Eq
       write(*,*) 'R0    :  ',i_R0
       write(*,*) 'IJK   :  ',i_IJK
       write(*,*) 'P     :  ',i_P
       write(*,*) 'Q     :  ',i_Q
       write(*,*) 'PQ    :  ',i_PQ
       write(*,*) 'ff    :  ',i_ff
       write(*,*) 'Rj    :  ',i_Rj

       stop

      end if

      call hf2mkr(Axyz,Bxyz,Cxyz,Dxyz,W0(i_ALPHAp),W0(i_ALPHAq),
     &            W0(i_R0),W0(i_IJK),W0(i_P),W0(i_Q),W0(i_PQ),
     &            W0(i_ff),W0(i_Rj),NPP,NPQ,Lr,Lr3)

#if defined(MSS_TIMER)
                     call second(t2)
                     if(itime.eq.1) time(3) = time(3) + (t2-t1)
#endif

c Compute the gradient ERIs.
c
c Currently, the 12 gradient integrals can be grouped into 3 sets. The 
c integrals in each set can be evaluated with a single call to the 
c routine "hfderi". The derivative indices specify the order of 
c differentiation for each of the 12 coordinates:
c
c     Set 1, inder1(i,n), i=1-3:   d/dPx, d/dPy, d/dPz
c                         i=4-6:   d/dQx, d/dQy, d/dQz
c     Set 2, inder2(i,n), i=1-3:   d/dRx, d/dRy, d/dRz
c     Set 3, inder3(i,n), i=1-3:   d/dSx, d/dSy, d/dSz
c
c N.B. The routine "hfderi" cannot simultaneously evaluate integrals from
c      different sets!
c
#if defined(MSS_TIMER)
                     call second(t1)
#endif

c Assign pointers to scratch space.

      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6
      Lq3 = ((Lq+1)*(Lq+2)*(Lq+3))/6

      i_E3   = i_IJK + (Lr+1)**3
      i_sum  = i_E3  + max(NPP,NPQ)
      i_top  = i_sum + NPR*6*Lq3 - 1

      if( i_top.gt.maxW0 )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_top
       write(*,*) '      allocated ',maxW0

       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'ERI   :  ',i_ERI
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'Eq    :  ',i_Eq
       write(*,*) 'R0    :  ',i_R0
       write(*,*) 'IJK   :  ',i_IJK
       write(*,*) 'E3    :  ',i_E3
       write(*,*) 'sum   :  ',i_sum

       stop

      end if

c Integral derivatives wrt P.

      i_dP = i_ERI

c Define the derivative indices.

      Nder1 = 3
      inder1(1,1) = 1
      inder1(2,1) = 0
      inder1(3,1) = 0
      inder1(4,1) = 0
      inder1(5,1) = 0
      inder1(6,1) = 0
      inder1(1,2) = 0
      inder1(2,2) = 1
      inder1(3,2) = 0
      inder1(4,2) = 0
      inder1(5,2) = 0
      inder1(6,2) = 0
      inder1(1,3) = 0
      inder1(2,3) = 0
      inder1(3,3) = 1
      inder1(4,3) = 0
      inder1(5,3) = 0
      inder1(6,3) = 0

      Nder2 = 1
      inder2(1,1) = 0
      inder2(2,1) = 0
      inder2(3,1) = 0

      Nder3 = 1
      inder3(1,1) = 0
      inder3(2,1) = 0
      inder3(3,1) = 0

      call hfderi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dP),
     &            W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &            MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &            canAB,canCD,canPQ)

c Integral derivatives wrt Q.

      i_dQ = i_dP + NPR*(Nint*3)

c Translation invariance requires d/dQ = -d/dP !

      do 200 n = 1,NPR*(Nint*3)
       W0(i_dQ+(n-1)) = -W0(i_dP+(n-1))
  200 continue

c Integral derivatives wrt R.

      i_dR = i_dQ + NPR*(Nint*3)

c Define the derivative indices.

      Nder1 = 1
      inder1(1,1) = 0
      inder1(2,1) = 0
      inder1(3,1) = 0
      inder1(4,1) = 0
      inder1(5,1) = 0
      inder1(6,1) = 0

      Nder2 = 3
      inder2(1,1) = 1
      inder2(2,1) = 0
      inder2(3,1) = 0
      inder2(1,2) = 0
      inder2(2,2) = 1
      inder2(3,2) = 0
      inder2(1,3) = 0
      inder2(2,3) = 0
      inder2(3,3) = 1

      Nder3 = 1
      inder3(1,1) = 0
      inder3(2,1) = 0
      inder3(3,1) = 0

      call hfderi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dR),
     &            W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &            MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &            canAB,canCD,canPQ)

c Integral derivatives wrt S.

      i_dS = i_dR + NPR*(Nint*3)

c Define the derivative indices.

      Nder1 = 1
      inder1(1,1) = 0
      inder1(2,1) = 0
      inder1(3,1) = 0
      inder1(4,1) = 0
      inder1(5,1) = 0
      inder1(6,1) = 0

      Nder2 = 1
      inder2(1,1) = 0
      inder2(2,1) = 0
      inder2(3,1) = 0

      Nder3 = 3
      inder3(1,1) = 1
      inder3(2,1) = 0
      inder3(3,1) = 0
      inder3(1,2) = 0
      inder3(2,2) = 1
      inder3(3,2) = 0
      inder3(1,3) = 0
      inder3(2,3) = 0
      inder3(3,3) = 1

      call hfderi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dS),
     &            W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &            MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &            canAB,canCD,canPQ)

#if defined(MSS_TIMER)
                     call second(t2)
                     if(itime.eq.1) time(4) = time(4) + (t2-t1)
#endif

c Contract primitive integrals.

#if defined(MSS_TIMER)
                     call second(t1)
#endif

c Assign pointers to scratch space.

      i_ff  = i_ERI + NPR*(Nint*12)
      i_top = i_ff  + NPR*2 - 1

      if( i_top.gt.maxW0 )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_top
       write(*,*) '      allocated ',maxW0

       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'ERI   :  ',i_ERI
       write(*,*) 'ff    :  ',i_ff

       stop

      end if

c Transform to obtain derivatives wrt (A,B,C,D).

#if defined(MSS_TIMER)
                     call second(t1)
#endif

      call hf2PRtoAB(W0(i_dP),W0(i_dR),bERI(1,1),bERI(1,4),
     &               W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ff),NPP,NPQ,Nint*3)

      call hf2QStoCD(W0(i_dQ),W0(i_dS),bERI(1,7),bERI(1,10),
     &               W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ff),NPP,NPQ,Nint*3)

#if defined(MSS_TIMER)
                     call second(t2)
                     if(itime.eq.1) time(5) = time(5) + (t2-t1)
#endif

      end
