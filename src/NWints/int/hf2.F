      Subroutine hf2(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,
     &               Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,
     &               bERI,Nint,canAB,canCD,canPQ,
     &               DryRun,W0,maxW0)
c $Id: hf2.F,v 1.10 1999-12-28 19:49:40 mg201 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical canAB,canCD,canPQ

      Logical GenCon,DryRun

      Logical doit_PQ
c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
      Dimension Cxyz(3),Cprims(NPC),Ccoefs(NPC,NCC)
      Dimension Dxyz(3),Dprims(NPD),Dcoefs(NPD,NCD)
c
c--> Block of Electron Repulsion Integrals

      Dimension bERI(Nint)

c--> Scratch Space

      Dimension W0(maxW0)
c
c Compute 4-ctr electron repulsion integrals (ERI) for four shells of 
c contracted Gaussian functions.
c
c******************************************************************************

      MXD = 0

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB
      NCQ = NCC*NCD

      GenCon = (NCP.ne.1) .or. (NCQ.ne.1)

      if( GenCon.and.(.not.dryrun) )then
       write(*,*) 'HF2:  Not prepared to handle general contractions.'
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
 
c Assign pointers to scratch space.

      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1

      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp

      if( i_left.ge.i_right .and. .not.DryRun )then

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

      if( DryRun )then

       MaxMem = i_left + (maxW0 - (i_right-1))
       NPP = NPA*NPB

      else

       call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &            Bxyz,Bprims,Bcoefs,NPB,NCB,
     &            GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)

      end if
 
c Define the prefactor of the overlap distribution "Q".

c Assign pointers to scratch space.

      i_ALPHAq = i_IPAIRp + 2*(NPA*NPB)
      i_IPAIRq = i_ALPHAq + 2*(NPC*NPD)
      i_left   = i_IPAIRq + 2*(NPC*NPD) - 1

      i_ESq   = i_right - 3*(NPC*NPD)
      i_right = i_ESq

      if( i_left.ge.i_right .and. .not.DryRun )then

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

      if( DryRun )then

       MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )
       NPQ = NPC*NPD

      else

       call hfset(Cxyz,Cprims,Ccoefs,NPC,NCC,
     &            Dxyz,Dprims,Dcoefs,NPD,NCD,
     &            GenCon,W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ESq),NPQ)

       end if

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      NPR = NPP*NPQ

      if (NPR.eq.0 .and. .not.DryRun) then
c... one of the charge distribution either P or Q has no overlap
c... therefore there are no integrals!!!!!!
        lszp = (la+1)*(la+2)*(lb+1)*(lb+2)/4
        lszq = (lc+1)*(lc+2)*(ld+1)*(ld+2)/4
        ninti = lszp*lszq
	call dfill(ninti,0.0d00,beri,1)
        return
      endif

      i_Ep   = i_IPAIRq + 2*(NPC*NPD)
      i_Eq   = i_Ep     + 3*NPP*(MXD+1)*((Lp+1)*(La+1)*(Lb+1))
      i_pf   = i_Eq     + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_left = i_pf     + 2*max(NPP,NPQ) - 1

      if( i_left.ge.i_right .and. .not.DryRun )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '      allocated ',maxW0

       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'Eq    :  ',i_Eq
       write(*,*) 'pf    :  ',i_pf
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp 
       write(*,*) 'ESq   :  ',i_ESq 

       stop

      end if

      if( DryRun )then

       MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )

      else

       do 100 nd = 0,MXD

        call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &             nd,NPP,MXD,La,Lb)

        call hfmke(Cxyz,Dxyz,W0(i_ALPHAq),W0(i_ESq),W0(i_Eq),W0(i_pf),
     &             nd,NPQ,MXD,Lc,Ld)

  100  continue

      end if

c Evaluate the auxiliary function integrals.

c Assign pointers to scratch space.

      i_R0  = i_Eq  + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_IJK = i_R0  + NPR*Lr3
      i_P   = i_IJK + (Lr+1)**3
      i_Q   = i_P   + 3*NPP
      i_PQ  = i_Q   + 3*NPQ
      i_ff  = i_PQ  + NPR*3
      i_Rj  = i_ff  + 2*NPR
      i_top = i_Rj  + NPR*(Lr+1)*Lr3 - 1

      if( i_top.gt.maxW0 .and. .not.DryRun )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_top
       write(*,*) '      allocated ',maxW0

       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
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

*rak:      doit_PQ = NPP.ge.NPQ

      if (Lp.eq.Lq) then
        doit_PQ = NPP.ge.NPQ
      else if (abs((NPP-NPQ)).ge.4) then
        doit_PQ = NPP.ge.NPQ
      else if (Lp.lt.Lq) then
        doit_PQ = .true.
      else
        doit_PQ = .false.
      endif

      if( DryRun )then

       MaxMem = max( MaxMem, i_top )

      else

        if (doit_PQ) then
          call hf2mkr(Axyz,Bxyz,Cxyz,Dxyz,W0(i_ALPHAp),W0(i_ALPHAq),
     &        W0(i_R0),W0(i_IJK),W0(i_P),W0(i_Q),W0(i_PQ),
     &        W0(i_ff),W0(i_Rj),
     &        NPP,NPQ,Lr,Lr3)
        else
          call hf2mkr(Cxyz,Dxyz,Axyz,Bxyz,W0(i_ALPHAq),W0(i_ALPHAp),
     &        W0(i_R0),W0(i_IJK),W0(i_Q),W0(i_P),W0(i_PQ),
     &        W0(i_ff),W0(i_Rj),
     &        NPQ,NPP,Lr,Lr3)
        endif

      end if

c Compute the ERI.

c Assign pointers to scratch space.

      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6
      Lq3 = ((Lq+1)*(Lq+2)*(Lq+3))/6
      Lpq3 = max(Lp3,Lq3)

      i_E3   = i_IJK + (Lr+1)**3
      i_sum  = i_E3  + max(NPP,NPQ)
      if (doit_PQ) then
        i_top  = i_sum + NPR*Lpq3 - 1
      else
        lszp = (la+1)*(la+2)*(lb+1)*(lb+2)/4
        lszq = (lc+1)*(lc+2)*(ld+1)*(ld+2)/4
        ninti = lszp*lszq
        i_beri = i_sum + NPR*Lpq3 
        i_top  = i_beri + ninti
      endif
      if( i_top.gt.maxW0 .and. .not.DryRun )then

       write(*,*) 'HF2:  Insufficient scratch space.'
       write(*,*) '      needed    ',i_top
       write(*,*) '      allocated ',maxW0

       write(*,*) 'ALPHAp:  ',i_ALPHAp 
       write(*,*) 'IPAIRp:  ',i_IPAIRp 
       write(*,*) 'ALPHAq:  ',i_ALPHAq 
       write(*,*) 'IPAIRq:  ',i_IPAIRq 
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'Eq    :  ',i_Eq
       write(*,*) 'R0    :  ',i_R0
       write(*,*) 'IJK   :  ',i_IJK
       write(*,*) 'E3    :  ',i_E3
       write(*,*) 'sum   :  ',i_sum

       stop

      end if

      if( DryRun )then

       MaxMem = max( MaxMem, i_top )

      else

       if (doit_PQ) then
         call hferi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),bERI,
     &       W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,MXD,
     &       canAB,canCD,canPQ)
       else
         call hferi(W0(i_Eq),W0(i_Ep),W0(i_R0),W0(i_IJK),W0(i_beri),
     &       W0(i_E3),W0(i_sum),NPQ,NPP,Ninti,Lc,Ld,La,Lb,Lr,MXD,
     &       canAB,canCD,canPQ)
         call hfreord_pq(bERI,W0(i_beri),lszp,lszq)
       endif
      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
00099 continue
      end
