c $Id$
      Subroutine hf2(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,
     &               Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,
     &               bERI,Nints,canAB,canCD,canPQ,
     &               DryRun,W0,maxW0)

      Implicit none
#include "stdio.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "case.fh"

      Integer NPA,NCA,La
      Integer NPB,NCB,Lb
      Integer NPC,NCC,Lc
      Integer NPD,NCD,Ld
      Integer Nints,maxW0

      Logical canAB,canCD,canPQ

      Logical GenCon,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Double Precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Double Precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
      Double Precision Cxyz(3),Cprims(NPC),Ccoefs(NPC,NCC)
      Double Precision Dxyz(3),Dprims(NPD),Dcoefs(NPD,NCD)
c
c--> Block of Electron Repulsion Integrals

      Double Precision bERI(Nints)

c--> Scratch Space

      Double Precision W0(maxW0)

c--> Local variables
      
      Integer NCP,NCQ,MXD,nd
      Integer La2,Lb2,Lc2,Ld2
      Integer lp,lq,lr,lp3,lq3,lr3,lpq3
      Integer i_ALPHAp,i_IPAIRp,i_ESp,i_Ep
      Integer i_ALPHAq,i_IPAIRq,i_ESq,i_Eq
      Integer i_left,i_right,i_pf
      Integer MaxMem,MaxAvail,MemPerPair,MaxPairs,NPass
      Integer maxpp,maxpq,ipp,ipq
      Integer MPP,MPQ,MPR,NPP,NPQ,NPR
      Integer lszp,lszq,ninti
      Integer i_R0,i_IJK,i_P,i_Q,i_PQ,i_ff,i_Rj
      Integer i_E3,i_t1,i_t2,i_t3,i_t4,i_eri,i_top
      Integer i,j,istep
      integer i_ffscr,l_ffscr,i_rscr,l_rscr,nvals

      Logical doit_PQ

c
c Compute 4-ctr electron repulsion integrals (ERI) for four shells of 
c contracted Gaussian functions.
c
c******************************************************************************

      MXD = 0
      istep = MA_sizeof(MT_DBL,1,MT_INT)
      if (istep .gt. 2) call errquit(
     &    'hf2:Too many integers per real for Ipair addressing',911,
     &       INT_ERR)
      istep = 2/istep

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB
      NCQ = NCC*NCD

      GenCon = (NCP.ne.1) .or. (NCQ.ne.1)

c Define the angular momentum of the overlap distributions.

      Lp = La + Lb
      Lq = Lc + Ld
      Lr = Lp + Lq

c Increment "Lr" to account for the order of differentiation.

      Lr = Lr + MXD

c Define the accumulated number of angular momentum functions <= Lr.

      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6
      Lq3 = ((Lq+1)*(Lq+2)*(Lq+3))/6
      Lpq3 = max(Lp3,Lq3)
      Lr3 = ((Lr+1)*(Lr+2)*(Lr+3))/6

c Define the prefactor of the overlap distribution "P".
 
c Assign pointers to scratch space.

      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp

      MaxMem = 1   ! take care of compiler warnings
      if (DryRun) then
        MaxMem = i_left + (maxW0 - (i_right-1))
        NPP = NPA*NPB
      else
        if (i_left.ge.i_right) then
          write(LuOut,*) 'HF2:  Insufficient scratch space.'
          write(LuOut,*) '      needed    ',i_left + (maxW0-(i_right-1))
          write(LuOut,*) '      allocated ',maxW0
          write(LuOut,*) 'From the left '
          write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
          write(LuOut,*) 'IPAIRp:  ',i_IPAIRp 
          write(LuOut,*) 'From the right '
          write(LuOut,*) 'ESp   :  ',i_ESp 
          call errquit('hf2: Not enough scratch space for hfset AB',911,
     &       MEM_ERR)
        end if
        call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,
     &      GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)
      end if
 
c Define the prefactor of the overlap distribution "Q".

c Assign pointers to scratch space.

      i_ALPHAq = i_IPAIRp + 2*(NPA*NPB)
      i_IPAIRq = i_ALPHAq + 2*(NPC*NPD)
      i_left   = i_IPAIRq + 2*(NPC*NPD) - 1
      i_ESq   = i_right - 3*(NPC*NPD)
      i_right = i_ESq

      if( DryRun )then
        MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )
        NPQ = NPC*NPD
      else
        if (i_left.ge.i_right) then
          write(LuOut,*) 'HF2:  Insufficient scratch space.'
          write(LuOut,*) '      needed    ',i_left + (maxW0-(i_right-1))
          write(LuOut,*) '      allocated ',maxW0
          write(LuOut,*) 'From the left '
          write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
          write(LuOut,*) 'IPAIRp:  ',i_IPAIRp 
          write(LuOut,*) 'ALPHAq:  ',i_ALPHAq 
          write(LuOut,*) 'IPAIRq:  ',i_IPAIRq 
          write(LuOut,*) 'From the right '
          write(LuOut,*) 'ESp   :  ',i_ESp 
          write(LuOut,*) 'ESq   :  ',i_ESq 
          call errquit('hf2: Not enough scratch space for hfset CD',911,
     &       MEM_ERR)
        end if
        call hfset(Cxyz,Cprims,Ccoefs,NPC,NCC,
     &      Dxyz,Dprims,Dcoefs,NPD,NCD,
     &      GenCon,W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ESq),NPQ)
      end if

c Zero out integral block. Return if NPR = 0, i.e. no integrals

      NPR = NPP*NPQ
      if (NPR .eq. 0) then 
!     if (.not.DryRun) call dfill(Nints,0.0d00,beri,1)
      if (.not.DryRun) call dcopy(Nints,0d0,0,beri,1) 
         go to 99
      endif

c Define cartesian components for shells, bra and ket length and 
c number of integrals.

      La2 = (la+1)*(la+2)/2
      Lb2 = (lb+1)*(lb+2)/2
      lszp = La2*Lb2*NCP

      Lc2 = (lc+1)*(lc+2)/2
      Ld2 = (ld+1)*(ld+2)/2
      lszq = Lc2*Ld2*NCQ

      ninti = lszp*lszq
      if (.not.DryRun) call dfill(Ninti,0.0d00,beri,1)

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      i_Ep   = i_IPAIRq + 2*(NPC*NPD)
      i_Eq   = i_Ep     + 3*NPP*(MXD+1)*((Lp+1)*(La+1)*(Lb+1))
      i_pf   = i_Eq     + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_left = i_pf     + 2*max(NPP,NPQ) - 1

      if( DryRun )then
        MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )
      else
        if (i_left.ge.i_right) then
          write(LuOut,*) 'HF2:  Insufficient scratch space.'
          write(LuOut,*) '      needed    ',i_left + (maxW0-(i_right-1))
          write(LuOut,*) '      allocated ',maxW0
          write(LuOut,*) 'From the left '
          write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
          write(LuOut,*) 'IPAIRp:  ',i_IPAIRp 
          write(LuOut,*) 'ALPHAq:  ',i_ALPHAq 
          write(LuOut,*) 'IPAIRq:  ',i_IPAIRq 
          write(LuOut,*) 'Ep    :  ',i_Ep
          write(LuOut,*) 'Eq    :  ',i_Eq
          write(LuOut,*) 'pf    :  ',i_pf
          write(LuOut,*) 'From the right '
          write(LuOut,*) 'ESp   :  ',i_ESp 
          write(LuOut,*) 'ESq   :  ',i_ESq 
          call errquit('hf2: Not enough scratch space for hfmke',911,
     &       MEM_ERR)
        end if

        do nd = 0,MXD
          call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &        nd,NPP,MXD,La,Lb)
          call hfmke(Cxyz,Dxyz,W0(i_ALPHAq),W0(i_ESq),W0(i_Eq),W0(i_pf),
     &        nd,NPQ,MXD,Lc,Ld)
        end do

      end if

c Evaluate the auxiliary function integrals.

c   Determine memory requirements and assign pointers to scratch space.
c   If insufficient memory, do multipassing.

      i_IJK = i_Eq  + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_R0  = i_IJK + (Lr+1)**3

      MaxAvail = maxW0-i_R0-3*(NPP+NPQ)+1
      MemPerPair = (Lr+2)*Lr3 + 5
      if (MaxAvail-MemPerPair .le. 0) then
        write(LuOut,*) 'HF2:  Insufficient scratch space.'
        write(LuOut,*) '      needed    ',i_R0+3*(NPP+NPQ)
     &      +MemPerPair*NPR
        write(LuOut,*) '      allocated ',maxW0
        write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
        write(LuOut,*) 'IPAIRp:  ',i_IPAIRp 
        write(LuOut,*) 'ALPHAq:  ',i_ALPHAq 
        write(LuOut,*) 'IPAIRq:  ',i_IPAIRq 
        write(LuOut,*) 'Ep    :  ',i_Ep
        write(LuOut,*) 'Eq    :  ',i_Eq
        write(LuOut,*) 'R0    :  ',i_R0
        write(LuOut,*) 'IJK   :  ',i_IJK
        call errquit('hf2: Not enough scratch space for hf2mkr',911,
     &       MEM_ERR)
      end if

      MaxPairs = Maxavail/MemPerPair
      
      if (MaxPairs .ge. NPR) then
        maxpp = NPP
        maxpq = NPQ
      else if (NPP .ge. NPQ) then
        if (MaxPairs .ge. NPP) then
          maxpp = NPP
          maxpq = MaxPairs/NPP
        else if (MaxPairs .ge. NPQ) then
          maxpq = NPQ
          maxpp = MaxPairs/NPQ
        else
          NPass = (NPP-1)/MaxPairs+1
          maxpp = min((NPP-1)/NPass+1,MaxPairs)
          maxpq = 1
        end if
      else
        if (MaxPairs .ge. NPQ) then
          maxpq = NPQ
          maxpp = MaxPairs/NPQ
        else if (MaxPairs .ge. NPP) then
          maxpp = NPP
          maxpq = MaxPairs/NPP
        else
          NPass = (NPP-1)/MaxPairs+1
          maxpp = min((NPP-1)/NPass+1,MaxPairs)
          maxpq = 1
        end if
      end if

      do ipp = 0,NPP-1,maxpp
        MPP = min(maxpp,NPP-ipp)
        do ipq = 0,NPQ-1,maxpq
          MPQ = min(maxpq,NPQ-ipq)
          MPR = MPP*MPQ
          i_P   = i_R0  + MPR*Lr3
          i_Q   = i_P   + 3*MPP
          i_PQ  = i_Q   + 3*MPQ
          i_ff  = i_PQ  + MPR*3
          i_Rj  = i_ff  + 2*MPR
          i_top = i_Rj  + MPR*(Lr+1)*Lr3 - 1

          if (Lp.eq.Lq) then
            doit_PQ = MPP.ge.MPQ
          else if (abs((MPP-MPQ)).ge.4) then
            doit_PQ = MPP.ge.MPQ
          else if (Lp.lt.Lq) then
            doit_PQ = .true.
          else
            doit_PQ = .false.
          endif

          if( DryRun )then
            MaxMem = max( MaxMem, i_top )
          else
             if(doscreen) then
                nvals=2*MPP*MPQ
                if (.not.MA_Push_get(mt_dbl,nvals,'ffscr',
     C               l_ffscr,i_ffscr)) call errquit(
     E               'hf2: ran out of MA memory',0,MEM_ERR)
                nvals=2*MPP*MPQ*lr3*(lr+1)
                if (.not.MA_Push_get(mt_dbl,nvals,'rscr',
     C               l_rscr,i_rscr)) call errquit(
     E               'hf2: ran out of MA memory',1,MEM_ERR)
                else
                   nvals=1
                   i_rscr=0
                   i_ffscr=0
                endif
            if (doit_PQ) then
              call hf2mkr(Axyz,Bxyz,Cxyz,Dxyz,
     &            W0(i_ALPHAp+2*ipp),W0(i_ALPHAq+2*ipq),
     &            W0(i_R0),W0(i_IJK),W0(i_P),W0(i_Q),W0(i_PQ),
     &            W0(i_ff),W0(i_Rj),
     &            dbl_mb(i_ffscr),dbl_mb(i_rscr),
     &            MPP,MPQ,Lr,Lr3)
            else
              call hf2mkr(Cxyz,Dxyz,Axyz,Bxyz,
     &            W0(i_ALPHAq+2*ipq),W0(i_ALPHAp+2*ipp),
     &            W0(i_R0),W0(i_IJK),W0(i_Q),W0(i_P),W0(i_PQ),
     &            W0(i_ff),W0(i_Rj),
     &            dbl_mb(i_ffscr),dbl_mb(i_rscr),
     &            MPQ,MPP,Lr,Lr3)
            end if
            if(doscreen) then
               if (.not.ma_chop_stack(l_ffscr)) call errquit(
     E              'hf2: machopstack failed',0,MEM_ERR)
            endif
          end if

c Compute the ERI.

c Assign pointers to scratch space.

          i_eri = i_P
          i_t2 = i_eri + ninti
          if (gencon) then
            if (doit_PQ) then
              i_E3 = i_t2 + MPQ*Lq3*max(MPP,NCP)
              i_t1 = i_E3
              i_t3 = i_E3 + MPQ*Lq3
              i_t4 = i_t3 + MPQ
              i_top = max(i_t1+MPQ*Lq3*NPB*NCA,i_E3+MPP,i_t4+NPD)-1
            else
              i_E3 = i_t2 + MPP*Lp3*max(MPQ,NCQ)
              i_t1 = i_E3
              i_t3 = i_E3 + MPP*Lp3
              i_t4 = i_t3 + MPP
              i_top = max(i_t1+MPP*Lp3*NPD*NCC,i_E3+MPQ,i_t4+NPB)-1
            end if
          else
            i_E3 = i_t2 + MPR*Lpq3
            i_top  = i_E3 + max(MPP,MPQ)-1
          end if
          if( DryRun )then
            MaxMem = max( MaxMem, i_top )
          else if (i_top.gt.maxW0) then
            write(LuOut,*) 'HF2:  Insufficient scratch space.'
            write(LuOut,*) '      needed    ',i_top
            write(LuOut,*) '      allocated ',maxW0
            write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
            write(LuOut,*) 'IPAIRp:  ',i_IPAIRp 
            write(LuOut,*) 'ALPHAq:  ',i_ALPHAq 
            write(LuOut,*) 'IPAIRq:  ',i_IPAIRq 
            write(LuOut,*) 'Ep    :  ',i_Ep
            write(LuOut,*) 'Eq    :  ',i_Eq
            write(LuOut,*) 'IJK   :  ',i_IJK
            write(LuOut,*) 'R0    :  ',i_R0
            write(LuOut,*) 'eri   :  ',i_eri
            write(LuOut,*) 't2    :  ',i_t2
            write(LuOut,*) 'E3    :  ',i_E3
            if (gencon) then
              write(LuOut,*) 't1  :  ',i_t1
              write(LuOut,*) 't3  :  ',i_t3
              write(LuOut,*) 't4  :  ',i_t4
            end if
            call errquit('hf2: Not enough scratch space for hferi',911,
     &       MEM_ERR)
          else if (doit_PQ) then
            if (GenCon) then
              call hferi_gen(W0(i_Ep+3*ipp),W0(i_Eq+3*ipq),
     &            W0(i_IPAIRp+istep*ipp),W0(i_IPAIRq+istep*ipq),
     &            W0(i_R0),W0(i_IJK),W0(i_eri),W0(i_E3),
     &            W0(i_t1),W0(i_t2),W0(i_t3),W0(i_t4),
     &            MPP,MPQ,NPP,NPQ,La,Lb,Lc,Ld,La2,Lb2,Lc2,Ld2,
     &            Lq,Lq3,Lr,Acoefs,Bcoefs,Ccoefs,Dcoefs,
     &            NPA,NPB,NPC,NPD,NCA,NCB,NCC,NCD,
     &            MXD,canAB,canCD,canPQ)
            else
              call hferi(W0(i_Ep+3*ipp),W0(i_Eq+3*ipq),
     &            W0(i_R0),W0(i_IJK),W0(i_eri),W0(i_E3),W0(i_t2),
     &            MPP,MPQ,NPP,NPQ,Nints,La,Lb,Lc,Ld,Lr,MXD,
     &            canAB,canCD,canPQ)
            end if
            call daxpy(ninti,1.0D0,W0(i_eri),1,beri,1)
          else
            if (GenCon) then
              call hferi_gen(W0(i_Eq+3*ipq),W0(i_Ep+3*ipp),
     &            W0(i_IPAIRq+istep*ipq),W0(i_IPAIRp+istep*ipp),
     &            W0(i_R0),W0(i_IJK),W0(i_eri),W0(i_E3),
     &            W0(i_t1),W0(i_t2),W0(i_t3),W0(i_t4),
     &            MPQ,MPP,NPQ,NPP,Lc,Ld,La,Lb,Lc2,Ld2,La2,Lb2,
     &            Lp,Lp3,Lr,Ccoefs,Dcoefs,Acoefs,Bcoefs,
     &            NPC,NPD,NPA,NPB,NCC,NCD,NCA,NCB,
     &            MXD,canCD,canAB,canPQ)
            else
              call hferi(W0(i_Eq+3*ipq),W0(i_Ep+3*ipp),
     &            W0(i_R0),W0(i_IJK),W0(i_eri),W0(i_E3),W0(i_t2),
     &            MPQ,MPP,NPQ,NPP,Ninti,Lc,Ld,La,Lb,Lr,MXD,
     &            canCD,canAB,canPQ)
            end if
            j = 1
            do i = 0,lszp-1
              call daxpy(lszq,1.0d0,W0(i_eri+i),lszp,beri(j),1)
              j = j+lszq
            end do
          end if

        end do
      end do

c Return the maximum amount of scratch space required by a "dry run".

   99 if( DryRun ) maxW0 = MaxMem
c
      end

