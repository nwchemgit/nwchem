      Subroutine hf2d(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0,maxW0)
      Implicit None
C $Id$
#include "global.fh"
#include "apiP.fh"
      Logical canAB,canCD,canPQ,dryrun
      
c--> Cartesian Coordinates, Primitives & Contraction Coefficients
      
      Integer La, Lb, Lc, Ld, Nint, MaxW0
      Integer NPA, NCA, NPB, NCB, NPC, NCC, NPD, NCD
      Double Precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Double Precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
      Double Precision Cxyz(3),Cprims(NPC),Ccoefs(NPC,NCC)
      Double Precision Dxyz(3),Dprims(NPD),Dcoefs(NPD,NCD)
      Integer ictra,ictrb,ictrc,ictrd
c--> Block of Electron Repulsion Integrals
      
      Double Precision bERI(Nint,*)
      
c--> Scratch Space
      
      Double Precision W0(maxW0)
c
      integer ithresh
      integer fa,fb,fc,fd,fall
      integer blocksz(4), index
c
      integer l2a,l2b,l2c,l2d,l2q
c
      logical okay
      logical split_all, split_some
      integer countsa, countss, countq, countok
      integer fallmax
      data fallmax /0/
      data countq  /0/
      data countss /0/
      data countsa /0/
      data countok /0/
      save fallmax
      save countq, countss, countsa, countok
c
      ithresh = intd_memthresh
*rak:      write(6,*)' ithresh = ',ithresh
c
c order of preference split2a>split2b>split3a>split2a2b>split3a2b
c
      l2a  = ((La+1)*(La+2)/2)
      l2b  = ((Lb+1)*(Lb+2)/2)
      l2c  = ((Lc+1)*(Lc+2)/2)
      l2d  = ((Ld+1)*(Ld+2)/2)
      l2q  = l2a*l2b*l2c*l2d
      blocksz(1) = npa
      blocksz(2) = npb
      blocksz(3) = npc
      blocksz(4) = npd
      split_all  = .false.
      split_some = .false.
*rak:      write(6,10000)countq,countss,countsa,countok,
*rak:     &    (countss+countsa+countok)
      countq = countq + 1
00001 continue
      fa   = blocksz(1)*l2a
      fb   = blocksz(2)*l2b
      fc   = blocksz(3)*l2c
      fd   = blocksz(4)*l2d
      fall = fa*fb*fc*fd
      fallmax = max(fallmax,fall)
      okay = fall.le.ithresh
      if (.not.okay) then
*rak:        write(6,*)'================================'
*rak:        write(6,*)' fa      :',fa
*rak:        write(6,*)' fb      :',fb
*rak:        write(6,*)' fc      :',fc
*rak:        write(6,*)' fd      :',fd
*rak:        write(6,*)' fall    :',fall
*rak:        write(6,*)' fallmax :',fallmax
*rak:        write(6,*)' blocksz :',blocksz
        call hf2d_getmax_index(blocksz,4,index)
        blocksz(index) = (blocksz(index)+1)/2
*rak:        write(6,*)' blocksz:',blocksz
        split_some = .true.
        if ((blocksz(1).eq.1).and.
     &      (blocksz(2).eq.1).and.
     &      (blocksz(3).eq.1).and.
     &      (blocksz(4).eq.1)) then
          split_all = .true.
          goto 00002
        else
          goto 00001
        endif
      endif
00002 continue
      if (split_all) then
        countsa = countsa + 1
        call hf2d_split_all(
     &      Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &      Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &      Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &      bERI,Nint,canAB,canCD,canPQ,dryrun,
     &      W0,maxW0)
      else if (split_some) then
        countss = countss + 1
        call hf2d_split_block(
     &      Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &      Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &      Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &      bERI,Nint,canAB,canCD,canPQ,dryrun,
     &      W0,maxW0,
     &      blocksz(1),blocksz(2),blocksz(3),blocksz(4))
      else
        countok = countok + 1
        call hf2dold(
     &      Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &      Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &      Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &      bERI,Nint,canAB,canCD,canPQ,dryrun,
     &      W0,maxW0)
      endif
10000 format(1x,'total=',i10,2x,'some=',i10,2x,'all=',i10,
     &    2x,'ok=',i10,2x,'sum=',i10)
      end
      Subroutine hf2dold(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0,maxW0)
      
      Implicit None
#include "mafdecls.fh"
#include "errquit.fh"
#include "case.fh"
      
      Logical canAB,canCD,canPQ,dryrun
      
      Logical GENCON
      
      
c--> Cartesian Coordinates, Primitives & Contraction Coefficients
      
      Integer La, Lb, Lc, Ld, Nint, MaxW0
      Integer NPA, NCA, NPB, NCB, NPC, NCC, NPD, NCD
      Double Precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Double Precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)
      Double Precision Cxyz(3),Cprims(NPC),Ccoefs(NPC,NCC)
      Double Precision Dxyz(3),Dprims(NPD),Dcoefs(NPD,NCD)
      Integer ictra,ictrb,ictrc,ictrd
c--> Block of Electron Repulsion Integrals
      
      Double Precision bERI(Nint,*)
      
c--> Scratch Space
      
      Double Precision W0(maxW0)
      
c--> Derivative Indices
      
      Integer mxder, Nder1, Nder2, Nder3
      Parameter (mxder=3)
      Integer inder1(6,mxder),inder2(3,mxder),inder3(3,mxder)
c:: local
      integer maxmem
      integer nd
      integer MXD, NCP, NCQ
      integer Lp, Lp3, Lq, Lq3, Lr, Lr3
      integer NPP, NPQ, NPR
      integer i_left, i_right, i_top
      integer i_ALPHAp, i_IPAIRp, i_ESp, i_ALPHAq, i_IPAIRq, i_ESq
      integer i_ERI, i_Eq, i_Ep, i_IJK, i_R0
      integer i_pf, i_P, i_Q, i_PQ, i_ff, i_Rj, i_E3, i_sum
      integer i_dP, i_dQ, i_dR, i_dS
      integer i_ffscr,l_ffscr,i_rscr,l_rscr,nvals
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
      
      
c Assign pointers to scratch space.
      
      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
      
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp
      
      if( i_left.ge.i_right .and. .not.dryrun )then
        
        write(*,*) 'HF2D:  Insufficient scratch space.'
        write(*,*) '      needed    ',i_left + (maxW0-(i_right-1))
        write(*,*) '      allocated ',maxW0
        
        write(*,*) 'From the left '
        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp 
        write(*,*) 'From the right '
        write(*,*) 'ESp   :  ',i_ESp 
        
        stop
        
      end if
      
      MaxMem = 1   ! take care of compiler warnings
      if (dryrun) then
        MaxMem = i_left + (maxW0-(i_right-1))
        NPP = NPA*NPB
      else
        call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,
     &      GENCON,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)
      endif
      
c Define the prefactor of the overlap distribution "Q".
      
c Assign pointers to scratch space.
      
      i_ALPHAq = i_IPAIRp + 2*(NPA*NPB)
      i_IPAIRq = i_ALPHAq + 2*(NPC*NPD)
      i_left   = i_IPAIRq + 2*(NPC*NPD) - 1
      
      i_ESq   = i_right - 3*(NPC*NPD)
      i_right = i_ESq
      
      if( i_left.ge.i_right .and. .not.dryrun )then
        
        write(*,*) 'HF2D:  Insufficient scratch space.'
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
      
      if( dryrun )then
        MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )
        NPQ = NPC*NPD
      else
        call hfset(Cxyz,Cprims,Ccoefs,NPC,NCC,
     &      Dxyz,Dprims,Dcoefs,NPD,NCD,
     &      GENCON,W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ESq),NPQ)
      endif
      
c Define the Hermite linear expansion coefficients.
      
c Assign pointers to scratch space.
      
      NPR = NPP*NPQ
      
      if (NPR.eq.0) then
c... no integrals because one charge distribution P or Q has no overlap
        call dfill((nint*12),0.0d00,bERI,1)
        return
      endif
      i_ERI  = i_IPAIRq + 2*(NPC*NPD)
      i_Ep   = i_ERI    + NPR*(Nint*6)
      i_Eq   = i_Ep     + 3*NPP*(MXD+1)*((Lp+1)*(La+1)*(Lb+1))
      i_pf   = i_Eq     + 3*NPQ*(MXD+1)*((Lq+1)*(Lc+1)*(Ld+1))
      i_left = i_pf     + 2*max(NPP,NPQ) - 1
      
      if( i_left.ge.i_right .and. .not.dryrun )then
        
        write(*,*) 'HF2D:  Insufficient scratch space.'
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
      
      if( dryrun )then
        MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )
      else
        do 100 nd = 0,MXD
          call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),
     &        W0(i_pf),nd,NPP,MXD,La,Lb)
          call hfmke(Cxyz,Dxyz,W0(i_ALPHAq),W0(i_ESq),W0(i_Eq),
     &        W0(i_pf),nd,NPQ,MXD,Lc,Ld)
100     continue
      endif
      
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
      
      if( i_top.gt.maxW0 .and. .not.dryrun )then
        
        write(*,*) 'HF2D:  Insufficient scratch space.'
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
      
      if( dryrun )then
        MaxMem = max( MaxMem, i_top)
      else
         nvals=1
         if(doscreen) nvals=2*NPP*NPQ
         if (.not.MA_Push_get(mt_dbl,nvals,'ffscr',
     C        l_ffscr,i_ffscr)) call errquit(
     E        'hf2d: ran out of MA memory',0,MEM_ERR)
         nvals=1
         if(doscreen) nvals=2*NPP*NPQ*lr3*(lr+1)
         if (.not.MA_Push_get(mt_dbl,nvals,'rscr',
     C        l_rscr,i_rscr)) call errquit(
     E        'hf2d: ran out of MA memory',1,MEM_ERR)
         
        call hf2mkr(
     &      Axyz,Bxyz,Cxyz,Dxyz,W0(i_ALPHAp),W0(i_ALPHAq),
     &      W0(i_R0),W0(i_IJK),W0(i_P),W0(i_Q),W0(i_PQ),
     &      W0(i_ff),W0(i_Rj),
     &        dbl_mb(i_ffscr),dbl_mb(i_rscr),
     N        NPP,NPQ,Lr,Lr3)
        if (.not.ma_chop_stack(l_ffscr)) call errquit(
     E       'hf2d: machopstack failed',0,MEM_ERR)
      endif
      
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
      
c Assign pointers to scratch space.
      
      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6
      Lq3 = ((Lq+1)*(Lq+2)*(Lq+3))/6
      
      i_E3   = i_IJK + (Lr+1)**3
      i_sum  = i_E3  + max(NPP,NPQ)
      i_top  = i_sum + NPR*6*Lq3      ! use i_sum for i_ff (size=NPR*2) to save memory 
      i_top  = i_top - 1              ! pointer offset of 1 removed
      
      if( i_top.gt.maxW0 .and. .not.dryrun )then
        
        write(*,*) 'HF2D:  Insufficient scratch space.'
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
      
      if( dryrun )then
        MaxMem = max( MaxMem, i_top)
      else
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
        
        call hfderi(
     &      W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dP),
     &      W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &      MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &      canAB,canCD,canPQ)
        
        i_dR  = i_dP + NPR*(Nint*3) ! allocated with i_eri
        
c Integral derivatives wrt R.
        
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
        
        if (ictra.ne.ictrb) then
          call hfderi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dR),
     &        W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &        MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &        canAB,canCD,canPQ)
        endif
        
c Contract primitive integrals. for a and b centers
        
c Transform to obtain derivatives wrt (A,B).
        
        i_ff  = i_sum        ! iff size=2*NPR (use i_sum instead)
        call hf2PRtoAB(
     &      W0(i_dP),W0(i_dR),bERI(1,1),bERI(1,4),
     &      W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ff),NPP,NPQ,Nint*3,
     &      ictra,ictrb)
c
c Integral derivatives wrt S.
        
        i_dS = i_dR 
        
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
        
        if (ictrc.ne.ictrd) then
          call hfderi(W0(i_Ep),W0(i_Eq),W0(i_R0),W0(i_IJK),W0(i_dS),
     &        W0(i_E3),W0(i_sum),NPP,NPQ,Nint,La,Lb,Lc,Ld,Lr,
     &        MXD,inder1,Nder1,inder2,Nder2,inder3,Nder3,
     &        canAB,canCD,canPQ)
        endif
        
c Contract primitive integrals. for centers c and d 

c Translation invariance requires d/dQ = -d/dP !
        
        i_dQ  = i_dP
        call dscal((NPR*(Nint*3)),-1.0d00,W0(i_dQ),1)
        call hf2QStoCD(
     &      W0(i_dQ),W0(i_dS),bERI(1,7),bERI(1,10),
     &      W0(i_ALPHAq),W0(i_IPAIRq),W0(i_ff),NPP,NPQ,Nint*3,
     &      ictrc,ictrd)
c
      endif
c
      if (dryrun) maxW0 = MaxMem
c
00099 continue
      end
      subroutine hf2d_getmax_index(b,nb,i)
      implicit none
      integer nb
      integer b(nb)
      integer i
c
      integer j, maxval
c
      i = 1
      maxval = b(i)
      do j = 2,nb
        if (b(j) .gt. maxval) then
          maxval = b(j)
          i = j
        endif
      enddo
      end
