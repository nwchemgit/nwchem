      Subroutine hf1(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,zan,exinv,ncenters,
     &               bO2I,bKEI,bNAI,Nint,O2I,KEI,NAI,canAB,
     &               DryRun,W0,maxW0)
c $Id: hf1.F,v 1.15 1999-10-26 21:52:25 mg201 Exp $
      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical O2I,KEI,NAI,canAB

      Logical GenCon,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates, Charges and Inverse Exponents

      Dimension Cxyz(3,ncenters),zan(ncenters),exinv(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      Dimension bO2I(Nint),bKEI(Nint),bNAI(Nint)

c--> Scratch Space.

      Dimension W0(maxW0)
#include "mafdecls.fh"
c
c Compute the overlap, kinetic energy, and nuclear attraction integrals for 
c two shells of contracted Gaussians functions. This driver is NOT capable of 
c evaluating integral derivatives.
c
c******************************************************************************
#if defined(INTDEBUG)
      if (.not.dryrun) then
        write(6,*)' inside hf1 '
        write(6,*)' npa,nca,la = ',npa,nca,la
        write(6,*)' npb,ncb,lb = ',npb,ncb,lb
        write(6,*)' ncenters   = ',ncenters
        write(6,*)' NINT       = ',nint
        write(6,*)' maxW0      = ',maxw0
        write(6,*)' <canAB:DryRun>-<',canab,':',dryrun,'>'
        write(6,*)' <o2i:kei:nai>-<',o2i,':',kei,':',nai,'>'
        write(6,'(a,3(2x,1pd20.10))')' Axyz =',Axyz
        write(6,'(a,3(2x,1pd20.10))')' Bxyz =',Bxyz
        write(6,'(a,100(3(2x,1pd20.10/)))')' Cxyz =',Cxyz
        do jjjj = 1,nca
        do iiii = 1,npa
          write(6,'(a,i3,a,2(2x,1pd20.10))')
     &        'Aprims:Acoeffs:(',iiii,') =',Aprims(iiii),
     &        Acoefs(iiii,jjjj)
        enddo
        enddo
        do jjjj = 1,ncb
        do iiii = 1,npb
          write(6,'(a,i3,a,2(2x,1pd20.10))')
     &        'Bprims:Bcoeffs:(',iiii,') =',Bprims(iiii),
     &        Bcoefs(iiii,jjjj)
        enddo
        enddo
      endif
      if (.not.dryrun) then
        write(6,*)' Li/Lj',La,'/',Lb
        write(6,*)' i_ngen ',nca
        write(6,*)' j_ngen ',ncb
        write(6,*)' int_hf1: lstv',Nint
        write(6,*)' int_hf1: lscr',maxW0
        call dfill(maxW0,0.0d00,W0,1)
        if (O2I) call dfill(Nint,0.0d00,bO2I,1)
        if (KEI) call dfill(Nint,0.0d00,bKEI,1)
        if (NAI) call dfill(Nint,0.0d00,bNAI,1)
*debug_ma:        call MA_summarize_allocated_blocks()
*debug_ma:        write(6,*)' int_hf1: ma verify 1-b4'
*debug_ma:        status = ma_verify_allocator_stuff()
*debug_ma:        write(6,*)' verstat = ',status
*debug_ma:        write(6,*)' int_hf1: ma verify 1-af'
        call util_flush(6)
      endif
      if (.not.dryrun) then
        call hf_print_set(1)
        call hf_print('hf1: a shell',axyz,aprims,acoefs,npa,nca,la)
        call hf_print('hf1: b shell',bxyz,bprims,bcoefs,npb,ncb,lb)
        call hf_print_set(0)
      endif
#endif
      MXD = 0

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB

      GenCon = NCP.ne.1

*      write(6,*)' gencon ',gencon,nca,ncb


c To determine all the Hermite expansion coefficients required to evaluate
c the kinetic energy integrals, increment the angular momenta by one.

      if( KEI )then
       Li = 1
      else
       Li = 0
      end if

c Define the angular momentum of the overlap distribution.

      Lp = La + Lb

c Increment "Lp" to account for the order of differentiation.

      Lp = Lp + MXD

c Define the accumulated number of angular momentum functions <= Lp.

      Lp3 = ((Lp+1)*(Lp+2)*(Lp+3))/6

c Define the prefactor of the overlap distribution "P".

c Assign pointers to scratch space.
 
      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
 
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp
 
      if( i_left.ge.i_right )then
 
       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0
       write(*,*) ' DryRun ',DryRun
 
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
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,
     &      GenCon,W0(i_ALPHAp),W0(i_IPAIRp),
     &      W0(i_ESp),NPP)

      end if

      if (npp.eq.0) then
          nintlocal = (La+1)*(La+2)/2*nca
          nintlocal = nintlocal*(Lb+1)*(Lb+2)/2*ncb
          if (O2I) call dfill(nintlocal,0.0d00,bO2I,1)
          if (KEI) call dfill(nintlocal,0.0d00,bKEI,1)
          if (NAI) call dfill(nintlocal,0.0d00,bNAI,1)
        return
      endif
c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_pf   = i_Ep     + 3*NPP*(MXD+1)*lprod
      iptr   = i_pf     + 2*NPP
      if (gencon) then
        La2 = (La+1)*(La+2)/2
        Lb2 = (Lb+1)*(Lb+2)/2
        i_gct_a    = iptr 
        i_gct_b    = i_gct_a   + NPP*NCA
        i_prim_ints = i_gct_b  + NPP*NCB
        i_half_ints = i_prim_ints + npp*La2*Lb2
        iptr        = i_half_ints +
     &        max((npp*La2*Lb2*nca),(ncb*La2*Lb2*nca))
      else
        i_gct_a     = -1
        i_gct_b     = -1
        i_prim_ints = -1
      endif
      i_left = iptr - 1
      if( i_left.ge.i_right )then

       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0
       write(*,*) ' DryRun ',DryRun

       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'pf    :  ',i_pf
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp

       stop

      end if
c form transformation matrix for 1 electron
      if (gencon.and.(.not.dryrun)) then
        call hf1_tran_gen(W0(i_gct_a),W0(i_gct_b),W0(i_IPAIRp),NPP,
     &      Acoefs,NPA,NCA,Bcoefs,NPB,NCB)
      endif                                        

      if( DryRun )then

       MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )

      else

       do 100 nd = 0,MXD
        call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &             nd,NPP,MXD,La+Li,Lb+Li)
  100  continue

      end if
       
c Compute the 2-center overlap integrals, <a|S|b>.

      if( O2I )then
        if( .not. DryRun )then
          if (gencon) then
            call hf2oi_gc(W0(i_Ep),bO2I,
     &          W0(i_prim_ints),W0(i_half_ints),
     &          (la2*lb2),
     &          NCA,NCB,NPP,La,Lb,Li,
     &          W0(i_gct_a),W0(i_gct_b),
     &          canAB)
          else
            call hf2oi(W0(i_Ep),bO2I,Nint,NPP,La,Lb,Li,canAB)
          endif
        end if
      end if

c Compute kinetic energy integrals, <a|T|b>.

      if( KEI )then

c Assign pointers to scratch space.

       if (gencon) then
         i_Ti  = iptr
         i_top = i_Ti + NPP - 1 
       else
         i_Ti  = i_Ep + 3*NPP*(MXD+1)*lprod
         i_top = i_Ti + NPP - 1
       endif
       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1:  Insufficient scratch space.'
        write(*,*) '       needed    ',i_top
        write(*,*) '       allocated ',maxW0
        write(*,*) ' DryRun ',DryRun

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'Ti    :  ',i_Ti

        stop

       end if

       if( DryRun )then

        MaxMem = max( MaxMem, i_top )

       else

         if (gencon) then
           call hfkei_gc(W0(i_ALPHAp),W0(i_Ep),bKEI,
     &         W0(i_prim_ints),W0(i_half_ints),
     &         W0(i_Ti),(la2*lb2),
     &         NCA,NCB,NPP,
     &         La,Lb,Li,
     &         W0(i_gct_a),W0(i_gct_b),
     &         canAB)

         else
           call hfkei(W0(i_ALPHAp),W0(i_Ep),bKEI,W0(i_Ti),
     &         Nint,NPP,La,Lb,Li,canAB)
         endif

       end if

      end if
       
c Compute nuclear attraction integrals, <a|V|b>.

      if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

       if (gencon) then
         i_R0  = iptr
       else
         i_R0  = i_Ep  + 3*NPP*(MXD+1)*lprod
       endif
       i_IJK = i_R0  + NPP*Lp3
       i_P   = i_IJK + (Lp+1)**3
       i_RS  = i_P   + NPP*3
       i_PC  = i_RS  + NPP
       i_ff  = i_PC  + NPP*3
       i_Rj  = i_ff  + NPP*2
       i_top = i_Rj  + NPP*(Lp+1)*Lp3 - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1:  Insufficient scratch space.'
        write(*,*) '       needed    ',i_top
        write(*,*) '       allocated ',maxW0
        write(*,*) ' DryRun ',DryRun

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'R0    :  ',i_R0
        write(*,*) 'IJK   :  ',i_IJK
        write(*,*) 'P     :  ',i_P
        write(*,*) 'RS    :  ',i_RS
        write(*,*) 'PC    :  ',i_PC
        write(*,*) 'ff    :  ',i_ff
        write(*,*) 'Rj    :  ',i_Rj

        stop

       end if

       if( DryRun )then

        MaxMem = max( MaxMem, i_top )

       else

        call hf1mkr(Axyz,Bxyz,Cxyz,zan,exinv,ncenters,
     &              W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &              W0(i_Rj),W0(i_R0),W0(i_R0),W0(i_IJK),
     &              NPP,Lp,Lp3,.FALSE.)

        if (gencon) then
          call hfnai_gc(W0(i_Ep),W0(i_R0),W0(i_IJK),bNAI,
     &        W0(i_prim_ints),W0(i_half_ints),
     &        (La2*Lb2),NCA,NCB,NPP,
     &        La,Lb,Li,Lp,Lp3,
     &        W0(i_gct_a),W0(i_gct_b),
     &        canAB)
        else
          call hfnai(W0(i_Ep),W0(i_R0),W0(i_IJK),bNAI,
     &        Nint,NPP,La,Lb,Li,Lp,Lp3,canAB)
        endif

       end if

      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
