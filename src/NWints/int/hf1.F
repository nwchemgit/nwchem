      Subroutine hf1(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,zan,exinv,ncenters,
     &               bO2I,bKEI,bNAI,Nints,O2I,KEI,NAI,canAB,
     &               DryRun,W0,maxW0)
c $Id$
      Implicit none
      integer NPA,NCA,La,NPB,NCB,Lb,ncenters,Nints,maxW0
      logical O2I,KEI,NAI,canAB
      logical GenCon,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      double precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      double precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates, Charges and Inverse Exponents

      double precision Cxyz(3,ncenters),zan(ncenters),exinv(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      double precision bO2I(Nints),bKEI(Nints),bNAI(Nints)

c--> Scratch Space.

      double precision W0(maxW0)

c--> Local variables

      integer MXD,NCP,NPP,La2,Lb2,Li,Lp,Lp3,MaxMem,nintlocal,lprod,nd
      integer i_ALPHAp,i_IPAIRp,i_ESp,i_left,i_right,i_top,i_Ep,i_ptr,
     &    i_pf,i_prim_ints,i_half_ints,i_Ti,i_R0,i_IJK,i_P,i_RS,i_PC,
     &    i_ff,i_Rj
#if defined(INTDEBUG)
      integer jjjj,iiii
#endif
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
c
c Compute the overlap, kinetic energy, and nuclear attraction integrals for 
c two shells of contracted Gaussians functions. This driver is NOT capable of 
c evaluating integral derivatives.
c
c******************************************************************************
#if defined(INTDEBUG)
      if (.not.dryrun) then
        write(LuOut,*)' inside hf1 '
        write(LuOut,*)' npa,nca,la = ',npa,nca,la
        write(LuOut,*)' npb,ncb,lb = ',npb,ncb,lb
        write(LuOut,*)' ncenters   = ',ncenters
        write(LuOut,*)' NINTS       = ',nints
        write(LuOut,*)' maxW0      = ',maxw0
        write(LuOut,*)' <canAB:DryRun>-<',canab,':',dryrun,'>'
        write(LuOut,*)' <o2i:kei:nai>-<',o2i,':',kei,':',nai,'>'
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
        write(LuOut,*)' Li/Lj',La,'/',Lb
        write(LuOut,*)' i_ngen ',nca
        write(LuOut,*)' j_ngen ',ncb
        write(LuOut,*)' int_hf1: lstv',Nints
        write(LuOut,*)' int_hf1: lscr',maxW0
        call dcopy(maxW0,0d0,0,W0,1)
        if (O2I) call dcopy(Nints,0d0,0,bO2I,1)
        if (KEI) call dcopy(Nints,0d0,0,bKEI,1)
        if (NAI) call dcopy(Nints,0d0,0,bNAI,1)
*debug_ma:        call MA_summarize_allocated_blocks()
*debug_ma:        write(LuOut,*)' int_hf1: ma verify 1-b4'
*debug_ma:        status = ma_verify_allocator_stuff()
*debug_ma:        write(LuOut,*)' verstat = ',status
*debug_ma:        write(LuOut,*)' int_hf1: ma verify 1-af'
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

      MaxMem = 1    ! take care of compiler warnings
      if( DryRun )then

        MaxMem = i_left + (maxW0 - (i_right-1))
        NPP = NPA*NPB

      else if (i_left.ge.i_right) then
 
        write(LuOut,*) 'HF1:  Insufficient scratch space.'
        write(LuOut,*) '       needed    ',i_left+(maxW0-(i_right-1))
        write(LuOut,*) '       allocated ',maxW0
        write(LuOut,*) ' DryRun ',DryRun
 
        write(LuOut,*) 'From the left '
        write(LuOut,*) 'ALPHAp:  ',i_ALPHAp
        write(LuOut,*) 'IPAIRp:  ',i_IPAIRp
        write(LuOut,*) 'From the right '
        write(LuOut,*) 'ESp   :  ',i_ESp
 
        call errquit('hf1: insufficient memory for hfset',911, MEM_ERR)
 
      else
 
        call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &      Bxyz,Bprims,Bcoefs,NPB,NCB,
     &      GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)

      end if

      La2 = (La+1)*(La+2)/2
      Lb2 = (Lb+1)*(Lb+2)/2

c Zero out the integrals. Return if screening test gives no pairs

      if (.not.DryRun) then
        nintlocal = La2*Lb2*nca*ncb
        if (O2I) call dcopy(nintlocal,0d0,0,bO2I,1)
        if (KEI) call dcopy(nintlocal,0d0,0,bKEI,1)
        if (NAI) call dcopy(nintlocal,0d0,0,bNAI,1)
        if (NPP.eq.0) return
      end if

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

C      write (LuOut,*) 'before hfmke'

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_ptr  = i_Ep     + 3*NPP*(MXD+1)*lprod

      i_prim_ints = i_ptr            ! take care of compiler warnings
      i_half_ints = i_prim_ints + NPP

      if (gencon) then
        i_prim_ints = i_ptr
        i_half_ints = i_prim_ints + NPP
        i_ptr       = i_half_ints + NPA*NCB
      end if

      i_pf   = i_ptr
      i_left = i_pf     + 2*NPP - 1

      if( DryRun )then

        MaxMem = max( MaxMem, i_left + (maxW0 - (i_right-1)) )

      else if( i_left.ge.i_right) then

        write(LuOut,*) 'HF1:  Insufficient scratch space.'
        write(LuOut,*) '       needed    ',i_left+(maxW0-(i_right-1))
        write(LuOut,*) '       allocated ',maxW0
        write(LuOut,*) ' DryRun ',DryRun

        write(LuOut,*) 'From the left '
        write(LuOut,*) 'ALPHAp:  ',i_ALPHAp
        write(LuOut,*) 'IPAIRp:  ',i_IPAIRp
        write(LuOut,*) 'Ep    :  ',i_Ep
        write(LuOut,*) 'pf    :  ',i_pf
        write(LuOut,*) 'From the right '
        write(LuOut,*) 'ESp   :  ',i_ESp

        call errquit('hf1: insufficient memory for hfmke',911, MEM_ERR)

      else

        do nd = 0,MXD
          call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),W0(i_pf),
     &        nd,NPP,MXD,La+Li,Lb+Li)
        end do

      end if
       
c Compute the 2-center overlap integrals, <a|S|b>.

      if( O2I )then
        if( .not. DryRun )then
          if (gencon) then
            call hf2oi_gc(W0(i_Ep),bO2I,W0(i_prim_ints),W0(i_half_ints),
     &          Acoefs,Bcoefs,W0(i_IPAIRp),NPA,NPB,NCA,NCB,NPP,
     &          La,Lb,La2,Lb2,Li,canAB)
          else
            call hf2oi(W0(i_Ep),bO2I,Nints,NPP,La,Lb,Li,canAB)
          endif
        end if
      end if

c Compute kinetic energy integrals, <a|T|b>.

      if (KEI) then

c Assign pointers to scratch space.

        i_Ti  = i_ptr
        i_top = i_Ti + NPP - 1

        if( DryRun )then

          MaxMem = max( MaxMem, i_top )

        else if( i_top.gt.maxW0 )then

          write(LuOut,*) 'HF1:  Insufficient scratch space.'
          write(LuOut,*) '       needed    ',i_top
          write(LuOut,*) '       allocated ',maxW0
          write(LuOut,*) ' DryRun ',DryRun

          write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
          write(LuOut,*) 'IPAIRp:  ',i_IPAIRp
          write(LuOut,*) 'Ep    :  ',i_Ep
          write(LuOut,*) 'Ti    :  ',i_Ti

          call errquit('hf1: insufficient memory for hfkei',911,
     &       MEM_ERR)

        else if (gencon) then

          call hfkei_gc(W0(i_ALPHAp),W0(i_Ep),bKEI,
     &        W0(i_prim_ints),W0(i_half_ints),W0(i_Ti),
     &        Acoefs,Bcoefs,W0(i_IPAIRp),
     &        NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,canAB)
        else

          call hfkei(W0(i_ALPHAp),W0(i_Ep),bKEI,W0(i_Ti),
     &        Nints,NPP,La,Lb,Li,canAB)

        end if

      end if
       
c Compute nuclear attraction integrals, <a|V|b>.

      if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

        i_R0  = i_ptr
        i_IJK = i_R0  + NPP*Lp3
        i_P   = i_IJK + (Lp+1)**3
        i_RS  = i_P   + NPP*3
        i_PC  = i_RS  + NPP
        i_ff  = i_PC  + NPP*3
        i_Rj  = i_ff  + NPP*2
        i_top = i_Rj  + NPP*(Lp+1)*Lp3 - 1

        if( DryRun )then

          MaxMem = max( MaxMem, i_top )
C          write (LuOut,*) MaxMem,maxW0

        else if( i_top.gt.maxW0 .and. .not.DryRun )then

          write(LuOut,*) 'HF1:  Insufficient scratch space.'
          write(LuOut,*) '       needed    ',i_top
          write(LuOut,*) '       allocated ',maxW0
          write(LuOut,*) ' DryRun ',DryRun

          write(LuOut,*) 'ALPHAp:  ',i_ALPHAp 
          write(LuOut,*) 'IPAIRp:  ',i_IPAIRp
          write(LuOut,*) 'Ep    :  ',i_Ep
          write(LuOut,*) 'R0    :  ',i_R0
          write(LuOut,*) 'IJK   :  ',i_IJK
          write(LuOut,*) 'P     :  ',i_P
          write(LuOut,*) 'RS    :  ',i_RS
          write(LuOut,*) 'PC    :  ',i_PC
          write(LuOut,*) 'ff    :  ',i_ff
          write(LuOut,*) 'Rj    :  ',i_Rj

          call errquit('hf1: insufficient memory for hfmkr',911,
     &       MEM_ERR)

        else

          call hf1mkr(Axyz,Bxyz,Cxyz,zan,exinv,ncenters,
     &        W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &        W0(i_Rj),W0(i_R0),W0(i_R0),W0(i_IJK),
     &        NPP,Lp,Lp3,.FALSE.)

          if (gencon) then
            call hfnai_gc(W0(i_Ep),W0(i_R0),W0(i_IJK),bNAI,
     &          W0(i_prim_ints),W0(i_half_ints),
     &          Acoefs,Bcoefs,W0(i_IPAIRp),
     &          NPA,NPB,NCA,NCB,NPP,La,Lb,La2,Lb2,Li,Lp,Lp3,canAB)
          else
            call hfnai(W0(i_Ep),W0(i_R0),W0(i_IJK),bNAI,
     &          Nints,NPP,La,Lb,Li,Lp,Lp3,canAB)
          end if

        end if

      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
