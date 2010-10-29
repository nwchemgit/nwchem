      Subroutine hf1_er(Axyz,Aprims,Acoefs,NPA,NCA,La,
     &               Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,
     &               Cxyz,zan,ncenters,
     &               bO2I,bKEI,bNAI,Nint,O2I,KEI,NAI,canAB,
     &               DryRun,W0,maxW0)
c $Id$
      implicit none

      integer NPA,NCA,NPB,NCB
      integer La,Lb,NInt
      integer ncenters,maxW0
      Logical O2I,KEI,NAI,canAB

      Logical GenCon,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      double precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      double precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates & Charges

      double precision Cxyz(3,ncenters),zan(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      double precision bO2I(Nint),bKEI(Nint),bNAI(ncenters,Nint)

c--> Scratch Space.

      double precision W0(maxW0)
C
C     local
C
      integer MXD,NCP,Li,Lp,Lp3,NPP
      integer i_ALPHAp,i_IPAIRp,i_left,i_ESp,i_right, i_exinv
      integer i_rs,i_p,i_ijk,i_ff,i_pc,i_r0c,i_pf,i_ep
      integer lprod,nd,i_r0,i_top,i_ti,i_rj
      integer MaxMem
c
c Compute the overlap, kinetic energy, and nuclear attraction integrals for 
c two shells of contracted Gaussians functions. This driver is NOT capable of 
c evaluating integral derivatives.
c
c******************************************************************************
      MXD = 0

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB

      GenCon = NCP.ne.1

      if( GenCon )then
       write(*,*) 'HF1: Not yet ready for general contraction.'
       stop
      end if

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
 
      i_exinv = 1
      i_ALPHAp = i_exinv+ncenters
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
 
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp
 
      if( i_left.ge.i_right )then
 
       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0
 
       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp
 
       stop
 
      end if
 
      MaxMem = 1       ! take care of compiler warnings

      if( DryRun )then

       MaxMem = i_left + (maxW0 - (i_right-1))
       NPP = NPA*NPB

      else

       call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &            Bxyz,Bprims,Bcoefs,NPB,NCB,
     &            GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)

      end if

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_pf   = i_Ep     + 3*NPP*(MXD+1)*lprod
      i_left = i_pf     + 2*NPP - 1

      if( i_left.ge.i_right )then

       write(*,*) 'HF1:  Insufficient scratch space.'
       write(*,*) '       needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '       allocated ',maxW0

       write(*,*) 'From the right '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'pf    :  ',i_pf
       write(*,*) 'From the left '
       write(*,*) 'ESp   :  ',i_ESp

       stop

      end if

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
        call hf2oi(W0(i_Ep),bO2I,Nint,NPP,La,Lb,Li,canAB)
       end if
      end if

c Compute kinetic energy integrals, <a|T|b>.

      if( KEI )then

c Assign pointers to scratch space.

       i_Ti  = i_Ep + 3*NPP*(MXD+1)*lprod
       i_top = i_Ti + NPP - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1:  Insufficient scratch space.'
        write(*,*) '       needed    ',i_top
        write(*,*) '       allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'Ti    :  ',i_Ti

        stop

       end if

       if( DryRun )then

        MaxMem = max( MaxMem, i_top )

       else

        call hfkei(W0(i_ALPHAp),W0(i_Ep),bKEI,W0(i_Ti),
     &             Nint,NPP,La,Lb,Li,canAB)
       end if

      end if
       
c Compute nuclear attraction integrals, <a|V|b>.

      if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

       i_R0  = i_Ep  + 3*NPP*(MXD+1)*lprod
       i_R0C = i_R0  + NPP*Lp3
       i_IJK = i_R0C  + NPP*Lp3*ncenters
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

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'R0    :  ',i_R0
        write(*,*) 'R0C   :  ',i_R0C
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

        call dfill (ncenters,0.0d0,W0(i_exinv),1)
        call hf1mkr(Axyz,Bxyz,Cxyz,zan,W0(i_exinv),ncenters,
     &              W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &              W0(i_Rj),W0(i_R0),W0(i_R0C),W0(i_IJK),
     &              NPP,Lp,Lp3,.TRUE.)

        call hfnai_er(ncenters,W0(i_Ep),W0(i_R0C),W0(i_IJK),bNAI,
     &             Nint,NPP,La,Lb,Li,Lp,Lp3,canAB)

       end if

      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
