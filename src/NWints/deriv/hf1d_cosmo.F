      Subroutine hf1d_cosmo(
     &       Axyz,Aprims,Acoefs,NPA,NCA,La,ictrA,
     &       Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrB,
     &       Cxyz,zan,
     &       bNAI,Nint,
     &       canAB,dryrun,W0,maxW0)
c $Id$
c
c     Calculate derivatives of an nuclear attraction type integral.
c     The integral is:
c
c       E = int(dr;Chi_A(r)Chi_B(r)Q_C/|r-R_C|)
c
c     We are interested in the derivatives with respect to the centers
c     A, B, and C. So upon return we want
c
c       bNAI(1:nint,1:3) = dE/dA                  (1)
c
c       bNAI(1:nint,4:6) = dE/dB                  (2)
c
c       bNAI(1:nint,7:9) = dE/dC = -dE/dA-dE/dB   (3)
c
c    Typically the calculation of derivative integrals relies on the
c    fact that the derivative of a Gaussian function with respect to
c    its center position is again a Gaussian function but of a higher
c    angular momentum. Hence the same code can be used to calculate
c    derivative integrals as well as energy integrals (I am not sure
c    whether the same would be true if directly calculating the 
c    derivative with respect to C as that would change the operator).
c
c    So we expect this code to evaluate eqs. (1) and (2) and use
c    translational invariance to evaluate eq. (3) (i.e. the second
c    equality).
c
c     Implicit real*8 (a-h,o-z)
c     Implicit integer (i-n)
      Implicit none

#include "errquit.fh"

      Logical dryrun
      Logical canAB,GenCon ! what is canAB???
      Integer Nint, maxW0

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Integer NPA, NCA, NPB, NCB, NCP
      Double precision Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Double precision Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates, Charges & Inverse Exponents

      Double precision Cxyz(3), zan

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      Double precision bNAI(Nint,9)

c--> Derivative Indices

      Integer inder1(3,3),inder2(6,6)

c--> Scratch Space.

      Double precision W0(maxW0)

c--> Other stuff

      Integer lprod, mxd, Lp3, Lp, Nder, NPP
      Integer La, Lb, Li, ictrA, ictrB
      Integer i_VA, i_VB, i_VP, i_VR, i_Ri, i_Rj, i_RS
      Integer MaxMem, nd
      Integer i_P, i_R0, i_top, i_pf, i_PC, i_R0C, i_left, i_right
      Integer i_IPAIRp, i_ESp, i_ALPHAp, i_Ep
      Integer i_IJK, i_ff

      Data inder1/ 1,0,0,
     &             0,1,0,
     &             0,0,1 /

      Data inder2/ 1,0,0,0,0,0,
     &             0,1,0,0,0,0,
     &             0,0,1,0,0,0,
     &             0,0,0,1,0,0,
     &             0,0,0,0,1,0,
     &             0,0,0,0,0,1 /
c
c Compute gradient of the overlap, kinetic energy, and nuclear attraction 
c integrals for 2 shells of contracted Gaussians functions.
c
c The derivative indices specify the order of differentiation for each 
c coordinate:
c
c      inder1(i,n), i=1-3:   d/dRx, d/dRy, d/dRz
c
c      inder2(i,n), i=1-3:   d/dPx, d/dPy, d/dPz
c                   i=4-6:   d/dRx, d/dRy, d/dRz
c
c******************************************************************************

#if defined(INTDEBUG)
      call hf_print_set(1)
      call hf_print('hf1d: a shell',axyz,aprims,acoefs,npa,nca,la)
      call hf_print('hf1d: b shell',bxyz,bprims,bcoefs,npb,ncb,lb)
      call hf_print_set(0)
#endif
      MXD = 1
      Li  = 0

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB

      GenCon = NCP.ne.1

      if( GenCon )then
        write(*,*) 'HF1D: Not yet ready for general contraction.'
        call errquit('HF1D: Not yet ready for general contraction.',
     &               0,UERR)
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
 
       write(*,*) 'HF1D 1:Insufficient scratch space.'
       write(*,*) '        needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '        allocated ',maxW0
 
       write(*,*) 'From the left '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'From the right '
       write(*,*) 'ESp   :  ',i_ESp
 
       call errquit('HF1D:  Insufficient scratch space.',i_right-i_left,
     &              UERR)
 
      end if
 
      MaxMem = 1    ! take care of compiler warnings
      if (dryrun) then
        MaxMem = i_left + (maxW0 - (i_right-1))
        NPP = NPA*NPB
      else
        call hfset(Axyz,Aprims,Acoefs,NPA,NCA,
     &             Bxyz,Bprims,Bcoefs,NPB,NCB,
     &             GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)
      endif

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_pf   = i_Ep     + 3*NPP*(MXD+1)*lprod
      i_left = i_pf     + 2*NPP - 1

      if( i_left.ge.i_right )then

       write(*,*) 'HF1D 2:Insufficient scratch space.'
       write(*,*) '        needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '        allocated ',maxW0

       write(*,*) 'From the right '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'pf    :  ',i_pf, npp,mxd,lprod
       write(*,*) 'From the left '
       write(*,*) 'ESp   :  ',i_ESp

       call errquit('HF1D:  Insufficient scratch space.',i_right-i_left,
     &              UERR)

      end if

      if (dryrun) then
        MaxMem = max(MaxMem, (i_left+(maxW0-(i_right-1))))
      else
        do 100 nd = 0,MXD

          call hfmke(Axyz,Bxyz,W0(i_ALPHAp),W0(i_ESp),W0(i_Ep),
     &           W0(i_pf),nd,NPP,MXD,La+Li,Lb+Li)

100     continue
      endif

       
       
c Compute nuclear attraction integrals, <a|V|b>.

c     if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

       i_R0  = i_Ep  + 3*NPP*(MXD+1)*lprod
       i_IJK = i_R0  + NPP*Lp3
       i_R0C = i_IJK + (Lp+1)**3
       i_P   = i_R0C + NPP*Lp3*1
       i_RS  = i_P   + NPP*3
       i_PC  = i_RS  + NPP
       i_ff  = i_PC  + NPP*3
       i_Rj  = i_ff  + NPP*2
       i_top = i_Rj  + NPP*(Lp+1)*Lp3 - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1D 3:Insufficient scratch space.'
        write(*,*) '        needed    ',i_top
        write(*,*) '        allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'R0    :  ',i_R0
        write(*,*) 'IJK   :  ',i_IJK
        write(*,*) 'R0C   :  ',i_R0C
        write(*,*) 'P     :  ',i_P
        write(*,*) 'RS    :  ',i_RS
        write(*,*) 'PC    :  ',i_PC
        write(*,*) 'ff    :  ',i_ff
        write(*,*) 'Rj    :  ',i_Rj

        call errquit('HF1D:  Insufficient scratch space.',maxW0-i_top,
     &               UERR)

       end if

       if (dryrun) then
         MaxMem = max(MaxMem, i_top)
       else
         call hf1mkr(Axyz,Bxyz,Cxyz,zan,0.0d0,1,
     &          W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &          W0(i_Rj),W0(i_R0),W0(i_R0C),W0(i_IJK),
     &          NPP,Lp,Lp3,.TRUE.)
c
c        Compute the 1-electron - nuclear attraction contributions of
c        the gradient wrt the nuclear coordinates
c
c        call hfefi(W0(i_Ep),W0(i_R0C),W0(i_IJK),bNAI,
c    &          NPP,Nint,La,Lb,Li,Lp,Lp3,1,
c    &          MXD,canAB,ictra,ictrb)

       endif
c Compute compute the derivative wrt to the centers of the Gaussian functions.

c Assign scratch space.

       i_VP  = i_IJK + (Lp+1)**3
       i_VR  = i_VP  + NPP*(Nint*3)
       i_VA  = i_VR  + NPP*(Nint*3)
       i_VB  = i_VA  + (Nint*3)
       i_ff  = i_VB  + (Nint*3)
       i_top = i_ff  + NPP*2 - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1D 4:Insufficient scratch space.'
        write(*,*) '        needed    ',i_top
        write(*,*) '        allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'R0    :  ',i_R0
        write(*,*) 'IJK   :  ',i_IJK
        write(*,*) 'VP    :  ',i_VP
        write(*,*) 'VR    :  ',i_VR
        write(*,*) 'VA    :  ',i_VA
        write(*,*) 'VB    :  ',i_VB
        write(*,*) 'ff    :  ',i_ff

        call errquit('HF1D:  Insufficient scratch space.',maxW0-i_top,
     &               UERR)

       end if

       if (dryrun) then
         MaxMem = max(MaxMem, i_top)
       else
c Compute the derivatives wrt to (P,R).

c Compute the derivatives of the primitive integrals.

         Nder = 6

         call hfdnai(W0(i_Ep),W0(i_R0),W0(i_IJK),W0(i_VP),
     &          NPP,Nint,La,Lb,Li,Lp,Lp3,1,
     &          MXD,inder2,Nder,canAB)
         
c Transform to obtain derivatives wrt (A,B).
         
         call hf1PRtoAB(W0(i_VP),W0(i_VR),W0(i_VA),W0(i_VB),
     &          W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ff),NPP,Nint*3,
     &          ictrA,ictrB)
         
c Combine the two components of the derivatives of NAIs.
c
c   (1) the derivative wrt to each of the nuclear attraction centers,
c   (2) the derivative wrt to the centers of the Gaussian functions.
         
         if( ictrA.eq.ictrB )then
           call daxpy(Nint*3,1.D0,W0(i_VA),1,bNAI(1,1),1) 
           call daxpy(Nint*3,-1.D0,W0(i_VA),1,bNAI(1,7),1) 
         else
           call daxpy(Nint*3,1.D0,W0(i_VA),1,bNAI(1,1),1)
           call daxpy(Nint*3,1.D0,W0(i_VB),1,bNAI(1,4),1)
           call daxpy(Nint*3,-1.D0,W0(i_VA),1,bNAI(1,7),1)
           call daxpy(Nint*3,-1.D0,W0(i_VB),1,bNAI(1,7),1)
         end if
         
       endif
c     end if
      
c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
