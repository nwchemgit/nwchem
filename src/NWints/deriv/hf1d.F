C> Old hf1d is now hf1dQ (below).  This hf1d is a wrapper
C> for backward compatibility.  hf1dQ has a new input argument,
C> hlfdrv, required for computing the half-derivative overlap matrix used
C> by COLUMBUS to compute non-adiabatic derivative couplings.  Routines
C> which need to calculate the HDOM should call hf1dQ directly,
C> specifying hlfdrv=.true.
C>
C> Lachlan Belcher 14 Feb 2023


      SUBROUTINE hf1d(
     &       Axyz,Aprims,Acoefs,NPA,NCA,La,ictrA,
     &       Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrB,
     &       Cxyz,zan,exinv,ncenters,
     &       bO2I,bKEI,bNAI,Nint,
     &       O2I,KEI,NAI,canAB,dryrun,W0,maxW0)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)
      logical hlfdrv    ! [input] if true, calculate same atom as well
                        ! as different atom integrals. This is necessary
                        ! when calculating the half-derivative overlap
                        ! matrix for NADCTs
      integer ictrA, ictrB !NB already used implicitly

      Logical dryrun
      Logical O2I,KEI,NAI,canAB,GenCon

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates, Charges & Inverse Exponents

      Dimension Cxyz(3,ncenters),zan(ncenters),exinv(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      Dimension bO2I(Nint,*),bKEI(Nint,*),bNAI(Nint,*)


      hlfdrv=.false.
      CALL hf1dQ(
     &       Axyz,Aprims,Acoefs,NPA,NCA,La,ictrA,
     &       Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrB,
     &       Cxyz,zan,exinv,ncenters,
     &       bO2I,bKEI,bNAI,Nint,
     &       O2I,KEI,NAI,canAB,dryrun,W0,maxW0,
     &       hlfdrv)



      end



C-------------------------------------------------------------------------
      Subroutine hf1dQ(
     &       Axyz,Aprims,Acoefs,NPA,NCA,La,ictrA,
     &       Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrB,
     &       Cxyz,zan,exinv,ncenters,
     &       bO2I,bKEI,bNAI,Nint,
     &       O2I,KEI,NAI,canAB,dryrun,W0,maxW0,
     &       hlfdrv)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

#include "errquit.fh"
#include "stdio.fh"
      LOGICAL smatm     ! if true, both integrals use the same atom
                        ! center.  Necessary for knowing when to
                        ! multiply by the 1/2 factor for half-derivative
                        ! integrals
      logical hlfdrv    ! [input] if true, calculate same atom as well
                        ! as different atom integrals. This is necessary
                        ! when calculating the half-derivative overlap
                        ! matrix for NADCTs
      integer ictrA, ictrB !NB already used implicitly

      Logical dryrun
      Logical O2I,KEI,NAI,canAB,GenCon

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoefs(NPA,NCA)
      Dimension Bxyz(3),Bprims(NPB),Bcoefs(NPB,NCB)

c--> Nuclear Cartesian Coordinates, Charges & Inverse Exponents

      Dimension Cxyz(3,ncenters),zan(ncenters),exinv(ncenters)

c--> Blocks of Overlap, Kinetic Energy & Nuclear Attraction Integrals

      Dimension bO2I(Nint,*),bKEI(Nint,*),bNAI(Nint,*)

c--> Derivative Indices

      Dimension inder1(3,3),inder2(6,6)

c--> Scratch Space.

      Dimension W0(maxW0)

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
      smatm=.FALSE.
      IF(ictrA.EQ.ictrB ) smatm = .TRUE.

c Determine whether general or segmented contraction is used.

      NCP = NCA*NCB

      GenCon = NCP.ne.1

      if( GenCon )then
        write(*,*) 'HF1D: Not yet ready for general contraction.'
        call errquit('HF1D: Not yet ready for general contraction.',
     &               0,UERR)
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
 
      i_ALPHAp = 1
      i_IPAIRp = i_ALPHAp + 2*(NPA*NPB)
      i_left   = i_IPAIRp + 2*(NPA*NPB) - 1
 
      i_ESp   = (maxW0+1) - 3*(NPA*NPB)
      i_right = i_ESp
 
      if( i_left.ge.i_right )then
 
       write(*,*) 'HF1D:  Insufficient scratch space.'
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
     &         Bxyz,Bprims,Bcoefs,NPB,NCB,
     &         GenCon,W0(i_ALPHAp),W0(i_IPAIRp),W0(i_ESp),NPP)
      endif

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      lprod = ((La+Li)+(Lb+Li)+1)*((La+Li)+1)*((Lb+Li)+1)

      i_Ep   = i_IPAIRp + 2*(NPA*NPB)
      i_pf   = i_Ep     + 3*NPP*(MXD+1)*lprod
      i_left = i_pf     + 2*NPP - 1

      if( i_left.ge.i_right )then

       write(*,*) 'HF1D:  Insufficient scratch space.'
       write(*,*) '        needed    ',i_left + (maxW0-(i_right-1))
       write(*,*) '        allocated ',maxW0

       write(*,*) 'From the right '
       write(*,*) 'ALPHAp:  ',i_ALPHAp
       write(*,*) 'IPAIRp:  ',i_IPAIRp
       write(*,*) 'Ep    :  ',i_Ep
       write(*,*) 'pf    :  ',i_pf
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

c Compute the 2-center overlap integrals, <a|S|b>.

      if( O2I.and.(.not.dryrun))then

       Nder = 3

       call hfd2oi(W0(i_Ep),bO2I,NPP,Nint,La,Lb,Li,MXD,inder1,Nder,
     &             canAB)

       do 200 n = 1,Nint
!        WRITE(*,*)"LB in hf1d, smatm=",smatm,"hlfdrv=",hlfdrv
!        WRITE(*,*)"ictrA=",ictrA,"ictrB=",ictrB
        IF ((.NOT. smatm) .AND. (hlfdrv)) THEN
          bO2I(n,1)=-5D-01*bO2I(n,1)
          bO2I(n,2)=-5D-01*bO2I(n,2)
          bO2I(n,3)=-5D-01*bO2I(n,3)
        ENDIF
        bO2I(n,4) = -bO2I(n,1)
        bO2I(n,5) = -bO2I(n,2)
        bO2I(n,6) = -bO2I(n,3)
  200  continue

#if defined(INTDEBUG)
       write(6,*)' ints from hf1d '
       call output(bO2i,1,nint,1,6,nint,6,1)
#endif
      end if
       
c Compute kinetic energy integrals, <a|T|b>.

      if( KEI )then

c Assign pointers to scratch space.

       i_Ti  = i_Ep + 3*NPP*(MXD+1)*lprod
       i_top = i_Ti + NPP - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1D:  Insufficient scratch space.'
        write(*,*) '        needed    ',i_top
        write(*,*) '        allocated ',maxW0

        write(*,*) 'ALPHAp:  ',i_ALPHAp 
        write(*,*) 'IPAIRp:  ',i_IPAIRp
        write(*,*) 'Ep    :  ',i_Ep
        write(*,*) 'Ti    :  ',i_Ti

        call errquit('HF1D:  Insufficient scratch space.',maxW0-i_top,
     &               UERR)

       end if

       if (dryrun) then
         MaxMem = max( MaxMem, i_top)
       else
         Nder = 3

         call hfdkei(W0(i_ALPHAp),W0(i_Ep),bKEI,W0(i_Ti),
     &          NPP,Nint,La,Lb,Li,MXD,inder1,Nder,canAB)

         do 300 n = 1,Nint
           bKEI(n,4) = -bKEI(n,1)
           bKEI(n,5) = -bKEI(n,2)
           bKEI(n,6) = -bKEI(n,3)
300      continue

       endif
      end if
       
c Compute nuclear attraction integrals, <a|V|b>.

      if( NAI )then

c Define the auxiliary function integrals.

c Assign scratch space.

       i_R0  = i_Ep  + 3*NPP*(MXD+1)*lprod
       i_IJK = i_R0  + NPP*Lp3
       i_R0C = i_IJK + (Lp+1)**3
       i_P   = i_R0C + NPP*Lp3*ncenters
       i_RS  = i_P   + NPP*3
       i_PC  = i_RS  + NPP
       i_ff  = i_PC  + NPP*3
       i_Rj  = i_ff  + NPP*2
       i_top = i_Rj  + NPP*(Lp+1)*Lp3 - 1

       if( i_top.gt.maxW0 )then

        write(*,*) 'HF1D:  Insufficient scratch space.'
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
         call hf1mkr(Axyz,Bxyz,Cxyz,zan,exinv,ncenters,
     &          W0(i_ALPHAp),W0(i_P),W0(i_RS),W0(i_PC),W0(i_ff),
     &          W0(i_Rj),W0(i_R0),W0(i_R0C),W0(i_IJK),
     &          NPP,Lp,Lp3,.TRUE.)

c Compute the derivative wrt to each of the nuclear attraction centers.
c This accounts for the operator dependence on the the nuclear coordinates.

         call hfefi(W0(i_Ep),W0(i_R0C),W0(i_IJK),bNAI,
     &          NPP,Nint,La,Lb,Li,Lp,Lp3,ncenters,
     &          MXD,canAB,ictra,ictrb)

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

        write(*,*) 'HF1D:  Insufficient scratch space.'
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
     &          NPP,Nint,La,Lb,Li,Lp,Lp3,ncenters,
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
           call daxpy(Nint*3,1.D0,W0(i_VA),1,bNAI(1,((ictrA-1)*3+1)),1)
         else
           call daxpy(Nint*3,1.D0,W0(i_VA),1,bNAI(1,((ictrA-1)*3+1)),1)
           call daxpy(Nint*3,1.D0,W0(i_VB),1,bNAI(1,((ictrB-1)*3+1)),1)
         end if
         
       endif
      end if
      
c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun ) maxW0 = MaxMem
c
      end
