      Subroutine hf3OIs(Axyz,Aprims,Acoef,NPA,La,
     &                  Bxyz,Bprims,Bcoef,NPB,Lb,
     &                  Cxyz,Cprims,Ccoef,NPC,Lc,
     &                  b3OI,Nint,TriDiag,
     &                  DryRun,W0,maxW0)
c $Id$

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Logical TriDiag,DryRun

c--> Cartesian Coordinates, Primitives & Contraction Coefficients

      Dimension Axyz(3),Aprims(NPA),Acoef(NPA)
      Dimension Bxyz(3),Bprims(NPB),Bcoef(NPB)
      Dimension Cxyz(3),Cprims(NPC),Ccoef(NPC)

c--> Block of 3-Center Overlap Integrals

      Dimension b3OI(Nint)

c--> Scratch Space.

      Dimension W0(maxW0)
c
c Compute 3-ctr overlap integrals (3OIs) for three shells of contracted 
c Gaussians functions.
c
c******************************************************************************

c Define the prefactor of the charge distribution.

c Assign pointers to scratch space.

      i_alpha = 1
      i_top   = i_alpha + (NPA*NPB*NPC)*4 - 1

      if((i_top.gt.maxW0).and.(.not.Dryrun))then

       write(*,*) 'HF3CTR:  Insufficient scratch space.'
       write(*,*) '         needed    ',i_top
       write(*,*) '         allocated ',maxW0

       write(*,*) 'alpha   :  ',i_alpha

       write(*,*) 'if you get this error doing higher multipoles with',
     & 'a small basis set, try modifying NWints/api/exact_mem.F:279' ! Jeff

       stop

      end if

      MaxMem = i_top    ! take care of compiler warnings

      if( DryRun )then

       MaxMem = i_top
       NABC = NPA*NPB*NPC

      else

       call hf1set3(Axyz,Aprims,Acoef,NPA,
     &              Bxyz,Bprims,Bcoef,NPB,
     &              Cxyz,Cprims,Ccoef,NPC,
     &              W0(i_alpha),NABC)

      end if

c Define the center of the charge distribution.

c Assign pointers to scratch space.

      i_E   = i_alpha + NABC*4
      i_G   = i_E     + NABC*3*(La+Lb+Lc+1)*(La+1)*(Lb+1)*(Lc+1)
      i_top = i_G     + NABC*3 - 1

      if((i_top.gt.maxW0).and.(.not.Dryrun)) then

       write(*,*) 'HF3CTR:  Insufficient scratch space.'
       write(*,*) '         needed    ',i_top
       write(*,*) '         allocated ',maxW0

       write(*,*) 'alpha   :  ',i_alpha
       write(*,*) 'E       :  ',i_E
       write(*,*) 'G       :  ',i_G

       write(*,*) 'if you get this error doing higher multipoles with',
     & 'a small basis set, try modifying NWints/api/exact_mem.F:279' ! Jeff

       stop

      end if

      if( DryRun )then

       MaxMem = max( MaxMem, i_top )

      else

       call hfctr3(Axyz,Bxyz,Cxyz,W0(i_alpha),W0(i_G),NABC)

      end if

c Define the Hermite linear expansion coefficients.

c Assign pointers to scratch space.

      i_GT    = i_G     + NABC*3
      i_ABC2I = i_GT    + NABC*3
      i_top   = i_ABC2I + NABC*3 - 1

      if( i_top .gt. maxW0 .and. .not.Dryrun)then

       write(*,*) 'HF3CTR:  Insufficient scratch space.'
       write(*,*) '         needed    ',i_top
       write(*,*) '         allocated ',maxW0

       write(*,*) 'alpha   :  ',i_alpha
       write(*,*) 'E       :  ',i_E
       write(*,*) 'G       :  ',i_G
       write(*,*) 'GT      :  ',i_GT
       write(*,*) 'ABC2I   :  ',i_ABC2I

       write(*,*) 'if you get this error doing higher multipoles with',
     & 'a small basis set, try modifying NWints/api/exact_mem.F:279' ! Jeff

       stop

      end if

      if( DryRun )then

       MaxMem = max( MaxMem, i_top )

      else

       call hf1mke3(Axyz,Bxyz,Cxyz,W0(i_alpha),W0(i_G),W0(i_GT),
     &              W0(i_ABC2I),W0(i_E),NABC,La,Lb,Lc)

      end if

c Return the maximum amount of scratch space required by a "dry run".

      if( DryRun )then
       maxW0 = MaxMem
       return
      end if

c Compute the 3-ctr OIs.

      call hfabc(W0(i_E),b3OI,NABC,La,Lb,Lc,TriDiag)

      end
