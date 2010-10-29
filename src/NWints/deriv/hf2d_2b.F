      subroutine hf2d_2b(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0,maxW0)
C $Id$
      Implicit None
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
      integer npb1, npb2
      integer b_exp, b_cof
      integer maxw0new, k0s, icont, pw0, kaddint
      integer mem1, mem1add
      integer mem2, mem2add
c
      npb1  = npb/2
      if (mod(npb,2).eq.1) npb1 = npb1 + 1
      npb2  = npb - npb1        ! size of second block
c... do first block
      b_exp = 1                 ! exponent block
      b_cof = b_exp + npb1      ! coeff block
      k0s  = b_cof + npb1*NCB   ! start of scratch
      maxw0new = maxw0 - k0s    ! new scratch size

      mem1add = k0s             ! take care of compiler warnings
      mem1    = mem1add
      mem2add = mem1

      if (dryrun) then
        mem1add = k0s
      else
        call dcopy(npb1,Bprims,1,W0(b_exp),1) ! copy exponents
        do 00100 icont = 1,NCB
          pw0 = b_cof+(icont-1)*npb1
          call dcopy(npb1,Bcoefs(1,icont),1,W0(pw0),1) ! copy coeffs
00100   continue
      endif
c
      call hf2dold(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,W0(b_exp),W0(b_cof),NPB1,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0(k0s),maxW0new)
      if (dryrun) mem1 = maxW0new + mem1add
c... do second block
      b_exp    = 1
      b_cof    = b_exp + npb2
      kaddint = b_cof + npb2*NCB
      k0s     = kaddint + 12*nint
      maxw0new = maxw0 - k0s
      if (dryrun) then
        mem2add = k0s
      else
        call dcopy(npb2,Bprims(npb1+1),1,W0(b_exp),1)
        do 00200 icont = 1,NCB
          pw0 = b_cof+(icont-1)*npb2
          call dcopy(npb2,Bcoefs((npb1+1),icont),1,W0(pw0),1)
00200   continue
      endif
      call hf2dold(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,W0(b_exp),W0(b_cof),NPB2,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    W0(kaddint),Nint,canAB,canCD,canPQ,dryrun,
     &    W0(k0s),maxW0new)
      if (dryrun) then
        mem2 = maxw0new + mem2add
        maxw0 = max(mem1,mem2)
      else
c... sum both blocks
        call daxpy((nint*12),1.0d00,w0(kaddint),1,beri,1)
      endif
      end
