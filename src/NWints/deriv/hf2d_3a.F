      subroutine hf2d_3a(
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
      integer npa1, npa2, npa3
      integer a_exp, a_cof
      integer maxw0new, k0s, icont, kaddint, pw0
      integer mem1add, mem1
      integer mem2add, mem2
      integer mem3add, mem3
c
      npa1  = npa/3
      npa2  = npa1
      if (mod(npa,3).eq.1) then
        npa1 = npa1 + 1
      elseif(mod(npa,3).eq.2) then
        npa1 = npa1 + 1
        npa2 = npa1
      endif
      npa3  = npa - npa1 - npa2        ! size of third block
c... do first block     (a1|*)
      a_exp = 1                 ! exponent block
      a_cof = a_exp + npa1      ! coeff block
      k0s  = a_cof + npa1*NCA   ! start of scratch
      maxw0new = maxw0 - k0s    ! new scratch size
      call dcopy(npa1,Aprims,1,W0(a_exp),1) ! copy exponents

      mem1    = k0s         ! take care of compiler warnings
      mem2    = k0s
      mem1add = mem1
      mem2add = mem2
      mem3add = mem2

      if (dryrun) then
        mem1add = k0s
      else
        do 00100 icont = 1,NCA
          pw0 = a_cof+(icont-1)*npa1
          call dcopy(npa1,Acoefs(1,icont),1,W0(pw0),1) ! copy coeffs
00100   continue
      endif
c
      call hf2dold(
     &    Axyz,W0(a_exp),W0(a_cof),NPA1,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0(k0s),maxW0new)
      if (dryrun) mem1 = maxW0new + mem1add
c... do second block  (a2|*)
      a_exp    = 1
      a_cof    = a_exp + npa2
      kaddint = a_cof + npa2*NCA
      k0s     = kaddint + 12*nint
      maxw0new = maxw0 - k0s
      if (dryrun) then
        mem2add = k0s
      else
        call dcopy(npa2,Aprims(npa1+1),1,W0(a_exp),1)
        do 00200 icont = 1,NCA
          pw0 = a_cof+(icont-1)*npa2
          call dcopy(npa2,Acoefs((npa1+1),icont),1,W0(pw0),1)
00200   continue
      endif
      call hf2dold(
     &    Axyz,W0(a_exp),W0(a_cof),NPA2,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    W0(kaddint),Nint,canAB,canCD,canPQ,dryrun,
     &    W0(k0s),maxW0new)
      if (dryrun) then
        mem2 = maxW0new + mem2add
      else
c... sum 1,2 blocks
        call daxpy((nint*12),1.0d00,w0(kaddint),1,beri,1)
      endif
c... do third block  (a2|*)
      a_exp    = 1
      a_cof    = a_exp + npa3
      kaddint = a_cof + npa3*NCA
      k0s     = kaddint + 12*nint
      maxw0new = maxw0 - k0s
      if (dryrun) then
        mem3add = k0s
      else
        call dcopy(npa3,Aprims(npa1+npa2+1),1,W0(a_exp),1)
        do 00300 icont = 1,NCA
          pw0 = a_cof+(icont-1)*npa3
          call dcopy(npa3,Acoefs((npa1+npa2+1),icont),1,W0(pw0),1)
00300   continue
      endif
      call hf2dold(
     &    Axyz,W0(a_exp),W0(a_cof),NPA3,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    W0(kaddint),Nint,canAB,canCD,canPQ,dryrun,
     &    W0(k0s),maxW0new)
      if (dryrun) then
        mem3 = maxw0new + mem3add
        maxw0 = max(mem1,mem2,mem3)
      else
c... sum both blocks
        call daxpy((nint*12),1.0d00,w0(kaddint),1,beri,1)
      endif
      end
