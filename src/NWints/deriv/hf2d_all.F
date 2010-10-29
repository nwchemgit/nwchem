      subroutine hf2d_split_all(
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
      integer maxw0new
      integer ngsize
      integer memscr
      integer mem_need
      integer a_p, b_p, c_p, d_p
      logical gencon
c
      gencon = nca.ne.1
      gencon = gencon.or.ncb.ne.1
      gencon = gencon.or.ncc.ne.1
      gencon = gencon.or.ncd.ne.1
      if (gencon) then
        write(*,*) 'HF2D:  Not prepared to handle general contractions.'
        stop
      endif
* this loop ignores the general contraction index.  Must change if 
* McMurchie-Davidson general contraction derivatives are ever done
*
      ngsize = nint*12
      memscr = ngsize + 1
      maxw0new = MaxW0 - ngsize
      mem_need = 0
      if (.not.dryrun) call dfill(ngsize,0.0d00,bERI,1)
      do a_p = 1, NPA
        do b_p = 1,NPB
          do c_p = 1,NPC
            do d_p = 1,NPD
              if (dryrun) then
                maxw0new = MaxW0 - ngsize
              else
                call dfill(ngsize,0.0d00,W0,1)
              endif
* compute primitive block
              call hf2dold(
     &            Axyz,Aprims(a_p),Acoefs(a_p,1),1,1,La,ictra,
     &            Bxyz,Bprims(b_p),Bcoefs(b_p,1),1,1,Lb,ictrb,
     &            Cxyz,Cprims(c_p),Ccoefs(c_p,1),1,1,Lc,ictrc,
     &            Dxyz,Dprims(d_p),Dcoefs(d_p,1),1,1,Ld,ictrd,
     &            w0,Nint,canAB,canCD,canPQ,dryrun,
     &            w0(memscr),maxW0new)
              if (dryrun) then
                mem_need = max(mem_need,(maxW0new+ngsize))
              else
* sum primitive block together
                call daxpy(ngsize,1.0d00,w0,1,beri,1)
              endif
            enddo
          enddo
        enddo
      enddo
c
      if (dryrun) then
        maxw0 = mem_need
      endif
      end
      subroutine hf2d_split_block(
     &    Axyz,Aprims,Acoefs,NPA,NCA,La,ictra,
     &    Bxyz,Bprims,Bcoefs,NPB,NCB,Lb,ictrb,
     &    Cxyz,Cprims,Ccoefs,NPC,NCC,Lc,ictrc,
     &    Dxyz,Dprims,Dcoefs,NPD,NCD,Ld,ictrd,
     &    bERI,Nint,canAB,canCD,canPQ,dryrun,
     &    W0,maxW0,ba,bb,bc,bd)
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
      Integer ba,bb,bc,bd ! block sizes for each index
c--> Block of Electron Repulsion Integrals
      
      Double Precision bERI(Nint,*)
      
c--> Scratch Space
      
      Double Precision W0(maxW0)
c
      integer maxw0new
      integer ngsize
      integer memscr
      integer mem_need
      integer apl, aph, bpl, bph, cpl, cph, dpl, dph
      integer dap, dbp, dcp, ddp
      logical gencon
c
      gencon = nca.ne.1
      gencon = gencon.or.ncb.ne.1
      gencon = gencon.or.ncc.ne.1
      gencon = gencon.or.ncd.ne.1
      if (gencon) then
        write(*,*) 'HF2D:  Not prepared to handle general contractions.'
        stop
      endif
*
* this loop ignores the general contraction index.  Must change if 
* McMurchie-Davidson general contraction derivatives are ever done
*
      ngsize = nint*12
      memscr = ngsize + 1
      maxw0new = MaxW0 - ngsize
      mem_need = 0
      if (.not.dryrun) call dfill(ngsize,0.0d00,bERI,1)
      do apl = 1, NPA, ba
        aph = min( NPA, (apl+ba-1))
        dap = aph - apl +1
        do bpl = 1, NPB, bb
          bph = min( NPB, (bpl+bb-1))
          dbp = bph - bpl + 1
          do cpl = 1, NPC, bc
            cph = min( NPC, (cpl+bc-1))
            dcp = cph - cpl + 1
            do dpl = 1, NPD, bd
              dph = min( NPD, (dpl+bd-1))
              ddp = dph - dpl + 1
              if (dryrun) then
                maxw0new = MaxW0 - ngsize
              else
                call dfill(ngsize,0.0d00,W0,1)
              endif
* compute primitive block
              call hf2dold(
     &            Axyz,Aprims(apl),Acoefs(apl,1),dap,NCA,La,ictra,
     &            Bxyz,Bprims(bpl),Bcoefs(bpl,1),dbp,NCB,Lb,ictrb,
     &            Cxyz,Cprims(cpl),Ccoefs(cpl,1),dcp,NCC,Lc,ictrc,
     &            Dxyz,Dprims(dpl),Dcoefs(dpl,1),ddp,NCD,Ld,ictrd,
     &            w0,Nint,canAB,canCD,canPQ,dryrun,
     &            w0(memscr),maxW0new)
              if (dryrun) then
                mem_need = max(mem_need,(maxW0new+ngsize))
              else
* sum primitive block together
                call daxpy(ngsize,1.0d00,w0,1,beri,1)
              endif
            enddo
          enddo
        enddo
      enddo
c
      if (dryrun) then
        maxw0 = mem_need
      endif
      end
