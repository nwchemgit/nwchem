      subroutine exactd_mem(bases,nbas)
c $Id: exactd_mem.F,v 1.6 1994-08-05 20:06:09 d3e129 Exp $
      implicit none
#include "bas.fh"
#include "apiP.fh"
c
c:: passed
      integer nbas
      integer bases(nbas)
c:: local
      integer lmem_2e4c, lmem_1e
      integer lsz_2e4c, lsz_1e
      integer ibas,jbas
      integer nbfmax,nbftest
c
      mem_2e4c = 0
      isz_2e4c = 0
      mem_1e   = 0
      isz_1e   = 0
c
      nbfmax = 0
      do 00100 ibas = 1,nbas
        do 00200 jbas = 1,ibas
c.... nbfmax for 2e3c and 2e4c
           if (.not.bas_nbf_cn_max(bases(ibas),nbftest))
     &          call errquit('exactd_mem:bas_nbf_cn_max error ',1)
           nbfmax = max(nbfmax,nbftest)
           if (.not.bas_nbf_cn_max(bases(jbas),nbftest))
     &          call errquit('exactd_mem:bas_nbf_cn_max error ',2)
           nbfmax = max(nbfmax,nbftest)
c.... 2e 4 center
           call emem_d2e4c(bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
           mem_2e4c = max(mem_2e4c, lmem_2e4c)
           isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 1e
           call emem_d1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
           mem_1e = max(mem_1e, lmem_1e)
           isz_1e = max(isz_1e, lsz_1e)
00200   continue
00100 continue
c     
c.. need extra memory for H1 call externalized by int_mem_h1 (Ta)
c
      mem_1e_h1 = mem_1e + nbfmax*nbfmax*3*2
c
c... 2e3c and 2e2c use same scratch space so the maxg is the only thing that changes      
c
*      isz_2e3c = nbfmax**3
*      isz_2e2c = nbfmax**2
*      mem_2e3c = mem_2e4c
*      mem_2e2c = mem_2e4c
c
      end
      subroutine emem_d1e(brain, ketin, max1e, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
c::functions
      integer int_unint
      external int_unint
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer max1e   ! [output] max size of 1e integrals buffer
      integer memsize ! [output] max size of scratch space for 1e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc, nint
      integer bra_type, ket_type, memdum, bra_geom, ket_geom
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket, nat
c
      memsize = 0
      max1e = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(6,*)'emem_d1e: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d1e: geom error ',911)
      endif
c
      nat = ncenter(bra_geom)
c
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
c
        if (bra.eq.ket) then
          ket_lim = ibc
        else
          ket_lim = ket_ucont
        endif
c
        do 00200 ikc = 1, ket_lim
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 9 000 000
          nint = int_unint(brain,ibc,ketin,ikc,0,0,0,0)
          call hf1d(
     &         exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,1,
     &         exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,2,
     &         exndcf,exndcf,nat,
     &         exndcf,exndcf,exndcf,nint,.true.,.true.,.true.,.false.,
     &         .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_d1e:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_d1e:bas_nbf_cn_max error ',2)
      max1e = max_nbf_bra*max_nbf_ket*nat*3
c
*     write(6,'(a,i10)')'emem_d1e: memsize =',memsize
*     write(6,'(a,i10)')'emem_d1e: max1e   =',max1e
c
      end
      subroutine emem_d2e4c(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geobasmapP.fh"
#include "geomP.fh"
c::functions
      integer int_unint
      external int_unint
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc, nint
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
      integer bra_geom, ket_geom, nat
c
      memsize = 0
      maxg = 0
c
      bra = brain + BASIS_HANDLE_OFFSET
      ket = ketin + BASIS_HANDLE_OFFSET
c
      bra_ucont = infbs_head(HEAD_NCONT,bra)
      ket_ucont = infbs_head(HEAD_NCONT,ket)
c
      bra_geom  = ibs_geom(bra)
      ket_geom  = ibs_geom(ket)
      if (bra_geom.ne.ket_geom) then
        write(6,*)'emem_d1e: two different geometries for',
     &         ' derivatives?'
        call errquit('emem_d2e4c: geom error ',911)
      endif
c
      nat = ncenter(bra_geom)
c
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
c
        if (bra.eq.ket) then
          ket_lim = ibc
        else
          ket_lim = ket_ucont
        endif
c
        do 00200 ikc = 1, ket_lim
          ket_type = infbs_cont(CONT_TYPE,ikc,ket)
          ket_prim = infbs_cont(CONT_NPRIM,ikc,ket)
          ket_gen  = infbs_cont(CONT_NGEN,ikc,ket)
          if (ket_type.lt.0) then
            ket_type = abs(ket_type)
            ket_gen  = 1
          endif
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
          call hf2d(
     &           exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,1,
     &           exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,2,
     &           exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,3,
     &           exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,4,
     &           exndcf,nint,.false.,.false.,.false.,
     &           .true.,exndcf,memdum)
*rak:-start-debug
*rak:          write(6,*)
*rak:          write(6,*)'===============================================',
*rak:     &        '================================='
*rak:          write(6,10000)ibc,ibc,ikc,ikc,
*rak:     &        bra_prim,bra_prim,ket_prim,ket_prim,
*rak:     &        nint,memdum
*rak:10000     format(1x,'<',4(i4),'>(',4(i3),') nint=',i10,' memory =',i10)
*rak:-end-debug
          memsize = max(memsize,memdum)
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
          nint = int_unint(brain,ibc,brain,ibc,ketin,ikc,ketin,ikc)
          call hf2d(
     &           exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,1,
     &           exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,2,
     &           exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,3,
     &           exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,4,
     &           exndcf,nint,.false.,.false.,.false.,
     &           .true.,exndcf,memdum)
*rak:-start-debug
*rak:          write(6,10000)ibc,ikc,ibc,ikc,
*rak:     &        bra_prim,ket_prim,bra_prim,ket_prim,
*rak:     &        nint,memdum
*rak:-end-debug
          memsize = max(memsize,memdum)
c
          if (bra.eq.ket) then
            memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
            nint = int_unint(brain,ibc,brain,ibc,brain,ibc,ketin,ikc)
            call hf2d(
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,1,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,2,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,3,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,4,
     &             exndcf,nint,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
*rak:-start-debug
*rak:          write(6,10000)ibc,ibc,ibc,ikc,
*rak:     &        bra_prim,bra_prim,bra_prim,ket_prim,
*rak:     &        nint,memdum
*rak:-end-debug
            memsize = max(memsize,memdum)
c
            memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bk|kk)
            nint = int_unint(brain,ibc,ketin,ikc,ketin,ikc,ketin,ikc)
            call hf2d(
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,1,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,2,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,3,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,4,
     &             exndcf,nint,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
*rak:-start-debug
*rak:          write(6,10000)ibc,ikc,ikc,ikc,
*rak:     &        bra_prim,ket_prim,ket_prim,ket_prim,
*rak:     &        nint,memdum
*rak:          write(6,*)'===============================================',
*rak:     &        '================================='
*rak:-end-debug
            memsize = max(memsize,memdum)
          endif
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_d2e4c:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_d2e4c:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
      maxg = maxg*12
      
c
*rak:      write(6,'(a,i10)')'emem_d2e4c: memsize =',memsize
*rak:      write(6,'(a,i10)')'emem_d2e4c: maxg    =',maxg
c
      end

