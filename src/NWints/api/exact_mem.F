      subroutine exact_mem(bases,nbas)
c $Id: exact_mem.F,v 1.5 1994-04-09 15:46:21 d3e129 Exp $
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
      integer ibas,jbas, kbas
      integer hang, tang, hang2
      integer lsz_3ov,lmem_3ov
c
      mem_2e4c = 0
      isz_2e4c = 0
      mem_1e   = 0
      isz_1e   = 0
c
      hang = 0
      do 00100 ibas = 1,nbas
        do 00200 jbas = 1,ibas
c.... angular momentum for 2e3c and 2e4c
          if(.not.bas_high_angular(bases(ibas),tang))
     &           call errquit('exact_mem:ang error ',0)
          hang = max(hang,tang)
          if(.not.bas_high_angular(bases(jbas),tang))
     &           call errquit('exact_mem:ang error ',0)
          hang = max(hang,tang)
c.... 2e 4 center
          call emem_2e4c(bases(ibas), bases(jbas), lsz_2e4c, lmem_2e4c)
          mem_2e4c = max(mem_2e4c, lmem_2e4c)
          isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 1e
          call emem_1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
          mem_1e = max(mem_1e, lmem_1e)
          isz_1e = max(isz_1e, lsz_1e)
00200   continue
00100 continue
c
c.. need extra memory for H1 call externalized by int_mem_h1
c
      mem_1e_h1 = mem_1e + isz_1e
c
c... 2e3c and 2e2c use same scratch space so the maxg is the only thing that changes      
c
      hang2 = (hang+1)*(hang+2)/2
      isz_2e3c = hang2**3
      isz_2e2c = hang2**2
      mem_2e3c = mem_2e4c
      mem_2e2c = mem_2e4c
c
c... now do loops for 3 center overlap integrals
c
c... add control mechanism for 3ov??
c
      isz_3ov = 0
      mem_3ov = 0
      do 00300 ibas = 1,nbas
        do 00400 jbas = 1,ibas
          do 00500 kbas = 1,jbas
            call emem_3ov(bases(ibas),bases(jbas),bases(kbas),
     &             lsz_3ov,lmem_3ov)
            isz_3ov = max(isz_3ov,lsz_3ov)
            mem_3ov = max(mem_3ov,lmem_3ov)
00500     continue
00400   continue
00300 continue
c
      end
      subroutine emem_3ov(ibasin,jbasin,kbasin,lsz_buf,memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::function
      integer nbf_from_ucont
      external nbf_from_ucont
c::passed
      integer ibasin, jbasin, kbasin, lsz_buf, memsize
c::local
      integer ibas, jbas, kbas
      integer it_ucont, jt_ucont, kt_ucont
      integer iic, ijc, ikc
      integer i_type, j_type, k_type
      integer i_prim, j_prim, k_prim
      integer memdum, j_lim, k_lim
      integer nbfimax, nbfjmax, nbfkmax, mynbf
c
      memsize = 0
      nbfimax = 0
      nbfjmax = 0
      nbfkmax = 0
c
      ibas = ibasin + BASIS_HANDLE_OFFSET
      jbas = jbasin + BASIS_HANDLE_OFFSET
      kbas = kbasin + BASIS_HANDLE_OFFSET
c
      it_ucont = infbs_head(HEAD_NCONT,ibas)
      jt_ucont = infbs_head(HEAD_NCONT,jbas)
      kt_ucont = infbs_head(HEAD_NCONT,kbas)
c
      do 00100 iic = 1,it_ucont
        mynbf = nbf_from_ucont(iic,ibasin)
        nbfimax = max(nbfimax,mynbf)
        i_type = infbs_cont(CONT_TYPE,iic,ibas)
        i_prim = infbs_cont(CONT_NPRIM,iic,ibas)
c
c... currently hf3ois does not allow for general contractions
c... so *_prim is not needed
c        i_gen  = infbs_cont(CONT_NGEN,iic,ibas)
c
        j_lim = jt_ucont
        if (ibas.eq.jbas) j_lim = iic
c
        do 00200 ijc = 1, j_lim
          mynbf = nbf_from_ucont(ijc,jbasin)
          nbfjmax = max(nbfjmax,mynbf)
          j_type = infbs_cont(CONT_TYPE,ijc,jbas)
          j_prim = infbs_cont(CONT_NPRIM,ijc,jbas)
c**          j_gen  = infbs_cont(CONT_NGEN,ijc,jbas)
c
          k_lim = kt_ucont
          if (kbas.eq.jbas) k_lim = ijc
          do 00300 ikc = 1, k_lim
            mynbf = nbf_from_ucont(ikc,kbasin)
            nbfkmax = max(nbfkmax,mynbf)
            k_type = infbs_cont(CONT_TYPE,ikc,kbas)
            k_prim = infbs_cont(CONT_NPRIM,ikc,kbas)
c**            k_gen  = infbs_cont(CONT_NGEN,ikc,kbas)
            
            memdum  = 9 000 000
            call hf3OIs(
     &             exndcf,exndcf,exndcf,i_prim,i_type,
     &             exndcf,exndcf,exndcf,j_prim,j_type,
     &             exndcf,exndcf,exndcf,k_prim,k_type,
     &             exndcf,1,.false.,
     &             .true.,exndcf,memdum)
            memsize = max(memsize,memdum)
00300     continue
00200   continue
00100 continue
c
      lsz_buf = nbfimax*nbfjmax*nbfkmax
      end
      subroutine emem_1e(brain, ketin, max1e, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer max1e   ! [output] max size of 1e integrals buffer
      integer memsize ! [output] max size of scratch space for 1e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen, bra_ang, ket_ang, angm
      logical status
c
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
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
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
c
          memdum = 9 000 000
          call hf1(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,1,
     &             exndcf,exndcf,exndcf,1,.true.,.true.,.true.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
00200   continue
00100 continue
c
      status = bas_high_angular(brain,bra_ang)
      status = status .and. bas_high_angular(ketin,ket_ang)
      if (.not.status) then
        write(6,*)' angular momentum not set for basis '
        write(6,*)' ket_ang = ',ket_ang
        write(6,*)' bra_ang = ',bra_ang
        stop ' error '
      endif
      angm = max(bra_ang,ket_ang)
      max1e = (angm+1)*(angm+2)/2
      max1e = max1e*max1e
c... insure scratch is large enough to hold 3*max1e
      memsize = max(max1e,memsize)
c
      write(6,'(a,i10)')'emem_1e: memsize =',memsize
      write(6,'(a,i10)')'emem_1e: max1e   =',max1e
c
      end
      subroutine emem_2e4c(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen, bra_ang, ket_ang, angm
      logical status
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
      do 00100 ibc = 1, bra_ucont
        bra_type = infbs_cont(CONT_TYPE,ibc,bra)
        bra_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        bra_gen  = infbs_cont(CONT_NGEN,ibc,bra)
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
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bk|kk)
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
00200   continue
00100 continue
c
      status = bas_high_angular(brain,bra_ang)
      status = status .and. bas_high_angular(ketin,ket_ang)
      if (.not.status) then
        write(6,*)' angular momentum not set for basis '
        write(6,*)' ket_ang = ',ket_ang
        write(6,*)' bra_ang = ',bra_ang
        stop ' error '
      endif
      angm = max(bra_ang,ket_ang)
      maxg = (angm+1)*(angm+2)/2
      maxg = maxg**4
c
      write(6,'(a,i10)')'emem_2e4c: memsize =',memsize
      write(6,'(a,i10)')'emem_2e4c: maxg    =',maxg
c
      end

