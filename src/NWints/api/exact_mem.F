      subroutine exact_mem(bases,nbas)
      implicit none
#include "apiP.fh"
c
c:: passed
      integer nbas
      integer bases(nbas)
c:: local
      integer lmem_2e4c, lmem_1e
      integer lsz_2e4c, lsz_1e
      integer ibas,jbas
c
      mem_2e4c = 0
      isz_2e4c = 0
      mem_1e   = 0
      isz_1e   = 0
c
      do 00100 ibas = 1,nbas
        do 00200 jbas = 1,ibas
c.... 2e 4 center
          call emem_2e4c(ibas, jbas, lsz_2e4c, lmem_2e4c)
          mem_2e4c = max(mem_2e4c, lmem_2e4c)
          isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 1e
          call emem_1e(ibas, jbas, lsz_1e, lmem_1e)
          mem_1e = max(mem_1e, lmem_1e)
          isz_1e = max(isz_1e, lsz_1e)
00200   continue
00100 continue
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
      
      status = bas_numcont(brain,bra_ucont,.true.)
      status = status .and. bas_numcont(ketin,ket_ucont,.true.)
      if (.not.status) then
        write(6,*)'emem_2e4c: number of unique contractions error '
        stop 'error'
      endif
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
      max1e = angm*(angm+1)/2
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
      
      status = bas_numcont(brain,bra_ucont,.true.)
      status = status .and. bas_numcont(ketin,ket_ucont,.true.)
      if (.not.status) then
        write(6,*)'emem_2e4c: number of unique contractions error '
        stop 'error'
      endif
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
      maxg = angm*(angm+1)/2
      maxg = maxg**4
c
      write(6,'(a,i10)')'emem_2e4c: memsize =',memsize
      write(6,'(a,i10)')'emem_2e4c: maxg    =',maxg
c
      end

