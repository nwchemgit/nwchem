      subroutine exact_mem(bases,nbas)
c $Id: exact_mem.F,v 1.18 1995-06-07 18:11:21 d3e129 Exp $
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
      integer lmem_2e3c, lsz_2e3c
      integer lmem_2e2c, lsz_2e2c
      integer ibas,jbas, kbas
*      integer nbfmax,nbftest
      integer lsz_3ov,lmem_3ov
c
c
c... ensure that mem_2e{432}c memory reqirements are nonzero
c
      mem_2e2c = 2
      isz_2e2c = 2
      mem_2e3c = 2
      isz_2e3c = 2
      mem_2e4c = 2
      isz_2e4c = 2
      mem_1e   = 2
      isz_1e   = 2
c
*      nbfmax = 0
      do 00100 ibas = 1,nbas
        do 00200 jbas = 1,ibas
*c.... nbfmax for 2e3c and 2e4c
*           if (.not.bas_nbf_cn_max(bases(ibas),nbftest))
*     &          call errquit('exact_mem:bas_nbf_cn_max error ',1)
*           nbfmax = max(nbfmax,nbftest)
*           if (.not.bas_nbf_cn_max(bases(jbas),nbftest))
*     &          call errquit('exact_mem:bas_nbf_cn_max error ',2)
*           nbfmax = max(nbfmax,nbftest)
c.... 2e 4 center
*debug:    call emem_2e4c_full
*debug:&           (bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
           call emem_2e4c(bases(ibas),bases(jbas),lsz_2e4c,lmem_2e4c)
           mem_2e4c = max(mem_2e4c, lmem_2e4c)
           isz_2e4c = max(isz_2e4c, lsz_2e4c)
c.... 2e 3 center
           call emem_2e3c(bases(ibas),bases(jbas),lsz_2e3c,lmem_2e3c)
           mem_2e3c = max(mem_2e3c, lmem_2e3c)
           isz_2e3c = max(isz_2e3c, lsz_2e3c)
c.... 2e 2 center
           call emem_2e2c(bases(ibas),bases(jbas),lsz_2e2c,lmem_2e2c)
           mem_2e2c = max(mem_2e2c, lmem_2e2c)
           isz_2e2c = max(isz_2e2c, lsz_2e2c)
c.... 1e
           call emem_1e(bases(ibas), bases(jbas), lsz_1e, lmem_1e)
           mem_1e = max(mem_1e, lmem_1e)
           isz_1e = max(isz_1e, lsz_1e)
00200   continue
00100 continue
c
c new label routines require additional isz_1e
c
      mem_1e = mem_1e + isz_1e
c     
c.. need extra memory for H1 call externalized by int_mem_h1
c
      mem_1e_h1 = mem_1e + isz_1e
c
c... now do loops for 3 center overlap integrals
c
c... add control mechanism for 3ov?
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
c new label routines require additional integral block size
c
      mem_2e4c  = mem_2e4c        +       isz_2e4c
      mem_2e3c  = mem_2e3c        +       isz_2e3c
      mem_2e2c  = mem_2e2c        +       isz_2e2c
      mem_1e    = mem_1e          +       isz_1e
      mem_1e_h1 = mem_1e_h1       +       2*isz_1e
      mem_3ov   = mem_3ov         +       isz_3ov
*obsolite due to new label routines:-start
*..: c
*..: c... sanity checks for lable routines       
*..: c
*..:       mem_2e4c  = max(mem_2e4c,  isz_2e4c)
*..:       mem_2e3c  = max(mem_2e3c,  isz_2e3c)
*..:       mem_2e2c  = max(mem_2e2c,  isz_2e2c)
*..:       mem_1e    = max(mem_1e,    isz_1e)
*..:       mem_1e_h1 = max(mem_1e_h1, (2*isz_1e))
*..:       mem_3ov   = max(mem_3ov,   isz_3ov)
*obsolite due to new label routines:-end
c.. debug write statements
*rak:      write(6,*) ' mem_2e4c =',mem_2e4c 
*rak:      write(6,*) ' mem_2e3c =',mem_2e3c 
*rak:      write(6,*) ' mem_2e2c =',mem_2e2c 
*rak:      write(6,*) ' mem_1e   =',mem_1e   
*rak:      write(6,*) ' mem_1e_h1=',mem_1e_h1
*rak:      write(6,*) ' mem_3ov  =',mem_3ov  
*rak:      write(6,*)' '
*rak:      write(6,*)' isz_2e4c =',isz_2e4c
*rak:      write(6,*)' isz_2e3c =',isz_2e3c
*rak:      write(6,*)' isz_2e2c =',isz_2e2c
*rak:      write(6,*)' isz_1e   =',isz_1e
*rak:      write(6,*)' isz_3ov  =',isz_3ov
*rak:      write(6,*)' '
*rak:      if(isz_2e2c.gt.0)stop 'tmp_stop'
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
      integer nbfimax, nbfjmax, nbfkmax
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
        i_type = abs(infbs_cont(CONT_TYPE,iic,ibas))
        i_prim = infbs_cont(CONT_NPRIM,iic,ibas)
c
c
c  ... not SP aware
c
c
c... currently hf3ois does not allow for general contractions
c... so *_gen is not needed
c        i_gen  = infbs_cont(CONT_NGEN,iic,ibas)
c
        j_lim = jt_ucont
        if (ibas.eq.jbas) j_lim = iic
c
        do 00200 ijc = 1, j_lim
          j_type = abs(infbs_cont(CONT_TYPE,ijc,jbas))
          j_prim = infbs_cont(CONT_NPRIM,ijc,jbas)
c**          j_gen  = infbs_cont(CONT_NGEN,ijc,jbas)
c
          k_lim = kt_ucont
          if (kbas.eq.jbas) k_lim = ijc
          do 00300 ikc = 1, k_lim
            k_type = abs(infbs_cont(CONT_TYPE,ikc,kbas))
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
      if (.not.bas_nbf_cn_max(ibasin,nbfimax))
     &       call errquit('emem_3ov:bas_nbf_cn_max error ibasin',2)
      if (.not.bas_nbf_cn_max(jbasin,nbfjmax))
     &       call errquit('emem_3ov:bas_nbf_cn_max error jbasin',2)
      if (.not.bas_nbf_cn_max(kbasin,nbfkmax))
     &       call errquit('emem_3ov:bas_nbf_cn_max error kbasin',2)
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
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
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
          call hf1(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,1,
     &             exndcf,exndcf,exndcf,1,.true.,.true.,.true.,.false.,
     &             .true.,exndcf,memdum)
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_1e:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_1e:bas_nbf_cn_max error ',2)
      max1e = max_nbf_bra*max_nbf_ket
c... insure scratch is large enough to hold 3*max1e
      memsize = max(3*max1e,memsize)
c
*     write(6,'(a,i10)')'emem_1e: memsize =',memsize
*     write(6,'(a,i10)')'emem_1e: max1e   =',max1e
c
      end
      subroutine emem_2e4c(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
      logical cando_sp_u
      external cando_sp_u
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra, max_nbf_ket
c
      logical bra_cando_sp
      logical ket_cando_sp
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
        if (bra_type.lt.0) then
          bra_type = abs(bra_type)
          bra_gen  = 1
        endif
        bra_cando_sp = cando_sp_u(brain,ibc,0)
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
          ket_cando_sp = cando_sp_u(ketin,ikc,0)
          if (bra_cando_sp.and.ket_cando_sp) goto 00201
c
          memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
*:debug   write(6,*)' coul ',ibc,ibc,ikc,ikc,memdum
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
*:debug   write(6,*)' exch ',ibc,ikc,ibc,ikc,memdum
          memsize = max(memsize,memdum)
c
          if (brain.eq.ketin) then
            memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . triple (bb|bk)
            call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
*:debug     write(6,*)' bbbk ',ibc,ibc,ibc,ikc,memdum
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
*:debug     write(6,*)' bkkk ',ibc,ikc,ikc,ikc,memdum
            memsize = max(memsize,memdum)
*:debug     write(6,*)' full ',memsize
          endif
00201     continue
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_2e4c:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_2e4c:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
c
c... add 10% for abbc and abcc class shells not checked for
      memsize = memsize + memsize/10
*      write(6,'(a,i10)')'emem_2e4c: memsize =',memsize
*      write(6,'(a,i10)')'emem_2e4c: maxg    =',maxg
c
*      stop 'tmp stop'
c
      end
      subroutine emem_2e3c(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
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
          call hf2(exndcf,exndcf,exndcf,1,1,0,
     &             exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
          memsize = max(memsize,memdum)
c
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_2e3c:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_2e3c:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_bra*max_nbf_ket
      maxg = max(maxg,(max_nbf_ket*max_nbf_ket*max_nbf_bra))
c
c... add 10% for safety
      memsize = memsize + memsize/10
*      write(6,'(a,i10)')'emem_2e3c: memsize =',memsize
*      write(6,'(a,i10)')'emem_2e3c: maxg    =',maxg
c
*      stop 'tmp stop'
c
      end
      subroutine emem_2e2c(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, ikc
      integer bra_type, ket_type, memdum
      integer bra_prim, ket_prim, ket_lim
      integer bra_gen, ket_gen
      integer max_nbf_bra,max_nbf_ket
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
          call hf2(exndcf,exndcf,exndcf,bra_prim,bra_gen,bra_type,
     &             exndcf,exndcf,exndcf,1,1,0,
     &             exndcf,exndcf,exndcf,ket_prim,ket_gen,ket_type,
     &             exndcf,exndcf,exndcf,1,1,0,
     &             exndcf,1,.false.,.false.,.false.,
     &             .true.,exndcf,memdum)
*          write(6,*)' emem_2e2c: memdum',memdum, memsize
          memsize = max(memsize,memdum)
c
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_2e2c:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_2e2c:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_ket
c
c... add 10% for safety
      memsize = memsize + memsize/10
*      write(6,'(a,i10)')'emem_2e2c: memsize =',memsize
*      write(6,'(a,i10)')'emem_2e2c: maxg    =',maxg
c
*      stop 'tmp stop'
c
      end
      subroutine emem_2e4c_full(brain, ketin, maxg, memsize)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
c::functions
      integer nbf_from_ucont
      external nbf_from_ucont
      logical cando_sp_u
      external cando_sp_u
c::passed
      integer brain   ! [input] bra basis set handle
      integer ketin   ! [input] ket basis set handle
      integer maxg    ! [output] max size of 2e integrals buffer
      integer memsize ! [output] max size of scratch space for 2e integral evaluation
c::local
      integer bra, ket, bra_ucont, ket_ucont, ibc, jbc, kbc, lbc
      integer memdum
      integer i_type, j_type, k_type, l_type
      integer i_prim, j_prim, k_prim, l_prim
      integer i_gen,  j_gen,  k_gen,  l_gen
      integer max_nbf_bra,max_nbf_ket
c
      logical bra_cando_sp
      logical ket_cando_sp
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
        i_type = infbs_cont(CONT_TYPE,ibc,bra)
        i_prim = infbs_cont(CONT_NPRIM,ibc,bra)
        i_gen  = infbs_cont(CONT_NGEN,ibc,bra)
        if (i_type.lt.0) then
          i_type = abs(i_type)
          i_gen  = 1
        endif
c
        do 00200 jbc = 1, ibc
          j_type = infbs_cont(CONT_TYPE,jbc,bra)
          j_prim = infbs_cont(CONT_NPRIM,jbc,bra)
          j_gen  = infbs_cont(CONT_NGEN,jbc,bra)
          if (j_type.lt.0) then
            j_type = abs(j_type)
            j_gen  = 1
          endif
          bra_cando_sp = cando_sp_u(brain,ibc,jbc)
c
          do 00300 kbc = 1, jbc
            k_type = infbs_cont(CONT_TYPE,kbc,ket)
            k_prim = infbs_cont(CONT_NPRIM,kbc,ket)
            k_gen  = infbs_cont(CONT_NGEN,kbc,ket)
            if (k_type.lt.0) then
              k_type = abs(k_type)
              k_gen  = 1
            endif
            do 00400 lbc = 1, kbc
              l_type = infbs_cont(CONT_TYPE,lbc,ket)
              l_prim = infbs_cont(CONT_NPRIM,lbc,ket)
              l_gen  = infbs_cont(CONT_NGEN,lbc,ket)
              if (l_type.lt.0) then
                l_type = abs(l_type)
                l_gen  = 1
              endif
              ket_cando_sp = cando_sp_u(ketin,kbc,jbc)
              if (bra_cando_sp.and.ket_cando_sp) goto 00201
c
              memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . coulumb
              call hf2(
     &            exndcf,exndcf,exndcf,i_prim,i_gen,i_type,
     &            exndcf,exndcf,exndcf,j_prim,j_gen,j_type,
     &            exndcf,exndcf,exndcf,k_prim,k_gen,k_type,
     &            exndcf,exndcf,exndcf,l_prim,l_gen,l_type,
     &            exndcf,1,.false.,.false.,.false.,
     &            .true.,exndcf,memdum)
              write(6,*)' coul ',ibc,jbc,kbc,lbc,memdum
              memsize = max(memsize,memdum)
c
              memdum = 9 000 000
c. . . . . . . . . . . . . . . . . . . . . . . . . . . . exchange
              call hf2(
     &            exndcf,exndcf,exndcf,i_prim,i_gen,i_type,
     &            exndcf,exndcf,exndcf,l_prim,l_gen,l_type,
     &            exndcf,exndcf,exndcf,k_prim,k_gen,k_type,
     &            exndcf,exndcf,exndcf,j_prim,j_gen,j_type,
     &            exndcf,1,.false.,.false.,.false.,
     &            .true.,exndcf,memdum)
              write(6,*)' exch ',ibc,jbc,kbc,lbc,memdum
              memsize = max(memsize,memdum)
              memsize = max(memsize,memdum)
              write(6,*)' full ',memsize
00201         continue
00400       continue
00300     continue
00200   continue
00100 continue
c
      if (.not.bas_nbf_cn_max(brain,max_nbf_bra))
     &       call errquit('emem_2e4c_full:bas_nbf_cn_max error ',2)
      if (.not.bas_nbf_cn_max(ketin,max_nbf_ket))
     &       call errquit('emem_2e4c_full:bas_nbf_cn_max error ',2)
c
      maxg = max_nbf_bra*max_nbf_bra
      maxg = maxg*max_nbf_ket*max_nbf_ket
c
      write(6,'(a,i10)')'emem_2e4c: memsize =',memsize
      write(6,'(a,i10)')'emem_2e4c: maxg    =',maxg
c
*      stop 'tmp stop'
c
      end
