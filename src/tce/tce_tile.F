      subroutine tce_tile
c
c $Id: tce_tile.F,v 1.4 2002-10-22 16:23:08 sohirata Exp $
c
c     Reorder and tile orbitals.  Minimal tiling scheme is
c     hole(alpha), hole(beta), particle(alpha), particle(beta)
c     with no symmetry.
c
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "tce.fh"
c
      integer ga_dble
      integer ma_dble
      integer available
      integer isize
      double precision size
      integer sym
      integer any, hole, particle
      integer i,j,k,nblocks
      logical nodezero
      integer spin_tmp(max_size)
      integer sym_tmp(max_size)
      integer range_tmp(max_size)
      character*4 irrepname
      character*5 spinname(2)
      data spinname/'alpha','beta '/
c
c     ====================
c     Retrieve Memory Info
c     ====================
c
      nodezero = (ga_nodeid().eq.0)
      ga_dble = dfloat(ga_memory_avail()) / dfloat(bytes)
      ma_dble=dfloat(ma_inquire_avail(mt_dbl))
      call ga_dgop(msg_tdmem,ga_dble,1,'+')
      if (ga_uses_ma()) then
        available = int(ga_dble)/2
      else
        available = int(ma_dble)
      endif
      if (nodezero.and.util_print('tile',
     1  print_default)) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'Memory Information',20,.true.)
        if (ga_uses_ma()) then
          write(LuOut,9000) int(ga_dble)
        else
          write(LuOut,9010) int(ga_dble),int(ma_dble)
c
c ====================================================================
c Warning! when ga_uses_ma()=.false., "shmmax" may become the limiting
c          memory size, rather than ga_memory_avail.  In that case,
c          despite the fact that ga_memory_avail gives a memory size
c          sufficient for the calculation, an out-of-memory error can
c          occur.  Check shmmax by "cat /proc/sys/kernel/shmmax".
c ====================================================================
c
        endif
      endif
 9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
 9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
     1       10x,'Available MA space size is    ',i10,' doubles')
c
c     ========================
c     Determine the block size
c     ========================
c
      size = dfloat(available)
      size = size / dfloat(nbf*(nbf+1)/2)
      size = size / 2.0d0
      size = dsqrt(size)
      isize = int(size)
      if (model.eq."ccd") then
        size = dfloat(available)
        size = size / 2.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."lccd") then
        size = dfloat(available)
        size = size / 2.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."ccsd") then
        size = dfloat(available)
        size = size / 2.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."lccsd") then
        size = dfloat(available)
        size = size / 2.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."ccsdt") then
        size = dfloat(available)
        size = size / 2.0d0
        size = size ** (1.0d0/6.0d0)
        isize = min(isize,int(size))
      else
        call errquit("tce_tile: the model NYI",0)
      endif
      if (nodezero.and.util_print('tile',print_default)) then
        write(LuOut,*)
        write(LuOut,9040) isize
      endif
 9040 format(1x,'Maximum block size',i10,' doubles')
c
c     ================================================
c     Reorder the orbitals in the ha, hb, pa, pb order
c     ================================================
c
      if (.not.ma_push_get(mt_dbl,nbf*(nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol)),"sorted MO coeffs",
     2  l_movecs_sorted,k_movecs_sorted))
     3  call errquit("tce_tile: MA problem",0)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted spins",
     2  l_spin_sorted,k_spin_sorted))
     3  call errquit("tce_tile: MA problem",0)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted irs",
     2  l_irs_sorted,k_irs_sorted))
     3  call errquit("tce_tile: MA problem",0)
      if (.not.ma_push_get(mt_dbl,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted evl",
     2  l_evl_sorted,k_evl_sorted))
     3  call errquit("tce_tile: MA problem",0)
      nirreps = sym_number_ops(geom) + 1
      if (nodezero.and.util_print('tile',print_debug)) then
        write(LuOut,*) "Number of irreps = ",nirreps
        do any = nfc(1)+1,nmo(1)-nfv(1)
          write(LuOut,*) any," spin=A, sym=",int_mb(k_irs(1)+any-1)
        enddo
        do any = nfc(ipol)+1,nmo(ipol)-nfv(ipol)
          write(LuOut,*) any," spin=B, sym=",int_mb(k_irs(ipol)+any-1)
        enddo
        call ga_print(g_movecs(1))
        call ga_print(g_movecs(ipol))
      endif
c
c     Hole Alpha
c
      i = 0
      noa = 0
      do sym = 0, nirreps
        j = 0
        do hole = nfc(1)+1, nocc(1)
          if (int_mb(k_irs(1)+hole-1).eq.sym) then          
             i = i + 1
             j = j + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+hole-1)
             call ga_get(g_movecs(1),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        if (j.gt.0) then
           nblocks = j/isize + 1
           do k = 1,nblocks
             noa = noa + 1
             spin_tmp(noa) = 1
             sym_tmp(noa) = sym
             if (k.eq.nblocks) then
               range_tmp(noa) = j-(j/nblocks)*(nblocks-1)
             else
               range_tmp(noa) = j/nblocks
             endif
           enddo
        endif
      enddo
c
c     Hole Beta
c
      nob = 0
      do sym = 0, nirreps
        j = 0
        do hole = nfc(ipol)+1, nocc(ipol)
          if (int_mb(k_irs(ipol)+hole-1).eq.sym) then          
             i = i + 1
             j = j + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+hole-1)
             call ga_get(g_movecs(ipol),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        if (j.gt.0) then
           nblocks = j/isize + 1
           do k = 1,nblocks
             nob = nob + 1
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             if (k.eq.nblocks) then
               range_tmp(noa+nob) = j-(j/nblocks)*(nblocks-1)
             else
               range_tmp(noa+nob) = j/nblocks
             endif
           enddo
        endif
      enddo
c
c     Particle Alpha
c
      nva = 0
      do sym = 0, nirreps
        j = 0
        do particle = nocc(1)+1, nmo(1)-nfv(1)
          if (int_mb(k_irs(1)+particle-1).eq.sym) then          
             i = i + 1
             j = j + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+particle-1)
             call ga_get(g_movecs(1),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        if (j.gt.0) then
           nblocks = j/isize + 1
           do k = 1,nblocks
             nva = nva + 1
             spin_tmp(noa+nob+nva) = 1
             sym_tmp(noa+nob+nva) = sym
             if (k.eq.nblocks) then
               range_tmp(noa+nob+nva) = j-(j/nblocks)*(nblocks-1)
             else
               range_tmp(noa+nob+nva) = j/nblocks
             endif
           enddo
        endif
      enddo
c
c     Particle Beta
c
      nvb = 0
      do sym = 0, nirreps
        j = 0
        do particle = nocc(ipol)+1, nmo(ipol)-nfv(ipol)
          if (int_mb(k_irs(ipol)+particle-1).eq.sym) then          
             i = i + 1
             j = j + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+particle-1)
             call ga_get(g_movecs(ipol),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        if (j.gt.0) then
           nblocks = j/isize + 1
           do k = 1,nblocks
             nvb = nvb + 1
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             if (k.eq.nblocks) then
               range_tmp(noa+nob+nva+nvb) = j-(j/nblocks)*(nblocks-1)
             else
               range_tmp(noa+nob+nva+nvb) = j/nblocks
             endif
           enddo
        endif
      enddo
c
c     Holes and particles
c
      noab = noa + nob
      nvab = nva + nvb
c
c     Allocate tiling related arrays
c
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Spin of blocks",
     1  l_spin,k_spin)) call errquit("tce_tile: MA problem",0)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Symmetry of blocks",
     1  l_sym,k_sym)) call errquit("tce_tile: MA problem",1)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Range of blocks",
     1  l_range,k_range)) call errquit("tce_tile: MA problem",2)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Offset of blocks",
     1  l_offset,k_offset)) call errquit("tce_tile: MA problem",3)
      j = 0
      do i = 1,noa+nob+nva+nvb
        int_mb(k_spin+i-1)=spin_tmp(i)
        int_mb(k_sym+i-1)=sym_tmp(i)
        int_mb(k_range+i-1)=range_tmp(i)
        int_mb(k_offset+i-1)=j
        j = j + range_tmp(i)
      enddo
c
      if (nodezero.and.util_print('tile',print_debug)) then
        write(LuOut,*) "Sorted"
        do any = 1,nmo(1)-nfv(1)-nfc(1)+nmo(ipol)-nfv(ipol)-nfc(ipol)
          write(LuOut,*) any," spin=",int_mb(k_spin_sorted+any-1)
     1                      ," irep=",int_mb(k_irs_sorted+any-1)
     2                      ," eval=",dbl_mb(k_evl_sorted+any-1)
        enddo
        call ma_print(dbl_mb(k_movecs_sorted),nbf,
     1  nmo(1)-nfv(1)-nfc(1)+nmo(ipol)-nfv(ipol)-nfc(ipol),"MO coeffs")
      endif
      if (nodezero.and.util_print('tile',print_default)) then
        write(LuOut,*)
        write(LuOut,9050)
        write(LuOut,*)
        write(LuOut,9070)
        do i = 1,noa+nob+nva+nvb
          call sym_irrepname(geom,int_mb(k_sym+i-1)+1,irrepname)
          write(LuOut,9060) i,spinname(int_mb(k_spin+i-1)),irrepname,
     1      int_mb(k_range+i-1),int_mb(k_offset+i-1)
        enddo
      endif
 9050 format(1x,'Suggested orbital blocking')
 9070 format(1x,'Block   Spin    Irrep     Size     Offset',/,
     1       1x,'-----------------------------------------')
 9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8)
c
c     ======
c     Return
c     ======
c
      do i = 1, ipol
        if (.not.ga_destroy(g_movecs(i))) 
     1    call errquit('tce_tile: GA problem',0)
      enddo
      return
      end             
