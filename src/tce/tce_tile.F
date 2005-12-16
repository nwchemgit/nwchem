      subroutine tce_tile(rtdb)
c
c $Id: tce_tile.F,v 1.16 2005-12-16 22:16:05 kowalski Exp $
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
#include "errquit.fh"
#include "tce.fh"
c
      integer rtdb
      double precision ga_dble
      double precision ma_dble
      integer isize,jsize
      double precision size
      integer sym
      integer any, hole, particle
      integer i,j,k,l,nblocks
      logical nodezero
      integer spin_tmp(max_size)
      integer sym_tmp(max_size)
      integer range_tmp(max_size)
      integer active_tiles(max_size)
c CCSDtq...
      integer j_a,j_ia
      integer l_a,l_ia
      integer nblocks_a,nblocks_ia 
c CCSDtq...
      character*4 irrepname
      character*5 spinname(2)
      data spinname/'alpha','beta '/
c
c     For orbital rotation test
c
c     integer l_unitary, k_unitary
c     integer l_movecs_rotated, k_movecs_rotated
c     integer n, orb1, orb2
c     double precision angle
c
c     ====================
c     Retrieve Memory Info
c     ====================
c
      nodezero = (ga_nodeid().eq.0)
      ga_dble = dfloat(ga_memory_avail()) / dfloat(bytes)
      ma_dble = dfloat(ma_inquire_avail(mt_dbl))
      call ga_dgop(msg_tcemem,ga_dble,1,'+')
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
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."lccd") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."ccsd") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."lccsd") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."ccsdta") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/6.0d0)
        isize = min(isize,int(size))
      else if (model.eq."ccsdt") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/6.0d0)
        isize = min(isize,int(size))
      else if (model.eq."ccsdtq") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else if (model.eq."qcisd") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."cisd") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."cisdt") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/6.0d0)
        isize = min(isize,int(size))
      else if (model.eq."cisdtq") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else if (model.eq."mbpt2") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."mbpt3") then
        size = dfloat(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."mbpt4") then
        size = dfloat(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else
        call errquit("tce_tile: the model NYI",0,CAPMIS_ERR)
      endif
      if (nodezero) write(LuOut,*)
      if (rtdb_get(rtdb,'tce:tilesize',mt_int,1,jsize)) then
        isize = jsize
        if (nodezero) write(LuOut,9050)
      endif
      if (nodezero.and.util_print('tile',print_default)) then
        write(LuOut,9040) isize
      endif
 9040 format(1x,'Maximum block size',i10,' doubles')
 9050 format(1x,'Maximum block size supplied by input')
c
c     ================================================
c     Reorder the orbitals in the ha, hb, pa, pb order
c     ================================================
c
      if (.not.ma_push_get(mt_dbl,nbf*(nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol)),"sorted MO coeffs",
     2  l_movecs_sorted,k_movecs_sorted))
     3  call errquit("tce_tile: MA problem",0,MA_ERR)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted spins",
     2  l_spin_sorted,k_spin_sorted))
     3  call errquit("tce_tile: MA problem",1,MA_ERR)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted irs",
     2  l_irs_sorted,k_irs_sorted))
     3  call errquit("tce_tile: MA problem",2,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted evl",
     2  l_evl_sorted,k_evl_sorted))
     3  call errquit("tce_tile: MA problem",3,MA_ERR)
      nirreps = sym_number_ops(geom) + 1
      if (nodezero.and.util_print('tile',print_debug)) then
        write(LuOut,*) "Number of irreps = ",nirreps
        do any = nfc(1)+1,nmo(1)-nfv(1)
          write(LuOut,*) any," spin=A, sym=",int_mb(k_irs(1)+any-1)
        enddo
        do any = nfc(ipol)+1,nmo(ipol)-nfv(ipol)
          write(LuOut,*) any," spin=B, sym=",int_mb(k_irs(ipol)+any-1)
        enddo
      endif
c
c     Hole Alpha
c
      i = 0
      noa = 0
      do sym = 0, nirreps
        j_ia = 0           !occ. inactive index
        j_a  = 0           !occ.   active index
        do hole = nfc(1)+1, nocc(1)-oactive(1)
          if (int_mb(k_irs(1)+hole-1).eq.sym) then          
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+hole-1)
             call ga_get(g_movecs(1),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        do hole =nocc(1)-oactive(1)+1,nocc(1)
          if (int_mb(k_irs(1)+hole-1).eq.sym) then
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+hole-1)
             call ga_get(g_movecs(1),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        j=j_ia+j_a
        if (j.gt.0) then
           nblocks_ia = j_ia/isize
           nblocks_a  =  j_a/isize
c           nblocks = nblocks_ia+nblocks_a
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           if (j_a  .gt. isize*nblocks_a)  nblocks_a  = nblocks_a+1
           l_ia = 0
           do k = 1,nblocks_ia
             noa = noa + 1
             active_tiles(noa) = 0
             spin_tmp(noa) = 1
             sym_tmp(noa) = sym
             range_tmp(noa) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa)
           enddo
           l_a = 0
           do k = 1,nblocks_a
             noa = noa + 1
             active_tiles(noa) = 1
             spin_tmp(noa) = 1
             sym_tmp(noa) = sym
             range_tmp(noa) = k*j_a/nblocks_a-l_a
             l_a = l_a + range_tmp(noa)
           enddo
        endif
      enddo
      oatiles(1)=nblocks_a
c
c     Hole Beta
c
      nob = 0
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        do hole = nfc(ipol)+1, nocc(ipol)-oactive(ipol)
          if (int_mb(k_irs(ipol)+hole-1).eq.sym) then          
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+hole-1)
             call ga_get(g_movecs(ipol),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        do hole=nocc(ipol)-oactive(ipol)+1,nocc(ipol)
          if (int_mb(k_irs(ipol)+hole-1).eq.sym) then
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+hole-1)
             call ga_get(g_movecs(ipol),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        j=j_ia+j_a
        if (j.gt.0) then
           nblocks_ia = j_ia/isize
           nblocks_a  =  j_a/isize
c           nblocks = nblocks_ia+nblocks_a
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           if (j_a  .gt. isize*nblocks_a)  nblocks_a  = nblocks_a+1
           l_ia = 0
           do k = 1,nblocks_ia
             nob = nob + 1
             active_tiles(noa+nob) = 0
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             range_tmp(noa+nob) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob)
           enddo
           l_a=0
           do k = 1,nblocks_a
             nob = nob + 1
             active_tiles(noa+nob) = 1
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             range_tmp(noa+nob) = k*j_a/nblocks_a-l_a
             l_a = l_a + range_tmp(noa+nob)
c *** debug ***
           write(6,*)'range active',k,nblocks_a,range_tmp(noa+nob)
           call util_flush(LuOut)  
c *************
           enddo
        endif
      enddo
      oatiles(ipol)=nblocks_a
c
c     Particle Alpha
c
      nva = 0
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        do particle = nocc(1)+1, nocc(1)+vactive(1)
          if (int_mb(k_irs(1)+particle-1).eq.sym) then          
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+particle-1)
             call ga_get(g_movecs(1),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        do particle = nocc(1)+vactive(1)+1, nmo(1)-nfv(1)
          if (int_mb(k_irs(1)+particle-1).eq.sym) then
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+particle-1)
             call ga_get(g_movecs(1),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        j=j_a+j_ia
        if (j.gt.0) then
           nblocks_a   = j_a/isize
           nblocks_ia  =  j_ia/isize
c           nblocks = nblocks_a+nblocks_ia
           if (j_a .gt. isize*nblocks_a) nblocks_a = nblocks_a+1
           if (j_ia  .gt. isize*nblocks_ia)  nblocks_ia  = nblocks_ia+1
           nblocks = nblocks_a+nblocks_ia
           l_a = 0
           do k = 1,nblocks_a
             nva = nva + 1
             active_tiles(noa+nob+nva) = 1
             spin_tmp(noa+nob+nva) = 1
             sym_tmp(noa+nob+nva) = sym
             range_tmp(noa+nob+nva) = k*j_a/nblocks_a-l_a
             l_a = l_a + range_tmp(noa+nob+nva)
           enddo
           l_ia = 0
           do k = 1,nblocks_ia
             nva = nva + 1
             active_tiles(noa+nob+nva) = 0 
             spin_tmp(noa+nob+nva) = 1
             sym_tmp(noa+nob+nva) = sym
             range_tmp(noa+nob+nva) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob+nva)
           enddo
        endif
      enddo
      vatiles(1)=nblocks_a
c
c     Particle Beta
c
      nvb = 0
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        do particle = nocc(ipol)+1, nocc(ipol)+vactive(ipol)
          if (int_mb(k_irs(ipol)+particle-1).eq.sym) then          
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+particle-1)
             call ga_get(g_movecs(ipol),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        do particle = nocc(ipol)+vactive(ipol)+1, nmo(ipol)-nfv(ipol)
          if (int_mb(k_irs(ipol)+particle-1).eq.sym) then
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+particle-1)
             call ga_get(g_movecs(ipol),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
          endif
        enddo
        j=j_a+j_ia
        if (j.gt.0) then
           nblocks_a   = j_a/isize
           nblocks_ia  =  j_ia/isize
c           nblocks = nblocks_a+nblocks_ia  
           if (j_a .gt. isize*nblocks_a) nblocks_a = nblocks_a+1
           if (j_ia  .gt. isize*nblocks_ia)  nblocks_ia  = nblocks_ia+1
           nblocks = nblocks_a+nblocks_ia 
           l_a = 0
           do k = 1,nblocks_a
             nvb = nvb + 1
             active_tiles(noa+nob+nva+nvb) = 1
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             range_tmp(noa+nob+nva+nvb) = k*j_a/nblocks_a-l_a
             l_a = l_a + range_tmp(noa+nob+nva+nvb)
           enddo
           l_ia = 0
           do k = 1,nblocks_ia
             nvb = nvb + 1
             active_tiles(noa+nob+nva+nvb) = 0
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             range_tmp(noa+nob+nva+nvb) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob+nva+nvb)
           enddo
        endif
      enddo
      vatiles(ipol)=nblocks_a
c
c     Holes and particles
c
      noab = noa + nob
      nvab = nva + nvb
c *** debug ***
c      write(6,*)'testkk'
c      write(6,*)'oactive(1)=',oactive(1)
c      write(6,*)'oactive(2)=',oactive(2)
c      write(6,*)'vactive(1)=',vactive(1)
c      write(6,*)'vactive(2)=',vactive(2)
c      write(6,*)'ipol=',ipol
c      write(6,*)'nvb=',nvb
c      write(6,*)'j_a=',j_a
c      write(6,*)'j_ia=',j_ia
c      write(6,*)'vatiles(1)=',vatiles(1)
c      write(6,*)'vatiles(2)=',vatiles(2)
c      write(6,*)'ACTIVE TILES'
c       do i=1,noab+nvab
c        write(6,*) active_tiles(i)
c       enddo
c      write(6,*)'        '
c      call util_flush(LuOut)
c *************
c
c     Allocate tiling related arrays
c
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Spin of blocks",
     1  l_spin,k_spin)) call errquit("tce_tile: MA problem",4,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Symmetry of blocks",
     1  l_sym,k_sym)) call errquit("tce_tile: MA problem",5,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Range of blocks",
     1  l_range,k_range)) call errquit("tce_tile: MA problem",6,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Offset of blocks",
     1  l_offset,k_offset)) call errquit("tce_tile: MA problem",7,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"alpha-izer",
     1  l_alpha,k_alpha)) call errquit("tce_tile: MA problem",8,
     2  MA_ERR)
      if(activecalc) THEN
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Active blocks ",
     1  l_active,k_active)) call errquit("tce_tile: MA problem",9,
     2  MA_ERR)
      end if
      j = 0
      do i = 1,noa+nob+nva+nvb
        int_mb(k_spin+i-1)=spin_tmp(i)
        int_mb(k_sym+i-1)=sym_tmp(i)
        int_mb(k_range+i-1)=range_tmp(i)
        int_mb(k_active+i-1)=active_tiles(i)
        int_mb(k_offset+i-1)=j
        j = j + range_tmp(i)
      enddo
      if (restricted) then
        do i = 1,noa
          int_mb(k_alpha+i-1) = i
        enddo
        do i = noa+1,noa+nob
          int_mb(k_alpha+i-1) = i - noa
        enddo
        do i = noa+nob+1,noa+nob+nva
          int_mb(k_alpha+i-1) = i
        enddo
        do i = noa+nob+nva+1,noa+nob+nva+nvb
          int_mb(k_alpha+i-1) = i - nva
        enddo
      else
        do i = 1,noa+nob+nva+nvb
          int_mb(k_alpha+i-1) = i
        enddo
      endif
c
      if (nodezero.and.util_print('tile',print_debug)) then
        write(LuOut,*) "Sorted"
        do any = 1,nmo(1)-nfv(1)-nfc(1)+nmo(ipol)-nfv(ipol)-nfc(ipol)
          write(LuOut,*) any," spin=",int_mb(k_spin_sorted+any-1),
     1                       " irep=",int_mb(k_irs_sorted+any-1),
     2                       " eval=",dbl_mb(k_evl_sorted+any-1)
        enddo
        call ma_print(dbl_mb(k_movecs_sorted),nbf,
     1  nmo(1)-nfv(1)-nfc(1)+nmo(ipol)-nfv(ipol)-nfc(ipol),"MO coeffs")
      endif
      if (nodezero.and.util_print('tile',print_default)) then
        write(LuOut,9070)
        do i = 1,noa+nob+nva+nvb
          call sym_irrepname(geom,int_mb(k_sym+i-1)+1,irrepname)
          write(LuOut,9060) i,spinname(int_mb(k_spin+i-1)),irrepname,
     1      int_mb(k_range+i-1),int_mb(k_offset+i-1),int_mb(k_alpha+i-1)
        enddo
      endif
 9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha',/,
     1         1x,'-------------------------------------------------')
 9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8)
c
c     =============================================
c     Test for a unitary transformation of orbitals
c     =============================================
c
c     n = nmo(1)-nfv(1)-nfc(1)+nmo(ipol)-nfv(ipol)-nfc(ipol)
c     if (nodezero) then
c       write(LuOut,*)
c       write(LuOut,*) " *************************************"
c       write(LuOut,*) " !!! CAUTION !!! Orbitals rotation !!!"
c       write(LuOut,*) " !!! CAUTION !!! Turn off symmetry !!!"
c       write(LuOut,*) " !!! CAUTION !!! Spin unrestricted !!!"
c       write(LuOut,*) " *************************************"
c       write(LuOut,*) " Jacobi can diverge even when the "
c       write(LuOut,*) " theory is intrinsically invariant"
c     endif
c     if (.not.ma_push_get(mt_dbl,n*n,"Unitary",
c    1  l_unitary,k_unitary)) call errquit("tce_tile: MA problem",-1,
c    2  MA_ERR)
c     if (.not.ma_push_get(mt_dbl,nbf*n,"rotated MO coeffs",
c    2  l_movecs_rotated,k_movecs_rotated))
c    3  call errquit("tce_tile: MA problem",-2,MA_ERR)
c     call ma_print(dbl_mb(k_movecs_sorted),nbf,n,
c    1  "Original MO coeffs")
c
c     Change from here ---------------------------------
c
c     orb1 = 4
c     orb2 = 5
c     angle = datan(1.0d0)*1.0d0
c
c     ------------------------------------------ to here
c
c     do i = 1,n
c       do j = 1,n
c         dbl_mb(k_unitary+(i-1)*(n)+j-1) = 0.0d0
c       enddo
c       dbl_mb(k_unitary+(i-1)*(n)+i-1) = 1.0d0
c     enddo
c     dbl_mb(k_unitary+(orb1-1)*n+orb1-1) = dcos(angle)
c     dbl_mb(k_unitary+(orb1-1)*n+orb2-1) = dsin(angle)
c     dbl_mb(k_unitary+(orb2-1)*n+orb1-1) = - dsin(angle)
c     dbl_mb(k_unitary+(orb2-1)*n+orb2-1) = dcos(angle)
c     call dgemm('N','N',nbf,n,n,1.0d0,
c    1  dbl_mb(k_movecs_sorted),nbf,dbl_mb(k_unitary),n,
c    2  0.0d0,dbl_mb(k_movecs_rotated),nbf)
c     do i = 1,nbf
c       do j = 1,n
c         dbl_mb(k_movecs_sorted+(i-1)*n+j-1) = 
c    1    dbl_mb(k_movecs_rotated+(i-1)*n+j-1)
c       enddo
c     enddo
c     write(LuOut,*) " Orbitals ",orb1," and ",orb2," are rotated by ",
c    1  angle
c     call ma_print(dbl_mb(k_unitary),n,n,
c    1  "Unitary rotation")
c     call ma_print(dbl_mb(k_movecs_sorted),nbf,n,
c    1  "Rotated MO coeffs")
c     if (.not.ma_pop_stack(l_movecs_rotated))
c    1  call errquit("tce_tile: MA problem",-3,MA_ERR)
c     if (.not.ma_pop_stack(l_unitary))
c    1  call errquit("tce_tile: MA problem",-4,MA_ERR)
c
c     ========
c     Test end
c     ========
c
      return
      end
