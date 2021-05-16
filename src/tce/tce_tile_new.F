      subroutine tce_tile_new(rtdb)
c
c $Id$
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
#include "tce_main.fh"
      integer rtdb
      double precision ga_dble
      double precision ma_dble
      integer isize,jsize
      double precision size
      integer any
      integer i,j
      logical nodezero
      integer spin_tmp(max_size)
      integer sym_tmp(max_size)
      integer range_tmp(max_size)
      integer active_tiles(max_size)
      integer b2am(max_size)
      integer spin_tmp_alpha(max_size)
      integer sym_tmp_alpha(max_size)
      integer range_tmp_alpha(max_size)
      character*4 irrepname
      character*5 spinname(2)
      data spinname/'alpha','beta '/
      integer ix1,ix2
      double precision xxx
      double precision vcut
      double precision bignumx
      integer max_tile
c
c     ====================
c     Retrieve Memory Info
c     ====================
c
      nodezero = (ga_nodeid().eq.0)
      ga_dble = dble(ga_memory_avail()) / dble(bytes)
      ma_dble = dble(ma_inquire_avail(mt_dbl))
      call ga_dgop(mt_dbl,ga_dble,1,'+')
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
        endif
      endif
c
c     ========================
c     Determine the block size
c     ========================
c
      size = dble(available)
      size = size / dble(nbf*(nbf+1)/2)
      size = size / 2.0d0
      size = dsqrt(size)
      isize = int(size)
c
      size = dble(available)
      size = size / 6.0d0
      if (model.eq."ccsd") then
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."ccsdt") then
        size = size ** (1.0d0/6.0d0)
        isize = min(isize,int(size))
      else if (model.eq."ccsdtq") then
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else
        call errquit("tce_tile_new: the model NYI",0,CAPMIS_ERR)
      endif
      if (nodezero) write(LuOut,*)
      if (rtdb_get(rtdb,'tce:tilesize',mt_int,1,jsize)) then
        tile_dim = jsize
        isize = jsize
        if (nodezero) write(LuOut,9050)
      endif
      if (nodezero.and.util_print('tile',print_default)) then
        write(LuOut,9040) isize
      endif
c
c     ================================================
c     Reorder the orbitals in the ha, hb, pa, pb order
c     ================================================
c
       bignumx = 10.0d+10
       if (.not.rtdb_get(rtdb,'tce:vcut',mt_dbl,1,vcut)) then
           vcut = bignumx
       end if
       if(vcut.lt.bignumx) then !vcut starts here
         ix1=0
         ix2=0
         do i=nocc(1)+1,nmo(1)
           xxx=dbl_mb(k_evl(1)+i-1)
           if(xxx.gt.vcut) then
               ix1=ix1+1
           end if
         enddo
c
         do i=nocc(ipol)+1,nmo(ipol)
           xxx=dbl_mb(k_evl(ipol)+i-1)
           if(xxx.gt.vcut) then
             ix2=ix2+1
           end if
         enddo
c
         if(ix1.eq.ix2) then
           nfv(1)=ix1
           nfv(ipol)=ix2
         else
           call errquit("tce_tile_new: dyn. frozen virtuals",1,MA_ERR)
       end if
c
       if(nodezero) then
        write(6,*)'redefined number of frozen virt. =',nfv(1)
        write(6,8778) dbl_mb(k_evl(1)+nmo(1)-1)
        call util_flush(6)
       end if
c link to QM/MM
        if (.not.rtdb_put(rtdb,'tce:eorbmax',mt_dbl,1,
     1           dbl_mb(k_evl(1)+nmo(1)-1)))
     2  call errquit('tce_energy: RTDB problem',0,MA_ERR)
c
       end if !vcut ends up here
c
      if (.not.ma_push_get(mt_dbl,nbf*(nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol)),"sorted MO coeffs",
     2  l_movecs_sorted,k_movecs_sorted))
     3  call errquit("tce_tile_new: MA problem",0,MA_ERR)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted spins",
     2  l_spin_sorted,k_spin_sorted))
     3  call errquit("tce_tile_new: MA problem",1,MA_ERR)
      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted irs",
     2  l_irs_sorted,k_irs_sorted))
     3  call errquit("tce_tile_new: MA problem",2,MA_ERR)
      if (.not.ma_push_get(mt_dbl,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted evl",
     2  l_evl_sorted,k_evl_sorted))
     3  call errquit("tce_tile_new: MA problem",3,MA_ERR)
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
c->d3p975
c
      if(.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1 +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted MO index",
     2 l_mo_index,k_mo_index))
     3 call errquit("tce_tile_new: MA problem",4,ma_err)
c
c<-d3p975
c
      if(intorb) then
        call tce_tile_new_intorb(isize,
     &           spin_tmp,sym_tmp,range_tmp,
     &           active_tiles,b2am,
     &           spin_tmp_alpha,sym_tmp_alpha,range_tmp_alpha)
      else
        call tce_tile_new_spinorb(isize,
     &           spin_tmp,sym_tmp,range_tmp,
     &           active_tiles)
      end if !intorb
c
c     Allocate tiling related arrays
c
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Spin of blocks",
     1  l_spin,k_spin)) call errquit("tce_tile_new: MA problem",4,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Symmetry of blocks",
     1  l_sym,k_sym)) call errquit("tce_tile_new: MA problem",5,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Range of blocks",
     1  l_range,k_range)) call errquit("tce_tile_new: MA problem",6,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Offset of blocks",
     1  l_offset,k_offset)) call errquit("tce_tile_new: MA problem",7,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"alpha-izer",
     1  l_alpha,k_alpha)) call errquit("tce_tile_new: MA problem",8,
     2  MA_ERR)
      if(activecalc) THEN
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Active blocks ",
     1  l_active,k_active)) call errquit("tce_tile_new: MA problem",9,
     2  MA_ERR)
      end if
      if(intorb) THEN
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"alpha-izer 2 ",
     1  l_b2am,k_b2am)) call errquit("tce_tile_new: MA problem",9,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha spins ",
     1  l_spin_alpha,k_spin_alpha))
     2  call errquit("tce_tile_new: MA problem",9,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha syms ",
     1  l_sym_alpha,k_sym_alpha))
     2  call errquit("tce_tile_new: MA problem",9,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha range ",
     1  l_range_alpha,k_range_alpha))
     2  call errquit("tce_tile_new: MA problem",9,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha range ",
     1  l_offset_alpha,k_offset_alpha))
     2  call errquit("tce_tile_new: MA problem",9,MA_ERR)
      end if
c
      if(.not.intorb) THEN
c
      j = 0
      do i = 1,noa+nob+nva+nvb
        int_mb(k_spin+i-1)=spin_tmp(i)
        int_mb(k_sym+i-1)=sym_tmp(i)
        int_mb(k_range+i-1)=range_tmp(i)
        int_mb(k_active+i-1)=active_tiles(i)
        int_mb(k_offset+i-1)=j
        j = j + range_tmp(i)
      enddo
c max_tile_size => tile_dim --------
        max_tile=0
        do i = 1,noa+nob+nva+nvb
          if(int_mb(k_range+i-1).gt.max_tile)
     &     max_tile=int_mb(k_range+i-1)
        enddo
        tile_dim=max_tile
c ----------------------------------
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
c
      end if !not intorb
c
      if(intorb) THEN
c
      j = 0
      do i = 1,noa+nob+nva+nvb
        int_mb(k_spin+i-1)=spin_tmp(i)
        int_mb(k_sym+i-1)=sym_tmp(i)
        int_mb(k_range+i-1)=range_tmp(i)
c        int_mb(k_active+i-1)=active_tiles(i) !l_active will not be used
        int_mb(k_offset+i-1)=j
        j = j + range_tmp(i)
      enddo
c max_tile_size => tile_dim --------
        max_tile=0
        do i = 1,noa+nob+nva+nvb
          if(int_mb(k_range+i-1).gt.max_tile)
     &     max_tile=int_mb(k_range+i-1)
        enddo
        tile_dim=max_tile
c ----------------------------------
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
        do i = 1,noa+nob+nva+nvb
          int_mb(k_b2am+i-1) = b2am(i)
        enddo
c
      j = 0
      do i = 1,noa+nva
        int_mb(k_spin_alpha+i-1)=spin_tmp_alpha(i)
        int_mb(k_sym_alpha+i-1)=sym_tmp_alpha(i)
        int_mb(k_range_alpha+i-1)=range_tmp_alpha(i)
c        int_mb(k_active+i-1)=active_tiles(i) !l_active will not be used
        int_mb(k_offset_alpha+i-1)=j
        j = j + range_tmp_alpha(i)
      enddo
c
      end if  ! intorb
c
c
      if(nodezero) then
       write(LuOut,9080) tile_dim
      end if
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
 8778  format('Orbital energy of the highest virtual orbital',2x,f17.8)
 9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
 9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
     1       10x,'Available MA space size is    ',i10,' doubles')
 9040 format(1x,'Maximum block size',i10,' doubles')
 9050 format(1x,'Maximum block size supplied by input')
 9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha',/,
     1         1x,'-------------------------------------------------')
 9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8)
 9080 format(/,1x,'tile_dim = ',2x,i4)
      return
      end
c
      subroutine tce_tile_new_spinorb(isize,
     &           spin_tmp,sym_tmp,range_tmp,
     &           active_tiles)
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer isize
      integer sym
      integer hole, particle
      integer i,j,k,nblocks
      integer spin_tmp(max_size)
      integer sym_tmp(max_size)
      integer range_tmp(max_size)
      integer active_tiles(max_size)
      integer j_a,j_ia
      integer l_a,l_ia
      integer nblocks_a,nblocks_ia
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*hole-1
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*hole-1
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*hole
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*hole
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*particle-1
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*particle-1
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*particle
c<-d3p975
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
c->d3p975
             int_mb(k_mo_index+i-1)=2*particle
c<-d3p975
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
c
c8778  format('Orbital energy of the highest virtual orbital',2x,f17.8)
c9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
c9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
c    1       10x,'Available MA space size is    ',i10,' doubles')
c9040 format(1x,'Maximum block size',i10,' doubles')
c9050 format(1x,'Maximum block size supplied by input')
c9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha',/,
c    1         1x,'-------------------------------------------------')
c9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8)
c9080 format(/,1x,'tile_dim = ',2x,i4)
      return
      end
c
      subroutine tce_tile_new_intorb(isize,
     &           spin_tmp,sym_tmp,range_tmp,
     &           active_tiles,b2am,
     &           spin_tmp_alpha,sym_tmp_alpha,range_tmp_alpha)
c
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer isize
      integer sym
      integer hole, particle
      integer i,j,k,nblocks
      integer spin_tmp(max_size)
      integer sym_tmp(max_size)
      integer range_tmp(max_size)
      integer active_tiles(max_size)
      integer b2am(max_size)
      integer spin_tmp_alpha(max_size)
      integer sym_tmp_alpha(max_size)
      integer range_tmp_alpha(max_size)
      integer j_a,j_ia
      integer l_a,l_ia
      integer nblocks_a,nblocks_ia
c
c     Hole Alpha
c
      i = 0
      noa = 0
      do sym = 0, nirreps
        j_ia = 0           !occ. inactive index
        j_a  = 0           !occ.   active index
        do hole = nfc(1)+1, nocc(2)
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
        do hole =nocc(2)+1,nocc(1) !FOR RHF this part does not
contribute
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
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           if (j_a  .gt. isize*nblocks_a)  nblocks_a  = nblocks_a+1
           l_ia = 0
           do k = 1,nblocks_ia
             noa = noa + 1
             active_tiles(noa) = 0
             spin_tmp(noa) = 1
              spin_tmp_alpha(noa) = 1
             sym_tmp(noa) = sym
              sym_tmp_alpha(noa) = sym
             range_tmp(noa) = k*j_ia/nblocks_ia-l_ia
              range_tmp_alpha(noa) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa)
           enddo
           l_a = 0
           do k = 1,nblocks_a
             noa = noa + 1
             active_tiles(noa) = 1
             spin_tmp(noa) = 1
              spin_tmp_alpha(noa) = 1
             sym_tmp(noa) = sym
              sym_tmp_alpha(noa) = sym
             range_tmp(noa) = k*j_a/nblocks_a-l_a
              range_tmp_alpha(noa) = k*j_a/nblocks_a-l_a
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
        do hole = nfc(ipol)+1, nocc(ipol) !beta electrons here no active part
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
        j=j_ia
        if (j.gt.0) then
           nblocks_ia = j_ia/isize
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           l_ia = 0
           do k = 1,nblocks_ia
             nob = nob + 1
             active_tiles(noa+nob) = 0
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             range_tmp(noa+nob) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob)
           enddo
        endif
      enddo
      oatiles(ipol)=0
c
c
c     Particle Alpha
c
c
      nva = 0
      do sym = 0, nirreps
        j_ia = 0
        do particle = nocc(1)+1, nmo(1)-nfv(1)
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
        j=j_ia
        if (j.gt.0) then
           nblocks_ia  =  j_ia/isize
           if (j_ia  .gt. isize*nblocks_ia)  nblocks_ia  = nblocks_ia+1
           nblocks = nblocks_ia
           l_ia = 0
           do k = 1,nblocks_ia
             nva = nva + 1
             active_tiles(noa+nob+nva) = 0
             spin_tmp(noa+nob+nva) = 1
              spin_tmp_alpha(noa+nva) = 1
             sym_tmp(noa+nob+nva) = sym
              sym_tmp_alpha(noa+nva) = sym
             range_tmp(noa+nob+nva) = k*j_ia/nblocks_ia-l_ia
              range_tmp_alpha(noa+nva) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob+nva)
           enddo
        endif
      enddo
      vatiles(1)=0
c
c
c     Particle Beta
c
      nvb = 0
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        do particle = nocc(ipol)+1, nocc(1) !for RHF this part does not
contribute
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
        do particle = nocc(1)+1, nmo(ipol)-nfv(ipol)
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
c
c forming b2am matrix
c
      do k=1,max_size
       b2am(k)=0
      enddo
c hole alpha
      do k=1,noa
       b2am(k)=k
      enddo
c hole beta
      j=1
      do k=1,noa
       if(active_tiles(k).eq.0) then
        b2am(noa+j)=k
        j=j+1
       end if
      enddo
c particle alpha
      do k=1,nva
       b2am(noab+k)=noa+k
      enddo
c particle beta
c  active part from hole alpha
      do k=1,nvb
       if(active_tiles(noab+nva+k).eq.1) then
         do j=1,noa
          if(active_tiles(j).eq.1) then
           b2am(noab+nva+k)=b2am(j)
ccc           b2am(j)=-1
           active_tiles(j)=-1
           go to 2900
          end if
         enddo
 2900    continue
       end if
      enddo
      do j=1,noa
       if(active_tiles(j).eq.-1) active_tiles(j)=1
      enddo
c inactive part from particle alpha
      do k=1,nvb
       if(active_tiles(noab+nva+k).eq.0) then
         do j=1,nva
          if(active_tiles(noab+j).eq.0) then
           b2am(noab+nva+k)=b2am(noab+j)
ccc           b2am(noab+j)=-1
           active_tiles(noab+j)=-1
           go to 2901
          end if
         enddo
 2901    continue
       end if
      enddo
      do j=1,nva
       if(active_tiles(noab+j).eq.-1) active_tiles(noab+j)=0
      enddo
c done with b2am matrix
c8778  format('Orbital energy of the highest virtual orbital',2x,f17.8)
c9000 format(10x,'Available GA+MA space size is ',i10,' doubles')
c9010 format(10x,'Available GA space size is    ',i10,' doubles',/,
c    1       10x,'Available MA space size is    ',i10,' doubles')
c9040 format(1x,'Maximum block size',i10,' doubles')
c9050 format(1x,'Maximum block size supplied by input')
c9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha',/,
c    1         1x,'-------------------------------------------------')
c9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8)
c9080 format(/,1x,'tile_dim = ',2x,i4)
      return
      end
