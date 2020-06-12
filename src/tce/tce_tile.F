      subroutine tce_tile(rtdb)
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
c
c
c THINGS TO DO 
c
c activate active_tiles, oatiles, and vatiles for intorb=.true.
c
c define b2a(max_size) matrix in tce_tile
c 
c declare in tce.fh:
c             l_b2am,k_b2am
c             l_spin_alpha,k_spin_alpha
c             l_sym_alpha,k_sym_alpha
c             l_range_alpha,k_range_alpha
c
c allocate: l_b2am,l_spin_alpha,l_sym_alpha,l_range_alpha
c de-allocate the above matrices in tce_energy.F
c
c
c
c
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
c ccsd_act/eomccsd_act
      integer active_tiles_n(max_size)
      integer j_an,nblocks_an,l_an
c INTORB
      integer b2am(max_size)
      integer spin_tmp_alpha(max_size)
      integer sym_tmp_alpha(max_size)
      integer range_tmp_alpha(max_size)
c CCSDtq...
      integer j_a,j_ia
      integer l_a,l_ia
      integer nblocks_a,nblocks_ia 
c CCSDtq...
      character*4 irrepname
      character*5 spinname(2)
      data spinname/'alpha','beta '/
c dynamically frozen virtuals
      integer ix1,ix2
      double precision xxx
      double precision vcut
      double precision bignumx
c for tile_dim
      integer max_tile
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
      ga_dble = dble(ga_memory_avail()) / dble(bytes)
      ma_dble = dble(ma_inquire_avail(mt_dbl))
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
c     previous estimation scheme led to huge tilesies in some cases
c     and users encountered segfaults
c     these estimates are what Jeff and Karol recommend anyways
c
      isize = 32
      if ((model.eq."ccd").or.(model.eq."lccd").or.(model.eq."cis")
     1.or.(model.eq."ccsd").or.(model.eq."lccsd")
     2.or.(model.eq."qcisd").or.(model.eq."cisd")
     4.or.(model.eq."mbpt2").or.(model.eq."mbpt3")
     5.or.(model.eq."ccsd_act")) then
        isize = 32
        if (rtdb_cget(rtdb,'tce:perturbative',1,perturbative)) then
          if (perturbative.eq.'(t)') isize = 24
          if (perturbative.eq.'[t]') isize = 24
          if (perturbative.eq.'lambda(t)') isize = 24
          if (perturbative.eq.'cr_(t)') isize = 24
          if (perturbative.eq.'cr_[t]') isize = 24
          if (perturbative.eq.'creom_(t)') isize = 20
          if (perturbative.eq.'emb1') isize = 24
          if (perturbative.eq.'emb2') isize = 24
          if (perturbative.eq.'2_t') isize = 24
          if (perturbative.eq.'2_tq') isize = 16
          if (perturbative.eq.'2_q') isize = 16
          if (perturbative.eq.'lr_(t)') isize = 24
          if (perturbative.eq.'lr_(tq1)') isize = 16
          if (perturbative.eq.'lr_(tq1p)') isize = 16
        endif
      else if (model.eq."ccsdta") then
        isize = 24
      elseif ((model.eq."ccsdt").or.(model.eq."cisdt")) then
        isize = 20
      elseif ((model.eq."ccsdtq").or.(model.eq."cisdtq")
     1   .or.(model.eq."mbpt4")) then
        isize = 12
c       special case of mbpt4sdq(t) can use CCSD(T)-like tilesizes
        if (rtdb_cget(rtdb,'tce:perturbative',1,perturbative)) then
          if (perturbative.eq.'(t)') isize = 24
        endif
      else
        isize = 32
        if (nodezero) then
          write(LuOut,*) 'Warning: unknown model (tce_tile.F)'
        endif
c        call errquit("tce_tile: the model NYI",0,CAPMIS_ERR)
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
 9040 format(1x,'Maximum block size',i10,' doubles')
 9050 format(1x,'Maximum block size supplied by input')
c
c ccsd_act/eomccsd_act
c orbital energies
ccx      if(model.eq.'ccsd_act') then
ccx      if(nodezero) then
ccx       write(6,*)'--------------------------------'
ccx       write(6,*)'ORBITAL ENERGIES OF CORRELATED ORBITALS'
ccx       do i=nfc(1)+1,nmo(1)-nfv(1)
ccx         write(6,200) i,dbl_mb(k_evl(1)+i-1)
ccx       enddo
ccx       write(6,*)'--------------------------------'
c
ccx      end if
ccx      end if
ccx 200  format(i6,2x,f20.8)
c 
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
        call errquit("tce_tile: dyn. frozen virtuals",1,MA_ERR)
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
 8778  format('Orbital energy of the highest virtual orbital',2x,f17.8)
c
c  initial zeroing
       do i=1,max_size
        active_tiles(i)=0
        active_tiles_n(i)=0
       enddo
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
c->d3p975
c
      if(.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1 +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted MO index",
     2 l_mo_index,k_mo_index))
     3 call errquit("tce_tile: MA problem",4,ma_err)
c
c<-d3p975
c 
      if(.not.intorb) then
c
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
c              int_mb(k_mo_index+i-1)=hole
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
c              int_mb(k_mo_index+i-1)=hole
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
c              int_mb(k_mo_index+i-1)=hole
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
c              int_mb(k_mo_index+i-1)= hole
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
c              int_mb(k_mo_index+i-1)=particle
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
c              int_mb(k_mo_index+i-1)=particle
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
c              int_mb(k_mo_index+i-1)=particle
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
c              int_mb(k_mo_index+i-1)=particle
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
      end if !not intorb
c
c
c
c
c
c
c
c
c
c
c
c INTORB can be used in CC calculation with 
c RHF or ROHF reference. Cannot be used with active-space CC
c approaches.
c In order to construct the RHF or ROHF tiles
c active space tiles algoritm was adopted
c
      if(intorb) then
c ccsd_act/eomccsd_act
      if(model.eq.'ccsd_act') then !----------------------------
c
c     Hole Alpha
c
      i = 0
      noa = 0
      do sym = 0, nirreps
        j_ia = 0           !occ. inactive index
        j_an = 0           !occ. semi-active
        j_a  = 0           !occ. active index
        do hole = nfc(1)+1, nocc(2)-oact
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
c              int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        do hole = nocc(2)-oact+1,nocc(2)
          if (int_mb(k_irs(1)+hole-1).eq.sym) then
             i = i + 1
             j_an = j_an + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+hole-1)
             call ga_get(g_movecs(1),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
c->d3p975
               int_mb(k_mo_index+i-1)=2*hole-1
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        do hole =nocc(2)+1,nocc(1) !FOR RHF this part does not contribute
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
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        j=j_ia+j_an+j_a
        if (j.gt.0) then
           nblocks_ia = j_ia/isize
           nblocks_an = j_an/isize
           nblocks_a  =  j_a/isize
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           if (j_an .gt. isize*nblocks_an) nblocks_an = nblocks_an+1
           if (j_a  .gt. isize*nblocks_a)  nblocks_a  = nblocks_a+1
           l_ia = 0
           do k = 1,nblocks_ia
             noa = noa + 1
             active_tiles(noa) = 0
             active_tiles_n(noa) = 0
             spin_tmp(noa) = 1
              spin_tmp_alpha(noa) = 1
             sym_tmp(noa) = sym
              sym_tmp_alpha(noa) = sym
             range_tmp(noa) = k*j_ia/nblocks_ia-l_ia
              range_tmp_alpha(noa) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa)
           enddo
           l_an = 0
           do k = 1,nblocks_an
             noa = noa + 1
             active_tiles(noa) = 0
             active_tiles_n(noa) = 1
             spin_tmp(noa) = 1
              spin_tmp_alpha(noa) = 1
             sym_tmp(noa) = sym
              sym_tmp_alpha(noa) = sym
             range_tmp(noa) = k*j_an/nblocks_an-l_an
              range_tmp_alpha(noa) = k*j_an/nblocks_an-l_an
             l_an = l_an + range_tmp(noa)
           enddo
           l_a = 0
           do k = 1,nblocks_a
             noa = noa + 1
             active_tiles(noa) = 1
             active_tiles_n(noa) = 1
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
        j_an = 0
        do hole = nfc(ipol)+1, nocc(ipol)-oact !beta electrons here no active part
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
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        do hole = nocc(ipol)-oact+1,nocc(ipol) !beta electrons semi active part
          if (int_mb(k_irs(ipol)+hole-1).eq.sym) then
             i = i + 1
             j_an = j_an + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+hole-1)
             call ga_get(g_movecs(ipol),1,nbf,hole,hole,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
c->d3p975
               int_mb(k_mo_index+i-1)=2*hole
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        j=j_ia+j_an
        if (j.gt.0) then
           nblocks_ia = j_ia/isize
           nblocks_an = j_an/isize
           if (j_ia .gt. isize*nblocks_ia) nblocks_ia = nblocks_ia+1
           if (j_an .gt. isize*nblocks_an) nblocks_an = nblocks_an+1
           l_ia = 0
           do k = 1,nblocks_ia
             nob = nob + 1
             active_tiles(noa+nob) = 0
             active_tiles_n(noa+nob) = 0
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             range_tmp(noa+nob) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob)
           enddo
           l_an = 0
           do k = 1,nblocks_an
             nob = nob + 1
             active_tiles(noa+nob) = 0
             active_tiles_n(noa+nob) = 1
             spin_tmp(noa+nob) = 2
             sym_tmp(noa+nob) = sym
             range_tmp(noa+nob) = k*j_an/nblocks_an-l_an
             l_an = l_an + range_tmp(noa+nob)
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
        j_an = 0 
        j_ia = 0
        do particle = nocc(1)+1,nocc(1)+uact 
          if (int_mb(k_irs(1)+particle-1).eq.sym) then
             i = i + 1
             j_an = j_an + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=1
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(1)+particle-1)
             call ga_get(g_movecs(1),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
c->d3p975
               int_mb(k_mo_index+i-1)=2*particle-1
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
          endif
        enddo
        do particle = nocc(1)+uact+1, nmo(1)-nfv(1)
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
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
          endif
        enddo
        j=j_an+j_ia
        if (j.gt.0) then
           nblocks_an  =  j_an/isize
           nblocks_ia  =  j_ia/isize
           if (j_an  .gt. isize*nblocks_an)  nblocks_an  = nblocks_an+1
           if (j_ia  .gt. isize*nblocks_ia)  nblocks_ia  = nblocks_ia+1
           nblocks = nblocks_ia+nblocks_an
           l_an = 0
           do k = 1,nblocks_an
             nva = nva + 1
             active_tiles(noa+nob+nva) = 0
             active_tiles_n(noa+nob+nva) = 1
             spin_tmp(noa+nob+nva) = 1
              spin_tmp_alpha(noa+nva) = 1
             sym_tmp(noa+nob+nva) = sym
              sym_tmp_alpha(noa+nva) = sym
             range_tmp(noa+nob+nva) = k*j_an/nblocks_an-l_an
              range_tmp_alpha(noa+nva) = k*j_an/nblocks_an-l_an
             l_an = l_an + range_tmp(noa+nob+nva)
           enddo
           l_ia = 0
           do k = 1,nblocks_ia
             nva = nva + 1
             active_tiles(noa+nob+nva) = 0
             active_tiles_n(noa+nob+nva) = 0
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
c     Particle Beta
c
      nvb = 0
      do sym = 0, nirreps
        j_ia = 0
        j_an = 0
        j_a  = 0
        do particle = nocc(ipol)+1, nocc(1) !for RHF this part does not contribute
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
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
          endif
        enddo
        do particle = nocc(1)+1, nocc(1)+uact
          if (int_mb(k_irs(ipol)+particle-1).eq.sym) then
             i = i + 1
             j_an = j_an + 1
             int_mb(k_irs_sorted+i-1)=sym
             int_mb(k_spin_sorted+i-1)=2
             dbl_mb(k_evl_sorted+i-1)=dbl_mb(k_evl(ipol)+particle-1)
             call ga_get(g_movecs(ipol),1,nbf,particle,particle,
     1         dbl_mb(k_movecs_sorted+(i-1)*nbf),nbf)
c->d3p975
               int_mb(k_mo_index+i-1)=2*particle
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
          endif
        enddo
        do particle = nocc(1)+uact+1, nmo(ipol)-nfv(ipol)
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
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
          endif
        enddo
        j=j_a+j_an+j_ia
        if (j.gt.0) then
           nblocks_a   = j_a/isize
           nblocks_an  = j_an/isize
           nblocks_ia  =  j_ia/isize
c           nblocks = nblocks_a+nblocks_ia
           if (j_a .gt. isize*nblocks_a) nblocks_a = nblocks_a+1
           if (j_an .gt. isize*nblocks_an) nblocks_an = nblocks_an+1
           if (j_ia  .gt. isize*nblocks_ia)  nblocks_ia  = nblocks_ia+1
           nblocks = nblocks_a+nblocks_an+nblocks_ia
           l_a = 0
           do k = 1,nblocks_a
             nvb = nvb + 1
             active_tiles(noa+nob+nva+nvb) = 1
             active_tiles_n(noa+nob+nva+nvb) = 1
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             range_tmp(noa+nob+nva+nvb) = k*j_a/nblocks_a-l_a
             l_a = l_a + range_tmp(noa+nob+nva+nvb)
           enddo
           l_an = 0
           do k = 1,nblocks_an
             nvb = nvb + 1
             active_tiles(noa+nob+nva+nvb) = 0
             active_tiles_n(noa+nob+nva+nvb) = 1
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             range_tmp(noa+nob+nva+nvb) = k*j_an/nblocks_an-l_an
             l_an = l_an + range_tmp(noa+nob+nva+nvb)
           enddo
           l_ia = 0
           do k = 1,nblocks_ia
             nvb = nvb + 1
             active_tiles(noa+nob+nva+nvb) = 0
             active_tiles_n(noa+nob+nva+nvb) = 0
             spin_tmp(noa+nob+nva+nvb) = 2
             sym_tmp(noa+nob+nva+nvb) = sym
             range_tmp(noa+nob+nva+nvb) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + range_tmp(noa+nob+nva+nvb)
           enddo
        endif
      enddo
      go to 1111
      end if !!!  (model.eq.'ccsd_act') !-----------------------
c
c
c
c still in intorb true (now general case)
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
c->d3p975
               int_mb(k_mo_index+i-1)=2*hole-1
c              int_mb(k_mo_index+i-1)=hole
c<-d3p975
          endif
        enddo
        do hole =nocc(2)+1,nocc(1) !FOR RHF this part does not contribute
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
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
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
c->d3p975
               int_mb(k_mo_index+i-1)=2*hole
c               int_mb(k_mo_index+i-1)=hole
c<-d3p975
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
c->d3p975
               int_mb(k_mo_index+i-1)=2*particle-1
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
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
        do particle = nocc(ipol)+1, nocc(1) !for RHF this part does not contribute
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
c               int_mb(k_mo_index+i-1)=particle
c<-d3p975
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
c->d3p975
               int_mb(k_mo_index+i-1)=2*particle
c               int_mb(k_mo_index+i-1)=particle
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
c
c ccsd_act/eomccsd_act
 1111  continue
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
c
c
c
      end if !intorb
c
c
c
c
c
c
c
c
c
c
c
c
c
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
      if(activecalc.or.(.not.intorb).or.(model.eq."ccsd_act")) THEN
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Active blocks ",
     1  l_active,k_active)) call errquit("tce_tile: MA problem",9,
     2  MA_ERR)
      end if
      if(model.eq."ccsd_act") THEN
      if (.not.ma_push_get(mt_int,noa+nva,"Active blocks alpha",
     1  l_active_o,k_active_o)) call errquit("tce_tile: MA problem",9,
     2  MA_ERR)
      end if
      if(intorb) THEN
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"alpha-izer 2 ",
     1  l_b2am,k_b2am)) call errquit("tce_tile: MA problem",9,
     2  MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha spins ",
     1  l_spin_alpha,k_spin_alpha)) 
     2  call errquit("tce_tile: MA problem",9,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha syms ",
     1  l_sym_alpha,k_sym_alpha))
     2  call errquit("tce_tile: MA problem",9,MA_ERR) 
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha range ",
     1  l_range_alpha,k_range_alpha))
     2  call errquit("tce_tile: MA problem",9,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nva,"test for alpha range ",
     1  l_offset_alpha,k_offset_alpha))
     2  call errquit("tce_tile: MA problem",9,MA_ERR)
      end if
c
c
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
c
c
c
c
c
c
      if(intorb) THEN
c
      j = 0
      do i = 1,noa+nob+nva+nvb
        int_mb(k_spin+i-1)=spin_tmp(i)
        int_mb(k_sym+i-1)=sym_tmp(i)
        int_mb(k_range+i-1)=range_tmp(i)
c        int_mb(k_active+i-1)=active_tiles(i) !l_active will not be used
c ccsd_act/eomccsd_act
        if(model.eq."ccsd_act") then 
        int_mb(k_active+i-1)=active_tiles_n(i) !l_active will not be used in stand.app.
        end if
c
        int_mb(k_offset+i-1)=j
        j = j + range_tmp(i)
      enddo
c ccsd_act/eomccsd_act
      if(model.eq."ccsd_act") then
       do i = 1,noa
         int_mb(k_active_o+i-1)=active_tiles_n(i)
       enddo
       do i = 1,nva
         int_mb(k_active_o+i+noa-1)=active_tiles_n(i+noa+nob)
       enddo
      end if
c
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
c --- debug ----
c      if(nodezero) then
c       write(6,*)'from tce_tile----------------'
c       write(6,*)'k_spin_alpha'
c       do i=1,noa+nva
c        write(6,*) spin_tmp_alpha(i),int_mb(k_spin_alpha+i-1)
c       enddo
c       write(6,*)'k_sym_alpha'
c       do i=1,noa+nva
c        write(6,*) sym_tmp_alpha(i),int_mb(k_sym_alpha+i-1)
c       enddo
c       write(6,*)'k_range_alpha'
c       do i=1,noa+nva
c        write(6,*) range_tmp_alpha(i),int_mb(k_range_alpha+i-1)
c       enddo
c       write(6,*)'k_offset'
c       do i=1,noa+nva
c        write(6,*) int_mb(k_offset_alpha+i-1)
c       enddo
c       call util_flush(6)
c      end if
c --------------
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
 9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha',/,
     1         1x,'-------------------------------------------------')
 9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8)
 9080 format(/,1x,'tile_dim = ',2x,i4)
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
