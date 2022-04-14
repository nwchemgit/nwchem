      subroutine texas_set_accy(pnl_thresh1)
      implicit none
      double precision pnl_thresh1 ! [input] accy threshold 1 for texas
c
      call setup_thresh(pnl_thresh1)
c
      end
      subroutine texas_terminate()
c $Id$
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
#include "mafdecls.fh"
#include "util.fh"
c----------------------------------------------------------------
      logical cent2,cent3,cent4
      common /what_was_calc/ cent2,cent3,cent4
      common /neglect/ eps,eps1,epsr
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_prep2 , k_txs_prep2 
c----------------------------------------------------------------------
c:added:rak:fix linxmem size to basis set
      integer texas_max_ang_val
      common /c_texas_max_ang_val/ texas_max_ang_val
c----------------------------------------------------------------------
      common /pnl_time/ time4texas_hf2_m,time4mul_quart
      common /txs_time/ timprep2,timcalc2,timblok2
c----------------------------------------------------------------
      common /pnl_nqrt/ ncall_pnl,ncall_41,nqrts_pnl,npart_pnl,nsize_pnl
c----------------------------------------------------------------
      common /pnl000/ xbluse,nbluse
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /ilepar/ lpar_exe,lpar_num
c----------------------------------------------------------------
      common /intgop/ ncache,maxprice,iprint,iblock
c----------------------------------------------------------------
      common /mem_pnl_scr/ nall_peak,mark_peak,mem_peak
c----------------------------------------------------------------
c 1998 (numerical instability in Tracy's recursive)
c
      common /count_primitive/ x_prim_all,
     *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
      common /lobsa_times/ xlobsa2,xlobsa4
C-----------------------------------------------------------------------
c     common /count_efficient/ xobsa2,xobsa2_ineff,xobsa2_more,
c    *                         xobsa4,xobsa4_ineff,xobsa4_more
c----------------------------------------------------------------
c
      if(ncall_pnl.gt.0) then 
c        average block-size in texas calculations for pnl 
c
         nvect_pnl=nsize_pnl/ncall_pnl
c
c        efficiency of using block's suructure :
c
         bleffiy=(xbluse/dble(ncall_pnl))*100.d0
      endif
c----------------------------------------------------------------
      texas_max_ang_val = -565   ! reset on termination
c----------------------------------------------------------------
c Information about integral's program performance :
c
c
      IF(iprint.gt.0) THEN
c
      write(8,*)'------------------------------------------------------'
      write(8,*)'Texas-int2 program has been called ',ncall_pnl,' times'
      if(ncall_pnl.eq.0) goto 00001
      write(8,*)'(including one-quartet-at-the-time ',ncall_41,' times)'
      write(8,*)'integrals have been calculated in  ',npart_pnl,' pass.'
      write(8,*)'number of processed shell quartets ',nqrts_pnl
      write(8,*)'average block-size (vector length) ',nvect_pnl
      write(8,18) bleffiy
   18 format(' blocking utilities of Texas used in ',f6.2,' %')
      if(cent2) write(8,*)'(type of integrals : 2-center two-electron)'
      if(cent3) write(8,*)'(type of integrals : 3-center two-electron)'
      if(cent4) write(8,*)'(type of integrals : 4-center two-electron)'
c----------------------------------------------------------------
      write(8,*)'------------------------------------------------------'
      write(8,*)'Pre-calculations for contracted shell pairs have been'
      write(8,*)'executed ',lpar_exe,' times for ',lpar_num,' pairs'
      write(8,*)'------------------------------------------------------'
      write(8,*)'Total number of contracted shell quartets ',ntotal
      write(8,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - - '
      write(8,*)'number of contracted quartets requested   ',nospec
      write(8,*)'reduced by the screening procedure to     ',nopres
      write(8,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - - '
      write(8,*)'number of primitive quartets requested    ',nrimtot
      write(8,*)'reduced by the screening procedure to     ',nrimret
      write(8,*)'------------------------------------------------------'
      write(8,*)'screening has been done with threshold  = ',eps     
      write(8,*)'------------------------------------------------------'
c
c CPU timing :
c
      write(8,*)'------------------------------------------'
      write(8,*)' CPU time of using TEXAS-integral program '
      write(8,*)' no of calls =',ncall_pnl,' no of quartets =',nqrts_pnl
      write(8,*)'------------------------------------------'
c
c..................................................................
c these below are timings for :
c
c     request_update + switch_scr + labels_udate (see texas_hf2_m)
c     time_pnl_2_txs=time4texas_hf2_m - time4mul_quart
c
c and txs_setup (see mul_quart)
c 
c     time_txs_setup=time4mul_quart   - timcalc2
c..................................................................
c
      time_interface=time4texas_hf2_m - timcalc2
c
      write(8,140) time_interface,timcalc2,time4texas_hf2_m
  140 format('  nwch-txs-nwch interface  time =',f9.1/
     *       '  txs integral calculation time =',f9.1/
     *       '  Total integral delivery  time =',f9.1)
      write(8,*)'------------------------------------------'
c
      call timepr(timblok2 + 0.001d0)
c
c----------------------------------------------------------------
c Check out memory status as r and and reserved in the TXS-BL() :
c
      call memstat(nreque_bl,nmark_bl,maxmem_bl,memtot_bl)
c
c take only memtot_bl from the call above, rest from 
c the common /mem_pnl_scr/
c
      write(8,*)' -------------------------------------------'
      write(8,*)'       LEAVING TEXAS-INTEGRAL PROGRAM'
      write(8,*)'       Statistics of memory in blscr:'
      write(8,*)'                                               '
ctest ?
c     write(8,*)' max number of allocations =',nreque_bl
c     write(8,*)' max number of marks       =',nmark_bl
c     write(8,*)' max memory used overall   =',maxmem_bl,' (high water)'
c     write(8,*)' total memory available    =',memtot_bl
ctest ?
      write(8,*)' max number of allocations =',nall_peak
      write(8,*)' max number of marks       =',mark_peak
      write(8,*)' max memory used overall   =',mem_peak,' (high water)'
      write(8,*)' total memory available    =',memtot_bl
      write(8,*)' -------------------------------------------'
c
c------------------------------------
c 1998 (numerical instability in Tracy's recursive)
c print out info about numerical instability in Tracy's recursive:
c
      write(8,9871) x_prim_all,x_prim_wr1,x_prim_wr2,x_prim_wr3
 9871 format(' total number of primitive quartets  =',f15.0/
     *       ' number of unstable quartets after 1 =',f15.0/
     *       ' number of unstable quartets after 2 =',f15.0/
     *       ' number of unstable quartets after 3 =',f15.0)
      write(8,9872) xlobsa2,xlobsa4
 9872 format(' number of shift cases: cent1 > cent3=',f15.0/
     *       ' number of shift cases: cent3 > cent1=',f15.0)
c------------------------------------
c     write(8,9872) xobsa2,xobsa2_ineff,xobsa2_more,
c    *              xobsa4,xobsa4_ineff,xobsa4_more
c9872 format(' number of cases (shift:cent1>cent3=',f15.0/
c    *       '       inefficient (for L_12<L_34) =',f15.0/
c    *       '       number of additional shifts =',f15.0/
c    *       ' number of cases (shift:cent3>cent1=',f15.0/
c    *       '       inefficient (for L_12>L_34) =',f15.0/
c    *       '       number of additional shifts =',f15.0)
c------------------------------------
00001 continue
c------------------------------------
      write(8,*)' -------------------------------------------'
      write(8,*)' TEXAS INTEGRAL PROGRAM HAS BEEN TERMINATED'
      write(8,*)' -------------------------------------------'
c
c     call util_flush(8)
c      close(8)
c------------------------------------
      ENDIF     !! end of iprint 
c----------------------------------------------------------------
c free up memory allocated by texas_init      
c
*rak:      write(6,*)'t h_txs_basnuc ',h_txs_basnuc 
*rak:      write(6,*)'t h_txs_mapp1  ',h_txs_mapp1
      if (h_txs_basnuc .ne. -1) then
         if (.not.MA_free_heap(h_txs_basnuc))
     &        call errquit('texas_terminate: failed ?', 911, INT_ERR)
      endif
      if (h_txs_mapp1 .ne. -1) then
         if (.not.MA_free_heap(h_txs_mapp1 ))
     &        call errquit('texas_terminate: failed ?', 911, INT_ERR)
      endif
      if (h_txs_prep2 .ne. -1) then
         if (.not.MA_free_heap(h_txs_prep2 ))
     &        call errquit('texas_terminate: failed ?', 911, INT_ERR)
      endif
c
      h_txs_basnuc = -1
      h_txs_mapp1  = -1 
      h_txs_prep2  = -1
      k_txs_basnuc =  0
      k_txs_mapp1  =  0
      k_txs_prep2  =  0
c
      end
c
c End of terminate 
c================================================================
      subroutine texas_init(rtdb,nbas,bases,nq_init,l_blsize,int_type)
c---------------------------------------------------------------
c l_blsize = size of the scratch array which will be used when
c             texas_hf2_m is called
c this is maximum size of the scratch needed to do one super-block
c
c int_type = type of a comming task : 
c            scfd_int,giao_int,der1_int,der2_int
c---------------------------------------------------------------
      implicit none
#include "errquit.fh"
c
c..  this routine is the current interface to the texas code
c
#include "bas.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "apiP.fh"
#include "global.fh"
c::functions
c::passed
      character*8 int_type
      integer rtdb ! rtdb handle
      integer nbas ! number of basis sets 
      integer bases(nbas) ! basis set handles
      integer nq_init ! [input] maximum number of quartets ever.
      integer l_blsize
c::local
      integer geom, basis 
      integer nat, natoms ! number of atoms
      integer npshells ! total number of primitive shells
      integer ncshells ! total number of contracted shells
      integer npfunctx ! total number of primitive functions
      integer ncfunctx ! total number of contracted functions
c
      integer npsh, ncsh 
      integer ibasis
c
      integer size_from_ma
      integer inuc,ibas
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_prep2 , k_txs_prep2 
c--------------------------------------------------------------------
c?    integer basis_init
c?    common /c_basis_init/ basis_init
c
      integer num_bas_1,num_bas_2,num_bas_3
      integer ncs_bas_1,ncs_bas_2,ncs_bas_3
      integer nps_bas_1,nps_bas_2,nps_bas_3
      integer nat_bas_1,nat_bas_2,nat_bas_3
      integer ncf_bas_1,ncf_bas_2,ncf_bas_3
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
c--------------------------------------------------------------------
c:added:rak:fix linxmem size to basis set
      integer texas_max_ang_val
      common /c_texas_max_ang_val/ texas_max_ang_val
c--------------------------------------------------------------------
      integer rtdb_copy
      common /c_rtdb_copy/ rtdb_copy
      integer bas
c--------------------------------------------------------------------
      integer txs_static, txs_dynam
c--------------------------------------------------------------------
c2002
      integer ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
c--------------------------------------------------------------------
c
      texas_max_ang_val = -565  ! ensure it is reset
      rtdb_copy = rtdb
c--------------------------------------------------------------------
c     if (nbas.gt.1 .and. ga_nodeid().eq.0) then
c       write(6,*)
c       write(6,*)' !! warning: texas can only handle one basis set'
c       write(6,*)
c       call util_flush(6)
c     endif
c
c     basis = bases(1) ! only one basis in texas
c     basis_init = basis
c
c     bas = basis + BASIS_HANDLE_OFFSET
c     geom = ibs_geom(bas)
c
c     nat = ncenter(geom)
c----------------------------------------------------------------------
c Send info about forthcomming task i.e. what type of integrals will be
c calculated. Int_type is used to estimate memory needed for that task.
c The memory requiment is taken into account when blocks of pairs are
c constructed. 
c     Whatever the value of int_type is passed here,the program will be
c able to handle ALL lower level tasks (e.g. if it is called here for
c int_type='der1_int', blocks will be not too big to handle ordinary
c 2-el.integrals).
c Memory needed for different kind of integrals changes as follows :
c
c   scfd_int < giao_int < der1_int < der2_int
c
c All lower than here, values of int_type parameter can be used when the
c texas_hf2_m routine is called for a particular type of calculations.

      call initialized_task(int_type)
c----------------------------------------------------------------------
c multi basis set : Total number of Atoms,Contracted & Primitive shells 
c
      num_bas_1=0
      ncs_bas_1=0
      nps_bas_1=0
      nat_bas_1=0
c
      num_bas_2=0
      ncs_bas_2=0
      nps_bas_2=0
      nat_bas_2=0
c
      num_bas_3=0
      ncs_bas_3=0
      nps_bas_3=0
      nat_bas_3=0
c
      natoms=0
      ncshells=0
      npshells=0
      do 100 ibasis=1,nbas
         basis = bases(ibasis) 
         bas = basis + BASIS_HANDLE_OFFSET
         geom  = ibs_geom(bas)
c
         natoms   = natoms   + ncenter(geom)
         ncshells = ncshells + ncont_tot_gb(bas)
         npshells = npshells + nprim_tot_gb(bas)
c
         if(ibasis.eq.1) then
            num_bas_1=basis
            ncs_bas_1=ncshells
            nps_bas_1=npshells
            nat_bas_1=natoms
         endif
         if(ibasis.eq.2) then
            num_bas_2=basis
            ncs_bas_2=ncshells
            nps_bas_2=npshells
            nat_bas_2=natoms
         endif
         if(ibasis.eq.3) then
            num_bas_3=basis
            ncs_bas_3=ncshells
            nps_bas_3=npshells
            nat_bas_3=natoms
         endif
         if(ibasis.gt.3) then
           if (ga_nodeid().eq.0) then
             write(6,*) ' too many basis sets : allowed=3, passed=',
     &           ibasis
             write(6,*) ' execution stopped in texas_init'
             call util_flush(6)
           endif
           call errquit('texas_face: execution stopped in texasinit',0,
     &       INT_ERR)
         endif
  100 continue
c----------------------------------------------------------------------
c Allocate memory for basis set & nuclear data :(5*natoms+13*npshells)
c Add one one more atom and one more shell for 2-&3-center integrals :
c
      size_from_ma = 5*(natoms+1) + 13*(npshells+1)
c
      if (.not. MA_alloc_get(MT_DBL, size_from_ma, 'texas basis/nuc ',
     &                       h_txs_basnuc, k_txs_basnuc))
     &    call errquit('texas_init: alloc for texas basis/nuc failed',
     &    911, MEM_ERR)
      call dfill(size_from_ma,0.0d00,dbl_mb(k_txs_basnuc),1)
      inuc = k_txs_basnuc
c>>>  ibas = inuc + 5*natoms
      ibas = inuc + 5*(natoms+1)
c----------------------------------------------------------------------
c make bl and inx arrays 
c
      call blinx(natoms,ncshells,npshells,inuc,ibas)
c----------------------------------------------------------------------
c make DATNUC , DATBAS and INX arrays :
c
      natoms=0
      ncshells=0
      npshells=0
      npfunctx=0  ! contr.func. counter , sum up inside
      ncfunctx=0  ! primi.func. counter , sum up inside
      do 200 ibasis=1,nbas
         basis = bases(ibasis) 
         bas   = basis + BASIS_HANDLE_OFFSET
         geom  = ibs_geom(bas)
         nat   = ncenter(geom)
c
         call makenuc(geom,natoms, nat,dbl_mb(inuc))
         natoms=natoms + nat
c
         ncsh = ncont_tot_gb(bas)
         npsh = nprim_tot_gb(bas)
         call makebas(ncshells,npshells, ncsh,bas, dbl_mb(ibas))
c
         call makeinx(npfunctx,ncfunctx,ncshells,ncsh,bas,basis)
c
         ncshells = ncshells + ncsh
         npshells = npshells + npsh
  200 continue
c----------------------------------------------------------------------
c2002
         ntotal=ncshells*(ncshells+1)/2
         ntotal=ntotal*(ntotal+1)/2
c----------------------------------------------------------------------
c test only : print INX,DATNUC,DATBAS arrays :
c
c     call print_inx_nuc_bas(nbas,dbl_mb(inuc),dbl_mb(ibas))
c
c----------------------------------------------------------------------
c Add one uncontracted s-shell with zero exponent for 2- and 3-center 
c integrals . 
*:rak Do so Even if there is only ONE basis set specified :
c
*:rak always add the zero s function
*:rak      if(nbas.gt.1) call add_s_zero( dbl_mb(ibas))
c2002 call add_s_zero( dbl_mb(ibas))
c
           if(nbas.gt.1) call add_s_zero( dbl_mb(ibas))
c----------------------------------------------------------------------
c Re-arrange a basis set according to the Texas convention
c 1) normalize and 2) re-order
c
c2002 call txsarrag(rtdb, nq_init, l_blsize)
      call txsarrag(rtdb, nq_init, l_blsize,int_type)
c
c----------------------------------------------------------------------
c print out basis sets info (in fort.8 if iprint is >0)
c
      call print_texas_nbas(nbas)
c----------------------------------------------------------------------
c... map texas memory to memory commons of the PNL/int API
c
*      write(6,*)' mem_2e4c :',mem_2e4c
*      write(6,*)' mem_2e3c :',mem_2e3c
*      write(6,*)' mem_2e2c :',mem_2e2c
*      write(6,*)' l_blsize :',l_blsize
*      mem_2e4c = max(mem_2e4c,l_blsize)
*      mem_2e3c = max(mem_2e3c,l_blsize)
*      mem_2e2c = max(mem_2e2c,l_blsize)
c tmp fix for memory split
c----------------------------------------------------------------------
      iszb_2e4c = isz_2e4c  ! use miminum of 1 quartet. 
      memb_2e4c = l_blsize
      call texas_memory (txs_static,1,txs_dynam)
      mem_2e4c = max(mem_2e4c,txs_static+txs_dynam)
      mem_2e3c = max(mem_2e3c,txs_static+txs_dynam)
      mem_2e2c = max(mem_2e2c,txs_static+txs_dynam)
c----------------------------------------------------------------------
c ONLY for TESTS : calculate two-electron integrals  
c
c
c        call test_cent2(num_bas_1,num_bas_1, l_blsize)
c        call test_cent2(num_bas_1,num_bas_2, l_blsize)
c        call test_cent2(num_bas_2,num_bas_1, l_blsize)
c        call test_cent2(num_bas_2,num_bas_2, l_blsize)
c
c         call test_cent2(num_bas_3,num_bas_3, l_blsize)
c         call test_cent2(num_bas_1,num_bas_3, l_blsize)
c         call test_cent2(num_bas_3,num_bas_1, l_blsize)
c
c        call test_cent3i_kl(num_bas_1,num_bas_1, l_blsize)
c        call test_cent3i_kl(num_bas_1,num_bas_2, l_blsize)
c        call test_cent3i_kl(num_bas_2,num_bas_2, l_blsize)
c
c        call test_cent3ij_k(num_bas_1,num_bas_1, l_blsize)
c        call test_cent3ij_k(num_bas_1,num_bas_2, l_blsize)
c        call test_cent3ij_k(num_bas_2,num_bas_2, l_blsize)
c
c        call test_cent4(num_bas_1,num_bas_1, l_blsize)
c        call test_cent4(num_bas_1,num_bas_2, l_blsize)
c        call test_cent4(num_bas_2,num_bas_1, l_blsize)
c        call test_cent4(num_bas_2,num_bas_2, l_blsize)
c
c End of tests; program will be stopped here
c
c       call texas_terminate()
c       STOP ' stopped by kw in texas_face '
c
c End of tests; program stopped here
c----------------------------------------------------------------------
*rak:      write(6,*)'i h_txs_basnuc ',h_txs_basnuc
*rak:      write(6,*)'i h_txs_mapp1  ',h_txs_mapp1
c----------------------------------------------------------------------
      end
c=================================================================
      subroutine blinx(natom,ncshell,npshell,inucx,ibasx)
      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c the common block /big/ bl(lcx) has been repleaced by pnl scratch 
c-----------------------------------------------------------------
      parameter(mxsh= 190 000)
c------------------------------------
      common /ctxs_index/ maxsh,ifp,inx(mxsh)
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c------------------------------------
      na=natom
      ncs=ncshell
      nsh=npshell
      inuc=inucx 
      ibas=ibasx 
c------------------------------------
      maxsh=mxsh
c------------------------------------
      call retinall
c------------------------------------
c reserve memory (in index) for inx
c
      call getint(12*(ncshell+1),inxx)
c------------------------------------
      end
c=================================================================
cccc  subroutine makenuc(geom,nat,datnuc)
      subroutine makenuc(geom,natoms,nat,datnuc)
      implicit real*8 (a-h,o-z)
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
ccccc integer geom, nat
cccc  double precision datnuc(5,nat)
c
      integer geom, natoms, nat
      double precision datnuc(5,*)
c
c     write(6,*)'entry MAKENUC : geom=',geom,' natoms,nat=',natoms,nat
c
      do 100 iat=1,nat
      iatoms=natoms+iat
      datnuc(1,iatoms)=charge(iat,geom)
      datnuc(2,iatoms)=coords(1,iat,geom)
      datnuc(3,iatoms)=coords(2,iat,geom)
      datnuc(4,iatoms)=coords(3,iat,geom)
      datnuc(5,iatoms)=iat
cccc  write(6,66) datnuc(5,iatoms),(datnuc(ii,iat),ii=1,4)
  100 continue
   66 format('atom Name=',i3,' charge=',f5.1,' xyz(au)=',3f10.5)
c
      end
c=================================================================
cccc  subroutine makebas(nshells,bas,datbas)
      subroutine makebas(ncshells,npshells,ncsh,bas,datbas)
c
c npshells = total number of primitive  shells up to this basis set
c ncshells = total number of contracted shells up to this basis set
c ncsh     = number of shells in this basis set (bas)
c
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
#include "nwc_const.fh"
#include "basP.fh"
      integer bas
      logical cart_2_sphe
      logical spcart_init
      external spcart_init
      common /ctxs_index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      dimension datbas(13,*)
c---------------------------------------------------------------
c test
c     write(6,*)' makebas entry :bas=',bas,
c    * '   ncshells,npshells,ncsh=', ncshells,npshells,ncsh
c---------------------------------------------------------------
c shell's type :
c
c ityp1= 1   2   3   4   5   6   7   8   9   10    !  texas type 
c        s   p   l   d5  d6  f7  f10 g15 h21 i28
c itype= 0   1   -1      2       3   4   5   6     !  pnl type
c---------------------------------------------------------------
c basis set type (cartesian or spherical harmonics)
c
      cart_2_sphe=bas_spherical(bas)
c---------------------------------------------------------------
c
      ityp_max=0
      ipshell = npshells
      do 100 ish=1,ncsh
         call basdat(bas,ncshells,ish,datbas,inx(inxx),ipshell,
     *               cart_2_sphe,itype)
         if(itype.gt.ityp_max) ityp_max=itype
  100 continue
c
c---------------------------------------------------------------
c
       if( cart_2_sphe ) then
*           write(6,*)' Basis Set no =',bas,' is Spherical Harmonics'
           if( spcart_init(ityp_max,.true.,.false.)) then
*              write(6,*)' transformation matrix initialized'
           else
              call errquit('Cart_Sphe.-init: failed ?', 911, INT_ERR)
           endif
       else
*           write(6,*)' Basis Set no =',bas,' is Cartesian '
       endif
c------------------------------------
c
      end
c==============1
      subroutine basdat(bas,ncshells,ish,datbas,inx, ipshell,
     *                  cart_2_sphe,itype)
      implicit real*8 (a-h,o-z)
****
**  This routine should be changed to be only the interface calls
**  RAK
****
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
c
c----------------------------------------------------------------------
c:added:rak:fix linxmem size to basis set
      integer texas_max_ang_val
      common /c_texas_max_ang_val/ texas_max_ang_val
c----------------------------------------------------------------------
      logical cart_2_sphe
      integer bas
      dimension datbas(13,*),inx(12,*)
      dimension isize_cart(10),isize_sphe(10)
c
      data isize_cart /1, 3, 4, 5, 6, 7, 10, 15, 21, 28 /
      data isize_sphe /1, 3, 4, 5, 5, 7,  7,  9, 11, 13 /
c
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"
c
c---------------------------------------------------------------
c shell's type :
c
c itype= 0   1   -1      2       3   4   5   6     !  pnl type
c        s   p   l   d5  d6  f7  f10 g15 h21 i28
c ityp1= 1   2   3   4   5   6   7   8   9   10    !  texas type 
c ityp1= 1   2   3   -   5   -   7   8   9   10    !  texas here
c only above TXS types can be used here
c---------------------------------------------------------------
c
      ics=ncshells+ish
c---------------------------------------------------------------
c test
c     write(6,*)' entry BASDAT bas=',bas,' ncshells,ish=',ncshells,ish,
c    *' icshell=',ics,' ipshell=',ipshell
c---------------------------------------------------------------
c
      inx(1,ics)=ipshell
c>>>> inx(2,ics)=(sf_ibs_cn2ce(ics,bas))
      inx(2,ics)=(sf_ibs_cn2ce(ish,bas))
c PNL types:
c>>>  itype = infbs_cont(CONT_TYPE,(sf_ibs_cn2ucn(ics,bas)),bas)
      itype = infbs_cont(CONT_TYPE,(sf_ibs_cn2ucn(ish,bas)),bas)
      texas_max_ang_val = max(texas_max_ang_val,itype)
C TXS types:
      ityp1 = itype + 1
      if(itype.eq.-1) ityp1=3
      if(itype.eq.2 ) ityp1=5
      if(itype.gt.2 ) ityp1=itype+4
c
c TXS shell sizes :
c
      if( cart_2_sphe) then
         inx(3,ics)=isize_sphe(ityp1)
      else
         inx(3,ics)=isize_cart(ityp1)
      endif
c
c>>>  igen = infbs_cont(CONT_NGEN,(sf_ibs_cn2ucn(ics,bas)),bas)
      igen = infbs_cont(CONT_NGEN,(sf_ibs_cn2ucn(ish,bas)),bas)
      inx(4,ics)=igen-1
      if(ityp1.eq.3) inx(4,ics)=0
c>>>  iprim = infbs_cont(CONT_NPRIM,(sf_ibs_cn2ucn(ics,bas)),bas)
      iprim = infbs_cont(CONT_NPRIM,(sf_ibs_cn2ucn(ish,bas)),bas)
      inx(5,ics)=ipshell+iprim
c------------------
c     inx(6,ics)= symm.
c     inx(7,ics)= symm.
c     inx(8,ics)= symm.
c     inx(9,ics)= ! not used
c------------------
c     inx(10,ics)= ! this +1 gives the first contr. function
c     inx(11,ics)= ! this    gives the last  contr. function
c------------------
      inx(12,ics)=ityp1
c
c transfer primitive shells info :
c
c>>>  iexppnl = infbs_cont(CONT_IEXP,(sf_ibs_cn2ucn(ics,bas)),bas) - 1
c>>>  icofpnl = infbs_cont(CONT_ICFP,(sf_ibs_cn2ucn(ics,bas)),bas) - 1
      iexppnl = infbs_cont(CONT_IEXP,(sf_ibs_cn2ucn(ish,bas)),bas) - 1
      icofpnl = infbs_cont(CONT_ICFP,(sf_ibs_cn2ucn(ish,bas)),bas) - 1
      do 100 ipr=1,iprim
        ips=ipr+ipshell
*       datbas(1,ips)=ext(ipr)
        datbas(1,ips)=sf_exndcf(iexppnl+ipr,bas)
        if(ityp1.eq.3) then
           datbas(2,ips)=sf_exndcf(icofpnl+ipr,bas)
           datbas(3,ips)=sf_exndcf(icofpnl+ipr+iprim,bas)
        else
           do 150 igc=1,igen
            datbas(1+igc,ips)=
     &             sf_exndcf(icofpnl+ipr+((igc-1)*iprim),bas)
 150       continue
        endif
c>>>>   iat = (sf_ibs_cn2ce(ics,bas))
        iat = (sf_ibs_cn2ce(ish,bas))
        datbas(11,ips)=coords(1,iat,ibs_geom(bas))
        datbas(12,ips)=coords(2,iat,ibs_geom(bas))
        datbas(13,ips)=coords(3,iat,ibs_geom(bas))
 100  continue
c
      ipshell = ipshell+iprim
      end
c=================================================================
      subroutine makeinx(nfunp,nfunc,ncshells,ncsh,bas,basis)
c-----------------------------------------------------------------
c nfunp   = total number of PRIM.basis funct. up to this basis set
c nfunc   = total number of CONT.basis funct. up to this basis set
c
c ncshells = total number of CONT.shells up to this basis set
c ncsh     = number of cont. shells in this basis set
c-----------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
#include "nwc_const.fh"
#include "basP.fh"
      integer bas,basis
      logical cart_2_sphe
      common /ctxs_index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c-----------------------------------------------------------------
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
c-----------------------------------------------------------------
c
      cart_2_sphe=bas_spherical(bas)
c
      call inxdat(nfunp,nfunc,ncshells,ncsh,inx(inxx),cart_2_sphe,ncf)
c
c save ncf for each basis set:
c
      if(basis.eq.num_bas_1) ncf_bas_1=ncf
      if(basis.eq.num_bas_2) ncf_bas_2=ncf
      if(basis.eq.num_bas_3) ncf_bas_3=ncf
c
      nbf=nfunp
c
      end
c==============1
      subroutine inxdat(nfunp,nfunc,ncshells,ncsh,inx,cart_2_sphe,ncf)
c
c nfunp   = total number of PRIM.basis funct. up to this basis set
c nfunc   = total number of CONT.basis funct. up to this basis set
c
c ncshells = total number of CONT.shells up to this basis set
c ncsh     = number of cont. shells in this basis set
c
c
      logical cart_2_sphe
      dimension inx(12,*)
      dimension isize_cart(10),isize_sphe(10)
      data isize_cart /1, 3, 4, 5, 6, 7, 10, 15, 21, 28 /
      data isize_sphe /1, 3, 4, 5, 5, 7,  7,  9, 11, 13 /
c
c>>>> nfunp=0
c>>>> nfunc=0
      do 100 icsh=1,ncsh
        ics=icsh + ncshells
c
        igenc=inx(4 ,ics)
        itype=inx(12,ics)
        if( cart_2_sphe) then
            ifunc=isize_sphe(itype)
        else
            ifunc=isize_cart(itype)
        endif
        inx(11,ics)=nfunc
        nfunc=nfunc+ifunc*(igenc+1)
        inx(10,ics)=nfunc
c contraction lenght :
c from  ia=inx(1,ics)+1
c to    ie=inx(5,ics)
        icont=inx(5,ics)-inx(1,ics)
        nfunp=nfunp+icont*ifunc
  100 continue
c
      ncf=inx(10,ncsh+ncshells)
c
      end
c=================================================================
c2002 subroutine txsarrag(rtdb, nq_init, l_blsize)
      subroutine txsarrag(rtdb, nq_init, l_blsize,int_type)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
#include "mafdecls.fh"
      logical firstd
      integer rtdb
      character*8 int_type
      character*11 scftype
      character*8 where
c----------------------------------------------------------------
      logical cent2,cent3,cent4
      common /what_was_calc/ cent2,cent3,cent4
c----------------------------------------------------------------
      integer h_txs_tempr, k_txs_tempr
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_prep2 , k_txs_prep2 
c----------------------------------------------------------------
      common /pnl_time/ time4texas_hf2_m,time4mul_quart
      common /txs_time/ timprep2,timcalc2,timblok2
c----------------------------------------------------------------
cccc  common /big/ bl(1) repleaced by pnl dbl_mb(ntxs_bl_scr) .
c
      common /bl_txs_add/ ntxs_bl_scr 
c----------------------------------------
      common /ctxs_index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /number/ acc
c------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /intgop/ ncache,maxprice,iprint,iblock
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /runtype/ scftype,where
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c------------------------------------
      common /pnl000/ xbluse,nbluse
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /ilepar/ lpar_exe,lpar_num
c------------------------------------
      common /pnl_nqrt/ ncall_pnl,ncall_41,nqrts_pnl,npart_pnl,nsize_pnl
      common /pnl008/ last_isupb,last_nparts
c------------------------------------
      common /mem_pnl_scr/ nall_peak,mark_peak,mem_peak
c------------------------------------
c 1998 (numerical instability in Tracy's recursive)
c
      common /count_primitive/ x_prim_all,
     *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
      common /lobsa_times/ xlobsa2,xlobsa4
C-----------------------------------------------------------------------
c     common /count_efficient/ xobsa2,xobsa2_ineff,xobsa2_more,
c    *                         xobsa4,xobsa4_ineff,xobsa4_more
c------------------------------------
c initialize commons number,ener,ganz,type :
c
      call txs_initial(rtdb)
c
c------------------------------------
c initialize common /count_primitive/ x_prim_all,
c    *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
c 1998 (numerical instability in Tracy's recursive)
c
c
      x_prim_all=0.d0
      x_prim_wr1=0.d0
      x_prim_wr2=0.d0
      x_prim_wr3=0.d0
c
      xlobsa2=0.d0
      xlobsa4=0.d0
c
c     xobsa2=0.d0
c     xobsa2_ineff=0.d0 
c     xobsa2_more =0.d0 
c     xobsa4=0.d0
c     xobsa4_ineff=0.d0
c     xobsa4_more =0.d0 
c------------------------------------
c initialize common /what_was_calc/ cent2,cent3,cent4
c
      cent2=.false.
      cent3=.false.
      cent4=.false.
c------------------------------------
c initialize common /mem_pnl_scr/ nall_peak,mark_peak,mem_peak
c (needed in texas_terminate for info :
c
      nall_peak=0
      mark_peak=0
      mem_peak=0
c------------------------------------
      ncall_41 =0
      ncall_pnl=0
      nqrts_pnl=0
      npart_pnl=0
      nsize_pnl=0
c
      last_isupb=0
      last_nparts=0
c------------------------------------
      xbluse=0.d0
      nbluse=0
c
      lpar_exe=0
      lpar_num=0
c
c2002 ntotal=0
      noijkl=0
      nopres=0
      nospec=0
c
      nrimtot=0
      nrimret=0
c------------------------------pnl time--------------------------
c zero out timings
c
      time4texas_hf2_m=0.0d0
      time4mul_quart  =0.0d0
c
      timprep2=0.0d0
      timcalc2=0.0d0
      timblok4=0.0d0
c
      call timepr(0.0d0)
c------------------------------pnl time--------------------------
c
      do 20 i=1,4
   20 nsy(i)=0
      nsym=0
c------------------------------------
c only texas correction to the pnl norm :
c
CC--> call txs_cor_norm(bl,inx(inxx),ibas,inuc,ncs,iprint0)
cccc  call txs_cor_norm(dbl_mb,inx(inxx),ibas-1,inuc-1,ncs,iprint0)
      call txs_cor_norm(dbl_mb(ibas),inx(inxx),ibas-1,inuc-1,ncs)
c------------------------------------
c reserve memory for: ncfunct,ncshell and iny :
c
      isize_from_ma = ncf+ncs+12*ncs 
c
      if (.not. MA_alloc_get(MT_DBL,isize_from_ma, 'txs-pnl mapping1',
     &                       h_txs_mapp1, k_txs_mapp1))
     &    call errquit('texas_init: alloc for texas mapping1 failed',
     &    911, MEM_ERR)
      call dfill(isize_from_ma,0.0d00,dbl_mb(k_txs_mapp1),1)
c
      ncfunct=k_txs_mapp1
      ncshell=ncfunct+ncf
      iny    =ncshell+ncs
c------------------------------------
c reorder the basis set according to TEXAS needs :
c
      call reordertxs(ncs,inx(inxx),
     *   dbl_mb(iny),dbl_mb(ncshell),dbl_mb(ncfunct),dbl_mb(inuc),
     *                                               dbl_mb(ibas))
c------------------------------------
c setup common /cpu/ intsize,iacc,icache,memreal
c
c machine characteristics - integer size, cache memory, real fast
c memory (as distinct from virtual memory)
#if defined(KSR)||defined(CRAY)||defined(EXT_INT)
      intsize=1
#else
cold      intsize=2
      intsize=8/ma_sizeof(MT_INT,1,MT_BYTE)
      if(intsize.ne.2) then
         write(6,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         write(6,*) ' 64-bit Architecture, you need to add '
         write(6,*) ' DEFINES+=-DEXT_INT in config/makefile.h'
         write(6,*) ' and recompile'
         write(6,*) ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         call errquit(' txsarrg: compilation error',0,CAPMIS_ERR)
      endif
#endif
c this above should not be used
      iacc=15
c assume 8Kb cache - change it for best efficiency
      icache=1024
c assume 32 M memory
*     memreal=8 388 608  ! RJH ... shurley not, ed.
      memreal=4 194 304
c these variables can be changed
c------------------------------------
c All two-electron integrals program options
c
c setup default values for 2-el.int.prog
c common /intgop/ ncache,maxprice,iprint,iblock
c common /intlim/ limxmem,limblks,limpair
c common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c common /route/ iroute
c
c
c2002 call texas_input
      call texas_input(int_type)
c------------------------------------
c make setup like in Twoint :
c common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c common /runtype/ scftype,where
c threshold is already setup by routine texas_set_accy(pnl_thresh1)
c--
      icache=icache*ncache
      nprint=iprint
      iblok=iblock
      ndirect=0
      nbeg=0
      nend=0
c
c the following variables are not used at all :
c       firstd,ndirect,nbeg,nend
c------------------------------------
c estimate memory requested in texas-scratch bl() by prepint2
c ( prepint2 calls BLOCKIN2 )
c There are 10 mem. allocations :
c 4 proportional to the NCS  : (3+nsymm)*ncs 
c 6 proportional to the NBL2 (no of block-pairs)
c on return size of bl() - ntxs_scr_size
c--------
c allocate memory needed temprorarly (ncs) in txs_scratch_siz1:
c
      if (.not. MA_alloc_get(MT_DBL,ncs+1,'for txs_scr_temp',
     &                       h_txs_tempr , k_txs_tempr ))
     &    call errquit('texas_init: alloc for for txs_scr_temp failed',
     &    911, MEM_ERR)
      call dfill((ncs+1),0.0d00,dbl_mb(k_txs_tempr),1)
c--------
c
      call txs_scratch_siz1(dbl_mb(ibas),dbl_mb(inuc),
     *                      dbl_mb(k_txs_tempr),
     *             iprint,ncs,inx(inxx),ntxs_scr_size)
c--------
c release memory allocated temprorarly in txs_scratch_siz1:
c
        if (.not.MA_free_heap(h_txs_tempr))
     &        call errquit('txs_scr_temp   : failed ?', 911, INT_ERR)
        h_txs_tempr  = -1
        k_txs_tempr  =  0
c--------
c
      lcore=ntxs_scr_size
c
      call retall
c
c------------------------------------
c allocate memory repleacing txs bl() according to the estimate above:
c
      if (.not. MA_alloc_get(MT_DBL,ntxs_scr_size, 'for txs scratch ',
     &                       h_txs_prep2 , k_txs_prep2 ))
     &    call errquit('texas_init: alloc for for txs scratch failed',
     &    911, MEM_ERR)
      call dfill(ntxs_scr_size,0.0d00,dbl_mb(k_txs_prep2),1)
c
      ntxs_bl_scr=k_txs_prep2 
c
c this is starting address in dbl_mb() which repleaces original
c texas scratch bl(). From hereafter texas-bl is eliminated.
c and lcore is now equal to ntxs_scr_size .
c
cccc  write(8,*)' Memory allocated to repleace txs BL :',ntxs_scr_size
c     write(6,*)'1 starting address is ntxs_bl_scr=',ntxs_bl_scr
c
c------------------------------------
c Preparation for the two-electron integrals calculation
c (blocking procedure)
c
      scftype='    '
      call prepint2(dbl_mb(ntxs_bl_scr),
     *              inuc,ibas,na,nbf,nsh,ncf,ncs,inx,lcore,
     *              maxprice,scftype)
c
c------------------------------------
c Check out maximum memory needed so far and reserved in the TXS-BL() :
c
      call memstat(nreque_bl,nmark_bl,maxmem_bl,memtot_bl)
c
c MAXMEM_BL is memory needed for blocking procedure
c
c     IF(IPRINT.gt.0) THEN
c     write(8,*)'-----------------------------------------------'
c     write(8,*)'   memory status after execution of prepint2   '
c     write(8,*)'   -----------------------------------------   '
c     write(8,*)' number of allocations =',nreque_bl
c     write(8,*)' number of marks       =',nmark_bl
c     write(8,*)' maximum memory used   =',maxmem_bl,'(high water)'
c     write(8,*)' total memory available=',memtot_bl
c     write(8,*)'-----------------------------------------------'
c     ENDIF
c
c------------------------------------
c Finally, blocking is done .
c------------------------------------
c MAXMEM_BL is memory needed for blocking procedure.
c To execute texas-integral program more memory is needed, for
c pnl-txs interface (in texas_hf2_m) and some for intermedites.
c The following call estimates this memory AND ADDS it to the
c memory requested for blocking (MAXMEM_BL). On return,
c              L_BLSIZE
c the TOTAL memory needed to call the texas-integral program 
c for NQUART_PNL number of shell quartets :
c
      nquart_pnl=nq_init ! max size of the PNL request (set in api/int_init)
c
cnew  call txs_scratch_siz2(iprint,nquart_pnl,maxmem_bl,l_blsize)
      call txs_scratch_siz2(iprint,ncs,inx(inxx),nquart_pnl,maxmem_bl,
     *                      l_blsize)
c
c------------------------------------
      end
c=================================================================
      subroutine txs_initial(rtdb)
      implicit real*8 (a-h,o-z)
c------------------------------------
c
c     .... set initial values
c
c------------------------------------
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      common /tx_unit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
      common /number/ acc
      common /ganz/ lcore,iov,last,nganz(67)
      common /ctxs_index/ maxsh,ifp,inx(1)
c------------------------------------
c for the case of several calls of texas_init()
c
      data n_open_8 /1/
      save n_open_8
c------------------------------------
c open files for texas
c
c     open(unit=11,status='scratch',form='unformatted')
c     open(unit=12,status='scratch',form='unformatted')
c     open(unit=13,status='scratch',form='unformatted')
c     open(unit=14,status='scratch',form='unformatted')
c     open(unit=15,status='scratch',form='unformatted')
c     open(unit=16,status='scratch',form='unformatted')
c     open(unit=30,file='fort.30',status='unknown',form='formatted')
c     open(unit=31,file='fort.31',status='unknown',form='formatted')
c     open(unit=32,file='fort.32',status='unknown',form='formatted')
c     open(unit=33,file='fort.33',status='unknown',form='formatted')
c------------------------------------
      if(.not.rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,iprint)) then
        iprint = 0
      endif
      if(iprint.gt.0 .and. n_open_8.gt.0) then
        open(unit=8,file='fort.8',status='unknown',form='formatted')
        n_open_8=0
      endif
c------------------------------------
      iov=1
c
      pi=3.1415926535898d0
      ix=5
c-----------
      iout=8
c
c     .... file directory dimensions
c
      ntap=9
      ltab=200
c
c     .... easy to change
c
      ang=1.889726342d0
      debye=0.39342658d0
      cbm=0.117946332d30
      ajoule=0.22936757d0
      evolt=0.036749026d0
      ckalmo=1.5936018d-3
      dkel=3.1667911d-6
      cmm1=4.5563352d-6
      hertz=1.51982983d-16
c
c     this is the required machine accuracy. it should be changed
c     according to the computer used
c
      acc=1.0d-14
c------------------------------------
      do 30 i=1,ntap
         nen(i)=1
         npl(i)=i+10
         nbl(i)=1
         lentry(i)=0
         if (i.gt.3) go to 30
c        np=npl(i)
c        if(np.ne.6) rewind np
   30 continue
c------------------------------------
      ij=0
      do 60 j=1,5
      do 60 i=1,40
         ij=ij+1
         icode(ij)=j*10000+i*100
   60 continue
c------------------------------------
c
      end
c=================================================================
cccc  subroutine txs_cor_norm(bl,inx,ibas,inuc,ncs,iprint)
      subroutine txs_cor_norm(bl,inx,ibas,inuc,ncs)
c bl is now replaced by pnl-array BASNUC .
      implicit real*8 (a-h,o-z)
      double precision pi
cold      dimension bl(*),inx(12,*)
      dimension bl(0:12,*),inx(12,*)
c---------------------------------------------------------------------
C (2/pi)**3/4 :
      pii=acos(-1d0)
      twopi=2d0/pii
      twopi3=twopi*twopi*twopi
      twopi75=sqrt(sqrt(twopi3))
c
      twopi75=1d0/twopi75
c---------------------------------------------------------------------
      do i=1,ncs
         ityp=inx(12,i)
c        begin and end of the contraction
         ig=inx(1,i)+1
         ie=inx(5,i)
         ngc=inx(4,i)
         do igc=1,ngc+1
            if(ig.eq.ie) then
               if(abs(bl(0,ig)).gt.1d-19) then
                  a75=(bl(0,ig)**(-0.75d0))*twopi75
               else
                  a75=twopi75
               endif
c-----------------------------
               bl(igc,ig)=bl(igc,ig)*a75
               if(ityp.eq.3) bl(2,ig)=bl(2,ig)* a75
            else
               do ii=ig,ie
c
c exponent may be zero (only once) to handle 3- and 2-center(index) integrals:
c
                  if(abs(bl(0,ii)).gt.1d-19) then
                     a75=(bl(0,ii)**(-0.75d0))*twopi75
                  else
                     a75=twopi75
                  endif
c-----------------------------
                  bl(igc,ii)=bl(igc,ii)*a75
                  if(ityp.eq.3) bl(2,ii)=bl(2,ii)* a75
c-----------------------------
               enddo
            endif
         enddo
      enddo
c
      end
c=================================================================
c2002 subroutine texas_input
      subroutine texas_input(int_type)
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      character*8 int_type
      logical firstd
      common /intgop/ ncache,maxprice,iprint,iblock
      common /intlim/ limxmem,limblks,limpair
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /mem_max_min/ ispblx,maxmem1,max_111,iforwhat
      common /route/ iroute
      common /allow_lost/ lost_allow   
c----------------------------------------------------------------------
c:added:rak:fix linxmem size to basis set
      integer texas_max_ang_val
      common /c_texas_max_ang_val/ texas_max_ang_val
c----------------------------------------------------------------------
      integer rtdb_copy
      common /c_rtdb_copy/ rtdb_copy
      integer int_tmp, rtdb
      integer rtdblim(3)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "global.fh"
c-----------------------------------------------------
c default values for options in two-el. integ. program 
c change these if needed by reading in new values from 
c some input (PNL must tell) :
c-----------------------------------------------------
c Defaults :
c
c
      rtdb = rtdb_copy
c
c
c
c:: read parameters from run time database use
c:: defaults if they are not set on the rtdb
c
c-----------------------------------------------------
c integral's threshould
c
c... now set by texas_set_accy
*      thresh1 = 1.0d-10
*      if (rtdb_get(rtdb,'int:txs:thre1',MT_DBL,1,dbl_tmp))then
*        thresh1 = dbl_tmp
*        write(luout,*)' texas integral default override: thre1 =',
*     &      thresh1
*      endif
*      thresh2 = 1.0d-10
*      if (rtdb_get(rtdb,'int:txs:thre2',MT_DBL,1,dbl_tmp))then
*        thresh2 = dbl_tmp
*        write(luout,*)' texas integral default override: thre2 =',
*     &      thresh2
*      endif
*
*      par(2)=thresh2
*      par(1)=thresh1
*c *** this value is transferred in par(26)
*      par(26)=thresh1
c-----------------------------------------------------
      icheck = 0
      if (rtdb_get(rtdb,'int:txs:icheck',MT_INT,1,int_tmp)) then
        icheck = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: icheck =',
     &      icheck
      endif
c-----------------------------------------------------
      ncache=32
      if (rtdb_get(rtdb,'int:txs:ncache',MT_INT,1,int_tmp)) then
        ncache = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: ncache =',
     &      ncache
      else if (rtdb_get(rtdb,'int:txs:ncac',MT_INT,1,int_tmp)) then
        ncache = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: ncache =',
     &      ncache
      endif
c-----------------------------------------------------
c
      maxprice=1    ! option name is IPAY
      if (rtdb_get(rtdb,'int:txs:maxprice',MT_INT,1,int_tmp)) then
        maxprice = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: maxprice =',
     &      maxprice
      else if(rtdb_get(rtdb,'int:txs:ipay',MT_INT,1,int_tmp)) then
        maxprice = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: maxprice =',
     &      maxprice
      endif
c-----------------------------------------------------
c
      iprint=0      ! option name is PRIN
      if (rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,int_tmp)) then
        iprint = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: iprint =',
     &      iprint
      else if(rtdb_get(rtdb,'int:txs:prin',MT_INT,1,int_tmp)) then
        iprint = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: iprint =',
     &      iprint
      endif
c-----------------------------------------------------
      iblock=0      ! option name is IBLO  not used
c                   ! option name is LIMI (3 numbers)
c...
c note: texas_max_ang_val is the PNL type for the maximum angular momentum in 
c       a basis set.
c----------------------------------------------------------------------
c                   ! option name is LIMI (3 numbers)
c default setup for limxmem,limblks,limpair :
c
      limxmem=4 000 000 ! limit for memory/block
      if (texas_max_ang_val.eq.6) then
        limxmem = 10 000 000    ! limit for memory/block (i functions)
      else if(texas_max_ang_val.eq.5) then
        limxmem = 6 000 000     ! limit for memory/block (g functions)
      else if(texas_max_ang_val.eq.4) then
        limxmem = 4 000 000     ! limit for memory/block (g functions)
      else if(texas_max_ang_val.eq.3) then
        limxmem = 3 000 000     ! limit for memory/block (f functions)
      else if(texas_max_ang_val.eq.2) then
        limxmem = 2 000 000     ! limit for memory/block (d functions)
      else if(texas_max_ang_val.eq.1) then
        limxmem = 2 000 000     ! limit for memory/block (p functions)
      else if(texas_max_ang_val.eq.0) then
        limxmem = 2 000 000     ! limit for memory/block (s functions)
      else
        write(6,*)
     &      'texas_input: guess for limxmem failed:using:',limxmem
      endif
      if (iforwhat.eq.4) limxmem=10 000 000
c.....................................................................
c2002 limblks=    5 000 ! limit for quartet-block-size
c2002 limpair=    1 000 ! limit for pair-block-size
c
c make these two dependent on a task and more like in pqs :
c
      if(int_type.eq.'scfd_int') then 
        limblks=0  !   limblks=750
        limpair=150
      endif
      if(int_type.eq.'giao_int') then
        limblks=0  !   limblks=500
        limpair=150
      endif
      if(int_type.eq.'der1_int') then 
        limblks=0  !   limblks=300
        limpair=100
      endif
      if(int_type.eq.'der2_int') then 
        limblks=0  !   limblks=150
        limpair=50
      endif
c
c2002 NOTE : for limblks=0 max block size will be determined 
c            automatically in the get_limit routine (spec_block.F)
c----------------------------------------------------------------------
      if (rtdb_get(rtdb,'int:txs:limi',MT_INT,3,rtdblim)) then
        if (rtdblim(1).gt.0) then
          limxmem = rtdblim(1)
          if (ga_nodeid().eq.0)
     &        write(luout,*)
     &        ' texas integral default override: limxmem =',
     &        limxmem
        endif
        if (rtdblim(2).gt.0) then
          limblks = rtdblim(2)
          if (ga_nodeid().eq.0)
     &        write(luout,*)
     &        ' texas integral default override: limblks =',
     &        limblks
        endif
        if (rtdblim(3).gt.0) then
          limpair = rtdblim(3)
          if (ga_nodeid().eq.0)
     &        write(luout,*)
     &        ' texas integral default override: limpair =',
     &        limpair
        endif
      endif
      if (rtdb_get(rtdb,'int:txs:limxmem',MT_INT,1,int_tmp)) then
        limxmem = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: limxmem =',
     &      limxmem
      endif
      if (rtdb_get(rtdb,'int:txs:limblks',MT_INT,1,int_tmp)) then
        limblks = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: limblks =',
     &      limblks
      endif
      if (rtdb_get(rtdb,'int:txs:limpair',MT_INT,1,int_tmp)) then
        limpair = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: limpair =',
     &      limpair
      endif
c-----------------------------------------------------
c default for iroute is 2 because of numerical stability
c
      iroute=0      ! option name is ROUT
c
      if (rtdb_get(rtdb,'int:txs:iroute',MT_INT,1,int_tmp)) then
        iroute = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: iroute =',
     &      iroute
      else if(rtdb_get(rtdb,'int:txs:rout',MT_INT,1,int_tmp)) then
        iroute = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: iroute =',
     &      iroute
      endif
c-----------------------------------------------------
c for numerical instability in texas :
c number of "lost digit" allowed to loose in texas integrals
c i.e. lost_allow=5 means accuracy for integrals is 10**-(14-5)
c
      lost_limit=0
      if(rtdb_get(rtdb,'int:txs:lost_limit',MT_INT,1,int_tmp)) then
        lost_limit = int_tmp
        if (ga_nodeid().eq.0)
     &      write(luout,*)
     &      ' texas integral default override: lost_limit=',lost_limit
      endif
      lost_allow=lost_limit
c
c-----------------------------------------------------
c
      end
c======================================================================
      subroutine txs_scratch_siz1(datbas,datnuc,ishell_blk,
     *                           iprint,ncs,inx,ntxs_scr_size)
      implicit real*8 (a-h,o-z)
c
c estimate memory requested in texas-scratch bl() by prepint2
c
c     common /cpu/ intsize,iacc,icache,memreal
c     common /intlim/ limxmem,limblks,limpair
cccc  common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
      common /ganz/ ll9(9)                           ,na,nn4(4)
      common /route/ iroute
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
      dimension ishell_blk(0:ncs)   ! blocks of shells 
c
c--------------------------------------------------------
c find out what blocking strategy will be used (txs93 or 95)
c
cold  if(iroute.eq.0) call whichblk(datnuc,na,inx,ncs,iroute)
c
      if(iroute.eq.0) then
         call whichblx(na,ncs,inx,iroute,datnuc,datbas,ishell_blk)
      endif
c
c output : iroute
c--------------------------------------------------------
c make blocks of shells and calculate number of bl1,bl2,bl4
c output :
c number of blocks of shells is : nbl1  
c number of blocks of pairs  is : nbl2  
c number of blocks of quarts is : nbl4  
c
      call blk_shells(ncs,inx,iroute,datnuc,datbas,
     *                ishell_blk,nbl1,nbl2,nbl4)
      call get_memory2(iprint,iroute,ncs,nbl1,nbl2,nbl4,memory2)
c--------------------------------------------------------
      ntxs_scr_size= memory2
c-------------------------------------------------------
      end
c======================================================================
cnew  subroutine txs_scratch_siz2(iprint,nquarts,maxmem_bl,l_blsize)
      subroutine txs_scratch_siz2(iprint,ncs,inx,nquarts,maxmem_bl,
     *                            l_blsize)
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------------
c estimate memory requested in texas-scratch bl() by texas_hf2_m
c when it is called for NQUARTS shell quartets
c----------------------------------------------------------------------
c INPUT : iprint - print level
c         nquarts - maximum number of quartets in the PNL request
c         maxmem_bl - memory allocated already in blocking procedure
c
c OUTPUT: l_blsize - total memory needed for texas-integrals
c                    maximum screatch size 
c----------------------------------------------------------------------
c maxme1 - max. memory needed for one block - calculated in blocin2
c maxme2 - max. memory needed in texas_setup (without uniq_pairs)
c maxme3 - max. memory needed in uniq_pairs when called from blockint
c maxme4 - max. memory needed in make_map2uniq & copy_block 
c                                both   called from blockint
c
c maxme2=4*nquartets_pnl + 2*nbl2 + 4*nbl4 + ngcd*size_ijkl
c
c the last term appears only for a request without labels, assume
c 3**4 * (pppp)=> 6561 or 1* ffff => 10000
c----------------------------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      common /mem_max_min/ ispblx,maxmem1,max_111,iforwhat
      common /memor1b/ nbl2
c transfered from blockin2 only to estimate memory here 
      common /max_pairs/ maxpar   
c this is maximum pairs allowed ; must be: maxpar.le.limpair
c
c memory common for use in TEXAS_MEMORY(m_static,nquarts,m_dynam) sub.
c
      common /txs_memory/ mem_stat
c
      dimension inx(12,*)
c----------------------------------------------------------------------
      nbl4=nbl2*(nbl2+1)/2
c------------------------------------------------
c max.memory needed for one block - calculated in blocin2
c
      maxme1=maxmem1
c------------------------------------------------
c max.memory needed in texas_setup (except uniq_pairs):
c
cnew  maxme2=4*nquarts + 2*nbl2 + 4*nbl4 + 10 000
c
      nolab=0
      do ics=1,ncs
         igcon=inx(4, ics) + 1
         isize=inx(3, ics)
         igcsize=igcon*isize
         if(igcsize.gt.nolab) nolab=igcsize
      enddo
      nolab4=nolab**4
c
      maxme2=6*nquarts + 2*nbl2 + 6*nbl4 + nolab4
      maxme2= maxme2/intsize + 11   ! 11 allocations
c
c dynamic memory :
c
      maxme2_dyna=4*nquarts
      maxme2_dyna= maxme2_dyna/intsize + 4
c------------------------------------------------
c memory reserved in uniq_pairs: two allocations
c
      maxme3=2*(maxpar/intsize + 1)
c------------------------------------------------
c memory needed in make_map2uniq: two allocations
c               &  copy_block   : one allocation
c
      maxme4=4*nquarts
      maxme4=maxme4/intsize + 4    ! 3 allocations
c
c dynamic memory :
c
      maxme4_dyna=maxme4
c------------------------------------------------
c final memory estimate :
c
c here:
c
      maxmem_scr = maxme1+maxme2+maxme3+maxme4
c
c total:
c
      maxmem_tot = maxmem_scr + maxmem_bl
c
cadded 10% extra for spin-orbit case failure
      l_blsize=(maxmem_tot*11)/10
c
c minimum memory (1-ij,1-kl,1-ijkl at the time)
c
      minmem_scr = max_111+maxme2+maxme3+maxme4
      minmem_tot = minmem_scr + maxmem_bl
c
c------------------------------------------------
c Separate memory indepent of the size of the PNL request
c
c dynamic (request dependent) :
c
      mem_dyna= maxme2_dyna + maxme4_dyna
c
c static (basis set dependent, request independent) :
c
      mem_stat=maxmem_tot - mem_dyna
c
c store MEM_STAT in a common block for latter use
c when TEXAS_MEMORY(m_static,nquarts,m_dynam) is called
c------------------------------------------------
c
      if(iprint.gt.1) then
        write(8,*)' -------------------------------------------'
        write(8,*)' PREDICTIONS FOR TEXAS MEMORY SCRATCH SIZE :'
        write(8,*)'           (in txs_scratch_siz2)            '
        write(8,*)' -------------------------------------------'
        write(8,*)' memory predicted for blockin4 :',maxme1
        write(8,*)' memory predicted for txs_setup:',maxme2
        write(8,*)' memory predicted for uniq_pair:',maxme3
        write(8,*)' memory predicted for map2uniq :',maxme4
        write(8,*)' MEMORY ESTIMATION FOR SCRATCH2:',maxmem_scr
        write(8,*)' '
        write(8,*)' -------------------------------------------'
        write(8,*)' TOTAL SCRATCH SIZE FOR TEXAS  :',maxmem_tot
        write(8,*)' ...........................................'
        write(8,*)' minimum scratch size for txs  :',minmem_tot
        write(8,*)' minimum memory for 1ij,1kl,1q :',max_111   
        write(8,*)' -------------------------------------------'
c       write(8,*)' memory static (basis set dep.):',mem_stat
c       write(8,*)' memory dynamic(request depen.):',mem_dyna
        write(8,*)' '
      endif
c------------------------------------------------
      end
c===============================================================
      subroutine texas_memory(m_static,nquarts,m_dynam)
      implicit real*8 (a-h,o-z)
      common /memor1b/ nbl2
      common /txs_memory/ mem_stat
c------------------------------------------------
c INPUT : nquarts 
c OUTPUT: m_static - memory independent of the request's size
c OUTPUT: m_dynam  - memory dependent on the request's size
c------------------------------------------------
c dynamic memory :
c
      maxme2_dyna=4*nquarts
      mem_dyna= maxme2_dyna 
c
      m_dynam=mem_dyna
c------------------------------------------------
c static memory :
c
      m_static=mem_stat
c------------------------------------------------
      end
c==================================================================
      subroutine get_int_price(icspnl,jcspnl,kcspnl,lcspnl,price)
c--------------------------------------------------------
c This routine gives back the price of one integral in the 
c PNL quartet of contracted shells
c
c Input  : icspnl,...,lcspnl 4 PNL shells
c Output : price ( relative to the most expensive integral)
c--------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c
      common /bl_txs_add/ ntxs_bl_scr 
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c--------------------------------------------------------
      common /ctxs_index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
c--------------------------------------------------------
      common /max_cost/ most_expensive
c--------------------------------------------------------
      nquart=1
c--------------------------------------------------------
c find texas-shells corresponding to the given set of pnl-shells.
c
      ics = icspnl             ! Copy shells so don't modify arguments
      jcs = jcspnl
      kcs = kcspnl
      lcs = lcspnl
      call txs_pnl(nquart,dbl_mb(ncshell),
     *             icstxs,jcstxs,kcstxs,lcstxs,
     *             ics,jcs,kcs,lcs)
c
c on return shells with texas numbers 
c--------------------------------------------------------
c reorder these shells in quartet according to texas convetion :
c
      call make_txs_quart(ics,jcs,kcs,lcs)
c--------------------------------------------------------
c get absolute price for ONE integral in this quartet 
c
      call get_price1(inx(inxx),ics,jcs,kcs,lcs,nfijkl,nprice)
c
c output : absolute price/int and number of integrals
c          nprice , nfijkl
c--------------------------------------------------------
c get absolute price of the most expensive integral overall:
c
      call get_price1(inx(inxx), 1 , 1 , 1 , 1 ,nf1111,nprmax)
c--------------------------------------------------------
c Finally, price relative to the most expensive integral's price
c
      price=dble(nprice)/dble(nprmax)
c--------------------------------------------------------
      end
c==================================================================
      subroutine add_s_zero(datbas)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------------------
c This routine adds one uncontracted s-orbital to the basis set:
c------------------------------------------------------------------
      common /ctxs_index/ maxsh,ifp,inx(1)    
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension datbas(13,*)
c------------------------------------
      ncs=ncs+1
      nsh=nsh+1
      nbf=nbf+1
      ncf=ncf+1
c
      call add_2_inxbas(inx(inxx), nbf,nsh,ncf,ncs, datbas)
c
c------------------------------------
      end
c==================================================================
      subroutine add_2_inxbas(inx,nbf,nsh,ncf,ncs,datbas)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*)
      dimension datbas(13,*)
c------------------------------------
c
c add to INX :
c
      inx(1,ncs)=nsh-1
      inx(2,ncs)=0        ! center at atom number 0 
      inx(3,ncs)=1        ! shell size 
      inx(4,ncs)=0        ! gen. cont.
      inx(5,ncs)=nsh
c------------------
c     inx(6,ncs)= symm.
c     inx(7,ncs)= symm.
c     inx(8,ncs)= symm.
c     inx(9,ncs)= not used
c------------------
c     inx(11,ncs)= ! this +1 gives the first contr. function
c     inx(10,ncs)= ! this    gives the last  contr. function
c------------------
      inx(10,ncs)=ncf
      inx(11,ncs)=ncf-1
      inx(12,ncs)=1
c------------------------------------
c
c add to DATBAS :
c
      datbas(1,nsh)=0.d0  !exponent
      datbas(2,nsh)=1.d0  !contrac. coef.
ccc   datbas(3,nsh) - (10,nbf) ! gen.contr. coeff.
      datbas(11,nsh)=0.d0   ! x-center
      datbas(12,nsh)=0.d0   ! y-center
      datbas(13,nsh)=0.d0   ! z-center
c
c------------------------------------
      end
c==================================================================
      subroutine print_texas_nbas(nbas)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /intgop/ ncache,maxprice,iprint,iblock
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
c-------------------------------------------------------
      if(iprint.eq.0) RETURN
c-------------------------------------------------------
      write(8,*)' -------------------------------------------'
      write(8,*)' The following basis sets can be used : '
      write(8,*)
c-------------------------------------------------------
      write(8,*)' ...........................................'
      write(8,*)' Basis  number of   primitive    contracted'
      write(8,*)'  Set    centers     shells     func. shells'
      write(8,*)' ...........................................'

      do ibasis=1,nbas
        if(ibasis.eq.1) then
           write(8,80) ibasis,nat_bas_1,
     *                        nps_bas_1,
     *                        ncf_bas_1,
     *                        ncs_bas_1
        endif
        if(ibasis.eq.2) then
           write(8,80) ibasis,nat_bas_2 - nat_bas_1,
     *                        nps_bas_2 - nps_bas_1,
     *                        ncf_bas_2 - ncf_bas_1,
     *                        ncs_bas_2 - ncs_bas_1
        endif
        if(ibasis.eq.3) then
           write(8,80) ibasis,nat_bas_3 - nat_bas_2,
     *                        nps_bas_3 - nps_bas_2,
     *                        ncf_bas_3 - ncf_bas_2,
     *                        ncs_bas_3 - ncs_bas_2
        endif
      enddo
   80 format('   no=',i1,3x,i4,9x,i4,5x,i5,3x,i4)
      write(8,*)' ...........................................'
c-------------------------------------------------------
      if(nbas.gt.1) then
         write(8,*)'                  plus             '
         write(8,*)'     one S-function with zero exponent : '
         write(8,*)'         COMBINED TEXAS BASIS SET  '
         write(8,*)
         write(8,80) nbas+1,na,nsh,ncf,ncs
         write(8,*)' ...........................................'
      endif
c-------------------------------------------------------
      write(8,*)' -------------------------------------------'
c-------------------------------------------------------
      end
c==================================================================
      subroutine print_inx_nuc_bas(nbas,datnuc,datbas)
      implicit real*8 (a-h,o-z)
#include "global.fh"
      common /ctxs_index/ maxsh,ifp,inx(1)    
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
      dimension datnuc(5,*), datbas(13,*)
c
      do ibasis=1,nbas
        if(ibasis.eq.1) then
          if (ga_nodeid().eq.0) then
            write(6,*)
            write(6,*)' PNL Basis Set no = 1'
            write(6,*)' natom=',nat_bas_1,' nsh=',nps_bas_1,
     *          ' ncf=',ncf_bas_1,' ncs=',ncs_bas_1
          endif
          ncs_b=1
          ncs_e=ncs_bas_1
          nps_b=1
          nps_e=nps_bas_1
          nat_b=1
          nat_e=nat_bas_1
        endif
        if(ibasis.eq.2) then
          if (ga_nodeid().eq.0) then
            write(6,*)
            write(6,*)' PNL Basis Set no = 2'
            write(6,*)' natom=',nat_bas_2,' nsh=',nps_bas_2,
     *          ' ncf=',ncf_bas_2,' ncs=',ncs_bas_2
          endif
          ncs_b=ncs_bas_1+1
          ncs_e=ncs_bas_2
          nps_b=nps_bas_1+1
          nps_e=nps_bas_2
          nat_b=nat_bas_1+1
          nat_e=nat_bas_2
        endif
        if(ibasis.eq.3) then
          if (ga_nodeid().eq.0) then
            write(6,*)
            write(6,*)' PNL Basis Set no = 3'
            write(6,*)' natom=',nat_bas_3,' nsh=',nps_bas_3,
     *          ' ncf=',ncf_bas_3,' ncs=',ncs_bas_3
          endif
          ncs_b=ncs_bas_2+1
          ncs_e=ncs_bas_3
          nps_b=nps_bas_2+1
          nps_e=nps_bas_3
          nat_b=nat_bas_2+1
          nat_e=nat_bas_3
        endif
        if (ga_nodeid().eq.0) then
          call do_print(inx(inxx),datnuc,datbas,ncs_b,ncs_e,
     *        nps_b,nps_e,
     *        nat_b,nat_e)
          write(6,*)
        endif
      enddo
c
      end
      subroutine do_print(inx,datnuc,datbas,ncs_b,ncs_e,
     *                                      nps_b,nps_e,
     *                                      nat_b,nat_e)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      dimension datnuc(5,*), datbas(13,*)
      dimension inx(12,*)
c
      write(6,*)' INX array :'
      do ics=ncs_b,ncs_e
      write(6,66) ics, (inx(ii,ics),ii=1,12)
   66 format('ics=',i3,2x,12(i5,1x))
      enddo
      write(6,*)
c
      write(6,*)' DATNUC array :'
      do iat=nat_b,nat_e
      write(6,67) iat, (datnuc(ii,iat),ii=1,5)
   67 format('iat=',i3,2x,5(f10.5,1x))
      enddo
      write(6,*)
c
      write(6,*)' DATBAS array :'
      do ips=nps_b,nps_e
      write(6,68) ips, (datbas(ii,ips),ii=1,13)
   68 format('ips=',i3,2x,5(f10.5,1x)/9x,5(f10.5,1x)/9x,5(f10.5,1x) )
      enddo
      write(6,*)
c
      end
c==================================================================
      subroutine initialized_task(int_type)
      character*8 int_type
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
      common /type_inited/ iforinit  ! use here & in texas_hf.F
c
      where='    '
      iforwhat=0
c
      if(int_type.eq.'scfd_int') then 
         where='buff'
         iforwhat=1
      endif
      if(int_type.eq.'giao_int') then
         where='shif'
         iforwhat=2
      endif
      if(int_type.eq.'der1_int') then 
         where='forc'
         iforwhat=3
      endif
      if(int_type.eq.'der2_int') then 
         where='hess'
         iforwhat=4
      endif
c
      if(iforwhat.eq.0) then
         write(6,*)'From texas_init (requested_task) '
         write(6,*)'int_type is undefined'
         call errquit('texas integral program cannot be initialized',0,
     &       INT_ERR)
      else
ctest
c        where='forc'
c        iforwhat=3
c        write(6,*)'overwritten int_type=',int_type,
c    *             ' where=',where,' iforwhat=',iforwhat
ctest
c
         iforinit=iforwhat
      endif
c
c only variable WHERE is used later on in texas
c IFORWHAT is used only in blksize1 & blksize2
c
      end
c==================================================================
      subroutine get_memory2(iprint,iroute,ncs,nbl1,nbl2,nbl4,
     *                       memory2)
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
ccc   common /memmax/ ispblx, maxme1,iforwhat
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
c-------------------------------------------------------------------
c Memory estimation for prepint2 (blocin2) : 5 memory allocations :
c
c allocated in blkpair :
c
      mem1=ncs+1     ! for nblock1(ncs+1)
      mem2=nbl2      ! for npar(nbl2)
      mem3=nbl4      ! for mxsize(nbl4)
      mem4=nbl4      ! for mxpair(nbl4)
      mem5=ncs       ! for nblock1_back(ncs)
c
c-------------------------------------------------------
      if(intsize.gt.1) then
         mem1 = mem1/intsize +1
         mem2 = mem2/intsize +1
         mem3 = mem3/intsize +1
         mem4 = mem4/intsize +1
         mem5 = mem5/intsize +1
      endif
c-------------------------------------------------------
c total memory :
c
      memory2= mem1+mem2+mem3+mem4+mem5
c-------------------------------------------------------
      if(iprint.gt.0) then
c 
      write(8,*)' -------------------------------------------'
      write(8,*)' TEXAS INTEGRAL PROGRAM HAS BEEN INITIALIZED'
      if(iforwhat.eq.1)
     &    write(8,*)'     for ordinary two-el.  integrals'
      if(iforwhat.eq.2)
     &    write(8,*)'     for first giao/nmr  derivatives'
      if(iforwhat.eq.3)
     &    write(8,*)'     for first  geometry derivatives'
      if(iforwhat.eq.4)
     &    write(8,*)'     for second geometry derivatives'
      write(8,*)' -------------------------------------------'
      write(8,*)' Maximum sizes of blocks upon restrictions :'
      write(8,*)'         cache memory < ',icache
      write(8,*)'         pairs limit  < ',limpair
      if(limblks.eq.0) then
         write(8,*)'         quart limit adjusted'
      else
         write(8,*)'         quart limit  < ',limblks
      endif
      write(8,*)'         memory limit < ',limxmem
      write(8,*)' -------------------------------------------'
c
        if(iroute.eq.1) then 
          write(8,*)' PREDICTIONS FOR TXS-93 BLOCKING PROCEDURE :'
        else
          write(8,*)' PREDICTIONS FOR TXS-95 BLOCKING PROCEDURE :'
        endif
c
c2002   write(8,*)'           (in txs_scratch_siz1)            '
        write(8,*)' -------------------------------------------'
        write(8,*)' number of shell-blocks        :',nbl1 
        write(8,*)' number of pair -blocks        :',nbl2 
        write(8,*)' number of quart-blocks        :',nbl4 
c2002   write(8,*)' MEMORY ESTIMATION FOR SCRATCH1:',memory2
        write(8,*)' -------------------------------------------'
c       write(8,*)' '
      endif
c-------------------------------------------------------
      end
c==================================================================
      subroutine make_txs_quart(icstxs,jcstxs,kcstxs,lcstxs)
c-----------------------------------------------------------------
c Input : icstxs,...,lcstxs - shells from one quartet with texas 
c Output: the same shells in texas order
c-----------------------------------------------------------------
        ics=icstxs
        jcs=jcstxs
        kcs=kcstxs
        lcs=lcstxs
c-----------------------------------------------------------------
        if(ics.lt.jcs) then
          ii=ics
          ics=jcs
          jcs=ii
        endif
        if(kcs.lt.lcs) then
          kk=kcs
          kcs=lcs
          lcs=kk
        endif
c
        ijcs=ics*(ics-1)/2+jcs
        klcs=kcs*(kcs-1)/2+lcs
c
        if(ijcs.lt.klcs) then
           ii=ics
           ics=kcs
           kcs=ii
           jj=jcs
           jcs=lcs
           lcs=jj
           ij=ijcs
c          ijcs=klcs
c          klcs=ij
        endif
c---------------------------------------------------------
        icstxs=ics
        jcstxs=jcs
        kcstxs=kcs
        lcstxs=lcs
c---------------------------------------------------------
      end
c==================================================================
