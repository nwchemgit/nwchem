      subroutine texas_set_accy(pnl_thresh1)
      implicit none
c
      double precision pnl_thresh1 ! [input] accy threshold 1 for texas
c
      call setup_thresh(pnl_thresh1)
      end
      subroutine texas_terminate()
c $Id: texas_face.F,v 1.10 1996-01-12 07:30:05 pg481 Exp $
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_mapp2 , k_txs_mapp2 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_mapp2 , k_txs_mapp2,
     *                      h_txs_prep2 , k_txs_prep2 
c----------------------------------------------------------------
      common /timing2/ timprep2,timstor2,timreca2,timblok2,timblok4
c------------------------------pnl time--------------------------
      common /pnl_time1/ time_pt,time_setup
      common /pnl_time2/ time_calc2,time_do2,time_do4,time_doqrt
      common /pnl_time3/ time_pr2sb,time_pr2sm
c
      common /pnl_nqrt/ ncall_pnl, nqrts_pnl
c------------------------------pnl time--------------------------
      common /pnl000/ xbluse,nbluse
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /ilepar/ lpartot,lpareal
c----------------------------------------------------------------
c efficiency of using block's suructure :
c
      if (nbluse.eq.0) then
        bieffiy = 0.0d00
      else
        ysmallbl=dble(nbluse)
        bleffiy=(xbluse/ysmallbl)*100.d0
      endif
c----------------------------------------------------------------
c Information about integral's program performance :
c
      write(8,*)'------------------------------------------------------'
      write(8,*)'Texas-int2 program has been called ',ncall_pnl,' times'
      if(ncall_pnl.eq.0) goto 00001
      write(8,*)'number of processed shell quartets ',nqrts_pnl
      write(8,18) bleffiy
   18 format(' blocking utilities of Texas used in ',f6.2,' %')
c----------------------------------------------------------------
c-print info :
      write(8,*)'------------------------------------------------------'
      write(8,*) 'Pre-calculations for pairs of contracted shells have'
      write(8,*) 'been executed ',lpareal,'(out of ',lpartot,') times '
      write(8,*)'------------------------------------------------------'
      write(8,*)'Total number of contracted shell quartets ',ntotal
      write(8,*)'has been reduced by molecular symmetry to ',noijkl
      write(8,*)'has been reduced by special request    to ',nospec
      write(8,*)'and due to the screening procedure to     ',nopres
      write(8,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - - '
      write(8,*)'number of primitive quartets after symmet.',nrimtot
      write(8,*)'reduced by the screening procedure to     ',nrimret
      write(8,*)'------------------------------------------------------'
c
c----------------------------------------------------------------
c CPU timing :
c
      write(8,*)'----------------------------------------'
      write(8,*)'CPU time of using TEXAS-integral program '
      write(8,*)'no of calls =',ncall_pnl,' no of quartets =',nqrts_pnl
      write(8,*)'----------------------------------------'
c
      tlocate= time_do2+time_do4+time_doqrt
      ttransf= time_setup+time_pt
      tpurcal=timreca2-tlocate + timprep2
c
      write(8,140) tpurcal,tlocate,ttransf,time_calc2+timprep2
  140 format('texas pure-calculations time  =',f10.2/
     *       'block & quartet locate  time  =',f10.2/
     *       'pnl-txs-pnl transfering time  =',f10.2/
     *       'Total integral delivery time  =',f10.2)
      write(8,*)'----------------------------------------'
c
      write(8,141) time_do2,time_do4,time_doqrt
  141 format('     doblock2 time =',f10.2/
     *       '     doblock4 time =',f10.2/
     *       '     doquarts time =',f10.2)
      write(8,*)'----------------------------------------'
      write(8,142) time_setup,time_pt
  142 format('     txs_setup time=',f10.2/   
     *       '     pnl_label time=',f10.2)
      write(8,*)'----------------------------------------'
c
      write(8,145)time_pr2sb,time_pr2sm
  145 format('split prec2  : sup-bl time=',f10.2,' small-bl=',f10.2/)
c
      call timepr(timblok2+timblok4)
c
c------------------------------pnl time--------------------------
c------------------------------pnl performance ------------------
c
          call txs_for_pnl( dbl_mb(k_txs_prep2) )
c
c----------------------------------------------------------------
c free up memory allocated by texas_init      
c
00001 continue
*rak:      write(6,*)'t h_txs_basnuc ',h_txs_basnuc 
*rak:      write(6,*)'t h_txs_mapp1  ',h_txs_mapp1
*rak:      write(6,*)'t h_txs_mapp2  ',h_txs_mapp2  
      if (.not.MA_free_heap(h_txs_basnuc))
     &        call errquit('texas_terminate: failed ?', 911)
      if (.not.MA_free_heap(h_txs_mapp1 ))
     &        call errquit('texas_terminate: failed ?', 911)
      if (.not.MA_free_heap(h_txs_mapp2 ))
     &        call errquit('texas_terminate: failed ?', 911)
      if (.not.MA_free_heap(h_txs_prep2 ))
     &        call errquit('texas_terminate: failed ?', 911)

c
      h_txs_basnuc = -1
      h_txs_mapp1  = -1 
      h_txs_mapp2  = -1
      h_txs_prep2  = -1
      k_txs_basnuc =  0
      k_txs_mapp1  =  0
      k_txs_mapp2  =  0
      k_txs_prep2  =  0
c
      end
c
c End of terminate 
c================================================================
      subroutine txs_for_pnl(bl)
c for performance estimation :
      implicit real*8 (a-h,o-z)
      common /pnlRJH/ isblcall
c------
      common /memor1a_R/ npard_R,ncost_R,nsupb_R, mxsize_R,nasize_R
c------
      common /memor1b/ nbl2,nbloks
c
      dimension bl(*)
c
         call rob_harr(bl(nasize_R),bl(isblcall),nbl2 )
c
      end
c---------------
      subroutine rob_harr(nasize,isbl_call,nbl2)
      dimension nasize(*),isbl_call(*)
c
      write(8,*) '===================================================='
      write(8,*) ' super-block  texas_size  pnl_calls average_pnl_size'
      write(8,*) '----------------------------------------------------'
c
      nbl4=nbl2*(nbl2+1)/2
c
      iaver_size_txs=0
      iaver_size_pnl=0
      nsbl_pnl=0
      do 10 isbl=1,nbl4
         isbl_called=isbl_call(isbl)
         if(isbl_called.eq.0) go to 10
c
         isbl_txs=nasize(isbl)
         isbl_pnl=isbl_txs/isbl_called
c
         iaver_size_txs=iaver_size_txs+isbl_txs
         iaver_size_pnl=iaver_size_pnl+isbl_pnl
c
         nsbl_pnl=nsbl_pnl+isbl_called
c
         write(8,66) isbl,isbl_txs,isbl_called, isbl_pnl
   10 continue
c
c average values :
c
c  it is like for pnl we would have more super-blocks :
c
         iaver_size_pn2=iaver_size_txs/nsbl_pnl 
c
         iaver_size_txs=iaver_size_txs/nbl4
         iaver_size_pnl=iaver_size_pnl/nbl4     
c
   66 format(i7,5x,3(i7,5x))
c
      write(8,*) '===================================================='
      write(8,*) '=========== average size of super-block ============'
      write(8,*) '                                                    '
      write(8,*) '--- texas ------------- pnl ------------- pnl2 -----'
      write(8,88) iaver_size_txs, iaver_size_pnl, iaver_size_pn2
   88 format(3(2x,i7,10x)) 
      write(8,*) '===================================================='
      end
c=================================================================
      subroutine texas_init(rtdb,nbas,bases,l_blsize)
c---------------------------------------------------------------
c l_blsize = size of the scratch array which will be used when
c             texas_hf2_m is called
c this is maximum size of the scratch needed to do one super-block
c---------------------------------------------------------------
      implicit none
c
c..  this routine is the current interface to the texas code
c
#include "bas.fh"
#include "geom.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "apiP.fh"
c::functions
c::passed
      integer rtdb ! rtdb handle
      integer nbas ! number of basis sets 
      integer bases(nbas) ! basis set handles
      integer l_blsize
c::local
      integer geom, basis 
      integer type, i, j, k, nprim, icoef, nshells, ishell, ngen
      integer sphcart, kend
      integer nat ! number of atoms
      integer center ! center for a given contr. shell
      integer npshells ! total number of primitive shells
      integer ipshell 
      integer size_from_ma
      integer inuc,ibas
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_mapp2 , k_txs_mapp2 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_mapp2 , k_txs_mapp2,
     *                      h_txs_prep2 , k_txs_prep2 
c----------------------------------------------------------------------
      integer basis_init
      common /c_basis_init/ basis_init
c----------------------------------------------------------------------
c
c only for tests :
c     integer ish,jsh,ksh,lsh,leri
c     integer ijsh,klsh,ijkl
c     double precision eri(614 656)  ! up to (ii/ii)
c     double precision ra(3),rb(3),rc(3),rd(3)
c only for tests - end
c----------------------------------------------------------------------
      integer rtdb_copy
      common /c_rtdb_copy/ rtdb_copy
      character*255 basis_name,trans_bname
      integer bas
c----------------------------------------------------------------------
c only for tests:
c
c     leri=614656
c----------------------------------------------------------------------
c
      rtdb_copy = rtdb
      if (nbas.gt.1) then
        write(6,*)' warning '
        write(6,*)
     &  ' texas can only handle 1 basis today/tommorrow/and Tuesday'
        write(6,*)' using the first basis passed as THE basis.'
        write(6,*)' '
        write(6,*)' Krys is working on the fix '
      endif
c
      basis = bases(1) ! only one basis in texas
      basis_init = basis
c
      bas = basis + BASIS_HANDLE_OFFSET
      geom = ibs_geom(bas)
c
      nat = ncenter(geom)
c------------------------------------
c total number of contracted shells
c
      nshells = ncont_tot_gb(bas)
c------------------------------------
c total number of primitive shells
c
      npshells = nprim_tot_gb(bas)
c------------------------------------
c For 3- and 2-center(index) integrals add one 
c uncontracted S-shell with exponent ZERO :
c If these integrals are going to be calculated latter
c number of basis SETS will be HERE greater than 1 :
c
      if(nbas.gt.1) then
      endif
c
c------------------------------------
c Put Basis set and nuclear data into the PNL location
c (remove it from bl() ) :
c 5*natoms
c 13*npshells
c
      size_from_ma = 5*nat + 13*npshells
c
      if (.not. MA_alloc_get(MT_DBL, size_from_ma, 'texas basis/nuc ',
     &                       h_txs_basnuc, k_txs_basnuc))
     &    call errquit('texas_init: alloc for texas basis/nuc failed',
     &    911)
      inuc = k_txs_basnuc
      ibas = inuc + 5*nat
c------------------------------------
c make bl and inx arrays 
c
cOLD  call blinx(nat,nshells,npshells)
      call blinx(nat,nshells,npshells,inuc,ibas)
c------------------------------------
c make datnuc array
c
      call makenuc(geom,nat,dbl_mb(inuc))
c------------------------------------
c
      call makebas(nshells,bas, dbl_mb(ibas))
c------------------------------------
c finish inx array
c
      call makeinx
c------------------------------------
c Re-arrange a basis set according to 
c the Texas convention
c 1) normalize and 2) re-order
c
      call txsarrag(dbl_mb, l_blsize)
c
c... map texas memory to memory commons of the PNL/int API
c
      write(6,*)' mem_2e4c :',mem_2e4c
      write(6,*)' mem_2e3c :',mem_2e3c
      write(6,*)' mem_2e2c :',mem_2e2c
      write(6,*)' l_blsize :',l_blsize
      mem_2e4c = max(mem_2e4c,l_blsize)
      mem_2e3c = max(mem_2e3c,l_blsize)
      mem_2e2c = max(mem_2e2c,l_blsize)
c
c------------------------------------
c ONLY for TESTS :
c  calculate two-electron integrals  
c  in one-quartet pices :
c
c---
c      call test_calc2(nshells)
c---
c End of tests . (program will be stoped in test_calc2) 
c------------------------------------
*rak:      write(6,*)'i h_txs_basnuc ',h_txs_basnuc
*rak:      write(6,*)'i h_txs_mapp1  ',h_txs_mapp1
*rak:      write(6,*)'i h_txs_mapp2  ',h_txs_mapp2
      end
c=================================================================
      subroutine blinx(natom,ncshell,npshell,inucx,ibasx)
      implicit real*8 (a-h,o-z)
c------------------------------------
c the common block /big/ bl(lcx) has been 
c repleaced by pnl allocated scratch 
c------------------------------------
      parameter(mxsh= 2 000 000)
c------------------------------------
      common /index/ maxsh,ifp,inx(mxsh)
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
      subroutine makenuc(geom,nat,datnuc)
      implicit real*8 (a-h,o-z)
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
      integer geom, nat
      double precision datnuc(5,nat)
c
      do 100 iat=1,nat
      datnuc(1,iat)=charge(iat,geom)
      datnuc(2,iat)=coords(1,iat,geom)
      datnuc(3,iat)=coords(2,iat,geom)
      datnuc(4,iat)=coords(3,iat,geom)
      datnuc(5,iat)=iat
cccc  write(6,66) datnuc(5,iat),(datnuc(ii,iat),ii=1,4)
  100 continue
   66 format('atom Name=',i3,' charge=',f5.1,' xyz(au)=',3f10.5)
c
      end
c=================================================================
      subroutine makebas(nshells,bas,datbas)
      implicit real*8 (a-h,o-z)
      integer bas
      common /index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      dimension datbas(13,*)
c
c------------------------------------------------
c shell's type :
c
c ityp1= 1   2   3   4   5   6   7   8   9   10    !  texas type 
c        s   p   l   d5  d6  f7  f10 g15 h21 i28
c itype= 0   1   -1      2       3   4   5   6     !  pnl type
c------------------------------------------------
c
      ipshell = 0
      do 00100 ics=1,nshells
ccccc    call basdat(bas, ics, bl(ibas+1),inx(inxx), ipshell)
         call basdat(bas, ics, datbas    ,inx(inxx), ipshell)
00100 continue
c
      end
c==============1
      subroutine basdat(bas, ics, datbas,inx, ipshell)
      implicit real*8 (a-h,o-z)
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
      integer bas
      dimension datbas(13,*),inx(12,*)
      dimension isize(10)
      data isize /1, 3, 4, 5, 6, 7, 10, 15, 21, 28 /
c---------------  s  p  l d5 d6  f7 f10 g15 h21 i28 -----------c
c
      inx(1,ics)=ipshell
      inx(2,ics)=ibs_cn2ce(ics,bas)
      itype = infbs_cont(CONT_TYPE,ibs_cn2ucn(ics,bas),bas)
      ityp1 = itype + 1
      if(itype.eq.-1) ityp1=3
      if(itype.eq.2 ) ityp1=5
      if(itype.gt.2 ) ityp1=itype+4
      inx(3,ics)=isize(ityp1)
      igen = infbs_cont(CONT_NGEN,ibs_cn2ucn(ics,bas),bas)
      inx(4,ics)=igen-1
      if(ityp1.eq.3) inx(4,ics)=0
      iprim = infbs_cont(CONT_NPRIM,ibs_cn2ucn(ics,bas),bas)
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
      iexppnl = infbs_cont(CONT_IEXP,ibs_cn2ucn(ics,bas),bas) - 1
      icofpnl = infbs_cont(CONT_ICFP,ibs_cn2ucn(ics,bas),bas) - 1
      do 100 ipr=1,iprim
        ips=ipr+ipshell
*       datbas(1,ips)=ext(ipr)
        datbas(1,ips)=exndcf(iexppnl+ipr,bas)
        if(ityp1.eq.3) then
           datbas(2,ips)=exndcf(icofpnl+ipr,bas)
           datbas(3,ips)=exndcf(icofpnl+ipr+iprim,bas)
        else
           do 150 igc=1,igen
            datbas(1+igc,ips)=exndcf(icofpnl+ipr+((igc-1)*iprim),bas)
 150       continue
        endif
        iat = ibs_cn2ce(ics,bas)
        datbas(11,ips)=coords(1,iat,ibs_geom(bas))
        datbas(12,ips)=coords(2,iat,ibs_geom(bas))
        datbas(13,ips)=coords(3,iat,ibs_geom(bas))
 100  continue
c
      ipshell = ipshell+iprim
      end
c=================================================================
      subroutine makeinx
      implicit real*8 (a-h,o-z)
      common /index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c
      call inxdat(ncs,inx(inxx),nbf,ncf)
c
      end
c==============1
      subroutine inxdat(ncs,inx,nbf,ncf)
      dimension inx(12,*)
      dimension isize(10)
      data isize /1, 3, 4, 5, 6, 7, 10, 15, 21, 28 /
c
      nfunc=0
      do 100 ics=1,ncs
        igenc=inx(4 ,ics)
        itype=inx(12,ics)
        ifunc=isize(itype)
        inx(11,ics)=nfunc
        nfunc=nfunc+ifunc*(igenc+1)
        inx(10,ics)=nfunc
  100 continue
c
      ncf=inx(10,ncs)
      nbf=inx(5, ncs)
c
      end
c=================================================================
      subroutine txsarrag(basnuc, l_blsize)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      character*11 scftype
      character*4 where
c----------------------------------------------------------------
      integer h_txs_basnuc, k_txs_basnuc
      integer h_txs_mapp1 , k_txs_mapp1 
      integer h_txs_mapp2 , k_txs_mapp2 
      integer h_txs_prep2 , k_txs_prep2 
      common /txs_ma_stuff/ h_txs_basnuc, k_txs_basnuc,
     *                      h_txs_mapp1 , k_txs_mapp1,
     *                      h_txs_mapp2 , k_txs_mapp2,
     *                      h_txs_prep2 , k_txs_prep2 
c----------------------------------------------------------------
c------------------------------pnl time--------------------------
      common /pnl_time1/ time_pt,time_setup
      common /pnl_time2/ time_calc2,time_do2,time_do4,time_doqrt
      common /pnl_time3/ time_pr2sb,time_pr2sm
c------------------------------pnl time--------------------------
cccc  common /big/ bl(1)
c repleaced by pnl dbl_mb(ntxs_bl_scr) .
c
      common /bl_txs_add/ ntxs_bl_scr 
c----------------------------------------
      common /index/ maxsh,ifp,inx(1)
      common /inde1/ inxx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /intgop/ ncache,maxprice,iprint,iblock
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /runtype/ scftype,where
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c------------------------------------
      common /pnl000/ xbluse,nbluse
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /ilepar/ lpartot,lpareal
c------------------------------------
      dimension basnuc(*)
c------------------------------------
c initialize commons number,ener,ganz,type :
c
      call txs_initial
c
c------------------------------------
c
      xbluse=0.d0
      nbluse=0
c
      lpartot=0
      lpareal=0
c
      ntotal=0
      noijkl=0
      nopres=0
      nospec=0
c
      nrimtot=0
      nrimret=0
c------------------------------pnl time--------------------------
      time_pt=0.d0
      time_setup=0.d0
c
      time_calc2=0.d0
      time_doqrt=0.d0
      time_do2  =0.d0
      time_do4  =0.d0
      time_pr2sb=0.d0
      time_pr2sm=0.d0
      call timepr(zero)
c------------------------------pnl time--------------------------
c
      do 20 i=1,4
   20 nsy(i)=0
      nsym=0
c------------------------------------
c only texas correction to the pnl norm :
c
      iprint0=0
CC--> call txs_cor_norm(bl,inx(inxx),ibas,inuc,ncs,iprint0)
      call txs_cor_norm(basnuc,inx(inxx),ibas-1,inuc-1,ncs,iprint0)
c------------------------------------
c reserve memory for: ncfunct,ncshell and iny :
c
      isize_from_ma = ncf+ncs+12*ncs 
c
      if (.not. MA_alloc_get(MT_DBL,isize_from_ma, 'txs-pnl mapping1',
     &                       h_txs_mapp1, k_txs_mapp1))
     &    call errquit('texas_init: alloc for texas basis/nuc failed',
     &    911)
c
      ncfunct=k_txs_mapp1
      ncshell=ncfunct+ncf
      iny    =ncshell+ncs
c------------------------------------
c reorder the basis set according to TEXAS needs :
c
      call reorder(ncs,inx(inxx),
     *     dbl_mb(iny),dbl_mb(ncshell),dbl_mb(ncfunct) )
c------------------------------------
C Old memory handling :
c Re-order a basis set according to Texas requirement :
c (remember the originally ordered basis set under bl(ncshell) )
c
c     call getmem(ncf,ncfunct)    ! cont.funct. : txs-->pnl
c     call getmem(ncs,ncshell)    ! cont.shells:  pnl-->txs
c     call getmem(ncs*12,iny)
c
c     call reorder(ncs,inx(inxx),bl(iny),bl(ncshell),bl(ncfunct) )
c
c     call retmem(1)
c------------------------------------
c Reserve memory for an array showing to which pair-block
c a given pair of shells belongs:
c
c     call getmem(ncs*(ncs+1)/2 , nblock2 )
c------
      jsize_from_ma = ncs*(ncs+1)/2
c
      if (.not. MA_alloc_get(MT_DBL,jsize_from_ma, 'txs-pnl mapping2',
     &                       h_txs_mapp2 , k_txs_mapp2 ))
     &    call errquit('texas_init: alloc for texas-mapping 2 failed',
     &    911)
c
      nblock2=k_txs_mapp2 
c------------------------------------
c setup common /cpu/ intsize,iacc,icache,memreal
c
c machine characteristics - integer size, cache memory, real fast
c memory (as distinct from virtual memory)
      intsize=2
c this above should not be used
      iacc=15
c assume 64K cache - change it for best efficiency
      icache=8192
c assume 32 M memory
      memreal=8 388 608
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
      call texas_input
c------------------------------------
c estimate memory requested in texas-scratch bl() by prepint2
c ( prepint2 calls BLOCKIN2 )
c There are 10 mem. allocations :
c 4 proportional to the NCS  : (3+nsymm)*ncs 
c 6 proportional to the NBL2 (no of block-pairs)
c on return size of bl() - ntxs_scr_size
c
      call txs_scratch_siz1(ncs,inx(inxx),nsym,ntxs_scr_size)
c
      lcore=ntxs_scr_size
c
      call retall
c
c------------------------------------
c allocate memory repleacing txs bl():
c
c
      if (.not. MA_alloc_get(MT_DBL,ntxs_scr_size, 'for txs scratch ',
     &                       h_txs_prep2 , k_txs_prep2 ))
     &    call errquit('texas_init: alloc for for txs scratch failed',
     &    911)
c
      ntxs_bl_scr=k_txs_prep2 
c
c this is starting address in dbl_mb() which repleaces original
c texas scratch bl(). Since hereafter texas-bl is eliminated.
c and lcore is now equal to ntxs_scr_size .
c
      write(8,*)' Memory allocated to repleace txs BL :',ntxs_scr_size
c
c------------------------------------
c make setup like in Twoint :
c common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
c common /runtype/ scftype,where
c
c     eps=par(1)
c threshold is already setup by subroutine texas_set_accy(pnl_thresh1)
c--
      eps=0.d0   
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
c Preparation for the two-electron integrals calculation
c (blocking procedure)
c
      scftype='    '
cold-
c     call prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,lcore,nsym,
c    *              maxprice,scftype)
cold-
c
      call prepint2(dbl_mb(ntxs_bl_scr),
     *              eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,lcore,nsym,
     *              maxprice,scftype)
c
cccc  write(6,100) scftype
c 100 format(
c    */'------------------------------------------------------'/
c    * '*  SCF will be performed in the ',a11,' fashion  *'/
c    * '------------------------------------------------------'/)
c
c------------------------------------
c Construct an array showing to which pair-block
c a given pair of shells belongs:
c
       call pair_in2(dbl_mb(ntxs_bl_scr),dbl_mb)
c
c on return nblock2(ijcs)--> ibl 
c placed in bl(nblock2)
c------------------------------------
c Check out maximum memory needed so far and reserved
c in the TXS-BL() :
c
      call memstat(nreque_bl,nmark_bl,maxmem_bl,memtot_bl)
c
      write(8,*)'-----------------------------------------------'
      write(8,*)' At the end of TXSARRAG in texas scratch BL() :'
      write(8,*)'                                              :'
      write(8,*)' number of allocations =',nreque_bl
      write(8,*)' number of marks       =',nmark_bl
      write(8,*)' maximum memory overall=',maxmem_bl
      write(8,*)' total memory available=',memtot_bl
      write(8,*)'-----------------------------------------------'
c
c------------------------------------
c estimate memory requested in texas-scratch bl() by texas_hf2_m :
c the size of the scratch bl_scr needed to call texas_hf2_m with :
c
      call txs_scratch_siz2(maxme1,maxme2)
c
      l_blsize=maxme1 + maxme2 + maxmem_bl
c
      write(8,*)' maximum size of bl_scr=',l_blsize  
      write(8,*)'-----------------------------------------------'
c
c maxme1 known from prepint2 execution (common /memmax/ )
c maxme2 estimated as needed in texas_setup
c------------------------------------
      end
c=================================================================
      subroutine txs_initial
      implicit real*8 (a-h,o-z)
c------------------------------------
c
c     .... set initial values
c
c------------------------------------
      common /tape/ inp,inp2,iout,ipun,ix,icond,itest,nentry,ltab,ntap,n
     1pl(9),nbl(9),nen(9),lentry(9),nam(200),num(200),irec(200),icode(20
     20),inpf,ioutf
      common /unit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /ganz/ lcore,iov,last,nganz(67)
      common /index/ maxsh,ifp,inx(1)
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
c-------------
      open(unit=8,file='fort.8',status='unknown',form='formatted')
ccccc open(unit=7,file='fort.7',status='unknown',form='formatted')
c--------------------------------
c ONLY FOR TESTS with test_calc2
      open(unit=81,file='integ.81',status='unknown',form='formatted')
      open(unit=82,file='integ.82',status='unknown',form='formatted')
c------------------------------------
      iov=1
c
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      five=5.0d0
      ten=10.0d0
      ten6=1.0d+6
       tenm8=1.0d-8
COLD  pi=3.1415926535d0
      pi=3.1415926535898d0
      ix=5
c-----------
      iout=8
c     ipun=7
c     inp=30
c     inp2=31
c     icond=6
c     inpf=32
c     ioutf=33
c
c     .... file directory dimensions
c
      ntap=9
      ltab=200
c
c     .... easy to change
c
      time=ten6
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
c     .... this is the required machine accuracy. it should be changed
c     according to the computer used
c
      acc=1.0d-14
c------------------------------------
c     do 10 i=1,4
c        nganz(i)=0
c  10 continue
c     do 20 i=1,67
c  20 nganz(i)=0
c    
c     nganz(18)=3070
c
c standard integral block-size
c
c------------------------------------
c     rewind inp
c     rewind inp2
c------------------------------------
      do 30 i=1,ntap
         nen(i)=1
         npl(i)=i+10
         nbl(i)=1
         lentry(i)=0
         if (i.gt.3) go to 30
         np=npl(i)
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
      subroutine txs_cor_norm(bl,inx,ibas,inuc,ncs,iprint)
c bl is now replaced by pnl-array BASNUC .
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension bl(*),inx(12,*)
c---------------------------------------------------------------------
      data pii /3.1415926535898d0/   
C (2/pi)**3/4 :
      twopi=two/pii
      twopi3=twopi*twopi*twopi
      twopi75=sqrt(sqrt(twopi3))
c
      twopi75=one/twopi75
c---------------------------------------------------------------------
      do 290 i=1,ncs
         ityp=inx(12,i)
c        begin and end of the contraction
         ig=inx(1,i)+1
         ie=inx(5,i)
         ngc=inx(4,i)
         do 285 igc=1,ngc+1
         do 280 ii=ig,ie
            ja=ibas+13*ii-12
            i1=inx(2,i)
            i2=inuc+i1*5
  270       a=bl(ja)
c
c exponent may be zero (only once) to handel 3- and 2-center(index) integrals:
c
            if(a.gt.zero) then
               a75=a**(-0.75d0)
            else
               a75=one
            endif
c-----------------------------
            bl(ja+igc)=bl(ja+igc)*twopi75 * a75
            if(ityp.eq.3) bl(ja+2)=bl(ja+2)*twopi75 * a75
c-----------------------------
  280    continue
c       loop over general contraction
  285    continue
  290 continue
c
      end
c=================================================================
      subroutine texas_input
      implicit real*8 (a-h,o-z)
      implicit integer (i-n)
      logical firstd
      common /intgop/ ncache,maxprice,iprint,iblock
      common /intlim/ limxmem,limblks,limpair
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /route/ iroute
      integer rtdb_copy
      common /c_rtdb_copy/ rtdb_copy
      integer int_tmp, rtdb
      integer rtdblim(3)
#include "mafdecls.fh"
#include "rtdb.fh"
#include "stdio.fh"
c-----------------------------------------------------
c default values for options in two-el. integ. program 
c change these if needed by reading in new values from 
c some input (PNL must tell) :
c-----------------------------------------------------
c Defaults :
c------------------------------------
c
      rtdb = rtdb_copy
c:: read parameters from run time database use
c:: defaults if they are not set on the rtdb
c
c------------------------------------
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
c------------------------------------
      icheck = 0
      if (rtdb_get(rtdb,'int:txs:icheck',MT_INT,1,int_tmp)) then
        icheck = int_tmp
        write(luout,*)' texas integral default override: icheck =',
     &      icheck
      endif
      ncache=10     ! option name is NCAC
      if (rtdb_get(rtdb,'int:txs:ncache',MT_INT,1,int_tmp)) then
        ncache = int_tmp
        write(luout,*)' texas integral default override: ncache =',
     &      ncache
      else if (rtdb_get(rtdb,'int:txs:ncac',MT_INT,1,int_tmp)) then
        ncache = int_tmp
        write(luout,*)' texas integral default override: ncache =',
     &      ncache
      endif
c
      maxprice=1    ! option name is IPAY
      if (rtdb_get(rtdb,'int:txs:maxprice',MT_INT,1,int_tmp)) then
        maxprice = int_tmp
        write(luout,*)' texas integral default override: maxprice =',
     &      maxprice
      else if(rtdb_get(rtdb,'int:txs:ipay',MT_INT,1,int_tmp)) then
        maxprice = int_tmp
        write(luout,*)' texas integral default override: maxprice =',
     &      maxprice
      endif
c
      iprint=0      ! option name is PRIN
      if (rtdb_get(rtdb,'int:txs:iprint',MT_INT,1,int_tmp)) then
        iprint = int_tmp
        write(luout,*)' texas integral default override: iprint =',
     &      iprint
      else if(rtdb_get(rtdb,'int:txs:prin',MT_INT,1,int_tmp)) then
        iprint = int_tmp
        write(luout,*)' texas integral default override: iprint =',
     &      iprint
      endif
      iblock=0      ! option name is IBLO  not used
c                   ! option name is LIMI (3 numbers)
c...
      limxmem=500 000 ! limit for memory/block
      limblks=  5 000 ! limit for quartet-block-size
      limpair=  1 000 ! limit for pair-block-size
c
      if (rtdb_get(rtdb,'int:txs:limi',MT_INT,3,rtdblim)) then
        if (rtdblim(1).gt.0) then
          limxmem = rtdblim(1)
          write(luout,*)' texas integral default override: limxmem =',
     &        limxmem
        endif
        if (rtdblim(2).gt.0) then
          limblks = rtdblim(2)
          write(luout,*)' texas integral default override: limblks =',
     &        limblks
        endif
        if (rtdblim(3).gt.0) then
          limpair = rtdblim(3)
          write(luout,*)' texas integral default override: limpair =',
     &        limpair
        endif
      endif
      if (rtdb_get(rtdb,'int:txs:limxmem',MT_INT,1,int_tmp)) then
        limxmem = int_tmp
          write(luout,*)' texas integral default override: limxmem =',
     &        limxmem
      endif
      if (rtdb_get(rtdb,'int:txs:limblks',MT_INT,1,int_tmp)) then
        limblks = int_tmp
          write(luout,*)' texas integral default override: limblks =',
     &        limblks
      endif
      if (rtdb_get(rtdb,'int:txs:limpair',MT_INT,1,int_tmp)) then
        limpair = int_tmp
          write(luout,*)' texas integral default override: limpair =',
     &        limpair
      endif
c
      iroute=0      ! option name is ROUT
      if (rtdb_get(rtdb,'int:txs:iroute',MT_INT,1,int_tmp)) then
        iroute = int_tmp
        write(luout,*)' texas integral default override: iroute =',
     &      iroute
      else if(rtdb_get(rtdb,'int:txs:rout',MT_INT,1,int_tmp)) then
        iroute = int_tmp
        write(luout,*)' texas integral default override: iroute =',
     &      iroute
      endif
c-----------------------------------------------------
c
      end
c======================================================================
      subroutine pair_in2(bl,bl_pnl)
      implicit real*8 (a-h,o-z)
      common /memor1/ iisd,jjsd,ijbld
cpnl: common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memor1b/ nbl2,nbloks
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /pnlRJH/ isblcall
      dimension bl(*)
      dimension bl_pnl(*)
c
      call whatbl2(nbl2,bl(ijbld),bl(npard),bl_pnl(nblock2))
c
c-----
c for performance estimation :
c (how many times a given s-bl was called)
c
      nsupblk=nbl2*(nbl2+1)/2 
      call getmem(nsupblk, isblcall)
      call zeroit(bl(isblcall),nsupblk)
c-----
c
      end
c---------------------
      subroutine whatbl2(nbl2,ijbl,npar, nblock2)
      dimension ijbl(nbl2,*), npar(nbl2), nblock2(*)
c
      do 10 ibl=1,nbl2
      ijpar=npar(ibl)
         do 20 ijp=1,ijpar
         ijcs=ijbl(ibl,ijp)
         nblock2(ijcs)=ibl 
   20    continue
   10 continue
c
      end
c======================================================================
      subroutine txs_scratch_siz1(ncs,inx,nsym,ntxs_scr_size)
c
c estimate memory requested in texas-scratch bl() by prepint2
c
      common /intlim/ limxmem,limblks,limpair
      dimension inx(12,*)
c
      iexch=1
      do 55 ics0=1,ncs-1
      icst0=inx(12,ics0)
      icsc0=inx(5,ics0)-inx(1,ics0)
      icsg0=inx(4,ics0)
      ics1=ics0+1
      icst1=inx(12,ics1)
      icsc1=inx(5,ics1)-inx(1,ics1)
      icsg1=inx(4,ics1)
      if(icst0.ne.icst1) then
        iexch=iexch+1
      else
        if(icsc0.ne.icsc1) iexch=iexch+1 
        if(icsg0.ne.icsg1) iexch=iexch+1 
      endif
   55 continue
c
c  number of different shells is : iexch 
c-------------------------------------------------------
c number of pairs :
c
      ncspair=ncs*(ncs+1)/2
c
c estimate number of pair-blocks :
c
      nbl2=iexch*(iexch+1)/2
      nbl2_size=ncspair/nbl2
c
c if above size is greater than the limit then 
c take 75% of the limpair as a size :
c
      if(nbl2_size.gt.limpair) nbl2=ncspair/(3*limpair/4)
c
c estimate number of quartet-blocks (super-blocks):
c 
      nbl4=nbl2*(nbl2+1)/2
c-------------------------------------------------------
c Memory estimation for prepint2 (blocin2) :
c 10 memory allocations :
c
      mem1=ncspair
      mem2=ncspair
      mem3=ncspair*nsym
      mem4=ncspair
      mem5=nbl4
c
      mem6=nbl2
      mem7=nbl2*limpair
      mem8=nbl4 
c
c       see blockin2 where 10K/intsize is used in preliminary res.
        mem678=15 003     ! 3* 5001
c
      mem9=nbl4
      mem0=nbl4
c
c-------------------------------------------------------
      ntxs_scr_size= mem1+mem2+mem3+mem4+mem5+mem6+mem7+mem8+mem9+mem0
      ntxs_scr_size= ntxs_scr_size + mem678
c
c there is one more call to getmem to allocate isblcall
c in PAIR_IN2 routine :
c
      ntxs_scr_size= ntxs_scr_size + nbl4
c-------------------------------------------------------
c test
c
c     write(8,*)' number of different shells is :',iexch
c     write(8,*)' number of pair-blocks         :',nbl2 
c     write(8,*)' number of quartet-blocks      :',nbl4 
c     write(8,*)' maximum memory in txs-scratch :',ntxs_scr_size
c-------------------------------------------------------
      end
c===============================================================
      subroutine txs_scratch_siz2(maxme1,maxme2)
c
c estimate memory requested in texas-scratch bl() by texas_hf2_m
c
      common /memmax/ ispblx, maxmem1,iforwhat
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c maxme1 - max. memory needed for one block - calculated in blocin2
c
c maxme2 - max. memory needed in texas_setup :
c          5*nquartets_pnl + 2*nbl2 + 2 *nbl4
c------------------------------------------------
c assumptin : nquarts_pnl <= 250 000          
c
      data nquarts /250 000/
c------------------------------------------------
      nbl4=nbl2*(nbl2+1)/2
c
      maxme1=maxmem1
      maxme2=5*nquarts + 2*nbl2 + 2*nbl4
c
c------------------------------------------------
      end
