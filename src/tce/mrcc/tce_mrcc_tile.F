      subroutine tce_mrcc_tile(rtdb)
c
c ------------------------
c Modified for the MR case
c by Jiri Brabec
c ------------------------
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
#include "tce_mrcc.fh"
#include "tce_main.fh"

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
cjb
c      logical mrccdebug
c      integer nref
      integer isact
      integer l_test,k_test
      integer iref
      logical active_tmp(max_size)
      logical lact
      integer iswap
      logical lswap
      integer ifrom,ito
      integer ispinfrom,ispinto
      integer iu
c      logical lusescffv,limprovet
c      logical limprovet
      integer imaxiact
cjb temporary arrays
      integer g_movecs_tmp(2)
      integer l_movecs_sortedmt,k_movecs_sortedmt
      logical lisa
c kk--
c      integer itype,isize1,isize2
c kk--
c
c     ======================
c     Initialize parameters
c     ======================
c
c       if (.not.rtdb_get(rtdb,'mrcc:debug',mt_log,1,mrccdebug))
c     1 mrccdebug = .false.
c      if (.not.rtdb_get(rtdb,'bwcc:nref', mt_int, 1, nref))
c     1  call ERRQUIT('tce_mrcc_readref',1,RTDB_ERR)
c       if (.not.rtdb_get(rtdb,'mrcc:usescffermiv',mt_log,1,lusescffv))
c     1 lusescffv = .false.
c       if (.not.rtdb_get(rtdb,'mrcc:improvetiling',mt_log,1,limprovet))
c     1 limprovet = .false.
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
      if (nodezero) then
        write(LuOut,*)
        call util_print_centered
     1    (LuOut,'Memory Information',20,.true.)
        if (ga_uses_ma()) then
          write(LuOut,9000) int(ga_dble)
        else
          write(LuOut,9010) int(ga_dble),int(ma_dble)
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
      size = dble(available)
      size = size / dble(nbf*(nbf+1)/2)
      size = size / 2.0d0
      size = dsqrt(size)
      isize = int(size)
      if (model.eq."bwccsd") then
        size = dble(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."bwccsdt") then
        size = dble(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else if (model.eq."mkccsd") then
        size = dble(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
c        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else if (model.eq."mkccsdt") then
        size = dble(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else if (model.eq."succsd") then
        size = dble(available)
        size = size / 6.0d0
        size = dsqrt(dsqrt(size))
        isize = min(isize,int(size))
      else if (model.eq."mbpt4") then
        size = dble(available)
        size = size / 6.0d0
        size = size ** (1.0d0/8.0d0)
        isize = min(isize,int(size))
      else
        call errquit("tce_tile: the model NYI",0,CAPMIS_ERR)
      endif
      if (nodezero) write(LuOut,*)
      if (rtdb_get(rtdb,'tce:tilesize',mt_int,1,jsize)) then
        tile_dim = jsize
        isize = jsize
        if (nodezero) write(LuOut,'(A,I5,A)')'Maximum block size',isize,
     +   ' supplied by input'
      endif
        if (nodezero) write(LuOut,'(A,A)')'Maximum block size',
     +   ' doubles'

      pgbsize = isize

c
c     ================================================
c     Reorder the orbitals in the ha, hb, pa, pb order
c     ================================================
c
       bignumx = 10.0d+10
       if (.not.rtdb_get(rtdb,'tce:vcut',mt_dbl,1,vcut)) then
         vcut = bignumx
       end if
c
c     =====================================
c     Alloc of internal arrays
c     =====================================
c
      do i=1,nref

      if (.not.ma_push_get(mt_dbl,nbf*(nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol)),"sorted MO coeffs",
     2  l_movecs_sortedm(i),k_movecs_sortedm(i)))
     3  call errquit("tce_tile: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted spins",
     2  l_spin_sortedm(i),k_spin_sortedm(i)))
     3  call errquit("tce_tile: MA problem",1,MA_ERR)

      if (.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted irs",
     2  l_irs_sortedm(i),k_irs_sortedm(i)))
     3  call errquit("tce_tile: MA problem",2,MA_ERR)

      if (.not.ma_push_get(mt_dbl,nmo(1)-nfv(1)-nfc(1)
     1  +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted evl",
     2  l_evl_sortedm(i),k_evl_sortedm(i)))
     3  call errquit("tce_tile: MA problem",3,MA_ERR)

      if(.not.ma_push_get(mt_int,nmo(1)-nfv(1)-nfc(1)
     1 +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted MO index",
     2 l_mo_indexm(i),k_mo_indexm(i)))
     3 call errquit("tce_tile: MA problem",4,ma_err)

      if(.not.ma_push_get(mt_log,nmo(1)-nfv(1)-nfc(1)
     1 +nmo(ipol)-nfv(ipol)-nfc(ipol),"sorted isactive",
     2 l_isactive(i),k_isactive(i)))
     3 call errquit("tce_tile: MA problem",5,ma_err)

      end do
c
c     =====================================
c     Alloc of internal arrays
c     =====================================
c
      do i=1,nref

      if (.not.ma_push_get(mt_int,nmo(1)+nmo(ipol),"spin iref",
     1  l_spin_tmpm(i),k_spin_tmpm(i)))
     2  call errquit("tce_mrcc_tile: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,nmo(1)+nmo(ipol),"sym iref",
     1  l_sym_tmpm(i),k_sym_tmpm(i)))
     2  call errquit("tce_mrcc_tile: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_int,nmo(1)+nmo(ipol),"range iref",
     1  l_range_tmpm(i),k_range_tmpm(i)))
     2  call errquit("tce_mrcc_tile: MA problem",0,MA_ERR)

      if (.not.ma_push_get(mt_log,nmo(1)+nmo(ipol),"active iref",
     1  l_active_tmpm(i),k_active_tmpm(i)))
     2  call errquit("tce_mrcc_tile: MA problem",0,MA_ERR)

      end do

      do i=1,ipol
        if (.not.ma_push_get(mt_dbl,nbf,'eigenvalues',
     1    l_evl_tile(i),k_evl_tile(i))) call errquit
     2    ('tce_mrcc_tile: failed to allocate evl',0,MA_ERR)
        if (.not.ma_push_get(mt_int,nbf,'irreps',
     1    l_irs_tile(i),k_irs_tile(i))) call errquit
     2    ('tce_mrcc_tile: failed to allocate irs',0,MA_ERR)
      enddo
c kk ---
ccx      if(.not.ga_create(mt_dbl,nbf,nmo,'alpha vec',nbf,0,
ccx     1 g_movecs_tmp(1)))
ccx     1 call errquit("tce_mrcc_tile: GA problem",1,MA_ERR)
ccx      if(.not.ga_create(mt_dbl,nbf,nmo,'beta vec',nbf,0,
ccx     1 g_movecs_tmp(2)))
ccx     1 call errquit("tce_mrcc_tile: GA problem",2,MA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'alpha vec',nbf,0,
     1 g_movecs_tmp(1)))
     1 call errquit("tce_mrcc_tile: GA problem",1,MA_ERR)
      if(.not.ga_create(mt_dbl,nbf,nbf,'beta vec',nbf,0,
     1 g_movecs_tmp(2)))
     1 call errquit("tce_mrcc_tile: GA problem",2,MA_ERR)
c kk ---
cc      do iref=1,nref
cc        if(.not.ga_create(mt_dbl,nbf,nmo,'alpha vec',nbf,0,
cc     1    g_movecsm(iref,1)))
cc     1    call errquit("tce_mrcc_tile: GA problem",1,MA_ERR)
cc        if(.not.ga_create(mt_dbl,nbf,nmo,'beta vec',nbf,0,
cc     1    g_movecsm(iref,2)))
cc     1    call errquit("tce_mrcc_tile: GA problem",2,MA_ERR)
cc      enddo
c
c
      do iref=1,nref
        if(.not.ga_create(mt_dbl,nbf,nbf,'alpha vec',nbf,0,
     1    g_movecsm(iref,1)))
     1    call errquit("tce_mrcc_tile: GA problem",1,MA_ERR)
        if(.not.ga_create(mt_dbl,nbf,nbf,'beta vec',nbf,0,
     1    g_movecsm(iref,2)))
     1    call errquit("tce_mrcc_tile: GA problem",2,MA_ERR)
      enddo

      nirreps = sym_number_ops(geom) + 1

      if (nodezero.and.mrccdebug) then
       write(LuOut,"('Number of irreps = ',I3)")nirreps
        write(LuOut,"(/,'-------------------------------------------')")
        write(LuOut,"('| ',A6,' |',A7,'|',A10,'|')")'Number',' Spin ',
     1' Symmetry '
        write(LuOut,"('-------------------------------------------',/)")
        do any = nfc(1)+1,nmo(1)-nfv(1)
          write(LuOut,"('| ',I6,' | alpha | ',I8,' |')") any,
     1int_mb(k_irs(1)+any-1)
        enddo
        do any = nfc(ipol)+1,nmo(ipol)-nfv(ipol)
          write(LuOut,"('| ',I6,' | beta  | ',I8,' |')") any,
     1int_mb(k_irs(ipol)+any-1)
        enddo
      endif
c
c --------------------------------------
c Warning. Code not ready for this vars 
c --------------------------------------
c
c      if((nfc(1).gt.0).or.(nfc(2).gt.0).or.(nfv(1).gt.0)
c     1.or.(nfv(2).gt.0).or.(oactive(1).gt.0).or.(oactive(2).gt.0)
c     2.or.(vactive(1).gt.0).or.(vactive(2).gt.0)) then
c        write(LuOut,"(/,'Tiling with frozen and explicit active orbs
c     1is not implemented',/)")
c        call errquit("tce_mrcc_tile: MA problem",4,ma_err)    
c      end if

      if(intorb.and.nodezero) then
        write(LuOut,"(/,'Intorb not available for the MR case.',/)")
        call errquit("tce_mrcc_tile: MA problem",200,ma_err)
      end if

      if(.not.intorb) then

      if(mrccdebug.and.nodezero)
     1  write(LuOut,"(/,'Starting tiles',/)")

        if(nodezero)
     1 write(LuOut,"(/,'nbf,nmo(1):',I6,I6,/)")nbf,nmo(1)


c
c ================
c    MR-tile loop
c ================
c
      do iref=1,nref
c
c ================
c     Ref swap
c ================
c for rhf ??
         do j=1,nbf
            do iu=1,ipol
              int_mb(k_irs_tile(iu)+j-1) = int_mb(k_irs(iu)+j-1)
              dbl_mb(k_evl_tile(iu)+j-1) = dbl_mb(k_evl(iu)+j-1)
              moindexes(j,iu,iref) = j
            enddo
         enddo
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step1'
c         write(6,*)'ipol = ',ipol
c         write(6,*)'nbf nmo(1) nmo(2) ',nbf,nmo(1),nmo(2)
c         call util_flush(6)
c       end if
c       call ga_inquire(g_movecs(1),itype,isize1,isize2)
c       if(nodezero) then
c        write(6,*)'g_movecs_1 type dim1 dim2 ',itype,isize1,isize2
c       end if
c       call ga_inquire(g_movecs_tmp(1),itype,isize1,isize2)
c       if(nodezero) then
c        write(6,*)'g_movecs_tmp_1 type dim1 dim2 ',itype,isize1,isize2
c       end if
c ------
         call ga_copy(g_movecs(1),g_movecs_tmp(1))
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step1a'
c         call util_flush(6)
c       end if
c ------
         if(ipol.eq.2)call ga_copy(g_movecs(ipol),g_movecs_tmp(ipol))
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step2'
c         call util_flush(6)
c       end if
c ------
         if(.not.lusescffv) then

         do j=1,ihefferank(1,iref)
            
            ispinfrom = 1
            ispinto = 1

            ifrom = mutualoper(j,1,iref)
            ito = mutualoper(2*ihefferank(1,iref)-j+1,1,iref)

            if(ifrom.lt.0) then
              ispinfrom = ipol
              ifrom = -ifrom
            endif

            if(ito.lt.0) then
              ispinto = ipol
              ito = -ito
            endif
c
c ----------------
c       Swap 
c ----------------
c
            if(ispinfrom.ne.ispinto) then
              write(LuOut,*)'ERROR: we are not doing spin-flip dude!'
            endif
            moindexes(ifrom,ispinfrom,iref) = ito
            moindexes(ito,ispinto,iref) = ifrom


            int_mb(k_irs_tile(ispinto)+ito-1) = 
     1 int_mb(k_irs(ispinfrom)+ifrom-1)
            int_mb(k_irs_tile(ispinfrom)+ifrom-1) = 
     1 int_mb(k_irs(ispinto)+ito-1)
     
            dbl_mb(k_evl_tile(ispinto)+ito-1) = 
     1 dbl_mb(k_evl(ispinfrom)+ifrom-1)
            dbl_mb(k_evl_tile(ispinfrom)+ifrom-1) =
     1 dbl_mb(k_evl(ispinto)+ito-1)
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step3'
c         call util_flush(6)
c       end if
c ------
            call ga_copy_patch('N',g_movecs(ispinfrom),1,nbf,
     1 ifrom,ifrom,g_movecs_tmp(ispinto),1,nbf,ito,ito)

c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step4'
c         call util_flush(6)
c       end if
c ------

            call ga_copy_patch('N',g_movecs(ispinto),1,nbf,
     1 ito,ito,g_movecs_tmp(ispinfrom),1,nbf,ifrom,ifrom)
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step5'
c         call util_flush(6)
c       end if
c ------

         enddo

         else

c         do j=1,ihefferank(1,iref)

c            ispinfrom = 1
c            ispinto = 1

c            ifrom = mutualoper(j,1,iref)
c            ito = mutualoper(2*ihefferank(1,iref)-j+1,1,iref)

c            if(ifrom.lt.0) then
c              ispinfrom = ipol
c              ifrom = -ifrom
c            endif

c            if(ito.lt.0) then
c              ispinto = ipol
c              ito = -ito
c            endif

c            if(ispinfrom.ne.ispinto) then
c              write(LuOut,*)'ERROR: Spin from .NE. to !!!'
c            endif

c            moindexes(ifrom,ispinfrom,iref) = ito
c            moindexes(ito,ispinto,iref) = ifrom
c            write(6,"(I4,'-->',I4)")ifrom,ito
c         enddo

         do j=1,ihefferank0(1,iref)

            ispinfrom = 1
            ispinto = 1

            ifrom = mutualoper0(j,iref,1)
            ito = mutualoper0(2*ihefferank0(1,iref)-j+1,iref,1)

            if(ifrom.lt.0) then
              ispinfrom = ipol
              ifrom = -ifrom
            endif

            if(ito.lt.0) then
              ispinto = ipol
              ito = -ito
            endif

c            write(6,"(I4,'-->',I4)")ifrom,ito
c
c ----------------
c       Swap 
c ----------------
c
            if(ispinfrom.ne.ispinto) then
              write(LuOut,*)'ERROR: Spin from .NE. to !!!'
            endif

            moindexes(ifrom,ispinfrom,iref) = ito
            moindexes(ito,ispinto,iref) = ifrom

            int_mb(k_irs_tile(ispinto)+ito-1) =
     1 int_mb(k_irs(ispinfrom)+ifrom-1)
            int_mb(k_irs_tile(ispinfrom)+ifrom-1) =
     1 int_mb(k_irs(ispinto)+ito-1)

            dbl_mb(k_evl_tile(ispinto)+ito-1) =
     1 dbl_mb(k_evl(ispinfrom)+ifrom-1)
            dbl_mb(k_evl_tile(ispinfrom)+ifrom-1) =
     1 dbl_mb(k_evl(ispinto)+ito-1)

            call ga_copy_patch('N',g_movecs(ispinfrom),1,nbf,
     1 ifrom,ifrom,g_movecs_tmp(ispinto),1,nbf,ito,ito)
            call ga_copy_patch('N',g_movecs(ispinto),1,nbf,
     1 ito,ito,g_movecs_tmp(ispinfrom),1,nbf,ifrom,ifrom)

         enddo

         endif


c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step6'
c         call util_flush(6)
c       end if
c ------
       call ga_copy(g_movecs_tmp(1),g_movecsm(iref,1))
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step7'
c         call util_flush(6)
c       end if
c ------
       call ga_copy(g_movecs_tmp(2),g_movecsm(iref,2))
c kk ---
c       if(nodezero) then
c         write(6,*)'tce_mrcc_tile step8'
c         call util_flush(6)
c       end if
c ------
c
c --------------------
c     Print debug
c --------------------
c
      if(mrccdebug.and..false.) then
         call ga_print(g_movecs_tmp(1))
         if(ipol.eq.2)
     1    call ga_print(g_movecs_tmp(ipol))
      endif


      if (.not.ma_push_get(mt_dbl,nbf*nbf*2,"coeffs",
     2  l_movecs_sortedmt,k_movecs_sortedmt))
     3  call errquit("tce_tile: MA problem",0,MA_ERR)

        if(nodezero) then
        call ga_get(g_movecs_tmp(1),1,nbf,1,nbf,
     1  dbl_mb(k_movecs_sortedmt),nbf)
        endif
       call ga_brdcst(65537, dbl_mb(k_movecs_sortedmt),
     G     MA_sizeof(MT_DBL,1,MT_BYTE)*nbf*nbf, 0)



        if(nodezero) then
        call ga_get(g_movecs_tmp(2),1,nbf,1,nbf,
     1  dbl_mb(k_movecs_sortedmt+nbf*nbf),nbf)
        endif
       call ga_brdcst(65538, dbl_mb(k_movecs_sortedmt+nbf*nbf),
     G     MA_sizeof(MT_DBL,1,MT_BYTE)*nbf*nbf, 0)

c        if(nodezero)
c     1 write(LuOut,"(/,'nbf,nmo(1):',I6,I6,/)")nbf,nmo(1)
c     1  call ma_print(dbl_mb(k_movecs_sortedmt),
c     2 nbf*(nmo(1)+nmo(2)),1,'MOVECS')

c        call ga_print(g_movecs_tmp(1))

c
c ================
c    Hole Alpha
c ================
c
      i = 0
      noa = 0
      imaxiact = 2
      if(limprovet)imaxiact=1
      do isact=1,imaxiact
        if(isact.eq.1)lact=.false.
        if(isact.eq.2)lact=.true.
      do sym = 0, nirreps
        j_ia = 0           !occ. inactive index
        j_a  = 0           !occ.   active index
        do hole = nfc(1)+1, nocc(1)-oactive(1)
          if((isactive(hole,1).eqv.lact).or.limprovet) then
          if (int_mb(k_irs_tile(1)+hole-1).eq.sym) then
          if((isactive(hole,1).eqv..true.).and.limprovet)lact=.true.
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=1
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(1)+hole-1)
c             call ga_get(g_movecs_tmp(1),1,nbf,hole,hole,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
             call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+(hole-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*hole-1
c<-d3p975
          endif
          endif ! isact
        enddo
        do hole =nocc(1)-oactive(1)+1,nocc(1)
          if (int_mb(k_irs_tile(1)+hole-1).eq.sym) then
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=1
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(1)+hole-1)
c             call ga_get(g_movecs_tmp(1),1,nbf,hole,hole,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
             call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+(hole-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)

c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*hole-1
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
             int_mb(k_spin_tmpm(iref)+noa-1) = 1
             int_mb(k_sym_tmpm(iref)+noa-1) = sym
             log_mb(k_active_tmpm(iref)+noa-1) = lact
             int_mb(k_range_tmpm(iref)+noa-1) = k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + int_mb(k_range_tmpm(iref)+noa-1)
c             write(LuOut,"('k=',I3,' range_tmp=',I8,' noa=',I4,
c     1' sym=',I4,' act=',L1)")
c     1 k,int_mb(k_range_tmpm(iref)+noa-1),
c     2 noa,int_mb(k_sym_tmpm(iref)+noa-1),lact
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
      enddo !isact
      oatiles(1)=nblocks_a

c
c     Hole Beta
c
      nob = 0
      imaxiact = 2
      if(limprovet)imaxiact=1
      do isact=1,imaxiact
        if(isact.eq.1)lact=.false.
        if(isact.eq.2)lact=.true.
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        do hole = nfc(ipol)+1, nocc(ipol)-oactive(ipol)
          if((isactive(hole,2).eqv.lact).or.limprovet) then
          if (int_mb(k_irs_tile(ipol)+hole-1).eq.sym) then
          if((isactive(hole,2).eqv..true.).and.limprovet)lact=.true.
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=2
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(ipol)+hole-1)
c             call ga_get(g_movecs_tmp(ipol),1,nbf,hole,hole,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
             call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+
     1 nbf*nbf+(hole-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)

c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*hole
c<-d3p975
          endif
          endif
        enddo
        do hole=nocc(ipol)-oactive(ipol)+1,nocc(ipol)
          if (int_mb(k_irs_tile(ipol)+hole-1).eq.sym) then
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=2
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(ipol)+hole-1)
c             call ga_get(g_movecs_tmp(ipol),1,nbf,hole,hole,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
             call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+
     1 nbf*nbf+(hole-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)

c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*hole
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
             int_mb(k_spin_tmpm(iref)+noa+nob-1) = 2
             int_mb(k_sym_tmpm(iref)+noa+nob-1) = sym
             log_mb(k_active_tmpm(iref)+noa+nob-1) = lact
             int_mb(k_range_tmpm(iref)+noa+nob-1) = 
     1 k*j_ia/nblocks_ia-l_ia
             l_ia = l_ia + int_mb(k_range_tmpm(iref)+noa+nob-1)
c             write(LuOut,"('k=',I3,' range_tmp=',I8,' nob=',I4,
c     1' sym=',I4,' act=',L1)")
c     1 k,int_mb(k_range_tmpm(iref)+noa+nob-1),
c     2 nob,int_mb(k_sym_tmpm(iref)+noa+nob-1),lact
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
      enddo ! isact
      oatiles(ipol)=nblocks_a
c
c     Particle Alpha
c
      nva = 0
      imaxiact = 2
      if(limprovet) imaxiact=1
      do isact=1,imaxiact
        if(isact.eq.1)lact=.true.
        if(isact.eq.2)lact=.false.
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        if(limprovet)lact=.false.
        do particle = nocc(1)+1, nocc(1)+vactive(1)
          if((isactive(particle,1).eqv.lact).or.limprovet) then
          if (int_mb(k_irs_tile(1)+particle-1).eq.sym) then
          if((isactive(particle,1).eqv..true.).and.limprovet)lact=.true.
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=1
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(1)+particle-1)
c             call ga_get(g_movecs_tmp(1),1,nbf,particle,particle,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
       call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+(particle-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)

c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*particle-1
c<-d3p975
          endif
          endif
        enddo
        do particle = nocc(1)+vactive(1)+1, nmo(1)-nfv(1)
          if((isactive(particle,1).eqv.lact).or.limprovet) then
          if (int_mb(k_irs_tile(1)+particle-1).eq.sym) then
          if((isactive(particle,1).eqv..true.).and.limprovet)lact=.true.
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=1
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(1)+particle-1)
c             call ga_get(g_movecs_tmp(1),1,nbf,particle,particle,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
       call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+(particle-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*particle-1
c<-d3p975
          endif
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
             log_mb(k_active_tmpm(iref)+noa+nob+nva-1) = lact
             int_mb(k_spin_tmpm(iref)+noa+nob+nva-1) = 1
             int_mb(k_sym_tmpm(iref)+noa+nob+nva-1) = sym
             int_mb(k_range_tmpm(iref)+noa+nob+nva-1) =
     1 k*j_ia/nblocks_ia-l_ia
              l_ia = l_ia + int_mb(k_range_tmpm(iref)+noa+nob+nva-1)
c            write(LuOut,"('k=',I3,' range_tmp=',I8,' nva=',I4,
c     1' sym=',I4,' act=',L1)")k,
c     2 int_mb(k_range_tmpm(iref)+noa+nob+nva-1),nva,sym,lact
           enddo
        endif
      enddo
      enddo
      vatiles(1)=nblocks_a
c
c     Particle Beta
c
      nvb = 0
      imaxiact = 2
      if(limprovet)imaxiact=1
      do isact=1,imaxiact
        if(isact.eq.1)lact=.true.
        if(isact.eq.2)lact=.false.
      do sym = 0, nirreps
        j_ia = 0
        j_a  = 0
        if(limprovet)lact=.false.
        do particle = nocc(ipol)+1, nocc(ipol)+vactive(ipol)
          if (int_mb(k_irs_tile(ipol)+particle-1).eq.sym) then
             i = i + 1
             j_a = j_a + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=2
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(ipol)+particle-1)
c             call ga_get(g_movecs_tmp(ipol),1,nbf,particle,particle,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
       call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+
     1 nbf*nbf+(particle-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)

c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*particle
c<-d3p975
          endif
        enddo
        do particle = nocc(ipol)+vactive(ipol)+1, nmo(ipol)-nfv(ipol)
          if((isactive(particle,2).eqv.lact).or.limprovet) then
          if (int_mb(k_irs_tile(ipol)+particle-1).eq.sym) then
          if((isactive(particle,2).eqv..true.).and.limprovet)lact=.true.
             i = i + 1
             j_ia = j_ia + 1
             int_mb(k_irs_sortedm(iref)+i-1)=sym
             int_mb(k_spin_sortedm(iref)+i-1)=2
             dbl_mb(k_evl_sortedm(iref)+i-1)=
     1 dbl_mb(k_evl_tile(ipol)+particle-1)
c             call ga_get(g_movecs_tmp(ipol),1,nbf,particle,particle,
c     1         dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
       call tce_mrcc_macpy(dbl_mb(k_movecs_sortedmt+
     1 nbf*nbf+(particle-1)*nbf),
     1 dbl_mb(k_movecs_sortedm(iref)+(i-1)*nbf),nbf)
c->d3p975
             int_mb(k_mo_indexm(iref)+i-1)=2*particle
c<-d3p975
          endif
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
             log_mb(k_active_tmpm(iref)+noa+nob+nva+nvb-1) = lact 
             int_mb(k_spin_tmpm(iref)+noa+nob+nva+nvb-1) = 2
             int_mb(k_sym_tmpm(iref)+noa+nob+nva+nvb-1) = sym 
             int_mb(k_range_tmpm(iref)+noa+nob+nva+nvb-1) =
     1 k*j_ia/nblocks_ia-l_ia 
             l_ia = l_ia + int_mb(k_range_tmpm(iref)+noa+nob+nva+nvb-1) 
c            write(LuOut,"('k=',I3,' range_tmp=',I8,' nvb=',I4,
c     1' sym=',I4,' act=',L1)")k,
c     2 int_mb(k_range_tmpm(iref)+noa+nob+nva+nvb-1),nvb,sym,lact
           enddo
        endif
      enddo
      enddo
      vatiles(ipol)=nblocks_a
c
c     Holes and particles
c
      noab = noa + nob
      nvab = nva + nvb

      nblcks(1,iref) = noa
      nblcks(2,iref) = nob
      nblcks(3,iref) = nva
      nblcks(4,iref) = nvb

      if (.not.ma_pop_stack(l_movecs_sortedmt))
     1  call errquit("tce_mrcc_tile: MA problem",10,MA_ERR)

      enddo

      if(mrccdebug.and.nodezero)
     1  write(LuOut,"(/,'Tiles created',/)")
c
      endif !not intorb
c
c
c     Allocate tiling related arrays for nref references
c
      do iref=1,nref

      noa = nblcks(1,iref)
      nob = nblcks(2,iref)
      nva = nblcks(3,iref)
      nvb = nblcks(4,iref)

      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Spin of blocks",
     1  l_spinm(iref),k_spinm(iref))) 
     2 call errquit("tce_mrcc_tile: MA problem",4,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Symmetry of blocks",
     1  l_symm(iref),k_symm(iref)))
     2 call errquit("tce_mrcc_tile: MA problem",5,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Range of blocks",
     1  l_rangem(iref),k_rangem(iref))) 
     2 call errquit("tce_mrcc_tile: MA problem",6,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"Offset of blocks",
     1  l_offsetm(iref),k_offsetm(iref)))
     2 call errquit("tce_tile: MA problem",7,MA_ERR)
      if (.not.ma_push_get(mt_int,noa+nob+nva+nvb,"alpha-izer",
     1  l_alpham(iref),k_alpham(iref)))
     2 call errquit("tce_tile: MA problem",8,MA_ERR)

      enddo
c
      if(.not.intorb) THEN
c

      do iref=1,nref
      j = 0
      do i = 1,nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+
     2 nblcks(4,iref)
        int_mb(k_spinm(iref)+i-1)=int_mb(k_spin_tmpm(iref)+i-1)
        int_mb(k_symm(iref)+i-1)=int_mb(k_sym_tmpm(iref)+i-1)
        int_mb(k_rangem(iref)+i-1)=int_mb(k_range_tmpm(iref)+i-1)
        int_mb(k_active+i-1)=active_tiles(i)
        int_mb(k_offsetm(iref)+i-1)=j
        log_mb(k_isactive(iref)+i-1)=log_mb(k_active_tmpm(iref)+i-1)
        j = j + int_mb(k_range_tmpm(iref)+i-1)
      enddo
      if (restricted) then
        do i = 1,nblcks(1,iref)
          int_mb(k_alpham(iref)+i-1) = i
        enddo
        do i = nblcks(1,iref)+1,nblcks(1,iref)+nblcks(2,iref)
          int_mb(k_alpham(iref)+i-1) = i - nblcks(1,iref)
        enddo
        do i = nblcks(1,iref)+nblcks(2,iref)+1,nblcks(1,iref)+
     2 nblcks(2,iref)+nblcks(3,iref)
          int_mb(k_alpham(iref)+i-1) = i
        enddo
        do i = nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+1,
     1 nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+nblcks(4,iref)
          int_mb(k_alpham(iref)+i-1) = i - nblcks(3,iref)
        enddo
      else
        do i = 1,nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+
     1 nblcks(4,iref)
          int_mb(k_alpham(iref)+i-1) = i
        enddo
      endif
      enddo ! iref
c
c
      end if !not intorb
c
      do iref=1,nref

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
        write(LuOut,"(/,'Blocks for reference ',I4)")iref
        write(LuOut,9070)
        do i = 1,nblcks(1,iref)+nblcks(2,iref)+nblcks(3,iref)+
     1 nblcks(4,iref)
          call sym_irrepname(geom,int_mb(k_symm(iref)+i-1)+1,irrepname)
          write(LuOut,9060) i,spinname(int_mb(k_spinm(iref)+i-1)),
     1 irrepname, int_mb(k_rangem(iref)+i-1),
     1 int_mb(k_offsetm(iref)+i-1),
     1 int_mb(k_alpham(iref)+i-1)
     1,log_mb(k_isactive(iref)+i-1)
        enddo
      endif
      enddo ! iref
 
 9070 format(/,1x,'Block   Spin    Irrep     Size     Offset   Alpha   
     1Active',/,
     1 1x,'-----------------------------------------------------------')
 9060 format(1x,i3,4x,a5,5x,a4,i4,' doubles',i8,i8,L8)

c
c     ===========================================
c     For testing purpose only, must be deleted 
c     ===========================================
c

      if (.not.ga_destroy(g_movecs_tmp(1)))
     1  call errquit('tce_mrcc_tile: GA problem',1,GA_ERR)
      if (.not.ga_destroy(g_movecs_tmp(2)))
     1  call errquit('tce_mrcc_tile: GA problem',1,GA_ERR)

      if (.not.ga_destroy(g_movecs(1)))
     1  call errquit('tce_mrcc_tile: GA problem',1,GA_ERR)
      if (.not.ga_destroy(g_movecs(2)))
     1  call errquit('tce_mrcc_tile: GA problem',1,GA_ERR)

      return
      end
c $Id$
