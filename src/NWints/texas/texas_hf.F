c $Id: texas_hf.F,v 1.15 1996-08-02 22:03:21 d3g681 Exp $
c=================================================================
      subroutine texas_hf2_m(ij_basis,ish,jsh,kl_basis,ksh,lsh,nquart,
     *                       q4, use_q4,
     *                       ra,rb,rc,rd,use_r,
     *                       eri,leri,icf,jcf,kcf,lcf,integ_n0,labels,
     *                       more_int,blscr,l_blscr)
c--------------------------------------------------------
c This subroutine delivers two-el.integrals 
c
c (1) all of them (zeros also) WITHOUT labels
c (2) only non-zero integrals WITH labels
c
c NQUART quartets of contr. shells is requested
c        ish(nquart), ..,lsh(nquart)
c Integrals return in ERI(leri) without labels or with labels
c in icf(leri),..,lcf(leri) .
c--------------------------------------------------------
      implicit none
#include "mafdecls.fh"
      logical use_r,labels, more_int
c PNL scratch for calculations replacing our bl() and its size :
      integer l_blscr
      double precision blscr(l_blscr)
c PNL requested set of contracted shell quartets :
      integer ij_basis, kl_basis
      integer nquart,ish(nquart),jsh(nquart),ksh(nquart),lsh(nquart)
      logical use_q4
      double precision q4(nquart)
c returning integral's indeces :
      integer icf(*),jcf(*),kcf(*),lcf(*), integ_n0
      integer leri
      double precision eri(leri)
      double precision ra(3),rb(3),rc(3),rd(3)
c----------------------------------------------------------------------
      integer ntxs_bl_scr
      common /bl_txs_add/ ntxs_bl_scr 
c----------------------------------------------------------------------
      integer basis_init
      common /c_basis_init/ basis_init
c----------------------------------------------------------------------
      integer num_bas_1,num_bas_2,num_bas_3
      integer ncs_bas_1,ncs_bas_2,ncs_bas_3
      integer nps_bas_1,nps_bas_2,nps_bas_3
      integer nat_bas_1,nat_bas_2,nat_bas_3
      integer ncf_bas_1,ncf_bas_2,ncf_bas_3
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3, ! Basis set handle
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3, ! Cummulative #shells
     *                     nps_bas_1,nps_bas_2,nps_bas_3, ! Cummulative #prims unused
     *                     nat_bas_1,nat_bas_2,nat_bas_3, ! Cummulative #atoms unused
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3  ! Cummulative #basis functions
c
c to see if any of these first shells from this request are ZERO
c
      integer ish_first,jsh_first,ksh_first,lsh_first
c----------------------------------------------------------------------
c... check for consistency the input basis is what was initialized.
c
      if(.not.more_int) then
c     ----------------------------------------------
         if(ij_basis.eq.num_bas_1) go to 1001
         if(ij_basis.eq.num_bas_2) go to 1001
         if(ij_basis.eq.num_bas_3) go to 1001
c
         write(6,*)' basis sets initialized :'
         if(num_bas_1.gt.0) write(6,*) num_bas_1
         if(num_bas_2.gt.0) write(6,*) num_bas_2
         if(num_bas_3.gt.0) write(6,*) num_bas_3
         write(6,*)' ij_basis   handle :',ij_basis
         call errquit
     &   ('texas_hf2_m: called with different basis set handle',911)
 1001    continue
c     ----------------------------------------------
         if(kl_basis.eq.num_bas_1) go to 1002
         if(kl_basis.eq.num_bas_2) go to 1002
         if(kl_basis.eq.num_bas_3) go to 1002
c
         write(6,*)' basis sets initialized :'
         if(num_bas_1.gt.0) write(6,*) num_bas_1
         if(num_bas_2.gt.0) write(6,*) num_bas_2
         if(num_bas_3.gt.0) write(6,*) num_bas_3
         write(6,*)' kl_basis   handle :',kl_basis
         call errquit
     &   ('texas_hf2_m: called with different basis set handle',911)
 1002    continue
c        --------------------------------------------------------------
c        remember the first shells in this request(if ZERO or not)
         ish_first=ish(1)
         jsh_first=jsh(1)
         ksh_first=ksh(1)
         lsh_first=lsh(1)
c        --------------------------------------------------------------
c        Check what basis sets are involved and what type of integrals
c        is requested (4-c , 3-c or 2-c ; c=center)
c
         call request_update(ish,jsh,ksh,lsh,nquart,ij_basis,kl_basis )
c        --------------------------------------------------------------
c        switch from txs bl() to pnl blscr()
c
         call switch_scr(dbl_mb(ntxs_bl_scr),blscr,l_blscr)
c        --------------------------------------------------------------
      endif
c----------------------------------------------------------------------
      call mul_quart(ish,jsh,ksh,lsh,nquart, q4,use_q4,
     *               ra,rb,rc,rd,use_r, blscr,l_blscr,eri,leri,
     *               icf,jcf,kcf,lcf,integ_n0,labels, more_int)
c----------------------------------------------------------------------
c In the case of different basis sets update the label's arrays
c
      if(labels) then
         call labels_update(icf,jcf,kcf,lcf,integ_n0,ij_basis,kl_basis,
     *                      ish_first,jsh_first,ksh_first,lsh_first)
      endif
c----------------------------------------------------------------------
      end
c=================================================================
      subroutine mul_quart(icspnl,jcspnl,kcspnl,lcspnl,nquart,q4,use_q4,
     *                     ra,rb,rc,rd,use_r, bl   ,l_blscr,eri,leri,
     *                     icf,jcf,kcf,lcf,integ_nx,labels, more_int)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical use_q4,use_r,labels,more_int
      double precision ra(3), rb(3), rc(3), rd(3)
      common /index/ maxsh,ifp,inx(1)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c----------------------------------------------------------------
      common /memmax/ ispblx, maxme1,iforwhat
c----------------------------------------------------------------
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c----------------------------------------------------------------
      common /pnl006/ nsplit,isplit,isbl_split, isbl_part
c----------------------------------------------------------------
      common /pnl_time1/ time_pt,time_setup
      common /pnl_time2/ time_calc2,time_do2,time_do4,time_doqrt
c
      common /pnl_nqrt/ ncall_pnl, nqrts_pnl
c------------------------------pnl time--------------------------
cpnl  common /big/ bl(1)
      dimension bl(l_blscr), q4(nquart)
c
c requested quartets of shells :
      dimension icspnl(nquart),jcspnl(nquart),
     *          kcspnl(nquart),lcspnl(nquart)
c returning integrals and indeces :
      dimension eri(leri)
      dimension icf(*),jcf(*),kcf(*),lcf(*)
c----------------------------------------------------------------
c eri - PNL buffer for returning integrals 
c----------------------------------------------------------------
             call txs_second(tim_beg)
c----------------------------------------------------------------
c memory checking ;
c
c     call getmem(0,last11)
c     call retmem(1)
c----------------------------------------------------------------
c MORE_INT=true means that it is called for the same request again 
c
  999 continue
c
      if(more_int) then
         nsplit=nsplit+1
         go to 1000
      endif
c----------------------------------------------------------------
      ncall_pnl=ncall_pnl+1
      nqrts_pnl=nqrts_pnl+nquart
ctest
c      write(6,*)' call no=',ncall_pnl,' scratch size=',l_blscr
c----------------------------------------------------------------
      iforwhat=1
c----------------------------------------------------------------
c if use_r is true replace current coordinates by Ra-Rd
c
c     if (use_r) then
c         call replacer(bl(inuc+1),ra,rb,rc,rd,inx,nquart,
c    *                  icspnl,jcspnl,kcspnl,lcspnl,'forw')
c     endif
c----------------------------------------------------------------
c find texas-shells corresponding to the given set of pnl-shells.
c re-oredr TXS-requested shells in quartets to TXS-work ordering :
c find wihch TXS pairs are present in a request (ijpres2,klpres2)
c find all pair-blocks to which requested shells  belong :
c
      call txs_second(tim_1setup)
c
      call txs_setup(bl,inx,ncs,leri,
     *               nquart,icspnl,jcspnl,kcspnl,lcspnl,labels,
     *               more_int)
c
c 11 calls to getmem 
c 12 calls to getmem
c if(.not.labels) +1 call to getmem
c
      call txs_second(tim_2setup)
      time_setup=time_setup + (tim_2setup-tim_1setup)
c----------------------------------------------------------------
c if this pnl-request is too big (too many integrals at once) 
c
      if( isplit.gt.0 ) then
          more_int=.true.
          nsplit=0
          go to 999
      endif
c----------------------------------------------------------------
 1000 continue
c----------------------------------------------------------------
c calculate two-electron integrals :
c
c----------------------------------------------------------------
      if(.not.labels) then
         ispec=1
         integ_n0=0
c
c        calculate integrals without labels :
c
*     ncfunct is not used on this call
         call calcint2(bl,inx, q4,use_q4, icf,jcf,kcf,lcf,eri,
     *        dbl_mb(ncfunct))
         integ_nx=integ_n0
c--------------------------------------------------------
c test only
c       write(6,*)' from texas_hf_m : q4 vector :nquart=',nquart
c       write(6,1234) (q4(iii),iii=1,nquart)
c       if(ncall_pnl.eq.1 .and. nsplit.eq.1) then
c       call print_int2(ncall_pnl,eri,leri,
c    *             icspnl,jcspnl,kcspnl,lcspnl,integ_n0,nquart)
c       endif
c--------------------------------------------------------
c
      else
c
         ispec=2
         integ_n0=0
c
c        calculate integrals with labels :
c
         call calcint2(bl,inx, q4,use_q4, icf,jcf,kcf,lcf,eri,
     $        dbl_mb(ncfunct))
         integ_nx=integ_n0
c
c        find pnl-labels corresponding to the given set of 
c        txs-labels , transfer labels from txs to pnl:
c
            call txs_second(tim_1pt)
c
C@@@@@   call pnl_txs(bl(ncfunct),icf,jcf,kcf,lcf,integ_n0 )
******    call pnl_txs(dbl_mb(ncfunct),icf,jcf,kcf,lcf,integ_n0 )
c----------------------------------------------------------------
c test only
c       write(6,*)' from texas_hf_m : q4 vector :nquart=',nquart
c       write(6,1234) (q4(iii),iii=1,nquart)
c1234 format(10(f5.3,1x))
c       if(ncall_pnl.eq.1 .and. nsplit.eq.1) then
c       call print_int1(ncall_pnl,eri,leri,icf,jcf,kcf,lcf,integ_n0)
c       endif
c--------------------------------------------------------
c
            call txs_second(tim_2pt)
            time_pt=time_pt+(tim_2pt - tim_1pt)
      endif
c
      if(isplit.gt.0) then
         if(nsplit.eq.isplit+1) more_int=.false.
      endif
c
      if(.not.more_int) then
         if(.not.labels) then
            call retmem(13)
         else
            call retmem(12)
         endif
c-------------------------
c if use_r is true return original coordinates 
c Probably NOT needed !
c---
c        if (use_r) then
c         call replacer(bl(inuc+1),ra,rb,rc,rd,inx,nquart,
c    *                  icspnl,jcspnl,kcspnl,lcspnl,'back')
c        endif
c-------------------------
c check and update memory status in pnl_scratch 
c after each pnl request is done; 
c
         call memstat_pnl
c-------------------------
      endif
c----------------------------------------------------------------
      call txs_second(tim_end)
      time_calc2=time_calc2+(tim_end -tim_beg)
c----------------------------------------------------------------
c test only :
c
c     integ_upto_sum=integ_upto_sum + integ_nx
c     write(6,*)'FROM MUL_QUART(end):LERI=',leri,' call no=',ncall_pnl,
c    *'req.size=',nquart,' split=',isplit+1,' repeat=',nsplit
c     write(6,*)
c    *'integrals: from this run=',integ_nx,',total=',integ_upto_sum,
c    *' more_int=',more_int,' labels=',labels
c-----
c----------------------------------------------------------------
      if(integ_nx.eq.0 .and. more_int) go to 999
c----------------------------------------------------------------
c-----
c          write(8,*)'FROM MUL_QUART (end): call no.=',ncall_pnl,
c    *               ' integ.buffer=',leri,' scratch size=',l_blscr
c----------------------------------------------------------------
c memory checking ;
c
c     call getmem(0,last12)
c     call retmem(1)
c        if(last11.ne.last12) then
c          write(6,*)'** Memory alocations in MUL_QUART          **'
c          write(6,*)' ncall_pnl=',ncall_pnl,'more_int=',more_int
c          write(6,*)' isplit=',isplit,' nsplit=',nsplit
c          write(6,*)' at the beginning of calc.=',last11
c          write(6,*)' at the end of cal.       =',last12
c        endif
c----------------------------------------------------------------
c
      end
c==================================================================
      subroutine replacer(datnuc,ra,rb,rc,rd,inx,nquart,
     *                    icspnl,jcspnl,kcspnl,lcspnl, direction)
      implicit real*8 (a-h,o-z)
      character*4 direction
      common /keepr/  ri(3),rj(3),rk(3),rl(3)
      dimension       ra(3),rb(3),rc(3),rd(3)
      dimension datnuc(5,*) 
      dimension inx(12,*)
      dimension icspnl(*),jcspnl(*),kcspnl(*),lcspnl(*)  ! dim=nquart
c
      if(direction.eq.'forw') then
        do 100 iq=1,nquart
        ics=icspnl(iq)
        jcs=jcspnl(iq)
        kcs=kcspnl(iq)
        lcs=lcspnl(iq)
c
        iat=inx(2,ics)
        jat=inx(2,jcs)
        kat=inx(2,kcs)
        lat=inx(2,lcs)
c
          do 10 ii=1,3
            if(iat.gt.0) ri(ii)=datnuc(ii+1,iat)
            if(jat.gt.0) rj(ii)=datnuc(ii+1,jat)
            if(kat.gt.0) rk(ii)=datnuc(ii+1,kat)
            if(lat.gt.0) rl(ii)=datnuc(ii+1,lat)
   10     continue
          do 20 ii=1,3
            if(iat.gt.0) datnuc(ii+1,iat)=ra(ii)
            if(jat.gt.0) datnuc(ii+1,jat)=rb(ii)
            if(kat.gt.0) datnuc(ii+1,kat)=rc(ii)
            if(lat.gt.0) datnuc(ii+1,lat)=rd(ii)
   20     continue
  100   continue
      endif
c  
      if(direction.eq.'back') then
        do 200 iq=1,nquart
        ics=icspnl(iq)
        jcs=jcspnl(iq)
        kcs=kcspnl(iq)
        lcs=lcspnl(iq)
c
        iat=inx(2,ics)
        jat=inx(2,jcs)
        kat=inx(2,kcs)
        lat=inx(2,lcs)
          do 30 ii=1,3
            if(iat.gt.0) datnuc(ii+1,iat)=ri(ii)
            if(jat.gt.0) datnuc(ii+1,jat)=rj(ii)
            if(kat.gt.0) datnuc(ii+1,kat)=rk(ii)
            if(lat.gt.0) datnuc(ii+1,lat)=rl(ii)
   30     continue
  200   continue
      endif
c  
      end
c==================================================================
      subroutine txs_setup(bl,inx,ncs,leri,
     *                     nquart,icspnl,jcspnl,kcspnl,lcspnl,labels,
     *               more_int)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical labels
c more_int is needed only to execute (or not) part of supblks_split
      logical more_int 
c
      common /bl_txs_add/ ntxs_bl_scr 
c
      common /memor1/ iisd,jjsd,ijbld
cpnl: common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memor1b/ nbl2,nbloks
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /pnl003/ nqrtpnl,icstxs,jcstxs,kcstxs,lcstxs
      common /pnl004/ isize,jsize,ksize,lsize,itxspnl
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /pnl006/ nsplit,isplit,isbl_split, isbl_part
      common /pnlRJH/ isblcall
      dimension icspnl(*),jcspnl(*),kcspnl(*),lcspnl(*)  ! dim=nquart
      dimension inx(12,*)
      dimension bl(*)
c--------------------------------------------------------
      nqrtpnl=nquart
      nsupblk=nbl2*(nbl2+1)/2
c--------------------------------------------------------
c find texas-shells corresponding to the given set of pnl-shells.
c
      call getmem(nquart,icstxs)
      call getmem(nquart,jcstxs)
      call getmem(nquart,kcstxs)
      call getmem(nquart,lcstxs)
c
C@@@@ call txs_pnl(nquart,bl(ncshell),
      call txs_pnl(nquart,dbl_mb(ncshell),
     *             bl(icstxs),bl(jcstxs),bl(kcstxs),bl(lcstxs),
     *             icspnl,jcspnl,kcspnl,lcspnl)
c
c on return TXS shells in icspnl... again
c
      call retmem(4)
c
c--------------------------------------------------------
c reserve memory 
c for pairs ij, kl which are present in a request
c
      call getmem(nquart, ijpres2)
      call getmem(nquart, klpres2)
c-
      call getmem(nquart, iqorder)
c--------------------------------------------------------
c reserve memory 
c for pair-blocks ij, kl which are present in a request
c
      call getmem(nbl2, ijblock)
      call getmem(nbl2, klblock)
c--------------------------------------------------------
c reserve memory for present super-blocks
c
      call getmem(nquart,isblqrts)
      call getmem(nsupblk,isblsize)
c--------------------------------------------------------
c  7 calls to getmem 
c--------------------------------------------------------
      call getmem(nquart,isupblk)
      call getmem(nsupblk,isblpoin)
      call swap_shells(nquart,icspnl,jcspnl,kcspnl,lcspnl,
     *     nbl2,
     *     bl(ijblock),bl(klblock),bl(ijpres2),bl(klpres2),
C@@@ *     bl(isupblk),bl(isblsize),bl(isblqrts),bl(nblock2),
     *     bl(isupblk),bl(isblsize),bl(isblqrts),dbl_mb(nblock2),
     *     bl(iqorder),bl(isblpoin), 
     *     dbl_mb(ntxs_bl_scr-1+isblcall) )
      call retmem(1)
      call retmem(1)
c  
      if(.not.labels) then
         call get_sizes(inx,icspnl,jcspnl,kcspnl,lcspnl,
     *                      isize,jsize,ksize,lsize,ngcont,ngcq)
ccc   write(6,*)'general contraction deep=',ngcont
         call getmem(ngcont*isize*jsize*ksize*lsize,itxspnl)
c
         call ordr_shells(icspnl,jcspnl,kcspnl,lcspnl,
     *                    isize,jsize,ksize,lsize,bl(itxspnl),
     *                    bl(iqorder),ngcq)
c       1 call to getmem
      endif
c--------------------------------------------------------
c update number of pairs from TXS to current PNL (NPAR) :
c
      call uniq_pairs(bl,nbl2,bl(npard),
     *     bl(ijblock),bl(klblock),bl(ijpres2),bl(klpres2),
     *     bl(isblsize),bl(isblqrts) )
c       2 call to getmem
c---------------------------------------------------------
c check current sizes of super-blocks against maximum sizes
c and available size of integral's buffer (leri):
c (maximum sizes calculated in blksizer (spec_block.f) )
c
      call getmem(nsupblk,isbl_part)
      call supblks_split(bl(mxsize),leri,inx,
     *                   bl(iisd),bl(jjsd),bl(ijbld),nbl2,
     *     bl(ijpres2),bl(klpres2),bl(isblsize),bl(isbl_part),
     *     more_int)
c
c---------------------------------------------------------
c calculate the size of the request and split it if needed
c
      call getmem(nsupblk+2,isbl_split)
      call request_split(leri,inx,bl(iisd),bl(jjsd),nbl2,
     *     bl(ijpres2),bl(klpres2),bl(isblsize),bl(isblqrts),
     *     bl(isbl_part),bl(mxsize),
     *     bl(isbl_split),isplit )
c
      call retmem(1)
      call getmem(isplit+2, isbl_split)
c---------------------------------------------------------
c calculate the pointer to the block's quartets 
c (to eliminate calcul. in the get_nbls_pnl routine)
c
      call getmem(nsupblk,isblpoint)
      call supblks_point(nsupblk,bl(isblsize),bl(isblpoint))
c--------------------------------------------------------
      end
c==================================================================
      subroutine get_sizes(inx, icstxs,jcstxs,kcstxs,lcstxs,
     *                     isize,jsize,ksize,lsize,ngcont,ngcq)
      dimension inx(12,*),icstxs(*),jcstxs(*), kcstxs(*),lcstxs(*)
c
      ics=icstxs(1)
      jcs=jcstxs(1)
      ijcs=(ics-1)*ics/2 +jcs
      kcs=kcstxs(1)
      lcs=lcstxs(1)
      klcs=(kcs-1)*kcs/2 +lcs
c
      igcon=inx(4, ics) + 1
      jgcon=inx(4, jcs) + 1
      kgcon=inx(4, kcs) + 1
      lgcon=inx(4, lcs) + 1
c
      ngcont=igcon*jgcon*kgcon*lgcon
c
      ngcq=1
      if(ngcont.gt.1) then
         ijpg=igcon*jgcon
         if(ics.eq.jcs) ijpg=igcon*(igcon+1)/2
         klpg=kgcon*lgcon
         if(kcs.eq.lcs) klpg=kgcon*(kgcon+1)/2
         ngcq=ijpg*klpg
         if(ijcs.eq.klcs) ngcq=ijpg*(ijpg+1)/2
      endif
c
      isize=inx(3, ics)
      jsize=inx(3, jcs)
      ksize=inx(3, kcs)
      lsize=inx(3, lcs)
c
      end
c==================================================================
      subroutine swap_shells(nquart,icstxs,jcstxs,kcstxs,lcstxs,
     *                                         nbl2,
     *                       ijblock,klblock,ijpres2,klpres2,
     *                       isupblk, isbl_s,isbl_q, nblock2,
     *                       iqorder,isbl_p,
     *                               isbl_c)  
      dimension icstxs(nquart),jcstxs(nquart),
     *          kcstxs(nquart),lcstxs(nquart)
      dimension ijblock(nbl2),klblock(nbl2)
      dimension ijpres2(nquart),klpres2(nquart)
      dimension iqorder(nquart)
      dimension isupblk(nquart)
      dimension isbl_s(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
      dimension isbl_q(nquart)
      dimension nblock2(*) !  dimension ncs*(ncs+1)/2
      dimension ipnl(4)
      dimension isbl_p(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
      dimension isbl_c(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
c
c
      isbl=0
      do 05 ibl=1,nbl2
      ijblock(ibl)=0
      klblock(ibl)=0
         do 05 kbl=1,ibl
         isbl=isbl+1
         isbl_s(isbl)=0
   05 continue
c
      do 10 iq=1,nquart
        ipnl(1)=1   
        ipnl(2)=2   
        ipnl(3)=3   
        ipnl(4)=4   
c
        ics=icstxs(iq)
        jcs=jcstxs(iq)
        kcs=kcstxs(iq)
        lcs=lcstxs(iq)
c
        if(ics.lt.jcs) then
          ii=ics
          ics=jcs
          jcs=ii
          ipnl(1)=2
          ipnl(2)=1
        endif
        if(kcs.lt.lcs) then
          kk=kcs
          kcs=lcs
          lcs=kk
          ipnl(3)=4
          ipnl(4)=3
        endif
c
c now find ijcs and klcs pairs and find to which 
c pair-blocks they belong :
c
        ijcs=ics*(ics-1)/2+jcs
        klcs=kcs*(kcs-1)/2+lcs
        ijblock1=nblock2(ijcs)
        klblock1=nblock2(klcs)
c
        iswitch=0
        if(ijblock1.lt.klblock1) iswitch=1
        if(ijblock1.eq.klblock1 .and. ijcs.lt.klcs) iswitch=1
c
        if(   iswitch.eq.1     ) then
           ii=ics
           ics=kcs
           kcs=ii
           jj=jcs
           jcs=lcs
           lcs=jj
           ij=ijcs
           ijcs=klcs
           klcs=ij
             i1=ipnl(1)
             j1=ipnl(2)
             ipnl(1)=ipnl(3)
             ipnl(2)=ipnl(4)
             ipnl(3)=i1
             ipnl(4)=j1
             ijblk=ijblock1
             ijblock1=klblock1
             klblock1=ijblk
        endif
c
        ijblock(ijblock1)=1
        klblock(klblock1)=1
c
        icstxs(iq)=ics
        jcstxs(iq)=jcs
        kcstxs(iq)=kcs
        lcstxs(iq)=lcs
c
        ijpres2(iq)=ijcs
        klpres2(iq)=klcs
c   
        isupblk(iq)=ijblock1*(ijblock1-1)/2+klblock1
c
        iorder=1000*ipnl(1)+100*ipnl(2)+10*ipnl(3)+ipnl(4)
        iqorder(iq)=iorder
c
   10 continue
c
c---------------------------------------------------------
c calculate how many super-blocks are matched by a given
c PNL request and to what degree :
c
c
      do 30 iq=1,nquart
      isblq=isupblk(iq)
      isbl_s(isblq)=isbl_s(isblq)+1
  30  continue
c
      ipoint=0
      isbl_p(1)=0
      do 40 ikbl=2,nbl2*(nbl2+1)/2
      ipoint=ipoint+isbl_s(ikbl-1)
      isbl_p(ikbl)=ipoint
  40  continue
c
      do 50 iq=1,nquart
      isblq=isupblk(iq)
      isbl_p(isblq)=isbl_p(isblq)+1
      iquart=isbl_p(isblq)
      isbl_q(iquart)=iq
  50  continue
c---------------------------------------------------------
c For performace estimation 
c(how many times a given s-bl was called)
c
      do 90 isbl=1,nbl2*(nbl2+1)/2 
      if(isbl_s(isbl).gt.0) isbl_c(isbl)=isbl_c(isbl)+1
   90 continue
c---------------------------------------------------------
c
      end
c==================================================================
      subroutine ordr_shells(icstxs,jcstxs,kcstxs,lcstxs,
     *                       isize,jsize,ksize,lsize,itxspnl,
     *                       iqorder,ngcq)
      dimension icstxs(*),jcstxs(*),kcstxs(*),lcstxs(*)
      dimension iqorder(*),itxspnl(*)
c
c all dimensions are = nquart
c---------------------------------------------------------
c Only for a run without labels :
c
      iorder=iqorder(1)
c
         if(iorder.eq.1234) then
c                  ijkl->ijkl
            ls=0
            ks=lsize
            js=ksize*ks
            is=jsize*js    
            im=0
            jm=0
            km=0
            lm=1
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.1243) then
c                  ijkl->ijlk
            ls=ksize
            ks=0
            js=lsize*ls
            is=jsize*js
            im=0
            jm=0
            km=1
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.2134) then 
c                  ijkl->jikl
            ls=0    
            ks=lsize
            is=ksize*ks
            js=isize*is
            im=0
            jm=0
            km=0
            lm=1
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.2143) then 
c                  ijkl->jilk
            ls=ksize
            is=lsize*ls
            js=isize*is
            ks=0
            im=0
            jm=0
            km=1
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
c
         if(iorder.eq.3412) then
c                  ijkl->klij
            is=jsize
            ls=isize*is
            ks=lsize*ls
            js=0
            im=0
            jm=1
            km=0
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.4312) then 
c                  ijkl->klji
            is=0
            js=isize
            ls=jsize*js
            ks=lsize*ls
            im=1
            jm=0
            km=0
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.3421) then 
c                  ijkl->lkij
            js=0
            is=jsize
            ks=isize*is
            ls=ksize*ks
            im=0
            jm=1
            km=0
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
         if(iorder.eq.4321) then
c                  ijkl->lkji
            is=0
            js=isize
            ks=jsize*js
            ls=ksize*ks
            im=1
            jm=0
            km=0
            lm=0
            call txspnl(is,im, js,jm, ks,km ,ls,lm, itxspnl,ngcq)
         endif
c
      end
c======================================================================
      subroutine txs_pnl(nquart,ncshell, 
     *                   icstxs,jcstxs,kcstxs,lcstxs,
     *                   icspnl,jcspnl,kcspnl,lcspnl)
c
      dimension ncshell(*)
      dimension 
     * icspnl(nquart),jcspnl(nquart),kcspnl(nquart),lcspnl(nquart),
     * icstxs(nquart),jcstxs(nquart),kcstxs(nquart),lcstxs(nquart)
c
c      txs-order      pnl-order
c
      do 10 iq=1,nquart
        icstxs(iq)=ncshell(icspnl(iq))
        jcstxs(iq)=ncshell(jcspnl(iq))
        kcstxs(iq)=ncshell(kcspnl(iq))
        lcstxs(iq)=ncshell(lcspnl(iq))
   10 continue
c
      do 20 iq=1,nquart
        icspnl(iq)=icstxs(iq)
        jcspnl(iq)=jcstxs(iq)
        kcspnl(iq)=kcstxs(iq)
        lcspnl(iq)=lcstxs(iq)
   20 continue
c
      end
c=================================================================
      subroutine pnl_txs(ncfunct, icf,jcf,kcf,lcf,integ_n0 )
      dimension ncfunct(*),icf(*),jcf(*),kcf(*),lcf(*)
c
      do 10 integ=1,integ_n0
      itxs=icf(integ)
      jtxs=jcf(integ)
      ktxs=kcf(integ)
      ltxs=lcf(integ)
c
      ipnl=ncfunct(itxs)
      jpnl=ncfunct(jtxs)
      kpnl=ncfunct(ktxs)
      lpnl=ncfunct(ltxs)
c
      icf(integ)=ipnl
      jcf(integ)=jpnl
      kcf(integ)=kpnl
      lcf(integ)=lpnl
   10 continue
c
      end
c=================================================================
      subroutine txspnl(is,im, js,jm, ks,km ,ls,lm,itxspnl,ngcq)
      common /pnl004/ isize,jsize,ksize,lsize,iiiiiii
      dimension itxspnl(*)
c
c Temprorarly re-order integrals txs-pnl
c
        increm=isize*jsize*ksize*lsize
        itxs=0
        do 10 iqu=1,ngcq
        integ=(iqu-1)*increm
           do 20 i=1,isize
           ii=(i-1)*is +i*im 
           do 20 j=1,jsize
           jj=(j-1)*js +j*jm 
           ij=ii+jj
           do 20 k=1,ksize
           kk=(k-1)*ks +k*km
           ijk=ij+kk
           do 20 l=1,lsize
           ll=(l-1)*ls +l*lm
           itxs=itxs+1
ccccc      ipnl=ijk+ll
           ipnl=ijk+ll + integ
           itxspnl(itxs)=ipnl
   20      continue
   10   continue
c
c     write(6 ,*)'itxs , ipnl ( itxspnl(itxs))'
c      do 1111 ii=1,isize*jsize*ksize*lsize
c     write(6 ,66)ii,itxspnl(ii)
c1111  continue
c  66 format('ii=',i3,1x,'itxspnl(ii)=',i3)
c
      end
c=====================================================================
      subroutine switch_scr(bltxs, blscr,l_blscr)
      implicit real*8 (a-h,o-z)
c
c copy  from local bl() as defined in texas_face3 :
c
      common /memor1_R/ iisd_R,jjsd_R,ijbld_R
      common /memor1a_R/ npard_R,ncost_R,nsupb_R, mxsize_R,nasize_R
      common /memors_R/      ijshp_R
c
c to pnl :
c
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memors/ nsym,ijshp,isymm
c
c and for corresp. sizes :
c
      common /memor1_S/ iisd_S,jjsd_S,ijbld_S
      common /memor1a_S/ npard_S,ncost_S,nsupb_S, mxsize_S,nasize_S
      common /memors_S/      ijshp_S
c
c
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nXXX,nganz(35),lopt(30)
c
      dimension  bltxs(*)
      dimension  blscr(l_blscr)
c-------------------------------
cccc  write(6,*)' in switch_scr ; l_blscr=',l_blscr
c-------------------------------
      lcore=l_blscr
c
c now get/ret operate on blscr .
c
      call retall
c-------------------------------
c
      call getmem(iisd_S ,  iisd)
      call getmem(jjsd_S ,  jjsd)
      call getmem(ijbld_S, ijbld)
c
      call getmem(npard_S ,npard)
      call getmem(ncost_S ,ncost)
      call getmem(nsupb_S ,nsupb)
      call getmem(mxsize_S,mxsize)
      call getmem(nasize_S,nasize)
c
      call getmem(ijshp_S, ijshp)
c-------------------------------
c copy data :
c
      call tfer_i(bltxs(iisd_R), blscr(iisd), iisd_S )
      call tfer_i(bltxs(jjsd_R), blscr(jjsd), jjsd_S )
      call tfer_i(bltxs(ijbld_R), blscr(ijbld), ijbld_S )
c
      call tfer_i(bltxs(npard_R), blscr(npard), npard_S )
      call tfer_i(bltxs(ncost_R), blscr(ncost), ncost_S )
      call tfer_i(bltxs(nsupb_R), blscr(nsupb), nsupb_S )
      call tfer_i(bltxs(mxsize_R),blscr(mxsize), mxsize_S )
      call tfer_i(bltxs(nasize_R),blscr(nasize), nasize_S )
c
      call tfer_i(bltxs(ijshp_R), blscr(ijshp), ijshp_S )
c
c-------------------------------
      end
c==========
      subroutine tfer_i(ia,ib,n)
      dimension ia(n),ib(n)
c
      do 10 ii=1,n
      ib(ii)=ia(ii)
   10 continue
c
      end
c=====================================================================
      subroutine supblks_split(maxsize,leri,inx,iis,jjs,ijbl,nbl2,
     *                         ijpres2,klpres2,isbl_s,isbl_part,
     *                         more_int)
c--------------------------------------------------------------------
c this routine is called for every new PNL request i.e when
c more_int=.false.
c--------------------------------------------------------------------
      logical more_int
      common /pnl_nqrt/ ncall_pnl, nqrts_pnl
      dimension maxsize(*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
      dimension isbl_s(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
      dimension ijpres2(*),klpres2(*)  ! dimension nquart
      dimension isbl_part(*)  ! dimension nbl2*(nbl2+1)/2
c--------------------------------------------------------------------
         isupb=0
         do 10 ibl=1,nbl2
         ijcs1=ijbl(ibl,1)
         ics1=iis(ijcs1)
         jcs1=jjs(ijcs1)
         icont_ij=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         integ_ij=inx(3,ics1)*inx(3,jcs1)*icont_ij
         do 10 kbl=1,ibl
            isupb=isupb+1 
            isbl_part(isupb)=1
            isbl_size=isbl_s(isupb)
            if(isbl_size.gt.0) then
               maxm_size=maxsize(isupb)
               klcs1=ijbl(kbl,1)
               kcs1=iis(klcs1)
               lcs1=jjs(klcs1)
               icont_kl=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
               integ_kl=inx(3,kcs1)*inx(3,lcs1)*icont_kl
               integ1=integ_ij*integ_kl
               integ=integ1*maxm_size
c
c first re-define maximum size according to leri :
               if(integ.gt.leri) then
                  nqrts=leri/integ1
                  if(nqrts.eq.0) nqrts=1
                  maxsize(isupb)=nqrts
               endif
c
c split each super-block in parts if its current pnl-size is greater
c than maximum-size allowed :
c
               maxm_size=maxsize(isupb)
               if(isbl_size.gt.maxm_size) then
                  nparts=isbl_size/maxm_size
                  nrem=mod(isbl_size,maxm_size)
                  if(nrem.gt.0) nparts=nparts+1
                  isbl_part(isupb)=nparts
               endif
            endif
   10    continue
c---------------------------------------------------------
      end
c=====================================================================
      subroutine request_split(leri,inx,iis,jjs,nbl2,ijpres2,klpres2,
     *                         isbl_s,isbl_q, isbl_part,maxsize,
     *                         isbl_split,isplit)
c--------------------------------------------------------------------
c if a given pnl-request is too big then it is devided into parts:
c--------------------------------------------------------------------
      common /pnl003/ nqrtpnl,icstxs,jcstxs,kcstxs,lcstxs
      common /pnl_nqrt/ ncall_pnl, nqrts_pnl
      dimension inx(12,*)
      dimension iis(*),jjs(*)
      dimension isbl_s(*)               ! dimension nbl2*(nbl2+1)/2  
      dimension isbl_q(*)               ! dimension nquart
      dimension ijpres2(*),klpres2(*)   ! dimension nquart
      dimension isbl_part(*),maxsize(*) ! dim. nbl2*(nbl2+1)/2
      dimension isbl_split(0:*)
c-----------------------------------
c
         itotal=0
         isplit=0
         integrals=0
         ipoint=0
c
         do 60 ibl=1,nbl2
         iibl=ibl*(ibl-1)/2
         do 60 kbl=1,ibl
         isupb=iibl+kbl
         isbl_size=isbl_s(isupb)
         maxm_size=maxsize(isupb)
         nparts=isbl_part(isupb)
         if(isbl_size.eq.0) go to 60
c
           iq1=isbl_q(ipoint+1)
           ijcs1=ijpres2(iq1)
           klcs1=klpres2(iq1)
           ics1=iis(ijcs1)
           jcs1=jjs(ijcs1)
           icont_ij=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
           kcs1=iis(klcs1)
           lcs1=jjs(klcs1)
           icont_kl=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
           integ1=inx(3,ics1)*inx(3,jcs1)*inx(3,kcs1)*inx(3,lcs1)
           integ1=integ1*icont_ij*icont_kl
c
c  splitting among super-blocks :
c
           if(nparts.eq.1) then
              integ=integ1*isbl_size
              if((integrals+integ).gt.leri) then
                 itotal=itotal+1
                 isplit=isplit+1
                 integrals=integ
                 isbl_split(isplit)=isupb
              else
                 integrals=integrals+integ
              endif
           else 
c  splitting inside of one super-block :
c force splitting :
              integrals=leri
              itotal=itotal+nparts
              isplit=isplit+1
              isbl_split(isplit)=isupb
           endif
c
           ipoint=ipoint+isbl_size
   60    continue
c
         isblast=nbl2*(nbl2+1)/2
         isbl_split(0)=1
         isbl_split(isplit+1)=isblast+1
c
c---------------------------------------------------------
ccc   if( isplit.gt.0 ) then
ccc    if(ncall_pnl.eq.3) then
c         write(6,*)'FROM REQUEST_SPLIT : isplit=',isplit+1
c         write(6,*)'pnl-request no=',ncall_pnl,' size=',nqrtpnl,
c    *              ' split into ',itotal+1,' pieces'
ccc  *              ' split into ',isplit+1,' pieces'
c         write(6,66)(isbl_split(ii),ii=0,isplit+1)
c 66      format(' isbl_split(ii)=',i3)
c
c        write(6,*)' blocks in the request no=3'
c        do isupb=1,nbl2*(nbl2+1)/2
c        maxm_size=maxsize(isupb)
c        isbl_size=isbl_s(isupb)
c        nparts=isbl_part(isupb)
c        write(6,77)isupb,maxm_size,isbl_size,nparts
c 77     format('block=',i3,' maxsize=',i5,' pnlsize=',i5,' parts=',i5)
c        enddo
ccc    endif
ccc   endif
c---------------------------------------------------------
      end
c=====================================================================
      subroutine supblks_point(nsupblk,isbl_size,isbl_point)
      dimension isbl_size(*), isbl_point(*)
c
c calculates the pointer to block's quartets 
c    - iq=isbl_q(IPOINT+iqp)
c
         ipoint=0
         isbl_point(1)=0
         do 10 isbl=2,nsupblk
         ipoint=ipoint+isbl_size(isbl-1)
         isbl_point(isbl)=ipoint
   10    continue
c
      end
c=====================================================================
      subroutine memstat_pnl
      common /mem_pnl_scr/ nall_peak,mark_peak,mem_peak
c
      call memstat(nalloc_bl,nmark_bl,maxmem_bl,memtot_bl)
c
      if(nalloc_bl.gt.nall_peak) nall_peak=nalloc_bl
      if(nmark_bl .gt.mark_peak) mark_peak=nmark_bl 
      if(maxmem_bl.gt. mem_peak) mem_peak=maxmem_bl
c
      end
c=====================================================================
      subroutine request_update(icspnl,jcspnl,kcspnl,lcspnl,nquart,
     *                          num_bas_ij,num_bas_kl)
c
c---------------------------------------------------------------------
c I,J shells might belong to different basis set than K,L shells
c---------------------------------------------------------------------
c num_bas_ij , num_bas_kl = basis sets for IJ,KL .
c For 4-center (ordinary) two-electron integrals ALL four vectors are
c non-zero (just numbers of contracted shells).
c For 3-center two-el. integrals ONE of these vectors should be ZERO.
c For 2-center two-el. integrals TWO of these vectors should be ZERO.
c
c If ZERO is found then it is replaced by NCS-the LAST contracted shell
c (in PNL order). This is the S-type shell (uncontracted & exp.=zero)
c which was added to the basis set. 
c---------------------------------------------------------------------
c
      logical cent2,cent3,cent4
      common /what_was_calc/ cent2,cent3,cent4
c---------------------------------------------------------------------
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c---------------------------------------------------------------------
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
c
      dimension icspnl(nquart),jcspnl(nquart),
     *          kcspnl(nquart),lcspnl(nquart)
c---------------------------------------------------------------------
ctest
c        write(6,*)'ncs_bas_1-3=',ncs_bas_1,ncs_bas_2,ncs_bas_3
c        write(6,*)' entr iscpnl :', icspnl
c        write(6,*)' entr jscpnl :', jcspnl
c        write(6,*)' entr kscpnl :', kcspnl
c        write(6,*)' entr lscpnl :', lcspnl
ctest
c
      ics=icspnl(1)
      jcs=jcspnl(1)
      kcs=kcspnl(1)
      lcs=lcspnl(1)
c
      ncenters=0
      if(ics.gt.0) ncenters=ncenters+1
      if(jcs.gt.0) ncenters=ncenters+1
      if(kcs.gt.0) ncenters=ncenters+1
      if(lcs.gt.0) ncenters=ncenters+1
c
c     write(6,*) ncenters,'-center integrals are requested'
c
      if(ncenters.eq.4) then
         cent4=.true.
c
c        all possible cases 
c        -----------------
c        if(num_bas_ij.eq.num_bas_1) then
c        endif
         if(num_bas_ij.eq.num_bas_2) then
            call ncshl_update(icspnl,nquart,ncs_bas_1)
            call ncshl_update(jcspnl,nquart,ncs_bas_1)
         endif
         if(num_bas_ij.eq.num_bas_3) then
            call ncshl_update(icspnl,nquart,ncs_bas_2)
            call ncshl_update(jcspnl,nquart,ncs_bas_2)
         endif
c
c        if(num_bas_kl.eq.num_bas_1) then
c        endif
         if(num_bas_kl.eq.num_bas_2) then
            call ncshl_update(kcspnl,nquart,ncs_bas_1)
            call ncshl_update(lcspnl,nquart,ncs_bas_1)
         endif
         if(num_bas_kl.eq.num_bas_3) then
            call ncshl_update(kcspnl,nquart,ncs_bas_2)
            call ncshl_update(lcspnl,nquart,ncs_bas_2)
         endif
      endif
c
      if(ncenters.eq.3) then
         cent3=.true.
c
c        all possible cases 
c        -----------------
         if(ics.eq.0) then
               do ii=1,nquart
                  icspnl(ii)=ncs
               enddo
            if(num_bas_ij.eq.num_bas_2) then
               call ncshl_update(jcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(jcspnl,nquart,ncs_bas_2)
            endif
            if(num_bas_kl.eq.num_bas_2) then
               call ncshl_update(kcspnl,nquart,ncs_bas_1)
               call ncshl_update(lcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(kcspnl,nquart,ncs_bas_2)
               call ncshl_update(lcspnl,nquart,ncs_bas_2)
            endif
         endif
c        -----------------
         if(jcs.eq.0) then
               do ii=1,nquart
                  jcspnl(ii)=ncs
               enddo
            if(num_bas_ij.eq.num_bas_2) then
               call ncshl_update(icspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(icspnl,nquart,ncs_bas_2)
            endif
            if(num_bas_kl.eq.num_bas_2) then
               call ncshl_update(kcspnl,nquart,ncs_bas_1)
               call ncshl_update(lcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(kcspnl,nquart,ncs_bas_2)
               call ncshl_update(lcspnl,nquart,ncs_bas_2)
            endif
         endif
c        -----------------
         if(kcs.eq.0) then
               do ii=1,nquart
                  kcspnl(ii)=ncs
               enddo
            if(num_bas_ij.eq.num_bas_2) then
               call ncshl_update(icspnl,nquart,ncs_bas_1)
               call ncshl_update(jcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(icspnl,nquart,ncs_bas_2)
               call ncshl_update(jcspnl,nquart,ncs_bas_2)
            endif
            if(num_bas_kl.eq.num_bas_2) then
               call ncshl_update(lcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(lcspnl,nquart,ncs_bas_2)
            endif
         endif
c        -----------------
         if(lcs.eq.0) then
               do ii=1,nquart
                  lcspnl(ii)=ncs
               enddo
            if(num_bas_ij.eq.num_bas_2) then
               call ncshl_update(icspnl,nquart,ncs_bas_1)
               call ncshl_update(jcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(icspnl,nquart,ncs_bas_2)
               call ncshl_update(jcspnl,nquart,ncs_bas_2)
            endif
            if(num_bas_kl.eq.num_bas_2) then
               call ncshl_update(kcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(kcspnl,nquart,ncs_bas_2)
            endif
         endif
      endif
c
      if(ncenters.eq.2) then
         cent2=.true.
c
c        only j=0 & l=0 case !!!
c        -----------------
               do ii=1,nquart
                  jcspnl(ii)=ncs
                  lcspnl(ii)=ncs
               enddo
            if(num_bas_ij.eq.num_bas_2) then
               call ncshl_update(icspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_ij.eq.num_bas_3) then
               call ncshl_update(icspnl,nquart,ncs_bas_2)
            endif
            if(num_bas_kl.eq.num_bas_2) then
               call ncshl_update(kcspnl,nquart,ncs_bas_1)
            endif
            if(num_bas_kl.eq.num_bas_3) then
               call ncshl_update(kcspnl,nquart,ncs_bas_2)
            endif
ctest
c        write(6,*)' exit iscpnl :', icspnl
c        write(6,*)' exit jscpnl :', jcspnl
c        write(6,*)' exit kscpnl :', kcspnl
c        write(6,*)' exit lscpnl :', lcspnl
ctest
      endif
c---------------------------------------------------------------------
      end
c=====================================================================
      subroutine ncshl_update(icspnl,nquart,ncs_bas_1)
      dimension icspnl(nquart)
c
            do ii=1,nquart
               icspnl(ii)=icspnl(ii)+ncs_bas_1
            enddo
c
      end
c=====================================================================
      subroutine labels_update(icf,jcf,kcf,lcf,integ,ij_basis,kl_basis,
     *                         ish_first,jsh_first,ksh_first,lsh_first)
c---------------------------------------------------------------------
      common /multi_basis/ num_bas_1,num_bas_2,num_bas_3,
     *                     ncs_bas_1,ncs_bas_2,ncs_bas_3,
     *                     nps_bas_1,nps_bas_2,nps_bas_3,
     *                     nat_bas_1,nat_bas_2,nat_bas_3,
     *                     ncf_bas_1,ncf_bas_2,ncf_bas_3 
      dimension icf(*),jcf(*),kcf(*),lcf(*)
c---------------------------------------------------------------------
ctest
c        write(6,*)'ncf_bas_1-3=',ncf_bas_1,ncf_bas_2,ncf_bas_3
c        write(6,*)' entr icf1-9 :',(icf(ii),ii=1,9)
c        write(6,*)' entr jcf1-9 :',(jcf(ii),ii=1,9)
c        write(6,*)' entr kcf1-9 :',(kcf(ii),ii=1,9)
c        write(6,*)' entr lcf1-9 :',(lcf(ii),ii=1,9)
ctest
c
      if(ij_basis.eq.num_bas_1 .and. kl_basis.eq.num_bas_1) RETURN
c
c ij
c
      if(ij_basis.eq.num_bas_2) then
        if(ish_first.gt.0) call ncfun_update(icf,integ,ncf_bas_1)
        if(jsh_first.gt.0) call ncfun_update(jcf,integ,ncf_bas_1)
      endif
c
      if(ij_basis.eq.num_bas_3) then
        if(ish_first.gt.0) call ncfun_update(icf,integ,ncf_bas_2)
        if(jsh_first.gt.0) call ncfun_update(jcf,integ,ncf_bas_2)
      endif
c
c kl
c
      if(kl_basis.eq.num_bas_2) then
        if(ksh_first.gt.0) call ncfun_update(kcf,integ,ncf_bas_1)
        if(lsh_first.gt.0) call ncfun_update(lcf,integ,ncf_bas_1)
      endif
c
      if(kl_basis.eq.num_bas_3) then
        if(ksh_first.gt.0) call ncfun_update(kcf,integ,ncf_bas_2)
        if(lsh_first.gt.0) call ncfun_update(lcf,integ,ncf_bas_2)
      endif
ctest
c        write(6,*)' exit icf1-9 :',(icf(ii),ii=1,9)
c        write(6,*)' exit jcf1-9 :',(jcf(ii),ii=1,9)
c        write(6,*)' exit kcf1-9 :',(kcf(ii),ii=1,9)
c        write(6,*)' exit lcf1-9 :',(lcf(ii),ii=1,9)
ctest
c
      end
c=====================================================================
      subroutine ncfun_update(icf,integ,ncf_bas_1)
      dimension icf(*)
c
         do ii=1,integ
         icf(ii)=icf(ii)-ncf_bas_1
         enddo
c
      end
c=====================================================================
