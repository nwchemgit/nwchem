      subroutine texas_hf2_m(ij_basis,ish,jsh,kl_basis,ksh,lsh,nquart,
     *                       q4, use_q4,
     *                       ra,rb,rc,rd,use_r,
     *                       eri,leri,icf,jcf,kcf,lcf,integ_n0,labels,
     *                       more_int,blscr,l_blscr,zerotol,int_type)
*
* $Id$
*
      implicit none
#include "mafdecls.fh"
#include "bas.fh"
c
      character*8 int_type
      integer ij_basis, ish(*), jsh(*), kl_basis, ksh(*), lsh(*)
      integer nquart
      double precision q4(*)
      logical use_q4
      double precision ra(*), rb(*), rc(*), rd(*)
      logical use_r
      double precision eri(*)
      integer leri
      integer icf(*), jcf(*), kcf(*), lcf(*)
      integer integ_n0
      logical labels
      logical more_int
      double precision blscr(*)
      integer l_blscr
      double precision zerotol
c
      double precision time_beg,time_end
      double precision time4texas_hf2_m,time4mul_quart
      common /pnl_time/ time4texas_hf2_m,time4mul_quart
c------------------------------------------------------------------
      call txs_second(time_beg)
c------------------------------------------------------------------
c send info about type of integrals requested :
c
      call requested_task(int_type)
c------------------------------------------------------------------
         call texas_hf2_m2(ij_basis,ish,jsh,kl_basis,ksh,lsh,nquart,
     *        q4, use_q4,
     *        ra,rb,rc,rd,use_r,
     *        eri,leri,icf,jcf,kcf,lcf,integ_n0,labels,
     *        more_int,blscr,l_blscr,zerotol)
c------------------------------------------------------------------
      call txs_second(time_end)
      time4texas_hf2_m=time4texas_hf2_m + (time_end-time_beg)
c------------------------------------------------------------------
      end
c=================================================================
      subroutine texas_hf2_m2(ij_basis,ish,jsh,kl_basis,ksh,lsh,nquart,
     *                       q4, use_q4,
     *                       ra,rb,rc,rd,use_r,
     *                       eri,leri,icf,jcf,kcf,lcf,integ_n0,labels,
     *                       more_int,blscr,l_blscr,zerotol)
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
#include "errquit.fh"
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
c screening threshold for output integrals
      double precision zerotol
c----------------------------------------------------------------------
      double precision savezerotol
      common /csavezerotol/ savezerotol ! Used in detbul
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
      savezerotol = zerotol
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
     &   ('texas_hf2_m: called with different basis set handle',911,
     &       BASIS_ERR)
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
     &   ('texas_hf2_m: called with different basis set handle',911,
     &       BASIS_ERR)
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
      common /ctxs_index/ maxsh,ifp,inx(1)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c----------------------------------------------------------------
cccc  common /memmax/ ispblx, maxme1,iforwhat
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
c----------------------------------------------------------------
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
c----------------------------------------------------------------
      common /pnl006/ nsplit,isplit,isbl_split, isbl_part
c----------------------------------------------------------------
      common /pnl_time/ time4texas_hf2_m,time4mul_quart
c----------------------------------------------------------------
      common /pnl_nqrt/ ncall_pnl,ncall_41,nqrts_pnl,npart_pnl,nsize_pnl
c----------------------------------------------------------------
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
             call txs_second(time_beg)
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
      if(nquart.eq.1) ncall_41=ncall_41 + 1
ctest
c     write(6,*)
c    * ' call no=',ncall_pnl,' nquart=',nquart,' nq_tot=',nqrts_pnl
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
      call txs_setup(bl,inx,ncs,leri,
     *               nquart,icspnl,jcspnl,kcspnl,lcspnl,labels)
ctest
c       call texas_terminate()
c       stop 'stopped by kw after txs_setup'
ctest
c
c 11 calls to getmem 
c 12 calls to getmem
c if(.not.labels) +1 call to getmem
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
      npart_pnl=npart_pnl+1
c
      if(.not.labels) then
         ispec=1
         integ_n0=0
c
c        calculate integrals without labels :
c
c        ncfunct is not used on this call
         call calcint2(bl,inx, q4,use_q4, icf,jcf,kcf,lcf,eri,
     *        dbl_mb(ncfunct))
         integ_nx=integ_n0
c
      else
c
         ispec=2
         integ_n0=0
c
c        calculate integrals with labels :
c
         call calcint2(bl,inx, q4,use_q4, icf,jcf,kcf,lcf,eri,
     *        dbl_mb(ncfunct))
         integ_nx=integ_n0
c
      endif
c
      if(isplit.gt.0) then
         if(nsplit.eq.isplit+1) more_int=.false.
      endif
c
      if(.not.more_int) then
         if(.not.labels) then
c10         call retmem(13)
            call retmem(11)
         else
c10         call retmem(12)
            call retmem(10)
         endif
c----------------------------------------------------------------
c if use_r is true return original coordinates 
c Probably NOT needed !
c---
c        if (use_r) then
c         call replacer(bl(inuc+1),ra,rb,rc,rd,inx,nquart,
c    *                  icspnl,jcspnl,kcspnl,lcspnl,'back')
c        endif
c----------------------------------------------------------------
c check and update memory status in pnl_scratch 
c after each pnl request is done; 
c
         call memstat_pnl
c
      endif
c----------------------------------------------------------------
      call txs_second(time_end)
      time4mul_quart=time4mul_quart  + (time_end -time_beg)
c----------------------------------------------------------------
      if(integ_nx.eq.0 .and. more_int) go to 999
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
     *                     nquart,icspnl,jcspnl,kcspnl,lcspnl,labels)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical labels
c
      common /bl_txs_add/ ntxs_bl_scr 
c
      common /memor1/ npard,mxsize,nblock1,nblock1_back
      common /memor11/ mxpair
      common /memor1b/ nbl2
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /pnl003/ nqrtpnl,icstxs,jcstxs,kcstxs,lcstxs
      common /pnl004/ isize,jsize,ksize,lsize,itxspnl
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /pnl006/ nsplit,isplit,isbl_split, isbl_part
      dimension icspnl(*),jcspnl(*),kcspnl(*),lcspnl(*)  ! dim=nquart
      dimension inx(12,*)
      dimension bl(*)
c--------------------------------------------------------
      nqrtpnl=nquart
      nsupblk=nbl2*(nbl2+1)/2
c--------------------------------------------------------
c find texas-shells corresponding to the given set of pnl-shells.
c
      call memo1_int(nquart,icstxs)
      call memo1_int(nquart,jcstxs)
      call memo1_int(nquart,kcstxs)
      call memo1_int(nquart,lcstxs)
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
      call memo1_int(nquart, ijpres2)
      call memo1_int(nquart, klpres2)
c-
      call memo1_int(nquart, iqorder)
c--------------------------------------------------------
c reserve memory 
c for pair-blocks ij, kl which are present in a request
c
      call memo1_int(nbl2, ijblock)
      call memo1_int(nbl2, klblock)
c--------------------------------------------------------
c reserve memory for present super-blocks
c
      call memo1_int(nquart,isblqrts)
      call memo1_int(nsupblk,isblsize)
c--------------------------------------------------------
c  7 calls to getmem 
c--------------------------------------------------------
      call memo1_int(nquart,isupblk)
      call memo1_int(nsupblk,isblpoin)
*         call txs_second(time_beg)
      call swap_shells(nquart,icspnl,jcspnl,kcspnl,lcspnl,nbl2,
     *     bl(ijblock),bl(klblock),bl(ijpres2),bl(klpres2),
     *     bl(isupblk),bl(isblsize),bl(isblqrts),bl(nblock1_back),
     *     bl(iqorder),bl(isblpoin) )
         call txs_second(time_end)
*         time_swap=time_swap +(time_end-time_beg)
      call retmem(1)
      call retmem(1)
c  
      if(.not.labels) then
         call get_sizes(inx,icspnl,jcspnl,kcspnl,lcspnl,
     *                      isize,jsize,ksize,lsize,ngcont)
ccc   write(6,*)'general contraction deep=',ngcont
c
         call memo1_int(ngcont*isize*jsize*ksize*lsize,itxspnl)
c
         call ordr_shells(icspnl,jcspnl,kcspnl,lcspnl,
     *                    isize,jsize,ksize,lsize,bl(itxspnl),
     *                    bl(iqorder),ngcont)
c       1 call to getmem
      endif
c---------------------------------------------------------
c calculate the pointer to the block's quartets 
c (to eliminate calcul. in the get_nbls_pnl routine)
c
      call memo1_int(nsupblk,isblpoint)
      call supblks_point(nsupblk,bl(isblsize),bl(isblpoint))
c--------------------------------------------------------
c check current sizes of super-blocks against maximum sizes
c and available size of integral's buffer (leri):
c (maximum sizes calculated in blksizer (spec_block.f) )
c
      call memo1_int(nsupblk,isbl_part)
      call supblks_split(bl,nquart,bl(mxsize),bl(mxpair),
     *                   bl(ijblock),bl(klblock),
     *                   leri,inx, bl(nblock1),nbl2,
     *                   bl(isblqrts),bl(isblpoint),
     *        bl(ijpres2),bl(klpres2),bl(isblsize),bl(isbl_part))
c---------------------------------------------------------
c calculate the size of the request and split it if needed
c
      call memo1_int(nsupblk+2,isbl_split)
      call request_split(leri,inx, nbl2,bl(ijblock),bl(klblock),
     *     bl(ijpres2),bl(klpres2),bl(isblsize),bl(isblqrts),
     *     bl(isbl_part),
     *     bl(isbl_split),isplit )
c
      call retmem(1)
      call memo1_int(isplit+2, isbl_split)
c--------------------------------------------------------
c2002 : efficeiency of using blocking capabilities of Texas:
c
      nquartot=ncs*(ncs+1)/2
      nquartot=nquartot*(nquartot+1)/2
      call blk_capab(nbl2,bl(npard),bl(isblsize),nquart,nquartot)
c--------------------------------------------------------
      end
c==================================================================
      subroutine get_sizes(inx, icstxs,jcstxs,kcstxs,lcstxs,
     *                     isize,jsize,ksize,lsize,ngcont)
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
     *                       isupblk, isbl_s,isbl_q, nblock1_back,
     *                       iqorder,isbl_p)
      dimension icstxs(nquart),jcstxs(nquart),
     *          kcstxs(nquart),lcstxs(nquart)
      dimension ijblock(nbl2),klblock(nbl2)
      dimension ijpres2(nquart),klpres2(nquart)
      dimension iqorder(nquart)
      dimension isupblk(nquart)
      dimension isbl_s(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
      dimension isbl_q(nquart)
      dimension nblock1_back(*) !  dimension ncs
      dimension ipnl(4)
      dimension isbl_p(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
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
cnew
        i_block1=nblock1_back(ics)
        j_block1=nblock1_back(jcs)
        ijblock1=i_block1*(i_block1-1)/2 + j_block1
        k_block1=nblock1_back(kcs)
        l_block1=nblock1_back(lcs)
        klblock1=k_block1*(k_block1-1)/2 + l_block1
cnew
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
ctest98
c       itxs=icstxs(iq)
c       ipnl=icspnl(iq)
c       jtxs=jcstxs(iq)
c       jpnl=jcspnl(iq)
c       ktxs=kcstxs(iq)
c       kpnl=kcspnl(iq)
c       ltxs=lcstxs(iq)
c       lpnl=ncshell(lcspnl(iq))
c       write(6,*)' PNL shells =',ipnl,jpnl,kpnl,lpnl
c       write(6,*)' TXS shells =',itxs,jtxs,ktxs,ltxs
ctest98
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
c copy  from local bl() as defined in texas_face :
c
      common /memor1_R/ npard_R,mxsize_R,nblock1_R,nblock1_back_R
      common /memor11_R/ mxpair_R
c
c to pnl :
c
      common /memor1/ npard,mxsize,nblock1,nblock1_back
      common /memor11/ mxpair
c
c and for corresp. sizes :
c
      common /memor1_S/ npard_S,mxsize_S,nblock1_S,nblock1_back_S
      common /memor11_S/ mxpair_S
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
      call memo1_int(npard_S ,npard)
      call memo1_int(mxsize_S,mxsize)
      call memo1_int(mxpair_S,mxpair)
      call memo1_int(nblock1_S,nblock1)
      call memo1_int(nblock1_back_S,nblock1_back)
c-------------------------------
c copy data :
c
      call tfer_i(bltxs(npard_R), blscr(npard), npard_S )
      call tfer_i(bltxs(mxsize_R),blscr(mxsize), mxsize_S )
      call tfer_i(bltxs(mxpair_R),blscr(mxpair), mxpair_S )
      call tfer_i(bltxs(nblock1_R),blscr(nblock1), nblock1_S )
      call tfer_i
     *    (bltxs(nblock1_back_R),blscr(nblock1_back), nblock1_back_S )
c-------------------------------
      end
c==========
      subroutine tfer_i(ia,ib,n)
      dimension ia(n),ib(n)
c
       do 10 ii=1,n
         ib(ii)=ia(ii)          ! let the compiler do the unrolling.
  10   continue
c
c$$$      m = mod(n,7)
c$$$      if( m .eq. 0 ) go to 20
c$$$        do 10 ii=1,n
c$$$          ib(ii)=ia(ii)
c$$$   10   continue
c$$$      if( n .lt. 7) return
c$$$   20 mp1=m+1
c$$$      do 30 i=mp1,n,7
c$$$         ib(i  ) = ia(i  )
c$$$         ib(i+1) = ia(i+1)
c$$$         ib(i+2) = ia(i+2)
c$$$         ib(i+3) = ia(i+3)
c$$$         ib(i+4) = ia(i+4)
c$$$         ib(i+5) = ia(i+5)
c$$$         ib(i+6) = ia(i+6)
c$$$   30 continue
c
      end
c=====================================================================
      subroutine supblks_split(bl,nquart,maxsize,maxpair,
     *                         ijblock,klblock,
     *                         leri,inx, nblock1,nbl2,
     *                         isbl_q,isbl_point,
     *                         ijpres2,klpres2,isbl_s,isbl_part)
c--------------------------------------------------------------------
c this routine is called for every new PNL request i.e when
c more_int=.false.
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
cccc  common /intlim/ limxmem,limblks,limpair
c
      common /pnl_nqrt/ ncall_pnl,ncall_41,nqrts_pnl,npart_pnl,nsize_pnl
c
c output from uniq_pairs_1 call:
c
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,    
     *                  nkl_uniqe, kl_uniqe_p, map_kl 
c
      dimension maxsize(*),maxpair(*)
      dimension inx(12,*),nblock1(*)
      dimension isbl_s(*) ! dimension nbl2*(nbl2+1)/2  num of super-blks
      dimension ijpres2(*),klpres2(*)  ! dimension nquart
      dimension isbl_part(*)  ! dimension nbl2*(nbl2+1)/2
      dimension isbl_q(*), isbl_point(*)
c
      dimension ijblock(nbl2),klblock(nbl2)
      dimension bl(*)
c---------------------------------------------------------------------
c for derivatives of all kinds :
cccc  common /memmax/ ispblx, maxme1,iforwhat
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
c---------------------------------------------------------------------
      n_times=1
      if(iforwhat.eq.2) n_times= 6 ! giao Ist derivatives
      if(iforwhat.eq.3) n_times=12 ! gradient derivatives
      if(iforwhat.eq.4) n_times=78 ! hessian  derivatives
c---------------------------------------------------------------------
ctest
c     write(6,*)'call=',ncall_pnl,' leri=',leri,' nquart=',nquart
c---------------------------------------------------------------------
         i_blks=0
         isupb=0
         do 10 ibl=1,nbl2
         if(ijblock(ibl).eq.0) then
            iibl=ibl*(ibl-1)/2
            iibl1=iibl+1
            iibl2=iibl+ibl
            do iiii=iibl1,iibl2
               isbl_part(iiii)=1
            enddo
            isupb=isupb+ibl
            go to 10
         endif
         call get_ics_jcs(nblock1,ibl,1,ics1,jcs1)
         icont_ij=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
         integ_ij=inx(3,ics1)*inx(3,jcs1)*icont_ij
            do 20 kbl=1,ibl
            isupb=isupb+1 
            isbl_part(isupb)=1
            isbl_size=isbl_s(isupb)
            if(isbl_size.gt.0) then
c
               call get_ics_jcs(nblock1,kbl,1,kcs1,lcs1)
               icont_kl=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
               integ_kl=inx(3,kcs1)*inx(3,lcs1)*icont_kl
               integ1=integ_ij*integ_kl
c
ccccccc        integ1=integ1*n_times  ! derivatives
c
c first re-define maximum size according to leri :
c
               maxm_size=maxsize(isupb)
               leri_size=leri/integ1
c
               if(leri_size.eq.0) then
c
                  write(6,*)' buffer size too small; leri=',leri,
     *            ' needed=',integ1,' =',n_times,'*',integ1/n_times
                  write(6,*)' shell sizes=',
     *            inx(3,ics1),inx(3,jcs1),inx(3,kcs1),inx(3,lcs1),
     *            ', shell g.con=',inx(4,ics1)+1,inx(4,jcs1)+1,
     *                            inx(4,kcs1)+1,inx(4,lcs1)+1
                  call errquit('texas: supblks_split',1016, INT_ERR)
               endif
c
               maxsize(isupb)=min(leri_size,maxm_size)
c
c split each super-block in parts if its current pnl-size is greater 
c than maximum-size allowed or number of pairs exceeds allowed limit
c
               maxm_size=maxsize(isupb)
               maxm_pair=maxpair(isupb)
c
               ipoint=isbl_point(isupb)
c
               call uniq_pairs_1('ks_split',ibl,kbl,bl,ijpres2,klpres2, 
     *                           isbl_size,isbl_q,ipoint)
               call retmem(2)    ! release allocations in uniq_pairs_1
c              output nij_uniqe & nkl_uniqe in  common /map4uniq/ 
c
              if(nij_uniqe.le.maxm_pair.and.nkl_uniqe.le.maxm_pair) then
                  if(isbl_size.gt.maxm_size) then
                     nparts=isbl_size/maxm_size
                     nrem  =mod(isbl_size,maxm_size)
                     if(nrem.gt.0) nparts=nparts+1
                     isbl_part(isupb)=nparts
                     maxsize(isupb)=maxm_size
                  endif
              else
                  maxm_size=min(maxm_size,maxm_pair)
                  nparts=isbl_size/maxm_size
                  nrem  =mod(isbl_size,maxm_size)
                  if(nrem.gt.0) nparts=nparts+1
                  isbl_part(isupb)=nparts
                  maxsize(isupb)=maxm_size
              endif
c
              i_blks=i_blks+isbl_part(isupb)
c
            endif
   20       continue
   10    continue
c---------------------------------------------------------
      if(i_blks.gt.0) nsize_pnl=nsize_pnl+ nquart/i_blks
c---------------------------------------------------------
      end
c=====================================================================
      subroutine request_split(leri,inx,nbl2, ijblock,klblock,
     *                         ijpres2,klpres2,
     *                         isbl_s,isbl_q, isbl_part,
     *                         isbl_split,isplit)
c--------------------------------------------------------------------
c if a given pnl-request is too big then it is divided into parts:
c--------------------------------------------------------------------
      dimension inx(12,*)
      dimension isbl_s(*)               ! dimension nbl2*(nbl2+1)/2  
      dimension isbl_q(*)               ! dimension nquart
      dimension ijpres2(*),klpres2(*)   ! dimension nquart
      dimension ijblock(nbl2),klblock(nbl2)
      dimension isbl_part(*)            ! dim. nbl2*(nbl2+1)/2
      dimension isbl_split(0:*)
c---------------------------------------------------------------------
c for derivatives of all kinds :
ccc   common /memmax/ ispblx, maxme1,iforwhat
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
c---------------------------------------------------------------------
      n_times=1
      if(iforwhat.eq.2) n_times= 6 ! giao Ist derivatives
      if(iforwhat.eq.3) n_times=12 ! gradient derivatives
      if(iforwhat.eq.4) n_times=78 ! hessian  derivatives
c---------------------------------------------------------------------
c
         isplit=0
         integrals=0
         ipoint=0
c
         do 50 ibl=1,nbl2
         if(ijblock(ibl).eq.0) go to 50
         iibl=ibl*(ibl-1)/2
            do 60 kbl=1,ibl
            isupb=iibl+kbl
            isbl_size=isbl_s(isupb)
            nparts=isbl_part(isupb)
            if(isbl_size.eq.0) go to 60
ctest
ccc      write(6,*)'block=',isupb,' size=',isbl_size,' nparts=',nparts
ctest
c
            iq1=isbl_q(ipoint+1)
            ijcs1=ijpres2(iq1)
            klcs1=klpres2(iq1)
            call get_ij_half(ijcs1,ics1,jcs1)
            icont_ij=(inx(4,ics1)+1)*(inx(4,jcs1)+1)
            call get_ij_half(klcs1,kcs1,lcs1)
            icont_kl=(inx(4,kcs1)+1)*(inx(4,lcs1)+1)
            integ1=inx(3,ics1)*inx(3,jcs1)*inx(3,kcs1)*inx(3,lcs1)
            integ1=integ1*icont_ij*icont_kl
c
cccccc      integ1=integ1*n_times ! derivatives
c
c  splitting among super-blocks :
c
            if(nparts.eq.1) then
               integ=integ1*isbl_size
               if((integrals+integ).gt.leri) then
                  isplit=isplit+1
                  integrals=integ
                  isbl_split(isplit)=isupb
               else
                  integrals=integrals+integ
               endif
            else 
c  splitting inside of one super-block : force splitting :
               integrals=leri
               isplit=isplit+1
               isbl_split(isplit)=isupb
            endif
c
            ipoint=ipoint+isbl_size
   60       continue
   50    continue
c
         isblast=nbl2*(nbl2+1)/2
         isbl_split(0)=1
         isbl_split(isplit+1)=isblast+1
c---------------------------------------------------------
      end
c=====================================================================
      subroutine supblks_point(nsupblk,isbl_size,isbl_point)
      dimension isbl_size(*), isbl_point(*)
c
c calculates the pointer to block's quartets : iq=isbl_q(IPOINT+iqp)
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
            if(num_bas_kl.eq.num_bas_3) then
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
            if(num_bas_kl.eq.num_bas_3) then
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
            if(num_bas_kl.eq.num_bas_3) then
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
            if(num_bas_kl.eq.num_bas_3) then
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
      subroutine requested_task(int_type)
      character*8 int_type
      character*27 type_i,type_r
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
      common /type_inited/ iforinit  ! use here & in texas_face.F
      common /mem_max_min/ ispblx,maxme1,max_111,iforwhat
c
      irequest = 1                   ! take care of compiler warnings
      if(int_type.eq.'scfd_int') then 
         where='buff'
         irequest=1
      endif
      if(int_type.eq.'giao_int') then
         where='shif'
         irequest=2
      endif
      if(int_type.eq.'der1_int') then 
         where='forc'
         irequest=3
      endif
      if(int_type.eq.'der2_int') then 
         where='hess'
         irequest=4
      endif
c
c check if the texas integral program has been initiated 
c for up to this requested task :
c
      if(irequest.gt.iforinit) then
         type_i =                 ' texas not initialized '
         if(iforinit.eq.1) type_i='ordinary 2-el.integrals    '
         if(iforinit.eq.2) type_i='giao integral derivatives  '
         if(iforinit.eq.3) type_i='first  geometry derivatives'
         if(iforinit.eq.4) type_i='second geometry derivatives'
c
         type_r =                 ' texas error on irequest'
         if(irequest.eq.1) type_r='ordinary 2-el.integrals    '
         if(irequest.eq.2) type_r='giao integral derivatives  '
         if(irequest.eq.3) type_r='first  geometry derivatives'
         if(irequest.eq.4) type_r='second geometry derivatives'
c
         write(6,66) type_i,type_r
   66    format('                   Texas integral program '/
     &       ' has been initialized only for :',a27,/,
     &       ' and can not be executed for   :',a27)
         call errquit('texas program stopped in requested_task',0,
     &       INT_ERR)
      else
         iforwhat=irequest
      endif
c
      end
c==================================================================
      subroutine get_ics_jcs(nblock1,ijblock,ijpar,ics,jcs)
c--------------------------------------------------
c ijblock - pair-block               ![input]
c ijpar   - a pair in that block     ![input]
c ics,jcs : shells                   ![output]
c--------------------------------------------------
      dimension nblock1(0:*   )      !    nblock1(0:nbl1)
c
      call get_ij_half(ijblock, iblock,jblock)
c
      ics_b= nblock1(iblock-1)+1     ! first shell 
      ics_e= nblock1(iblock)         ! last shell
      jcs_b= nblock1(jblock-1)+1     ! first shell 
      jcs_e= nblock1(jblock)         ! last shell
c
c     increm=ics_e-ics_b+1
c     jncrem=jcs_e-jcs_b+1
c
      if(jblock.eq.iblock) then
         call get_ij_half(ijpar,ishell,jshell)
      else
         jncrem=jcs_e-jcs_b+1
         call get_ij_full(ijpar,jncrem,ishell,jshell)
      endif
c
c ishell,jshell - shells in the ijblock (numbered from 1 to increm
c                                                      1 to jncrem)
c
      ics=ishell+ics_b -1
      jcs=jshell+jcs_b -1
c
      end
c==================================================================
      subroutine blk_capab(nbl2,npar,isbl_s,nquart_pnl,nquart_txs)
      implicit real*8 (a-h,o-z)
      common /pnl000/ xbluse,nbluse
      dimension isbl_s(*) , npar(*)
c
      xquart_pnl=dble(nquart_pnl)
      xquart_txs=dble(nquart_txs)
c
      ibluse=0
      ikbl=0
      do ibl=1,nbl2
         do kbl=1,ibl
            ikbl=ikbl+1
            iquart_pnl=isbl_s(ikbl)
            if(iquart_pnl.gt.0) ibluse=ibluse+1
         enddo
      enddo
c
      nbluse=nbluse+ibluse
c
      xblk_pnl=dble(ibluse)
      xblk_txs=dble(ikbl)
c
      size_pnl=xquart_pnl/xblk_pnl
      size_txs=xquart_txs/xblk_txs
      bluse=size_pnl/size_txs
      
      xbluse=xbluse+bluse
c
      end
c==================================================================
