c===============================================================
      subroutine blockin2(bl,lcore,inx,nbl2,nbloks,iret_limit)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
c
      common /route/ iroute
c
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /memor1/ iisd,jjsd,ijbld
cpnl: common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memors/ nsym,ijshp,isymm
c 
c save sizes :
      common /memor1_S/ iisd_S,jjsd_S,ijbld_S
      common /memor1a_S/ npard_S,ncost_S,nsupb_S, mxsize_S,nasize_S
      common /memors_S/      ijshp_S
c
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
c begining of the blocking procedure for pairs
c
      call memo1(ncs,ijchd)
c
c------------------------------------------------
c Check out how many atoms of the same type the 
c molecule contains and decide which blocking
c procedure will be used :
c
      if(iroute.eq.0) call whichblk(dbl_mb(inuc),na,inx,ncs,iroute)
c 
c output iroute (=1 or =2) ( stored in the common /route/
c
      if(nprint.gt.0) write(8,*)' blocking procedure route=',iroute
c------------------------------------------------
c Constructe blocks of contracted shell pairs
c ( calculates number of block-pairs nbl2 and 
c  reserves memory for npar, ijbl and nsupb 
c  and sets up these arrays )
c
      call blkpair(bl,ncs,inx,bl(iisd),bl(jjsd), bl(ijchd),
     *             nbl2,nsupb,npard,ijbld, mxsize,nasize,lcore,nsym,
     *             dbl_mb(inuc),dbl_mb(ibas),na,iret_limit)
c
c  output : arrays : nsupb, npar ijbl mxsize )
c------------------------------------------------
c returm in the case when LIMPAIR was changed in the BLKSIZER:
c
      if(iret_limit.gt.0) RETURN
c------------------------------------------------
c Memory checking in former texas bl() , now repleaced by PNL dbl_mb:
c
      call getmem(0,last1)
      call retmem(1)
      last2=nsym*ncs*(ncs+1)/2  +  nbl2*(nbl2+1)/2
      memneed=last1+last2
c--------
      if(memneed .gt. lcore )then
         write(6,*) ' Memory problem in txs-bl-replacing scratch : '
         write(6,*) ' requested=',memneed,' ; available=',lcore  
         write(6,*) ' check in TXSARRAG : ntxs_scr_size too small '
         write(6,*) ' nbl2=',nbl2,' last1=',last1,' last2=',last2
         if(icheck.eq.0)  call errquit(
     $        'texas:Memory problem in Blockin2', 0)
*stop 'Memory problem in Blockin2'
      endif
c--------
c------------------------------------------------
      call getmem(nbl2*(nbl2+1)/2, ncost)
      call price(inx,bl(iisd),bl(jjsd),bl(ijbld),bl(npard),
     *           bl(nsupb),nbl2,bl(ncost),nprint)
c------------------------------------------------
c save sizes:
c 
      ncspairs=ncs*(ncs+1)/2
      iisd_S=ncspairs
      jjsd_S=ncspairs
      ijshp_S=nsym*ncspairs
      ncost_S=nbl2*(nbl2+1)/2
c------------------------------------------------
c save addresses :
c
      call data_save
c------------------------------------------------
c This is the end of blocking procedure for PAIRS 
c
c The following is known at this point :
c
c  pairs ijcs are given by ics=iis(ijcs) and jcs=jjs(ijcs)
c  number of pair-blocks = NBL2
c  number of pairs in the block ibl is NPAR(ibl)
c  which pairs belong to this block : ijcs=ijbl(ibl,1-npar)
c---
c The number of contracted shell Quartets would be
c NSUPB=NBL2*(NBL2+1)/2 if they are not too big. I call
c them Super-blocks of c.s.quartets.The maximum size for 
c each of them is known and kept in array MAXSIZE(i=1,NSUPB) 
c (bl(mxsize)) (no of quartets). According to this max.size
c super-blocks will be split into a number of smaller blocks 
c which is also known (NBLOKS) but is NOT transfered.
c------------------------------------------------
      end
c===============================================================
      subroutine blkpair(bl,ncs,inx,iis,jjs,ijcheck, nbl2,nsupb,npard,
     *                   ijbld,mxsize,nasize,lcore,nsym,
     *                   datnuc,datbas,natoms,iret_limit)
      implicit real*8 (a-h,o-z)
c May96 :to be used in txs_scratch_siz2 routine to estimate memory :
      common /max_pairs/ maxpar
c save sizes :
      common /memor1_S/ iisd_S,jjsd_S,ijbld_S
      common /memor1a_S/ npard_S,ncost_S,nsupb_S, mxsize_S,nasize_S
      common /memors_S/      ijshp_S
c
      common /route/ iroute
      common /cpu/ intsize,iacc,icache,memreal
      dimension bl(*)
      dimension inx(12,*)
      dimension iis(*),jjs(*),ijcheck(*)
      dimension datnuc(5,*),datbas(13,*)
c------------------------------------------------
c Constructe pairs of contracted shells ( I>=J )
c
      call pairs(ncs,inx,iis,jjs)
c
c------------------------------------------------
c Calculate number of blocks of pairs (NBL2)
c and setup the IJBL(nbl2,*) and NPAR(*) arrays
c
c Preliminary memory reservation for 
c NPAR , IJBL , MXSIZE
c
      ndim=10000
      call memo1a(ndim,ndim,ndim,npard,ijbld,mxsize)
  100 continue
c
      call blockpnr(nbl2,maxpar,ncs,bl(ijbld),bl(npard),
     *              ijcheck,iis,jjs,inx,  ndim,ireturn,
     *              datnuc,datbas,natoms,iroute)
      if(ireturn.eq.1) then 
          ndim=2*ndim
          call retmem(3)
          call memo1a(ndim,ndim,ndim,npard,ijbld,mxsize)
          go to 100
      endif
c
c Final memory reservation for NPAR and IJBL and for NSUPB
c         for further uses
c
      npardim=nbl2
      ijbldim=nbl2*maxpar
      msizdim=nbl2*(nbl2+1)/2
      call retmem(3)
      call memo1a(npardim,ijbldim,msizdim,npard,ijbld,mxsize)
c
c------------------------------------------------
c save sizes:
      npard_S= npardim
      ijbld_S=ijbldim
      mxsize_S=msizdim
c------------------------------------------------
c Select the blocks of pairs of this same type c
c Type means this same TYPE and CONTRACTION LENGHT 
c        and general contraction deep 
c        
c Setup the IJBL(nbl2,*) and NPAR(*) arrays
c
c Constructe these blocks of contracted pairs :
c
      call blockpar(nbl2,ncs,bl(ijbld),bl(npard),ijcheck,iis,jjs,inx,
     *              datnuc,datbas,natoms,iroute)
c
c------------------------------------------------
c  Set up the vectors :
c  MXSIZE(super-block)-  maximum size of small blocks in a super-block
c  NSUPB (super-block)-  number of small blocks in this sup-block 
c and for PNL :
c  NASIZE(super-block)-  actual size of a super block (no. of quartets
c                                                     without spliting)
c
      nosb=nbl2*(nbl2+1)/2
      call getmem(nosb,nsupb)
      call getmem(nosb,nasize)
c
      call blksizer(bl(nsupb),nbl2,bl(ijbld),bl(npard), iis,jjs,inx,
     *              bl(mxsize),bl(nasize),ncs,iroute,iret_limit)
c
c iret_limit shows if everything is OK with LIMXMEM and LIMPAIR:
c
c iret_limit=0 OK
c iret_limitM .NE.0 :
c                LIMPAIR was to high for the current value of LIMXMEM :
c                LIMPAIR was decreased to iret_limpit in blksizer.
c                Whole blocking procedure has to be repeated .
c                begining from memory allocation for scratch1 
c                in texas_face.F. This info is sent to texas_face 
c                through common block /max_pairs/
c
      if(iret_limit.gt.0) RETURN
c------------------------------------------------
c save sizes:
      nsupb_S=nosb
      nasize_S=nosb
c------------------------------------------------
      end
c===============================================================
      subroutine pairs(ncs,inx,iis,jjs)
      dimension inx(12,*)
      dimension iis(*),jjs(*)
c
      ijcs=0
      do 400 ics=1,ncs
      do 400 jcs=1,ics
         ijcs=ijcs+1
         iis(ijcs)=ics
         jjs(ijcs)=jcs
  400 continue
c-------------------------------------------------------------
      end
c===============================================================
      subroutine blockpnr(nbl2,maxpar,ncs,ijbl1,npar,
     *                    ijcheck,iis,jjs,inx,
     *                    ndim,ireturn,
     *                    datnuc,datbas,natoms,iroute)
c-----------------------------------------------------------
c*  calculates a Number of blocks of contracted pairs (nbl2)
c*  and a number of pairs in a block ( npar(nbl2) )
c*  and selects a block with a maximum number of pairs
c-----
c upper limits for : maximum memory in bl() for ONE block  
c                    number of quart in a quartet-blocks
c                    number of pairs in a pair-block
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical cond1,cond2,cond3,cond4,cond5
      common /intgop/ ncache,maxprice,iprint,iblock
      common /intlim/ limxmem,limblks,limpair
      dimension datnuc(5,*),datbas(13,*)
      dimension npar(*),ijbl1(*)
      dimension ijcheck(*)
      dimension inx(12,*),iis(*),jjs(*)
c
      maxpair=limpair
c
      ireturn=0
c
      ncsp=ncs*(ncs+1)/2
      do 100 ijcs=1,ncsp
      ijcheck(ijcs)=0
  100 continue
c
c    select the blocks of pairs of this same type c
c    / type means this same TYPE and CONTRACTION LENGHT /
c    /                 and general contraction         /
c
      ibl=1
      ijbl1(1)=1
      npar(1)=1
      ijcheck(1)=1
      ijcs=1
cc
 1000 continue
      if(ijcheck(ijcs).eq.0) then
cc
         ics0=iis(ijcs)
         jcs0=jjs(ijcs)
c atoms and charges:
         iat0=inx(2,ics0)
         jat0=inx(2,jcs0)
c charge
         if(iat0.eq.0) then
            nzi0=0
         else
            nzi0=datnuc(1,iat0)
         endif
         if(jat0.eq.0) then
            nzj0=0
         else
            nzj0=datnuc(1,jat0)
         endif
c types:
         it0=inx(12,ics0)
         jt0=inx(12,jcs0)
c contraction:
         ic0=inx(5,ics0)-inx(1,ics0)
         jc0=inx(5,jcs0)-inx(1,jcs0)
c first primitive exponent:
         icb0=inx(1,ics0)+1
         jcb0=inx(1,jcs0)+1
         exi0=datbas(1,icb0)
         exj0=datbas(1,jcb0)
c gen.contr.
         ngci0=inx(4,ics0)
         ngcj0=inx(4,jcs0)
           do 1100 ib=1,ibl
           ijcs1=ijbl1(ib)
           ijp=npar(ib)
ccc
           ics1=iis(ijcs1)
           jcs1=jjs(ijcs1)
           iat1=inx(2,ics1)
           jat1=inx(2,jcs1)
           if(iat1.eq.0) then
              nzi1=0
           else
              nzi1=datnuc(1,iat1)
           endif
           if(jat1.eq.0) then
              nzj1=0
           else
              nzj1=datnuc(1,jat1)
           endif
           it1=inx(12,ics1)
           jt1=inx(12,jcs1)
           ic1=inx(5,ics1)-inx(1,ics1)
           jc1=inx(5,jcs1)-inx(1,jcs1)
           icb1=inx(1,ics1)+1
           jcb1=inx(1,jcs1)+1
           exi1=datbas(1,icb1)
           exj1=datbas(1,jcb1)
           ngci1=inx(4,ics1)
           ngcj1=inx(4,jcs1)
ccc
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           if(it0  .eq.  it1 .and. jt0  .eq.  jt1) cond1=.true.
           if(ic0  .eq.  ic1 .and. jc0  .eq.  jc1) cond2=.true.
           if(ngci0.eq.ngci1 .and. ngcj0.eq.ngcj1) cond3=.true.
           if(nzi0 .eq. nzi1 .and.  nzj0.eq. nzj1) cond4=.true.
           if(exi0 .eq. exi1 .and.  exj0.eq. exj1) cond5=.true.
ccc
      IF( iroute.eq.1 ) THEN
        if(cond1) then
          if(cond2) then
            if(cond3) then
                    ijp=ijp+1
c                   limit a number of pairs in a block :
                    if(ijp.gt.maxpair) then
                      go to 1100
                    else
                      npar(ib)=ijp
                      ijcheck(ijcs)=1
                      if(ijcs.eq.ncsp) go to 2000
                      ijcs=ijcs+1
                      go to 1000
                    endif
            endif
          endif
        endif
      ELSE
        if(cond1) then
          if(cond2) then
            if(cond3) then
              if(cond4) then
                if(cond5) then
                    ijp=ijp+1
c                   limit a number of pairs in a block :
                    if(ijp.gt.maxpair) then
                      go to 1100
                    else
                      npar(ib)=ijp
                      ijcheck(ijcs)=1
                      if(ijcs.eq.ncsp) go to 2000
                      ijcs=ijcs+1
                      go to 1000
                    endif
                endif
              endif
            endif
          endif
        endif
      ENDIF
 1100      continue
c
         ibl=ibl+1
         if(ibl.ge.ndim) then
            ireturn=1
            return
         endif
         ijbl1(ibl)=ijcs
         npar(ibl)=1
         ijcheck(ijcs)=1
         if(ijcs.eq.ncsp) go to 2000
         ijcs=ijcs+1
         go to 1000
      else
         if(ijcs.eq.ncsp) go to 2000
         ijcs=ijcs+1
         go to 1000
      endif
 2000 nbl2=ibl
c
c--------------------------------
      maxpar=0
      iblmax=1
      do 2500 ibl=1,nbl2
      ijpar=npar(ibl)
      if(ijpar.gt.maxpar) then
           maxpar=ijpar
           iblmax=ibl
      endif
 2500 continue
c---------------------------------------------------------------
      if(iprint.gt.0) write(8,495) nbl2,ncsp,maxpar,iblmax
  495 format(/' Number of Blocks of Contracted Pairs =',i7/
     *        ' total number of contracted pairs     =',i7/ 
     *        ' maximum number of pairs =',i6,' in block=',i5/)
c-------------------
      nbloks=nbl2*(nbl2+1)/2
      ntotq=ncsp*(ncsp+1)/2
      if(iprint.gt.0) write(8,496) nbloks,ntotq
  496 format(/' Number of Blocks of Contracted Quartets =',i15/
     *        ' total number of contracted quartets     =',i15/)
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blockpar(nbl2,ncs,ijbl,npar,ijcheck,iis,jjs,inx,
     *                   datnuc,datbas,natoms,iroute)
c------------------------------------------------------------
c*    select the blocks of pairs of this same type c
c*      type means :
c*    1. the same TYPE 
c*    2. the same CONTRACTION LENGHT 
c*    3. the same GENERAL CONTRACTION DEEP
c*    4. the same TYPE OF ATOMS
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical cond1,cond2,cond3,cond4,cond5
      common /intlim/ limxmem,limblks,limpair
      dimension datnuc(5,*),datbas(13,*)
      dimension npar(*),ijbl(nbl2,*)
      dimension ijcheck(*)
      dimension inx(12,*),iis(*),jjs(*)
c
      maxpair=limpair
c
      ncsp=ncs*(ncs+1)/2
      do 10 ijcs=1,ncsp
      ijcheck(ijcs)=0
   10 continue
c
      ibl=1
      ijbl(1,1)=1
      npar(1)=1
      ijcheck(1)=1
      ijcs=1
cc
 1000 continue
      if(ijcheck(ijcs).eq.0) then
cc
         ics0=iis(ijcs)
         jcs0=jjs(ijcs)
c atoms and charges:
         iat0=inx(2,ics0)
         jat0=inx(2,jcs0)
         if(iat0.eq.0) then
            nzi0=0
         else
            nzi0=datnuc(1,iat0)
         endif
         if(jat0.eq.0) then
            nzj0=0
         else
            nzj0=datnuc(1,jat0)
         endif
c types and contractions :
         it0=inx(12,ics0)
         jt0=inx(12,jcs0)
         ic0=inx(5,ics0)-inx(1,ics0)
         jc0=inx(5,jcs0)-inx(1,jcs0)
c first primitive exponent:
         icb0=inx(1,ics0)+1
         jcb0=inx(1,jcs0)+1
         exi0=datbas(1,icb0)
         exj0=datbas(1,jcb0)
c gen.contr.
         ngci0=inx(4,ics0)
         ngcj0=inx(4,jcs0)
cc
           do 1100 ib=1,ibl
           ijcs1=ijbl(ib,1)
           ijp=npar(ib)
cc
           ics1=iis(ijcs1)
           jcs1=jjs(ijcs1)
           iat1=inx(2,ics1)
           jat1=inx(2,jcs1)
           if(iat1.eq.0) then
              nzi1=0
           else
              nzi1=datnuc(1,iat1)
           endif
           if(jat1.eq.0) then
              nzj1=0
           else
              nzj1=datnuc(1,jat1)
           endif
           it1=inx(12,ics1)
           jt1=inx(12,jcs1)
           ic1=inx(5,ics1)-inx(1,ics1)
           jc1=inx(5,jcs1)-inx(1,jcs1)
           icb1=inx(1,ics1)+1
           jcb1=inx(1,jcs1)+1
           exi1=datbas(1,icb1)
           exj1=datbas(1,jcb1)
c gen.contr.
           ngci1=inx(4,ics1)
           ngcj1=inx(4,jcs1)
cc
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           if(it0  .eq.  it1 .and. jt0  .eq.  jt1) cond1=.true.
           if(ic0  .eq.  ic1 .and. jc0  .eq.  jc1) cond2=.true.
           if(ngci0.eq.ngci1 .and. ngcj0.eq.ngcj1) cond3=.true.
           if(nzi0 .eq. nzi1 .and.  nzj0.eq. nzj1) cond4=.true.
           if(exi0 .eq. exi1 .and.  exj0.eq. exj1) cond5=.true.
cc
      IF( iroute.eq.1 ) THEN
        if(cond1) then
          if(cond2) then
            if(cond3) then
                    ijp=ijp+1
c                   limit a number of pairs in a block :
                    if(ijp.gt.maxpair) then
                       go to 1100
                    else
                       npar(ib)=ijp
                       ijbl(ib,ijp)=ijcs
                       ijcheck(ijcs)=1
                       if(ijcs.eq.ncsp) go to 2000
                       ijcs=ijcs+1
                       go to 1000
                    endif
            endif
          endif
        endif
      ELSE
        if(cond1) then
          if(cond2) then
            if(cond3) then
              if(cond4) then
                if(cond5) then
                    ijp=ijp+1
c                   limit a number of pairs in a block :
                    if(ijp.gt.maxpair) then
                       go to 1100
                    else
                       npar(ib)=ijp
                       ijbl(ib,ijp)=ijcs
                       ijcheck(ijcs)=1
                       if(ijcs.eq.ncsp) go to 2000
                       ijcs=ijcs+1
                       go to 1000
                    endif
                endif
              endif
            endif
          endif
        endif
      ENDIF
 1100      continue
c
         ibl=ibl+1
         ijbl(ibl,1)=ijcs
         npar(ibl)=1
         ijcheck(ijcs)=1
         if(ijcs.eq.ncsp) go to 2000
         ijcs=ijcs+1
         go to 1000
      else
         if(ijcs.eq.ncsp) go to 2000
         ijcs=ijcs+1
         go to 1000
      endif
 2000 continue
c-------------
c checking :
      if(ibl.ne.nbl2) then
         write(6,*)'number of pair-blocks :'
         write(6,*)'  from blockpNr =',nbl2
         write(6,*)'  from blockpAr =',ibl 
         call errquit('texas: pair conflict', 0)
*         stop 'pair conflict'
      endif
c---------------------------------------------------------------
cMay15 
c     do 2500 ibl=1,nbl2
c     ijcs1=ijbl(ibl,1)
c     ics1=iis(ijcs1)
c     jcs1=jjs(ijcs1)
c     it1=inx(12,ics1)
c     jt1=inx(12,jcs1)
c     ic1=inx(5,ics1)-inx(1,ics1)
c     jc1=inx(5,jcs1)-inx(1,jcs1)
c     ijpar=npar(ibl)
c     write(8,502) ibl,ijpar,it1,jt1,ic1,jc1
c2500 continue
c 502 format(/'Block=',i3,2x,'contains ',i3,' pairs of csh.',
c    *' Itype=',i2,' Jtype=',i2,' Icon=',i2,' Jcon=',i2/)
cMay15 
c---------------------------------------------------------------
      end
c===============================================================
      subroutine blksizer(nsupb ,nbl2,ijbl,npar,iis,jjs,inx,
     *                    mxsize,nasize,ncs,iroute,iret_limit)
c-----------------------------------------------------------
c This subroutine determines the block-size for two-electron
c integrals calculations.
c The memory requimant is also estimated for ordinary two-el.
c integrals as well as for GIAO derivatives.
c-----------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical firstd
c
c absolute limits 
      common /intlim/ limxmem,limblks,limpair
c
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /cpu/ intsize,iacc,icache,memreal
      common /memmax/ ispblx, maxme1,iforwhat
c
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
c
      dimension nasize (*)
      dimension mxsize(*), nsupb(*)
      dimension npar(*),ijbl(nbl2,*)
      dimension inx(12,*),iis(*),jjs(*)
c
c-----------------------------------------------------
      iret_limit=0
c-----------------------------------------------------
c  iforwhat shows what the Blockin2 routins are called for :
c  1) for ordinary two-el.integrals (ifor=1)
c  2) for GIAO two-el. derivatives  (ifor=2)
c  3) for gradient derivatives      (ifor=3)
c  4) for second derivatives        (ifor=4)
c
c  This is used to send info about memory to the Calcint2
c-----------------------------------------------------
c  Calculates the total number of blocks of contracted
c  shells quartets NBLOKS (only for info).
c  Sets up the MXSIZE(*) and NSUPB(*) arrays for each
c  Super-block (number of them is NBL2*(NBL2+1)/2 ) .
c
      mem2max=0
      mem2ijx=0
      mem2klx=0
      mem4max=0
c
      nblsmax=0
      maxme1=0
      nbloks=0     ! counter of all blocks
      nbl12=0      ! counter of super-blocks
      do 100 ibl=1,nbl2
      ijpar=npar(ibl)
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
c
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
      itype1=itype
      jtype1=jtype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      nfij=lenn(itype1)*lenn(jtype1)
c
      lci=inx(5,ics1)-inx(1,ics1)
      lcj=inx(5,jcs1)-inx(1,jcs1)
c
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
c
      do 100 kbl=1,ibl     
      nbl12=nbl12+1
c
      klpar=npar(kbl)
      klcs1=ijbl(kbl,1)
      kcs1=iis(klcs1)
      lcs1=jjs(klcs1)
      ktype=inx(12,kcs1)
      ltype=inx(12,lcs1)
c
      ktype1=ktype
      ltype1=ltype
      if(ktype.gt.4) ktype1=ktype-1
      if(ltype.gt.4) ltype1=ltype-1
      if(ktype1.gt.5) ktype1=ktype1-1
      if(ltype1.gt.5) ltype1=ltype1-1
      nfkl=lenn(ktype1)*lenn(ltype1)
      nfijkl=nfij*nfkl
c
      lck=inx(5,kcs1)-inx(1,kcs1)
      lcl=inx(5,lcs1)-inx(1,lcs1)
c
      ngck=inx(4,kcs1)
      ngcl=inx(4,lcs1)
c
c calculate maximum memory needed for one block :
c 1) for ordinary two-el.integrals (ifor=1)
c 2) for GIAO two-el. derivatives  (ifor=2)
c 3) for gradient derivatives      (ifor=3)
c 4) for second derivatives        (ifor=4)
c
ctest
c     write(8,*)'    '
c     write(8,*)'Super-block=',nbl12,' Ibl=',ibl,' Kbl=',kbl
c
c----
      ifor=1
c----
      IF( iroute.eq.1 ) THEN
      call blksize1(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *              lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *              memor2,memor2ij,memor2kl,memor4,ifor)
      ELSE
      call blksize2(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *              lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *              memor2,memor2ij,memor2kl,memor4,ifor)
      ENDIF
c     write(8,*)' memor2,memor2ij,memor2kl,memor4'
c     write(8,*) memor2,memor2ij,memor2kl,memor4
c----
      if(memor2.gt.mem2max) then
         mem2max=memor2
         ispbl2=nbl12
      endif
      if(memor2ij.gt.mem2ijx) then
         mem2ijx=memor2ij
         ispbl2ij=nbl12
      endif
      if(memor2kl.gt.mem2klx) then
         mem2klx=memor2kl
         ispbl2kl=nbl12
      endif
      if(memor4.gt.mem4max) then
         mem4max=memor4
         ispbl4=nbl12
      endif
c-------------------------------------------------
      maxsize=       icache /nfijkl
      if(maxsize.eq.0) maxsize=1
      if(maxsize.gt.limblks) maxsize=limblks
c-------------------------------------------------
 1234 continue
c
      ikbl=0       ! counter of blocks belonging to one super-block
      maxmem=0
      maxqrt=0
      mxsize(nbl12)=maxsize
c
      nquart=ijpar*klpar
      if(ibl.eq.kbl) nquart=ijpar*(ijpar+1)/2
cpnl:
      nasize(nbl12)=nquart       ! texas-size of a super-block
c
c  calculate the number of blocks of contracted quartets
c
          ikbl=ikbl+1
          ijsize=ijpar
          klsize=klpar
c----------------------------------------------
c put the upper limit on the maxmem :
c-->  mqrt1=memor2ij+memor2kl+memor4
c memory for ONE quartet and one pair in block2 :
c
      mqrt1=memor2ij+memor2kl+ memor4
      if(mqrt1.ge.limxmem) then
        write(6,*)'A block with only one quartet needs more memory than'
        write(6,*)'the upper limit for a whole block given in the input'
        write(6,*)'input data =',limxmem,' request for 1 quart.=',mqrt1
        write(6,*)'     Change the value of LIMXMEM '
        write(6,*)'program stoped in the Blksizer routine'
        call errquit('texas:blksizer', 0)
cccc     stop 'blksizer'
      endif
c----------------------------------------------
c put the upper limit on the maxmem :
c     mqrt1=ijpar*memor2ij+klpar*memor2kl+ memor4
c memory for ONE quartet but all pairs :
c
      mqrt1=ijpar*memor2ij+klpar*memor2kl+ memor4
c
      if(mqrt1.ge.limxmem) then
        write(6,*)'       For the current value of LIMXMEM             '
        write(6,*)'            pair-blocks are too big.                '
c
c npar_limit=[ sqrt( 4*m2 +4*m4*limxmem) - 2*m2 ]/2*m4
c 
        m2=min(memor2ij,memor2kl)
        m4=memor4
        npar_lim=sqrt( dble(m2*m2 + m4*limxmem) ) -m2
        npar_lim=npar_lim/m4
        if(npar_lim.eq.0) then
           write(6,*)'program stoped in the Blksizer routine'
           write(6,*)'LIMPAIR can not be adjusted to LIMXMEM'
           write(6,*)'     Change the value of LIMXMEM '
           call errquit('texas:blksizer', 0)
cccc       stop 'blksizer'
        endif
        write(6,*)' LIMPAIR is decreased from ',limpair,' to',npar_lim
        limpair=npar_lim
        iret_limit=1
        RETURN
      endif
c----------------------------------------------
c for PNL calculate what SHOULD be maximum size of a block due to the
c memory limit :
c
      memorx=ijpar*memor2ij+klpar*memor2kl+ mxsize(nbl12)*memor4
      if(memorx.gt.limxmem) then
         mxqrts=(limxmem-ijpar*memor2ij-klpar*memor2kl)/memor4
         if(mxqrts.eq.0) mxqrts=1
         if(mxqrts.gt.limblks) mxqrts=limblks
c        write(6,*)' from sizer : sup-blk=',nbl12,' txs_size=',
c    *   nasize(nbl12),' memory=',memorx,' limxmem=',limxmem
c        write(6,*)' max.size changed from ',mxsize(nbl12),' to',mxqrts
         mxsize(nbl12)=mxqrts
      endif
c----------------------------------------------
      if(mxsize(nbl12).gt.maxqrt) maxqrt=mxsize(nbl12)
c----------------------------------------------
      memory=ijpar*memor2ij+klpar*memor2kl+ mxsize(nbl12)*memor4
c----------------------------------------------
c MEMORY is the memory needed for pairs and quartets in a block
c it is also necessary to take into account a memory reserved
c already for prec2ij and prec2kl as well as memory in Blockin4:
c
      if(intsize.ne.1) then
        memin4=7*(ikbl/intsize+1) +  maxqrt+2*( maxqrt/intsize +1 )
      else
        memin4=7*ikbl+3*maxqrt
      endif
c----------------------------------------------
      memory=memory+memor2 + memin4
c----------------------------------------------
      if(memory.gt.maxmem) then
         ijpax=ijpar
         klpax=klpar
         nblsx=mxsize(nbl12)
         maxmem=memory
      endif
c----------------------------------------------
c select a super-block with maximum memory requirement for ONE block
c
      nsupb(nbl12)=ikbl
      nbloks=nbloks+ikbl
      if(maxmem.gt.maxme1) then
         maxme1=maxmem
         ispblx=nbl12
c        remember also : 
         iblx=ibl
         kblx=kbl
c
         itypx=itype1
         jtypx=jtype1
         ktypx=ktype1
         ltypx=ltype1
c
         lcix=lci
         lcjx=lcj
         lckx=lck
         lclx=lcl
c
         ngcix=ngci
         ngcjx=ngcj
         ngckx=ngck
         ngclx=ngcl
c
         m2  =memor2 
         min4=memin4
         m2ij=memor2ij
         m2kl=memor2kl
         m4  =memor4
         ijpay=ijpax
         klpay=klpax
         nblsy=nblsx
      endif
c--------------
c Select a super-block with the maximum 
c block-size for ONE block :
c
      if(maxqrt.gt.nblsmax) then
         nblsmax=maxqrt
         ispblq=nbl12
      endif
c--------------
  100 continue
c-----------------------------------------------------------
c Print out info about maximum sizes of super-blocks (PNL):
c
      IF(nprint.gt.0) THEN
c
      write(8,*)'--------------------------------------------------'
      write(8,*)' Maximum sizes of super-blocks upon restrictions :'
      write(8,*)'         cache memory < ',icache
      write(8,*)'         block size   < ',limblks
      write(8,*)'         memory limit < ',limxmem
      write(8,*)'--------------------------------------------------'
c     do 200 isupb=1,nbl12
c     write(8,*)' super-block=',isupb,' txs-size=',nasize(isupb),
c    *                                ' max-size=',mxsize(isupb)
c 200 continue
c     write(8,*)'--------------------------------------------------'
c
      ENDIF
c-----------------------------------------------------------
c This is memory requested for one block
c and it should be the same like maxme1 .
c        mblok1=ijpay*m2ij+klpay*m2kl+nblsy*m4
c        mblok1=mblok1+m2 + min4
c-----------------------------------------------------------
c Now calculate how much memory will be needed for the NMR
c two-el. int. derivatives with the SAME blok-size structure
c like for ordinary integrals :
c
      ifor=2
      IF( iroute.eq.1 ) THEN
      call blksize1(iblx,kblx,ijpay,klpay,itypx,jtypx,ktypx,ltypx,
     *              lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *              memor2,memor2ij,memor2kl,memor4,ifor)
      ELSE
      call blksize2(iblx,kblx,ijpay,klpay,itypx,jtypx,ktypx,ltypx,
     *              lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *              memor2,memor2ij,memor2kl,memor4,ifor)
      ENDIF
c
         maxme2=ijpay*memor2ij+klpay*memor2kl+nblsy*memor4
         maxme2=maxme2+memor2 + min4
c-----------------------------------------------------------
c Print out info concerning memory requiments :
c-----------------------------------------------------------
      IF(nprint.eq.0) RETURN
c-----------------------------------------------------------
      write(8,*)
      write(8,497) nbl12,nblsmax,ispblq,nbloks
  497 format(/' Number of Super-Blocks of Contracted Quartets =',i20/
     *        ' maximum block-size =',i7,' in S-block=',i6/
     *        ' Total number of Blocks of Contracted Quartets =',i20/)
c------------pnl -------------------------------------------
c
c print it only for the CHECK run :
c
      IF(icheck.ne.0) THEN
      write(8,498) nblsmax,ispblq, mem2max,ispbl2, mem2ijx,ispbl2ij,
     *             mem2klx,ispbl2kl,mem4max,ispbl4, maxme1,ispblx 
  498 format( ' ***** Maximum values required for Super-blocks *****'/
     *        ' maximum block-size    =',i10,'   in S-block=',i6/
     *        ' memory for pair-precal=',i10,'   in S-block=',i6/
     *        ' memory for one pair IJ=',i10,'   in S-block=',i6/
     *        ' memory for one pair KL=',i10,'   in S-block=',i6/
     *        ' memory for one quartet=',i10,'   in S-block=',i6/
     *        ' memory for ONE block  =',i10,'   in S-block=',i6/)
c
      ENDIF
c
c print it for the CHECK run and for ordinary two-el.integrals :
c
      IF(icheck.ne.0 .or. iforwhat.eq.1) THEN
c     write(6,499) ispblx, iblx,kblx,m2,
c    *  m2ij,ijpay, m2kl,klpay, m4,nblsy, maxme1,maxme1+1
c
      write(8,499) ispblx, iblx,kblx,m2,
     *  m2ij,ijpay, m2kl,klpay, m4,nblsy, maxme1,maxme1+1
  499 format( ' ** Most memory consuming block in the Super-block  **'/
     *        ' number =',i6,' made out of pair-blocks =',2i6//
     *        '   for ordinary two-electron integrals requires :'//
     *        '   ',i7,' d.w. of memory for pair-precalc.'/
     *        '   ',i7,' d.w. of memory for one pair ij (',i5,' pairs)'/
     *        '   ',i7,' d.w. of memory for one pair kl (',i5,' pairs)'/
     *        '   ',i7,' d.w. of memory for one quartet (',i5,' quart)'/
     *        '   '  ,i9,' d.w. of memory for the whole block  '/
     *        ' first free address in bl will be=',i10/)
      ENDIF
c
c NMR :
c print it for the CHECK run and for the GIAO two-el.integrals :
c
cPNL----------
c     IF(icheck.ne.0 .or. iforwhat.eq.2) THEN
c     write(8,501) ispblx, iblx,kblx,memor2,
c    *  memor2ij,ijpay, memor2kl,klpay, memor4,nblsy, maxme2
c 501 format( ' ** Most memory consuming block in the Super-block  **'/
c    *        ' number =',i6,' made out of pair-blocks =',2i6//
c    *        '   for the GIAO two-electron integrals requires :'//
c    *        '   ',i7,' d.w. of memory for pair-precalc.'/
c    *        '   ',i7,' d.w. of memory for one pair ij (',i5,' pairs)'/
c    *        '   ',i7,' d.w. of memory for one pair kl (',i5,' pairs)'/
c    *        '   ',i7,' d.w. of memory for one quartet (',i5,' quart)'/
c    *        '   '  ,i9,' d.w. of memory for the whole block  '/)
c transmit the maxme2 to the Calcint2 through the common /memmax/
c
c        if(iforwhat.eq.2) maxme1=maxme2
c     ENDIF
cPNL----------
c--------
      end
c===============================================================
      subroutine blksize1(ibl,kbl,ijpar,klpar,ityp,jtyp,ktyp,ltyp,
     *                    lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *                    memory2,memor2ij,memor2kl,memory4,ifor)
      implicit real*8 (a-h,o-z)
      character*4 where
c----------------------------------------------------------------------
c  Calculates a memory needed in a Super-block for :
c   1) prec2ij and prec2kl (called for each super-block)
c                                     - memory2
c
c  for each block in a super-block calculates :
c
c   2) memory needed for ONE pair ij  - memor2ij
c   3) memory needed for ONE pair kl  - memor2kl
c   4) memory needed for ONE quartet  - memory4
c   
c----------------------------------------------------------------------
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c----------------------------------------------------------------------
      common /logic4/ nfu(1)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
c----------------------------------------------------------------------
c ifor=1   for ordinary two-el.integrals
c ifor=2   for NMR GIAO integral derivatives
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
c----------------------------------------------------------------------
      lci=lcix
      lcj=lcjx
      lck=lckx
      lcl=lclx
      ngci=ngcix
      ngcj=ngcjx
      ngck=ngckx
      ngcl=ngclx
c----------------------------------------------------------------------
c set up commons /obarai/ and /shell/
c
      NQI=NDEGE(ITYP)
      NQJ=NDEGE(JTYP)
      NQK=NDEGE(KTYP)
      NQL=NDEGE(LTYP)
      NSIJ=NQI+NQJ-1
      NSKL=NQK+NQL-1
c
      if(where.eq.'shif') then
        nsij=nsij+1
        nskl=nskl+1
      endif
C
      MMAX=NSIJ+NSKL-1
      mmax1=mmax-1
C
      LNI=LEN(ITYP)
      LNJ=LEN(JTYP)
      LNK=LEN(KTYP)
      LNL=LEN(LTYP)
      LNIJKL=LNI*LNJ*LNK*LNL
C
      LNIJ=LENSM(NSIJ)
      LNKL=LENSM(NSKL)
C
       NSIJ1=NSIJ+1
       NSKL1=NSKL+1
       NQIJ=NQI
       IF(NQJ.GT.NQI) NQIJ=NQJ
       NQIJ1=NQIJ+1
       NQKL=NQK
       IF(NQL.GT.NQK) NQKL=NQL
       NQKL1=NQKL+1
c
       LSHELLT=0
       IF(ITYP.EQ.3) LSHELLT=LSHELLT+1
       IF(JTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(KTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(LTYP.EQ.3) LSHELLT=LSHELLT+1
c----------------------------------------------------------------------
       lcij=lci*lcj
       lckl=lck*lcl
c----------------------------------------------------------------------
c for new general contraction handling :
c
      ngcij=(ngci+1)*(ngcj+1)
      ngckl=(ngck+1)*(ngcl+1)
      ngcd =ngcij*ngckl
c----------------------------------------------------------------------
c Memory reserved for the WHOLE super-block in prec2ij,prec2kl :
c
      mem2ij=3*ijpar*lcij
      mem2kl=3*klpar*lckl
c
      if(kbl.ne.ibl) then 
         memory2=mem2ij+mem2kl
      else
         memory2=mem2ij
      endif
c
c----------------------------------------------------------------------
      nfumax=nfu(mmax)
      if(nsij.ge.nskl) then
         nfha=nfumax*klpar*lckl
      else
         nfha=nfumax*ijpar*lcij
      endif
c----------------------------------------------------------------------
c as allocated in memo5c_1 :
c
      memory2=memory2 + 3*nfha 
c----------------------------------------------------------------------
c Memory for the ONE pair ij,kl and ONE quartet in a Super-block is :
c
      ijpar1=1
      klpar1=1
      nbls1=1
c------
c pair-precalculations:
c
      call in5a(ijpar1,mmax1, memprij)
      call in5b(klpar1,mmax1, memprkl)
c------------------------
      memor2ij=memprij
      memor2kl=memprkl
c
      if(where.eq.'shif') then
c      add memory reserved in memo6 routine :
         memor2ij=memor2ij+3
         memor2kl=memor2kl+3
      endif
c------------------------
c quartet-calculations:
c
      call in5c(nbls1,mmax1,nfumax,mempre4)
c------------------------
c quartet-trobsa,assemble:
c
      call in4a(nbls1,memasse,memtrob,where)
c------------------------
c quartet-amshift:
c
      call in4b(nbls1,memamsh,where)
c------------------------
c
      memasse=memasse+3   ! convert in assemble
cMay96memamsh=memamsh+3   ! convert in amshift 
      memamsh=memamsh+6   ! convert in amshift 
      if(ngcd.gt.1) then 
         memgenc=(ngcij+ngckl+ngcd*lnij*lnkl)
         memasse=memasse+memgenc
      endif
c----------------------------------------------------------------------
c What is needed at once ? Memory for :
c
c    1. Precalc4 + Trobsa+Assemble     proportional to NBLS
c    2. Precalc4 + Assemble+Amshift    proportional to NBLS
c-----
      mem1 = mempre4+memtrob+memasse
      mem2 = mempre4+memasse+memamsh
c
c------------------------
      memory4 =max(mem1,mem2)
c------------------------
      end
c===============================================================
      subroutine blksize2(ibl,kbl,ijpar,klpar,ityp,jtyp,ktyp,ltyp,
     *                    lcix,lcjx,lckx,lclx, ngcix,ngcjx,ngckx,ngclx,
     *                    memory2,memor2ij,memor2kl,memory4,ifor)
      implicit real*8 (a-h,o-z)
      character*4 where
c----------------------------------------------------------------------
c  Calculates a memory needed in a Super-block for :
c   1) prec2ij and prec2kl (called for each super-block)
c                                     - memory2
c
c  for each block in a super-block calculates :
c
c   2) memory needed for ONE pair ij  - memor2ij
c   3) memory needed for ONE pair kl  - memor2kl
c   4) memory needed for ONE quartet  - memory4
c   
c----------------------------------------------------------------------
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c----------------------------------------------------------------------
      common /logic4/ nfu(1)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
c----------------------------------------------------------------------
c ifor=1   for ordinary two-el.integrals
c ifor=2   for NMR GIAO integral derivatives
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
c----------------------------------------------------------------------
      lci=lcix
      lcj=lcjx
      lck=lckx
      lcl=lclx
      lcij=lci*lcj
      lckl=lck*lcl
      ngci=ngcix
      ngcj=ngcjx
      ngck=ngckx
      ngcl=ngclx
c----------------------------------------------------------------------
c set up commons /obarai/ and /shell/
c for subroutines in4a,in4b 
c
      NQI=NDEGE(ITYP)
      NQJ=NDEGE(JTYP)
      NQK=NDEGE(KTYP)
      NQL=NDEGE(LTYP)
      NSIJ=NQI+NQJ-1
      NSKL=NQK+NQL-1
c
      if(where.eq.'shif') then
        nsij=nsij+1
        nskl=nskl+1
      endif
C
      MMAX=NSIJ+NSKL-1
      mmax1=mmax-1
C
      LNI=LEN(ITYP)
      LNJ=LEN(JTYP)
      LNK=LEN(KTYP)
      LNL=LEN(LTYP)
      LNIJKL=LNI*LNJ*LNK*LNL
C
      LNIJ=LENSM(NSIJ)
      LNKL=LENSM(NSKL)
C
       NSIJ1=NSIJ+1
       NSKL1=NSKL+1
       NQIJ=NQI
       IF(NQJ.GT.NQI) NQIJ=NQJ
       NQIJ1=NQIJ+1
       NQKL=NQK
       IF(NQL.GT.NQK) NQKL=NQL
       NQKL1=NQKL+1
c
       LSHELLT=0
       IF(ITYP.EQ.3) LSHELLT=LSHELLT+1
       IF(JTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(KTYP.EQ.3) LSHELLT=LSHELLT+1
       IF(LTYP.EQ.3) LSHELLT=LSHELLT+1
c----------------------------------------------------------------------
c for new general contraction handling :
c
      ngcij=(ngci+1)*(ngcj+1)
      ngckl=(ngck+1)*(ngcl+1)
      ngcd =ngcij*ngckl
c----------------------------------------------------------------------
c Memory reserved for the WHOLE super-block in prec2ij,prec2kl :
c
      mem2ij=(ijpar+2)*lcij
      mem2kl=(klpar+2)*lckl
c
      if(kbl.ne.ibl) then 
         memory2=mem2ij+mem2kl
      else
         memory2=mem2ij
      endif
c
c----------------------------------------------------------------------
c Memory for the ONE pair ij,kl and ONE quartet in a Super-block is :
c----------------------------------------------------------------------
c pair-precalculations (as reserved in memo5a,b):
c
c proportional to the number of pairs :
      memprij=3+14*lcij
      memprkl=3+14*lckl
c
c constat part (not proportional to pairs):
c
      memcoij=lci+lcj+lcij + mmax1*lcij
      memcokl=lck+lcl+lckl + mmax1*lckl
      if(ngcd.gt.1) then
        memcoij=memcoij+ngcij*lcij
        memcokl=memcokl+ngckl*lckl
      endif
      memcons=memcoij+memcokl
c------------------------
      memor2ij=memprij
      memor2kl=memprkl
c
      if(where.eq.'shif') then
c      add memory reserved in memo6 routine :
         memor2ij=memor2ij+3
         memor2kl=memor2kl+3
      endif
c------------------------
c quartet-calculations (reserved in memo5c :
c
c proportional to block-szie (nbls):
c
      mempre4=20
      if(ngcd.gt.1) mempre4=20+ngcd+1
c
c constant (not prop. to nbls):
c
      nfumax=nfu(mmax)
      nfha=nfumax*max(lcij,lckl)
      memcons=memcons+(4*(ngcd+1)+3*nfha)
c------------------------
c add constant part of requested memory to memory2
c
      memory2=memory2 + memcons
c------------------------
c quartet-trobsa,assemble:
c
      nbls1=1
c
      call in4a(nbls1,memasse,memtrob,where)
c------------------------
c quartet-amshift:
c
      call in4b(nbls1,memamsh,where)
c------------------------
c
cMay96memamsh=memamsh+3   ! convert in amshift 
      memamsh=memamsh+6   ! convert in amshift 
      if(ngcd.gt.1) then 
         memgenc=(ngcij+ngckl+ngcd*lnij*lnkl)
         memasse=memasse+memgenc
      endif
c----------------------------------------------------------------------
c What is needed at once ? Memory for :
c
c    1. Precalc4 + Trobsa+Assemble     proportional to NBLS
c    2. Precalc4 + Assemble+Amshift    proportional to NBLS
c-----
      mem1 = mempre4+memtrob+memasse
      mem2 = mempre4+memasse+memamsh
c
c------------------------
      memory4 =max(mem1,mem2)
c------------------------
      end
c===============================================================
c
c This set of routines named in... is called from
c the blksize1 subroutine in oreder to estimate
c the memory request for a given Super-block  
c and use it latter to set up the block-size .
c==============================================================
      subroutine in4a(nbls,memasse,memtrob,where)
      character*4 where
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /logic1/ ndege(1)
      common /logic2/ len(1)
      common /logic3/ lensm(1)
      common /logic4/ nfu(1)
c------------------------------------------
c Memory requested in the Memo4a subroutine
c------------------------------------------
      memtrob=0
      memasse=0
c------------------------------------------
c for trobsa and  for assemble :
c
      mem0=lnij*lnkl
c ------------------------------------------
c       Memory for "assemble"
c ------------------------------------------
      ngcijkl=(ngci+1)*(ngcj+1)*(ngck+1)*(ngcl+1)
      nblsg=nbls*ngcijkl
c
      if(where.ne.'shif') then
c       ----------------------------
        memasse=nblsg*(lnijkl + mem0)
c       ----------------------------
      else
c       ----------------------------
        memasse=nblsg*(7*lnijkl + mem0+6*nfu(nsij)*nfu(nskl)  )
c       ----------------------------
      endif
c
      if(mmax.le.2) then
         memasse=memasse+2*nbls
         return
      endif
c
        IF(LSHELLT.GT.0) THEN
          if(where.eq.'shif' .or. where.eq.'forc') then
           mbfkl12=lnij*nfu(nqkl1+1)*nbls + 6*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12=nfu(nqij1+1)*lnkl*nbls + 6*nfu(nqij+1)*nfu(nskl)*nbls
          else
           mbfkl12=lnij*nfu(nqkl+1)*nbls 
           mbfij12=nfu(nqij+1)*lnkl*nbls
          endif
c
          if(lshellt.gt.1) then
c           ----------------------
            memasse=memasse+2*(mbfij12+mbfkl12)
c           ----------------------
          else
c           ----------------------
            memasse=memasse+(mbfij12+mbfkl12)
c           ----------------------
          endif
c     
        IF( LSHELLT.GT.1 ) THEN
c
          if(where.eq.'shif' .or. where.eq.'forc') then
            mbf2l=nfu(nqij1+1)*nfu(nqkl1+1)*nbls 
     *         +6*nfu(nqij +1)*nfu(nqkl +1)*nbls
            mbfkl3=lnij*4*nbls + 6*nfu(nsij)*nbls
            mbfij3=4*lnkl*nbls + 6*nfu(nskl)*nbls
          else
            mbf2l=nfu(nqij+1)*nfu(nqkl+1)*nbls 
            mbfkl3=lnij*nbls
            mbfij3=lnkl*nbls
          endif
c
          if(lshellt.gt.2) then
c           ----------------------
            memasse=memasse+4*mbf2l
c           ----------------------
          else
c           ----------------------
            memasse=memasse+2*mbf2l
c           ----------------------
          endif
c           ----------------------
            memasse=memasse+(mbfij3+mbfkl3)
c           ----------------------
c
        IF( LSHELLT.GT.2 ) THEN
c
          if(where.eq.'shif' .or. where.eq.'forc') then
            mbf3l1=max( nfu(nqij1+1),nfu(nqkl1+1) )
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=4*mbf3l1*nbls + 6*mbf3l0*nbls
          else
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=mbf3l0*nbls
          endif
c
          if(lshellt.gt.3) then
c           ----------------------
            memasse=memasse+4*mbf3l
c           ----------------------
           else
c           ----------------------
            memasse=memasse+2*mbf3l
c           ----------------------
           endif
c
        IF( LSHELLT.GT.3 ) then
c
          if(where.eq.'shif' .or. where.eq.'forc') then
            i4s =16*nbls + 6*nbls
          else
            i4s =nbls
          endif
c           ----------------------
            memasse=memasse+i4s
c           ----------------------
        ENDIF
        ENDIF
        ENDIF
        ENDIF
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Memory handling for Obara-Saika-Tracy method
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
c     l11=mmax
c     l12=lensm(mmax)     
c     mem1=l11*l12
      mem1=mmax*lensm(mmax)
cc
      mem2=0
      if(nsij.ge.nskl) then
        klstep=0
        do 10 ijstep=mmax,nsij,-1
        klstep=klstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2=mem2+ijkld
   10   continue
      else
        ijstep=0
        do 11 klstep=mmax,nskl,-1
        ijstep=ijstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2=mem2+ijkld
   11   continue
      endif
c           ----------------------
            memtrob=nbls*(mem0+mem1+mem2)
c           ----------------------
c
      end
c
c==============================================================
      subroutine in4b(nbls,memamsh,where)
      character*4 where
c--
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
c------------------------------------------
c Memory requested in the Memo4b subroutine
c------------------------------------------
c       Memory for amshift 
c
            mwvus=max(lnij,lnkl)*max(nfu(nqj+1),nfu(nql+1))
            mxij=nfu(nqi+1)*nfu(nqij+1)*lnkl
c
            mwij=mwvus
            mwij=mwij*nbls
            mxij=mxij*nbls
        if(where.eq.'shif') then
            mwij=6*mwij
            mxij=6*mxij
        endif
c           ----------------------
            memamsh=mwij+mxij
c           ----------------------
        IF(LSHELLT.GT.0) THEN
c
            mvus=mwvus
            myz=nfu(nqi+1)*nfu(nqj+1)*nfu(nqkl+1)
            mvus=mvus*nbls
            myz=myz*nbls
        if(where.eq.'shif') then
            mvus=6*mvus
            myz =6*myz 
        endif
c           ----------------------
            memamsh=memamsh+(mvus+myz)
c           ----------------------
c
        IF( LSHELLT.GT.1 ) THEN
            mbf2l=nfu(nqij+1)*nfu(nqkl+1) *nbls
            if(where.eq.'shif') then
               mbf2l=6*mbf2l
            endif
c           ----------------------
            memamsh=memamsh+(2*mvus+myz)
c           ----------------------
c
          if(lshellt.gt.2) then
c           ----------------------
            memamsh=memamsh+4*mbf2l 
c           ----------------------
          else
c           ----------------------
            memamsh=memamsh+mbf2l
c           ----------------------
          endif
c
        IF( LSHELLT.GT.2 ) THEN
c
         mnbls=nbls
         if(where.eq.'shif') mnbls=6*nbls
c
         if(lshellt.gt.3) then
c           ----------------------
            memamsh=memamsh+4*mnbls
c           ----------------------
          else
c           ----------------------
            memamsh=memamsh+2*mnbls
c           ----------------------
          endif
c
        ENDIF
        ENDIF
        ENDIF
c
      end
c==============================================================
      subroutine in5a(npij,mmax1, memory)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory requested in the Memo5a subroutine
c------------------------------------------
c Memory handling for left-hand pairs:
c Total number of calls of Getmem is 12 or 14 (if gen.con.)
c------------------------------------------
      ijpar=npij
c------------------------------------------
       ndi=   ijpar*lci
       ndj=   ijpar*lcj
c     ---------------------------------
      memory=2*(ndi+ndj) + 3*ijpar
c     ---------------------------------
       ndij =ndi*lcj
       ndij3=ndij*3
c     ---------------------             
      memory=memory+3*ndij3
c     ---------------------                  
CMay96 memory=memory+2*ndij
      memory=memory+3*ndij
c     ---------------------             
      memory=memory+ndij3
c     ---------------------             
      ndijm=ndij*mmax1
c     ---------------------             
      memory=memory+ndijm
c     ---------------------                 
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
      if(ngcd.gt.1) then
        ndig=ndi*ngci1
        ndjg=ndj*ngcj1
c       -----------------------
        memory=memory+ndig+ndjg
c       -----------------------
      endif
c
      end
c==============================================================
      subroutine in5b(npkl,mmax1, memory)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory requested in the Memo5b subroutine
c------------------------------------------
c Total number of calls of Getmem is 12 or 14 (if gen.con.)
c------------------------------------------
      klpar=npkl
c------------------------------------------
       ndk=   klpar*lck
       ndl=   klpar*lcl
c     ---------------------------------
      memory=2*(ndk+ndl) + 3*klpar
c     ---------------------------------
       ndkl=ndk*lcl
       ndkl3=ndkl*3
c     ---------------------------------
      memory=memory+3*ndkl3
c     ---------------------------------
CMay96 memory=memory+2*ndkl
      memory=memory+3*ndkl
c     ---------------------------------
      memory=memory+ndkl3
c     ---------------------------------
      ndklm=ndkl*mmax1
c     ---------------------------------
      memory=memory+ndklm
c     ---------------------------------
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
      if(ngcd.gt.1) then
        ndkg=ndk*ngck1
        ndlg=ndl*ngcl1
c       -------------------------------
        memory=memory+ndkg+ndlg
c       -------------------------------
      endif
c
      end
c==============================================================
      subroutine in5c(nbls,mmax1,nfumax,memor4)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c memory requested in the Memo5c subroutine
c------------------------------------------
c memory handling 
c Total number of calls of Getmem is 24 or 26 (if gen.cont)
c reserve memory for quartets ijkl
c------------------------------------------
      nblsi=nbls
      if(intsize.ne.1) nblsi=nbls/intsize+1
c     ------------------------------------
      memor4=4*nblsi
c     ------------------------
      memor4=memor4 + 4*nbls
c     ------------------------
      nbmx=nbls*mmax1
c     ------------------------
      memor4=memor4 + 2*nbmx
c     ------------------------
      nbls3=nbls*3
c     ------------------------
      memor4=memor4 + 5*nbls3 + 3*nbls
c     ------------------------
      memor4=memor4 + nbls3*nfumax
c     ------------------------
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c     ------------------------
cccc  memor4=memor4 + 4*ngcd
c     ------------------------
      if(ngcd.gt.1) then
c       ----------------------
        memor4=memor4 + nbls*(1+ngcd)
c       ----------------------
      endif
c
      end
c===============================================================
c The new price routine :
c=====
      subroutine price(inx,iis,jjs,ijbl,npar,nsupb,nbl2,ncost,nprint)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
      common /logic4/ nfu(1)
c
      dimension nsupb(*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*), ncost(*),npar(*)
c-------------------------
c*
      if(nprint.gt.1) then
      write(8,*  )' FIRST 1000 BLOCKS '
      write(8,505)
  505 format(/
     *' Block  MMAX   type   contrac.   nfun/q      Nqrt   Nfunc',2x,
     * 'Price/i')
      endif
c
c contracted integrals counter
c
      maxpr=0
      minpr=100000
      maxst=0
      ntotq=0
c
      ikbl=0
      do 100 ibl=1,nbl2
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
c     int=inx(3,ics1)
c     jnt=inx(3,jcs1)
c     ijnt=int*jnt
      itype1=itype
      jtype1=jtype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      icont=inx(5,iis(ijcs1))-inx(1,iis(ijcs1))
      jcont=inx(5,jjs(ijcs1))-inx(1,jjs(ijcs1))
      ijcont=icont*jcont
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngci1=ngci+1
      ngcj1=ngcj+1
      nparij=npar(ibl)
      nfij  =lenn(itype1)*lenn(jtype1)
      ngcij =ngci1*ngcj1
c-
         do 100 kbl=1,ibl
         ikbl=ikbl+1
         klcs1=ijbl(kbl,1)
         kcs1=iis(klcs1)
         lcs1=jjs(klcs1)
         nparkl=npar(kbl)
c
         nquart=nparij*nparkl
         if(kbl.eq.ibl) nquart=nparij*(nparij+1)/2
c
         ntotq=ntotq+nquart
c
c        knt=inx(3,kcs1)
c        lnt=inx(3,lcs1)
c        klnt=knt*lnt
         ktype=inx(12,kcs1)
         ltype=inx(12,lcs1)
         ktype1=ktype
         ltype1=ltype
         if(ktype.gt.4) ktype1=ktype-1
         if(ltype.gt.4) ltype1=ltype-1
         if(ktype1.gt.5) ktype1=ktype1-1
         if(ltype1.gt.5) ltype1=ltype1-1
c
         nfijkl=nfij*lenn(ktype1)*lenn(ltype1)
         nfun=nquart*nfijkl
c
         kcont=inx(5,iis(klcs1))-inx(1,iis(klcs1))
         lcont=inx(5,jjs(klcs1))-inx(1,jjs(klcs1))
         ijklc=ijcont*kcont*lcont
c
c number of general contractions
c
         ngck=inx(4,kcs1)
         ngcl=inx(4,lcs1)
         ngck1=ngck+1
         ngcl1=ngcl+1
cccccc   ngcd=ngci1*ngcj1*ngck1*ngcl1
         ngcd=ngcij*ngck1*ngcl1
c-
         nqi=ndege(itype1)
         nqj=ndege(jtype1)
         nsij=nqi+nqj-1
         lnij=lensm(nsij)
         nqij=nqi
         if(nqj.gt.nqi) nqij=nqj
c
         nqk=ndege(ktype1)
         nql=ndege(ltype1)
         nskl=nqk+nql-1
         lnkl=lensm(nskl)
         nqkl=nqk
         if(nql.gt.nqk) nqkl=nql
c
         mmax=nsij+nskl-1
c
c estimate the price for a block
c
         if ( mmax.le.2 ) then
           nprice=39
           if(mmax.eq.2) nprice=54
           if( max(itype,jtype,ktype,ltype).eq.3) nprice=56
         else
c from prec4neg(36) and XWPQ(15) :
           nprice=51
           nprice=nprice+4*mmax -3
c     
c from Obasai :
c
           isum=0
           do 851 im=1,mmax-2
           isum=isum + (mmax-1-im)*( nfu(im+1)-nfu(im) )
  851      continue
           nfact=30-mmax
           if(mmax.ge.9) nfact=18
           nprice=nprice+nfact*isum
c
           nsxx=nskl
           nqxx=nqij
           if(nskl.gt.nsij) then
             nsxx=nsij
             nqxx=nqkl
           endif
c
c from Tracy :
c
           isum=0
           do 852 kp=2,nsxx
           ndix=nqxx-nsxx+kp
           if(ndix.le.0) ndix=1
           isum=isum+(nfu(mmax+1-kp)-nfu( ndix ))*(nfu(kp+1)-nfu(kp))
  852      continue
           nprice=nprice+7*isum
c
c from assemble ;
c
           nassem=(lnkl-nfu(nqkl))*(lnij-nfu(nqij))
c
           nprice=nprice+nassem
c        
           if(itype.eq.3) nprice=nprice+2*nassem
           if(jtype.eq.3) nprice=nprice+2*nassem
           if(ktype.eq.3) nprice=nprice+2*nassem
           if(ltype.eq.3) nprice=nprice+2*nassem
c
         endif
c
ccccccc  nprice=nprice*icont*jcont*kcont*lcont
         nprice=nprice*ijklc
c
c from shifting of angular momentum :
c
         if(nqij.eq.nsij .and. nqkl.eq.nskl) then
         else if (mmax.gt.2) then
            nq1=nqj
            nq2=nqi
               if(nqj.gt.nqi) then
                  nq1=nqi
                  nq2=nqj
               endif
            isum=0
            do 853 j=2,nq1
            do 853 i=nsij+1-j,nq2,-1
            isum=isum+(nfu(j+1)-nfu(j))*(nfu(i+1)-nfu(i))
  853       continue
            isum=2*isum
c from tfer
            isum=isum*(nfu(nskl+1)-nfu(nqkl))
c
            nq1=nql
            nq2=nqk
               if(nql.gt.nqk) then
                  nq1=nqk
                  nq2=nql
               endif
            ksum=0
            do 854 l=2,nq1
            do 854 k=nskl+1-l,nq2,-1
            ksum=ksum+(nfu(l+1)-nfu(l))*(nfu(k+1)-nfu(k))
  854       continue
            ksum=2*ksum
c from tfer
            ksum=ksum*(nfu(nqi+1)-nfu(nqi))*(nfu(nqj+1)-nfu(nqj))
c
            nprice=nprice + isum+ksum
c
         endif
c
        nprice=nprice/nfijkl
c
        ngen1=16*ngcd
        ngen2=(ngci1+1)*(ngcj1+1)*(ngck1+1)*(ngcl1+1)
        ngen12=ngen1/ngen2
c
        nprice=nprice*ngen12
        nprice=nprice/ngcd
c
c-----> npricet=nprice*nfun
c
        ncost(ikbl)=nprice
c
        if(nprice.gt.maxpr) maxpr=nprice
        if(nprice.lt.minpr) minpr=nprice
c
c  end of price
c
       if(nprint.gt.1 .and. ikbl.le.1000) then
         if(mmax.ne.maxst) then
c* for printing purpose only
           write(8,*)'   '
           maxst=mmax
         endif
         write(8,506) ikbl,mmax,itype,jtype,ktype,ltype,
     *                icont,jcont,kcont,lcont,nfijkl,nquart,nfun,
ccc  *                ngci1,ngcj1,ngck1,ngcl1,nfijkl,nquart,nfun,
     *                nprice
         if(nsupb(ikbl).gt.1) then
            write(8,507) nsupb(ikbl)
         endif
       endif
  100 continue
c
c 506 format(2x,i4,2x,i2,2x,4i2,2x,4i2,3x,i5,4x,i4,1x,i7,3x,i5,2x,i8)
  506 format(2x,i4,2x,i2,2x,4i2,2x,4i2,3x,i5,4x,i6,1x,i8,3x,i6      )
  507 format(5x,' this super-block is split into ',i5,' blocks')
c
c---------------------------------------------------------------
c print statistics : number of 2e-int. and its price disribution
c 
cpnl  if(scftype.eq.'    ') then
cpnl    call statint(inx,iis,jjs,ijbl,npar,nbl2,ncost,minpr,maxpr)
cpnl  endif
c----------------------------------------------------------------
c
      end
c================================================
      subroutine statint(inx,iis,jjs,ijbl,npar,nbl2,ncost,
     *                   minpr,maxpr)
      implicit real*8 (a-h,o-z)
c
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*), ncost(*),npar(*)
      dimension intnu(9)
      dimension xintxx(9)
      data half /0.5d0/
c-------------------------------------------
c number of integrals more expensive than :
c     10% of maximum price
c     25% of maximum price 
c     50% of maximum price
c     75% of maximum price
c    100% of maximum price
c-------------------------------------------
      do 15 ii=1,9
      xintxx(ii)=0
   15 continue
      xprice=dble(maxpr)
      x001=0.001d0*xprice
      x005=0.005d0*xprice
      x010=0.010d0*xprice
      x050=0.050d0*xprice
      x100=0.10d0*xprice
      x250=0.25d0*xprice
      x500=0.50d0*xprice
      x750=0.75d0*xprice
      x990=0.99d0*xprice
c
      intnu(1)=int(x001)
      intnu(2)=int(x005)
      intnu(3)=int(x010)
      intnu(4)=int(x050)
      intnu(5)=int(x100)
      intnu(6)=int(x250)
      intnu(7)=int(x500)
      intnu(8)=int(x750)
      intnu(9)=int(x990)
c-------------------------
c
c contracted integrals counter
c
      xinteg=0
c
      ikbl=0
      do 100 ibl=1,nbl2
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
      inte=inx(3,ics1)
      jnte=inx(3,jcs1)
      ijnte=inte*jnte
      nparij=npar(ibl)
c-
      call parinteg(ibl,nparij,ijbl,nbl2,iis,jjs,ijnte,inte,ijxint)
c-
         do 100 kbl=1,ibl
         ikbl=ikbl+1
         klcs1=ijbl(kbl,1)
         kcs1=iis(klcs1)
         lcs1=jjs(klcs1)
         knte=inx(3,kcs1)
         lnte=inx(3,lcs1)
         klnte=knte*lnte
         nparkl=npar(kbl)
c--
         if(kbl.ne.ibl) then
           call parinteg(kbl,nparkl,ijbl,nbl2,iis,jjs,klnte,knte,klxint)
           xinteg1=dble(ijxint)*dble(klxint)
         else
           xinteg1=dble(ijxint)*dble(ijxint+1)
           xinteg1=xinteg1*half
         endif
c
c  integrals in one super-block = xinteg1
c  total number of integrals:
c
         xinteg=xinteg + xinteg1
c--
         nprice=ncost(ikbl)
c
c count number of integrals more expensive than ..
c
         do 150 ii=1,9
          if(nprice.gt.intnu(ii)) xintxx(ii)=xintxx(ii)+xinteg1
  150    continue
c
  100 continue
c
        call statint1(xintxx,intnu,minpr,maxpr,xinteg)
c
      end
c================================================
      subroutine parinteg(ibl,nparij,ijbl,nbl2,iis,jjs,
     *                    ijnt,int, ijxint)
      implicit real*8 (a-h,o-z)
      dimension iis(*),jjs(*),ijbl(nbl2,*)
c-
       ijeq=0
       do 1001 ijpar=1,nparij
        ijcs=ijbl(ibl,ijpar)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        if(jcs.eq.ics) ijeq=ijeq+1
 1001  continue
       ijxint=(nparij-ijeq)*ijnt + ijeq*int*(int+1)/2
c-
      end
c================================================
      subroutine statint1(xintxx,intnu,minpr,maxpr,xinteg)
      implicit real*8 (a-h,o-z)
      dimension xintxx(9), intnu(9)
c-----------------------------------
c print statistics : number of 2e-int. and its price disribution
c 
c--
      x1024=1024.d0
      xmega=1.d0/(x1024*x1024)
      if(xinteg.gt.0.d0) xint1=1.d0/xinteg
c--
c     write(6,500)
      write(8,500)
  500 format(70('.'))
c     write(6,501)
      write(8,501)
  501 format(' Integral statistics (w/o symmetry and neglect)'/)
c     write(6,502) minpr,maxpr
      write(8,502) minpr,maxpr
  502 format('        computational cost per integral :'/
     *       '        minimum = ',i3,' maximum = ',i6/)
c--
      xbyte=xinteg*8.d0*xmega
      part=100.d0
      ipay=0
c     write(6,504) xinteg,part,ipay     ,xbyte
      write(8,504) xinteg,part,ipay     ,xbyte
c--
c     write(6,500)
      write(8,500)
      part=0.d0
      do 20 ii=1,9
        if(xintxx(ii).gt.0.d0) then
           xintx=xintxx(ii)
           xbyte=xintx*8.d0*xmega
           if(xinteg.gt.0.d0) part=xintx*xint1
           part=part*100.d0
c          write(6,504) xintx,part,intnu(ii),xbyte
           write(8,504) xintx,part,intnu(ii),xbyte
        endif
   20 continue
c     write(6,500)
      write(8,500)
  504 format
     *(f17.0,' (',f5.1,'%) int. stored for ipay=',i6,' take ',
     *                                               f10.2,'MB')
c---
      end
c================================================================
c    subroutines for blocking quartets :
c================================================================
      subroutine blockin4(isbl,ibl,kbl,nbloks,nbls_pnl,
     *                    bl,inx,npard,nbl2)
c---------------------------------------------
c this routine is called for each super-block
c and NBLOKS is a number of blocks belonging
c to the given super-block ISBL .
c (whatever the "super-block" means )
c---------------------------------------------
      implicit real*8 (a-h,o-z)
      logical firstd
c
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /memor1/ iisd,jjsd,ijbld
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      common /memor3/ nblok1d
      common /memors/ nsym,ijshp,isymm
c
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
c
c*  constructe blocks of quartets of contracted shells :
c*  set up vectors : nqrt,nibl,nkbl, nijb,nije, nklb,nkle
c
      call memo2(nbloks)  ! res.mem. for 7 arrays (above)
c
      call blockqur(isbl,ibl,kbl,nbloks,maxqrt,nbl2,bl(ijbld),bl(npard),
     *bl(nqrtd),bl(nibld),bl(nkbld),bl(nijbd),bl(nijed),
     *bl(nklbd),bl(nkled),bl(iisd),bl(jjsd),inx)
c
CCCC  call memo3(maxqrt)
      call memo3(nbls_pnl)
c
c end of the blocking procedure
c
      end
c===============================================================
      subroutine blockqur(isbl,ibl,kbl,nbloks,maxqrt,nbl2,ijbl,npar,
     * nqrt,nibl,nkbl, nijb,nije,nklb,nkle,iis,jjs,inx)
c
      dimension npar(*),ijbl(nbl2,*)
      dimension inx(12,*),iis(*),jjs(*)
c*
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nqrt(*)
c*
c-----------------------------------------------------------
c This is the super-block ISBL  and it is constructed 
c from pair-blocks ibl and kbl .
c
c Constructe NBLOKS blocks of quartets of contracted shells
c belonging to the super-block ISBL
c
      ikbl=0
c
      if(ibl.gt.kbl) then
         call nondiax(nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      else
         call diagonx(nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      endif
c
c It should be  nbloks=ikbl
c
c NOT true for PNL :
cpnl  if(ikbl.ne.nbloks) then
cpnl     write(8,*)' In the Super-block no.=',isbl
cpnl     write(8,*)' number of blocks of contracted quartets'
cpnl     write(8,*)' from Blksizer and Blockqur subroutines'
cpnl     write(8,*)'          is not the same             '
cpnl :   write(8,*)'  execution stoped '
cpnl :   stop600
cpnl  endif
c
cpnl
      nbloks=ikbl
c
      maxqrt=0
      ikblmx=1
      ntotq=0
      do 650 ikbl=1,nbloks
      ntotq=ntotq+nqrt(ikbl)
          if(nqrt(ikbl).gt.maxqrt) then
             maxqrt=nqrt(ikbl)
             ikblmx=ikbl
          endif
  650 continue
c------------------------------------------------------------------
c     write(8,497) isbl,nbloks,ntotq,maxqrt,ikblmx
c 497 format(/'     In the Super-block no.=',i6/ 
c    *        '  number of Blocks of Contracted Quartets =',i6/
c    *        '  total number of contracted quartets =',i7/
c    *        '  maximum number of quartets =',i6,' in block=',i6/)
c
c------------------------------------------------------------------
      end
c==================================================================
      subroutine nondiax(nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
c
      dimension npar(*),ijbl(nbl2,*)
      dimension nqrt(*)
c
      dimension nibl(*), nkbl(*)
      dimension nijb(*),nije(*),nklb(*),nkle(*)
c
      ijpar=npar(ibl)
      klpar=npar(kbl)
      nquart=ijpar*klpar
c
          ijsize=ijpar
          klsize=klpar
          ijdev=1
          kldev=1
          ijrem=0
          klrem=0
c
CNOSPL.   if(nquart.le.maxsize) go to 99
CNOSPLIT : CODE has been removed !!!!
   99 continue
c
          nqrt1=ijsize*klsize
          nqrt2=ijrem*klsize
          nqrt3=ijsize*klrem
          nqrt4=ijrem*klrem
          ijds=ijdev*ijsize
          klds=kldev*klsize
c
          do 100 ij=1,ijdev
          ij1=(ij-1)*ijsize
             do 100 kl=1,kldev
             kl1=(kl-1)*klsize
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ij1+1
             nije(ikbl)=ij1+ijsize
             nklb(ikbl)=kl1+1
             nkle(ikbl)=kl1+klsize
             nqrt(ikbl)=nqrt1
  100      continue
c
           if(ijrem.gt.0) then
             do 200 kl=1,kldev
             kl1=(kl-1)*klsize
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ijds+1
             nije(ikbl)=ijpar
             nklb(ikbl)=kl1+1
             nkle(ikbl)=kl1+klsize
             nqrt(ikbl)=nqrt2
  200        continue
           endif
c
           if(klrem.gt.0) then
             do 300 ij=1,ijdev
             ij1=(ij-1)*ijsize
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ij1+1
             nije(ikbl)=ij1+ijsize
             nklb(ikbl)=klds+1
             nkle(ikbl)=klpar
             nqrt(ikbl)=nqrt3
  300        continue
           endif
c
           if(nqrt4.gt.0) then
             ikbl=ikbl+1
             nibl(ikbl)=ibl
             nkbl(ikbl)=kbl
             nijb(ikbl)=ijds+1
             nije(ikbl)=ijpar
             nklb(ikbl)=klds+1
             nkle(ikbl)=klpar
             nqrt(ikbl)=nqrt4
           endif
c
      end
c=====================================================================
      subroutine diagonx(nqrt,ikbl,ibl,kbl,ijbl,nbl2,npar,
     *                   nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
      dimension npar(*),ijbl(nbl2,*)
      dimension nqrt(*),nibl(*),nkbl(*)
      dimension nijb(*),nije(*),nklb(*),nkle(*)
c
      ijpar=npar(ibl)
      nquart=ijpar*(ijpar+1)/2
c***
cNOSPLIT :
C     if(nquart.le.nbls) then
c
         ikbl=ikbl+1
c
         nibl(ikbl)=ibl
         nkbl(ikbl)=kbl
         nijb(ikbl)=1
         nije(ikbl)=ijpar
         nklb(ikbl)=1
         nkle(ikbl)=0
         nqrt(ikbl)=nquart
c
c****
CNOS  else
CNOSPLIT    CODE has been removed !
c****
CNOS  endif
c---------------
      end
c======================================================================
      subroutine data_save
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memors/ nsym,ijshp,isymm
c
c copy it to :
c
c refer to local bl() as defined in texas_face3 :
c
      common /memor1_R/ iisd_R,jjsd_R,ijbld_R
      common /memor1a_R/ npard_R,ncost_R,nsupb_R, mxsize_R,nasize_R
      common /memors_R/      ijshp_R
c
      iisd_R  =iisd
      jjsd_R  =jjsd
      ijbld_R =ijbld
      npard_R =npard
      ncost_R =ncost
      nsupb_R =nsupb 
      mxsize_R=mxsize
      nasize_R=nasize
      ijshp_R =ijshp
c
      end
c======================================================================
      subroutine whichblk(datnuc,natoms,inx,ncs,iroute)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
      dimension n(0:104)        ! rjh
      integer inx(12,*)
      integer ncs
c     
c     RJH. I think that this routine attempts to decide if texas
c     93 would be faster than 95.  It breaks if the charge
c     on a dummy center is less than one or greater than 104.
c     
c     For now extend n() to start from 0 and ignore atoms with 
c     negative or >104 charge since they probably won't have 
c     any basis functions on them.
c     
      iroute=1
      if(natoms.eq.1) return
c
      iqmax=0
      do 10 iat=1,natoms
         iq=int( datnuc(1,iat) )
         if (iq.ge.0 .and. iq.le.104) then ! rjh
            if(iq.gt.iqmax) iqmax=iq
         endif                  ! rjh
 10   continue
c     
      do 15 iel=0,iqmax         ! rjh
         n(iel)=0
 15   continue
c     
      do 20 iat=1,natoms
         iq=int( datnuc(1,iat) )
         if (iq.ge.0 .and. iq.le.104) then ! rjh
            do 25 iel=0,iqmax   ! rjh
               if(iq.eq.iel) n(iel)=n(iel)+1
 25         continue
         endif                  ! rjh
 20   continue
c     
      ntypes=0
      do 30 iel=0,iqmax         ! rjh
         if(n(iel).gt.0) ntypes=ntypes+1
 30   continue
c     
      if( (natoms/ntypes) .ge.2 ) then
*         write(6,*) ' TEXAS 95 because of many similar atoms'
         iroute=2
      endif
c
c     If there are G or higher functions force the use of texas 93
c     since it will use much less memory
c
      maxtype = 0
      do ics = 1, ncs
         maxtype = max(maxtype,inx(12,ics))
      enddo
      if (maxtype .ge. 8) then
*         write(6,*) ' TEXAS 93 because of high angular momentum'
         iroute = 1
      endif
c     
      end
c======================================================================
