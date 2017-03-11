* $Id$
c===================================================================
      subroutine calcint2(bl,inx, q4,use_q4,
     *                    lab1,lab2,lab3,lab4,fockmat,
     $     ncfunct)
      implicit real*8 (a-h,o-z)
      integer ncfunct(*)        ! Map from txs to pnl basis functions
      logical firstd
      logical use_q4
      character*11 scftype
      character*8 where
c----------------------------------------------------------------
c The variable where shows the type of calculated integrals:
c     where.eq.'buff'  - for ordinary two-el. integrals
c     where.eq.'shif'  - for the GIAO/NMR integral derivativs
c     where.eq.'forc'  - for first geometry derivatives
c     where.eq.'hess'  - for second geometry derivatives
c It can be extanded for other targets.
c----------------------------------------------------------------
      common /runtype/ scftype,where
      common /txs_time/ timprep2,timcalc2,timblok2
      common /number/ acc
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /memor1/ npard,mxsize,nblock1,nblock1_back
      common /memor1b/ nbl2
c
      dimension bl(*)
      dimension inx(12,*)
      dimension fockmat(*)
      dimension lab1(*),lab2(*),lab3(*),lab4(*)
      dimension q4(*)   ! pnl-symmetry factors
c----------------------------------------------------------------
      call txs_second(time_beg)
c----------------------------------------------------------------
c calculations of two-electron integrals in blocks
c Nbl2   - number of blocks of pairs
c
        call blockint(bl,nbl2,inx,bl(nblock1),bl(npard),
     *                          bl(mxsize),
     *                q4,use_q4,
     *                lab1,lab2,lab3,lab4,fockmat,ncfunct)
c
c----------------------------------------------------------------
       call txs_second(time_end)
       timcalc2=timcalc2 + time_end-time_beg
c----------------------------------------------------------------
      end
c===================================================================
      subroutine blockint(bl,nbl2,inx,nblock1,npar,
     *                          mxsize, 
     *                    q4,use_q4,
     *                    lab1,lab2,lab3,lab4,fockmat, ncfunct)
c
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      integer ncfunct(*)        ! txs to pnl basis map
      logical firstd
      logical use_q4
      character*11 scftype
      character*8 where
c----------------------------------------------------------------
c where is 'disk' or 'fock' for SCF : non -,semi- and full-direct
c It can be extanded for special purposes :
c it is 'shif' or 'sf10' for nmr shift calculations
c it is 'forc' for gradient
c----------------------------------------------------------------
c ispec - indicates a special use of the two-el.int.program (pnl)
      common /pnl000/ xbluse,nbluse
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /pnl006/ nsplit,isplit,isbl_split, isbl_part
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c------
      common /superbl/ isupb,ibl,kbl
      common /runtype/ scftype,where
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      dimension bl(*), inx(12,*)
cnew  dimension ncost(*),mxsize(*)
      dimension          mxsize(*)
      dimension nblock1(*),npar(*)
      dimension fockmat(*)
      dimension lab1(*),lab2(*),lab3(*),lab4(*)
      dimension q4(*)
c----------------------------------------------------------------
c setup info needed if SPEC (spec.use of two-el.prog) is ON
c
      if(ispec.ne.0) scftype='**  PNL  **'
c----------------------------------------------------------------
c for checking only :
      call getmem(1,last0)
c     write(8,*)'From Blockint (entry): first free address in BL=',last0
      call retmem(1)
c------------------------------------------------------------------c
c   begining of the loop over blocks of contracted quartets        c
c------------------------------------------------------------------c
c loop over quartet-blocks : nbl4  =nbl2*(nbl2+1)/2
c
      isupb=0
      do 100 ibl=1,nbl2
      call doblock2(ibl,bl(ijblock),idoit)
      if(idoit.eq.0) go to 100
      ibl12=ibl*(ibl-1)/2
        do 200 kbl=1,ibl
        call doblock2(kbl,bl(klblock),jdoit)
        if(jdoit.eq.0) go to 200
        isupb=ibl12+kbl
c
        call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                    isbl_size,isbl_point)
c
        call dosupblk(isupb,bl(isblsize),bl(isblqrts),isbl_point,
     *                isbdo,bl(isbl_split) )
        if(isbdo.eq.0) go to 200
c
        if   (where.eq.'buff'                    .or.
     *        where.eq.'shif'                    .or.
     *        where.eq.'forc'                    .or.
     *        where.eq.'hess'                  ) then
c
           maxsize=mxsize(isupb)
c
c for each super-block do update to the current pnl-sizes :
c
           call get_nparts(isupb,bl(isbl_part),nparts)
c
c check in quartets for a super-block which is too big :
c
           if(nparts.gt.1) then
c-----------> call getmem(isbl_size,isbl_copy) ! <-----------------c
              call memo1_int(isbl_size,isbl_copy)
              call copy_block(isbl_size,bl(isblqrts),isbl_point,
     *                        bl(isbl_copy) )
              call check_in(isupb,isbl_size,bl(isblqrts),isbl_point,
     *                       bl(isbl_copy),maxsize)
           else
              isbl_copy = 0
           endif
c
           call uniq_pairs_1('blockint',
     *                       ibl,kbl,bl,bl(ijpres2),bl(klpres2),
     *                       isbl_size,bl(isblqrts),isbl_point)
ctest
c             WRITE(6,*)
c    *        'SUPER=',isupb,' IBL,KBL=',ibl,kbl,' PNL-SIZE=',isbl_size,
c    *        'MAX-SIZE=',maxsize,' PARTS=',nparts
ctest
c                                              ! two calls to getmem
          call prec2ij(ibl,nblock1,bl,inx)     ! one call to getmem
          call prec2kl(kbl,nblock1,bl,inx)     ! one call to getmem
c
c
c make map from present quartets to unique pairs in a given block
c
           call make_map2uniq(bl,ibl,kbl,
     *                        isbl_size,bl(isblqrts),isbl_point,
     *                        bl(ijpres2),bl(klpres2) )
c
c                                              ! two calls to getmem
c
c in above call there are two memory allocations (getmem=memo1_int)
c output from this call is in the map4uniq common block
c
c Both, isbl_size and maxsize are needed in the following call
c for the case when a block is split into parts. Memory in memo5c should
c be reserved according to min(isbl_size,maxsize) while in memo3 ALWAYS
c according to isbl_size :
c
           call onesuper(isupb,ibl,kbl, isbl_size,maxsize,
     *                   bl,inx,npar,nbl2, where, q4,use_q4,
     *                   lab1,lab2,lab3,lab4 ,fockmat, ncfunct)
c
           call retmem(2)  ! release mem.alloc. in  make_map2uniq
c
c check out calculated quartests :
c
           call check_out(isupb,bl(isblsize),bl(isblqrts),isbl_point,
     *                    bl(isbl_split),bl(isbl_copy),nparts,maxsize,
     *                    ibl,bl(ijblock),igo_back )
c
           call retmem(1)                  ! from prec2kl
           call retmem(1)                  ! from prec2ij
           call retmem(2)                  ! from uniq_pairs_1
c
           if(nparts.gt.1) call retmem(1)     ! release isbl_copy
c
           if(igo_back.eq.1) go to 300
c
      endif
  200 continue
  100 continue
c
  300 continue
c
c------------------------------------------------------------------c
c******* end of the loop over blocks of contracted quartets *******c
c------------------------------------------------------------------c
c    
c for checking only :
      call getmem(1,last1)
      call retmem(1)
        if(last1.ne.last0) then 
           write(8,*)'**        WRONG memory alocations          **'
           write(8,*)' at the beginning of calc.=',last0
           write(8,*)' at the end of cal.       =',last1
           write(8,*)'** Execution has been stoped in BLOCKINT   **'
*           stop100
           call errquit('texas:blockint',0, INT_ERR)
        endif
c*
      end
c===================================================================
      subroutine onesuper(isupb,ibl,kbl, nbls_pnl,maxsize,
     *                    bl,inx,npar,nbl2, where, q4,use_q4,
     *                    lab1,lab2,lab3,lab4 ,fockmat,ncfunct)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      integer ncfunct(*)        ! txs to pnl basis map
      logical use_q4
      character*8 where
c----------------------------------------------------------------
c nbls_size = full pnl-size
c maxsize   = maximum size (can be smaller than above)
c----------------------------------------------------------------
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
      common /memor3/ nblok1d
      dimension bl(*), inx(12,*)
      dimension npar(nbl2)
      dimension fockmat(*)
      dimension lab1(*),lab2(*),lab3(*),lab4(*)
      dimension q4(*)    !  pnl-symmetry factors
c----------------------------------------------------------------
      call mmark
c----------------------------------------------------------------
c Construct blocks of contracted shell quartets for a given
c super-block (isupb) ; set up common /memor2/ and /memor3/ .
c
c----
       call blockin4(isupb,ibl,kbl, nbls_pnl,bl,npar)
c----
c----------------------------------------------------------------
c memory checking :
         call getmem(0,last11)
         call retmem(1)
c----------------------------------------------------------------
c
ccccc    do 200 ikbl=1,nblokx     ! always =1 for pnl
c
           ikbl=1    
           call oneblock(isupb,ikbl,bl, bl(nibld),bl(nkbld),
     *                   bl(nijbd),bl(nijed),bl(nklbd),bl(nkled),
     *                   bl(nblok1d),bl(nqrtd), nbl2,inx,
     *                        where, q4,use_q4,maxsize,
     *                   lab1,lab2,lab3,lab4,fockmat,ncfunct)
c
c----------------------------------------------------------------
c memory checking :
         call getmem(0,last12)
         call retmem(1)
         if(last11.ne.last12) then
           write(6,*)'**        WRONG memory alocations          **'
           write(6,*)' sup-block no=',isupb,'small-block=',ikbl
           write(6,*)'       on entry  =',last11
           write(6,*)'       on exit   =',last12
           write(6,*)'** Execution has been stopped in ONESUPER   **'
*           stop200
           call errquit('texas:onsuper',0, INT_ERR)
         endif
c
cc200    continue
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------
      end
c===================================================================
      subroutine oneblock(isupb,ikbl,bl,nibl,nkbl, nijb,nije,nklb,nkle,
     *                    nblok1,nqrt, nbl2,inx,
     *                         where, q4,use_q4,maxsize,
     *                    lab1,lab2,lab3,lab4,fockmat,map_txs_pnl)
c
      implicit real*8 (a-h,o-z)
      integer map_txs_pnl(*)        ! txs to pnl basis map = ncfunct
      logical first,firstd
      logical use_q4
      character*8 where
c------------------------------------------------------
      common /route/ iroute
c------------------------------------------------------
      common /pnl000/ xbluse,nbluse
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl002/ ncshell,ncfunct,nblock2,integ_n0
      common /pnl003/ nqrtpnl,icstx,jcstx,kcstx,lcstx
      common /pnl004/ isize,jsize,ksize,lsize,itxspnl
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c------------------------------------------------------
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /number/ acc
c
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time0/ tprec2
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time2/ terint,ttrans
      common /time3/ tzeroi,tspeci
ctime
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
#include "texas_lpar.fh"
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
c
c nmr derivatives :
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
      common /memor1/ npard,mxsize,nblock1,nblock1_back
c
      dimension bl(*), inx(12,*)
      dimension nibl(*),nkbl(*),nijb(*),nije(*),nklb(*),nkle(*)
      dimension nblok1(2,*),nqrt(*)
c
      dimension fockmat(*)
      dimension lab1(*),lab2(*),lab3(*),lab4(*)
      dimension q4(*)
c-----------------------------------------------------------
c     if(icheck.gt.0) return
c-----------------------------------------------------------
c  number of ij and kl pairs (npij,npkl)
c  (zero for npkl means that the block is diagonal)
      ibl=nibl(ikbl)
      kbl=nkbl(ikbl)
c
      nijbeg=nijb(ikbl)
      nijend=nije(ikbl)
      npij=nijend-nijbeg+1
c
      nklbeg=nklb(ikbl)
      nklend=nkle(ikbl)
      npkl=nklend-nklbeg+1
      if(nklend.eq.0) npkl=0
      npklx=npkl
      if(npkl.eq.0) npklx=npij
c-----------------------------------------------------------
c Get pnl_size of the ISUPB block and a pointer to its qrts. 
c
      call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                  isbl_size,isbl_point)
c-----------------------------------------------------------
c get the NBLOK1 array and pnl-block-size nbls_pnl=nqrt(ikbl)
c (for pnl the ikbl is ALWAYS equal 1). The ISYMM(ii) is zero out
c which is needed if a sup-block is split
c
      call block1_pnl(isupb,isbl_size,bl(isblqrts),isbl_point,maxsize,
     *                bl(ijpres2),bl(klpres2),bl(isymm),nblok1,nqrt)
ccc   nbls_pnl=nqrt(ikbl)
      nbls_pnl=nqrt( 1  )
c
c
      nbls=nbls_pnl
c-----------------------------------------------------------
c When the SPEC option is on
c
      nblspec=0
      if(ispec.ne.0) then
c
c2002    ntotal=ncs*(ncs+1)/2
c        ntotal=ntotal*(ntotal+1)/2
c
c        call txs_second(tim_1doq)
c
c        needed when a given super-block is split in parts :
c
         call doquarts(isupb,isbl_size,bl(isblqrts),isbl_point,
     *                                         bl(isymm),nblspec)
ccc  *                 bl(ijpres2),bl(klpres2),bl(isymm),nblspec)
c
c nblspec= could be the same or smaller than nbls(pnl)
c is used ONLY for performance estimation 
c-----------------------------------------------------------
c
c        call txs_second(tim_2doq)
c        time_doqrt=time_doqrt + tim_2doq-tim_1doq
c
c efficiency :
c
         nospec=nospec+nblspec
c
c2002
c        xbluse=xbluse+dble(nblspec)/dble(nbls)
c        nbluse=nbluse+1
c
      endif
c-----------------------------------------------------------
c return PNL-check-run 
c     if(icheck.gt.0) return
c-----------------------------------------------------------
c  first quartet of shells :
cnew
      call get_ics_jcs(bl(nblock1),ibl,  1  ,ics1,jcs1)
      call get_ics_jcs(bl(nblock1),kbl,  1  ,kcs1,lcs1)
cnew
      ijcs1=ics1*(ics1-1)/2+jcs1
      klcs1=kcs1*(kcs1-1)/2+lcs1
cnew
c-----------------------------------------------------------
c  set up the type (common types) and length of shells and
c  information concerning contractions (commons contr,gcont) :
c
      call txs_shells(inx,ics1,jcs1,kcs1,lcs1)
c-----------------------------------------------------------
c (setup the obarai and shell commons) :
c
c
      call iobara(itype1,jtype1,ktype1,ltype1,where)
c-----------------------------------------------------------
c for memory handling 
c
      nfumax=nfu(mmax)
      mmax1=mmax-1
      if(mmax1.eq.0) mmax1=1
c
c98   if(nsij.ge.nskl) then
c98      nhabcd=nkl_uniqe
c98      nfha=nfumax*nhabcd*lckl
c98   else
c98      nhabcd=nij_uniqe
c98      nfha=nfumax*nhabcd*lcij
c98   endif
c98
c
c this is needed ONLY for TRACY's recursive 
c and NOW do it ALWAYS as follows :
c
         nhabcd=nkl_uniqe
         nfha=nfumax*nhabcd*lckl
c-----------------------------------------------------------
c memory handling 
c
c 1: for pairs precalculations (20 quantities)
c (for whole block of contracted quartets and 
c  all primitive quartets belonging to this block)
c 
c 2: and quartets precalculations (12 quantities)
c (for whole block of contracted quartets and 
c        one primitive quartet )
c 
c 3: and for individual shells (8 quantities)
c ( aa, bb, cc, dd - exponents )
c ( cis,cjs,cks,cls -contr coef.)
c
c nbls   = a full size of the pnl-requested super-block
c-----------------------------------------------------------
      IF( iroute.eq.1 ) THEN
         call memo5a_1(nij_uniqe ,mmax1)
         call memo5b_1(nkl_uniqe,mmax1)
         call memo5c_1(bl, nbls,mmax1,nij_uniqe,nkl_uniqe,nfha,nfumax)
c
c        50 or 56 calls of getmem
         ncalls=50
         if(ngcd.gt.1) ncalls=56
      ELSE
         call memo5a_2(nij_uniqe ,mmax1)
         call memo5b_2(nkl_uniqe,mmax1)
         call memo5c_2(nbls,mmax1,nij_uniqe,nkl_uniqe,nfumax)
c
c        43 or 47 calls of getmem plus 3 calls for where='forc'or'hess'
c         
         ncalls=43
         if(ngcd.gt.1) ncalls=47
         if(where.eq.'forc' .or. where.eq.'hess') ncalls=ncalls+3
      ENDIF
c-----------------------------------------------------------
c Perform pre-calculations for pairs of contracted sells
c
      call txs_second(tprec2b)
c
      IF( iroute.eq.1 ) THEN
         call precalc2_1(isupb,bl,mmax,mmax1,nhabcd,nfumax,nbl2,nbls,
     *                   inx,ibl,kbl,npkl)
      ELSE
         call precalc2_2(isupb,bl,mmax,mmax1, nfumax, nbl2,nbls,
     *                   inx,ibl,kbl,npkl)
      ENDIF
c-----------------------------------------------------------
c for NMR derivatives 
c
c   NOT READY for PNL : inuc, ibas refer to dbl_mb()
c
cnmr  if(where.eq.'shif') then
cnmr     call memo6(npij,npklx)        !   2 calls of getmem
cnmr     call getmem(npij ,ijcent)
cnmr     call getmem(npklx,klcent)
cnmr     call precal2d(bl(inuc+1),iis,jjs,inx, npij,npklx,npkl,
cnmr *           ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
cnmr *           bl(ixyab),bl(ixycd),
cnmr *           bl(isymm),bl(ijcent),bl(klcent) )
cnmr     call retmem(2)
cnmr  endif
c-----------------------------------------------------------
c neglect some of contracted quartets on the cont. shell level :
c isym(ijkl)=0 or 1 , 0 means that ijkl quartet is neglected.
c
c (neglect on the contected shell level is not performed for PNL)
c
c  NO call to the CSHNEG routine !!!!
c-----------------------------------------------------------
c At this moment some of c.s.quartets may not be present if this run
c is with MORE_INT=.true. which means that  this super-block was split
c and some of its quartets were calculated in a previous call.
c
c Thus, reduce the block size from  nbls ---> nblsp (present)
c (only present c.s.quartets are considered)
c        
c define idx1 (ijpar <--- ijkl quartets mapping)
c    and idx2 (klpar <--- ijkl quartets mapping)
c
      call indexp_pnl(isupb,isbl_size,bl(isblqrts),isbl_point,
     *                bl(map_ij),bl(map_kl),
     *                bl(idx1),bl(idx2),bl(isymm),bl(indxp),nblsp)
c
c-----------------------------------------------------------
         call txs_second(tprec2e)
         tprec2=tprec2+tprec2e-tprec2b
c-----------------------------------------------------------
c FROM hereafter the block-size is reduced nbls-->nblsp
c
      if(nblsp.eq.0) then 
c       write(8,*)' From oneblock block=',ikbl,'ngcd=',ngcd
c       write(8,*)' nbls--->nblsp ',nbls,nblsp
c ncalls-1 comes from the fact that no memory has been res.
c for buf in this case :
        ncalls=ncalls-1
        go to 110
      endif
c
      nopres=nopres+nblsp
cccc  nblsrem=nbls
      nbls=nblsp
c-----------------------------------------------------------
c calculate special integrals with mmax<=2 
c (ssss), (psss),(lsss) etc.
c
      if(mmax.le.2) then
         call txs_second(tspec1)
       if(icheck.eq.0) then
         IF( iroute.eq.1 ) THEN
cnew       call erintsp_1(isupb,bl,first,nbls,acc, ikbl,npij,npklx,npkl)
           call erintsp_1(isupb,bl,first,nbls,acc, ikbl,           npkl)
         ELSE
cnew       call erintsp_2(isupb,bl,first,nbls,acc, ikbl,npij,npklx,npkl)
           call erintsp_2(isupb,bl,first,nbls,acc, ikbl,           npkl)
         ENDIF
       endif
         call txs_second(tspec2)
         tspeci=tspeci+tspec2-tspec1
         if( first ) go to 110
         go to 120
      endif
c
c-----------------------------------------------------------
c calculate integrals with MMAX > 2
c
ctime
      call txs_second(teri)
c
      if(icheck.eq.0) then
         IF( iroute.eq.1 ) THEN
cnew     call erinteg_1(isupb,bl,first,nbls,acc,ikbl,npij,npklx,npkl,
c98      call erinteg_1(isupb,bl,first,nbls,acc,ikbl,           npkl,
c98  *                  lobsa,immax,kmmax, where)
         call erinteg_1(isupb,bl,first,nbls,acc,ikbl,           npkl,
     *                                     where)
         ELSE
cnew     call erinteg_2(isupb,bl,first,nbls,acc,ikbl,npij,npklx,npkl,
c98      call erinteg_2(isupb,bl,first,nbls,acc,ikbl,           npkl,
c98  *                  lobsa,immax,kmmax, where)
         call erinteg_2(isupb,bl,first,nbls,acc,ikbl,           npkl,
     *                                     where)
         ENDIF
      endif
c
c-----------------------------------------------------------
      call txs_second(terie)
      terint=terint+(terie-teri)
c
      if( first ) go to 110
c
c-----------------------------------------------------------
c transformatin d6->d5, f10->f7, g15->g9, h21->h11, i28->i13
c
c do we need to transform ?
c
      lenall=ilen+jlen+klen+llen
      lenall1=ilen1+jlen1+klen1+llen1
      if(lenall1-lenall.ne.0 ) then
c>>>> if( cart_2_sphe ) then
         call txs_second(trans1)
c        ------------------------------------------------------
c         get spherical harmonics PNL---types : ityps - ltyps :
c         using Texas cartesian types itype1-ltype1
c
c         cartesian basis set shell sizes : ilen1 - llen1 :
c         spherical harmonics shell sizes : ilen  - llen  :
c         spherical harmonics PNL---types : ityps - ltyps :
c
          call get_spher_pnl_type(itype1,ityps)
          call get_spher_pnl_type(jtype1,jtyps)
          call get_spher_pnl_type(ktype1,ktyps)
          call get_spher_pnl_type(ltype1,ltyps)
c        -------------------------------------------------------
c
c
         mnbls=nbls
         nbuf=ibuf
         if(where.eq.'shif') then
c        --- nmr derivatives ----
           mnbls=6*nbls
           nbuf=ibuf+ngcd*nbls*lnijkl
         endif
         if(where.eq.'forc') then
c        --- gradient derivatives ----
           mnbls=9*nbls
           nbuf=ibuf
         endif
         if(where.eq.'hess') then
c        --- hessian  derivatives ----
           mnbls=45*nbls   ! for second derivatives only
cccc       mnbls=54*nbls   ! for second and first derivatives together
           nbuf=ibuf
         endif
c
         incrt=mnbls*lnijkl
c
           do 130 iqu=1,ngcd                   ! over gen.cont.
           jbuf=nbuf+(iqu-1)*incrt   
           call transfor(bl,mnbls,jbuf,
     *                   ityps , jtyps, ktyps, ltyps,
     *                   ilen1 , jlen1, klen1, llen1,
     *                   ilen  , jlen , klen , llen )
  130      continue
  125    continue
         call txs_second(trans2)
         ttrans=ttrans+trans2-trans1
      endif
c
c end of transformation
c-----------------------------------------------------------
c
  120 continue
c-----------------------------------------------------------
      if(icheck.gt.0) go to 110
c-----------------------------------------------------------
        call txs_second(tdest)
c--------------------------------------------------------------
c Final use of integrals : 
c
c  For Battelle PNL : special request (SPEC option in the input)
c
c  5. put integrals from a selected quartet of shells into
c     a buffer 
c        (when where='buff'  - call destbuf or destbul)
c        (when where='forc'  - call destduf or destdul)
c
c--------------------------------------------------------------
c
      if ( where.eq.'buff') then
c-5>  put two-el.int. into the buffer - only a set of integrals 
c         fockmat is now the buffer for integrals
c
       if(ispec.eq.1) then
c       integrals without labels (all of them):
        call destbuf(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat,bl(itxspnl), q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm))
       endif
       if(ispec.eq.2) then
c       integrals (non-zeros) with labels 
        call destbul(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat, lab1,lab2,lab3,lab4, q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm),bl(iqorder),map_txs_pnl)
       endif
      endif
c 
      if ( where.eq.'forc') then
c         gradient derivatives :
c
       if(ispec.eq.1) then
c       integrals without labels (all of them):
        call destduf(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat,bl(itxspnl), q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm), bl(iqorder))
       endif
       if(ispec.eq.2) then
c       integrals (non-zeros) with labels 
        call destdul(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat, lab1,lab2,lab3,lab4, q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm),bl(iqorder),map_txs_pnl)
       endif
      endif
c 
      if ( where.eq.'hess') then
c         hessian  derivatives :
c
       if(ispec.eq.1) then
c       integrals without labels (all of them):
        call desthuf(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat,bl(itxspnl), q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm), bl(iqorder))
       endif
       if(ispec.eq.2) then
c       integrals (non-zeros) with labels 
        call desthul(ikbl,nbls,nblok1,        ncs,inx,
     *       bl(ibuf),fockmat, lab1,lab2,lab3,lab4, q4,use_q4,
     *       bl(icfg),bl(jcfg),bl(kcfg),bl(lcfg),ngcd,lnijkl,
     *       bl(indxp) ,bl(isymm),bl(iqorder),map_txs_pnl)
       endif
      endif
c-----------------------------------------------------------
c
        call txs_second(tdeste)
        tdesti=tdesti+(tdeste-tdest)
c
  110 continue
c
c-----------------------------------------------------------
c release memory at the end of a given block :
c
c nmr deriv:
c
      if(where.eq.'shif') ncalls=ncalls+2
c
      if(icheck.eq.0) then
        call retmem(ncalls+1)   ! +1 to release ibuf
      else
        call retmem(ncalls)
      endif
c----------------------------------------
      end
c===================================================================
      subroutine erinteg_1(isupb,bl,first,nbls,acc,ikbl,npkl, where)
c 
c character variable WHERE is needed for derivatives only
c
      implicit real*8 (a-h,o-z)
      logical first
      character*8 where
c
      common /primij/ iabprim, ijdim ,ijpar1 
      common /primkl/ kabprim, kldim ,klpar1 
c
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
ctime
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
#include "texas_lpar.fh"
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor3/ nblok1d
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
c nmr der.
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
c-----------------------------------------------------------------
c  memory handling for a block 
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(bl, nbls, l11,l12,mem2,igmcnt)
c
      nfumax=nfu(mmax)
      mmax1=mmax-1
      nfumax1=nfumax
      narray=1
      nqijr=nqij
      nqklr=nqkl
      if(where.eq.'shif') nfumax1=nfu(mmax1)
      if(where.eq.'forc') then 
         nfumax1=nfu(mmax1)
         narray=4
         nqij=nqij-1
         nqkl=nqkl-1
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
      if(where.eq.'hess') then 
         nfumax1=nfu(mmax1-1)
         narray=10
         nqij=nqij-2
         nqkl=nqkl-2
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
c-----------------------------------------------------------------
      first=.true.
c-----------------------------------------------------------------
      igcoet=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcij)
        call getmem(ngckl*nbls,igckl)
        call getmem(ngcd*nbls          ,igcoet)
      endif
c
      if(where.eq.'forc' .or. where.eq.'hess') then
        call getmem(nbls,iaax)
        call getmem(nbls,ibbx)
        call getmem(nbls,iccx)
      endif
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
      call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                  isbl_size,isbl_point)
c
      call getmem(isbl_size,map_n0)
      call get_nbls_n0(isbl_size,isbl_point,bl(isblqrts),
     *                 bl(map_n0), nbls_n0)
c-----------------------------------------------------------------
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      if(ngcd.gt.1) then
        call gcparij(nbls, bl(idx1),nij_uniqe,
     *               l1,l2,ngci1,ngcj1,ngcij,
     *               bl(igci),bl(igcj), bl(igcij))
      endif
      lcij=lcij+1
         lckl=0
         nblsr=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         if(ngcd.gt.1) then
           call gcparkl(nbls, bl(idx2),nkl_uniqe,
     *                  l3,l4,ngck1,ngcl1,ngckl,
     *                  bl(igck),bl(igcl), bl(igckl))
         endif
         lckl=lckl+1
c
c calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
ctime
      call txs_second(t4b)
c
c2002 call prec4neg_1(isbl_size,isbl_point,bl(isblqrts),npkl,
c
      call prec4neg_1(nbls_n0,bl(map_n0), lcij,lckl,lc12,lc34, 
     *     bl(ieab),bl(iecd),bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *                bl(ixp),bl(ixq),
     *                bl(map_ij),bl(map_kl),nij_uniqe,nkl_uniqe,
     *     bl(irppq),bl(irhoapb),bl(irhocpd),bl(irys),bl(iconst),
     *                          nbls1,bl(indx) )
c
       call txs_second(t4e)
       tpre4n=tpre4n+(t4e-t4b)
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c  note :
c from here a given block (nbls) was reduced to (nbls1) ***
c              for given l1,l2,l3,l4 !!!
c                      and
c rho, rppq, rhoapb, rhocpd, rys, and const
c  have dimension nbls1 (not  (nbls) )
c-----------------------------------------------------------------
c substitute zero for all needed bufors for quartets
c which do not appear first time after neglect
c
         if(first) then
ctime
              call txs_second(tzer)
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
ckw      call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd)
         call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd,narray)
            endif
              call txs_second(tzero)
              tzeroi=tzeroi+(tzero-tzer)
         endif
c-----------------------------------------------------------------
c calculate products of contraction coefficients 
c for general contracted basis set
c
      if(ngcd.GT.1) then
         call gcqijkl(nbls,nbls1, bl(indx),
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(indgc),bl(igcoef),
     *                bl(igcij),ngcij, bl(igckl),ngckl)
      endif
c-----------------------------------------------------------------
         call txs_second(txwpqb)
c
         call xwpq_1(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     *               lcij,lckl,bl(idx1),bl(idx2),bl(indx),
     *               nij_uniqe,nkl_uniqe,
     *               bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
     *               bl(itxab),bl(itxcd),
     *               bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) )
c
c-----------------------------------------------------------------
c calculate abcd(i)=apb/cpd or vice versa depending on stability :
c
         call abcd_1(nbls1, l1,l2, lcij,lckl,
     *               bl(idx1),bl(idx2),bl(indx),nij_uniqe,nkl_uniqe,
     *               bl(iabcd), bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb),
     *               bl(iaa),bl(ibb),bl(iconst),bl(igcoef),ngcd)
c-----------------------------------------------------------------
c convert:  abnia,cdnia, xpn,xqn, habcd 
c from pair to quartet (nbls1) quantities 
c for trobsa 
c
cxxxxij
ctime
         call txs_second(tc1xb)
         txwpq=txwpq+(tc1xb-txwpqb)
c
      if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
         call conv1x_1(nbls1,mmax1,nij_uniqe ,lcij, bl(idx1),bl(indx),
     *                 bl(iabnia), bl(ixpn), bl(iabnix),bl(ixpnx)  )
        if(where.eq.'forc' .or. where.eq.'hess') then
           call conv1der_1(nbls1,nij_uniqe,l1,bl(idx1),bl(indx),
     *                     bl(iaa),bl(iaax) )
           call conv1der_1(nbls1,nij_uniqe,l2,bl(idx1),bl(indx),
     *                     bl(ibb),bl(ibbx) )
        endif
      endif
cxxxxkl
         call conv1x_1(nbls1,mmax1,nkl_uniqe,lckl, bl(idx2),bl(indx),
     *                 bl(icdnia), bl(ixqn), bl(icdnix),bl(ixqnx)  )
        if(where.eq.'forc' .or. where.eq.'hess') then
           call conv1der_1(nbls1,nkl_uniqe,l3,bl(idx2),bl(indx),
     *                     bl(icc),bl(iccx) )
        endif
ctime
         call txs_second(tc1xe)
         tconv1=tconv1+(tc1xe-tc1xb)
cxxxx2
c98    if(nsij.ge.nskl) then
         call conv2x(nbls1,nfumax1,nkl_uniqe,lckl, bl(idx2),bl(indx),
     *               bl(ihabcd),nfumax, bl(ihabcdx)  )
c98    else
c98      if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
c98      call conv2x(nbls1,nfumax1,nij_uniqe,lcij, bl(idx1),bl(indx),
c98  *               bl(ihabcd),nfumax, bl(ihabcdx)  )
c98      endif
c98    endif
ctime
         call txs_second(tc2xe)
         tconv2=tconv2+(tc2xe-tc1xe)
c
c98      call trobsa(bl,nbls1,l11,l12,mem2,immax,kmmax,lobsa)
         call trobsa_1(bl,nbls1,l11,l12,mem2)
c
ctime
        call txs_second(ttrob)
        ttrobs=ttrobs+(ttrob-tc2xe)
cxxxxxx
      if(ngcd.eq.1) then
         call assemblx(bl,first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,lcij,lckl,nij_uniqe,nkl_uniqe)
cnew *                 l1,l2,l3,l4,lcij,lckl,npij,npklx)
      else
c98      call gcqijkl(nbls,nbls1, bl(indx),
c98  *                     ngci1,ngcj1,ngck1,ngcl1,ngcd,
c98  *          bl(indgc),bl(igcoef),
c98  *          bl(igcij),ngcij, bl(igckl),ngckl)
c
c transpose gcoef(ngcd,nbls) into gcoet(nbls,ngcd)
         call trspmo(bl(igcoef),ngcd,bl(igcoet),nbls)
c
         call assemblg(bl,first,nbls,nbls1,lnij,lnkl,ngcd,igcoet)
      endif
ctime
         call txs_second(tasse)
         tassem=tassem+(tasse-ttrob)
cxxxxx
c
         nblsr=nbls1
c
   34    continue        
   12 continue
c
c grad  &  hess derivatives:
c
      nqij=nqijr
      nqkl=nqklr
c
c-----------------------------------------------------------------
c release memory :
c
      call retmem(1)    ! allocated for map_n0
c
      if(where.eq.'forc' .or. where.eq.'hess') call retmem(3)
c
      if(ngcd.gt.1) call retmem(3)
c
      if( first ) then
          call retmem(igmcnt)                         
          return
      endif
c
c release memory from trobsa
c
      call retmem(3)
c
c------------------------------------
c for NMR, GRADIENT and HESSIAN derivatives :
c
      lnijr=lnij
      lnklr=lnkl
c------------------------------------
c reserves memory for amshift
c
      call memo4b(bl, nbls,igmcnt1)
c
c------------------------------------
c for NMR derivatives :
c
      if(where.eq.'shif') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=1
        klderiv=1
        call iobarb(ijderiv,klderiv)
c
c construct nmr/giao derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call shift_der(nbls,lnijr,lnklr,nij_uniqe,nkl_uniqe,ngcd,
     *                idx1,idx2, ixab,ixcd,ixyab,ixycd)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c for GRADIENT derivatives :
c
      if(where.eq.'forc') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=1
        klderiv=1
        call iobarb(ijderiv,klderiv)
c
c construct gradient derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call force_der(bl,nbls,lnijr,lnklr,nij_uniqe,ngcd,idx1,ixab)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c for HESSIAN  derivatives :
c
      if(where.eq.'hess') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=2
        klderiv=2
        call iobarb(ijderiv,klderiv)
c
c construct hessian derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call hessian_der(bl,nbls,lnijr,lnklr,nij_uniqe,ngcd,idx1,ixab)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c*  shift angular momentum (a->b, c->d) :
ctime   
           call txs_second(tamsb)
cnew    call amshift(bl,nbls,lnij,lnkl,npij,npklx,ngcd)
        call amshift(bl,nbls,lnij,lnkl,nij_uniqe,nkl_uniqe,ngcd)
           call txs_second(tamse)
           tamshf=tamshf+(tamse-tamsb)
c-----
      call retmem(igmcnt+igmcnt1-3)                         
c-----
c
      end
c===================================================================
      subroutine erintsp_1(isupb,bl,first,nbls,acc,ikbl,npkl)
      implicit real*8 (a-h,o-z)
      logical first
c
      common /primij/ iabprim, ijdim ,ijpar1 
      common /primkl/ kabprim, kldim ,klpar1 
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
#include "texas_lpar.fh"
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common /memor3/ nblok1d
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memors/ nsym,ijshp,isymm
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
c-----------------------------------------------------------------
c  memory handling for a block 
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(bl, nbls, l11,l12,mem2,igmcnt)
c-----------------------------------------------------------------
      first=.true.
c-----------------------------------------------------------------
      igcoet=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcij)
        call getmem(ngckl*nbls,igckl)
        call getmem(ngcd*nbls       ,igcoet)
      endif
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
      call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                  isbl_size,isbl_point)
c
      call getmem(isbl_size,map_n0)
      call get_nbls_n0(isbl_size,isbl_point,bl(isblqrts),
     *                 bl(map_n0), nbls_n0)
c-----------------------------------------------------------------
c
c* loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      if(ngcd.gt.1) then
        call gcparij(nbls, bl(idx1),nij_uniqe,
     *               l1,l2,ngci1,ngcj1,ngcij,
     *               bl(igci),bl(igcj), bl(igcij))
      endif
      lcij=lcij+1
         lckl=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         if(ngcd.gt.1) then
           call gcparkl(nbls, bl(idx2),nkl_uniqe,
     *                  l3,l4,ngck1,ngcl1,ngckl,
     *                  bl(igck),bl(igcl), bl(igckl))
      endif
         lckl=lckl+1
c
c* calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
      call precspec_1(nbls_n0,bl(map_n0), lcij,lckl,lc12,lc34, 
     *                bl(idx1),bl(idx2),
     *     bl(ieab),bl(iecd), bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq),bl(i1apb),bl(i1cpd),bl(itxab),bl(itxcd),
     *     bl(map_ij),bl(map_kl),nij_uniqe,nkl_uniqe,
     *     bl(irys),bl(iconst),bl(ixwp),bl(ixwq), nbls1,bl(indx) )
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c  note :
c* from here a given block (nbls) was reduced to (nbls1) ***
c*              for given l1,l2,l3,l4 !!!
c*                      and
c* rho, rppq, rhoapb, rhocpd, rys, and const
c  have dimension nbls1 (not  (nbls) )
c******************************************************
c* substitute zero for all needed bufors for quartets
c* which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
ckw      call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd)
         call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd,1)
            endif
         endif
c
c********
c
c*    special code for (ss ss) and (xs ss), (sx ss), (ss xs), (ss sx)
c*          itegrals ( x= p or l )
c
      if(ngcd.eq.1) then
          call specase_1(bl,first,nbls,nbls1,bl(indx),bl(idx1),bl(idx2),
     *          nij_uniqe,nkl_uniqe,l1,l2,l3,l4,
     *          bl(icis),bl(icjs),bl(icks),bl(icls),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),bl(irr1) )
      else
         call gcqijkl(nbls,nbls1, bl(indx),
     *                     ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *          bl(indgc),bl(igcoef),
     *          bl(igcij),ngcij, bl(igckl),ngckl)
c
c transpose gcoef(ngcd,nbls) into gcoet(nbls,ngcd)
         call trspmo(bl(igcoef),ngcd,bl(igcoet),nbls)
c
          call specasg(bl,first,nbls,nbls1, bl(indx),bl(idx1),bl(idx2),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),
     *          ngcd,bl(indgc),bl(igcoet),lnijkl)
      endif
c
   34    continue        
   12 continue
c
c-----------------------------------------------------------------
c release memory 
c
      call retmem(1)    ! allocated for map_n0
      if(ngcd.gt.1) call retmem(3)
c----
      call retmem(igmcnt)                         
c
      end
c===================================================================
      subroutine erinteg_2(isupb,bl,first,nbls,acc,ikbl,npkl, where)
c 
c character variable WHERE is needed for derivatives only
c
      implicit real*8 (a-h,o-z)
      logical first
      character*8 where
c
      common /primij/ iabprim, ijdim ,ijpar1 
      common /primkl/ kabprim, kldim ,klpar1 
c
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
ctime
      common /timex/ tconv1,tconv2,ttrobs
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
ctime
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
#include "texas_lpar.fh"
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor3/ nblok1d
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
c nmr der.
      common /memor6/ ixyab,ixycd
c
      common /memors/ nsym,ijshp,isymm
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
        call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                    isbl_size,isbl_point)
c-----------------------------------------------------------------
c  memory handling for a block 
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(bl, nbls, l11,l12,mem2,igmcnt)
c-----------------------------------------------------------------
      nfumax=nfu(mmax)
      mmax1=mmax-1
      nfumax1=nfumax
      narray=1
      nqijr=nqij
      nqklr=nqkl
      if(where.eq.'shif') nfumax1=nfu(mmax1)
      if(where.eq.'forc') then 
         nfumax1=nfu(mmax1)
         narray=4
         nqij=nqij-1
         nqkl=nqkl-1
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
      if(where.eq.'hess') then 
         nfumax1=nfu(mmax1-1)
         narray=10
         nqij=nqij-2
         nqkl=nqkl-2
         if(nqij.le.0) nqij=1
         if(nqkl.le.0) nqkl=1
      endif
c-----------------------------------------------------------------
      first=.true.
c-----------------------------------------------------------------
      igcoet=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcijx)
        call getmem(ngckl*nbls,igcklx)
        call getmem(ngcd*nbls          ,igcoet)
      endif
      if(where.eq.'forc' .or. where.eq.'hess') then
        call getmem(nbls,iaax)
        call getmem(nbls,ibbx)
        call getmem(nbls,iccx)
      endif
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
      call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                  isbl_size,isbl_point)
c
      call getmem(isbl_size,map_n0)
      call get_nbls_n0(isbl_size,isbl_point,bl(isblqrts),
     *                 bl(map_n0), nbls_n0)
c-----------------------------------------------------------------
c
c loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      lcij=lcij+1
      if(ngcd.gt.1) then
        call gcpairs(nbls,lcij,ngcij,bl(igcij),bl(igcijx) )
      endif
         lckl=0
         nblsr=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         lckl=lckl+1
         if(ngcd.gt.1) then
          call gcpairs(nbls,lckl,ngckl,bl(igckl),bl(igcklx))
         endif
c
c* calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
ctime
      call txs_second(t4b)
c
c2002 call prec4neg_2(isbl_size,isbl_point,bl(isblqrts),npkl,
c
      call prec4neg_2(nbls_n0,bl(map_n0),
     *                lcij,lckl, lc12,lc34, 
     *     bl(ieab),bl(iecd),bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq),
     *     bl(map_ij),bl(map_kl),nij_uniqe,nkl_uniqe,
     *     bl(irppq),bl(irhoapb),bl(irhocpd),bl(irys),bl(iconst),
     *                          nbls1,bl(indx) )
c
ctime
       call txs_second(t4e)
       tpre4n=tpre4n+(t4e-t4b)
c
c 
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c  note :
c* from here a given block (nbls) was reduced to (nbls1) ***
c*              for given l1,l2,l3,l4 !!!
c*                      and
c* rho, rppq, rhoapb, rhocpd, rys, and const
c  have dimension nbls1 (not  (nbls) )
c******************************************************
c* substitute zero for all needed bufors for quartets
c* which do not appear first time after neglect
c
         if(first) then
ctime
              call txs_second(tzer)
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
         call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd,narray)
            endif
              call txs_second(tzero)
              tzeroi=tzeroi+(tzero-tzer)
         endif
c----------------------------------------------------------------
c calculate products of contraction coefficients 
c for general contracted basis set
c
      if(ngcd.GT.1) then
         call gcquart(nbls,nbls1, bl(indx),
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(igcijx),ngcij, bl(igcklx),ngckl,
     *                bl(indgc),bl(igcoef)  )  
      endif
c----------------------------------------------------------------
c needed for obara-saika recursives :
c calculate XWP, XWQ and p1234 where 
c p1234= - [ b_exp*(A-B) + d_exp*(C-D) ]/(c_exp+d_exp) 
c
      call txs_second(txwpqb)
c
      call xwpq_2(nbls1,bl(ixwp),bl(ixwq),bl(ip1234),
     *            lc12,lc34,lcij,lckl,
     *            bl(idx1),bl(idx2),bl(indx),
     *            nij_uniqe,nkl_uniqe,
     *            bl(irppq),bl(ixp),bl(ixq),bl(ixpp),bl(ixqq),
c98  *            bl(itxab),bl(itxcd),bl(iabcd),
     *            bl(itxab),bl(itxcd),
     *            bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb) )
c----------------------------------------------------------------
c needed for tracy's recursives :
c depending on numerical stability calculate 
c abcd=(a+b)/(c+d) or (c+d)/(a+b) 
c
      call abcd_2(nbls1,lcij,lckl,bl(indx),
     *            bl(iabcd), bl(iapb),bl(i1cpd),bl(icpd),bl(i1apb),
     *            bl(iconst),bl(igcoef),ngcd)
c----------------------------------------------------------------
c convert:  abnia,cdnia, xpn,xqn, habcd 
c from pair to quartet (nbls1) quantities 
c for trobsa 
c
c  ij
ctime
         call txs_second(tc1xb)
         txwpq=txwpq+(tc1xb-txwpqb)
c
      if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
         call conv1x_2(nbls1,mmax1,nij_uniqe,lcij, bl(idx1),bl(indx),
     *                 bl(ixpn), bl(ixpnx)  )
                       iabnix=iabnia+(lcij-1)*mmax1
        if(where.eq.'forc' .or. where.eq.'hess') then
           call conv1der_2(nbls1,nij_uniqe,l1,bl(idx1),bl(indx),
     *                     bl(iaa),bl(iaax) )
           call conv1der_2(nbls1,nij_uniqe,l2,bl(idx1),bl(indx),
     *                     bl(ibb),bl(ibbx) )
        endif
      endif
cxxxxkl
         call conv1x_2(nbls1,mmax1,nkl_uniqe,lckl, bl(idx2),bl(indx),
     *                 bl(ixqn), bl(ixqnx)  )
                       icdnix=icdnia+(lckl-1)*mmax1
        if(where.eq.'forc' .or. where.eq.'hess') then
           call conv1der_2(nbls1,nkl_uniqe,l3,bl(idx2),bl(indx),
     *                     bl(icc),bl(iccx) )
        endif
ctime
         call txs_second(tc1xe)
         tconv1=tconv1+(tc1xe-tc1xb)
cxxxx2
c98... if(nsij.ge.nskl) then
ccccc    call conv2x(nbls1,nfumax1,npklx,lckl, bl(idx2),bl(indx),
ccccc*               bl(ihabcd),nfumax, bl(ihabcdx)  )
                     ihabcdx=ihabcd+(lckl-1)*nfumax*3
c98... else
cccc     if(nbls1.ne.nbls .or. nblsr.ne.nbls) then
cccc     call conv2x(nbls1,nfumax1,npij ,lcij, bl(idx1),bl(indx),
cccc *               bl(ihabcd),nfumax, bl(ihabcdx)  )
c98                  ihabcdx=ihabcd+(lcij-1)*nfumax*3
cccc     endif
c98... endif
ctime
         call txs_second(tc2xe)
         tconv2=tconv2+(tc2xe-tc1xe)
c
c98      call trobsa(bl,nbls1,l11,l12,mem2,immax,kmmax,lobsa)
         call trobsa_2(bl,nbls1,l11,l12,mem2)
c
ctime
        call txs_second(ttrob)
        ttrobs=ttrobs+(ttrob-tc2xe)
cxxxxxx
      if(ngcd.eq.1) then
         call assemblx(bl,first,nbls,nbls1,lnij,lnkl,
     *                 l1,l2,l3,l4,lcij,lckl,nij_uniqe,nkl_uniqe)
cnew *                 l1,l2,l3,l4,lcij,lckl,npij,npklx)
      else
c98      call gcquart(nbls,nbls1, bl(indx),
c98  *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
c98  *                bl(igcijx),ngcij, bl(igcklx),ngckl,
c98  *                bl(indgc),bl(igcoef)  )  
c
c transpose gcoef(ngcd,nbls) into gcoet(nbls,ngcd)
         call trspmo(bl(igcoef),ngcd,bl(igcoet),nbls)
c
         call assemblg(bl,first,nbls,nbls1,lnij,lnkl,ngcd,igcoet)
c ,npkl) .... spurious argument?
      endif
ctime
         call txs_second(tasse)
         tassem=tassem+(tasse-ttrob)
cxxxxx
c
         nblsr=nbls1
c
   34    continue        
   12 continue
c
c grad derivatives:
c
      nqij=nqijr
      nqkl=nqklr
c
c---------------------------------------------
c release memory :
c
      call retmem(1)  ! allocation for map_n0
c
      if(where.eq.'forc' .or. where.eq.'hess') call retmem(3)
c
c transpose but2 into buf2
      if(ngcd.gt.1) call retmem(3)
c
      if( first ) then
          call retmem(igmcnt)                         
          return
      endif
c
c release memory from trobsa
c
      call retmem(3)
c
c---------------------------------------------
c for NMR , GRADIENT and HESSIAN derivatives :
c
      lnijr=lnij
      lnklr=lnkl
c---------------------------------------------
c reserves memory for amshift
c
      call memo4b(bl, nbls,igmcnt1)
c---------------------------------------------
c for NMR derivatives :
c
c
      if(where.eq.'shif') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=1
        klderiv=1
        call iobarb(ijderiv,klderiv)
c
c construct nmr/giao derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call shift_der(nbls,lnijr,lnklr,nij_uniqe,nkl_uniqe,ngcd,
     *                idx1,idx2, ixab,ixcd,ixyab,ixycd)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c for GRADIENT derivatives :
c
      if(where.eq.'forc') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=1
        klderiv=1
        call iobarb(ijderiv,klderiv)
c
c construct gradient derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call force_der(bl,nbls,lnijr,lnklr,nij_uniqe,ngcd,idx1,ixab)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c for HESSIAN  derivatives :
c
      if(where.eq.'hess') then
        call txs_second(tderb)
c
c       return to the original values of nsij,nskl and mmax :
c
        ijderiv=2
        klderiv=2
        call iobarb(ijderiv,klderiv)
c
c construct hessian derivatives : (i+j,s|k+l,s)(x,y,z)
c
        call hessian_der(bl,nbls,lnijr,lnklr,nij_uniqe,ngcd,idx1,ixab)
c
        call txs_second(tdere)
        tderiv=tderiv+tdere-tderb
      endif
c------------------------------------
c
c*  shift angular momentum (a->b, c->d) :
ctime   
           call txs_second(tamsb)
cnew    call amshift(bl,nbls,lnij,lnkl,npij,npklx,ngcd)
        call amshift(bl,nbls,lnij,lnkl,nij_uniqe,nkl_uniqe,ngcd)
           call txs_second(tamse)
           tamshf=tamshf+(tamse-tamsb)
c-----
      call retmem(igmcnt+igmcnt1-3)                         
c-----
c
      end
c===================================================================
      subroutine erintsp_2(isupb,bl,first,nbls,acc,ikbl,npkl)
      implicit real*8 (a-h,o-z)
      logical first
c
      common /primij/ iabprim, ijdim ,ijpar1 
      common /primkl/ kabprim, kldim ,klpar1 
      common /howmany/ ntotal,noijkl,nopres,nospec,nrimtot,nrimret
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
#include "texas_lpar.fh"
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common /memor3/ nblok1d
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
      common /memor5f/ indxp
      common /memors/ nsym,ijshp,isymm
c
      common /pnl001/ ispec,ijpres2,klpres2, ijblock,klblock,iqorder
      common /pnl005/ isblsize,isblqrts,isblpoint
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
        call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                    isbl_size,isbl_point)
c-----------------------------------------------------------------
c  memory handling for a block 
c (independent of contractions)
c reserves memory for trobsa and assemble
c
      call memo4a(bl, nbls, l11,l12,mem2,igmcnt)
c-----------------------------------------------------------------
      first=.true.
c-----------------------------------------------------------------
      igcoet=1
      if(ngcd.gt.1) then
        ngcij=ngci1*ngcj1
        ngckl=ngck1*ngcl1
        call getmem(ngcij*nbls,igcijx)
        call getmem(ngckl*nbls,igcklx)
        call getmem(ngcd*nbls       ,igcoet)
      endif
c-----------------------------------------------------------------
c get pnl-size of the super-block and the pointer to its quartets:
c
      call get_nbls_pnl(isupb,bl(isblsize),bl(isblpoint),
     *                  isbl_size,isbl_point)
c
      call getmem(isbl_size,map_n0)
      call get_nbls_n0(isbl_size,isbl_point,bl(isblqrts),
     *                 bl(map_n0), nbls_n0)
c-----------------------------------------------------------------
c
c* loop over contraction :
c
      lcij=0
      do 12 l1=1,lci
      do 12 l2=1,lcj
      lcij=lcij+1
      if(ngcd.gt.1) then
        call gcpairs(nbls,lcij,ngcij,bl(igcij),bl(igcijx) )
      endif
         lckl=0
         do 34 l3=1,lck
         do 34 l4=1,lcl
         lckl=lckl+1
         if(ngcd.gt.1) then
          call gcpairs(nbls,lckl,ngckl,bl(igckl),bl(igcklx))
      endif
c
c* calculate :  rppq,rho,    rysx,rhoapb,rhocpd :
c
c2002 call precspec_2(isbl_size,isbl_point,bl(isblqrts),npkl,
c
      call precspec_2(nbls_n0,bl(map_n0),
     *                lcij,lckl, lc12,lc34,bl(idx1),bl(idx2),
     *                bl(ieab),bl(iecd),
     *                bl(iapb),bl(icpd),bl(icij),bl(ickl),
     *     bl(ixp),bl(ixq),bl(i1apb),bl(i1cpd),bl(itxab),bl(itxcd),
     *     bl(map_ij),bl(map_kl),nij_uniqe,nkl_uniqe,
     *     bl(irys),bl(iconst),bl(ixwp),bl(ixwq), nbls1,bl(indx) )
c
         nrimret=nrimret+nbls1
         nrimtot=nrimtot+nbls
         if(nbls1.eq.0) go to 34
c
c  note :
c* from here a given block (nbls) was reduced to (nbls1) ***
c*              for given l1,l2,l3,l4 !!!
c*                      and
c* rho, rppq, rhoapb, rhocpd, rys, and const
c  have dimension nbls1 (not  (nbls) )
c******************************************************
c* substitute zero for all needed bufors for quartets
c* which do not appear first time after neglect
c
         if(first) then
            nblsnot=nbls-nbls1
            if(nblsnot.ne.0) then
ckw      call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd)
         call zeroint(bl,nbls,nbls1,nblsnot,lnij,lnkl,indx,ngcd,1)
            endif
         endif
c
c********
c
c*    special code for (ss ss) and (xs ss), (sx ss), (ss xs), (ss sx)
c*          itegrals ( x= p or l )
c
      if(ngcd.eq.1) then
          call specase_2(bl,first,nbls,nbls1, bl(indx),
     *          nij_uniqe,nkl_uniqe,l1,l2,l3,l4,
     *          bl(icis),bl(icjs),bl(icks),bl(icls),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),bl(irr1) )
      else
         call gcquart(nbls,nbls1, bl(indx),
     *                ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                bl(igcijx),ngcij, bl(igcklx),ngckl,
     *                bl(indgc),bl(igcoef)  )  
c
c transpose gcoef(ngcd,nbls) into gcoet(nbls,ngcd)
         call trspmo(bl(igcoef),ngcd,bl(igcoet),nbls)
c
         call specasg(bl,first,nbls,nbls1, bl(indx),bl(idx1),bl(idx2),
     *          bl(ibuf),bl(ibuf2),
     *          bl(iconst),bl(irys),bl(ixwp),bl(ixwq),
     *          ngcd,bl(indgc),bl(igcoet),lnijkl)
      endif
c
   34    continue        
   12 continue
c
c release memory :
c
      call retmem(1)    ! allocated for map_n0
      if(ngcd.gt.1) call retmem(3)
c----
      call retmem(igmcnt)                         
c
      end
c===================================================================
c instead of the dgetmo routine of blas :
c transpse A(la,lb) into B(lb,la)
c
      subroutine trspmo(amat     ,lda,        bmat     ,ldb)
      implicit real*8 (a-h,o-z)
      dimension amat(lda,ldb), bmat(ldb,lda)
c
      do 10 i=1,lda
      do 10 j=1,ldb
      bmat(j,i)=amat(i,j)
   10 continue
c
      end
c===============================================================
      subroutine block1_pnl(isupb, isbl_size,isbl_q, ipoint,maxsize,
     *                      ijpres2,klpres2,ijklpres,nblok1,nqrt)
      dimension ijpres2(*),klpres2(*),isbl_q(*) ! dim=nquart
      dimension ijklpres(*)
      dimension nblok1(2,*), nqrt(*)
c------------------------------------------------
         ijkl=0
         do 100 iqp=1,isbl_size
             ijkl=ijkl+1
         ijklpres(iqp)=0    ! needed if sup-block was split 
         iq=isbl_q(ipoint+iqp)
         if(iq.eq.0) go to 100
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
c
             nblok1(1,ijkl)=ijcsq
             nblok1(2,ijkl)=klcsq
c
  100    continue
c
         nqrt(1)=min(isbl_size,maxsize)
c
      end
c===============================================================
      subroutine doquarts(isupb,isbl_size,isbl_q,ipoint,
     *                                    ijklpres,nblspec)
      implicit real*8 (a-h,o-z)
      dimension isbl_q(*) ! dim=nquart
      dimension ijklpres(*) ! dim.=nbls
c------------------------------------------------
c For PNL ijklpres(ijkl) is set up to zeros at the very begining
c------------------------------------------------
         ijklspec=0
!DIR$ NEXTSCALAR
         do iqp=1,isbl_size
         iq=isbl_q(ipoint+iqp)
         if(iq.ne.0) then
            ijklspec=ijklspec+1
            ijkl=iqp
check it in:
            ijklpres(ijkl)=iq
         endif
         enddo
c
         nblspec=ijklspec
c
      end
c===============================================================
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
      subroutine indexp_pnl(isupb,isbl_size,isbl_q,ipoint,map_ij,map_kl,
     *                      indxij,indxkl,ipres,indxp,nblsp)
      dimension isbl_q(*) ! dim=nquart
      dimension map_ij(isbl_size), map_kl(isbl_size)
      dimension indxij(*),indxkl(*),ipres(*),indxp(*)
c-----------------------------------------------------------------
c  setup relation between PRESENT c.s.quartets and pairs
c            ijklp and ijpar=1,..., and klpar=1,...
c-----------------------------------------------------------------
         ijkl=0
         ijklp=0
         do 100 iqp=1,isbl_size
             ijkl=ijkl+1
         iq=isbl_q(ipoint+iqp)
         if(iq.eq.0) go to 100
c
c----------old-----------------
c           ijcsq=ijpres2(iq)
c           klcsq=klpres2(iq)
c           ijpar=map_ij(ijcsq)
c           klpar=map_kl(klcsq)
c----------old-----------------
cnew : map to uniqe 
c
            ijpar=map_ij(iqp)
            klpar=map_kl(iqp)
c
c-----?      ijkl=ijkl+1
c
             if(ipres(ijkl).ne.0) then
                ijklp=ijklp+1
                indxij(ijklp)=ijpar
                indxkl(ijklp)=klpar
                indxp(ijklp) =ijkl
             endif
c
  100    continue
c
c
      nblsp=ijklp
c
c     write(82,*)' from indexp_pnl : indxIJ pairs   :'
c     write(82,88) (indxij(ii),ii=1,nblsp)    
c     write(82,*)' from indexp_pnl : indxKL pairs   :'
c     write(82,88) (indxkl(ii),ii=1,nblsp)   
c 88  format(15(i3,1x))
c
      end
c====================================================================
      subroutine get_nbls_pnl(isupb,isbl_s,isbl_p, isbl_size,isbl_poin)
      dimension isbl_s(*),isbl_p(*)
c
c returns size of the pnl-block - isbl_size
c and pointer to its quartets   - iq=isbl_q(IPOINT+iqp)
c
         isbl_size=isbl_s(isupb)
         isbl_poin=isbl_p(isupb)
c
      end
c====================================================================
      subroutine doblock2(ibl,ijblock,idoit)
      dimension ijblock(*)
c
      idoit=ijblock(ibl)
c
      end
c====================================================================
      subroutine dosupblk(isupb,isbl_s,isbl_q,ipoint,idoit,isbl_split)
      common /pnl006/ nsplit,isplit,isblsplit, isblpart
      dimension isbl_s(*),isbl_q(*),isbl_split(0:*)
ctest
c     write(6,*)'FROM DOSUPBLK :isplit=',isplit,' rep.no=',nsplit
c     write(6,*)'superblock=',isupb,' its size=',isbl_s(isupb)
ctest
c
      isbl_size=isbl_s(isupb)
      idoit=isbl_size
      if(idoit.eq.0) return
c
      icount=0
      do 10 iqp=1,isbl_size
      iq=isbl_q(ipoint+iqp)
      if(iq.gt.0) icount=icount+1
   10 continue
c
      if(icount.eq.0) then
         idoit=0
         return
      endif
c
      if(isplit.gt.0) then
         ibeg=isbl_split(nsplit-1)
         iend=isbl_split(nsplit)-1
ccc   write(6,*)'superblock=',isupb,' in range=',ibeg,iend
         if(ibeg.le.isupb .and. isupb.le.iend) then
         else
            idoit=0
         endif
      endif
c
      end
c===============================================================
      subroutine copy_block(isbl_size,isbl_q,isbl_point,isbl_copy )
      dimension isbl_q(*),isbl_copy(*)    
c
      do 10 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      isbl_copy(iqp)=iq
   10 continue
c
      end
c===============================================================
      subroutine check_in (isupb,isbl_size,isbl_q,isbl_point,isbl_copy,
     *                     maxsize)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------------
      common /time5/ tcheck,tuniq2,tmap4u
c------------------------------------------------------------
      dimension isbl_q(*),isbl_copy(*)    
c-----
      call txs_second(time1)
c-----
c
c copy back first MAXSIZE non-zero quartets from isbl_copy to isbl_q :
c
      nonzero=0
      do 10 iqp=1,isbl_size
      isbl_q(isbl_point+iqp)=0
c
      iq=isbl_copy(iqp)
      if(iq.gt.0) then
         nonzero=nonzero+1
         if(nonzero.le.maxsize) then
            isbl_q(isbl_point+iqp)=iq
         endif
      endif
   10 continue
c
c-----
      call txs_second(time2)
      tcheck=tcheck+time2-time1
c-----
c
      end
c===============================================================
      subroutine check_out(isupb,isbl_s,isbl_q,isbl_point,isbl_split,
     *                     isbl_copy,nparts,maxsize,
     *                     ibl,ijblock,igo_back )
      implicit real*8 (a-h,o-z)
c------------------------------------------------------------
      common /time5/ tcheck,tuniq2,tmap4u
c------------------------------------------------------------
      common /pnl006/ nsplit,isplit,isblsplit, isblpart
      dimension isbl_s(*),isbl_q(*), isbl_copy(*)
      dimension isbl_split(0:*)
      dimension ijblock(*)
c--------------------------------------------------------
      call txs_second(time1)
c--------------------------------------------------------
      igo_back=0
c
      isbl_size=isbl_s(isupb)
c
      IF(nparts.eq.1) THEN
         do 10 iqp=1,isbl_size
         isbl_q(isbl_point+iqp)=0
   10    continue
c
         if(isplit.gt.0) then
            ibeg=isbl_split(nsplit-1)
            iend=isbl_split(nsplit)-1
            if(ibeg.le.isupb .and. isupb.le.iend) then
               if(isupb.gt.1) ijblock(ibl-1)=0
               if(isupb.eq.iend) igo_back=1
            endif
         endif
      ENDIF
c
      IF(nparts.gt.1) THEN
c        zero out first MAXSIZE non-zero quartets in isbl_copy  :
         nonzero=0
         do 20 iqp=1,isbl_size
         iq=isbl_copy(iqp)
         if(iq.gt.0) then
            nonzero=nonzero+1
            if(nonzero.le.maxsize) then
               isbl_copy(iqp)=0
            endif
         endif
   20    continue
c
c copy from isbl_copy into isbl_q :
c
         non0_left=0
         do 30 iqp=1,isbl_size
         iq=isbl_copy(iqp)
         if(iq.gt.0) non0_left=non0_left+1
         isbl_q(isbl_point+iqp)=iq
   30    continue
c
         if(isplit.eq.0) write(6,*)' NPATRS=',nparts,' ISPLIT=0 !!'
c
         ibeg=isbl_split(nsplit-1)
         iend=isbl_split(nsplit)-1
         if(isplit.gt.0) then
            if(ibeg.le.isupb .and. isupb.le.iend) then
ckwol ?        if(isupb.gt.1) ijblock(ibl-1)=0
               igo_back=1
ckwol ?        if(isupb.eq.iend) igo_back=1
            endif
         endif
c
         if(non0_left.gt.0) then
            nsplit=nsplit-1
         else
ctry        if(isupb.lt.iend) nsplit=nsplit-1
         endif
c
      ENDIF
c
c--------------------------------------------------------
      call txs_second(time2)
      tcheck=tcheck+time2-time1
c--------------------------------------------------------
      end
c===============================================================
      subroutine get_nparts(isupb,isbl_part,nparts)
      dimension isbl_part(*)
c
      nparts=isbl_part(isupb)
c
      end
c=======================================================================
      subroutine make_map2uniq(bl,ibl,kbl,isbl_size,isbl_q,isbl_point,
     *                         ijpres2,klpres2)
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
c------------------------------------------------------------
      common /time5/ tcheck,tuniq2,tmap4u
c------------------------------------------------------------
c output :
c
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c------------------------------------------------------------
      dimension bl(*)
      dimension isbl_q(*), ijpres2(*),klpres2(*)
c--------------------------------------------------------
      call txs_second(time1)
c--------------------------------------------------------
c
c number of uniqe pairs in IBL & Kbl (nij_uniqe,nkl_uniqe) and pointers
c to the ijcs shell-pairs are obtained in prec2ij & prec2kl
c
      call memo1_int(isbl_size,map_ij)
      call memo1_int(isbl_size,map_kl)
c
      call map_2_uniq(isbl_size,isbl_q,isbl_point,ijpres2,klpres2,
     *                nij_uniqe, bl(ij_uniqe_p), 
     *                nkl_uniqe, bl(kl_uniqe_p), 
     *                bl(map_ij),bl(map_kl) )
c--------------------------------------------------------
      call txs_second(time2)
      tmap4u=tmap4u+time2-time1
c--------------------------------------------------------
c
      end
c===============================================================
      subroutine map_2_uniq(isbl_size,isbl_q,isbl_point,
     *                      ijpres2,klpres2,
     *                      nij_uniqe, ij_uniqe, 
     *                      nkl_uniqe, kl_uniqe, 
     *                      map_ij,map_kl)
      dimension isbl_q(*), ijpres2(*),klpres2(*)
      dimension ij_uniqe(nij_uniqe),kl_uniqe(nkl_uniqe)
      dimension map_ij(isbl_size),map_kl(isbl_size)
c
      do 10 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      if(iq.eq.0) go to 10
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
c
         call seartch_4_uniqe(ijcsq,ij_uniqe,nij_uniqe,ij_u)
         call seartch_4_uniqe(klcsq,kl_uniqe,nkl_uniqe,kl_u)
c
         map_ij(iqp)=ij_u
         map_kl(iqp)=kl_u
   10 continue
c
      end
c===============================================================
      subroutine        seartch_4_uniqe(ijcsq,ij_uniqe,nij_uniqe,ij_u)
      dimension ij_uniqe(nij_uniqe)
c
      if(nij_uniqe.eq.0) then
           write(6,*)' number of uniqe pair is ZERO ; wrong'
           call errquit('s_4_u:execution stopped in seartch_4_uniqe',0,
     &       INT_ERR)
      endif
      if(nij_uniqe.eq.1) then
         ij_u=1
         return
      endif
c
      ijhalf=nij_uniqe/2 
      ijdelta=ijhalf/2
c
  100 continue
      ijcsh =ij_uniqe(ijhalf)
c              write(6,*)' ijhalf=',ijhalf,' ijcsh=',ijcsh
      if(ijcsq.eq.ijcsh) then
         ij_u=ijhalf
         return
      endif
c
      if(ijcsq.gt.ijcsh) then
         ijhalf=ijhalf + ijdelta
         if(ijhalf.gt.nij_uniqe) ijhalf=nij_uniqe
         ijdelta=ijdelta/2
         if(ijdelta.eq.0) ijdelta=1
         go to 100
      endif
      if(ijcsq.lt.ijcsh) then
         ijhalf=ijhalf - ijdelta
         if(ijhalf.le. 0       ) ijhalf=1
         ijdelta=ijdelta/2
         if(ijdelta.eq.0) ijdelta=1
         go to 100
      endif
c
      end
c=======================================================================
      subroutine get_nbls_n0(isbl_size,isbl_point,isbl_q,map_n0,nbls_n0)
      dimension isbl_q(*)
      dimension map_n0(*)
c
      ijklp=0
      do 100 iqp=1,isbl_size
         iq=isbl_q(isbl_point+iqp)
         if(iq.eq.0) go to 100
         ijklp=ijklp+1
         map_n0(ijklp)=iqp
  100 continue
      nbls_n0=ijklp
c
      end
c=======================================================================
