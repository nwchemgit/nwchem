c===============================================================
c kw 22.05.97 : nsupb(*) has been eleminated.(not used for PNNL)
c kw 22.05.97 : nasize(*)has been eleminated. 
c===============================================================
      subroutine blockin2(bl,lcore,inx,nbl2)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /cpu/ intsize,iacc,icache,memreal
c
      common /route/ iroute
c
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /memor1/ npard,mxsize,nblock1,nblock1_back
c
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
c begining of the blocking procedure for pairs
c------------------------------------------------
c Constructe blocks of contracted shells and shell-pairs
c calculates number of blocks of :
c shells (nbl1), pairs-shell (nbl2) and quartets-shell (nbl4)
c allocates memory for 
c nblock1(), npar(), mxsize(*)
c (and sets up these arrays )
c
      call blkpair(bl,ncs,inx,lcore,
     *             dbl_mb(inuc),dbl_mb(ibas),na,
     *             nblock1,nblock1_back,npard,mxsize,       ! output
     *             nbl1,nbl2,nbl4)                          ! output
c
c  output : array's addresses :
c  nblock1, nblock1_back,npard, mxsize
c------------------------------------------------
      if(nprint.gt.1) then
          call block_info(inx,bl(nblock1),bl(npard),nbl2)
      endif
c------------------------------------------------
c save addresses :
c
      call data_save
c------------------------------------------------
c This is the end of blocking procedure for PAIRS 
c The following is known at this point :
c
c  number of shell-blocks = NBL1
c  number of pair-blocks  = NBL2
c  number of quart-blocks = NBL4
c  number of pairs in the block ibl is NPAR(ibl)
c---
c The number of contracted shell Quartets is
c NBL4=NBL2*(NBL2+1)/2. The maximum size for each of 
c them is known and kept in array MAXSIZE(i=1,NSUPB) 
c (bl(mxsize)) (no of quartets). According to this max.size
c super-blocks will be split into a number of smaller blocks 
c------------------------------------------------
      end
c===============================================================
      subroutine blkpair(bl,ncs,inx,lcore, datnuc,datbas,natoms,
     *           nblock1,nblock1_back,npard,mxsize,    ! output
     *                   nbl1,nbl2,nbl4)               ! output
      implicit real*8 (a-h,o-z)
c May96 :to be used in txs_scratch_siz2 routine to estimate memory :
      common /max_pairs/ maxpar
      common /route/ iroute
c save sizes :
      common /memor1_S/ npard_S,mxsize_S,nblock1_S,nblock1_back_S
c
      dimension bl(*)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c------------------------------------------------
      call memo1_int(ncs+1,nblock1)
      call memo1_int(ncs  ,nblock1_back)
c------------------------------------------------
c make blocks of single shells :
c
      call blk_shells(ncs,inx,iroute,datnuc,datbas,
     *                bl(nblock1),nbl1,nbl2,nbl4,maxpar)
c
c output : nblock1(i-1)+1  : first shell in i-blck1   
c output : nblock1(i)      : last  shell in i-blck1   
c output : nbl1 - number of blocks of shells
c output : nbl2 - number of blocks of shell-pairs 
c output : nbl4 - number of blocks of shell-quartets
c output : maxpar-maximum number of pairs in blokcs of pairs
c------------------------------------------------
c make nblock1_back(*) array
c
      call make_nblock1_back(ncs,nbl1,bl(nblock1),bl(nblock1_back))
c
c output : nblock1_back(ics) => block of shells ics belong to
c------------------------------------------------
c allocate memory for arrays : NPAR(nbl2) and MXSIZE(NBL4)
c
      call memo1_int(nbl2,npard)
      call memo1_int(nbl4,mxsize)
c------------------------------------------------
c save sizes:
c
      npard_S        =nbl2
      mxsize_S       =nbl4
      nblock1_S      =ncs+1
      nblock1_back_S =ncs
c------------------------------------------------
c Constructe blocks of contracted shell pairs :
c
      call blk_pairs(nbl1,nbl2,bl(nblock1),bl(npard))
c
c output: NPAR(*)
c------------------------------------------------
c  Set up the vector :
c MXSIZE(super-block)-  maximum size of a super-block
c
      call blksizer(inx,ncs,nbl2,iroute,bl(nblock1),bl(npard),
     *              bl(mxsize))
c
c------------------------------------------------
      end
c===============================================================
      subroutine blksizer(inx,ncs,nbl2,iroute,nblock1,npar,
     *                    mxsize)
c-----------------------------------------------------------
c This routine determines maximum for the block-size for
c integrals calculations. The memory requimant is also 
c estimated for ordinary two-el. as well as for GIAO, GRADIENT
c and HESSIAN  derivatives.
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
      dimension inx(12,*),nblock1(0:*)
      dimension npar(*)
      dimension mxsize(*)
c
c-----------------------------------------------------
c  iforwhat shows what the Blockin2 routins are called for :
c  1) for ordinary two-el.integrals (ifor=1)
c  2) for GIAO two-el. derivatives  (ifor=2)
c  3) for gradient derivatives      (ifor=3)
c  4) for second derivatives        (ifor=4)
c
c  This is used to send info about memory to the Calcint2
c-----------------------------------------------------
c  Sets up the MXSIZE(*) array 
c
c
      nblsmax=0
      nblsmin=1 000 000   
      maxme1=0
      nbloks=0     ! counter of all blocks
      nbl12=0      ! counter of super-blocks
      do 100 ibl=1,nbl2
      ijpar=npar(ibl)
cnew  only first pair in a block:
      call get_ics_jcs(nblock1,ibl,  1  ,ics1,jcs1)
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
cnew
      call get_ics_jcs(nblock1,kbl,  1  ,kcs1,lcs1)
c
      ktype=inx(12,kcs1)
      ltype=inx(12,lcs1)
c
      ktype1=ktype
      ltype1=ltype
      if(ktype.gt.4)  ktype1=ktype-1
      if(ltype.gt.4)  ltype1=ltype-1
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
c-------------------------------------------------
      nquart_all=ijpar*klpar
      if(ibl.eq.kbl) nquart_all=ijpar*(ijpar+1)/2
c
c maximum pairs allowed by the input limit 
c for pairs/pair-block
c
      ijpar_lim=ijpar
      klpar_lim=klpar
      if(ijpar.gt.limpair) ijpar_lim=limpair
      if(klpar.gt.limpair) klpar_lim=limpair
c
c     if(ibl.eq.kbl) nquart_lim=ijpar_lim*(ijpar_lim+1)/2
c
      nquart_lim=ijpar_lim*klpar_lim
c
      nquart=min(nquart_all,nquart_lim)
c
      ijpar=ijpar_lim
      klpar=klpar_lim
c
ctest
c     write(6,*) 'blksizer:block=',nbl12,
c    *  ' nquart_all=',nquart_all,' nquart_lim=',nquart_lim
ctest
c-------------------------------------------------
c maximum block-size allowed by the cache size
c
      maxsize=       icache /nfijkl
      if(maxsize.eq.0) maxsize=1
c
c maximum block-size allowed by the input limit 
c for quartets/block
c
      if(maxsize.gt.limblks) maxsize=limblks
c
c maximum block-size according to the current task
c i.e. to the number of quartets of the same type
c
      mxsize(nbl12)=min(maxsize,nquart)
c------------------------------------------------
c calculate maximum memory needed for one block :
c 1) for ordinary two-el.integrals (ifor=1)
c 2) for GIAO two-el. derivatives  (ifor=2)
c 3) for gradient derivatives      (ifor=3)
c 4) for second derivatives        (ifor=4)
c------------------------------------------------
c
      IF( iroute.eq.1 ) THEN
        call blksize1(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                memor2,memor2ij,memor2kl,memor4,iforwhat)
      ELSE
        call blksize2(ibl,kbl,ijpar,klpar,itype1,jtype1,ktype1,ltype1,
     *                lci,lcj,lck,lcl, ngci,ngcj,ngck,ngcl,
     *                memor2,memor2ij,memor2kl,memor4,iforwhat)
      ENDIF
c     write(6,*)'block=',nbl12,'  memor2,memor2ij,memor2kl,memor4 =',
c    *                            memor2,memor2ij,memor2kl,memor4
c-------------------------------------------------
      ikbl=0       !counter of blocks belonging to one super-block
      maxmem=0
      maxqrt=0
      ikbl=ikbl+1  !always =1 for pnl
c-------------------------------------------------
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
c redefine maximum block-szie according to the memory limit :
c
      memorx=ijpar*memor2ij+klpar*memor2kl+ mxsize(nbl12)*memor4
      if(memorx.gt.limxmem) then
         mxqrts=(limxmem-ijpar*memor2ij-klpar*memor2kl)/memor4
         if(mxqrts.le.0) mxqrts=1
         if(mxqrts.gt.limblks) mxqrts=limblks
c        write(6,*)' blksizer : sup-blk=',nbl12,' full size=',
c    *   nquart_all   ,' memory needed=',memorx,' memory limit=',limxmem
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
      nbloks=nbloks+ikbl
      if(maxmem.gt.maxme1) then
         maxme1=maxmem
         ispblx=nbl12
c        remember also : 
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
c Select a super-block with the maximum & minimum block-size :
c
      if(maxqrt.gt.nblsmax) then
         nblsmax=maxqrt
         ispblq=nbl12
      endif
      if(maxqrt.lt.nblsmin) then
         nblsmin=maxqrt
         ispblm=nbl12
      endif
c--------------
  100 continue
c-----------------------------------------------------------
      IF(nprint.gt.0) THEN
c
      write(8,*)'--------------------------------------------------'
      write(8,*)' Maximum sizes of super-blocks upon restrictions :'
      write(8,*)'         cache memory < ',icache
      write(8,*)'         pairs limit  < ',limpair
      write(8,*)'         quart limit  < ',limblks
      write(8,*)'         memory limit < ',limxmem
      write(8,*)'--------------------------------------------------'
c
c
c Print out info concerning memory requiments :
c
      npairs=ncs*(ncs+1)/2
      pairs=dble(npairs)
      quart=pairs*(pairs+1)/2
      write(8,*)
      write(8,497) nbl12,nbloks,quart,nblsmax,ispblq
     *                                    ,nblsmin,ispblm
  497 format(/' Number of Super-Blocks of Contracted Quartets =',i15/
     *        ' Total number of Blocks of Contracted Quartets =',i15/
     *        ' Total number of contracted shell quartets     =',f16.0/
     *        ' maximum block-size =',i7,'  in S-block =',i6/ 
     *        ' minimum block-size =',i7,'  in S-block =',i6/)
c
      ENDIF
c-----------------------------------------------------------
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
c ifor=3   for gradient integral derivatives (first  derivatives)
c ifor=4   for hessian  integral derivatives (second derivatives)
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
      if(ifor.eq.4) where='hess'
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
      if(where.eq.'shif' .or. where.eq.'forc') then
        nsij=nsij+1
        nskl=nskl+1
      endif
      if(where.eq.'hess') then
        nsij=nsij+2
        nskl=nskl+2
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
c ifor=3   for gradient integral derivatives (first  derivatives)
c ifor=4   for hessian  integral derivatives (second derivatives)
c
      where='    '
      if(ifor.eq.2) where='shif'
      if(ifor.eq.3) where='forc'
      if(ifor.eq.4) where='hess'
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
      if(where.eq.'shif' .or. where.eq.'forc') then
        nsij=nsij+1
        nskl=nskl+1
      endif
      if(where.eq.'hess') then
        nsij=nsij+2
        nskl=nskl+2
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
      if(where.eq.'forc') then
         memor2ij=memor2ij+lci+lcj
         memor2kl=memor2kl+lck
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
ccccc if(where.ne.'shif' .or. where.ne.'forc') then
c       ----------------------------
        memasse=nblsg*(lnijkl + mem0)
c       ----------------------------
ccccc endif
      if(where.eq.'shif') then
c       ----------------------------
        memasse=nblsg*(7*lnijkl + mem0+6*nfu(nsij)*nfu(nskl)  )
c       ----------------------------
      endif
      if(where.eq.'forc') then
c       ----------------------------
change  memasse=nblsg*(9*lnijkl +4*mem0+10*nfu(nsij)*nfu(nskl)  )
        memasse=nblsg*max(9*lnijkl,4*mem0) 
        memasse=memasse + 10*nfu(nsij)*nfu(nskl)
c       ----------------------------
      endif
      if(where.eq.'hess') then
c       ----------------------------
c       belove 54 if grad. are returned to
ccc     memasse=nblsg*max(45*lnijkl,10*mem0)  ! only second returned
        memasse=nblsg*max(54*lnijkl,10*mem0)  ! first & second returned
        memasse=memasse + 55*nfu(nsij)*nfu(nskl)
c       ----------------------------
      endif
c
      if(mmax.le.2) then
         memasse=memasse+2*nbls
         return
      endif
c
        IF(LSHELLT.GT.0) THEN
           mbfkl12=lnij*nfu(nqkl+1)*nbls 
           mbfij12=nfu(nqij+1)*lnkl*nbls
c
          if(where.eq.'shif') then
           mbfkl12=lnij*nfu(nqkl1+1)*nbls + 6*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12=nfu(nqij1+1)*lnkl*nbls + 6*nfu(nqij+1)*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
           mbfkl12= 4*lnij*nfu(nqkl1+1)*nbls 
     *            +10*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12= 4*nfu(nqij1+1)*lnkl*nbls 
     *            +10*nfu(nqij+1)*nfu(nskl)*nbls
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
          mbf2l =nfu(nqij+1)*nfu(nqkl+1)*nbls 
          mbfkl3=lnij*nbls
          mbfij3=lnkl*nbls
          if(where.eq.'shif') then
            mbf2l=nfu(nqij1+1)*nfu(nqkl1+1)*nbls 
     *         +6*nfu(nqij +1)*nfu(nqkl +1)*nbls
            mbfkl3=lnij*4*nbls + 6*nfu(nsij)*nbls
            mbfij3=4*lnkl*nbls + 6*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
            mbf2l=4*nfu(nqij1+1)*nfu(nqkl1+1)*nbls 
     *         +10*nfu(nqij +1)*nfu(nqkl +1)*nbls
            mbfkl3=4*lnij*4*nbls +10*nfu(nsij)*nbls
            mbfij3=4*4*lnkl*nbls +10*nfu(nskl)*nbls
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
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=mbf3l0*nbls
          if(where.eq.'shif') then
            mbf3l1=max( nfu(nqij1+1),nfu(nqkl1+1) )
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=4*mbf3l1*nbls + 6*mbf3l0*nbls
          endif
          if(where.eq.'forc') then
            mbf3l1=max( nfu(nqij1+1),nfu(nqkl1+1) )
            mbf3l0=max( nfu(nqij +1),nfu(nqkl +1) )
            mbf3l=(4*4*mbf3l1 +10*mbf3l0)*nbls
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
          i4s =nbls
          if(where.eq.'shif') i4s =16*nbls + 6*nbls
          if(where.eq.'forc') i4s =4*16*nbls +10*nbls
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
        if(where.eq.'forc') then
            mwij=10*mwij
            mxij=10*mxij
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
        if(where.eq.'forc') then
            mvus=10*mvus
            myz =10*myz 
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
            if(where.eq.'forc') then
               mbf2l=10*mbf2l
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
c=======================================================================
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
      subroutine block_info(inx,nblock1,npar,nbl2)
      implicit real*8 (a-h,o-z)
      common /local_info/ mmax,itype,jtype,ktype,ltype,
     *                         icont,jcont,kcont,lcont
      dimension inx(12,*),nblock1(*),npar(*)
c-------------------------------------------------------------
c this is called only for print purposes (iprint>1)
c-------------------------------------------------------------
      write(8,*  )' FIRST 1000 BLOCKS '
      write(8,505)
  505 format(/
     *' Block  MMAX   type   contrac.   nfun/q      Nqrt   Nfunc',2x,
     * 'Price/i')
c
      maxst=0
      ikbl=0
      do 100 ibl=1,nbl2
      nparij=npar(ibl)
      call get_ics_jcs(nblock1,ibl,  1  ,ics1,jcs1)
c
         do 200 kbl=1,ibl
         nparkl=npar(kbl)
         ikbl=ikbl+1
c
         if(ikbl.gt.1000) RETURN
c
         call get_ics_jcs(nblock1,kbl,  1  ,kcs1,lcs1)
c
         nquart=nparij*nparkl
         if(kbl.eq.ibl) nquart=nparij*(nparij+1)/2
c
c---------------------------------------------------------------
c estimate the price of one integral in a block
c
         call get_price1(inx,ics1,jcs1,kcs1,lcs1,nfijkl,nprice)
c
c output:  nprice and nfijkl
c and mmax, types, contractions in the local common /local_info/
c---------------------------------------------------------------
c
         if(mmax.ne.maxst) then
c           for printing purpose only
            write(8,*)'   '
            maxst=mmax
         endif
c
         nfun=nquart*nfijkl
         write(8,506) ikbl,mmax,itype,jtype,ktype,ltype,
     *                icont,jcont,kcont,lcont,nfijkl,nquart,nfun,nprice
ccc  *                ngci1,ngcj1,ngck1,ngcl1,nfijkl,nquart,nfun,nprice
  200    continue
  100 continue
c
  506 format(2x,i4,2x,i2,2x,4i2,2x,4i2,3x,i5,4x,i6,1x,i8,3x,i6      )
c
c---------------------------------------------------------------
c
      end
c=======================================================================
      subroutine get_price1(inx,ics1,jcs1,kcs1,lcs1,nfijkl,nprice)
      common /local_info/ mmax,itype,jtype,ktype,ltype,
     *                         icont,jcont,kcont,lcont
      common /logic1/ ndege(1)
      common /logic2/ lenn(1)
      common /logic3/ lensm(1)
      common /logic4/ nfu(1)
      dimension inx(12,*)
c------------------------- ij part -----------
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
      itype1=itype
      jtype1=jtype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      icont=inx(5,ics1)-inx(1,ics1)
      jcont=inx(5,jcs1)-inx(1,jcs1)
      ijcont=icont*jcont
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngci1=ngci+1
      ngcj1=ngcj+1
      nfij  =lenn(itype1)*lenn(jtype1)
      ngcij =ngci1*ngcj1
c
         nqi=ndege(itype1)
         nqj=ndege(jtype1)
         nsij=nqi+nqj-1
         lnij=lensm(nsij)
         nqij=nqi
         if(nqj.gt.nqi) nqij=nqj
c
c------------------------- kl part -----------
c
         ktype=inx(12,kcs1)
         ltype=inx(12,lcs1)
         ktype1=ktype
         ltype1=ltype
         if(ktype.gt.4) ktype1=ktype-1
         if(ltype.gt.4) ltype1=ltype-1
         if(ktype1.gt.5) ktype1=ktype1-1
         if(ltype1.gt.5) ltype1=ltype1-1
         kcont=inx(5,kcs1)-inx(1,kcs1)
         lcont=inx(5,lcs1)-inx(1,lcs1)
         ngck=inx(4,kcs1)
         ngcl=inx(4,lcs1)
         ngck1=ngck+1
         ngcl1=ngcl+1
c
         nqk=ndege(ktype1)
         nql=ndege(ltype1)
         nskl=nqk+nql-1
         lnkl=lensm(nskl)
         nqkl=nqk
         if(nql.gt.nqk) nqkl=nql
c
c----------------------  ijkl part -----------
c
         nfijkl=nfij*lenn(ktype1)*lenn(ltype1)
         mmax=nsij+nskl-1
c
c---------------------------------------------
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
         nprice=nprice*icont*jcont*kcont*lcont
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
c general contraction factor :
c
        ngcfact=(ngci1+1)*(ngcj1+1)*(ngck1+1)*(ngcl1+1)
        ngcfact=ngcfact/16
c
        nprice=nprice/ngcfact
c---------------------------------------------
        nprice=nprice/nfijkl
c---------------------------------------------
      end
c=======================================================================
c================================================================
c    subroutines for blocking quartets :
c================================================================
      subroutine blockin4(isbl,ibl,kbl,nbls_pnl,bl,npar)
c---------------------------------------------
c this routine is called for each super-block
c and NBLOKS is a number of blocks belonging
c to the given super-block ISBL .
c (whatever the "super-block" means )
c---------------------------------------------
      implicit real*8 (a-h,o-z)
c
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
c
      dimension bl(*)
      dimension npar(*)      !    npar(nbl2)
c------------------------------------------------
c
c*  constructe blocks of quartets of contracted shells :
c*  set up vectors : nqrt,nibl,nkbl, nijb,nije, nklb,nkle
c
      nbloks=1            ! always for pnl
      call memo2(nbloks)  ! res.mem. for 7 arrays (above)
c
      call blockqur(isbl,ibl,kbl,maxqrt,npar,
     *              bl(nqrtd),bl(nibld),bl(nkbld),bl(nijbd),bl(nijed),
     *              bl(nklbd),bl(nkled))
c
      call memo3(nbls_pnl)
c
c end of the blocking procedure
c
      end
c===============================================================
      subroutine blockqur(isbl,ibl,kbl,maxqrt,npar,
     *                    nqrt,nibl,nkbl, nijb,nije,nklb,nkle)
c
      dimension npar(*)
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
         call nondiax(nqrt,ikbl,ibl,kbl,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      else
         call diagonx(nqrt,ikbl,ibl,kbl,npar,
     *                nibl,nkbl, nijb,nije,nklb,nkle)
      endif
c
cpnl
      nbloks=ikbl
      maxqrt=nqrt(1) 
c------------------------------------------------------------------
      end
c==================================================================
      subroutine nondiax(nqrt,ikbl,ibl,kbl,npar,
     *                   nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
c
      dimension npar(*)
      dimension nqrt(*)
c
      dimension nibl(*), nkbl(*)
      dimension nijb(*),nije(*),nklb(*),nkle(*)
c
      ijpar=npar(ibl)
      klpar=npar(kbl)
cccc  nquart=ijpar*klpar
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
      subroutine diagonx(nqrt,ikbl,ibl,kbl,npar,
     *                   nibl,nkbl, nijb,nije,nklb,nkle)
      implicit real*8 (a-h,o-z)
      dimension npar(*)
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
      common /memor1/ npard,mxsize,nblock1,nblock1_back
c
c addresses refer to local bl() as defined in texas_face :
c
      common /memor1_R/ npard_R,mxsize_R,nblock1_R,nblock1_back_R
c
      npard_R        =npard
      mxsize_R       =mxsize
      nblock1_R      =nblock1
      nblock1_back_R =nblock1_back
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
c
c     NEW routines 1997
c
c
c===============================================================
      subroutine blk_shells(ncs,inx,iroute,datnuc,datbas,
     *                      nblock1,nbl1,nbl2,nbl4,maxpar)
      implicit real*8 (a-h,o-z)
      logical txs93,txs95, cond1,cond2,cond3,cond4,cond5
      common /intlim/ limxmem,limblks,limpair
      dimension nblock1(0:ncs)
      dimension inx(12,*)
      dimension datnuc(5,*),datbas(13,*)
c
c calculate number of different shells (i.e. different blocks)
c according to two different criterions (txs93 or txs95)
c
      iexch=0
      do 55 ics0=1,ncs-1
         ist0=inx(12,ics0)                 ! type of shell
         iat0=inx(2,ics0)                  ! atom
         nzi0=0                            ! charge
         if(iat0.gt.0) nzi0=datnuc(1,iat0)
         isc0=inx(5,ics0)-inx(1,ics0)      ! contraction lenght
         isp0=inx(1,ics0)+1                ! first primitive
         exi0=datbas(1,isp0)               ! exponent of first primitive
         ngci0=inx(4,ics0)                 ! general contraction deep
c
c next shell :
c
         ics1=ics0+1
         ist1=inx(12,ics1)     
         iat1=inx(2,ics1)       
         nzi1=0                
         if(iat1.gt.0) nzi1=datnuc(1,iat1)
         isc1=inx(5,ics1)-inx(1,ics1)  
         isp1=inx(1,ics1)+1           
         exi1=datbas(1,isp1)         
         ngci1=inx(4,ics1)          
c
c---------------------------------------------
c check against criterion (3 or 5 conditions) :
c
           cond1=.false.
           cond2=.false.
           cond3=.false.
           cond4=.false.
           cond5=.false.
           if(ist0 .eq. ist1 ) cond1=.true.
           if(isc0 .eq. isc1 ) cond2=.true.
           if(ngci0.eq.ngci1 ) cond3=.true.
           if(nzi0 .eq. nzi1 ) cond4=.true.
           if(exi0 .eq. exi1 ) cond5=.true.
c
           txs93=.false.
           txs95=.false.
           if(cond1.and.cond2.and.cond3) txs93=.true.
           if(cond4.and.cond5.and.txs93) txs95=.true.
c
           if(iroute.eq.1) then
              if(.not.txs93) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
              endif
           else
              if(.not.txs95) then
                 iexch=iexch+1
                 nblock1(iexch)=ics0
              endif
           endif
c
   55 continue
c
      nblock1(0)=0
      nblock1(iexch+1)=ncs
c
c---------------------------------------------
c number of shell-blocks:
c
      nbl1=iexch+1
c---------------------------------------------
c find maximum number of pairs per block of pairs :
c
      maxshell=0
      do ibl=1,nbl1
         icsb=nblock1(ibl-1)+1
         icse=nblock1(ibl)
         ishell=icse-icsb+1
         do jbl=1,ibl
            jcsb=nblock1(jbl-1)+1
            jcse=nblock1(jbl)
            jshell=jcse-jcsb+1
            ijshell=ishell*jshell
            if(jbl.eq.ibl) ijshell=ishell*(ishell+1)/2
            if(ijshell.gt.maxshell) maxshell=ijshell
         enddo
      enddo
c
      nbl2=nbl1*(nbl1+1)/2
      nbl4=nbl2*(nbl2+1)/2
      maxpar=maxshell
c---------------------------------------------
ctest only
c     write(6,*)' nbl1=',nbl1
c     write(6,*)' nbl2=',nbl2
c     write(6,*)' nbl4=',nbl4
c---------------------------------------------
c
      end
c===============================================================
      subroutine make_nblock1_back(ncs,nbl1,nblock1,nblock1_back)
      dimension nblock1(0:ncs),nblock1_back(ncs)
c construct nblock1_back() array
      do ibl=1,nbl1
         icsb=nblock1(ibl-1)+1
         icse=nblock1(ibl)
         do ics=icsb,icse
            nblock1_back(ics)=ibl
         enddo
      enddo
      end
c===============================================================
      subroutine blk_pairs(nbl1,nbl2,nblock1, npar)
      implicit real*8 (a-h,o-z)
      dimension nblock1(0:*)                  ! nblock1(0:ncs)
      dimension npar(nbl2)                    ! output 
c
c nbl1 - number of single shell blocks 
c nbl2 - number of shell-pairs  blocks 
c    
c
c  constructe pairs of contracted shells
c
      ijbl=0             !    blocks counter (no limit)
      do 100 ibl=1,nbl1     
         ibeg=nblock1(ibl-1)+1
         iend=nblock1(ibl)
         do 200 jbl=1,ibl
            ijbl=ijbl+1
            jbeg=nblock1(jbl-1)+1
            jend=nblock1(jbl)
            ijpar=0
            do 300 ics=ibeg,iend
               jenx=jend
               if(jbl.eq.ibl) jenx=ics
               do 400 jcs=jbeg,jenx
                  ijpar=ijpar+1
                  npar(ijbl)=ijpar
  400          continue
  300       continue
  200    continue
  100 continue
c
      end
c===============================================================
