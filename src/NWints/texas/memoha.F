* $Id$
c====================================================================
c kw Feb. 18,1994
c there is the new subroutine memo5 (memory handling for pairs)
c
c====================================================================
c    Memory handling subroutines for 2-electron integrals program
c
c====================================================================
      subroutine memo1_int(namount,iaddress)
      common /cpu/ intsize,iacc,icache,memreal
c
      needed=namount
      if(intsize.ne.1) needed=namount/intsize+1 
      call getmem(needed,iaddress)
c
      end
c====================================================================
      subroutine memo2(nbloks)
      common /cpu/ intsize,iacc,icache,memreal
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
c
      ndim=nbloks
      if(intsize.ne.1) ndim=ndim/intsize+1
c
      call getmem(ndim,nqrtd)     ! for nqrt array
      call getmem(ndim,nibld)     ! for nibl array
      call getmem(ndim,nkbld)     ! for nkbl array
      call getmem(ndim,nijbd)     ! for nijb array
      call getmem(ndim,nijed)     ! for nije array
      call getmem(ndim,nklbd)     ! for nklb array
      call getmem(ndim,nkled)     ! for nkle array
c
      return
      end
c====================================================================
      subroutine memo3(maxqrt)
      common /cpu/ intsize,iacc,icache,memreal
      common /memor3/ nblok1d
      common /memors/ nsym,ijshp,isymm
c
c--------------------------------------------------
      ndim=maxqrt*2
      if(intsize.ne.1) ndim=ndim/intsize+1
c
      call getmem(ndim,nblok1d)      ! for nblok1(2,*)
      call getmem(maxqrt,isymm)      ! for isymm(*)
c--------------------------------------------------
c     call memo1_int(maxqrt*2, nblok1d)  ! for nblok1(2*maxqrt) 
c     call memo1_int(maxqrt  , nsymm  )  ! for symm(maxqrt) 
c--------------------------------------------------
      end
c********
      subroutine memo4a(bl, nbls, l11,l12,mem2,igmcnt)
      double precision bl(*)
c nmr deriv
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
c--
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
#include "texas_lpar.fh"     
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
c dimensions for assembling :
      common /dimasse/ lqij,lqkl,lqmx,lij3,lkl3,l3l,lsss
c dimensions for a.m.shifting :
c     common /dimamsh/ 
c
C************************************************************
cxxx  DATA LENSM/1,4,10,20,35,56,84,120,165,220,286,364,455,560,680/
C*******  UP TO: S P D F G H I J K L M N O P Q *******
C     LENSM(NSIJ)=TOTAL NUMBER OF FUNCTIONS UP TO GIVEN NSIJ
C************************************************************
c---------------------------------------------------------------------
c  dimensions for assembling :
c  buf2(nbls,lnij,lnkl), bfij1(nbls,lqij,lnkl), bfkl1(nbls,lnij,lqkl)
c                        bfij2(nbls,lqij,lnkl), bfkl2(nbls,lnij,lqkl)
c                        bfij3(nbls,lij3,lnkl), bfkl3(nbls,lnij,lkl3)
c
c                        bf2l1(nbls,lqij,lqkl), bf2l2(nbls,lqij,lqkl)
c                        bf2l3(nbls,lqij,lqkl), bf2l4(nbls,lqij,lqkl)
c
c                        bf3l1(nbls,l3l ,lqmx), bf3l2(nbls,l3l ,lqmx)
c                        bf3l3(nbls,lqmx,l3l ), bf3l4(nbls,lqmx,l3l )
c
c                         ssss(nbls,lsss,lsss)
c---------------------------------------------------------------------
c
       lqij=nfu(nqij +1)
       lqkl=nfu(nqkl +1)
       lij3=1
       lkl3=1
       l3l =1
       lsss=1
       if(where.eq.'shif' .or. where.eq.'forc') then
          lqij=nfu(nqij1+1)
          lqkl=nfu(nqkl1+1)
          if(lshellt.gt.1) then
            lij3=4
            lkl3=4
          endif
          if(lshellt.gt.2) l3l =4
          if(lshellt.gt.3) lsss=4
       endif
       lqmx=max( lqij,lqkl )
c
c---------------------------------------------------------------------
c l11,l12,mem2 are not used for mmax.le.2 (psss)
c
       l11=1
       l12=1
       mem2=1
c---------------------------------------------------------------------
c
c* initiate all addresses :
c for trobsa :
       iwt0=1
       iwt1=1
       iwt2=1
c for assemble :
       ibuf=1
       ibuf2=1
       ibfij1=1
       ibfij2=1
       ibfkl1=1
       ibfkl2=1
       ibf2l1=1
       ibf2l2=1
       ibf2l3=1
       ibf2l4=1
       ibfij3=1
       ibfkl3=1
       ibf3l=1
c
c      ibf3l1=ibf3l
c
       ibf3l1=1
       ibf3l2=1
       ibf3l3=1
       ibf3l4=1
c
       issss=1
c
      mem0=lnij*lnkl
c
C******************************************************
c       Memory for "assemble"
c
c ------------------------------------------
c
c gen.contr.
      ngcijkl=(ngci+1)*(ngcj+1)*(ngck+1)*(ngcl+1)
      nblsg=nbls*ngcijkl
c
ccccc if(where.ne.'shif' .and. where.ne.'forc') then
      if(where.eq.'buff') then
        call getmem_zero(bl,nblsg*lnijkl,ibuf)  ! for buf(nbls,lnijkl)    ZERO
        call getmem_zero(bl,nblsg*mem0,ibuf2)  ! for buf2(nbls,lnij,lnkl) ZERO
      endif
      if(where.eq.'shif') then
c     - for nmr derivatives -
        call getmem(7*nblsg*lnijkl,ibuf)  ! for buf(nbls,lnijkl)
        ixxx=nblsg*mem0 + 6*nblsg*nfu(nsij)*nfu(nskl)
        call getmem(ixxx      ,ibuf2)  ! for buf2(nbls,lnij,lnkl)
      endif
      if(where.eq.'forc') then
c     memory allocated for ibuf will be used twice : first for
c     assembling (instead of buf2) and then for final derivatives.
c     For ibuf allocate maximum of :
        iyyy=nblsg*max(9*lnijkl,4*mem0) 
c     and for ibuf2 :
        ixxx=               10*nblsg*nfu(nsij)*nfu(nskl)
c     instead of ixxx=4*nblsg*mem0 + 10*nblsg*nfu(nsij)*nfu(nskl)
c
c 4*nblsg*mem0 is probably ALWAYS greater than 9*nblsg*lnijkl
c
c 4 comes from : ordinary contraction 
c              + rescaled contrac. with 2*expA
c              + rescaled contrac. with 2*expB
c              + rescaled contrac. with 2*expC
c 10 comes from 9 different derivatives with respect to 
c Ax,y,z , Bx,y,z and Cx,y,z (center positions)
c     plus 1 location for ordinary integrals.
c
        call getmem(iyyy  ,ibuf )  ! for buf (nbls,lnijkl)
        call getmem(ixxx  ,ibuf2)  ! for buf2(nbls,lnij,lnkl)
      endif
c
      if(where.eq.'hess') then
        iyyy=nblsg*max(54*lnijkl,10*mem0) 
        ixxx=55*nblsg*nfu(nsij)*nfu(nskl)
c
c 10 comes from : ordinary contraction 
c               + rescaled contrac. with 2*expA
c               + rescaled contrac. with 2*expB
c               + rescaled contrac. with 2*expC
c               + rescaled contrac. with 2*expA*2expB
c               + rescaled contrac. with 2*expA*2expC
c               + rescaled contrac. with 2*expB*2expC
c               + rescaled contrac. with (2*expA)**2    
c               + rescaled contrac. with (2*expB)**2
c               + rescaled contrac. with (2*expC)**2
c 54 comes from :  9 first derivatives 
c                +45 second derivatives
c
c 55 comes from :  1 ordinary integrals
c                  9 first derivatives 
c                +45 second derivatives
c
        call getmem(iyyy  ,ibuf )  ! for buf (nbls,lnijkl)
        call getmem(ixxx  ,ibuf2)  ! for buf2(nbls,lnij,lnkl)
      endif
c
c
c  count calls of getmem :
c
change  igmcnt=2     !  to save ibuf
        igmcnt=1
c
      if(mmax.le.2) return
c
        IF(LSHELLT.GT.0) THEN
c for ordinary integrals:
c
           mbfkl12=lnij*nfu(nqkl+1)*nbls 
           mbfij12=nfu(nqij+1)*lnkl*nbls
c
          if(where.eq.'shif') then
           mbfkl12=lnij*nfu(nqkl1+1)*nbls + 6*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12=nfu(nqij1+1)*lnkl*nbls + 6*nfu(nqij+1)*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
           mbfkl12=4*lnij*nfu(nqkl1+1)*nbls
     *            +10*nfu(nsij)*nfu(nqkl+1)*nbls
           mbfij12=4*nfu(nqij1+1)*lnkl*nbls
     *            +10*nfu(nqij+1)*nfu(nskl)*nbls
          endif
c
          if(lshellt.gt.1) then
            call getmem_zero(bl,mbfij12,ibfij1)  ! for bfij1 ZERO
            call getmem_zero(bl,mbfij12,ibfij2)  ! for bfij2 ZERO
            call getmem_zero(bl,mbfkl12,ibfkl1)  ! for bfkl1 ZERO
            call getmem_zero(bl,mbfkl12,ibfkl2)  ! for bfkl2 ZERO
            igmcnt=igmcnt+4
          else
            call getmem_zero(bl,mbfij12,ibfij1)  ! for bfij1 ZERO
            ibfij2=ibfij1
            call getmem_zero(bl,mbfkl12,ibfkl1)  ! for bfkl1 ZERO
            ibfkl2=ibfkl1
            igmcnt=igmcnt+2
          endif
c     
        IF( LSHELLT.GT.1 ) THEN
c
            mbf2l=nfu(nqij+1)*nfu(nqkl+1)*nbls 
            mbfkl3=lnij*nbls
            mbfij3=lnkl*nbls
c
          if(where.eq.'shif') then
            mbf2l=nfu(nqij1+1)*nfu(nqkl1+1)*nbls 
     *         +6*nfu(nqij +1)*nfu(nqkl +1)*nbls
c
            mbfkl3=lnij*4*nbls + 6*nfu(nsij)*nbls
            mbfij3=4*lnkl*nbls + 6*nfu(nskl)*nbls
          endif
          if(where.eq.'forc') then
            mbf2l=4*nfu(nqij1+1)*nfu(nqkl1+1)*nbls 
     *           +10*nfu(nqij +1)*nfu(nqkl +1)*nbls
c
            mbfkl3=4*(lnij*4*nbls) + 10*nfu(nsij)*nbls
            mbfij3=4*(4*lnkl*nbls) + 10*nfu(nskl)*nbls
          endif
c
          if(lshellt.gt.2) then
            call getmem_zero(bl,mbf2l,ibf2l1)   ! for bf2l1 ZERO
            call getmem_zero(bl,mbf2l,ibf2l2)   ! for bf2l2 ZERO
            call getmem_zero(bl,mbf2l,ibf2l3)   ! for bf2l3 ZERO
            call getmem_zero(bl,mbf2l,ibf2l4)   ! for bf2l4 ZERO
            igmcnt=igmcnt+4
          else
            call getmem_zero(bl,mbf2l,ibf2l1)   ! for bf2l1 ZERO
            ibf2l2=ibf2l1
            call getmem_zero(bl,mbf2l,ibf2l3)   ! for bf2l3 ZERO
            ibf2l4=ibf2l3
            igmcnt=igmcnt+2
          endif
c
            call getmem_zero(bl,mbfij3,ibfij3)  ! for bfij3 ZERO
            call getmem_zero(bl,mbfkl3,ibfkl3)  ! for bfkl3 ZERO
            igmcnt=igmcnt+2
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
            mbf3l=4*(4*mbf3l1*nbls) + 10*mbf3l0*nbls
          endif
c
          if(lshellt.gt.3) then
            call getmem(mbf3l,ibf3l1)  ! for bf3l1
            call getmem(mbf3l,ibf3l2)  ! for bf3l2
            call getmem(mbf3l,ibf3l3)  ! for bf3l3
            call getmem(mbf3l,ibf3l4)  ! for bf3l4
            igmcnt=igmcnt+4
           else
            call getmem(mbf3l,ibf3l1)  ! for bf3l1
            ibf3l2=ibf3l1
            call getmem(mbf3l,ibf3l3)  ! for bf3l3
            ibf3l4=ibf3l3
            igmcnt=igmcnt+2
           endif
c
        IF( LSHELLT.GT.3 ) then
c
            i4s =nbls
          if(where.eq.'shif') then
            i4s =16*nbls + 6*nbls
          endif
          if(where.eq.'forc') then
            i4s =4*16*nbls + 10*nbls
          endif
c
           call getmem_zero(bl,i4s  ,issss)  ! for ssss ZERO
c
            igmcnt=igmcnt+1
        ENDIF
        ENDIF
        ENDIF
        ENDIF
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         Memory handling for Obara-Saika-Tracy method
c     
c     0) for target classes WT0 or XT0(nbls,lnij,lnkl)
c
c     1) for recursive formulas in Obara-Saika:
c
c         WT1 or XT1( mmax, nbls, lensm(mmax) )
c
c     2) for recursive formulas in Tracy :
c        WT2(nbls,mem2)  where mem2 is a sum of all matrices  
c        from xt1(lensm(mmax),1) to  xt1(lensm(nsij),lensm(nskl))
c  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
c  for target classes
c
cc
c  for Obara-Saika
c
      l11=mmax
      l12=lensm(mmax)     
      mem1=l11*l12
cc
c  for Tracy 
c
      mem2_1=0
c98   if(nsij.ge.nskl) then
        klstep=0
        do 10 ijstep=mmax,nsij,-1
        klstep=klstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2_1=mem2_1+ijkld
   10   continue
c98   else
      mem2_2=0
        ijstep=0
        do 11 klstep=mmax,nskl,-1
        ijstep=ijstep+1
        ijdim=lensm(ijstep)
        kldim=lensm(klstep)
        ijkld=ijdim*kldim
        mem2_2=mem2_2+ijkld
   11   continue
c98   endif
c98
      mem2=max(mem2_1,mem2_2)
c
ccc   write(6,*)' memoha: mem2_1,mem2_2,mem2=',mem2_1,mem2_2,mem2
c
      call getmem_zero(bl,nbls*mem0,iwt0)   ! for wt0(nbls,lnij,lnkl) ZERO
      call getmem_zero(bl,nbls*mem1,iwt1)   ! for wt1(l11,nbls,l12) ZERO
      call getmem_zero(bl,nbls*mem2,iwt2)      ! for wt2(nbls,mem2) ZERO
c
      igmcnt=igmcnt+3
c
      return
      end
c
c********
      subroutine memo4b(bl,nbls,igmcnt)
      double precision bl(*)
c nmr deriv
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
c--
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
#include "texas_lpar.fh"
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
C
C************************************************************
c
c* initiate all addresses :
c
c for amshift :
       ix2l1=1
       ix2l2=1
       ix2l3=1
       ix2l4=1
       ix3l1=1
       ix3l2=1
       ix3l3=1
       ix3l4=1
       ixij=1
       iyij=1
       izij=1
       iwij=1
       ivij=1
       iuij=1
       isij=1
c
c------------------------------------------------
c       Memory for "shifts"
c
c* for wij and xij :
c
c---new----
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
        if(where.eq.'hess') then
            mwij=55*mwij
            mxij=55*mxij
        endif
c---new----
c
            call getmem(mwij,iwij)    ! for wij
            call getmem_zero(bl,mxij,ixij)    ! for xij ZERO
c
c  count calls of getmem :
c
            igmcnt=2
c
        IF(LSHELLT.GT.0) THEN
c
c* for vij10:
c
c--new--    mvus=lnij2
            mvus=mwvus
            myz=nfu(nqi+1)*nfu(nqj+1)*nfu(nqkl+1)
            mvus=mvus*nbls
            myz=myz*nbls
c
        if(where.eq.'shif') then
            mvus=6*mvus
            myz =6*myz 
        endif
        if(where.eq.'forc') then
            mvus=10*mvus
            myz =10*myz 
        endif
c
            call getmem(mvus,ivij)      ! for vij
            call getmem(myz ,iyij)      ! for yij
c 
           igmcnt=igmcnt+2
c     
        IF( LSHELLT.GT.1 ) THEN
            mbf2l=nfu(nqij+1)*nfu(nqkl+1) *nbls
            if(where.eq.'shif') then
               mbf2l=6*mbf2l
            endif
            if(where.eq.'forc') then
               mbf2l=10*mbf2l
            endif
c
c* for x2l1-4, uij and sij:
c
            call getmem(mvus,iuij)      ! for uij
            call getmem(mvus,isij)      ! for sij
            call getmem(myz ,izij)      ! for zij
            igmcnt=igmcnt+3
cc
          if(lshellt.gt.2) then
            call getmem(mbf2l,ix2l1)    ! for x2l1
            call getmem(mbf2l,ix2l2)    ! for x2l2
            call getmem(mbf2l,ix2l3)    ! for x2l3
            call getmem(mbf2l,ix2l4)    ! for x2l4
            igmcnt=igmcnt+4
          else
            call getmem(mbf2l,ix2l1)    ! for x2l1
            ix2l2=ix2l1                 ! for x2l2
            ix2l3=ix2l1                 ! for x2l3
            ix2l4=ix2l1                 ! for x2l4
            igmcnt=igmcnt+1
          endif
c
        IF( LSHELLT.GT.2 ) THEN
c
         mnbls=nbls
         if(where.eq.'shif') mnbls=6*nbls
         if(where.eq.'forc') mnbls=10*nbls
c
         if(lshellt.gt.3) then
            call getmem(mnbls*nfu(nqkl+1), ix3l1) ! for x3l1
            call getmem(mnbls*nfu(nqkl+1), ix3l2) ! for x3l2
            call getmem(mnbls*nfu(nqij+1), ix3l3) ! for x3l3
            call getmem(mnbls*nfu(nqij+1), ix3l4) ! for x3l4
            igmcnt=igmcnt+4
          else
            call getmem(mnbls*nfu(nqkl+1), ix3l1) ! for x3l1
            ix3l2=ix3l1
            call getmem(mnbls*nfu(nqij+1), ix3l3) ! for x3l3
            ix3l4=ix3l3
            igmcnt=igmcnt+2
          endif
c
        ENDIF
        ENDIF
        ENDIF
c
      return
      end
c
c================================================================
      subroutine memo5a_2(npij,mmax1)
c------------------------------------------
c Memory handling for left-hand pairs:
c
c 1: for individual shells (2 quantities)
c   cis,cjs - contr coef. dimensions are (lci), (lcj)
c
c 2: for : xab(ijpar,3) and xp, xpn, xpp all (ijpar,3,lcij)
c
c 3: for : apb, rapb, factij, (lcij)
c          ceofij and sij all (ijpar,lcij)
c
c 4. for : txab(ijpar,3,lcij)
c
c Total number of calls of Getmem is 11 or 12 (if gen.con.)
c OR 13 or 14 if where='forc'
c------------------------------------------
c for gradient derivatives:
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
c
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c------------------------------------------
      ijpar=npij
c------------------------------------------
c reserve memory for left-hand pairs IJ :
c
       ndi=   ijpar*lci
       ndj=   ijpar*lcj
c
      call getmem(lci,icis)       ! for cis(lci)                 1
      call getmem(lcj,icjs)       ! for cjs(lcj)                 2
      call getmem(ijpar*3,ixab)   ! for xab(ijpar,3)             3
c
       ndij=ndi*lcj
       ndij3=ndij*3
c
ckw Do not change this order
      call getmem(ndij3,ixp)     ! for xp(ijpar,3,lcij)          4
      call getmem(ndij3,ixpn)    ! for xpn(ijpar,3,lcij)         5
      call getmem(ndij3,ixpp)    ! for xpp(ijpar,3,lcij)         6
ckw up to here.
c
      call getmem(lcij,ifij)     ! for factij(lcij)              7 
      call getmem(ndij,icij)     ! for coefij(ijpar,lcij)        8
      call getmem(ndij,ieab)     ! for eab(ijpar,lcij)           9
      call getmem(ndij3,itxab)   ! for txab(ijpar,3,lcij)       10
c
      ndijm=lcij*mmax1
      call getmem(ndijm,iabnia)  ! for abnia(mmax-1,lcij)       11
c
c------------------------------------------
c for general contraction on IJ-pairs
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c-----
c
      igcij=1
      if(ngcd.gt.1) then
        ndijg=lcij*ngci1*ngcj1
        call getmem(ndijg,igcij)              !               12
      endif
c
      iaa=1
      ibb=1
      if(where.eq.'forc' .or. where.eq.'hess') then
         call getmem(ndi,iaa)     ! for  aa(ijpar,lci)        13  
         call getmem(ndj,ibb)     ! for  bb(ijpar,lcj)        14   
      endif
c------------------------------------------
      end
c================================================================
      subroutine memo5b_2(npkl,mmax1)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory handling for right-hand pairs:
c------------------------------------------
c for gradient derivatives:
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c------------------------------------------
      klpar=npkl
c------------------------------------------
c reserve memory for right-hand pairs KL :
c
       ndk=   klpar*lck
cccc   ndl=   klpar*lcl
c
      call getmem(lck,icks)       ! for cks(lck)                1
      call getmem(lcl,icls)       ! for cls(lcl)                2
      call getmem(klpar*3,ixcd)  ! for xcd(klpar,3)             3
c
       ndkl=ndk*lcl
       ndkl3=ndkl*3
c
ckw Do not change this order
      call getmem(ndkl3,ixq)     ! for xq(klpar,3,lckl)         4
      call getmem(ndkl3,ixqn)    ! for xqn(klpar,3,lckl)        5
      call getmem(ndkl3,ixqq)    ! for xqq(klpar,3,lckl)        6
ckw up to here.
c
      call getmem(ndkl,ifkl)     ! for factkl(klapr,lckl)       7
      call getmem(ndkl,ickl)     ! for coefkl(klapr,lckl)       8
      call getmem(ndkl,iecd)     ! for ecd(klapr,lckl)          9
      call getmem(ndkl3,itxcd)   ! for txcd(klpar,3,lckl)      10
c
      ndklm=lckl*mmax1
      call getmem(ndklm,icdnia)  ! for cdnia(mmax-1,lckl)      11
c------------------------------------------
c for general contraction on KL-pairs
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c-----
      igckl=1
      if(ngcd.gt.1) then
        ndklg=lckl*ngck1*ngcl1
        call getmem(ndklg,igckl)      !               12
      endif
c------------------------------------------
      icc=1
      if(where.eq.'forc' .or. where.eq.'hess') then
         call getmem(ndk,icc)   ! for  cc(klpar,lck) 13
      endif
c------------------------------------------
      end
c================================================================
      subroutine memo5c_2(nbls,mmax1,npij,npkl,nfumax)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory handling 
c
c 3: and quartets precalculations (12 quantities)
c (for whole block of contracted quartets and 
c        one primitive quartet )
c
c Total number of calls of Getmem is 21 or 23 (if gen.cont)
c------------------------------------------
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
c------------------------------------------
c reserve memory for quartets ijkl
c------------------------------------------
      nblsi=nbls
      if(intsize.ne.1) nblsi=nbls/intsize+1
c------------------------------------------
      call getmem(nblsi,indxp)   !                    1
      call getmem(nblsi,idx1)    ! for indxij         2
      call getmem(nblsi,idx2)    ! for indxkj         3
      call getmem(nblsi,indx)    ! for index          4
c
      call getmem(1   ,irppq)    ! for rppq(1   )     5
      call getmem(nbls,irr1)     ! for rr1(nbls)      6  
c
      call getmem(1   ,irhoapb)  ! for rhoapb(1   )   7
      call getmem(1   ,irhocpd)  ! for rhocpd(1   )   8
c
      nbls3=nbls*3
      call getmem(nbls3,ixpnx)   !                    9
      call getmem(nbls3,ixwp)    ! for xwp(nbls,3)   10
      call getmem(nbls3,ixqnx)   !                   11
      call getmem(nbls3,ixwq)    ! for xwq(nbls,3)   12
      call getmem(nbls3,ip1234)  ! for p1234(nbls,3) 13
      call getmem(1   ,iabcd)    ! for abcd(1   )    14
      call getmem(nbls,iconst)   ! for const(nbls)   15
      call getmem(nbls,irys)     ! for rys(nbls)     16
c
      nfha=3*nfumax*max(lcij,lckl)
      call getmem(nfha,ihabcd)    !                  17
c------------------------------------------
c for general contraction
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c------------------------------------------
c for both gen.contr. and segmented basis sets
c because of the common Destiny
c
      call getmem(ngcd,icfg)        !               18
      call getmem(ngcd,jcfg)        !               19
      call getmem(ngcd,kcfg)        !               20
      call getmem(ngcd,lcfg)        !               21
c
c------------------------------------------
c for general contraction
c
      indgc=1
      igcoef=1
c
      if(ngcd.gt.1) then
        call getmem(nbls,indgc)       !             22
        call getmem(nbls*ngcd,igcoef) !             23
      endif
c
c------------------------------------------
      end
c====================================================================
      subroutine memo6(npij,npkl)
      common /memor6/ ixyab,ixycd
c**************
c
c Memory handling for NMR derivatives
c reserve memory for pair quantities :
c
c  ( Xa*Yb - Ya*Xb ) = xyab(ijpar,3)  - contributes to Z deriv.
c  (-Xa*Zb + Za*Xb ) = xyab(ijpar,2)  - contributes to Y deriv.
c  ( Ya*Zb + Za*Yb ) = xyab(ijpar,1)  - contributes to X deriv.
c
c  ( Xc*Yd - Yc*Xd ) = xycd(klpar,3)  - contributes to Z deriv.
c  (-Xc*Zd + Zc*Xd ) = xycd(klpar,2)  - contributes to Y deriv.
c  ( Yc*Zd + Zc*Yd ) = xycd(klpar,1)  - contributes to X deriv.
c
c**************
c
      npij3=3*npij
      npkl3=3*npkl
c
      call getmem(npij3,ixyab)
      call getmem(npkl3,ixycd)
c
      end
c================================================================
c used when iroute=1 (old) :
c
      subroutine memo5a_1(npij,mmax1)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory handling for left-hand pairs:
c
c 1: for individual shells (4 quantities)
c ( aa, bb  - exponents ) and  ( cis,cjs - contr coef.)
c  dimensions are (ijpar,lcij)
c
c 2: for : xab(ijpar,3) and xp, xpn, xpp all (ijpar,3,lcij)
c
c 3: for : apb, rapb, factij, ceofij and sij all (ijpar,lcij)
c
c 4. for : txab(ijpar,3,lcij)
c
c Total number of calls of Getmem is 13 or 15 (if gen.con.)
c------------------------------------------
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
c------------------------------------------
      ijpar=npij
c------------------------------------------
c reserve memory for left-hand pairs IJ :
c
       ndi=   ijpar*lci
       ndj=   ijpar*lcj
c
      call getmem(ndi,iaa)        ! for  aa(ijpar,lci)           1
      call getmem(ndj,ibb)        ! for  bb(ijpar,lcj)           2
      call getmem(ndi,icis)       ! for cis(ijpar,lci)           3
      call getmem(ndj,icjs)       ! for cjs(ijpar,lcj)           4
      call getmem(ijpar*3,ixab)   ! for xab(ijpar,3)              5
c
       ndij=ndi*lcj
       ndij3=ndij*3
c
ckw Do not change this order
      call getmem(ndij3,ixp)     ! for xp(ijpar,3,lcij)          6
      call getmem(ndij3,ixpn)    ! for xpn(ijpar,3,lcij)         7
      call getmem(ndij3,ixpp)    ! for xpp(ijpar,3,lcij)         8
ckw up to here.
c
c     call getmem(ndij,iapb)     ! for apb(ijpar,lcij)             
c     call getmem(ndij,i1apb)    ! for rapb(ijpar,lcij)         
      call getmem(ndij,ifij)     ! for factij(ijpar,lcij)        9 
      call getmem(ndij,icij)     ! for coefij(ijpar,lcij)       10
      call getmem(ndij,ieab)     ! for eab(ijpar,lcij)          
c
      call getmem(ndij3,itxab)   ! for txab(ijpar,3,lcij)       11
c
      ndijm=ndij*mmax1
      call getmem(ndijm,iabnia)  ! for abnia(ijpar,mmax-1,lcij) 12
c
c------------------------------------------
c for general contraction on IJ-pairs
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c-----
c
      igci=1
      igcj=1
c
      if(ngcd.gt.1) then
        ndig=ndi*ngci1
        ndjg=ndj*ngcj1
        call getmem(ndig,igci)        !               13
        call getmem(ndjg,igcj)        !               14
      endif
c
c------------------------------------------
      end
c================================================================
      subroutine memo5b_1(npkl,mmax1)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory handling for right-hand pairs:
c
c 1: for individual shells (4 quantities)
c ( cc, dd  - exponents ) and  ( cks,cls - contr coef.)
c  dimensions are (klpar,lcij)
c
c 2: for : xcd(ijpar,3) and xq, xqn, xqq all (klpar,3,lckl)
c
c 3: for : cpd, rcpd, factkl, coefkl and skl all (klpar,lckl)
c
c 4. for : txcd(klpar,3,lckl)
c
c Total number of calls of Getmem is 13 or 15 (if gen.con.)
c------------------------------------------
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
c------------------------------------------
      klpar=npkl
c------------------------------------------
c reserve memory for right-hand pairs KL :
c
       ndk=   klpar*lck
       ndl=   klpar*lcl
c
      call getmem(ndk,icc)        ! for  cc(klpar,lck)           1
      call getmem(ndl,idd)        ! for  dd(klpar,lcl)           2
      call getmem(ndk,icks)       ! for cks(klpar,lck)           3
      call getmem(ndl,icls)       ! for cls(klpar,lcl)           4
      call getmem(klpar*3,ixcd)  ! for xcd(klpar,3)              5
c
       ndkl=ndk*lcl
       ndkl3=ndkl*3
c
ckw Do not change this order
      call getmem(ndkl3,ixq)     ! for xq(klpar,3,lckl)          6
      call getmem(ndkl3,ixqn)    ! for xqn(klpar,3,lckl)         7
      call getmem(ndkl3,ixqq)    ! for xqq(klpar,3,lckl)         8
ckw up to here.
c
c     call getmem(ndkl,icpd)     ! for cpd(klapr,lckl)     
c     call getmem(ndkl,i1cpd)    ! for rcpd(klapr,lckl)   
      call getmem(ndkl,ifkl)     ! for factkl(klapr,lckl)        9
      call getmem(ndkl,ickl)     ! for coefkl(klapr,lckl)       10
      call getmem(ndkl,iecd)     ! for ecd(klapr,lckl)  
c
      call getmem(ndkl3,itxcd)   ! for txcd(klpar,3,lckl)       11
c
      ndklm=ndkl*mmax1
      call getmem(ndklm,icdnia)  ! for cdnia(klpar,mmax-1,lckl) 12
c
c------------------------------------------
c for general contraction on KL-pairs
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c-----
c
      igck=1
      igcl=1
c
      if(ngcd.gt.1) then
        ndkg=ndk*ngck1
        ndlg=ndl*ngcl1
        call getmem(ndkg,igck)        !               13
        call getmem(ndlg,igcl)        !               14
      endif
c------------------------------------------
      end
c================================================================
      subroutine memo5c_1(bl,nbls,mmax1,npij,npkl,nfha,nfumax)
      double precision bl(*)
      common /cpu/ intsize,iacc,icache,memreal
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
c------------------------------------------
c Memory handling 
c
c 3: and quartets precalculations (12 quantities)
c (for whole block of contracted quartets and 
c        one primitive quartet )
c
c Total number of calls of Getmem is 24 or 26 (if gen.cont)
c------------------------------------------
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
c
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
      common /memor5f/ indxp
c------------------------------------------
c reserve memory for quartets ijkl
c------------------------------------------
      nblsi=nbls
      if(intsize.ne.1) nblsi=nbls/intsize+1
c------------------------------------------
      call getmem(nblsi,indxp)   !                    3
c------------------------------------------
c
      call getmem(nblsi,idx1)    ! for indxij         4
      call getmem(nblsi,idx2)    ! for indxkj         5
      call getmem(nblsi,indx)    ! for index          6
c
      call getmem(nbls,irppq)    ! for rppq(nbls)     7    
cNOT  call getmem(nbls,irho)     ! for rho(nbls)      8   
      call getmem(nbls,irr1)     ! for rr1(nbls)      9  
c       
c
      call getmem(nbls,irhoapb)  ! for rhoapb(nbls)   10
      call getmem(nbls,irhocpd)  ! for rhocpd(nbls)   11
c
      nbmx=nbls*mmax1
      call getmem(nbmx,iabnix)   !                    12
      call getmem(nbmx,icdnix)   !                    13
c
      nbls3=nbls*3
      call getmem(nbls3,ixpnx)   !                    14
      call getmem(nbls3,ixwp)    ! for xwp(nbls,3)    15
      call getmem(nbls3,ixqnx)   !                    16
      call getmem(nbls3,ixwq)    ! for xwq(nbls,3)    17
      call getmem(nbls3,ip1234)  ! for p1234(nbls,3)  18
      call getmem(nbls,iabcd)    ! for abcd(nbls)     19
      call getmem(nbls,iconst)   ! for const(nbls)    20
      call getmem(nbls,irys)     ! for rys(nbls)      21
c
      call getmem(nfha*3,ihabcd) !                    22
      call getmem_zero(bl,nbls3*nfumax,ihabcdx)  !            23 ZERO
c
c------------------------------------------
c for general contraction
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c------------------------------------------
c for both gen.contr. and segmented basis sets
c because of the common Destiny
c
      call getmem(ngcd,icfg)          !               24
      call getmem(ngcd,jcfg)          !               25
      call getmem(ngcd,kcfg)          !               26
      call getmem(ngcd,lcfg)          !               27
c
c------------------------------------------
c for general contraction
c
      indgc=1
      igcoef=1
c
      if(ngcd.gt.1) then
        call getmem(nbls,indgc)       !               32
        call getmem(nbls*ngcd,igcoef) !               33
      endif
c
c------------------------------------------
      end
c====================================================================
