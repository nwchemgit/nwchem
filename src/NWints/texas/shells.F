* $Id: shells.f 19696 2010-10-29 16:53:42Z d3y133 $
c--------------------------------------------------------------------
      subroutine txs_shells(inx,ics1,jcs1,kcs1,lcs1)
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /lengt/ ilen,jlen,klen,llen, ilen1,jlen1,klen1,llen1
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
#include "texas_lpar.fh"      
      dimension inx(12,*)
c
c This subroutine sets up TYPE and LENGTH of shells and
c contraction information including general contraction.
c
c  type of shells :
c
      itype=inx(12,ics1)
      jtype=inx(12,jcs1)
      ktype=inx(12,kcs1)
      ltype=inx(12,lcs1)
c
      itype1=itype
      jtype1=jtype
      ktype1=ktype
      ltype1=ltype
      if(itype.gt.4) itype1=itype-1
      if(jtype.gt.4) jtype1=jtype-1
      if(ktype.gt.4) ktype1=ktype-1
      if(ltype.gt.4) ltype1=ltype-1
c
      if(itype1.gt.5) itype1=itype1-1
      if(jtype1.gt.5) jtype1=jtype1-1
      if(ktype1.gt.5) ktype1=ktype1-1
      if(ltype1.gt.5) ltype1=ltype1-1
c
c needed for transformation i.e. d6->d5, f10->f7  :
c
      ilen=inx(3,ics1)
      jlen=inx(3,jcs1)
      klen=inx(3,kcs1)
      llen=inx(3,lcs1)
c
      ilen1=len(itype1)
      jlen1=len(jtype1)
      klen1=len(ktype1)
      llen1=len(ltype1)
c
c number of general contractions
c
      ngci=inx(4,ics1)
      ngcj=inx(4,jcs1)
      ngck=inx(4,kcs1)
      ngcl=inx(4,lcs1)
c
      ngci1=ngci+1
      ngcj1=ngcj+1
      ngck1=ngck+1
      ngcl1=ngcl+1
      ngcd=ngci1*ngcj1*ngck1*ngcl1
c
c contraction begining and end :
c
      ia=inx(1,ics1)+1
      ie=inx(5,ics1)
      ja=inx(1,jcs1)+1
      je=inx(5,jcs1)
      ka=inx(1,kcs1)+1
      ke=inx(5,kcs1)
      la=inx(1,lcs1)+1
      le=inx(5,lcs1)
c
c length of contractions
c
      lci=ie-ia+1
      lcj=je-ja+1
      lck=ke-ka+1
      lcl=le-la+1
      lcij=lci*lcj
      lckl=lck*lcl
c
      return
      end
c-----------------------------------------------
      subroutine iobara(ityp,jtyp,ktyp,ltyp,where)
      character*8 where
c
c Set up the OBARAI, SHELL and LCASES commons :
c
ctest
cccc  common /derivat/ ijderiv,klderiv
ctest
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
#include "texas_lpar.fh"      
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
C
c     DIMENSION NDEGE(16),LEN(16),LENSM(15)
c     DATA NDEGE/1,2,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
c     DATA LEN  /1,3,4,6,10,15,21,28,36,45,55,66,78,91,105,120/
c     DATA LENSM/1,4,10,20,35,56,84,120,165,220,286,364,455,560,680/
C*******  UP TO: S P D F G H I J K L M N O P Q *******
C     LENSM(NSIJ)=TOTAL NUMBER OF FUNCTIONS UP TO GIVEN NSIJ
C************************************************************
      NQI=NDEGE(ITYP)
      NQJ=NDEGE(JTYP)
      NQK=NDEGE(KTYP)
      NQL=NDEGE(LTYP)
      NSIJ=NQI+NQJ-1
      NSKL=NQK+NQL-1
c--------
c for NMR derivatives
c     nsij=nsij+ijderiv
c     nskl=nskl+klderiv
c
      if (where.eq.'shif' .or. where.eq.'forc') then
        nsij=nsij+1
        nskl=nskl+1
      endif
      if (where.eq.'hess') then
        nsij=nsij+2
        nskl=nskl+2
      endif
c--------
C
      MMAX=NSIJ+NSKL-1
c     MMAX1=MMAX-1
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
cccccccccccccccccccccccccccccccc
c    handle the L-shell cases
c
       IJBEG=NFU(NQIJ)+1
       KLBEG=NFU(NQKL)+1
ccccc
       LHELP=0
       LSHELIJ=0
       LSHELKL=0
       LSHELLT=0
C
       LCAS2(1)=0
       LCAS2(2)=0
       LCAS2(3)=0
       LCAS2(4)=0
C
       LCAS3(1)=0
       LCAS3(2)=0
       LCAS3(3)=0
       LCAS3(4)=0
C
       IF(ITYP.EQ.3) THEN
         LSHELLT=LSHELLT+1
         LSHELIJ=LSHELIJ+1
         IF(JTYP.LE.3) IJBEG=1
         IF(KTYP.EQ.3) LCAS2(1)=1
         IF(LTYP.EQ.3) LCAS2(2)=1
       ENDIF
       IF(JTYP.EQ.3) THEN
         LSHELLT=LSHELLT+1
         LSHELIJ=LSHELIJ+2
         IF(ITYP.LE.3) IJBEG=1
         IF(KTYP.EQ.3) LCAS2(3)=1
         IF(LTYP.EQ.3) LCAS2(4)=1
       ENDIF
CCC
       IF(KTYP.EQ.3) THEN
         LSHELLT=LSHELLT+1
         LSHELKL=LSHELKL+1
         IF(LTYP.LE.3) KLBEG=1
         IF(LSHELIJ.EQ.3)  LCAS3(1)=1
       ENDIF
       IF(LTYP.EQ.3) THEN
         LSHELLT=LSHELLT+1
         LSHELKL=LSHELKL+2
         IF(KTYP.LE.3) KLBEG=1
         IF(LSHELIJ.EQ.3)  LCAS3(2)=1
         IF(LCAS2(1).EQ.1) LCAS3(3)=1
         IF(LCAS2(3).EQ.1) LCAS3(4)=1
       ENDIF
c
                           lcase=1        ! no l-shells
       if(lshellt.eq.1) then
         if(lshelij.eq.1)  lcase=2        ! lx/yz
         if(lshelij.eq.2)  lcase=3        ! xl/yz
         if(lshelkl.eq.1)  lcase=4        ! xy/lz
         if(lshelkl.eq.2)  lcase=5        ! xy/zl
       endif
       if(lshellt.eq.2) then        
         if(lshelij.eq.3)  lcase=6        ! ll/xy
         if(lshelkl.eq.3)  lcase=7        ! xy/ll
         if(lcas2(1).eq.1) lcase=8        ! lx/ly
         if(lcas2(2).eq.1) lcase=9        ! lx/yl
         if(lcas2(3).eq.1) lcase=10       ! xl/ly
         if(lcas2(4).eq.1) lcase=11       ! xl/yl
       endif
       if(lshellt.eq.3) then
         if(lcas3(1).eq.1) lcase=12       ! ll/lx
         if(lcas3(2).eq.1) lcase=13       ! ll/xl
         if(lcas3(3).eq.1) lcase=14       ! lx/ll
         if(lcas3(4).eq.1) lcase=15       ! xl/ll
       endif
       if(lshellt.eq.4)    lcase=16       ! ll/ll
c
      return
      end
c***************************
c-----------------------------------------------
      subroutine iobarb(ijderiv,klderiv)
c
c return original value of nsij,nskl and lnij,lnkl
c
c     common /derivat/ ijderiv,klderiv
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     *     NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
#include "texas_lpar.fh"      
c
C************************************************************
c
c---
      nsij=nsij-ijderiv
      nskl=nskl-klderiv
c---
C
      MMAX=NSIJ+NSKL-1
C
      LNIJ=LENSM(NSIJ)
      LNKL=LENSM(NSKL)
C
      NSIJ1=NSIJ+1
      NSKL1=NSKL+1
      end
