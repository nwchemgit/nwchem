c $Id: assemblx.f,v 1.1 1995-10-30 20:57:33 d3e129 Exp $
c----------------------------------------------------------------
c      ASSEMBLY OF THE 2-EL. INTEGRALS (I+J,0|K+L,0)
c----------------------------------------------------------------
      subroutine assemblx(bl,firstc,nbls,nbls1,l01,l02,
     *                    lci,lcj,lck,lcl,lcij,lckl,npij,npkl,ndiag)
      implicit real*8 (a-h,o-z)
      logical firstc,firstx
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension bl(*)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
c dimensions for assembling :
      common /dimasse/ lqij,lqkl,lqmx,lij3,lkl3,l3l,lsss
c
c----------------------------------------------------------------
c
        nqijr=nqij
        nqklr=nqkl
      if(where.eq.'shif') then
c-    - for nmr derivatives -
        nqijr=nqij1
        nqklr=nqkl1
      endif
c-----------------------------------------------------------------------
c
      firstx=firstc
c
c-----------------------------------------------------------------------
c-                  --- for buf2  ---
c
        call conbuf2(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibuf2),bl(indx))
c
      IF(lshellt.eq.0) go to 100
c
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-                     --- for bfij1 s from -> lx/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfij1),lqij,lnkl, bl(icis),npij,lci,
     *           bl(idx1),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-                        --- for bfij2 s from xl/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfij2),lqij,lnkl, bl(icjs),npij,lcj,
     *           bl(idx1),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-                          --- for bfkl1 s from xy/lz ---
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfkl1),lnij,lqkl, bl(icks),npkl,lck,
     *           bl(idx2),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-                          --- for bfkl2 s from xy/zl ---
c
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then 
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfkl2),lnij,lqkl,  bl(icls),npkl,lcl,
     *           bl(idx2),bl(indx), ijenx,klenx)
      endif
c
      IF(lshellt.eq.1) go to 100
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
c-
      if(lshelij.eq.3) then
c-                          --- for bfij3  ss from ll/xy ---
c
                ij3b=1
                kl3b=klbeg
        firstx=firstc
        call conijkl3 (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibfij3),lij3,lnkl, bl(ifij),npij,lcij,
     *                 bl(idx1),bl(indx),ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.3) then
c-                         --- for bfkl3 ss from xy/ll ---
c
                ij3b=ijbeg
                kl3b=1
        firstx=firstc
        call conijkl3 (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibfkl3),lnij,lkl3, bl(ifkl),npkl,lckl,
     *                 bl(idx2),bl(indx), ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(1).eq.1) then
c-                          --- for bf2l1  ss from lx/ly ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l1),lqij,lqkl,
     *               bl(icis),bl(icks),npij,npkl,lci,lck,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(2).eq.1) then
c-                          --- for bf2l2  ss from lx/yl ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l2),lqij,lqkl,
     *               bl(icis),bl(icls),npij,npkl,lci,lcl,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(3).eq.1) then
c-                          --- for bf2l3  ss from xl/ly ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l3),lqij,lqkl,
     *               bl(icjs),bl(icks),npij,npkl,lcj,lck,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(4).eq.1) then
c-                          --- for bf2l4  ss from xl/yl ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l4),lqij,lqkl,
     *               bl(icjs),bl(icls),npij,npkl,lcj,lcl,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c
      IF(lshellt.eq.2) go to 100
c-----------------------------------------------------------------------
c-                      --- for bf3l  ---
c
c     if(lcase.eq.12.or.lcase.eq.16) then
c-
      if(lcas3(1).eq.1) then
c-            --- for bf3l1  sss from ll/lx ---
c
        firstx=firstc
        call conb3la(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icks),bl(ifij),  bl(ibf3l1),l3l,lqmx,
     *               lck,lcij,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq.13.or.lcase.eq.16) then
c-
      if(lcas3(2).eq.1) then
c-            --- for bf3l2  sss from ll/xl ---
c
        firstx=firstc
        call conb3la(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icls),bl(ifij),  bl(ibf3l2),l3l,lqmx,
     *               lcl,lcij,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas3(3).eq.1) then
c-            --- for bf3l3  sss from lx/ll ---
c
        firstx=firstc
        call conb3lb(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icis),bl(ifkl),  bl(ibf3l3),lqmx,l3l,
     *               lci,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c
c     if(lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas3(4).eq.1) then
c-            --- for bf3l4  sss from xl/ll ---
c
        firstx=firstc
        call conb3lb(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icjs),bl(ifkl),  bl(ibf3l4),lqmx,l3l,
     *               lcj,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c
      IF(lshellt.eq.3) go to 100
c-----------------------------------------------------------------------
c
      if(lcase.eq.16) then
c-    -- for ssss(nbls)  ssss from ll/ll --
c
        firstx=firstc
        call conssss (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                bl(ifij),bl(ifkl), bl(issss),lsss ,
     *                lcij,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c
c-----------------------------------------------------------------------
c
  100 CONTINUE
c
      firstc=firstx
c
      end
c===============================================================
      subroutine conbuf2(firstc,nbls,nbls1,xt1,lt1,lt2, buf2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
cnmr
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
cnmr
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension buf2(nbls,lt1,lt2)
C
c-------
c
      IF(where.ne.'shif') THEN
c
        ijs=nfu(nqij)+1
        IF (FIRSTC) THEN
           DO 501 KL=nfu(nqkl)+1,LNKL
           DO 501 ij=ijs,LNIJ
           do 501 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  501      CONTINUE
C
           FIRSTC=.FALSE.
        ELSE
C
           DO 601 KL=nfu(nqkl)+1,LNKL
           DO 601 ij=ijs,LNIJ
           do 601 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=BUF2(ijkl,IJ,KL)+XT1(i,IJ,KL)
  601      CONTINUE
        ENDIF
c
      ELSE
c
        ijs=nfu(nqij)+1
        lnijx=nfu(nsij)
        lnklx=nfu(nskl)
        IF (FIRSTC) THEN
           DO 551 KL=nfu(nqkl)+1,nfu(nskl)
           DO 551 ij=ijs,nfu(nsij)
           do 551 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  551      CONTINUE
           DO 552 KL=nfu(nqkl)+1,nfu(nskl)
           DO 552 ij=nfu(nsij)+1,nfu(nsij+1)
           do 552 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  552      CONTINUE
           DO 553 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 553 ij=nfu(nqij)+1,nfu(nsij)
           do 553 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  553      CONTINUE
           FIRSTC=.FALSE.
c
        ELSE
C
           DO 651 KL=nfu(nqkl)+1,lnklx
           DO 651 ij=ijs,lnijx
           do 651 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  651      CONTINUE
           DO 652 KL=nfu(nqkl)+1,lnklx
           DO 652 ij=nfu(nsij)+1,nfu(nsij+1)
           do 652 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  652      CONTINUE
           DO 653 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 653 ij=nfu(nqij)+1,nfu(nsij)
           do 653 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  653      CONTINUE
        ENDIF
      ENDIF
c
      return
      end
c===============================================================
      subroutine conijkl1(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij1,lt3,lt4, facti,npij,lci,
     *                    idx1,indx, ijenx,klenx)
      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
      dimension bl(*)
C
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
      dimension facti(npij,*)
      dimension idx1(*),indx(*)
c
c********
      call convr1(bl,nbls1,ifni,facti,lci,npij,idx1,indx)
c  
      call assel1(firstc,xt1,lt1,lt2,nbls,
     *            bl(ifni),  bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
c
      call retmem(1)
c
      return
      end
c===============================================================
      subroutine assel1(firstc,xt1,lt1,lt2,nbls,
     *            facti, bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
c
c***
      implicit real*8 (a-h,o-z)
      logical firstc
C***
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
      dimension facti(*)
C
      IF (FIRSTC) THEN
C
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
              do 504 i=1,nbls1
              ijkl=indx(i)
              BFIJ1(ijkl,IJ,KL)=XT1(i,IJ,KL)*FACTI(i)
  504         CONTINUE
C
           FIRSTC=.FALSE.
      ELSE
C
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
              do 604 i=1,nbls1
              ijkl=indx(i)
              BFIJ1(ijkl,IJ,KL)=BFIJ1(ijkl,IJ,KL)+XT1(i,IJ,KL)*FACTI(i)
  604         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conijkl3(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij3,lt3,lt4,  factij,npij,lcij,
     *                    idx1,indx, ij3b,kl3b)
      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
      dimension factij(npij,*)
      dimension idx1(*),indx(*)
c********
c
      call convr1(bl,nbls1,ifnij,factij,lcij, npij,idx1,indx)
c  
      call assel2a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bfij3,lt3,lt4, bl(ifnij), indx, ij3b,kl3b)
c
      call retmem(1)
c
      return
      end
c***************
      subroutine assel2a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij3,lt3,lt4, factij, indx, ij3b,kl3b)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
      dimension factij(*)
C
      IF (FIRSTC) THEN
c
              do 502 kl=kl3b,lt4 
              do 502 ij=ij3b,lt3
              do 502 i=1,nbls1
                 ijkl=indx(i)
                 BFIJ3(ijkl,ij,kl)=XT1(i,ij,KL)*FACTIJ(i)
  502         CONTINUE
c
          FIRSTC=.FALSE.
      ELSE
C
              DO 602 KL=kl3b,lt4 
              do 602 ij=ij3b,lt3
              do 602 i=1,nbls1
              ijkl=indx(i)
              BFIJ3(ijkl,ij,kl)=bfij3(ijkl,ij,kl)+XT1(i,ij,KL)*FACTIJ(i)
  602         CONTINUE
C
      ENDIF
      return
      end
c===============================================================
      subroutine conb2ln(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bf2l1, lt3,lt4,
     *                    facti,factk,npij,npkl,
     *                    lci,lck,idx1,idx2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
      dimension facti(npij,*),factk(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c********
c
      call convr2(bl,nbls1,ifni ,facti ,lci ,npij,idx1,
     *                  ifnk ,factk ,lck ,npkl,idx2, indx)
c  
      call assel2b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bf2l1,lt3,lt4, bl(ifni),bl(ifnk),indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel2b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   bf2l1,lt3,lt4,facti,factk,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
      dimension facti(*),factk(*)
c********
C
             if(where.eq.'shif') then
c            -- for nmr derivatives -- 
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
                klenx=nfu(nqkl1+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c
      IF (FIRSTC) THEN
c
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
                 do 504 i=1,nbls1
                 ijkl=indx(i)
                 xIJ1=XT1(i,IJ,KL)*FACTI(i)
                 BF2L1(ijkl,ij,kl)=xij1*FACTK(i)
  504         CONTINUE
c*
          FIRSTC=.FALSE.
      ELSE
C
c-----------> DO 604 KL=KLBEG,NFU(NQKL+1)
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
                 do 604 i=1,nbls1
                 ijkl=indx(i)
                 xIJ1=XT1(i,IJ,KL)*FACTI(i)
                 BF2L1(ijkl,ij,kl)=BF2l1(ijkl,ij,kl)+xij1*FACTK(i)
  604         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conb3la(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     factk,factij, bf3l,lt5,lt6,
     *                     lck,lcij,npij,npkl,idx1,idx2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      dimension bl(*)
cc
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension factk(npkl,*), factij(npij,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifnk ,factk ,lck ,npkl,idx2,
     *                  ifnij,factij,lcij,npij,idx1,indx)
c
      call assel3a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bl(ifnk),bl(ifnij), bf3l,lt5,lt6,  indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel3a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   factk,factij, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension factk(*),factij(*)
c
c********
C
             if(where.eq.'shif') then
c            -- for nmr derivatives -- 
                ijenx=nfu(nqij1+1)
c--------nie->  if(nqij1.eq.nsij) ijenx=1
                klenx=nfu(nqkl1+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqij+1)
             else
                ijenx=nfu(nqij+1)
c--------nie->  if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c--
      IF (FIRSTC) THEN
c
              DO 503 KL=KLBEG,klenx
              do 503 ij=1,lt5
              do 503 i=1,nbls1
                 ijkl=indx(i)
                 XIJ3=XT1(i,ij,KL)*FACTIJ(i)
                 BF3L(ijkl,ij,KL)=xij3*FACTK(i)
  503         CONTINUE
c
          FIRSTC=.FALSE.
      ELSE
c
            DO 603 KL=KLBEG,klenx
            do 603 ij=1,lt5 
               do 603 i=1,nbls1
               ijkl=indx(i)
               XIJ3=XT1(i,ij,KL)*FACTIJ(i)
               BF3L(ijkl,ij,KL)=bf3l(ijkl,ij,kl)+xij3*FACTK(i)
  603       CONTINUE
      endif
      return
      end
ccccccc
c===============================================================
      subroutine conb3lb(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     facti,factkl, bf3l,lt5,lt6,
     *                     lci,lckl,npij,npkl,idx1,idx2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      dimension bl(*)
cc
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension facti(npij,*),factkl(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifni ,facti ,lci ,npij,idx1,
     *                  ifnkl,factkl,lckl,npkl,idx2,indx)
c
      call assel3b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bl(ifni),bl(ifnkl), bf3l,lt5,lt6, indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel3b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   facti,factkl, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2), bf3l(nbls,lt5,lt6)
      dimension facti(*), factkl(*)
c
c********
c--
             if(where.eq.'shif') then
c            -- for nmr derivatives -- 
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
                klenx=nfu(nqkl1+1)
c--------nie--> if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
c--------nie--> if(nqkl.eq.nskl) klenx=1
             endif
c--
      IF (FIRSTC) THEN
              DO 511 kl=1,lt6
              DO 511 IJ=IJBEG,ijenx
              do 511 i=1,nbls1
                 ijkl=indx(i)
                 xKL3=XT1(i,IJ,kl)*FACTKL(i)
                 BF3L(ijkl,IJ,kl)=xKL3*FACTI(i)
  511         CONTINUE
c*
          FIRSTC=.FALSE.
      ELSE
C
              DO 611 kl=1,lt6
              DO 611 IJ=IJBEG,ijenx
              do 611 i=1,nbls1
              ijkl=indx(i)
                 xKL3=XT1(i,IJ,kl)*FACTKL(i)
                 BF3L(ijkl,IJ,kl)=BF3L(ijkl,IJ,kl)+xkl3*FACTI(i)
  611         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conssss (bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    factij,factkl, ssss,isdim,
     *                    lcij,lckl,npij,npkl,idx1,idx2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
c--------------------------------------------------------------
c
      dimension factij(npij,*),factkl(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifnij,factij,lcij,npij,idx1,
     *                  ifnkl,factkl,lckl,npkl,idx2,indx)
c
      call assel4(firstc,nbls,nbls1,xt1,lt1,lt2,
     *            bl(ifnij),bl(ifnkl), ssss,isdim, indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel4(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                  factij,factkl, ssss,isdim,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
      dimension factij(*),factkl(*)
C
      IF (FIRSTC) THEN
c
              do 507 kl=1,isdim 
              do 507 ij=1,isdim 
              do 507 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)= XT1(i,ij,kl)*FACTIJ(i)*FACTKL(i)
  507         continue
c
           FIRSTC=.FALSE.
      ELSE
C
c------------------------
              do 607 kl=1,isdim 
              do 607 ij=1,isdim 
              do 607 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)=ssss(ijkl,ij,kl)+
     *                         XT1(i,ij,kl)*FACTIJ(i)*FACTKL(i)
  607         continue
c------------------------
C
      ENDIF
      return
      end
c===============================================================
C***
C***   assembling for gen. contr.
c***  It uses the Transpose buf2 array which is called BUT2 here.
C***
      subroutine assemblg(bl,firstc,nbls,nbls1,l01,l02,ngcd,ibut2)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      dimension bl(*)
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg
C******************************************************
C**    ASSEMBLY OF THE 2-EL. INTEGRALS (I+J,0|K+L,0)
C**
c********
c
        call asselg(firstc,bl(iwt0),l01,l02,nbls, bl(ibut2),
     *                bl(indx),nbls1, ngcd,bl(indgc),bl(igcoef) )
c
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine asselg(firstc,xt1,lt1,lt2,nbls,but2,indx,nbls1,
     *                  ngcd,indgc,gcoef)
      implicit real*8 (a-h,o-z)
      logical firstc
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
cc
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
cccc  dimension buf2(nbls,lt1,lt2,ngcd)
      dimension but2(ngcd,nbls,lt1,lt2)
c
      dimension indgc(nbls) 
      dimension gcoef(ngcd,nbls)
C
      ijs=nfu(nqij)+1
      IF (FIRSTC) THEN
        DO 501 KL=nfu(nqkl)+1,LNKL
        DO 501 ij=ijs,LNIJ
        do 501 i=1,nbls1
        ijkl=indx(i)
        ngcq=indgc(ijkl)
        xint=xt1(i,ij,kl)
        if(abs(xint).gt.0.d0) then
          do 502 iqu=1,ngcq
          but2(iqu,ijkl,ij,kl)=xint*gcoef(iqu,ijkl)
  502     continue
        else
          do 503 iqu=1,ngcq
          but2(iqu,ijkl,ij,kl)=0.d0
  503     continue
        endif
  501   continue
C
           FIRSTC=.FALSE.
      ELSE
C
        DO 601 KL=nfu(nqkl)+1,LNKL
        DO 601 ij=ijs,LNIJ
        do 601 i=1,nbls1
        ijkl=indx(i)
        ngcq=indgc(ijkl)
        xint=xt1(i,ij,kl)
        if(abs(xint).gt.0.d0) then
          do 602 iqu=1,ngcq
          but2(iqu,ijkl,ij,kl)=but2(iqu,ijkl,ij,kl)+xint*gcoef(iqu,ijkl)
  602     continue
        endif
  601   continue
c
      ENDIF
      end
c===============================================================
      subroutine convr1(bl,nbls,ifni,facti,lci,npij,idx1,indx)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension idx1(*),indx(*)
      dimension facti(npij,*)
c
          call getmem(nbls,ifni)
c
       ifni1=ifni-1
       do 100 i=1,nbls
       ijkl=indx(i)
       ijpar=idx1(ijkl)
       bl(ifni1+i)=facti(ijpar,lci)
  100  continue
c
      return
      end
c===============================================================
      subroutine convr2(bl,nbls,ifnk ,factk ,lck ,npkl,idx2,
     *                       ifnij,factij,lcij,npij,idx1,indx)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension idx1(*),idx2(*),indx(*)
      dimension factk(npkl,*),factij(npij,*)
c
          call getmem(nbls,ifnk)
          call getmem(nbls,ifnij)
c
       ifnk1=ifnk-1
       ifnij1=ifnij-1
c
       do 100 i=1,nbls
       ijkl=indx(i)
       ijpar=idx1(ijkl)
       klpar=idx2(ijkl)
c
       bl(ifnk1+i)=factk(klpar,lck)
       bl(ifnij1+i)=factij(ijpar,lcij)
  100  continue
c
      return
      end
c===============================================================
