* $Id$
c===================================================================
c kw, 1998 : the indgc(nbls) array is not used anymore 
c===================================================================
c          FOR GENERAL CONTRACTED SHELLS
c     subroutine gcparij(nbls, indxij,npij,
c     subroutine gcparkl(nbls,indxkl,npkl,
c     subroutine gcqijkl(nbls,nbls1, index,indxij,indxkl,npij,npkl,
c     subroutine gcpairs(nbls, 
c     subroutine gcquart(nbls,nbls1, index,
c     subroutine specasg(bl,first,nbls,nbls1, index,indxij,indxkl,
c     subroutine assemblg(bl,firstc,nbls,nbls1,l01,l02,ngcd,
c     subroutine asselg(firstc,xt1,lt1,lt2,nbls,indx,nbls1,
c===================================================================
c          FOR GENERAL CONTRACTED SHELLS
c
c          Used when iroute=1 (old) : tx93 
c
c NOTE :
c npij and npkl denote now UNIQE pairs , not all pairs in a block.
c===================================================================
      subroutine gcparij(nbls, indxij,npij, ii,jj,ngci1,ngcj1,ngcij,
     *                   gci,gcj, gcij)
      implicit real*8 (a-h,o-z)
c
      dimension indxij(*)
      dimension gci(npij,ngci1,*),gcj(npij,ngcj1,*)
      dimension gcij(ngcij,nbls)
c-------------------------------------------------------------------
c      This is called from Erinteg and Erintsp 
c             From contraction loops
c         FOR GENERAL CONTRACTED SHELLS
c-------------------------------------------------------------------
      do 204 ijkl=1,nbls
      ijpar=indxij(ijkl)
             ijpg=0
             do 2041 igc=1,ngci1
             coefi=gci(ijpar,igc,ii)
             do 2041 jgc=1,ngcj1
             coefj=gcj(ijpar,jgc,jj)*coefi
             ijpg=ijpg+1
             gcij(ijpg,ijkl)=coefj
 2041        continue
  204 continue
c
      end
c====================================================================
      subroutine gcparkl(nbls,indxkl,npkl, kk,ll,ngck1,ngcl1,ngckl,
     *                   gck,gcl,gckl)
      implicit real*8 (a-h,o-z)
      dimension indxkl(*)
      dimension gck(npkl,ngck1,*),gcl(npkl,ngcl1,*)
      dimension gckl(ngckl,nbls)
c-------------------------------------------------------------------
c      This is called from Erinteg and Erintsp 
c             From contraction loops (kl)
c          FOR GENERAL CONTRACTED SHELLS
c-------------------------------------------------------------------
c
      do 204 ijkl=1,nbls
      klpar=indxkl(ijkl)
             klpg=0
             do 2042 kgc=1,ngck1
             coefk=gck(klpar,kgc,kk)
             do 2042 lgc=1,ngcl1
             coefl=gcl(klpar,lgc,ll)*coefk
             klpg=klpg+1
             gckl(klpg,ijkl)=coefl
 2042        continue
c
  204 continue
c
      end
c====================================================================
      subroutine gcqijkl(nbls,nbls1, index,
     *                   ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                   indgc,gcoef,
     *                   gcij,ngcij, gckl,ngckl)
      implicit real*8 (a-h,o-z)
      dimension index(*)
c
      dimension indgc(nbls) 
      dimension gcoef(ngcd,nbls)
      dimension gcij(ngcij,nbls),gckl(ngckl,nbls)
c-------------------------------------------------------------------
c      This is called from Erinteg and Erintsp 
c             From contraction loops
c       FOR GENERAL CONTRACTED SHELLS
c-------------------------------------------------------------------
      ijpg=ngci1*ngcj1
      klpg=ngck1*ngcl1
c
      do 204 i=1,nbls1
      ijkl=index(i)
             ijklg=0
             do 2043 ijp1=1,ijpg
             gcoefij=gcij(ijp1,ijkl)
             do 2043 klp1=1,klpg
             ijklg=ijklg+1
             gcoef(ijklg,ijkl)=gcoefij*gckl(klp1,ijkl)
 2043        continue
c98       indgc(ijkl)=ijklg
  204 continue
c
      end
c===================================================================
c===================================================================
c     Used when iroute=2 (new) : tx95
c
c===================================================================
      subroutine gcpairs(nbls, lcij,ngcij, gcij,gcijx) 
      implicit real*8 (a-h,o-z)
c------------------------------------------------------
      dimension gcij(ngcij,*) ! (ngcj1,ngci1,*)
      dimension gcijx(ngcij,nbls)
c------------------------------------------------------
c     This is called from Erinteg and Erintsp 
c     From contraction loops
c------------------------------------------------------
c     
      if (ngcij .eq. 1) then
         call dfill(nbls,gcij(1,lcij),gcijx,1)
      else
         do ijkl=1,nbls
            do ijpg = 1, ngcij
               gcijx(ijpg,ijkl)=gcij(ijpg,lcij)
            enddo
         enddo
      endif
c     
      end
c====================================================================
      subroutine gcquart(nbls,nbls1, index,
     *                   ngci1,ngcj1,ngck1,ngcl1,ngcd,
     *                   gcij,ngcij,  gckl,ngckl,
ccc   output :
     *                   indgc,gcoef)
      implicit real*8 (a-h,o-z)
c------------------------------------------------------
      dimension index(*)
      dimension indgc(nbls) 
      dimension gcoef(ngcd,nbls), gcij(ngcij,nbls),gckl(ngckl,nbls)
c------------------------------------------------------
c     This is called from Erinteg and Erintsp 
c     From contraction loops
c     
c     FOR GENERAL CONTRACTED SHELLS
c------------------------------------------------------
c     
      ijpg = ngci1*ngcj1
      klpg = ngck1*ngcl1
c
      if ( klpg      .eq.1) then
c
c     Either have ijcs.ne.klcs OR both are not generally contracted 
c     (should not happen here) ... in both cases don't need to worry 
c     about off-diagonals
c
         do i=1,nbls1
            ijkl=index(i)
            do ijp1=1,ijpg
               gcoef(ijp1,ijkl)=gcij(ijp1,ijkl)*gckl(1,ijkl)
            enddo
c98             indgc(ijkl)=ijpg
         enddo
      else
         do i=1,nbls1
            ijkl=index(i)
c     
            ijklg=0
               do ijp1=1,ijpg
                  gcoefij=gcij(ijp1,ijkl)
                  do klp1=1,klpg
                     gcoef(ijklg+klp1,ijkl)=gcoefij*gckl(klp1,ijkl)
                  enddo
                  ijklg=ijklg+klpg
               enddo
c98             indgc(ijkl)=ijklg
         enddo
      endif
c     
      end
c====================================================================
      subroutine specasg(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                   buf,buf1, const,rysx,xpqr,txxr,
     *                   ngcd,indgc,gcoef,ijkln)
      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c-------------------------------------------------------------------
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension const(*),rysx(*)
cccc  dimension buf(ngcd,nbls,ijkln),buf1(nbls1,*)
      dimension buf(nbls,ijkln,ngcd),buf1(nbls1,*)
      dimension indgc(nbls)
cccc  dimension gcoef(ngcd,nbls)
      dimension gcoef(nbls,ngcd)
c***************************************************************
c**     FOR GENERAL CONTRACTED SHELLS
c**  this routine constitues the special code for
c**  two types of integrals over nbls quartets of primitive shells
c**  1. (ss|ss)
c**  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p 
c**  these integrals are also contracted here.
c**
c**  input
c**  ------
c**  all precalculated values for whole block :
c**
c**  const(nbls) - contains consts=pi3*sabcd/(pq*sqrt(ppq)) for all int.
c**  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c**  xp,xp      - geometry for p,q
c**
c**  output
c**  ------
c**  buf(ngcd,nbls,ijkln) - contains final integrals
c
c***************************************************************
c
c* memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
      do 100 i=1,nbls1
      xrys=rysx(i)
      call ft0
      bl(if00+i)=f0
      bl(if11+i)=f1
  100 continue
c
c *** special code for (ss ss) integrals
c
      if(mmax.eq.1) then
          do 2031 i=1,nbls1
          buf1(i,1)=const(i)*bl(if00+i)
 2031     continue
      go to 203
      endif
c
c *** special code for (ps ss), (sp ss) (ss ps) and (ss sp)
c
cxxxx if(mmax.eq.2) then
      if (ityp.gt.1) then
        do 101 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 101    continue
      else if (jtyp.gt.1) then
        do 102 i=1,nbls1
        xpqr(1,i)=xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=xpqr(3,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 102    continue
      else if (ktyp.gt.1) then
        do 103 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) - txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) - txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) - txxr(3,i)*bl(if00+i)
 103    continue
      else
        do 104 i=1,nbls1
        xpqr(1,i)=-xpqr(1,i)*bl(if11+i) + txxr(1,i)*bl(if00+i)
        xpqr(2,i)=-xpqr(2,i)*bl(if11+i) + txxr(2,i)*bl(if00+i)
        xpqr(3,i)=-xpqr(3,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
 104    continue
      endif
c*************
        do 106 i=1,nbls1
          buf1(i,1)=-xpqr(1,i)*const(i)
          buf1(i,2)=-xpqr(2,i)*const(i)
          buf1(i,3)=-xpqr(3,i)*const(i)
106     continue
c
c**********************************************************
c
  203 continue
c
c non-diagonal and diagonal cases are not distinguished anymore
c ngcq=ngcd (always)
c
         if(first) then
              do 304 iqu=1,ngcd
              do 304 icx=1,lnijkl
              do 304 i=1,nbls1
              xint=buf1(i,icx)
              ijkl=index(i)
                buf(ijkl,icx,iqu)=xint*gcoef(ijkl,iqu)
  304       continue
            first=.false.
         else
              do 305 iqu=1,ngcd
              do 305 icx=1,lnijkl
              do 305 i=1,nbls1
              xint=buf1(i,icx)
              ijkl=index(i)
                buf(ijkl,icx,iqu)=buf(ijkl,icx,iqu)+xint*gcoef(ijkl,iqu)
  305         continue
         endif
c
c release memory
c
      call retmem(2)
c
      end
c====================================================================
C
C     Assembling of primitive integrals for general contraction
C     
C      It uses the transposed gen.cont. coefficient matrix :
C      gcoef(ngcd,nbls) --transpose into--> gcoef(nbls,ngcd)
C....................................................................
      subroutine assemblg(bl,firstc,nbls,nbls1,l01,l02,ngcd,igcoet)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*8 where
      common /runtype/ scftype,where
c
      logical firstc
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
c
      dimension bl(*)
c--------------------------------------------------------
c for ordinary scf integrals:
c
      if(where.eq.'buff' .or. where.eq.'shif') then
           call asselg(firstc,bl(iwt0),l01,l02,nbls,bl(ibuf2),
     *                   bl(indx),nbls1, ngcd,bl(igcoet) )
      endif
c
c--------------------------------------------------------
c for gradient integral derivatives:
c
      if(where.eq.'forc') then
         ibut2=ibuf
c
         call getmem(ngcd*nbls,igc_ax)
         call getmem(ngcd*nbls,igc_bx)
         call getmem(ngcd*nbls,igc_cx)
         call coef_x_exp(nbls,nbls1,bl(indx),ngcd,
     *                   bl(igcoet),bl(iaax),bl(ibbx),bl(iccx),
     *                   bl(igc_ax),bl(igc_bx),bl(igc_cx) )
c
         call asselg_der1(firstc,bl(iwt0),l01,l02,nbls, bl(ibut2),
     *                    bl(indx),nbls1, ngcd,bl(igcoet),
     *                    bl(igc_ax),bl(igc_bx),bl(igc_cx) )
         call retmem(3)
cold
cold     call asselg_der(firstc,bl(iwt0),l01,l02,nbls,bl(ibut2),
cold *                   bl(indx),nbls1, ngcd,bl(igcoet) ,
cold *                   bl(iaax),bl(ibbx),bl(iccx))
cold
      endif
c--------------------------------------------------------
      end
c===============================================================
      subroutine asselg(firstc,xt1,lt1,lt2,nbls,buf2,
     *                    indx,nbls1,ngcd,gcoef)
      implicit real*8 (a-h,o-z)
      logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,ngcd)
      dimension gcoef(nbls,ngcd)
c-------------------------------------------------------------
      ijs=nfu(nqij)+1
      kls=nfu(nqkl)+1
c-------------------------------------------------------------
c
      IF (FIRSTC) THEN
        do 501 iqu=1,ngcd
        do 501 kl=kls,lnkl
        do 501 ij=ijs,lnij
        do 501 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
          buf2(ijkl,ij,kl,iqu)=xint*gcoef(ijkl,iqu)
  501   continue
        firstc=.false.
      ELSE
        do 601 iqu=1,ngcd
        do 601 kl=kls,lnkl
        DO 601 ij=ijs,lnij
        do 601 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
          buf2(ijkl,ij,kl,iqu)=buf2(ijkl,ij,kl,iqu)+xint*gcoef(ijkl,iqu)
  601   continue
      ENDIF
c-------------------------------------------------------------
      end
c=======================================================================
      subroutine coef_x_exp(nbls,nbls1,indx,ngcd,
     *                      gcoef,aax,bbx,ccx,
     *                      gc_ax,gc_bx,gc_cx )
      implicit real*8 (a-h,o-z)
      common /route/ iroute
      dimension indx(*)
      dimension gcoef(nbls,ngcd)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
      dimension gc_ax(nbls,ngcd),gc_bx(nbls,ngcd),gc_cx(nbls,ngcd)
c
      if(iroute.eq.1) then
         do iqu=1,ngcd
            do i=1,nbls1
               ijkl=indx(i)
               gc_ax(ijkl,iqu)=gcoef(ijkl,iqu)*aax(i)
               gc_bx(ijkl,iqu)=gcoef(ijkl,iqu)*bbx(i)
               gc_cx(ijkl,iqu)=gcoef(ijkl,iqu)*ccx(i)
            enddo
         enddo
      else
         aaxi=aax(1)
         bbxi=bbx(1)
         ccxi=ccx(1)
         do iqu=1,ngcd
            do i=1,nbls1
               ijkl=indx(i)
               gc_ax(ijkl,iqu)=gcoef(ijkl,iqu)*aaxi
               gc_bx(ijkl,iqu)=gcoef(ijkl,iqu)*bbxi
               gc_cx(ijkl,iqu)=gcoef(ijkl,iqu)*ccxi
            enddo
         enddo
      endif 
c
      end
c=======================================================================
      subroutine asselg_der1(firstc,xt1,lt1,lt2,nbls,buf2,indx,nbls1,
     *                       ngcd,gcoef,gc_ax,gc_bx,gc_cx)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
c2002 dimension buf2(4,nbls,lt1,lt2,ngcd)
      dimension buf2(nbls,lt1,lt2,ngcd,4)
      dimension gcoef(nbls,ngcd)
      dimension gc_ax(nbls,ngcd),gc_bx(nbls,ngcd),gc_cx(nbls,ngcd)
c-----------------------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------------------
c
      ijs=nfu(nqij)+1
      kls=nfu(nqkl)+1
c
      IF(FIRSTC) THEN
        do iqu=1,ngcd
           do kl=kls,lnkl
              do ij=ijs,lnij
                 do i=1,nbls1
                    ijkl=indx(i)
                    xint=xt1(i,ij,kl)
                    if(abs(xint).gt.0.0d0) then
                       buf2(ijkl,ij,kl,iqu,1)=xint*gcoef(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,2)=xint*gc_ax(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,3)=xint*gc_bx(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,4)=xint*gc_cx(ijkl,iqu)
                    else
                       buf2(ijkl,ij,kl,iqu,1)=0.0d0
                       buf2(ijkl,ij,kl,iqu,2)=0.0d0
                       buf2(ijkl,ij,kl,iqu,3)=0.0d0
                       buf2(ijkl,ij,kl,iqu,4)=0.0d0
                    endif
                 enddo
              enddo
           enddo
        enddo
        FIRSTC=.FALSE.
      ELSE
c...... region I
        do iqu=1,ngcd
           do kl=kls,nfu(nskl)
              do ij=ijs,nfu(nsij)
                 do i=1,nbls1
                    ijkl=indx(i)
                    xint=xt1(i,ij,kl)
                    if(abs(xint).gt.0.0d0) then
                       buf2(ijkl,ij,kl,iqu,1)=buf2(ijkl,ij,kl,iqu,1)
     *                                        + xint*gcoef(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,2)=buf2(ijkl,ij,kl,iqu,2)
     *                                        + xint*gc_ax(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,3)=buf2(ijkl,ij,kl,iqu,3)
     *                                        + xint*gc_bx(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,4)=buf2(ijkl,ij,kl,iqu,4)
     *                                        + xint*gc_cx(ijkl,iqu)
                    endif
                 enddo
              enddo
           enddo
        enddo
c...... region II
        do iqu=1,ngcd
           do kl=nfu(nskl)+1,nfu(nskl+1)
              do ij=ijs,nfu(nsij)
                 do i=1,nbls1
                    ijkl=indx(i)
                    xint=xt1(i,ij,kl)
                    if(abs(xint).gt.0.0d0) then
                       buf2(ijkl,ij,kl,iqu,4)=buf2(ijkl,ij,kl,iqu,4)
     *                                        + xint*gc_cx(ijkl,iqu)
                    endif
                 enddo
              enddo
           enddo
        enddo
c...... region III
        do iqu=1,ngcd
           do kl=kls,nfu(nskl)
              do ij=nfu(nsij)+1,nfu(nsij+1)
                 do i=1,nbls1
                    ijkl=indx(i)
                    xint=xt1(i,ij,kl)
                    if(abs(xint).gt.0.0d0) then
                       buf2(ijkl,ij,kl,iqu,2)=buf2(ijkl,ij,kl,iqu,2)
     *                                        + xint*gc_ax(ijkl,iqu)
                       buf2(ijkl,ij,kl,iqu,3)=buf2(ijkl,ij,kl,iqu,3)
     *                                        + xint*gc_bx(ijkl,iqu)
                    endif
                 enddo
              enddo
           enddo
        enddo
      ENDIF
c
      end
c=======================================================================
c old routine used until 1998 :
c
c  for the gradient integral derivatives
c
      subroutine asselg_der(firstc,xt1,lt1,lt2,nbls,buf2,
     *                        indx,nbls1,ngcd,gcoef,
     *                        aax,bbx,ccx)
      implicit real*8 (a-h,o-z)
      logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension gcoef(nbls,ngcd)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
C
      dimension buf2(4,nbls,lt1,lt2,ngcd)
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-------------------------------------------------------------
      ijs=nfu(nqij)+1
      kls=nfu(nqkl)+1
c-------------------------------------------------------------
c
      IF (FIRSTC) THEN
        do 501 iqu=1,ngcd
        do 501 kl=kls,lnkl
        do 501 ij=ijs,lnij
        do 501 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
          buf2(1,ijkl,ij,kl,iqu)=xint*gcoef(ijkl,iqu)
          buf2(2,ijkl,ij,kl,iqu)=xint*gcoef(ijkl,iqu)*aax(i)
          buf2(3,ijkl,ij,kl,iqu)=xint*gcoef(ijkl,iqu)*bbx(i)
          buf2(4,ijkl,ij,kl,iqu)=xint*gcoef(ijkl,iqu)*ccx(i)
  501   continue
        firstc=.false.
      ELSE
        do 601 iqu=1,ngcd
        do 601 kl=kls,lnkl
        DO 601 ij=ijs,lnij
        do 601 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
          buf2(1,ijkl,ij,kl,iqu)=buf2(1,ijkl,ij,kl,iqu)
     *                          +xint*gcoef(ijkl,iqu)
          buf2(2,ijkl,ij,kl,iqu)=buf2(2,ijkl,ij,kl,iqu)
     *                          +xint*gcoef(ijkl,iqu)*aax(i)
          buf2(3,ijkl,ij,kl,iqu)=buf2(3,ijkl,ij,kl,iqu)
     *                          +xint*gcoef(ijkl,iqu)*bbx(i)
          buf2(4,ijkl,ij,kl,iqu)=buf2(4,ijkl,ij,kl,iqu)
     *                          +xint*gcoef(ijkl,iqu)*ccx(i)
  601   continue
      ENDIF
c-------------------------------------------------------------
      end
