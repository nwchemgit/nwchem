      subroutine prec2ij(ibl, bl,inx,npar,nbl2, iis,jjs,ijbl)
c $Id: precalc2.F,v 1.8 1996-02-12 19:44:35 d3g681 Exp $
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /route/ iroute
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primij/ iabprim, ijdim ,ijpar1 
      common /time0/ tprec2
c
      common /map4uniq/ map_ij,ijcs_b,ijcs_e, map_kl,klcs_b,klcs_e
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
c     if(icheck.gt.0) then
c        call getmem(0,l0)
c        return
c     endif
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
      nparij=npar(ibl)
c---------------------------------------
c get : nij_uniqe=no of unique pairs and pointer to the list
c
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
c---------------------------------------
      call dimenij(ibl,inx,nparij,nbl2, iis,jjs,ijbl,ijdim,ijcont)
c---------------------------------------
      ijpar1=nparij
c---------------------------------------
      IF( iroute.eq.1 ) THEN
         call getmem(3*ijdim,iabprim)
         call ab_prim_1(ibl,nparij,ijbl,nbl2,
     *                  dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                  bl(iabprim),ijcont,
     *                  nij_uniqe,bl(ij_uniqe_p), 
     *                  bl(map_ij),ijcs_b,ijcs_e )
      ELSE
         call getmem(2*ijcont+ijdim,iabprim)
         iapb =iabprim
         i1apb=iabprim+ijcont
         isab =iabprim+ijcont*2
         call ab_prim_2(ibl,nparij,ijbl,nbl2,
     *                  dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                  bl(iapb),bl(i1apb),bl(isab),ijcont,
     *                  nij_uniqe,bl(ij_uniqe_p), 
     *                  bl(map_ij),ijcs_b,ijcs_e )
      ENDIF
c---------------------------------------
      call txs_second(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
      subroutine prec2kl(ibl, bl,inx,npar,nbl2, iis,jjs,ijbl)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /route/ iroute
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primkl/ kabprim, kldim ,klpar1 
      common /time0/ tprec2
c
      common /map4uniq/ map_ij,ijcs_b,ijcs_e, map_kl,klcs_b,klcs_e
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
c     if(icheck.gt.0) then
c        call getmem(0,l0)
c        return
c     endif
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
      nparkl=npar(ibl)
c---------------------------------------
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
c---------------------------------------
      call dimenij(ibl,inx,nparkl,nbl2, iis,jjs,ijbl,kldim,klcont)
c---------------------------------------
      klpar1=nparkl
c---------------------------------------
      IF( iroute.eq.1 ) THEN
         call getmem(3*kldim,kabprim)
         call ab_prim_1(ibl,nparkl,ijbl,nbl2,
     *                  dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                  bl(kabprim),klcont,
     *                  nij_uniqe,bl(ij_uniqe_p), 
     *                  bl(map_kl),klcs_b,klcs_e )
      ELSE
         call getmem(2*klcont+kldim,kabprim)
         icpd =kabprim
         i1cpd=kabprim+klcont
         iscd =kabprim+klcont*2
         call ab_prim_2(ibl,nparkl,ijbl,nbl2,
     *                  dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                  bl(icpd),bl(i1cpd),bl(iscd),klcont,
     *                  nij_uniqe,bl(ij_uniqe_p), 
     *                  bl(map_kl),klcs_b,klcs_e )
      ENDIF
c---------------------------------------
      call txs_second(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c====================================================================
      subroutine dimenij(ibl,inx,nparij,nbl2, iis,jjs,ijbl,ijdim,ijcont)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
c
      ijcs1=ijbl(ibl,1)
      ics1=iis(ijcs1)
      jcs1=jjs(ijcs1)
      icont=inx(5,ics1)-inx(1,ics1)
      jcont=inx(5,jcs1)-inx(1,jcs1)
      ijcont=icont*jcont
      ijdim=nparij*ijcont
      end
c====================================================================
      subroutine precdiag
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /route/ iroute
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c---------------------------------------------------------------
c* for a diagonal case : npklx=0
c
c* since this is for a diagonal block
c* pairs IJ and KL are the same
c
      ixcd=ixab
      ixq =ixp
      ixqn=ixpn
      ixqq=ixpp
c   
      icpd=iapb
      i1cpd=i1apb
      ickl=icij
      ifkl=ifij
      iscd=isab
c
      icc=iaa
      idd=ibb
c
      icks=icis
      icls=icjs
c
      iecd=ieab
c
      itxcd=itxab
c
      icdnia=iabnia
c-----------------------------------
      IF( iroute.eq.1 ) THEN
        igck=igci
        igcl=igcj
      ELSE
        igckl=igcij
      ENDIF
c-----------------------------------
      end
c====================================================================
      subroutine specasg(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                   buf,buf1, const,rysx,xpqr,txxr,
     *                   ngcd,indgc,gcoef,ijkln)
c****
c Works with Transpose BUF array 
c****
      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
c******************************************************
c
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension const(*),rysx(*)
      dimension buf(ngcd,nbls,ijkln),buf1(nbls1,*)
      dimension indgc(nbls)
      dimension gcoef(ngcd,nbls)
c
c***************************************************************
c**     FOR GENERAL CONTRACTED SHELLS
c**  this subroutine constitues the special code for
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
c*******
c
c*
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
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
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
      if(first) then
c
           do 204 icx=1,lnijkl
           do 204 i=1,nbls1
           xint=buf1(i,icx)
             ijkl=index(i)
             ngcq=indgc(ijkl)
             do 2041 iqu=1,ngcq
             buf(iqu,ijkl,icx)=xint*gcoef(iqu,ijkl)
 2041        continue
  204    continue
c
         first=.false.
      else
c
           do 205 icx=1,lnijkl
           do 205 i=1,nbls1
           xint=buf1(i,icx)
             ijkl=index(i)
             ngcq=indgc(ijkl)
             do 2051 iqu=1,ngcq
             buf(iqu,ijkl,icx)=buf(iqu,ijkl,icx)+xint*gcoef(iqu,ijkl)
 2051        continue
  205    continue
c
      endif
c
c* relise memory
c
      call retmem(2)
c*
      return
      end
c====================================================================
      subroutine make_map(ibl,ijbl,nbl2,npar,bl,map_ij,ijcs_b,ijcs_e)
      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension ijbl(nbl2,*), npar(nbl2)
c------------------------------------------------------------------
c find dimensions :
c
      nparij=npar(ibl)
c
      ijcs1=ijbl(ibl,1)
      ijcs2=ijbl(ibl,nparij)
c
      ijdim=ijcs2-ijcs1+1
      call getmem(ijdim,map_ij)
      ijcs_b=ijcs1
      ijcs_e=ijcs2
      call map_pairs(ibl,nparij,ijbl,nbl2,ijcs1,ijcs2,bl(map_ij))
c
c     write(6,*)' make_map: npairs=',npairs,' ijcs1,2=',ijcs1,ijcs2
c
      end
c================================================================
      subroutine map_pairs(ibl,nparij,ijbl,nbl2,ijcs1,ijcs2,map_ij)
      dimension ijbl(nbl2,*), map_ij(ijcs1:ijcs2)
c----------------------------------------------------------------
c construct a mapping array that shows map(ijcs)--> ijpar (pnl):
c needed only for pnl
c----------------------------------------------------------------
      do 50 ijpar=1,nparij
      ijcs=ijbl(ibl,ijpar)
      map_ij(ijcs)=ijpar
   50 continue
c
      end
c================================================================
      subroutine ab_prim_1(ibl,nparij,ijbl,nbl2,
     *                     datnuc,datbas,inx,iis,jjs,
     *                     abprim,ijcont,
     *                     nij_uniqe,ij_uniqe,map,ij1,ij2)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension datbas(13,*),datnuc(5,*)
      dimension abprim(nparij,ijcont,3)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension ij_uniqe(nij_uniqe), map(ij1:ij2)
c
c     data pi23 / 0.258012275487719489d0 /
c---------------------------------------------------------------
ctest
c     write(6,*)' from ab_prim : ibl=',ibl,' nij_uniqe=',nij_uniqe
c     call prt_uniqe(ij_uniqe,nij_uniqe)
c     write(6,*)' from ab_prim : ij1=',ij1,' ij2=',ij2   
ctest
c
cccc  do 100 ijpar=1,nparij
cccc  ijcs=ijbl(ibl,ijpar)
      do 100 ij_u =1,nij_uniqe
      ijcs=ij_uniqe(ij_u)
      ijpar=map(ijcs)
c
      ics=iis(ijcs)
      jcs=jjs(ijcs)
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
        xab=datnuc(2,iat)-datnuc(2,jat)
        yab=datnuc(3,iat)-datnuc(3,jat)
        zab=datnuc(4,iat)-datnuc(4,jat)
        rr=xab*xab+yab*yab+zab*zab
        ij=0
           do 200 is=ia,ie
           aa=datbas(1,is)
c----------------------------------
             do 200 js=ja,je
             bb=datbas(1,js)
             ij=ij+1
c----------------------------------
c Sab=(a*b)**3/4 * (a+b)**-1 * exp( -ab/(a+b) * Rab**2)
c
             if(aa.gt.zero .and. bb.gt.zero) then
                axb=aa*bb
                apb=aa+bb
                apb1=one/apb
                e=axb*apb1
                abprim(ijpar,ij,1)=apb
                abprim(ijpar,ij,2)=apb1
                abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp(-e*rr)
             else
c center23:
                if(aa.le.zero .and. bb.gt.zero) then
                   write(8,*)' exp_a=zero ; ics=',ics
                   apb=aa+bb
                   apb1=one/apb
                   abprim(ijpar,ij,1)=apb
                   abprim(ijpar,ij,2)=apb1
                   abprim(ijpar,ij,3)=sqrt(half*sqrt(bb)) ! sqrt(Sbb)
                endif
                if(aa.gt.zero .and. bb.le.zero) then
                   write(8,*)' exp_b=zero ; jcs=',jcs
                   apb=aa+bb
                   apb1=one/apb
                   abprim(ijpar,ij,1)=apb
                   abprim(ijpar,ij,2)=apb1
                   abprim(ijpar,ij,3)=sqrt(half*sqrt(aa)) ! sqrt(Saa)
                endif
                if(aa.le.zero .and. bb.le.zero) then
                   abprim(ijpar,ij,1)=zero
                   abprim(ijpar,ij,2)=zero
                   abprim(ijpar,ij,3)=zero
                endif
             endif
c----------------------------------
  200     continue
  100 continue
      end
c============================================================
      subroutine precalc2_1(isupb,bl,mmax,mmax1,nhabcd,nfumax,nbl2,nbls,
     *                      inx,iis,jjs,ijbl,nblok1,
     *                      ibl,nijbeg,nijend,npij,
     *                      kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c-----------------------------------------------------------
      common /begin/ ijbegin,klbegin
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpartot,lpareal
c
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
      common /map4uniq/ map_ij,ijcs_b,ijcs_e, map_kl,klcs_b,klcs_e
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nblok1(2,*)
c---------------------------------------------------------------
      common /pnl008/ last_isupb,last_nparts
c---------------------------------------------------------------
c     if(isupb.eq.last_isupb) then
c        lpartot=lpartot+1
c        if(kbl.eq.ibl) then
c            kabprim=iabprim
c            kldim  =ijdim
c            klpar1 =ijpar1
c            call precdiag
c        endif
c        RETURN
c     endif
c---------------------------------------------------------------
c get unique pairs :
c
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
      call get_uniq(kbl,nkl_uniqe,kl_uniqe_p,bl)
c---------------------------------------------------------------
c For PNL applications there are NO SMALL BLOCKS thus, no re-use !
CPNL  call reuseij(ibl,nbls,mmax,nblok1,reuse)
c-----------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c-----------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      lpareal=lpareal+1
      call precal2a_1(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npij,
     *              ibl,ijbl,nbl2,nijbeg,nijend,
     *              bl(iabprim      ),ijpar1 ,
     *              lcij,bl(iaa),bl(ibb),bl(ieab),bl(icis),bl(icjs),
     *              bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *              bl(igci),bl(igcj),ngci1,ngcj1,'left ',
     *              nij_uniqe,bl(ij_uniqe_p),
     *              bl(map_ij),ijcs_b,ijcs_e)
c----------
c for abnia
      if(mmax.gt.2) then
         call precal2b_1(mmax1,lcij,npij, bl(i1apb),
     *                 ijpar1,ijbegin, bl(iabnia),
     *                 nij_uniqe,bl(ij_uniqe_p),
     *                 bl(map_ij),ijcs_b,ijcs_e)
      endif
c----------------------------------------
c precalculations for pairs kl
c
      if(kbl.eq.ibl) then
          kabprim=iabprim
          kldim  =ijdim
          klpar1 =ijpar1
      endif
      if(npkl.ne.0) then
          call precal2x_1(kabprim,kldim,icpd,i1cpd,iscd)
c
          lpartot=lpartot+1
          lpareal=lpareal+1
          call precal2a_1(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npkl,
     *                  kbl,ijbl,nbl2,nklbeg,nklend,
     *                  bl(kabprim      ),klpar1 ,
     *                  lckl,bl(icc),bl(idd),bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igck),bl(igcl),ngck1,ngcl1,'right',
     *                  nkl_uniqe,bl(kl_uniqe_p),
     *                  bl(map_kl),klcs_b,klcs_e)
c--
c for cdnia
         if(mmax.gt.2) then
            call precal2b_1(mmax1,lckl,npklx,bl(i1cpd),
     *                    klpar1,klbegin, bl(icdnia),
     *                    nkl_uniqe,bl(kl_uniqe_p),
     *                    bl(map_kl),klcs_b,klcs_e)
         endif
      else
          call precdiag
      endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_1(npij,npklx,npkl, bl(i1apb),bl(i1cpd),
     *                 ijpar1,klpar1, bl(ihabcd),nhabcd,nfumax ,
     *               nij_uniqe,bl(ij_uniqe_p),bl(map_ij),ijcs_b,ijcs_e,
     *               nkl_uniqe,bl(kl_uniqe_p),bl(map_kl),klcs_b,klcs_e)
      endif
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_1(datbas,datnuc,iis,jjs,inx,npij,
     * ibl,ijbl,nbl2,nijbeg,nijend,
     * abprim,ijpar1,lcij, aaa,bbb,estab,cis,cjs,
     * xab,xparij,coefij,factij,txab, gci,gcj,ngci1,ngcj1,which ,
     *                nij_uniqe,ij_uniqe,map_ij,ijcs_b,ijcs_e)
c
      implicit real*8 (a-h,o-z)
      character*5 which
c
      COMMON /types/itype,jtype,ktype,ltype,ityp,jtyp,ktyp,ltyp
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /neglect/ eps,eps1,epsr
c
      dimension abprim(ijpar1,lcij,3)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension aaa(npij,*),bbb(npij,*), cis(npij,*),cjs(npij,*)
      dimension xab(npij,3), xparij(npij,3,lcij,3)
      dimension estab(npij,lcij)
      dimension coefij(npij,lcij), factij(npij,lcij)
      dimension xa(3),xb(3),txab(npij,3,*)
      dimension gci(npij,ngci1,*),gcj(npij,ngcj1,*)
c
      dimension ij_uniqe(nij_uniqe), map_ij(ijcs_b:ijcs_e)
c
c---------------------------------------------------------------
cnot  par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
c
      if(which.eq.'left ') then
          itypp=ityp
          jtypp=jtyp
          ngcii=ngci
          ngcjj=ngcj
          nqii=nqi
          nqjj=nqj
          lcii=lci
          lcjj=lcj
      else
          itypp=ktyp
          jtypp=ltyp
          ngcii=ngck
          ngcjj=ngcl
          nqii=nqk
          nqjj=nql
          lcii=lck
          lcjj=lcl
      endif
c
c-------------------------------------------------------------
c for integral neglect (to massage contraction coef.)
c
      contrij=dble(lcii*lcjj)
      contrij=sqrt(contrij)
c-------------------------------------------------------------
c* for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c---------------
c precalculations for the pairs IJ :
c
cccc  ijpar=0
cccc  do 100 ijp=nijbeg,nijend
cccc  ijpar=ijpar+1
cccc    ijcs=ijbl(ibl,ijp)
      do 100 ij_u=1,nij_uniqe
        ijcs=ij_uniqe(ij_u)
        ijp  =map_ij(ijcs)
        ijpar=ijp
c
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        fact1=one
cpnl    if(ics.eq.jcs) then
cpnl       if(ngctot.eq.0) fact1=fact1*half
cpnl    endif
c
c      starting contr
        ia=inx(1,ics)+1
        ja=inx(1,jcs)+1
c      last contr
        ie=inx(5,ics)
        je=inx(5,jcs)
c       number of general contr.
c       ngci=inx(4,ics)
c       ngcj=inx(4,jcs)
c in the common block contr
c
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iatom.eq.0) iatom=jatom
        if(jatom.eq.0) jatom=iatom
c
C?????    rab2=0.d0
          do 150 i=1,3
          if(iatom.eq.0 .and. jatom.eq.0) then
            xa(i)=zero
            xb(i)=zero
          else
            xa(i)=datnuc(1+i,iatom)
            xb(i)=datnuc(1+i,jatom)
          endif
          xab(ijpar,i)=xa(i)-xb(i)
C?????    rab2=rab2+xab(ijpar,i)*xab(ijpar,i)
  150     continue
c
         ji=0
         is1=0
         do 200 is=ia,ie
            is1=is1+1
            aa=datbas(1,is)
            aaa(ijpar,is1)=aa
               if(ngctot.eq.0) then
                  csi=datbas(2,is)
                  cpi=datbas(3,is)
                  coefi=csi
                  est_is=csi
                  est_ip=1.d0
                  if(itypp.eq.3) then
                     est_ip=cpi
                     coefi=cpi
                     csi=csi/cpi
                     facti=csi
                  endif
                  coefi=coefi*fact1
                  if(which.eq.'right') coefi=coefi*par268
                  cis(ijpar,is1)=csi
               else              
c* gen.contr. shell is somewhere
                  est_is=one
                  do 210 ig=0,ngcii
                  gci(ijpar,ig+1,is1)=datbas(ig+2,is)
  210             continue
               endif
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas(1,js)
            bbb(ijpar,js1)=bb
               if(ngctot.eq.0) then
                  csj=datbas(2,js)
                  cpj=datbas(3,js)
                  coefj=csj
                  est_js=csj
                  est_jp=1.d0
                  if(jtypp.eq.3) est_jp=cpj
                  if(jtypp.eq.3) coefj=cpj
c
                  coefij(ijpar,ji )=coefi*coefj
c
                  if(jtypp.eq.3) then
                     csj=csj/cpj
                     factj=csj
                     if(itypp.eq.3) then
                       factij(ijpar,ji  )=facti*factj
                     endif
                  endif
                     cjs(ijpar,js1)=csj
               else              
c* gen.contr.
                  coefij(ijpar,ji )=one
                  est_js=one
                  if(which.eq.'right') coefij(ijpar,ji)=par268
c
                  do 220 jg=0,ngcjj
                  gcj(ijpar,jg+1,js1)=datbas(jg+2,js)
  220             continue
               endif
c
            rapb=abprim(ijp,ji,2)
            sab =abprim(ijp,ji,3)
c
            aa1=aa*abprim(ijp,ji,2)
            bb1=bb*abprim(ijp,ji,2)
            coefij(ijpar,ji)=coefij(ijpar,ji)     *sab
c
            do 230 l=1,3
            xparij(ijpar,l,ji,1)=aa1*xa(l)+bb1*xb(l) ! xp(ijpar,l,ji)
            xxl=xa(l)
            if(nqii.lt.nqjj) xxl=xb(l)
            xparij(ijpar,l,ji,2)=xparij(ijpar,l,ji,1)-xxl  ! xpn
            xparij(ijpar,l,ji,3)=aa*xa(l)+bb*xb(l)         ! xpp
c center23:
            if(aa.le.zero .or. bb.le.zero) then
               xparij(ijpar,l,ji,2)=zero
            endif
  230       continue
c------------------------------------------------------------------
c For integral neglect :
c no l-shells
                if(itypp.ne.3 .and. jtypp.ne.3) then
                   factab=abs( est_is*est_js )
                endif
c l-shells
                if(itypp.eq.3 .and. jtypp.ne.3) then
                   fcts= abs( est_is*est_js )
                   fctp= abs( est_ip*est_js )
                   factab=max( fcts , fctp )
                endif
                if(itypp.ne.3 .and. jtypp.eq.3) then
                   fcts= abs( est_is*est_js )
                   fctp= abs( est_is*est_jp )
                   factab=max( fcts , fctp )
                endif
                if(itypp.eq.3 .and. jtypp.eq.3) then
                   fctss= abs( est_is*est_js )
                   fctsp= abs( est_is*est_jp )
                   fctps= abs( est_ip*est_js )
                   fctpp= abs( est_ip*est_jp )
                   factab=max(fctss,fctsp,fctps,fctpp)
                endif
c
c----
            factab=factab*contrij
            if(factab.gt.1.d0) factab=1.d0
c----
            estab(ijpar,ji)=factab*sab
            estab(ijpar,ji)=estab(ijpar,ji)*estab(ijpar,ji)
c
c------------------------------------------------------------------
c
  200    continue
  100 continue
c
ctxab
c
      if(nqii.ge.nqjj) then
         ijs1=0
         do 151 is1=1,lcii
         do 151 js1=1,lcjj
         ijs1=ijs1+1
         do 151 ijpar=1,npij
            txab(ijpar,1,ijs1)=-bbb(ijpar,js1)*xab(ijpar,1)
            txab(ijpar,2,ijs1)=-bbb(ijpar,js1)*xab(ijpar,2)
            txab(ijpar,3,ijs1)=-bbb(ijpar,js1)*xab(ijpar,3)
c center23:
            aa=aaa(ijpar,is1)
            if(aa.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  151    continue
      else
         ijs1=0
         do 152 is1=1,lcii
         do 152 js1=1,lcjj
         ijs1=ijs1+1
         do 152 ijpar=1,npij
            txab(ijpar,1,ijs1)= aaa(ijpar,is1)*xab(ijpar,1)
            txab(ijpar,2,ijs1)= aaa(ijpar,is1)*xab(ijpar,2)
            txab(ijpar,3,ijs1)= aaa(ijpar,is1)*xab(ijpar,3)
c center23:
            bb=bbb(ijpar,js1)
            if(bb.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  152    continue
      endif
c 
c for the case with mmax=2 - special cases
c
      if ( mmax.eq.2 ) then
         if (itypp.gt.1 .or. jtypp.gt.1 ) then
           do 153 ijs1=1,lcij
           do 153 ijpar=1,npij
            rapb1=abprim(nijbeg-1+ijpar,ijs1,2)
c
            txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
            txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
            txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c--------------
      end
c====================================================================
      subroutine precal2b_1(mmax1,lcij,npij,rapb,ijpar1,ijbegin,abnia,
     *                      nij_uniqe,ij_uniqe,map_ij,ijcs_b,ijcs_e)
c-------------------------------------------------------------------
c  OUTPUT
c  -------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  1. ABNIA(IJ,L,ij) -   L*( 0.5/(a+b) )  with L=1,2,...MMAX-1
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---> dimension rapb(npij,*), abnia(npij,mmax1,*)
      dimension rapb(ijpar1,*), abnia(npij,mmax1,*)
      dimension ij_uniqe(nij_uniqe), map_ij(ijcs_b:ijcs_e)
      data half /0.5d0/
c--------------------------
      ijstart=ijbegin-1
c--------------------------
      do 150 ij=1,lcij
         do 200 ij_u=1,nij_uniqe
         ijcs =ij_uniqe(ij_u)
         ijpar=map_ij(ijcs)
            apb2=half*rapb(ijstart+ijpar,ij)
            abnia(ijpar,1,ij)=apb2
               do 250 i=2,mmax1
               abnia(ijpar,i,ij)=abnia(ijpar,i-1,ij)+apb2
  250          continue
  200    continue
  150  continue
c
      end
c====================================================================
      subroutine precal2c_1(npij,npkl,npklx, rapb,rcpd,ijpar1,klpar1,
     *                      habcd,nhabcd,nfumax, 
     *                      nij_uniqe,ij_uniqe,map_ij, ijcs_b,ijcs_e,
     *                      nkl_uniqe,kl_uniqe,map_kl, klcs_b,klcs_e )
c-----------------------------------------------------------------------
c  OUTPUT
c  ------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  3. HABCD(IJ,lx,ifu,ij)  or  HABCD(KL,lx,ifu,kl) 
c
c  where IFU denotes number of elementary function from 1 up to the
c  total number of functions corresponding to the MMAX-1 /for example
c  for mmax=3 it is from 1 to 10 - s,x,y,z,xx,yy,zz,xy,xz,yz /
c  The second index lx stays for x,y or z and is used to find the power
c  of an elementary function in these directions from matrix HNIA which
c  is constant and is set up in BLOCK DATA logobsa. For example :
c
c     habcd(klpar,1=x,ifu,kl)=hnia(1,ifu)*rcpd(klpar,kl)
c
c  HABCD is used only in TRACY's recursive in subroutines TRACIJ, TRACKL
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /begin/ ijbegin,klbegin
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /logicd/ hnia(3,1)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
c---> dimension rapb(npij,*),rcpd(npkl,*)
      dimension rapb(ijpar1,*),rcpd(klpar1,*)
      dimension habcd(nhabcd,3,nfumax,*)
c
      dimension ij_uniqe(nij_uniqe), map_ij(ijcs_b:ijcs_e)
      dimension kl_uniqe(nkl_uniqe), map_kl(klcs_b:klcs_e)
c
      data half /0.5d0/
c
c-----------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c-----------------------
      if(nsij.ge.nskl) then
         do 170 kl=1,lckl
           do 180 ifu=1,nfumax
           hx=hnia(1,ifu)
           hy=hnia(2,ifu)
           hz=hnia(3,ifu)
cccccccc     do 190 klpar=1,npkl
             do 190 kl_u=1,nkl_uniqe
             klcs =kl_uniqe(kl_u)
             klpar=map_kl(klcs)
             habcd(klpar,1,ifu,kl)=hx*rcpd(klstart+klpar,kl)
             habcd(klpar,2,ifu,kl)=hy*rcpd(klstart+klpar,kl)
             habcd(klpar,3,ifu,kl)=hz*rcpd(klstart+klpar,kl)
  190        continue
  180      continue
  170    continue
      else
         do 270 ij=1,lcij
           do 280 ifu=1,nfumax
           hx=hnia(1,ifu)
           hy=hnia(2,ifu)
           hz=hnia(3,ifu)
ccccccc      do 290 ijpar=1,npij
             do 290 ij_u=1,nij_uniqe
             ijcs=ij_uniqe(ij_u)
             ijpar=map_ij(ijcs)
             habcd(ijpar,1,ifu,ij)=hx*rapb(ijstart+ijpar,ij)
             habcd(ijpar,2,ifu,ij)=hy*rapb(ijstart+ijpar,ij)
             habcd(ijpar,3,ifu,ij)=hz*rapb(ijstart+ijpar,ij)
  290        continue
  280      continue
  270    continue
      endif
c
      end
c====================================================================
      subroutine prec4neg_1(isbl_size,isbl_point,isbl_q,ijpres2,klpres2,
     *                      npij,npkl,ndiag,ij,kl,
     *                      ijpar1,lc12, klpar1,lc34,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq, 
     *                      map_ij,ijcs_b,ijcs_e,map_kl,klcs_b,klcs_e,
coutput
     *                      rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
      common /begin/ ijbegin,klbegin
c-------
      dimension map_ij(ijcs_b:ijcs_e),map_kl(klcs_b:klcs_e)
      dimension isbl_q(*),ijpres2(*),klpres2(*)
      dimension index(*)
c-------
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension xp(npij,3,*),xq(npkl,3,*)
c----output
      dimension rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
c not par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      sqrpold=one
      rpold=one
      abpcd1r=zero
      abpcdr=zero
c-----------------------------------------------------------------
      ijkl1=0
      ijklp=0
      do 100 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      if(iq.eq.0) go to 100
c
         ijklp=ijklp+1
c
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
         ijpar=map_ij(ijcsq)
         klpar=map_kl(klcsq)
c
         apb1=apb(ijstart+ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
c
         abpcd1=apb1+cpd1
         if(abpcd1.ne.abpcd1r) then
            abpcd1r=abpcd1
            abpcdr=one/abpcd1r
         endif
c
         estim=esti1*esti2*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            if(ndiag.eq.0) coef12=coef12*par268
            ijkl1=ijkl1+1
            index(ijkl1)=ijklp
c
            rppq(ijkl1)=abpcdr
            rhoapb(ijkl1)=abpcdr*cpd1
            rhocpd(ijkl1)=abpcdr*apb1
            abxcd=apb1*cpd1
            rho1=abxcd*abpcdr
c
            x1= xp(ijpar,1,ij) - xq(klpar,1,kl)
            x2= xp(ijpar,2,ij) - xq(klpar,2,kl)
            x3= xp(ijpar,3,ij) - xq(klpar,3,kl)
c
            rr2=x1*x1 + x2*x2 + x3*x3
            rys(ijkl1)=rr2*rho1
c
              if(abpcdr.ne.rpold) then
                rpold=abpcdr
                sqrpold=sqrt(rpold)
              endif
c
            const(ijkl1)=coef12*sqrpold
cpnl        if(ndiag.eq.0.and.ijpar.eq.klpar) then
cpnl           const(ijkl1)=const(ijkl1)*0.5d0
cpnl        endif
         endif
 100  continue
c
      nbls1=ijkl1
c
      end
c====================================================================
      subroutine precspec_1(isbl_size,isbl_point,
     *                      isbl_q,ijpres2,klpres2,
     *                      npij,npkl,ndiag, ij,kl,
     *                      ijpar1,lc12, klpar1,lc34, indxij,indxkl,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq,
     *                      rapb,rcpd,txab,txcd,
     *                      map_ij,ijcs_b,ijcs_e,map_kl,klcs_b,klcs_e,
c output
     *                      rys,const,xpqr,txxr,nbls1,index)
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
      common /begin/ ijbegin,klbegin
c---------------
      dimension map_ij(ijcs_b:ijcs_e),map_kl(klcs_b:klcs_e)
      dimension isbl_q(*),ijpres2(*),klpres2(*)
      dimension indxij(*),indxkl(*),index(*)
c---------------
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension rapb(ijpar1,*),rcpd(klpar1,*)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
c-------
      dimension txab(npij,3,*),txcd(npkl,3,*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
      par268=     8.d0
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      sqrpold=one
      rpold=one
c
      abpcd1r=zero
      abpcdr=zero
c-----
      ijkl1=0
      ijklp=0
      do 100 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      if(iq.eq.0) go to 100
c
         ijklp=ijklp+1
c
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
         ijpar=map_ij(ijcsq)
         klpar=map_kl(klcsq)
c
         apb1=apb(ijstart+ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
c
         abpcd1=apb1+cpd1
         if(abpcd1.ne.abpcd1r) then
            abpcd1r=abpcd1
            abpcdr=one/abpcd1r
         endif
c
         estim=esti1*esti2*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            if(ndiag.eq.0) coef12=coef12*par268
            ijkl1=ijkl1+1
            index(ijkl1)=ijklp
c
            abxcd=apb1*cpd1
            rho1=abxcd*abpcdr
c
            x1= xp(ijpar,1,ij) - xq(klpar,1,kl)
            x2= xp(ijpar,2,ij) - xq(klpar,2,kl)
            x3= xp(ijpar,3,ij) - xq(klpar,3,kl)
c
            xpqr(1,ijkl1)=x1*rho1
            xpqr(2,ijkl1)=x2*rho1
            xpqr(3,ijkl1)=x3*rho1
c
            rr2=x1*x1 + x2*x2 + x3*x3
            rys(ijkl1)=rr2*rho1
c
              if(abpcdr.ne.rpold) then
                rpold=abpcdr
                sqrpold=sqrt(rpold)
              endif
c
            const(ijkl1)=coef12*sqrpold
cpnl        if(ndiag.eq.0.and.ijpar.eq.klpar) then
cpnl           const(ijkl1)=const(ijkl1)*0.5d0
cpnl        endif
c
         endif
 100  continue
c
      nbls1=ijkl1
      if(nbls1.eq.0) return
c
cnew
      if (ityp.gt.1 .or. jtyp.gt.1 ) then
        do 210 i=1,nbls1
        ijkl=index(i)
        ijpar=indxij(ijkl)
c
        xpqr(1,i)=xpqr(1,i)*rapb(ijstart+ijpar,ij)
        xpqr(2,i)=xpqr(2,i)*rapb(ijstart+ijpar,ij)
        xpqr(3,i)=xpqr(3,i)*rapb(ijstart+ijpar,ij)
c
        txxr(1,i)=txab(ijpar,1,ij)
        txxr(2,i)=txab(ijpar,2,ij)
        txxr(3,i)=txab(ijpar,3,ij)
  210   continue
      endif
c
      if (ktyp.gt.1 .or. ltyp.gt.1 ) then
        do 220 i=1,nbls1
        ijkl=index(i)
        klpar=indxkl(ijkl)
        xpqr(1,i)=xpqr(1,i)*rcpd(klstart+klpar,kl)
        xpqr(2,i)=xpqr(2,i)*rcpd(klstart+klpar,kl)
        xpqr(3,i)=xpqr(3,i)*rcpd(klstart+klpar,kl)
c
        txxr(1,i)=txcd(klpar,1,kl)
        txxr(2,i)=txcd(klpar,2,kl)
        txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
c
      end
c====================================================================
      subroutine xwpq_1(nbls1,xwp,xwq,p1234,ijpar1,lc12, klpar1,lc34,
     *                  lcij,lckl,npij,npkl,indxij,indxkl,index,
     *                  rppq,xp,xq,xpp,xqq,
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb)
c-----------------------------------------------
c  OUTPUT
c
c  XWP(ijkl1) - coordinates of W-P,  where W=(XPP+XQQ)/(a+b+c+d)
c  XWQ(ijkl1) - coordinates of W-Q
c  P1234(ijkl1)=(txab+txcd)*(rcpd OR rapb)
c  ABCD(ijkl1)= apb*rcpd   OR  cpd*rapb
c-----------------------------------------------
c
      implicit real*8 (a-h,o-z)
cflops
      common /flops/ iflop(20)
cnowy
      common /begin/ ijbegin,klbegin
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
c
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension rppq(*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy dimension apb(npij,*),rapb(npij,*),cpd(npkl,*),rcpd(npkl,*)
      dimension apb(ijpar1,*),rapb(ijpar1,*)
      dimension cpd(klpar1,*),rcpd(klpar1,*)
cnowy
      dimension abcd(nbls1)
C
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
C
      do 100 i=1,nbls1
      ijkl=index(i)
      ijpar=indxij(ijkl)
      klpar=indxkl(ijkl)
c
      rppq1=rppq(i)
c
      xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
      ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
      zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
      XWP(i,1)=xwl-XP(ijpar,1,lcij)
      XWP(i,2)=ywl-XP(ijpar,2,lcij)
      XWP(i,3)=zwl-XP(ijpar,3,lcij)
c
      XWQ(i,1)=XWL-XQ(klpar,1,lckl)
      XWQ(i,2)=YWL-XQ(klpar,2,lckl)
      XWQ(i,3)=ZWL-XQ(klpar,3,lckl)
c
  100 CONTINUE
c
c
c  for tracy1 :
c
      if(nsij.ge.nskl) then
c
          do 210 i=1,nbls1
          ijkl=index(i)
          ijpar=indxij(ijkl)
          klpar=indxkl(ijkl)
          rcpd1=rcpd(klstart+klpar,lckl)
c
          abcd(i)=apb(ijstart+ijpar,lcij)*rcpd1
c
          p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
          p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
          p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
  210     continue
c
      else
c
          do 310 i=1,nbls1
          ijkl=index(i)
          ijpar=indxij(ijkl)
          klpar=indxkl(ijkl)
          rapb1=rapb(ijstart+ijpar,lcij)
c
          abcd(i)=cpd(klstart+klpar,lckl)*rapb1
c
          p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rapb1
          p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rapb1
          p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rapb1
  310     continue
c
      endif
cflops
cxxx   iflop(7)=iflop(7)+18*nbls1
c
      return
      end
c====================================================================
      subroutine specase_1(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                     npij,npkl,ii,jj,kk,ll,
     *                     cis,cjs,cks,cls,
     *                     buf,buf1, const,rysx,xpqr,txxr,concoe)
c****
      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
c******************************************************
c
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension cis(npij,*),cjs(npij,*),cks(npkl,*),cls(npkl,*)
      dimension const(*),rysx(*),concoe(*)
      dimension buf(nbls,*),buf1(nbls1,4)
c
c***************************************************************
c**  this subroutine constitues the special code for
c**  two types of integrals over nbls quartets of primitive shells
c**  1. (ss|ss)
c**  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c**  these integrals are also contracted here.
c**
c**  this routine is called from the twoe subroutine.
c**
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
c**  buf(nbls,*) -contains final contracted integrals
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
c*******
c
c*
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
      intct=0
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
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
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
        if( max(ityp,jtyp,ktyp,ltyp).eq.3) then
          if(ityp.eq.3) then
            do 1051 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)
            concoe(i)=const(i)*cis(ijpar,ii)
 1051       continue
          endif
          if(jtyp.eq.3) then
            do 1052 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)
            concoe(i)=const(i)*cjs(ijpar,jj)
 1052       continue
          endif
          if(ktyp.eq.3) then
            do 1053 i=1,nbls1
            ijkl=index(i)
            klpar=indxkl(ijkl)
            concoe(i)=const(i)*cks(klpar,kk)
 1053       continue
          endif
          if(ltyp.eq.3) then
            do 1054 i=1,nbls1
            ijkl=index(i)
            klpar=indxkl(ijkl)
            concoe(i)=const(i)*cls(klpar,ll)
 1054       continue
          endif
c*
           do 105 i=1,nbls1
           buf1(i,1)=concoe(i)*bl(if00+i)
  105      continue
           intct=1
        endif
c*************
        do 106 i=1,nbls1
          buf1(i,intct+1)=-xpqr(1,i)*const(i)
          buf1(i,intct+2)=-xpqr(2,i)*const(i)
          buf1(i,intct+3)=-xpqr(3,i)*const(i)
106     continue
c
c**********************************************************
c
  203 continue
      if(first) then
         do 204 icx=1,lnijkl
         do 204 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf1(i,icx)
  204    continue
         first=.false.
      else
         do 205 icx=1,lnijkl
         do 205 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf(ijkl,icx)+buf1(i,icx)
  205    continue
      endif
c
c
c* relise memory
c
      call retmem(2)
c*
      return
      end
c====================================================================
      subroutine precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c-------------------------------------------------------------
c* precalculations for the pairs IJ :
c
      iapb =iabprim
      i1apb=iabprim   + ijdim
      isab =iabprim   + ijdim*2
c--------------
      end
c====================================================================
      subroutine ab_prim_2(ibl,nparij,ijbl,nbl2,
     *                     datnuc,datbas,inx,iis,jjs,
     *                     apb,rapb,sab,ijcont,
     *                     nij_uniqe,ij_uniqe,map,ij1,ij2)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension datbas(13,*),datnuc(5,*)
cccc  dimension abprim(nparij,ijcont,3)
      dimension apb(ijcont),rapb(ijcont),sab(nparij,ijcont)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension ij_uniqe(nij_uniqe), map(ij1:ij2)
c
      dimension sqrtx(400),eexx(400)
c---------------------------------------------------------------
cccc  do 100 ijpar=1,nparij
cccc  ijcs=ijbl(ibl,ijpar)
cnew
      do 100 ij_u =1,nij_uniqe
      ijcs=ij_uniqe(ij_u)
      ijpar=map(ijcs)
c
      ics=iis(ijcs)
      jcs=jjs(ijcs)
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
        xab=datnuc(2,iat)-datnuc(2,jat)
        yab=datnuc(3,iat)-datnuc(3,jat)
        zab=datnuc(4,iat)-datnuc(4,jat)
        rr=xab*xab+yab*yab+zab*zab
           ij=0
           do 200 is=ia,ie
           aa=datbas(1,is)
c----------------------------------
             do 200 js=ja,je
             bb=datbas(1,js)
             ij=ij+1
c----------------------------------
             if(aa.gt.zero .and. bb.gt.zero) then
cnew            if(ijpar.eq.1) then
                if(ij_u .eq.1) then
                  axb=aa*bb
                  apb(ij)=aa+bb
                  rapb(ij)=one/apb(ij)
                  e=axb*rapb(ij)
                  eexx(ij)=e
                  sqrtx(ij)=sqrt(sqrt(axb))**3
                  sqrt3=sqrtx(ij)
                else
                  sqrt3=sqrtx(ij)
                  e=eexx(ij)
                endif
                sab(ijpar,ij)=rapb(ij)*sqrt3*exp(-e*rr)
             else
c center23:
                if(aa.le.zero .and. bb.gt.zero) then
                   write(8,*)' exp_a=zero ; ics=',ics
cnew               if(ijpar.eq.1) then
                   if(ij_u .eq.1) then
                      apb(ij)=aa+bb
                      rapb(ij)=one/apb(ij)
                      sqrtx(ij)=sqrt(half*sqrt(bb)) ! sqrt(Sbb)
                      sqrt3=sqrtx(ij)
                   else
                      sqrt3=sqrtx(ij)
                   endif
                   sab(ijpar,ij)=sqrt3                ! sqrt(Sbb)
                endif
                if(aa.gt.zero .and. bb.le.zero) then
                   write(8,*)' exp_b=zero ; jcs=',jcs
cnew               if(ijpar.eq.1) then
                   if(ij_u .eq.1) then
                      apb(ij)=aa+bb
                      rapb(ij)=one/apb(ij)
                      sqrtx(ij)=sqrt(half*sqrt(aa)) ! sqrt(Saa)
                      sqrt3=sqrtx(ij)
                   else
                      sqrt3=sqrtx(ij)
                   endif
                   sab(ijpar,ij)=sqrt3              ! sqrt(Saa)
                endif
                if(aa.le.zero .and. bb.le.zero) then
                   apb(ij)=zero
                   rapb(ij)=zero
                   sab(ijpar,ij)=zero
                endif
             endif
c
c----------------------------------
  200     continue
  100 continue
      end
c============================================================
cccc  subroutine precalc2_2(isupb,bl,mmax,mmax1,nhabcd,nfumax,nbl2,nbls,
      subroutine precalc2_2(isupb,bl,mmax,mmax1,       nfumax,nbl2,nbls,
     *                      inx,iis,jjs,ijbl,nblok1,
     *                      ibl,nijbeg,nijend,npij,
     *                      kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c-----------------------------------------------------------
      common /begin/ ijbegin,klbegin
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpartot,lpareal
c
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
c
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5d/ iabnix,icdnix,ixpnx,ixqnx,ihabcdx
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c
      common /map4uniq/ map_ij,ijcs_b,ijcs_e, map_kl,klcs_b,klcs_e
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nblok1(2,*)
c---------------------------------------------------------------
      common /pnl008/ last_isupb,last_nparts
c---------------------------------------------------------------
c     if(isupb.eq.last_isupb) then
c        lpartot=lpartot+1
c        if(kbl.eq.ibl) then
c            kabprim=iabprim
c            kldim  =ijdim
c            klpar1 =ijpar1
c            call precdiag
c        endif
c        RETURN
c     endif
c---------------------------------------------------------------
c get unique pairs :
c
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
      call get_uniq(kbl,nkl_uniqe,kl_uniqe_p,bl)
c---------------------------------------------------------------
c For PNL applications there are NO SMALL BLOCKS thus, no re-use !
CPNL  call reuseij(ibl,nbls,mmax,nblok1,reuse)
c-----------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c-----------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_2(iabprim,lcij,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      lpareal=lpareal+1
      call precal2a_2(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npij,
     *              ibl,ijbl,nbl2,nijbeg,nijend,
     *                  bl(iabprim+lcij),bl(iabprim+2*lcij),ijpar1,
     *                  lcij,bl(ieab), bl(icis),bl(icjs),
     *              bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *                  bl(igcij),ngci1,ngcj1,'left ',
     *              nij_uniqe,bl(ij_uniqe_p),
     *              bl(map_ij),ijcs_b,ijcs_e)
c----------
c for abnia
      if(mmax.gt.2) then
         call precal2b_2(mmax1,lcij, bl(i1apb),bl(iabnia))
      endif
c----------------------------------------
c precalculations for pairs kl
c
      if(kbl.eq.ibl) then
          kabprim=iabprim
          kldim  =ijdim
          klpar1 =ijpar1
      endif
      if(npkl.ne.0) then
          call precal2x_2(kabprim,lckl,icpd,i1cpd,iscd)
c
          lpartot=lpartot+1
          lpareal=lpareal+1
          call precal2a_2(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npkl,
     *                  kbl,ijbl,nbl2,nklbeg,nklend,
     *                  bl(kabprim+lckl),bl(kabprim+2*lckl),klpar1,
     *                  lckl,bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igckl),ngck1,ngcl1,'right',
     *                  nkl_uniqe,bl(kl_uniqe_p),
     *                  bl(map_kl),klcs_b,klcs_e)
c--
c for cdnia
         if(mmax.gt.2) then
            call precal2b_2(mmax1,lckl, bl(i1cpd),bl(icdnia))
         endif
      else
          call precdiag
      endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_2(lcij,lckl,bl(i1apb),bl(i1cpd),
     *                   bl(ihabcd),nfumax )
      endif
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_2(datbas,datnuc,iis,jjs,inx,npij,
     *                    ibl,ijbl,nbl2,nijbeg,nijend,
     *                    rapb,sab,ijpar1, lcij,estab, cis,cjs,
     *                    xab,xparij,coefij,factij,txab,
     *                    gcij,ngci1,ngcj1,which ,
     *                nij_uniqe,ij_uniqe,map_ij,ijcs_b,ijcs_e)
c
      implicit real*8 (a-h,o-z)
      character*5 which
c
      COMMON /types/itype,jtype,ktype,ltype,ityp,jtyp,ktyp,ltyp
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lc12,lc34
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /neglect/ eps,eps1,epsr
c---------------------------------------------------------------
      dimension rapb(lcij),sab(ijpar1,lcij) 
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension cis(*),cjs(*)
      dimension xab(npij,3), xparij(npij,3,lcij,3)
      dimension estab (npij,lcij)
      dimension coefij(npij,lcij), factij(lcij)
      dimension xa(3),xb(3),txab(npij,3,*)
      dimension gcij(ngci1,ngcj1,lcij)
      dimension aexp(100),bexp(100)
c
      dimension ij_uniqe(nij_uniqe), map_ij(ijcs_b:ijcs_e)
c---------------------------------------------------------------
cnot  par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
c
      if(which.eq.'left ') then
          itypp=ityp
          jtypp=jtyp
          ngcii=ngci
          ngcjj=ngcj
          nqii=nqi
          nqjj=nqj
          lcii=lci
          lcjj=lcj
      else
          itypp=ktyp
          jtypp=ltyp
          ngcii=ngck
          ngcjj=ngcl
          nqii=nqk
          nqjj=nql
          lcii=lck
          lcjj=lcl
      endif
c
c-------------------------------------------------------------
c for integral neglect (to massage contraction coef.)
c
      contrij=dble(lcii*lcjj)
      contrij=sqrt(contrij)
c-------------------------------------------------------------
c* for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c---------------
c precalculations for the pairs IJ :
c
      do 100 ij_u=1,nij_uniqe
        ijcs=ij_uniqe(ij_u)
        ijp  =map_ij(ijcs)
        ijpar=ijp
c
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        fact1=one
cpnl    if(ics.eq.jcs) then
cpnl       if(ngctot.eq.0) fact1=fact1*half
cpnl    endif
c
c      starting contr
        ia=inx(1,ics)+1
        ja=inx(1,jcs)+1
c      last contr
        ie=inx(5,ics)
        je=inx(5,jcs)
c       number of general contr.
c       ngci=inx(4,ics)
c       ngcj=inx(4,jcs)
c in the common block contr
c
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iatom.eq.0) iatom=jatom
        if(jatom.eq.0) jatom=iatom
c
          do 150 i=1,3
          if(iatom.eq.0 .and. jatom.eq.0) then
            xa(i)=zero
            xb(i)=zero
          else
            xa(i)=datnuc(1+i,iatom)
            xb(i)=datnuc(1+i,jatom)
          endif
          xab(ijpar,i)=xa(i)-xb(i)
  150     continue
c
         ji=0
         is1=0
         do 200 is=ia,ie
            is1=is1+1
            aa=datbas(1,is)
            aexp(is1)=aa
               if(ngctot.eq.0) then
                  csi=datbas(2,is)
                  cpi=datbas(3,is)
                  coefi=csi
                  est_is=csi
                  est_ip=one
                  if(itypp.eq.3) then
                     est_ip=cpi
                     coefi=cpi
                     csi=csi/cpi
                     facti=csi
                  endif
                  coefi=coefi*fact1
                  if(which.eq.'right') coefi=coefi*par268
                  cis(is1)=csi
               endif
c
         js1=0
         do 200 js=ja,je
            js1=js1+1
            ji=ji+1
            bb=datbas(1,js)
            bexp(js1)=bb
               if(ngctot.eq.0) then
                  csj=datbas(2,js)
                  cpj=datbas(3,js)
                  coefj=csj
                  est_js=csj
                  est_jp=one
                  if(jtypp.eq.3) est_jp=cpj
                  if(jtypp.eq.3) coefj=cpj
                  coefij(ijpar,ji )=coefi*coefj
                  if(jtypp.eq.3) then
                     csj=csj/cpj
                     factj=csj
                     if(itypp.eq.3) then
                       factij(ji)=facti*factj
                     endif
                  endif
                  cjs(js1)=csj
               else              
c                 gen.contr.
                  est_is=one
                  est_js=one
                  coefij(ijpar,ji )=one
                  if(which.eq.'right') coefij(ijpar,ji)=par268
               endif
c
c--------------------------------
c for general contractions :
               if(ngctot.gt.0 .and. ijpar.eq.1) then
                  do 210 ig=0,ngcii
                  gci=datbas(ig+2,is)
                     do 210 jg=0,ngcjj
                       gcj=gci*datbas(jg+2,js)
                       gcij(ig+1,jg+1,ji)=gcj
  210                continue
               endif
c--------------------------------
c
            rapb1=rapb(ji)
            sab1 =sab(ijp,ji)
            aa1=aa*rapb1
            bb1=bb*rapb1
c
c-overlap   coefij(ijpar,ji)=coefij(ijpar,ji)*rapb1*sab1
            coefij(ijpar,ji)=coefij(ijpar,ji)      *sab1
c
            do 230 l=1,3
            xparij(ijpar,l,ji,1)=aa1*xa(l) + bb1*xb(l)  ! xp(ijpar,l,ji
                xxl=xa(l)
                if(nqii.lt.nqjj) xxl=xb(l)
            xparij(ijpar,l,ji,2)=xparij(ijpar,l,ji,1)-xxl  ! xpn
            xparij(ijpar,l,ji,3)=aa*xa(l)+bb*xb(l)         ! xpp
c center23:
            if(aa.le.zero .or. bb.le.zero) then
               xparij(ijpar,l,ji,2)=zero
            endif
c
  230       continue
c------------------------------------------------------------------
c For integral neglect :
c no l-shells
                if(itypp.ne.3 .and. jtypp.ne.3) then
                   factab=abs( est_is*est_js )
                endif
c l-shells
                if(itypp.eq.3 .and. jtypp.ne.3) then
                   fcts= abs( est_is*est_js )
                   fctp= abs( est_ip*est_js )
                   factab=max( fcts , fctp )
                endif
                if(itypp.ne.3 .and. jtypp.eq.3) then
                   fcts= abs( est_is*est_js )
                   fctp= abs( est_is*est_jp )
                   factab=max( fcts , fctp )
                endif
                if(itypp.eq.3 .and. jtypp.eq.3) then
                   fctss= abs( est_is*est_js )
                   fctsp= abs( est_is*est_jp )
                   fctps= abs( est_ip*est_js )
                   fctpp= abs( est_ip*est_jp )
                   factab=max(fctss,fctsp,fctps,fctpp)
                endif
c
c----
            factab=factab*contrij
            if(factab.gt.1.d0) factab=1.d0
c----
            estab(ijpar,ji)=factab*sab1
            estab(ijpar,ji)=estab(ijpar,ji)*estab(ijpar,ji)
c
c------------------------------------------------------------------
c
  200    continue
  100 continue
c
ctxab
c
      if(nqii.ge.nqjj) then
         ijs1=0
         do 151 is1=1,lcii
         aa1=aexp(is1)
         do 151 js1=1,lcjj
         bb1=bexp(js1)
         ijs1=ijs1+1
cold     do 151 ijpar=1,npij
         do 151 ij_u=1,nij_uniqe
         ijcs=ij_uniqe(ij_u)
         ijpar=map_ij(ijcs)
            txab(ijpar,1,ijs1)=-bb1*xab(ijpar,1)
            txab(ijpar,2,ijs1)=-bb1*xab(ijpar,2)
            txab(ijpar,3,ijs1)=-bb1*xab(ijpar,3)
c center23:
            if(aa1.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  151    continue
      else
         ijs1=0
         do 152 is1=1,lcii
         aa1=aexp(is1)
         do 152 js1=1,lcjj
         bb1=bexp(js1)
         ijs1=ijs1+1
cold     do 152 ijpar=1,npij
         do 152 ij_u=1,nij_uniqe
         ijcs=ij_uniqe(ij_u)
         ijpar=map_ij(ijcs)
            txab(ijpar,1,ijs1)= aa1*xab(ijpar,1)
            txab(ijpar,2,ijs1)= aa1*xab(ijpar,2)
            txab(ijpar,3,ijs1)= aa1*xab(ijpar,3)
c center23:
            if(bb1.le.zero) then
               txab(ijpar,1,ijs1)=zero
               txab(ijpar,2,ijs1)=zero
               txab(ijpar,3,ijs1)=zero
            endif
  152    continue
      endif
c 
c for the case with mmax=2 - special cases
c
      if ( mmax.eq.2 ) then
         if (itypp.gt.1 .or. jtypp.gt.1 ) then
           do 153 ijs1=1,lcij
           rapb1=rapb(ijs1)
cold       do 153 ijpar=1,npij
           do 153 ij_u=1,nij_uniqe
           ijcs=ij_uniqe(ij_u)
           ijpar=map_ij(ijcs)
            txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
            txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
            txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c--------------
      end
c====================================================================
      subroutine precal2b_2(mmax1,lcij, rapb, abnia)
      implicit real*8 (a-h,o-z)
      dimension rapb(*), abnia(mmax1,*)
      data half /0.5d0/
c--------------------------
c     ijstart=ijbegin-1
c--------------------------
      do 150 ij=1,lcij
            apb2=half*rapb(ij)
            abnia(1,ij)=apb2
               do 250 i=2,mmax1
               abnia(i,ij)=abnia(i-1,ij)+apb2
  250          continue
  150 continue
c
      end
c====================================================================
      subroutine precal2c_2(lcij,lckl,rapb,rcpd,habcd,nfumax)
      implicit real*8 (a-h,o-z)
c     common /begin/ ijbegin,klbegin
c     common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /logicd/ hnia(3,1)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
cccc  dimension rapb(ijpar1,*),rcpd(klpar1,*)
cccc  dimension habcd(nhabcd,3,nfumax,*)
      dimension rapb(*),rcpd(*), habcd(3,nfumax,*)
c
      data half /0.5d0/
c
c-----------------------
c     ijstart=ijbegin-1
c     klstart=klbegin-1
c-----------------------
c Do it only for one pair (the first one)
c
      if(nsij.ge.nskl) then
c-----
         do 170 kl=1,lckl
         rcpdkl=rcpd(kl)
           do 180 ifu=1,nfumax
           hx=hnia(1,ifu)
           hy=hnia(2,ifu)
           hz=hnia(3,ifu)
             habcd(1,ifu,kl)=hx*rcpdkl
             habcd(2,ifu,kl)=hy*rcpdkl
             habcd(3,ifu,kl)=hz*rcpdkl
  180      continue
  170    continue
c-----
      else
c-------
         do 270 ij=1,lcij
         rapbij=rapb(ij)
           do 280 ifu=1,nfumax
           hx=hnia(1,ifu)
           hy=hnia(2,ifu)
           hz=hnia(3,ifu)
             habcd(1,ifu,ij)=hx*rapbij
             habcd(2,ifu,ij)=hy*rapbij
             habcd(3,ifu,ij)=hz*rapbij
  280      continue
  270    continue
c-------
      endif
c
      end
c====================================================================
      subroutine prec4neg_2(isbl_size,isbl_point,isbl_q,ijpres2,klpres2,
     *                      npij,npkl,ndiag,ij,kl, lc12, lc34,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq, 
     *                      map_ij,ijcs_b,ijcs_e,map_kl,klcs_b,klcs_e,
coutput
     *                      rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
ccc   common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
c-------
      dimension map_ij(ijcs_b:ijcs_e),map_kl(klcs_b:klcs_e)
      dimension isbl_q(*),ijpres2(*),klpres2(*)
      dimension index(*)
c-------
cold  dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension apb(lc12),cpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
      dimension xp(npij,3,*),xq(npkl,3,*)
c----output
ccc   dimension rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c now scalars:  rppq,rhoapb,rhocpd  
      dimension rys(*),const(*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
c not par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
c NEWER
         apb1=apb(ij)
         cpd1=cpd(kl)
c
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
         abpcdr=one/abpcd1
         rho1=abxcd*abpcdr
c        
         sqrpold=sqrt(abpcdr)
c
         rppq  =abpcdr
         rhoapb=abpcdr*cpd1
         rhocpd=abpcdr*apb1
c-----------------------------------------------------------------
      ijkl1=0
      ijklp=0
      do 100 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      if(iq.eq.0) go to 100
c
         ijklp=ijklp+1
c
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
         ijpar=map_ij(ijcsq)
         klpar=map_kl(klcsq)
c
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
c
         estim=esti1*esti2*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            if(ndiag.eq.0) coef12=coef12*par268
            ijkl1=ijkl1+1
            index(ijkl1)=ijklp
c
            x1= xp(ijpar,1,ij) - xq(klpar,1,kl)
            x2= xp(ijpar,2,ij) - xq(klpar,2,kl)
            x3= xp(ijpar,3,ij) - xq(klpar,3,kl)
            rr2=x1*x1 + x2*x2 + x3*x3
c
            rys(ijkl1)=rr2*rho1
c
            const(ijkl1)=coef12*sqrpold
cpnl        if(ndiag.eq.0.and.ijpar.eq.klpar) then
cpnl           const(ijkl1)=const(ijkl1)*0.5d0
cpnl        endif
         endif
 100  continue
c
      nbls1=ijkl1
c
      end
c====================================================================
      subroutine precspec_2(isbl_size,isbl_point,
     *                      isbl_q,ijpres2,klpres2,
     *                      npij,npkl,ndiag, ij,kl,
     *                      lc12, lc34, indxij,indxkl,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq,
     *                      rapb,rcpd,txab,txcd,
     *                      map_ij,ijcs_b,ijcs_e,map_kl,klcs_b,klcs_e,
c output
     *                      rys,const,xpqr,txxr,nbls1,index)
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
      common /begin/ ijbegin,klbegin
c---------------
      dimension map_ij(ijcs_b:ijcs_e),map_kl(klcs_b:klcs_e)
      dimension isbl_q(*),ijpres2(*),klpres2(*)
      dimension indxij(*),indxkl(*),index(*)
c---------------
      dimension apb(lc12),cpd(lc34)
      dimension rapb(lc12),rcpd(lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
c-------
      dimension txab(npij,3,*),txcd(npkl,3,*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
      par268=     8.d0
c---------------------------------------------------------------
c NEWER
         apb1=apb(ij)
         cpd1=cpd(kl)
c
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
         abpcd1r=abpcd1
         abpcdr=one/abpcd1r
         rho1=abxcd*abpcdr
c        
         sqrpold=sqrt(abpcdr)
c------------
      ijkl1=0
      ijklp=0
      do 100 iqp=1,isbl_size
      iq=isbl_q(isbl_point+iqp)
      if(iq.eq.0) go to 100
c
         ijklp=ijklp+1
c
         ijcsq=ijpres2(iq)
         klcsq=klpres2(iq)
         ijpar=map_ij(ijcsq)
         klpar=map_kl(klcsq)
c
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
c
         estim=esti1*esti2*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            if(ndiag.eq.0) coef12=coef12*par268
            ijkl1=ijkl1+1
            index(ijkl1)=ijklp
c
            x1= xp(ijpar,1,ij) - xq(klpar,1,kl)
            x2= xp(ijpar,2,ij) - xq(klpar,2,kl)
            x3= xp(ijpar,3,ij) - xq(klpar,3,kl)
c
            xpqr(1,ijkl1)=x1*rho1
            xpqr(2,ijkl1)=x2*rho1
            xpqr(3,ijkl1)=x3*rho1
c
            rr2=x1*x1 + x2*x2 + x3*x3
            rys(ijkl1)=rr2*rho1
c
            const(ijkl1)=coef12*sqrpold
cpnl        if(ndiag.eq.0.and.ijpar.eq.klpar) then
cpnl           const(ijkl1)=const(ijkl1)*0.5d0
cpnl        endif
c
         endif
 100  continue
c
      nbls1=ijkl1
      if(nbls1.eq.0) return
c
cnew
      if (ityp.gt.1 .or. jtyp.gt.1 ) then
          rapbxxx=rapb(ij)
        do 210 i=1,nbls1
        ijkl=index(i)
        ijpar=indxij(ijkl)
        xpqr(1,i)=xpqr(1,i)*rapbxxx 
        xpqr(2,i)=xpqr(2,i)*rapbxxx 
        xpqr(3,i)=xpqr(3,i)*rapbxxx
c
        txxr(1,i)=txab(ijpar,1,ij)
        txxr(2,i)=txab(ijpar,2,ij)
        txxr(3,i)=txab(ijpar,3,ij)
  210   continue
      endif
c
      if (ktyp.gt.1 .or. ltyp.gt.1 ) then
          rcpdxxx=rcpd(kl)
        do 220 i=1,nbls1
        ijkl=index(i)
        klpar=indxkl(ijkl)
        xpqr(1,i)=xpqr(1,i)*rcpdxxx     
        xpqr(2,i)=xpqr(2,i)*rcpdxxx    
        xpqr(3,i)=xpqr(3,i)*rcpdxxx   
c
        txxr(1,i)=txcd(klpar,1,kl)
        txxr(2,i)=txcd(klpar,2,kl)
        txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
c
      end
c====================================================================
      subroutine xwpq_2(nbls1,xwp,xwq,p1234,
     1        ijpar1,lc12, klpar1,lc34,
     *        lcij,lckl,npij,npkl,
     *        indxij,indxkl,index,
     *        rppq,xp,xq,xpp,xqq,
     *        txab,txcd,abcd,apb,rcpd,cpd,rapb)
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /flops/ iflop(20)
      common /begin/ ijbegin,klbegin
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
      dimension indxij(*),indxkl(*),index(*)
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension xpp(npij,3,*),xqq(npkl,3,*)
      dimension txab(npij,3,*),txcd(npkl,3,*)
cnowy
      dimension apb(lc12),rapb(lc12)
      dimension cpd(lc34),rcpd(lc34)
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
      rppq1=rppq
C
      do 100 i=1,nbls1
      ijkl=index(i)
      ijpar=indxij(ijkl)
      klpar=indxkl(ijkl)
c
      xwl=( xpp(ijpar,1,lcij) + xqq(klpar,1,lckl) )*rppq1
      ywl=( xpp(ijpar,2,lcij) + xqq(klpar,2,lckl) )*rppq1
      zwl=( xpp(ijpar,3,lcij) + xqq(klpar,3,lckl) )*rppq1
c
      XWP(i,1)=xwl-XP(ijpar,1,lcij)
      XWP(i,2)=ywl-XP(ijpar,2,lcij)
      XWP(i,3)=zwl-XP(ijpar,3,lcij)
c
      XWQ(i,1)=XWL-XQ(klpar,1,lckl)
      XWQ(i,2)=YWL-XQ(klpar,2,lckl)
      XWQ(i,3)=ZWL-XQ(klpar,3,lckl)
c
  100 CONTINUE
c
c
c  for tracy1 :
c
      if(nsij.ge.nskl) then
c-new
          rcpd1=rcpd(lckl)
          abcd=apb(lcij)*rcpd1
c
          do 210 i=1,nbls1
          ijkl=index(i)
          ijpar=indxij(ijkl)
          klpar=indxkl(ijkl)
c
          p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
          p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
          p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
  210     continue
c
      else
c-new
          rapb1 =rapb(lcij)
          abcd  =cpd(lckl)*rapb1
c
          do 310 i=1,nbls1
          ijkl=index(i)
          ijpar=indxij(ijkl)
          klpar=indxkl(ijkl)
c
          p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rapb1
          p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rapb1
          p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rapb1
  310     continue
c
      endif
c
      end
c====================================================================
      subroutine specase_2(bl,first,nbls,nbls1, index,
     *                     npij,npkl,ii,jj,kk,ll,
     *                     cis,cjs,cks,cls,
     *                     buf,buf1, const,rysx,xpqr,txxr,concoe)
      implicit real*8 (a-h,o-z)
      logical first
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      dimension bl(*)
      dimension index(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension cis(*),cjs(*),cks(*),cls(*)
      dimension const(*),rysx(*),concoe(*)
      dimension buf(nbls,*),buf1(nbls1,4)
c---------------------------------------------------------------
c**  this subroutine constitues the special code for
c**  two types of integrals over nbls quartets of primitive shells
c**  1. (ss|ss)
c**  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c**  these integrals are also contracted here.
c**
c**  this routine is called from the twoe subroutine.
c**
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
c**  buf(nbls,*) -contains final contracted integrals
c---------------------------------------------------------------
c
c* memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
c*******
c
c*
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
      intct=0
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
        xpqr(3,i)=xpqr(2,i)*bl(if11+i) + txxr(3,i)*bl(if00+i)
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
c------------------- NEW ------------------------
        if( max(ityp,jtyp,ktyp,ltyp).eq.3) then
c
          if(ityp.eq.3) then
            coeff=cis(ii)
            do 1051 i=1,nbls1
            concoe(i)=const(i)*coeff 
 1051       continue
          endif
          if(jtyp.eq.3) then
            coeff=cjs(jj)
            do 1052 i=1,nbls1
            concoe(i)=const(i)*coeff
 1052       continue
          endif
          if(ktyp.eq.3) then
            coeff=cks(kk)
            do 1053 i=1,nbls1
            concoe(i)=const(i)*coeff 
 1053       continue
          endif
          if(ltyp.eq.3) then
            coeff=cls(ll)
            do 1054 i=1,nbls1
            concoe(i)=const(i)*coeff
 1054       continue
          endif
           do 105 i=1,nbls1
           buf1(i,1)=concoe(i)*bl(if00+i)
  105      continue
           intct=1
        endif
c------------------- NEW end --------------------
        do 106 i=1,nbls1
          buf1(i,intct+1)=-xpqr(1,i)*const(i)
          buf1(i,intct+2)=-xpqr(2,i)*const(i)
          buf1(i,intct+3)=-xpqr(3,i)*const(i)
106     continue
c------------------------------------------------
  203 continue
      if(first) then
         do 204 icx=1,lnijkl
         do 204 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf1(i,icx)
  204    continue
         first=.false.
      else
         do 205 icx=1,lnijkl
         do 205 i=1,nbls1
         ijkl=index(i)
         buf(ijkl,icx)=buf(ijkl,icx)+buf1(i,icx)
  205    continue
      endif
c
c relise memory
c
      call retmem(2)
c
      end
c====================================================================
      subroutine precal2x_2(iabprim,lcij, iapb,i1apb,isab)
c-------------------------------------------------------------
c* precalculations for the pairs IJ :
c
c     iapb =iabprim
c     i1apb=iabprim   + ijdim
c     isab =iabprim   + ijdim*2
      iapb =iabprim
      i1apb=iabprim+lcij
      isab =iabprim+2*lcij
c--------------
      end
c====================================================================
