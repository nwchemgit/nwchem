      subroutine prec2ij(ibl, bl,inx,npar,nbl2, iis,jjs,ijbl)
c $Id: precalc2.F,v 1.5 1996-01-19 07:30:44 d3g681 Exp $
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primij/ iabprim, ijdim ,ijpar1 
      common /time0/ tprec2
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
c 
c on return : nij_uniqe=no of unique pairs and pointer to the list
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
      nparij=npar(ibl)
c---------------------------------------
      call dimenij(ibl,inx,nparij,nbl2, iis,jjs,ijbl,ijdim,ijcont)
c---------------------------------------
      ijpar1=nparij
      call getmem(3*ijdim,iabprim)
c---------------------------------------
         call ab_prim(ibl,nparij,ijbl,nbl2,
     *                dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                bl(iabprim),ijcont,
     *                nij_uniqe,bl(ij_uniqe_p) )
c
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
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primkl/ kabprim, kldim ,klpar1 
      common /time0/ tprec2
c
      dimension bl(*)
      dimension inx(12,*),iis(*),jjs(*), ijbl(nbl2,*), npar(*)
c---------------------------------------
      call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
c---------------------------------------
      if(icheck.gt.0) then
         call getmem(0,l0)
         return
      endif
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
      nparkl=npar(ibl)
c---------------------------------------
      call dimenij(ibl,inx,nparkl,nbl2, iis,jjs,ijbl,kldim,klcont)
c---------------------------------------
      klpar1=nparkl
      call getmem(3*kldim,kabprim)
c---------------------------------------
         call ab_prim(ibl,nparkl,ijbl,nbl2,
     *                dbl_mb(inuc),dbl_mb(ibas),inx,iis,jjs,
     *                bl(kabprim),klcont,
     *                nij_uniqe,bl(ij_uniqe_p) )
c
c---------------------------------------
      call txs_second(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
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
c=======================
      subroutine ab_prim(ibl,nparij,ijbl,nbl2,
     *                   datnuc,datbas,inx,iis,jjs,
     *                   abprim,ijcont,
     *                   nij_uniqe,ij_uniqe)
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datbas(13,*),datnuc(5,*)
      dimension abprim(nparij,ijcont,3)
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
c
      dimension ij_uniqe(nij_uniqe)
ctest ?
      dimension map(2000) ! from uniqe to ijpar
c
c     data pi23 / 0.258012275487719489d0 /
c---------------------------------------------------------------
ctest
c     write(6,*)' from ab_prim : ibl=',ibl,' nij_uniqe=',nij_uniqe
c     call prt_uniqe(ij_uniqe,nij_uniqe)
c     write(6,*)' from ab_prim : end    nparij=',nparij
ctest
c
      do 50 ijpar=1,nparij
      ijcs=ijbl(ibl,ijpar)
      map(ijcs)=ijpar
   50 continue
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
      subroutine precalc2(isupb,bl,mmax,mmax1,nhabcd,nfumax, nbl2,nbls,
     *                       inx,iis,jjs,ijbl,nblok1,
     *                       ibl,nijbeg,nijend,npij,
     *                       kbl,nklbeg,nklend,npklx,npkl)
      implicit real*8 (a-h,o-z)
      character*3 reuse
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
     *                 icfg,jcfg,kcfg,lcfg
c
      dimension bl(*), inx(12,*),iis(*),jjs(*), ijbl(nbl2,*)
      dimension nblok1(2,*)
c---------------------------------------------------------------
      call reuseij(ibl,nbls,mmax,nblok1,reuse)
c     write(6,*)'-------------------------------------'
c     write(6,*)' from precalc2 reuse=',reuse,' nbls=',nbls
c     write(6,*)' ibl,kbl=',ibl,kbl,' ngcd=',ngcd
c     write(6,*)' address igci-l=',igci,igcj,igck,igcl
c     write(6,*)' npij,lcij=',npij,lcij,' npkl,lckl=',nplkx,lckl
c     write(6,*)' lci,lcj=',lci,lcj
c     write(6,*)' lck,lcl=',lck,lcl
c-----------------------------------------------------------
      ijbegin=nijbeg
      klbegin=nklbeg
c-----------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x(iabprim,ijdim,iapb,i1apb,isab)
c--------
      lpartot=lpartot+1
      if(reuse.eq.'no ') then
          lpareal=lpareal+1
          call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
          call precal2a(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npij,
     *                  ibl,ijbl,nbl2,nijbeg,nijend,
     *                  bl(iabprim      ),ijpar1 ,
     *                  lcij,bl(iaa),bl(ibb),bl(ieab),bl(icis),bl(icjs),
     *                  bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *                  bl(igci),bl(igcj),ngci1,ngcj1,'left ',
     *                nij_uniqe,bl(ij_uniqe_p) )
c----------
c for abnia
          if(mmax.gt.2) then
             call precal2b(mmax1,lcij,npij, bl(i1apb),
     *                     ijpar1,ijbegin, bl(iabnia))
          endif
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
          call precal2x(kabprim,kldim,icpd,i1cpd,iscd)
c
          lpartot=lpartot+1
          lpareal=lpareal+1
          call get_uniq(kbl,nkl_uniqe,kl_uniqe_p,bl)
          call precal2a(dbl_mb(ibas),dbl_mb(inuc),iis,jjs,inx, npkl,
     *                  kbl,ijbl,nbl2,nklbeg,nklend,
     *                  bl(kabprim      ),klpar1 ,
     *                  lckl,bl(icc),bl(idd),bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igck),bl(igcl),ngck1,ngcl1,'right',
     *                nkl_uniqe,bl(kl_uniqe_p) )
c--
c for cdnia
         if(mmax.gt.2) then
            call precal2b(mmax1,lckl,npklx,bl(i1cpd),
     *                    klpar1,klbegin, bl(icdnia))
         endif
      else
          call precdiag
      endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c(npij,npklx,npkl, bl(i1apb),bl(i1cpd),
     *                 ijpar1,klpar1,
     *                 bl(ihabcd),nhabcd,nfumax )
      endif
c---------------------------------------------------------------
      end
c============================================================
      subroutine index4(npij,npkl,npklx,nijbeg,nijend, nklbeg,nklend,
     *                  indxij,indxkl)
      dimension indxij(*),indxkl(*)
c-------
c   INPUT :
c
c   For a given block of contracted Quartets it is needed to 
c   consider here :
c  
c    NPIJ pairs <IJ  -  from NIJBEG to NIJEND in IJBL(ibl,nijbeg-nijend)
c    NPKL pairs KL>  -  from NKLBEG to NKLEND in IJBL(kbl,nklbeg-nklend)
c
c    NPKLX  - is 0 if this subroutine is called for a diagonal block
c                  (IJ,IJ), or it is NPKL otherwise.
c-------
c  OUTPUT :
c
c     INDXIJ(ijkl)  - gives the first  pair IJ for IJKL quartet
c     INDXKL(ijkl)  - gives the second pair KL for IJKL quartet
c
c-------
c
c* setup relation between ijkl and ijpar=1,..., and klpar=1,...
c
      ijkl=0
      ijpar=0
      do 10 ijp=nijbeg,nijend
      ijpar=ijpar+1
         klpar=0
         nklendx=nklend
         if(npklx.eq.0) nklendx=ijp
         do 10 klp=nklbeg,nklendx
         klpar=klpar+1
         ijkl=ijkl+1
         indxij(ijkl)=ijpar
         indxkl(ijkl)=klpar
   10 continue
c
c         write(8,88) (indxij(ii),ii=1,ijkl)
c         write(8,88) (indxkl(ii),ii=1,ijkl)
c 88      format('from index4 ;indxij,kl=',5i3)
c
      end
c====================================================================
      subroutine indexp(npij,npkl,npklx,nijbeg,nijend, nklbeg,nklend,
     *                  indxij,indxkl,ipres,indxp,nblsp)
      dimension indxij(*),indxkl(*),ipres(*),indxp(*)
c-------
c  setup relation between PRESENT c.s.quartets and pairs
c            ijklp and ijpar=1,..., and klpar=1,...
c
      ijklp=0
      ijkl=0
      ijpar=0
      do 10 ijp=nijbeg,nijend
      ijpar=ijpar+1
         klpar=0
         nklendx=nklend
         if(npklx.eq.0) nklendx=ijp
         do 10 klp=nklbeg,nklendx
         klpar=klpar+1
         ijkl=ijkl+1
             if(ipres(ijkl).ne.0) then
                ijklp=ijklp+1
                indxij(ijklp)=ijpar
                indxkl(ijklp)=klpar
                indxp(ijklp) =ijkl
             endif
   10 continue
c
      nblsp=ijklp
c
c         write(8,88) (indxij(ii),ii=1,nblsp)    
c         write(8,88) (indxkl(ii),ii=1,nblsp)   
c 88      format('from index4 ;indxij,kl=',5i3)
c         write(8,89) (indxp(ii),ii=1,nblsp)
c 89      format('from index4 ;indxp    =',5i3)
c
      end
c====================================================================
c********************************************************************
c**       Description of the PRECAL2A subroutine
c**  This subroutine is used to calculate some pair-quantities
c**  for left-hand as well as for right-hand pairs
c**
c**  This subroutine calculates 26 quantities for contracted
c**  shell pairs which belong to a given block of contracted
c**  quartets. There are NPIJ pairs <IJ and NPKL pairs KL>  
c**  which form a given block of contracted quartets. The 
c**  relationship between IJKL quartet and these pairs is
c**  through vectors INDXIJ and INDXKL :
c**    <IJ = INDXIJ( ijkl )  and      KL> = INDXKL( ijkl ) 
c
c  INPUT :
c
c     INX(12,ICS)     - informations about ICS contracted shell
c     DATNUC(5,iatom) - charge,position,.. for atom IATOM
c     DATBAS(13,ips)   - basis set data for IPS primitive shell
c    
c     LCIJ, LCKL      - lci*lcj  and lck*lcl where lci,lcj,lck,lcl are
c                       contraction lengths for ICS,JCS,KCS,LCS shells
c
c     NGCI1,NGCJ1,NGCK1,NGCL1 - depth of General Contraction
c                               for shells ICS,JCS,KCS,LCS
c     
c
c   Information about PAIRS of contracted shells :
c      ( comes from Blkprep routine )
c
c     NBL2  - number of blocks of contracted PAIRS 
c     IJBL(ibl,ij) - IJCS pair as the IJth pair from IBL block
c
c     IIS(ijcs) -> ICS - gives first (ICS) contracted shell in IJCS pair
c     JJS(ijcs) -> JCS - gives last  (JCS) contracted shell in IJCS pair
c
c   For a given block of contracted Quartets it is needed to 
c   consider here :
c  
c    NPIJ pairs <IJ  -  from NIJBEG to NIJEND in IJBL(ibl,nijbeg-nijend)
c    NPKL pairs KL>  -  from NKLBEG to NKLEND in IJBL(kbl,nklbeg-nklend)
c
c    NPKLX  - is 0 if this subroutine is called for a diagonal block
c                  (IJ,IJ), or it is NPKL otherwise.
c*************
c
c  OUTPUT :
c
c     INDXIJ(ijkl)  - gives the first  pair IJ for IJKL quartet
c     INDXKL(ijkl)  - gives the second pair KL for IJKL quartet
c
c     AAA, CIS(ij,ips) - exponent and re-scaled contracion coefficient
c                        of each IPS primitive shell from ICS contracted
c                        shell belonging to the pair IJ
c     BBB, CJS(ij,jps) - like above for JPS primitive shell from
c                        JCS contracted shell
c     CCC, CKS(kl,kps) - like above for KPS, KCS and pair KL
c     DDD, CLS(kl,lps) - like above for LPS, LCS 
c
c     XAB(ij,3) - coordinaties x,y,z  of (A-B)  for each pair IJ 
c     XCD(kl,3) - coordinaties x,y,z  of (C-D)  for each pair KL 
c
c* Matrices XPARIJ, XPARKL and YPARIJ, YPARKL 
c* have the following meaning :
c* (a,b,c,d- exponents,  A,B,C,D- orbital centers ,
c*     csi,cpi-  s- and p-contraction coefficients )
c
c  XPARIJ(npij,3,lcij,1) =  XP(npij,3,lcij)   -   ( a*A + b*B )/(a+b)
c  XPARIJ(npij,3,lcij,2) = XPN(npij,3,lcij)   -   XP - A  or XP -B
c  XPARIJ(npij,3,lcij,3) = XPP(npij,3,lcij)   -   a*A  +  b*B
c
c  XPARKL(npkl,3,lckl,1) =  XQ(npkl,3,lckl)   - like above for centers
c  XPARKL(npkl,3,lckl,2) = XQN(npkl,3,lckl)        C,D --> Q, QN, QQ
c  XPARKL(npkl,3,lckl,3) = XQQ(npkl,3,lckl)      
c
c  YPARIJ(npij,lcij,1)  =    APB(npij,lcij)   -      a+b
c  YPARIJ(npij,lcij,2)  =   RAPB(npij,lcij)   -    1/(a+b)
c  YPARIJ(npij,lcij,3)  = FACTIJ(npij,lcij)   - only when ICS and JCS
c                                                are L-shells :
c                			       csi/cpi * csj/cpj
c
c  YPARIJ(npij,lcij,4) = COEFIJ(npij,lcij)   - (csi or cpi)*(csj or cpj)
c                                                   *sab/(a+b)
c  YPARIJ(npij,lcij,5) =    SAB(npij,lcij)   - pure overlap sab
c
c  YPARKL(npkl,lckl,1)  =   CPD(npkl,lckl)   -      c+d         
c  YPARKL(npkl,lckl,2)  =  RCPD(npkl,lckl)   -   1/(c+d)
c  YPARKL(npkl,lckl,3)  = FACTKL(npkl,lckl)  - only when KCS and LCS
c                                                are L-shells :
c                			       csk/cpk * csl/cpl
c  YPARKL(npkl,lckl,4)  = COEFKL(npkl,lckl)  - (csk or cpk)*(csl or cpl)
c                                                 8*scd/(c+d)
c  YPARKL(npkl,lckl,5)  =    SCD(npkl,lckl)  -  scd *(2/pi)**3
c
c  TXAB(IJ,3,lcij) - coordinates XAB re-scaled by exponents a or b
c  TXCD(KL,3,lckl) - coordinates XAB re-scaled by exponents c or d
c                   Txab and Txcd are used in XWPQ routine to calculate
c                   P1234(ijkl,3) needed in Trobsa / ijkl-quartet /
c
c  GCI(npij,ngci1,lci)- general contraction coef. for ICS 
c  GCJ(npij,ngcj1,lcj)- general contraction coef. for JCS 
c  GCk(npkl,ngck1,lck)- general contraction coef. for KCS 
c  GCL(npkl,ngcl1,lcl)- general contraction coef. for LCS 
c
c*******************
c     NOTE :
c
c* Here SAB, SCD as well as Coefij, Coefkl are not "symmetrical" :
c*
c*      SAB=sab  without any factors (pure overlap)
c*      SCD=scd  multiplied by (2/pi)**3
c* These two above are used as a prodauct in Negelect only.
c*
c* Coefij and Coefkl as a product are used in Neglect only
c* to calculate Const. They are :
c* Coefij  is    sab * contraction ij  * 1/(a+b)
c* Coefkl  is  8*scd * contraction kl  * 1/(c+d)
c
c********************************************************************
      subroutine precal2a(datbas,datnuc,iis,jjs,inx,npij,
     * ibl,ijbl,nbl2,nijbeg,nijend,
     * abprim,ijpar1,lcij, aaa,bbb,estab,cis,cjs,
     * xab,xparij,coefij,factij,txab, gci,gcj,ngci1,ngcj1,which ,
     *                nij_uniqe,ij_uniqe )
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
      dimension ij_uniqe(nij_uniqe)
ctest ?
      dimension map1(2000),map2(2000) ! for unique pairs 
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
c mapping :
      ijpar=0
      do 50 ijp=nijbeg,nijend
      ijpar=ijpar+1
      ijcs=ijbl(ibl,ijp)
      map1(ijcs)=ijp
      map2(ijcs)=ijpar
   50 continue
c
c---------------
c
c* precalculations for the pairs IJ :
c
cccc  ijpar=0
cccc  do 100 ijp=nijbeg,nijend
cccc  ijpar=ijpar+1
cccc    ijcs=ijbl(ibl,ijp)
      do 100 ij_u=1,nij_uniqe
        ijcs=ij_uniqe(ij_u)
        ijp  =map1(ijcs)
        ijpar=map2(ijcs)
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
      subroutine precdiag
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------
      common /memor5x/ ieab,iecd
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5c/ itxab,itxcd,iabcd,ihabcd
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg
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
c--
      iecd=ieab
c
      icks=icis
      icls=icjs
c
      itxcd=itxab
c
      igck=igci
      igcl=igcj
c
      icdnia=iabnia
c-----------------------------------
      end
c====================================================================
      subroutine precal2b(mmax1,lcij,npij, rapb,ijpar1,ijbegin, abnia)
c-------------------------------------------------------------------
c  Description :
c
c  This subroutine calculates 1 quantity for contracted shell pairs 
c  which belong to a given block of contracted shell quartets .
c    This subroutine is called (from Blockint) only for blocks
c  with total angular momentum MMAX > 2 ,where MMAX is the sum 
c  of angular momentum of four contracted shells (+1) /for exemple :
c  (ss ss):mmax=1, (pp,ps):mmax=4, (dd,ll):mmax=7, (ff,ff):mmax=13 /
c
c  INPUT 
c  ------
c  MMAX1 = MMAX-1 , where MMAX is the sum of angular momentum 
c  NPIJ  - number of pairs <IJ in a given block of contracted quartets
c  RAPB(ij,lcij) - 1/(a+b) : reversed sum of exponents of each
c                     pair of primitive shells from IJth pair of  
c                     contracted shells.
c  OUTPUT
c  -------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  1. ABNIA(IJ,L,ij) -   L*( 0.5/(a+b) )  with L=1,2,...MMAX-1
c
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c---> dimension rapb(npij,*), abnia(npij,mmax1,*)
      dimension rapb(ijpar1,*), abnia(npij,mmax1,*)
      data half /0.5d0/
c--------------------------
      ijstart=ijbegin-1
c--------------------------
c
      do 150 ij=1,lcij
         do 200 ijpar=1,npij
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
      subroutine precal2c(npij,npkl,npklx, rapb,rcpd,ijpar1,klpar1,
     *                    habcd,nhabcd,nfumax)
c-----------------------------------------
c  Description of the PRECAL2b subroutine
c
c  This subroutine calculates 1 quantity for contracted shell pairs 
c  which belong to a given block of contracted shell quartets .
c    This subroutine is called (from Blockint) only for blocks
c  with total angular momentum MMAX > 2 ,where MMAX is the sum 
c  of angular momentum of four contracted shells (+1) /for exemple :
c  (ss ss):mmax=1, (pp,ps):mmax=4, (dd,ll):mmax=7, (ff,ff):mmax=13 /
c
c  INPUT 
c  ------
c
c    NPIJ, NPKL  - number of pairs <IJ  nd KL> in a given block
c                  of contracted quartets
c    NPKLX  - is 0 if this subroutine is called for a diagonal block
c                  (IJ,IJ), or it is NPKL otherwise.
c    RAPB(ij,lcij) - 1/(a+b) : reversed sum of exponents of each
c                     pair of primitive shells from IJth pair of  
c                     contracted shells.
c    RCPD(kl,lckl) - 1/(c+d)  : like above for pairs KL>
c
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
c*************
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
             do 190 klpar=1,npkl
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
             do 290 ijpar=1,npij
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
      subroutine prec4neg(nbls,npij,npkl,ndiag,ij,kl,
     1     ijpar1,lc12, klpar1,lc34,indxij,indxkl,
     2     estab,estcd,
     *     apb,cpd,coefij,coefkl,xp,xq, nsym,isymm,rnsym,
coutput
     *     rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c*************
c  Description of the PRECNEG4 subroutine
c
c  For two primitive pairs ij and kl which contribute 
c  to all NBLS contracted quartets IJKL in a given block
c  an estimation is made and some of quartets are neglected.
c  The block size NBLS is reduced to NBLS1. Indecies of remaining
c  quartets are stored in vector INDEX and for these the CONST vector
c  is constructed.
c
c  INPUT
c********
c  
c  NBLS - block size (number of quartets)
c  NPIJ, NPKL - number of contracted pairs forming a given block
c  IJ, KL     - primitive shell pairs
c  EPSR- thresold for neglect
c
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c
c  COEFIJ(ijpar,lcij) - re-scaled contraction coef.(see Precal2a)
c  COEFKL(klpar,lckl) - re-scaled contraction coef.(see Precal2a)
c
c
c  OUTPUT
c*********
c
c  NBLS1 - reduced block-size
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows 
c                 the ijkl  quartet in Original block
c  RHO(ijkl1)   =(a+b)*(c+d)/(a+b+c+d)
c  RPPQ(ijkl1)  = 1/(a+b+c+d)      where a,b,c,d, are exponents
c  CONST(ijkl1) = 8*rapb*rcpd*SQRT(rppq)
c                  *contri*contrj*contrk*contrl * sab*scd
c                 /sab,scd denote "pure" overlaps /
c*************
c
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cflops
      common /flops/ iflop(20)
      common /neglect/ eps,eps1,epsr
cnowy
      common /superbl/ isupb,ibl,kbl
      common /begin/ ijbegin,klbegin
c
      dimension indxij(*),indxkl(*),index(*)
c-------
      dimension apb(ijpar1,lc12),cpd(klpar1,lc34)
      dimension coefij(npij,lc12),coefkl(npkl,lc34)
      dimension  estab(npij,lc12), estcd(npkl,lc34)
c-------
      dimension xp(npij,3,*),xq(npkl,3,*)
      dimension        rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c symmetry
      dimension isymm(*)
c
      data zero /0.d0/
      data one /1.d0/
c---------------------------------------------------------------
c not par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
c
      sqrpold=one
      rpold=one
      abpcd1r=zero
      abpcdr=zero
c
      ijkl=0
      ijkl1=0
      ijklp=0
      do 100 ijpar=1,npij
         apb1=apb(ijstart+ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
         if(ndiag.eq.0) coef1=coef1*par268
         xp1= xp(ijpar,1,ij)
         xp2= xp(ijpar,2,ij)
         xp3= xp(ijpar,3,ij)
      npx=npkl
      if(ndiag.eq.0) npx=ijpar
      do 100 klpar=1,npx
      ijkl=ijkl+1
      ijklsm=isymm(ijkl)
      if(ijklsm.eq.0) go to 100
         if(nsym.gt.0)  symfac=rnsym*dble(ijklsm)
         ijklp=ijklp+1
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
c
      if(abpcd1.ne.abpcd1r) then
         abpcd1r=abpcd1
         abpcdr=one/abpcd1r
      endif
         rho1=abxcd*abpcdr
c
         estim=esti1*esti2
         estim=estim*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            ijkl1=ijkl1+1
c------->   index(ijkl1)=ijkl
            index(ijkl1)=ijklp
c
            rppq(ijkl1)=abpcdr
            rhoapb(ijkl1)=abpcdr*cpd1
            rhocpd(ijkl1)=abpcdr*apb1
c
            x1= xp1 - xq(klpar,1,kl)
            x2= xp2 - xq(klpar,2,kl)
            x3= xp3 - xq(klpar,3,kl)
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
            if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac
         endif
 100  continue
c
      nbls1=ijkl1
c
cflops
cxxxxx iflop(6)=iflop(6)+3*nbls+2*nbls1
c
      end
c====================================================================
      subroutine precspec(nbls,npij,npkl,ndiag, ij,kl,
     1     ijpar1,lc12, klpar1,lc34, indxij,indxkl,
     2     estab,estcd,
     *     apb,cpd,coefij,coefkl,xp,xq, rapb,rcpd,txab,txcd,
     *                   nsym,isymm,rnsym,
c output
     *     rys,const,xpqr,txxr,nbls1,index)
c*************
c*
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cflops
      common /flops/ iflop(20)
      common /neglect/ eps,eps1,epsr
cnowy
      common /begin/ ijbegin,klbegin
      common /superbl/ isupb,ibl,kbl
c
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
      dimension isymm(*)
c
      data zero /0.d0/
      data one /1.d0/
c---------------------------------------------------------------
c not par268=eps1*8.d0
      par268=     8.d0
c---------------------------------------------------------------
      ijstart=ijbegin-1
      klstart=klbegin-1
c---------------------------------------------------------------
c
      sqrpold=one
      rpold=one
c
      abpcd1r=zero
      abpcdr=zero
c
      ijkl=0
      ijkl1=0
      ijklp=0
      do 100 ijpar=1,npij
         apb1=apb(ijstart+ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         if(ndiag.eq.0) coef1=coef1*par268
         xp1= xp(ijpar,1,ij)
         xp2= xp(ijpar,2,ij)
         xp3= xp(ijpar,3,ij)
      npx=npkl
      if(ndiag.eq.0) npx=ijpar
      do 100 klpar=1,npx
      ijkl=ijkl+1
c
      ijklsm=isymm(ijkl)
      if(ijklsm.eq.0) go to 100
         if(nsym.gt.0) symfac=rnsym*dble(ijklsm)
c-----------
         ijklp=ijklp+1
c
         cpd1=cpd(klstart+klpar,kl)
         coef2=coefkl(klpar,kl)
         esti2=estcd (klpar,kl)
         abpcd1=apb1+cpd1
         abxcd=apb1*cpd1
cnew
c
      if(abpcd1.ne.abpcd1r) then
         abpcd1r=abpcd1
         abpcdr=one/abpcd1r
      endif
         rho1=abxcd*abpcdr
c
         estim=esti1*esti2
         estim=estim*abpcdr
         if(estim.gt.epsr) then
            coef12=coef1*coef2
            ijkl1=ijkl1+1
c---->      index(ijkl1)=ijkl
            index(ijkl1)=ijklp
c
            x1= xp1 - xq(klpar,1,kl)
            x2= xp2 - xq(klpar,2,kl)
            x3= xp3 - xq(klpar,3,kl)
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
            if(nsym.gt.0) const(ijkl1)=const(ijkl1)*symfac
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
cflops
cxxxxx iflop(6)=iflop(6)+3*nbls+2*nbls1
c
      end
c====================================================================
      subroutine xwpq(nbls1,xwp,xwq,p1234,
     1     ijpar1,lc12, klpar1,lc34,
     *        lcij,lckl,npij,npkl,
     *        indxij,indxkl,index,
     *        rppq,xp,xq,xpp,xqq,
     *        txab,txcd,abcd,apb,rcpd,cpd,rapb)
c-----------------------------------------------
c  Description of the XWPQ     subroutine
c
c  For two primitive pairs lcij and lckl which contribute 
c  to all NBLS1 contracted quartets IJKL in a given block
c  4 quatities are calculated for each quartet (xwp,xwq,p1234,abcd)
c
c  INPUT :
c*******
c 
c  NBLS1 - blocksize
c  lcij, lckl - pairs of primitive shells
c  NPIJ,NPKL  - number of contracted pairs
c  INDEX(ijkl1) - for ijkl1 quartet in Reduced block shows 
c                 the ijkl  quartet in Original block
c  INDXIJ(ijkl) - shows the contracted pair IJPAR for ijkl quartet
c  INDXKL(ijkl) - shows the contracted pair KLPAR for ijkl quartet
c  RPPQ(ijkl1)  = 1/(a+b+c+d)         where a,b,c,d, are exponents
c  XP(ijpar,3,lcij) - coordinates of P
c  XQ(klpar,3,lckl) - coordinates of Q
c  XPP(ijpar,3,lcij)- coordinates of ( a*A + b*B)
c  XQQ(klpar,3,lckl)- coordinates of ( c*C + d*D)
c  TXAB(ijpar,3,  )  - coordinates XAB re-scaled by exponents a or b
c  TXCD(klpar,3,  )  - coordinates XAB re-scaled by exponents c or d
c  APB(ijpar,lcij) =a+b
c  CPD(klpar,lckl) =c+d
c  RAPB(ijpar,lcij)=1/(a+b)
c  RCPD(klpar,lckl)=1/(c+d)
c
c  OUTPUT
c********
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
      subroutine specase(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                   npij,npkl,ii,jj,kk,ll,
     *                   cis,cjs,cks,cls,
     *                   buf,buf1, const,rysx,xpqr,txxr,concoe)
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
      subroutine precal2x(iabprim,ijdim,
     *                  iapb,i1apb,isab)
c-------------------------------------------------------------
c* precalculations for the pairs IJ :
c
      iapb =iabprim
      i1apb=iabprim   + ijdim
      isab =iabprim   + ijdim*2
c--------------
      end
c====================================================================
      subroutine reuseij(ibl,nbls,mmax,nblok1,reuse)
      character*3 reuse
      logical cond1,cond2
      common /ijcsfl/ ijblokp,ijprevf,ijprevl,ijtprev,maxprev,ngcprev
      common /types/itype,jtype,ktype,ltype,itype1,jtype1,ktype1,ltype1 
      common /gcont/ ngci1,ngcj1,ngck1,ngcl1,ngcd
      dimension nblok1(2,*)
c-----------------------------------------------------------
      reuse='no '
c---------------------------------------------------------------
cnew in order to remember left-hand pairs IJ in a block ikbl
c keep it in common /ijcsfl/ and use it 
c
      ijcurrf=nblok1(1,1)
      ijcurrl=nblok1(1,nbls)
      ijtcurr=0
      if(mmax.eq.2) then
        if(itype1.gt.1 .or. jtype1.gt.1) ijtcurr=2
      endif
      maxcurr=mmax
c     ngccurr=ngcd
      ngccurr=1000*ngci1+100*ngcj1+10*ngck1+ngcl1
c----
      if(ibl.ne.ijblokp) go to 1000
c----------------------------------------
      cond1=.false.
      cond2=.false.
c-----
        if(ijcurrf.eq.ijprevf .and. ijcurrl.eq.ijprevl) cond1=.true.
ctest.............
c       if(ngccurr.eq.ngcprev. or. 
c    *  (ngccurr.gt.1.and.ngcprev.gt.1)) cond2=.true.
ctest.............
        if(ngccurr.eq.ngcprev ) cond2=.true.
c----------------------------------------
c  check if the current block has the same left-hand pairs IJ
c  like the previous one. If it has DO NOT precalculate 
c  corresponding IJ-pair quantities. Instead, re-use old ones 
c  ( the memory location is the same )
c
c ijtprev and ijtcurr may have value of 0 or 2. They are needed only
c for the mmax=2 and give information about re-scaling (2) or not (0)
c precalculated TXAB 
c       
      if(cond1 .and. cond2) then
c------>  if(maxprev.ne.2) go to 9999
          if(maxcurr.eq.maxprev) go to 9999
          if(ijtprev.eq.2 .and. ijtcurr.ne.2) then
          else
            go to 9999
          endif
      endif
          go to 8888
 9999 continue
      reuse='yes'
 8888 continue
c
c----------------------------------------
cnew in order to remember left-hand pairs IJ in a block ikbl
c keep it in common /ijcsfl/ and use it in precalc...
c
 1000 continue
      ijblokp=ibl
      ijprevf=ijcurrf
      ijprevl=ijcurrl
      ijtprev=ijtcurr
      maxprev=maxcurr
      ngcprev=ngccurr
c----------------------------------------
      end
c================================================================
