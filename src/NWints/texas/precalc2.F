c--------------------------------------------------------------------
      subroutine prec2ij(ibl,nblock1, bl,inx)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /route/ iroute
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primij/ iabprim, ijdim ,ijpar1 
      common /time0/ tprec2
c
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
      dimension inx(12,*)
      dimension nblock1(*)
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
c get : nij_uniqe=no of unique pairs and the pointer to the list
c
c10   call get_uniq(ibl,nij_uniqe,ij_uniqe_p,bl)
c
c nij_uniqe = number of uniqe pairs in ibl block of pairs
c (no of present pairs)
c ij_uniqe_p = pointer to the list
c---------------------------------------
      call dimenij(ibl,inx,nij_uniqe,nblock1,ijdim,ijcont)
c---------------------------------------
      IF( iroute.eq.1 ) THEN
         call getmem(3*ijdim,iabprim)
         call ab_prim_1(ibl, dbl_mb(inuc),dbl_mb(ibas),inx,
     *                  bl(iabprim),ijcont,
     *                  nij_uniqe,bl(ij_uniqe_p) )
      ELSE
         call getmem(2*ijcont+ijdim,iabprim)
         iapb =iabprim
         i1apb=iabprim+ijcont
         isab =iabprim+ijcont*2
         call ab_prim_2(ibl, dbl_mb(inuc),dbl_mb(ibas),inx,
     *                  bl(iapb),bl(i1apb),bl(isab),ijcont,
     *                  nij_uniqe,bl(ij_uniqe_p) ) 
      ENDIF
c---------------------------------------
      call txs_second(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c=======================
      subroutine prec2kl(kbl,nblock1, bl,inx)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
      logical firstd
      common /route/ iroute
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /primkl/ kabprim, kldim ,klpar1 
      common /time0/ tprec2
c
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*)
      dimension inx(12,*)
      dimension nblock1(*)
c---------------------------------------
      call txs_second(tprec2b)
c---------------------------------------
c10   call get_uniq(kbl,nkl_uniqe,kl_uniqe_p,bl)
c---------------------------------------
      call dimenij(kbl,inx,nkl_uniqe,nblock1,kldim,klcont)
c---------------------------------------
      IF( iroute.eq.1 ) THEN
         call getmem(3*kldim,kabprim)
         call ab_prim_1(kbl, dbl_mb(inuc),dbl_mb(ibas),inx,
     *                  bl(kabprim),klcont,
     *                  nkl_uniqe,bl(kl_uniqe_p) )
      ELSE
         call getmem(2*klcont+kldim,kabprim)
         icpd =kabprim
         i1cpd=kabprim+klcont
         iscd =kabprim+klcont*2
         call ab_prim_2(kbl,dbl_mb(inuc),dbl_mb(ibas),inx,
     *                  bl(icpd),bl(i1cpd),bl(iscd),klcont,
     *                  nkl_uniqe,bl(kl_uniqe_p) )
      ENDIF
c---------------------------------------
      call txs_second(tprec2e)
      tprec2=tprec2+tprec2e-tprec2b
c---------------------------------------
      end
c====================================================================
      subroutine dimenij(ibl,inx,nparij,nblock1,ijdim,ijcont)
      dimension inx(12,*),nblock1(*)   ! (0:ncs)
c
      call get_ics_jcs(nblock1,ibl,1,ics1,jcs1)
c
      icont=inx(5,ics1)-inx(1,ics1)
      jcont=inx(5,jcs1)-inx(1,jcs1)
      ijcont=icont*jcont
      ijdim=nparij*ijcont
c
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
      subroutine ab_prim_1(ibl, datnuc,datbas,inx, abprim,ijcont,
     *                     nij_uniqe,ij_uniqe )
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension inx(12,*)
      dimension datbas(13,*),datnuc(5,*)
      dimension abprim(nij_uniqe,ijcont,3)
      dimension ij_uniqe(nij_uniqe) 
c---------------------------------------------------------------
ctest
c     write(6,*)' from ab_prim : ibl=',ibl,' nij_uniqe=',nij_uniqe
c     call prt_uniqe(ij_uniqe,nij_uniqe)
c     call prt_uniqe(ibl,ij_uniqe,nij_uniqe,1111     )
c     write(6,*)' from ab_prim : ij1=',ij1,' ij2=',ij2   
ctest
c
cccc  do 100 ijpar=1,nparij
cccc  ijcs=ijbl(ibl,ijpar)
cnew 
c
      do 100 ijpar =1,nij_uniqe   ! over uniqe pairs 
c
      ijcs=ij_uniqe(ijpar)
      call get_ij_half(ijcs,ics,jcs)
c
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iat.eq.0) iat=jat
        if(jat.eq.0) jat=iat
c
        if(iat.eq.0 .and. jat.eq.0) then
           xab=zero
           yab=zero
           zab=zero
        else
           xab=datnuc(2,iat)-datnuc(2,jat)
           yab=datnuc(3,iat)-datnuc(3,jat)
           zab=datnuc(4,iat)-datnuc(4,jat)
        endif
c
        rr=xab*xab+yab*yab+zab*zab
c
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
cold            axb=aa*bb
cold            apb=aa+bb
cold            apb1=one/apb
cold            e=axb*apb1
cold            abprim(ijpar,ij,1)=apb
cold            abprim(ijpar,ij,2)=apb1
cold            abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp(-e*rr)
cnew98.......................................................
c
                axb=aa*bb
                apb=aa+bb
                apb1=one/apb
                e=axb*apb1
c98
                err=e*rr
                if(err.gt.32) then
                   exp_err=0.d0         ! exp(-32)=1.27*10**-14
                else
                   exp_err=exp(-err)
                endif
c98
                abprim(ijpar,ij,1)=apb
                abprim(ijpar,ij,2)=apb1
                abprim(ijpar,ij,3)=apb1*sqrt(sqrt(axb))**3*exp_err
cnew98.......................................................
             else
c center23:
                if(aa.le.zero .and. bb.gt.zero) then
c>>>>>>>>>         write(8,*)' exp_a=zero ; ics=',ics
                   apb=aa+bb
                   apb1=one/apb
                   abprim(ijpar,ij,1)=apb
                   abprim(ijpar,ij,2)=apb1
c>>*a-3/4          abprim(ijpar,ij,3)=     sqrt(bb)             !Sbb
                   abprim(ijpar,ij,3)=one/sqrt(sqrt(bb))   !sqrt(Sbb)
                endif
                if(aa.gt.zero .and. bb.le.zero) then
c>>>>>>>>          write(8,*)' exp_b=zero ; jcs=',jcs
                   apb=aa+bb
                   apb1=one/apb
                   abprim(ijpar,ij,1)=apb
                   abprim(ijpar,ij,2)=apb1
c>>>>>>>>>>>>      abprim(ijpar,ij,3)=sqrt(sqrt(aa)) ! sqrt(Saa)
c>>*a-3/4          abprim(ijpar,ij,3)=     sqrt(aa)             !Saa
                   abprim(ijpar,ij,3)=one/sqrt(sqrt(aa))   !sqrt(Saa)
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
      subroutine precalc2_1(isupb,bl,mmax,mmax1,nhabcd,nfumax,
     *                      nbl2,nbls, inx,ibl,kbl,npkl) 
c
c  npkl=ndiag 
c
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c-----------------------------------------------------------
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpar_exe,lpar_num
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
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*), inx(12,*)
c---------------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c--------
      lpar_exe=lpar_exe+1
      lpar_num=lpar_num+nij_uniqe
      call precal2a_1(dbl_mb(ibas),dbl_mb(inuc),inx, ibl, nbl2,
     *              bl(iabprim      ),
     *              lcij,bl(iaa),bl(ibb),bl(ieab),bl(icis),bl(icjs),
     *              bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *              bl(igci),bl(igcj),ngci1,ngcj1,'left ',
     *              nij_uniqe,bl(ij_uniqe_p) )
c----------
c for abnia
      if(mmax.gt.2) then
         call precal2b_1(mmax1,lcij, bl(i1apb), bl(iabnia),
     *                   nij_uniqe,bl(ij_uniqe_p) )
      endif
c----------------------------------------
c precalculations for pairs kl
c
c10   if(kbl.eq.ibl) then
c10       kabprim=iabprim
c10       kldim  =ijdim
c10   endif
c10   if(npkl.ne.0) then
          call precal2x_1(kabprim,kldim,icpd,i1cpd,iscd)
c
          lpar_exe=lpar_exe+1
          lpar_num=lpar_num+nkl_uniqe
          call precal2a_1(dbl_mb(ibas),dbl_mb(inuc),inx,kbl,nbl2,
     *                  bl(kabprim      ),
     *                  lckl,bl(icc),bl(idd),bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igck),bl(igcl),ngck1,ngcl1,'right',
     *                  nkl_uniqe,bl(kl_uniqe_p) )
c--
c for cdnia
         if(mmax.gt.2) then
            call precal2b_1(mmax1,lckl,bl(i1cpd), bl(icdnia),
     *                      nkl_uniqe,bl(kl_uniqe_p))
         endif
c10   else
c10       call precdiag
c10   endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_1(bl(i1apb),bl(i1cpd),bl(ihabcd),nhabcd,nfumax,
     *               nij_uniqe,bl(ij_uniqe_p),nkl_uniqe,bl(kl_uniqe_p))
      endif
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_1(datbas,datnuc,inx, ibl,nbl2, abprim,lcij,
     *                      aaa,bbb,estab,cis,cjs,
     *                      xab,xparij,coefij,factij,txab, 
     *                      gci,gcj,ngci1,ngcj1,which ,
     *                      nij_uniqe,ij_uniqe )
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
c
      dimension inx(12,*)
      dimension abprim(nij_uniqe,lcij,3)
      dimension datbas(13,*),datnuc(5,*)
c
      dimension aaa(nij_uniqe,*),bbb(nij_uniqe,*)
      dimension cis(nij_uniqe,*),cjs(nij_uniqe,*)
      dimension xab(nij_uniqe,3), xparij(nij_uniqe,3,lcij,3)
      dimension estab(nij_uniqe,lcij)
      dimension coefij(nij_uniqe,lcij), factij(nij_uniqe,lcij)
      dimension txab(nij_uniqe,3,*)
      dimension gci(nij_uniqe,ngci1,*),gcj(nij_uniqe,ngcj1,*)
c
      dimension ij_uniqe(nij_uniqe)
c
      dimension xa(3),xb(3)
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
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c-------------------------------------------------------------
c precalculations for the pairs IJ :
c
c-------------------------------old
cccc  ijpar=0
cccc  do 100 ijp=nijbeg,nijend
cccc  ijpar=ijpar+1
cccc    ijcs=ijbl(ibl,ijp)
c-------------------------------old
c
      do 100 ijpar=1,nij_uniqe
        ijcs=ij_uniqe(ijpar)
        call get_ij_half(ijcs,ics,jcs)
c
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
c
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
c gen.contr. shell is somewhere
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
c gen.contr.
                  coefij(ijpar,ji )=one
                  est_js=one
                  if(which.eq.'right') coefij(ijpar,ji)=par268
c
                  do 220 jg=0,ngcjj
                  gcj(ijpar,jg+1,js1)=datbas(jg+2,js)
  220             continue
               endif
c
            sab=   abprim(ijpar,ji,3)
            aa1=aa*abprim(ijpar,ji,2)
            bb1=bb*abprim(ijpar,ji,2)
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
         do 151 ijpar=1,nij_uniqe
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
         do 152 ijpar=1,nij_uniqe
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
           do 153 ijpar=1,nij_uniqe
              rapb1=abprim(ijpar,ijs1,2)
              txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
              txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
              txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c------------------------------------------------------------------
c
      end
c====================================================================
      subroutine precal2b_1(mmax1,lcij,rapb,abnia, nij_uniqe,ij_uniqe )
c-------------------------------------------------------------------
c  OUTPUT
c  -------
c  For each pair of primitive shells ij (kl) from IJ (KL) contracted
c  pair belonging to a given block of quartets with total ang.mom. MMAX
c
c  1. ABNIA(IJ,L,ij) -   L*( 0.5/(a+b) )  with L=1,2,...MMAX-1
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension rapb(nij_uniqe,*),abnia(nij_uniqe,mmax1,*)
      dimension ij_uniqe(nij_uniqe)
      data half /0.5d0/
c--------------------------
      do 150 ij=1,lcij
         do 200 ijpar=1,nij_uniqe
            apb2=half*rapb(ijpar,ij)
            abnia(ijpar,1,ij)=apb2
               do 250 i=2,mmax1
                  abnia(ijpar,i,ij)=abnia(ijpar,i-1,ij)+apb2
  250          continue
  200    continue
  150  continue
c
      end
c====================================================================
      subroutine precal2c_1(rapb,rcpd, habcd,nhabcd,nfumax, 
     *                      nij_uniqe,ij_uniqe, nkl_uniqe,kl_uniqe)
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
c  HABCD is used only in TRACY's recursive in routines TRACIJ, TRACKL
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /contr/ ngci,ngcj,ngck,ngcl,lci,lcj,lck,lcl,lcij,lckl
      common /logicd/ hnia(3,1)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension rapb(nij_uniqe,*),rcpd(nkl_uniqe,*)
      dimension habcd(nhabcd,3,nfumax,*)
      dimension ij_uniqe(nij_uniqe),kl_uniqe(nkl_uniqe)
c
      data half /0.5d0/
c-----------------------
         do 170 kl=1,lckl
           do 180 ifu=1,nfumax
           hx=hnia(1,ifu)
           hy=hnia(2,ifu)
           hz=hnia(3,ifu)
cccccccc     do 190 klpar=1,npkl
             do 190 klpar=1,nkl_uniqe
             habcd(klpar,1,ifu,kl)=hx*rcpd(klpar,kl)
             habcd(klpar,2,ifu,kl)=hy*rcpd(klpar,kl)
             habcd(klpar,3,ifu,kl)=hz*rcpd(klpar,kl)
  190        continue
  180      continue
  170    continue
c
      end
c====================================================================
      subroutine prec4neg_1(isbl_size,isbl_point,isbl_q,ndiag,
     *                      ij,kl,lc12,lc34,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq, 
     *                      map_ij, map_kl,
     *                      nij_uniqe,nkl_uniqe,
coutput
     *                      rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
c-------
      dimension map_ij(isbl_size),map_kl(isbl_size)
      dimension isbl_q(*)
      dimension index(*)
c-------
      dimension apb(nij_uniqe,lc12),cpd(nkl_uniqe,lc34)
      dimension coefij(nij_uniqe,lc12),coefkl(nkl_uniqe,lc34)
      dimension  estab(nij_uniqe,lc12), estcd(nkl_uniqe,lc34)
      dimension xp(nij_uniqe,3,*),xq(nkl_uniqe,3,*)
c----output
      dimension rppq(*),rhoapb(*),rhocpd(*),rys(*),const(*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
c not par268=eps1*8.d0
      par268=     8.d0
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
c-----------old-------------
c        ijcsq=ijpres2(iq)
c        klcsq=klpres2(iq)
c        ijpar=map_ij(ijcsq)
c        klpar=map_kl(klcsq)
c-----------old-------------
cnew: map to uniqe pairs 
c
         ijpar=map_ij(iqp)
         klpar=map_kl(iqp)
c
         apb1=apb(ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         cpd1=cpd(klpar,kl)
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
c10         if(ndiag.eq.0) coef12=coef12*par268
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
      subroutine precspec_1(isbl_size,isbl_point, isbl_q, ndiag,
     *                      ij,kl,lc12,lc34, indxij,indxkl,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq,
     *                      rapb,rcpd,txab,txcd,
     *                      map_ij,map_kl,
     *                      nij_uniqe,nkl_uniqe,
c output
     *                      rys,const,xpqr,txxr,nbls1,index)
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
c---------------
      dimension map_ij(isbl_size),map_kl(isbl_size)
      dimension isbl_q(*)
      dimension indxij(*),indxkl(*),index(*)
c---------------
      dimension apb(nij_uniqe,lc12),cpd(nkl_uniqe,lc34)
      dimension rapb(nij_uniqe,*),rcpd(nkl_uniqe,*)
      dimension coefij(nij_uniqe,lc12),coefkl(nkl_uniqe,lc34)
      dimension  estab(nij_uniqe,lc12), estcd(nkl_uniqe,lc34)
      dimension txab(nij_uniqe,3,*),txcd(nkl_uniqe,3,*)
      dimension xp(nij_uniqe,3,*),xq(nkl_uniqe,3,*)
c
      dimension rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
      par268=     8.d0
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
c-----------------old-------------
c        ijcsq=ijpres2(iq)
c        klcsq=klpres2(iq)
c        ijpar=map_ij(ijcsq)
c        klpar=map_kl(klcsq)
c-----------------old-------------
cnew : map to uniqe pairs
c
         ijpar=map_ij(iqp)
         klpar=map_kl(iqp)
c
         apb1=apb(ijpar,ij)
         coef1=coefij(ijpar,ij)
         esti1=estab (ijpar,ij)
c
         cpd1=cpd(klpar,kl)
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
c10         if(ndiag.eq.0) coef12=coef12*par268
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
        ijpar=indxij(ijkl)                   ! uniqe pair 
c
        xpqr(1,i)=xpqr(1,i)*rapb(ijpar,ij)
        xpqr(2,i)=xpqr(2,i)*rapb(ijpar,ij)
        xpqr(3,i)=xpqr(3,i)*rapb(ijpar,ij)
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
        xpqr(1,i)=xpqr(1,i)*rcpd(klpar,kl)
        xpqr(2,i)=xpqr(2,i)*rcpd(klpar,kl)
        xpqr(3,i)=xpqr(3,i)*rcpd(klpar,kl)
c
        txxr(1,i)=txcd(klpar,1,kl)
        txxr(2,i)=txcd(klpar,2,kl)
        txxr(3,i)=txcd(klpar,3,kl)
  220   continue
      endif
c
      end
c====================================================================
      subroutine xwpq_1(nbls1,xwp,xwq,p1234,lc12,lc34, lcij,lckl,
     *                  indxij,indxkl,index, nij_uniqe,nkl_uniqe,
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
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
c
      dimension indxij(*),indxkl(*),index(*)
c
      dimension  xp(nij_uniqe,3,*), xq(nkl_uniqe,3,*)
      dimension xpp(nij_uniqe,3,*),xqq(nkl_uniqe,3,*)
      dimension apb(nij_uniqe,*),rapb(nij_uniqe,*)
      dimension cpd(nkl_uniqe,*),rcpd(nkl_uniqe,*)
      dimension txab(nij_uniqe,3,*),txcd(nkl_uniqe,3,*)
      dimension rppq(*)
      dimension abcd(nbls1)
c
      logical stable 
      common /tracy_stability/ stable(1)
c---------------------------------------------------------------
         do 100 i=1,nbls1
c
            ijkl=index(i)
            ijpar=indxij(ijkl)          ! uniqe pair now !!!!
            klpar=indxkl(ijkl)          ! uniqe pair now !!!!
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
c for Tracy's recursive :
c
            rcpd1=rcpd(klpar,lckl)
c
            p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
            p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
            p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
c
  100    continue
c
c---------------------------------------------------------------
c because of possible numerical instablity in Tracy's recursive 
c it is needed to check it out (see xwpq_2) 
c
c FOR now IT IS NOT DONE yet for route=1 . It will go as stable
c if nsij.GE.nskl and as unstable for nsij.LT.nskl 
c (just like it used to be) 
c---------------------------------------------------------------
c old fashion way regardless of stability or instability (to be changed)
c
         if(nsij.ge.nskl) then
            stable(1)=.true.
            do i=1,nbls1
               ijkl=index(i)
               ijpar=indxij(ijkl)      
               klpar=indxkl(ijkl)     
               abcd(i)=apb(ijpar,lcij)*rcpd(klpar,lckl)
            enddo
         else
            stable(1)=.false.
            do i=1,nbls1
               ijkl=index(i)
               ijpar=indxij(ijkl)    
               klpar=indxkl(ijkl)   
               abcd(i)=rapb(ijpar,lcij)*cpd(klpar,lckl)
            enddo
         endif
c---------------------------------------------------------------
c
      end
c====================================================================
      subroutine specase_1(bl,first,nbls,nbls1, index,indxij,indxkl,
     *                     npij,npkl,ii,jj,kk,ll,
     *                     cis,cjs,cks,cls,
     *                     buf,buf1, const,rysx,xpqr,txxr,concoe)
c*
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
c-------------------------------------------------------------------
c
      dimension bl(*)
      dimension index(*),indxij(*),indxkl(*)
      dimension xpqr(3,*),txxr(3,*)
      dimension cis(npij,*),cjs(npij,*),cks(npkl,*),cls(npkl,*)
      dimension const(*),rysx(*),concoe(*)
      dimension buf(nbls,*),buf1(nbls1,4)
c
c-------------------------------------------------------------------
c  this routine constitues the special code for
c  two types of integrals over nbls quartets of primitive shells
c  1. (ss|ss)
c  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c  these integrals are also contracted here.
c
c
c  input
c  ------
c  all precalculated values for whole block :
c
c  const(nbls) - contains consts=pi3sabcd/(pqsqrt(ppq)) for all int.
c  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c  xp,xp      - geometry for p,q
c
c  output
c  ------
c  buf(nbls,) -contains final contracted integrals
c-------------------------------------------------------------------
c
c memory for f00,f11 :
c
      call getmem(nbls1,if00)
      call getmem(nbls1,if11)
c
      if00=if00-1
      if11=if11-1
c
c
      do 100 i=1,nbls1
      xrys=rysx(i)
      call ft0
      bl(if00+i)=f0
      bl(if11+i)=f1
  100 continue
c
c special code for (ss ss) integrals
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
c
        if( max(ityp,jtyp,ktyp,ltyp).eq.3) then
          if(ityp.eq.3) then
            do 1051 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)                 ! uniqe pair now !!
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
c
           do 105 i=1,nbls1
           buf1(i,1)=concoe(i)*bl(if00+i)
  105      continue
           intct=1
        endif
c
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
c relise memory
c
      call retmem(2)
c
      end
c====================================================================
      subroutine precal2x_1(iabprim,ijdim,iapb,i1apb,isab)
c-------------------------------------------------------------
c precalculations for the pairs IJ :
c
      iapb =iabprim
      i1apb=iabprim   + ijdim
      isab =iabprim   + ijdim*2
c
      end
c====================================================================
      subroutine ab_prim_2(ibl, datnuc,datbas,inx,apb,rapb,sab,ijcont,
     *                     nij_uniqe,ij_uniqe )
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      dimension inx(12,*)
      dimension datbas(13,*),datnuc(5,*)
cccc  dimension abprim(nparij,ijcont,3)
      dimension apb(ijcont),rapb(ijcont),sab(nij_uniqe,ijcont)
      dimension ij_uniqe(nij_uniqe)
c
      dimension sqrtx(400),eexx(400)
c---------------------------------------------------------------
cccc  do 100 ijpar=1,nparij
cccc  ijcs=ijbl(ibl,ijpar)
cnew
      do 100 ijpar =1,nij_uniqe   ! over uniqe pairs
c
      ijcs=ij_uniqe(ijpar)
      call get_ij_half(ijcs,ics,jcs)
c
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      iat=inx(2,ics)
        ja=inx(1,jcs)+1
        je=inx(5,jcs)
        jat=inx(2,jcs)
c
c center23 (atom number 0 is the one with s-orbital/zero-exp.)
c
        if(iat.eq.0) iat=jat
        if(jat.eq.0) jat=iat
c
        if(iat.eq.0 .and. jat.eq.0) then
           xab=zero
           yab=zero
           zab=zero
        else
           xab=datnuc(2,iat)-datnuc(2,jat)
           yab=datnuc(3,iat)-datnuc(3,jat)
           zab=datnuc(4,iat)-datnuc(4,jat)
        endif
c
        rr=xab*xab+yab*yab+zab*zab
c
           ij=0
           do 200 is=ia,ie
           aa=datbas(1,is)
c----------------------------------
             do 200 js=ja,je
             bb=datbas(1,js)
             ij=ij+1
c----------------------------------
             if(aa.gt.zero .and. bb.gt.zero) then
                if(ijpar.eq.1) then
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
cold            sab(ijpar,ij)=rapb(ij)*sqrt3*exp(-e*rr)
cnew98......
                err=e*rr
                if(err.gt.32.d0) then
                   sab(ijpar,ij)=0.d0
                else
                   sab(ijpar,ij)=rapb(ij)*sqrt3*exp(-err)
                endif
cnew98......
             else
c center23:
                if(aa.le.zero .and. bb.gt.zero) then
c>>>>>>>>>         write(8,*)' exp_a=zero ; ics=',ics
                   if(ijpar.eq.1) then
                      apb(ij)=aa+bb
                      rapb(ij)=one/apb(ij)
c>>*b-3/4             sqrtx(ij)=     sqrt(bb)       ! Sbb 
                      sqrtx(ij)=one/sqrt(sqrt(bb))  ! sqrt(Sbb)
                      sqrt3=sqrtx(ij)
                   else
                      sqrt3=sqrtx(ij)
                   endif
                   sab(ijpar,ij)=sqrt3               ! Sbb
                endif
                if(aa.gt.zero .and. bb.le.zero) then
c>>>>>>>           write(8,*)' exp_b=zero ; jcs=',jcs
                   if(ijpar.eq.1) then
                      apb(ij)=aa+bb
                      rapb(ij)=one/apb(ij)
c>>>*a-3/4            sqrtx(ij)=     sqrt(aa)       ! Saa 
                      sqrtx(ij)=one/sqrt(sqrt(aa))  ! sqrt(Saa)
                      sqrt3=sqrtx(ij)
                   else
                      sqrt3=sqrtx(ij)
                   endif
                   sab(ijpar,ij)=sqrt3              ! Saa
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
      subroutine precalc2_2(isupb,bl,mmax,mmax1,nfumax,
     *                      nbl2,nbls, inx,ibl,kbl,npkl)
      implicit real*8 (a-h,o-z)
#include "mafdecls.fh"
c-----------------------------------------------------------
      common /primij/ iabprim, ijdim ,ijpar1
      common /primkl/ kabprim, kldim ,klpar1
      common /ilepar/ lpar_exe,lpar_num
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
      common /map4uniq/ nij_uniqe, ij_uniqe_p, map_ij,
     *                  nkl_uniqe, kl_uniqe_p, map_kl
c
      dimension bl(*), inx(12,*)
c---------------------------------------------------------------
c precalculations for pairs ij
c
      call precal2x_2(iabprim,lcij,iapb,i1apb,isab)
c----------
      lpar_exe=lpar_exe+1
      lpar_num=lpar_num+nij_uniqe
      call precal2a_2(dbl_mb(ibas),dbl_mb(inuc),inx,ibl,nbl2,
     *               bl(iabprim+lcij),bl(iabprim+2*lcij),
     *               lcij,bl(iaa),bl(ibb),bl(ieab),bl(icis),bl(icjs),
     *               bl(ixab),bl(ixp),bl(icij),bl(ifij), bl(itxab),
     *               bl(igcij),ngci1,ngcj1,'left ',
     *               nij_uniqe,bl(ij_uniqe_p) )
c----------
c for abnia
      if(mmax.gt.2) then
         call precal2b_2(mmax1,lcij, bl(i1apb),bl(iabnia))
      endif
c----------------------------------------
c precalculations for pairs kl
c
c10   if(kbl.eq.ibl) then
c10       kabprim=iabprim
c10       kldim  =ijdim
c10   endif
c10   if(npkl.ne.0) then
          call precal2x_2(kabprim,lckl,icpd,i1cpd,iscd)
c
          lpar_exe=lpar_exe+1
          lpar_num=lpar_num+nkl_uniqe
          call precal2a_2(dbl_mb(ibas),dbl_mb(inuc),inx,kbl,nbl2,
     *                  bl(kabprim+lckl),bl(kabprim+2*lckl),
c not an error, icc should be twice :
     *                  lckl,bl(icc),bl(icc),bl(iecd),bl(icks),bl(icls),
     *                  bl(ixcd),bl(ixq),bl(ickl),bl(ifkl), bl(itxcd),
     *                  bl(igckl),ngck1,ngcl1,'right',
     *                  nkl_uniqe,bl(kl_uniqe_p) )
c--
c for cdnia
         if(mmax.gt.2) then
            call precal2b_2(mmax1,lckl, bl(i1cpd),bl(icdnia))
         endif
c10   else
c10       call precdiag
c10   endif
c----------------------------------------
c for habcd
c
      if(mmax.gt.2) then
         call precal2c_2(lcij,lckl,bl(i1apb),bl(i1cpd),
     *                   bl(ihabcd),nfumax )
      endif
c---------------------------------------------------------------
c check for numerical stablility :
c
      call num_stability_2(dbl_mb(ibas),inx,
     *                     bl(ij_uniqe_p),bl(kl_uniqe_p))
c---------------------------------------------------------------
      end
c====================================================================
      subroutine precal2a_2(datbas,datnuc,inx,ibl,nbl2,
     *                      rapb,sab,lcij,
     *                      aaa,bbb,estab, cis,cjs,
     *                      xab,xparij,coefij,factij,txab,
     *                      gcij,ngci1,ngcj1,which ,
     *                      nij_uniqe,ij_uniqe )
c
      implicit real*8 (a-h,o-z)
      character*5 which
c
c for gradient derivatives only :
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
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
      dimension rapb(lcij),sab(nij_uniqe,lcij) 
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,*)
c
      dimension cis(*),cjs(*)
      dimension xab(nij_uniqe,3), xparij(nij_uniqe,3,lcij,3)
      dimension estab (nij_uniqe,lcij)
      dimension coefij(nij_uniqe,lcij), factij(lcij)
      dimension txab(nij_uniqe,3,*)
      dimension gcij(ngcj1,ngci1,lcij)
      dimension aexp(100),bexp(100)
c only for where='forc' or 'hess'
      dimension aaa(nij_uniqe,*),bbb(nij_uniqe,*)
c
      dimension ij_uniqe(nij_uniqe)
c
      dimension xa(3),xb(3)
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
c for gen.contr
c
      ngctot=ngci+ngcj+ngck+ngcl
c-------------------------------------------------------------
c precalculations for the pairs IJ :
c
      do 100 ijpar=1,nij_uniqe
        ijcs=ij_uniqe(ijpar)
        call get_ij_half(ijcs,ics,jcs)
c
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
            if(where.eq.'forc') aaa(ijpar,is1)=aa
            if(where.eq.'hess') aaa(ijpar,is1)=aa
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
            if(where.eq.'forc' .and. which.eq.'left ') bbb(ijpar,js1)=bb
            if(where.eq.'hess' .and. which.eq.'left ') bbb(ijpar,js1)=bb
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
                       gcij(jg+1,ig+1,ji)=gcj
  210                continue
               endif
c--------------------------------
c
            rapb1=rapb(ji)
            sab1 =sab(ijpar,ji)
            aa1=aa*rapb1
            bb1=bb*rapb1
c
c-overlap   coefij(ijpar,ji)=coefij(ijpar,ji)*rapb1*sab1
            coefij(ijpar,ji)=coefij(ijpar,ji)      *sab1
c--------------------------------
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
         do 151 ijpar=1,nij_uniqe
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
         do 152 ijpar=1,nij_uniqe
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
           do 153 ijpar=1,nij_uniqe
              txab(ijpar,1,ijs1)=txab(ijpar,1,ijs1)*rapb1
              txab(ijpar,2,ijs1)=txab(ijpar,2,ijs1)*rapb1
              txab(ijpar,3,ijs1)=txab(ijpar,3,ijs1)*rapb1
  153      continue
         endif
      endif
c------------------------------------------------------------------
      end
c====================================================================
      subroutine precal2b_2(mmax1,lcij, rapb, abnia)
      implicit real*8 (a-h,o-z)
      dimension rapb(*), abnia(mmax1,*)
      data half /0.5d0/
c
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
      common /logicd/ hnia(3,1)
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension rapb(*),rcpd(*), habcd(3,nfumax,*)
c
      data half /0.5d0/
c
c-----------------------
c Do it only for one pair (the first one)
c
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
c
      end
c====================================================================
      subroutine prec4neg_2(isbl_size,isbl_point,isbl_q,ndiag,
     *                      ij,kl, lc12, lc34,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq, 
     *                      map_ij,map_kl,
     *                      nij_uniqe,nkl_uniqe,
coutput
     *                      rppq,rhoapb,rhocpd,rys,const,nbls1,index)
c--------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /neglect/ eps,eps1,epsr
c-------
      dimension map_ij(isbl_size),map_kl(isbl_size)
      dimension isbl_q(*)
      dimension index(*)
c-------
      dimension apb(lc12),cpd(lc34)
      dimension coefij(nij_uniqe,lc12),coefkl(nkl_uniqe,lc34)
      dimension  estab(nij_uniqe,lc12), estcd(nkl_uniqe,lc34)
      dimension xp(nij_uniqe,3,*),xq(nkl_uniqe,3,*)
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
c-----------old--------------
c        ijcsq=ijpres2(iq)
c        klcsq=klpres2(iq)
c        ijpar=map_ij(ijcsq)
c        klpar=map_kl(klcsq)
c-----------old--------------
cnew : map to uniqe pairs 
c
         ijpar=map_ij(iqp)
         klpar=map_kl(iqp)
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
c10         if(ndiag.eq.0) coef12=coef12*par268
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
      subroutine precspec_2(isbl_size,isbl_point, isbl_q,ndiag,
     *                      ij,kl, lc12, lc34, indxij,indxkl,
     *                      estab,estcd, apb,cpd,coefij,coefkl,xp,xq,
     *                      rapb,rcpd,txab,txcd,
     *                      map_ij,map_kl,
     *                      nij_uniqe,nkl_uniqe,
c     output
     *                      rys,const,xpqr,txxr,nbls1,index)
c-------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /neglect/ eps,eps1,epsr
c---------------
      dimension map_ij(isbl_size),map_kl(isbl_size)
      dimension isbl_q(*)
      dimension indxij(*),indxkl(*),index(*)
c---------------
      dimension apb(lc12),cpd(lc34)
      dimension rapb(lc12),rcpd(lc34)
      dimension coefij(nij_uniqe,lc12),coefkl(nkl_uniqe,lc34)
      dimension  estab(nij_uniqe,lc12), estcd(nkl_uniqe,lc34)
      dimension txab(nij_uniqe,3,*),txcd(nkl_uniqe,3,*)
      dimension xp(nij_uniqe,3,*),xq(nkl_uniqe,3,*)
c
      dimension rys(*),const(*)
      dimension xpqr(3,*),txxr(3,*)
c     
      data zero,one /0.d0 , 1.d0/
c---------------------------------------------------------------
      par268=     8.d0
c---------------------------------------------------------------
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
c----------------old-------------
c        ijcsq=ijpres2(iq)
c        klcsq=klpres2(iq)
c        ijpar=map_ij(ijcsq)
c        klpar=map_kl(klcsq)
c----------------old-------------
cnew : map to uniqe pairs
c
         ijpar=map_ij(iqp)
         klpar=map_kl(iqp)
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
c10         if(ndiag.eq.0) coef12=coef12*par268
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
c     pnl        if(ndiag.eq.0.and.ijpar.eq.klpar) then
c     pnl           const(ijkl1)=const(ijkl1)*0.5d0
c     pnl        endif
c     
         endif
 100  continue
c     
      nbls1=ijkl1
      if(nbls1.eq.0) return
c     
c     new
      if (ityp.gt.1 .or. jtyp.gt.1 ) then
         rapbxxx=rapb(ij)
         do 210 i=1,nbls1
            ijkl=index(i)
            ijpar=indxij(ijkl)                  ! uniqe pair
            xpqr(1,i)=xpqr(1,i)*rapbxxx 
            xpqr(2,i)=xpqr(2,i)*rapbxxx 
            xpqr(3,i)=xpqr(3,i)*rapbxxx
c     
            txxr(1,i)=txab(ijpar,1,ij)
            txxr(2,i)=txab(ijpar,2,ij)
            txxr(3,i)=txab(ijpar,3,ij)
 210     continue
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
 220     continue
      endif
c     
      end
c====================================================================
      subroutine xwpq_2(nbls1,xwp,xwq,p1234,lc12, lc34,lcij,lckl,
     *                  indxij,indxkl,index, nij_uniqe,nkl_uniqe,
     *                  rppq,xp,xq,xpp,xqq,
c98  *                  txab,txcd,abcd,apb,rcpd,cpd,rapb)
     *                  txab,txcd,abcd,apb,rcpd,cpd,rapb,
     *                  const_max)
c---------------------------------------------------------------
c const_max is maximum of abs values of const(nbls1) 
c OUTPUT :  xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3), abcd
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /flops/ iflop(20)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      dimension xwp(nbls1,3),xwq(nbls1,3),p1234(nbls1,3)
      dimension indxij(*),indxkl(*),index(*)
c
      dimension xp(nij_uniqe,3,*),xq(nkl_uniqe,3,*)
      dimension xpp(nij_uniqe,3,*),xqq(nkl_uniqe,3,*)
      dimension txab(nij_uniqe,3,*),txcd(nkl_uniqe,3,*)
c
      dimension apb(lc12),rapb(lc12)
      dimension cpd(lc34),rcpd(lc34)
c---------------------------------------------------------------
c accuracy :
      logical stable_2
      common /stability_2/ stable_2      
      common /stability_3/ b_over_a(1000) 
c
c used only here and in trobsa.f :
c
      logical stable 
      common /tracy_stability/ stable(1)
c
c set up in set_lost_allowed() in texas_hf 
c
      common /allow_lost/ lost_allow
c---------------------------------------------------------------
c For test only : counter of all & unstable blocks of primitive quartets
c used here and in raktest.F
c
c     common /count_primitive/ x_prim_all,
c    *                         x_prim_wr1, x_prim_wr2, x_prim_wr3
c     common /count_efficient/ xobsa2,xobsa2_ineff,xobsa2_more,
c    *                         xobsa4,xobsa4_ineff,xobsa4_more
c---------------------------------------------------------------
      rppq1=rppq
      rcpd1=rcpd(lckl)
c
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
c for Tarcy's recursive :
c
         p1234(i,1)=(txab(ijpar,1,lcij)+txcd(klpar,1,lckl))*rcpd1
         p1234(i,2)=(txab(ijpar,2,lcij)+txcd(klpar,2,lckl))*rcpd1
         p1234(i,3)=(txab(ijpar,3,lcij)+txcd(klpar,3,lckl))*rcpd1
c
  100 continue
c
c---------------------------------------------------------------
c for accuracy in Tracy's recursive : 
c
      stable(1)=.true.
      abcd=apb(lcij)*rcpd1
c
      x_prim_all=x_prim_all + 1.d0
c
        lost_digit_1=0
        lost_digit=0
        if( .not.stable_2 ) then
ctest      x_prim_wr1=x_prim_wr1 + 1.d0
           boa=b_over_a(lcij)
           if( boa .gt. 0.d0 ) then
ctest         x_prim_wr2=x_prim_wr2 + 1.d0
c
c............ boa=0.d0 means that it was ok ......
c
              lost_const=int( log10(const_max) )
c 
              xost_digit_1= log10(boa*rcpd(lckl))
              lost_digit_1= int( xost_digit_1 )
              lost_digit  = (nskl-1)*lost_digit_1 + lost_const
c
              if(lost_digit.gt.lost_allow) then
                 stable(1)=.false.
                 abcd=rapb(lcij)*cpd(lckl)
ctest
ctest            x_prim_wr3 =x_prim_wr3 + 1.d0
ctest
c                write(6,*)'const()=',(const(ii),ii=1,nbls1)
c                write(6,*)'a+b=',apb(lcij),' c+d=',cpd(lckl)
c                write(6,*)'Lost : 1 step=',lost_digit_1,
c    *              ' const=',lost_const,
c    *              ' total=',lost_digit,'(allowed=',lost_allow,')',
c    *              '=> stable=',stable(1),
c    *              ' an.mom. (',nqi-1,nqj-1,'|',nqk-1,nql-1,')'
              endif
           endif
        endif
ctest
c
c---------------------------------------------------------------
c calculate how many times the program goes according to :
c (1) lobsa=2 (when it is stable)
c
c     if(stable(1)) then
c        xobsa2=xobsa2+1.0
c        if(nsij.LT.nskl) then
c           xobsa2_ineff=xobsa2_ineff + 1.0
c           xobsa2_more =xobsa2_more + dble(nskl-nsij)
c        endif
c     else
c        xobsa4=xobsa4+1.0
c        if(nsij.GT.nskl) then
c           xobsa4_ineff=xobsa4_ineff + 1.0
c           xobsa4_more =xobsa4_more + dble(nsij-nskl)
c        endif
c     endif
c---------------------------------------------------------------
c old fashion way :
c
c     if(nsij.ge.nskl) then
c        stable(1)=.true.
c        abcd=apb(lcij)*rcpd1
c     else
c        stable(1)=.false.
c        abcd=rapb(lcij)*cpd(lckl)
c     endif
c---------------------------------------------------------------
c     write(6,*)' nsij,nskl=',nsij,nskl,
c    *              ' an.mom. (',nqi,nqj,'|',nqk,nql,')',
c    *                                  ' stable=',stable(1)
c
c for tests :
c
c     write(6,*) ' forced to go by lobsa=2'
c        stable(1)=.true.
c        abcd=apb(lcij)*rcpd1
c
c     write(6,*) ' forced to go by lobsa=4'
c        stable(1)=.false.
c        abcd=rapb(lcij)*cpd(lckl)
c
c---------------------------------------------------------------
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
c  this subroutine constitues the special code for
c  two types of integrals over nbls quartets of primitive shells
c  1. (ss|ss)
c  2. (xs|ss),(sx|ss),(ss|xs),(ss|sx) where x=p or l
c  these integrals are also contracted here.
c
c  this routine is called from the twoe subroutine.
c
c
c  input
c  ------
c  all precalculated values for whole block :
c
c  const(nbls) - contains consts=pi3sabcd/(pqsqrt(ppq)) for all int.
c  rysx(nbls) - contains  xrys i.e. arguments for fm,ft0 routines
c  xp,xp      - geometry for p,q
c
c  output
c  ------
c  buf(nbls,) -contains final contracted integrals
c---------------------------------------------------------------
c
c memory for f00,f11 :
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
c special code for (ss ss) integrals
c
      if(mmax.eq.1) then
          do 2031 i=1,nbls1
          buf1(i,1)=const(i)*bl(if00+i)
 2031     continue
      go to 203
      endif
c
c  special code for (ps ss), (sp ss) (ss ps) and (ss sp)
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
      subroutine get_min_exp(datbas,inx,nij_uniqe,ij_uniqe,
     *                       a_exp_min, b_exp_min )
c---------------------------------------------------------------
c Because of numerical instability in Tracy's recursive :
c select the minimum a_ & b_exponent
c this routine is called for right side pairs ONLY
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension datbas(13,*)
      dimension inx(12,*)
      dimension ij_uniqe(nij_uniqe)
c---------------------------------------------------------------
c
              ijcs=ij_uniqe(1)    
              call get_ij_half(ijcs,ics,jcs)
              ie=inx(5,ics)
              je=inx(5,jcs)
              a_exp_min=datbas(1,ie)
              b_exp_min=datbas(1,je)
c---------------------------------------------------------------
      end
c====================================================================
      subroutine num_stability_2(datbas,inx,ij_uniqe,kl_uniqe)
      implicit real*8 (a-h,o-z)
      logical stable_2
      common /stability_2/ stable_2       ! use here and in xwpq_2
      common /stability_3/ b_over_a(1000) ! use here and in xwpq_2
c 1000 above is enough for contractions as long as 31 (31*31=961)
c
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c number of digit allowed to loose (calculated in prepint2)
      common /allow_lost/ lost_allow
c
      dimension datbas(13,*)
      dimension inx(12,*)
      dimension ij_uniqe(*)  ! ij_uniqe(nij_uniqe)
      dimension kl_uniqe(*)  ! kl_uniqe(nkl_uniqe)
c-------------------------------------------------------------
c June 1998
c Numerical stability of Tracy's recursive depends on the relative
c value of exponents involved. For (ab/cd) instability occurs when
c the second exponent (b) is much bigger than others.
c In such a case one has to proceed as follows:
c construct (s,s|i+j+k+l,s) integrals and then shift ang.mom.
c (s,s|i+j+k+l,s) ---> (i+j,s|k+l,s)
c from center 3 to center 1 (tracy's recursive) 
c 
c Before numerical instability in Tracy's recursive has been
c realized the program ALWAYS calculated (i+j+k+l,s|s,s) and
c then shifted ang. mom. from 1 to 3.
c-------------------------------------------------------------
c calculate maximum of 2*b/max(1,a) and maximum of 1/(c+d)
c these two above will determine numerical stability in
c Tracy's recursive 
c-------------------------------------------------------------
c Lost digit :
c 
c in one recursive step :
c
c
c                    2*b
c          log[-------------------]
c                max(1,a) * (c+d)
c
c in all (nskl-1) steps :
c
c                    2*b
c          log[-------------------] * (nskl-1) + log(const_max)
c                max(1,a) * (c+d)
c
c
c where const_max is max of const(ijkl) constracted in prec4neg_ .
c
c Number of lost digit given above gives the final accuracy for 
c calculated integrals :
c
c          acc_obatained = 10**-(14-lost_digit)
c
c The requested accuracy is given by the user's threshold
c
c          acc_requested = 10**-10 =10**-(14-4)  for example.
c
c (14 comes from double precission)
c-------------------------------------------------------------
c
           stable_2=.true.
c
           if(mmax.gt.2 .and. nskl.gt.1) then
c
c left side pairs:
c
              ijcs=ij_uniqe(1)    
              call get_ij_half(ijcs,ics,jcs)
              ie=inx(5,ics)
              ia=inx(1,ics)+1
              je=inx(5,jcs)
              ja=inx(1,jcs)+1
                 b_exp_max=datbas(1,ja)
                 a_exp_min=datbas(1,ie)
                 a_1_max=max(1.d0,a_exp_min)
c
c two b over max(1,a)  :
c
                 boa_max=2.d0*b_exp_max/a_1_max
c
              if(boa_max .lt. 10.d0) return  ! stable_2=.true. 
c
c right side pairs:
c
              klcs=kl_uniqe(1)    
              call get_ij_half(klcs,kcs,lcs)
              ke=inx(5,kcs)
              le=inx(5,lcs)
                 c_exp_min=datbas(1,ke)
                 d_exp_min=datbas(1,le)
                 cpd_min=c_exp_min+d_exp_min
                 cpd_max=1.d0/cpd_min
c
c Lost allow is -log(10**-14)=14 (double precision)
c                    minus (for example)
c               -log(10**-10)=10 (given integral's threshold)
c
c Lost commited :
c
              xost_digit_1= log10( boa_max*cpd_max )
              lost_digit_1= int( xost_digit_1 )
c
              lost_digit  = (nskl-1)*lost_digit_1
c
c this is very conservative because the contraction coeff. are
c not taken into account. It is done later in XWPQ_2 .
c
c     write(6,*)' lost_digits: in 1 step=',lost_digit_1,
c    *          ' ; in all steps=',lost_digit
c
              if(lost_digit.gt.lost_allow) stable_2=.false.
c
c     write(6,*)'num_stability_2 : stable_2=',stable_2
c-------------------------------------------------------------
              if(stable_2) return
c-------------------------------------------------------------
c calculate and save 2*b/max(1,a) for each left side contraction 
c
              ij=0
              do is=ia,ie
                 aa=datbas(1,is)
                 a_1_max=max(1.d0,aa)
                 a2=2.d0/a_1_max
                 do js=ja,je
                    bb=datbas(1,js)
                    ij=ij+1
c
                    b_over_a(ij)=0.d0
                    boa=bb*a2
                    boa_max_cpd=boa*cpd_max
                    if(boa_max_cpd .ge. 10.d0) then
                       xost_digit_1= log10( boa_max_cpd )
                       lost_digit_1= int( xost_digit_1 ) 
                       lost_digit  = (nskl-1)*lost_digit_1
                       if(lost_digit.gt.lost_allow) then
                          b_over_a(ij)=boa
                       endif
                    endif
                 enddo ! js=ja,je
              enddo    ! is=ia,ie
c-------------------------------------------------------------
           endif   !  if(mmax.gt.2 .and. nskl.gt.1) then
c
      end
c===================================================================
