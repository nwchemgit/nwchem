c $Id: prepint.f,v 1.7 1996-07-12 18:25:46 pg481 Exp $
      subroutine prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,
     *                    lcore,nsym,maxprice,scftype,iret_limit)
c---------------------------------------------------------------
c EPS is not used here, threshold is now set by calling 
c texas_set_accy before texas_init is called.
c---------------------------------------------------------------
c dimensions for logical matrices in commons logicd,logic1-11 :
c
c up to ff,ff :
c     parameter (lpar1=14,lpar2= 455,lpar3= 364,lpar4=5,lpar5=13)
c up to gg,gg :
c     parameter (lpar1=18,lpar2= 969,lpar3= 816,lpar4=6,lpar5=17)
c up to hh,hh :
c     parameter (lpar1=22,lpar2=1771,lpar3=1540,lpar4=7,lpar5=21)
c up to ii,ii :
c     parameter (lpar1=26,lpar2=2925,lpar3=2600,lpar4=8,lpar5=25)
c up to jj,jj :
c     parameter (lpar1=30,lpar2=4495,lpar3=4060,lpar4=9,lpar5=29)
c up to kk,kk (which is ii,ii + 2 )
      parameter (lpar1=34,lpar2=6545,lpar3=5984,lpar4=10,lpar5=33)
c up to ll,ll :
c     parameter (lpar1=38,lpar2=9139,lpar3=8436,lpar4=11,lpar5=37)
c up to mm,mm :
c     parameter (lpar1=42,lpar2=12341,lpar3=11480,lpar4=12,lpar5=41)
c---------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical first,firstd
      character*11 scftype
c
      common /timing2/ timprep2,timstor2,timreca2,timblok2,timblok4
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common /logicd/ hnia(3,lpar2)
      common /logic1/ ndege(lpar4)
      common /logic2/ len(lpar4)
      common /logic3/ lensm(lpar5)
c
      common /logic4/ nfu(lpar1)
      common /logic5/ icoor(lpar2)
      common /logic6/ icool(lpar2)
      common /logic7/ ifrst(lpar2)
      common /logic8/ ilast(lpar2)
      common /logic9/ nia(3,lpar2)
      common /logic10/ nmxyz(3,lpar2)
      common /logic11/ npxyz(3,lpar3)
cxx
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inucx,ibasx,nax,nbfx,nshx,ncfx,ncsx
c
      common /memor1/ iisd,jjsd,ijbld
cpnl: common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1a/ npard,ncost,nsupb, mxsize,nasize
      common /memor1b/ nbl2,nbloks
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
c
      common /memors/ nsymx,ijshp,isymm
      common /symmet/ rnsym
c----------------------------------------------------------------
c
      dimension bl(*)
      dimension inx(12,*)
c
c----------------------------------------------------------
c
      call txs_second(tim1)
c
c----------------------------------------------------------
c set up integral threshold :
c
c      call setup_thres(eps) 
c now done in texas_set_accy which calls setup_thresh() routine 
c----------------------------------------------------------
c zero out timings
c
      timprep2=zero
      timstor2=zero
      timreca2=zero
      timblok4=zero
c
c***************************
c* setup the infob common :
c
      inucx=inuc
      ibasx=ibas 
      nax  =na
      nbfx =nbf
      nshx =nsh
      ncfx =ncf
      ncsx =ncs
c
      nsymx=nsym
      rnsym=one/(one+float(nsym))
c***************************
c
       call fprep
c
c***************************
c*  setup logical matrices - old logobsa; common  /logicd/, logic1-11/
c
       call datlog(inx,ncs,lpar1,lpar2,lpar3,lpar4,lpar5)
c
c***************************
c* blocking procedure for pairs of contracted shells
c
       call txs_second(tim3)
       call blockin2(bl,lcore,inx,nbl2,nbloks,iret_limit)
c
       if(iret_limit.gt.0) return
c
       call txs_second(tim4)
       timblok2=tim4-tim3
c
c***************************
C DO NOT DO IT for PNL
c
c decide how to run scf :
c nbloks=nbl2*(nbl2+1)/2
c
CPNL   nbloks=nbl2*(nbl2+1)/2
CPNL   call howtorun(nbloks,maxprice,bl(ncost),scftype)
c
c  For each super-block :
c   ncost(ikbl) is now 0 - re-calculate int. from the ikbl block
c                      1 - store int. from the ikbl block
c***************************
c
      call txs_second(tim2)
      timprep2=tim2-tim1
c
      return
      end
c----------------------------------
      subroutine howtorun(nbloks,maxprice,ncost,scftype)
      character*11 scftype
      dimension ncost(*)
c---------------------------------
c This is loop over super-blocks :
c
      i1=0
      do 100 ikbl=1,nbloks
      iprice=ncost(ikbl)
         if(iprice.le.maxprice) then
           ncost(ikbl)=0
         else
           ncost(ikbl)=1
           i1=i1+1
         endif
  100 continue
c
      scftype='semi-direct'
      if(i1.eq.nbloks) scftype='non -direct'
      if(i1.eq.0     ) scftype='full-direct'
c
c-------------------------------------------------
c
      if(scftype.eq.'semi-direct') then
c--->    write(6,210) nbloks,i1,nbloks-i1,maxprice
         write(8,210) nbloks,i1,nbloks-i1,maxprice
  210 format(
     *' Total number of Super-blocks of contracted quartet =',i15/
     *' integrals from ',i12,' Super-blocks will be stored on a disk'/
     *' integrals from ',i12,' Super-blocks will be re-calculated'/
     *' with the maximum price of calculations accepted to pay ',i8/)
      endif
c
      end
c================================================================
c
      subroutine timepr(tblock)
      implicit real*8 (a-h,o-z)
c
      common /timex/ tconv1,tconv2,ttrobs
      common /time0/ tprec2
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time2/ terint,ttrans
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
c
      data zero /0.d0/
c
      if(tblock.eq.zero) then
ctime
        tprec2=zero
        tconv1=zero
        tconv2=zero
        ttrobs=zero
        tpre4n=zero
        txwpq =zero
        tassem=zero
        tamshf=zero
        tdesti=zero
        terint=zero
        tzeroi=zero
        tspeci=zero
        ttrans=zero
        tderiv=zero
ctime
      else
c
        write(8,120) tblock,tprec2,tpre4n,tzeroi,txwpq,tconv1,tconv2,
     *          ttrobs,tassem,tderiv,tamshf,terint,tspeci,ttrans,tdesti
c
  120   format(
     *  /'***************************',2x,'***************************'/
     *        '   S U B R O U T I N E       T I M I N G  (in sec)     '/
     *   '---------------------------',2x,'---------------------------'/
     *        'time for BLOCKIN=',f10.2,'  time for PRECL2 =',f10.2/
     *   '---------------------------',2x,'---------------------------'/
     *   'time for prec4n =',f10.2/
     *   'time for zeroin =',f10.2/
     *   'time for xwpq   =',f10.2/
     *   'time for conv1x =',f10.2/
     *   'time for conv2x =',f10.2/
     *   'time for trobsa =',f10.2/
     *   'time for assemb =',f10.2/
     *   'time for derivx =',f10.2/
     *   'time for amshif =',f10.2/
     *   '---------------------------',2x,'---------------------------'/
     *   'time for ERINTEG=',f10.2,'  time for ERINTSP=',f10.2/
     *   '---------------------------',2x,'---------------------------'/
     *    'time for TRANSFO=',f10.2,'  time for DESTINY=',f10.2/
     *   '***************************',2x,'***************************')
c
        write(8,125) tblock+tprec2+terint+tspeci+ttrans+tdesti
  125   format('BLOCKIN+PRECL2+ERINTEG+ERINTSP+TRANSF+DESTINY=',f10.2/)
      endif
c
      end
c================================================================
      subroutine setup_thresh(thres)
      implicit real*8 (a-h,o-z)
c----------------------------------------------------------------
c Set up the integral threshold :
c----------------------------------------------------------------
c Estimator1=16/(pi)**1/2 * Sab*Scd * 1/(a+b+c+d)**1/2
c
c  where Sab=(a+b)**-1 * (a*b)**+3/4 *exp(-ab/(a+b) * Rab**2 )
c
c
c    Estimator1 > eps    do integrals .
c
c It is done in a quadratic form by  Es1*Es1> eps**2
c
c    Sab**2  *  Scd**2 > pi/256 * eps**2  * (a+b+c+d)
c
c  epsr= pi/256 * eps**2
c
c  eps =eps -original thresh. for integrals,
c  eps1=1/eps  used to rescale integr. in precalc2a subroutine
c  epsr=pi/256 *eps**2-used in prec4neg and precspec subroutines.
c----------------------------------------------------------------
c pi256=pi/256
c
      common /neglect/ eps,eps1,epsr
      data one,pi, pi256 /1.d0 , 3.1415926535d0 , 0.012271846d0 /
c
      eps =thres
      eps1=one/eps
      epsr=pi256*eps*eps
c
c----------------------------------------------------------
      end
