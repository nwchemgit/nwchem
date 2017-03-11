c $Id$
cccc  subroutine prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,
      subroutine prepint2(bl,    inuc,ibas,na,nbf,nsh,ncf,ncs,inx,
     *                    lcore,maxprice,scftype)
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
cold      parameter (lpar1=34,lpar2=6545,lpar3=5984,lpar4=10,lpar5=33)
c up to ll,ll :
c     parameter (lpar1=38,lpar2=9139,lpar3=8436,lpar4=11,lpar5=37)
c up to mm,mm :
c     parameter (lpar1=42,lpar2=12341,lpar3=11480,lpar4=12,lpar5=41)
c---------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical firstd
      character*11 scftype
c
#include "texas_lpar.fh"
cxx
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inucx,ibasx,nax,nbfx,nshx,ncfx,ncsx
c
      common /memor1b/ nbl2
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
c
      common /memors/ nsymx,ijshp,isymm
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
c setup the infob common :
c
      inucx=inuc
      ibasx=ibas 
      nax  =na
      nbfx =nbf
      nshx =nsh
      ncfx =ncf
      ncsx =ncs
c
cnew  nsymx=nsym
      nsymx=0
c----------------------------------------------------------
       call fprep
c
c----------------------------------------------------------
c setup logical matrices - old logobsa; common  /logicd/, logic1-11/
c
       call datlog(inx,ncs,lpar1,lpar2,lpar3,lpar4,lpar5)
c
c----------------------------------------------------------
c blocking procedure for pairs of contracted shells
c
       call txs_second(tim3)
       call blockin2(bl,lcore,inx,nbl2)
c
       call txs_second(tim4)
       timblok2=tim4-tim3
c
c----------------------------------------------------------
c
      call txs_second(tim2)
      timprep2=tim2-tim1
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
      common /time5/ tcheck,tuniq2,tmap4u
c
      if(tblock.eq.0.0d0) then
ctime
        tprec2=0.0d0
        tconv1=0.0d0
        tconv2=0.0d0
        ttrobs=0.0d0
        tpre4n=0.0d0
        txwpq =0.0d0
        tassem=0.0d0
        tamshf=0.0d0
        tdesti=0.0d0
        terint=0.0d0
        tzeroi=0.0d0
        tspeci=0.0d0
        ttrans=0.0d0
        tderiv=0.0d0
c
        tcheck=0.0d0
        tuniq2=0.0d0
        tmap4u=0.0d0
ctime
      else
c
        write(8,120) tblock,tprec2,
     *                             tpre4n,tzeroi,txwpq,tconv1,tconv2,
     *                             ttrobs,tassem,tderiv,tamshf,ttrans,
     *               terint,tspeci,
     *               tblock+tprec2+terint+tspeci,      ! calculation
     *                             tuniq2,tmap4u,tcheck,tdesti,
     *               tuniq2+tmap4u+tcheck+tdesti,      ! interface
     *       tblock+tprec2+terint+tspeci + tuniq2+tmap4u+tcheck+tdesti 
c     call util_flush(8)
  120 format(
     *       '========================================================'/
     *       '                                                        '/
     *       'T I M I N G  O F  T E X A S  S U B R O U T I N E S (sec)'/
     *       '                                                        '/
     *       '========================================================'/
     *       'time for BLOCKIN=',f10.1/
     *       '--------------------------------------------------------'/
     *       'time for PRECALC=',f10.1/
     *       '--------------------------------------------------------'/
     *       '                             time for prec4n =',f10.1/   
     *       '                             time for zeroin =',f10.1/   
     *       '                             time for xwpq   =',f10.1/   
     *       '                             time for conv1x =',f10.1/   
     *       '                             time for conv2x =',f10.1/   
     *       '                             time for trobsa =',f10.1/   
     *       '                             time for assemb =',f10.1/   
     *       '                             time for derivx =',f10.1/   
     *       '                             time for amshif =',f10.1/   
     *       '                             time for transf =',f10.1/   
     *       '                             ---------------------------'/
     *       'time for ERINTEG=',f10.1/
     *       '--------------------------------------------------------'/
     *       'time for ERINTSP=',f10.1/
     *       '========================================================'/
     *       'CALCULATION TIME=',f10.1,'   (pure texas calculations)  '/
     *       '========================================================'/
     *       '                             time for uniq_pa=',f10.1/   
     *       '                             time for map2uni=',f10.1/   
     *       '                             time for checkio=',f10.1/   
     *       '                             time for destiny=',f10.1/   
     *       '                             ---------------------------'/
     *       'INTERFACE   TIME=',f10.1,'   (internal interface)       '/
     *       '========================================================'/
     *       'TOTAL TEXAS TIME=',f10.1/
     *       '========================================================')
c
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
      data pi256 / 0.012271846d0 /
c
      eps =thres
      eps1=1.0d0/eps
      epsr=pi256*eps*eps
c
ccc   write(6,*)' integ thresh from setup_thresh()=',eps
c----------------------------------------------------------
      end
