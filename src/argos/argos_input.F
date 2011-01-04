      subroutine argos_input(irtdb)
c
c $Id$
c
      implicit none
c
      real*8 one
      parameter(one=1.0d0)
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item
      character*80 project,server
c
      integer mdalgo,iset,iforw,msdit,mcgit,ncgcy,mrun,ibinar,ipme
      integer maxlam,npgdec,issscl,mequi,mdacq,ldacq,ipscal,itscal
      integer ipolt,mshitw,mshits,npstep,npener,npdist,nptmng,nfpair
      integer nfrdf,nfcntr,nfoutp,nfstat,nfrest,keepr,nfcoor,nfscoo
      integer nfvelo,nfsvel,nfprop,nffree,nfsync,iformt,madbox,loadb
      integer nfforc,nfwfor,nfsfor,npstat,nfesp,impfr,impto,nftri
      integer nflong,nfgaus,mropt,nfauto,nfconv,nffet,idebug,mpolit
      integer npatom,nptopw,nptops,npx,npy,npz,npg,nbx,nby,nbz,itrack
      integer ngx,ngy,ngz,numfix,ifix,noshak,nfslow,nfefld,nfsfld
      integer ngl,isetp1,isetp2,ifield,iapprx,nacfl,npxpct,imfft
      integer ireact,itemp,nscb,idscb(10),memlim,morder,isolvo,lfout6
      integer mwmreq,msareq,itest,nodpme,lbpair,ireset,npmf,nopack
      integer ndistr,ndaver,idevel,itime(500),nftime,idtime,nfdrss
      integer itload,ihess,latom,icomm,nfanal,iport,iguide,icentr
      integer mbbreq,nfwcoo,nfwvel,ioload,iprof,isload,icntrl,istart
      integer ifidi,ncoll,ipopt,ivopt,ilambd,mlambd,includ,npforc
      integer iprop,iprpmf,nfhop,ipbtyp,iradgy,idifco,nfnewf,ncolt
      integer icmopt,imembr,icrit,nbget,nprec,nfload,ncp
      real*8 dx0sd,dxsdmx,dxmsd,dx0cg,dxcgmx,dxmcg,edacq,fdacq,ddacq
      real*8 delta,stime,tstep,prsext,prsrlx,compr,tmpext,tmprlx
      real*8 tmsrlx,tgauss,rshort,rlong,ptol,tlwsha,tlssha,factld
      real*8 rrdf,xfield,xfvect(3),xffreq,weight,dielec,ealpha
      real*8 frgaus,rbox,disrlx,drsscl,eqm,tguide,fguide,rfix,fcoll
      real*8 tmpext2,tann1,tann2,scaleq,rsgm,facpmf,rhop,thop
      character*12 afix
c
      integer i,mtimes,len
      character*255 filcmd
      integer lfncmd,numcmd,ifirst,ilast
      logical expert
c
#include "bitops.fh"
c
      lfncmd=69
      numcmd=0
c
      mtimes=250
c
c     set the defaults
c
      ncp=1
      nbx=0
      nby=0
      nbz=0
      npx=0
      npy=0
      npz=0
      npg=1
      includ=0
      itrack=0
      mdalgo=1
      iset=1
      isetp1=0
      isetp2=0
      iforw=1
      ilambd=0
      mlambd=0
      msdit=100
      mcgit=0
      ncgcy=10
      mrun=0
      maxlam=21
      npgdec=0
      issscl=0
      mequi=0
      mdacq=1000
      ldacq=500
      ipscal=0
      ipopt=0
      ivopt=0
      itscal=0
      nfgaus=0
      ipolt=0
      ifidi=0
      mpolit=10
      mshitw=100
      mshits=100
      npstep=0
      npstat=0
      npatom=0
      nptopw=0
      nptops=0
      npener=0
      npforc=0
      npdist=0
      nptmng=0
      nfpair=1
      nfesp=0
      nfrdf=0
      impfr=0
      impto=0
      nftri=0
      ngl=1000
      nflong=1
      nfcntr=0
      nfanal=0
      nfoutp=10
      nfstat=100
      nfrest=0
      keepr=0
      nfslow=0
      nfcoor=-1
      nfwcoo=-1
      nfscoo=-1
      nfvelo=-1
      nfwvel=-1
      nfsvel=-1
      nfforc=-1
      nfwfor=-1
      nfsfor=-1
      nfprop=0
      nffree=1
      nfsync=0
      nfauto=0
      nfconv=0
      nffet=0
      nfefld=0
      nfsfld=0
      nfnewf=0
      ibinar=0
      iformt=1
      madbox=6
      loadb=3
      nfload=0
      ireset=0
      mropt=1
      idebug=0
      icntrl=0
      ifield=0
      iapprx=0
      nacfl=1000
      npxpct=0
      ipme=0
      morder=4
      ngx=nbx
      ngy=nby
      ngz=nbz
      len=0
      nbget=0
      nprec=0
      iprop=0
      nfhop=0
      rhop=0.35d0
      thop=0.020d0
c      ifix=0
c      iunfix=0
      numfix=0
      noshak=0
      nscb=1
      idscb(1)=1
      ireact=0
      memlim=0
      isolvo=0
      lfout6=0
      iprpmf=0
      imfft=1
      mwmreq=0
      msareq=0
      mbbreq=0
      itest=0
      nodpme=0
      lbpair=10
      ndistr=0
      npmf=0
      ndaver=0
      idevel=0
      nftime=0
      nfdrss=0
      itload=0
      ioload=0
      isload=0
      iprof=0
      ihess=0
      latom=0
      icomm=0
      iport=0
      iguide=0
      icmopt=0
      imembr=0
      nopack=0
      istart=0
      ncoll=0
      ncolt=0
      ipbtyp=1
      iradgy=0
      idifco=0
      icrit=0
      fcoll=1.0d0
      rsgm=-1.0d0
c
      server=' '
c
      dx0sd=0.01d0
      dxsdmx=0.00001d0
      dxmsd=0.05d0
      dx0cg=0.001d0
      dxcgmx=0.00001d0
      dxmcg=0.05d0
      edacq=5.0d0
      ddacq=5.0d0
      fdacq=0.75d0
      delta=0.075d0
      stime=0.0d0
      tstep=0.001d0
      prsext=1.025e5
      prsrlx=0.5d0
      compr=4.53e-10
      tmpext=298.15d0
      tmpext2=298.15d0
      tmprlx=0.1d0
      tmsrlx=0.1d0
      tann1=0.0d0
      tann2=0.0d0
      tgauss=298.15d0
      rlong=0.9d0
      rshort=0.9d0
      ptol=0.001d0
      tlwsha=0.001d0
      tlssha=0.001d0
      factld=0.75d0
      rrdf=0.0d0
      xfvect(1)=0.0d0
      xfvect(2)=0.0d0
      xfvect(3)=1.0d0
      xffreq=0.0d0
      weight=0.0d0
      ealpha=0.0d0
      dielec=80.0d0
      frgaus=1.0d0
      rbox=0.0d0
      disrlx=0.1d0
      drsscl=one
      eqm=0.0d0
      tguide=0.0d0
      fguide=0.0d0
      scaleq=-1.0d0
c      do 4 i=1,msf
c      lsffix(i)=1
c    4 continue
      do 7 i=1,mtimes
      itime(i)=0
    7 continue
      idtime=0
      project=' '
      expert=.false.
c
c     read the nwchem input file
c
    1 continue
c
      if(.not.inp_read()) call md_abort('argos: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
      if(inp_compare(.false.,'system',item)) then
      if(len.gt.0) call md_abort('argos: duplicate system',0)
      if(.not.inp_a(project)) call md_abort('argos: error system',0)
      len=index(project,' ')-1
      filcmd=project(1:len)//'.cmd '
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='unknown',err=1999)
      rewind(unit=lfncmd)
      numcmd=0
      goto 2
      endif
c
      if(inp_compare(.false.,'expert',item)) then
      expert=.true.
      goto 2
      endif
c
      if(inp_compare(.false.,'finish',item)) then
      istart=1
      goto 2
      endif
c
      if(inp_compare(.false.,'resume',item)) then
      istart=2
      goto 2
      endif
c
      if(inp_compare(.false.,'prefetch',item)) then
      if(.not.inp_i(nbget)) nbget=-1
      goto 2
      endif
c
      if(inp_compare(.false.,'server',item)) then
      if(.not.inp_a(server)) call md_abort('argos: error server',0)
      if(.not.inp_i(iport)) call md_abort('argos: error port',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'boxmin',item)) then
      if(.not.inp_f(rbox)) call md_abort('argos: error minbox',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'segmentsize',item)) then
      if(.not.inp_f(rsgm))
     + call md_abort('argos: error segmentsize',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'boxes',item)) then
      if(.not.inp_i(nbx)) call md_abort('argos: error nbx',0)
      if(.not.inp_i(nby)) call md_abort('argos: error nby',0)
      if(.not.inp_i(nbz)) call md_abort('argos: error nbz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'cells',item)) then
      if(.not.inp_i(nbx)) call md_abort('argos: error nbx',0)
      if(.not.inp_i(nby)) call md_abort('argos: error nby',0)
      if(.not.inp_i(nbz)) call md_abort('argos: error nbz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'nodes',item)) then
      if(.not.inp_i(npx)) call md_abort('argos: error npx',0)
      if(.not.inp_i(npy)) call md_abort('argos: error npy',0)
      if(.not.inp_i(npz)) call md_abort('argos: error npz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'procs',item)) then
      if(.not.inp_i(npx)) call md_abort('argos: error npx',0)
      if(.not.inp_i(npy)) call md_abort('argos: error npy',0)
      if(.not.inp_i(npz)) call md_abort('argos: error npz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'groups',item)) then
      if(.not.inp_i(npg)) call md_abort('argos: error npg',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'cores',item)) then
      if(.not.inp_i(ncp)) call md_abort('argos: error ncp',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(iset)) call md_abort('argos: error set',2)
      if(iset.le.0) call md_abort('argos: illegal iset',2)
      if(iset.gt.3) call md_abort('argos: illegal iset',2)
      goto 2
      endif
c
      if(inp_compare(.false.,'pset',item)) then
      if(.not.inp_i(isetp1)) call md_abort('argos: error pset',2)
      if(isetp1.le.1) call md_abort('argos: illegal ipset1',2)
      if(isetp1.gt.3) call md_abort('argos: illegal ipset1',2)
      if(.not.inp_i(isetp2)) goto 2
      if(isetp1.eq.3) call md_abort('argos: illegal ipset1',2)
      if(isetp2.ne.3) call md_abort('argos: illegal ipset2',2)
      goto 2
      endif
c
      if(inp_compare(.false.,'lambda',item)) then
      if(.not.inp_i(ilambd)) call md_abort('argos: error lambda',0)
      if(.not.inp_i(mlambd)) call md_abort('argos: error lambda',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'numerical',item)) then
      ifidi=1
      goto 2
      endif
c
      if(inp_compare(.false.,'pbc',item)) then
      if(.not.inp_a(item)) call md_abort('argos: error pbc',1)
      if(inp_compare(.false.,'atom',item)) then
      ipbtyp=1
      elseif(inp_compare(.false.,'residue',item)) then
      ipbtyp=2
      elseif(inp_compare(.false.,'molecule',item)) then
      ipbtyp=3
      else
      call md_abort('argos: error pbc',2)
      endif
      goto 2
      endif
c
c     energy minimization
c
c     em [[sd] [<i>] [init <r>] [max <i>] [min <i>]]
c     [cg [<i>] [cy <i>] [init <r>] [max <r>] [min <i>]]
c
      if(inp_compare(.false.,'sd',item)) then
      if(.not.inp_i(msdit)) call  md_abort('argos: error em sd',2)
  101 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'init',item)) then
      if(.not.inp_f(dx0sd)) call md_abort('argos: error sd init',2)
      goto 101
      endif
      if(inp_compare(.false.,'min',item)) then
      if(.not.inp_f(dxsdmx)) call md_abort('argos: error sd min',2)
      goto 101
      endif
      if(inp_compare(.false.,'max',item)) then
      if(.not.inp_f(dxmsd)) call md_abort('argos: error sd max',2)
      goto 101
      endif
      if(inp_compare(.false.,'criterion',item)) then
      if(.not.inp_i(icrit)) call md_abort('argos: error sd crit',2)
      goto 101
      endif
      endif
c
c     conjugate gradient
c
      if(inp_compare(.false.,'cg',item)) then
      if(.not.inp_i(mcgit)) call  md_abort('argos: error em cg',2)
  103 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'init',item)) then
      if(.not.inp_f(dx0cg)) call md_abort('argos: error cg init',2)
      goto 103
      endif
      if(inp_compare(.false.,'min',item)) then
      if(.not.inp_f(dxcgmx)) call md_abort('argos: error cg min',2)
      goto 103
      endif
      if(inp_compare(.false.,'max',item)) then
      if(.not.inp_f(dxmcg)) call md_abort('argos: error cg max',2)
      goto 103
      endif
      if(inp_compare(.false.,'cy',item)) then
      if(.not.inp_i(ncgcy)) call md_abort('argos: error cg cy',2)
      goto 103
      endif
      if(inp_compare(.false.,'criterion',item)) then
      if(.not.inp_i(icrit)) call md_abort('argos: error cg crit',2)
      goto 103
      endif
      endif
c
c     forward mcti
c
      if(inp_compare(.false.,'forward',item)) then
      iforw=1
      if(inp_i(mrun)) then
      maxlam=mrun
      if(.not.inp_a(item)) goto 1
      if(.not.inp_compare(.false.,'of',item)) goto 3
      if(.not.inp_i(maxlam)) call md_abort('argos: error mcti',7)
      endif
      goto 2
      endif
c
c
c     reverse mcti
c
      if(inp_compare(.false.,'reverse',item)) then
      iforw=-1
      if(inp_i(mrun)) then
      maxlam=mrun
      if(.not.inp_a(item)) goto 1
      if(.not.inp_compare(.false.,'of',item)) goto 3
      if(.not.inp_i(maxlam)) call md_abort('argos: error mcti',7)
      endif
      goto 2
      endif
c
c     error
c
      if(inp_compare(.false.,'error',item)) then
      if(.not.inp_f(edacq)) call md_abort('argos: error error',8)
      goto 2
      endif
c
c     drift
c
      if(inp_compare(.false.,'drift',item)) then
      if(.not.inp_f(ddacq)) call md_abort('argos: error drift',9)
      goto 2
      endif
c
c     factor
c
      if(inp_compare(.false.,'factor',item)) then
      if(.not.inp_f(fdacq)) call md_abort('argos: error factor',10)
      goto 2
      endif
c
c     decomposition
c
      if(inp_compare(.false.,'decomp',item)) then
      if(.not.inp_i(npgdec)) npgdec=1
      goto 2
      endif
c
c     separation shifted scaling
c
      if(inp_compare(.false.,'sss',item)) then
      issscl=1
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'delta',item)) then
      if(.not.inp_f(delta)) call md_abort('argos: error delta',11)
      goto 2
      endif
      endif
c
c     new mcti
c
      if(inp_compare(.false.,'new',item)) then
      mropt=1
      goto 2
      endif
c
c     renew mcti
c
      if(inp_compare(.false.,'renew',item)) then
      mropt=2
      goto 2
      endif
c
c     extend mcti
c
      if(inp_compare(.false.,'extend',item)) then
      mropt=3
      goto 2
      endif
c
c     input for autocorrelation function error analysis
c
      if(inp_compare(.false.,'auto',item)) then
      if(.not.inp_i(nacfl)) call md_abort('argos: error auto',3)
  202 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'fit',item)) then
      if(.not.inp_i(iapprx)) call md_abort('argos: error fit',3)
      if(iapprx.lt.0) iapprx=0
      if(iapprx.gt.15) iapprx=15
      goto 202
      endif
      if(inp_compare(.false.,'weight',item)) then
      if(.not.inp_f(weight)) call md_abort('argos: error weight',3)
      if(weight.lt.0.0) weight=0.0
      goto 202
      endif
      goto 3
      endif
c
c     input for md and mcti
c
      if(inp_compare(.false.,'time',item)) then
      if(.not.inp_f(stime)) call md_abort('argos: error start',3)
      goto 2
      endif
c
      if(inp_compare(.false.,'step',item)) then
      if(.not.inp_f(tstep)) call md_abort('argos: error step',3)
      goto 2
      endif
c
      if(inp_compare(.false.,'leapfrog',item)) then
      mdalgo=1
      goto 2
      endif
c
      if(inp_compare(.false.,'leapfrog_bc',item)) then
      mdalgo=2
      goto 2
      endif
c
      if(inp_compare(.false.,'equil',item)) then
      if(.not.inp_i(mequi)) call md_abort('argos: error equil',3)
      goto 2
      endif
c
      if(inp_compare(.false.,'data',item)) then
      if(.not.inp_i(mdacq)) call md_abort('argos: error data',4)
      ldacq=mdacq
      if(.not.inp_a(item)) goto 1
      if(.not.inp_compare(.false.,'over',item)) goto 3
      if(.not.inp_i(ldacq)) call md_abort('argos: error over',5)
      goto 2
      endif
c
c     constant pressure
c
c     isobar [<pres_ext>] [trelax <prsrlx>] [compress <compr>]
c
      if(inp_compare(.false.,'isobar',item)) then
      ipscal=1
      if(inp_f(prsext)) ipscal=1
  250 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trelax',item)) then
      if(.not.inp_f(prsrlx))
     + call md_abort('argos: error p_trelax',12)
      goto 250
      endif
      if(inp_compare(.false.,'compress',item)) then
      if(.not.inp_f(compr)) call md_abort('argos: error compress',13)
      goto 250
      endif
      if(inp_compare(.false.,'anisotropic',item)) then
      ipscal=2
      goto 250
      endif
      if(inp_compare(.false.,'xy',item)) then
      ipopt=12
      goto 250
      endif
      if(inp_compare(.false.,'z',item)) then
      ipopt=3
      goto 250
      endif
      if(inp_compare(.false.,'xy-z',item)) then
      ipopt=123
      goto 250
      endif
      goto 3
      endif
c
c     constant temperature
c
c     isotherm [<tmpext> [tmpext2>]] [trelax <tmprlx> [<tmsrlx>]] \
c              [anneal <tanneal1> [<tanneal2>]]
c
      if(inp_compare(.false.,'isotherm',item)) then
      itscal=1
      if(inp_f(tmpext)) itscal=1
      if(.not.inp_f(tmpext2)) tmpext2=tmpext
  252 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trelax',item)) then
      if(.not.inp_f(tmprlx))
     + call md_abort('argos: error t_trelax',14)
      if(inp_f(tmsrlx)) itscal=2
      goto 252
      endif
      if(inp_compare(.false.,'anneal',item)) then
      if(.not.inp_f(tann1))
     + call md_abort('argos: error t_tanneal',14)
      if(.not.inp_f(tann2)) then
      tann2=tann1
      tann1=0.0d0
      endif
      goto 252
      endif
      goto 3
      endif
c
c     qhop proton hopping
c
      if(inp_compare(.false.,'qhop',item)) then
      if(.not.inp_i(nfhop)) nfhop=10
      if(.not.inp_f(rhop)) rhop=0.35d0
      if(.not.inp_f(thop)) thop=0.02d0
      goto 2
      endif
c
c     guided dynamics
c
      if(inp_compare(.false.,'guided',item)) then
      iguide=1
      if(.not.inp_f(fguide)) fguide=0.2d0
      if(.not.inp_f(tguide)) tguide=0.2d0
      goto 2
      endif
c
c     membrane setup options
c
      if(inp_compare(.false.,'membrane',item)) then
      imembr=1
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'rotations',item)) then
      imembr=2
      goto 2
      endif
      goto 3
      endif
c
c     solute center of mass options
c
      if(inp_compare(.false.,'scmxy',item)) then
      if(.not.inp_i(icmopt)) icmopt=1
      goto 2
      endif
c
c     reaction field
c
      if(inp_compare(.false.,'react',item)) then
      ireact=1
      if(inp_f(dielec)) ireact=1
      goto 2
      endif
c
c     external field
c
      if(inp_compare(.false.,'field',item)) then
      if(.not.inp_f(xfield)) call md_abort('argos: error xfield',13)
      ifield=1
  260 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'freq',item)) then
      if(.not.inp_f(xffreq)) call md_abort('argos: error xffreq',13)
      goto 260
      endif
      if(inp_compare(.false.,'vector',item)) then
      if(.not.inp_f(xfvect(1)))
     + call md_abort('argos: error xfvector x',13)
      if(.not.inp_f(xfvect(2)))
     + call md_abort('argos: error xfvector y',13)
      if(.not.inp_f(xfvect(3)))
     + call md_abort('argos: error xfvector z',13)
      goto 260
      endif
      goto 3
      endif
c
c     pmf
c
      if(inp_compare(.false.,'pmf',item)) then
      npmf=1
      facpmf=1.0d0
  288 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'equilharm',item)) then
      if(.not.inp_i(npmf)) call md_abort('argos: error equilharm',0)
      if(npmf.le.1) call md_abort('argos: error equilharm',npmf)
      if(.not.inp_a(item)) goto 1
      elseif(inp_compare(.false.,'scale',item).or.
     + inp_compare(.false.,'scale_factor',item)) then
      if(.not.inp_f(facpmf))
     + call md_abort('argos: error pmf scale',0)
      endif
      goto 2
      endif
c
c     collapse
c
      if(inp_compare(.false.,'collapse',item)) then
      ncolt=0
      ncoll=3
      fcoll=1.0d+01
  278 continue
      if(.not.inp_f(fcoll)) then
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'segment',item)) then
      ncolt=3
      ncoll=6
      goto 278
      endif
      if(inp_compare(.false.,'z',item)) then
      ncoll=ncolt+1
      goto 278
      endif
      if(inp_compare(.false.,'xy',item)) then
      ncoll=ncolt+2
      goto 278
      endif
      endif
      goto 2
      endif
c
c     include
c
      if(inp_compare(.false.,'include',item)) then
      if(.not.inp_a(item)) call md_abort('argos: error include',0)
      if(inp_compare(.false.,'fixed',item)) then
      includ=1
      else
      call md_abort('argos: error include option',0)
      endif
      goto 2
      endif
c
c     distance restraints
c
      if(inp_compare(.false.,'distar',item)) then
      ndistr=1
  241 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'draver',item)) then
      ndaver=1
      if(inp_f(disrlx)) ndaver=1
      goto 241
      elseif(inp_compare(.false.,'scale',item)) then
      if(.not.inp_f(drsscl))
     + call md_abort('argos: error distar scale',13)
      goto 241
      elseif(inp_compare(.false.,'after',item)) then
      if(.not.inp_i(nfdrss))
     + call md_abort('argos: error distar scale',13)
      goto 241
      endif
      endif
c
c     development flag
c
      if(inp_compare(.false.,'develop',item)) then
      if(.not.inp_i(idevel)) call md_abort('argos: error develop',3)
      goto 2
      endif
c
c     particle-mesh Ewald
c
c     ewald [ grid <i ngx> [ <i ngy> [ <i ngz> ] ] ]
c
      if(inp_compare(.false.,'ewald',item).or.
     + inp_compare(.false.,'pme',item)) then
      ipme=1
  251 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'grid',item)) then
      if(.not.inp_i(ngx)) call md_abort('argos: error grid',13)
      ngy=ngx
      ngz=ngx
      if(inp_i(ngy)) then
      if(.not.inp_i(ngz)) ngz=ngy
      endif
      goto 251
      endif
      if(inp_compare(.false.,'alpha',item)) then
      if(.not.inp_f(ealpha)) call md_abort('argos: error alpha',13)
      goto 251
      endif
      if(inp_compare(.false.,'order',item)) then
      if(.not.inp_i(morder)) call md_abort('argos: error order',13)
      goto 251
      endif
      if(inp_compare(.false.,'fft',item)) then
      if(.not.inp_i(imfft)) call md_abort('argos: error fft',13)
      goto 251
      endif
      if(inp_compare(.false.,'procs',item)) then
      if(.not.inp_i(nodpme)) call md_abort('argos: error procs',13)
      goto 251
      endif
      if(inp_compare(.false.,'nodes',item)) then
      if(.not.inp_i(nodpme)) call md_abort('argos: error nodes',13)
      goto 251
      endif
      if(inp_compare(.false.,'solvent',item)) then
      isolvo=1
      goto 251
      endif
      goto 3
      endif
c
c     fix [solvent] [solute]
c
      if(inp_compare(.false.,'fix',item).or.
     + inp_compare(.false.,'free',item)) then
c
      ifix=1
      afix=' '
      rfix=0.0d0
      if(inp_compare(.false.,'free',item)) ifix=0
c
      if(.not.inp_a(item)) call md_abort('argos: error fix',1)
c
      if(inp_compare(.false.,'solvent',item)) then
      if(.not.inp_i(ifirst)) then
      if(.not.inp_a(item)) then
      ifirst=0
      ilast=0
      else
      ifirst=-1
      if(inp_compare(.false.,'within',item)) then
      if(.not.inp_f(rfix)) call md_abort('argos: error fix',2)
      if(.not.inp_a(afix)) call md_abort('argos: error fix',3)
      ilast=-1
      elseif(inp_compare(.false.,'beyond',item)) then
      if(.not.inp_f(rfix)) call md_abort('argos: error fix',4)
      if(.not.inp_a(afix)) call md_abort('argos: error fix',5)
      ilast=-2
      else
      call md_abort('argos: error fix',6)
      endif
      endif
      else
      if(.not.inp_i(ilast)) ilast=ifirst
      endif
      write(lfncmd,1001) ifix,ifirst,ilast,rfix,afix
 1001 format('fix w ',i2,2i10,f12.6,a)
      numfix=numfix+1
      goto 2
      endif
c
      if(inp_compare(.false.,'solute',item)) then
      if(.not.inp_i(ifirst)) then
      if(.not.inp_a(item)) then
      ifirst=0
      ilast=0
      else
      ifirst=-1
      if(inp_compare(.false.,'within',item)) then
      if(.not.inp_f(rfix)) call md_abort('argos: error fix',7)
      if(.not.inp_a(afix)) call md_abort('argos: error fix',8)
      ilast=-1
      elseif(inp_compare(.false.,'beyond',item)) then
      if(.not.inp_f(rfix)) call md_abort('argos: error fix',9)
      if(.not.inp_a(afix)) call md_abort('argos: error fix',10)
      ilast=-2
      else
      call md_abort('argos: error fix',11)
      endif
      endif
      else
      if(.not.inp_i(ilast)) ilast=ifirst
      endif
      if(.not.inp_a(item)) then
      write(lfncmd,1002) ifix,ifirst,ilast,rfix,afix
 1002 format('fix s ',i2,2i10,f12.6,a)
      numfix=numfix+1
      goto 2
      else
      if(inp_compare(.false.,'heavy',item)) then
      write(lfncmd,1003) ifix,ifirst,ilast
 1003 format('fix x ',i2,2i10)
      numfix=numfix+1
      goto 2
      else
 1113 continue
      write(lfncmd,1004) ifix,ifirst,ilast,item(1:index(item,' ')-1)
 1004 format('fix s ',i2,2i10,12x,a)
      numfix=numfix+1
      if(inp_a(item)) goto 1113
      goto 2
      endif
      endif
      endif
c
      if(inp_compare(.false.,'permanent',item)) then
      write(lfncmd,1005)
 1005 format('fix p ')
      numfix=numfix+1
      goto 2
      endif
c
      if(inp_compare(.false.,'select',item)) then
      if(.not.inp_a(item)) call md_abort('argos_input fix select',14)
      write(lfncmd,1006) ifix,0.0d0,item
 1006 format('fix S ',i2,20x,f12.6,a)
      numfix=numfix+1
      goto 2
      endif
c
      call md_abort('argos: unrecognized fix/free command',13)
      endif
c
c      if(inp_compare(.false.,'fix',item)) then
c  261 continue
c      if(.not.inp_a(item)) goto 1
c  361 continue
c      if(inp_compare(.false.,'all',item)) then
c      ifix=ior(ifix,1)
c      goto 261
c      endif
c      if(inp_compare(.false.,'solvent',item)) then
c      ifix=ior(ifix,2)
c      goto 261
c      endif
c      if(inp_compare(.false.,'solute',item)) then
c      ifix=ior(ifix,4)
c      if(.not.inp_a(item)) goto 1
c      if(inp_compare(.false.,'fraction',item)) then
cc      do 5 i=1,msf
cc      lsffix(i)=0
cc    5 continue
cc      do 6 i=1,msf
cc      if(.not.inp_i(n)) goto 261
cc      if(n.gt.0.and.n.le.msf) lsffix(n)=1
cc    6 continue
c      goto 261
c      endif
c      goto 361
c      endif
c      if(inp_compare(.false.,'non-H',item)) then
c      ifix=ior(ifix,8)
c      goto 261
c      endif
c      goto 3
c      endif
c
c     fix [solvent] [solute]
c
c      if(inp_compare(.false.,'unfix',item)) then
c  262 continue
c      if(.not.inp_a(item)) goto 1
c      if(inp_compare(.false.,'all',item)) then
c      iunfix=ior(iunfix,1)
c      goto 262
c      endif
c      if(inp_compare(.false.,'solvent',item)) then
c      iunfix=ior(iunfix,2)
c      goto 262
c      endif
c      if(inp_compare(.false.,'solute',item)) then
c      iunfix=ior(iunfix,4)
c      goto 262
c      endif
c      if(inp_compare(.false.,'non-H',item)) then
c      iunfix=ior(iunfix,8)
c      goto 262
c      endif
c      goto 3
c      endif
c
c     noshake [solvent] [solute]
c
      if(inp_compare(.false.,'noshake',item)) then
c      noshak=3
  263 continue
      if(.not.inp_a(item)) then
        if(noshak.eq.0) noshak=3
        goto 1
      end if
      if(inp_compare(.false.,'solvent',item)) then
      if(noshak.eq.0) noshak=1
      if(noshak.eq.2) noshak=3
      goto 263
      endif
      if(inp_compare(.false.,'solute',item)) then
      if(noshak.eq.0) noshak=2
      if(noshak.eq.1) noshak=3
      goto 263
      endif
      goto 3
      endif
c
c     hessian
c
      if(inp_compare(.false.,'profile',item)) then
      if(.not.inp_i(iprof)) iprof=1
      goto 2
      endif
c
c     hessian
c
      if(inp_compare(.false.,'hessian',item)) then
      ihess=1
      goto 2
      endif
c
c     radius of gyration
c
      if(inp_compare(.false.,'radius_gyration',item)) then
      iradgy=1
      goto 2
      endif
c
c     diffusion coefficients
c
      if(inp_compare(.false.,'diffusion',item)) then
      idifco=1
      goto 2
      endif
c
c     velocity reassignment
c
      if(inp_compare(.false.,'vreass',item).or.
     + inp_compare(.false.,'reassign',item)) then
      if(.not.inp_i(nfgaus)) call md_abort('argos: error vreass',13)
      if(nfgaus.lt.0) nfgaus=-1
      if(.not.inp_f(tgauss)) call md_abort('argos: error vreass',13)
  900 continue
      if(.not.inp_a(item)) goto 1
c
      if(inp_compare(.false.,'fraction',item)) then
      if(.not.inp_f(frgaus)) then
      frgaus=0.5
      goto 900
      endif
      if(frgaus.lt.0.0) frgaus=0.0
      if(frgaus.gt.1.0) frgaus=1.0
      goto 900
      endif
c
      if(inp_compare(.false.,'once',item)) then
      ivopt=ior(ivopt,1)
      goto 900
      endif
c
      if(inp_compare(.false.,'first',item).or.
     + inp_compare(.false.,'initial',item)) then
      ivopt=ior(ivopt,2)
      goto 900
      endif
c
      if(inp_compare(.false.,'last',item).or.
     + inp_compare(.false.,'final',item)) then
      ivopt=ior(ivopt,4)
      goto 900
      endif
c
      goto 2
      endif
c
c     QM zero energy
c
      if(inp_compare(.false.,'eqm',item)) then
      if(.not.inp_f(eqm)) call md_abort('argos: error eqm',0)
      goto 2
      endif
c
c     input for all calculations
c
      if(inp_compare(.false.,'cutoff',item)) then
  300 continue
      if(inp_f(rshort)) then
      rlong=rshort
      goto 300
      endif
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'short',item)) then
      if(.not.inp_f(rshort)) call md_abort('argos: error short',14)
      goto 300
      endif
      if(inp_compare(.false.,'long',item)) then
      if(.not.inp_f(rlong)) call md_abort('argos: error long',14)
      goto 300
      endif
      goto 3
      endif
c
      if(inp_compare(.false.,'atomlist',item)) then
      latom=1
      goto 2
      endif
c
      if(inp_compare(.false.,'polar',item)) then
      ipolt=1
  400 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'first',item)) then
      ipolt=1
      goto 400
      endif
      if(inp_compare(.false.,'scf',item)) then
      ipolt=3
      if(.not.inp_i(mpolit)) mpolit=10
      if(inp_f(ptol)) goto 400
      goto 400
      endif
      goto 3
      endif
c
      if(inp_compare(.false.,'shake',item)) then
      if(inp_i(mshitw)) then
      mshits=mshitw
      if(inp_i(mshits)) then
      endif
      endif
      if(inp_f(tlwsha)) then
      tlssha=tlwsha
      if(inp_f(tlssha)) then
      endif
      endif
      goto 2
      endif
c
c     communication option
c
      if(inp_compare(.false.,'comlim',item)) then
      if(.not.inp_a(item)) call md_abort('argos: error comlim',0)
      if(inp_compare(.false.,'on',item)) then
      icomm=1
      goto 2
      elseif(inp_compare(.false.,'off',item)) then
      icomm=0
      goto 2
      else
      call md_abort('argos: error comlim',0)
      endif
      endif
c
      if(inp_compare(.false.,'detail',item)) then
      if(.not.inp_i(npstat)) npstat=1
      goto 2
      endif
c
c     print flags
c
      if(inp_compare(.false.,'print',item)) then
  500 continue
      if(.not.inp_a(item)) goto 1
  501 continue
      if(inp_compare(.false.,'step',item)) then
      if(.not.inp_i(nfoutp)) call md_abort('argos: error step',20)
      goto 500
      endif
      if(inp_compare(.false.,'stat',item)) then
      if(.not.inp_i(nfstat)) call md_abort('argos: error stat',20)
      goto 500
      endif
      if(inp_compare(.false.,'dayout',item)) then
      itrack=1
      goto 500
      endif
      if(inp_compare(.false.,'out6',item)) then
      lfout6=1
      goto 500
      endif
      if(inp_compare(.false.,'pmf',item)) then
      if(.not.inp_i(iprpmf)) iprpmf=-1
      goto 500
      endif
      if(inp_compare(.false.,'extra',item)) then
      npstep=1
      goto 500
      endif
      if(inp_compare(.false.,'energy',item)) then
      if(.not.inp_i(npener)) npener=1
      goto 500
      endif
      if(inp_compare(.false.,'energies',item)) then
      if(.not.inp_i(npener)) npener=1
      goto 500
      endif
      if(inp_compare(.false.,'forces',item)) then
      if(.not.inp_i(npforc)) npforc=1
      goto 500
      endif
      if(inp_compare(.false.,'matrix',item)) then
      npdist=1
      goto 500
      endif
      if(inp_compare(.false.,'expect',item)) then
      if(.not.inp_i(npxpct)) call md_abort('argos: error expect',20)
      goto 500
      endif
      if(inp_compare(.false.,'timing',item)) then
      nptmng=ga_nnodes()-1
      goto 500
      endif
      if(inp_compare(.false.,'topol',item)) then
  502 continue
      if(.not.inp_a(item)) then
      nptopw=1
      nptops=1
      npatom=1
      goto 1
      endif
  503 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'solvent',item)) then
      nptopw=1
      goto 503
      endif
      if(inp_compare(.false.,'solute',item)) then
      nptops=1
      goto 503
      endif
      if(inp_compare(.false.,'nonbond',item)) then
      npatom=1
      goto 503
      endif
      goto 501
      endif
      goto 3
      endif
c
c     update
c
      if(inp_compare(.false.,'update',item)) then
  800 continue
      if(.not.inp_a(item)) goto 1
  801 continue
      if(inp_compare(.false.,'pairs',item)) then
      if(.not.inp_i(nfpair)) call md_abort('argos: error pair',20)
      goto 800
      endif
      if(inp_compare(.false.,'esp',item)) then
      if(.not.inp_i(nfesp)) call md_abort('argos: error esp',20)
      goto 800
      endif
      if(inp_compare(.false.,'long',item)) then
      if(.not.inp_i(nflong)) call md_abort('argos: error long',20)
      goto 800
      endif
      if(inp_compare(.false.,'analysis',item).or.
     + inp_compare(.false.,'analyze',item)) then
      if(.not.inp_i(nfanal)) call md_abort('argos: error anal',20)
      goto 800
      endif
      if(inp_compare(.false.,'center',item)) then
      if(.not.inp_i(nfcntr)) call md_abort('argos: error center',20)
      icentr=3
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'zonly',item)) then
      icentr=1
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'xyonly',item)) then
      icentr=2
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'fraction',item)) then
      nscb=0
  802 continue
      if(inp_i(itemp)) then
      nscb=nscb+1
      idscb(nscb)=itemp
      if(nscb.lt.10) goto 802
      goto 800
      endif
      goto 800
      else
      nscb=1
      idscb(1)=1
      goto 801
      endif
      goto 800
      endif
      if(inp_compare(.false.,'motion',item)) then
      if(.not.inp_i(nfslow)) call md_abort('argos: error motion',20)
      goto 800
      endif
      if(inp_compare(.false.,'rdf',item)) then
      if(.not.inp_i(nfrdf)) call md_abort('argos: error rdf',20)
  803 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'range',item)) then
      if(.not.inp_f(rrdf)) call md_abort('argos: error range',20)
      goto 803
      endif
      if(inp_compare(.false.,'bins',item)) then
      if(.not.inp_i(ngl)) call md_abort('argos: error bins',20)
      goto 803
      endif
      goto 801
      endif
      goto 2
      endif
c
c     recording
c
      if(inp_compare(.false.,'record',item)) then
  600 continue
      if(.not.inp_a(item)) goto 1
  601 continue
      if(inp_compare(.false.,'rest',item)) then
      if(.not.inp_i(nfrest)) call md_abort('argos: error restart',20)
      if(.not.inp_a(item)) goto 1
      if(.not.inp_compare(.false.,'keep',item)) goto 601
      keepr=1
      goto 600
      endif
      if(inp_compare(.false.,'acf',item)) then
      nfauto=1
      goto 600
      endif
      if(inp_compare(.false.,'cnv',item)) then
      nfconv=1
      goto 600
      endif
      if(inp_compare(.false.,'fet',item)) then
      nffet=1
      goto 600
      endif
      if(inp_compare(.false.,'precision',item)) then
      if(.not.inp_i(nprec)) call md_abort('argos: error coord',20)
      goto 600
      endif
      if(inp_compare(.false.,'coord',item)) then
      if(.not.inp_i(nfcoor)) call md_abort('argos: error coord',20)
      goto 600
      endif
      if(inp_compare(.false.,'wcoor',item)) then
      if(.not.inp_i(nfwcoo)) call md_abort('argos: error scoor',20)
      goto 600
      endif
      if(inp_compare(.false.,'scoor',item)) then
      if(.not.inp_i(nfscoo)) call md_abort('argos: error scoor',20)
      goto 600
      endif
      if(inp_compare(.false.,'veloc',item)) then
      if(.not.inp_i(nfvelo)) call md_abort('argos: error veloc',20)
      goto 600
      endif
      if(inp_compare(.false.,'wvelo',item)) then
      if(.not.inp_i(nfwvel)) call md_abort('argos: error svelo',20)
      goto 600
      endif
      if(inp_compare(.false.,'svelo',item)) then
      if(.not.inp_i(nfsvel)) call md_abort('argos: error svelo',20)
      goto 600
      endif
      if(inp_compare(.false.,'force',item)) then
      if(.not.inp_i(nfforc)) call md_abort('argos: error force',20)
      goto 600
      endif
      if(inp_compare(.false.,'wforc',item)) then
      if(.not.inp_i(nfwfor)) call md_abort('argos: error sforc',20)
      goto 600
      endif
      if(inp_compare(.false.,'sforc',item)) then
      if(.not.inp_i(nfsfor)) call md_abort('argos: error sforc',20)
      goto 600
      endif
      if(inp_compare(.false.,'prop',item)) then
      if(.not.inp_i(nfprop)) call md_abort('argos: error prop',20)
      iprop=0
      if(nfprop.lt.0) then
      iprop=1
      nfprop=-nfprop
      endif
      goto 600
      endif
      if(inp_compare(.false.,'prop_average',item)) then
      if(.not.inp_i(nfprop)) call md_abort('argos: error prop',20)
      iprop=1
      goto 600
      endif
      if(inp_compare(.false.,'free',item)) then
      if(.not.inp_i(nffree)) call md_abort('argos: error free',20)
      goto 600
      endif
      if(inp_compare(.false.,'sync',item)) then
      if(.not.inp_i(nfsync)) call md_abort('argos: error sync',20)
      goto 600
      endif
      if(inp_compare(.false.,'binary',item)) then
      ibinar=1
      goto 600
      endif
      if(inp_compare(.false.,'ascii',item)) then
      ibinar=0
      goto 600
      endif
      if(inp_compare(.false.,'ecce',item)) then
      iformt=1
      goto 600
      endif
      if(inp_compare(.false.,'argos',item)) then
      iformt=0
      goto 600
      endif
      if(inp_compare(.false.,'times',item)) then
      if(inp_i(nftime)) then
  602 continue
      if(inp_i(i)) then
      if(i.gt.0) then
      idtime=idtime+1
      else
      i=-i
      endif
      itime(i)=idtime
      goto 602
      endif
      goto 600
      endif
      endif
      goto 2
      endif
c
c     batch
c
      if(inp_compare(.false.,'batch',item)) then
      if(.not.inp_i(nfnewf)) call md_abort('argos: error batch',15)
      goto 2
      endif
c
c     importing
c
      if(inp_compare(.false.,'import',item)) then
      if(.not.inp_i(impfr)) then
      impfr=1
      impto=1
      nftri=1
      goto 2
      endif
      if(.not.inp_i(impto)) then
      impto=impfr
      impfr=1
      nftri=1
      goto 2
      endif
      if(.not.inp_i(nftri)) nftri=1
      goto 2
      endif
c
c     other stuff
c
      if(inp_compare(.false.,'memory',item)) then
      if(.not.inp_i(memlim)) call md_abort('argos: error memory',15)
      goto 2
      endif
      if(inp_compare(.false.,'extra',item)) then
      if(.not.inp_i(madbox)) call md_abort('argos: error extra',15)
      goto 2
      endif
      if(inp_compare(.false.,'mwm',item)) then
      if(.not.inp_i(mwmreq)) call md_abort('argos: error mwm',15)
      goto 2
      endif
      if(inp_compare(.false.,'msa',item)) then
      if(.not.inp_i(msareq)) call md_abort('argos: error msa',15)
      goto 2
      endif
      if(inp_compare(.false.,'mbb',item)) then
      if(.not.inp_i(mbbreq)) call md_abort('argos: error mbb',15)
      goto 2
      endif
      if(inp_compare(.false.,'mcells',item)) then
      if(.not.inp_i(mbbreq)) call md_abort('argos: error mcells',15)
      goto 2
      endif
      if(inp_compare(.false.,'load',item)) then
      loadb=0
      if(.not.inp_i(nfload)) nfload=1
  700 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'none',item)) then
      loadb=0
      goto 700
      endif
      if(inp_compare(.false.,'reset',item)) then
      ireset=1
      goto 700
      endif
      if(inp_compare(.false.,'last',item)) then
      itload=0
      goto 700
      endif
      if(inp_compare(.false.,'minimum',item)) then
      itload=1
      goto 700
      endif
      if(inp_compare(.false.,'average',item)) then
      itload=2
      goto 700
      endif
      if(inp_compare(.false.,'combination',item)) then
      itload=3
      goto 700
      endif
      if(inp_compare(.false.,'iotime',item)) then
      ioload=1
      goto 700
      endif
      if(inp_compare(.false.,'experimental',item)) then
      ioload=2
      goto 700
      endif
      if(inp_compare(.false.,'sizez',item)) then
      isload=1
      if(loadb.eq.0) loadb=1
      if(loadb.eq.2) loadb=3
      if(inp_f(factld)) goto 700
      goto 700
      endif
      if(inp_compare(.false.,'size',item)) then
      if(loadb.eq.0) loadb=1
      if(loadb.eq.2) loadb=3
      if(inp_f(factld)) goto 700
      goto 700
      endif
      if(inp_compare(.false.,'pairs',item)) then
      if(loadb.eq.0) loadb=2
      if(loadb.eq.1) loadb=3
      if(.not.inp_i(lbpair)) goto 700
      lbpair=min(lbpair,50)
      goto 700
      endif
      goto 2
      endif
      if(inp_compare(.false.,'debug',item)) then
      if(.not.inp_i(i)) call md_abort('argos: error debug',15)
  555 continue
      if(i.le.0) then
      idebug=0
      else
      idebug=ior(idebug,2**(i-1))
      endif
      if(inp_i(i)) goto 555
      goto 2
      endif
      if(inp_compare(.false.,'test',item)) then
      if(.not.inp_i(itest)) call md_abort('argos: error test',15)
      goto 2
      endif
      if(inp_compare(.false.,'nopack',item)) then
      nopack=1
      goto 2
      endif
      if(inp_compare(.false.,'pack',item)) then
      nopack=0
      goto 2
      endif
      if(inp_compare(.false.,'control',item)) then
      if(.not.inp_i(icntrl)) call md_abort('argos: error control',15)
      goto 2
      endif
      if(inp_compare(.false.,'scale',item)) then
      if(.not.inp_f(scaleq)) call md_abort('argos: error scale',15)
      goto 2
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,9999) item
 9999 format('argos: unknown field ',a)
      call md_abort('argos: error reading input file',0)
      endif
c
      if(.not.expert) then
      if(nflong.gt.1.and.ipme.gt.0)
     + call md_abort('argos: PME nflong>1 only in expert mode',0)
      endif
c
      if(nfscoo.lt.0) nfscoo=nfcoor
      if(nfsvel.lt.0) nfsvel=nfvelo
      if(nfsfor.lt.0) nfsfor=nfforc
      if(nfwcoo.lt.0) nfwcoo=nfcoor
      if(nfwvel.lt.0) nfwvel=nfvelo
      if(nfwfor.lt.0) nfwfor=nfforc
      if(nfscoo.lt.0) nfscoo=0
      if(nfsvel.lt.0) nfsvel=0
      if(nfsfor.lt.0) nfsfor=0
      if(nfwcoo.lt.0) nfwcoo=0
      if(nfwvel.lt.0) nfwvel=0
      if(nfwfor.lt.0) nfwfor=0
c
      numcmd=numfix
c
c     put input in rtdb
c
      if(project(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'md:project',1,project))
     + call md_abort('argos: rtdb_put failed',0)
      endif
      if(.not.rtdb_put(irtdb,'md:port',mt_int,1,iport))
     + call md_abort('argos: rtdb_put failed',0)
      if(iport.gt.0) then
      if(.not.rtdb_cput(irtdb,'md:server',1,server))
     + call md_abort('argos: rtdb_put failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'md:istart',mt_int,1,istart))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nbx',mt_int,1,nbx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nby',mt_int,1,nby))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nbz',mt_int,1,nbz))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npx',mt_int,1,npx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npy',mt_int,1,npy))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npz',mt_int,1,npz))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npg',mt_int,1,npg))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ncp',mt_int,1,ncp))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mdalgo',mt_int,1,mdalgo))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iset',mt_int,1,iset))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isetp1',mt_int,1,isetp1))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isetp2',mt_int,1,isetp2))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iforw',mt_int,1,iforw))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:msdit',mt_int,1,msdit))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mcgit',mt_int,1,mcgit))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icrit',mt_int,1,icrit))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ncgcy',mt_int,1,ncgcy))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mrun',mt_int,1,mrun))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:maxlam',mt_int,1,maxlam))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npgdec',mt_int,1,npgdec))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:issscl',mt_int,1,issscl))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mequi',mt_int,1,mequi))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mdacq',mt_int,1,mdacq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ldacq',mt_int,1,ldacq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nacfl',mt_int,1,nacfl))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iapprx',mt_int,1,iapprx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipscal',mt_int,1,ipscal))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipopt',mt_int,1,ipopt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itscal',mt_int,1,itscal))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfgaus',mt_int,1,nfgaus))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ivopt',mt_int,1,ivopt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipolt',mt_int,1,ipolt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mpolit',mt_int,1,mpolit))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mshitw',mt_int,1,mshitw))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mshits',mt_int,1,mshits))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itrack',mt_int,1,itrack))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npstep',mt_int,1,npstep))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npstat',mt_int,1,npstat))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npener',mt_int,1,npener))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npforc',mt_int,1,npforc))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npdist',mt_int,1,npdist))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptmng',mt_int,1,nptmng))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npatom',mt_int,1,npatom))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptopw',mt_int,1,nptopw))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptops',mt_int,1,nptops))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npxpct',mt_int,1,npxpct))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfpair',mt_int,1,nfpair))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfesp',mt_int,1,nfesp))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfrdf',mt_int,1,nfrdf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nflong',mt_int,1,nflong))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfcntr',mt_int,1,nfcntr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icentr',mt_int,1,icentr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfanal',mt_int,1,nfanal))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfslow',mt_int,1,nfslow))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfoutp',mt_int,1,nfoutp))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfstat',mt_int,1,nfstat))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfrest',mt_int,1,nfrest))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:keepr',mt_int,1,keepr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfcoor',mt_int,1,nfwcoo))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfscoo',mt_int,1,nfscoo))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfvelo',mt_int,1,nfwvel))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsvel',mt_int,1,nfsvel))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfforc',mt_int,1,nfwfor))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsfor',mt_int,1,nfsfor))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfprop',mt_int,1,nfprop))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iprop',mt_int,1,iprop))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nffree',mt_int,1,nffree))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsync',mt_int,1,nfsync))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfauto',mt_int,1,nfauto))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfconv',mt_int,1,nfconv))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nffet',mt_int,1,nffet))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:impfr',mt_int,1,impfr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:impto',mt_int,1,impto))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nftri',mt_int,1,nftri))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iformt',mt_int,1,iformt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:madbox',mt_int,1,madbox))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ireset',mt_int,1,ireset))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:loadb',mt_int,1,loadb))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mropt',mt_int,1,mropt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idebug',mt_int,1,idebug))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icntrl',mt_int,1,icntrl))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ifidi',mt_int,1,ifidi))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipbtyp',mt_int,1,ipbtyp))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngl',mt_int,1,ngl))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ifield',mt_int,1,ifield))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipme',mt_int,1,ipme))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngx',mt_int,1,ngx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngy',mt_int,1,ngy))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngz',mt_int,1,ngz))
     + call md_abort('argos: rtdb_put failed',0)
c      if(.not.rtdb_put(irtdb,'md:ifix',mt_int,1,ifix))
c     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:numfix',mt_int,1,numfix))
     + call md_abort('argos: rtdb_put failed',0)
c      if(.not.rtdb_put(irtdb,'md:lsffix',mt_int,msf,lsffix))
c     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:noshak',mt_int,1,noshak))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfefld',mt_int,1,nfefld))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsfld',mt_int,1,nfsfld))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nscb',mt_int,1,nscb))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idscb',mt_int,nscb,idscb))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ireact',mt_int,1,ireact))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:memlim',mt_int,1,memlim))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:morder',mt_int,1,morder))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isolvo',mt_int,1,isolvo))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:lfout6',mt_int,1,lfout6))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iprpmf',mt_int,1,iprpmf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:imfft',mt_int,1,imfft))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mwmreq',mt_int,1,mwmreq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:msareq',mt_int,1,msareq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mbbreq',mt_int,1,mbbreq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itest',mt_int,1,itest))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nodpme',mt_int,1,nodpme))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:lbpair',mt_int,1,lbpair))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ndistr',mt_int,1,ndistr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npmf',mt_int,1,npmf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:facpmf',mt_dbl,1,facpmf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfdrss',mt_int,1,nfdrss))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ndaver',mt_int,1,ndaver))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idevel',mt_int,1,idevel))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itime',mt_int,mtimes,itime))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nftime',mt_int,1,nftime))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfload',mt_int,1,nfload))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itload',mt_int,1,itload))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ioload',mt_int,1,ioload))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isload',mt_int,1,isload))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ihess',mt_int,1,ihess))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:latom',mt_int,1,latom))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icomm',mt_int,1,icomm))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfhop',mt_int,1,nfhop))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rhop',mt_dbl,1,rhop))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:thop',mt_dbl,1,thop))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iguide',mt_int,1,iguide))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icmopt',mt_int,1,icmopt))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:imembr',mt_int,1,imembr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nopack',mt_int,1,nopack))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iprof',mt_int,1,iprof))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:numcmd',mt_int,1,numcmd))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ncoll',mt_int,1,ncoll))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ilambd',mt_int,1,ilambd))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mlambd',mt_int,1,mlambd))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:includ',mt_int,1,includ))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iradgy',mt_int,1,iradgy))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idifco',mt_int,1,idifco))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfnewf',mt_int,1,nfnewf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nbget',mt_int,1,nbget))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nprec',mt_int,1,nprec))
     + call md_abort('argos: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'md:fguide',mt_dbl,1,fguide))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tguide',mt_dbl,1,tguide))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dx0sd',mt_dbl,1,dx0sd))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxsdmx',mt_dbl,1,dxsdmx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxmsd',mt_dbl,1,dxmsd))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dx0cg',mt_dbl,1,dx0cg))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxcgmx',mt_dbl,1,dxcgmx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxmcg',mt_dbl,1,dxmcg))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:edacq',mt_dbl,1,edacq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ddacq',mt_dbl,1,ddacq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:fdacq',mt_dbl,1,fdacq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:delta',mt_dbl,1,delta))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:stime',mt_dbl,1,stime))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tstep',mt_dbl,1,tstep))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:prsext',mt_dbl,1,prsext))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:prsrlx',mt_dbl,1,prsrlx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:compr',mt_dbl,1,compr))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmpext',mt_dbl,1,tmpext))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmpext2',mt_dbl,1,tmpext2))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmprlx',mt_dbl,1,tmprlx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmsrlx',mt_dbl,1,tmsrlx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tann1',mt_dbl,1,tann1))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tann2',mt_dbl,1,tann2))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tgauss',mt_dbl,1,tgauss))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:frgaus',mt_dbl,1,frgaus))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rlong',mt_dbl,1,rlong))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rshort',mt_dbl,1,rshort))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ptol',mt_dbl,1,ptol))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tlwsha',mt_dbl,1,tlwsha))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tlssha',mt_dbl,1,tlssha))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:factld',mt_dbl,1,factld))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rrdf',mt_dbl,1,rrdf))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xfield',mt_dbl,1,xfield))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xffreq',mt_dbl,1,xffreq))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xfvect',mt_dbl,3,xfvect))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:weight',mt_dbl,1,weight))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dielec',mt_dbl,1,dielec))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ealpha',mt_dbl,1,ealpha))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rbox',mt_dbl,1,rbox))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rsgm',mt_dbl,1,rsgm))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:disrlx',mt_dbl,1,disrlx))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:drsscl',mt_dbl,1,drsscl))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:uqmatm',mt_dbl,1,eqm))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:fcoll',mt_dbl,1,fcoll))
     + call md_abort('argos: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:scaleq',mt_dbl,1,scaleq))
     + call md_abort('argos: rtdb_put failed',0)
c
      close(unit=lfncmd)
c
      return
c
 1999 continue
      call md_abort('argos: unable to open cmd file',13)
c
      return
      end
