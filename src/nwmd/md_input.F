      subroutine md_input(irtdb)
c
c $Id: md_input.F,v 1.33 2000-10-06 19:43:55 d3j191 Exp $
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
      integer nflong,nfgaus,mropt,nfauto,nfconv,nffet,idebug,mpolit
      integer npatom,nptopw,nptops,npx,npy,npz,nbx,nby,nbz,itrack
      integer ngx,ngy,ngz,numfix,ifix,noshak,nfslow,nfefld,nfsfld
      integer ngl,isetp1,isetp2,ifield,iapprx,nacfl,npxpct,imfft
      integer ireact,itemp,nscb,idscb(10),memlim,morder,isolvo,lfout6
      integer mwmreq,msareq,itest,nodpme,lbpair,ireset,npmf,nopack
      integer ndistr,ndaver,idevel,itime(500),nftime,idtime,nfdrss
      integer itload,ihess,latom,icomm,nfanal,iport,iguide
      integer mbbreq,nfwcoo,nfwvel,ioload,iprof,isload,icntrl,istart
      real*8 dx0sd,dxsdmx,dxmsd,dx0cg,dxcgmx,dxmcg,edacq,fdacq,ddacq
      real*8 delta,stime,tstep,prsext,prsrlx,compr,tmpext,tmprlx
      real*8 tmsrlx,tgauss,rshort,rlong,ptol,tlwsha,tlssha,factld
      real*8 rrdf,xfield,xfvect(3),xffreq,weight,rqmmm,dielec,ealpha
      real*8 frgaus,rbox,disrlx,drsscl,eqm,tguide,fguide
c
      integer i,mtimes,mdold,len
      character*255 filcmd
      integer lfncmd,numcmd,ifirst,ilast
c
#include "bitops.fh"
c
      mdold=0
      if(.not.rtdb_get(irtdb,'md:mdold',mt_int,1,mdold)) mdold=0
      if(mdold.ne.0) then
      call arg_input(irtdb)
      return
      endif
c
      lfncmd=69
      numcmd=0
c
      mtimes=250
c
c     set the defaults
c
      nbx=0
      nby=0
      nbz=0
      npx=0
      npy=0
      npz=0
      itrack=0
      mdalgo=1
      iset=1
      isetp1=0
      isetp2=0
      iforw=1
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
      itscal=0
      nfgaus=0
      ipolt=0
      mpolit=10
      mshitw=100
      mshits=100
      npstep=0
      npatom=0
      nptopw=0
      nptops=0
      npener=0
      npdist=0
      nptmng=0
      nfpair=1
      nfrdf=0
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
      nfprop=0
      nffree=1
      nfsync=0
      nfauto=0
      nfconv=0
      nffet=0
      nfefld=0
      nfsfld=0
      ibinar=0
      iformt=1
      madbox=6
      loadb=0
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
      imfft=1
      mwmreq=0
      msareq=0
      mbbreq=0
      itest=0
      nodpme=0
      lbpair=5
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
      nopack=0
      istart=0
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
      tmprlx=0.1d0
      tmsrlx=0.1d0
      tgauss=298.15d0
      rlong=0.9d0
      rshort=0.9d0
      rqmmm=0.9d0
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
c      do 4 i=1,msf
c      lsffix(i)=1
c    4 continue
      do 7 i=1,mtimes
      itime(i)=0
    7 continue
      idtime=0
      project=' '
c
c     read the nwchem input file
c
    1 continue
c
      if(.not.inp_read()) call errquit('md_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
      if(inp_compare(.false.,'system',item)) then
      if(len.gt.0) call errquit('md_input: duplicate system',0)
      if(.not.inp_a(project)) call errquit('md_input: error system',0)
      len=index(project,' ')-1
      filcmd=project(1:len)//'.cmd '
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='unknown',err=1999)
      rewind(unit=lfncmd)
      numcmd=0
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
      if(inp_compare(.false.,'server',item)) then
      if(.not.inp_a(server)) call errquit('md_input: error server',0)
      if(.not.inp_i(iport)) call errquit('md_input: error port',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'boxmin',item)) then
      if(.not.inp_f(rbox)) call errquit('md_input: error minbox',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'boxes',item)) then
      if(.not.inp_i(nbx)) call errquit('md_input: error nbx',0)
      if(.not.inp_i(nby)) call errquit('md_input: error nby',0)
      if(.not.inp_i(nbz)) call errquit('md_input: error nbz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'nodes',item)) then
      if(.not.inp_i(npx)) call errquit('md_input: error npx',0)
      if(.not.inp_i(npy)) call errquit('md_input: error npy',0)
      if(.not.inp_i(npz)) call errquit('md_input: error npz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'procs',item)) then
      if(.not.inp_i(npx)) call errquit('md_input: error npx',0)
      if(.not.inp_i(npy)) call errquit('md_input: error npy',0)
      if(.not.inp_i(npz)) call errquit('md_input: error npz',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(iset)) call errquit('md_input: error set',2)
      if(iset.le.0) call errquit('md_input: illegal iset',2)
      if(iset.gt.3) call errquit('md_input: illegal iset',2)
      goto 2
      endif
c
      if(inp_compare(.false.,'pset',item)) then
      if(.not.inp_i(isetp1)) call errquit('md_input: error pset',2)
      if(isetp1.le.1) call errquit('md_input: illegal ipset1',2)
      if(isetp1.gt.3) call errquit('md_input: illegal ipset1',2)
      if(.not.inp_i(isetp2)) goto 2
      if(isetp1.eq.3) call errquit('md_input: illegal ipset1',2)
      if(isetp2.ne.3) call errquit('md_input: illegal ipset2',2)
      goto 2
      endif
c
c     energy minimization
c
c     em [[sd] [<i>] [init <r>] [max <i>] [min <i>]]
c     [cg [<i>] [cy <i>] [init <r>] [max <r>] [min <i>]]
c
      if(inp_compare(.false.,'sd',item)) then
      if(.not.inp_i(msdit)) call  errquit('md_input: error em sd',2)
  101 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'init',item)) then
      if(.not.inp_f(dx0sd)) call errquit('md_input: error sd init',2)
      goto 101
      endif
      if(inp_compare(.false.,'min',item)) then
      if(.not.inp_f(dxsdmx)) call errquit('md_input: error sd min',2)
      goto 101
      endif
      if(inp_compare(.false.,'max',item)) then
      if(.not.inp_f(dxmsd)) call errquit('md_input: error sd max',2)
      goto 101
      endif
      endif
c
c     conjugate gradient
c
      if(inp_compare(.false.,'cg',item)) then
      if(.not.inp_i(mcgit)) call  errquit('md_input: error em cg',2)
  103 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'init',item)) then
      if(.not.inp_f(dx0cg)) call errquit('md_input: error cg init',2)
      goto 103
      endif
      if(inp_compare(.false.,'min',item)) then
      if(.not.inp_f(dxcgmx)) call errquit('md_input: error cg min',2)
      goto 103
      endif
      if(inp_compare(.false.,'max',item)) then
      if(.not.inp_f(dxmcg)) call errquit('md_input: error cg max',2)
      goto 103
      endif
      if(inp_compare(.false.,'cy',item)) then
      if(.not.inp_i(ncgcy)) call errquit('md_input: error cg cy',2)
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
      if(.not.inp_i(maxlam)) call errquit('md_input: error mcti',7)
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
      if(.not.inp_i(maxlam)) call errquit('md_input: error mcti',7)
      endif
      goto 2
      endif
c
c     error
c
      if(inp_compare(.false.,'error',item)) then
      if(.not.inp_f(edacq)) call errquit('md_input: error error',8)
      goto 2
      endif
c
c     drift
c
      if(inp_compare(.false.,'drift',item)) then
      if(.not.inp_f(ddacq)) call errquit('md_input: error drift',9)
      goto 2
      endif
c
c     factor
c
      if(inp_compare(.false.,'factor',item)) then
      if(.not.inp_f(fdacq)) call errquit('md_input: error factor',10)
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
      if(.not.inp_f(delta)) call errquit('md_input: error delta',11)
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
      if(.not.inp_i(nacfl)) call errquit('md_input: error auto',3)
  202 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'fit',item)) then
      if(.not.inp_i(iapprx)) call errquit('md_input: error fit',3)
      if(iapprx.lt.0) iapprx=0
      if(iapprx.gt.15) iapprx=15
      goto 202
      endif
      if(inp_compare(.false.,'weight',item)) then
      if(.not.inp_f(weight)) call errquit('md_input: error weight',3)
      if(weight.lt.0.0) weight=0.0
      goto 202
      endif
      goto 3
      endif
c
c     input for md and mcti
c
      if(inp_compare(.false.,'time',item)) then
      if(.not.inp_f(stime)) call errquit('md_input: error start',3)
      goto 2
      endif
c
      if(inp_compare(.false.,'step',item)) then
      if(.not.inp_f(tstep)) call errquit('md_input: error step',3)
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
      if(.not.inp_i(mequi)) call errquit('md_input: error equil',3)
      goto 2
      endif
c
      if(inp_compare(.false.,'data',item)) then
      if(.not.inp_i(mdacq)) call errquit('md_input: error data',4)
      ldacq=mdacq
      if(.not.inp_a(item)) goto 1
      if(.not.inp_compare(.false.,'over',item)) goto 3
      if(.not.inp_i(ldacq)) call errquit('md_input: error over',5)
      goto 2
      endif
c
c     constant pressure
c
c     isobar [<pres_ext>] [trelax <prsrlx>] [compres <compr>]
c
      if(inp_compare(.false.,'isobar',item)) then
      ipscal=1
      if(inp_f(prsext)) ipscal=1
  250 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trelax',item)) then
      if(.not.inp_f(prsrlx))
     + call errquit('md_input: error p_trelax',12)
      goto 250
      endif
      if(inp_compare(.false.,'compres',item)) then
      if(.not.inp_f(compr)) call errquit('md_input: error compres',13)
      goto 250
      endif
      if(inp_compare(.false.,'anisotropic',item)) then
      ipscal=2
      goto 250
      endif
      if(inp_compare(.false.,'noz',item)) then
      ipscal=-iabs(ipscal)
      goto 250
      endif
      goto 3
      endif
c
c     constant temperature
c
c     isotherm [<temp_ext>] [trelax <tmprlx> [<tmsrlx>]]
c
      if(inp_compare(.false.,'isotherm',item)) then
      itscal=1
      if(inp_f(tmpext)) itscal=1
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trelax',item)) then
      if(.not.inp_f(tmprlx))
     + call errquit('md_input: error t_trelax',14)
      if(inp_f(tmsrlx)) itscal=2
      goto 2
      endif
      goto 3
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
      if(.not.inp_f(xfield)) call errquit('md_input: error xfield',13)
      ifield=1
  260 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'freq',item)) then
      if(.not.inp_f(xffreq)) call errquit('md_input: error xffreq',13)
      goto 260
      endif
      if(inp_compare(.false.,'vector',item)) then
      if(.not.inp_f(xfvect(1)))
     + call errquit('md_input: error xfvector x',13)
      if(.not.inp_f(xfvect(2)))
     + call errquit('md_input: error xfvector y',13)
      if(.not.inp_f(xfvect(3)))
     + call errquit('md_input: error xfvector z',13)
      goto 260
      endif
      goto 3
      endif
c
c     pmf
c
      if(inp_compare(.false.,'pmf',item)) then
      npmf=1
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'constraint',item)) then
      npmf=2
      if(.not.inp_a(item)) goto 1
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
     + call errquit('md_input: error distar scale',13)
      goto 241
      elseif(inp_compare(.false.,'after',item)) then
      if(.not.inp_i(nfdrss))
     + call errquit('md_input: error distar scale',13)
      goto 241
      endif
      endif
c
c     development flag
c
      if(inp_compare(.false.,'develop',item)) then
      if(.not.inp_i(idevel)) call errquit('md_input: error develop',3)
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
      if(.not.inp_i(ngx)) call errquit('md_input : error grid',13)
      ngy=ngx
      ngz=ngx
      if(inp_i(ngy)) then
      if(.not.inp_i(ngz)) ngz=ngy
      endif
      goto 251
      endif
      if(inp_compare(.false.,'alpha',item)) then
      if(.not.inp_f(ealpha)) call errquit('md_input : error alpha',13)
      goto 251
      endif
      if(inp_compare(.false.,'order',item)) then
      if(.not.inp_i(morder)) call errquit('md_input : error order',13)
      goto 251
      endif
      if(inp_compare(.false.,'fft',item)) then
      if(.not.inp_i(imfft)) call errquit('md_input : error fft',13)
      goto 251
      endif
      if(inp_compare(.false.,'procs',item)) then
      if(.not.inp_i(nodpme)) call errquit('md_input : error procs',13)
      goto 251
      endif
      if(inp_compare(.false.,'nodes',item)) then
      if(.not.inp_i(nodpme)) call errquit('md_input : error nodes',13)
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
      if(inp_compare(.false.,'free',item)) ifix=0
c
      if(.not.inp_a(item)) call errquit('md_input : error fix',13)
c
      if(inp_compare(.false.,'solvent',item)) then
      if(.not.inp_i(ifirst)) then
      ifirst=0
      ilast=0
      else
      if(.not.inp_i(ilast)) ilast=ifirst
      endif
      write(lfncmd,1001) ifix,ifirst,ilast
 1001 format('fix w ',i2,2i10)
      numfix=numfix+1
      goto 2
      endif
c
      if(inp_compare(.false.,'solute',item)) then
      if(.not.inp_i(ifirst)) then
      ifirst=0
      ilast=0
      else
      if(.not.inp_i(ilast)) ilast=ifirst
      endif
      if(.not.inp_a(item)) then
      write(lfncmd,1002) ifix,ifirst,ilast
 1002 format('fix s ',i2,2i10)
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
 1004 format('fix s ',i2,2i10,a)
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
      call errquit('md_input : unrecognized fix/free command',13)
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
      noshak=3
  263 continue
      if(.not.inp_a(item)) goto 1
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
c     velocity reassignment
c
      if(inp_compare(.false.,'vreass',item)) then
      if(.not.inp_i(nfgaus)) call errquit('md_input: error vreass',13)
      if(.not.inp_f(tgauss)) call errquit('md_input: error vreass',13)
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'fraction',item)) then
      if(.not.inp_f(frgaus)) then
      frgaus=0.5
      goto 2
      endif
      if(frgaus.lt.0.0) frgaus=0.0
      if(frgaus.gt.1.0) frgaus=1.0
      else
      goto 3
      endif
      goto 2
      endif
c
c     QM zero energy
c
      if(inp_compare(.false.,'eqm',item)) then
      if(.not.inp_f(eqm)) call errquit('md_input: error eqm',0)
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
      if(.not.inp_f(rshort)) call errquit('md_input: error short',14)
      goto 300
      endif
      if(inp_compare(.false.,'long',item)) then
      if(.not.inp_f(rlong)) call errquit('md_input: error long',14)
      goto 300
      endif
      if(inp_compare(.false.,'qmmm',item)) then
      if(.not.inp_f(rqmmm)) call errquit('md_input: error qmmm',14)
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
      if(.not.inp_a(item)) call errquit('md_input: error comlim',0)
      if(inp_compare(.false.,'on',item)) then
      icomm=1
      goto 2
      elseif(inp_compare(.false.,'off',item)) then
      icomm=0
      goto 2
      else
      call errquit('md_input: error comlim',0)
      endif
      endif
c
c     print flags
c
      if(inp_compare(.false.,'print',item)) then
  500 continue
      if(.not.inp_a(item)) goto 1
  501 continue
      if(inp_compare(.false.,'step',item)) then
      if(.not.inp_i(nfoutp)) call errquit('md_input: error step',20)
      goto 500
      endif
      if(inp_compare(.false.,'stat',item)) then
      if(.not.inp_i(nfstat)) call errquit('md_input: error stat',20)
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
      if(inp_compare(.false.,'extra',item)) then
      npstep=1
      goto 500
      endif
      if(inp_compare(.false.,'energy',item)) then
      if(.not.inp_i(npener)) npener=1
      goto 500
      endif
      if(inp_compare(.false.,'matrix',item)) then
      npdist=1
      goto 500
      endif
      if(inp_compare(.false.,'expect',item)) then
      if(.not.inp_i(npxpct)) call errquit('md_input: error expect',20)
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
      if(inp_compare(.false.,'solvent',item)) then
      nptopw=1
      goto 502
      endif
      if(inp_compare(.false.,'solute',item)) then
      nptops=1
      goto 502
      endif
      if(inp_compare(.false.,'nonbond',item)) then
      npatom=1
      goto 502
      endif
      nptopw=1
      nptops=1
      npatom=1
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
      if(.not.inp_i(nfpair)) call errquit('md_input: error pair',20)
      goto 800
      endif
      if(inp_compare(.false.,'long',item)) then
      if(.not.inp_i(nflong)) call errquit('md_input: error long',20)
      goto 800
      endif
      if(inp_compare(.false.,'analysis',item).or.
     + inp_compare(.false.,'analyze',item)) then
      if(.not.inp_i(nfanal)) call errquit('md_input: error anal',20)
      goto 800
      endif
      if(inp_compare(.false.,'center',item)) then
      if(.not.inp_i(nfcntr)) call errquit('md_input: error center',20)
      if(.not.inp_a(item)) goto 1
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
      if(.not.inp_i(nfslow)) call errquit('md_input: error motion',20)
      goto 800
      endif
      if(inp_compare(.false.,'rdf',item)) then
      if(.not.inp_i(nfrdf)) call errquit('md_input: error rdf',20)
  803 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'range',item)) then
      if(.not.inp_f(rrdf)) call errquit('md_input: error range',20)
      goto 803
      endif
      if(inp_compare(.false.,'bins',item)) then
      if(.not.inp_i(ngl)) call errquit('md_input: error bins',20)
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
      if(.not.inp_i(nfrest)) call errquit('md_input: error restart',20)
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
      if(inp_compare(.false.,'coord',item)) then
      if(.not.inp_i(nfcoor)) call errquit('md_input: error coord',20)
      goto 600
      endif
      if(inp_compare(.false.,'wcoor',item)) then
      if(.not.inp_i(nfwcoo)) call errquit('md_input: error scoor',20)
      goto 600
      endif
      if(inp_compare(.false.,'scoor',item)) then
      if(.not.inp_i(nfscoo)) call errquit('md_input: error scoor',20)
      goto 600
      endif
      if(inp_compare(.false.,'veloc',item)) then
      if(.not.inp_i(nfvelo)) call errquit('md_input: error veloc',20)
      goto 600
      endif
      if(inp_compare(.false.,'wvelo',item)) then
      if(.not.inp_i(nfwvel)) call errquit('md_input: error svelo',20)
      goto 600
      endif
      if(inp_compare(.false.,'svelo',item)) then
      if(.not.inp_i(nfsvel)) call errquit('md_input: error svelo',20)
      goto 600
      endif
      if(inp_compare(.false.,'prop',item)) then
      if(.not.inp_i(nfprop)) call errquit('md_input: error prop',20)
      goto 600
      endif
      if(inp_compare(.false.,'free',item)) then
      if(.not.inp_i(nffree)) call errquit('md_input: error free',20)
      goto 600
      endif
      if(inp_compare(.false.,'sync',item)) then
      if(.not.inp_i(nfsync)) call errquit('md_input: error sync',20)
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
c     other stuff
c
      if(inp_compare(.false.,'memory',item)) then
      if(.not.inp_i(memlim)) call errquit('md_input: error memory',15)
      goto 2
      endif
      if(inp_compare(.false.,'extra',item)) then
      if(.not.inp_i(madbox)) call errquit('md_input: error extra',15)
      goto 2
      endif
      if(inp_compare(.false.,'mwm',item)) then
      if(.not.inp_i(mwmreq)) call errquit('md_input: error mwm',15)
      goto 2
      endif
      if(inp_compare(.false.,'msa',item)) then
      if(.not.inp_i(msareq)) call errquit('md_input: error msa',15)
      goto 2
      endif
      if(inp_compare(.false.,'mbb',item)) then
      if(.not.inp_i(mbbreq)) call errquit('md_input: error mbb',15)
      goto 2
      endif
      if(inp_compare(.false.,'load',item)) then
      loadb=0
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
      if(inp_i(lbpair)) goto 700
      goto 700
      endif
      goto 2
      endif
      if(inp_compare(.false.,'debug',item)) then
      if(.not.inp_i(idebug)) call errquit('md_input: error debug',15)
      goto 2
      endif
      if(inp_compare(.false.,'test',item)) then
      if(.not.inp_i(itest)) call errquit('md_input: error test',15)
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
      if(.not.inp_i(icntrl)) call errquit('md_input: error control',15)
      goto 2
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,9999) item
 9999 format('md_input: unknown field ',a)
      call errquit('md_input: error reading input file',0)
      endif
c
      if(nfscoo.lt.0) nfscoo=nfcoor
      if(nfsvel.lt.0) nfsvel=nfvelo
      if(nfwcoo.lt.0) nfwcoo=nfcoor
      if(nfwvel.lt.0) nfwvel=nfvelo
      if(nfscoo.lt.0) nfscoo=0
      if(nfsvel.lt.0) nfsvel=0
      if(nfwcoo.lt.0) nfwcoo=0
      if(nfwvel.lt.0) nfwvel=0
c
      numcmd=numfix
c
c     put input in rtdb
c
      if(project(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'md:project',1,project))
     + call errquit('md_input: rtdb_put failed',0)
      endif
      if(.not.rtdb_put(irtdb,'md:port',mt_int,1,iport))
     + call errquit('md_input: rtdb_put failed',0)
      if(iport.gt.0) then
      if(.not.rtdb_cput(irtdb,'md:server',1,server))
     + call errquit('md_input: rtdb_put failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'md:istart',mt_int,1,istart))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nbx',mt_int,1,nbx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nby',mt_int,1,nby))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nbz',mt_int,1,nbz))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npx',mt_int,1,npx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npy',mt_int,1,npy))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npz',mt_int,1,npz))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mdalgo',mt_int,1,mdalgo))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iset',mt_int,1,iset))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isetp1',mt_int,1,isetp1))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isetp2',mt_int,1,isetp2))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iforw',mt_int,1,iforw))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:msdit',mt_int,1,msdit))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mcgit',mt_int,1,mcgit))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ncgcy',mt_int,1,ncgcy))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mrun',mt_int,1,mrun))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:maxlam',mt_int,1,maxlam))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npgdec',mt_int,1,npgdec))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:issscl',mt_int,1,issscl))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mequi',mt_int,1,mequi))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mdacq',mt_int,1,mdacq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ldacq',mt_int,1,ldacq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nacfl',mt_int,1,nacfl))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iapprx',mt_int,1,iapprx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipscal',mt_int,1,ipscal))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itscal',mt_int,1,itscal))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfgaus',mt_int,1,nfgaus))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipolt',mt_int,1,ipolt))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mpolit',mt_int,1,mpolit))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mshitw',mt_int,1,mshitw))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mshits',mt_int,1,mshits))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itrack',mt_int,1,itrack))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npstep',mt_int,1,npstep))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npener',mt_int,1,npener))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npdist',mt_int,1,npdist))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptmng',mt_int,1,nptmng))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npatom',mt_int,1,npatom))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptopw',mt_int,1,nptopw))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nptops',mt_int,1,nptops))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npxpct',mt_int,1,npxpct))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfpair',mt_int,1,nfpair))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfrdf',mt_int,1,nfrdf))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nflong',mt_int,1,nflong))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfcntr',mt_int,1,nfcntr))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfanal',mt_int,1,nfanal))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfslow',mt_int,1,nfslow))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfoutp',mt_int,1,nfoutp))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfstat',mt_int,1,nfstat))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfrest',mt_int,1,nfrest))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:keepr',mt_int,1,keepr))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfcoor',mt_int,1,nfwcoo))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfscoo',mt_int,1,nfscoo))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfvelo',mt_int,1,nfwvel))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsvel',mt_int,1,nfsvel))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfprop',mt_int,1,nfprop))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nffree',mt_int,1,nffree))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsync',mt_int,1,nfsync))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfauto',mt_int,1,nfauto))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfconv',mt_int,1,nfconv))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nffet',mt_int,1,nffet))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iformt',mt_int,1,iformt))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:madbox',mt_int,1,madbox))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ireset',mt_int,1,ireset))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:loadb',mt_int,1,loadb))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mropt',mt_int,1,mropt))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idebug',mt_int,1,idebug))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icntrl',mt_int,1,icntrl))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngl',mt_int,1,ngl))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ifield',mt_int,1,ifield))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ipme',mt_int,1,ipme))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngx',mt_int,1,ngx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngy',mt_int,1,ngy))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ngz',mt_int,1,ngz))
     + call errquit('md_input: rtdb_put failed',0)
c      if(.not.rtdb_put(irtdb,'md:ifix',mt_int,1,ifix))
c     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:numfix',mt_int,1,numfix))
     + call errquit('md_input: rtdb_put failed',0)
c      if(.not.rtdb_put(irtdb,'md:lsffix',mt_int,msf,lsffix))
c     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:noshak',mt_int,1,noshak))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfefld',mt_int,1,nfefld))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfsfld',mt_int,1,nfsfld))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nscb',mt_int,1,nscb))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idscb',mt_int,nscb,idscb))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ireact',mt_int,1,ireact))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:memlim',mt_int,1,memlim))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:morder',mt_int,1,morder))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isolvo',mt_int,1,isolvo))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:lfout6',mt_int,1,lfout6))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:imfft',mt_int,1,imfft))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mwmreq',mt_int,1,mwmreq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:msareq',mt_int,1,msareq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:mbbreq',mt_int,1,mbbreq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itest',mt_int,1,itest))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nodpme',mt_int,1,nodpme))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:lbpair',mt_int,1,lbpair))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ndistr',mt_int,1,ndistr))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:npmf',mt_int,1,npmf))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nfdrss',mt_int,1,nfdrss))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ndaver',mt_int,1,ndaver))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:idevel',mt_int,1,idevel))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itime',mt_int,mtimes,itime))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nftime',mt_int,1,nftime))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:itload',mt_int,1,itload))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ioload',mt_int,1,ioload))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:isload',mt_int,1,isload))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ihess',mt_int,1,ihess))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:latom',mt_int,1,latom))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:icomm',mt_int,1,icomm))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iguide',mt_int,1,iguide))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:nopack',mt_int,1,nopack))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:iprof',mt_int,1,iprof))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:numcmd',mt_int,1,numcmd))
     + call errquit('md_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'md:fguide',mt_dbl,1,fguide))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tguide',mt_dbl,1,tguide))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dx0sd',mt_dbl,1,dx0sd))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxsdmx',mt_dbl,1,dxsdmx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxmsd',mt_dbl,1,dxmsd))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dx0cg',mt_dbl,1,dx0cg))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxcgmx',mt_dbl,1,dxcgmx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dxmcg',mt_dbl,1,dxmcg))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:edacq',mt_dbl,1,edacq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ddacq',mt_dbl,1,ddacq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:fdacq',mt_dbl,1,fdacq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:delta',mt_dbl,1,delta))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:stime',mt_dbl,1,stime))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tstep',mt_dbl,1,tstep))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:prsext',mt_dbl,1,prsext))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:prsrlx',mt_dbl,1,prsrlx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:compr',mt_dbl,1,compr))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmpext',mt_dbl,1,tmpext))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmprlx',mt_dbl,1,tmprlx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tmsrlx',mt_dbl,1,tmsrlx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tgauss',mt_dbl,1,tgauss))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:frgaus',mt_dbl,1,frgaus))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rlong',mt_dbl,1,rlong))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rshort',mt_dbl,1,rshort))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rqmmm',mt_dbl,1,rqmmm))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ptol',mt_dbl,1,ptol))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tlwsha',mt_dbl,1,tlwsha))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:tlssha',mt_dbl,1,tlssha))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:factld',mt_dbl,1,factld))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rrdf',mt_dbl,1,rrdf))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xfield',mt_dbl,1,xfield))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xffreq',mt_dbl,1,xffreq))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:xfvect',mt_dbl,3,xfvect))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:weight',mt_dbl,1,weight))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:dielec',mt_dbl,1,dielec))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:ealpha',mt_dbl,1,ealpha))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:rbox',mt_dbl,1,rbox))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:disrlx',mt_dbl,1,disrlx))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:drsscl',mt_dbl,1,drsscl))
     + call errquit('md_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'md:uqmatm',mt_dbl,1,eqm))
     + call errquit('md_input: rtdb_put failed',0)
c
      close(unit=lfncmd)
c
      return
c
 1999 continue
      call errquit('md_input : unable to open cmd file',13)
c
      return
      end
