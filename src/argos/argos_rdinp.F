      subroutine argos_rdinp()
c
      implicit none
c
#include "argos.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "argos_constants.fh"
#if defined(NEED_LOC)
      external loc
      integer loc
#endif
c
      character*20 operat
      character*32 theory
      logical lstate,lstart,lrstrt,lcont
      integer igt,igr,nbytes,ibl,iun
      integer niperw,nbperw,istart
      character*5 gid
c
c     set process id and number of processes
c     --------------------------------------
c
      me=ga_nodeid()
      np=ga_nnodes()
c
c     set logical file numbers
c
      lfninp=10
      lfnout=11
      lfntop=12
      lfnrst=13
      lfntrj=14
      lfnprp=15
      lfngib=16
      lfnqrs=17
      lfndbg=18
      lfnmri=19
      lfnmro=20
      lfntim=21
      lfnsyn=22
c
      lfnmrr=23
      lfnarg=24
      lfnsum=25
      lfnlog=26
      lfnind=27
      lfnsin=28
      lfnrdi=29
      lfnfld=30
      lfnsfl=31
      lfntst=32
      lfnrdf=33
      lfnppd=34
      lfndip=35
      lfndef=36
      lfnhis=37
      lfnhbo=38
      lfnkrk=39
      lfnqsc=40
      lfnacf=41
      lfncnv=42
      lfnfet=43
      lfnuse=44
      lfnmsg=45
      lfnday=46
      lfncmd=47
      lfnpmf=48
      lfntri=49
      lfnhop=50
c
c     set defaults
c
      titinp(1)='NWChem:MD input'
      titinp(2)=' '
      titinp(3)=' '
c
      call swatch(datinp,timinp)
c
c     print flags
c
      nptopw=0
      nptops=0
      nptopt=0
      npstep=0
      npstat=0
      npener=0
      npforc=0
      nptmng=0
      npmemo=0
c
c     default task is single point energy set 1
c
      ntype=0
      mdtype=1
      iset=1
      mdordr=0
      nserie=0
c
c     polarization defaults
c
      ipolt=0
      mpolit=1
      ptol=1.0d-5
c
c     default cutoff radii
c
      iswtch=0
      rshort=0.9d0
      rlong=0.9d0
c
c     md defaults
c
      stime=zero
      mdstep=0
      tstep=0.001d0
      kequi=0
      mequi=0
      kdacq=0
      mdacq=100
      npgdrv=0
      npgdec=0
      iffdrv=0
      isaltb=0
c
c     mcti defaults
c
      maxlam=21
      elam=one
      ddacq=5.0d0
      edacq=5.0d0
      fdacq=0.75d0
      macfl=1
      ixcl=0
      iapprx=0
      weight=-2.5d-2
      facapp=two
c
      dgscl=one
      dgref=5.0d-4
      ddgscl=5.0d-4
c
      ssshft=0.075d0
c
c     multiple run defaults
c
      mrun=0
      mropt=0
      ldacq=0
      mdopt=0
      msplit=0
c
c     em defaults
c
      mintyp=0
      nem=0
      nemcrt=3
      nfqrs=25
c
      msdit=100
      dx0sd=0.01d0
      dxmsd=0.05d0
      dxsdmx=0.00001d0
c
      mcgit=0
      ncgcy=0
      dx0cg=0.01d0
      dxcgmx=0.00001d0
c
c     shake defaults
c
      mshitw=100
      mshits=100
      tlwsha=0.001d0
      tlssha=0.001d0
      ifss=0
c
c     pairlist defaults
c
      lwtype=0
c
c     frequency defaults
c
      nfoutp=-1
      nfstat=-1
      nfrest=0
      nffree=-1
c
      nfcoor=0
      nfscoo=0
      nfprop=0
      nfvelo=0
      nfsvel=0
      nfforc=0
      nfsfor=0
      nfpold=0
      nrwrec=0
      nfindu=0
      nfsind=0
c
      nfcntr=0
      nfslow=0
      nfshkw=0
      nfshks=0
c
      nfrdf=0
      nfdip=0
      nfkirk=0
      nkirk=0
      nfhbo=0
      drdf=0.0d0
c
      lenhis=0
      lendis=0
c
      nfsync=0
      iprop=0
c
      lserver=.false.
c
c     distributions
c
      ngl=0
      rrdf=zero
      ngr=0
      ngc=0
      ngrww=0
      ngrsw=0
      ngrss=0
c
      ndip=0
      rdip=zero
      rkirk=zero
c
      numhis=0
      numdis=0
      lnghis=0
      lngdis=0
c
      nfacfa=0
      nfauto=0
      nfconv=0
c
      nfdip=0
      nfkirk=0
c
      issscl=0
c
      nfnewf=0
      ibatch=0
c
c     constant p
c
      ipscal=0
      ipopt=0
      prsext=1.025d5
      prsrlx=0.5d0
      compr=4.53d-10
c
c     constant t
c
      itscal=0
      tmpext=298.15d0
      tmpext1=298.15d0
      tmpext2=298.15d0
      tmprlx=0.4d0
      tmsrlx=0.4d0
      tann1=0.0d00
      tann2=1.0d02
c
c     velocity reassignment
c
      ivreas=0
      ivopt=0
      tgauss=298.15d0
      iseed=12345
c
c     reaction field
c
      ireact=0
      dielec=one
c
c     external electric field
c
      ifield=0
      field=zero
c
c     solute centering
c
      nfcntr=0
c
c     print default
c
      nprint=3
c
c     recording defaults
c
      ibinar=0
c
c     load balancing
c
      factld=zero
c
c     other defaults (some are inactive)
c
      icbw=1
      icbs=1
      irr=0
      sil=0.001d0
      idipol=0
      iwarn=0
      irlow=0
      idebug=0
      iumbr=0
      igmass=0
      lowmcr=0
      noone=0
      npolit=0
c
      nwarn=0
      nform=0
      nopack=0
c
      urlow=tiny
c
      ignore=0
      mdo=6
      mwork=-1
      hbdis1=pthree
      hbdis2=pfour
      hbang1=two*atan(one)
      hbang2=two*atan(one)
      verinp=zero
c
      iformt=1
c
c     parallel execution defaults
c
c     many of the above defaults are superceded by defaults in the rtdb
c
c     general initializations
c
c     nbperw : number of bytes per word
c     niperw : number of integers per word
c
      nbperw=ma_sizeof(mt_dbl,1,mt_byte)
      niperw=nbperw/ma_sizeof(mt_int,1,mt_byte)
c
c     get input variables on process 0 only
c
      if(me.eq.0) then
c
c     sequential access to rtdb
c
      lstate=rtdb_parallel(.false.)
c
c     retrieve project name
c
      if(.not.rtdb_cget(irtdb,'md:project',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'prep:sysnam',1,filnam)) then
      if(.not.rtdb_cget(irtdb,'file_prefix',1,filnam))
     + call md_abort('argos_rdname: error rtdb_cget project',0)
      endif
      endif
      projct=filnam
c
      ibl=index(filnam,' ')-1
      if(ibl.lt.1) then
      filnam='nwmd'
      ibl=4
      endif
      iun=index(filnam,'_')-1
      if(iun.lt.0) iun=ibl
      if(iun.eq.0) then
      ibl=index(filnam,' ')-1
      iun=5
      endif
c
c     set file names
c
      root=filnam(1:ibl)
      filday=filnam(1:ibl)//'.day'
      filinp=filnam(1:ibl)//'.inp'
      filrst=filnam(1:ibl)//'.rst'
      filtop=filnam(1:iun)//'.top'
      filtrj=filnam(1:ibl)//'.trj'
      filprp=filnam(1:ibl)//'.prp'
      filgib=filnam(1:ibl)//'.gib'
      filppd=filnam(1:ibl)//'.ppd'
      filhis=filnam(1:ibl)//'.his'
      filhbo=filnam(1:ibl)//'.hbo'
      filmro=filnam(1:ibl)//'.mro'
      filmrr=filnam(1:ibl)//'.mrr'
      filmri=filnam(1:ibl)//'.mri'
      filout=filnam(1:ibl)//'.out'
      filmem=filnam(1:ibl)//'.mem'
      fildef=filnam(1:ibl)//'.def'
      filqsc=filnam(1:ibl)//'.qsc'
      filsur=filnam(1:ibl)//'.sur'
      filqrs=filnam(1:ibl)//'.qrs'
      filrdf=filnam(1:ibl)//'.rdf'
      filddf=filnam(1:ibl)//'.ddf'
      filkrk=filnam(1:ibl)//'.krk'
      fildbg=filnam(1:ibl)//'.dbg'
      filacf=filnam(1:ibl)//'.acf'
      filcnv=filnam(1:ibl)//'.cnv'
      filfet=filnam(1:ibl)//'.fet'
      filind=filnam(1:ibl)//'.ind'
      filsin=filnam(1:ibl)//'.sin'
      filtim=filnam(1:ibl)//'.tim'
      filsyn=filnam(1:ibl)//'.syn'
      filrdi=filnam(1:ibl)//'.rdi'
      filfld=filnam(1:ibl)//'.fld'
      filsfl=filnam(1:ibl)//'.sfl'
      filtst=filnam(1:ibl)//'.tst'
      filcmd=filnam(1:ibl)//'.cmd'
      filtri=filnam(1:ibl)//'.tri'
      filpmf=filnam(1:ibl)//'.pmf'
      filhop=filnam(1:ibl)//'.hop'
      rfile=filrst
c
c     input from rtdb
c
      if(.not.rtdb_get(irtdb,'md:port',mt_int,1,iport))
     + call md_abort('argos_rtdbin: rtdb_get failed',0)
      if(iport.gt.0) then
      if(.not.rtdb_cget(irtdb,'md:server',1,server))
     + call md_abort('argos_rtdbin: rtdb_get failed',0)
      endif
      if(.not.rtdb_get(irtdb,'md:istart',mt_int,1,istart))
     + call md_abort('argos_rtdbin: rtdb_get failed',1)
      if(.not.rtdb_get(irtdb,'md:nbx',mt_int,1,nbx))
     + call md_abort('argos_rtdbin: rtdb_get failed',1)
      if(.not.rtdb_get(irtdb,'md:nby',mt_int,1,nby))
     + call md_abort('argos_rtdbin: rtdb_get failed',2)
      if(.not.rtdb_get(irtdb,'md:nbz',mt_int,1,nbz))
     + call md_abort('argos_rtdbin: rtdb_get failed',3)
      if(.not.rtdb_get(irtdb,'md:npx',mt_int,1,npx))
     + call md_abort('argos_rtdbin: rtdb_get failed',4)
      if(.not.rtdb_get(irtdb,'md:npy',mt_int,1,npy))
     + call md_abort('argos_rtdbin: rtdb_get failed',5)
      if(.not.rtdb_get(irtdb,'md:npz',mt_int,1,npz))
     + call md_abort('argos_rtdbin: rtdb_get failed',6)
      if(.not.rtdb_get(irtdb,'md:npg',mt_int,1,npg))
     + call md_abort('argos_rtdbin: rtdb_get failed',6)
      if(.not.rtdb_get(irtdb,'md:ncp',mt_int,1,ncp))
     + call md_abort('argos_rtdbin: rtdb_get failed',6)
      if(.not.rtdb_get(irtdb,'md:mdalgo',mt_int,1,mdalgo))
     + call md_abort('argos_rtdbin: rtdb_get failed',7)
      if(.not.rtdb_get(irtdb,'md:iset',mt_int,1,iset))
     + call md_abort('argos_rtdbin: rtdb_get failed',8)
      if(.not.rtdb_get(irtdb,'md:isetp1',mt_int,1,isetp1))
     + call md_abort('argos_rtdbin: rtdb_get failed',9)
      if(.not.rtdb_get(irtdb,'md:isetp2',mt_int,1,isetp2))
     + call md_abort('argos_rtdbin: rtdb_get failed',10)
      if(.not.rtdb_get(irtdb,'md:iforw',mt_int,1,lamtyp))
     + call md_abort('argos_rtdbin: rtdb_get failed',11)
      if(.not.rtdb_get(irtdb,'md:msdit',mt_int,1,msdit))
     + call md_abort('argos_rtdbin: rtdb_get failed',12)
      if(.not.rtdb_get(irtdb,'md:mcgit',mt_int,1,mcgit))
     + call md_abort('argos_rtdbin: rtdb_get failed',13)
      if(.not.rtdb_get(irtdb,'md:icrit',mt_int,1,icrit))
     + call md_abort('argos_rtdbin: rtdb_get failed',13)
      if(.not.rtdb_get(irtdb,'md:ncgcy',mt_int,1,ncgcy))
     + call md_abort('argos_rtdbin: rtdb_get failed',14)
      if(.not.rtdb_get(irtdb,'md:mrun',mt_int,1,mrun))
     + call md_abort('argos_rtdbin: rtdb_get failed',15)
      if(.not.rtdb_get(irtdb,'md:maxlam',mt_int,1,maxlam))
     + call md_abort('argos_rtdbin: rtdb_get failed',16)
      if(.not.rtdb_get(irtdb,'md:npgdec',mt_int,1,npgdec))
     + call md_abort('argos_rtdbin: rtdb_get failed',17)
      if(.not.rtdb_get(irtdb,'md:issscl',mt_int,1,issscl))
     + call md_abort('argos_rtdbin: rtdb_get failed',18)
      if(.not.rtdb_get(irtdb,'md:mequi',mt_int,1,mequi))
     + call md_abort('argos_rtdbin: rtdb_get failed',19)
      if(.not.rtdb_get(irtdb,'md:mdacq',mt_int,1,mdacq))
     + call md_abort('argos_rtdbin: rtdb_get failed',20)
      if(.not.rtdb_get(irtdb,'md:ldacq',mt_int,1,ldacq))
     + call md_abort('argos_rtdbin: rtdb_get failed',21)
      if(.not.rtdb_get(irtdb,'md:iapprx',mt_int,1,iapprx))
     + call md_abort('argos_rtdbin: rtdb_get failed',22)
      if(.not.rtdb_get(irtdb,'md:nacfl',mt_int,1,nacfl))
     + call md_abort('argos_rtdbin: rtdb_get failed',23)
      if(.not.rtdb_get(irtdb,'md:ipscal',mt_int,1,ipscal))
     + call md_abort('argos_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:ipopt',mt_int,1,ipopt))
     + call md_abort('argos_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:ivopt',mt_int,1,ivopt))
     + call md_abort('argos_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:itscal',mt_int,1,itscal))
     + call md_abort('argos_rtdbin: rtdb_get failed',25)
      if(.not.rtdb_get(irtdb,'md:nfgaus',mt_int,1,nfgaus))
     + call md_abort('argos_rtdbin: rtdb_get failed',26)
      if(.not.rtdb_get(irtdb,'md:ipolt',mt_int,1,ipolt))
     + call md_abort('argos_rtdbin: rtdb_get failed',27)
      if(.not.rtdb_get(irtdb,'md:mpolit',mt_int,1,mpolit))
     + call md_abort('argos_rtdbin: rtdb_get failed',28)
      if(.not.rtdb_get(irtdb,'md:mshitw',mt_int,1,mshitw))
     + call md_abort('argos_rtdbin: rtdb_get failed',29)
      if(.not.rtdb_get(irtdb,'md:mshits',mt_int,1,mshits))
     + call md_abort('argos_rtdbin: rtdb_get failed',30)
      if(.not.rtdb_get(irtdb,'md:itrack',mt_int,1,itrack))
     + call md_abort('argos_rtdbin: rtdb_get failed',31)
      if(.not.rtdb_get(irtdb,'md:npstep',mt_int,1,npstep))
     + call md_abort('argos_rtdbin: rtdb_get failed',32)
      if(.not.rtdb_get(irtdb,'md:npstat',mt_int,1,npstat))
     + call md_abort('argos_rtdbin: rtdb_get failed',32)
      if(.not.rtdb_get(irtdb,'md:npener',mt_int,1,npener))
     + call md_abort('argos_rtdbin: rtdb_get failed',33)
      if(.not.rtdb_get(irtdb,'md:npforc',mt_int,1,npforc))
     + call md_abort('argos_rtdbin: rtdb_get failed',33)
      if(.not.rtdb_get(irtdb,'md:npdist',mt_int,1,npdist))
     + call md_abort('argos_rtdbin: rtdb_get failed',34)
      if(.not.rtdb_get(irtdb,'md:nptmng',mt_int,1,nptmng))
     + call md_abort('argos_rtdbin: rtdb_get failed',35)
      if(.not.rtdb_get(irtdb,'md:npatom',mt_int,1,npatom))
     + call md_abort('argos_rtdbin: rtdb_get failed',36)
      if(.not.rtdb_get(irtdb,'md:nptopw',mt_int,1,nptopw))
     + call md_abort('argos_rtdbin: rtdb_get failed',37)
      if(.not.rtdb_get(irtdb,'md:nptops',mt_int,1,nptops))
     + call md_abort('argos_rtdbin: rtdb_get failed',38)
      if(.not.rtdb_get(irtdb,'md:npxpct',mt_int,1,npxpct))
     + call md_abort('argos_rtdbin: rtdb_get failed',39)
      if(.not.rtdb_get(irtdb,'md:nfpair',mt_int,1,nfpair))
     + call md_abort('argos_rtdbin: rtdb_get failed',40)
      if(.not.rtdb_get(irtdb,'md:nfesp',mt_int,1,nfesp))
     + call md_abort('argos_rtdbin: rtdb_get failed',40)
      if(.not.rtdb_get(irtdb,'md:nfrdf',mt_int,1,nfrdf))
     + call md_abort('argos_rtdbin: rtdb_get failed',41)
      if(.not.rtdb_get(irtdb,'md:nflong',mt_int,1,nflong))
     + call md_abort('argos_rtdbin: rtdb_get failed',42)
      if(.not.rtdb_get(irtdb,'md:nfcntr',mt_int,1,nfcntr))
     + call md_abort('argos_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:icentr',mt_int,1,icentr))
     + call md_abort('argos_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:nfanal',mt_int,1,nfanal))
     + call md_abort('argos_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:nfslow',mt_int,1,nfslow))
     + call md_abort('argos_rtdbin: rtdb_get failed',44)
      if(.not.rtdb_get(irtdb,'md:nfoutp',mt_int,1,nfoutp))
     + call md_abort('argos_rtdbin: rtdb_get failed',45)
      if(.not.rtdb_get(irtdb,'md:nfstat',mt_int,1,nfstat))
     + call md_abort('argos_rtdbin: rtdb_get failed',46)
      if(.not.rtdb_get(irtdb,'md:nfrest',mt_int,1,nfrest))
     + call md_abort('argos_rtdbin: rtdb_get failed',47)
      if(.not.rtdb_get(irtdb,'md:keepr',mt_int,1,keepr))
     + call md_abort('argos_rtdbin: rtdb_get failed',48)
      if(.not.rtdb_get(irtdb,'md:nfcoor',mt_int,1,nfcoor))
     + call md_abort('argos_rtdbin: rtdb_get failed',50)
      if(.not.rtdb_get(irtdb,'md:nfscoo',mt_int,1,nfscoo))
     + call md_abort('argos_rtdbin: rtdb_get failed',51)
      if(.not.rtdb_get(irtdb,'md:nfvelo',mt_int,1,nfvelo))
     + call md_abort('argos_rtdbin: rtdb_get failed',52)
      if(.not.rtdb_get(irtdb,'md:nfsvel',mt_int,1,nfsvel))
     + call md_abort('argos_rtdbin: rtdb_get failed',53)
      if(.not.rtdb_get(irtdb,'md:nfforc',mt_int,1,nfforc))
     + call md_abort('argos_rtdbin: rtdb_get failed',50)
      if(.not.rtdb_get(irtdb,'md:nfsfor',mt_int,1,nfsfor))
     + call md_abort('argos_rtdbin: rtdb_get failed',51)
      if(.not.rtdb_get(irtdb,'md:nfprop',mt_int,1,nfprop))
     + call md_abort('argos_rtdbin: rtdb_get failed',54)
      if(.not.rtdb_get(irtdb,'md:iprop',mt_int,1,iprop))
     + call md_abort('argos_rtdbin: rtdb_get failed',54)
      if(.not.rtdb_get(irtdb,'md:nffree',mt_int,1,nffree))
     + call md_abort('argos_rtdbin: rtdb_get failed',55)
      if(.not.rtdb_get(irtdb,'md:nfsync',mt_int,1,nfsync))
     + call md_abort('argos_rtdbin: rtdb_get failed',56)
      if(.not.rtdb_get(irtdb,'md:nfauto',mt_int,1,nfauto))
     + call md_abort('argos_rtdbin: rtdb_get failed',57)
      if(.not.rtdb_get(irtdb,'md:nfconv',mt_int,1,nfconv))
     + call md_abort('argos_rtdbin: rtdb_get failed',58)
      if(.not.rtdb_get(irtdb,'md:nffet',mt_int,1,nffet))
     + call md_abort('argos_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:impfr',mt_int,1,impfr))
     + call md_abort('argos_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:impto',mt_int,1,impto))
     + call md_abort('argos_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:nftri',mt_int,1,nftri))
     + call md_abort('argos_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:iformt',mt_int,1,iformt))
     + call md_abort('argos_rtdbin: rtdb_get failed',60)
      if(.not.rtdb_get(irtdb,'md:madbox',mt_int,1,madbox))
     + call md_abort('argos_rtdbin: rtdb_get failed',61)
      if(.not.rtdb_get(irtdb,'md:loadb',mt_int,1,loadb))
     + call md_abort('argos_rtdbin: rtdb_get failed',62)
      if(.not.rtdb_get(irtdb,'md:ireset',mt_int,1,ireset))
     + call md_abort('argos_rtdbin: rtdb_get failed',62)
      if(.not.rtdb_get(irtdb,'md:mropt',mt_int,1,mropt))
     + call md_abort('argos_rtdbin: rtdb_get failed',63)
      if(.not.rtdb_get(irtdb,'md:idebug',mt_int,1,idebug))
     + call md_abort('argos_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:icntrl',mt_int,1,icntrl))
     + call md_abort('argos_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:ifidi',mt_int,1,ifidi))
     + call md_abort('argos_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:ipbtyp',mt_int,1,ipbtyp))
     + call md_abort('argos_rtdbin: rtdb_get ipbtyp failed',64)
      if(.not.rtdb_get(irtdb,'md:ngl',mt_int,1,ngl))
     + call md_abort('argos_rtdbin: rtdb_get failed',65)
      if(.not.rtdb_get(irtdb,'md:ifield',mt_int,1,ifield))
     + call md_abort('argos_rtdbin: rtdb_get failed',66)
      if(.not.rtdb_get(irtdb,'md:ipme',mt_int,1,ipme))
     + call md_abort('argos_rtdbin: rtdb_get failed',67)
      if(.not.rtdb_get(irtdb,'md:ngx',mt_int,1,ngx))
     + call md_abort('argos_rtdbin: rtdb_get failed',68)
      if(.not.rtdb_get(irtdb,'md:ngy',mt_int,1,ngy))
     + call md_abort('argos_rtdbin: rtdb_get failed',69)
      if(.not.rtdb_get(irtdb,'md:ngz',mt_int,1,ngz))
     + call md_abort('argos_rtdbin: rtdb_get failed',70)
      if(.not.rtdb_get(irtdb,'md:numfix',mt_int,1,numfix))
     + call md_abort('argos_rtdbin: rtdb_get failed',71)
c      if(.not.rtdb_get(irtdb,'md:iunfix',mt_int,1,iunfix))
c     + call md_abort('argos_rtdbin: rtdb_get failed',72)
c      if(.not.rtdb_get(irtdb,'md:lsffix',mt_int,msf,lsffix))
c     + call md_abort('argos_rtdbin: rtdb_get failed',72)
      if(.not.rtdb_get(irtdb,'md:noshak',mt_int,1,noshak))
     + call md_abort('argos_rtdbin: rtdb_get failed',73)
      if(.not.rtdb_get(irtdb,'md:nfefld',mt_int,1,nfefld))
     + call md_abort('argos_rtdbin: rtdb_get failed',74)
      if(.not.rtdb_get(irtdb,'md:nfsfld',mt_int,1,nfsfld))
     + call md_abort('argos_rtdbin: rtdb_get failed',75)
      if(.not.rtdb_get(irtdb,'md:nscb',mt_int,1,nscb))
     + call md_abort('argos_rtdbin: rtdb_get failed',76)
      if(nscb.gt.10) call md_abort('Increase dimension of idscb',0)
      if(.not.rtdb_get(irtdb,'md:idscb',mt_int,nscb,idscb))
     + call md_abort('argos_rtdbin: rtdb_get failed',77)
      if(.not.rtdb_get(irtdb,'md:ireact',mt_int,1,ireact))
     + call md_abort('argos_rtdbin: rtdb_get failed',78)
      if(.not.rtdb_get(irtdb,'md:memlim',mt_int,1,memlim))
     + call md_abort('argos_rtdbin: rtdb_get failed',79)
      if(.not.rtdb_get(irtdb,'md:morder',mt_int,1,morder))
     + call md_abort('argos_rtdbin: rtdb_get failed',80)
      if(.not.rtdb_get(irtdb,'md:isolvo',mt_int,1,isolvo))
     + call md_abort('argos_rtdbin: rtdb_get failed',81)
      if(.not.rtdb_get(irtdb,'md:lfout6',mt_int,1,lfout6))
     + call md_abort('argos_rtdbin: rtdb_get failed',82)
      if(.not.rtdb_get(irtdb,'md:iprpmf',mt_int,1,iprpmf))
     + call md_abort('argos_rtdbin: rtdb_get failed',82)
      if(.not.rtdb_get(irtdb,'md:imfft',mt_int,1,imfft))
     + call md_abort('argos_rtdbin: rtdb_get failed',83)
      if(.not.rtdb_get(irtdb,'md:mwmreq',mt_int,1,mwmreq))
     + call md_abort('argos_rtdbin: rtdb_get failed',84)
      if(.not.rtdb_get(irtdb,'md:msareq',mt_int,1,msareq))
     + call md_abort('argos_rtdbin: rtdb_get failed',85)
      if(.not.rtdb_get(irtdb,'md:mbbreq',mt_int,1,mbbreq))
     + call md_abort('argos_rtdbin: rtdb_get failed',85)
      if(.not.rtdb_get(irtdb,'md:itest',mt_int,1,itest))
     + call md_abort('argos_rtdbin: rtdb_get failed',86)
      if(.not.rtdb_get(irtdb,'md:nodpme',mt_int,1,nodpme))
     + call md_abort('argos_rtdbin: rtdb_get failed',87)
      if(.not.rtdb_get(irtdb,'md:lbpair',mt_int,1,lbpair))
     + call md_abort('argos_rtdbin: rtdb_get failed',88)
      if(.not.rtdb_get(irtdb,'md:ndistr',mt_int,1,ndistr))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:npmf',mt_int,1,npmf))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:facpmf',mt_dbl,1,facpmf))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:ndaver',mt_int,1,ndaver))
     + call md_abort('argos_rtdbin: rtdb_get failed',90)
      if(.not.rtdb_get(irtdb,'md:idevel',mt_int,1,idevel))
     + call md_abort('argos_rtdbin: rtdb_get failed',91)
c      if(.not.rtdb_get(irtdb,'md:itime',mt_int,mtimes,itime))
c     + call md_abort('argos_rtdbin: rtdb_get failed',92)
      if(.not.rtdb_get(irtdb,'md:nftime',mt_int,1,nftime))
     + call md_abort('argos_rtdbin: rtdb_get failed',92)
      if(.not.rtdb_get(irtdb,'md:nfdrss',mt_int,1,nfdrss))
     + call md_abort('argos_rtdbin: rtdb_get failed',93)
      if(.not.rtdb_get(irtdb,'md:nfload',mt_int,1,nfload))
     + call md_abort('argos_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:itload',mt_int,1,itload))
     + call md_abort('argos_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:ioload',mt_int,1,ioload))
     + call md_abort('argos_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:isload',mt_int,1,isload))
     + call md_abort('argos_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:ihess',mt_int,1,ihess))
     + call md_abort('argos_rtdbin: rtdb_get failed',95)
      if(.not.rtdb_get(irtdb,'md:latom',mt_int,1,latom))
     + call md_abort('argos_rtdbin: rtdb_get failed',96)
      if(.not.rtdb_get(irtdb,'md:icomm',mt_int,1,icomm))
     + call md_abort('argos_rtdbin: rtdb_get failed',97)
      if(.not.rtdb_get(irtdb,'md:nfhop',mt_int,1,nfhop))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:rhop',mt_dbl,1,rhop))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:thop',mt_dbl,1,thop))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:iguide',mt_int,1,iguide))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:icmopt',mt_int,1,icmopt))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:imembr',mt_int,1,imembr))
     + call md_abort('argos_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:nopack',mt_int,1,nopack))
     + call md_abort('argos_rtdbin: rtdb_get failed',99)
      if(.not.rtdb_get(irtdb,'md:iprof',mt_int,1,iprof))
     + call md_abort('argos_rtdbin: rtdb_get failed',99)
      if(.not.rtdb_get(irtdb,'md:ncoll',mt_int,1,ncoll))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:ilambd',mt_int,1,ilambd))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:mlambd',mt_int,1,mlambd))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:includ',mt_int,1,includ))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:iradgy',mt_int,1,iradgy))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:idifco',mt_int,1,idifco))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nfnewf',mt_int,1,nfnewf))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nbget',mt_int,1,nbget))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nprec',mt_int,1,nprec))
     + call md_abort('argos_rtdbin: rtdb_get failed',89)
c
      if(.not.rtdb_get(irtdb,'md:fguide',mt_dbl,1,fguide))
     + call md_abort('argos_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:tguide',mt_dbl,1,tguide))
     + call md_abort('argos_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:dx0sd',mt_dbl,1,dx0sd))
     + call md_abort('argos_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:dxsdmx',mt_dbl,1,dxsdmx))
     + call md_abort('argos_rtdbin: rtdb_get failed',101)
      if(.not.rtdb_get(irtdb,'md:dxmsd',mt_dbl,1,dxmsd))
     + call md_abort('argos_rtdbin: rtdb_get failed',102)
      if(.not.rtdb_get(irtdb,'md:dx0cg',mt_dbl,1,dx0cg))
     + call md_abort('argos_rtdbin: rtdb_get failed',103)
      if(.not.rtdb_get(irtdb,'md:dxcgmx',mt_dbl,1,dxcgmx))
     + call md_abort('argos_rtdbin: rtdb_get failed',104)
      if(.not.rtdb_get(irtdb,'md:dxmcg',mt_dbl,1,dxmcg))
     + call md_abort('argos_rtdbin: rtdb_get failed',105)
      if(.not.rtdb_get(irtdb,'md:edacq',mt_dbl,1,edacq))
     + call md_abort('argos_rtdbin: rtdb_get failed',106)
      if(.not.rtdb_get(irtdb,'md:ddacq',mt_dbl,1,ddacq))
     + call md_abort('argos_rtdbin: rtdb_get failed',107)
      if(.not.rtdb_get(irtdb,'md:fdacq',mt_dbl,1,fdacq))
     + call md_abort('argos_rtdbin: rtdb_get failed',108)
      if(.not.rtdb_get(irtdb,'md:delta',mt_dbl,1,delta))
     + call md_abort('argos_rtdbin: rtdb_get failed',109)
      if(.not.rtdb_get(irtdb,'md:stime',mt_dbl,1,stime))
     + call md_abort('argos_rtdbin: rtdb_get failed',110)
      if(.not.rtdb_get(irtdb,'md:tstep',mt_dbl,1,tstep))
     + call md_abort('argos_rtdbin: rtdb_get failed',111)
      if(.not.rtdb_get(irtdb,'md:prsext',mt_dbl,1,prsext))
     + call md_abort('argos_rtdbin: rtdb_get failed',112)
      if(.not.rtdb_get(irtdb,'md:prsrlx',mt_dbl,1,prsrlx))
     + call md_abort('argos_rtdbin: rtdb_get failed',113)
      if(.not.rtdb_get(irtdb,'md:compr',mt_dbl,1,compr))
     + call md_abort('argos_rtdbin: rtdb_get failed',114)
      if(.not.rtdb_get(irtdb,'md:tmpext',mt_dbl,1,tmpext1))
     + call md_abort('argos_rtdbin: rtdb_get failed',115)
      tmpext=tmpext1
      if(.not.rtdb_get(irtdb,'md:tmpext2',mt_dbl,1,tmpext2))
     + call md_abort('argos_rtdbin: rtdb_get failed',115)
      if(.not.rtdb_get(irtdb,'md:tmprlx',mt_dbl,1,tmprlx))
     + call md_abort('argos_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tmsrlx',mt_dbl,1,tmsrlx))
     + call md_abort('argos_rtdbin: rtdb_get failed',117)
      if(.not.rtdb_get(irtdb,'md:tann1',mt_dbl,1,tann1))
     + call md_abort('argos_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tann2',mt_dbl,1,tann2))
     + call md_abort('argos_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tgauss',mt_dbl,1,tgauss))
     + call md_abort('argos_rtdbin: rtdb_get failed',118)
      if(.not.rtdb_get(irtdb,'md:frgaus',mt_dbl,1,frgaus))
     + call md_abort('argos_rtdbin: rtdb_get failed',119)
      if(.not.rtdb_get(irtdb,'md:rlong',mt_dbl,1,rlong))
     + call md_abort('argos_rtdbin: rtdb_get failed',120)
      if(.not.rtdb_get(irtdb,'md:rshort',mt_dbl,1,rshort))
     + call md_abort('argos_rtdbin: rtdb_get failed',121)
      if(.not.rtdb_get(irtdb,'md:ptol',mt_dbl,1,ptol))
     + call md_abort('argos_rtdbin: rtdb_get failed',123)
      if(.not.rtdb_get(irtdb,'md:tlwsha',mt_dbl,1,tlwsha))
     + call md_abort('argos_rtdbin: rtdb_get failed',124)
      if(.not.rtdb_get(irtdb,'md:tlssha',mt_dbl,1,tlssha))
     + call md_abort('argos_rtdbin: rtdb_get failed',125)
      if(.not.rtdb_get(irtdb,'md:factld',mt_dbl,1,factld))
     + call md_abort('argos_rtdbin: rtdb_get failed',126)
      if(.not.rtdb_get(irtdb,'md:rrdf',mt_dbl,1,rrdf))
     + call md_abort('argos_rtdbin: rtdb_get failed',127)
      if(.not.rtdb_get(irtdb,'md:xfield',mt_dbl,1,xfield))
     + call md_abort('argos_rtdbin: rtdb_get failed',128)
      if(.not.rtdb_get(irtdb,'md:xffreq',mt_dbl,1,xffreq))
     + call md_abort('argos_rtdbin: rtdb_get failed',129)
      if(.not.rtdb_get(irtdb,'md:xfvect',mt_dbl,3,xfvect))
     + call md_abort('argos_rtdbin: rtdb_get failed',130)
      if(.not.rtdb_get(irtdb,'md:weight',mt_dbl,1,weight))
     + call md_abort('argos_rtdbin: rtdb_get failed',131)
      if(.not.rtdb_get(irtdb,'md:dielec',mt_dbl,1,dielec))
     + call md_abort('argos_rtdbin: rtdb_get failed',132)
      if(.not.rtdb_get(irtdb,'md:ealpha',mt_dbl,1,ealpha))
     + call md_abort('argos_rtdbin: rtdb_get failed',133)
      if(.not.rtdb_get(irtdb,'md:rbox',mt_dbl,1,rbox))
     + call md_abort('argos_rtdbin: rtdb_get failed',134)
      if(.not.rtdb_get(irtdb,'md:rsgm',mt_dbl,1,rsgm))
     + call md_abort('argos_rtdbin: rtdb_get failed',134)
      if(.not.rtdb_get(irtdb,'md:disrlx',mt_dbl,1,disrlx))
     + call md_abort('argos_rtdbin: rtdb_get failed',135)
      if(.not.rtdb_get(irtdb,'md:drsscl',mt_dbl,1,drsscl))
     + call md_abort('argos_rtdbin: rtdb_get failed',136)
      if(.not.rtdb_get(irtdb,'md:pmetol',mt_dbl,1,pmetol))
     + pmetol=1.0d-05
      if(.not.rtdb_get(irtdb,'md:fcoll',mt_dbl,1,fcoll))
     + fcoll=1.0d+01
      if(.not.rtdb_get(irtdb,'md:scaleq',mt_dbl,1,scaleq))
     + scaleq=-1.0d+00
c
      if(.not.rtdb_cget(irtdb,'task:theory',1,theory))
     + call md_abort('argos_rtdbin: rtdg_cget failed',137)
c
      if(.not.rtdb_cget(irtdb,'task:operation',1,operat))
c     + call md_abort('argos_rtdbin: rtdg_cget failed',138)
     + operat='energy'
c
      if(operat.eq.'energy'.or.operat.eq.'ENERGY') then
      ntype=0
      endif
c
      if(operat.eq.'optimize'.or.operat.eq.'OPTIMIZE') then
      ntype=1
      if(nfrest.gt.0) nfqrs=nfrest
      nfvelo=0
      nfsvel=0
      endif
c
      if(operat.eq.'dynamics'.or.operat.eq.'DYNAMICS') then
      ntype=2
      mdtype=iset
      if(isetp1.eq.2) mdtype=4
      if(isetp1.eq.3) mdtype=5
      if(isetp1.eq.2.and.isetp2.eq.3) mdtype=6
      if(mdtype.gt.3) iset=1
      endif
c
      if(operat.eq.'thermodynamics'.or.operat.eq.'THERMODYNAMICS') then
      ntype=3
      endif
c
      ssshft=delta
c
      if(nodpme.gt.ngz) nodpme=ngz
      if(nodpme.gt.np) nodpme=np
      if(nodpme.eq.0) nodpme=np
      if(nodpme.gt.ngz) nodpme=ngz
c
c     determine if this is a start, restart or continue
c
      call util_get_rtdb_state(irtdb,lstart,lrstrt,lcont)
c
      if(istart.gt.0) then
      lstart=istart.eq.0
      lrstrt=istart.eq.1
      lcont=istart.eq.2
      endif
c
      if(lstart) then
      nserie=0
      else
      if(lrstrt) then
      nserie=1
      else
      if(.not.lcont)
     + call md_abort('argos_rtdbin: failed to determine rtdb state',0)
      nserie=2
      endif
      endif
c
c     change lfnout if output to unit 6 is requested
c
      if(lfout6.ne.0) lfnout=6
c
c     pairlist type
c
      lstype=1
      if(latom.eq.1) lstype=0
c
c     Count radial distribution functions contributions
c     -------------------------------------------------
c
      ngc=0
      ngr=0
      ngrww=0
      ngrsw=0
      ngrss=0
      if(nfrdf.gt.0) then
      open(unit=lfnrdi,file=filrdi,form='formatted',status='old')
      rewind(unit=lfnrdi)
    1 continue
      read(lfnrdi,*,end=9999) igt,igr
      ngc=ngc+1
      if(igr.gt.ngr) ngr=igr
      if(igt.eq.1) ngrww=ngrww+1
      if(igt.eq.2) ngrsw=ngrsw+1
      if(igt.eq.3) ngrss=ngrss+1
      goto 1
 9999 continue
      close(unit=lfnrdi)
c
      if(ngc.gt.0) then
      if(ngl.eq.0) ngl=1000
      if(rrdf.lt.small) then
      if(ngrww.gt.0.and.rrdf.lt.rshort) rrdf=rshort
      if(ngrsw.gt.0.and.rrdf.lt.rshort) rrdf=rshort
      if(ngrss.gt.0.and.rrdf.lt.rshort) rrdf=rshort
      endif
      else
      ngl=1
      nfrdf=0
      endif
c
      endif
c     routine argos_chckin performs a check of selected input
c
c     MCTI initializations
c
      if(mdtype.ge.7.and.ntype.eq.3.and.mrun.eq.0) mrun=maxlam
      if(mdtype.ge.7.and.ntype.eq.3.and.nffree.lt.0) nffree=1
      if(mdtype.ge.7.and.ntype.eq.3.and.nfoutp.lt.0) nfoutp=1
      if(mdtype.ge.7.and.ntype.eq.3.and.nfstat.lt.0) nfstat=1
      if(mdtype.ge.7.and.ntype.eq.3.and.ldacq.eq.0) ldacq=mdacq
      if(mdtype.ge.7.and.ntype.eq.3.and.macfl.le.mdacq) macfl=mdacq+1
      if(mdtype.ge.7.and.maxlam.le.1)
     + call md_abort('Number of MCTI ensembles too small:',maxlam)
c
c     MD initializations
c
      if(nfoutp.lt.0) nfoutp=mdacq/100
      if(nfstat.lt.0) nfstat=mdacq/10
c
c     square SHAKE tolerances
c
      tlwsha=tlwsha*tlwsha
      tlssha=tlssha*tlssha
c
c     handle zero frequencies
c
      iint=0
      if(nfprop.lt.0) then
      iint=1
      nfprop=-nfprop
      endif
      if(nserie.ne.0) irr=1
c
      if(mrun.lt.0) then
      mrun=-mrun
      noone=mrun-1
      msplit=3
      endif
c
      if(nfdip.eq.0) then
      nfdip=1
      ndip=0
      endif
c
      if(nfkirk.eq.0) then
      nfkirk=1
      nkirk=0
      endif
c
      ivreas=nfgaus
      if(nfgaus.lt.0) nfgaus=-nfgaus
c
c     only 3D-periodic boundary systems can be at constant pressure
c
c      if(npbtyp.ne.1) ipscal=0
c
      lstate=rtdb_parallel(.true.)
c
c     get variables from restart file if not initial run
c
      if(nserie.gt.0) call argos_rdrest(lfnrst,filrst)
c
      endif
c
c     broadcast to all nodes
c
      if(np.gt.1) then
      nbytes=loc(inp_ptr)-loc(ptol)
      call ga_brdcst(mrg_d36,ptol,nbytes,0)
      nbytes=7*ma_sizeof(mt_int,1,mt_byte)
      call ga_brdcst(mrg_d37,npx,nbytes,0)
      call util_char_ga_brdcst(mrg_d38,filnam,0)
      endif
c
      lesp=.false.
      lpert2=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.2.or.isetp2.eq.2))
      lpert3=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.3.or.isetp2.eq.3))
c
      lpme=ipme.gt.0
      ltwin=rlong.gt.rshort.or.lpme
c
      lpack=nopack.eq.0
c
      msa=msareq
      mwm=mwmreq
c
      factgf=one
      if(iguide.gt.0) factgf=tstep/tguide
      factgg=one-factgf
c
      costio=zero
      lpmfc=.true.
c
      return
      end
c $Id$
