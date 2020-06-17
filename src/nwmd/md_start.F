      subroutine md_start
c
c $Id$
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      character*5 gid
      integer iun,ibl
      logical lexist
c
      if(npg.gt.1) call md_partition
c
      if(me.eq.0) then
      if(npg.gt.1) then
      write(gid,3300) ipg
 3300 format(i5.5)
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
      filtop=filnam(1:iun)//'.top'
      filrst=filnam(1:ibl)//'.rst'
      filout=filnam(1:ibl)//gid//'.out'
      filtrj=filnam(1:ibl)//gid//'.trj'
      filprp=filnam(1:ibl)//gid//'.prp'
      filmro=filnam(1:ibl)//gid//'.mro'
      filgib=filnam(1:ibl)//gid//'.gib'
      if(npg.gt.1) then
      filrst=filnam(1:ibl)//gid//'.rst'
      inquire(file=filrst(1:index(filrst,' ')-1),exist=lexist)
      if(.not.lexist) filrst=filnam(1:ibl)//'.rst'
      endif
      endif
      endif
c
c     Open debug file if requested
c     ----------------------------
c
      if(idebug.gt.0) then
      lfndbg=18
      write(fildbg,'(a,i5.5,a)') 'nwchem_',me,'.dbg'
      open(unit=lfndbg,file=fildbg,form='formatted',status='unknown')
      endif
c
c     open recording file
c
      if(me.eq.0) then
      call md_fopen(.false.)
c
      if(ntype.eq.3) then
      open(unit=lfngib,file=filgib(1:index(filgib,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(ntype.eq.0.and.nftri.gt.0) then
      open(unit=lfntri,file=filtri(1:index(filtri,' ')-1),
     + form='formatted',status='unknown')
      endif
c
      if(itest.gt.0) then
      open(unit=lfntst,file=filtst(1:index(filtst,' ')-1),
     + form='formatted',status='unknown')
      endif
      endif
c
      if(lqmd) call qmd_start()
c
c     print input information
c     -----------------------
c
      call md_print()
c
c     start the spacial decomposition API
c     -----------------------------------
c
      if(.not.lqmd) then
      call sp_start(lfnout,lfntop,filtop,lfnrst,filrst,lfnsyn,filsyn,
     + nfsync,rshort,rlong,zero,rsgm,
     + npx,npy,npz,nbx,nby,nbz,
     + npbtyp,nbxtyp,box,vlat,lpbc,
     + nwm,mwm,nwa,mwa,nsf,msf,nsm,msm,nsa,msa,
     + loadb,lbpair,factld,ipolt.ne.0,.false.,temp,tempw,temps,lqmd,
     + iguide,lfndbg,idebug,projct,mbbreq,nserie,isload,ireset,icntrl,
     + nseq,i_lseq,ndums,nbget,nprec,madbox)
      endif
c
c     start the analysis API
c     ----------------------
c
c      if(nfanal.gt.0) call ana_init(nsa,msa,.false.)
c
c     start the classical forces API
c     ------------------------------
c
      call cf_start(irtdb,lqmd,lfnout,lfntop,filtop,ndistr,npmf,npmfi,
     + nwm,mwm,nwa,mwa,nsf,msf,nsm,msm,nsa,msa,
     + mdalgo,npbtyp,nbxtyp,rshort,rlong,rqmmm,box,
     + ipme,morder,ngx,ngy,ngz,nodpme,pmetol,
     + tstep,tlwsha,mshitw,mshits,noshak,
     + iqmmm,ipolt,itscal,ipscal,ipopt,tmpext,prsext,
     + tmprlx,tmsrlx,prsrlx,scaleq,facpmf,
     + 0,temp,tempw,temps,compr,ntype,iset,isetp1,isetp2,
     + issscl,delta,nfanal,lpbc,npgdec,xfield,xfvect,xffreq,
     + npener,icntrl,nbias,mropt,includ,ltwin,
     + nseq,i_lseq,nfhop,rhop,thop,ndums,ipbtyp,lfnhop,iradgy,
     + 0,nprec)
c     + nbget)
c
      if(mlambd.gt.0.and.ilambd.gt.0.and.ilambd.le.mlambd)
     + call cf_lambda(lamtyp,irun,maxlam,elam,lfnout,lfnpmf,
     + rlambd,dlambd,filnam)
c
c     print topology information
c     --------------------------
c
      call cf_print_top(lfnout,npatom,nptopw,nptops)
c
      msm=max(1,nsm)
      msf=max(1,nsf)
      mst=max(msm,nseq)
c
c     start the property API
c     ----------------------
c
      call prp_start(nserie,ntype,nftri,lfnrst,filrst,lfnout,lfnprp,
     + lfngib,nfoutp,nfstat,nfprop,iprop,
     + .true.,.true.,ltwin,ipolt.ne.0,ipme.ne.0,
     + npstep.ne.0,npener.ne.0,npstat,
     + nwm,msf,nsf,mpe,mdacq,mrun,iset,isetp1,isetp2,tstep,msm,nsm,
     + nsa,ddacq,edacq,iprof,npmf,npener,ndistr,lpbc,nbias,nodpme,npmfi,
     + iguide.ne.0,iqmmm.ne.0,lqmd,iradgy,idifco,nbget,ipg,npg)
c
c     allocate memory for coordinates, velocities, etc.
c     -------------------------------------------------
c
      if(.not.ma_push_get(mt_int,mwm*miw2,'iw',l_iw,i_iw))
     + call md_abort('Failed to allocate memory for iw',0)
      if(.not.ma_push_get(mt_int,msa*mis2,'is',l_is,i_is))
     + call md_abort('Failed to allocate memory for is',0)
      if(.not.ma_push_get(mt_int,mwm,'iwz',l_iwz,i_iwz))
     + call md_abort('Failed to allocate memory for iw',0)
      if(.not.ma_push_get(mt_int,msa,'isz',l_isz,i_isz))
     + call md_abort('Failed to allocate memory for is',0)
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'xw',l_xw,i_xw))
     + call md_abort('Failed to allocate memory for xw',0)
      if(.not.ma_push_get(mt_dbl,3*mwm,'xwm',l_xwm,i_xwm))
     + call md_abort('Failed to allocate memory for xwm',0)
      if(.not.ma_push_get(mt_dbl,mwm,'rtos',l_rtos,i_rtos))
     + call md_abort('Failed to allocate memory for rtos',0)
      if(.not.ma_push_get(mt_dbl,3*msa,'xs',l_xs,i_xs))
     + call md_abort('Failed to allocate memory for xs',0)
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'yw',l_yw,i_yw))
     + call md_abort('Failed to allocate memory for yw',0)
      if(.not.ma_push_get(mt_dbl,3*msa,'ys',l_ys,i_ys))
     + call md_abort('Failed to allocate memory for ys',0)
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'vw',l_vw,i_vw))
     + call md_abort('Failed to allocate memory for vw',0)
      if(.not.ma_push_get(mt_dbl,3*msa,'vs',l_vs,i_vs))
     + call md_abort('Failed to allocate memory for vs',0)
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'vwt',l_vwt,i_vwt))
     + call md_abort('Failed to allocate memory for vwt',0)
      if(.not.ma_push_get(mt_dbl,3*msa,'vst',l_vst,i_vst))
     + call md_abort('Failed to allocate memory for vst',0)
      if(.not.ma_push_get(mt_dbl,6*mwa*mwm,'fw',l_fw,i_fw))
     + call md_abort('Failed to allocate memory for fw',0)
      if(.not.ma_push_get(mt_dbl,6*msa,'fs',l_fs,i_fs))
     + call md_abort('Failed to allocate memory for fs',0)
      if(.not.ma_push_get(mt_dbl,3*mwm,'xwcr',l_xwcr,i_xwcr))
     + call md_abort('Failed to allocate memory for xwcr',0)
      if(.not.ma_push_get(mt_dbl,3*msm,'xsm',l_xsm,i_xsm))
     + call md_abort('Failed to allocate memory for xsm',0)
      if(.not.ma_push_get(mt_dbl,4*mst,'tsm',l_tsm,i_tsm))
     + call md_abort('Failed to allocate memory for tsm',0)
      if(.not.ma_push_get(mt_dbl,3*msm,'xsm',l_xsmp,i_xsmp))
     + call md_abort('Failed to allocate memory for xsmp',0)
      if(.not.ma_push_get(mt_dbl,8*msm,'gsm',l_gsm,i_gsm))
     + call md_abort('Failed to allocate memory for gsm',0)
      if(.not.ma_push_get(mt_dbl,3*msm,'xscr',l_xscr,i_xscr))
     + call md_abort('Failed to allocate memory for xscr',0)
      if(.not.ma_push_get(mt_dbl,msm,'dsr',l_dsr,i_dsr))
     + call md_abort('Failed to allocate memory for dsr',0)
      if(.not.ma_push_get(mt_dbl,18*msm,'zs',l_zs,i_zs))
     + call md_abort('Failed to allocate memory for zs',0)
      if(.not.ma_push_get(mt_dbl,2*mpe*msf,'esw',l_esw,i_esw))
     + call md_abort('Failed to allocate memory for esw',0)
      if(.not.ma_push_get(mt_dbl,2*mpe*msf*msf,'ess',l_ess,i_ess))
     + call md_abort('Failed to allocate memory for ess',0)
      if(.not.ma_push_get(mt_dbl,6*msf*msf,'fss',l_fss,i_fss))
     + call md_abort('Failed to allocate memory for fss',0)
      if(.not.ma_push_get(mt_dbl,msf,'esk',l_esk,i_esk))
     + call md_abort('Failed to allocate memory for esk',0)
      if(.not.ma_push_get(mt_dbl,mwa+msa,'wws',l_wws,i_wws))
     + call md_abort('Failed to allocate memory for wws',me)
      if(npener.eq.0) then
      if(.not.ma_push_get(mt_dbl,1,'esa',l_esa,i_esa))
     + call md_abort('Failed to allocate memory for esa',0)
      else
      if(.not.ma_push_get(mt_dbl,2*nsa,'esa',l_esa,i_esa))
     + call md_abort('Failed to allocate memory for esa',0)
      endif
      if(ipolt.gt.0) then
      if(.not.ma_push_get(mt_dbl,6*mwa*mwm,'pw',l_pw,i_pw))
     + call md_abort('Failed to allocate memory for pw',0)
      if(.not.ma_push_get(mt_dbl,6*msa,'ps',l_ps,i_ps))
     + call md_abort('Failed to allocate memory for ps',0)
      if(lpert2.or.lpert3) then
      if(.not.ma_push_get(mt_dbl,12*mwa*mwm,'pwp',l_pwp,i_pwp))
     + call md_abort('Failed to allocate memory for pwp',0)
      if(.not.ma_push_get(mt_dbl,12*msa,'psp',l_psp,i_psp))
     + call md_abort('Failed to allocate memory for psp',0)
      endif
      endif
      if(iguide.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'gw',l_gw,i_gw))
     + call md_abort('Failed to allocate memory for gw',0)
      if(.not.ma_push_get(mt_dbl,3*msa,'gs',l_gs,i_gs))
     + call md_abort('Failed to allocate memory for gs',0)
      endif
      if(icmopt.gt.0) then
      if(.not.ma_push_get(mt_dbl,5*msm,'fcm',l_fcm,i_fcm))
     + call md_abort('Failed to allocate memory for fcm',0)
      endif
c
      if(imembr.gt.0) then
      if(.not.ma_push_get(mt_int,2*msa,'mm',l_mm,i_mm))
     + call md_abort('Failed to allocate memory for mm',me)
      if(.not.ma_push_get(mt_dbl,7*msm,'fm',l_fm,i_fm))
     + call md_abort('Failed to allocate memory for fm',me)
      endif
c
c     retrieve current coordinates for this node
c     ------------------------------------------
c
      if(lqmd) then
      call qmd_setup(nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),dbl_mb(i_gs),nsaloc)
      else
      call sp_setup(me,int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_xwcr),
     + dbl_mb(i_vw),dbl_mb(i_gw),nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_xscr),dbl_mb(i_vs),
     + dbl_mb(i_gs),nsaloc,lpack)
      endif
      if(.not.lqmd)
     + call sp_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
c
c     initialize packing
c     ------------------
c
      if(lpack)
     + call sp_pack_init(int_mb(i_is),nsaloc,int_mb(i_iw),nwmloc)
c
c     spacial decomposition
c     ---------------------
c
      if(.not.lqmd) then
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
      endif
c
c     calculate mass factors
c     ----------------------
c
      call cf_weight(nwmloc,nsaloc,int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lshop-1)*msa),wbox)
c
c     calculate centers of mass
c     -------------------------
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     fix
c     ------------
c
      if(me.eq.0.and.numfix.gt.0) then
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='old')
      rewind(lfncmd)
      endif
      call cf_fix(lfnout,lfncmd,numfix,int_mb(i_iw+(lwgmn-1)*mwm),
     + int_mb(i_iw+(lwdyn-1)*mwm),nwmloc,
     + int_mb(i_is+(lsgan-1)*msa),int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lssgm-1)*msa),nsaloc,
     + dbl_mb(i_xwm),dbl_mb(i_xs))
      if(me.eq.0.and.numfix.gt.0) then
      close(unit=lfncmd)
      endif
      if(.not.lqmd)
     + call sp_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
c
c     print decomposition information
c     -------------------------------
c
      if(.not.lqmd) call sp_print()
c
      call prp_setup(wbox)
c
c     write file headers
c     ------------------
c
      if(me.eq.0.and.ntype.ne.3) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0)
     + call cf_trjhdr(lfntrj)
      endif
c
      if(npg.gt.1.and.me.eq.0) then
      root=filnam(1:ibl)//gid
      filrst=filnam(1:ibl)//gid//'.rst'
      rfile=filrst
      endif
c
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('ma problems at end of md_start',me)
      return
      end
      subroutine md_rdinp()
c
      implicit none
c
#include "md_common.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwmd_constants.fh"
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
      uqmd=zero
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
     + call md_abort('md_rdname: error rtdb_cget project',0)
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
     + call md_abort('md_rtdbin: rtdb_get failed',0)
      if(iport.gt.0) then
      if(.not.rtdb_cget(irtdb,'md:server',1,server))
     + call md_abort('md_rtdbin: rtdb_get failed',0)
      endif
      if(.not.rtdb_get(irtdb,'md:istart',mt_int,1,istart))
     + call md_abort('md_rtdbin: rtdb_get failed',1)
      if(.not.rtdb_get(irtdb,'md:nbx',mt_int,1,nbx))
     + call md_abort('md_rtdbin: rtdb_get failed',1)
      if(.not.rtdb_get(irtdb,'md:nby',mt_int,1,nby))
     + call md_abort('md_rtdbin: rtdb_get failed',2)
      if(.not.rtdb_get(irtdb,'md:nbz',mt_int,1,nbz))
     + call md_abort('md_rtdbin: rtdb_get failed',3)
      if(.not.rtdb_get(irtdb,'md:npx',mt_int,1,npx))
     + call md_abort('md_rtdbin: rtdb_get failed',4)
      if(.not.rtdb_get(irtdb,'md:npy',mt_int,1,npy))
     + call md_abort('md_rtdbin: rtdb_get failed',5)
      if(.not.rtdb_get(irtdb,'md:npz',mt_int,1,npz))
     + call md_abort('md_rtdbin: rtdb_get failed',6)
      if(.not.rtdb_get(irtdb,'md:npg',mt_int,1,npg))
     + call md_abort('md_rtdbin: rtdb_get failed',6)
      if(.not.rtdb_get(irtdb,'md:mdalgo',mt_int,1,mdalgo))
     + call md_abort('md_rtdbin: rtdb_get failed',7)
      if(.not.rtdb_get(irtdb,'md:iset',mt_int,1,iset))
     + call md_abort('md_rtdbin: rtdb_get failed',8)
      if(.not.rtdb_get(irtdb,'md:isetp1',mt_int,1,isetp1))
     + call md_abort('md_rtdbin: rtdb_get failed',9)
      if(.not.rtdb_get(irtdb,'md:isetp2',mt_int,1,isetp2))
     + call md_abort('md_rtdbin: rtdb_get failed',10)
      if(.not.rtdb_get(irtdb,'md:iforw',mt_int,1,lamtyp))
     + call md_abort('md_rtdbin: rtdb_get failed',11)
      if(.not.rtdb_get(irtdb,'md:msdit',mt_int,1,msdit))
     + call md_abort('md_rtdbin: rtdb_get failed',12)
      if(.not.rtdb_get(irtdb,'md:mcgit',mt_int,1,mcgit))
     + call md_abort('md_rtdbin: rtdb_get failed',13)
      if(.not.rtdb_get(irtdb,'md:icrit',mt_int,1,icrit))
     + call md_abort('md_rtdbin: rtdb_get failed',13)
      if(.not.rtdb_get(irtdb,'md:ncgcy',mt_int,1,ncgcy))
     + call md_abort('md_rtdbin: rtdb_get failed',14)
      if(.not.rtdb_get(irtdb,'md:mrun',mt_int,1,mrun))
     + call md_abort('md_rtdbin: rtdb_get failed',15)
      if(.not.rtdb_get(irtdb,'md:maxlam',mt_int,1,maxlam))
     + call md_abort('md_rtdbin: rtdb_get failed',16)
      if(.not.rtdb_get(irtdb,'md:npgdec',mt_int,1,npgdec))
     + call md_abort('md_rtdbin: rtdb_get failed',17)
      if(.not.rtdb_get(irtdb,'md:issscl',mt_int,1,issscl))
     + call md_abort('md_rtdbin: rtdb_get failed',18)
      if(.not.rtdb_get(irtdb,'md:mequi',mt_int,1,mequi))
     + call md_abort('md_rtdbin: rtdb_get failed',19)
      if(.not.rtdb_get(irtdb,'md:mdacq',mt_int,1,mdacq))
     + call md_abort('md_rtdbin: rtdb_get failed',20)
      if(.not.rtdb_get(irtdb,'md:ldacq',mt_int,1,ldacq))
     + call md_abort('md_rtdbin: rtdb_get failed',21)
      if(.not.rtdb_get(irtdb,'md:iapprx',mt_int,1,iapprx))
     + call md_abort('md_rtdbin: rtdb_get failed',22)
      if(.not.rtdb_get(irtdb,'md:nacfl',mt_int,1,nacfl))
     + call md_abort('md_rtdbin: rtdb_get failed',23)
      if(.not.rtdb_get(irtdb,'md:ipscal',mt_int,1,ipscal))
     + call md_abort('md_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:ipopt',mt_int,1,ipopt))
     + call md_abort('md_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:ivopt',mt_int,1,ivopt))
     + call md_abort('md_rtdbin: rtdb_get failed',24)
      if(.not.rtdb_get(irtdb,'md:itscal',mt_int,1,itscal))
     + call md_abort('md_rtdbin: rtdb_get failed',25)
      if(.not.rtdb_get(irtdb,'md:nfgaus',mt_int,1,nfgaus))
     + call md_abort('md_rtdbin: rtdb_get failed',26)
      if(.not.rtdb_get(irtdb,'md:ipolt',mt_int,1,ipolt))
     + call md_abort('md_rtdbin: rtdb_get failed',27)
      if(.not.rtdb_get(irtdb,'md:mpolit',mt_int,1,mpolit))
     + call md_abort('md_rtdbin: rtdb_get failed',28)
      if(.not.rtdb_get(irtdb,'md:mshitw',mt_int,1,mshitw))
     + call md_abort('md_rtdbin: rtdb_get failed',29)
      if(.not.rtdb_get(irtdb,'md:mshits',mt_int,1,mshits))
     + call md_abort('md_rtdbin: rtdb_get failed',30)
      if(.not.rtdb_get(irtdb,'md:itrack',mt_int,1,itrack))
     + call md_abort('md_rtdbin: rtdb_get failed',31)
      if(.not.rtdb_get(irtdb,'md:npstep',mt_int,1,npstep))
     + call md_abort('md_rtdbin: rtdb_get failed',32)
      if(.not.rtdb_get(irtdb,'md:npstat',mt_int,1,npstat))
     + call md_abort('md_rtdbin: rtdb_get failed',32)
      if(.not.rtdb_get(irtdb,'md:npener',mt_int,1,npener))
     + call md_abort('md_rtdbin: rtdb_get failed',33)
      if(.not.rtdb_get(irtdb,'md:npforc',mt_int,1,npforc))
     + call md_abort('md_rtdbin: rtdb_get failed',33)
      if(.not.rtdb_get(irtdb,'md:npdist',mt_int,1,npdist))
     + call md_abort('md_rtdbin: rtdb_get failed',34)
      if(.not.rtdb_get(irtdb,'md:nptmng',mt_int,1,nptmng))
     + call md_abort('md_rtdbin: rtdb_get failed',35)
      if(.not.rtdb_get(irtdb,'md:npatom',mt_int,1,npatom))
     + call md_abort('md_rtdbin: rtdb_get failed',36)
      if(.not.rtdb_get(irtdb,'md:nptopw',mt_int,1,nptopw))
     + call md_abort('md_rtdbin: rtdb_get failed',37)
      if(.not.rtdb_get(irtdb,'md:nptops',mt_int,1,nptops))
     + call md_abort('md_rtdbin: rtdb_get failed',38)
      if(.not.rtdb_get(irtdb,'md:npxpct',mt_int,1,npxpct))
     + call md_abort('md_rtdbin: rtdb_get failed',39)
      if(.not.rtdb_get(irtdb,'md:nfpair',mt_int,1,nfpair))
     + call md_abort('md_rtdbin: rtdb_get failed',40)
      if(.not.rtdb_get(irtdb,'md:nfesp',mt_int,1,nfesp))
     + call md_abort('md_rtdbin: rtdb_get failed',40)
      if(.not.rtdb_get(irtdb,'md:nfrdf',mt_int,1,nfrdf))
     + call md_abort('md_rtdbin: rtdb_get failed',41)
      if(.not.rtdb_get(irtdb,'md:nflong',mt_int,1,nflong))
     + call md_abort('md_rtdbin: rtdb_get failed',42)
      if(.not.rtdb_get(irtdb,'md:nfcntr',mt_int,1,nfcntr))
     + call md_abort('md_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:icentr',mt_int,1,icentr))
     + call md_abort('md_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:nfanal',mt_int,1,nfanal))
     + call md_abort('md_rtdbin: rtdb_get failed',43)
      if(.not.rtdb_get(irtdb,'md:nfslow',mt_int,1,nfslow))
     + call md_abort('md_rtdbin: rtdb_get failed',44)
      if(.not.rtdb_get(irtdb,'md:nfoutp',mt_int,1,nfoutp))
     + call md_abort('md_rtdbin: rtdb_get failed',45)
      if(.not.rtdb_get(irtdb,'md:nfstat',mt_int,1,nfstat))
     + call md_abort('md_rtdbin: rtdb_get failed',46)
      if(.not.rtdb_get(irtdb,'md:nfrest',mt_int,1,nfrest))
     + call md_abort('md_rtdbin: rtdb_get failed',47)
      if(.not.rtdb_get(irtdb,'md:keepr',mt_int,1,keepr))
     + call md_abort('md_rtdbin: rtdb_get failed',48)
      if(.not.rtdb_get(irtdb,'md:nfcoor',mt_int,1,nfcoor))
     + call md_abort('md_rtdbin: rtdb_get failed',50)
      if(.not.rtdb_get(irtdb,'md:nfscoo',mt_int,1,nfscoo))
     + call md_abort('md_rtdbin: rtdb_get failed',51)
      if(.not.rtdb_get(irtdb,'md:nfvelo',mt_int,1,nfvelo))
     + call md_abort('md_rtdbin: rtdb_get failed',52)
      if(.not.rtdb_get(irtdb,'md:nfsvel',mt_int,1,nfsvel))
     + call md_abort('md_rtdbin: rtdb_get failed',53)
      if(.not.rtdb_get(irtdb,'md:nfforc',mt_int,1,nfforc))
     + call md_abort('md_rtdbin: rtdb_get failed',50)
      if(.not.rtdb_get(irtdb,'md:nfsfor',mt_int,1,nfsfor))
     + call md_abort('md_rtdbin: rtdb_get failed',51)
      if(.not.rtdb_get(irtdb,'md:nfprop',mt_int,1,nfprop))
     + call md_abort('md_rtdbin: rtdb_get failed',54)
      if(.not.rtdb_get(irtdb,'md:iprop',mt_int,1,iprop))
     + call md_abort('md_rtdbin: rtdb_get failed',54)
      if(.not.rtdb_get(irtdb,'md:nffree',mt_int,1,nffree))
     + call md_abort('md_rtdbin: rtdb_get failed',55)
      if(.not.rtdb_get(irtdb,'md:nfsync',mt_int,1,nfsync))
     + call md_abort('md_rtdbin: rtdb_get failed',56)
      if(.not.rtdb_get(irtdb,'md:nfauto',mt_int,1,nfauto))
     + call md_abort('md_rtdbin: rtdb_get failed',57)
      if(.not.rtdb_get(irtdb,'md:nfconv',mt_int,1,nfconv))
     + call md_abort('md_rtdbin: rtdb_get failed',58)
      if(.not.rtdb_get(irtdb,'md:nffet',mt_int,1,nffet))
     + call md_abort('md_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:impfr',mt_int,1,impfr))
     + call md_abort('md_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:impto',mt_int,1,impto))
     + call md_abort('md_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:nftri',mt_int,1,nftri))
     + call md_abort('md_rtdbin: rtdb_get failed',59)
      if(.not.rtdb_get(irtdb,'md:iformt',mt_int,1,iformt))
     + call md_abort('md_rtdbin: rtdb_get failed',60)
      if(.not.rtdb_get(irtdb,'md:madbox',mt_int,1,madbox))
     + call md_abort('md_rtdbin: rtdb_get failed',61)
      if(.not.rtdb_get(irtdb,'md:loadb',mt_int,1,loadb))
     + call md_abort('md_rtdbin: rtdb_get failed',62)
      if(.not.rtdb_get(irtdb,'md:ireset',mt_int,1,ireset))
     + call md_abort('md_rtdbin: rtdb_get failed',62)
      if(.not.rtdb_get(irtdb,'md:mropt',mt_int,1,mropt))
     + call md_abort('md_rtdbin: rtdb_get failed',63)
      if(.not.rtdb_get(irtdb,'md:idebug',mt_int,1,idebug))
     + call md_abort('md_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:icntrl',mt_int,1,icntrl))
     + call md_abort('md_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:ifidi',mt_int,1,ifidi))
     + call md_abort('md_rtdbin: rtdb_get failed',64)
      if(.not.rtdb_get(irtdb,'md:ipbtyp',mt_int,1,ipbtyp))
     + call md_abort('md_rtdbin: rtdb_get ipbtyp failed',64)
      if(.not.rtdb_get(irtdb,'md:ngl',mt_int,1,ngl))
     + call md_abort('md_rtdbin: rtdb_get failed',65)
      if(.not.rtdb_get(irtdb,'md:ifield',mt_int,1,ifield))
     + call md_abort('md_rtdbin: rtdb_get failed',66)
      if(.not.rtdb_get(irtdb,'md:ipme',mt_int,1,ipme))
     + call md_abort('md_rtdbin: rtdb_get failed',67)
      if(.not.rtdb_get(irtdb,'md:ngx',mt_int,1,ngx))
     + call md_abort('md_rtdbin: rtdb_get failed',68)
      if(.not.rtdb_get(irtdb,'md:ngy',mt_int,1,ngy))
     + call md_abort('md_rtdbin: rtdb_get failed',69)
      if(.not.rtdb_get(irtdb,'md:ngz',mt_int,1,ngz))
     + call md_abort('md_rtdbin: rtdb_get failed',70)
      if(.not.rtdb_get(irtdb,'md:numfix',mt_int,1,numfix))
     + call md_abort('md_rtdbin: rtdb_get failed',71)
c      if(.not.rtdb_get(irtdb,'md:iunfix',mt_int,1,iunfix))
c     + call md_abort('md_rtdbin: rtdb_get failed',72)
c      if(.not.rtdb_get(irtdb,'md:lsffix',mt_int,msf,lsffix))
c     + call md_abort('md_rtdbin: rtdb_get failed',72)
      if(.not.rtdb_get(irtdb,'md:noshak',mt_int,1,noshak))
     + call md_abort('md_rtdbin: rtdb_get failed',73)
      if(.not.rtdb_get(irtdb,'md:nfefld',mt_int,1,nfefld))
     + call md_abort('md_rtdbin: rtdb_get failed',74)
      if(.not.rtdb_get(irtdb,'md:nfsfld',mt_int,1,nfsfld))
     + call md_abort('md_rtdbin: rtdb_get failed',75)
      if(.not.rtdb_get(irtdb,'md:nscb',mt_int,1,nscb))
     + call md_abort('md_rtdbin: rtdb_get failed',76)
      if(nscb.gt.10) call md_abort('Increase dimension of idscb',0)
      if(.not.rtdb_get(irtdb,'md:idscb',mt_int,nscb,idscb))
     + call md_abort('md_rtdbin: rtdb_get failed',77)
      if(.not.rtdb_get(irtdb,'md:ireact',mt_int,1,ireact))
     + call md_abort('md_rtdbin: rtdb_get failed',78)
      if(.not.rtdb_get(irtdb,'md:memlim',mt_int,1,memlim))
     + call md_abort('md_rtdbin: rtdb_get failed',79)
      if(.not.rtdb_get(irtdb,'md:morder',mt_int,1,morder))
     + call md_abort('md_rtdbin: rtdb_get failed',80)
      if(.not.rtdb_get(irtdb,'md:isolvo',mt_int,1,isolvo))
     + call md_abort('md_rtdbin: rtdb_get failed',81)
      if(.not.rtdb_get(irtdb,'md:lfout6',mt_int,1,lfout6))
     + call md_abort('md_rtdbin: rtdb_get failed',82)
      if(.not.rtdb_get(irtdb,'md:iprpmf',mt_int,1,iprpmf))
     + call md_abort('md_rtdbin: rtdb_get failed',82)
      if(.not.rtdb_get(irtdb,'md:imfft',mt_int,1,imfft))
     + call md_abort('md_rtdbin: rtdb_get failed',83)
      if(.not.rtdb_get(irtdb,'md:mwmreq',mt_int,1,mwmreq))
     + call md_abort('md_rtdbin: rtdb_get failed',84)
      if(.not.rtdb_get(irtdb,'md:msareq',mt_int,1,msareq))
     + call md_abort('md_rtdbin: rtdb_get failed',85)
      if(.not.rtdb_get(irtdb,'md:mbbreq',mt_int,1,mbbreq))
     + call md_abort('md_rtdbin: rtdb_get failed',85)
      if(.not.rtdb_get(irtdb,'md:itest',mt_int,1,itest))
     + call md_abort('md_rtdbin: rtdb_get failed',86)
      if(.not.rtdb_get(irtdb,'md:nodpme',mt_int,1,nodpme))
     + call md_abort('md_rtdbin: rtdb_get failed',87)
      if(.not.rtdb_get(irtdb,'md:lbpair',mt_int,1,lbpair))
     + call md_abort('md_rtdbin: rtdb_get failed',88)
      if(.not.rtdb_get(irtdb,'md:ndistr',mt_int,1,ndistr))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:npmf',mt_int,1,npmf))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:facpmf',mt_dbl,1,facpmf))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:ndaver',mt_int,1,ndaver))
     + call md_abort('md_rtdbin: rtdb_get failed',90)
      if(.not.rtdb_get(irtdb,'md:idevel',mt_int,1,idevel))
     + call md_abort('md_rtdbin: rtdb_get failed',91)
c      if(.not.rtdb_get(irtdb,'md:itime',mt_int,mtimes,itime))
c     + call md_abort('md_rtdbin: rtdb_get failed',92)
      if(.not.rtdb_get(irtdb,'md:nftime',mt_int,1,nftime))
     + call md_abort('md_rtdbin: rtdb_get failed',92)
      if(.not.rtdb_get(irtdb,'md:nfdrss',mt_int,1,nfdrss))
     + call md_abort('md_rtdbin: rtdb_get failed',93)
      if(.not.rtdb_get(irtdb,'md:nfload',mt_int,1,nfload))
     + call md_abort('md_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:itload',mt_int,1,itload))
     + call md_abort('md_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:ioload',mt_int,1,ioload))
     + call md_abort('md_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:isload',mt_int,1,isload))
     + call md_abort('md_rtdbin: rtdb_get failed',94)
      if(.not.rtdb_get(irtdb,'md:ihess',mt_int,1,ihess))
     + call md_abort('md_rtdbin: rtdb_get failed',95)
      if(.not.rtdb_get(irtdb,'md:latom',mt_int,1,latom))
     + call md_abort('md_rtdbin: rtdb_get failed',96)
      if(.not.rtdb_get(irtdb,'md:icomm',mt_int,1,icomm))
     + call md_abort('md_rtdbin: rtdb_get failed',97)
      if(.not.rtdb_get(irtdb,'md:nfhop',mt_int,1,nfhop))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:rhop',mt_dbl,1,rhop))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:thop',mt_dbl,1,thop))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:iguide',mt_int,1,iguide))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:icmopt',mt_int,1,icmopt))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:imembr',mt_int,1,imembr))
     + call md_abort('md_rtdbin: rtdb_get failed',98)
      if(.not.rtdb_get(irtdb,'md:nopack',mt_int,1,nopack))
     + call md_abort('md_rtdbin: rtdb_get failed',99)
      if(.not.rtdb_get(irtdb,'md:iprof',mt_int,1,iprof))
     + call md_abort('md_rtdbin: rtdb_get failed',99)
      if(.not.rtdb_get(irtdb,'md:ncoll',mt_int,1,ncoll))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:ilambd',mt_int,1,ilambd))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:mlambd',mt_int,1,mlambd))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:includ',mt_int,1,includ))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:iradgy',mt_int,1,iradgy))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:idifco',mt_int,1,idifco))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nfnewf',mt_int,1,nfnewf))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nbget',mt_int,1,nbget))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
      if(.not.rtdb_get(irtdb,'md:nprec',mt_int,1,nprec))
     + call md_abort('md_rtdbin: rtdb_get failed',89)
c
      if(.not.rtdb_get(irtdb,'md:fguide',mt_dbl,1,fguide))
     + call md_abort('md_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:tguide',mt_dbl,1,tguide))
     + call md_abort('md_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:dx0sd',mt_dbl,1,dx0sd))
     + call md_abort('md_rtdbin: rtdb_get failed',100)
      if(.not.rtdb_get(irtdb,'md:dxsdmx',mt_dbl,1,dxsdmx))
     + call md_abort('md_rtdbin: rtdb_get failed',101)
      if(.not.rtdb_get(irtdb,'md:dxmsd',mt_dbl,1,dxmsd))
     + call md_abort('md_rtdbin: rtdb_get failed',102)
      if(.not.rtdb_get(irtdb,'md:dx0cg',mt_dbl,1,dx0cg))
     + call md_abort('md_rtdbin: rtdb_get failed',103)
      if(.not.rtdb_get(irtdb,'md:dxcgmx',mt_dbl,1,dxcgmx))
     + call md_abort('md_rtdbin: rtdb_get failed',104)
      if(.not.rtdb_get(irtdb,'md:dxmcg',mt_dbl,1,dxmcg))
     + call md_abort('md_rtdbin: rtdb_get failed',105)
      if(.not.rtdb_get(irtdb,'md:edacq',mt_dbl,1,edacq))
     + call md_abort('md_rtdbin: rtdb_get failed',106)
      if(.not.rtdb_get(irtdb,'md:ddacq',mt_dbl,1,ddacq))
     + call md_abort('md_rtdbin: rtdb_get failed',107)
      if(.not.rtdb_get(irtdb,'md:fdacq',mt_dbl,1,fdacq))
     + call md_abort('md_rtdbin: rtdb_get failed',108)
      if(.not.rtdb_get(irtdb,'md:delta',mt_dbl,1,delta))
     + call md_abort('md_rtdbin: rtdb_get failed',109)
      if(.not.rtdb_get(irtdb,'md:stime',mt_dbl,1,stime))
     + call md_abort('md_rtdbin: rtdb_get failed',110)
      if(.not.rtdb_get(irtdb,'md:tstep',mt_dbl,1,tstep))
     + call md_abort('md_rtdbin: rtdb_get failed',111)
      if(.not.rtdb_get(irtdb,'md:prsext',mt_dbl,1,prsext))
     + call md_abort('md_rtdbin: rtdb_get failed',112)
      if(.not.rtdb_get(irtdb,'md:prsrlx',mt_dbl,1,prsrlx))
     + call md_abort('md_rtdbin: rtdb_get failed',113)
      if(.not.rtdb_get(irtdb,'md:compr',mt_dbl,1,compr))
     + call md_abort('md_rtdbin: rtdb_get failed',114)
      if(.not.rtdb_get(irtdb,'md:tmpext',mt_dbl,1,tmpext1))
     + call md_abort('md_rtdbin: rtdb_get failed',115)
      tmpext=tmpext1
      if(.not.rtdb_get(irtdb,'md:tmpext2',mt_dbl,1,tmpext2))
     + call md_abort('md_rtdbin: rtdb_get failed',115)
      if(.not.rtdb_get(irtdb,'md:tmprlx',mt_dbl,1,tmprlx))
     + call md_abort('md_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tmsrlx',mt_dbl,1,tmsrlx))
     + call md_abort('md_rtdbin: rtdb_get failed',117)
      if(.not.rtdb_get(irtdb,'md:tann1',mt_dbl,1,tann1))
     + call md_abort('md_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tann2',mt_dbl,1,tann2))
     + call md_abort('md_rtdbin: rtdb_get failed',116)
      if(.not.rtdb_get(irtdb,'md:tgauss',mt_dbl,1,tgauss))
     + call md_abort('md_rtdbin: rtdb_get failed',118)
      if(.not.rtdb_get(irtdb,'md:frgaus',mt_dbl,1,frgaus))
     + call md_abort('md_rtdbin: rtdb_get failed',119)
      if(.not.rtdb_get(irtdb,'md:rlong',mt_dbl,1,rlong))
     + call md_abort('md_rtdbin: rtdb_get failed',120)
      if(.not.rtdb_get(irtdb,'md:rshort',mt_dbl,1,rshort))
     + call md_abort('md_rtdbin: rtdb_get failed',121)
      if(.not.rtdb_get(irtdb,'md:rqmmm',mt_dbl,1,rqmmm))
     + call md_abort('md_rtdbin: rtdb_get failed',122)
      if(.not.rtdb_get(irtdb,'md:ptol',mt_dbl,1,ptol))
     + call md_abort('md_rtdbin: rtdb_get failed',123)
      if(.not.rtdb_get(irtdb,'md:tlwsha',mt_dbl,1,tlwsha))
     + call md_abort('md_rtdbin: rtdb_get failed',124)
      if(.not.rtdb_get(irtdb,'md:tlssha',mt_dbl,1,tlssha))
     + call md_abort('md_rtdbin: rtdb_get failed',125)
      if(.not.rtdb_get(irtdb,'md:factld',mt_dbl,1,factld))
     + call md_abort('md_rtdbin: rtdb_get failed',126)
      if(.not.rtdb_get(irtdb,'md:rrdf',mt_dbl,1,rrdf))
     + call md_abort('md_rtdbin: rtdb_get failed',127)
      if(.not.rtdb_get(irtdb,'md:xfield',mt_dbl,1,xfield))
     + call md_abort('md_rtdbin: rtdb_get failed',128)
      if(.not.rtdb_get(irtdb,'md:xffreq',mt_dbl,1,xffreq))
     + call md_abort('md_rtdbin: rtdb_get failed',129)
      if(.not.rtdb_get(irtdb,'md:xfvect',mt_dbl,3,xfvect))
     + call md_abort('md_rtdbin: rtdb_get failed',130)
      if(.not.rtdb_get(irtdb,'md:weight',mt_dbl,1,weight))
     + call md_abort('md_rtdbin: rtdb_get failed',131)
      if(.not.rtdb_get(irtdb,'md:dielec',mt_dbl,1,dielec))
     + call md_abort('md_rtdbin: rtdb_get failed',132)
      if(.not.rtdb_get(irtdb,'md:ealpha',mt_dbl,1,ealpha))
     + call md_abort('md_rtdbin: rtdb_get failed',133)
      if(.not.rtdb_get(irtdb,'md:rbox',mt_dbl,1,rbox))
     + call md_abort('md_rtdbin: rtdb_get failed',134)
      if(.not.rtdb_get(irtdb,'md:rsgm',mt_dbl,1,rsgm))
     + call md_abort('md_rtdbin: rtdb_get failed',134)
      if(.not.rtdb_get(irtdb,'md:disrlx',mt_dbl,1,disrlx))
     + call md_abort('md_rtdbin: rtdb_get failed',135)
      if(.not.rtdb_get(irtdb,'md:drsscl',mt_dbl,1,drsscl))
     + call md_abort('md_rtdbin: rtdb_get failed',136)
      if(.not.rtdb_get(irtdb,'md:pmetol',mt_dbl,1,pmetol))
     + pmetol=1.0d-05
      if(.not.rtdb_get(irtdb,'md:fcoll',mt_dbl,1,fcoll))
     + fcoll=1.0d+01
      if(.not.rtdb_get(irtdb,'md:scaleq',mt_dbl,1,scaleq))
     + scaleq=-1.0d+00
c
      if(.not.rtdb_cget(irtdb,'task:theory',1,theory))
     + call md_abort('md_rtdbin: rtdg_cget failed',137)
      if(theory.eq.'md') then
      iquant=0
      iqmmm=0
      else
      iqmmm=1
      if(.not.rtdb_get(irtdb,'task:QMMM',mt_log,1,lqmmm)) iqmmm=0
      if(iqmmm.eq.1) then
      if(.not.rtdb_get(irtdb,'qmmm:uqmatm',mt_dbl,1,uqmatm))
     + uqmatm=0.0d0
      if(.not.rtdb_get(irtdb,'qmmm:linkatm',mt_int,1,linkatm))
     + linkatm=0
      if(.not.rtdb_get(irtdb,'qmmm:nobq',mt_int,1,nobq))
     + nobq=1
      uqmatm=cau2kj*uqmatm
      iquant=0
      else
      iquant=1
      endif
      endif
c
      if(.not.rtdb_cget(irtdb,'task:operation',1,operat))
c     + call md_abort('md_rtdbin: rtdg_cget failed',138)
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
     + call md_abort('md_rtdbin: failed to determine rtdb state',0)
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
c     routine md_chckin performs a check of selected input
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
      if(nserie.gt.0) call md_rdrest(lfnrst,filrst)
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
      lqmd=iquant.ne.0
      lqmmm=iqmmm.ne.0
      lpert2=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.2.or.isetp2.eq.2))
      lpert3=ntype.eq.3.or.(iset.eq.1.and.(isetp1.eq.3.or.isetp2.eq.3))
c
      lpme=ipme.gt.0
      ltwin=rlong.gt.rshort.or.lpme
c
      lpack=nopack.eq.0
      if(lqmd) lpack=.false.
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
      subroutine md_print()
c
      implicit none
c
#include "md_common.fh"
c
      character*22 ctype
c
      if(me.eq.0) then
      if(lfnout.ne.6)
     + open(unit=lfnout,file=filout(1:index(filout,' ')-1),
     + form='formatted',status='unknown')
c
      if(ntype.eq.0) then
      ctype='ENERGY EVALUATION     '
      elseif(ntype.eq.1) then
      ctype='ENERGY MINIMIZATION   '
      elseif(ntype.eq.2) then
      ctype='MOLECULAR DYNAMICS    '
      elseif(ntype.eq.3) then
      ctype='FREE ENERGY EVALUATION'
      else
      call md_abort('Unknown calculation type',ntype)
      endif
c
      call swatch(today,now)
      write(lfnout,1000) ctype,today,now
 1000 format(/,1x,a22,t110,2a10)
c
      write(lfnout,1001) titinp(1),datinp,timinp,titinp(2),titinp(3)
 1001 format(/,' Title ',t10,a,t110,2a10,/,t10,a,/,t10,a)
c
      write(lfnout,1002) filnam(1:index(filnam,' ')-1)
 1002 format(/,' System ',a)
c
      if(ntype.le.2) then
      write(lfnout,1003) iset
 1003 format(/,' Force field parameter set ',i4)
      endif
c
      if(nserie.eq.0) then
      write(lfnout,1026)
 1026 format(/,' Initial simulation',/)
      elseif(nserie.eq.1) then
      write(lfnout,1027)
 1027 format(/,' Restart simulation',/)
      else
      write(lfnout,1028)
 1028 format(/,' Continuation simulation',/)
      endif
c
      if(ntype.eq.3) then
      if(mropt.eq.0) then
      write(lfnout,1032)
 1032 format(' Multiple run initial simulation')
      elseif(mropt.eq.1) then
      write(lfnout,1033)
 1033 format(' Multiple run simulation with initial conditions')
      else
      write(lfnout,1034)
 1034 format(' Multiple run extension simulation')
      endif
      write(lfnout,1025) mrun
 1025 format(' Number of runs ',i7)
      endif
      if(ntype.ge.2) then
      write(lfnout,1024) mequi,mdacq
 1024 format(' Number of equilibration steps ',i7,/
     + ' Number of data gathering steps',i7) 
      endif
c
      if(ntype.eq.3.and.issscl.gt.0) then
      write(lfnout,1029) ssshft
 1029 format(/,' Separation shifted scaling, delta ',f12.6,' nm**2') 
      endif
c
      if(ipme.gt.0) then
      write(lfnout,1030) morder,ngx,ngy,ngz,nodpme,pmetol
 1030 format(/,' Particle-mesh Ewald, spline to order ',i5,/,
     + ' Grid size ',i5,'x',i5,'x',i5,/,
     + ' Number of processors for p-FFT ',i5,/,
     + ' Tolerance at cutoff ',1pe12.5)
      if(imfft.le.1) then
      write(lfnout,1036)
 1036 format(' Using NWChem 3D-pFFT')
      else
      write(lfnout,1037)
#if defined(ESSL)
 1037 format(' Using PESSL 3D-pFFT')
#else
 1037 format(' Using system specific 3D-pFFT')
#endif
      endif
      endif
c
      if(ipolt.eq.1) then
      write(lfnout,1031) mpolit,ptol
 1031 format(/,' Polarization model, maximum iterations ',i5,/,
     + ' Tolerance ',1pe12.5)
      endif
c
      if(.not.ltwin) then
      write(lfnout,1004) rshort
 1004 format(/,' Cutoff radius ',f12.6, ' nm')
      else
      write(lfnout,1005) rshort,rlong
 1005 format(/,' Cutoff radius short range ',f12.6, ' nm',
     +       /,' Cutoff radius long range  ',f12.6, ' nm')
      endif
c
      if(ifield.gt.0) then
      write(lfnout,1035) xfield,xffreq,xfvect
 1035 format(/,' External electrostatic field strength ',f12.6,/,
     + 25x,'frequency ',f12.6,/,
     + 25x,'field vector ',3f12.6)
      endif
c
      if(ntype.ge.2) then
      if(ivreas.ne.0) then
      if(ivreas.gt.0) then
      write(lfnout,1022) tgauss,ivreas
 1022 format(/,' Velocity reassignment temperature ',f12.6,' K',/,
     + ' Velocity reassignment frequency   ',i5,/)
      else
      write(lfnout,1122) tgauss
 1122 format(/,' Velocity reassignment temperature ',f12.6,' K',/,
     + ' Velocity reassignment in first step only',/)
      endif
      endif
      if(itscal.gt.0) then
      write(lfnout,1021) tmpext1,tmpext2,tmprlx,tmsrlx,tann1,tann2
 1021 format(/,' Isothermal ensemble external temperature ',
     + 2f12.6,' K',/,
     + ' Temperature relaxation time solvent      ',f12.6,' ps',/,
     + ' Temperature relaxation time solute       ',f12.6,' ps',/,
     + ' Temperature annealing between times      ',2f12.6,' ps')
      endif
      if(ipscal.ne.0) then
      write(lfnout,1023) prsext,prsrlx,compr
 1023 format(' Isobaric ensemble external pressure ',1pe12.5,' Pa',/,
     + ' Pressure relaxation time ',0pf12.6,' ps',/,
     + ' Compressebility ',1pe12.5,' m**2/N')
      if(ipopt.eq.12) write(lfnout,2023)
 2023 format(' Pressure scaling in x and y only')
      if(ipopt.eq.3) write(lfnout,3023)
 3023 format(' Pressure scaling in z only')
      if(ipopt.eq.123) write(lfnout,3024)
 3024 format(' Pressure scaling in x and y separate from z')
      endif
      if(iguide.gt.0) then
      write(lfnout,1038) fguide,tguide
 1038 format(/,' Self-guided molecular dynamics with force coupling of',
     + f8.5,' and relaxation time of',f8.5,' ps')
      endif
      if(mdalgo.eq.1) then
      write(lfnout,1039)
 1039 format(/,' Leap-frog integration')
      else
      write(lfnout,1040)
 1040 format(/,' Leap-frog integration (Brown-Clarke)')
      endif
      endif
c
      if(scaleq.ge.zero) then
      write(lfnout,1041) scaleq
 1041 format(/,' Solute charges scaled by factor ',f8.5)
      endif
c
      if(ntype.eq.1) then
      if(msdit.gt.0) then
      write(lfnout,1006) msdit,dx0sd,dxmsd,dxsdmx
 1006 format(/,' Maximum number of steepest descent steps ',i5,/,
     + ' Initial step size ',f12.6,' nm',/,
     + ' Maximum step size ',f12.6,' nm',/,
     + ' Minimum step size ',f12.6,' nm')
      endif
      if(mcgit.gt.0) then
      write(lfnout,1007) mcgit,ncgcy,dx0cg,dxsdmx
 1007 format(/,' Maximum number of conjugate gradient steps ',i5,/,
     + ' Number of conjugate gradient steps per cycle ',i5,/,
     + ' Initial interval size ',f12.6,' nm',/,
     + ' Minimum step size ',f12.6,' nm')
      endif
      endif
c
      write(lfnout,1008) mshitw,tlwsha,mshits,tlssha
 1008 format(/,' Maximum number of solvent SHAKE iterations ',i5,
     + ', solvent SHAKE tolerance ',f12.6,' nm',/,
     + ' Maximum number of solute SHAKE iterations  ',i5,
     + ', solute SHAKE tolerance  ',f12.6,' nm',/)
c
      if(ntype.ge.2) then
      write(lfnout,1019) nfpair
 1019 format(' Frequency update pair lists ',t55,i5)
      if(ltwin) write(lfnout,1020) nflong
 1020 format(' Frequency update long range forces ',t55,i5)
      write(lfnout,1017) nfslow
 1017 format(' Frequency removal overall motion ',t55,i5)
      write(lfnout,1018) nfcntr
 1018 format(' Frequency solute centering ',t55,i5)
      if(icentr.eq.1) write(lfnout,1818)
 1818 format(' Centering in z-direction only')
      if(icentr.eq.2) write(lfnout,1819)
 1819 format(' Centering in xyplane only')
      write(lfnout,1009) nfoutp
 1009 format(' Frequency printing step information ',t55,i5)
      write(lfnout,1010) nfstat
 1010 format(' Frequency printing statistical information ',t55,i5)
      endif
      if(ntype.eq.1) then
      write(lfnout,1011) nfqrs
 1011 format(' Frequency recording minimum energy structure ',t55,i5)
      endif
      if(ntype.ge.2) then
      write(lfnout,1012) nfrest
 1012 format(' Frequency recording restart file ',t55,i5)
      write(lfnout,1013) nfcoor
 1013 format(' Frequency recording system coordinates ',t55,i5)
      write(lfnout,1014) nfscoo
 1014 format(' Frequency recording solute coordinates ',t55,i5)
      write(lfnout,1015) nfvelo
 1015 format(' Frequency recording system velocities ',t55,i5)
      write(lfnout,1016) nfsvel
 1016 format(' Frequency recording solute velocities ',t55,i5)
      write(lfnout,1115) nfforc
 1115 format(' Frequency recording system forces ',t55,i5)
      write(lfnout,1116) nfsfor
 1116 format(' Frequency recording solute forces ',t55,i5)
      endif
c
      write(lfnout,2000)
 2000 format(/,' LOAD BALANCING',/)
c

      if(loadb.eq.0) write(lfnout,2001)
 2001 format(' None')
      if(loadb.eq.2) write(lfnout,2002)
 2002 format(' Redistribution of inter-processor box pairs')
      if(loadb.eq.1) write(lfnout,2003) factld
 2003 format(' Resizing of processor domains: smallest scaling ',
     + t55,f12.6)
      if(loadb.eq.3) write(lfnout,2004) lbpair,factld
 2004 format(' Redistribution of inter-processor box pairs: retry ',
     + t55,i5,/,' Resizing of processor domains: smallest scaling ',
     + t55,f12.6)
      if(isload.gt.0) write(lfnout,2104)
 2104 format(' Resizing in z-dimension only')
c
      if(ireset.gt.0) write(lfnout,2005)
 2005 format(/,' Load balancing reset')
c
      write(lfnout,2012) nfload
 2012 format(/,' Load balancing frequency ',i5)
c
      if(itload.eq.0) write(lfnout,2006)
 2006 format(/,' Load balancing based on last synchronization time')
      if(itload.eq.1) write(lfnout,2007)
 2007 format(/,' Load balancing based on minimum synchronization time')
      if(itload.eq.2) write(lfnout,2008)
 2008 format(/,' Load balancing based on average synchronization time')
      if(itload.eq.3) write(lfnout,2009)
 2009 format(/,' Load balancing based on average synchronization time',
     + ' with minimum synchronization time on processor 0')
c
      if(ioload.eq.1) write(lfnout,2010)
 2010 format(/,' Time for I/O counted at each step')
      if(ioload.eq.1) write(lfnout,2011)
 2011 format(/,' Experimental load balancing')
c
      if(nfhop.gt.0) then
c
      write(lfnout,3000)
 3000 format(/,' QHOP PROTON HOPPING',/)
c
      write(lfnout,3001) nfhop
 3001 format(' Frequency hopping attempts ',t55,i5)
      write(lfnout,3002) rhop
 3002 format(' Donor-acceptor cutoff distance',t55,f12.6,' nm')
      write(lfnout,3003) thop
 3003 format(' Minimum time before backhop',t55,f12.6,' ps')
c
      endif
c
      endif
c
      return
      end
      subroutine md_rdrest(lfn,fil)
c
      implicit none
c
#include "md_common.fh"
c
      integer lfn
      character*255 fil
c
      character*13 string
c
      if(me.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + status='old',form='formatted',err=9999)
      rewind(lfn)
c
    1 continue
      read(lfn,1000,end=9997) string
 1000 format(a13)
      if(string.ne.'restart input') goto 1
      read(lfn,1001,end=9998,err=9998) ntype,mdtype
 1001 format(11i7)
      read(lfn,1001) nfpair,nflong
      read(lfn,1001) lwtype,lstype,nfrest,keepr
      if(nserie.eq.1) then
      read(lfn,1002) krun,kequi,kdacq,mrun,mequi,mdacq,ldacq
 1002 format(7i7)
      else
      read(lfn,1002) krun,kequi,kdacq
      endif
      read(lfn,1003) stime,tstep
 1003 format(2f12.6)
      read(lfn,1004) rshort,rlong
 1004 format(2f12.6)
      read(lfn,1005) mshitw,tlwsha
      read(lfn,1005) mshits,tlssha
 1005 format(i7,f12.6)
      read(lfn,1006) ipscal,prsext,prsrlx,compr,ipopt
 1006 format(i5,e12.5,f12.6,e12.5,i5)
      read(lfn,1007) itscal,tmpext1,tmprlx,tmsrlx,tmpext2,
     + tann1,tann2
 1007 format(i5,6f12.6)
      read(lfn,1008) nfgaus,ivopt,tgauss,iseed
 1008 format(2i7,f12.6,i12)
      read(lfn,1009) nfoutp,nfstat,nfprop,nfnewf,ibatch
 1009 format(11i7)
      read(lfn,1009) ibinar,iformt
      read(lfn,1009) nfcoor,nfscoo,nfvelo,nfsvel,nfforc,nfsfor
      read(lfn,1009) nffree
      read(lfn,1009) nfcntr,nfslow
      read(lfn,1010) nfrdf,numrdf,ngc,ngr,ngl,ngrww,ngrsw,ngrss
 1010 format(8i7)
      read(lfn,1011) rrdf,drdf
 1011 format(2f12.6)
      read(lfn,1012) numdis,lendis
 1012 format(2i7)
      read(lfn,1013) numhis,lenhis
 1013 format(11i7)
      read(lfn,1014) nfdip,ndip,rdip
 1014 format(2i7,f12.6)
      read(lfn,1015) nfkirk,nkirk,rkirk
 1015 format(2i7,f12.6)
      endif
c
 9997 continue
      if(me.eq.0) close(unit=lfn)
      return
c
 9998 continue
      call md_abort('Unable to read restart file in md_rest ',me)
      return
 9999 continue
      call md_abort('Unable to open restart file in md_rest ',me)
      return
      end
      subroutine md_wtrest(lfn)
c
      implicit none
c
#include "md_common.fh"
c
      integer lfn
c
      if(me.ne.0) return
c
      write(lfn,1000)
 1000 format('restart input')
      write(lfn,1001) ntype,mdtype
 1001 format(11i7)
      write(lfn,1001) nfpair,nflong
      write(lfn,1001) lwtype,lstype,nfrest,keepr
      write(lfn,1002) irun,iequi,idacq,mrun,mequi,mdacq,ldacq
 1002 format(11i7)
      write(lfn,1003) stime,tstep
 1003 format(2f12.6)
      write(lfn,1004) rshort,rlong
 1004 format(2f12.6)
      write(lfn,1005) mshitw,tlwsha
      write(lfn,1005) mshits,tlssha
 1005 format(i7,f12.6)
      write(lfn,1006) ipscal,prsext,prsrlx,compr,ipopt
 1006 format(i5,e12.5,f12.6,e12.5,i5)
      write(lfn,1007) itscal,tmpext1,tmprlx,tmsrlx,tmpext2,
     + tann1,tann2
 1007 format(i5,6f12.6)
      tmpext=tmpext1
      write(lfn,1008) nfgaus,ivopt,tgauss,iseed
 1008 format(2i7,f12.6,i12)
      write(lfn,1009) nfoutp,nfstat,nfprop,nfnewf,ibatch
 1009 format(11i7)
      write(lfn,1009) ibinar,iformt
      write(lfn,1009) nfcoor,nfscoo,nfvelo,nfsvel,nfforc,nfsfor
      write(lfn,1009) nffree
      write(lfn,1009) nfcntr,nfslow
      write(lfn,1010) nfrdf,numrdf,ngc,ngr,ngl,ngrww,ngrsw,ngrss
 1010 format(8i7)
      write(lfn,1011) rrdf,drdf
 1011 format(2f12.6)
      write(lfn,1012) numdis,lendis
 1012 format(2i7)
      write(lfn,1013) numhis,lenhis
 1013 format(11i7)
      write(lfn,1014) nfdip,ndip,rdip
 1014 format(2i7,f12.6)
      write(lfn,1015) nfkirk,nkirk,rkirk
 1015 format(2i7,f12.6)
c
      return
      end
      subroutine md_fopen(lclose)
c
      implicit none
c
#include "md_common.fh"
c
      logical lclose
c
      if(nfnewf.gt.0) then
      ibatch=ibatch+1
      write(filtrj,1000) root(1:index(root,' ')-1),
     + ibatch,'.trj'
      write(filhop,1000) root(1:index(root,' ')-1),
     + ibatch,'.hop'
      write(filprp,1000) root(1:index(root,' ')-1),
     + ibatch,'.prp'
      write(filpmf,1000) root(1:index(root,' ')-1),
     + ibatch,'.pmf'
      write(filtim,1000) root(1:index(root,' ')-1),
     + ibatch,'.tim'
      write(rfile,1000) root(1:index(root,' ')-1),
     + ibatch,'.rst'
 1000 format(a,i3.3,a)
      endif
c
      if(ntype.ne.3) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      if(lclose) close(lfntrj)
      open(unit=lfntrj,file=filtrj(1:index(filtrj,' ')-1),
     + form='formatted',status='unknown')
      call cf_trjhdr(lfntrj)
      endif
      if(nfhop.gt.0) then
      if(lclose) close(lfnhop)
      open(unit=lfnhop,file=filhop(1:index(filhop,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nfprop.gt.0) then
      if(lclose) close(lfnprp)
      open(unit=lfnprp,file=filprp(1:index(filprp,' ')-1),
     + form='formatted',status='unknown')
      call prp_header()
      endif
      endif
      if(ntype.eq.2.and.iprpmf.ne.0) then
      if(lclose) close(lfnpmf)
      open(unit=iabs(lfnpmf),file=filpmf(1:index(filpmf,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nftime.gt.0) then
      if(lclose) close(lfntim)
      open(unit=lfntim,file=filtim(1:index(filtim,' ')-1),
     + form='formatted',status='unknown')
      call md_hdrtim()
      endif
c
      return
      end
      subroutine md_hdrtim()
c
      implicit none
c
#include "md_common.fh"
c
      write(lfntim,1200) np,56
 1200 format(2i5)
      write(lfntim,1201)
      write(lfntim,1202)
      write(lfntim,1203)
      write(lfntim,1204)
      write(lfntim,1205)
      write(lfntim,1206)
      write(lfntim,1207)
      write(lfntim,1208)
      write(lfntim,1209)
      write(lfntim,1210)
      write(lfntim,1211)
      write(lfntim,1212)
      write(lfntim,1213)
      write(lfntim,1214)
      write(lfntim,1215)
      write(lfntim,1216)
      write(lfntim,1217)
      write(lfntim,1218)
      write(lfntim,1219)
      write(lfntim,1220)
      write(lfntim,1221)
      write(lfntim,1222)
      write(lfntim,1223)
      write(lfntim,1224)
      write(lfntim,1225)
      write(lfntim,1226)
      write(lfntim,1227)
      write(lfntim,1228)
      write(lfntim,1229)
      write(lfntim,1230)
      write(lfntim,1231)
      write(lfntim,1232)
      write(lfntim,1233)
      write(lfntim,1234)
      write(lfntim,1235)
      write(lfntim,1236)
      write(lfntim,1237)
      write(lfntim,1238)
      write(lfntim,1239)
      write(lfntim,1240)
      write(lfntim,1241)
      write(lfntim,1242)
      write(lfntim,1243)
      write(lfntim,1244)
      write(lfntim,1245)
      write(lfntim,1246)
      write(lfntim,1247)
      write(lfntim,1248)
      write(lfntim,1249)
      write(lfntim,1250)
      write(lfntim,1251)
      write(lfntim,1252)
      write(lfntim,1253)
      write(lfntim,1254)
      write(lfntim,1255)
      write(lfntim,1256)
 1201 format('  1 : Initialization velocities')
 1202 format('  2 : Periodic Boundary Conditions')
 1203 format('  3 : Dynamic Load Balancing')
 1204 format('  4 : Atom Redistribution')
 1205 format('  5 : Center of Mass Coordinates')
 1206 format('  6 : Recording Trajectory')
 1207 format('  7 : CAFE Initialization')
 1208 format('  8 : QMD Forces Evaluation')
 1209 format('  9 : Synchronization')
 1210 format(' 10 : SPACE Initialization')
 1211 format(' 11 : Synchronization')
 1212 format(' 12 : Collapse Option')
 1213 format(' 13 : External Fields')
 1214 format(' 14 : Multiprocs Data Initialization')
 1215 format(' 15 : Multiprocs Data Communication')
 1216 format(' 16 : Multiprocs Data Forces Evaluation')
 1217 format(' 17 : Restraints Data Initialization')
 1218 format(' 18 : Restraints Data Communication')
 1219 format(' 19 : Restraints Data Forces Evaluation')
 1220 format(' 20 : PMF Data Initialization')
 1221 format(' 21 : PMF Data Communication')
 1222 format(' 22 : PMF Data Forces Evaluation')
 1223 format(' 23 : Induced Dipoles Evaluation')
 1224 format(' 24 : PME Initialization')
 1225 format(' 25 : PME Charge Grid Evaluation')
 1226 format(' 26 : PME Charge Grid Communication')
 1227 format(' 27 : PME Charge Grid Retrieval')
 1228 format(' 28 : PME Reverse 3D-pFFT')
 1229 format(' 29 : PME Energy Evaluation')
 1230 format(' 30 : PME Synchronization')
 1231 format(' 31 : PME Forward 3D-pFFT')
 1232 format(' 32 : PME Synchronization')
 1233 format(' 33 : Local Coordinates Retrieval')
 1234 format(' 34 : Local Forces Evaluation')
 1235 format(' 35 : Local Force Accumulation')
 1236 format(' 36 : Non-Local Coordinates Retrieval')
 1237 format(' 37 : Non-Local Forces Evaluation')
 1238 format(' 38 : Non-Local Force Accumulation')
 1239 format(' 39 : PME Wait')
 1240 format(' 40 : PME Forces Evaluation')
 1241 format(' 41 : PME Forces Communication')
 1242 format(' 42 : PME Node Synchronization')
 1243 format(' 43 : PME Flag')
 1244 format(' 44 : Synchronization')
 1245 format(' 45 : SPACE Finalization')
 1246 format(' 46 : QM/MM Forces Evaluation')
 1247 format(' 47 : Forces Normalization')
 1248 format(' 48 : Guided MD Forces')
 1249 format(' 49 : MD Time Step')
 1250 format(' 50 : SHAKE')
 1251 format(' 51 : CAFE Finalization')
 1252 format(' 52 : Center of Mass Motion Removal')
 1253 format(' 53 : SPACE Update')
 1254 format(' 54 : Property Evaluation')
 1255 format(' 55 : Restart File Update')
 1256 format(' 56 : TOTAL TIME PER STEP')
c          
      return
      end
      subroutine md_partition
c
      implicit none
c
#include "md_common.fh"
#include "global.fh"
#include "util.fh"
c  size should equal maximum number of processors that could be in
c  a group
      integer MAX_GRP_SIZE
      parameter (MAX_GRP_SIZE = 1000)
      integer grp_size, nproc
      integer list(MAX_GRP_SIZE)
      integer i, grp_handle
      real*8 rantmp
c
c      write(*,1000) npg
c 1000 format(' The number of processor groups requested is ',i5)
c      write(*,1001) ga_pgroup_get_world()
c 1001 format(' The world group handle is ',i5)
c      write(*,1002) ga_pgroup_get_default()
c 1002 format(' The default group handle is ',i5)
c
      nproc = ga_nnodes()
      me = ga_nodeid()
      rantmp=util_random(iseed*(me+1))
      if (mod(nproc,npg).ne.0) then
      call md_abort('Number of processor groups is not a multiple',
     + ' of number of processors',ga_nodeid())
      endif
      grp_size = nproc/npg
      ipg = (me - mod(me,grp_size))/grp_size
      do i = 1, grp_size
        list(i) =ipg*grp_size + i - 1
      end do
      grp_handle = ga_pgroup_create(list,grp_size)
      call ga_pgroup_set_default(grp_handle)
      np = ga_nnodes()
      me = ga_nodeid()
      return
      end
c
      subroutine md_partition_end
c
      implicit none
c
#include "md_common.fh"
#include "global.fh"
      integer grp_handle
      logical status
c
      grp_handle = ga_pgroup_get_default()
      call ga_pgroup_set_default(ga_pgroup_get_world())
      status = ga_pgroup_destroy(grp_handle)
c
      call ga_sync()
c
      np = ga_nnodes()
      me = ga_nodeid()
      return
      end
