      subroutine argos_start
c
c $Id$
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
c     determine the partition of processes, process ids, etc
c     ------------------------------------------------------
c
      call argos_partition
c
c     generate file names and open all files
c     --------------------------------------
c
      call argos_fopen(.false.)
c
c     print input information
c     -----------------------
c
      call argos_print()
c
c     start the spacial decomposition API
c     -----------------------------------
c
      call argos_space_start(lfnout,lfntop,filtop,lfnrst,filrst,lfnsyn,
     + filsyn,
     + nfsync,rshort,rlong,zero,rsgm,
     + npx,npy,npz,nbx,nby,nbz,
     + npbtyp,nbxtyp,box,vlat,lpbc,
     + nwm,mwm,nwa,mwa,nsf,msf,nsm,msm,nsa,msa,
     + loadb,lbpair,factld,ipolt.ne.0,.false.,temp,tempw,temps,
     + iguide,lfndbg,idebug,projct,mbbreq,nserie,isload,ireset,icntrl,
     + nseq,i_lseq,ndums,nbget,nprec,madbox)
c
c     start the analysis API
c     ----------------------
c
c      if(nfanal.gt.0) call ana_init(nsa,msa,.false.)
c
c     start the classical forces API
c     ------------------------------
c
      call argos_cafe_start(irtdb,lfnout,lfntop,filtop,ndistr,npmf,
     + npmfi,
     + nwm,mwm,nwa,mwa,nsf,msf,nsm,msm,nsa,msa,
     + mdalgo,npbtyp,nbxtyp,rshort,rlong,box,
     + ipme,morder,ngx,ngy,ngz,nodpme,pmetol,
     + tstep,tlwsha,mshitw,mshits,noshak,
     + ipolt,itscal,ipscal,ipopt,tmpext,prsext,
     + tmprlx,tmsrlx,prsrlx,scaleq,facpmf,
     + 0,temp,tempw,temps,compr,ntype,iset,isetp1,isetp2,
     + issscl,delta,nfanal,lpbc,npgdec,xfield,xfvect,xffreq,
     + npener,icntrl,nbias,mropt,includ,ltwin,
     + nseq,i_lseq,nfhop,rhop,thop,ndums,ipbtyp,lfnhop,iradgy,
     + 0,nprec)
c     + nbget)
c
      if(mlambd.gt.0.and.ilambd.gt.0.and.ilambd.le.mlambd)
     + call argos_cafe_lambda(lamtyp,irun,maxlam,elam,lfnout,lfnpmf,
     + rlambd,dlambd,filnam)
c
c     print topology information
c     --------------------------
c
      call argos_cafe_print_top(lfnout,npatom,nptopw,nptops)
c
      msm=max(1,nsm)
      msf=max(1,nsf)
      mst=max(msm,nseq)
c
c     start the property API
c     ----------------------
c
      call argos_prop_start(nserie,ntype,nftri,lfnrst,filrst,lfnout,
     + lfnprp,lfngib,nfoutp,nfstat,nfprop,iprop,
     + .true.,.true.,ltwin,ipolt.ne.0,ipme.ne.0,
     + npstep.ne.0,npener.ne.0,npstat,
     + nwm,msf,nsf,mpe,mdacq,mrun,iset,isetp1,isetp2,tstep,msm,nsm,
     + nsa,ddacq,edacq,iprof,npmf,npener,ndistr,lpbc,nbias,nodpme,npmfi,
     + iguide.ne.0,iradgy,idifco,nbget,meg,npg)
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
      call argos_space_setup(me,int_mb(i_iw),dbl_mb(i_xw),
     + dbl_mb(i_xwcr),
     + dbl_mb(i_vw),dbl_mb(i_gw),nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_xscr),dbl_mb(i_vs),
     + dbl_mb(i_gs),nsaloc,lpack)
      call argos_space_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
c
c     initialize packing
c     ------------------
c
      if(lpack)
     + call argos_space_pack_init(int_mb(i_is),nsaloc,int_mb(i_iw),
     + nwmloc)
c
c     spacial decomposition
c     ---------------------
c
      call argos_space_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),
     + dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
c     calculate mass factors
c     ----------------------
c
      call argos_cafe_weight(nwmloc,nsaloc,int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lshop-1)*msa),wbox)
c
c     calculate centers of mass
c     -------------------------
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
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
      call argos_cafe_fix(lfnout,lfncmd,numfix,
     + int_mb(i_iw+(lwgmn-1)*mwm),
     + int_mb(i_iw+(lwdyn-1)*mwm),nwmloc,
     + int_mb(i_is+(lsgan-1)*msa),int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lssgm-1)*msa),nsaloc,
     + dbl_mb(i_xwm),dbl_mb(i_xs))
      if(me.eq.0.and.numfix.gt.0) then
      close(unit=lfncmd)
      endif
      call argos_space_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
c
c     print decomposition information
c     -------------------------------
c
      call argos_space_print()
c
      call argos_prop_setup(wbox)
c
c     write file headers
c     ------------------
c
      if(me.eq.0.and.ntype.ne.3) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0)
     + call argos_cafe_trjhdr(lfntrj)
      endif
c
      if(npg.gt.1.and.me.eq.0) then
c      root=filnam(1:ibl)//gid
c      filrst=filnam(1:ibl)//gid//'.rst'
c      rfile=filrst
      filrst=rfile
      endif
c
      if(.not.ma_verify_allocator_stuff())
     + call md_abort('ma problems at end of argos_start',me)
      return
      end
