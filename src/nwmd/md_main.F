      subroutine md_main()
c
      implicit none
c
#include "md_common.fh"
c
      if(ntype.eq.0) then
c
c     single energy
c     -------------
c
      call md_sp()
c
      elseif(ntype.eq.1) then
c     energy minimization
c     -------------------
c
      call md_em()
c
      elseif(ntype.eq.2) then
c     molecular dynamics
c     ------------------
c
      call md_md()
c
      elseif(ntype.eq.3) then
c
c     free energy simulation
c     ----------------------
c
      call md_ti()
c
      else
      call errquit('Unknown calculation type',ntype)
      endif
c
      return
      end
      subroutine md_sp()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      return
      end
      subroutine md_em()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      integer i_pcgw,l_pcgw,i_pcgs,l_pcgs
c
      llong=ltwin
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     periodic boundary conditions
c
      call md_fold(int_mb(i_iw),int_mb(i_is),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     atom redistribution
c
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + int_mb(i_is),nsaloc)
c
      if(msdit.gt.0) call md_stdesc()
c
      if(mcgit.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'pcgw',l_pcgw,i_pcgw))
     + call errquit('Failed to allocate memory for pcgw',me)
      if(.not.ma_push_get(mt_dbl,3*msa,'pcgs',l_pcgs,i_pcgs))
     + call errquit('Failed to allocate memory for pcgs',me)
      call md_congra(dbl_mb(i_pcgw),dbl_mb(i_pcgs))
      if(.not.ma_pop_stack(l_pcgs))
     + call errquit('Failed to deallocate memory for pcgs',me)
      if(.not.ma_pop_stack(l_pcgw))
     + call errquit('Failed to deallocate memory for pcgw',me)
      endif
c
      return
      end
      subroutine md_stdesc()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      real*8 epot,edif,epsd,dx,dxf,dxmax
      character*1 cqrs
c
      isdit=1
c
c     initialize cf
c
      call cf_init(lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),eww,
     + dbl_mb(i_esw),dbl_mb(i_ess))
c
c     atomic forces and potential energies
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),epme,epot)
c
      epsd=epot
      dx=dx0sd
      mdstep=0
c
    1 continue
c
      dxf=dx/fmax
c
      call cf_emstep(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_fw),dbl_mb(i_vw),nwmloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_fs),dbl_mb(i_vs),nsaloc,dxf)
c
c     shake
c
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      isdit=isdit+1
c
      lpair=mod(isdit,nfpair).eq.0
      llong=(mod(isdit,nflong).eq.0.or.lpair).and.ltwin
c
c     initialize cf
c
      call cf_init(lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),eww,
     + dbl_mb(i_esw),dbl_mb(i_ess))
c
c     atomic forces and potential energies
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),epme,epot)
c
c     if energy goes up restore coordinates and forces
c
      if(epot.gt.epsd.and.isdit.lt.msdit.and.dxmax.gt.dxsdmx) then
c
      call cf_restor(dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_fw),
     + dbl_mb(i_vw),nwmloc,dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_fs),
     + dbl_mb(i_vs),nsaloc)
c
      dx=half*dx
      goto 1
      endif
c
      edif=epot-epsd
      epsd=epot
c
      cqrs=' '
      if(mod(isdit,nfqrs).eq.0) then
      call sp_wrtrst(lfnqrs,filqrs,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      cqrs='*'
      endif
c
      if(me.eq.0) then
      write(lfnout,600) isdit,cqrs,edif,epsd,dxmax
  600 format(' sd',i7,1x,a1,3x,2(1pe13.5),0pf12.8)
      endif
c
      dx=min(1.2d0*dx,dxmsd)
c
      if(dxmax.gt.dxsdmx.and.isdit.lt.msdit) goto 1

      return
      end
      subroutine md_congra(pcgw,pcgs)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      real*8 pcgw(mwm,3,mwa),pcgs(msa,3)
c
      real*8 epot,edif,epcg,dx,dxf,dxmax,beta1
      character*1 cqrs
c
      icgit=1
c
c     initialize cf
c
      call cf_init(lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),eww,
     + dbl_mb(i_esw),dbl_mb(i_ess))
c
c     atomic forces and potential energies
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),epme,epot)
c
      epcg=epot
      dx=dx0cg
      beta1=zero
      mdstep=0
c
      return
      end
      subroutine md_md()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer iequi
c
c     equilibration
c
      do 1 iequi=kequi+1,mequi
      call md_newton()
    1 continue
c
c     data gathering
c
      do 2 iequi=kdacq+1,mdacq
      call timer_start(1)
c
      call md_newton()
c
      call prp_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme)
c
      if(mod(mdstep,nfrest).eq.0)
     + call sp_wrtrst(lfnrst,filrst,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
c
      call timer_stop(1)
    2 continue
c
      return
      end
      subroutine md_ti()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      return
      end
      subroutine md_newton()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      external timer_wall
      real*8 timer_wall
c
      real*8 dxmax
c
      mdstep=mdstep+1
      stime=stime+tstep
c
c     reassign velocities
c
c     if(mod(mdstep,nfgaus).eq.0)
c    + call cf_gauss(nwmloc,nsaloc,vw,vs,tgaus)
c
      lpair=mdstep.eq.1.or.mod(mdstep,nfpair).eq.0
      llong=(mdstep.eq.1.or.mod(mdstep,nflong).eq.0.or.lpair).and.ltwin
c
      if(lpair) then
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     periodic boundary conditions
c
      call md_fold(int_mb(i_iw),int_mb(i_is),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     dynamic load balancing
c
      call sp_balanc(timer_wall(2))
c
c     atom redistribution
c
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + int_mb(i_is),nsaloc)
c
      endif
c
c     initialize cf
c
      call cf_init(lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),eww,
     + dbl_mb(i_esw),dbl_mb(i_ess))
c
c     atomic forces and potential energies
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     time step
c
      call cf_mdstep(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_vwt),dbl_mb(i_fw),nwmloc,
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs),dbl_mb(i_vst),
     + dbl_mb(i_fs),nsaloc,int_mb(i_iw+(lwgmn-1)*mwm),
     + int_mb(i_is+(lsgan-1)*msa))
c
c     shake
c
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     coordinate scaling
c
      call cf_final(dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_yw),
     + dbl_mb(i_vw),dbl_mb(i_vwt),nwmloc,
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_ys),dbl_mb(i_vs),
     + dbl_mb(i_vst),int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + int_mb(i_is+(lsfrc-1)*msa),dbl_mb(i_zs),
     + dbl_mb(i_esk),nsaloc,box,pres,temp,tempw,temps)
c
c     update decomposition module
c
      call sp_update(me,vlat,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_xwcr),dbl_mb(i_vw),nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),nsaloc)
c
      return
      end
      subroutine md_shake(xw,yw,iwl,xs,ys,isl,dmax)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),xs(msa,3),ys(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 dmax
c
      integer i,j,lhandl,ibbl,iwfr,iwto,isfr,isto,nbbl
      logical lself
c
      dmax=zero
c
      if(nwmloc.gt.0) then
      call cf_shakew(xw,yw,iwl(1,lwgmn),iwl(1,lwdyn),nwmloc)
      do 1 j=1,nwa
      do 2 i=1,nwmloc
      dmax=max(dmax,(xw(i,1,j)-yw(i,1,j))**2+
     + (xw(i,2,j)-yw(i,2,j))**2+(xw(i,3,j)-yw(i,3,j))**2)
    2 continue
    1 continue
      endif
c
      if(nsaloc.gt.0) then
      call sp_nbbl(nbbl)
      do 3 ibbl=1,nbbl
      call sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
      if(lself) call cf_shakes(lhandl,xs,ys,isl(1,lsgan),isl(1,lsatt),
     + isl(1,lsdyn),isfr,isto)
    3 continue
      do 4 i=1,nsaloc
      dmax=max(dmax,(xs(i,1)-ys(i,1))**2+(xs(i,2)-ys(i,2))**2+
     + (xs(i,3)-ys(i,3))**2)
    4 continue
      endif
c
      dmax=sqrt(dmax)
      call ga_dgop(mrg_d45,dmax,1,'max')
c
      return
      end
      subroutine md_fold(iwl,isl,xw,xwm,xs,xsm)
c
      implicit none
c
#include "md_common.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3),xsm(msm,3)
c
      call cf_fold(nwmloc,xw,xwm,nsaloc,isl(1,lsatt),isl(1,lsmol),
     + xs,xsm)
c
      return
      end
      subroutine md_forces(iwl,isl,xw,xwm,xs,fw,fs)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
c
      integer i,j,k
c
      call ga_sync()
c
      call sp_initf(fw,fs,llong)
c
      call sp_putix(me,iwl,xw,nwmloc,isl,xs,nsaloc)
c
      call ga_sync()
c
      call md_fclass(iwl,isl,xw,xwm,xs,fw,fs)
c
      call timer_start(2)
      call ga_sync()
      call timer_stop(2)
c
      call sp_final(fw,fs)
c
      if(ltwin) then
      do 1 j=1,3
      do 2 k=1,nwa
      do 3 i=1,nwmloc
      fw(i,j,k,1)=fw(i,j,k,1)+fw(i,j,k,2)
    3 continue
    2 continue
      do 4 i=1,nsaloc
      fs(i,j,1)=fs(i,j,1)+fs(i,j,2)
    4 continue
    1 continue
      endif
c
      call cf_fnorm(fw,nwmloc,isl(1,lsatt),fs,nsaloc,fnorm,fmax)
c
      return
      end
      subroutine md_fclass(iwl,isl,xw,xwm,xs,fw,fs)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lself
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer nbbl
c
      call cf_multi(nsaloc,dbl_mb(i_xs),dbl_mb(i_fs),isl(1,lsgan),
     + isl(1,lsfrc),isl(1,lsdyn),isl(1,lsct1),dbl_mb(i_ess))
c
      if(iewald.gt.0) then
      call pme_chgrid(dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),isl(1,lsct1),
     + nsaloc,epme)
      endif
c
      call sp_nbbl(nbbl)
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + isl(1,lsatt),isl(1,lsmol),dbl_mb(i_xs),dbl_mb(i_xsm))
c
      lpbc=.true.
c
      do 1 ibbl=1,nbbl
c
      call sp_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,lself)
c
      call cf_comw(xw,xwm,jwfr,jwto)
c
      call forces(lself,xw,xwm,fw,zw,dbl_mb(i_rtos),iwl(1,lwdyn),
     + int_mb(i_iwz),
     + iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),fs,dbl_mb(i_zs),
     + isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),isl(1,lsfrc),
     + isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),isl(1,lsct3),
     + int_mb(i_isz),isfr,isto,jsfr,jsto,lpbc,lhandl,
     + .true.,eww,dbl_mb(i_esw),dbl_mb(i_ess))
c
      call sp_accfbl(ibbl,lhandl,fw,fs)
c
    1 continue
c
      if(iewald.gt.0) then
      call pme_forces(fw,nwmloc,fs,isl(1,lsct1),nsaloc)
      endif
c
      return
      end

