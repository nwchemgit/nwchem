      subroutine argos_newton()
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
c
      logical frequency
      external frequency
c
      external timer_wall,timer_wall_average,timer_wall_minimum
      real*8 timer_wall,timer_wall_average,timer_wall_minimum
      external timer_wall_total
      real*8 timer_wall_total
c
      real*8 dxmax,wallt,syntim
c
c     start timer
c
      call timer_start(1)
c
c     debug code
c
      if(idebug.gt.0) then
      write(lfndbg,'(a,i5,f12.6)') 'Timestep ',mdstep,stime
      endif
c
c     end debug code
c
      if(tann2.gt.0.and.tmpext1.ne.tmpext2.and.tann1.lt.tann2)then
      if(stime.gt.tann1.and.stime.lt.tann2) then
      tmpext=((tann2-stime)*tmpext1+(stime-tann1)*tmpext2)/(tann2-tann1)
      else
      tmpext=tmpext1
      if(stime.ge.tann2) tmpext=tmpext2
      endif
      endif
c
      wallt=timer_wall(203)
c
c     dynamic load balancing data
c
      if(itload.eq.0) then
      wallt=timer_wall(203)
      syntim=timer_wall(204)
      elseif(itload.eq.1) then
      wallt=timer_wall_minimum(203)
      syntim=timer_wall_minimum(204)
      elseif(itload.eq.2) then
      wallt=timer_wall_average(203)
      syntim=timer_wall_average(204)
      else
      wallt=timer_wall_average(203)
      syntim=timer_wall_average(204)
      if(me.eq.0) then
      wallt=timer_wall_minimum(203)
      syntim=timer_wall_minimum(204)
      endif
      endif
c
      if(me.eq.0.and.ioload.eq.1) then
      if(corrio.ge.syntim) then
      syntim=zero
      else
      syntim=syntim-corrio
      endif
      endif
c
      if(lpair) call timer_reset(203)
      call timer_start(203)
c
c     reassign velocities
c
      if(frequency(mdstep,nfgaus)) then
      call argos_cafe_gauss(tgauss,frgaus,nwmloc,nsaloc,
     + dbl_mb(i_vw),dbl_mb(i_vs),int_mb(i_iw+(lwdyn-1)*mwm),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      endif
c
      lpair=lfirst.or.frequency(mdstep,nfpair)
      lload=lfirst.or.frequency(mdstep,nfload)
      lhop=frequency(mdstep,nfhop)
      llong=(lfirst.or.frequency(mdstep,nflong).or.lpair).and.ltwin
c
      call timer_stop(1)
c
      if(lpair) then
c
      call timer_start(2)
c
c     center of mass coordinates
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     periodic boundary conditions
c
      call argos_fold(int_mb(i_iw),int_mb(i_is),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_xsm))
c
      call timer_stop(2)
c
      if(lload) then
c
      call timer_start(3)
c
      if(me.eq.0.and.ioload.eq.2) then
      if(corrio.ge.syntim.and.
     + ((nfcoor.gt.nfpair.and.nfcoor.gt.np).or.
     +  (nfcoor.eq.0.and.nfscoo.gt.nfpair.and.nfscoo.gt.np))) then
      syntim=zero
      else
      syntim=syntim-corrio
      endif
      endif
c
      call argos_space_balanc(stime,syntim,wallt,frequency(mdstep,
     + nfsync))
      call timer_reset(204)
c
      call timer_stop(3)
      call timer_start(4)
c
c     atom redistribution
c
      call argos_space_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),
     + dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
      call timer_stop(4)
c
      endif
c
      endif
c
      call timer_start(5)
c
c     center of mass coordinates
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     subtract center of mass coordinates from reference coordinates
c
      if(idifco.gt.0)
     + call argos_addref(.false.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
      call timer_stop(5)
      call timer_start(6)
c
      if(lfw.or.lfs) then
      call argos_space_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),
     + nsaloc)
      call argos_space_wrttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
c
      call timer_stop(6)
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
c
      corrio=costio
      if(me.eq.0.and..not.lequi.and.(lxw.or.lvw.or.lxs.or.lvs)) then
c
      call timer_start(6)
c
      if(.not.lfw.and..not.lfs) then
      call argos_space_wrttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
      call timer_stop(6)
c
      costio=max(costio,timer_wall(6))
      corrio=costio-timer_wall(6)
c
      endif
c
c     atomic forces and potential energies
c
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     self-guided forces
c
      if(iguide.gt.0) then
      call timer_start(48)
      call argos_guided(dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_gw),dbl_mb(i_gs))
      call timer_stop(48)
      endif
c
c     center of mass options
c
      if(icmopt.gt.0) then
      call timer_start(48)
      call argos_cmopt(dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_fcm),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      call timer_stop(48)
      endif
c
      if(imembr.ne.0) call argos_membrane_forces(int_mb(i_mm),
     + dbl_mb(i_fm),dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_fs),
     + dbl_mb(i_wws))
c
c     time step
c
      call timer_start(49)
      call argos_cafe_mdstep(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_vwt),dbl_mb(i_fw),nwmloc,
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs),dbl_mb(i_vst),
     + dbl_mb(i_fs),nsaloc,int_mb(i_iw+(lwgmn-1)*mwm),
     + int_mb(i_is+(lsgan-1)*msa),int_mb(i_is+(lssgm-1)*msa),tmpext,
     + int_mb(i_is+(lshop-1)*msa))
      call timer_stop(49)
c
c     shake
c
      call timer_start(50)
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
      call timer_stop(50)
c
c     recalculate velocities
c
      call argos_cafe_veloc(nwmloc,dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_vw),
     + nsaloc,dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs))
c
c     velocity scaling to preset temperature
c
      if(frequency(mdstep,nfgaus)) then
      call argos_cafe_vscale(tgauss,nwmloc,nsaloc,
     + dbl_mb(i_vw),dbl_mb(i_vwt),dbl_mb(i_vs),dbl_mb(i_vst),
     + int_mb(i_iw+(lwdyn-1)*mwm),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      if(iand(ivopt,1).eq.1) nfgaus=0
      endif
c
c     coordinate scaling
c
      call timer_start(51)
      call argos_cafe_final(dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_yw),
     + dbl_mb(i_vw),dbl_mb(i_vwt),nwmloc,
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_ys),dbl_mb(i_vs),
     + dbl_mb(i_vst),int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + int_mb(i_is+(lsfrc-1)*msa),int_mb(i_is+(lshop-1)*msa),
     + dbl_mb(i_zs),
     + dbl_mb(i_esk),nsaloc,box,vlat,pres,temp,tempw,temps)
      call timer_stop(51)
c
c     center of mass coordinates
c
      call timer_start(52)
c
      if(idifco.gt.0.or.
     + frequency(mdstep,nfcntr).or.frequency(mdstep,nfslow)) then
c
      call argos_cafe_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     add center of mass coordinates from reference coordinates
c
      if(idifco.gt.0)
     + call argos_addref(.true.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
c     center solute in box
c
      if(frequency(mdstep,nfcntr)) then
      call argos_cafe_center(dbl_mb(i_xw),nwmloc,
     + int_mb(i_is+(lsfrc-1)*msa),
     + dbl_mb(i_xs),nsaloc,idscb,nscb,icentr)
      endif
c
c     remove overall translational motion
c
      if(frequency(mdstep,nfslow)) then
      call argos_cafe_slow(dbl_mb(i_xw),dbl_mb(i_vw),nwmloc,
     + dbl_mb(i_xs),
     + dbl_mb(i_vs),int_mb(i_is+(lsatt-1)*msa),nsaloc)
      endif
c
      endif
      call timer_stop(52)
      call timer_start(53)
c
c     update decomposition module
c
      call argos_space_update(me,vlat,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_xwcr),dbl_mb(i_vw),nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),nsaloc)
c
      if(itest.gt.0) call argos_test()
c
      lfirst=.false.
c
      call timer_stop(53)
      call timer_stop(203)
c
      return
      end
c $Id$
