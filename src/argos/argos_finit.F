      subroutine argos_finit(iwl,isl,xw,xwm,xs,fw,fs,xsm,xsmp)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      real*8 xsm(msm,3),xsmp(msm,3)
c
      integer i
c
      call timer_start(7)
c
      do 1 i=1,nsm
      xsmp(i,1)=xsm(i,1)
      xsmp(i,2)=xsm(i,2)
      xsmp(i,3)=xsm(i,3)
    1 continue
c
c     initialize cafe
c
      call argos_cafe_init(stime,lpair,llong,box,vlat,vlati,zw,
     + dbl_mb(i_zs),
     + eww,dbl_mb(i_esw),dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esa))
c
      if(lpair) call argos_zinit(int_mb(i_iwz),int_mb(i_isz))
c
      call timer_stop(7)
c
      call timer_start(9)
      call ga_sync()
      call timer_stop(9)
c
      call timer_start(10)
c
      call argos_space_initf(fw,fs,llong,int_mb(i_iwz),int_mb(i_isz),
     + lpair)
c
      call argos_space_putix(me,iwl,xw,nwmloc,isl,xs,nsaloc)
c
      call timer_stop(10)
c
      call timer_start(11)
      call ga_sync()
      call timer_stop(11)
c
      if(ncoll.gt.0) then
      call timer_start(12)
      call argos_cafe_collapse(ncoll,fcoll,nsaloc,nwmloc,isl(1,lsmol),
     + isl(1,lssgm),
     + dbl_mb(i_xs),dbl_mb(i_xsm),mst,dbl_mb(i_tsm),dbl_mb(i_fs),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_fw))
      call timer_stop(12)
      endif
      if(ifield.gt.0) then
      call timer_start(13)
      call argos_cafe_extern(stime,nsaloc,dbl_mb(i_fs),
     + isl(1,lsct1),nwmloc,dbl_mb(i_fw))
      call timer_stop(13)
      endif
      call argos_cafe_multi(nsaloc,dbl_mb(i_xs),dbl_mb(i_fs),isl(1,
     + lsgan),
     + isl(1,lsatt),isl(1,lsfrc),isl(1,lsdyn),isl(1,lsct1),
     + dbl_mb(i_ess),dbl_mb(i_fss),lfnpmf,lprpmf,iprpmf,
     + npmf.eq.1.or.lpmfc)
c
      if(ipolt.gt.0) then
      call timer_start(23)
      call argos_induce_start(iwl,isl,xw,xwm,xs,dbl_mb(i_pw),
     + dbl_mb(i_pwp),dbl_mb(i_ps),dbl_mb(i_psp))
      call timer_stop(23)
      endif
c
      if(lpme.and.llong) then
      if(lpert2) then
      call argos_pme_energy(2,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(2))
      call timer_start(24)
      call argos_pme_init()
      call timer_stop(24)
      endif
      if(lpert3) then
      call argos_pme_energy(3,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(3))
      call timer_start(24)
      call argos_pme_init()
      call timer_stop(24)
      endif
      call argos_pme_chgrid(iset,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(iset))
      endif
c
      return
      end
c $Id$
