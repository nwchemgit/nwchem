      subroutine argos_fclass(iwl,isl,xw,xwm,xs,fw,fs)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical argos_space_local,argos_cafe_hopping
      external argos_space_local,argos_cafe_hopping
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      logical lself,local,lpbcs
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer nbbl
      integer itcom1,itforc,itcom2
c
      logical lforce,ldone,lnew
c
      local=.true.
      itcom1=33
      itforc=34
      itcom2=35
c
      if(nbget.eq.0) then
      call argos_space_nbbl(nbbl)
      else
      itforc=35
      call timer_start(33)
      if(ipolt.eq.0) then
      if(nbget.lt.0) then
      call argos_space_prefetch_all(nbbl,iwl,xw,isl,xs)
      else
      call argos_space_prefetch(nbbl,iwl,xw,isl,xs)
      endif
      else
      call argos_space_prefetch_p(nbbl,iwl,xw,dbl_mb(i_pw),
     + dbl_mb(i_pwp),
     + isl,xs,dbl_mb(i_ps),dbl_mb(i_psp))
      endif
      call timer_stop(33)
      endif
c
      lforce=.true.
      ldone=.true.
c
      if(lpair.and.lhop) then
      lforce=.false.
      ldone=.false.
      endif
c
  100 continue
c
      do 1 ibbl=1,nbbl
c
      if(local.and..not.argos_space_local(ibbl)) then
      itcom1=36
      itforc=37
      itcom2=38
      local=.false.
      if(nbget.ne.0) itforc=36
      endif
c
      if(nbget.eq.0) then
      call timer_start(itcom1)
      if(ipolt.eq.0) then
      call argos_space_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lself,lpbcs)
      else
      call argos_space_getxpbl(ibbl,lhandl,
     + iwl,xw,dbl_mb(i_pw),dbl_mb(i_pwp),iwfr,iwto,jwfr,jwto,
     + isl,xs,dbl_mb(i_ps),dbl_mb(i_psp),isfr,isto,jsfr,jsto,
     + lself,lpbcs)
      endif
      call timer_stop(itcom1)
      else
      call timer_start(34)
      call argos_space_nbwait(ibbl,lnew,lhandl,lself,lpbcs,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,isl)
      call timer_stop(34)
      if(nbget.gt.0.and.lnew) then
      call timer_start(33)
      call argos_space_prefetch_next(iwl,xw,isl,xs)
      call timer_stop(33)
      endif
      endif
c
      call timer_start(itforc)
      call argos_cafe_comw(xw,xwm,jwfr,jwto)
c
      if(ipolt.eq.0) then
      call argos_cafe(lself,lpbcs,xw,xwm,fw,zw,dbl_mb(i_rtos),
     + iwl(1,lwdyn),
     + int_mb(i_iwz),iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),fs,
     + dbl_mb(i_zs),isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),
     + isl(1,lsfrc),isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),
     + isl(1,lsct3),isl(1,lssgm),isl(1,lshop),int_mb(i_isz),
     + isfr,isto,jsfr,jsto,lpbc,lhandl,
     + .true.,eww,dbl_mb(i_esw),dbl_mb(i_ess),dbl_mb(i_fss),
     + dbl_mb(i_esa),int_mb(i_lseq),lforce)
      else
      call argos_cafe_polar(lself,lpbcs,xw,xwm,fw,dbl_mb(i_pw),
     + dbl_mb(i_pwp),zw,
     + dbl_mb(i_rtos),iwl(1,lwdyn),int_mb(i_iwz),iwfr,iwto,jwfr,jwto,xs,
     + dbl_mb(i_xsm),fs,dbl_mb(i_ps),dbl_mb(i_psp),dbl_mb(i_zs),
     + isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),isl(1,lsfrc),
     + isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),isl(1,lsct3),
     + isl(1,lssgm),isl(1,lshop),int_mb(i_isz),
     + isfr,isto,jsfr,jsto,lpbc,lhandl,
     + .true.,eww,dbl_mb(i_esw),dbl_mb(i_ess),dbl_mb(i_fss),
     + dbl_mb(i_esa),int_mb(i_lseq))
      endif
      call timer_stop(itforc)
c
      if(nbget.eq.0) then
      call timer_start(itcom2)
      call argos_space_accfbl(ibbl,lhandl,fw,fs,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(itcom2)
      else
      call timer_start(37)
      call argos_space_nbaccfbl(ibbl,lhandl,fw,fs,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(37)
      endif
c
    1 continue
c
      if(nbget.ne.0) then
      call timer_start(38)
      call argos_space_nbwaitf()
      call timer_stop(38)
      endif
c
      if(.not.ldone) then
      if(argos_cafe_hopping(lpbc,lpbcs,stime,isl,isl(1,lssgm),isl(1,
     + lsgan),
     + isl(1,lsct3),isl(1,lshop),xs,nsaloc,lfnhop)) then
      if(lpair) then
      lpair=.false.
      lload=.false.
      lforce=.true.
      ldone=.false.
      else
      lpair=.true.
      lload=.true.
      lforce=.true.
      ldone=.true.
      call argos_space_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
      endif
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      goto 100
      endif
      endif
c
      if(lpme.and.llong) then
      call argos_pme_forces(fw(1,1,1,2),nwmloc,fs(1,1,2),
     + isl(1,lsct1),isl(1,lssgm),nsaloc)
      endif
c
      return
      end
c $Id$
