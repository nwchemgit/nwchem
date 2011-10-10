      subroutine argos_md()
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical frequency
      real*8 timer_wall_total
      external frequency,timer_wall_total
c
      character*1 mdt
      logical ltmp
c
      lfirst=.true.
      lxw=.false.
      lvw=.false.
      lxs=.false.
      lvs=.false.
      lesp=.false.
c
c     equilibration
c
      ltmp=lhop
      lhop=.false.
      lequi=.true.
      lprpmf=.false.
      do 1 iequi=kequi+1,mequi
c
      mdstep=mdstep+1
      stime=stime+tstep
      lpmfc=npmf.gt.1.and.iequi.gt.npmf
      call timer_start(201)
      call argos_newton()
      call argos_prop_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call argos_server
      call timer_stop(201)
    1 continue
      lpmfc=.true.
c
c     data gathering
c
      lhop=ltmp
      lequi=.false.
      lprpmf=iprpmf.ne.0
      if(lprpmf) lfnpmf=-iabs(lfnpmf)
      call timer_start(205)
      do 2 idacq=kdacq+1,mdacq
c
      mdstep=mdstep+1
      stime=stime+tstep
c
      lxw=frequency(mdstep,nfcoor)
      lvw=frequency(mdstep,nfvelo)
      lfw=frequency(mdstep,nfforc)
      lxs=frequency(mdstep,nfscoo)
      lvs=frequency(mdstep,nfsvel)
      lfs=frequency(mdstep,nfsfor)
      lesp=frequency(mdstep,nfesp)
c
      call argos_timer_init()
c
      call timer_start(201)
c
      call argos_newton()
c
      call timer_start(6)
c
      if(lfw.or.lfs) then
      call argos_space_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),
     + nsaloc)
      endif
c
      call timer_stop(6)
c
      call timer_start(55)
      mdt=' '
      if(iguide.gt.0) mdt='g'
      write(projct,4000) nserie,mdt,mequi,idacq,tmpext,prsext,
     + filnam(1:38)
 4000 format(i2,' md',a1,i7,' + ',i7,' @',f7.2,e9.2,' ',a)
      if(frequency(mdstep,nfrest)) then
      call argos_wrtrst(lfnrst,rfile,.true.)
      endif
      if(frequency(mdstep,nftime)) call argos_wrtime
c
      call timer_stop(55)
c
      call argos_prop_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call argos_server
      call argos_prop_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme)
c
      call timer_stop(201)
c
      if(frequency(mdstep,nfnewf).and.idacq.ne.mdacq) then
      if(me.eq.0) call argos_fopen(.true.)
      endif
c
c
      if(lstop) then
      if(me.eq.0) write(*,1000) tleft,tneed
 1000 format(///,' Time left (',f12.3,' s) is less than twice the ',
     + ' time needed to reach writing the next restart file (',
     + f12.3,' s)',//,' Simulation aborted',//)
      return
      endif
c
    2 continue
c
      return
      end
c $Id$
