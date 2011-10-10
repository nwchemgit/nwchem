      subroutine argos_congra(iwdt,xw,yw,vw,fw,pcgw,
     + isdt,xs,ys,vs,fs,pcgs,ww,ws)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      logical frequency
      external frequency
c
      integer iwdt(mwm),isdt(msa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),fs(msa,3)
      real*8 pcgw(mwm,3,mwa),pcgs(msa,3),ww(mwa),ws(msa)
c
      integer iwa,iwm,isa,ix,inner
      logical ldone
      real*8 alpha,beta1,beta2,beta3,beta4,beta5,gamma,zeta
      real*8 dx,ecgnew,eqrs,ecgold,ecgdif,edt,dxf
      real*8 dxsmax,dxwmax,fnorm1,fnorm2,ypa,ypb,pnorm,dxstep
      real*8 dxmax,dxfi,ecg0,ecg1,ecg2,epcgw,epcgsw,epcgs
      character*1 cqrs
c
      if(me.eq.0) write(lfnout,1000)
 1000 format(/,' CONJUGATE GRADIENT MINIMIZATION',//,
     + '   Step File     Energy       Energy       Energy   ',
     + '    Energy       Energy     Largest  ',/,
     + '        wrt     gradient       Total      solvent   ',
     + '   slv-sol       solute  displacement',/,
     + '                 kJ/mol       kJ/mol       kJ/mol   ',
     + '    kJ/mol       kJ/mol        nm',/)
c
      dx=dx0cg
      beta1=zero
      icgit=0
      lpair=.true.
      lload=.true.
      lhop=.false.
c
      call timer_start(201)
c
c     atomic forces and potential energies
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call argos_prop_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      ecgnew=epot
      eqrs=epot
      ecgold=ecgnew
      ecgdif=zero
      edt=zero
c
c     copy initial coordinates into vw and vs
c     initialize search direction vectors pcgw and pcgs to zero
c
      if(nwmloc.gt.0) then
      do 1 iwa=1,nwa
      do 2 ix=1,3
      do 3 iwm=1,nwmloc
      vw(iwm,ix,iwa)=xw(iwm,ix,iwa)
      pcgw(iwm,ix,iwa)=zero
    3 continue
    2 continue
    1 continue
      endif
      if(nsaloc.gt.0) then
      do 4 ix=1,3
      do 5 isa=1,nsaloc
      vs(isa,ix)=xs(isa,ix)
      pcgs(isa,ix)=zero
    5 continue
    4 continue
      endif
c
c     take one steepest descent step to get initial search direction
c
      dxf=half*dx/fmax
      if(nwmloc.gt.0) then
      do 6 iwa=1,nwa
      do 7 ix=1,3
      do 8 iwm=1,nwmloc
      yw(iwm,ix,iwa)=xw(iwm,ix,iwa)
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      xw(iwm,ix,iwa)=xw(iwm,ix,iwa)+dxf*fw(iwm,ix,iwa)/ww(iwa)
      endif
    8 continue
    7 continue
    6 continue
      endif
      if(nsaloc.gt.0) then
      do 9 ix=1,3
      do 10 isa=1,nsaloc
      ys(isa,ix)=xs(isa,ix)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      xs(isa,ix)=xs(isa,ix)+dxf*fs(isa,ix)/ws(isa)
      endif
   10 continue
    9 continue
      endif
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      fnorm1=zero
      dxfi=one/dxf
      if(nwmloc.gt.0) then
      do 11 iwa=1,nwa
      do 12 ix=1,3
      do 13 iwm=1,nwmloc
      fw(iwm,ix,iwa)=(xw(iwm,ix,iwa)-yw(iwm,ix,iwa))*dxfi*ww(iwa)
      fnorm1=fnorm1+fw(iwm,ix,iwa)**2
   13 continue
   12 continue
   11 continue
      endif
      if(nsaloc.gt.0) then
      do 14 ix=1,3
      do 15 isa=1,nsaloc
      fs(isa,ix)=(xs(isa,ix)-ys(isa,ix))*dxfi*ws(isa)
      fnorm1=fnorm1+fs(isa,ix)**2
   15 continue
   14 continue
      endif
c
c     global sum fnorm1
c
      call ga_dgop(mrg_d08,fnorm1,1,'+')
c
      ecg0=ecgnew
      ecg1=ecgnew
      ecg2=ecgnew
      icgit=0
c
      call timer_stop(201)
c
c     outer loop
c
  100 continue
c
      if(icgit.eq.(icgit/ncgcy)*ncgcy) beta1=zero
      icgit=icgit+1
      lpair=frequency(icgit,nfpair)
      lload=frequency(icgit,nfload)
      lhop=.false.
c
      ypa=zero
      pnorm=zero
      if(nwmloc.gt.0) then
      do 16 iwa=1,nwa
      do 17 ix=1,3
      do 18 iwm=1,nwmloc
      pcgw(iwm,ix,iwa)=fw(iwm,ix,iwa)+beta1*pcgw(iwm,ix,iwa)
      pnorm=pnorm+pcgw(iwm,ix,iwa)**2
      ypa=ypa+pcgw(iwm,ix,iwa)*fw(iwm,ix,iwa)
   18 continue
   17 continue
   16 continue
      endif
      if(nsaloc.gt.0) then
      do 19 ix=1,3
      do 20 isa=1,nsaloc
      pcgs(isa,ix)=fs(isa,ix)+beta1*pcgs(isa,ix)
      pnorm=pnorm+pcgs(isa,ix)**2
      ypa=ypa+pcgs(isa,ix)*fs(isa,ix)
   20 continue
   19 continue
      endif
c
c     accumulate pnorm
c
      call ga_dgop(mrg_d09,ypa,1,'+')
      call ga_dgop(mrg_d10,pnorm,1,'+')
c
      if(pnorm.lt.zero) call md_abort('congra: pnorm<zero',0)
      pnorm=sqrt(pnorm)
c
      alpha=zero
      ecg1=ecg2
      beta2=dx/pnorm
      inner=0
c
c     inner loop
c
  200 continue
      call timer_start(201)
c
      inner=inner+1
c
      if(nwmloc.gt.0) then
      do 21 iwa=1,nwa
      do 22 ix=1,3
      do 23 iwm=1,nwmloc
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      xw(iwm,ix,iwa)=vw(iwm,ix,iwa)+beta2*pcgw(iwm,ix,iwa)/ww(iwa)
      endif
   23 continue
   22 continue
   21 continue
      endif
      if(nsaloc.gt.0) then
      do 24 ix=1,3
      do 25 isa=1,nsaloc
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      xs(isa,ix)=vs(isa,ix)+beta2*pcgs(isa,ix)/ws(isa)
      endif
   25 continue
   24 continue
      endif
c
c     atomic forces and potential energies
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call argos_em_fnorm(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      lpair=.false.
      lload=.false.
      lhop=.false.
c
      beta3=beta2*pnorm/fnorm
c
      if(nwmloc.gt.0) then
      do 26 iwa=1,nwa
      do 27 ix=1,3
      do 28 iwm=1,nwmloc
      yw(iwm,ix,iwa)=xw(iwm,ix,iwa)
      dxstep=beta3*fw(iwm,ix,iwa)/ww(iwa)
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      xw(iwm,ix,iwa)=yw(iwm,ix,iwa)+dxstep
      endif
   28 continue
   27 continue
   26 continue
      endif
      if(nsaloc.gt.0) then
      do 29 ix=1,3
      do 30 isa=1,nsaloc
      ys(isa,ix)=xs(isa,ix)
      dxstep=beta3*fs(isa,ix)/ws(isa)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      xs(isa,ix)=ys(isa,ix)+dxstep
      endif
   30 continue
   29 continue
      endif
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
c     find constrained forces
c
      ypb=zero
      if(nwmloc.gt.0) then
      do 34 iwa=1,nwa
      do 35 ix=1,3
      do 36 iwm=1,nwmloc
      fw(iwm,ix,iwa)=(xw(iwm,ix,iwa)-yw(iwm,ix,iwa))*ww(iwa)/beta3
      ypb=ypb+pcgw(iwm,ix,iwa)*fw(iwm,ix,iwa)
   36 continue
   35 continue
   34 continue
      endif
      if(nsaloc.gt.0) then
      do 37 ix=1,3
      do 38 isa=1,nsaloc
      fs(isa,ix)=(xs(isa,ix)-ys(isa,ix))*ws(isa)/beta3
      ypb=ypb+pcgs(isa,ix)*fs(isa,ix)
   38 continue
   37 continue
      endif
      call ga_dgop(mrg_d11,ypb,1,'+')
c
      call argos_prop_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      ecg2=epot
      epcgw=epotw
      epcgsw=epotsw
      epcgs=epots
      ecgdif=ecg2-ecgold
      edt=ecg2-ecg0
      call timer_stop(201)
c
c     check if interval is appropriate
c
      if(ypb.ge.zero.and.ecg2.lt.ecg1) then
      alpha=beta2
      ecg1=ecg2
      ypa=ypb
      beta2=two*beta2
      goto 200
      endif
      call timer_start(201)
c
c     find interpolation in interval
c
      zeta=three*(ecg1-ecg2)/(beta2-alpha)-ypa-ypb
      gamma=zeta**2-ypa*ypb
c
      if(gamma.lt.zero) then
      gamma=zero
      else
      gamma=sqrt(gamma)
      endif
c
      beta4=beta2-(gamma-zeta-ypb)*(beta2-alpha)/(ypa-ypb+two*gamma)
c
c     advance coordinates to interpolated point
c
      dxmax=zero
      if(nwmloc.gt.0) then
      do 39 iwa=1,nwa
      do 40 ix=1,3
      do 41 iwm=1,nwmloc
      yw(iwm,ix,iwa)=vw(iwm,ix,iwa)
      dxstep=beta4*pcgw(iwm,ix,iwa)/ww(iwa)
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      xw(iwm,ix,iwa)=vw(iwm,ix,iwa)+dxstep
      if(abs(dxstep).gt.dxmax) dxmax=abs(dxstep)
      endif
   41 continue
   40 continue
   39 continue
      endif
      if(nsaloc.gt.0) then
      do 42 ix=1,3
      do 43 isa=1,nsaloc
      ys(isa,ix)=vs(isa,ix)
      dxstep=beta4*pcgs(isa,ix)/ws(isa)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      xs(isa,ix)=vs(isa,ix)+dxstep
      if(abs(dxstep).gt.dxmax) dxmax=abs(dxstep)
      endif
   43 continue
   42 continue
      endif
      call ga_dgop(mrg_d12,dxmax,1,'max')
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
c     atomic forces and potential energies
c
      call argos_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call argos_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call argos_prop_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      dxmax=zero
      dxwmax=zero
      dxsmax=zero
c
      beta5=beta4*pnorm/fnorm
c
c     advance coordinates with these forces
c
      if(nwmloc.gt.0) then
      do 44 iwa=1,nwa
      do 45 ix=1,3
      do 46 iwm=1,nwmloc
      yw(iwm,ix,iwa)=xw(iwm,ix,iwa)
      dxstep=beta5*fw(iwm,ix,iwa)/ww(iwa)
      if(iand(iwdt(iwm),mfixed).ne.lfixed) then
      xw(iwm,ix,iwa)=yw(iwm,ix,iwa)+dxstep
      if(abs(dxstep).gt.dxwmax) dxwmax=abs(dxstep)
      endif
   46 continue
   45 continue
   44 continue
      if(dxwmax.gt.dxmax) dxmax=dxwmax
      endif
      if(nsaloc.gt.0) then
      do 47 ix=1,3
      do 48 isa=1,nsaloc
      ys(isa,ix)=xs(isa,ix)
      dxstep=beta5*fs(isa,ix)/ws(isa)
      if(iand(isdt(isa),mfixed).ne.lfixed) then
      xs(isa,ix)=ys(isa,ix)+dxstep
      if(abs(dxstep).gt.dxsmax) dxsmax=abs(dxstep)
      endif
   48 continue
   47 continue
      if(dxsmax.gt.dxmax) dxmax=dxsmax
      endif
      call ga_dgop(mrg_d12,dxmax,1,'max')
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      fnorm2=zero
c
      if(nwmloc.gt.0) then
      do 52 iwa=1,nwa
      do 53 ix=1,3
      do 54 iwm=1,nwmloc
      fw(iwm,ix,iwa)=(xw(iwm,ix,iwa)-yw(iwm,ix,iwa))*ww(iwa)/beta5
      pcgw(iwm,ix,iwa)=(yw(iwm,ix,iwa)-vw(iwm,ix,iwa))*ww(iwa)/beta4
      vw(iwm,ix,iwa)=yw(iwm,ix,iwa)
      fnorm2=fnorm2+fw(iwm,ix,iwa)**2
   54 continue
   53 continue
   52 continue
      endif
      if(nsaloc.gt.0) then
      do 55 ix=1,3
      do 56 isa=1,nsaloc
      fs(isa,ix)=(xs(isa,ix)-ys(isa,ix))*ws(isa)/beta5
      pcgs(isa,ix)=(ys(isa,ix)-vs(isa,ix))*ws(isa)/beta4
      vs(isa,ix)=ys(isa,ix)
      fnorm2=fnorm2+fs(isa,ix)**2
   56 continue
   55 continue
      endif
c
c     global sum fnorm2
c
      call ga_dgop(mrg_d13,fnorm2,1,'+')
c
      beta1=sqrt(fnorm2/fnorm1)
      fnorm1=fnorm2
c
      ecg2=epot
      epcgw=epotw
      epcgsw=epotsw
      epcgs=epots
c
      ecgdif=ecg2-ecgold
      ecgold=ecg2
      edt=ecg2-ecg0
      ecg1=ecg2
      ecgnew=ecg2
c
c     record to lfntrj
c
      lxw=frequency(mdstep,nfcoor)
      lxs=frequency(mdstep,nfscoo)
c
      if(lxw.or.lxs) then
      call argos_space_wrttrj(lfntrj,lxw,.false.,.false.,lxs,.false.,
     + .false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
c
c     write restart file
c
      ldone=.not.(icgit.lt.mcgit.and.dxmax.gt.dxcgmx.and.ecgdif.lt.zero)
c
      cqrs=' '
      if(frequency(icgit,nfqrs).or.(ldone.and.eqrs.gt.epot)) then
      write(projct,4000) nserie,msdit,icgit,filnam(1:56)
 4000 format(i2,' em ',i7,' + ',i7,' ',a)
      call argos_wrtrst(lfnqrs,filqrs,.false.)
      cqrs='*'
      eqrs=ecg1
      endif
      if(me.eq.0.and.frequency(icgit,nfprop)) call argos_prop_record()
c
c     print minimization step data to output
c
      if(me.eq.0) then
      write(lfnout,600) icgit,cqrs,ecgdif,ecg1,epcgw,epcgsw,epcgs,dxmax
  600 format(i7,1x,a1,3x,5(1pe13.5),0pf12.8)
      endif
c
      if(itest.gt.0) call argos_test()
c
      call timer_stop(201)
      if(.not.ldone) goto 100
c
      ecgdif=edt
c
      return
      end
c $Id$
