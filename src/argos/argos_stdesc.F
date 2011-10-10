      subroutine argos_stdesc(iwdt,xw,yw,vw,fw,isdt,ismol,
     + xs,ys,vs,fs,ww,ws,mm,fm,xsm)
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical frequency
      external frequency
c
      integer iwdt(mwm),isdt(msa),ismol(msa),mm(msa,2)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),fs(msa,3)
      real*8 ww(mwa),ws(msa),xsm(msm,3),fm(msm,7)
c
      integer i,j
      logical ldone
      real*8 edif,epsd,epsdw,epsdsw,epsds,ecrit
      real*8 da,damsd,dx,dxf,dxmax,factor,dxstep,eqrs
      character*1 cqrs
      real*8 angle,o(3),p(3),x(3),y(3)
c
      double precision grms(2)
      integer nrms(2)
c
      damsd=-1.1d0
      if(imembr.gt.0) call argos_membrane_init(ismol,mm,xs,xsm,fm)
c
      if(me.eq.0) write(lfnout,1000)
 1000 format(/,' STEEPEST DESCENT MINIMIZATION',//,
     + '   Step File     Energy       Energy       Energy   ',
     + '    Energy       Energy     Largest  ',/,
     + '        wrt     gradient       Total      solvent   ',
     + '   slv-sol       solute  displacement',/,
     + '                 kJ/mol       kJ/mol       kJ/mol   ',
     + '    kJ/mol       kJ/mol        nm',/)
c
      isdit=0
c
      lpair=.true.
      lload=.true.
      lhop=.false.
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
      epsd=epot
      epsdw=epotw
      epsdsw=epotsw
      epsds=epots
      eqrs=epot
      dx=dx0sd
      mdstep=0
      da=0.01
c
    1 continue
      call timer_start(201)
c
      dxf=dx/fmax
c
      if(nwmloc.gt.0) then
      do 2 j=1,nwa
      do 3 i=1,nwmloc
      yw(i,1,j)=xw(i,1,j)
      yw(i,2,j)=xw(i,2,j)
      yw(i,3,j)=xw(i,3,j)
      vw(i,1,j)=fw(i,1,j)
      vw(i,2,j)=fw(i,2,j)
      vw(i,3,j)=fw(i,3,j)
    3 continue
      factor=one/ww(j)
      do 4 i=1,nwmloc
      if(iand(iwdt(i),mfixed).ne.lfixed) then
      dxstep=factor*dxf*fw(i,1,j)
      xw(i,1,j)=xw(i,1,j)+dxstep
      dxstep=factor*dxf*fw(i,2,j)
      xw(i,2,j)=xw(i,2,j)+dxstep
      dxstep=factor*dxf*fw(i,3,j)
      xw(i,3,j)=xw(i,3,j)+dxstep
      endif
    4 continue
    2 continue
      endif
c
      if(nsaloc.gt.0) then
      do 5 i=1,nsaloc
      ys(i,1)=xs(i,1)
      ys(i,2)=xs(i,2)
      ys(i,3)=xs(i,3)
      vs(i,1)=fs(i,1)
      vs(i,2)=fs(i,2)
      vs(i,3)=fs(i,3)
    5 continue
      if(imembr.eq.0) then
      do 6 i=1,nsaloc
      if(iand(isdt(i),mfixed).ne.lfixed) then
      factor=one/ws(i)
      dxstep=factor*dxf*fs(i,1)
      xs(i,1)=xs(i,1)+dxstep
      dxstep=factor*dxf*fs(i,2)
      xs(i,2)=xs(i,2)+dxstep
      dxstep=factor*dxf*fs(i,3)
      xs(i,3)=xs(i,3)+dxstep
      endif
    6 continue
      else
      do 16 i=1,msm
      fm(i,1)=zero
      fm(i,2)=zero
      fm(i,3)=zero
      fm(i,4)=zero
      fm(i,5)=zero
      fm(i,6)=zero
      fm(i,7)=zero
   16 continue
      do 17 i=1,nsaloc
      factor=one/ws(i)
      fm(mm(i,2),1)=fm(mm(i,2),1)+factor*fs(i,1)
      fm(mm(i,2),2)=fm(mm(i,2),2)+factor*fs(i,2)
      fm(mm(i,2),3)=fm(mm(i,2),3)+factor*
     + ((xs(i,1)-xsm(mm(i,2),1))*fs(i,2)-
     +  (xs(i,2)-xsm(mm(i,2),2))*fs(i,1))
   17 continue
      if(np.gt.1) call ga_dgop(mrg_d50,fm,3*msm,'+')
      do 18 i=1,nsm
      fm(i,4)=fm(i,3)
      write(*,'(a,i5,a,i5)') 'Mol ',i,', Angle ',da*fm(i,4)
   18 continue
      do 19 i=1,nsaloc
      if(iand(isdt(i),mfixed).ne.lfixed) then
c
c     rotations
c
      o(1)=xsm(mm(i,2),1)
      o(2)=xsm(mm(i,2),1)
      o(3)=xsm(mm(i,2),1)
      p(1)=o(1)
      p(2)=o(2)
      p(3)=o(3)+1.0d0
      x(1)=xs(i,1)
      x(2)=xs(i,2)
      x(3)=xs(i,3)
      angle=da*fm(mm(i,2),4)
      call rotate(o,p,angle,x,y)
      xs(i,1)=y(1)
      xs(i,2)=y(2)
      xs(i,3)=y(3)
c
c     translations
c
      dxstep=dxf*fm(mm(i,2),1)
      xs(i,1)=xs(i,1)+dxstep
      dxstep=dxf*fm(mm(i,2),2)
      xs(i,2)=xs(i,2)+dxstep
c
      endif
   19 continue
      endif
      endif
c
c     shake
c
      call argos_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      isdit=isdit+1
c
      lpair=frequency(isdit,nfpair)
      lload=lpair
      lhop=.false.
      llong=(frequency(isdit,nflong).or.lpair).and.ltwin
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
c     if energy goes up restore coordinates and forces
c
      ecrit=epot
      if(icrit.eq.1) ecrit=epot
c
c     changing the logic here so that bad coordinates
c     are not preserved in the last md iteration (M.V.)
c     ------------------------------------------------
      if(ecrit.gt.epsd) then
c
      call argos_cafe_restor(dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_fw),
     + dbl_mb(i_vw),nwmloc,dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_fs),
     + dbl_mb(i_vs),nsaloc)
c
      if(isdit.lt.msdit.and.dxmax.gt.dxsdmx) then
      dx=half*dx
      da=half*da
      call timer_stop(201)
      if(dx.gt.dxsdmx) goto 1
      else
c
c       this insures the global restoration of coordinates
c       so that the right restart file written out
c       --------------------------------------------------
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
      endif
      endif
c
      edif=ecrit-epsd
      epsd=ecrit
      epsdw=epotw
      epsdsw=epotsw
      epsds=epots
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
      ldone=.not.(dxmax.gt.dxsdmx.and.isdit.lt.msdit)
c
      cqrs=' '
      if(frequency(isdit,nfqrs).or.(ldone.and.eqrs.gt.epot)) then
      write(projct,4000) nserie,isdit,0,filnam(1:56)
 4000 format(i2,' em ',i7,' + ',i7,' ',a)
      call argos_wrtrst(lfnqrs,filqrs,.false.)
      cqrs='*'
      eqrs=epot
      endif
      if(me.eq.0.and.frequency(isdit,nfprop)) call argos_prop_record()
c
      if(me.eq.0) then
      write(lfnout,600) isdit,cqrs,edif,epsd,epsdw,epsdsw,epsds,dxmax
  600 format(i7,1x,a1,3x,5(1pe13.5),0pf12.8)
      call util_flush(lfnout)
      endif
c
      if(itest.gt.0) call argos_test()
c
      dx=min(1.2d0*dx,dxmsd)
      da=min(1.2d0*da,damsd)
c
      call timer_stop(201)
      if(.not.ldone) goto 1
c
      return
      end
c $Id$
