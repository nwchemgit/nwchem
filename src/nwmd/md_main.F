      subroutine md_main()
c
c $Id: md_main.F,v 1.71 2000-04-13 20:22:54 d3j191 Exp $
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
c
c     energy minimization
c     -------------------
c
      call md_em()
c
      elseif(ntype.eq.2) then
c
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
      integer i_wws,l_wws
c
      lpair=.true.
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
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
      if(.not.ma_push_get(mt_dbl,mwa+msa,'wws',l_wws,i_wws))
     + call errquit('Failed to allocate memory for wws',me)
c
      call cf_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call md_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
c
      call prp_print()
c
      call rtdb_put(irtdb,'md:energy',mt_dbl,1,epot)
c
      if(.not.ma_pop_stack(l_wws))
     + call errquit('Failed to deallocate memory for wws',me)
c
c     print energies
c
      call cf_print_energy(lfnout)
c
      call sp_printf(lfnout,filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs))
c
      if(itest.gt.0) call md_test()
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
      integer i_pcgw,l_pcgw,i_pcgs,l_pcgs,i_wws,l_wws
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
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
      if(.not.ma_push_get(mt_dbl,mwa+msa,'wws',l_wws,i_wws))
     + call errquit('Failed to allocate memory for wws',me)
c
      call cf_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call md_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
      if(msdit.gt.0) call md_stdesc(int_mb(i_iw+(lwdyn-1)*mwm),
     + dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),
     + int_mb(i_is+(lsdyn-1)*msa),dbl_mb(i_xs),dbl_mb(i_ys),
     + dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_wws),
     + dbl_mb(i_wws+mwa))
c
      if(mcgit.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'pcgw',l_pcgw,i_pcgw))
     + call errquit('Failed to allocate memory for pcgw',me)
      if(.not.ma_push_get(mt_dbl,3*msa,'pcgs',l_pcgs,i_pcgs))
     + call errquit('Failed to allocate memory for pcgs',me)
      call md_congra(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),dbl_mb(i_pcgw),
     + int_mb(i_is+(lsdyn-1)*msa),dbl_mb(i_xs),dbl_mb(i_ys),
     + dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_pcgs),dbl_mb(i_wws),
     + dbl_mb(i_wws+mwa))
      if(.not.ma_pop_stack(l_pcgs))
     + call errquit('Failed to deallocate memory for pcgs',me)
      if(.not.ma_pop_stack(l_pcgw))
     + call errquit('Failed to deallocate memory for pcgw',me)
      endif
c
      if(.not.ma_pop_stack(l_wws))
     + call errquit('Failed to deallocate memory for wws',me)
c
      call cf_print_energy(lfnout)
c
      call sp_printf(lfnout,filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs))
c
      return
      end
      subroutine md_eminit(xw,yw,xs,ys)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3)
c
      integer i,j
      real*8 dxmax
c
      if(nwmloc.gt.0) then
      do 1 j=1,nwa
      do 2 i=1,nwmloc
      yw(i,1,j)=xw(i,1,j)
      yw(i,2,j)=xw(i,2,j)
      yw(i,3,j)=xw(i,3,j)
    2 continue
    1 continue
      endif
c
      if(nsaloc.gt.0) then
      do 3 i=1,nsaloc
      ys(i,1)=xs(i,1)
      ys(i,2)=xs(i,2)
      ys(i,3)=xs(i,3)
    3 continue
      endif
c
      lpair=.true.
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     shake
c
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      return
      end
      subroutine md_stdesc(iwdt,xw,yw,vw,fw,isdt,xs,ys,vs,fs,ww,ws)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      logical frequency
      external frequency
c
      integer iwdt(mwm),isdt(msa)
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),vw(mwm,3,mwa),fw(mwm,3,mwa)
      real*8 xs(msa,3),ys(msa,3),vs(msa,3),fs(msa,3)
      real*8 ww(mwa),ws(msa)
c
      integer i,j
      logical ltr,ldone
      real*8 edif,epsd,epsdw,epsdsw,epsds
      real*8 dx,dxf,dxmax,factor,dxstep,eqrs
      character*1 cqrs
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
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
c
      epsd=epot
      epsdw=epotw
      epsdsw=epotsw
      epsds=epots
      eqrs=epot
      dx=dx0sd
      mdstep=0
c
    1 continue
      call timer_start(1)
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
      endif
c
c      call cf_sdstep(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
c     + dbl_mb(i_yw),dbl_mb(i_fw),dbl_mb(i_vw),nwmloc,
c     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
c     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_fs),dbl_mb(i_vs),nsaloc,dxf)
c
c     shake
c
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
      isdit=isdit+1
c
      lpair=frequency(isdit,nfpair)
      llong=(frequency(isdit,nflong).or.lpair).and.ltwin
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
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
      call timer_stop(1)
      goto 1
      endif
c
      edif=epot-epsd
      epsd=epot
      epsdw=epotw
      epsdsw=epotsw
      epsds=epots
c
      ltr=frequency(isdit,nfem)
c
      if(ltr) then
      if(lqmd) then
      call qmd_wrttrj(lfnemt,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfnemt,.true.,.false.,.true.,.false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
      endif
c
      ldone=.not.(dxmax.gt.dxsdmx.and.isdit.lt.msdit)
c
      cqrs=' '
      if(frequency(isdit,nfqrs).or.(ldone.and.eqrs.gt.epot)) then
      write(projct,4000) isdit,0,filnam(1:55)
 4000 format('em sd ',i7,' cg ',i7,' ',a)
      call md_wrtrst(lfnqrs,filqrs,.false.)
      cqrs='*'
      eqrs=epot
      endif
      if(me.eq.0.and.frequency(isdit,nfprop)) call prp_record()
c
      if(me.eq.0) then
      write(lfnout,600) isdit,cqrs,edif,epsd,epsdw,epsdsw,epsds,dxmax
  600 format(i7,1x,a1,3x,5(1pe13.5),0pf12.8)
      endif
c
      if(itest.gt.0) call md_test()
c
      dx=min(1.2d0*dx,dxmsd)
c
      call timer_stop(1)
      if(.not.ldone) goto 1

      return
      end
      subroutine md_congra(iwdt,xw,yw,vw,fw,pcgw,
     + isdt,xs,ys,vs,fs,pcgs,ww,ws)
c
      implicit none
c
#include "md_common.fh"
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
      logical ltr,ldone
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
c
      call timer_start(1)
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
c
      ecgnew=epot
      eqrs=epot
      ecgold=ecgnew
      ecgdif=zero
      edt=zero
c
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
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
      call timer_stop(1)
c
c     outer loop
c
  100 continue
c
      if(icgit.eq.(icgit/ncgcy)*ncgcy) beta1=zero
      icgit=icgit+1
      lpair=frequency(icgit,nfpair)
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
      if(pnorm.lt.zero) call errquit('congra: pnorm<zero',0)
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
      call timer_start(1)
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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      lpair=.false.
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
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
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
c
      ecg2=epot
      epcgw=epotw
      epcgsw=epotsw
      epcgs=epots
      ecgdif=ecg2-ecgold
      edt=ecg2-ecg0
      call timer_stop(1)
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
      call timer_start(1)
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
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
c     record to lfnemt
c
      ltr=frequency(icgit,nfem)
c
      if(ltr) then
      if(lqmd) then
      call qmd_wrttrj(lfnemt,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfnemt,.true.,.false.,.true.,.false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
      endif
c
c     write restart file
c
      ldone=.not.(icgit.lt.mcgit.and.dxmax.gt.dxcgmx.and.ecgdif.lt.zero)
c
      cqrs=' '
      if(frequency(icgit,nfqrs).or.(ldone.and.eqrs.gt.epot)) then
      write(projct,4000) msdit,icgit,filnam(1:55)
 4000 format('em sd ',i7,' cg ',i7,' ',a)
      call md_wrtrst(lfnqrs,filqrs,.false.)
      cqrs='*'
      eqrs=ecg1
      endif
      if(me.eq.0.and.frequency(icgit,nfprop)) call prp_record()
c
c     print minimization step data to output
c
      if(me.eq.0) then
      write(lfnout,600) icgit,cqrs,ecgdif,ecg1,epcgw,epcgsw,epcgs,dxmax
  600 format(i7,1x,a1,3x,5(1pe13.5),0pf12.8)
      endif
c
      if(itest.gt.0) call md_test()
c
      call timer_stop(1)
      if(.not.ldone) goto 100
c
      ecgdif=edt
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
      logical frequency
      external frequency
c
      lfirst=.true.
      lxw=.false.
      lvw=.false.
      lxs=.false.
      lvs=.false.
c
c     equilibration
c
      lequi=.true.
      do 1 iequi=kequi+1,mequi
      call timer_start(1)
      call md_newton()
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
      call md_server
      call timer_stop(1)
    1 continue
c
c     data gathering
c
      lequi=.false.
      do 2 idacq=kdacq+1,mdacq
      call timer_start(1)
c
      call md_newton()
c
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
      call md_server
      call prp_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm)
c
      lxw=frequency(mdstep,nfcoor)
      lvw=frequency(mdstep,nfvelo)
      lxs=frequency(mdstep,nfscoo)
      lvs=frequency(mdstep,nfsvel)
c
      if(me.eq.0.and.idacq.eq.mdacq.and.(lxw.or.lvw.or.lxs.or.lvs)) then
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
      endif
c
      write(projct,4000) mequi,idacq,tmpext,filnam(1:48)
 4000 format('md ',i7,'+ ',i7,' @ ',f7.2,' K ',a)
      if(frequency(mdstep,nfrest)) call md_wrtrst(lfnrst,filrst,.true.)
      if(frequency(mdstep,nftime)) call md_wrtime
      if(frequency(mdstep,nfanal)) call ana_data(stime,nsaloc,
     + int_mb(i_is+(lsgan-1)*msa),dbl_mb(i_xs))
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
#include "global.fh"
#include "msgids.fh"
c
      logical frequency,prp_mcti_step
      external frequency,prp_mcti_step
c
      logical done
      integer npp,irunp
c
      if(me.eq.0) then
      if(mropt.ge.2) then
      open(unit=lfnmri,file=filmri(1:index(filmri,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmri) npp
      if(npp.ne.np) call errquit('Number of nodes changed',npp)
      krun=0
      endif
      open(unit=lfnmro,file=filmro(1:index(filmro,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      write(lfnmro) np,mrun,mequi,mdacq
      endif
      call ga_brdcst(mrg_d30,krun,ma_sizeof(mt_int,1,mt_byte),0)
c
      do 1 irun=krun+1,mrun
c
      lfirst=.true.
c
c     initialize parameters
c
      call cf_lambda(lamtyp,irun,maxlam,elam,lfnout,rlambd,dlambd)
c
c     property initialization
c
      call prp_init()
c
      if(mropt.ge.2) then
      if(me.eq.0) then
      read(lfnmri) irunp,kequi,kdacq
      if(irunp.ne.irun) call errquit('Number of run changed',irunp)
      if(kequi.lt.mequi) kdacq=0
      call prp_rdmri(lfnmri)
      endif
      call sp_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      call ga_brdcst(mrg_d31,kequi,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mrg_d32,kdacq,ma_sizeof(mt_int,1,mt_byte),0)
      endif
c
      lxw=.false.
      lvw=.false.
      lxs=.false.
      lvs=.false.
c
c     equilibration
c
      lequi=.true.
      do 2 iequi=kequi+1,mequi
      call timer_start(1)
      call md_newton()
      call timer_stop(1)
    2 continue
c
c     data gathering
c
      mdstep=kdacq
      if(kdacq.eq.0) stime=zero
c
      lequi=.false.
      do 3 idacq=kdacq+1,mdacq
c
      call timer_start(1)
c
      call md_newton()
c
      done=prp_mcti_step(idacq)
c
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot)
      call prp_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_esk),epme,uqmd,uqmmm)
c
      lxw=frequency(mdstep,nfcoor)
      lvw=frequency(mdstep,nfvelo)
      lxs=frequency(mdstep,nfscoo)
      lvs=frequency(mdstep,nfsvel)
c
      if(me.eq.0.and.idacq.eq.mdacq.and.(lxw.or.lvw.or.lxs.or.lvs)) then
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
      endif
c
      write(projct,4000) irun,mrun,mequi,idacq,tmpext,filnam(1:35)
 4000 format('ti ',i5,'/',i5,' :',i7,'+ ',i7,' @ ',f7.2,' K ',a)
      if(frequency(mdstep,nfrest)) call md_wrtrst(lfnrst,filrst,.true.)
      if(frequency(mdstep,nftime)) call md_wrtime
c
      call timer_stop(1)
c
      if(done) goto 4
c
    3 continue
    4 continue
c
      if(me.eq.0) then
      write(lfnmro) irun,mequi,mdacq
      call prp_wrtmro(lfnmro)
      call sp_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
c
      call prp_mcti_run(rlambd,dlambd)
c
    1 continue
c
      call prp_mcti()
c
      return
c
 9999 continue
      call errquit('Failed to open file mro',0)
      return
      end
      subroutine md_newton()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      logical frequency
      external frequency
c
      external timer_wall,timer_wall_average,timer_wall_minimum
      real*8 timer_wall,timer_wall_average,timer_wall_minimum
c
      real*8 dxmax,wallt
c
      mdstep=mdstep+1
      stime=stime+tstep
c
      wallt=timer_wall(3)
      call timer_start(3)
c
c     reassign velocities
c
      if(frequency(mdstep,nfgaus)) then
      call cf_gauss(tgauss,frgaus,nwmloc,nsaloc,
     + dbl_mb(i_vw),dbl_mb(i_vs),int_mb(i_iw+(lwdyn-1)*mwm),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      endif
c
      lpair=lfirst.or.frequency(mdstep,nfpair)
      llong=(lfirst.or.frequency(mdstep,nflong).or.lpair).and.ltwin
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
      if(.not.lqmd) then
c
c     dynamic load balancing
c
      if(itload.eq.0) then
      if(me.eq.0) then
      call sp_balanc(stime,timer_wall_minimum(4),wallt,
     + frequency(mdstep,nfsync))
      else
      call sp_balanc(stime,timer_wall(4),wallt,frequency(mdstep,nfsync))
      endif
      else
      call sp_balanc(stime,timer_wall_average(4),wallt,
     + frequency(mdstep,nfsync))
      endif
      call timer_reset(4)
c
c     atom redistribution
c
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
      endif
c
      endif
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     subtract center of mass coordinates from reference coordinates
c
      call md_addref(.false.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      if(me.eq.0.and..not.lequi.and.(lxw.or.lvw.or.lxs.or.lvs)) then
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfntrj,lxw,lvw,lxs,lvs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))
      endif
      endif
c
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     self-guided forces
c
      if(iguide.gt.0) then
      call md_guided(dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_gw),dbl_mb(i_gs))
      endif
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
c     coordinate scaling
c
      call cf_final(dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_yw),
     + dbl_mb(i_vw),dbl_mb(i_vwt),nwmloc,
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_ys),dbl_mb(i_vs),
     + dbl_mb(i_vst),int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + int_mb(i_is+(lsfrc-1)*msa),dbl_mb(i_zs),
     + dbl_mb(i_esk),nsaloc,box,vlat,pres,temp,tempw,temps)
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
c
c     add center of mass coordinates from reference coordinates
c
      call md_addref(.true.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
c     center solute in box
c
      if(frequency(mdstep,nfcntr)) then
      call cf_center(dbl_mb(i_xw),nwmloc,int_mb(i_is+(lsfrc-1)*msa),
     + dbl_mb(i_xs),nsaloc,idscb,nscb)
      endif
c
c     remove overall translational motion
c
      if(frequency(mdstep,nfslow)) then
      call cf_slow(dbl_mb(i_vw),nwmloc,dbl_mb(i_vs),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
      endif
c
c     velocity scaling to preset temperature
c
      if(frequency(mdstep,nfgaus)) then
      call cf_vscale(tgauss,nwmloc,nsaloc,
     + dbl_mb(i_vw),dbl_mb(i_vs),int_mb(i_iw+(lwdyn-1)*mwm),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      endif
c
c     update decomposition module
c
      if(.not.lqmd) call sp_update(me,vlat,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_xwcr),dbl_mb(i_vw),nwmloc,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),nsaloc)
c
      if(itest.gt.0) call md_test()
c
      lfirst=.false.
c
      call timer_stop(3)
c
      return
      end
      subroutine md_addref(ladd,xwm,xwcr,xsm,xscr,dsr)
c
      implicit none
c
#include "md_common.fh"
c
      logical ladd
      real*8 xwm(mwm,3),xwcr(mwm,3),xsm(msm,3),xscr(msm,3),dsr(msm)
c
      integer i
c
c      return
      if(ladd) then
      dwr=zero
      do 1 i=1,nwmloc
      xwcr(i,1)=xwcr(i,1)+xwm(i,1)
      xwcr(i,2)=xwcr(i,2)+xwm(i,2)
      xwcr(i,3)=xwcr(i,3)+xwm(i,3)
      dwr=dwr+xwcr(i,1)*xwcr(i,1)+xwcr(i,2)*xwcr(i,2)+
     + xwcr(i,3)*xwcr(i,3)
    1 continue
      do 2 i=1,nsm
      xscr(i,1)=xscr(i,1)+xsm(i,1)
      xscr(i,2)=xscr(i,2)+xsm(i,2)
      xscr(i,3)=xscr(i,3)+xsm(i,3)
      dsr(i)=xscr(i,1)*xscr(i,1)+xscr(i,2)*xscr(i,2)+xscr(i,3)*xscr(i,3)
    2 continue
      else
      do 3 i=1,nwmloc
      xwcr(i,1)=xwcr(i,1)-xwm(i,1)
      xwcr(i,2)=xwcr(i,2)-xwm(i,2)
      xwcr(i,3)=xwcr(i,3)-xwm(i,3)
    3 continue
      do 4 i=1,nsm
      xscr(i,1)=xscr(i,1)-xsm(i,1)
      xscr(i,2)=xscr(i,2)-xsm(i,2)
      xscr(i,3)=xscr(i,3)-xsm(i,3)
    4 continue
      endif
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
      logical cf_shakep
      external cf_shakep
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
    3 continue
      do 4 ibbl=1,nbbl
      call sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
      if(lself) call cf_shakes(lhandl,xs,ys,isl(1,lsgan),isl(1,lsatt),
     + isl(1,lsdyn),isfr,isto)
    4 continue
      do 5 i=1,nsaloc
      dmax=max(dmax,(xs(i,1)-ys(i,1))**2+(xs(i,2)-ys(i,2))**2+
     + (xs(i,3)-ys(i,3))**2)
    5 continue
      if(.not.cf_shakep(xs,ys,isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),
     + nsaloc)) goto 3
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
      subroutine md_finit(iwl,isl,xw,xwm,xs,fw,fs)
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
c     initialize cf
c
      call cf_init(lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),eww,
     + dbl_mb(i_esw),dbl_mb(i_ess))
c
      if(lpair) call md_zinit(int_mb(i_iwz),int_mb(i_isz))
c
      if(lqmd) then
      call qmd_forces(irtdb,int_mb(i_is+(lsgan-1)*msa),
     + int_mb(i_is+(lsatt-1)*msa),dbl_mb(i_xs),dbl_mb(i_fs),msa,
     + nsaloc,uqmd)
      uqmd=uqmd-uqmatm
      return
      endif
c
      call ga_sync()
c
      call sp_initf(fw,fs,llong,int_mb(i_iwz),int_mb(i_isz),lpair)
c
      call sp_putix(me,iwl,xw,nwmloc,isl,xs,nsaloc)
c
      call ga_sync()
c
      call timer_start(10)
      call cf_multi(nsaloc,dbl_mb(i_xs),dbl_mb(i_fs),isl(1,lsgan),
     + isl(1,lsfrc),isl(1,lsdyn),isl(1,lsct1),dbl_mb(i_ess))
      call timer_stop(10)
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
      if(.not.lqmd) then
c
      call md_fclass(iwl,isl,xw,xwm,xs,fw,fs)
c
      call timer_start(2)
      call timer_start(4)
      call ga_sync()
      call timer_stop(4)
      call timer_stop(2)
c
      call sp_final(fw,fs,lpair,int_mb(i_iwz),int_mb(i_isz))
c
      if(lqmmm) then
      call qmmm_forces(irtdb,mwm,nwmloc,mwa,nwa,int_mb(i_iwz),xw,fw,
     + msa,nsaloc,int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + int_mb(i_is+(lsct1-1)*msa),int_mb(i_isz),xs,fs,uqmmm)
      uqmmm=uqmmm-uqmatm
      endif
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
      endif
c
      call cf_fnorm(iwl(1,lwdyn),fw,nwmloc,
     + isl(1,lsatt),isl(1,lsdyn),fs,nsaloc,fnorm,fmax)
c
      return
      end
      subroutine md_guided(fw,fs,gw,gs)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 fw(mwm,3,mwa,2),fs(msa,3,2)
      real*8 gw(mwm,3,mwa),gs(msa,3)
c
      integer i,j,k
c
      do 1 k=1,nwa
      do 2 j=1,3
      do 3 i=1,nwmloc
      gw(i,j,k)=factgg*gw(i,j,k)+factgf*(fw(i,j,k,1)+fguide*gw(i,j,k))
      fw(i,j,k,1)=fw(i,j,k,1)+fguide*gw(i,j,k)
    3 continue
    2 continue
    1 continue
c
      do 4 j=1,3
      do 5 i=1,nsaloc
      gs(i,j)=factgg*gs(i,j)+factgf*(fs(i,j,1)+fguide*gs(i,j))
      fs(i,j,1)=fs(i,j,1)+fguide*gs(i,j)
    5 continue
    4 continue
c
      return
      end
      subroutine md_zinit(iwz,isz)
c
      implicit none
c
#include "md_common.fh"
c
      integer iwz(mwm),isz(msa)
c
      integer i
c
      do 1 i=1,mwm
      iwz(i)=0
    1 continue
c
      do 2 i=1,msa
      isz(i)=0
    2 continue
c
      return
      end
      subroutine md_fclass(iwl,isl,xw,xwm,xs,fw,fs)
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
      logical lself,local,locali,lpbcs
c
      integer ibbl,lhandl,iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer nbbl
c
      if(ipme.gt.0) then
      if(lpert2) then
      call pme_energy(2,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),nsaloc,epme(2))
      call pme_init()
      endif
      if(lpert3) then
      call pme_energy(3,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),nsaloc,epme(3))
      call pme_init()
      endif
      call pme_chgrid(iset,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),nsaloc,epme(iset))
      endif
c
      locali=.true.
      call timer_start(13)
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + isl(1,lsatt),isl(1,lsmol),dbl_mb(i_xs),dbl_mb(i_xsm))
c
      if(ipolt.gt.0) then
      call md_induce(iwl,isl,xw,xwm,xs,dbl_mb(i_pw),
     + dbl_mb(i_pwp),dbl_mb(i_ps),dbl_mb(i_psp))
      endif
c
      call sp_nbbl(nbbl)
c
      do 1 ibbl=1,nbbl
c
      call sp_getxbl(ibbl,lhandl,local,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lself,lpbcs)
c
      if(locali.and..not.local) then
      call timer_stop(13)
      call timer_start(14)
      locali=.false.
      endif
c
      call cf_comw(xw,xwm,jwfr,jwto)
c
      call forces(lself,lpbcs,xw,xwm,fw,zw,dbl_mb(i_rtos),
     + iwl(1,lwdyn),int_mb(i_iwz),
     + iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),fs,dbl_mb(i_zs),
     + isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),isl(1,lsfrc),
     + isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),isl(1,lsct3),
     + int_mb(i_isz),isfr,isto,jsfr,jsto,lpbc,lhandl,
     + .true.,eww,dbl_mb(i_esw),dbl_mb(i_ess))
c
      call sp_accfbl(ibbl,lhandl,fw,fs,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
c
    1 continue
c
      if(locali) then
      call timer_stop(13)
      else
      call timer_stop(14)
      endif
c
      if(ipme.gt.0) then
      call pme_forces(fw,nwmloc,fs,isl(1,lsct1),nsaloc)
      endif
c
      return
      end
      subroutine md_induce(iwl,isl,xw,xwm,xs,pw,pwp,ps,psp)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      integer iwl(mwm,miw2),isl(msa,mis2)
      real*8 xw(mwm,3,mwa),xwm(mwm,3),xs(msa,3)
      real*8 pw(mwm,3,mwa,2),ps(msa,3,2)
      real*8 pwp(mwm,3,mwa,2,2),psp(msa,3,2,2)
c
      logical lself,local,lpbcs
      integer npolit,ibbl,nbbl,lhandl
      integer iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto
      integer i,j,k
      real*8 pmax
c
c     initialize induced fields to zero for first order polarization
c
      if(ipolt.eq.1) then
      if(nwmloc.gt.0) then
      do 1 k=1,nwa
      do 2 j=1,3
      do 3 i=1,nwmloc
      pw(i,j,k,1)=zero
      pw(i,j,k,2)=zero
    3 continue
    2 continue
    1 continue
      endif
      if(nsaloc.gt.0) then
      do 4 j=1,3
      do 5 i=1,nsaloc
      ps(i,j,1)=zero
      ps(i,j,2)=zero
    5 continue
    4 continue
      endif
      if(lpert2.or.lpert3) then
      if(nwmloc.gt.0) then
      do 6 k=1,nwa
      do 7 j=1,3
      do 8 i=1,nwmloc
      pwp(i,j,k,1,1)=zero
      pwp(i,j,k,1,2)=zero
      pwp(i,j,k,2,1)=zero
      pwp(i,j,k,2,2)=zero
    8 continue
    7 continue
    6 continue
      endif
      if(nsaloc.gt.0) then
      do 9 j=1,3
      do 10 i=1,nsaloc
      psp(i,j,1,1)=zero
      psp(i,j,1,2)=zero
      psp(i,j,2,1)=zero
      psp(i,j,2,2)=zero
   10 continue
    9 continue
      endif
      endif
      endif
c
c     iterative cycle to generate induced fields
c     ------------------------------------------
c
      npolit=0
      call sp_nbbl(nbbl)
   11 continue
c
c     copy fields from previous iteration
c     -----------------------------------
c
      if(nwmloc.gt.0) then
      do 12 k=1,nwa
      do 13 j=1,3
      do 14 i=1,nwmloc
      pw(i,j,k,2)=pw(i,j,k,1)
      pw(i,j,k,1)=zero
   14 continue
   13 continue
   12 continue
      endif
      if(nsaloc.gt.0) then
      do 15 j=1,3
      do 16 i=1,nsaloc
      ps(i,j,2)=ps(i,j,1)
      ps(i,j,1)=zero
   16 continue
   15 continue
      endif
      if(mdtype.gt.3) then
      if(nwmloc.gt.0) then
      do 17 k=1,nwa
      do 18 j=1,3
      do 19 i=1,nwmloc
      pwp(i,j,k,1,2)=pwp(i,j,k,1,1)
      pwp(i,j,k,1,1)=zero
      pwp(i,j,k,2,2)=pwp(i,j,k,2,1)
      pwp(i,j,k,2,1)=zero
   19 continue
   18 continue
   17 continue
      endif
      if(nsaloc.gt.0) then
      do 20 j=1,3
      do 21 i=1,nsaloc
      psp(i,j,1,2)=psp(i,j,1,1)
      psp(i,j,1,1)=zero
      psp(i,j,2,2)=psp(i,j,2,1)
      psp(i,j,2,1)=zero
   21 continue
   20 continue
      endif
      endif
c
c     copy current fields into local global array
c     -------------------------------------------
c
      call sp_putp(me,pw,pwp,nwmloc,ps,psp,nsaloc,lpert2.or.lpert3)
c
c     synchronize to ensure induced fields are available
c     --------------------------------------------------
c
      call ga_sync()
c
      do 22 ibbl=1,nbbl
c
      call sp_getxpbl(ibbl,lhandl,local,
     + iwl,xw,pw,pwp,iwfr,iwto,jwfr,jwto,
     + isl,xs,ps,psp,isfr,isto,jsfr,jsto,lself,lpbcs)
c
      call cf_comw(xw,xwm,jwfr,jwto)
c
      call induce(lself,lpbcs,xw,xwm,pw,pwp,iwl(1,lwdyn),int_mb(i_iwz),
     + iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),ps,psp,
     + isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),isl(1,lsfrc),
     + isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),isl(1,lsct3),
     + int_mb(i_isz),isfr,isto,jsfr,jsto,lpbc,lhandl)
c
      call sp_accpbl(ibbl,lhandl,pw,pwp,ps,psp,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
c
   22 continue
c
      pmax=0.0d0
      do 23 k=1,nwa
      do 24 j=1,3
      do 25 i=1,nwmloc
      pmax=max(pmax,abs(pw(i,j,k,2)-pw(i,j,k,1)))
   25 continue
   24 continue
   23 continue
c
      if(np.gt.0) call ga_dgop(mrg_d06,pmax,1,'max')
c
      write(*,'(7f12.6)') pmax,(pw(1,j,1,1),j=1,3),(pw(1,j,1,2),j=1,3)
c
      if(pmax.gt.ptol) goto 11
c
c     if not converged goto 11
c
c
      return
      end
      subroutine md_wrtrst(lfn,fil,lveloc)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      integer lfn
      character*255 fil
      character*7 cdr
      logical lveloc
c
      if(me.eq.0) then
      if(keepr.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown',err=9999)
      else
      write(cdr,'(i7.7)') keepr
      open(unit=lfn,file=fil(1:index(fil,' ')-1)//cdr,
     + form='formatted',status='unknown',err=9999)
      endif
      endif
c
      call sp_wrtrst(lfn,fil,lveloc,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_gw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),projct)
c
      call md_wtrest(lfn)
c
      call prp_wtrest(lfn)
c
      call sp_wtrest(lfn)
c
      if(me.eq.0) then
      close(unit=lfn)
      endif
c
      if(keepr.gt.0) keepr=keepr+1
c
      return
c
 9999 continue
      call errquit('Unable to open restart for writing',me)
      return
      end
      subroutine md_server
c
      implicit none
c
#include "md_common.fh"
c
      integer*4 create_client_socket
      integer client_socket_write
      external create_client_socket,client_socket_write
c
      integer*4 iiport,lens
c
      character*255 string
      integer numbyt
c
      if(iport.le.0.or.me.gt.0) return
c
c     open socket to server
c
      iiport=iport
      if(.not.lserver) then
      write(*,1) server(1:index(server,' ')-1),iiport
    1 format('Attempt to open socket to ',a,' port ',i5)
      isocket=create_client_socket(server,iiport)
      lserver=isocket.gt.0
      endif
c
      if(lserver) then
      print*,'server socket open'
      else
      print*,'server socket error'
      endif
c
      if(.not.lserver) return
c
      write(string,1000) stime,etot,temp,pres*1.0125d-5,volume
 1000 format("tETpV",4f12.3,f12.6)
c
c      write(string,1000) stime,temp,pres*1.0125d-5,volume
c 1000 format("TIM",f12.6,"TMP",f12.6,"PRS",f12.6,"VOL",f12.6)
c      write(*,'(a)') string
      lens=66
      string(66:66)=char(13)
c
      numbyt=client_socket_write(isocket,string,lens)
c
      print*,'Bytes written to socket is ',numbyt
c
      return
      end
      subroutine md_wrtime
c
      implicit none
c
#include "md_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      external timer_wall
      real*8 timer_wall
c
      integer i,j
      real*8 tim(8,1024)
c
      do 1 i=1,np
      do 2 j=1,8
      tim(j,i)=zero
    2 continue
    1 continue
c
      tim(1,me+1)=timer_wall(10)
      tim(2,me+1)=timer_wall(11)
      tim(3,me+1)=timer_wall(12)
      tim(4,me+1)=timer_wall(13)
      tim(5,me+1)=timer_wall(14)
      tim(6,me+1)=timer_wall(15)
      tim(7,me+1)=timer_wall(16)
      tim(8,me+1)=timer_wall(2)
c
      call ga_dgop(mrg_d04,tim,8*np,'+')
c
      if(me.eq.0) then
      write(lfntim,1000) stime
 1000 format('timings',/,f12.6)
      write(lfntim,1001) ((tim(i,j),i=1,8),j=1,np)
 1001 format(8f10.3)
      endif
c
      return
      end
      subroutine md_test
c
c $Id: md_main.F,v 1.71 2000-04-13 20:22:54 d3j191 Exp $
c
      implicit none
c
#include "md_common.fh"
c
      itest=itest-1
c
      if(me.ne.0) return
c
      if(iquant.eq.0) then
c
      if(ntype.eq.1) then
      write(lfntst,1000) etot
 1000 format('Energy           = ',1pe12.3)
      endif
c
      if(ntype.eq.2) then
      write(lfntst,1100) stime,temp,volume,pres,etot
 1100 format('Time             = ',0pf9.3,/,
     + 'Temperature      = ',0pf8.2,/,
     + 'Volume           = ',0pf8.2,/,
     + 'Pressure         = ',1pe12.2,/,
     + 'Energy           = ',1pe12.3,/)
      endif
c
      if(ntype.eq.3) write(lfntst,1200) isdit+icgit,etot
 1200 format('Iteration        = ',i10,/,
     + 'Energy           = ',1pe12.3,/)
c
      else
c
      if(ntype.eq.1) then
      write(lfntst,2000) etot
 2000 format('Energy           = ',1pe12.3)
      endif
c
      if(ntype.eq.2) then
      write(lfntst,2100) stime,temp,etot
 2100 format('Time             = ',0pf9.3,/,
     + 'Temperature      = ',0pf8.2,/,
     + 'Energy           = ',1pe12.3,/)
      endif
c
      endif
c
      if(itest.eq.0) close(lfntst)
c
      return
      end
