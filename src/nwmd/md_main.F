      subroutine md_main()
c
c $Id$
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
      if(nftri.eq.0) then
      call md_sp()
      else
      call md_spi()
      endif
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
      call md_abort('Unknown calculation type',ntype)
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
      lpair=.true.
      lload=.false.
      lhop=.false.
      llong=ltwin
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
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
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
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
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      if(.not.lqmmm) call prp_print()
c
      call rtdb_put(irtdb,'md:energy',mt_dbl,1,epot)
c
c     print energies
c
      if(.not.lqmmm) then
      call cf_print_energy(lfnout)
      call sp_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
      endif
c
      if(ifidi.ne.0) then
      call md_fd(int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),dbl_mb(i_fs),
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw))
      endif
c
      if(itest.eq.1) call md_test()
c
      return
      end
      subroutine md_sp_qmmm()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
      lpair=.false.
      lload=.false.

c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))

      call rtdb_put(irtdb,'md:energy',mt_dbl,1,epot)

      return
      end
      subroutine md_spi()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      external sp_rdtrj,sp_skip,frequency
      logical sp_rdtrj,sp_skip,frequency
c
      integer import
c
      if(impfr.gt.1) then
      if(.not.sp_skip(lfntri,impfr-1)) return
      endif
c
      do 1 import=impfr,impto
c
      if(.not.frequency(import+1-impfr,nftri)) then
      if(.not.sp_skip(lfntri,1)) return
      goto 1
      endif
c
      lpair=.true.
      lload=.true.
      lhop=.false.
      llong=ltwin
c
c     read frame from trajectory file
c
      if(.not.sp_rdtrj(lfntri,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs),nwmloc,nsaloc)) return
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
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
      call cf_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call md_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(import,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      call prp_step(import,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm)
c
cx      call prp_print()
c
      call rtdb_put(irtdb,'md:energy',mt_dbl,1,epot)
c
c     print energies
c
      if(me.eq.0.and.frequency(import,npener))
     + call cf_print_energy(lfnout)
c      if(me.eq.0.and.frequency(import,nfprop)) call prp_record()
c
      if(me.eq.0.and.frequency(import,npforc))
     + call sp_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
c
      if(ifidi.ne.0) then
      call md_fd(int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),dbl_mb(i_fs),
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw))
      endif
c
      if(itest.eq.1) call md_test()
c
    1 continue
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
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
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
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
      call cf_mass(dbl_mb(i_wws),dbl_mb(i_wws+mwa),
     + int_mb(i_is+(lsatt-1)*msa),nsaloc)
c
      call md_eminit(dbl_mb(i_xw),dbl_mb(i_yw),
     + dbl_mb(i_xs),dbl_mb(i_ys))
c
      if(msdit.gt.0) call md_stdesc(int_mb(i_iw+(lwdyn-1)*mwm),
     + dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs),dbl_mb(i_fs),
     + dbl_mb(i_wws),dbl_mb(i_wws+mwa),int_mb(i_mm),dbl_mb(i_fm),
     + dbl_mb(i_xsm))
c
      if(mcgit.gt.0) then
      if(.not.ma_push_get(mt_dbl,3*mwa*mwm,'pcgw',l_pcgw,i_pcgw))
     + call md_abort('Failed to allocate memory for pcgw',me)
      if(.not.ma_push_get(mt_dbl,3*msa,'pcgs',l_pcgs,i_pcgs))
     + call md_abort('Failed to allocate memory for pcgs',me)
      call md_congra(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
     + dbl_mb(i_yw),dbl_mb(i_vw),dbl_mb(i_fw),dbl_mb(i_pcgw),
     + int_mb(i_is+(lsdyn-1)*msa),dbl_mb(i_xs),dbl_mb(i_ys),
     + dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_pcgs),dbl_mb(i_wws),
     + dbl_mb(i_wws+mwa))
      if(.not.ma_pop_stack(l_pcgs))
     + call md_abort('Failed to deallocate memory for pcgs',me)
      if(.not.ma_pop_stack(l_pcgw))
     + call md_abort('Failed to deallocate memory for pcgw',me)
      endif
c
      call cf_print_energy(lfnout)
c
      call sp_printf(filtop,lfntop,
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_fs),npener,dbl_mb(i_esa))
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
      return
c
      lpair=.true.
      lload=.true.
      lhop=.false.
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
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
      subroutine md_stdesc(iwdt,xw,yw,vw,fw,isdt,ismol,
     + xs,ys,vs,fs,ww,ws,mm,fm,xsm)
c
      implicit none
c
#include "md_common.fh"
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
      if(imembr.gt.0) call md_membrane_init(ismol,mm,xs,xsm,fm)
c
      if(lqmmm) then
        if(me.eq.0) write(6,601) "@","@"
      end if
601   format(
     $     /,a1,' Step       Energy      Delta E   SGrms',
     $     '     WGrms     Xmax',
     $     /,a1,' ---- ---------------- -------- --------',
     $     ' -------- -------- ')

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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      if(lqmmm) call qmmm_print_energy1(irtdb,epot,uqmmm)
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
      if(lqmmm) call qmmm_print_energy1(irtdb,epot,uqmmm)
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
      call cf_restor(dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_fw),
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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
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
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,box,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfntrj,lxw,.false.,.false.,lxs,.false.,.false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
      endif
c
      ldone=.not.(dxmax.gt.dxsdmx.and.isdit.lt.msdit)
c
      cqrs=' '
      if(frequency(isdit,nfqrs).or.(ldone.and.eqrs.gt.epot)) then
      write(projct,4000) nserie,isdit,0,filnam(1:56)
 4000 format(i2,' em ',i7,' + ',i7,' ',a)
      call md_wrtrst(lfnqrs,filqrs,.false.)
      cqrs='*'
      eqrs=epot
      endif
      if(me.eq.0.and.frequency(isdit,nfprop)) call prp_record()
c
      if(me.eq.0) then
      write(lfnout,600) isdit,cqrs,edif,epsd,epsdw,epsdsw,epsds,dxmax
  600 format(i7,1x,a1,3x,5(1pe13.5),0pf12.8)
      call util_flush(lfnout)
      endif
c
      if(lqmmm) then
      nrms(2) = 0
      grms(2) = 0.0d0 
      if(nwmloc.gt.0) then
        do j=1,nwa
          do i=1,nwmloc
            if(iand(iwdt(i),mfixed).ne.lfixed) then
              nrms(2) = nrms(2) + 1
              grms(2) = grms(2) + 
     +               fw(i,1,j)*fw(i,1,j) +
     +               fw(i,2,j)*fw(i,2,j) + 
     +               fw(i,3,j)*fw(i,3,j) 
            endif
          end do
        end do
      endif
c
      nrms(1) = 0
      grms(1) = 0.0d0
      do i=1,nsaloc
        if(iand(isdt(i),mfixed).ne.lfixed) then
          nrms(1) = nrms(1) + 1
          grms(1) = grms(1) + 
     +           fs(i,1)*fs(i,1)+
     +           fs(i,2)*fs(i,2)+
     +           fs(i,3)*fs(i,3)
        endif
      end do
      call ga_dgop(msg_qmmm_force,grms,2,'+')
      call ga_igop(msg_qmmm_nat1,nrms,2,'+')

      do i=1,2
      if(nrms(i).ne.0) then
        grms(i) = sqrt(grms(i)/dble(nrms(i)))
        grms(i) = grms(i)*5.29177249d-02/2.625499962d+03
      end if
      end do
      if(me.eq.0)    
     $  write(6,602) "@", isdit, epsd/2.625499962d+03,
     $  edif/2.625499962d+03,
     $  grms(1), grms(2), 
     $  dxmax/5.29177249d-02
602   format(a1,i5,f17.8,1p,d9.1,0p,3f9.5)

      end if
c

      if(itest.gt.0) call md_test()
c
      dx=min(1.2d0*dx,dxmsd)
      da=min(1.2d0*da,damsd)
c
      call timer_stop(201)
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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
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
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
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
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(isdit,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
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
c     record to lfntrj
c
      lxw=frequency(mdstep,nfcoor)
      lxs=frequency(mdstep,nfscoo)
c
      if(lxw.or.lxs) then
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,box,stime,pres,temp,tempw,temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      call sp_wrttrj(lfntrj,lxw,.false.,.false.,lxs,.false.,.false.,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
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
      call timer_stop(201)
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
      call md_newton()
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call md_server
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
      call md_timer_init()
c
      call timer_start(201)
c
      call md_newton()
c
      call timer_start(6)
c
      if(lfw.or.lfs) then
      call sp_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),nsaloc)
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
      call md_wrtrst(lfnrst,rfile,.true.)
      endif
      if(frequency(mdstep,nftime)) call md_wrtime
c
      call timer_stop(55)
c
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call md_server
      call prp_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm)
c
      call timer_stop(201)
c
      if(frequency(mdstep,nfnewf).and.idacq.ne.mdacq) then
      if(me.eq.0) call md_fopen(.true.)
      endif
c
c
      if(lstop.and.(.not.lqmmm)) then
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
      subroutine md_ti()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      logical frequency,prp_mcti_step,prp_rdmri,sp_rdmri
      integer prp_dfr
      external frequency,prp_mcti_step,prp_rdmri,sp_rdmri,prp_dfr
c
      logical done
      integer npp,i,irunp,jequi,jdacq,ndec,ndum
      real*8 rdum
      character*256 filrun
c
      if(nserie.eq.0) then
      if(mropt.ge.2) then
      if(me.eq.0) then
      open(unit=lfnmri,file=filmri(1:index(filmri,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmri) npp
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      endif
      krun=0
      endif
      if(me.eq.0) then
      open(unit=lfnmro,file=filmro(1:index(filmro,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      write(lfnmro) np,mrun,mequi,mdacq
      endif
      else
c
      if(me.eq.0) then
      open(unit=lfnmro,file=filmro(1:index(filmro,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmro) npp,mrun,mequi,mdacq
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      endif
      krun=0
      do 11 irun=1,mrun
      if(me.eq.0) then
      read(lfnmro,end=12,err=12) irunp,kequi,kdacq
      endif
      if(.not.prp_rdmri(lfnmro,ndec,mropt)) goto 12
      if(.not.sp_rdmri(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))) goto 12
      krun=irun
   11 continue
   12 continue
c
      if(me.eq.0) then
      open(unit=lfnmri,file=filmri(1:index(filmri,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmri) npp
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      rewind(lfnmro)
      read(lfnmro) npp,mrun,mequi,mdacq
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      rewind(lfngib)
      endif
c
      do 13 irun=1,mrun
      if(me.eq.0) then
      read(lfnmri) irunp,kequi,kdacq
      endif
      if(.not.prp_rdmri(lfnmri,ndec,mropt))
     + call md_abort('Error in mri file',0)
      if(.not.sp_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))) 
     + call md_abort('Error in mri file',0)
      read(lfnmro) irunp,kequi,kdacq
      if(.not.prp_rdmri(lfnmro,ndec,mropt))
     + call md_abort('Error in mro file',0)
      if(.not.sp_rdmri(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))) 
     + call md_abort('Error in mro file',0)
      if(me.eq.0) then
      read(lfngib,2000) ndum
 2000 format(7x,i7)
      read(lfngib,2001) (rdum,i=1,24)
 2001 format(4e20.12)
      read(lfngib,2001) (rdum,i=1,ndum)
      read(lfngib,2001) (rdum,i=1,ndum)
      read(lfngib,2002) ndum
 2002 format(i10)
      read(lfngib,2003) rdum
      read(lfngib,2003) rdum
 2003 format(e20.12)
      if(ndec.gt.0) read(lfngib,2001) (rdum,i=1,nsa)
      endif
   13 continue
c
      endif
c
      call ga_brdcst(mrg_d30,krun,ma_sizeof(mt_int,1,mt_byte),0)
c
      do 1 irun=krun+1,mrun
c
      if(npg.gt.1) then
      if(irun.ne.ipg+1) goto 1
      endif
c
      if(me.eq.0) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filtrj(1:index(filtrj,'.trj')-1),'-',
     + irun,'.trj '
      open(unit=lfntrj,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      call cf_trjhdr(lfntrj)
      endif
      if(nfprop.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filprp(1:index(filprp,'.prp')-1),'-',
     + irun,'.prp '
      open(unit=lfnprp,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nfrdf.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filrdf(1:index(filprp,'.prp')-1),'-',
     + irun,'.rdf '
      open(unit=lfnrdf,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      endif
      endif
c
      if(irun.eq.1.and.iand(ivopt,2).eq.2) nfgaus=ivreas
      if(irun.eq.maxlam.and.iand(ivopt,4).eq.4) nfgaus=ivreas
c
      lfirst=.true.
c
c     initialize parameters
c
      call cf_lambda(lamtyp,irun,maxlam,elam,lfnout,lfnpmf,
     + rlambd,dlambd,filnam)
c
c     property initialization
c
      call prp_init()
c
      if(mropt.ge.2) then
      if(me.eq.0) then
      read(lfnmri,end=399,err=399) irunp,kequi,kdacq
      if(irunp.ne.irun) call md_abort('Number of run changed',irunp)
      if(.not.prp_rdmri(lfnmri,ndec,mropt))
     + call md_abort('Error in mri',0)
      if(kequi+kdacq.lt.mequi) then
      kequi=kequi+kdacq
      kdacq=0
      elseif(kequi.lt.mequi) then
      if(kdacq+kequi-mequi.gt.0) then
      kdacq=prp_dfr(kdacq+kequi-mequi)
      else
      kdacq=0
      endif
      kequi=mequi
      endif
      goto 398
  399 continue
      kequi=0
      kdacq=0
  398 continue
      endif
      call ga_brdcst(mrg_d31,kequi,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mrg_d32,kdacq,ma_sizeof(mt_int,1,mt_byte),0)
c
c     kequi will be 0 if no records could be read from lfnmri, and
c     the coordinates and velocities will be used from the previous 
c     lambda run.
c     if this is the first run there must be something wrong
c
      if(kequi.eq.0.and.kdacq.eq.0) then
      if(irun.eq.1) call md_abort('No records found on mri',me)
      else
      if(.not.sp_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs)))
     + call md_abort('Error reading mri in sp_rdmri',0)
      endif
      endif
c
      lxw=.false.
      lvw=.false.
      lxs=.false.
      lvs=.false.
      lesp=.false.
c
      if(mropt.eq.2) then
      kequi=0
      kdacq=0
      call prp_init
      endif
c
c     equilibration
c
      lequi=.true.
      jequi=kequi
      lprpmf=.false.
      iprpmf=0
      do 2 iequi=kequi+1,mequi
c
      mdstep=mdstep+1
      lpmfc=npmf.gt.1.and.iequi.gt.npmf
      call timer_start(201)
      call md_newton()
      call timer_stop(201)
      jequi=iequi
      stime=stime+tstep
    2 continue
      lpmfc=.true.
c
c     data gathering
c
      mdstep=kdacq
      if(kdacq.eq.0) stime=zero
c
      ndec=0
      if(npgdec.gt.1) call cf_dera_init()
c
      lequi=.false.
      jdacq=kdacq
      lprpmf=.true.
      iprpmf=-1
      do 3 idacq=kdacq+1,mdacq
c
      mdstep=mdstep+1
c
      lxw=frequency(mdstep,nfcoor)
      lvw=frequency(mdstep,nfvelo)
      lfw=frequency(mdstep,nfforc)
      lxs=frequency(mdstep,nfscoo)
      lvs=frequency(mdstep,nfsvel)
      lfs=frequency(mdstep,nfsfor)
c
      call md_timer_init()
c
      call timer_start(201)
c
      call md_newton()
c
      if(npgdec.gt.1) ndec=ndec+1
c
      call cf_mcti_kin(int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsgan-1)*msa),dbl_mb(i_vs),nsaloc)
c
      done=prp_mcti_step(idacq,ldacq)
c
      call prp_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call prp_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm)
c
      if(lfw.or.lfs) then
      call sp_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),nsaloc)
      endif
c
      write(projct,4000) nserie,irun,mrun,mequi,idacq,tmpext,
     + filnam(1:32)
 4000 format(i2,' ti ',i5,'/',i5,' :',i7,'+ ',i7,' @ ',f7.2,' K ',a)
      if(frequency(mdstep,nfrest)) call md_wrtrst(lfnrst,rfile,.true.)
      if(frequency(mdstep,nftime)) call md_wrtime
c
      call timer_stop(201)
c
      jdacq=idacq
      if(idacq.gt.ldacq.and.done) goto 4
c
      stime=stime+tstep
c
    3 continue
    4 continue
c
      if(me.eq.0) then
      write(lfnmro) irun,jequi,jdacq
      endif
      call prp_wrtmro(lfnmro,ndec)
      if(me.eq.0) then
      call sp_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),projct)
      endif
c
      call prp_mcti_run(rlambd,dlambd,ndec)
c
      if(ndec.gt.2) call cf_print_dera(lfnout,ndec)
c
      if(me.eq.0) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      close(unit=lfntrj,status='keep')
      endif
      if(nfprop.gt.0) then
      close(unit=lfnprp,status='keep')
      endif
      endif
c
    1 continue
c
      call prp_mcti(npgdec,filnam)
c
      if(me.eq.0) then
      if(nserie.eq.0) then
      if(mropt.ge.2) then
      close(unit=lfnmri,status='keep')
      endif
      close(unit=lfnmro,status='keep')
      else
      close(unit=lfnmro,status='keep')
      endif
      endif
c
      return
c
 9999 continue
      call md_abort('Failed to open file mro',0)
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
      call cf_gauss(tgauss,frgaus,nwmloc,nsaloc,
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
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     periodic boundary conditions
c
      call md_fold(int_mb(i_iw),int_mb(i_is),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_xsm))
c
      call timer_stop(2)
c
      if(lload) then
      if(.not.lqmd) then
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
      call sp_balanc(stime,syntim,wallt,frequency(mdstep,nfsync))
      call timer_reset(204)
c
      call timer_stop(3)
      call timer_start(4)
c
c     atom redistribution
c
      call sp_travel(box,dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + dbl_mb(i_gw),int_mb(i_iw),nwmloc,dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),int_mb(i_is),nsaloc)
c
      call timer_stop(4)
c
      endif
      endif
c
      endif
c
      call timer_start(5)
c
c     center of mass coordinates
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     subtract center of mass coordinates from reference coordinates
c
      if(idifco.gt.0)
     + call md_addref(.false.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
      call timer_stop(5)
      call timer_start(6)
c
      if(lfw.or.lfs) then
      call sp_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),nsaloc)
      call sp_wrttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
c
      call timer_stop(6)
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
c
      corrio=costio
      if(me.eq.0.and..not.lequi.and.(lxw.or.lvw.or.lxs.or.lvs)) then
c
      call timer_start(6)
c
      if(lqmd) then
      call qmd_wrttrj(lfntrj,mwm,nwmloc,mwa,nwa,msa,nsaloc,
     + .true.,.false.,.true.,.false.,box,stime-tstep,pres,temp,tempw,
     + temps,
     + dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xs),dbl_mb(i_vs))
      else
      if(.not.lfw.and..not.lfs) then
      call sp_wrttrj(lfntrj,lxw,lvw,lfw,lxs,lvs,lfs,
     + stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_fw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_fs))
      endif
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
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
c     self-guided forces
c
      if(iguide.gt.0) then
      call timer_start(48)
      call md_guided(dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_gw),dbl_mb(i_gs))
      call timer_stop(48)
      endif
c
c     center of mass options
c
      if(icmopt.gt.0) then
      call timer_start(48)
      call md_cmopt(dbl_mb(i_vs),dbl_mb(i_fs),dbl_mb(i_fcm),
     + int_mb(i_is+(lsmol-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      call timer_stop(48)
      endif
c
      if(imembr.ne.0) call md_membrane_forces(int_mb(i_mm),
     + dbl_mb(i_fm),dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_fs),
     + dbl_mb(i_wws))
c
c     time step
c
      call timer_start(49)
      call cf_mdstep(int_mb(i_iw+(lwdyn-1)*mwm),dbl_mb(i_xw),
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
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is),dxmax)
      call timer_stop(50)
c
c     recalculate velocities
c
      call cf_veloc(nwmloc,dbl_mb(i_xw),dbl_mb(i_yw),dbl_mb(i_vw),
     + nsaloc,dbl_mb(i_xs),dbl_mb(i_ys),dbl_mb(i_vs))
c
c     velocity scaling to preset temperature
c
      if(frequency(mdstep,nfgaus)) then
      call cf_vscale(tgauss,nwmloc,nsaloc,
     + dbl_mb(i_vw),dbl_mb(i_vwt),dbl_mb(i_vs),dbl_mb(i_vst),
     + int_mb(i_iw+(lwdyn-1)*mwm),
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsatt-1)*msa))
      if(iand(ivopt,1).eq.1) nfgaus=0
      endif
c
c     coordinate scaling
c
      call timer_start(51)
      call cf_final(dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_yw),
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
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm),dbl_mb(i_gsm))
c
c     add center of mass coordinates from reference coordinates
c
      if(idifco.gt.0)
     + call md_addref(.true.,dbl_mb(i_xwm),dbl_mb(i_xwcr),
     + dbl_mb(i_xsm),dbl_mb(i_xscr),dbl_mb(i_dsr))
c
c     center solute in box
c
      if(frequency(mdstep,nfcntr)) then
      call cf_center(dbl_mb(i_xw),nwmloc,int_mb(i_is+(lsfrc-1)*msa),
     + dbl_mb(i_xs),nsaloc,idscb,nscb,icentr)
      endif
c
c     remove overall translational motion
c
      if(frequency(mdstep,nfslow)) then
      call cf_slow(dbl_mb(i_xw),dbl_mb(i_vw),nwmloc,dbl_mb(i_xs),
     + dbl_mb(i_vs),int_mb(i_is+(lsatt-1)*msa),nsaloc)
      endif
c
      endif
      call timer_stop(52)
      call timer_start(53)
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
      call timer_stop(53)
      call timer_stop(203)
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
    3 continue
      if(nsaloc.gt.0) then
      call sp_nbbl(nbbl)
      do 4 ibbl=1,nbbl
      call sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
      if(lself) call cf_shakes(lhandl,xs,ys,isl(1,lsgan),isl(1,lsatt),
     + isl(1,lssgm),isl(1,lsdyn),isl(1,lshop),isfr,isto)
    4 continue
      do 5 i=1,nsaloc
      dmax=max(dmax,(xs(i,1)-ys(i,1))**2+(xs(i,2)-ys(i,2))**2+
     + (xs(i,3)-ys(i,3))**2)
    5 continue
      endif
      if(npmf.eq.1.or.lpmfc) then
      if(.not.cf_shakep(xs,ys,isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),
     + isl(1,lshop),nsaloc)) goto 3
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
      subroutine md_finit(iwl,isl,xw,xwm,xs,fw,fs,xsm,xsmp)
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
      call cf_init(stime,lpair,llong,box,vlat,vlati,zw,dbl_mb(i_zs),
     + eww,dbl_mb(i_esw),dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esa))
c
      if(lpair) call md_zinit(int_mb(i_iwz),int_mb(i_isz))
c
      call timer_stop(7)
c
      if(lqmd) then
      call timer_start(8)
      call qmd_forces(irtdb,int_mb(i_is+(lsgan-1)*msa),
     + int_mb(i_is+(lsatt-1)*msa),dbl_mb(i_xs),dbl_mb(i_fs),msa,
     + nsaloc,uqmd,lesp)
      uqmd=uqmd-uqmatm
      call timer_stop(8)
      return
      endif
c
      call timer_start(9)
      call ga_sync()
      call timer_stop(9)
c
      call timer_start(10)
c
      call sp_initf(fw,fs,llong,int_mb(i_iwz),int_mb(i_isz),lpair)
c
      call sp_putix(me,iwl,xw,nwmloc,isl,xs,nsaloc)
c
      call timer_stop(10)
c
      call timer_start(11)
      call ga_sync()
      call timer_stop(11)
c
      if(ncoll.gt.0) then
      call timer_start(12)
      call cf_collapse(ncoll,fcoll,nsaloc,nwmloc,isl(1,lsmol),
     + isl(1,lssgm),
     + dbl_mb(i_xs),dbl_mb(i_xsm),mst,dbl_mb(i_tsm),dbl_mb(i_fs),
     + dbl_mb(i_xw),dbl_mb(i_xwm),dbl_mb(i_fw))
      call timer_stop(12)
      endif
      if(ifield.gt.0) then
      call timer_start(13)
      call cf_extern(stime,nsaloc,dbl_mb(i_fs),
     + isl(1,lsct1),nwmloc,dbl_mb(i_fw))
      call timer_stop(13)
      endif
      call cf_multi(nsaloc,dbl_mb(i_xs),dbl_mb(i_fs),isl(1,lsgan),
     + isl(1,lsatt),isl(1,lsfrc),isl(1,lsdyn),isl(1,lsct1),
     + dbl_mb(i_ess),dbl_mb(i_fss),lfnpmf,lprpmf,iprpmf,
     + npmf.eq.1.or.lpmfc)
c
      if(ipolt.gt.0) then
      call timer_start(23)
      call md_induce(iwl,isl,xw,xwm,xs,dbl_mb(i_pw),
     + dbl_mb(i_pwp),dbl_mb(i_ps),dbl_mb(i_psp))
      call timer_stop(23)
      endif
c
      if(lpme.and.llong) then
      if(lpert2) then
      call pme_energy(2,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(2))
      call timer_start(24)
      call pme_init()
      call timer_stop(24)
      endif
      if(lpert3) then
      call pme_energy(3,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(3))
      call timer_start(24)
      call pme_init()
      call timer_stop(24)
      endif
      call pme_chgrid(iset,dbl_mb(i_xw),nwmloc,dbl_mb(i_xs),
     + isl(1,lsct1),isl(1,lssgm),nsaloc,epme(iset))
      endif
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
      call timer_start(44)
      call timer_start(202)
      call timer_start(204)
      call ga_sync()
      call timer_stop(204)
      call timer_stop(202)
      call timer_stop(44)
c
      call timer_start(45)
      call sp_final(fw,fs,lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(45)

      if(lqmmm) then
      call timer_start(46)
      call qmmm_forces(irtdb,mwm,nwmloc,mwa,nwa,int_mb(i_iwz),xw,fw,
     + msa,nsaloc,int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsdyn-1)*msa),
     + int_mb(i_is+(lsct1-1)*msa),int_mb(i_isz),xs,fs,uqmmm)
c      uqmmm=uqmmm-uqmatm
      call timer_stop(46)
      endif

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

      call timer_start(47)
      call cf_fnorm(iwl(1,lwdyn),fw,nwmloc,
     + isl(1,lsatt),isl(1,lsdyn),fs,nsaloc,fnorm,fmax)
      call timer_stop(47)
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
      call sp_copyg(dbl_mb(i_gw),dbl_mb(i_gs))
      call ga_sync()
c
      return
      end
      subroutine md_cmopt(vs,fs,fcm,ismol,isatt)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      real*8 vs(msa,3),fs(msa,3,2),fcm(msm,5)
      integer ismol(msa),isatt(msa)
c
      call cf_scmfor(icmopt,ismol,isatt,vs,fs,nsaloc,fcm)
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
      logical sp_local,cf_hopping
      external sp_local,cf_hopping
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
      call sp_nbbl(nbbl)
      else
      itforc=35
      call timer_start(33)
      if(ipolt.eq.0) then
      if(nbget.lt.0) then
      call sp_prefetch_all(nbbl,iwl,xw,isl,xs)
      else
      call sp_prefetch(nbbl,iwl,xw,isl,xs)
      endif
      else
      call sp_prefetch_p(nbbl,iwl,xw,dbl_mb(i_pw),dbl_mb(i_pwp),
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
      if(local.and..not.sp_local(ibbl)) then
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
      call sp_getxbl(ibbl,lhandl,
     + iwl,xw,iwfr,iwto,jwfr,jwto,isl,xs,isfr,isto,jsfr,jsto,
     + lself,lpbcs)
      else
      call sp_getxpbl(ibbl,lhandl,
     + iwl,xw,dbl_mb(i_pw),dbl_mb(i_pwp),iwfr,iwto,jwfr,jwto,
     + isl,xs,dbl_mb(i_ps),dbl_mb(i_psp),isfr,isto,jsfr,jsto,
     + lself,lpbcs)
      endif
      call timer_stop(itcom1)
      else
      call timer_start(34)
      call sp_nbwait(ibbl,lnew,lhandl,lself,lpbcs,
     + iwfr,iwto,jwfr,jwto,isfr,isto,jsfr,jsto,iwl,isl)
      call timer_stop(34)
      if(nbget.gt.0.and.lnew) then
      call timer_start(33)
      call sp_prefetch_next(iwl,xw,isl,xs)
      call timer_stop(33)
      endif
      endif
c
      call timer_start(itforc)
      call cf_comw(xw,xwm,jwfr,jwto)
c
      if(ipolt.eq.0) then
      call forces(lself,lpbcs,xw,xwm,fw,zw,dbl_mb(i_rtos),iwl(1,lwdyn),
     + int_mb(i_iwz),iwfr,iwto,jwfr,jwto,xs,dbl_mb(i_xsm),fs,
     + dbl_mb(i_zs),isl(1,lsgan),isl(1,lsatt),isl(1,lsdyn),isl(1,lsgrp),
     + isl(1,lsfrc),isl(1,lsmol),isl(1,lssss),isl(1,lsct1),isl(1,lsct2),
     + isl(1,lsct3),isl(1,lssgm),isl(1,lshop),int_mb(i_isz),
     + isfr,isto,jsfr,jsto,lpbc,lhandl,
     + .true.,eww,dbl_mb(i_esw),dbl_mb(i_ess),dbl_mb(i_fss),
     + dbl_mb(i_esa),int_mb(i_lseq),lforce)
      else
      call forcep(lself,lpbcs,xw,xwm,fw,dbl_mb(i_pw),dbl_mb(i_pwp),zw,
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
      call sp_accfbl(ibbl,lhandl,fw,fs,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(itcom2)
      else
      call timer_start(37)
      call sp_nbaccfbl(ibbl,lhandl,fw,fs,
     + lpair,int_mb(i_iwz),int_mb(i_isz))
      call timer_stop(37)
      endif
c
    1 continue
c
      if(nbget.ne.0) then
      call timer_start(38)
      call sp_nbwaitf()
      call timer_stop(38)
      endif
c
      if(.not.ldone) then
      if(cf_hopping(lpbc,lpbcs,stime,isl,isl(1,lssgm),isl(1,lsgan),
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
      call sp_update_i(nsaloc,int_mb(i_is),nwmloc,int_mb(i_iw))
      endif
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      goto 100
      endif
      endif
c
      if(lpme.and.llong) then
      call pme_forces(fw(1,1,1,2),nwmloc,fs(1,1,2),
     + isl(1,lsct1),isl(1,lssgm),nsaloc)
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
      logical lself,lpbcs
      integer ibbl,nbbl,lhandl
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
      npolit=npolit+1
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
      call sp_getxpbl(ibbl,lhandl,
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
      lpair=.false.
      lload=.false.
c
      if(np.gt.0) call ga_sync()
c
      if(nwmloc.gt.0) then
      do 26 k=1,nwa
      do 27 j=1,3
      do 28 i=1,nwmloc
      pw(i,j,k,2)=pw(i,j,k,1)
      pw(i,j,k,1)=zero
   28 continue
   27 continue
   26 continue
      endif
      if(nsaloc.gt.0) then
      do 29 j=1,3
      do 30 i=1,nsaloc
      ps(i,j,2)=ps(i,j,1)
      ps(i,j,1)=zero
   30 continue
   29 continue
      endif
c
      call sp_getp(me,pw,pwp,nwmloc,ps,psp,nsaloc,lpert2.or.lpert3,1)
c
      if(nwmloc.gt.0) then
      do 31 k=1,nwa
      do 32 j=1,3
      do 33 i=1,nwmloc
      pw(i,j,k,1)=pw(i,j,k,1)+pw(i,j,k,2)
      pw(i,j,k,2)=zero
   33 continue
   32 continue
   31 continue
      endif
      if(nsaloc.gt.0) then
      do 34 j=1,3
      do 35 i=1,nsaloc
      ps(i,j,1)=ps(i,j,1)+ps(i,j,2)
      ps(i,j,2)=zero
   35 continue
   34 continue
      endif
c
      call sp_getp(me,pw,pwp,nwmloc,ps,psp,nsaloc,lpert2.or.lpert3,2)
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
      if(np.gt.1) call ga_dgop(mrg_d06,pmax,1,'max')
c
      if(pmax.gt.ptol.and.npolit.le.mpolit.and.ipolt.gt.1) goto 11
c
c     copy current fields into local global array
c     -------------------------------------------
c
      call sp_putp(me,pw,pwp,nwmloc,ps,psp,nsaloc,lpert2.or.lpert3)
c
      return
      end
      subroutine md_wrtrst(lfn,fil,lveloc)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "util.fh"
c
      real*8 timer_wall
      external timer_wall
c
      integer lfn
      character*255 fil
      logical lveloc
      integer i,left
      character*255 filn
c
      if(me.eq.0) then
      if(keepr.eq.0) then
      open(unit=lfn,file=fil(1:index(fil,' ')-1),
     + form='formatted',status='unknown',err=9999)
      else
      if(ntype.eq.3.and.npg.le.1) then
      print*,npg,ntype,irun
      print*,fil(1:index(fil,' ')-1)
      print*,fil(1:index(fil,'.rst')-1)
      write(filn,'(a,i5.5,a)') fil(1:index(fil,'.rst')-1),
     + irun,'.rst '
      else
      write(filn,'(a,a,i5.5,a)') fil(1:index(fil,'.rst')-1),'-',
     + keepr,'.rst '
      endif
      open(unit=lfn,file=filn(1:index(filn,' ')-1),
     + form='formatted',status='unknown',err=9999)
      endif
      endif
c
      call sp_wrtrst(lfn,fil,lveloc,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_gw),
     + dbl_mb(i_xwcr),int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),
     + dbl_mb(i_gs),dbl_mb(i_xscr),projct,int_mb(i_lseq))
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
      call timer_stop(205)
      tneed=timer_wall(205)
      call timer_reset(205)
      call timer_start(205)
      left=util_time_remaining(irtdb)
      tleft=dble(left)
      i=0
      if(left.lt.0) then
      i=1
      elseif(tleft.gt.two*tneed) then
      i=1
      endif
      call ga_brdcst(mrg_d48,i,ma_sizeof(mt_int,1,mt_byte),0)
      lstop=i.eq.0
c
      return
c
 9999 continue
      call md_abort('Unable to open restart for writing',me)
      return
      end
      subroutine md_server
c
      implicit none
c
#include "md_common.fh"
c
#if !defined(WIN32)
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
      lens=66
      string(66:66)=char(13)
c
      numbyt=client_socket_write(isocket,string,lens)
c
      print*,'Bytes written to socket is ',numbyt
c
#endif
      return
      end
      subroutine md_timer_init()
c
      implicit none
c
#include "md_common.fh"
c
      integer i
c
      do 1 i=1,200
      call timer_reset(i)
    1 continue
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
      external timer_wall,timer_wall_total
      real*8 timer_wall,timer_wall_total
c
      integer i,j
      real*8 tim(56,1024)
c
      do 1 i=1,np
      do 2 j=1,56
      tim(j,i)=zero
    2 continue
    1 continue
c
      do 3 i=1,55
      tim(i,me+1)=timer_wall_total(i)
      tim(56,me+1)=tim(56,me+1)+tim(i,me+1)
    3 continue
c
      call ga_dgop(mrg_d04,tim,56*np,'+')
c
      if(me.eq.0) then
      write(lfntim,1000) stime
 1000 format('timings',/,f12.6)
      do 1002 j=1,np
      write(lfntim,1001) (tim(i,j),i=1,56)
 1001 format(10f7.3)
 1002 continue
      endif
c
      return
      end
      subroutine md_test
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
      subroutine md_fd(isg,xs,fst,fs,iwg,xw,fwt,fw)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer isg(msa),iwg(mwm)
      real*8 xs(msa,3),fst(msa,3),fs(msa,3)
      real*8 xw(mwm,3,mwa),fwt(mwm,3,mwa),fw(mwm,3,mwa)
c
      integer i,j,k,lsa
      real*8 etott,xsk(3),ft(3,3),dx,xdv
c
      if(me.eq.0) write(lfnout,1000)
 1000 format(//,' FINITE DIFFERENCE SOLUTE FORCES',//,
     + '  Atom',
     + '                 Analytic forces       ',
     + '             Finite difference forces  ',
     + '             Deviation                 ',/,
     + '      ',
     + '             fx          fy          fz',
     + '             fx          fy          fz',
     + '            dfx         dfy         dfz',/)
c
      lpair=.false.
      lload=.false.
c      xdv=0.000001
      xdv=dx0sd
c
      do 1 i=1,nsaloc
      fst(i,1)=fs(i,1)
      fst(i,2)=fs(i,2)
      fst(i,3)=fs(i,3)
    1 continue
      do 2 i=1,nwmloc
      do 3 j=1,nwa
      fwt(i,1,j)=fw(i,1,j)
      fwt(i,2,j)=fw(i,2,j)
      fwt(i,3,j)=fw(i,3,j)
    3 continue
    2 continue
      etott=etot
c
      do 4 i=1,nsa
      lsa=0
      ft(1,1)=zero
      ft(2,1)=zero
      ft(3,1)=zero
      ft(1,2)=zero
      ft(2,2)=zero
      ft(3,2)=zero
      ft(1,3)=zero
      ft(2,3)=zero
      ft(3,3)=zero
      do 5 j=1,nsaloc
      if(isg(j).eq.i) then
      lsa=j
      xsk(1)=xs(j,1)
      xsk(2)=xs(j,2)
      xsk(3)=xs(j,3)
      ft(1,1)=fst(j,1)
      ft(2,1)=fst(j,2)
      ft(3,1)=fst(j,3)
      goto 6
      endif
    5 continue
    6 continue
c
      do 7 j=2,3
      do 8 k=1,3
c
      dx=-xdv
      if(j.eq.3) dx=xdv
c
      if(lsa.gt.0) xs(lsa,k)=xs(lsa,k)+dx
c
c     atomic forces and potential energies
c
      call md_finit(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs),
     + dbl_mb(i_xsm),dbl_mb(i_xsmp))
      call md_forces(int_mb(i_iw),int_mb(i_is),dbl_mb(i_xw),
     + dbl_mb(i_xwm),dbl_mb(i_xs),dbl_mb(i_fw),dbl_mb(i_fs))
c
      call prp_proper(0,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,uqmd,uqmmm,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsm))
c
c
      if(lsa.gt.0) then
      xs(lsa,k)=xsk(k)
      ft(k,j)=etott-etot
      endif
c
    8 continue
    7 continue
c
      if(np.gt.1) call ga_dgop(mrg_d14,ft,6,'+')
c
      if(me.eq.0) write(lfnout,1001) i,(ft(j,1),j=1,3),
     + ((ft(j,3)-ft(j,2))/(two*xdv),j=1,3),
     + (ft(j,1)-(ft(j,3)-ft(j,2))/(two*xdv),j=1,3)
 1001 format(i7,5x,3f12.3,3x,3f12.3,3x,3E12.3)
c
    4 continue
c
      return
      end
      subroutine md_membrane_init(ismol,mm,xs,xsm,fm)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer ismol(msa),mm(msa,2)
      real*8 xs(msa,3),xsm(msm,3),fm(msm,7)
c
      integer i,j,k,numg
      real*8 dx,dk
c
      numg=0
c
      do 1 i=1,msm
      mm(i,1)=0
    1 continue
      do 2 i=1,nsaloc
      mm(ismol(i),1)=mm(ismol(i),1)+1
      mm(i,2)=ismol(i)
    2 continue
      if(np.gt.1) call ga_dgop(mrg_d49,mm,msm,'+')
      do 3 i=1,nsaloc
      if(mm(ismol(i),1).eq.1) then
      k=0
      dk=zero
      do 4 j=1,nsm
      if(ismol(i).ne.j) then
      dx=(xs(i,1)-xsm(j,1))**2+(xs(i,2)-xsm(j,2))**2+
     + (xs(i,3)-xsm(j,3))**2
      if(k.eq.0.or.dx.lt.dk) then
      dk=dx
      k=j
      endif
      endif
    4 continue
      mm(i,2)=k
      numg=numg+1
      endif
    3 continue
      if(numg.gt.0.and.me.eq.0) then
      write(*,2000) numg
 2000 format(' Regrouping of',i5,' atoms',/)
      endif
c
      return
      end
      subroutine md_membrane_forces(mm,fm,xs,xsm,fs,ws)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer mm(msa,2)
      real*8 fm(msm,7),xs(msa,3),xsm(msm,3),fs(msa,3),ws(msa)
c
      integer i,j
      real*8 factor
c
      do 1 i=1,msm
      fm(i,1)=zero
      fm(i,2)=zero
      fm(i,3)=zero
      fm(i,4)=zero
      fm(i,5)=zero
      fm(i,6)=zero
      fm(i,7)=zero
    1 continue
c
      do 2 i=1,nsaloc
      factor=one/ws(i)
      fm(mm(i,2),1)=fm(mm(i,2),1)+factor*fs(i,1)
      fm(mm(i,2),2)=fm(mm(i,2),2)+factor*fs(i,2)
      fm(mm(i,2),3)=fm(mm(i,2),3)+factor*fs(i,3)
      fm(mm(i,2),4)=fm(mm(i,2),4)+factor*
     + ((xs(i,1)-xsm(mm(i,2),1))*fs(i,2)-
     +  (xs(i,2)-xsm(mm(i,2),2))*fs(i,1))
    2 continue
      if(np.gt.1) call ga_dgop(mrg_d50,fm,4*msm,'+')
c
      do 3 i=1,nsm
      fm(i,4)=fm(i,3)
      if(me.eq.0) write(lfnout,1000) i,(fm(i,j)/dble(mm(i,1)),j=1,3)
 1000 format(i5,3f12.3)
    3 continue
c
c     molecular rotations only
c
      if(imembr.eq.2) then
      do 4 i=1,nsaloc
      fs(i,1)=fs(i,1)-ws(i)*fm(mm(i,2),1)/dble(mm(i,1))
      fs(i,2)=fs(i,2)-ws(i)*fm(mm(i,2),2)/dble(mm(i,1))
      fs(i,3)=fs(i,3)-ws(i)*fm(mm(i,2),3)/dble(mm(i,1))
    4 continue
      endif
c
      return
      end
