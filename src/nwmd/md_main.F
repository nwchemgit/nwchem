      subroutine md_main()
c
      implicit none
c
#include "md_common.fh"
c
c     single energy
c     -------------
c
c     energy minimization
c     -------------------
c
c     molecular dynamics
c     ------------------
c
      call md_md()
c
c     free energy simulation
c     ----------------------
c
      return
      end
      subroutine md_md()
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
c     equilibration
c     -------------
c
      do 1 iequi=kequi+1,mequi
      if(.not.ma_verify_allocator_stuff()) stop '0000'
      call md_newton()
c
      call prp_step(mdstep,temp,pres,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess))
c
    1 continue
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
      mdstep=mdstep+1
c
c     reassign velocities
c
c     if(mod(mdstep,nfgaus).eq.0)
c    + call cf_gauss(nwmloc,nsaloc,vw,vs,tgaus)
c
      llong=(mdstep.eq.1.or.mod(mdstep,nflong).eq.0).and.ltwin
      lpair=mdstep.eq.1.or.mod(mdstep,nfpair).eq.0
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
     + dbl_mb(i_fs),nsaloc)
c
c     shake
c
      call md_shake(dbl_mb(i_xw),dbl_mb(i_yw),int_mb(i_iw),
     + dbl_mb(i_xs),dbl_mb(i_ys),int_mb(i_is))
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
     + int_mb(i_is+(lsdyn-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_zs),nsaloc,box,pres,prsscl,temp)
c
      return
      end
      subroutine md_shake(xw,yw,iwl,xs,ys,isl)
c
      implicit none
c
#include "md_common.fh"
#include "mafdecls.fh"
c
      real*8 xw(mwm,3,mwa),yw(mwm,3,mwa),xs(msa,3),ys(msa,3)
      integer iwl(mwm,miw2),isl(msa,mis2)
c
      integer lhandl,ibbl,iwfr,iwto,isfr,isto
      logical lself
c
      if(nwmloc.gt.0)
     + call cf_shakew(xw,yw,iwl(1,lwgmn),iwl(1,lwdyn),nwmloc)
c
      if(nsaloc.gt.0) then
      call sp_nbbl(nbbl)
      do 1 ibbl=1,nbbl
c
      call sp_gethdl(ibbl,lhandl,lself,iwfr,iwto,isfr,isto)
c
      if(lself) call cf_shakes(lhandl,xs,ys,isl(1,lsgan),isl(1,lsatt),
     + isl(1,lsdyn),isfr,isto)
c
    1 continue
c
      endif
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
      call ga_sync()
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
c
      call sp_nbbl(nbbl)
c
      call cf_cenmas(nwmloc,dbl_mb(i_xw),dbl_mb(i_xwm),nsaloc,
     + int_mb(i_is+(lsatt-1)*msa),int_mb(i_is+(lsmol-1)*msa),
     + dbl_mb(i_xs),dbl_mb(i_xsm))
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
      return
      end
