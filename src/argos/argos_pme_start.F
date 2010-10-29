      subroutine argos_pme_start(a,m,imffti,nodpmi,
     + ngxi,ngyi,ngzi,nwm,nwa,nsa,ictrl,nbgeti)
c $Id$
      implicit none
c
      integer nodpmi,ngxi,ngyi,ngzi,imffti,nwm,nwa,nsa,m
      integer ictrl,nbgeti
      real*8 a
c
      integer mem
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      me=ga_nodeid()
      np=ga_nnodes()
c
      pi=four*atan(one)
      twopi=two*pi
      alpha=a
      morder=m
      nodpme=nodpmi
      icntrl=ictrl
      nbget=nbgeti
c
      imfft=imffti
      lcorr(1)=.false.
      lcorr(2)=.false.
      lcorr(3)=.false.
c
      mwa=nwa
      mwm=nwm
      msa=nsa
c
      ngx=ngxi
      ngy=ngyi
      ngz=ngzi
      ngmax=max(ngx,ngy,ngz)
      mgx=ngx
      mgy=ngy
      mgz=(ngz/nodpme)+1
      ngrx=ngx+morder
      ngry=ngy+morder
      ngrz=ngz
      mgrx=ngrx
      mgry=ngry
      mgrz=mgz
      if(nodpme*mgz.lt.ngz)
     + call md_abort('Inconsistent number of pme nodes',0)
c
c     allocate memory for spline coefficients
c
      if(.not.ma_push_get(mt_dbl,3*ngmax,'bmod',l_bmod,i_bmod))
     + call md_abort('Failed to allocate bmod',0)
      if(.not.ma_push_get(mt_dbl,ngmax,'barr',l_barr,i_barr))
     + call md_abort('Failed to allocate barr',0)
c
c     allocate memory for the grids
c
      mem=2*mgx*mgy*mgz
      if(.not.ma_push_get(mt_dbl,mem,'grid1',l_grd1,i_grd1))
     + call md_abort('Failed to allocate grid1',0)
      if(.not.ma_push_get(mt_dbl,mem,'grid2',l_grd2,i_grd2))
     + call md_abort('Failed to allocate grid2',0)
      if(.not.ma_push_get(mt_dbl,mgrx*mgry*mgz,'grid',l_grd,i_grd))
     + call md_abort('Failed to allocate grid',0)
      if(nbget.gt.0) then
      if(.not.ma_push_get(mt_dbl,mgrx*mgry*mgz,'gridt',l_grdt,i_grdt))
     + call md_abort('Failed to allocate gridt',0)
      if(.not.ma_push_get(mt_int,mgz,'gridh',l_grdh,i_grdh))
     + call md_abort('Failed to allocate gridh',0)
      else
      if(.not.ma_push_get(mt_dbl,1,'gridt',l_grdt,i_grdt))
     + call md_abort('Failed to allocate gridt',0)
      if(.not.ma_push_get(mt_int,1,'gridh',l_grdh,i_grdh))
     + call md_abort('Failed to allocate gridh',0)
      endif
c
c     allocate memory for node arrays
c
      if(.not.ma_push_get(mt_int,np,'lnode',l_lnod,i_lnod))
     + call md_abort('Failed to allocate lnode',0)
      if(.not.ma_push_get(mt_int,np,'kfrom',l_kfr,i_kfr))
     + call md_abort('Failed to allocate kfrom',0)
      if(.not.ma_push_get(mt_int,np,'kto',l_kto,i_kto))
     + call md_abort('Failed to allocate kto',0)
      if(.not.ma_push_get(mt_int,np,'nodep',l_nodp,i_nodp))
     + call md_abort('Failed to allocate nodep',0)
c
      if(.not.ma_push_get(mt_int,ngz,'knode',l_knod,i_knod))
     + call md_abort('Failed to allocate knode',0)
c
      if(.not.ma_push_get(mt_log,mgz*np,'lsldo',l_lsld,i_lsld))
     + call md_abort('Failed to allocate ldldo',0)
      if(.not.ma_push_get(mt_int,4*mgz*np,'slmax',l_slmax,i_slmax))
     + call md_abort('Failed to allocate slmax',0)
c
c     allocate memory for fft node maps
c
      mem=3*max(ngx,ngy,ngz)
      if(.not.ma_push_get(mt_int,mem,'map1',l_map1,i_map1))
     + call md_abort('Failed to allocate map1',0)
      if(.not.ma_push_get(mt_int,mem,'map2',l_map2,i_map2))
     + call md_abort('Failed to allocate map2',0)
c
c     allocate memory for slab arrays
c
      mem=mgz*np
      if(.not.ma_push_get(mt_int,mem,'slab',l_slab,i_slab))
     + call md_abort('Failed to allocate slab',0)
c
c     create the global arrays for the slabs
c
      call ga_create_list(np*mgrx,mgry,'sl',mgrx,mgry,mgz,lslab)
c
c     create the flag array
c
      if(.not.ga_create(mt_int,np,5,'flg',1,5,ga_flg))
     + call md_abort('Failed to allocate global array flg',0)
      call argos_pme_flag(0,0,me)
      call argos_pme_flag(1,0,me)
c
c     select the fft
c
      call argos_pme_select(zero)
c
c     test the fft
c
      if(lpnode)
     + call argos_pme_fftest(dbl_mb(i_grd1),dbl_mb(i_grd2),
     + int_mb(i_nodp))
c
      call argos_pme_coeff(dbl_mb(i_bmod),dbl_mb(i_barr))
c
      return
      end
