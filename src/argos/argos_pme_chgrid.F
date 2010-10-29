      subroutine argos_pme_chgrid(ipset,xw,nwmloc,xs,isq,
     + isgm,nsaloc,epme)
c $Id$
      implicit none
c
#include "argos_pme_common.fh"
#include "mafdecls.fh"
c
      integer ipset,nwmloc,nsaloc
      real*8 xw(mwm,3,mwa),xs(msa,3),epme
      integer isq(msa),isgm(msa)
c
      call timer_start(24)
c
      call argos_pme_init()
c
      if(.not.lcorr(ipset)) then
      call argos_cafe_pmecor(ipset,xw,nwmloc,xs,isq,isgm,
     + nsaloc,ecorr(ipset))
      lcorr(ipset)=.true.
      endif
      call timer_stop(24)
c
      call argos_pme_grid(ipset,xw,nwmloc,xs,isq,isgm,nsaloc,
     + int_mb(i_nodp),int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + int_mb(i_lnod),log_mb(i_lsld),int_mb(i_slmax),
     + dbl_mb(i_grd1),dbl_mb(i_grd2),dbl_mb(i_grd),
     + int_mb(i_grdh),dbl_mb(i_grdt))
c
      epme=ecorr(ipset)
c
      if(lpnode) then
c
      call timer_start(26)
c
      call argos_pme_fft(ngx,ngy,ngz,-2,dbl_mb(i_grd2),dbl_mb(i_grd1),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(26)
      call timer_start(27)
c
      call argos_cafe_epme(int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + dbl_mb(i_bmod),dbl_mb(i_grd1),mgz,epme)
c
      epme=epme+ecorr(ipset)
c
      call timer_stop(27)
      call timer_start(28)
c
      call argos_pme_sync(nodpme,int_mb(i_nodp))
c
      call timer_stop(28)
      call timer_start(31)
c
      call argos_pme_fft(ngx,ngy,ngz,1,dbl_mb(i_grd1),dbl_mb(i_grd2),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(31)
      call timer_start(32)
c
      call argos_pme_putgrd(dbl_mb(i_grd2),dbl_mb(i_grd))
c
      call timer_stop(32)
c
c     synchronize pme nodes
c
      call timer_start(33)
c
      call argos_pme_sync(nodpme,int_mb(i_nodp))
c
      call timer_stop(33)
c
c     set flag on all nodes to indicate that potential
c     grid calculation has completed
c
      call timer_start(34)
c
      if(me.eq.int_mb(i_nodp)) call argos_pme_flag(0,1,-1)
c
      call timer_stop(34)
      endif
c
      call timer_start(34)
      if(icntrl.eq.1) call ga_sync()
c
      call timer_stop(34)
c
      return
      end
