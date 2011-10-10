      subroutine argos_pme_energy(ipset,xw,nwmloc,xs,
     + isq,isgm,nsaloc,epme)
c
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
      if(.not.lpnode) return
c
      call timer_start(28)
c
      call argos_pme_fft(ngx,ngy,ngz,-2,dbl_mb(i_grd2),dbl_mb(i_grd1),
     + int_mb(i_map1),int_mb(i_map2),int_mb(i_knod),imfft)
c
      call timer_stop(28)
      call timer_start(29)
c
      call argos_cafe_epme(int_mb(i_knod),int_mb(i_kfr),int_mb(i_kto),
     + dbl_mb(i_bmod),dbl_mb(i_grd1),mgz,epme)
c
      epme=epme+ecorr(ipset)
c
      call timer_stop(29)
c
      return
      end
c $Id$
