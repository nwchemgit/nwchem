      subroutine grid_setspac(xyz, xyzw, Rij, rq, p, ictr, nq,
     ,     expo, rad)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c
#include "cdft.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
#include "dftps.fh"
c
      integer ictr
      integer nq
      double precision xyzw(4,nq), xyz(3,ncenters)
      double precision Rij((ncenters*(ncenters+1))/2)
      double precision rq(nq,ncenters)
      double precision p(nq,ncenters)
      double precision rad
      double precision expo(*)
c
      integer iscr,lscr,lgth_scr,ixyzm,lxyzm
      integer ncenters_scr ! no. of centers after screening
      integer ictr_scr  ! grid center in the screened reduced list
      integer k_indx,l_indx,i
      integer basis
      integer ga_nodeid
      external ga_nodeid
c
c
c     set scratch space for grid compression
c
      lgth_scr=nq
      if(lssw) lgth_scr=max(lgth_scr,ncenters)
      if (.not.MA_Push_get(mt_dbl,3*ncenters,' xyzm',lxyzm,ixyzm))
     &   call errquit('grid_setspac: cannot allocate xyzm',1, MA_ERR)
      if (.not.MA_Push_get(mt_int,lgth_scr,' scr gridcmp',lscr,iscr))
     &   call errquit('grid_setspac: cannot allocate lscr',2, MA_ERR)
      call dcopy(3*ncenters,xyz,1,dbl_mb(ixyzm),1)

c
c     should have all possible spatial weight possibilities here ... for now 
c     try Becke weights.
c
      ncenters_scr=ncenters
      ictr_scr=ictr
      if (.not.MA_Push_get(mt_int,ncenters,' indx',l_indx,k_indx))
     &   call errquit('grid_setspac: cannot allocate indx',1, MA_ERR)
      do i=1,ncenters
         int_mb(k_indx+i-1)=i
      enddo
      if(.not.lssw) then
c
c     becke
c
         if(trunc_neigh) 
     .        call grid_neighb_list(xyz(1,ictr), dbl_mb(ixyzm), 
     ,        ncenters_scr, ictr_scr,int_mb(k_indx))
         call grid_beckew(xyzw, dbl_mb(ixyzm),  rq, p, ictr_scr,
     ,        nq,int_mb(iscr),ncenters_scr,int_mb(k_indx))
      else
c
c     call to ssw routine
c
         if (odftps) call pstat_on(ps_ssw)
         basis=ao_bas_han
         if (adft) basis=cd_bas_han
         if(trunc_neigh)
     .        call grid_signf(basis,  ncenters_scr, 
     &        xyz,dbl_mb(ixyzm),ictr_scr,
     &        rad, expo, int_mb(k_indx),iatype_pt_chg,
     &        nq,xyzw)
         if(whichssw.eq.'erf1') then
            call grid_ssw(xyzw, dbl_mb(ixyzm), rij,rq, 
     &           p, ictr_scr, ncenters_scr, nq,
     ,           int_mb(iscr),int_mb(k_indx), 'erf1')
         else
            call grid_sswold(xyzw, dbl_mb(ixyzm), rij,rq, 
     &           p, ictr_scr, ncenters_scr, nq,
     ,           int_mb(iscr),int_mb(k_indx),whichssw)
         endif
         if (odftps) call pstat_off(ps_ssw)
      endif
         


c
      if (.not. MA_chop_Stack(lxyzm))
     &   call errquit('grid_setspac: pop stack failed.',1, MA_ERR)
      return
      end 
