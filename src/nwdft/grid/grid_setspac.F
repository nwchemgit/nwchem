      subroutine grid_setspac(rtdb, xyz, xyzw, Rij, rq, p, ictr, nq,
     ,     lssw, whichf, expo, rad)
c
c$Id: grid_setspac.F,v 1.3 1999-05-01 17:28:11 edo Exp $
c
      implicit none
      integer rtdb
c
#include "cdft.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
c
      integer ictr
      integer nq
      double precision xyzw(4,nq), xyz(3,ncenters)
      double precision Rij((ncenters*(ncenters+1))/2)
      double precision rq(nq,ncenters)
      double precision p(nq,ncenters)
      double precision rad
      double precision expo(*)
      logical lssw
      character*4 whichf
c
      integer iscr,lscr,lgth_scr,ixyzm,lxyzm
      integer ncenters_scr ! no. of centers after screening
      integer ictr_scr  ! grid center in the screened reduced list
c
c     set scratch space for grid compression
c
      lgth_scr=nq
      if(lssw) lgth_scr=max(lgth_scr,ncenters)
      if (.not.MA_Push_get(mt_dbl,3*ncenters,' xyzm',lxyzm,ixyzm))
     &   call errquit('grid_setspac: cannot allocate xyzm',1)
      if (.not.MA_Push_get(mt_int,lgth_scr,' scr gridcmp',lscr,iscr))
     &   call errquit('grid_setspac: cannot allocate lscr',2)
      call dcopy(3*ncenters,xyz,1,dbl_mb(ixyzm),1)

c
c     should have all possible spatial weight possibilities here ... for now 
c     try Becke weights.
c
      ncenters_scr=ncenters
      ictr_scr=ictr
      if(.not.lssw) then
c
c     becke
c
         call grid_neighb_list(xyz(1,ictr), dbl_mb(ixyzm), 
     ,        ncenters_scr, ictr_scr)
         call grid_beckew(xyzw, dbl_mb(ixyzm),  rq, p, ictr_scr,
     ,        nq,int_mb(iscr),ncenters_scr)
      else
c
c     call to ssw routine
c
         call grid_signf(AO_bas_han,  nbf_ao, ncenters_scr, 
     &        xyz,dbl_mb(ixyzm),ictr_scr,
     &        rad, expo, int_mb(iscr))
         call grid_ssw(xyzw, dbl_mb(ixyzm),rij, rq, 
     &        p, ictr_scr, ncenters_scr, nq,
     ,        int_mb(iscr),whichf)
      endif
         


c
      if (.not. MA_Pop_Stack(lscr))
     &   call errquit('grid_setspac: pop stack failed.',2)
      if (.not. MA_Pop_Stack(lxyzm))
     &   call errquit('grid_setspac: pop stack failed.',1)
      return
      end 
