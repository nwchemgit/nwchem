      subroutine tce_mrcc_force_orben()
        implicit none
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_mrcc.fh"
#include "tce_main.fh"
#include "geom.fh"

      logical nodezero
      double precision dsumh,dsump
      integer i,iref
      double precision p,h
      integer g6b
      

      nodezero = (ga_nodeid().eq.0)

      do iref=1,nref

         noa = nblcks(1,iref)
         nob = nblcks(2,iref)
         nva = nblcks(3,iref)
         nvb = nblcks(4,iref)

         noab = noa+nob
         nvab = nva+nvb

      dsumh = 0.0d0
      dsump = 0.0d0
      p = 0.0d0
      h = 0.0d0

      DO g6b = 1,noab+nvab
        do i=1,int_mb(k_rangem(iref)+g6b-1)

           h = h + 1.0d0
           dsumh = dsumh + dbl_mb(k_evl_sortedm(iref)+i-1)

        enddo

      orbdegenenergy(iref,g6b) = (dsumh/h)
      if(nodezero)write(6,"('Denom forced: ',4F16.12)")
     1 orbdegenenergy(iref,g6b)

       enddo
      enddo ! iref

      return
      end

c $Id$
