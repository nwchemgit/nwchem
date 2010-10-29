      subroutine argos_diana_rmsd
c
c $Id$
      implicit none
c
#include "argos_diana_params.fh"
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer i_rms,l_rms,ga_rms
      integer ilr,ihr,jlr,jhr,ndxr,ldr
c
      filrms=card(14:80)
c
      if(me.eq.0) then
      if(.not.ma_push_get(mt_dbl,nfrtot*3,'trms',l_rms,i_rms))
     + call md_abort('Failed to allocate trms',0)
      endif
c
      call ga_distribution(ga_trj,me,il,ih,jl,jh)
      call ga_access(ga_trj,il,ih,jl,jh,ndx,ld)
c
c     allocate global array to hold the rms data
c
      if(.not.ga_create(mt_dbl,2,nfrtot,'rms',2,nfrdim,ga_rms))
     + call md_abort('Failed to create ga_rms',0)
      call ga_distribution(ga_rms,me,ilr,ihr,jlr,jhr)
      call ga_access(ga_rms,ilr,ihr,jlr,jhr,ndxr,ldr)
c
      if(card(8:13).eq.'averag') then
      elseif(card(8:13).eq.'refere') then
      call argos_diana_rmsdev(dbl_mb(ndx),dbl_mb(i_xref),dbl_mb(i_swt),
     + dbl_mb(i_wt),dbl_mb(i_xrms),dbl_mb(ndxr),dbl_mb(i_rms),ga_rms,
     + int_mb(i_isel),byte_mb(i_snam))
      else
      call md_abort('Error in rmsdev directive',0)
      endif
c
      if(.not.ga_destroy(ga_rms))
     + call md_abort('Could not destroy ga_rms',0)
c
      if(me.eq.0) then
      if(.not.ma_pop_stack(l_rms))
     + call md_abort('Failed to deallocate trms',0)
      endif
c
      return
      end
