      subroutine argos_diana_edfinal()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
c     deallocate memory
c
      if(.not.ma_pop_stack(l_ord))
     + call md_abort('Failed to deallocate memory for ord',me)
      if(.not.ma_pop_stack(l_cov))
     + call md_abort('Failed to deallocate memory for cov',me)
      if(.not.ma_pop_stack(l_ndx))
     + call md_abort('Failed to deallocate memory for ndx',me)
c
c     destroy global array covariance matrix
c
      if(.not.ga_destroy(ga_vec))
     + call md_abort('Failed to destroy global array vector',me)
      if(.not.ga_destroy(ga_cov))
     + call md_abort('Failed to destroy global array covar',me)
c
      return
      end
c $Id$
