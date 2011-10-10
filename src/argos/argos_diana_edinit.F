      subroutine argos_diana_edinit()
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      if(me.eq.0) then
      write(*,1000) (ma_sizeof(mt_dbl,3*nsel,mt_byte)+
     + ma_sizeof(mt_int,6*nsel,mt_byte))/104857,
     + ma_inquire_avail(mt_byte)/1048576,
     + ma_sizeof(mt_dbl,18*nsel*nsel/np,mt_byte)/1048576,
     + ga_memory_avail()/1048576
 1000 format(/,' Memory   Required    Available',//,
     + '   ma   ',2i10,' MB',/,
     + '   ga   ',2i10,' MB',/)
      endif
c      print*,'nsel=',nsel
c
c     create nsel x nsel global array covariance matrix
c
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'covar',0,3*nsel,ga_cov))
     + call md_abort('Failed to create global array covar',me)
      call ga_zero(ga_cov)
      if(.not.ga_create(mt_dbl,3*nsel,3*nsel,'vector',0,3*nsel,ga_vec))
     + call md_abort('Failed to create global array vector',me)
c
c     allocate memory
c
      if(.not.ma_push_get(mt_int,3*nsel,'ndx',l_ndx,i_ndx))
     + call md_abort('Failed to allocate memory for ndx',me)
      if(.not.ma_push_get(mt_dbl,3*nsel,'cov',l_cov,i_cov))
     + call md_abort('Failed to allocate memory for cov',me)
      if(.not.ma_push_get(mt_int,3*nsel,'ord',l_ord,i_ord))
     + call md_abort('Failed to allocate memory for ord',me)
c
      if(me.eq.0) then
      write(*,1002) 3*nsel
 1002 format(' Memory allocated for covariance ',i5)
      endif
c
      return
      end
c $Id$
