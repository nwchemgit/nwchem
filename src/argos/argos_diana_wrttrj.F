      subroutine argos_diana_wrttrj(ga_t,fname,fmt,nbatch,logw)
c
c     read the trajectory
c     -------------------
c
      implicit none
c
#include "argos_diana_common.fh"
#include "argos_diana_params.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "msgids.fh"
#include "util.fh"
c
      integer ga_t,nbatch
      character*255 fname
      real*8 x(msa,3)
      logical logw
c
      integer ifil,ibatch
      character*3 fmt
      integer i,ipd,ifr
c
      if(nbatch.le.0) return
      if(me.ne.0) return
c
      ibatch=0
      ifil=0
      ipd=index(fname,'.')
      if(ipd.eq.0) ipd=index(fname,' ')
      lsonly=.true.
c
      do 1 ifr=1,nfrtot
c
      if(ibatch.eq.0) then
      ifil=ifil+1
      write(filcop,'(a,i5.5,a,a)') fname(1:ipd-1),ifil,'.',fmt
      call argos_diana_wthdr(fmt,byte_mb(i_snam),byte_mb(i_tag),
     + int_mb(i_isel),
     + logw)
      endif
      ibatch=ibatch+1
c
c     write frame
c
      call ga_get(ga_t,1,3*nsa,ifr,ifr,dbl_mb(i_xdat),nsa)
      call argos_diana_wtfram(fmt,byte_mb(i_snam),dbl_mb(i_xdat),
     + dbl_mb(i_wdat),
     + int_mb(i_isel),byte_mb(i_tag),dbl_mb(i_val),int_mb(i_wsel),
     + int_mb(i_ndxw))
c
      if(ibatch.eq.nbatch.or.ifr.eq.nfrtot) then
      ibatch=0
      close(unit=lfncop)
      endif
c
c
    1 continue
c
      return
      end
c $Id$
