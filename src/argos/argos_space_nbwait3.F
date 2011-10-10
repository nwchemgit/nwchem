      subroutine argos_space_nbwait3(lbbl)
c
      implicit none
c
#include "argos_space_common.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer ibbl
      integer lbbl(mbbl,mbb2)
      integer handle_iw,handle_fw,handle_is,handle_fs
      integer handle_fwl,handle_fsl
c
      do 1 ibbl=1,nbbl
c
      handle_iw=lbbl(ibbl,15)
      handle_fw=lbbl(ibbl,16)
      handle_is=lbbl(ibbl,17)
      handle_fs=lbbl(ibbl,18)
      handle_fwl=lbbl(ibbl,19)
      handle_fsl=lbbl(ibbl,20)
c
      if(handle_iw.ne.0) then
      call ga_nbwait(handle_iw)
      endif
      if(handle_fw.ne.0) then
      call ga_nbwait(handle_fw)
      endif
      if(handle_is.ne.0) then
      call ga_nbwait(handle_is)
      endif
      if(handle_fs.ne.0) then
      call ga_nbwait(handle_fs)
      endif
      if(handle_fwl.ne.0) then
      call ga_nbwait(handle_fwl)
      endif
      if(handle_fsl.ne.0) then
      call ga_nbwait(handle_fsl)
      endif
c
    1 continue
c
      return
      end
c $Id$
