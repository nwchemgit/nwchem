      subroutine argos_pme_flag(iflag,icolor,inode)
c
      implicit none
c
#include "argos_pme_common.fh"
c
      integer iflag,icolor,inode
      integer il,ih,jl,jh,node,icolr
c
c     set flag iflag on node inode to color icolor:
c
c     icolor>=0 set flag to icolor
c     icolor<0  accumulate iabs(icolor) to flag
c
c     inode>=0  set/accumulate flag on node inode
c     inode<0   set/accumulate flag on all nodes
c
      icolr=iabs(icolor)
c
      if(inode.lt.0) then
      do 1 node=0,np-1
      call ga_distribution(ga_flg,node,il,ih,jl,jh)
      if(icolor.ge.0) then
      call ga_put(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1,1)
      endif
    1 continue
      else
c
c     set flag on single node
c
      call ga_distribution(ga_flg,inode,il,ih,jl,jh)
      if(icolor.ge.0) then
      call ga_put(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1)
      else
      call ga_acc(ga_flg,il,ih,jl+iflag,jl+iflag,icolr,1,1)
      endif
      endif
c
      return
      end
c $Id$
