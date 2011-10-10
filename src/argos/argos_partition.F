      subroutine argos_partition
c
      implicit none
c
#include "argos.fh"
#include "global.fh"
#include "util.fh"
c
c  size should equal maximum number of processors that could be in
c  a group
c
      integer MAX_GRP_SIZE
      parameter (MAX_GRP_SIZE = 1000)
      integer grp_size
      integer list(MAX_GRP_SIZE)
      integer i, grp_handle
      real*8 rantmp
c
      if(npg.eq.0) npg=1
c
      npw=ga_nnodes()
      mew=ga_nodeid()
c
      if(npg.ne.1) then
      rantmp=util_random(iseed*(me+1))
      if (mod(npw,npg).ne.0) then
      call md_abort('Number of processes is not a multiple',
     + ' of number of process groups',mew)
      endif
      grp_size=npw/npg
      meg =(mew-mod(mew,grp_size))/grp_size
      do 1 i=1,grp_size
      list(i)=meg*grp_size+i-1
    1 continue
      grp_handle=ga_pgroup_create(list,grp_size)
      call ga_pgroup_set_default(grp_handle)
      endif
c
      np=ga_nnodes()
      me=ga_nodeid()
c
      if(me.eq.0) then
      if(util_print('distribution', print_high)) then
      write(*,1000) npg
 1000 format(' The number of processor groups requested is ',i5)
      write(*,1001) npw
 1001 format(' The number of processors is ',i5)
c
      if(npg.ne.1) then
      write(*,1002) ga_pgroup_get_world()
 1002 format(' The world group handle is ',i5)
      write(*,1003) ga_pgroup_get_default()
 1003 format(' The default group handle is ',i5)
      endif
c
      endif
      endif
c
      return
      end
c
c $Id$
