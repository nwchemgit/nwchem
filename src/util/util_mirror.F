      subroutine util_mirrmat(nmat,g_a,g_amirr)
c $Id: util_mirror.F,v 1.1 2005-12-15 20:15:17 edo Exp $
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nmat,g_a(*),g_amirr(*)
c
      integer i
      do i=1,nmat
         call util_mirrclone(g_a(i),g_amirr(i))
         call ga_mask_sync(.true.,.false.)
         call ga_copy(g_a(i),g_amirr(i))

      enddo
      return
      end
      subroutine util_mirrclone(g_a,g_amirr)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nmat,g_a,g_amirr
c
      integer maxdim
      parameter(maxdim=10)
      integer gnblock(maxdim),plist
      character*32 array_name
      integer type, dims(maxdim),i,ndim
c
      do  i=1,maxdim
         gnblock(i)=0
      enddo
      plist=ga_pgroup_get_mirror()
      call nga_inquire(g_a, type, ndim, dims)
      write(0,123) ga_nodeid(), ' ndim ',ndim,' dims ',
     D     dims(1),dims(2),dims(3)
 123  format(I2,A,I2,A,3I4)
      if(ndim.gt.maxdim) call errquit(
     Q     ' mirrclone: increase maxdim  to ',ndim,0)
      if(.not.nga_create_config(MT_DBL, ndim, dims,
     .     'mirrored GA',
     .     gnblock,plist, g_amirr)) then
         call ga_inquire_name(g_a, array_name)
         write(0,*) ga_nodeid(), ' rep  failed on ',
     A        array_name,' ndim ',ndim,' dims ',
     D     dims(1),dims(2),dims(3)
         call errquit(' mirrclone failed ',0,0)
      endif

      return
      end
      subroutine util_mirrstop(ndens,g_densrep)
      implicit none
#include "global.fh"
      integer ndens,g_densrep(*)
c
      integer i
      do i=1,ndens
         if(.not.ga_destroy(g_densrep(i))) call errquit
     Q        (' giarepdm destroy ',0,0)
      enddo
      return
      end
      subroutine util_mirrmerge(g_densrep,g_dens)
      implicit none
      integer g_densrep,g_dens
c
      call ga_mask_sync(.true.,.false.)
      call ga_merge_mirrored(g_densrep)
      call ga_mask_sync(.true.,.false.)
      call ga_copy(g_densrep,g_dens)
      return
      end
