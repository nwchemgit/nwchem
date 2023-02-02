      logical function util_mirrmat(nmat,g_a,g_amirr,
     D     docopy, dozero)
c $Id$
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nmat,g_a(*),g_amirr(*)
      logical docopy, dozero
      logical util_mirrclone,util_mirrmat_chunk
      external util_mirrclone,util_mirrmat_chunk
c
      integer maxdim
      parameter(maxdim=10)
      integer gnblock(maxdim)
      integer type, dims(maxdim)
c
      integer i,ii
      integer ndim
      util_mirrmat=.false.
#ifdef DISABLE_GAMIRROR
#warning disabling gamirror
      return
#endif
      call nga_inquire(g_a, type, ndim, dims)
      if(ndim.gt.maxdim) call errquit(
     Q     ' mirrmat: increase maxdim  to ',ndim,0)
c
      do  i=1,maxdim
         gnblock(i)=0
      enddo
      util_mirrmat=util_mirrmat_chunk(nmat,g_a,g_amirr,
     D     docopy, dozero, gnblock)
      return
      end
      logical function util_mirrmat_chunk(nmat,g_a,g_amirr,
     D     docopy, dozero, gnblock)
c $Id$
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer nmat,g_a(*),g_amirr(*)
      logical docopy, dozero
      integer maxdim
      integer gnblock(*)
      logical util_mirrclone
      external util_mirrclone
c
      integer i,ii
      util_mirrmat_chunk=.false.
#ifdef DISABLE_GAMIRROR
#warning disabling gamirror
      return
#endif
c
c     need to add sanity check on docopy vs dozero
      do i=1,nmat
         if(util_mirrclone(g_a(i),g_amirr(i),gnblock)) then
            call ga_mask_sync(.true.,.false.)
            if(docopy) call ga_copy(g_a(i),g_amirr(i))
            if(dozero) call ga_zero(g_amirr(i))
         else
            if(i.gt.1) then
               do ii=1,i-1
                  call util_mirrstop(g_amirr(ii))
               enddo
            endif
            return
         endif
      enddo
      util_mirrmat_chunk=.true.
      return
      end
      logical function util_mirrclone(g_a,g_amirr,gnblock)
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
      plist=ga_pgroup_get_mirror()
      call nga_inquire(g_a, type, ndim, dims)
#ifdef DEBUG
      write(0,123) ga_nodeid(), ' ndim ',ndim,' dims ',
     D     dims(1),dims(2),dims(3)
 123  format(I2,A,I2,A,3I4)
#endif
      if(ndim.gt.maxdim) call errquit(
     Q     ' mirrclone: increase maxdim  to ',ndim,0)
      if(.not.nga_create_config(MT_DBL, ndim, dims,
     .     'mirrored GA',
     .     gnblock,plist, g_amirr)) then
         call ga_inquire_name(g_a, array_name)
         write(0,*) ga_nodeid(), ' rep  failed on ',
     A        array_name,' ndim ',ndim,' dims ',
     D     dims(1),dims(2),dims(3)
         util_mirrclone=.false.
      else
         util_mirrclone=.true.
      endif

      return
      end
      subroutine util_mirrstop(g_densrep)
      implicit none
#include "global.fh"
      integer ndens,g_densrep
c
         if(.not.ga_destroy(g_densrep)) call errquit
     Q        (' giarepdm destroy ',0,0)
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
      subroutine util_mirrmrgadd(g_densrep,g_dens)
      implicit none
      integer g_densrep,g_dens,g_densdis
#include "global.fh"
c
      if(.not.ga_duplicate(g_dens, g_densdis, 'GA distr'))
     C     call errquit(' utilmirr: gaduplicate failed ',0,0)
      call util_mirrmerge(g_densrep,g_densdis)
      call util_mirrstop(g_densrep)
      call ga_add(1d0,g_densdis,1d0,g_dens,g_dens)
      call util_mirrstop(g_densdis)
      return
      end
      integer function util_mirrmem(nmat,g_a)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
      integer g_a(*)
      integer nmat
c
      integer ppn
      integer nbytes,i,imat
      integer maxdim
      parameter(maxdim=10)
      integer type, dims(maxdim),ndim
      util_mirrmem=0
      call util_getppn(ppn)
      do imat=1,nmat
         call nga_inquire(g_a(imat), type, ndim, dims)
      if(ndim.gt.maxdim) call errquit(
     Q     ' mirrmem: increase maxdim  to ',ndim,0)
#if 1
      if(ga_nodeid().eq.0) then
         write(0,123)  ' ndim ',ndim,' dims ',
     D        (dims(i),i=1,ndim)
      endif
 123     format(A,I2,A,3I10)
#endif
         nbytes=1
         do i=1,ndim
            nbytes=nbytes*dims(i)
         enddo
         nbytes=ma_sizeof(mt_dbl,1,mt_byte)*nbytes
         util_mirrmem=util_mirrmem+nbytes
      enddo
      if(ga_nnodes().ne.1) then
c     divide by ppn & round it by 30%
         util_mirrmem=(util_mirrmem*130)/(ppn*100)
      endif
      return
      end
