c
c $Id$
c
c  This file holds the DRA write and read routines for the Sx matrices
c
      subroutine sx_write(g_sx, nmat, nsize, dra_handle)
c
      implicit none
#include "errquit.fh"
c
#include "dra.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer nmat              ! [input] number of Sx matrices
      integer g_sx              ! [input] GA handles for Sx matrices
      integer nsize             ! [input] dimension of Sx (nbf)
      integer dra_handle        ! [input] handle for the DRA
c
      integer loa(3), hia(3), lob(2), hib(2)
      integer ind, status, request, g_temp
c
c Loop over each ga in g_sx and put out to disk, the dra has dimensions
c i = 1, nsize*nmat and j = 1, nsize
c
      loa(2) = 1
      hia(2) = nsize
      loa(3) = 1
      hia(3) = nsize
      lob(1) = 1
      hib(1) = nsize
      lob(2) = 1
      hib(2) = nsize
      if (.not.ga_create(MT_DBL, nsize, nsize, 'g_temp', -1, -1, 
     &     g_temp)) call errquit('sx_dra: ga_create failed',555,
     &       GA_ERR)
c
      do ind = 1, nmat
        loa(1) = ind
        hia(1) = ind
        call nga_copy_patch('N', g_sx, loa, hia, g_temp, lob, hib)
        status = dra_write_section(.false.,g_temp,1,nsize,1,nsize,
     &              dra_handle,nsize*(ind-1)+1,nsize*ind,1,nsize,
     &              request)
        if (status.ne.0) call errquit
     &     ('sx_write: problem with write for index ',ind, DISK_ERR)
c
c Make sure that the write is finished
c
        if (dra_wait(request).ne.0) call errquit
     &     ('sx_write: dra_wait failed for index',ind, DISK_ERR)
      enddo
c
      if (.not. ga_destroy(g_temp))
     &  call errquit('sx_write: could not destroy g_temp',555, GA_ERR)
c
      return
      end
c
c  The "reverse" of the previous routine
c
      subroutine sx_read(g_sx, nmat, nsize, dra_handle)
c
      implicit none
#include "errquit.fh"
c
#include "dra.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer nmat              ! [input] number of Sx matrices
      integer g_sx(nmat)        ! [input] GA handles for Sx matrices
      integer nsize             ! [input] dimension of Sx (nbf)
      integer dra_handle        ! [input] handle for the DRA
c
      integer loa(3), hia(3), lob(2), hib(2)
      integer ind, status, request, g_temp
c
c Loop over each ga in g_sx and put out to disk, the dra has dimensions
c i = 1, nsize*nmat and j = 1, nsize
c
      loa(2) = 1
      hia(2) = nsize
      loa(3) = 1
      hia(3) = nsize
      lob(1) = 1
      hib(1) = nsize
      lob(2) = 1
      hib(2) = nsize
      if (.not.ga_create(MT_DBL, nsize, nsize, 'g_temp', -1, -1,
     &     g_temp)) call errquit('sx_dra: ga_create failed',555, GA_ERR)
c
      do ind = 1, nmat
        loa(1) = ind
        hia(1) = ind
        status = dra_read_section(.false.,g_temp,1,nsize,1,nsize,
     &              dra_handle,nsize*(ind-1)+1,nsize*ind,1,nsize,
     &              request)
c
c Make sure that the write is finished
c
        if (dra_wait(request).ne.0) call errquit
     &     ('sx_read: dra_wait failed for index',ind, CALC_ERR)
c
        call nga_copy_patch('N', g_temp, lob, hib, g_sx, loa, hia)
        if (status.ne.0) call errquit
     &     ('sx_read: problem with read for index ',ind, GA_ERR)
      enddo
c
      if (.not. ga_destroy(g_temp))
     &  call errquit('sx_read: could not destroy g_temp',555, GA_ERR)
c
      return
      end
c
c Original read without 3D GA
c
      subroutine sx_read_orig(g_sx, nmat, nsize, dra_handle)
c
      implicit none
#include "errquit.fh"
c
#include "dra.fh"
c
      integer nmat              ! [input] number of Sx matrices
      integer g_sx(nmat)        ! [input] GA handles for Sx matrices
      integer nsize             ! [input] dimension of Sx (nbf)
      integer dra_handle        ! [input] handle for the DRA
c
      integer ind, status, request
c
c Loop over each ga in g_sx and put out to disk, the dra has dimensions
c i = 1, nsize*nmat and j = 1, nsize
c
      do ind = 1, nmat
        status = dra_read_section(.false.,g_sx(ind),1,nsize,1,nsize,
     &              dra_handle,nsize*(ind-1)+1,nsize*ind,1,nsize,
     &              request)
        if (status.ne.0) call errquit
     &     ('sx_read: problem with read for index ',ind, UNKNOWN_ERR)
c
c Make sure that the write is finished
c
        if (dra_wait(request).ne.0) call errquit
     &     ('sx_read: dra_wait failed for index',ind, CALC_ERR)
      enddo
c
      return
      end
