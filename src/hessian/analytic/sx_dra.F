c
c $Id: sx_dra.F,v 1.1 2000-05-30 14:47:30 windus Exp $
c
c  This file holds the DRA write and read routines for the Sx matrices
c
      subroutine sx_write(g_sx, nmat, nsize, dra_handle)
c
      implicit none
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
        status = dra_write_section(.false.,g_sx(ind),1,nsize,1,nsize,
     &              dra_handle,nsize*(ind-1)+1,nsize*ind,1,nsize,
     &              request)
        if (status.ne.0) call errquit
     &     ('sx_write: problem with write for index ',ind)
c
c Make sure that the write is finished
c
        if (dra_wait(request).ne.0) call errquit
     &     ('sx_write: dra_wait failed for index',ind)
      enddo
c
      return
      end
c
c  The "reverse" of the previous routine
c
      subroutine sx_read(g_sx, nmat, nsize, dra_handle)
c
      implicit none
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
     &     ('sx_read: problem with read for index ',ind)
c
c Make sure that the write is finished
c
        if (dra_wait(request).ne.0) call errquit
     &     ('sx_read: dra_wait failed for index',ind)
      enddo
c
      return
      end
