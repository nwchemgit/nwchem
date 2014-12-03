
c $Id$

*     ************************************************
*     *                                              *
*     *             borb_projector_init              *
*     *                                              *
*     ************************************************
      subroutine borb_projector_init(npsp0)
      implicit none
      integer npsp0

#include "bafdecls.fh"
#include "borb_projector.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value

*     **** external functions ****
      logical  control_brillioun_ondisk
      external control_brillioun_ondisk

      brill_ondisk = control_brillioun_ondisk()
      curr_nb       = 1

      npsp   = npsp0
      tagmax = 0

      value = BA_alloc_get(mt_int,npsp,'nsize_prj',
     >                     nsize_prj(2),nsize_prj(1))
      value = value.and.BA_alloc_get(mt_int,npsp,'nprj_prj',
     >                               nprj_prj(2),nprj_prj(1))
      value = value.and.BA_alloc_get(mt_int,npsp,'nbrill_prj',
     >                               nbrill_prj(2),nbrill_prj(1))
      value = value.and.BA_alloc_get(mt_int,npsp,'indx_prj',
     >                               indx_prj(2),indx_prj(1))
      value = value.and.BA_alloc_get(mt_int,npsp,'hndl_prj',
     >                               hndl_prj(2),hndl_prj(1))
      if (.not.value) 
     > call errquit('borb_projector_init: out of heap memory',0,MA_ERR)

      return
      end

*     ************************************************
*     *                                              *
*     *             borb_projector_alloc             *
*     *                                              *
*     ************************************************
*     
      integer function borb_projector_alloc(nbrill,nprj,nsize)
      implicit none
      integer nbrill,nprj,nsize

#include "bafdecls.fh"
#include "borb_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,hndl,indx

      tag    = tagmax 
      int_mb(nbrill_prj(1)+tag) = nbrill
      int_mb(nprj_prj(1)+tag)   = nprj
      int_mb(nsize_prj(1)+tag)  = nsize

      value = BA_alloc_get(mt_dbl,nbrill*nprj*nsize,
     >                     'prj_alloc',hndl,indx)
      if (.not.value) 
     > call errquit('borb_projector_alloc: out of heap memory',0,MA_ERR)

      int_mb(indx_prj(1)+tag)  = indx
      int_mb(hndl_prj(1)+tag)  = hndl

      tagmax = tagmax + 1

      borb_projector_alloc = tag
      return
      end

      
*     ************************************************
*     *                                              *
*     *             borb_projector_end               *
*     *                                              *
*     ************************************************
*     
      subroutine borb_projector_end()
      implicit none

#include "bafdecls.fh"
#include "borb_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,hndl,indx

      value = .true.
      do tag=0,tagmax-1
         hndl = int_mb(hndl_prj(1)+tag)
         value = value.and.BA_free_heap(hndl)
      end do
      value = value.and.BA_free_heap(nbrill_prj(2))
      value = value.and.BA_free_heap(nprj_prj(2))
      value = value.and.BA_free_heap(nsize_prj(2))
      value = value.and.BA_free_heap(indx_prj(2))
      value = value.and.BA_free_heap(hndl_prj(2))
      if (.not.value) 
     > call errquit('borb_projector_end:error freeing heap',0,MA_ERR)

      return
      end 

*     ************************************************
*     *                                              *
*     *             borb_projector_add               *
*     *                                              *
*     ************************************************
*     
      subroutine borb_projector_add(tag,nb,n,proj)
      implicit none
      integer tag,nb,n
      real*8 proj(*)

#include "bafdecls.fh"
#include "borb_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer indx,nprj,nsize,shift

      indx  = int_mb(indx_prj(1)+tag)
      nprj  = int_mb(nprj_prj(1)+tag)
      nsize = int_mb(nsize_prj(1)+tag)

      curr_nb = nb
c      if (brill_ondisk) then
c         shift = indx + ((nb-1)*nprj + (n-1))*nsize
c         .... add to local diskfile ***
c      else
c         shift = indx + ((n-1))*nsize
c      end if

      shift = indx + ((nb-1)*nprj + (n-1))*nsize
      call dcopy(nsize,proj,1,dbl_mb(shift),1)

      return
      end 

*     ************************************************
*     *                                              *
*     *             borb_projector_get_ptr           *
*     *                                              *
*     ************************************************
*     
      integer function borb_projector_get_ptr(tag,nb,n)
      implicit none
      integer tag,nb,n

#include "bafdecls.fh"
#include "borb_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer indx,nprj,nsize,shift

      indx  = int_mb(indx_prj(1)+tag)
      nprj  = int_mb(nprj_prj(1)+tag)
      nsize = int_mb(nsize_prj(1)+tag)
      borb_projector_get_ptr = indx + ((nb-1)*nprj + (n-1))*nsize
      return
      end
