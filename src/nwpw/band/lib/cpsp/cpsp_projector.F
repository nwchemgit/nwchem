
c $Id$

*     ************************************************
*     *                                              *
*     *             cpsp_projector_init              *
*     *                                              *
*     ************************************************
      subroutine cpsp_projector_init(npsp0)
      implicit none
      integer npsp0

#include "bafdecls.fh"
#include "cpsp_projector.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer ltag

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
      value = value.and.BA_alloc_get(mt_int,npsp,'list_prj',
     >                               list_prj(2),list_prj(1))
      if (.not.value) 
     > call errquit('cpsp_projector_init: out of heap memory',0,MA_ERR)

      do ltag=0,npsp-1
         int_mb(list_prj(1)+ltag) = ltag
      end do

      return
      end

*     ************************************************
*     *                                              *
*     *             cpsp_projector_alloc             *
*     *                                              *
*     ************************************************
*     
      integer function cpsp_projector_alloc(nbrill,nprj,nsize)
      implicit none
      integer nbrill,nprj,nsize

#include "bafdecls.fh"
#include "cpsp_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,hndl,indx

      tag    = int_mb(list_prj(1)+tagmax)
      int_mb(nbrill_prj(1)+tag) = nbrill
      int_mb(nprj_prj(1)+tag)   = nprj
      int_mb(nsize_prj(1)+tag)  = nsize

      value = BA_alloc_get(mt_dbl,nbrill*nprj*nsize,
     >                     'prj_alloc',hndl,indx)
      if (.not.value) 
     > call errquit('cpsp_projector_alloc: out of heap memory',0,MA_ERR)

      int_mb(indx_prj(1)+tag)  = indx
      int_mb(hndl_prj(1)+tag)  = hndl

      tagmax = tagmax + 1

      cpsp_projector_alloc = tag
      return
      end

*     ************************************************
*     *                                              *
*     *             cpsp_projector_dealloc           *
*     *                                              *
*     ************************************************

      subroutine cpsp_projector_dealloc(tagin)
      implicit none
      integer tagin

#include "bafdecls.fh"
#include "cpsp_projector.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ltagin,taglast,hndl
               
      hndl = int_mb(hndl_prj(1)+tagin)
      if (.not.BA_free_heap(hndl))
     >call errquit('cpsp_projector_dealloc:error freeing heap',0,MA_ERR)
               
      tagmax  = tagmax - 1
      taglast = int_mb(list_prj(1)+tagmax)
               
*     **** swap list_prj(tagin) with list_prj(tagmax) ***
      if (tagin.ne.taglast) then
         !*** find where tagin is in the list ****
         ltagin = 0 
         do while (int_mb(list_prj(1)+ltagin).ne.tagin)
             ltagin = ltagin + 1
         end do
         if (ltagin.gt.tagmax)
     > call errquit('cpsp_projector_dealloc:cannot find tagin',tagin,0)
         int_mb(list_prj(1) + ltagin) = int_mb(list_prj(1) + tagmax)
         int_mb(list_prj(1) + tagmax) = tagin
      end if
      return
      end

      
*     ************************************************
*     *                                              *
*     *             cpsp_projector_end               *
*     *                                              *
*     ************************************************
*     
      subroutine cpsp_projector_end()
      implicit none

#include "bafdecls.fh"
#include "cpsp_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,ltag,hndl,indx

      value = .true.
      do ltag=0,tagmax-1
         tag  = int_mb(list_prj(1)+ltag)
         hndl = int_mb(hndl_prj(1)+tag)
         value = value.and.BA_free_heap(hndl)
      end do
      value = value.and.BA_free_heap(nbrill_prj(2))
      value = value.and.BA_free_heap(nprj_prj(2))
      value = value.and.BA_free_heap(nsize_prj(2))
      value = value.and.BA_free_heap(indx_prj(2))
      value = value.and.BA_free_heap(hndl_prj(2))
      value = value.and.BA_free_heap(list_prj(2))
      if (.not.value) 
     > call errquit('cpsp_projector_end:error freeing heap',0,MA_ERR)

      return
      end 

*     ************************************************
*     *                                              *
*     *             cpsp_projector_add               *
*     *                                              *
*     ************************************************
*     
      subroutine cpsp_projector_add(tag,nb,n,proj)
      implicit none
      integer tag,nb,n
      real*8 proj(*)

#include "bafdecls.fh"
#include "cpsp_projector.fh"
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
*     *             cpsp_projector_get_ptr           *
*     *                                              *
*     ************************************************
*     
      integer function cpsp_projector_get_ptr(tag,nb,n)
      implicit none
      integer tag,nb,n

#include "bafdecls.fh"
#include "cpsp_projector.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer indx,nprj,nsize,shift

      indx  = int_mb(indx_prj(1)+tag)
      nprj  = int_mb(nprj_prj(1)+tag)
      nsize = int_mb(nsize_prj(1)+tag)
      cpsp_projector_get_ptr = indx + ((nb-1)*nprj + (n-1))*nsize
      return
      end
