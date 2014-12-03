*
* $Id$
*

*     ************************************************
*     *                                              *
*     *                cpsi_data_init                *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_init(npsi0)
      implicit none
      integer npsi0

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value
      integer ltag

*     **** external functions ****
      logical  control_brillioun_ondisk
      external control_brillioun_ondisk

      brill_ondisk = control_brillioun_ondisk()

      npsi   = 2*npsi0
      tagmax = 0
      tagmin = npsi-1

      value = BA_alloc_get(mt_int,npsi,'nsize_psi',
     >                     nsize_psi(2),nsize_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'neall_psi',
     >                               neall_psi(2),neall_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'nbrill_psi',
     >                               nbrill_psi(2),nbrill_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'indx_psi',
     >                               indx_psi(2),indx_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'hndl_psi',
     >                               hndl_psi(2),hndl_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'curnb_psi',
     >                               curnb_psi(2),curnb_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'update_psi',
     >                               update_psi(2),update_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'list_psi',
     >                               list_psi(2),list_psi(1))
      value = value.and.BA_alloc_get(mt_int,npsi,'next_psi',
     >                               next_psi(2),next_psi(1))

      call icopy(npsi,0,0,int_mb(update_psi(1)),1)
      do ltag=0,npsi-1
         int_mb(list_psi(1)+ltag) = ltag
      end do

      if (.not.value) 
     > call errquit('cpsi_data_init: out of heap memory',0,MA_ERR)

      return
      end
*     ************************************************
*     *                                              *
*     *                cpsi_data_end                 *
*     *                                              *
*     ************************************************
*     
      subroutine cpsi_data_end()
      implicit none

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,ltag,hndl

      if (tagmax.gt.0) write(*,*) "CPSI_DATA_END, tagmax=",tagmax
      value = .true.
      do ltag=0,tagmax-1
         tag  = int_mb(list_psi(1)+ltag)
         hndl = int_mb(hndl_psi(1)+tag)
         value = value.and.BA_free_heap(hndl)
      end do
      value = value.and.BA_free_heap(nbrill_psi(2))
      value = value.and.BA_free_heap(neall_psi(2))
      value = value.and.BA_free_heap(nsize_psi(2))
      value = value.and.BA_free_heap(indx_psi(2))
      value = value.and.BA_free_heap(hndl_psi(2))
      value = value.and.BA_free_heap(curnb_psi(2))
      value = value.and.BA_free_heap(update_psi(2))
      value = value.and.BA_free_heap(list_psi(2))
      value = value.and.BA_free_heap(next_psi(2))
      if (.not.value)
     > call errquit('cpsi_data_end:error freeing heap',0,MA_ERR)

      return
      end 




*     ************************************************
*     *                                              *
*     *                cpsi_data_alloc               *
*     *                                              *
*     ************************************************
*     
      integer function cpsi_data_alloc(nbrill,neall,nsize)
      implicit none
      integer nbrill,neall,nsize

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,hndl,indx


      tag = int_mb(list_psi(1)+tagmax)
      int_mb(nbrill_psi(1)+tag) = nbrill
      int_mb(neall_psi(1) +tag) = neall
      int_mb(nsize_psi(1) +tag) = nsize
      int_mb(curnb_psi(1) +tag) = 1
      int_mb(update_psi(1)+tag) = 0
      int_mb(next_psi(1)+tag)   = -1

      value = BA_alloc_get(mt_dbl,(nbrill*neall*nsize),
     >                     'cpsi_alloc',hndl,indx)
      if (.not.value) 
     > call errquit('cpsi_data_alloc: out of heap memory',0,MA_ERR)


      int_mb(indx_psi(1)+tag)  = indx
      int_mb(hndl_psi(1)+tag)  = hndl

      tagmax = tagmax + 1

      cpsi_data_alloc = tag
      return
      end

*     ************************************************
*     *                                              *
*     *                cpsi_data_dealloc             *
*     *                                              *
*     ************************************************
*     
      subroutine cpsi_data_dealloc(tagin)
      implicit none
      integer tagin

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer ltagin,taglast,hndl

      hndl = int_mb(hndl_psi(1)+tagin)
      if (.not. BA_free_heap(hndl))
     > call errquit('cpsi_data_dealloc: error freeing heap',0,MA_ERR)

      tagmax  = tagmax - 1
      taglast = int_mb(list_psi(1)+tagmax)

*     **** swap list_psi(tagin) with list_psi(tagmax) ***
      if (tagin.ne.taglast) then
         !*** find where tagin is in the list ****
         ltagin = 0
         do while (int_mb(list_psi(1)+ltagin).ne.tagin)
             ltagin = ltagin + 1
         end do
         if (ltagin.gt.tagmax)
     > call errquit('cpsi_data_dealloc:cannot find tagin',tagin,0)
         int_mb(list_psi(1) + ltagin) = int_mb(list_psi(1) + tagmax)
         int_mb(list_psi(1) + tagmax) = tagin
      end if
      return
      end



*     ************************************************
*     *                                              *
*     *                cpsi_data_push_stack          *
*     *                                              *
*     ************************************************
*     
      integer function cpsi_data_push_stack(nbrill,neall,nsize)
      implicit none
      integer nbrill,neall,nsize

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      logical value
      integer tag,hndl,indx

      tag    = tagmin
      int_mb(nbrill_psi(1)+tag) = nbrill
      int_mb(neall_psi(1) +tag) = neall
      int_mb(nsize_psi(1) +tag) = nsize
      int_mb(curnb_psi(1) +tag) = 1
      int_mb(update_psi(1)+tag) = 0
      int_mb(next_psi(1)+tag)   = -1

      value = BA_push_get(mt_dbl,nbrill*neall*nsize,
     >                     'cpsi_stack',hndl,indx)
      if (.not.value)
     > call errquit('cpsi_data_push_stack: out of stack',0,MA_ERR)

      int_mb(indx_psi(1)+tag)  = indx
      int_mb(hndl_psi(1)+tag)  = hndl

      tagmin = tagmin - 1

      cpsi_data_push_stack = tag
      return
      end

*     ************************************************
*     *                                              *
*     *                cpsi_data_pop_stack           *
*     *                                              *
*     ************************************************
*     
      subroutine cpsi_data_pop_stack(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer hndl

      if (tag.ne.(tagmin+1))
     > call errquit('cpsi_data_push_stack:tag!=tagmin+1',0,MA_ERR)

      tagmin = tagmin + 1
      hndl = int_mb(hndl_psi(1)+tag) 
      if (.not.BA_pop_stack(hndl))
     > call errquit('cpsi_data_pop_stack: error popping stack',1,MA_ERR)

      return
      end

      
*     ************************************************
*     *                                              *
*     *                 cpsi_data_add                *
*     *                                              *
*     ************************************************
*     
      subroutine cpsi_data_add(tag,nb,n,psi)
      implicit none
      integer tag,nb,n
      real*8 psi(*)

#include "bafdecls.fh"
#include "cpsi_data.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer indx,neall,nsize,shift,curr_nb

      indx  = int_mb(indx_psi(1) +tag)
      neall = int_mb(neall_psi(1)+tag)
      nsize = int_mb(nsize_psi(1)+tag)

      curr_nb = int_mb(curnb_psi(1)+tag)
c      if (brill_ondisk) then
c         shift = indx + ((nb-1)*nprj + (n-1))*nsize
c         .... add to local diskfile ***
c      else
c         shift = indx + ((n-1))*nsize
c      end if

      shift = indx + ((nb-1)*neall + (n-1))*nsize
      call dcopy(nsize,psi,1,dbl_mb(shift),1)

      return
      end 

*     ************************************************
*     *                                              *
*     *               cpsi_data_nsize                *
*     *                                              *
*     ************************************************
      integer function cpsi_data_nsize(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      cpsi_data_nsize = int_mb(nsize_psi(1)+tag)
      return
      end

*     ************************************************
*     *                                              *
*     *               cpsi_data_update               *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_update(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      int_mb(update_psi(1)+tag) = 1
      return
      end
*     ************************************************
*     *                                              *
*     *               cpsi_data_noupdate             *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_noupdate(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      int_mb(update_psi(1)+tag) = 0
      return
      end


*     ************************************************
*     *                                              *
*     *               cpsi_data_set_next             *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_set_next(tag,nexttag)
      implicit none
      integer tag,nexttag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      int_mb(next_psi(1)+tag) = nexttag
      return
      end



*     ************************************************
*     *                                              *
*     *               cpsi_data_get_next             *
*     *                                              *
*     ************************************************
      integer function cpsi_data_get_next(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      cpsi_data_get_next = int_mb(next_psi(1)+tag)
      return
      end



*     ************************************************
*     *                                              *
*     *               cpsi_data_get_ptr              *
*     *                                              *
*     ************************************************
*     
      integer function cpsi_data_get_ptr(tag,nb,n)
      implicit none
      integer tag,nb,n

#include "bafdecls.fh"
#include "cpsi_data.fh"

*     ***** local variables ****
      integer indx,neall,nsize

      indx  = int_mb(indx_psi(1) +tag)
      neall = int_mb(neall_psi(1)+tag)
      nsize = int_mb(nsize_psi(1)+tag)

      cpsi_data_get_ptr = indx +((nb-1)*neall+(n-1))*nsize
      return
      end


*     ************************************************
*     *                                              *
*     *               cpsi_data_get_chnk             *
*     *                                              *
*     ************************************************
*
      integer function cpsi_data_get_chnk(tag,nb)
      implicit none
      integer tag,nb

#include "bafdecls.fh"
#include "cpsi_data.fh"

*     ***** local variables ****
      integer indx,neall,nsize

      indx  = int_mb(indx_psi(1) +tag)
      neall = int_mb(neall_psi(1)+tag)
      nsize = int_mb(nsize_psi(1)+tag)

      cpsi_data_get_chnk = indx +(nb-1)*neall*nsize
      return
      end


*     ************************************************
*     *                                              *
*     *               cpsi_data_get_allptr           *
*     *                                              *
*     ************************************************
*
      integer function cpsi_data_get_allptr(tag)
      implicit none
      integer tag

#include "bafdecls.fh"
#include "cpsi_data.fh"

      cpsi_data_get_allptr = int_mb(indx_psi(1) +tag)
      return
      end


*     ************************************************
*     *                                              *
*     *               cpsi_data_copy_brill           *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_copy_brill(tag,nb1,nb2)
      implicit none
      integer tag,nb1,nb2

#include "bafdecls.fh"
#include "cpsi_data.fh"
      
*     ***** local variables ****
      integer indx,neall,nsize,shift1,shift2
      
      indx  = int_mb(indx_psi(1) +tag)
      neall = int_mb(neall_psi(1)+tag)
      nsize = int_mb(nsize_psi(1)+tag)

      shift1 = indx + (nb1-1)*neall*nsize
      shift2 = indx + (nb2-1)*neall*nsize
      call dcopy(neall*nsize,dbl_mb(shift1),1,dbl_mb(shift2),1)

      return
      end 

*     ************************************************
*     *                                              *
*     *               cpsi_data_copyall              *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_copyall(tag1,tag2)
      implicit none
      integer tag1,tag2

#include "bafdecls.fh"
#include "cpsi_data.fh"

*     ***** local variables ****
      integer nbrill,neall,nsize
      integer indx1,neall1,nsize1,nbrill1
      integer indx2,neall2,nsize2,nbrill2
      
      indx1   = int_mb(indx_psi(1)  +tag1)
      indx2   = int_mb(indx_psi(1)  +tag2)
      neall1  = int_mb(neall_psi(1) +tag1)
      neall2  = int_mb(neall_psi(1) +tag2)
      nsize1  = int_mb(nsize_psi(1) +tag1)
      nsize2  = int_mb(nsize_psi(1) +tag2)
      nbrill1 = int_mb(nbrill_psi(1)+tag1)
      nbrill2 = int_mb(nbrill_psi(1)+tag2)

      nbrill = nbrill1
      nsize  = nsize1
      neall  = neall1
      call dcopy(nbrill*neall*nsize,dbl_mb(indx1),1,dbl_mb(indx2),1)

      return
      end

*     ************************************************
*     *                                              *
*     *               cpsi_data_zeroall              *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_zeroall(tag1)
      implicit none
      integer tag1

#include "bafdecls.fh"
#include "cpsi_data.fh"

*     ***** local variables ****
      integer nbrill,neall,nsize
      integer indx1,neall1,nsize1,nbrill1

      indx1   = int_mb(indx_psi(1)  +tag1)
      neall1  = int_mb(neall_psi(1) +tag1)
      nsize1  = int_mb(nsize_psi(1) +tag1)
      nbrill1 = int_mb(nbrill_psi(1)+tag1)

      nbrill = nbrill1
      nsize  = nsize1
      neall  = neall1
      call dcopy(nbrill*neall*nsize,0.0d0,0,dbl_mb(indx1),1)
      return
      end


*     ************************************************
*     *                                              *
*     *               cpsi_data_scalall              *
*     *                                              *
*     ************************************************
      subroutine cpsi_data_scalall(alpha,tag1)
      implicit none
      real*8 alpha
      integer tag1

#include "bafdecls.fh"
#include "cpsi_data.fh"

*     ***** local variables ****
      integer nbrill,neall,nsize
      integer indx1,neall1,nsize1,nbrill1

      indx1   = int_mb(indx_psi(1)  +tag1)
      neall1  = int_mb(neall_psi(1) +tag1)
      nsize1  = int_mb(nsize_psi(1) +tag1)
      nbrill1 = int_mb(nbrill_psi(1)+tag1)

      nbrill = nbrill1
      nsize  = nsize1
      neall  = neall1
      call dscal(nbrill*neall*nsize,alpha,dbl_mb(indx1),1)
      return
      end

