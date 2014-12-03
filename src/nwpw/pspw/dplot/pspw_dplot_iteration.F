*
* $Id$
*
      subroutine dplot_iteration_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "btdb.fh"

*     **** dplot_iteration_block ****
      logical exist
      integer iteration(2),iteration_count
      common / dplot_iteration_Block / iteration,iteration_count,exist

*     **** local variables ****
      integer    rtdb,ma_type

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb
      
      rtdb = control_rtdb()

*     **** get dplot_iteration_lists have been imposed
      if (rtdb_ma_get(rtdb, 'pspw_dplot:iteration_list', ma_type,
     >        iteration_count, iteration(2))) then
        exist = .true.

         if (.not.BA_get_index(iteration(2),iteration(1)))
     >      call errquit(
     >       'pspw_dplot_iteration_init: ma_get_index failed',0, MA_ERR)
      else
        exist = .false.
      end if
      return
      end


      subroutine dplot_iteration_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** dplot_iteration_block ****
      logical exist
      integer iteration(2),iteration_count
      common / dplot_iteration_Block / iteration,iteration_count,exist

      if (exist) then

*        **** free heap ****
         if (.not. BA_free_heap(iteration(2))) then
            call errquit('dplot_iteration_end:freeing heap',0, MA_ERR)
         end if

      end if

      return
      end



      logical function dplot_iteration_check(it)
      implicit none
      integer it

#include "bafdecls.fh"

*     **** dplot_iteration_block ****
      logical exist
      integer iteration(2),iteration_count
      common / dplot_iteration_Block / iteration,iteration_count,exist

*     **** local variables ****
      logical value
      integer ii

      
      value = .false.

      if (exist) then

         ii=1
         do while ((ii.le.iteration_count).and.(.not.value))
           if (it.eq.int_mb(iteration(1)+ii-1)) value = .true.
           ii = ii + 1
         end do

      end if

      dplot_iteration_check = value
      return
      end

      subroutine dplot_iteration(it,ispin,ne,psi,dn,psi_r)
      implicit none
      integer    it
      integer    ispin,ne(2)
      complex*16 psi(*)
      real*8     dn(*)
      real*8     psi_r(*)

#include "bafdecls.fh"
#include "btdb.fh"

*     **** local variables ****
      integer nfft3d,npack1,nemax,rtdb
      character*8 tag

*     **** external functions ****
      integer  control_rtdb
      external control_rtdb

c
      if      (it .lt. 10) then
         write(tag,'(I1)') it
      else if (it .lt. 100) then
         write(tag,'(I2)') it
      else if (it .lt. 1000) then
         write(tag,'(I3)') it
      else if (it .lt. 10000) then
         write(tag,'(I4)') it
      else if (it .lt. 100000) then
         write(tag,'(I5)') it
      else
         write(tag,'(I6)') it
      end if


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(1,npack1)
      nemax = ne(1)+ne(2)
      rtdb  = control_rtdb()

      call dplot_loop(rtdb,
     >                ispin,ne,
     >                npack1,nfft3d,nemax,
     >                psi,dn,psi_r,
     >                .true.,tag)

      return
      end
