c $Id: bq_input.F,v 1.2 2005-03-30 23:34:28 marat Exp $
      subroutine bq_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "bq_params.fh"
#include "bq.fh"
c     
      integer rtdb
c
      character*32 pname
      character*255 field
      character*255 namespace
      integer max_nbq
      logical status
      integer h_cbq,i_cbq
      integer h_qbq,i_qbq
      integer ncenter
      double precision x,y,z,q
      integer i
      integer handle
c
      pname = "bq_input"
c
      write(*,*) "in "//pname
c
      max_nbq = min( 4*1000, ma_inquire_avail(MT_INT))
      max_nbq = max_nbq/4
c
      if ( .not. ma_push_get( MT_DBL, 3*max_nbq, 
     &      'cbq scratch',
     &      h_cbq, i_cbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      max_nbq, MA_ERR)
c
      if ( .not. ma_push_get( MT_DBL, max_nbq, 
     &      'qbq scratch',
     &      h_qbq, i_qbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      max_nbq, MA_ERR)

c
c     set namespace if any
c     --------------------
      if(inp_n_field().eq.1) then
        namespace=bq_default
      else if(inp_n_field().eq.2) then
        call inp_set_field(1)
        if (.not. inp_a(field)) then
          call errquit('cons_input: inp_read failed',0,0)
       else
          namespace=field
       end if
      end if 
      
      ncenter = 0
 20   if (inp_read()) then
       if(inp_a(field)) then
         if (inp_compare(.false., 'end', field)) then
           goto 30
         end if
         call inp_set_field(0)
         status = .true.
         status = status .and. inp_f(q)
         status = status .and. inp_f(x)
         status = status .and. inp_f(y)         
         status = status .and. inp_f(z)         
         write(*,*) q,x,y,z
          if (.not. status) 
     &       call errquit('bq_input: inp_read failed',0,0)
         dbl_mb(i_qbq+ncenter) = q
         dbl_mb(i_cbq+3*ncenter)   = x
         dbl_mb(i_cbq+3*ncenter+1) = y
         dbl_mb(i_cbq+3*ncenter+2) = z
         ncenter = ncenter + 1
         status = ncenter .le. max_nbq 
         if (.not. status) 
     $       call errquit('bq_input: too many centers',0,0)
         go to 20
       end if
      end if
 30   continue

      if(.not.bq_create(namespace,handle))   
     &       call errquit('bq_create failed',0,0)

      if(.not.bq_set(handle,ncenter,dbl_mb(i_qbq),dbl_mb(i_cbq)))
     &       call errquit('bq_set failed',0,0)

      call bq_print_info(handle)
      
      do i=1,ncenter
         dbl_mb(i_qbq+i) = 0
         dbl_mb(i_cbq+3*i)   = 0
         dbl_mb(i_cbq+3*i+1) = 0
         dbl_mb(i_cbq+3*i+2) = 0
      end do

      if(.not.bq_get(handle,ncenter,dbl_mb(i_qbq),dbl_mb(i_cbq)))
     &       call errquit('bq_set failed',0,0)

      write(*,*) "testing"
      do i=0,ncenter-1
         write(*,*) dbl_mb(i_qbq+i) 
         write(*,*) dbl_mb(i_cbq+3*i)   
         write(*,*) dbl_mb(i_cbq+3*i+1) 
         write(*,*) dbl_mb(i_cbq+3*i+2) 
      end do

      ncenter = 0
      if(.not.bq_ncenter(handle,ncenter))
     &       call errquit('bq_ncenter failed',0,0)
      
      write(*,*) "bq_ncenter",ncenter
      write(*,*) "testing"

      if(.not.bq_rtdb_store(rtdb,handle))
     &       call errquit('bq_store failed',0,0)

      if(.not.rtdb_print(rtdb))
     &       call errquit('rtdb print failed',0,0)

      if(.not.bq_destroy(handle))
     &       call errquit('bq_destroy failed',0,0)

      write(*,*) "before load"
      call bq_print_info(handle)

      if(.not.bq_create(namespace,handle))   
     &       call errquit('bq_create failed',0,0)

      if(.not.bq_rtdb_load(rtdb,handle))   
     &       call errquit('bq_create failed',0,0)

      write(*,*) "after load"
      call bq_print_info(handle)
  
      return
      end
c

