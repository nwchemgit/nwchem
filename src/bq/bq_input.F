c $Id: bq_input.F,v 1.5 2005-04-04 22:58:58 marat Exp $
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
      double precision  scale 
c
      pname = "bq_input"
c
      write(*,*) "in "//pname
c
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

      if(inp_n_field().gt.2) 
     >  call errquit(pname//'too many fields',0,0)

      namespace = bq_default
      scale = 1.88972598858d0

      do i=1,inp_n_field()-1
        call inp_set_field(i)
        if (.not. inp_a(field))   
     &    call errquit(pname//': inp_read failed',0,0)
        if (inp_compare(.false.,field,'au')) then
          scale = 1.0d0       
        else if (inp_compare(.false.,field,'angstroms')) then
          scale = 1.88972598858d0
        else
         namespace=field 
        end if 
      end do

     
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

      call dscal(3*ncenter, scale, dbl_mb(i_cbq), 1)

      if(.not.bq_create(namespace,handle))   
     &       call errquit('bq_create failed',0,0)

      if(.not.bq_set(handle,ncenter,dbl_mb(i_qbq),dbl_mb(i_cbq)))
     &       call errquit('bq_set failed',0,0)

      call bq_print_info(handle)

      if(.not.bq_rtdb_store(rtdb,handle))
     &       call errquit('bq_store failed',0,0)

      if(.not.bq_destroy(handle))
     &       call errquit('bq_destroy failed',0,0)

c
c     activate default namespace if found
c     -----------------------------------
      if(ncenter.ne.0 .and. namespace.eq.bq_default) then
      if(.not. rtdb_cput(rtdb,"bq" , 1,bq_default)) then
         call errquit(pname//'failed to set default namespace',
     >                0,RTDB_ERR)
      end if
      end if

      if (.not. ma_pop_stack(h_qbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      max_nbq, MA_ERR)

      if (.not. ma_pop_stack(h_cbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      max_nbq, MA_ERR)
 
      return
      end
c

