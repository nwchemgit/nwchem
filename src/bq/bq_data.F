c $Id: bq_data.F,v 1.1 2005-03-30 02:21:41 marat Exp $
      block data initial_bq_data
C$Id: bq_data.F,v 1.1 2005-03-30 02:21:41 marat Exp $
      implicit none
#include "bq_data.fh"
c     
      integer i             
      data bq_active /max_bq*.false./
      data bq_ncent /max_bq*0/

      end block data initial_bq_data

      function bq_create(namespace,handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       character*(*) namespace
       logical bq_create
c      local variables
       integer i
       integer handle

       bq_create = .false.
       do i=1,max_bq
        if(bq_ncent(i).eq.0) then
          bq_create = .true.
          bq_name(i) = namespace
          handle = i
          return
        end if
       end do
       
       return
       end

      function bq_set(handle,n,q,c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       integer n
       double precision q(n)
       double precision c(3*n)
       logical bq_set
c      local variables
       integer i
       integer h_c,i_c
       integer h_q,i_q
       character(32) pname

       pname = "bq_set"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_set = .false.
         return
       else
         bq_set = .true.
       end if

       if(.not.ma_alloc_get(MT_DBL, 3*n, 'bqdata c',
     &      h_c, i_c) ) call errquit(
     &      pname//' unable to allocate heap space',
     &      3*n, MA_ERR)

       if(.not.ma_alloc_get(MT_DBL, n, 'bqdata q',
     &      h_q, i_q) ) call errquit(
     &      pname//' unable to allocate heap space',
     &      n, MA_ERR)

      
       do i=1,n
         dbl_mb(i_q+i-1) = q(i)
       end do 

       do i=1,3*n
         dbl_mb(i_c+i-1) = c(i)
       end do 

       bq_ncent(handle)  = n
       bq_charge(handle) = h_q
       bq_coord(handle)  = h_c

       return
       end

      function bq_get(handle,n,q,c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       integer n
       double precision q(n)
       double precision c(3*n)
       logical bq_get
c      local variables
       integer i
       integer h_c,i_c
       integer h_q,i_q
       character(32) pname

       pname = "bq_get"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_get = .false.
         return
       else
         bq_get = .true.
       end if

       if(n.ne.bq_ncent(handle)) then
         bq_get = .false.
         return
       end if


       h_q = bq_charge(handle)
       h_c = bq_coord(handle)

       if(.not.ma_get_index( h_c, i_c) ) call errquit(
     &      pname//' unable to locate coord handle',
     &      0, MA_ERR)


       if(.not.ma_get_index( h_q, i_q) ) call errquit(
     &      pname//' unable to locate charge handle',
     &      0, MA_ERR)

       do i=1,n
         q(i) = dbl_mb(i_q+i-1) 
       end do 

       do i=1,3*n
         c(i) = dbl_mb(i_c+i-1) 
       end do 

       return
       end

      subroutine bq_print_info(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"

       integer handle
c      local variables
       integer i
       integer k
       integer h_c,i_c
       integer h_q,i_q
       character(32) pname

       pname = "bq_print_info"

       if(handle .lt.0 .or. handle .gt. max_bq) then
           call errquit(
     &      pname//' unable to locate handle ',
     &      0,0)
       end if

       h_q = bq_charge(handle)
       h_c = bq_coord(handle)

       if(.not.ma_get_index( h_c, i_c) ) call errquit(
     &      pname//' unable to locate coord handle',
     &      0, MA_ERR)


       if(.not.ma_get_index( h_q, i_q) ) call errquit(
     &      pname//' unable to locate charge handle',
     &      0, MA_ERR)

       if(ga_nodeid().eq.0) then
      
         write(*,*) "bq information"
         write(*,*) "  name", bq_name(handle)
         write(*,*) "  number of charges",bq_ncent(handle)
         do i=1,bq_ncent(handle)
           write(*,'(4F12.6)') dbl_mb(i_q+i-1),
     >          (dbl_mb(i_c+3*(i-1)+k-1),k=1,3) 
         end do

       end if

       return
       end

