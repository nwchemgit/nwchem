c $Id: bq_data.F,v 1.12 2008-03-18 18:10:11 marat Exp $
      block data initial_bq_data
      implicit none
#include "bq_data.fh"
c     
      data bq_active /max_bq*.false./
      data bq_mem /max_bq*.false./
      data bq_ncent /max_bq*0/

      end

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

      function bq_get_handle(namespace,handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       character*(*) namespace
       logical bq_get_handle
c      local variables
       integer i
       integer handle

       bq_get_handle = .false.
       do i=1,max_bq
        if(bq_name(i).eq.namespace) then
          bq_get_handle = .true.
          handle = i
          return
        end if
       end do
       
       return
       end

      function bq_activate(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_activate
c      local variables
       integer i

       bq_activate = .true.
       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_activate = .false.
         write(*,*) "bq handle is out of bounds"
         return
       end if
       if(bq_ncent(handle).eq.0) then
         bq_activate = .false.
         write(*,*) "bq ncent is zero"
         return
       end if

       do i=1,max_bq
        bq_active(i)=.false.
       end do
        
       bq_active(handle) = .true.
       
       return
       end

      function bq_deactivate(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_deactivate
c      local variables

       bq_deactivate = .true.
       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_deactivate = .false.
         return
       end if
       if(bq_ncent(handle).eq.0) then
         bq_deactivate = .false.
         return
       end if

       bq_active(handle) = .false.
       
       return
       end

      function bq_get_active(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_get_active
c      local variables
       integer i

       do i=1,max_bq
        if(bq_active(i)) then
          bq_get_active = .true.
          handle = i
          return
        end if
       end do
        
       bq_get_active = .false.
       
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
       character*(32) pname

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
       bq_mem(handle)    = .true.

       return
       end

      function bq_alloc(handle,n)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       integer n
       logical bq_alloc
c      local variables
       integer i
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname

       pname = "bq_alloc"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_alloc = .false.
         return
       else
         bq_alloc = .true.
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
         dbl_mb(i_q+i-1) = 0.0d0
       end do 

       do i=1,3*n
         dbl_mb(i_c+i-1) = 0.0d0
       end do 

       bq_ncent(handle)  = n
       bq_charge(handle) = h_q
       bq_coord(handle)  = h_c
       bq_mem(handle)    = .true.

       return
       end

      function bq_pset(handle,n,h_q,h_c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       integer n
       integer h_q,h_c
       logical bq_pset
c      local variables
       character*(32) pname

       pname = "bq_pset"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_pset = .false.
         return
       else
         bq_pset = .true.
       end if

       bq_ncent(handle)  = n
       bq_charge(handle) = h_q
       bq_coord(handle)  = h_c
       bq_mem(handle)    = .false.

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
       character*(32) pname

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

      function bq_ncenter(handle,n)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       integer n
       logical bq_ncenter
c      local variables
       character*(32) pname

       pname = "bq_ncenter"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_ncenter = .false.
         return
       else
         bq_ncenter = .true.
       end if

       n = bq_ncent(handle)

       return
       end

      function bq_namespace(handle,namespace)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       character*(*) namespace
       logical bq_namespace
c      local variables
       character*(32) pname

       pname = "bq_namespace"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_namespace = .false.
         return
       else
         bq_namespace = .true.
       end if

       namespace = bq_name(handle)

       return
       end

      function bq_check_handle(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_check_handle
c      local variables

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_check_handle = .false.
       else
         bq_check_handle = .true.
       end if

       return
       end

      subroutine bq_print_info(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"

       integer handle
c      local variables
       integer i,printlevel
       integer k
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname
       logical status,oprint

       pname = "bq_print_info"

       call util_print_get_level(printlevel)
       oprint = (ga_nodeid().eq.0).and.(printlevel.le.print_medium)

       if(handle .lt.0 .or. handle .gt. max_bq) then
           call errquit(
     &      pname//' unable to locate handle ',
     &      0,0)
       end if

       if(bq_ncent(handle).eq.0) then
          write(*,*) "No charges are found"
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

       if (oprint) then
         call util_print_centered(6,
     >      "Bq Structure Information (Angstroms)",
     >              36, .true.)
      
         write(*,*) "Name: ", bq_name(handle)
         write(*,*) "Number of centers: ",bq_ncent(handle)
         do i=1,bq_ncent(handle)
           write(6,FMT=9000)
     >           i,(dbl_mb(i_c+3*(i-1)+k-1)*0.529177249d00,k=1,3),
     >           dbl_mb(i_q+i-1)

         end do
         write(*,*)
         write(*,*)

       end if
9000   format(i5,2x,"Bq",4x,3f15.8,3x,"charge",3x,f15.8)

       return
       end

      function bq_destroy(handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_destroy
c      local variables
       integer h_c
       integer h_q
       character*(32) pname

       pname = "bq_destroy"

       bq_destroy = .true.

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_destroy = .false.
         return
       else
         bq_destroy = .true.
       end if

       bq_name(handle)=" "
       bq_active(handle)=.false.

       if(bq_ncent(handle).eq.0) return
       bq_ncent(handle) = 0 

       if(.not.bq_mem(handle)) return
       bq_mem(handle) = .false.
        
       h_q = bq_charge(handle)
       h_c = bq_coord(handle)

       if(.not.ma_free_heap(h_q))
     &      call errquit(
     &      pname//' unable to deallocate heap space',
     &      0, MA_ERR)

       if(.not.ma_free_heap(h_c))
     &      call errquit(
     &      pname//' unable to deallocate heap space',
     &      0, MA_ERR)


       return
       end

      function bq_destroy_all()
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle
       logical bq_destroy_all
c      local variables
       character*(32) pname
        
       logical bq_destroy
       external bq_destroy

       pname = "bq_destroy_all"

       bq_destroy_all = .false.

       do handle=1,max_bq
         if(.not.bq_destroy(handle)) return
       end do
       bq_destroy_all = .true.

       return
       end

      function bq_index_coord(handle,i_c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
       integer handle
       logical bq_index_coord
c      local variables
       integer h_c,i_c
       character*(32) pname

       pname = "bq_index_coord"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_index_coord = .false.
         return
       else
         bq_index_coord = .true.
       end if

       h_c = bq_coord(handle)

       if(.not.ma_get_index( h_c, i_c) ) call errquit(
     &      pname//' uunable to locate coord handle',
     &      0, MA_ERR)

       return
       end

      function bq_index_charge(handle,i_q)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
       integer handle
       logical bq_index_charge
c      local variables
       integer h_q,i_q
       character*(32) pname

       pname = "bq_index_charge"

       if(handle .lt.0 .or. handle .gt. max_bq) then
         bq_index_charge = .false.
         return
       else
         bq_index_charge = .true.
       end if

       h_q = bq_charge(handle)

       if(.not.ma_get_index( h_q, i_q) ) call errquit(
     &      pname//' unable to locate charge handle',
     &      0, MA_ERR)

       return
       end

      function bq_on()
      implicit none
#include "bq_data.fh"
#include "errquit.fh"
       logical bq_on
c      local variables
       integer i

       bq_on = .false.

       do i=1,max_bq
        if(bq_active(i).and.bq_ncent(i).gt.0) then
          bq_on = .true.
          return
        end if
       end do
        
       return
       end

