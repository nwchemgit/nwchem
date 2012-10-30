c $Id$
      block data initial_bq_data
      implicit none
#include "bq_data.fh"
c     
      data bq_active /max_bq*.false./
      data bq_mem /max_bq*.false./
      data bq_ncent /max_bq*-1/

      end
C>
C> \defgroup bq Point charges
C>
C> The Bq module plays a role that is similar to the role of the
C> geometry. I.e. it deals with centers distributed in space where
C> the centers may have some properties. 
C>
C> Because Bq centers are often used to generate embedding potentials
C> in QM/MM calculations Bq instances have to be able to contain
C> thousands of centers. The geometry instance typically contain only
C> atoms and hence are restricted in the number of centers they can
C> contain.
C
C> \ingroup bq
C> @{
C>
C> \brief Create a Bq instance
C>
      function bq_create(namespace,handle)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       character*(*) namespace !< [Input] The Bq instance name
       logical bq_create
c      local variables
       integer i
       integer handle !< [Output] The Bq instance handle

       bq_create = .false.
       do i=1,max_bq
        if(bq_ncent(i).eq.-1) then
          bq_create = .true.
          bq_name(i) = namespace
          bq_ncent(i) = 0
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

       logical bq_check_handle
       external bq_check_handle

       bq_get_handle = .false.
       do i=1,max_bq
        if(bq_name(i).eq.namespace.and.bq_check_handle(i)) then
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

       logical bq_check_handle
       external bq_check_handle

       bq_activate = .true.
       if(.not.bq_check_handle(handle)) then
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

       logical bq_check_handle
       external bq_check_handle

       bq_deactivate = .true.
       if(.not.bq_check_handle(handle)) then
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

       logical bq_check_handle
       external bq_check_handle

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
C>
C> \brief Set the coordinates and charges for a Bq instance
C>
C> Allocates memory and associates it with the Bq instance and 
C> initializes that memory with the charges and coordinates provided.
C> The Bq instance assumes responsibility for the memory. I.e.
C> the memory will be deallocated when the Bq instance is eventually
C> destroyed.
C>
      function bq_set(handle,n,q,c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle !< [Input] The Bq instance handle
       integer n      !< [Input] The number of centers
       double precision q(n)   !< [Input] The charges
       double precision c(3*n) !< [Input] The coordinates
       logical bq_set
c      local variables
       integer i
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_set"

       if(.not.bq_check_handle(handle)) then
         bq_set = .false.
         return
       else
         bq_set = .true.
       end if
c
c      If there is
c        - already memory associated with this Bq instance
c        - but it is not enough to hold the new data and
c        - the Bq instance is responsible for this memory 
c      then 
c        free the memory before allocating new memory to avoid
c        memory leaks
c
       if (bq_ncent(handle).ne.0.and.bq_ncent(handle).lt.n) then
         h_c = bq_coord(handle)
         h_q = bq_charge(handle)
         if (bq_mem(handle)) then
           if (.not.ma_free_heap(h_c))
     &       call errquit(pname//' unable to free heap space',
     &                    h_c,MA_ERR)
           if (.not.ma_free_heap(h_q))
     &       call errquit(pname//' unable to free heap space',
     &                    h_q,MA_ERR)
         endif
         bq_ncent(handle) = 0
       endif
c
c      If this Bq instance has no memory then
c        allocate some
c      else
c        look up the offsets
c
       if (bq_ncent(handle).eq.0) then
         if(.not.ma_alloc_get(MT_DBL, 3*n, 'bqdata c',
     &        h_c, i_c) ) call errquit(
     &        pname//' unable to allocate heap space',
     &        3*n, MA_ERR)

         if(.not.ma_alloc_get(MT_DBL, n, 'bqdata q',
     &        h_q, i_q) ) call errquit(
     &        pname//' unable to allocate heap space',
     &        n, MA_ERR)
       else
         h_c = bq_coord(handle)
         h_q = bq_charge(handle)
         if(.not.ma_get_index( h_c, i_c) ) call errquit( 
     &        pname//' unable to locate coord handle', 
     &        0, MA_ERR)

         if(.not.ma_get_index( h_q, i_q) ) call errquit(
     &        pname//' unable to locate charge handle',
     &        0, MA_ERR)
       endif

      
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
C>
C> \brief Allocate and initialize space for a Bq instance
C>
C> This routine always allocates and initializes memory for this Bq
C> instance. If this instance already has memory associated with it
C> and if this memory is the responsibility of the Bq instance it will
C> deallocated.
C>
C> \return Return .true. if successful and .false. otherwise.
C>
      function bq_alloc(handle,n)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle !< [Input] The Bq instance handle
       integer n      !< [Input] The number of centers
       logical bq_alloc
c      local variables
       integer i
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_alloc"

       if(.not.bq_check_handle(handle)) then
         bq_alloc = .false.
         return
       else
         bq_alloc = .true.
       end if

       if(bq_ncent(handle).ne.0.and.bq_mem(handle)) then
         h_c = bq_coord(handle)
         h_q = bq_charge(handle)
         if (.not.ma_free_heap(h_c))
     &     call errquit(pname//' unable to free heap space',
     &                  h_c,MA_ERR)
         if (.not.ma_free_heap(h_q))
     &     call errquit(pname//' unable to free heap space',
     &                  h_q,MA_ERR)
         bq_ncent(handle) = 0
       endif

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
C>
C> \brief Associate memory handles with a Bq instance
C>
C> This function associates chunks of memory containing the charges
C> and the coordinates with a Bq instance. The memory is supposed to
C> be allocated on the heap in the calling routine. The memory remains
C> the responsibility of the application. I.e. if the Bq instance is
C> destroyed the memory chunks will not be deallocated. 
C> To set the memory chunks and transfer the associated responsibility
C> to the Bq instance use bq_pset_mem instead.
C>
C> Any memory that was associated with this Bq instance and was the
C> responsibility of the Bq instance will be deallocated before the
C> new memory gets associated.
C>
C> \return Returns .true. if successfull and .false. otherwise.
C>
      function bq_pset(handle,n,h_q,h_c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle !< [Input] The Bq instance handle
       integer n      !< [Input] The number of centers
       integer h_q    !< [Input] The memory handle for the charges
       integer h_c    !< [Input] The memory handle for the coordinates
       logical bq_pset
c      local variables
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_pset"

       if(.not.bq_check_handle(handle)) then
         bq_pset = .false.
         return
       else
         bq_pset = .true.
       end if

       if(bq_ncent(handle).ne.0.and.bq_mem(handle)) then
         h_c = bq_coord(handle)
         h_q = bq_charge(handle)
         if (.not.ma_free_heap(h_c))
     &     call errquit(pname//' unable to free heap space',
     &                  h_c,MA_ERR)
         if (.not.ma_free_heap(h_q))
     &     call errquit(pname//' unable to free heap space',
     &                  h_q,MA_ERR)
         bq_ncent(handle) = 0
       endif

       bq_ncent(handle)  = n
       bq_charge(handle) = h_q
       bq_coord(handle)  = h_c
       bq_mem(handle)    = .false.

       return
       end
C>
C> \brief Transfer memory handles and associated responsibility to a Bq
C> instance
C>
C> This function associates chunks of memory containing the charges
C> and the coordinates with a Bq instance, and transfers the 
C> responsibility for managing this memory as well.
C> The memory is supposed to be allocated on the heap in the calling
C> routine.
C> To set the memory chunks and not transfer the associated
C> responsibility to the Bq instance use bq_pset instead.
C>
C> Any memory that was associated with this Bq instance and was the
C> responsibility of the Bq instance will be deallocated before the
C> new memory gets associated.
C>
C> \return Returns .true. if successfull and .false. otherwise.
C>
      function bq_pset_mem(handle,n,h_q,h_c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
       integer handle !< [Input] The Bq instance handle
       integer n      !< [Input] The number of centers
       integer h_q    !< [Input] The memory handle for the charges
       integer h_c    !< [Input] The memory handle for the coordinates
       logical bq_pset_mem
c      local variables
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_pset_mem"

       if(.not.bq_check_handle(handle)) then
         bq_pset_mem = .false.
         return
       else
         bq_pset_mem = .true.
       end if

       if(bq_ncent(handle).ne.0.and.bq_mem(handle)) then
         h_c = bq_coord(handle)
         h_q = bq_charge(handle)
         if (.not.ma_free_heap(h_c))
     &     call errquit(pname//' unable to free heap space',
     &                  h_c,MA_ERR)
         if (.not.ma_free_heap(h_q))
     &     call errquit(pname//' unable to free heap space',
     &                  h_q,MA_ERR)
         bq_ncent(handle) = 0
       endif

       bq_ncent(handle)  = n
       bq_charge(handle) = h_q
       bq_coord(handle)  = h_c
       bq_mem(handle)    = .true.

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

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_get"

       if(.not.bq_check_handle(handle)) then
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
C>
C> \brief Get the number of centers of a Bq instance
C>
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

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_ncenter"

       if(.not.bq_check_handle(handle)) then
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

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_namespace"

       if(.not.bq_check_handle(handle)) then
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

       if(handle .lt.1 .or. handle .gt. max_bq) then
         bq_check_handle = .false.
       else if (bq_ncent(handle).eq.-1) then
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
       double precision bq_charge_total

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_print_info"

       call util_print_get_level(printlevel)
       oprint = (ga_nodeid().eq.0).and.(printlevel.le.print_medium)

       if(.not.bq_check_handle(handle)) then
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

c        == tally up bq charges ==
         bq_charge_total = 0.d0
         do i=1,bq_ncent(handle)
           write(6,FMT=9000)
c     >           i,(dbl_mb(i_c+3*(i-1)+k-1),k=1,3),
     >           i,(dbl_mb(i_c+3*(i-1)+k-1)*0.529177249d00,k=1,3),
     >           dbl_mb(i_q+i-1)
           bq_charge_total = bq_charge_total + dbl_mb(i_q+i-1)
         end do
         write(*,*) "Total Bq charge: ",bq_charge_total
         write(*,*)
         write(*,*)
c
       end if
9000   format(i5,2x,"Bq",4x,3f15.8,3x,"charge",3x,f15.8)

       return
       end
C>
C> \brief Destroy a Bq instance
C>
C> Destroys a Bq instance. Whether the memory for the coordinates
C> and the charges is deallocated depends on how this memory was
C> associated with the Bq instance.
C>
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

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_destroy"

       bq_destroy = .true.

       if(.not.bq_check_handle(handle)) then
         bq_destroy = .false.
         return
       else
         bq_destroy = .true.
       end if

       bq_name(handle)=" "
       bq_active(handle)=.false.

       if(bq_ncent(handle).eq.0) then
         bq_ncent(handle) = -1
         return
       else
         bq_ncent(handle) = -1
       endif

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
C>
C> \brief Destroy all Bq instances
C>
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

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_destroy_all"

       bq_destroy_all = .false.

       do handle=1,max_bq
         if(bq_check_handle(handle)) then
           if(.not.bq_destroy(handle)) return
         endif
       end do
       bq_destroy_all = .true.

       return
       end
C>
C> \brief Retrieve the memory index for the coordinates of a Bq instance
C>
      function bq_index_coord(handle,i_c)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
       integer handle !< [Input] The Bq instance handle
       logical bq_index_coord
c      local variables
       integer h_c
       integer i_c !< [Output] The coordinates memory index
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_index_coord"

       if(.not.bq_check_handle(handle)) then
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
C>
C> \brief Retrieve the memory index for the charges of a Bq instance
C>
      function bq_index_charge(handle,i_q)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
       integer handle !< [Input] The Bq instance handle
       logical bq_index_charge
c      local variables
       integer h_q
       integer i_q !< [Output] The charges memory index
       character*(32) pname

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_index_charge"

       if(.not.bq_check_handle(handle)) then
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

      subroutine bq_force_status(rtdb,ostatus)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "rtdb.fh"
#include "errquit.fh"
       integer rtdb
       logical ostatus
c       
       external bq_on
       logical bq_on
c      local variables
       integer i

       logical bq_check_handle
       external bq_check_handle

       ostatus = .false.

       if(.not.bq_on()) return
       if(.not. rtdb_get(rtdb,"bq:force",mt_log,1,ostatus)) 
     >        ostatus = .false.
        
       return
       end


      function bq_on()
      implicit none
#include "bq_data.fh"
#include "errquit.fh"
       logical bq_on
c      local variables
       integer i

       logical bq_check_handle
       external bq_check_handle

       bq_on = .false.

       do i=1,max_bq
        if(bq_active(i).and.bq_ncent(i).gt.0) then
          bq_on = .true.
          return
        end if
       end do
        
       return
       end

      subroutine bq_print_xyz(handle,un)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"

       integer handle
       integer un
c      local variables
       integer i,printlevel
       integer k
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname
       logical status,oprint

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_print_info"

       call util_print_get_level(printlevel)
       oprint = (ga_nodeid().eq.0).and.(printlevel.le.print_medium)

       if(.not.bq_check_handle(handle)) then
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

         write(un,*) bq_ncent(handle)
         write(un,*)
         do i=1,bq_ncent(handle)
           write(un,FMT=9000)
     >           i,(dbl_mb(i_c+3*(i-1)+k-1)*0.529177249d00,k=1,3),
     >              dbl_mb(i_q+i-1)

         end do
         write(*,*)
         write(*,*)

9000   format(i5,2x,"Bq",4x,4f15.8)

       return
       end

      subroutine bq_print_xyzq(handle,un)
      implicit none
#include "mafdecls.fh"
#include "bq_data.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"

       integer handle
       integer un
c      local variables
       integer i,printlevel
       integer k
       integer h_c,i_c
       integer h_q,i_q
       character*(32) pname
       logical status,oprint

       logical bq_check_handle
       external bq_check_handle

       pname = "bq_print_info"

       call util_print_get_level(printlevel)
       oprint = (ga_nodeid().eq.0).and.(printlevel.le.print_medium)

       if(.not.bq_check_handle(handle)) then
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

         write(un,*) bq_ncent(handle)
         write(un,*)
         do i=1,bq_ncent(handle)
           write(un,FMT=9000)
     >           i,(dbl_mb(i_c+3*(i-1)+k-1)*0.529177249d00,k=1,3),
     >             dbl_mb(i_q+i-1)

         end do
         write(*,*)
         write(*,*)

9000   format(i5,2x,"Bq",4x,4f15.8)

       return
       end
C> @}
