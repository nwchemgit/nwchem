*
* $Id$
*

      block data initial_smd_data
      implicit none
#include "smd_data.fh"
c     
      data smd_n /0/
      data smd_nd /smd_max_n*0/

      end block data

      subroutine smd_data_get_sindex(sname,ns,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      integer ns
      logical result
c
      integer i

      result = .false.
      ns = 0
      do i=1,smd_n
       if(smd_name(i).eq.sname) then
         ns = i
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_get_dindex(sname,dname,nd,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer nd
      logical result
c
      integer i,ns
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_dindex"    
 
      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         nd = i
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_get_index(sname,dname,ind,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ind
      logical result
c
      integer i,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_index"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         ind = smd_index(ns,i)
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_destroy(sname,dname,result)
      implicit none
#include "smd_data.fh"
#include "mafdecls.fh"
      character*(*) sname
      character*(*) dname
      integer ind
      logical result
c
      integer i,j,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_index"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result) return

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         ind = smd_index(ns,i)
         result = .true.
         goto 10
       end if
      end do
10    continue
      if(.not.result) return
c
      if(.not.ma_free_heap(smd_handle(ns,i)))
     +     call errquit(pname//'memory deallocation'//smd_dname(ns,i),
     +           smd_dim(ns,i),
     +           0)

c
      do j=i,smd_nd(ns)-1
        smd_dname(ns,j)  = smd_dname(ns,j+1)
        smd_handle(ns,j) = smd_handle(ns,j+1)
        smd_index(ns,j)  = smd_index(ns,j+1)
        smd_type(ns,j)   = smd_type(ns,j+1)
        smd_dim(ns,j)    = smd_dim(ns,j+1)
      end do
      smd_dname(ns,smd_nd(ns))  = " "
      smd_handle(ns,smd_nd(ns)) = 0
      smd_index(ns,smd_nd(ns))  = 0
      smd_type(ns,smd_nd(ns))   = 0
      smd_dim(ns,smd_nd(ns))    = 0
      smd_nd(ns) = smd_nd(ns)-1 
      return
      end

      subroutine smd_data_get_size(sname,dname,nsize,result)
      implicit none
#include "smd_data.fh"
#include "smd_const_data.fh"
      character*(*) sname
      character*(*) dname
      integer nsize
      logical result
c
      integer i,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_size"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         if(smd_type(ns,i).eq.MT_STR) then
           nsize = smd_dim(ns,i)*smd_string_size
         else
           nsize = smd_dim(ns,i)
         end if
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_get_ind_dim(sname,dname,ind,nq,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ind,nq
      logical result
c
      integer i,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_index"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         ind = smd_index(ns,i)
         nq  = smd_dim(ns,i)
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_get_attrib(sname,dname,ind,nq,dtype,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ind,nq,dtype
      logical result
c
      integer i,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_index"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         ind = smd_index(ns,i)
         nq  = smd_dim(ns,i)
         dtype = smd_type(ns,i)
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_data_get_dim(sname,dname,ndim,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      logical result
c
      integer i,ns,nd
      character*72 buffer
      character*30 pname

      pname = "smd_data_get_dim"    
 
      buffer = sname

      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      result = .false.
      nd = 0
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         ndim = smd_dim(ns,i)
         result = .true.
         return
       end if
      end do

      return
      end

      subroutine smd_namespace_create(sname)
      implicit none
#include "smd_data.fh"
      character*(*) sname
c
      character*(30) pname
      pname = "smd_namespace_create"

      call smd_namespace_destroy(sname)
      smd_n = smd_n + 1
      if(smd_n.gt.smd_max_n) 
     +   call errquit(pname//'increase namespace size to',
     +         smd_n,
     +         0)
      smd_name(smd_n) = sname

      end

      subroutine smd_namespace_destroy(sname)
      implicit none
#include "smd_data.fh"
#include "errquit.fh"
#include "mafdecls.fh"
      character*(*) sname
c
      integer i,j
      integer ns
      logical result
      character*72 buffer
      character*(30) pname
      pname = "smd_namespace_destroy"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result) return

      do i=smd_nd(ns),1,-1
        if(.not.ma_free_heap(smd_handle(ns,i)))
     +   call errquit(pname//'memory deallocation'//smd_dname(ns,i),
     +         smd_dim(ns,i),
     +         MA_ERR)

        smd_dname(ns,i) = " "
        smd_dim(ns,i) = 0

      end do

      smd_nd(ns) = 0
      smd_name(ns) = " "

      do i=ns,smd_n-1
       smd_name(i) = smd_name(i+1)
       smd_nd(i)   = smd_nd(i+1)
       do j=1,smd_nd(i)
         smd_dname(i,j)  = " "
         smd_handle(i,j) = 0
         smd_index(i,j)  = 0
         smd_type(i,j)   = 0
         smd_dim(i,j)    = 0
       end do
       do j=1,smd_nd(i+1)
         smd_dname(i,j)  = smd_dname(i+1,j) 
         smd_handle(i,j) = smd_handle(i+1,j) 
         smd_index(i,j)  = smd_index(i+1,j) 
         smd_type(i,j)   = smd_type(i+1,j) 
         smd_dim(i,j)    = smd_dim(i+1,j) 
       end do
      end do
      smd_n = smd_n - 1
      end

      subroutine smd_data_destroy_all()
      implicit none
#include "smd_data.fh"
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer i
      integer ns
      character*(30) pname
      pname = "smd_data_destroy_all"

      do ns=smd_n,1,-1

        do i=smd_nd(ns),1,-1
          if(.not.ma_free_heap(smd_handle(ns,i)))
     +     call errquit(pname//'memory deallocation'//smd_dname(ns,i),
     +           smd_dim(ns,i),
     +           MA_ERR)

          smd_dname(ns,i) = " "
          smd_dim(ns,i) = 0

        end do

        smd_nd(ns) = 0
        smd_name(i) = " "
      end do
      smd_n = 0
     
      return
      end

      subroutine smd_data_print(un)
      implicit none
#include "smd_data.fh"
#include "errquit.fh"
#include "mafdecls.fh"
c
      integer i,un
      integer ns
      character*(30) pname
      pname = "smd_namespace_destroy"

      do ns=1,smd_n
        call smd_data_namespace_print(smd_name(ns),un)
      end do
     
      return
      end

      subroutine smd_data_create(sname,dname,ndim,dtype)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ma_type,ma_ndim
      integer ndim0,dtype0
      logical okeep,oexists

      call util_io_enable()
      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

c     see if the data already exists and is of the same nature
      call smd_data_inquire(sname,dname,ndim0,dtype0,oexists)
      if(oexists) 
     +  call smd_data_get_dindex(sname,dname,nd,result)
      okeep = oexists.and.
     +        ndim.eq.ndim0.and.
     +        dtype.eq.dtype0 
      if(okeep) goto 10
      
      if(oexists) then
      if(.not.ma_free_heap(smd_handle(ns,nd)))
     +     call errquit(pname//'memory deallocation'//smd_dname(ns,nd),
     +           smd_dim(ns,nd),
     +           0)
      else 
        smd_nd(ns) = smd_nd(ns) + 1
        if(smd_nd(ns).gt.smd_max_nd) 
     +     call errquit(pname//'increase dataspace size to',
     +           smd_nd(ns),
     +           0)

        nd = smd_nd(ns)
      end if
      smd_dim(ns,nd)=ndim
      smd_type(ns,nd)=dtype
      smd_dname(ns,nd)=dname

      if(dtype.eq.MT_STR) then
        ma_type = MT_BYTE
        ma_ndim = ndim*smd_string_size
      else
        ma_type = dtype
        ma_ndim = ndim
      end if

      if(.not.ma_alloc_get(ma_type,
     +   ma_ndim,
     +   smd_dname(ns,nd),
     +   smd_handle(ns,nd),
     +   smd_index(ns,nd)))
     +   call errquit(pname//'memory allocation'//smd_dname(ns,nd),
     +         smd_dim(ns,nd),
     +         MA_ERR)

10    continue      
      if(smd_type(ns,nd).eq.MT_DBL) then
        call dfill(smd_dim(ns,nd),
     +             0.d0,
     +             dbl_mb(smd_index(ns,nd)),
     +             1)
      end if

      if(smd_type(ns,nd).eq.MT_INT) then
        call ifill(smd_dim(ns,nd),
     +             0,
     +             int_mb(smd_index(ns,nd)),
     +             1)
      end if
      if(smd_type(ns,nd).eq.MT_LOG) then
        do i=1,smd_dim(ns,nd)
          log_mb(smd_index(ns,nd)+i-1) = .false.
        end do
      end if
      if(smd_type(ns,nd).eq.MT_STR) then
        call smd_wipe_ma_char(smd_dim(ns,nd),
     +                        byte_mb(smd_index(ns,nd)))
      end if
cc     should probably fill logical and char arrays
      call util_io_disable()

      end

      subroutine smd_data_create_get(sname,dname,ndim,dtype,ind)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
      integer ind
c
      logical result
      character*30 pname,buffer

      pname = "smd_data_create_get"

      call smd_data_create(sname,dname,ndim,dtype)
      call smd_data_get_index(sname,dname,ind,result)
      if(.not.result) 
     +  call errquit(pname//" error 1",0,0)
      end

      subroutine smd_data_create_get0(sname,dname,ndim,dtype,ind)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
      integer ind
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ma_type,ma_ndim
      integer dtype0,ndim0
      logical okeep,oexists

      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

c     see if the data already exists and is of the same nature
      call smd_data_inquire(sname,dname,ndim0,dtype0,result)
      okeep = result.and.
     +        ndim.eq.ndim0.and.
     +        dtype.eq.dtype0 
      if(okeep) then
        call smd_data_get_dindex(sname,dname,nd,result)
        if(.not.result) 
     +   call errquit(pname//"no index",0,0)
        write(*,*) "found same data"
      else
        write(*,*) "creating new data"
        call smd_data_destroy(sname,dname,result)
        smd_nd(ns) = smd_nd(ns) + 1
        if(smd_nd(ns).gt.smd_max_nd) 
     +     call errquit(pname//'increase dataspace size to',
     +           smd_nd(ns),
     +           0)

        nd = smd_nd(ns)
        smd_dim(ns,nd)=ndim
        smd_type(ns,nd)=dtype
        smd_dname(ns,nd)=dname

        if(dtype.eq.MT_STR) then
          ma_type = MT_BYTE
          ma_ndim = ndim*smd_string_size
        else
          ma_type = dtype
          ma_ndim = ndim
        end if

        if(.not.ma_alloc_get(ma_type,
     +     ma_ndim,
     +     smd_dname(ns,nd),
     +     smd_handle(ns,nd),
     +     smd_index(ns,nd)))
     +     call errquit(pname//'memory allocation'//smd_dname(ns,nd),
     +           smd_dim(ns,nd),
     +           MA_ERR)


      end if
      ind = smd_index(ns,nd)
c
      if(smd_type(ns,nd).eq.MT_DBL) then
        call dfill(smd_dim(ns,nd),
     +             0.d0,
     +             dbl_mb(smd_index(ns,nd)),
     +             1)
      end if

      if(smd_type(ns,nd).eq.MT_INT) then
        call ifill(smd_dim(ns,nd),
     +             0,
     +             int_mb(smd_index(ns,nd)),
     +             1)
      end if
      if(smd_type(ns,nd).eq.MT_LOG) then
        do i=1,smd_dim(ns,nd)
          log_mb(smd_index(ns,nd)+i-1) = .false.
        end do
      end if
cc     should probably fill logical and char arrays

      end

      subroutine smd_data_inquire(sname,dname,ndim,dtype,result)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
      logical result
c
      integer i,ns
      character*72 buffer
      character*30 pname

      pname = "smd_data_inquire"    
 
      buffer = sname

      ndim  = 0
      dtype = 0 
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result) goto 100

      result = .false.
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         result = .true.
         ndim = smd_dim(ns,i)
         dtype=smd_type(ns,i)
       end if
      end do

100   continue
      return
      end

      subroutine smd_data_inquire_type(sname,dname,dtype)
      implicit none
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer dtype
c      
      integer ndim
      logical result
      integer i,ns
      character*72 buffer
      character*30 pname

      pname = "smd_data_inquire_type"    
 
      buffer = sname

      ndim  = 0
      dtype = 0 
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result) goto 100

      result = .false.
      buffer = dname
      do i=1,smd_nd(ns)
       if(smd_dname(ns,i).eq.buffer) then
         result = .true.
         dtype=smd_type(ns,i)
         goto 100
       end if
      end do

100   continue
      if(.not.result)
     +   call errquit(pname//'data not found'//buffer,
     +         0,
     +         0)
      return
      end

      subroutine smd_data_namespace_brdcst(sname,root)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "msgids.fh"
#include "global.fh"
#include "smd_const_data.fh"
      character*(*) sname
      integer root
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim,dtype,ind

      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      nd = smd_nd(ns)
      do i=1,nd
        ndim  = smd_dim(ns,i)
        dtype = smd_type(ns,i)
        ind   = smd_index(ns,i) 
        buffer = smd_dname(ns,i)
        if(dtype.eq. MT_DBL) then
          call ga_brdcst(msg_smd,dbl_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_INT) then
          call ga_brdcst(msg_smd,int_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
c could be problems here!!!!
        else if(dtype.eq.MT_STR) then
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     smd_string_size*ndim*ma_sizeof(MT_BYTE,1,mt_byte),root)
         else if(dtype.eq.MT_BYTE) then
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_LOG) then
          call ga_brdcst(msg_smd,log_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else
         call errquit(pname//'unknown type',0,0)
        endif
      end do 

      return
      end

      subroutine smd_data_namespace_print(sname,un)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "msgids.fh"
#include "global.fh"
#include "smd_const_data.fh"
      character*(*) sname
      integer un
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,j,k
      integer ndim,dtype,ind

      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      nd = smd_nd(ns)
      write(un,*) "NAMESPACE",sname
      do k=1,nd
        ndim  = smd_dim(ns,k)
        dtype = smd_type(ns,k)
        ind   = smd_index(ns,k) 
        buffer = smd_dname(ns,k)
        if(dtype.eq. MT_DBL) then
          write(un,*) "DATA NAME ",buffer
          write(un,*) "DATA TYPE"," double"
          write(un,*) "DATA SIZE",ndim
          do i=1,ndim
          write(un,*) dbl_mb(ind+i-1)
          end do
        else if(dtype.eq.MT_LOG) then
          write(un,*) "DATA NAME ",buffer
          write(un,*) "DATA TYPE"," logical"
          write(un,*) "DATA SIZE",ndim
          do i=1,ndim
          write(un,*) log_mb(ind+i-1)
          end do
         else if(dtype.eq.MT_INT) then
          write(un,*) "DATA NAME ",buffer
          write(un,*) "DATA TYPE"," integer"
          write(un,*) "DATA SIZE",ndim
          do i=1,ndim
          write(un,*) int_mb(ind+i-1)
          end do
        else if(dtype.eq.MT_STR) then
          write(un,*) "DATA NAME ",buffer
          write(un,*) "DATA TYPE"," string"
          write(un,*) "DATA SIZE",ndim
c          call smd_ma_char_print(ndim*smd_string_size,byte_mb(ind),un)
          do i=1,ndim
          write(un,'(72A1)') 
     >      (byte_mb(ind+smd_string_size*(i-1)+j-1),j=1,smd_string_size)
          end do
        else
         call errquit(pname//'unknown type',0,0)
        endif
      end do 

      return
      end

      subroutine smd_data_namespace_brdcst_create(sname,root)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "msgids.fh"
#include "global.fh"
#include "smd_const_data.fh"
      character*(*) sname
      integer root
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim,dtype,ind

      pname = "smd_data_create"

      if(ga_nodeid().eq.root) then
      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)
       nd = smd_nd(ns)
      else
       nd = 0
      end if

      if(ga_nodeid().ne.root) then
        call smd_namespace_destroy(sname)
        call smd_namespace_create(sname)
      end if
c
      call ga_brdcst(msg_smd,nd,
     >     ma_sizeof(mt_int,1,mt_byte),root)
c
      call ga_sync()
      do i=1,nd
        if(ga_nodeid().eq.root) then
          ndim  = smd_dim(ns,i)
          dtype = smd_type(ns,i)
          buffer = smd_dname(ns,i)
        end if

        call ga_brdcst(msg_smd,ndim,
     >       ma_sizeof(mt_int,1,mt_byte),root)
        call ga_brdcst(msg_smd,dtype,
     >       ma_sizeof(mt_int,1,mt_byte),root)
        call ga_brdcst(msg_smd,buffer,
     >     72*ma_sizeof(mt_byte,1,mt_byte),root)
     
        call ga_sync()
     
        if(ga_nodeid().ne.root) then
          call smd_data_create_get(sname,buffer,ndim,dtype,ind)
        else
          ind = smd_index(ns,i)
        end if
        call ga_sync()

        if(dtype.eq. MT_DBL) then
          call ga_brdcst(msg_smd,dbl_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_INT) then
          call ga_brdcst(msg_smd,int_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_STR) then
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     smd_string_size*ndim*ma_sizeof(MT_BYTE,1,mt_byte),root)
       else if(dtype.eq.MT_BYTE) then
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_LOG) then
          call ga_brdcst(msg_smd,log_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else
         call errquit(pname//'unknown type',0,0)
        endif

      end do 

      return
      end

      subroutine smd_data_brdcst_create(root)
      implicit none
#include "smd_data.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      integer root
c
      integer ns,nsmax
      character*(30) pname
      character*(72) buffer
      pname = "smd_namespace_destroy"

      
      if(ga_nodeid().eq.root) then
        nsmax = smd_n
      end if
      call ga_brdcst(msg_smd,nsmax,
     >     ma_sizeof(mt_int,1,mt_byte),root)
      do ns=1,nsmax
      if(ga_nodeid().eq.root) then
        buffer = smd_name(ns)
      end if
      call ga_sync()
      call ga_brdcst(msg_smd,buffer,
     >     72*ma_sizeof(mt_byte,1,mt_byte),root)
       call smd_data_namespace_brdcst_create(buffer,root)   
      end do
     
      return
      end

      subroutine smd_data_put_log(sname,dname,ndim,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      logical datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i

      pname = "smd_data_put_int"

      call smd_data_create_get(sname,dname,ndim,MT_LOG,ind)
      do i=1,ndim
        log_mb(ind+i-1) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_log(sname,dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      logical datum(ndim)
c
      integer ind,dtype
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_get_int"

      call smd_data_get_attrib(sname,dname,ind,ndim1,dtype,result)
      if(.not.result) return
      if(dtype.ne.MT_LOG)
     >     call errquit(pname//"incompatible types",0,0)
      if(ndim1.ne.ndim)
     >     call errquit(pname//"incompatible sizes",0,0)

      do i=1,ndim
        datum(i) = log_mb(ind+i-1)
      end do
      return
      end

      subroutine smd_data_put_int1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      integer datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_int1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
       do i=1,ndim
        int_mb(ind+ilo+i-2) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_int1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      integer datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_int1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
       do i=1,ndim
        datum(i) = int_mb(ind+ilo+i-2)
      end do
      return
      end

      subroutine smd_data_put_dbl1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      double precision datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_dbl1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      do i=1,ndim
        dbl_mb(ind+ilo+i-2) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_dbl1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      double precision datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_dbl1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
       do i=1,ndim
        datum(i) = dbl_mb(ind+ilo+i-2)
      end do
      return
      end

      subroutine smd_data_put_log1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      logical datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_int1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
       do i=1,ndim
        log_mb(ind+ilo+i-2) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_log1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      logical datum(ihi-ilo+1)
c
      integer ndim
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_int1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      do i=1,ndim
        datum(i) = log_mb(ind+ilo+i-2) 
      end do
      return
      end

      subroutine smd_data_get_log2(sname,dname,ilo,ihi,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      logical datum(ihi-ilo+1)
      logical result
c
      integer ndim
      integer ind
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_put_int1"

      ndim = ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           return
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      do i=1,ndim
        datum(i) = log_mb(ind+ilo+i-2) 
      end do
      return
      end

      subroutine smd_data_put_int(sname,dname,ndim,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i

      pname = "smd_data_put_int"

      call smd_data_create_get(sname,dname,ndim,MT_INT,ind)
      do i=1,ndim
        int_mb(ind+i-1) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_int(sname,dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer datum(ndim)
c
      integer ind,dtype
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_get_int"

      call smd_data_get_attrib(sname,dname,ind,ndim1,dtype,result)
      if(.not.result) return
      if(dtype.ne.MT_INT)
     >     call errquit(pname//"incompatible types",0,0)
      if(ndim1.ne.ndim)
     >     call errquit(pname//"incompatible sizes",0,0)

      do i=1,ndim
        datum(i) = int_mb(ind+i-1)
      end do
      return
      end

      subroutine smd_data_put_dbl(sname,dname,ndim,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      double precision datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i

      pname = "smd_data_put_int"

      call smd_data_create_get(sname,dname,ndim,MT_DBL,ind)
      do i=1,ndim
        dbl_mb(ind+i-1) = datum(i)
      end do
      return
      end

      subroutine smd_data_get_dbl(sname,dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      double precision datum(ndim)
c
      integer ind,dtype
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_get_int"

      call smd_data_get_attrib(sname,dname,ind,ndim1,dtype,result)
      if(.not.result) return
      if(dtype.ne.MT_DBL)
     >     call errquit(pname//"incompatible types",0,0)

      if(ndim1.ne.ndim)
     >     call errquit(pname//"incompatible sizes",0,0)

      do i=1,ndim
        datum(i) = dbl_mb(ind+i-1)
      end do
      return
      end

      subroutine smd_data_put_string(sname,dname,ndim,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
      character*(*) datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,l
      integer ndim1

      pname = "smd_data_create"

c      ndim1 = smd_tag_size*ndim

      do i=1,ndim
       l = inp_strlen(datum(i))
       if(l.gt.smd_string_size)
     >     call errquit(pname//"too big size",l,0) 
      end do
      call smd_data_create_get(sname,dname,ndim,MT_STR,ind)
      call smd_set_char_to_ma(ndim,datum,byte_mb(ind))
      return
      end

      subroutine smd_data_put_string1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      integer dtype
      character*(*) datum(ihi-ilo+1)
c
      integer ind
      integer ndim
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,l
      integer ndim1
      integer offset

      pname = "smd_data_create"

c      ndim1 = smd_tag_size*ndim

      ndim=ihi-ilo+1
      do i=1,ndim
       l = inp_strlen(datum(i))
       if(l.gt.smd_string_size)
     >     call errquit(pname//"too big size",l,0) 
      end do
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      offset = smd_string_size*(ilo-1)
      call smd_set_char_to_ma1(datum,ilo,ihi,ndim1,byte_mb(ind))
      return
      end

      subroutine smd_data_get_string1(sname,dname,ilo,ihi,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      integer dtype
      character*(*) datum(ihi-ilo+1)
c
      integer ind
      integer ndim
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,l
      integer ndim1
      integer offset

      pname = "smd_data_get_string1"

      ndim=ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           buffer = dname
           call errquit(pname//"cannot find"//buffer,0,0) 
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      offset = smd_string_size*(ilo-1)
      call smd_set_ma_to_char1(datum,ilo,ihi,ndim1,byte_mb(ind))
      return
      end

      subroutine smd_data_get_string2(sname,dname,ilo,ihi,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) sname
      character*(*) dname
      integer ilo,ihi
      integer dtype
      character*(*) datum(ihi-ilo+1)
c
      integer ind
      integer ndim
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i,l
      integer ndim1
      integer offset

      pname = "smd_data_get_string1"

      ndim=ihi-ilo+1
      call smd_data_get_ind_dim(sname,dname,ind,ndim1,result)
      if(.not.result) then
           return
      end if
      if(ihi.gt.ndim1) then
          buffer = dname
          call errquit(pname//"out of bounds"//buffer,0,0)
      end if
      offset = smd_string_size*(ilo-1)
      call smd_set_ma_to_char1(datum,ilo,ihi,ndim1,byte_mb(ind))
      return
      end

      subroutine smd_data_get_string(sname,dname,ndim,datum,result)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "smd_const_data.fh"
#include "inp.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      character*(*) datum(ndim)
c
      integer ind
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1,dtype
      integer chunk
      logical result

      pname = "smd_data_get_string"

      call smd_data_get_attrib(sname,dname,ind,ndim1,dtype,result)
      if(.not.result) return
c      write(*,*) "LEN",LEN(datum(1))
      if(dtype.ne.MT_STR)
     >     call errquit(pname//"incompatible types",0,0)
      if(ndim1.ne.ndim)
     >     call errquit(pname//"incompatible sizes",0,0)
      call smd_set_ma_to_char(ndim1,byte_mb(ind),datum)
      return
      end

      subroutine smd_data_namespace_rtdb_get(rtdb,sname,prefix)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
#include "msgids.fh"
#include "global.fh"
#include "smd_const_data.fh"
#include "inp.fh"
#include "rtdb.fh"
      integer rtdb
      character*(*) sname
      character*(*) prefix
c
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      character*256 big_buffer
      integer i
      integer istart,iend
      integer ndim,dtype,ind
      character*(smd_string_size) token
      character*256 mtoken(100)
      integer  itoken(100)
      double precision  ftoken(100)
      logical ltoken(100)
      pname = "smd_data_namespace_rtdb_get"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)
      nd = smd_nd(ns)
     
      if(.not.rtdb_first(rtdb,token))
     +  call errquit(pname//'rtdb is empty',0,0)
c
300   continue
      if(index(token,prefix).gt.0) then
      if(.not.rtdb_get_info(rtdb, token, dtype, ndim, buffer))
     >  call errquit("rtdb_get_info"//token,0,0)

      if(dtype.eq. MT_DBL) then
        if(.not.rtdb_get(rtdb,token,mt_dbl,ndim,ftoken))
     >     call errquit("rtdb_get"//token,0,0)
        call smd_data_put_dbl(sname,token,ndim,ftoken)
      else if(dtype.eq.MT_INT) then
        if(.not.rtdb_get(rtdb,token,mt_int,ndim,itoken))
     >     call errquit("rtdb_get"//token,0,0)
        call smd_data_put_int(sname,token,ndim,itoken)
      else if(dtype.eq.MT_C_CHAR) then
        if(.not.rtdb_cget_size(rtdb,token,ndim))
     >     call errquit("rtdb_get_size"//token,0,0)
        if(.not.rtdb_cget(rtdb,token,ndim,mtoken))
     >     call errquit("rtdb_get"//token,0,0)
        call smd_data_put_string(sname,token,ndim,mtoken)
      else if(dtype.eq.MT_LOG) then
        if(.not.rtdb_get(rtdb,token,mt_log,ndim,ltoken))
     >     call errquit("rtdb_get"//token,0,0)
        call smd_data_put_log(sname,token,ndim,ltoken)
      else
        call errquit(pname//'unknown type',0,0)
      endif

      end if

      result = rtdb_next(rtdb,big_buffer)
      token = big_buffer
      if(result) goto 300
c      call smd_data_put_int(sname,"smd:rtdb_handle",1,rtdb)

c      if(.not. rtdb_getfname(rtdb, big_buffer)) call
c     *        errquit('rtdb_getfname failed',0,0)
c      istart = 0
c   10 if (inp_strtok(big_buffer, '/', istart, iend)) then
c      buffer = big_buffer(istart:iend)
c      goto 10
c      endif
c      call smd_data_put_string(sname,"smd:rtdb_name",1,
c     >                         buffer)

      return
      end

      subroutine smd_data_test()
C$Id$
      implicit none
c
#include "mafdecls.fh"
#include "smd_const_data.fh"
c
      integer rtdb, ma_handle, ma_index
      integer itest(3), ibuf(3)
      integer itest1(4)
      real ftest(4), fbuf(4)
      double precision dtest(4), dbuf(5)
      double precision dtest1(5)
      character *20 cbuf(4), ccbuf(4)
      character *20 name,rtdb_fname
      character *26 date
      logical status
      integer type, nelem, i
      character*(smd_string_size) sname,tag
c
      data itest/1, 2, 3/
      data itest1/10, 20, 30,40/
      data ftest/1.0e0, 2.0e0, 3.0e0, 4.0e0/
      data dtest/1.0d0, 2.0d0, 3.0d0, 4.0d0/
      data dtest1/10.0d0, 20.0d0, 30.0d0, 40.0d0, 50.0d0/
c
      cbuf(1) = 'Have'
      cbuf(2) = 'a'
      cbuf(3) = 'nice'
      cbuf(4) = 'day, Robert!'
c
      write(*,*) "Doing a test aren't you :)"
      sname="marat"

      call smd_namespace_create(sname)
      tag = "itest"
      call smd_data_put_int(sname,tag,3,itest)
c      
      tag = "dtest"
      call smd_data_put_dbl(sname,tag,4,dtest)
c
      call smd_data_namespace_print(sname,6) 

      tag = "itest"
      call smd_data_put_int(sname,tag,4,itest1)
      call smd_data_namespace_print(sname,6) 
      end
