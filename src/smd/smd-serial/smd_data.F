*
* $Id: smd_data.F,v 1.9 2008-05-05 21:15:59 marat Exp $
*

      block data initial_smd_data
      implicit none
#include "smd_data.fh"
c     
      data smd_n /0/
      data smd_nd /smd_max_n*0/

      end

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

      subroutine smd_data_get_index_size(sname,dname,ind,nq,result)
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

      subroutine smd_data_get_size(sname,dname,ndim,result)
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
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

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
      pname = "smd_namespace_destroy"

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

        write(un,*) "namespace number and name",ns,smd_name(ns)
        write(un,*) "contents"
        do i=1,smd_nd(ns)
          write(un,*) "data name", smd_dname(ns,i)
          write(un,*) "data handle", smd_handle(ns,i)
          write(un,*) "data index", smd_index(ns,i)
          write(un,*) "data type", smd_type(ns,i)
          write(un,*) "data dims", smd_dim(ns,i)

        end do
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

      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      
      smd_nd(ns) = smd_nd(ns) + 1
      if(smd_nd(ns).gt.smd_max_nd) 
     +   call errquit(pname//'increase dataspace size to',
     +         smd_nd(ns),
     +         0)

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
     +   ma_ndim,
     +   smd_dname(ns,nd),
     +   smd_handle(ns,nd),
     +   smd_index(ns,nd)))
     +   call errquit(pname//'memory allocation'//smd_dname(ns,nd),
     +         smd_dim(ns,nd),
     +         MA_ERR)

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
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ma_type,ma_ndim

      pname = "smd_data_create"

      buffer = sname
      call smd_data_get_sindex(buffer,ns,result)
      if(.not.result)
     +   call errquit(pname//'no namespace'//buffer,
     +         0,
     +         0)

      
      smd_nd(ns) = smd_nd(ns) + 1
      if(smd_nd(ns).gt.smd_max_nd) 
     +   call errquit(pname//'increase dataspace size to',
     +         smd_nd(ns),
     +         0)

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
     +   ma_ndim,
     +   smd_dname(ns,nd),
     +   smd_handle(ns,nd),
     +   smd_index(ns,nd)))
     +   call errquit(pname//'memory allocation'//smd_dname(ns,nd),
     +         smd_dim(ns,nd),
     +         MA_ERR)


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
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,dbl_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_INT) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,int_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_STR) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     smd_string_size*ndim*ma_sizeof(dtype,1,mt_byte),root)
         else if(dtype.eq.MT_BYTE) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_LOG) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,log_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
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
     >     72*ma_sizeof(mt_int,1,mt_byte),root)
     
        call ga_sync()
     
        if(ga_nodeid().ne.root) then
          call smd_data_create_get(sname,buffer,ndim,dtype,ind)
        else
          ind = smd_index(ns,i)
        end if
        call ga_sync()

        if(dtype.eq. MT_DBL) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,dbl_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_INT) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,int_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_STR) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     smd_string_size*ndim*ma_sizeof(dtype,1,mt_byte),root)
       else if(dtype.eq.MT_BYTE) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
          call ga_brdcst(msg_smd,byte_mb(ind),
     >     ndim*ma_sizeof(dtype,1,mt_byte),root)
        else if(dtype.eq.MT_LOG) then
          write(*,*) "brdcst",ndim,dtype,buffer
          call util_flush(6)
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
c
      integer root
c
      integer ns
      character*(30) pname
      character*(72) buffer
      pname = "smd_namespace_destroy"

      do ns=1,smd_n
       call smd_data_namespace_brdcst_create(smd_name(ns),root)   
      end do
     
      return
      end

      subroutine smd_data_put_int(sname,dname,dtype,ndim,datum)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "smd_data.fh"
      character*(*) sname
      character*(*) dname
      integer ndim
      integer dtype
      integer datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i

      pname = "smd_data_create"

      call smd_data_create_get(sname,dname,ndim,dtype,ind)
      do i=1,ndim
        int_mb(ind+i-1) = datum(i)
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
      character*(smd_tag_size) datum(ndim)
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_create"

      ndim1 = smd_tag_size*ndim

      call smd_data_create_get(sname,dname,ndim1,MT_BYTE,ind)
      call smd_set_ma_char(ndim1,datum,byte_mb(ind))
      return
      end

      subroutine smd_data_get_string(sname,dname,ndim,datum)
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
      character*(smd_tag_size) datum(ndim)
c
      integer ind
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1
      integer chunk
      logical result

      pname = "smd_data_create"

      call smd_data_get_index_size(sname,dname,ind,ndim1,result)
      call smd_set_ma_char(ndim1,byte_mb(ind),datum)
      return
      end

      subroutine smd_data_put_string1(sname,dname,datum)
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
      character*(*) datum
c
      integer ind
      logical result
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1

      pname = "smd_data_create"

      ndim1 = inp_strlen(datum)

      call smd_data_create_get(sname,dname,ndim1,MT_BYTE,ind)
      call smd_set_ma_char(ndim1,datum,byte_mb(ind))
      return
      end

      subroutine smd_data_get_string1(sname,dname,datum)
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
      character*(*) datum
c
      integer ind
      integer ns,nd
      character*30 pname
      character*72 buffer
      integer i
      integer ndim1
      logical result

      pname = "smd_data_create"

      call smd_data_get_index_size(sname,dname,ind,ndim1,result)
      call smd_set_ma_char(ndim1,byte_mb(ind),datum)
      return
      end
