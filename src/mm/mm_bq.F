      block data initial_mm_bq_data
      implicit none
#include "mm_bq_data.fh"
      
      data nbq /0/   
   
      end

      subroutine mm_bq_init(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mm_bq_data.fh"
#include "mm_coords_data.fh"

      integer rtdb

      character*180 filename
      integer n, nbq0
      integer fn
      character*180 tag
      character*30 pname
      character*32 prefix

      pname = "mm_bq_init"
c     write(*,*) pname
c     prefix = "mm:coords:load"

c     tag=trim(prefix)//":name"
c     if (.not. rtdb_cget(rtdb, tag, 1, filename))
c    &     goto 911

c     if(.not.rtdb_cget(rtdb,"mm:bq",1,filename))
c    & goto 911
c     write(*,*) "filename in mm_bq_init= ", filename

c     call mm_open_file(filename,fn)
c     call util_xyz_nlines(nbq0,fn)
c     close(fn)

c     deallocate all previous allocated arrays just in case
      call mm_bq_end()

      nbq0 = nmm

c      
      if(nbq0.ne.nbq) then
        call mm_bq_end()
        if(.not.ma_alloc_get(mt_int,nbq0,'mm bq ind',
     &                       h_ibq,i_ibq))
     &   call errquit(pname//'Failed to allocate heap',nbq0,
     &         MA_ERR)

c       if(.not.ma_alloc_get(mt_int,nbq0,'mm bq chg',
c    &                       h_qbq,i_qbq))
c    &   call errquit(pname//'Failed to allocate heap',nbq0,
c    &         MA_ERR)
c     
c       if(.not.ma_alloc_get(mt_int,3*nbq0,'mm bq coord',
c    &                       h_cbq,i_cbq))
c    &   call errquit(pname//'Failed to allocate heap',nbq0,
c    &         MA_ERR)
      end if
c     
      nbq = nbq0
      e_bq = 0.0d0

      call mm_bq_load(rtdb)
      call mm_bq_create(rtdb)

c    force separate calculation of bq energy
c    ---------------------------------------
      if (.not. rtdb_put(rtdb,'dft:bq_energy',mt_dbl,1,0.0d0))
     + call errquit(pname//'setting dft:bq_energy',0,CALC_ERR)
  
      return

911   call errquit("error "//trim(pname),0,
     &        -1)

      end

      subroutine mm_bq_end()
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "mm_bq_data.fh"
#include "bq.fh"
      
      integer rtdb
        
      logical bq_end
      external bq_end

      character*30 pname
      logical status

      pname = "mm_bq_end"

      if(nbq.gt.0) then
        if (.not.ma_free_heap(h_ibq))   goto 911

        status = bq_end()

        if ((ga_nodeid().eq.0) .and. (.not. status)) then
          call errquit("failed to uninitialize bq", 0, CALC_ERR)
        end if

        nbq = 0
      end if

      return

911   call errquit("error "//trim(pname),0,
     &        -1)

      end

      subroutine mm_bq_load(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mm_bq_data.fh"
#include "mm_coords_data.fh"
      integer rtdb

      character*180 filename
      integer n
      integer fn
      character*180 tag
      character*30 pname
      character*32 prefix

      pname = "mm_bq_load"
c     write(*,*) pname    
 
      prefix = "mm:coords:load"

c     tag=trim(prefix)//":name"
c     if (.not. rtdb_cget(rtdb, tag, 1, filename))
c    $     goto 911

c     if(.not.rtdb_cget(rtdb,"mm:bq",1,filename))
c    & goto 911
c     write(*,*) "filename in mm_bq_load= ", filename

c     call mm_open_file(filename,fn)
c     call inp_init(fn,LuOut)

c     n = 0
c     do
c       if(.not.inp_read()) exit
c       if(inp_n_field().ge.1) then
c         n = n + 1
c         if(.not.inp_i(int_mb(i_ibq+n-1)))    goto 911
c         write(*,*) 'ibq = ', int_mb(i_ibq+n-1)
c       end if
c     end do

c     close(fn)
c     call inp_init(LuIn,LuOut)

c     checking bound
c     if (n.ne.nbq) goto 911

      do n=1,nbq
        int_mb(i_ibq+n-1) = int_mb(i_imm+n-1)
      end do

      return

911   call errquit("error "//trim(pname),0,
     &        -1)

      end

      subroutine mm_bq_create(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "global.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mm_coords_data.fh"
#include "bq_params.fh"
#include "bq.fh"

      integer rtdb


      logical bq_init
      external bq_init
      logical status
      integer handle
      integer ncenter
      integer i, j
      integer h_cbq, i_cbq
      integer h_qbq, i_qbq
      double precision scale
      character*255 aunit

      character*255 namespace
      character*180 pname 
     
      pname = "mm_bq_create"
c     write(*,*) pname

      namespace = bq_default
c     write(*,*) 'namespace = ', namespace

      aunit = "ang"

      ncenter = nmm

      call util_length_scale(aunit,scale)

      if ( .not. ma_push_get( MT_DBL, 3*ncenter, 
     &      'cbq scratch',
     &      h_cbq, i_cbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      ncenter, MA_ERR)
c
      call dfill(3*ncenter,0.0d0,dbl_mb(i_cbq),1)
c
      if ( .not. ma_push_get( MT_DBL, ncenter, 
     &      'qbq scratch',
     &      h_qbq, i_qbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      ncenter, MA_ERR)
c
      call dfill(ncenter,0.0d0,dbl_mb(i_qbq),1)
  
      do i=1,ncenter
         dbl_mb(i_cbq+3*(i-1))   = dbl_mb(i_rmm+3*(i-1))*scale
         dbl_mb(i_cbq+3*(i-1)+1) = dbl_mb(i_rmm+3*(i-1)+1)*scale
         dbl_mb(i_cbq+3*(i-1)+2) = dbl_mb(i_rmm+3*(i-1)+2)*scale
         dbl_mb(i_qbq+i-1)       = dbl_mb(i_chgmm+i-1)
         
c        write(*,*) 'bq coords= ', (dbl_mb(i_cbq+3*(i-1)+j),j=0,2),
c    >                              dbl_mb(i_qbq+i-1)
      end do

      call bq_rtdb_delete(rtdb,namespace)
 
c     write(*,*) 'bq_create'
      if(.not.bq_create(namespace,handle))   
     &       call errquit('bq_create failed',0,0)
c     write(*,*) 'handle= ', handle

      if(.not.bq_set(handle,ncenter,dbl_mb(i_qbq),dbl_mb(i_cbq)))
     &       call errquit('bq_set failed',0,0)
      
c     call bq_print_info(handle)

      if(.not.bq_rtdb_store(rtdb,handle))
     &       call errquit('bq_store failed',0,0)

      if(.not.bq_destroy(handle))
     &       call errquit('bq_destroy failed',0,0)
      
c     write(*,*) 'trying put bq into rtdb'
      if(ncenter.ne.0 .and. namespace.eq.bq_default) then
c       write(*,*) 'rtdb_cput'
        if(.not. rtdb_cput(rtdb,"bq" , 1,bq_default)) then
          call errquit(pname//'failed to set default namespace',
     >                 0,RTDB_ERR)
        end if
c       write(*,*) 'calling bq_init'
        status = bq_init(rtdb)
c       write(*,*) 'status in bq create= ', status
        if ((ga_nodeid().eq.0) .and. (.not. status)) then
          call errquit("failed to initialize bq", 0, CALC_ERR)
        end if
      end if

      if (.not. ma_pop_stack(h_qbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      ncenter, MA_ERR)

      if (.not. ma_pop_stack(h_cbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      ncenter, MA_ERR)
          
c     write(*,*) 'out', pname  

      end
