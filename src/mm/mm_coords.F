c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_coords_data
      implicit none
#include "mm_coords_data.fh"
c     
      data ntot /0/
      data nmm /0/
      data nqm /0/

      end

      subroutine  mm_coords_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
c
      character*180 filename
      integer       nmax,n
c
      character*180 message
      character*30 pname
      integer nqm,nmm
c
      pname = "mm_coords_init"
c
      call mm_coords_nqm(nqm,rtdb)
      call mm_coords_nmm(nmm,rtdb)
c
      call mm_coords_allocate(nqm,nmm)
c     
      call mm_coords_mm_load(rtdb)
      call mm_coords_qm_load(rtdb)
c
      call mm_coords_test()
      return
911   call errquit("error "//trim(message),0,-1)
      return
      end

      subroutine  mm_coords_allocate(nqm0,nmm0)
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer nqm0,nmm0
c
      character*180 message
      character*30 pname
c
      integer ntot0
      logical ignore
c
      pname = "mm_coords_allocate"
      ntot0 = nqm0 + nmm0
      if(ntot0.ne.ntot) then
        call mm_coords_end()
        if ( .not. ma_allocate_heap(mt_dbl,3*ntot0, 
     +                          "mm coords",
     +                           h_r)) goto 911
        if(.not.ma_allocate_heap(mt_byte,16*ntot0,'mm tags',h_t))
     +        goto 911
      end if
c
      if(.not.ma_get_index(h_r,i_r)) goto 911
      if(.not.ma_get_index(h_t,i_t)) goto 911
      ntot = ntot0
      nqm  = nqm0
      nmm  = nmm0
c     
      call dfill(3*ntot,0.0d0,dbl_mb(i_r),1)
      call mm_coords_tags_reset(ntot,byte_mb(i_t))
c      
      return
911   call errquit("error "//trim(message),0,-1)
      return
      end

      subroutine mm_coords_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer rtdb
c
      character*30 pname
      pname = "mm_coords_end"
c
      if(ntot.gt.0) then
          if (.not.ma_free_heap(h_t)) goto 911
          if (.not.ma_free_heap(h_r)) goto 911
          ntot = 0
          nqm  = 0
          nmm = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_coords_nmm(nmm,rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer nmm
      integer rtdb
c
      character*180 filename
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
      integer iloc(3)
      integer i
      character*32 prefix
      character*32 aunit
c
      pname = "mm_coords_nmm"
      prefix = "mm:coords:load"
c
      tag=trim(prefix)//":name"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      call mm_open_file(filename,fn)
c     ----------------------
c     count numbers of lines
c     ----------------------
      call util_xyz_nlines(nmm,fn)
      close(fn)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_coords_nqm(nqm,rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"

      integer nqm,rtdb
      external geom_rtdb_ncent
      logical geom_rtdb_ncent
      character*30 message
c

      message = "calling geom_rtdb_ncent"
      if (.not. geom_rtdb_ncent(rtdb, 'geometry', nqm)) goto 911

      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_coords_mm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      integer rtdb
c
      character*180 filename
c
      integer n
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
      integer iloc(3)
      integer i
      character*32 prefix
      character*32 aunit
      double precision scale
      integer i_rmm,i_tmm
      character*16 mmtag
c      
      i_rmm = i_r+3*nqm
      i_tmm = i_t+16*nqm
c
      pname = "mm_coords_nmm"
      prefix = "mm:coords:load"
c
      tag=trim(prefix)//":name"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      tag=trim(prefix)//":format"
      if(.not. rtdb_get(rtdb,tag,mt_int,3,iloc))   
     $   iloc = (/ 1, 2, 3 /) 
c
      call mm_open_file(filename,fn)
c      
      call util_xyz_read_coords(nmm,n,iloc,dbl_mb(i_rmm),fn)
      close(fn)
c      
      tag=trim(prefix)//":units"
      if (.not.rtdb_cget(rtdb, tag, 1,aunit)) aunit="ang"
      call util_convert_units(aunit,"angstrom",scale)   
      if (scale.ne.1.0) then
        call dscal(3*nmm, scale,dbl_mb(i_rmm),1)
      end if
      mmtag = "MM"
      call mm_coords_tags_set(nmm,"MM",byte_mb(i_tmm))
c      
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_coords_qm_load(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "mm_coords_data.fh"

      integer rtdb
c
      integer geom
      integer nat
      logical status
      character*30 message
      double precision scale
      logical geom_cart_get1
      external geom_cart_get1
c
c     load geometry 
c     -------------
      if (.not. geom_create(geom, 'geometry'))
     &     call errquit('cons_load_geom',0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     &     call errquit('cons_load_geom',0, RTDB_ERR)
c
c     get cart coordinates
c     --------------------
      status=geom_ncent(geom,nat)
      message = "nat not eauals nmax"
      if(nat.ne.nqm) then
          write(*,*) nat,nqm
          goto 911
      end if
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if(.not.geom_cart_get1(geom,nqm,byte_mb(i_t),dbl_mb(i_r)))
     &    goto 911

      if(.not.geom_destroy(geom)) 
     &    goto 911

      call util_convert_units("au","angstrom",scale)
      call dscal(3*nat, scale,dbl_mb(i_r),1)
      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_coords_test()
      implicit none
#include "mm_coords_data.fh"
#include "mafdecls.fh"

      integer i,j

      do i=1,ntot
         write(6,'(16A1,1X,3(1X,F10.6))') 
     $        (byte_mb(i_t+16*(i-1)+j),j=0,15),
     $        (dbl_mb(i_r+3*(i-1)+j),j=0,2)
      end do
      
      end

      subroutine mm_coords_tags_set(n,tag,t)
      implicit none

      integer n
      character*(*) tag
      character*16 t(n)
      integer i

      do i=1,n
         t(i) = trim(tag)
      end do
      
      end

      subroutine mm_coords_tags_reset(n,t)
      implicit none

      integer n
      character*16 t(n)
      integer i

      do i=1,n
         t(i) = "X"
      end do
      
      end
