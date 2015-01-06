c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_coords_data
      implicit none
#include "mm_data.fh"
c     
      data ntot /0/
      data nmm /0/
      data nqm /0/

      end

      function  mm_coords_end(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      logical mm_end
      integer rtdb
c
      character*30 pname
      pname = "mm_end"
c
      if(ntot.gt.0) then
          if (.not.ma_free_heap(h_t)) goto 911
          if (.not.ma_free_heap(h_r)) goto 911
          ntot = 0
          nqm  = 0
          nmm = 0
      end if
      mm_end = .true.
      return
911   mm_end = .false.
      return
      end
      subroutine  mm_coords_allocate(nqm0,nmm0)
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
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
      double precision scale
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
      integer nqm1,nmm1
c
      pname = "mm_init"
      ntot0 = nqm0 + nmm0
      if(ntot0.ne.ntot) then
        if ( .not. ma_alloc_get(MT_DBL,3*ntot, 
       $                        "mm coords",
       $                         h_r,i_r)) goto 911
         if(.not.ma_alloc_get(mt_byte,16*ntot,'mm tags',h_t,i_t))
       +      goto 911

      end if
c
      call mm_coords_nqm(nqm,rtdb)
      call mm_coords_nmm(nmm,rtdb)
      ntot = nqm + nmm
c     
      message = "allocating mm coords"
      return
911   call errquit("error "//trim(message),0,-1)
      return
      end

      subroutine  mm_coords_init(rtdb)
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
      integer       nmax,n
c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*180 tag
      character*30 pname,atag
      double precision scale
c
      logical util_io_unit
      external util_io_unit
      integer max_vdw
      integer nqm1,nmm1
c
      pname = "mm_init"
c
      call mm_coords_nqm(nqm,rtdb)
      call mm_coords_nmm(nmm,rtdb)
      ntot = nqm + nmm
c     
      message = "allocating mm coords"
      if ( .not. ma_alloc_get(MT_DBL,3*ntot, 
     $                        "mm coords",
     $                         h_r,i_r)) goto 911
       if(.not.ma_alloc_get(mt_byte,16*ntot,'mm tags',h_t,i_t))
     +      goto 911
      call mm_coords_mm_load(rtdb,nmm,byte_mb(i_t+16*nqm),
     $                        dbl_mb(i_r+3*nqm))
      call mm_coords_qm_load(rtdb,nqm,byte_mb(i_t),dbl_mb(i_r))
      call mm_test(ntot,byte_mb(i_t),dbl_mb(i_r))
      return
911   call errquit("error "//trim(message),0,-1)
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

      subroutine mm_coords_mm_load(rtdb,nmax,t,c)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer nmax
      character*(16) t(nmax)
      double precision c(3,nmax)
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
      call util_xyz_read_coords(nmax,n,iloc,c,fn)
      close(fn)
c      
      tag=trim(prefix)//":units"
      if (.not.rtdb_cget(rtdb, tag, 1,aunit)) aunit="ang"
      call util_convert_units(aunit,"angstrom",scale)   
      if (scale.ne.1.0) then
        call dscal(3*n, scale,c,1)
      end if
      do i=1,nmax
        t(i) = "MM"
      enddo
c      
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

      subroutine mm_coords_qm_load(rtdb,nmax,t,c)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "mm_data.fh"

      integer rtdb
      integer nmax
      double precision c(3,nmax)
      character*(16) t(nmax)
c
      integer i_c,h_c
      integer geom
      integer nat
      logical status
      character*30 message
      integer i
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
      if(nat.ne.nmax) then
          write(*,*) nat,nmax
          goto 911
      end if
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if(.not.geom_cart_get1(geom,nat,t,c))
     &    goto 911
 
      if(.not.geom_destroy(geom)) 
     &    goto 911

      call util_convert_units("au","angstrom",scale)
      call dscal(3*nat, scale,c,1)
      return
911   call errquit("error "//trim(message),0,-1)
      end

