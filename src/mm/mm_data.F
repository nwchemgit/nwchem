c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_data
      implicit none
#include "mm_data.fh"
c     
      data nvdw /0/
      data nmm /0/
      data nqm /0/

      end

      function  mm_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
      logical mm_init
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
c
      pname = "mm_init"
c
      if (nvdw .eq. 0) then
          call mm_data_load_vdw(rtdb)
      else
          goto 911
      end if
      if (nmm .eq. 0) then
          call mm_data_load_coords(rtdb)
      else
          goto 911
      end if
      call mm_load_qm_coords(rtdb)
      mm_init = .true.
      return
911   mm_init = .false.
      return
      end

      function  mm_end(rtdb)
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
      if (nmm .gt. 0) then
        if (.not.ma_free_heap(h_rmm)) goto 911
        nmm = 0
      end if
      if (nvdw .gt. 0) then
        if (.not.ma_free_heap(h_bvdw)) goto 911
        if (.not.ma_free_heap(h_avdw)) goto 911
        if (.not.ma_free_heap(h_jvdw)) goto 911
        if (.not.ma_free_heap(h_ivdw)) goto 911
        nvdw = 0
      end if
      if(nqm.gt.0) then
          if (.not.ma_free_heap(h_tqm)) goto 911
          if (.not.ma_free_heap(h_rqm)) goto 911
          nqm = 0
      end if
      mm_end = .true.
      return
911   mm_end = .false.
      return
      end

      subroutine mm_test()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"

      write(*,*) "handles",nvdw,h_ivdw,h_jvdw,h_avdw,h_bvdw
      end

      subroutine mm_data_load_coords(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
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
      double precision scale
c
      pname = "mm_input_load_coords"
      prefix = "mm:coords:load"
c
      tag=trim(prefix)//":name"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      tag=trim(prefix)//":format"
      if(.not. rtdb_get(rtdb,tag,mt_int,3,iloc))   
     $   iloc = (/ 1, 2, 3 /) 

      call mm_open_file(filename,fn)
c     ----------------------
c     count numbers of lines
c     ----------------------
      call util_xyz_nlines(nmm,fn)
      write(*,*) "nmm=",nmm

      if ( .not. ma_alloc_get( MT_DBL, 3*nmm, 
     &      'mm coords',
     &      h_rmm, i_rmm) ) call errquit(
     &      'failed memory allocation '//pname,
     &      nmm, MA_ERR)
      call util_xyz_read_coords(nmm,n,iloc,dbl_mb(i_rmm),fn)
      close(fn)
c      
      tag=trim(prefix)//":units"
      if (.not.rtdb_cget(rtdb, tag, 1,aunit)) aunit="ang"
      call util_convert_units(aunit,"angstrom",scale)   
      if (scale.ne.1.0) then
        call dscal(3*nmm, scale,dbl_mb(i_rmm),1)
      end if
      do i=1,nmm
         write(54,*) dbl_mb(i_rmm+i-1)
      enddo
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_data_load_vdw(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_data.fh"
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
c
      pname = "mm_input_load_vdw"
c
      tag="mm:vdw:load"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      call mm_open_file(filename,fn)
c     ------------------------- 
c     count numbers of vdw pairs
c     ------------------------- 
      call inp_init(fn,LuOut)
      n = 0
      do
        if (.not.inp_read()) exit
        n=n+1
      end do
      nvdw=n
      rewind(fn)
      call inp_init(fn,LuOut)
c
      max_vdw = 2500
      call util_xyz_nlines(max_vdw,fn)
      write(*,*) "nvdw,max_vdw",nvdw,max_vdw
c
      if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &      'vdw scratch i',
     &      h_ivdw, i_ivdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &      'vdw scratch j',
     &      h_jvdw, i_jvdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &      'vdw scratch A',
     &      h_avdw, i_avdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &      'vdw scratch B',
     &      h_bvdw, i_bvdw) ) call errquit(
     &      'failed scratch space '//pname,
     &      nvdw, MA_ERR)

      call inp_init(fn,LuOut)
      n = 0
10    continue
      if (.not.inp_read()) goto 30
      message = "checking for number of fields"
      if(inp_n_field().ne.4) goto 911
      message = "checking for bounds"
      if(n.gt.nvdw) goto 911
      message = "reading vdw params"
      n = n+1
      write(*,*) "n=",n
      if (.not.inp_i(int_mb(i_ivdw+n-1))) goto 911
      if (.not.inp_i(int_mb(i_jvdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_avdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_bvdw+n-1))) goto 911
      goto 10
30    continue
      close(fn)
      call inp_init(LuIn,LuOut)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_open_file(filename,fn)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       fn
c
      character*180 buffer
      character*180 message
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
c
      pname = "mm_open_file"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
c     first try to open file in the run directory
      buffer = filename
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=10)
      goto 800
10    continue
c     now try perm directory
      call util_file_name_resolve(buffer, .false.)
      message = "opening file "//buffer
      open(unit=fn,file=buffer,status='old',form="formatted",ERR=911)
800   continue
      write(luout,*) "Successfully "//trim(message)
      write(luout,*)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mm_data.fh"
#include "mafdecls.fh"

      integer rtdb
      double precision e
      double precision scale_energy,scale_length
      double precision scale_force
      double precision evdw
      character*30 message
      integer i,j

      call mm_load_qm_coords(rtdb)
      call mm_add_vdw_energy(nqm,dbl_mb(i_rqm),
     &                      nmm,dbl_mb(i_rmm),
     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
     &                      evdw)

      call util_convert_units("kcal","au",scale_energy)
 
      evdw = evdw*scale_energy
C      e=0
      e = e + evdw
      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mm_data.fh"
#include "mafdecls.fh"

      integer rtdb
      double precision e
      integer n
      integer i_g,h_g
      double precision g(3,n)
      double precision scale_energy,scale_length
      double precision scale_force
      double precision evdw
      character*30 message
      character*16 thetag
      double precision crd(3),q
      integer i,j

      call mm_load_qm_coords(rtdb)
      if ( .not. ma_push_get( MT_DBL, 3*nqm, 
     &      'vdw scratch B',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*nqm,0.0d0,dbl_mb(i_g),1)

      call mm_add_vdw_egrad(nqm,dbl_mb(i_rqm),
     &                      nmm,dbl_mb(i_rmm),
     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
C     &                      n,e,g)
     &                      n,evdw,dbl_mb(i_g))

      call util_convert_units("kcal","au",scale_energy)
      call util_convert_units("ang","au",scale_length)
 
      evdw = evdw*scale_energy
      e = e + evdw
      scale_force = scale_energy/scale_length
C      call dscal(3*nqm, scale_force,dbl_mb(i_g),1)
C      call dfill(3*nqm,0.0d0,g,1)
      call daxpy(3*nqm,scale_force,dbl_mb(i_g),1,g,1)

      write(6,1000) "MM/VDW",
     $     'x','y','z','x','y','z'
      do i = 1, nqm
         write(6,2000) i, byte_mb(i_tqm+16*(i-1)),
     $        (dbl_mb(i_rqm+3*(i-1)+j),j=0,2),
     $        (scale_force*dbl_mb(i_g +3*(i-1)+j),j=0,2)
      enddo
      write(6,*)


      write(6,1000) "QM + MM/VDW",
     $     'x','y','z','x','y','z'
      do i = 1, nqm
         write(6,2000) i, byte_mb(i_tqm+16*(i-1)),
     $        (dbl_mb(i_rqm+3*(i-1)+j),j=0,2),
c     $        (g(j,i),j=1,3)
     $        g(1,i),g(2,i),g(3,i)
      enddo
      write(6,*)



 1000 format(/,/,25X,A,' ENERGY GRADIENTS',/,/,4X,'atom',15X,
     $     'coordinates',
     $     24X,'gradient',/,6X,2(1X,(3(10X,A1))))
 2000 format(1X,I3,1X,A4,2(1X,3(1X,F10.6)))
      write(6,*)
      call util_flush(6)

      message = "memory deallocation"
      if(.not.ma_pop_stack(h_g)) goto 911
      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_add_vdw_egrad(nqm,rqm,
     &                            nmm,rmm,
     &                            np,p1,p2,
     &                            a,b,
     &                            ng,evdw,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer nqm
      double precision rqm(3,nqm)
      integer nmm
      double precision rmm(3,nmm)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      integer ng
      double precision e
      double precision g(3,ng)
c
      integer n
      integer i1,i2
      double precision evdw
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 message

      evdw = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)-nqm
        message = "checking bounds"
        if(i1.gt.nqm .or. i2.gt.nmm) then
            write(*,*) i1,i2
            goto 911
        end if
        dx = rmm(1,i2) - rqm(1,i1)
        dy = rmm(2,i2) - rqm(2,i1)
        dz = rmm(3,i2) - rqm(3,i1)
      
        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw = evdw + f12 - f6
        
        df = (12.0d0*f12 - 6.0d0*f6)*r2inv
        g(1,i1) = g(1,i1) + df*dx
        g(2,i1) = g(2,i1) + df*dy
        g(3,i1) = g(3,i1) + df*dz
   
C        write(*,*) "qm coord=",rqm(:,i1),i1
C        write(*,*) "mm coord=",rmm(:,i2),i2
C        write(*,*) "force =",g(:,i1),i1
      end do
      write(*,*) "vdw energy",evdw
      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_add_vdw_energy(nqm,rqm,
     &                            nmm,rmm,
     &                            np,p1,p2,
     &                            a,b,
     &                            evdw)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer nqm
      double precision rqm(3,nqm)
      integer nmm
      double precision rmm(3,nmm)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      double precision e
c
      integer n
      integer i1,i2
      double precision evdw
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 message

      evdw = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)-nqm
        message = "checking bounds"
        if(i1.gt.nqm .or. i2.gt.nmm) then
            write(*,*) i1,i2
            goto 911
        end if
        dx = rmm(1,i2) - rqm(1,i1)
        dy = rmm(2,i2) - rqm(2,i1)
        dz = rmm(3,i2) - rqm(3,i1)
      
        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw = evdw + f12 - f6
        
      end do
      write(*,*) "vdw energy",evdw
      return
911   call errquit("error "//trim(message),0,-1)
      end
 
      subroutine mm_load_qm_coords(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "mm_data.fh"

      integer rtdb
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
      if(.not.status) 
     & call errquit('cons_init: geom_create?',70, GEOM_ERR)
      if (nat.ne.nqm) then
        if(nqm.gt.0) then
          if (.not.ma_free_heap(h_tqm)) goto 911
          if (.not.ma_free_heap(h_rqm)) goto 911
        end if
        if ( .not. ma_alloc_get( MT_DBL, nat*3, 'qm coord',
     &        h_rqm, i_rqm) ) goto 911
        if(.not.ma_alloc_get(mt_byte,16*nat,'qm tags',h_tqm,i_tqm))
     +       goto 911

      end if
      nqm = nat
      if(.not.geom_cart_get1(geom, nqm,byte_mb(i_tqm),dbl_mb(i_rqm)))
     &    goto 911
 
      if(.not.geom_destroy(geom)) 
     &    goto 911

      call util_convert_units("au","angstrom",scale)
      call dscal(3*nqm, scale,dbl_mb(i_rqm),1)
      do i=1,3*nqm
        write(99,*) dbl_mb(i_rqm+i-1)
      end do
      return
911   call errquit("error "//trim(message),0,-1)
      end
 
 
