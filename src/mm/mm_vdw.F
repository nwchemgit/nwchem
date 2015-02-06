c $Id: bq_data.F 23019 2012-10-30 00:59:12Z d3y133 $
      block data initial_mm_vdw_data
      implicit none
#include "mm_vdw_data.fh"
c     
      data nvdw /0/

      end

      subroutine mm_vdw_end()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"
c
      character*30 pname
      pname = "mm_vdw_end"
c
      if (nvdw .gt. 0) then
        if (.not.ma_free_heap(h_bvdw)) goto 911
        if (.not.ma_free_heap(h_avdw)) goto 911
        if (.not.ma_free_heap(h_jvdw)) goto 911
        if (.not.ma_free_heap(h_ivdw)) goto 911
        nvdw = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_vdw_allocate(n)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"
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
      pname = "mm_vdw_allocate"

      if(n.ne.nvdw) then 
        call mm_vdw_end()
        nvdw = n
        write(*,*) "allocating vdw data"
c
        if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &        'vdw scratch i',
     &        h_ivdw, i_ivdw) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw, MA_ERR)

        if ( .not. ma_alloc_get( MT_INT, nvdw, 
     &        'vdw scratch j',
     &        h_jvdw, i_jvdw) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw, MA_ERR)

        if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &        'vdw scratch A',
     &        h_avdw, i_avdw) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw, MA_ERR)

        if ( .not. ma_alloc_get( MT_DBL, nvdw, 
     &        'vdw scratch B',
     &        h_bvdw, i_bvdw) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw, MA_ERR)

      end if
      call ifill(n,0.0,int_mb(i_ivdw),1)
      call ifill(n,0.0,int_mb(i_jvdw),1)
      call dfill(n,0.0d0,dbl_mb(i_avdw),1)
      call dfill(n,0.0d0,dbl_mb(i_bvdw),1)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_vdw_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"
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
      pname = "mm_vdw_load"
c
      tag="mm:vdw:load"
      if (.not. rtdb_cget(rtdb, tag, 1, filename))
     $     goto 911
c
      call mm_open_file(filename,fn)
c
      call util_xyz_nlines(n,fn)
      call mm_vdw_allocate(n)
      rewind(fn)
c      
c     loading parameters here
c
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
      if (.not.inp_i(int_mb(i_ivdw+n-1))) goto 911
      if (.not.inp_i(int_mb(i_jvdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_avdw+n-1))) goto 911
      if (.not.inp_f(dbl_mb(i_bvdw+n-1))) goto 911
      goto 10
30    continue
      close(fn)
      call inp_init(LuIn,LuOut)
      call mm_vdw_test()
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_vdw_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mm_vdw_data.fh"
#include "mm_coords_data.fh"
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

      call mm_coords_qm_load(rtdb,nqm,byte_mb(i_t),dbl_mb(i_r))
      if ( .not. ma_push_get( MT_DBL, 3*nqm, 
     &      'vdw scratch B',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*nqm,0.0d0,dbl_mb(i_g),1)

      call mm_vdw_egrad(ntot,dbl_mb(i_r),
     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
     &                      n,evdw,dbl_mb(i_g))

C      call mm_vdw_egrad(nqm,dbl_mb(i_r),
C     &                      nmm,dbl_mb(i_r+3*nqm),
C     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
C     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
C     &                      n,evdw,dbl_mb(i_g))
C
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
         write(6,2000) i, byte_mb(i_t+16*(i-1)),
     $        (dbl_mb(i_r+3*(i-1)+j),j=0,2),
     $        (scale_force*dbl_mb(i_g +3*(i-1)+j),j=0,2)
      enddo
      write(6,*)


      write(6,1000) "QM + MM/VDW",
     $     'x','y','z','x','y','z'
      do i = 1, nqm
         write(6,2000) i, byte_mb(i_t+16*(i-1)),
     $        (dbl_mb(i_r+3*(i-1)+j),j=0,2),
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

      subroutine mm_vdw_add_energy(rtdb,e)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"
#include "mm_coords_data.fh"

      integer rtdb
      double precision e
      double precision scale_energy,scale_length
      double precision scale_force
      double precision evdw
      character*30 message
      integer i,j

      call mm_coords_qm_load(rtdb,nqm,byte_mb(i_t),dbl_mb(i_r))
      call mm_vdw_energy(ntot,dbl_mb(i_r),
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

      subroutine mm_vdw_egrad(nt,r,
     &                            np,p1,p2,
     &                            a,b,
     &                            ng,evdw,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer nt
      double precision r(3,nt)
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
        i2 = p2(n)
        message = "checking bounds"
        if(i1.gt.nt .or. i2.gt.nt) then
            write(*,*) i1,i2
            goto 911
        end if
        dx = r(1,i2) - r(1,i1)
        dy = r(2,i2) - r(2,i1)
        dz = r(3,i2) - r(3,i1)
      
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
   
        write(*,*) "qm coord=",r(:,i1),i1
        write(*,*) "mm coord=",r(:,i2),i2
C        write(*,*) "force =",g(:,i1),i1
      end do
      write(*,*) "vdw energy",evdw
      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_vdw_energy(nt,r,
     &                            np,p1,p2,
     &                            a,b,
     &                            evdw)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer nt
      double precision r(3,nt)
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
      character*30 pname

      pname = "mm_vdw_energy1"
      evdw = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)
        message = "checking bounds"
        if(i1.gt.nt .or. i2.gt.nt) then
            write(*,*) i1,i2
            goto 911
        end if
        dx = r(1,i2) - r(1,i1)
        dy = r(2,i2) - r(2,i1)
        dz = r(3,i2) - r(3,i1)
      
        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw = evdw + f12 - f6
        
        write(*,*) "r=",r(:,i1)
      end do
      write(*,*) "vdw energy",evdw
      return
911   call errquit("error "//trim(pname)//trim(message),0,-1)
      end

      subroutine mm_vdw_test()
      implicit none
#include "mm_vdw_data.fh"
#include "mafdecls.fh"

      integer i,j

      do i=1,nvdw
         write(6,'(2(1X,I3),1X,2(1X,E12.6))') 
     $        (int_mb(i_ivdw+i-1)),
     $        (int_mb(i_jvdw+i-1)),
     $        (dbl_mb(i_avdw+i-1)),
     $        (dbl_mb(i_bvdw+i-1))
      end do
      
      end
