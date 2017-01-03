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

      subroutine mm_vdw14_end()
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
      pname = "mm_vdw14_end"
c
      if (nvdw14 .gt. 0) then
        if (.not.ma_free_heap(h_bvdw14)) goto 911
        if (.not.ma_free_heap(h_avdw14)) goto 911
        if (.not.ma_free_heap(h_jvdw14)) goto 911
        if (.not.ma_free_heap(h_ivdw14)) goto 911
        if (.not.ma_free_heap(h_scnb)) goto 911
c       if (.not.ma_free_heap(h_scee)) goto 911
        nvdw14 = 0
      end if
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_vdw_deallocate()
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"
      call mm_vdw_end()    
      call mm_vdw14_end()    

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
c     write(*,*) pname

      if(n.ne.nvdw) then 
        call mm_vdw_end()
        nvdw = n
c       write(*,*) "allocating vdw data"
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

      subroutine mm_vdw14_allocate(n)
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
      pname = "mm_vdw14_allocate"

      if(n.ne.nvdw14) then 
        call mm_vdw14_end()
        nvdw14 = n
c       write(*,*) "allocating vdw14 data"
c
        if ( .not. ma_alloc_get( MT_INT, nvdw14, 
     &        'vdw scratch i',
     &        h_ivdw14, i_ivdw14) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw14, MA_ERR)

        if ( .not. ma_alloc_get( MT_INT, nvdw14, 
     &        'vdw scratch j',
     &        h_jvdw14, i_jvdw14) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw14, MA_ERR)

        if ( .not. ma_alloc_get( MT_DBL, nvdw14, 
     &        'vdw scratch A',
     &        h_avdw14, i_avdw14) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw14, MA_ERR)

        if ( .not. ma_alloc_get( MT_DBL, nvdw14, 
     &        'vdw scratch B',
     &        h_bvdw14, i_bvdw14) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw14, MA_ERR)

        if ( .not. ma_alloc_get( MT_DBL, nvdw14, 
     &        'vdw scratch scnb',
     &        h_scnb, i_scnb) ) call errquit(
     &        'failed scratch space '//pname,
     &        nvdw14, MA_ERR)

c       if ( .not. ma_alloc_get( MT_DBL, nvdw14, 
c    &        'vdw scratch scee',
c    &        h_scee, i_scee) ) call errquit(
c    &        'failed scratch space '//pname,
c    &        nvdw14, MA_ERR)
      end if
      call ifill(n,0.0,int_mb(i_ivdw14),1)
      call ifill(n,0.0,int_mb(i_jvdw14),1)
      call dfill(n,0.0d0,dbl_mb(i_avdw14),1)
      call dfill(n,0.0d0,dbl_mb(i_bvdw14),1)
      call dfill(n,0.0d0,dbl_mb(i_scnb),1)
c     call dfill(n,0.0d0,dbl_mb(i_scee),1)

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
      integer       nmax,n1,n2
      integer       nvdw0, n14
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
c     write(*,*) pname

      nvdw0 = 0
      n14   = 0
      call mm_get_nvdw(rtdb,nvdw0)
      call mm_get_nvdw14(rtdb,n14)

      call mm_vdw_allocate(nvdw0)
      call mm_vdw14_allocate(n14)

      call mm_vdw_parm_load(rtdb)
      call mm_vdw14_parm_load(rtdb)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_get_nvdw(rtdb,nvdw0)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb 
      integer nvdw0

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
        
      pname = "mm_get_nvdw"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      nvdw0 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'vdw',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.4) then
                nvdw0 = nvdw0 + 1
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_get_nvdw14(rtdb,n14)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb 
      integer n14

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
        
      pname = "mm_get_nvdw14"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      n14 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'scaled_vdw',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
c             else if(inp_n_field().eq.6) then
                n14 = n14 + 1
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_vdw_parm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"

      integer rtdb 
      integer n1

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
        
      pname = "mm_vdw_parm_load"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      n1 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'vdw',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.4) then
                message = "checking for bounds"
                if(n1.gt.nvdw) goto 911
                message = "reading vdw params"
                n1 = n1 + 1
                if (.not.inp_i(int_mb(i_ivdw+n1-1))) goto 911
                if (.not.inp_i(int_mb(i_jvdw+n1-1))) goto 911
                if (.not.inp_f(dbl_mb(i_avdw+n1-1))) goto 911
                if (.not.inp_f(dbl_mb(i_bvdw+n1-1))) goto 911
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine mm_vdw14_parm_load(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "mm_vdw_data.fh"

      integer rtdb 
      integer n1

      character*180 filename
      character*255 field
      character*180 message
      character*30  pname
      integer       fn
        
      pname = "mm_vdw14_parm_load"   
   
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911

      call mm_open_file(filename,fn)
      
      call inp_init(fn,LuOut)
      n1 = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'scaled_vdw',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
c             else if(inp_n_field().eq.6) then
                message = "checking for bounds"
                if(n1.gt.nvdw14) goto 911
                message = "reading vdw14 params"
                n1 = n1 + 1
                if (.not.inp_i(int_mb(i_ivdw14+n1-1))) goto 911
                if (.not.inp_i(int_mb(i_jvdw14+n1-1))) goto 911
                if (.not.inp_f(dbl_mb(i_avdw14+n1-1))) goto 911
                if (.not.inp_f(dbl_mb(i_bvdw14+n1-1))) goto 911
                if (.not.inp_f(dbl_mb(i_scnb+n1-1))) goto 911
c               if (.not.inp_f(dbl_mb(i_scee+n1-1))) goto 911
                goto 10
              end if
            end if
            goto 10
          end if
        end if
        goto 5
      end if
20    continue
      close(fn)
      call inp_init(LuIn,LuOut)

      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end
c
c     adjust scnb between QM-MMlnk to be 1.0d0
c     
      subroutine mm_vdw_adj_scnb()
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_data.fh"
#include "mm_link_data.fh"
c     local variables
      integer i, j, k
      integer i_indx, j_indx
      integer i1, i2

      character*30 pname
      pname = "mm_vdw_adj_scnb"
c     write(*,*) pname     

c     TP: 1-4 VDW between QM and MMlnk is not scaled
      do i = 1, nvdw14
        i1 = int_mb(i_ivdw14+i-1)
        i2 = int_mb(i_jvdw14+i-1)
        do j = 1, nqm 
          i_indx = int_mb(i_iqm+j-1)
          do k = 1,nlink
            j_indx = int_mb(i_lb+nlink+k-1)
            if(i1.eq.i_indx.and.i2.eq.j_indx) then 
              dbl_mb(i_scnb+i-1) = 1.0d0
            else if(i1.eq.j_indx.and.i2.eq.i_indx) then
              dbl_mb(i_scnb+i-1) = 1.0d0
            end if       
          end do
        end do
      end do
      end
c
      subroutine mm_vdw_add_egrad(rtdb,e,n,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_data.fh"
#include "mm_vdw_coords_data.fh"

      integer rtdb
      double precision e
      integer n
      integer i_g,h_g
      double precision g(3,n)
      double precision scale_energy,scale_length
      double precision scale_force
      double precision evdw,evdw14
      character*30 message
      character*16 thetag
      double precision crd(3),q
      integer i,j
      character*30 pname

      pname = "mm_vdw_add_egrad"
c     write(*,*) pname

      if ( .not. ma_push_get( MT_DBL, 3*n, 
     &      'vdw scratch B',
     &      h_g, i_g) ) call errquit(
     &      'failed scratch space ',
     &      nqm, MA_ERR)

      call dfill(3*n,0.0d0,dbl_mb(i_g),1)

      call mm_vdw_egrad(nvtot,dbl_mb(i_rvdw),
     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
     &                      n,evdw,dbl_mb(i_g))

      call mm_vdw14_egrad(nv14tot,dbl_mb(i_rvdw14),
     &                      nvdw14,int_mb(i_ivdw14),int_mb(i_jvdw14),
     &                      dbl_mb(i_avdw14),dbl_mb(i_bvdw14),
     &                      dbl_mb(i_scnb),n,evdw14,dbl_mb(i_g))
C
      call util_convert_units("kcal","au",scale_energy)
      call util_convert_units("ang","au",scale_length)

      evdw = evdw*scale_energy
      evdw14 = evdw14*scale_energy
      e = e + evdw + evdw14

      scale_force = scale_energy/scale_length
      call daxpy(3*n,scale_force,dbl_mb(i_g),1,g,1)


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
#include "mm_vdw_coords_data.fh"
#include "mm_coords_data.fh"

      integer rtdb
      double precision e
      double precision scale_energy,scale_length
      double precision scale_force
      double precision evdw, evdw14
      character*30 message
      integer i,j
      logical qmcoords

      call mm_vdw_energy(nvtot,dbl_mb(i_rvdw),
     &                      nvdw,int_mb(i_ivdw),int_mb(i_jvdw),
     &                      dbl_mb(i_avdw),dbl_mb(i_bvdw),
     &                      evdw)

      call mm_vdw14_energy(nv14tot,dbl_mb(i_rvdw14),
     &                      nvdw14,int_mb(i_ivdw14),int_mb(i_jvdw14),
     &                      dbl_mb(i_avdw14),dbl_mb(i_bvdw14),
     &                      dbl_mb(i_scnb),evdw14)

      call util_convert_units("kcal","au",scale_energy)
 
      evdw = evdw*scale_energy
      evdw14 = evdw14*scale_energy

      e = e + evdw + evdw14

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
#include "mafdecls.fh"
#include "global.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_coords_data.fh"
      integer nt
      double precision r(3,nt)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      integer ng
      double precision e
      double precision evdw
      double precision g(3,ng)
c
      integer n, m
      integer i1,i2
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 pname
      character*30 message
      integer i_i, i_j
      integer i_indx

      pname = "mm_vdw_egrad"
c     write(*,*) pname

      evdw = 0.0d0

      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)
        i_i = 0
        i_j = 0
        do m = 1,nvtot
          i_indx = int_mb(i_icvdw+m-1)
          if(i1.eq.i_indx) i_i = m 
          if(i2.eq.i_indx) i_j = m 
        end do
 
        message = "checking bounds"
        if(i_i.eq.0 .or. i_j.eq.0) then
            goto 911
        end if

        dx = r(1,i_j) - r(1,i_i)
        dy = r(2,i_j) - r(2,i_i)
        dz = r(3,i_j) - r(3,i_i)
      
        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw = evdw + f12 - f6
        
        df = (12.0d0*f12 - 6.0d0*f6)*r2inv

        do m=1,nact
          i_indx = int_mb(i_iact+m-1)
          if(i1.eq.i_indx) then 
            g(1,m) = g(1,m) + df*dx
            g(2,m) = g(2,m) + df*dy
            g(3,m) = g(3,m) + df*dz
          end if
          if(i2.eq.i_indx) then
            g(1,m) = g(1,m) - df*dx
            g(2,m) = g(2,m) - df*dy
            g(3,m) = g(3,m) - df*dz
          end if
        end do
      end do
c     write(*,*) "evdw = ",evdw, " kcal/mol"

      return
911   call errquit("error "//trim(message),0,-1)
      end

      subroutine mm_vdw14_egrad(nt,r,
     &                            np,p1,p2,
     &                            a,b,scnb,
     &                            ng,evdw14,g)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "mm_geom_data.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_coords_data.fh"
#include "mm_link_data.fh"
      integer nt
      double precision r(3,nt)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      double precision scnb(np)
      integer ng
      double precision e
      double precision evdw14
      double precision g(3,ng)
c
      integer n, m, k
      integer i1,i2
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 message
      character*30 pname
      integer i_i, i_j
      integer i_indx, j_indx

      pname = "mm_vdw14_egrad"
c     write(*,*) pname

      evdw14 = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)
        i_i = 0
        i_j = 0
        do m = 1,nv14tot
          i_indx = int_mb(i_icvdw14+m-1)
          if(i1.eq.i_indx) i_i = m 
          if(i2.eq.i_indx) i_j = m 
        end do

c       TP: 1-4 VDW between QM and MMlnk is not scaled
c       do k = 1, nqm 
c         i_indx = int_mb(i_iqm+k-1)
c         do m = 1,nlink
c           j_indx = int_mb(i_lb+nlink+m-1)
c           write(*,*) 'TP: i j =', i_indx, j_indx
c           if(i1.eq.i_indx.and.i2.eq.j_indx) then 
c             write(*,*) 'i j scnb bf=', i_indx, j_indx, scnb(n)
c             scnb(n) = 1.0d0
c             write(*,*) 'i j scnb af=', i_indx, j_indx, scnb(n) 
c           else if(i1.eq.j_indx.and.i2.eq.i_indx) then
c             write(*,*) 'i j scnb bf=', i_indx, j_indx, scnb(n)  
c             scnb(n) = 1.0d0
c             write(*,*) 'i j scnb af=', i_indx, j_indx, scnb(n)  
c           end if       
c         end do
c       end do

        message = "checking bounds"
        if(i_i.eq.0 .or. i_j.eq.0) goto 911

        dx = r(1,i_j) - r(1,i_i)
        dy = r(2,i_j) - r(2,i_i)
        dz = r(3,i_j) - r(3,i_i)
      
        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw14 = evdw14 + (f12 - f6)*scnb(n)
       
        df = (12.0d0*f12 - 6.0d0*f6)*r2inv*scnb(n)

        do m=1,nact
          i_indx = int_mb(i_iact+m-1)
          if(i1.eq.i_indx) then 
            g(1,m) = g(1,m) + df*dx
            g(2,m) = g(2,m) + df*dy
            g(3,m) = g(3,m) + df*dz
          end if
          if(i2.eq.i_indx) then
            g(1,m) = g(1,m) - df*dx
            g(2,m) = g(2,m) - df*dy
            g(3,m) = g(3,m) - df*dz
          end if
        end do
      end do
c     write(*,*) "e_vdw14 = ",evdw14, " kcal/mol"

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
#include "mafdecls.fh"
#include "global.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_coords_data.fh"
      integer nt
      double precision r(3,nt)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      double precision e
      double precision evdw
c
      integer n,m
      integer i1,i2
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 message
      character*30 pname
      integer i_i, i_j
      integer i_indx

      pname = "mm_vdw_energy"
c     write(*,*) pname      
 
      evdw = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)
        i_i = 0
        i_j = 0
        do m = 1,nvtot
          i_indx = int_mb(i_icvdw+m-1)
          if(i1.eq.i_indx) i_i = m 
          if(i2.eq.i_indx) i_j = m 
        end do

        message = "checking bounds"
        if(i_i.eq.0 .or. i_j.eq.0) then
            goto 911
        end if

        dx = r(1,i_j) - r(1,i_i)
        dy = r(2,i_j) - r(2,i_i)
        dz = r(3,i_j) - r(3,i_i)

        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw = evdw + f12 - f6
        
      end do
c     write(*,*) "evdw = ",evdw, " kcal/mol"

      return
911   call errquit("error "//trim(pname)//trim(message),0,-1)
      end

      subroutine mm_vdw14_energy(nt,r,
     &                            np,p1,p2,
     &                            a,b,
     &                            scnb,evdw14)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "mm_coords_data.fh"
#include "mm_vdw_coords_data.fh"
      integer nt
      double precision r(3,nt)
      integer np
      integer p1(np)
      integer p2(np)
      double precision a(np)
      double precision b(np)
      double precision scnb(np)
      double precision e
      double precision evdw14
c
      integer n,m
      integer i1,i2
      double precision dx,dy,dz,df
      double precision rr
      double precision rinv
      double precision r2inv
      double precision r6,r12
      double precision f6,f12
      character*30 message
      character*30 pname
      integer i_i, i_j
      integer i_indx

      pname = "mm_vdw14_energy"
c     write(*,*) pname

      evdw14 = 0.0d0
      do n = 1, np
        i1 = p1(n)
        i2 = p2(n)
        i_i = 0
        i_j = 0
        do m = 1,nv14tot
          i_indx = int_mb(i_icvdw14+m-1)
          if(i1.eq.i_indx) i_i = m 
          if(i2.eq.i_indx) i_j = m 
        end do
        message = "checking bounds"
        if(i_i.eq.0 .or. i_j.eq.0) then
            goto 911
        end if

        dx = r(1,i_j) - r(1,i_i)
        dy = r(2,i_j) - r(2,i_i)
        dz = r(3,i_j) - r(3,i_i)

        rr = dx*dx + dy*dy + dz*dz
        rinv = sqrt(1.0d0/rr)
        r2inv = rinv*rinv     
        r6 = r2inv*r2inv*r2inv
        r12 = r6*r6
    
        f12 = a(n)*r12
        f6  = b(n)*r6
        evdw14 = evdw14 + (f12 - f6)*scnb(n)
        
      end do

c     write(*,*) "evdw14 = ",evdw14, " kcal/mol"
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

      subroutine mm_write_vdw(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_vdw_data.fh"

      integer lunit
      character*30 pname 

      integer i
      integer indx_i, indx_j
      double precision acoef
      double precision bcoef

      pname = "mm_write_vdw"
c     write(*,*) pname

      write(lunit,'(a)') 'vdw'
      write(lunit,'(a)') '#      i       j    A_coeff         B_coeff'
      write(lunit,'(a)') '#'

      do i=1,nvdw
        indx_i = int_mb(i_ivdw+i-1)
        indx_j = int_mb(i_jvdw+i-1)
        acoef = dbl_mb(i_avdw+i-1)
        bcoef = dbl_mb(i_bvdw+i-1)
        write(lunit,'(2(x,i7),2ES16.8)') indx_i, indx_j, acoef, bcoef
      end do

      write(lunit,'(a,/)') 'end'

      end

      subroutine mm_write_vdw14(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_vdw_data.fh"

      integer lunit
      character*30 pname 

      integer i
      integer indx_i, indx_j
      double precision acoef
      double precision bcoef
      double precision scnb
      double precision scee

      pname = "mm_write_vdw14"
c     write(*,*) pname

      write(lunit,'(a)') 'scaled_vdw'
      write(lunit,'(a)') '#      i       j    A_coeff         B_coeff'//
     >                   '          one_scnb'
c    >                   '          one_scnb        one_scee'
      write(lunit,'(a)') '#'

c     write(*,*) 'nvdw14 =', nvdw14
      do i=1,nvdw14
        indx_i = int_mb(i_ivdw14+i-1)
        indx_j = int_mb(i_jvdw14+i-1)
        acoef = dbl_mb(i_avdw14+i-1)
        bcoef = dbl_mb(i_bvdw14+i-1)
        scnb = dbl_mb(i_scnb+i-1)
c       scee = dbl_mb(i_scee+i-1)
        write(lunit,'(2(x,i7),3ES16.8)') indx_i, indx_j, acoef, bcoef, 
     >                                   scnb
c       write(lunit,'(2(x,i7),4ES16.8)') indx_i, indx_j, acoef, bcoef, 
c    >                                   scnb, scee
      end do

      write(lunit,'(a,/)') 'end'

      end
