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
      integer nqm0,nmm0
c
      pname = "mm_coords_init"
c
c     deallocate previous allocated just in case
      call mm_coords_deallocate()

      call mm_coords_nqm(rtdb,nqm0)
      call mm_coords_nmm(rtdb,nmm0)
c
      call mm_coords_allocate(nqm0,nmm0)
c     
      call mm_coords_qm_load(rtdb)
      call mm_coords_mm_load(rtdb)

c     write(*,*) 'out ', pname
c
c     call mm_coords_test()
c     call mm_coords_test2()
      return
911   call errquit("error "//trim(message),0,-1)
      return
      end

      subroutine mm_coords_qm_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb    

      character*180 message
      character*30 pname
      integer nqm0

      pname = "mm_coords_qm_init"
c     write(*,*) pname

      call mm_coords_nqm(rtdb,nqm0)
c     call mm_coords_qm_allocate(nqm0)
      call mm_coords_qm_load(rtdb)

      end

      subroutine mm_coords_mm_init(rtdb)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb    

      character*180 message
      character*30 pname
      integer nmm0

      pname = "mm_coords_qm_init"
c     write(*,*) pname

      call mm_coords_nmm(rtdb,nmm0)
c     call mm_coords_mm_allocate(nmm0)
      call mm_coords_mm_load(rtdb)

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
c TP: 

      pname = "mm_coords_allocate"
       
      ntot0 = nqm0 + nmm0
      if(ntot0.ne.ntot) then
        call mm_coords_end()
        if ( .not. ma_alloc_get(mt_dbl,3*ntot0, 
     +                          "total coords",
     +                           h_r,i_r)) goto 911
        if(.not.ma_alloc_get(mt_byte,16*ntot0,'mm tags',h_t,i_t))
     +        goto 911
      end if
    
      if(nqm0.ne.nqm) then
        call mm_coords_qm_end()
        if(.not.ma_alloc_get(mt_dbl,3*nqm0,"qm coords",h_rqm,i_rqm))
     &       goto 911

        if(.not.ma_alloc_get(mt_int,nqm0,"qm indices",h_iqm,i_iqm))
     &       goto 911

        if(.not.ma_alloc_get(mt_byte,16*nqm0,"qm tag",h_tqm,i_tqm))
     &       goto 911
      end if
 
      if(nmm0.ne.nmm) then
        call mm_coords_mm_end()
        if(.not.ma_alloc_get(mt_dbl,3*nmm0,"mm coords",h_rmm,i_rmm))
     &       goto 911
        if(.not.ma_alloc_get(mt_int,nmm0,"mm indices",h_imm,i_imm))
     &       goto 911
        if(.not.ma_alloc_get(mt_byte,16*nmm0,"mm tag",h_tmm,i_tmm))
     &       goto 911
        if(.not.ma_alloc_get(mt_dbl,nmm0,"mm charge",h_chgmm,
     &    i_chgmm)) goto 911
      end if
c
c     write(*,*) 'get_index'
c     if(.not.ma_get_index(h_r,i_r)) goto 911
c     if(.not.ma_get_index(h_t,i_t)) goto 911
c     if(.not.ma_get_index(h_rqm,i_rqm)) goto 911
c     if(.not.ma_get_index(h_rmm,i_rmm)) goto 911
c     if(.not.ma_get_index(h_iqm,i_iqm)) goto 911
c     if(.not.ma_get_index(h_imm,i_imm)) goto 911
c     if(.not.ma_get_index(h_tqm,i_tqm)) goto 911
c     if(.not.ma_get_index(h_tmm,i_tmm)) goto 911
c     if(.not.ma_get_index(h_chgmm,i_chgmm)) goto 911
      ntot = ntot0
      nqm  = nqm0
      nmm  = nmm0
c     write(*,*) 'nqm =', nqm, 'nmm= ', nmm
c     
      call dfill(3*ntot,0.0d0,dbl_mb(i_r),1)
      call dfill(3*nqm,0.0d0,dbl_mb(i_rqm),1)
      call dfill(3*nmm,0.0d0,dbl_mb(i_rmm),1)
      call ifill(nqm,0.0,int_mb(i_iqm),1)
      call ifill(nmm,0.0,int_mb(i_imm),1)
      call mm_coords_tags_reset(ntot,byte_mb(i_t))
c     
c     write(*,*) "out ", pname
      return
911   call errquit("error "//trim(message),0,-1)
      return
      end

      subroutine  mm_coords_deallocate()
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "mm_coords_data.fh"
      
        call mm_coords_end()
        call mm_coords_qm_end()
        call mm_coords_mm_end()

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
c     write(*,*) pname
c     write(*,*) ntot
      if(ntot.gt.0) then
          if (.not.ma_free_heap(h_t)) goto 911
c         write(*,*) 'done free h_t'
          if (.not.ma_free_heap(h_r)) goto 911
c         write(*,*) 'done free h_r'
          ntot = 0
c         nqm  = 0
c         nmm = 0
      end if

      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_coords_qm_end()
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
      pname = "mm_coords_qm_end"
c
c     write(*,*) pname
c     write(*,*) 'nqm = ', nqm
      if(nqm.gt.0) then
c         write(*,*) 'start free heaps'
          if (.not.ma_free_heap(h_iqm)) goto 911
c         write(*,*) 'done free iqm'
          if (.not.ma_free_heap(h_rqm)) goto 911
c         write(*,*) 'done free rqm'
          if (.not.ma_free_heap(h_tqm)) goto 911
          nqm = 0
      end if 

      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end
      
      subroutine mm_coords_mm_end()
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
      pname = "mm_coords_mm_end"
c
c     write(*,*) pname
      if(nmm.gt.0) then
          if (.not.ma_free_heap(h_rmm)) goto 911
c         write(*,*) 'done free rmm'
          if (.not.ma_free_heap(h_imm)) goto 911
c         write(*,*) 'done free imm'
          if (.not.ma_free_heap(h_tmm)) goto 911
          if (.not.ma_free_heap(h_chgmm)) goto 911
          nmm = 0
      end if 
 
      return
911   call errquit("error "//trim(pname),0,-1)
      return
      end

      subroutine mm_coords_nmm(rtdb,nmm)
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
      character*255 filename
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
c     prefix = "mm:coords:load"
c
      tag = "mm:crdparms:load"
c     if (.not. rtdb_cget(rtdb, tag, 1, filename))
c    $     goto 911
c
c     filename = "nwchem.mmcoords"
      if(.not.rtdb_cget(rtdb,tag,1,filename))
     & goto 911
c     write(*,*) "filename in mm_coords_nmm= ", filename

      call mm_open_file(filename,fn)
c     ----------------------
c     count numbers of lines
c     ----------------------
c     call util_xyz_nlines(nmm,fn)
      call inp_init(fn,LuOut)
      nmm = 0

5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'mm',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.6) then
                nmm = nmm + 1
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

c     write(*,*) 'nmm from crdparms file =', nmm

      return
911   call errquit("error "//trim(pname),0,
     >        -1)
      end

      subroutine mm_coords_nqm(rtdb,nqm)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "geom.fh"
#include "inp.fh"
#include "stdio.fh"

      integer nqm,rtdb
      external geom_rtdb_ncent
      logical geom_rtdb_ncent
      character*30 message
c
c     Marat's code commented out for now
c     ================================
c     message = "calling geom_rtdb_ncent"
c     if (.not. geom_rtdb_ncent(rtdb, 'geometry', nqm)) goto 911

c  TP : read qm from qmcoord file
      
      character*180 filename
      integer       nmax,n
c
      integer fn,k
      character*180 buffer
c     character*180 message
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
      character*255 field
c
      pname = "mm_coords_nqm"
c     filename = "nwchem.qmcoords"
      
      tag = "mm:crdparms:load"
      if(.not.rtdb_cget(rtdb,tag,1,filename))
     & goto 911
c     write(*,*) "filename in mm_coords_nqm= ", filename
     
      call mm_open_file(filename,fn)
c     ----------------------
c     count numbers of lines
c     ----------------------
c     call util_xyz_nlines(nqm,fn)
c     close(fn)  
      call inp_init(fn,LuOut)
      nqm = 0

5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'qm',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
                nqm = nqm + 1
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

c     write(*,*) 'nqm from qmcoord file =', nqm
c     write(*,*) 'nqm from crdparms file =', nqm
    
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
c     integer i_tmm
c     integer i_rmm,i_tmm
      character*16 mmtag
c      
c     i_rmm = i_r+3*nqm
c     i_tmm = i_t+16*nqm
c
      pname = "mm_coords_nmm"
c     write(*,*) pname     

c     prefix = "mm:coords:load"
      tag = "mm:crdparms:load"
c     filename = "nwchem.mmcoords"
      if(.not.rtdb_cget(rtdb,tag,1,filename))
     & goto 911
c
c     tag=trim(prefix)//":name"
c     if (.not. rtdb_cget(rtdb, tag, 1, filename))
c    $     goto 911
c
c     tag=trim(prefix)//":format"
c     if(.not. rtdb_get(rtdb,tag,mt_int,3,iloc))   
c    $   iloc = (/ 1, 2, 3 /) 
c
c     write(*,*) pname
      call mm_open_file(filename,fn)
      
c      
c     call util_xyz_read_coords(nmm,n,iloc,dbl_mb(i_rmm),fn)
      call inp_init(fn,LuOut)
      n = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
          if(inp_compare(.false.,'mm',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
                if(inp_a(field)) then
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.6) then
                n = n + 1
c               if(n > nqm) goto 911
                if(.not.inp_i(int_mb(i_imm+n-1)))       goto 911
                if(.not.inp_a(mmtag))                   goto 911
c               write(*,*) 'qmtag = ', qmtag
                call mm_coords_tag_set(mmtag,byte_mb(i_tmm+(n-1)*16))
c               write(*,*) 'tag qm = ', byte_mb(i_tqm+(n-1)*16)
                if(.not.inp_f(dbl_mb(i_rmm+(n-1)*3)))   goto 911
                if(.not.inp_f(dbl_mb(i_rmm+(n-1)*3+1))) goto 911
                if(.not.inp_f(dbl_mb(i_rmm+(n-1)*3+2))) goto 911
                if(.not.inp_f(dbl_mb(i_chgmm+n-1)))     goto 911
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
c      
c     tag=trim(prefix)//":units"
c     if (.not.rtdb_cget(rtdb, tag, 1,aunit)) aunit="ang"
c     call util_convert_units(aunit,"angstrom",scale)   
c     if (scale.ne.1.0) then
c       call dscal(3*nmm, scale,dbl_mb(i_rmm),1)
c     end if
c      mmtag = "MM"
c      call mm_coords_tags_set(nmm,"MM",byte_mb(i_tmm))
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
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
#include "geom.fh"
#include "mm_coords_data.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"

      integer rtdb
c
      integer n,fn
      integer geom
      integer nat
      logical status
      character*30 message
      character*30 pname
      double precision scale
      logical geom_cart_get1
      external geom_cart_get1
      character*180 filename
      character*255 field
      character*16 qmtag
c
c     load geometry 
c     -------------
c     if (.not. geom_create(geom, 'geometry'))
c    &     call errquit('cons_load_geom',0, GEOM_ERR)
c     if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
c    &     call errquit('cons_load_geom',0, RTDB_ERR)
c
c     get cart coordinates
c     --------------------
c     status=geom_ncent(geom,nat)
c     message = "nat not equals nmax"
c     write(*,*) 'aux_geom =', aux_geom
c     write(*,*) 'nat = ', nat
c     write(*,*) 'nqm = ', nqm
c     write(*,*) 'nlink = ', nlink
c     if(.not.aux_geom) 
c    &   nat = nat - nlink
c     if(nat.ne.nqm) then
c         write(*,*) nat,nqm
c         goto 911
c     end if
c     if(.not.status) 
c    & call errquit('cons_init: geom_create?',70, GEOM_ERR)
c     if(.not.geom_cart_get1(geom,nqm,byte_mb(i_t),dbl_mb(i_r)))
c    &    goto 911

c     if(.not.geom_destroy(geom)) 
c    &    goto 911
      pname = "mm_coords_qm_load"
c     write(*,*) pname

c     filename = "nwchem.qmcoords"
      if(.not.rtdb_cget(rtdb,"mm:crdparms:load",1,filename))
     & goto 911
c     write(*,*) 'filename in qm_load=', filename

      call mm_open_file(filename,fn)

      call inp_init(fn,LuOut)
      n = 0
5     continue
      if(inp_read()) then
        if(inp_a(field)) then
c         write(*,*) field
          if(inp_compare(.false.,'qm',field)) then
10          continue
            if(inp_read()) then
              if(inp_n_field().eq.1) then
c               write(*,*) 'nfield = 1'
                if(inp_a(field)) then
c                 write(*,*) field
                  if(inp_compare(.false.,'end',field)) goto 20
                end if
              else if(inp_n_field().eq.5) then
c               write(*,*) 'nfield = 5'
                n = n + 1
c               if(n > nqm) goto 911
                if(.not.inp_i(int_mb(i_iqm+n-1)))       goto 911
                if(.not.inp_a(qmtag))                   goto 911
c               write(*,*) 'qmtag = ', qmtag
                call mm_coords_tag_set(qmtag,byte_mb(i_tqm+(n-1)*16))
c               write(*,*) 'tag qm = ', byte_mb(i_tqm+(n-1)*16)
                if(.not.inp_f(dbl_mb(i_rqm+(n-1)*3)))   goto 911
                if(.not.inp_f(dbl_mb(i_rqm+(n-1)*3+1))) goto 911
                if(.not.inp_f(dbl_mb(i_rqm+(n-1)*3+2))) goto 911
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
     
c     write(*,*) 'out', pname
c     call util_convert_units("au","angstrom",scale)
c     call dscal(3*nat, scale,dbl_mb(i_r),1)
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

      subroutine mm_coords_test2()
      implicit none
#include "mm_coords_data.fh"
#include "mafdecls.fh"

      integer i,j
      
      do i=1,nqm
        write(6,'(1X,I5,3(1X,F10.6))') 
     &       int_mb(i_iqm+i-1), (dbl_mb(i_rqm+3*(i-1)+j),j=0,2)
      end do

      do i=1,nmm
        write(6,'(1X,I5,3(1X,F10.6))') 
     &       int_mb(i_imm+i-1), (dbl_mb(i_rmm+3*(i-1)+j),j=0,2)
      end do

      end 

      subroutine mm_coords_tags_set(n,tag,t)
      implicit none

      integer n
      character*(*) tag
      character*16 t(n)
      integer i

      do i=1,n
         t(i) = tag
c        t(i) = trim(tag)
      end do
      
      end

      subroutine mm_coords_tag_set(tag,t)
      implicit none

      character*16 tag
      character*16 t
      integer i
c        write(*,*) 'tag = ', tag
         t = tag
c        t = trim(tag)
c        write(*,*) 't = ', t
          
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

      subroutine mm_write_qmcoords(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_coords_data.fh"

      integer lunit
      integer i
      integer indx
      character*16 symbol
      double precision coord(3)

      character*30 pname
  
      pname = "mm_write_qmcoords"
c     write(*,*) pname

      write(lunit,'(a)') 'qm'
     
      do i=1,nqm
        indx = int_mb(i_iqm+i-1)
        coord(1) = dbl_mb(i_rqm+3*(i-1))
        coord(2) = dbl_mb(i_rqm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rqm+3*(i-1)+2)
        call mm_coords_tag_set(byte_mb(i_tqm+(i-1)*16),symbol)
        write(lunit,'(x,i7,a4,3F25.16)') indx, trim(symbol), coord(1:3)
c       write(lunit,'(x,i7,2x,a4,3F25.16)') indx, symbol, coord(1:3)
      end do

      write(lunit,'(a,/)') 'end'

      end

      subroutine mm_write_mmcoords(lunit)
      implicit none
#include "mafdecls.fh"
#include "mm_coords_data.fh"

      integer lunit
      integer i
      integer indx
      character*16 symbol
      double precision coord(3)
      double precision mmchg

      character*30 pname
  
      pname = "mm_write_mmcoords"
c     write(*,*) pname

      write(lunit,'(a)') 'mm'
     
      do i=1,nmm
        indx = int_mb(i_imm+i-1)
        coord(1) = dbl_mb(i_rmm+3*(i-1))
        coord(2) = dbl_mb(i_rmm+3*(i-1)+1)
        coord(3) = dbl_mb(i_rmm+3*(i-1)+2)
        mmchg = dbl_mb(i_chgmm+i-1)
        call mm_coords_tag_set(byte_mb(i_tmm+(i-1)*16),symbol)
c       write(lunit,'(x,i7,2x,a4,4F25.16)') indx, symbol, coord(1:3),
        write(lunit,'(x,i7,a4,4F25.16)') indx, trim(symbol), coord(1:3),
     >                                      mmchg
      end do

      write(lunit,'(a,/)') 'end'

      end
  
      subroutine mm_write_active_region_xyz(iunit, filename)
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "util.fh"
#include "msgids.fh"
#include "mm_link_data.fh"
#include "mm_geom_data.fh"

      integer iunit
      character*255 filename

      integer i, k
      character*16 tag

      if (ga_nodeid().eq.0) then
        open(iunit,file=filename,form='formatted')
        write(iunit,'(1x,i5)') nact
        write(iunit,'(a)')
        do i=1,nact
          call mm_coords_tag_set(byte_mb(i_tqml+16*(i-1)), tag)
          write(iunit,'(1x,a16,1x,3f15.8)') tag, 
     >                              (dbl_mb(i_rqml+3*(i-1)+k), k=0,2)
        end do
        close(iunit)
      end if

      end
