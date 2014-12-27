      subroutine mm_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "mm_params.fh"
#include "mm.fh"
#include "stdio.fh"
c     
      integer rtdb !< [Input] The RTDB handle
c
      logical hprint
      integer print_level
      character*32 pname
      character*255 field
      character*255 namespace
      character*255 message
      character*255 aunit
      character*255 filename,rtdb_charge_tag,rtdb_coord_tag
      integer max_nbq
      logical status
      integer h_cbq,i_cbq
      integer h_qbq,i_qbq
      integer ncenter,ncenter_tmp
      double precision x,y,z,q
      integer i,ipos
      integer handle
      double precision  scale 
      double precision  charge_scale
      integer ip,np
      character*32 mtoken(10),tag
c
      pname = "mm_input"

      call util_print_get_level(print_level)
      hprint = (ga_nodeid().eq.0).and.(print_high.le.print_level)
c
c     ------------------
c     parse main bq line 
c     set namespace
c     ------------------
      namespace = mm_default
      call inp_set_field(1)
      if (inp_a(field)) then 
         namespace = field
      end if
      write(*,*) "namespace=",namespace

c    ---------------------------
c    begin input processing loop
c    ---------------------------
 20   if (inp_read()) then
       if(inp_a(field)) then
         if (inp_compare(.false., 'end', field)) then
            goto 30
C         else if (inp_compare(.false., 'clear', field)) then
C            call mm_rtdb_delete(rtdb,namespace)
C            goto 20
         else if (inp_compare(.false., 'load', field)) then
            np = inp_n_field()
            if(np.ne.3) goto 911
            do ip = 1,np-1
              if(.not.inp_a(mtoken(ip))) goto 911
            end do
            if (inp_compare(.false., 'vdw', mtoken(1))) then
              tag="mm:load:vdw"
              if (.not.rtdb_cput(rtdb,tag,1,mtoken(2))) 
     >              call errquit(pname//'failed to store'//tag,0,
     >               RTDB_ERR)
                call  mm_data_load_vdw(rtdb)
            else
               message = "unknown keyword "
               goto 911
            end if
            goto 20
         end if

       end if
      end if
 30   continue
c    -------------------------
c    end input processing loop
c    -------------------------

      return
911   call errquit(message,0,-1)
      end
C
C      subroutine mm_input_load_vdw(filename)
C      implicit none
C#include "util.fh"
C#include "errquit.fh"
C#include "inp.fh"
C#include "stdio.fh"
C#include "mafdecls.fh"
C      character*(*) filename
C      integer       nmax,n
Cc
C      integer fn,k
C      character*180 buffer,field
C      character*180 message
C      character*30 pname,atag
Cc
C      logical util_io_unit
C      external util_io_unit
C      integer h_i,i_i
C      integer h_j,i_j
C      integer h_a,i_a
C      integer h_b,i_b
C      integer max_vdw
Cc
C      pname = "mm_input_load_vdw"
Cc
C      max_vdw = 2500
Cc
C      if ( .not. ma_push_get( MT_INT, max_vdw, 
C     &      'vdw scratch i',
C     &      h_i, i_i) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      max_vdw, MA_ERR)
C
C      if ( .not. ma_push_get( MT_INT, max_vdw, 
C     &      'vdw scratch j',
C     &      h_j, i_j) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      max_vdw, MA_ERR)
C
C      if ( .not. ma_push_get( MT_DBL, max_vdw, 
C     &      'vdw scratch A',
C     &      h_a, i_a) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      max_vdw, MA_ERR)
C
C      if ( .not. ma_push_get( MT_DBL, max_vdw, 
C     &      'vdw scratch B',
C     &      h_b, i_b) ) call errquit(
C     &      'failed scratch space '//pname,
C     &      max_vdw, MA_ERR)
C
Cc
C      call mm_open_file(filename,fn)
Cc      
C      call inp_init(fn,LuOut)
C      n = 0
C10    continue
C      if (.not.inp_read()) goto 30
C      message = "checking for number of fields"
C      if(inp_n_field().ne.4) goto 911
C      message = "checking for bounds"
C      if(n.gt.max_vdw) goto 911
C      message = "reading vdw params"
C      n = n+1
C      write(*,*) "n=",n
C      if (.not.inp_i(int_mb(i_i+n-1))) goto 911
C      if (.not.inp_i(int_mb(i_j+n-1))) goto 911
C      if (.not.inp_f(dbl_mb(i_a+n-1))) goto 911
C      if (.not.inp_f(dbl_mb(i_b+n-1))) goto 911
C      goto 10
C30    continue
C      close(fn)
C      call inp_init(LuIn,LuOut)
C      if (.not. ma_pop_stack(h_b)) call errquit(
C     &      'failed to free scratch space '//pname,
C     &      max_vdw, MA_ERR)
C      if (.not. ma_pop_stack(h_a)) call errquit(
C     &      'failed to free scratch space '//pname,
C     &      max_vdw, MA_ERR)
C      if (.not. ma_pop_stack(h_j)) call errquit(
C     &      'failed to free scratch space '//pname,
C     &      max_vdw, MA_ERR)
C      if (.not. ma_pop_stack(h_i)) call errquit(
C     &      'failed to free scratch space '//pname,
C     &      max_vdw, MA_ERR)
C      return
C911   call errquit("error "//trim(message),0,
C     >        -1)
C      end
C
C      subroutine mm_open_file(filename,fn)
C      implicit none
C#include "util.fh"
C#include "errquit.fh"
C#include "inp.fh"
C#include "stdio.fh"
C      character*(*) filename
C      integer       fn
Cc
C      character*180 buffer
C      character*180 message
C      character*30 pname,atag
Cc
C      logical util_io_unit
C      external util_io_unit
Cc
C      pname = "mm_open_file"
Cc
C      if(.not.util_io_unit(80,90,fn))
C     +  call errquit(pname//"cannot get io unit",0,0)
Cc     first try to open file in the run directory
C      buffer = filename
C      message = "loading vdw information from "//buffer
C      open(unit=fn,file=buffer,status='old',form="formatted",ERR=10)
C      goto 800
C10    continue
Cc     now try perm directory
C      call util_file_name_resolve(buffer, .false.)
C      message = "loading vdw information from "//buffer
C      open(unit=fn,file=buffer,status='old',form="formatted",ERR=911)
C800   continue
C      write(luout,*) "Successfully "//trim(message)
C      write(luout,*)
C      return
C911   call errquit("error "//trim(message),0,
C     >        -1)
C      end
C
C
