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
c      namespace = mm_default
C      call inp_set_field(1)
C      if (inp_a(field)) then 
C         namespace = field
C      end if
C      write(*,*) "namespace=",namespace

c    ---------------------------
c    begin input processing loop
c    ---------------------------
 20   if (inp_read()) then
       if(inp_a(field)) then
         if (inp_compare(.false., 'end', field)) then
            goto 30
         else if (inp_compare(.false., 'coords', field)) then
             call mm_input_coords(rtdb)
             goto 20
         else if (inp_compare(.false., 'vdw', field)) then
             call mm_input_vdw(rtdb)
             goto 20
         else if (inp_compare(.false., 'load', field)) then
            np = inp_n_field()
            if(np.lt.3) goto 911
            do ip = 1,np-1
              if(.not.inp_a(mtoken(ip))) goto 911
            end do
            if (inp_compare(.false., 'vdw', mtoken(1))) then
              tag="mm:load:vdw"
              if (.not.rtdb_cput(rtdb,tag,1,mtoken(2))) 
     >              call errquit(pname//'failed to store'//tag,0,
     >               RTDB_ERR)
            else if (inp_compare(.false., 'positions', mtoken(1))) then
              tag="mm:load:positions"
              if (.not.rtdb_cput(rtdb,tag,1,mtoken(2))) 
     >              call errquit(pname//'failed to store'//tag,0,
     >               RTDB_ERR)
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

      subroutine mm_input_coords(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c     
      integer nmax,n
      character*255 message
      character*255 field
      character*255 svalue
c      character*255 afile_coords 
      character*255 afile_charges
      character*255 afile
      logical oscale
      integer iloc(4)
      double precision charge_scale
      integer i,k,n0,ipos
      integer dn
      double precision scale
      character*30 pname
      character*32 prefix
      character*84 tag

      pname = "mm_input_vdw" 
      prefix = "mm:coords:load"
      iloc = [2,3,4,5]
c      afile_coords  = " "
      afile_charges = " "
      afile = " "
      oscale = .false.
      charge_scale = 1.0

10    continue
      if (inp_a(field)) then
        message = "processing field "//trim(field)
        if (inp_contains(.false.,'load',field,ipos)) then
          message = "reading load"
          if (.not. inp_a(svalue)) goto 911
          tag=trim(prefix)//":name"
          if (.not.rtdb_cput(rtdb,tag,1,svalue)) goto 911
          goto 10
        else if (inp_compare(.false., 'format', field)) then
           do k=1,3
             if (.not.inp_i(iloc(k))) goto 911
           end do
           tag=trim(prefix)//":format"
           if(.not. rtdb_put(rtdb,tag,mt_int,3,iloc)) goto 911 
           goto 10
        else if (inp_compare(.false., 'units', field)) then
          message = "reading units "
          if (.not. inp_a(svalue)) goto 911
          tag=trim(prefix)//":units"
          if (.not.rtdb_cput(rtdb,tag,1,svalue)) goto 911
          goto 10
        else
          message = "uknown field "//trim(field)
          goto 911
        end if
        goto 10
      end if
      return

911   call errquit("error "//message,0,
     >        -1)

      end

      subroutine mm_input_vdw(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c     
      integer nmax,n
      character*255 message
      character*255 field
      character*255 svalue
c      character*255 afile_coords 
      character*255 afile_charges
      character*255 afile
      logical oscale
      integer iloc(4)
      double precision charge_scale
      integer i,k,n0,ipos
      integer dn
      double precision scale
      character*30 pname
      character*84 tag

      pname = "mm_input_vdw" 
      iloc = [2,3,4,5]
c      afile_coords  = " "
      afile_charges = " "
      afile = " "
      oscale = .false.
      charge_scale = 1.0

10    continue
      if (inp_a(field)) then
        message = "processing field "//trim(field)
c        if (inp_compare(.false., 'coords', field)) then
c           if (.not.inp_a(afile)) goto 911
c             goto 10
        if (inp_contains(.false.,'load',field,ipos)) then
          message = "reading load"
          if (.not. inp_a(svalue)) goto 911
          tag="mm:vdw:load"
          if (.not.rtdb_cput(rtdb,tag,1,svalue)) 
     >          call errquit(pname//'failed to store'//tag,0,
     >           RTDB_ERR)
           goto 10
        else if (inp_compare(.false., 'charges', field)) then
           if (.not.inp_a(afile_charges)) goto 911
             goto 10
        else if (inp_compare(.false., 'format', field)) then
           do k=1,4
             if (.not.inp_i(iloc(k))) goto 911
           end do
           goto 10
        else if (inp_compare(.false., 'scale', field)) then
           if (.not.inp_f(charge_scale)) goto 911
             oscale = .true.
             goto 10
        else
          if(afile.eq. " " ) afile = field 
        end if
        goto 10
      end if
      return

911   call errquit("error "//message,0,
     >        -1)

      end
C
C      subroutine mm_input_vdw(rtdb)
C      implicit none
C#include "errquit.fh"
C#include "inp.fh"
C#include "mafdecls.fh"
C#include "rtdb.fh"
C#include "util.fh"
C#include "global.fh"
C      integer rtdb
Cc     
C      integer nmax,n
C      character*255 message
C      character*255 field
C      character*255 svalue
Cc      character*255 afile_coords 
C      character*255 afile_charges
C      character*255 afile
C      logical oscale
C      integer iloc(4)
C      double precision charge_scale
C      integer i,k,n0,ipos
C      integer dn
C      double precision scale
C      character*30 pname
C      character*84 tag
C
C      pname = "mm_input_vdw" 
C      iloc = [2,3,4,5]
Cc      afile_coords  = " "
C      afile_charges = " "
C      afile = " "
C      oscale = .false.
C      charge_scale = 1.0
C
C10    continue
C      if (inp_a(field)) then
C        message = "processing field "//trim(field)
Cc        if (inp_compare(.false., 'coords', field)) then
Cc           if (.not.inp_a(afile)) goto 911
Cc             goto 10
C        if (inp_contains(.false.,'load',field,ipos)) then
C          message = "reading load"
C          if (.not. inp_a(svalue)) goto 911
C          tag="mm:load:vdw"
C          if (.not.rtdb_cput(rtdb,tag,1,svalue)) 
C     >          call errquit(pname//'failed to store'//tag,0,
C     >           RTDB_ERR)
C           goto 10
C        else if (inp_compare(.false., 'charges', field)) then
C           if (.not.inp_a(afile_charges)) goto 911
C             goto 10
C        else if (inp_compare(.false., 'format', field)) then
C           do k=1,4
C             if (.not.inp_i(iloc(k))) goto 911
C           end do
C           goto 10
C        else if (inp_compare(.false., 'scale', field)) then
C           if (.not.inp_f(charge_scale)) goto 911
C             oscale = .true.
C             goto 10
C        else
C          if(afile.eq. " " ) afile = field 
C        end if
C        goto 10
C      end if
C      return
C
C911   call errquit("error "//message,0,
C     >        -1)
C
C      end
C
