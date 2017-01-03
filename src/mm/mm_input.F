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
c     namespace = mm_default
C     call inp_set_field(1)
C     if (inp_a(field)) then 
C        namespace = field
C     end if
C     write(*,*) "namespace=",namespace

c     ---------------------------
c     begin input processing loop
c     ---------------------------
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
c        TP: added
         else if (inp_compare(.false., 'crdparms', field)) then
             call mm_input_crdparms(rtdb)
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
c     character*255 afile_coords 
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
c     afile_coords  = " "
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
c     character*255 afile_coords 
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
c     afile_coords  = " "
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

      subroutine mm_input_crdparms(rtdb)
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
      integer i,k,n0,ipos
      character*30 pname
      character*84 tag

      pname = "mm_input_crdparms" 

10    continue
      if (inp_a(field)) then
        message = "processing field "//trim(field)
        if (inp_contains(.false.,'load',field,ipos)) then
          message = "reading load"
          if (.not. inp_a(svalue)) goto 911
          tag="mm:crdparms:load"
          if (.not.rtdb_cput(rtdb,tag,1,svalue)) 
     >          call errquit(pname//'failed to store'//tag,0,
     >           RTDB_ERR)
           goto 10
        end if
        goto 10
      end if
      return

911   call errquit("error "//message,0,
     >        -1)

      end
