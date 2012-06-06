c $Id$

      subroutine bq_input(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "bq_params.fh"
#include "bq.fh"
c     
      integer rtdb
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

c
      pname = "bq_input"

      call util_print_get_level(print_level)
      hprint = (ga_nodeid().eq.0).and.(print_high.le.print_level)
c
      max_nbq = 25000  ! temporary limit for now
c
      if ( .not. ma_push_get( MT_DBL, 3*max_nbq, 
     &      'cbq scratch',
     &      h_cbq, i_cbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      max_nbq, MA_ERR)
c
      call dfill(3*max_nbq,0.0d0,dbl_mb(i_cbq),1)
c
      if ( .not. ma_push_get( MT_DBL, max_nbq, 
     &      'qbq scratch',
     &      h_qbq, i_qbq) ) call errquit(
     &      'failed scratch space '//pname,
     &      max_nbq, MA_ERR)
c
      call dfill(max_nbq,0.0d0,dbl_mb(i_qbq),1)


c     ------------------
c     parse main bq line 
c     ------------------
      namespace = bq_default
      aunit = "ang"
      call inp_set_field(1)
10    continue
      if (inp_a(field)) then 
        if (inp_contains(.false.,'unit',field,ipos)) then
          message = "reading units"
          if (.not. inp_a(aunit)) goto 911  
          goto 10
        else
          if(namespace.eq.bq_default) namespace = field
          goto 10
        end if 
      end if

      call util_length_scale(aunit,scale)
     
      ncenter = 0
 20   if (inp_read()) then
       if(inp_a(field)) then

         if (inp_compare(.false., 'end', field)) then
           goto 30

         else if (inp_compare(.false., 'load', field)) then
            call  bq_input_load(rtdb,max_nbq,ncenter,
     &                       dbl_mb(i_qbq),
     &                       dbl_mb(i_cbq))
            goto 20
         else if (inp_compare(.false., 'bq_rtdb', field)) then
           if (.not.inp_a(rtdb_charge_tag)) rtdb_charge_tag='cosmo:efcz'
           if (.not.inp_a(rtdb_coord_tag))  rtdb_coord_tag ='cosmo:efcc'

           call bq_fromrtdb(rtdb,rtdb_charge_tag,rtdb_coord_tag,
     &                      ncenter_tmp,
     &                      dbl_mb(i_qbq+ncenter),
     &                      dbl_mb(i_cbq+3*ncenter))
           ncenter = ncenter + ncenter_tmp
           goto 20

         else
           call inp_set_field(0)
           status = .true.
           status = status .and. inp_f(x)
           status = status .and. inp_f(y)         
           status = status .and. inp_f(z)         
           status = status .and. inp_f(q)
           if (hprint) write(*,*) q,x,y,z
            if (.not. status) 
     &         call errquit('bq_input: inp_read failed',0,0)
           dbl_mb(i_qbq+ncenter) = q
           dbl_mb(i_cbq+3*ncenter)   = x
           dbl_mb(i_cbq+3*ncenter+1) = y
           dbl_mb(i_cbq+3*ncenter+2) = z
           ncenter = ncenter + 1
           status = ncenter .le. max_nbq 
           if (.not. status) 
     $         call errquit('bq_input: too many centers',0,0)
           go to 20
       end if
      end if
      end if
 30   continue

      call dscal(3*ncenter, scale, dbl_mb(i_cbq), 1)

      if(.not.bq_create(namespace,handle))   
     &       call errquit('bq_create failed',0,0)

      if(.not.bq_set(handle,ncenter,dbl_mb(i_qbq),dbl_mb(i_cbq)))
     &       call errquit('bq_set failed',0,0)

      call bq_print_info(handle)

      if(.not.bq_rtdb_store(rtdb,handle))
     &       call errquit('bq_store failed',0,0)

      if(.not.bq_destroy(handle))
     &       call errquit('bq_destroy failed',0,0)

c
c     activate default namespace if found
c     -----------------------------------
      if(ncenter.ne.0 .and. namespace.eq.bq_default) then
      if(.not. rtdb_cput(rtdb,"bq" , 1,bq_default)) then
         call errquit(pname//'failed to set default namespace',
     >                0,RTDB_ERR)
      end if
      end if

      if (.not. ma_pop_stack(h_qbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      max_nbq, MA_ERR)

      if (.not. ma_pop_stack(h_cbq)) call errquit(
     &      'failed to free scratch space '//pname,
     &      max_nbq, MA_ERR)
 
      stop
      return
998   call errquit('no token found '//field,0,
     >        RTDB_ERR)
911   call errquit(message,0,-1)
      end

      subroutine bq_input_load(rtdb,nmax,n,q,c)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "util.fh"
#include "global.fh"
#include "bq_params.fh"
#include "bq.fh"
c     
      integer rtdb
      integer nmax,n
      double precision q(nmax),c(3,nmax)
      character*255 message
      character*255 field
      character*255 afile_coords 
      character*255 afile_charges
      character*255 afile
      logical oscale
      integer iloc(4)
      double precision charge_scale
      integer i,k,n0

      
      iloc = [2,3,4,5]
      afile_coords  = " "
      afile_charges = " "
      afile = " "
      oscale = .false.
      charge_scale = 1.0
10    continue
      if (inp_a(field)) then
        message = "processing field "//trim(field)
        if (inp_compare(.false., 'coords', field)) then
           if (.not.inp_a(afile_coords)) goto 998
             goto 10
        else if (inp_compare(.false., 'charges', field)) then
           if (.not.inp_a(afile_charges)) goto 998
             goto 10
        else if (inp_compare(.false., 'format', field)) then
           do k=1,4
             if (.not.inp_i(iloc(k))) goto 998
           end do
           goto 10
        else if (inp_compare(.false., 'scale', field)) then
           if (.not.inp_f(charge_scale)) goto 998
             oscale = .true.
             goto 10
        else
          afile = field 
        end if
      end if

      if(afile_coords.ne." "
     >   .and. 
     >   afile_charges.ne." ") then
        n0=n
        call bq_load_xyz(afile_coords,nmax,n,iloc(1:3),c)
        n = n0
        call bq_load_q(afile_charges,nmax,n,iloc(4),q)
      else if (afile.ne." ") then
        call bq_load(afile,nmax,n,iloc,q,c)
      else
        message = "give me some files now"
        goto 998
      end if

      if(oscale) then
        call dscal(n, charge_scale, q, 1)
      end if

c      do i=1,n
c       write(*,*) c(1,i),c(2,i),c(3,i),q(i)
c      end do
c      write(*,*) "oscale ",oscale
c      write(*,*) "scale ",charge_scale
c      write(*,*) "afile ",afile
c      write(*,*) "afile_coords ",afile_coords
c      write(*,*) "afile_charges ",afile_charges

      return

998   call errquit("error "//message,0,
     >        -1)

      end

*     *************************************************
*     *                                               *
*     *               bq_get_scalefrominput           *
*     *                                               *
*     *************************************************
      subroutine bq_get_scalefrominput(scale)
      implicit none
      double precision scale

#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ***
      double precision bohr,ang2au,newang2au
      parameter (bohr=0.529177249d0)
      parameter (ang2au = 1.0d0/bohr)   ! ang2au =1.88972598858d0

      integer ind,unit_ind
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 3)
      character*22 dirs(num_dirs)
      character*255 test

      integer nunits
      parameter (nunits=9)
      character*10 units(nunits)

      data dirs / 'units',
     >            'angstrom_to_au','ang2au' /

      data units / 'au','atomic','bohr',
     >             'angstroms', 'an',
     >             'nanometers','nm',
     >             'picometers','pm' /




      unit_ind  = 4 !*** default scale is angstroms ***
      newang2au = ang2au
      call inp_set_field(1)
 40   if (inp_a(test)) then

         if (inp_match(num_dirs, .false., test, dirs, ind)) then
            goto (41,42,42) ind

*           **** units ****
  41        if (.not. inp_a(test)) call errquit
     >       ('simulation_cell_input: error reading units', 0,
     &       INPUT_ERR)
            if (.not.inp_match(nunits, .false., test, units, unit_ind))
     >       call errquit('simulation_cell_input: error reading units',
     $       0,INPUT_ERR)
            goto 40

*           **** angstrom_to_au ****
  42        if (.not. inp_f(newang2au)) newang2au = ang2au
            goto 40
         end if
      end if

*     **** set the scale ****
      goto ( 52,52,52,53,53,54,54,55,55 ) unit_ind
  52  scale = 1.0d0 !**** atomic units ****
      goto 60
  53  scale = 1.0d0*newang2au !**** angstrom units ****
      goto 60
  54  scale = 10.0d0*newang2au !**** nanometer units ****
      goto 60
  55  scale = 1.0d0*newang2au/100.0d0 !**** picometer units ****
      goto 60

  60  continue
      return
      end

      subroutine bq_get_length_scale(aunit,scale)
      implicit none
      double precision scale
      character*(*) aunit

#include "inp.fh"
#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ***
      double precision bohr,ang2au,newang2au
      parameter (bohr=0.529177249d0)
      parameter (ang2au = 1.0d0/bohr)   ! ang2au =1.88972598858d0

      integer ind,unit_ind
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 3)
      character*22 dirs(num_dirs)
      character*10 test
      character*255 message

      integer nunits
      parameter (nunits=9)
      character*10 units(nunits)

      data units / 'au','atomic','bohr',
     >             'angstroms', 'an',
     >             'nanometers','nm',
     >             'picometers','pm' /


       test = aunit

       message = "matching unit "//test
       if (.not.inp_match(nunits, .false., test, units, unit_ind))
     >    goto 59

*     **** set the scale ****
      goto ( 52,52,52,53,53,54,54,55,55 ) unit_ind
  52  scale = 1.0d0 !**** atomic units ****
      goto 60
  53  scale = 1.0d0*newang2au !**** angstrom units ****
      goto 60
  54  scale = 10.0d0*newang2au !**** nanometer units ****
      goto 60
  55  scale = 1.0d0*newang2au/100.0d0 !**** picometer units ****
      goto 60

  59  call errquit(trim(message),0,-1) 
  60  continue
      return
      end

      subroutine bq_load_xyz(filename,nmax,n,iloc,c)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       nmax,n
      integer iloc(3)
      double precision   c(3,nmax)

c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*30 pname,atag
      integer iloc_max
c
      logical util_io_unit
      external util_io_unit
      double precision f3(3)
c
      pname = "bq_load_xyz"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
      call util_file_name_resolve(filename, .false.)
      open(unit=fn,file=filename,status='old',form="formatted")
c      
      iloc_max = maxval(iloc)
c      
      call inp_init(fn,LuOut)
10    continue
      if (.not.inp_read()) goto 30
      if(inp_n_field().lt.iloc_max) goto 10
      message = "checking for bounds"
      if(n.gt.nmax) goto 911
      message = "reading coordinates"
      do k=1,3
        call inp_set_field(iloc(k)-1) 
        if (.not.inp_f(f3(k))) goto 10
      end do
      n = n+1
      c(:,n) = f3
      goto 10
30    continue
      close(fn)
      call inp_init(LuIn,LuOut)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine bq_load_q(filename,nmax,n,iloc,q)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       nmax,n
      integer iloc
      double precision   q(nmax)
      double precision   c(3,nmax)

c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*30 pname,atag
      integer iloc_max
      double precision f
c
      logical util_io_unit
      external util_io_unit
c
      pname = "bq_load_q"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
      call util_file_name_resolve(filename, .false.)
      open(unit=fn,file=filename,status='old',form="formatted")
c      
      iloc_max = iloc
c      
      call inp_init(fn,LuOut)
10    continue
      if (.not.inp_read()) goto 30
      if(inp_n_field().lt.iloc_max) goto 10
      message = "checking for bounds"
      if(n.gt.nmax) goto 911
      message = "reading charges"
      call inp_set_field(iloc-1) 
      if (.not.inp_f(f)) goto 10
      n=n+1
      q(n) = f
      goto 10
30    continue
      call inp_init(LuIn,LuOut)
      close(fn)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end

      subroutine bq_load(filename,nmax,n,iloc,q,c)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      integer       nmax,n
      integer iloc(4)
      double precision   q(nmax)
      double precision   c(3,nmax)

c
      integer fn,k
      character*180 buffer,field
      character*180 message
      character*30 pname,atag
      integer iloc_max
      double precision f4(4)
c
      logical util_io_unit
      external util_io_unit
c
      pname = "bq_load"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
      call util_file_name_resolve(filename, .false.)
      open(unit=fn,file=filename,status='old',form="formatted")
c      
      iloc_max = maxval(iloc)
c      
      call inp_init(fn,LuOut)
10    continue
      if (.not.inp_read()) goto 30
      if(inp_n_field().lt.iloc_max) goto 10
      message = "checking for bounds"
      if(n.gt.nmax) goto 911
      message = "reading coordinates"
      do k=1,3
        call inp_set_field(iloc(k)-1) 
        if (.not.inp_f(f4(k))) goto 10
      end do
      k = 4
      call inp_set_field(iloc(k)-1) 
      if (.not.inp_f(f4(k))) goto 10
      n=n+1
      c(:,n) = f4(1:3)
      q(n) = f4(4)
      goto 10
30    continue
      call inp_init(LuIn,LuOut)
      close(fn)
      return
911   call errquit("error "//trim(message),0,
     >        -1)
      end
c      
      subroutine bq_readfile(filename,ncenter,q,c)
      implicit none
#include "util.fh"
#include "errquit.fh"
#include "inp.fh"
      character*(*) filename
      integer       ncenter
      double precision   q(ncenter)
      double precision   c(3,ncenter)

c
      integer fn,fn1,k,i,ipos
      character*180 buffer
      character*30 pname,atag
c
      logical util_io_unit
      external util_io_unit
c
      pname = "bq_readfile"
c
      if(.not.util_io_unit(80,90,fn))
     +  call errquit(pname//"cannot get io unit",0,0)
      open(unit=fn,file=filename,status='old',form="formatted")
      i=0
10    continue        
       read(fn,'(A180)',ERR=30,END=30) buffer
       if(inp_contains(.false.,"bq",buffer,ipos)) then
         i=i+1
         if(i.gt.ncenter)
     +     call errquit(pname//"increase bq stack",i,0)
         read(buffer,*) atag,(c(k,i),k=1,3),q(i)
       end if
      goto 10
  30  close(fn)
      ncenter = i
      return
      end
*     *************************************************
*     *                                               *
*     *               bq_fromrtdb                     *
*     *                                               *
*     *************************************************
      subroutine bq_fromrtdb(rtdb,rtdb_charge_tag,rtdb_coord_tag,
     &                      ncenter,charge,coord)
      implicit none
      integer       rtdb
      character*(*) rtdb_charge_tag
      character*(*) rtdb_coord_tag
      integer       ncenter
      real*8        charge(*)
      real*8        coord(*)

#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"

*     **** local variables ****
      integer      ncenter_tmp,ncenter3_tmp,matype1,matype2
      character*26 date
      
      ncenter = 0
      if (rtdb_get_info(rtdb,rtdb_charge_tag,
     &                  matype1,ncenter_tmp,date).and.
     &   rtdb_get_info(rtdb,rtdb_coord_tag,
     &            matype2,ncenter3_tmp,date)) then

        if ((ncenter3_tmp .eq. (3*ncenter_tmp)) .and.
     &      (matype1.eq.matype2)                .and.
     &      (matype1.eq.mt_dbl))                then

          !*** get charges ***
          if (.not.rtdb_get(rtdb,rtdb_charge_tag,
     &                      mt_dbl,ncenter_tmp,charge))
     &      goto 35

          !*** get coordinates ***
          if (.not.rtdb_get(rtdb,rtdb_coord_tag,
     &                      mt_dbl,ncenter3_tmp,coord))
     &      goto 35

            ncenter = ncenter + ncenter_tmp
        end if
      end if

  35  continue
      return
      end

