c $Id: bq_input.F,v 1.8 2006-10-06 06:56:59 marat Exp $



c     ***  The (current) format of the bq input is as folllows:
c     ***
c     ***   bq [units <string units default angstroms> ] 
c     ***      [(angstrom_to_au || ang2au) <real angstrom_to_au default 1.8897265>]
c     ***
c     ***       [bq_filename <string filename default bqcharges>]
c     ***       [bq_rtdb <string rtdb_charge_tag default cosmo:efcz> <string rtdb_coord_tag default cosmo:efcc>]
c     ***
c     ***
c     ***       { <real charge> <real x y z>
c     ***         ... }
c     ***
c     ***   end
c     ***
c     *** units - keyword specifying that a value will be entered by the user
c     *** for the string variable <units>. The default units for the geometry 
c     *** input are Angstrøms (Note: atomic units or Bohr are used within the
c     *** code, regardless of the option specified for the input units.  The 
c     *** code recognizes the following possible values  for the string variable <units>:         
c     ***   angstroms or an -- Angstroms (Å), the default (converts to A.U. using the Åto A.U. conversion factor)
c     ***   au or atomic or bohr -- Atomic units (A.U.)
c     ***   nm or nanometers -- nanometers (converts to A.U. using a conversion factor computed as  times the Å to A.U. conversion factor)
c     ***   pm or picometers -- picometers (converts to A.U. using a conversion factor computed as  times the Å to A.U. conversion factor)


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
      character*255 filename,rtdb_charge_tag,rtdb_coord_tag
      integer max_nbq
      logical status
      integer h_cbq,i_cbq
      integer h_qbq,i_qbq
      integer ncenter,ncenter_tmp
      double precision x,y,z,q
      integer i
      integer handle
      double precision  scale 
      double precision  charge_scale

c
      pname = "bq_input"

      call util_print_get_level(print_level)
      hprint = (ga_nodeid().eq.0).and.(print_high.le.print_level)

c
c      write(*,*) "in "//pname
c
c
      max_nbq = min( 4*7000, ma_inquire_avail(MT_INT))
      max_nbq = max_nbq/4
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

c      if(inp_n_field().gt.2) 
c     >  call errquit(pname//'too many fields',0,0)

      namespace = bq_default




c     *** get scale using the form used in geometry ***
      call bq_get_scalefrominput(scale)

c      scale = 1.88972598858d0
c      do i=1,inp_n_field()-1
c        call inp_set_field(i)
c        if (.not. inp_a(field))   
c     &    call errquit(pname//': inp_read failed',0,0)
c        if (inp_compare(.false.,field,'au')) then
c          scale = 1.0d0       
c        else if (inp_compare(.false.,field,'angstroms')) then
c          scale = 1.88972598858d0
c        else
c         namespace=field 
c        end if 
c      end do

     
      ncenter = 0
 20   if (inp_read()) then
       if(inp_a(field)) then

         if (inp_compare(.false., 'end', field)) then
           goto 30
         end if

         !*** read charges from file ***
         if (inp_compare(.false., 'load', field)) then
            filename = 'bq_charges'
            charge_scale = 1.0d0
            ncenter = 0
 200        continue
            if (inp_a(field)) then
              if (inp_compare(.false., 'filename', field)) then
                if (.not.inp_a(filename)) goto 998
                goto 200
              end if
              if (inp_compare(.false., 'scale', field)) then
                if (.not.inp_f(charge_scale)) goto 998
                goto 200
              end if
            end if
            call util_file_name_resolve(filename, .false.)

            call bq_readfile(filename,ncenter,
     &                       charge_scale,
     &                       dbl_mb(i_qbq),
     &                       dbl_mb(i_cbq))
            goto 20
         end if

         !*** read charges from rtdb ***
         if (inp_compare(.false., 'bq_rtdb', field)) then
           if (.not.inp_a(rtdb_charge_tag)) rtdb_charge_tag='cosmo:efcz'
           if (.not.inp_a(rtdb_coord_tag))  rtdb_coord_tag ='cosmo:efcc'

           call bq_fromrtdb(rtdb,rtdb_charge_tag,rtdb_coord_tag,
     &                      ncenter_tmp,
     &                      dbl_mb(i_qbq+ncenter),
     &                      dbl_mb(i_cbq+3*ncenter))
           ncenter = ncenter + ncenter_tmp
           goto 20
         end if

         call inp_set_field(0)
         status = .true.
         status = status .and. inp_f(q)
         status = status .and. inp_f(x)
         status = status .and. inp_f(y)         
         status = status .and. inp_f(z)         
         if (hprint) write(*,*) q,x,y,z
          if (.not. status) 
     &       call errquit('bq_input: inp_read failed',0,0)
         dbl_mb(i_qbq+ncenter) = q
         dbl_mb(i_cbq+3*ncenter)   = x
         dbl_mb(i_cbq+3*ncenter+1) = y
         dbl_mb(i_cbq+3*ncenter+2) = z
         ncenter = ncenter + 1
         status = ncenter .le. max_nbq 
         if (.not. status) 
     $       call errquit('bq_input: too many centers',0,0)
         go to 20
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
 
      return
998   call errquit('no token found '//field,0,
     >        RTDB_ERR)
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


*     *************************************************
*     *                                               *
*     *               bq_readfile                     *
*     *                                               *
*     *************************************************
      subroutine bq_readfile(filename,ncenter,s,charge,coord)
      implicit none
      character*(*) filename
      integer       ncenter
      double precision s
      real*8        charge(*)
      real*8        coord(3,*)

c
      double precision q,c(3)

      ncenter = 0
      open(unit=99,file=filename,status='old')
        do while(.true.)
           read(99,*,ERR=30,END=30) q,c(1),c(2),c(3)
           charge(ncenter+1)  = charge(ncenter+1)+s*q
           coord(1,ncenter+1) = c(1)
           coord(2,ncenter+1) = c(2)
           coord(3,ncenter+1) = c(3)
           ncenter = ncenter + 1
        end do
  30  close(99)
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

