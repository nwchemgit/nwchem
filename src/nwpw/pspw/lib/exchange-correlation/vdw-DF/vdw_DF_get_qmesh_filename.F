
*     *****************************************
*     *                                       *
*     *         vdw_DF_get_qmesh_filename     *
*     *                                       *
*     *****************************************
*
*     This function returns the filename of the qmesh datafile.
*
*     *** Order of precedence for choosing name                     ***
*     *** 1) value of NWCHEM_QMESH_DATA environment variable        ***
*     *** 2) value of NWCHEM_QMESH_DATA set in $HOME/.nwchemrc file ***
*     *** 3) value of the compiled in library name                  ***
*
*     This is a serial io routine
*
      subroutine vdw_DF_get_qmesh_filename(qmesh_data_name)
      implicit none
      character*(*) qmesh_data_name

#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"

*     **** local variables ****
      logical mprint,hprint,debug,does_it_exist
      logical from_environment,from_compile,from_nwchemrc
      integer iop,lgth,unitf,print_level,i,j
      character*255 qmesh_library

*     **** external functions ****
      logical  util_find_dir
      external util_find_dir

      call util_print_get_level(print_level)
      mprint = print_medium.le.print_level
      hprint = print_high  .le.print_level
      debug  = print_debug .le.print_level
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.

*     **** Try to get from NWCHEM_QMESH_DATA environment variable ****
      call util_getenv('NWCHEM_QMESH_DATA',qmesh_data_name)
      lgth=inp_strlen(qmesh_data_name)
      if (lgth.gt.0) then
         if (util_find_dir(qmesh_data_name)) then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &       qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         end if
      end if


*     **** Try to get from NWCHEM_QMESH_DATA nwchemrc ****
*2:   check for NWCHEM_QMESH_DATA defined in users .nwchemrc file
*     assumed structure in .nwchemrc file is variable [whitespace] value
*     one setting per line
*
      qmesh_library='nwchem_qmesh_data'
      call inp_save_state() ! save state of any inp unit
      if(.not.util_nwchemrc_get(qmesh_library,qmesh_data_name)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        does_it_exist=util_find_dir(qmesh_data_name)
        if (does_it_exist)then
          from_nwchemrc = .true.
          call inp_restore_state() ! restore state of any inp unit
          goto 99
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_QMESH_DATA set to: <',
     &     qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      call inp_restore_state() ! restore state of any inp unit



*     **** Try to get from compile ****
      from_compile = .true.
      call util_nwchem_srcdir(qmesh_data_name)
      qmesh_data_name
     >     =qmesh_data_name(1:inp_strlen(qmesh_data_name))
     >     //"/nwpw/pspw/lib/exchange-correlation/vdw-DF/vdw_qmesh.dat"


 99   continue

      if (from_environment) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: environment'
      else if (from_nwchemrc) then
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: .nwchemrc'
      else
          if (debug)
     >     write(luout,*)
     >     ' nwchem_qmesh_data name resolved from: compiled reference'
      endif
      if (debug) then
         write(luout,*) ' NWCHEM_QMESH_DATA set to: <',
     >    qmesh_data_name(1:inp_strlen(qmesh_data_name)),'>'
      end if

      return
      end

