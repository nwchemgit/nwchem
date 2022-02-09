*
* $Id$
*

*
* order of precidence for choosing name
* 1) value of NWCHEM_NWPW_LIBRARY environment variable
* 2) value of the compiled in library name
*      
c the string pspw_default is appended by default for pspw
c the string paw_default is appended by default for paw
c
      subroutine nwpw_libfile(flag,filename,element,library_type)
      implicit none

#include "inp.fh"
#include "btdb.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "util.fh"
#include "bafdecls.fh"


      integer       flag
      character*(*) filename !RESULT
      character*(*) element
      integer       library_type
c
      logical mprint,hprint,debug,lexist
      logical does_it_exist,foundit,is_it_open,noslash
      logical nwrcopen,from_environment,from_nwchemrc,from_compile

      character*(nw_max_path_len) libname,rtdb_name,library,library_name
      character*128 test
      character*24  basis_library
      integer lgth,l,ind,rtdb,unitrc

*     **** external functions ****
      logical  util_find_dir
      external util_find_dir
      logical  control_print
      external control_print
      integer  control_rtdb
      external control_rtdb
c

*
* order of precedence for choosing name
* 1) value of NWCHEM_NWPW_LIBRARY environment variable
* 2) value of NWCHEM_NWPW_LIBRARY set in $HOME/.nwchemrc file
* 3) value of the compiled in library name
*
      mprint = control_print(print_medium)
      hprint = control_print(print_high)
      debug  = control_print(print_debug)
      nwrcopen         = .false.
      from_environment = .false.
      from_nwchemrc    = .false.
      from_compile     = .false.
      noslash          = .false.
      unitrc           = -1

*     **** Try to get from NWCHEM_NWPW_LIBRARY environment variable ****
      call util_getenv('NWCHEM_NWPW_LIBRARY',libname)
      lgth=inp_strlen(libname)
      if (lgth.gt.0) then
         if (util_find_dir(libname)) then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_NWPW_LIBRARY set to: <',
     &           libname(1:inp_strlen(libname)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         end if
      end if

*     **** Try to get from NWCHEM_NWPW_LIBRARY nwchemrc ****
*2: check for NWCHEM_NWPW_LIBRARY defined in users .nwchemrc file
*   assumed structure in .nwchemrc file is variable [whitespace] value
*   one setting per line
*
      basis_library='nwchem_nwpw_library'
      call inp_save_state() ! save state of any inp unit
      if (.not.util_nwchemrc_get(basis_library,libname)) then
        if (debug) then
          write(luout,*)'util_nwchemrc_get failed'
        endif
      else
        lgth = inp_strlen(libname)
        does_it_exist = .false.
        if (libname(lgth:lgth).eq.'/') then
          does_it_exist=util_find_dir(libname)
        else
          inquire(file=libname,exist=does_it_exist)
          noslash=does_it_exist
        endif
        if (does_it_exist)then
          from_nwchemrc = .true.
          call inp_restore_state() ! restore state of any inp unit
          goto 99
        else
          write(luout,*)' warning:::::::::::::: from_nwchemrc'
          write(luout,*)' NWCHEM_NWPW_LIBRARY set to: <',
     &      libname(1:inp_strlen(libname)),'>'
          write(luout,*)' but file does not exist !'
          write(luout,*)' using compiled in library'
        endif
      endif
      call inp_restore_state() ! restore state of any inp unit

*     **** Try to get from compile ****
      from_compile = .true.
      call util_nwchem_srcdir(libname)
c     libname = 
c    &NWPW_LIBRARY
c
c     add here nwpw/libraryps bit
c
      libname=libname(1:inp_strlen(libname))
     &     //"/nwpw/libraryps/"
      if (util_find_dir(libname)) then
         goto 99
      else
         write(luout,*)' warning:::::::::::::: from_compile'
         write(luout,*)' NWCHEM_NWPW_LIBRARY is: <',
     &        libname(1:inp_strlen(libname)),'>'
         write(luout,*)' but file does ','not exist or you ',
     &        'do not have ','access to it !'
      endif
      call errquit('nwpwlibfile: no nwpw library found',0, DISK_ERR)

 99   continue

      if (from_environment) then
          if (mprint)
     >     write(luout,*)' library name resolved from: environment'
      else if (from_nwchemrc) then
          if (mprint)
     >     write(luout,*)' library name resolved from: .nwchemrc'
      else
          if (mprint)
     >     write(luout,*)
     >      ' library name resolved from: compiled reference'
      endif
      if (mprint) then 
         write(luout,*) ' NWCHEM_NWPW_LIBRARY set to: <',
     >                  libname(1:inp_strlen(libname)),'>'
      end if

      if(noslash) then
         write(luout,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(luout,*) "! DANGEROUS: nwchem_nwpw_library set equal  !"
         write(luout,*) "! to a filename and not to a directory since!"
         write(luout,*) "! since you did not use the ending slash.   !"
         write(luout,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      end if
      if (nwrcopen) then
        close(unit=unitrc,status='keep') ! close file after library name found
      end if

*     **** library for paw ****
      if (flag.eq.2) then

        library_type = 0
        rtdb = control_rtdb()
        l = index(element,' ')  -1
        if (l.lt.1) l = 2
        rtdb_name = element(1:l)//':paw_library'
        if (.not.btdb_cget(rtdb,rtdb_name,1,library_name)) then
           rtdb_name = '*'//':paw_library'
           if (.not.btdb_cget(rtdb,rtdb_name,1,library_name)) then
              library_name = 'paw_default'
           end if
         end if

        ind = index(library_name,' ')  -1
        library = library_name(1:ind)//'/'
        ind = index(library,' ')  -1

c        filename=
c     >   libname(1:inp_strlen(libname))//library(1:ind)//element(1:l)
c     >   //'_basis'
        if (library(1:1).eq.'/') then
             filename=library(1:ind)//element(1:l)
        else
           filename=
     >   libname(1:inp_strlen(libname))//library(1:ind)//element(1:l)
       end if

        if (hprint) then
        write(luout,*)
     >   ' library file name is: ',filename(1:inp_strlen(filename))
        write(luout,*)
        end if


*     **** library for pspw and band ****
      else

        rtdb = control_rtdb()
        l = index(element,' ')  -1
        if (l.lt.1) l = 2
        rtdb_name = element(1:l)//':pspw_library'
        if (.not.btdb_cget(rtdb,rtdb_name,1,library_name)) then
           rtdb_name = '*'//':pspw_library'
           if (.not.btdb_cget(rtdb,rtdb_name,1,library_name))
     >       library_name = 'pspw_default'
        end if

        rtdb_name = element(1:l)//':pspw_library_type'
        if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,library_type)) then
           rtdb_name = '*'//':pspw_library_type'
           if (.not.btdb_get(rtdb,rtdb_name,mt_int,1,library_type))
     >        library_type = 0
        end if


        !*** library directory ***
        if (library_type.eq.0) then
          ind = index(library_name,' ')  -1
          library = library_name(1:ind)//'/'
          ind = index(library,' ')  -1

          if (library(1:1).eq.'/') then
             filename=library(1:ind)//element(1:l)
          else
             filename=
     >     libname(1:inp_strlen(libname))//library(1:ind)//element(1:l)
          end if

        !*** cpi filename ***
        else if (library_type.eq.2) then
          ind = index(library_name,' ')  -1
          filename= library_name(1:ind)
          call util_file_name_resolve(filename, .false.)
        !*** teter filename ***
        else if (library_type.eq.3) then
          ind = index(library_name,' ')  -1
          filename= library_name(1:ind)
          call util_file_name_resolve(filename, .false.)
        end if


        if (hprint) then
        write(luout,*)
     >   ' library file name is: ',filename(1:inp_strlen(filename))
        write(luout,*) 
        end if

      end if

c
c     check if file exists
c
      inquire(file=filename,exist=lexist)
      if(.not.lexist) then
         write(luout,*) ' could not find library for ',element
         call errquit(
     .     'nwpw_libfile: missing library file ',0, DISK_ERR)
      endif
      return
      end
