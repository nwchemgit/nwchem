*
* order of precidence for choosing name
* 1) value of NWCHEM_NWPW_LIBRARY environment variable
* 2) value of the compiled in library name
*      
c the string library1 is appended by default for now
c
C $Id: nwpw_libfile.F,v 1.2 2002-01-09 01:19:42 edo Exp $
      subroutine nwpw_libfile(filename,element)
      implicit none
#include "inp.fh"
#include "stdio.fh"
      character*(*) filename
      character*(*) element
c
      character*255 libname
      integer lgth
      logical lexist,from_environment
c
      call util_getenv('NWCHEM_NWPW_LIBRARY',libname)
      lgth=inp_strlen(libname)
      from_environment = .false.
      if (lgth.gt.0) then
         lexist = .false.
         inquire(file=libname,exist=lexist)
         if (lexist)then
            from_environment = .true.
            goto 99
         else
            write(luout,*)' warning:::::::::::::: from_environment'
            write(luout,*)' NWCHEM_NWPW_LIBRARY set to: <',
     &           libname(1:inp_strlen(libname)),'>'
            write(luout,*)' but file does not exist !'
            write(luout,*)' using compiled library'
         endif
      endif
      libname = 
     &NWPW_LIBRARY
         lexist = .false.
         inquire(file=libname(1:inp_strlen(libname)),exist=lexist)
         if (.not.lexist)then
        write(luout,*)' warning:::::::::::::: from_compile'
        write(luout,*)' NWPW_BASIS_LIBRARY is: <',
     &      libname(1:inp_strlen(libname)),'>'
        write(luout,*)' but file does ','not exist or you ',
     &      'do not have ','access to it !'
         endif
 99   continue
      if (from_environment) then
          write(luout,*)' library name resolved from: environment'
       else
           write(luout,*)
     &      ' library name resolved from: compiled reference'
       endif
      write(luout,*)' library file name is: <',
     &     libname(1:inp_strlen(libname)),'>'
      write(luout,*) ' '

      filename=libname(1:inp_strlen(libname))//'/library1/'//element
c
c     check if file exists
c
      inquire(file=filename,exist=lexist)
      if(.not.lexist) then
         write(luout,*) ' could not find library for ',element
         call errquit(
     .     'nwpw_libfile: missing library file ',0)
      endif
      return
      end
