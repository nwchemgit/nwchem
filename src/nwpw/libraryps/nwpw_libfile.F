*
* $Id: nwpw_libfile.F,v 1.7 2003-10-17 22:56:19 carlfahl Exp $
*

*
* order of precidence for choosing name
* 1) value of NWCHEM_NWPW_LIBRARY environment variable
* 2) value of the compiled in library name
*      
c the string library1 is appended by default for pspw
c the string library2 is appended by default for paw
c
      subroutine nwpw_libfile(flag,filename,element)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"
      integer       flag
      character*(*) filename !RESULT
      character*(*) element
c
      character*255 libname
      external util_find_dir
      logical  util_find_dir,lexist
      integer lgth,l
      logical from_environment
c
      call util_getenv('NWCHEM_NWPW_LIBRARY',libname)
      lgth=inp_strlen(libname)
      from_environment = .false.
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
         endif
      endif
      libname = 
     &NWPW_LIBRARY
      if (util_find_dir(libname)) then
         goto 99
      else
         write(luout,*)' warning:::::::::::::: from_compile'
         write(luout,*)' NWPW_BASIS_LIBRARY is: <',
     &        libname(1:inp_strlen(libname)),'>'
         write(luout,*)' but file does ','not exist or you ',
     &        'do not have ','access to it !'
      endif
      call errquit('nwpwlibfile: no nwpw library found',0, DISK_ERR)
 99   continue
      if (from_environment) then
          write(luout,*)' library name resolved from: environment'
       else
           write(luout,*)
     &      ' library name resolved from: compiled reference'
       endif


*     **** library for paw ****
      if (flag.eq.2) then
        write(luout,*)' library file name is: <',
     &     libname(1:inp_strlen(libname)),
     &     'library2/', element, '_basis','>'
        write(luout,*) ' '

        l = index(element,' ')  -1
        if (l.lt.1) l = 2
        filename=
     > libname(1:inp_strlen(libname))
     > //'/library2/'//element(1:l)//'_basis'

*     **** library for pspw and band ****
      else
        write(luout,*)' library file name is: <',
     &     libname(1:inp_strlen(libname)),
     &     'library1/', element, '>'
        write(luout,*) ' '

        filename=libname(1:inp_strlen(libname))//'/library1/'//element
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
