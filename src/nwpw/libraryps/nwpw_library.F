C $Id: nwpw_library.F,v 1.1 2001-10-31 02:07:52 edo Exp $
	subroutine nwpw_libgeninp(element,filename)
	implicit none
#include "inp.fh"
#include "stdio.fh"
	character*(*) element ! [in]
	character*(*) filename
c
	character*255 elem_fname
	integer lef,unitrc
	logical isopen
c
	call nwpw_libfile(elem_fname,element)
	lef=inp_strlen(elem_fname)
	write(luout,*) ' got this library ',elem_fname(1:lef)
c
c       open the library file
c
	unitrc = 68
 1	continue
	unitrc = unitrc + 1
	if (unitrc.gt.80) then
	   write(luout,*)' could not find a free unit 69->80'
	   call errquit('bas_input: no free units for library',911)
	endif
	inquire(unit=unitrc,opened=isopen)
	if (isopen) goto 1
	open (unit=unitrc,  file=elem_fname,
     &        status='old', err=98,
     &        form='formatted',
     &        access='sequential')
	open (unit=unitrc+1,  file=filename,
     &        status='new', err=99,
     &        form='formatted',
     &        access='sequential')
	close(unitrc)
	close(unitrc+1)
	return
 98	continue
	write(luout,*) 'problems opening library file'
	call errquit('nwpw failed ',0)
 99	continue
	write(luout,*) 'problems opening psp input file'
	call errquit('nwpw failed ',0)
	return
	end
	integer function nwpw_libglmax(atom)
	implicit none
	character*(*) atom
	nwpw_libglmax= -1
	return
	end
	integer function nwpw_libglocp(atom)
	implicit none
	character*(*) atom
	nwpw_libglocp= -1
	return
	end
	double precision function nwpw_libgrlocal(atom)
	implicit none
	character*(*) atom
	nwpw_libgrlocal = 1.0d0
	return
	end
