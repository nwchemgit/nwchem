C $Id: nwpw_library.F,v 1.2 2001-11-03 01:12:50 edo Exp $
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
	call nwpw_libparse(unitrc,unitrc+1)
	close(unitrc)
	close(unitrc+1)
	return
 98	continue
	write(luout,*) 'problems opening library file ',elem_fname
	call errquit('nwpw_libgeninp: nwpw failed ',0)
 99	continue
	write(luout,*) 'problems opening psp input file ',filename
	call errquit('nwpw failed ',0)
	return
	end
	integer function nwpw_libglmax(atom)
	implicit none
#include "nwpw_inp.fh"
 	character*(*) atom
	nwpw_libglmax= lmax_read
	return
	end
	integer function nwpw_libglocp(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libglocp= locp_read
	return
	end
	double precision function nwpw_libgrlocal(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libgrlocal = rlocal_read
	return
	end
	subroutine nwpw_libparse(u_library,u_psp_inp)
	implicit none
#include "inp.fh"
#include "nwpw_inp.fh"
	integer u_library,u_psp_inp
	character*1024 line
	character*12 keyread
c
	integer l_lgth,i,sqptr_l,sqptr_r,k_lgth
	integer read_int
	double precision read_float
	logical gotsq
c
c       square parenthesis [] are used for lmax locp rlocal
c
	lmax_read=-1
	locp_read=-1
	rlocal_read=1.d0
 1	continue
	read(u_library,'(A)',end=2001,err=2002) line
	l_lgth=inp_strlen(line)
c
c       scan for [
c
	gotsq=.false.
	do i=1,l_lgth
	   if(line(i:i).eq.'[') then
	      gotsq=.true.
	      sqptr_l=i
	   endif
	enddo
	if(gotsq) then
	   do i=sqptr_l,l_lgth
	      if(line(i:i).eq.']') then
		 sqptr_r=i
	      endif
	   enddo
	   keyread=line(sqptr_l+1:sqptr_r-1)
	   k_lgth=sqptr_r-sqptr_l-1
c       keyword recognition
c       lmax
	   if(keyread(1:k_lgth).eq.'lmax') then
	      read(line(sqptr_r+1:l_lgth),*) lmax_read
	   elseif(keyread(1:k_lgth).eq.'locp') then
	      read(line(sqptr_r+1:l_lgth),*) locp_read
	   elseif(keyread(1:k_lgth).eq.'rlocal') then
	      read(line(sqptr_r+1:l_lgth),*) rlocal_read
	   else
	      write(0,*) ' failed keyword recognition for',
     .  keyread(1:k_lgth)
c       should be put warning for non recognized keyword?
	   endif
	else
	   write(u_psp_inp,*) line(1:l_lgth)
	endif
	goto 1
 2001	continue
	return
 2002	call errquit('read error in nwpwlibparse',0)
	end
