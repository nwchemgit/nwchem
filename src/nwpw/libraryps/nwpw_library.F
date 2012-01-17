
c
c $Id$
c
	subroutine nwpw_libgeninp(iopt,element,filename)
	implicit none
        integer iopt
	character*(*) element ! [in]
	character*(*) filename

#include "errquit.fh"
#include "inp.fh"
#include "stdio.fh"


c
	character*255 elem_fname
	integer lef,unitrc,library_type
	logical isopen
c
	call nwpw_libfile(iopt,elem_fname,element,library_type)
	lef=inp_strlen(elem_fname)
c
c       open the library file
c
	unitrc = 68
 1	continue
	unitrc = unitrc + 1
	if (unitrc.gt.80) then
	   write(luout,*)' could not find a free unit 69->80'
	   call errquit('bas_input: no free units for library',911,
     &       BASIS_ERR)
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
	call nwpw_libparse(library_type,unitrc,unitrc+1)
	close(unitrc)
	close(unitrc+1)

        call util_file_copy(elem_fname,filename)

	return
 98	continue
	write(luout,*) 'problems opening library file ',elem_fname
	call errquit('nwpw_libgeninp: nwpw failed ',0, DISK_ERR)
 99	continue
	write(luout,*) 'problems opening psp input file ',filename
	call errquit('nwpw failed ',0, DISK_ERR)
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

	integer function nwpw_libgpsp_type(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libgpsp_type= psp_type
	return
	end

	integer function nwpw_libgefg_type(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libgefg_type= efg_type
	return
	end

	double precision function nwpw_libgrlocal(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libgrlocal = rlocal_read
	return
	end

	double precision function nwpw_libgcutoff(atom)
	implicit none
#include "nwpw_inp.fh"
	character*(*) atom
	nwpw_libgcutoff = cutoff_read
	return
	end

	subroutine nwpw_libparse(library_type,u_library,u_psp_inp)
	implicit none
#include "errquit.fh"
#include "inp.fh"
#include "nwpw_inp.fh"
        integer library_type
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
        psp_type=0
        efg_type=0
        cutoff_read = 30.0d0

        if      (library_type.eq.2) then
          psp_type = 2
        else if (library_type.eq.3) then
          psp_type = 3
        end if

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
	   elseif(keyread(1:k_lgth).eq.'cutoff') then
	      read(line(sqptr_r+1:l_lgth),*) cutoff_read
	   elseif(keyread(1:k_lgth).eq.'psp_type') then
	      read(line(sqptr_r+1:l_lgth),*) psp_type
	   elseif(keyread(1:k_lgth).eq.'efg_type') then
	      read(line(sqptr_r+1:l_lgth),*) efg_type
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
 2002	call errquit('read error in nwpwlibparse',0, DISK_ERR)
	end
