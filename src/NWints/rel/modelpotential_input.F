c
c     process model potential input
c
      subroutine modelpotential_input(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "util.fh"
#include "inp.fh"
#include "modelpotential.fh"
c
      integer rtdb ! [input] handle to database
c    
      character*255 libname
      character*32 pname
      character*2 atsym
      integer i,j
      integer unitno
      integer necp,igau,ngau
      integer mpoption
c
      integer ok
c
      character*16 tagi      ! tag of atom i
      character*2 symi
      character*16 elemi
      double precision ci(3),chgi
      integer atni
      logical status
      logical fileexists
      logical does_it_exist
c
      logical util_io_unit
      external util_io_unit
c
      logical util_find_dir
      external util_find_dir
c
      pname = "modelpotential_input"
c      
c     get model potential file name
      if(.not.rtdb_get(rtdb,'modelpotential:option',mt_int,1,mpoption)) 
     &  call errquit(pname//'unable to access mpotion',0,RTDB_ERR)
c
c     first: look for the environment variable
      does_it_exist=.false.
      call util_getenv('NWCHEM_BASIS_LIBRARY',libname)
      if (inp_strlen(libname).gt.0) then
        mpdir=libname(1:inp_strlen(libname))//"rel-modelpotentials/"
        does_it_exist = util_find_dir(mpdir) 
        if (does_it_exist) goto 10
      end if
c       
c     second: look for the src directory
      does_it_exist=.false.
      call util_nwchem_srcdir(libname)
      if (inp_strlen(libname).gt.0) then
        mpdir=libname(1:inp_strlen(libname))
     $     //"/basis/libraries/rel-modelpotentials/" 
        does_it_exist = util_find_dir(mpdir) 
        if (does_it_exist) goto 10
      end if
c
c     third: look in the local or permanent directory
      does_it_exist=.false.
      call util_directory_name(mpdir, .false., ga_nodeid())
      if (inp_strlen(mpdir).gt.0) then
        mpdir=mpdir(1:inp_strlen(mpdir))//"/"
        does_it_exist = util_find_dir(mpdir) 
        if (does_it_exist) goto 10
      end if
c
c     abort if the directory does not exist
      if (.not.does_it_exist)
     $  call errquit(pname//'modelpotential not found',0,INPUT_ERR)
c
  10  continue
c
c     build mpfile name
      if (mpoption.eq.1) then
        mpfile=mpdir(1:inp_strlen(mpdir))//"modbas.4c"
      else if (mpoption.eq.2) then
        mpfile=mpdir(1:inp_strlen(mpdir))//"modbas.2c"
      else
        call errquit(pname//'unknown option',0,RTDB_ERR)
      end if
c
c     initialize parameters
      necp = 0
      ngau = 0
      do i = 1,120
        mpelem(i) = ""
        do igau = 1,50
           mpexpo(i,igau) = 0.d0
           mpcoef(i,igau) = 0.d0
        end do
      end do
c
c     get a file unit
      if(.not.util_io_unit(80,90,unitno))
     +  call errquit(pname//"util_io_unit failed",0,0)
c
c     check if model potential file is present
      fileexists = .false.
      inquire(file=mpfile,exist=fileexists)
      if (.not.fileexists) 
     $  call errquit(pname//'modelpotential not found',0,INPUT_ERR)
c
c     open and read the file
      open(unit=unitno, file=mpfile,status='unknown')
      i=0
 100  continue
      read(unitno,*,end=200) atsym,necp,ngau
      i=i+1
      mpelem(i) = atsym
c      write(luout,*) i,mpelem(i)
      read(unitno,*,end=200)
      do igau = 1,ngau
       read(unitno,*,end=200) mpexpo(i,igau),mpcoef(i,igau)
c       write(luout,*) i,igau,mpexpo(i,igau),mpcoef(i,igau)
      end do
      goto 100
 200  continue
      close(unitno)
c
      mpmaxelem = i  ! maximum number of elements found in file
      if(.not.rtdb_put(rtdb,'modelpotential:mpmaxelem',mt_int,1,
     &   mpmaxelem))
     $    call errquit(pname//'rtdb_put failed', 0, RTDB_ERR)
c
      return
      end
