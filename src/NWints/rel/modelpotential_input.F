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
      character*255 mpfile
      character*32 pname
      character*2 atsym
      integer i,j
      integer unitno
      integer necp,igau,ngau
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
c
      logical util_io_unit
      external util_io_unit
c
      pname = "modelpotential_input"
c      
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
      if(.not.util_io_unit(80,90,unitno))
     +  call errquit(pname//"util_io_unit failed",0,0)
c
c     get model potential file
      if (.not. rtdb_cget(rtdb,'modelpotential:file', 1,mpfile))
     $    call errquit(pname//'rtdb_cget failed', 0, RTDB_ERR)
c
c     print model potential information
      write(luout,*)
      write(luout,1000)
      write(luout,*) "Reading model potential file:",
     &            mpfile(1:inp_strlen(mpfile))
      write(luout,*)
c
c     check if model potential file is present
      fileexists = .false.
      inquire(file=mpfile,exist=fileexists)
      if (.not.fileexists) 
     $  call errquit(pname//'modelpotential file not found',0,INPUT_ERR)
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
c
 1000 format(' Using Relativistic Model Potential Approach'/
     & ' C. van Wullen, J. Chem. Phys. 109, 392 (1998)'/
     & ' C. van Wullen, C. Michauk, J. Chem. Phys. 123, 204113 (2005)'/)
c
      return
      end
