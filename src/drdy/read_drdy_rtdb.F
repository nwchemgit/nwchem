C $Id$
c
c ***************************************************************************
c     read5 -> read_drdy_rtdb
c ***************************************************************************
      subroutine read_drdy_rtdb(rtdb)
      implicit none 
#include "errquit.fh"
* set up local input file from master nwchem input on rtdb.
#include "drdyP.fh"
#include "stdio.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
*::passed:
      Integer rtdb ! [input] RTDB handle
*::local:
      character*255 filename
*
      character * 80 string
      integer istrt, isect, iend
      integer i, j
      character*255 rtdbname
      character*255 drdystring
      integer len_drdy
      integer count
c
c read and write file
c
      lu_drdy_rline = 55
      call util_file_name('drdy_input',.false.,.false.,filename)
      open(unit=lu_drdy_rline,file=filename,status='unknown')
      close(unit=lu_drdy_rline,status='delete')
      call util_file_unlink(filename)
      open(unit=lu_drdy_rline,file=filename,
     &      form='formatted',status='new', err=99999)
c
      count = 0
      rtdbname = 'drdy:nwchem:line:count'
      if (.not.rtdb_get(rtdb,rtdbname,mt_int,1,count)) then
        write(fu6,*)' no input for drdy '
        call errquit('read_drdy_rtdb: fatal error ',911, RTDB_ERR)
      endif
c
      do i = 1,count
        write(rtdbname,10000)i
        if (.not.rtdb_cget(rtdb,rtdbname,1,drdystring)) then
          write(fu6,*)' could not rtdb_cget ',
     &          rtdbname(1:(inp_strlen(rtdbname)))
          call errquit('read_drdy_rtdb: fatal error ',911, RTDB_ERR)
        endif
        len_drdy = inp_strlen(drdystring)
        write(lu_drdy_rline,'(a)') drdystring(1:len_drdy)
      enddo
      rewind lu_drdy_rline
c
c     start off by finding the first section given in the input deck
c
      call drdy_rline(string,istrt,isect,iend)
c
c     error check - check that the file wasn't empty and that the 
c     first line found was the start of a section
c
      if (iend.eq.1) then
         write(fu6,1000)
1000     format(1x,'error:  input file was empty')
         call errquit('read_drdy_rtdb: fatal error',911, RTDB_ERR)
      end if
c
      if (isect.ne.1) then
         write(fu6,1100)
1100     format(1x,'error:  first non-comment, non-blank line must be'
     *             ,' a section header that begins with a *')
         call errquit('read_drdy_rtdb: fatal error',911, RTDB_ERR)
      end if
c
c     Call the correct routine to read in the data in this section.  
c     control will come back here once another section header has
c     been found.  This will continue until an end-of-file mark
c     has been encoutered.
c
      do while (iend.eq.0)
         j = istrt + 1
         if (string(j:j+5).eq.'genera') then
               call drdy_rgener(string,iend,istrt)
         else if (string(j:j+5).eq.'react1') then
               irepr(1) = 1
               call drdy_rstat(string,iend,istrt,1)
         else if (string(j:j+5).eq.'react2') then
               irepr(2) = 1
               call drdy_rstat(string,iend,istrt,2)
         else if (string(j:j+4).eq.'prod1') then
               irepr(3) = 1
               call drdy_rstat(string,iend,istrt,3)
         else if (string(j:j+4).eq.'prod2') then
               irepr(4) = 1
               call drdy_rstat(string,iend,istrt,4)
         else if (string(j:j+4).eq.'start') then
               irepr(5) = 1
               call drdy_rstat(string,iend,istrt,5)
         else if (string(j:j+3).eq.'path') then
               call drdy_rpath(string,iend,istrt)
         else
               write(fu6,1200) string(j-1:80)
               call errquit('read_drdy_rtdb: fatal error',911, RTDB_ERR)
         end if
      end do
c
1200  format(3x,'Error:  the following is not a valid section name',
     *         /A80)
c
      return
99999 continue
      write(fu6,*)' read_drdy_rtdb: could not open: ',filename
      call errquit('read_drdy_rtdb: fatal error',911, RTDB_ERR)
10000 format('drdy:',i10)
      end

