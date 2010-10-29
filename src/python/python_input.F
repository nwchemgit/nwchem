      subroutine python_input(rtdb)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
      integer rtdb
c
      integer k, line, leading
      character*255 buf, key
      logical oprint
c
      oprint = .true.
c
c     Process any keys on the PYTHON directive
c
 1    if (inp_a(key)) then
         if (inp_compare(.false.,key,'print')) then
            oprint = .true.
         else if (inp_compare(.false.,key,'noprint')) then
            oprint = .false.
         endif
         goto 1
      endif
      call inp_clear_err
c
      if (oprint) then
         write(6,*)
         call util_print_centered(6,'NWChem Python program',40,.true.)
         write(6,*)
      endif
c
      open(1,file="nwchem.py",form="formatted",status="unknown")
c
c     Need to circumvent NWChem's processing of special characters
c
c     Also, python insists that the initial indent level is zero
c     so to permit NWChem python code to be indented within the
c     python block quietly strip off from all input lines
c     the amount of white space present on the first line
c
      line = 1
      leading = 0
 10   if (inp_read_physical_line(buf)) then
         if (inp_compare(.false.,buf,'end')) goto 20
         if (line .eq. 1 .and. buf.ne.' ') then
            do k = 1, len(buf)
               if (buf(k:k) .ne. ' ') goto 7
            enddo
 7          leading = k - 1
         endif
         line = line + 1
c
         do k = 1, leading
            if (buf(k:k) .ne. ' ') call errquit
     $           ('python_input: indentation must be >= that'//
     $           ' of first line',leading, INPUT_ERR)
         enddo
         write(1,15) (buf(k:k),k=leading+1,inp_strlen(buf))
         if (oprint)
     $        write(6,15) (buf(k:k),k=leading+1,inp_strlen(buf))
 15      format(255a1) 
         goto 10
      endif
c
 20   close(1,status='keep')
      if (oprint) then
         write(6,*)
         call util_flush(6)
      endif
c
      end
