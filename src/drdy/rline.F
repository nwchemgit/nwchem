C $Id$
c ***************************************************************************
c     rline
c ***************************************************************************
c
      subroutine drdy_rline(string,istrt,isect,iend)
c
c      This subroutine finds the first non-comment and non-blank line 
c      and the location of the first character on that line.  The contents
c      of the line are placed in the variable string.  It will 
c      then change it to all lower case using the routine 'case'.
c      A flag is set if the line is the start of a new section.
c
c      istrt:   location in 80 character string that is the first non-blank
c      isect:   a 1 is returned if the first non-blank character is a * 
c               signalling the start of a section
c      iend:    a 1 is returned if the end of file 5 is encountered
c      string:  80 character string that is the contents of the current line
c
      implicit none ! double precision (a-h,o-z)
c
#include "drdyP.fh"
#include "inp.fh"
c
*::passed:      
      character * 80 string
      integer istrt
      integer isect
      integer iend
*::local:
      integer iblnk
      integer i, j
c
1000  format(A80)
c
c       initialize variables
c
      isect = 0
      iend = 0
      iblnk = 1
c
c     read in next line in file 5 - find first non-blank character and
c     determine if it is a comment (# sign signals a comment).  Once a
c     comment is found the rest of the line is skipped (by setting the
c     counter to the last character [i=80])
c
      do while (iblnk.eq.1)
        read(lu_drdy_rline,1000,end=9999) string
        i = 1
        do while (iblnk.eq.1.and.i.le.80)
          if (string(i:i).ne.' ') then
            if (string(i:i).ne.'#') then
              iblnk = 0
            else
              i = 80
            end if
          end if
          i = i+1
        end do
      end do
c
c      set location of first non-blank character
      istrt = i-1
c
c      check if it is a section header (first character a *)
c
      if (string(istrt:istrt).eq.'*') isect=1
c
c      strip of any trailing comments on the line
c
      do  j=istrt,80
        if (string(j:j).eq.'#') string(j:80) = ' '
      enddo
c
c      change to lower case
c
*       string = case(string)
      call inp_lcase(string)
*      write(fu6,*)' rline string:<',string,'>'
      return
9999  iend=1
      return
      end
