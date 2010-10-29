C $Id$
c
c ***************************************************************************
c     rword
c ***************************************************************************
c
      subroutine drdy_rword(string,istrt,ierr)
c
c     This subroutine is a utility routine that takes a character string
c     and searches from the current character (istrt) for the first
c     character that occurs after at least one blank (ie for the next 
c     word on a line).  It returns the location of that character in
c     istrt.  ierr is set to 1 if no other characters were found in 
c     the string.
c
      implicit double precision (a-h,o-z)
      character * 80 string
c
c     clear flags
      jblank = 0
      ierr = 0
c
c     search for the first blank - signaling the end of the current word
c
      do while (istrt.le.79.and.jblank.eq.0)
         istrt = istrt+1
         if (string(istrt:istrt).eq.' ') then
c                                               find first letter of next word
            jblank = 1
            do while (string(istrt:istrt).eq.' '.and.ierr.eq.0)
               istrt = istrt + 1
               if (istrt.gt.80) ierr=1
            end do
         end if
         if (istrt.eq.80) ierr=1
      end do
      return
      end
