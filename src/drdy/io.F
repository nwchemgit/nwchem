C $Id$
      Character*80 Function drdy_case (string)
c
c   Function which takes a string of 80 characters and converts the 
c   upper case letters in the string to lower case letters
c
      Character*80 string
C
      Integer i,itry
      Character * 1 xlett
c
      do i = 1, 80
         xlett = string(i:i)
         itry = ichar (xlett)
         if (xlett .ge. 'A' .and. xlett .le. 'Z') then 
            itry = itry + 32
            string (i:i) = char (itry)
         endif
      enddo
c
      drdy_case  = string 
c
      return
      end
      Double Precision Function drdy_cfloat(string)
C
      Implicit None
#include "errquit.fh"
C
      Character*80 string,number
      Character ch
      Logical lexp,ldec
      Integer length,i,ibeg,iend,idec,iue,ile,iud,ild,isum,iexp,n,j
      Double Precision sign,f1,f2,f3,esign
      Integer fu6
      Data fu6/6/
C
      lexp = .false.
      ldec = .false.
      length = len(string)
      if (length .eq. 0) then
         drdy_cfloat = 0.0d0
         return
      endif
c     write(fu6,*) length,string
C
C     Find the first nonblank character
C
      i = 1
      do while (string(i:i) .eq. ' ' .and. i .le. length)
         i = i + 1
      enddo
C
C     If it is a blank string set Function to zero
C
      if (i .gt. length) then
         drdy_cfloat = 0.0d0
         return
      endif
      ibeg = i
c
c     find the first blank character after the number
c
      i = ibeg+1
      do while (string(i:i) .ne. ' ' .and. i .le. length)
         i = i + 1
      enddo
      iend = i-1
C
C     Stripe the blanks before and after the number
C
      number = string(ibeg:iend)
      length = iend - ibeg + 1
C   
C     Make sure there is no blank left
C
      if (index(number,' ') .le. length) then
         write(fu6,6000) string
         call errquit ('drdy_cfloat:1: fatal error',911, UNKNOWN_ERR)
      endif
C
C     Find the decimal point
C
      idec = index(number,'.')
      if (idec .ne. 0) ldec = .true.
C
C     Find the exponential symbol
C
      iue = index(number,'E')
      ile = index(number,'e')
      iud = index(number,'D')
      ild = index(number,'d')
      isum = iue + ile + iud + ild
      iexp = max0(iue,ile,iud,ild)
      if (isum .gt. iexp) then
         write(fu6,6000) string
         call errquit('drdy_cfloat:2: fatal error',911, UNKNOWN_ERR)
      endif
      if (iexp .ne. 0) then
         lexp = .true.
      else
         iexp = length + 1
      endif
c
      if (.not. ldec) idec = iexp
C
C     Get the number before decimal
C
      ibeg = 2
      if (number(1:1) .eq. '+') then
         sign = 1.0d0
      elseif(number(1:1) .eq. '-') then
         sign = -1.0d0
      else
         sign = 1.0d0
         ibeg = 1
      endif
      if (ibeg .eq. iexp) then
         f1 = 1.0d0
      else
         f1 = 0.0d0
      endif
      do i = ibeg,idec-1
         ch = number(i:i)
         if (ch .ge. '0' .and. ch .le. '9') then
            n = ichar(ch) - ichar('0')
            f1 = f1 * 10.0d0 + dble(n)
         else
            write(fu6,6000) string
            call errquit('drdy_cfloat:3: fatal error',911,
     &       UNKNOWN_ERR)
         endif
      enddo
C
C     Get the number after decimal 
C
      f2 = 0.0d0
      if (ldec) then
         j = 0
         do i = idec+1,iexp-1
            ch = number(i:i)
            if (ch .ge. '0' .and. ch .le. '9') then
               n = ichar(ch) - ichar('0')
               f2 = f2 * 10.0d0 + dble(n)
               j = j + 1
            else
               write(fu6,6000) string
               call errquit('drdy_cfloat:4: fatal error',911,
     &       UNKNOWN_ERR)
            endif
         enddo
         f2 = f2 / 10.0d0 ** dble(j)
      endif
C
C    Get the exponent
C
      esign = 1.0d0
      f3 = 0.0d0
      if (lexp) then 
         ibeg = iexp + 2
         if (number(iexp+1:iexp+1) .eq. '+') then
            esign = 1.0d0
         elseif(number(iexp+1:iexp+1) .eq. '-') then
            esign = -1.0d0
         else
            esign = 1.0d0
            ibeg = iexp + 1
         endif
         do i = ibeg,length
            ch = number(i:i)
            if (ch .ge. '0' .and. ch .le. '9') then
               n = ichar(ch) - ichar('0')
               f3 = f3 * 10.0d0 + dble(n)
            else
               write(fu6,6000) string
               call errquit('drdy_cfloat:5: fatal error',911,
     &       UNKNOWN_ERR)
            endif
         enddo
      endif 
C
      drdy_cfloat = (sign * (f1 + f2)) * 10.0d0 ** (esign*f3)
c
      return
C
6000  Format(/1X,'Illegal number: ',A80)
C
      end
      Subroutine clen (line,ibeg,nlen,ndim)
C Get length of line and location of first nonblank character
      Character*80 line
      Integer ibeg,nlen,ndim
C
      Integer n
C
      n = ndim
      do while (n.gt.0 .and. line(n:n).eq.' ')
         n = n - 1
      enddo
      nlen = n
      ibeg = 1
      do while (ibeg.lt.nlen .and. line(ibeg:ibeg).eq.' ')
         ibeg = ibeg + 1
      enddo
      return
      end
      Integer Function drdy_icint(svalue)
c
c   Function which takes a character string which is really a integer and
c   converts it to an integer.  Error checking not possible using ichar is
c   done.  Program will stop if an illegal character is passed via the 
c   string 'svalue'.
c
      Implicit None
#include "errquit.fh"
C
      Character * 80 svalue
c
      Integer istart,isign,j,k
      Integer fu6
      Data fu6/6/
C
c     check for negative or positive sign (- or +)
c
      drdy_icint = 0
      istart = 1
      isign = +1
      if (svalue(1:1).eq.'+') then
         isign = +1
         istart = 2
      else if(svalue(1:1).eq.'-') then
         isign = -1
         istart = 2
      end if
c
c     now convert each character to an integer - and store in drdy_icint
c
      k = istart
      do while (svalue(k:k).ne.' ')
         j = ichar(svalue(k:k)) - ichar('0')
         if (j.lt.0.or.j.gt.9) then
            write(fu6,6000)svalue
             call errquit('drdy_icint: fatal error',911,
     &       UNKNOWN_ERR)
         end if
c
         drdy_icint = 10.d00*drdy_icint + float(j)
         k = k+1
      end do
      drdy_icint = drdy_icint * isign
c
      return
6000  Format(3x,'ERROR: illegal digit in an integer: ',A80)
      end
      Subroutine drdy_rline(string,istrt,isect,iend)
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
      Implicit none
      Character*80 string,drdy_case
      External drdy_case
      Integer istrt,isect,iend
#include "drdyP.fh"
C
      Integer i,iblnk,j
C
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
         read(lu_drdy_rline,6000,end=9999,err=9999) string
c         read(fu5,6000,end=9999,err=9999) string
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
      string = drdy_case(string)
      return
9999  iend=1
      return
6000  format(A80)
      end
      Subroutine drdy_rword(string,istrt,ierr)
c
c     This subroutine is a utility routine that takes a character string
c     and searches from the current character (istrt) for the first
c     character that occurs after at least one blank (ie for the next 
c     word on a line).  It returns the location of that character in
c     istrt.  ierr is set to 1 if no other characters were found in 
c     the string.
c
      Implicit None
      Character * 80 string
      Integer istrt,ierr
c
      Integer jblank
C
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
