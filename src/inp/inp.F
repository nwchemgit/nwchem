      block data inp_data
C$Id$
      implicit none
#include "inpP.fh"
c
      data iread   /5/
      data iwrite  /6/
      data jrec    /-1/
      data jump    /0/
      data oswit   /.false. /
      data nerr    /999/
      data nline   /0/
      data noline  /0/
      data ierrpos /-1/
      data errmsg  /' '/
      data input_line /0/
      data xblnk /' '/
#if defined(CRAY) || defined(HPUX) || defined(WIN32) || defined(PSCALE) || defined(__FLANG) || ( __GNUC__ >= 4) || defined(FUJITSU)
      data xtab  /'	'/      ! Tab ... no backslash necessary
#elif (defined(LINUX) || defined(MACX)) && !defined(PGLINUX) && !defined(XLFLINUX) &&!( __GNUC__ >= 4)
      data xtab  /9/            ! Tab ... g77 has trouble with escape sequence
#else
      data xtab  /'\	'/      ! Tab ... note backslash for cpp
#endif
      data xsplit/';'/
      data xcomm /'#'/
      data xback /'\\'/         ! Backslash ... note backslash for cpp
      data xquote/'"'/
c
      data save_level /0/            ! Depth of nesting in save/restore
      data include_level /0/
c
      end
C> \ingroup inp
C> @{
C>
C> \brief Initializes free format input routines
C>
C> This function initializes free format input routines to take input
C> from Fortran unit `ir` and send output to fortran unit `iw`. The
C> input file is processed from the current location.
C>
C> Function `inp_init()` should be invoked each time the input file is
C> repositioned using other than `inp_*()` routines (e.g. by using
C> rewind). It should also be called when the focus of the input
C> processing transfers to another file. This happens for example
C> when the main input file includes another file.
c
      subroutine inp_init(ir,iw)
      implicit none
#include "inpP.fh"
      integer ir !< [Input] The standard input unit number
      integer iw !< [Input] The standard output unit number
      external inp_data  ! For T3D linker
c     
      iread = ir
      iwrite = iw
      jrec = -1
      jump = 0
      oswit = .false. 
      nerr = 999
      nline = 0
      noline = 0
      ierrpos = -1
      errmsg = ' '
      input_line = 0
c     
      end
C
C> \brief Returns the number of fields in a line of input
C>
C> For a given line of input this function returns the number of
C> white-space separated fields on that line. If the number of fields
C> equals 0 then the end of the input file has been reached.
C>
C> \return Return the number of fields in a line of input
c
      integer function inp_n_field()
      implicit none
#include "inpP.fh"
c     
c     return no. of fields in the input line ... 0 = EOF
c     
      inp_n_field = jump
c     
      end
C
C> \brief Returns the number of fields processed
C>
C> This function returns the number of fields of the current input line
C> processed so far. The return value is in the range from 0 to
C> inp_n_field().
C>
C> \return Return the number of fields processed
c
      integer function inp_cur_field()
      implicit none
#include "inpP.fh"
c     
c     return no. of fields processed so far (0,...,inp_n_field())
c     
      inp_cur_field = jrec
c     
      end
C
C> \brief Set the number of the field to be read next
C>
C> This subroutine moves the read position on the current input line
C> so that the field with the specified number is read next. Valid
C> arguments are 1 to inp_n_field().
c
      subroutine inp_set_field(ivalue)
      implicit none
#include "errquit.fh"
#include "inpP.fh"
      integer ivalue !< [Input] The number of the next field
      integer inp_n_field
      external inp_n_field
c     
c     set field to be read next (ivalue=0,...,inp_n_field())
c     
      if (ivalue.lt.0 .or. ivalue.gt.inp_n_field())
     $     call errquit('inp_set_field: stupid field value',ivalue,
     &       INPUT_ERR)
      jrec = ivalue
c
      ierrpos = -1
      errmsg = ' '
c     
      end
C
C> \brief Retrieve as much of the current input line as the buffer can
C> hold
C>
C> Simply copy the current input line to the buffer. The Fortran 
C> standard specifies that as much will be copied as will fit. No 
C> tests are performed to check whether the data will fit. Also, whether
C> the data fits does not affect the return value.
C>
C> \return Return .true. if there is a current input line, and .false.
C> otherwise
c
      logical function inp_line(z)
      implicit none
#include "inpP.fh"
      character*(*) z !< [Output] (part of) the current input line
c     
c     set the variable z to be as much of the current input line
c     that it can hold
c     
      if (jump .gt. 0) then
         z = ia
         inp_line = .true.
         ierrpos = -1
         errmsg = ' '
      else
         errmsg = 'no input line available'
         ierrpos = -1
         inp_line = .false.
      endif
c     
      end
C
C> \brief Read a physical line of input straight into buffer
C>
C> This routine is intended for programms embedded into NWChem
C> that need to read from the NWChem input file bypassing any
C> processing that NWChem does on input data, but while still
C> maintaining a correct count of the input line number for
C> error reporting.
C>
C> \return Return .true. if there is an input line to be read, and
C> .false. otherwise.

      logical function inp_read_physical_line(buf)
      implicit none
#include "inpP.fh"
      character*(*) buf
c
c     Read a physical line of input into buf() WITHOUT any
c     tokenizing, scanning, etc. 
c
c     This routine is intended for programms embedded into NWChem
c     that need to read from the NWChem input file bypassing any
c     processing that NWChem does on input data, but while still
c     maintaining a correct count of the input line number for
c     error reporting.
c
c     First void out any info about the current input line
c
      ierrpos = -1
      jump = 0
      jrec = 0
c
      buf = ' '
      read(5,'(a)',end=10,err=10) buf
      input_line = input_line + 1
c
      inp_read_physical_line = .true.
      return
c
 10   inp_read_physical_line = .false.
      errmsg  = 'unexpected end of data file'
c
      end
C
C> \brief Read and parse a line of input
C>
C> This function reads a line of input and tokenizes it into fields.
C> The data this produces is stored internally, but it can subsequently
C> be queried using other routines in this module. This routine also 
C> deals with special cases such as include files and comment lines.
C>
C> The lines from the input are split into white space (blank or
C> tab) separated fields. White space may be incorporated into a field
C> by enclosing it in quotes (for example, "new name"). The case of
C> the input is preserved.
C>
C> Blank lines are ignored and text beginning with a pound or hash
C> symbol (`#`) is treated as a comment. A backslash (`\`) at the end
C> of a line (followed only by white space) can be used to concatentate
C> physical input lines into one logical input line.  A semicolon
C> (`;`) may be used to split a single physical input line into
C> multiple logical input lines. The special command characters hash
C> (`#`), semicolon (`;`) and quotation mark (`"`) will be treated
C> simply as characters only if prefaced by a backslash. (NOTE: This
C> must be done even when the character appears within a character
C> string enclosed in quotes.)
C>
C>The number of fields read is initially set to 0, there being a total
C> of `inp_n_field()` fields in the line.
C>
C> If a non-blank line is successfully parsed then `.true.` is returned;
C> otherwise an internal error message is set and `.false.` is returned.
C> Possible errors include such actions as detection of unexpected EOF
C> (which can be checked for with `function inp_eof()`), or failure to
C> parse the line (e.g., a character string without a terminating
C> quote).
C>
C> End of file (EOF) is usually indicated by reaching the actual
C> end of the physical input file.  Alternatively, the user can specify
C> the end of file location at any point by inserting a
C> physical input line that begins with an asterisk (*), or a period,
C> or the letters EOF (which may be in upper or lower case), and is
C> followed only by trailing white space.
C>
C> The maximum input line width is 1024 characters.
C>
C> \return Return .true. if another line of input was read successfully,
C> and .false. otherwise.
C
      logical function inp_read()
      implicit none
#include "errquit.fh"
#include "inpP.fh"
c     
c     this routine reads a data card and scans it for non - space fields
c     the number of fields is stored in jump, the starting point of a
c     field in istrt(i) and the number of characters in that field
c     in inumb(i).
c     
      character*1 xprev
      integer lenja, i, k, mark, j, nbegin, nfini
      integer inp_strlen, jwidth
      logical inp_compare
      external inp_strlen, inp_compare
      integer ncol(max_field)
      logical ois_ws            ! Inline funtion  
      character*1 xtest
c
      ois_ws(xtest) = (xtest.eq.xblnk .or. xtest.eq.xtab)
c
      inp_read = .false.  ! First assume things will not work
c
c     
 1    nline=nline+1
      if(nline.le.noline)go to 150
c
      if (oswit) then
         ierrpos = -1
         errmsg = 'unexpected end of data file'
         inp_read = .false.
         jump=0
         jrec=0
         return
      else
         ierrpos = -1
         errmsg = ' '
      endif
c
c     read next physical input line
c
 101  lenja = 0
 100  read(iread,'(a)',end=300)ja(lenja+1:max_width)
      input_line = input_line + 1
      lenja = inp_strlen(ja)
c
c     Check for . * eof at beginning of line to indicate EOF
c
      if (lenja.eq.1 .and. (ja(1:1).eq.'.' .or. ja(1:1).eq.'*'))
     $     goto 300
      if (lenja.eq.3 .and. inp_compare(.false., 'eof', ja(1:3))) 
     $     goto 300
c
c     Handle include statement
c
      if (inp_compare(.false.,ja(1:7),'include')) then
         if (include_level .eq. max_include_level) call errquit
     $        ('inp_read: include nested too deep ', include_level,
     &       INPUT_ERR)
         include_level = include_level + 1
         include_file_name(include_level) = ja(9:)
         write(6,*) ' include: start of ', 
     $        include_file_name(include_level)
     $        (1:inp_strlen(include_file_name(include_level)))
         open(80+include_level,file=include_file_name(include_level),
     $        form='formatted', status='old', err=105)
         call inp_save_state()
         call inp_init(80+include_level,6)
         goto 1
      endif
c
c     handle blank lines and concatenation using backslash
c
      if (lenja.eq.0) then
         goto 100
      else
         if (ja(lenja:lenja) .eq. xback) then
            ja(lenja:lenja) = xblnk
            goto 100
         endif
      endif
      jwidth = inp_strlen(ja)
c
c     handle comments from # to eol ... allow for backslash quoting
c
      xprev = xblnk
      do i=1, jwidth
 91      if (ja(i:i) .eq. xcomm .and. xprev.ne.xback) then
            lenja = inp_strlen(ja)
*            write(iwrite,90) ja(i+1:lenja)
* 90         format(/' comment :-',1x,a)
            ja(i:max_width) = xblnk
            goto 80
         else if (ja(i:i) .eq. xcomm .and. xprev.eq.xback) then
c     Shuffle string down to overwrite quoting backslash
            tmp = ja
            ja(i-1:jwidth) = tmp(i:jwidth)
            xprev = xblnk
            goto 91
         endif
         xprev = ja(i:i)
      enddo
c
 80   if (inp_strlen(ja(1:i)) .eq. 0) goto 101 ! All line comments
c
c     figure out where ; splits physical line into multiple logical lines
c     again handling quoted backslash
c
      k=jwidth
      mark=0
      xprev = xblnk
      do i=1,jwidth
 81      if(ja(i:i).eq.xsplit .and. xprev.ne.xback) then
            mark=mark+1
            ncol(mark)=i
         else if(ja(i:i).eq.xsplit .and. xprev.eq.xback) then
c     Shuffle string down to overwrite quoting backslash
            tmp = ja
            ja(i-1:jwidth) = tmp(i:jwidth)
            xprev = xblnk
            goto 81
         endif
         xprev = ja(i:i)
      enddo
      noline=1
      if(mark.eq.0) then
         nstart(noline)=1
         nend(noline)=jwidth
      else
         i=ncol(mark)+1
         if(i.le.jwidth) then
            do j=i,jwidth
               if(.not. ois_ws(ja(j:j))) go to 170
            enddo
         endif
         k=ncol(mark)-1
         mark=mark-1
c     
 170     noline=mark+1
         nstart(1)=1
         do i=1,mark
            j=ncol(i)
            nend(i)=j-1
            nstart(i+1)=j+1
         enddo
         nend(noline)=k
      endif
      nline=1
c     
c     Start processing next logical input line (put into ia(1:iwidth))
c
 150  jump=0
      jrec=0
      nbegin = nstart(nline)
      nfini  = nend(nline)
      iwidth = nfini-nbegin+1
      ia = xblnk
      ia(1:iwidth)=ja(nbegin:nfini)
c     
c     partition input line into strings inside double quotes or
c     white space separated fields
c
      i = 1
 151  continue
      do j = i, iwidth
         if (.not. ois_ws(ia(j:j))) goto 152
      enddo
      goto 155                  ! Done
 152  i = j
c
      jump = jump + 1
      istrt(jump) = i
      if (ia(i:i) .eq. xquote) then
c
c     Quoted string ... look for closing quote
c
         do j = i+1, iwidth
 154        if (ia(j:j).eq.xquote .and. ia(j-1:j-1).ne.xback) then
               goto 153
            else if (ia(j:j).eq.xquote .and. ia(j-1:j-1).eq.xback) then
               tmp = ia
               ia(j-1:max_width) = tmp(j:max_width)
               goto 154
            endif
         enddo
         ierrpos = j
         errmsg = 'no terminating quote for string'
         inp_read = .false.
         oswit = .false.
         return
 153     continue
      else
c
c     Simple field ... look for next ws
c
         do j = i+1, iwidth
            if (ois_ws(ia(j:j))) goto 157
         enddo
 157     j = j - 1
      endif
c
      inumb(jump) = j - istrt(jump) + 1
      i = j + 1
      goto 151
c
 155  continue                  ! Finished tokenizing
      if (jump .gt. 0) iwidth = istrt(jump)+inumb(jump)-1
      inp_read = .true.
      return
c
 300  if (include_level .gt. 0) then ! End of file detected
         close(80+include_level)
         write(6,*) ' include: end of ', 
     $        include_file_name(include_level)
     $        (1:inp_strlen(include_file_name(include_level)))
         include_level = include_level - 1
         call inp_restore_state()
         goto 1
      else
         oswit = .true.
         ierrpos = -1
         errmsg = 'unexpected end of data file'
         inp_read = .false.
         jump=0
         jrec=0
      endif
      return
c
 105  call errquit('inp_read: failed to open include file',0, INPUT_ERR)
c
      end
C
C> \brief Return .true. if the end-of-file was reached
C>
C> \return Return .true. if the end-of-file was reached, and .false. 
C> otherwise
c
      logical function inp_eof()
      implicit none
#include "inpP.fh"
c      
      inp_eof = oswit
c
      end
c
C> \brief Clear all error conditions and error messages
c
      subroutine inp_clear_err()
      implicit none
#include "inpP.fh"
c
c     Clear error conditions and messages
c
      ierrpos = -1
      errmsg = ' '
c
      end
c
C> \brief Print the current input line and read position if an error
C> has occurred.
c
      subroutine inp_errout
      implicit none
#include "inpP.fh"
      integer length, inp_strlen, i
      external inp_strlen
      character*1 xpt, xstp
      data xpt,xstp/'*', '.'/
c     
c     If an error has occured print out the error message
c     and the position in the current input line
c     
      if (include_level .gt. 0) then
         write(6,*) ' Include file stack '
         do i = 1, include_level
            write(6,321) i, include_file_name(i)
     $        (1:inp_strlen(include_file_name(i)))
 321        format(1x,i5,2x,a)
         enddo
      endif
      if (errmsg .ne. ' ') then
         length = inp_strlen(errmsg)
         write(iwrite, 40) input_line, errmsg(1:length)
 40      format(' input error at line', i5,': ', a)
c
****  jrec=-1 ! Why was this being set?  Seems undesirable
         write(iwrite,50)ia(1:iwidth)
 50      format(1x,a)
         if (ierrpos .gt. 0) then
            do 60 i=1,iwidth
               tmp(i:i)=xstp
 60         continue
            tmp(ierrpos:ierrpos)=xpt
            write(iwrite,50)tmp(1:iwidth)
         endif
      endif
c     
      end
c
C> \brief Write the current input line to standard output
c
      subroutine inp_outrec
      implicit none
#include "inpP.fh"
c
c     Write out the current input line
c
      write(iwrite,50) input_line, ia(1:iwidth)
50    format(1x,i5,': ',a)
c
      end
c
C> \brief Retrieve the next field as a character string
C>
C> Retrieves the next token as a text string. Any enclosing quotes are
C> removed from the field. If the value does not fit the buffer provided
C> an error condition is raised. The value of the buffer is changed
C> only if the function is successfull.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function inp_a(a)
      implicit none
#include "inpP.fh"
      integer i1, i2, length
      character*(*) a !< [Output] The text value of the next field
c     
c     Return field as character string, minus any enclosing quotes
c     with an error if it does not fit
c     
      ierrpos = -1
      errmsg = ' '
      if(jrec .ge. jump) then
         a = xblnk
         inp_a = .false.
         ierrpos = 0
         errmsg = 'at end of line looking for character string'
         return
      endif
      i1 = istrt(jrec+1)
      i2 = istrt(jrec+1)+inumb(jrec+1)-1
      if (ia(i1:i1).eq.xquote .and. ia(i2:i2).eq.xquote) then
         i1 = i1+1
         length = inumb(jrec+1)-2
      else
         length = inumb(jrec+1)
      endif
      if (len(a) .lt. length) then
         a = xblnk
         inp_a = .false.
         ierrpos = 0
         errmsg = 'inp_a: string is too large for argument'
         return
      else
         jrec = jrec + 1
         a = ia(i1:i1+length-1)
         inp_a = .true.
         return
      endif
c     
      end 
c
C> \brief Retrieve the next field as a character string
C>
C> Retrieves the next token as a text string. Any enclosing quotes are
C> removed from the field. If the value does not fit the buffer provided
C> the value is silently truncated to as much as will fit.
C>
C> \return Return .true. if successfull, and .false. otherwise.
c
      logical function inp_a_trunc(a)
      implicit none
#include "inpP.fh"
      integer i1, i2, length
      character*(*) a !< [Output] The text value of the next field
c     
c     Return field as character string, minus any enclosing quotes
c     quietly truncating if it does not fit
c     
      ierrpos = -1
      errmsg = ' '
      if(jrec .ge. jump) then
         a = xblnk
         inp_a_trunc = .false.
         ierrpos = 0
         errmsg = 'at end of line looking for character string'
         return
      endif
      i1 = istrt(jrec+1)
      i2 = istrt(jrec+1)+inumb(jrec+1)-1
      if (ia(i1:i1).eq.xquote .and. ia(i2:i2).eq.xquote) then
         i1 = i1+1
         length = inumb(jrec+1)-2
      else
         length = inumb(jrec+1)
      endif
      jrec = jrec + 1
      a = ia(i1:i1+length-1)
      inp_a_trunc = .true.
      return
c     
      end
c
C> \brief Retrieve the floating point value of the next field
C>
C> Interpret the next field as a floating point value and assign
C> its value to the buffer provided. If the next field does not 
C> represent a valid floating point value an error condition is
C> raised. The buffer is not changed unless the function is
C> successfull.
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function inp_f (buf)
      implicit none
#include "inpP.fh"
      double precision buf !< [Output] The floating point value of the
                           !< next field
      double precision ten, dtmp
      integer i1, i2, ie2, isign, ie, iexp, ie1, itmp, i, j
      logical orep
      character*1 xchar(17)
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1     ,'+','-','.','e','d','E','D'/
      data ten/10.0d0/
c     
      ierrpos = -1
      errmsg = ' '
      dtmp=0.0d0
      if (jrec.ge.jump) then
         inp_f = .false.
         errmsg = 'at end of line looking for floating point number'
         ierrpos=-1
         return
      endif
      jrec=jrec+1
      i1=istrt(jrec)
      i2=i1+inumb(jrec)-1
      ie2=i2
c...  sign
      isign=1
      if (ia(i1:i1).eq.xchar(12))isign=-1
      if (ia(i1:i1).eq.xchar(12).or.ia(i1:i1).eq.xchar(11)) i1=i1+1
c...  exponent
      do ie=i1+1,i2
         if (ia(ie:ie).eq.xchar(14) .or. ia(ie:ie).eq.xchar(15) .OR.
     $      ia(ie:ie).eq.xchar(16) .or. ia(ie:ie).eq.xchar(17)) goto 20
      enddo
      iexp=0
      go to 50
 20   i2=ie-1
      iexp=1
      ie1=ie+1
      if (ia(ie1:ie1).eq.xchar(12))iexp=-1
      if (ia(ie1:ie1).eq.xchar(12).or.ia(ie1:ie1).eq.xchar(11))
     *     ie1=ie1+1
      itmp=0
      do i=ie1,ie2
         do j=1,10
            if (ia(i:i).eq.xchar(j)) go to 41
         enddo
         goto 100
 41      itmp=itmp*10+j-1
      enddo
      iexp=iexp*itmp
c.... the number itself
 50   orep=.false.
      do i=i1,i2
         if(ia(i:i).ne.xchar(13)) then
            do j=1,10
               if (ia(i:i).eq.xchar(j)) go to 70
            enddo
            goto 100
 70         dtmp=dtmp*ten+ dble(j-1)
         else
            if(orep)go to 100
            iexp=iexp+i-i2
            orep=.true.
         endif
      enddo
      dtmp = dtmp * dble(isign) * ten**iexp
      inp_f = .true.
      buf = dtmp
      return
c
 100  inp_f = .false.
      jrec = jrec-1             ! Position to re-read the field
      ierrpos = i
      errmsg = 'illegal character reading floating point number'
c     
      end
C
C> \brief Mark an input error at the beginning of the current field
C>
C> When an input error is detected we want to generate an error message
C> that tells the user which field caused a problem. For this purpose
C> the program needs to be able to specify the field for which an 
C> error was detected. The routine provides that function, in addition
C> it allows an error message to be specified as well.
c
      subroutine inp_mark_err(message)
      implicit none
#include "inpP.fh"
      character*(*) message !< [Input] The error message
c
c     Mark an input error at the beginning of the current input field
c
      ierrpos = istrt(min(max_field,max(1,jrec)))
      errmsg  = message
c
      end
c
C> \brief Retrieve the integer of the next field
C>
C> Interpret the next field as an integer value and assign
C> its value to the buffer provided. If the next field does not 
C> represent a valid integer value an error condition is
C> raised. The buffer is not changed unless the function is
C> successfull.
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function inp_i(jbuf)
      implicit none
#include "inpP.fh"
      character*1 xchar(12)
      integer n, ifact, ist, nstrt, i, j
      character*1 xtemp
      integer jbuf !< [Output] The integer value read
      integer jtmp
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1     ,'+','-'/
c     
c     subroutine for reading integers from the array ia,
c     starting at ia(istrt(jrec)) and going on for inumb(jrec))
c     elements. plus signs are ignored, the answer is accumulated
c     in jtmp
c     
      ierrpos = -1
      errmsg = ' '
      jtmp = 0
      if(jrec.ge.jump) then
         inp_i = .false.
         ierrpos = -1
         errmsg = 'at end of line looking for integer'
         return
      endif
      jrec = jrec + 1
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      do i = 1,n
         xtemp = ia(nstrt:nstrt)
         do j=1,12
            if(xchar(j).eq.xtemp)go to 130
         enddo
         goto 120
c
 130     if(j.ge.11) then
            if(nstrt.ne.ist)go to 120
            if(j.ge.12)jtmp=-jtmp
            go to 160
         endif
         jtmp=jtmp+(j-1)*ifact
         ifact = ifact * 10
         nstrt=nstrt-1
      enddo
 160  continue
      inp_i = .true.
      jbuf = jtmp
      return
c
 120  ierrpos = nstrt
      errmsg  = 'illegal character when reading integer'
      inp_i = .false.
      jrec = jrec-1
      return
c     
      end
c
C> \brief Retrieve the logical value of the next field
C>
C> Interpret the next field as a logical value and assign
C> its value to the buffer provided. If the next field does not 
C> represent a valid logical value an error condition is
C> raised. Valid string representations of the logical values
C> true and false are given in the table below:
C>
C>    | True   | False   |
C>    |:------:|:-------:|
C>    | .true. | .false. |
C>    | true   | false   |
C>    | t      | f       |
C>    | yes    | no      |
C>    | y      | n       |
C>    | on     | off     |
C>
C> The value of the buffer is changed only if the function is
C> successfull.
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function inp_l(lbuf)
      implicit none
#include "inpP.fh"
#include "errquit.fh"
      integer maxlen
      parameter (maxlen = 7)
      character*7 xtoken
      integer n, ifact, ist, nstrt, i, j
      logical lbuf !< [Output] The logical value read
c     
      if (maxlen.ne.len(xtoken))
     +  call errquit('inp_l: coding error',0,UERR)
c
      inp_l = .false.
      ierrpos = -1
      errmsg = ' '
      if(jrec.ge.jump) then
         inp_l = .false.
         ierrpos = -1
         errmsg = 'at end of line looking for logical'
         return
      endif
      jrec = jrec + 1
      n    = inumb(jrec)
      ist  = istrt(jrec)
      if (n.gt.maxlen) then
c
c       the token is longer than the longest valid value of logical
c       so this cannot match
c
        inp_l = .false.
        jrec = jrec - 1   ! position cursor to reread this token
        return
      endif
      xtoken = ia(ist:ist+n-1)
      call inp_lcase(xtoken)
c
c     Check whether the string represents a valid value of .true.
c
      if (xtoken(1:n).eq.'.true.') then
        inp_l = .true.
        lbuf  = .true.
      else if (xtoken(1:n).eq.'true') then
        inp_l = .true.
        lbuf  = .true.
      else if (xtoken(1:n).eq.'t') then
        inp_l = .true.
        lbuf  = .true.
      else if (xtoken(1:n).eq.'yes') then
        inp_l = .true.
        lbuf  = .true.
      else if (xtoken(1:n).eq.'y') then
        inp_l = .true.
        lbuf  = .true.
      else if (xtoken(1:n).eq.'on') then
        inp_l = .true.
        lbuf  = .true.
c
c     Check whether the string represents a valid value of .false.
c
      else if (xtoken(1:n).eq.'.false.') then
        inp_l = .true.
        lbuf  = .false.
      else if (xtoken(1:n).eq.'false') then
        inp_l = .true.
        lbuf  = .false.
      else if (xtoken(1:n).eq.'f') then
        inp_l = .true.
        lbuf  = .false.
      else if (xtoken(1:n).eq.'no') then
        inp_l = .true.
        lbuf  = .false.
      else if (xtoken(1:n).eq.'n') then
        inp_l = .true.
        lbuf  = .false.
      else if (xtoken(1:n).eq.'off') then
        inp_l = .true.
        lbuf  = .false.
c
c     Else we cannot interpret this field as a logical
c
      else
        inp_l = .false.
        jrec = jrec-1   ! position cursor to reread this token
      endif
      return
c     
      end
c
C> \brief Compare two strings
C>
C> This function compares two strings, either case sensitively or
C> case insensitively. For the two strings to match they have to be
C> of equal length at least.
C>
C> \return Return .true. if the two strings match, and .false. otherwise
C
      logical function inp_compare(ocase, a, b)
      implicit none
      logical ocase   !< [Input] If .true. do case sensitive comparison
      character*(*) a !< [Input] String A to match B
      character*(*) b !< [Input] String B to match A
      integer la, lb, i
      character*1 atest, btest
      integer inp_strlen
      external inp_strlen
c
      inp_compare = .false.
      la = inp_strlen(a)
      lb = inp_strlen(b)
      if (la .ne. lb) then      ! use .gt. for short match
         return
      else if (ocase) then
         inp_compare = a(1:la) .eq. b(1:lb)
         return
      else
         do i = 1, la
            atest = a(i:i)
            btest = b(i:i)
            call inp_lcase(atest)
            call inp_lcase(btest)
            if (atest.ne.btest) return
         enddo
         inp_compare = .true.
         return
      endif
c
      end
c
C> \brief Matches a string to the elements of an array of strings
C>
C> This function looks for a string in an array of strings. A match
C> is found if the search string appears once and only once in the 
C> array. The index in the array where the match was found is retrieved.
C> If the function fails to find a unique match then
C>
C> - ind = -1 if no occurence of string was found in array
C>
C> - ind = 0 if multiple occurences of string were found in array
C>
C> \return Return .true. if one match was found, and .false. otherwise
c
      logical function inp_match(nrec, ocase, test, array, ind)
      implicit none
      integer nrec           !< [Input] The length of the array of
                             !< strings
      logical ocase          !< [Input] If .true. do case sensitive
                             !< comparison
      logical inp_compare
      character*(*) test     !< [Input] String to find in array
      character*(*) array(*) !< [Input] Array of strings to search
      integer ind            !< [Output] The index of the element of
                             !< array that matches test
      integer i, j, l, inp_strlen
      external inp_compare, inp_strlen
c
      l = inp_strlen(test)
      inp_match = .false.
      ind = -1
c     
      do i=1,nrec
         if (inp_compare(ocase, test(1:l), array(i))) then
            if (inp_match) then
               inp_match = .false. ! Ambiguity
               ind = 0
c
               write(6,1) test(1:l), (array(j),j=1,nrec)
 1             format('inp: ambiguous match for ', a,', in:'/
     $              100(1x,a/))
c
               return
            else
               inp_match = .true. ! First match
               ind = i
            endif
         endif
      enddo
c     
      end
C
C> \brief Checks if string A is contained in string B
C>
C> Checks whether string A is contained in string B. If so, then the 
C> position where A was found in B is retrieved. The strings can be
C> compared case sensitively or case insensitive. If no matches were
C> found the position retrieved is -1.
C>
C> \return Return .true. if string A is contained in B, and .false.
C> otherwise
c
      logical function inp_contains(ocase, a, b,ipos)
c
c  check if string a is contained in b. return starting
c  location of string a in b.
c
      implicit none
      logical ocase   !< [Input] If .true. do case sensitive comparison
      character*(*) a !< [Input] String A to be found in B
      character*(*) b !< [Input] String B to be searched
      integer ipos    !< [Output] The position where A was found in B
      integer la, lb, i, j
      character*1 atest, btest
      integer inp_strlen
      external inp_strlen
c
      ipos = -1
      inp_contains = .false.
      la = inp_strlen(a)
      lb = inp_strlen(b)
      if (la .gt. lb) then   
         return
      else if (ocase) then
         do i = 1, lb - la + 1
            inp_contains = a(1:la) .eq. b(i:i+la)
            if (inp_contains) then
               ipos = i
               return
            endif
         enddo
         return
      else
         do j = 0, lb - la 
            do i = 1, la
               atest = a(i:i)
               btest = b(j+i:j+i)
               call inp_lcase(atest)
               call inp_lcase(btest)
               if (atest.ne.btest) goto 00011
            enddo
            inp_contains = .true.
            ipos = j + 1
            return
00011    continue
         enddo
         return
      endif
c
      end
C
C> \brief Move the cursor one field back
C>
C> Moves the cursor one field back so that a field can be read again
C> with the next read call.
c
      subroutine inp_prev_field()
      implicit none
#include "inp.fh"
c
      call inp_set_field(max(0,inp_cur_field()-1))
c
      end
c
C> \brief Returns the length of the contents of a string
C>
C> In Fortran the length of a string is determined by the amount of 
C> memory allocated for it. Often one just wants to know how long the
C> value of the string is without any trailing spaces. Which is what
C> this function returns.
C>
C> \return The length of the value of a string without any trailing 
C> spaces
c
      integer function inp_strlen(a)
      implicit none
#include "inpP.fh"
      character*(*) a !< [Input] The string
      integer i
      logical ois_ws
      intrinsic len
      character*1 xtest
      ois_ws(xtest) = (xtest.eq.xblnk .or. xtest.eq.xtab)
c
      do i = len(a),1,-1
         if (.not. ois_ws(a(i:i))) goto 10
      enddo
c
 10   inp_strlen = i
c
      end
c
C> \brief Convert a string to lower case
C>
C> Replaces the contents of a string with the same value but in lower
C> case characters only.
c
      subroutine inp_lcase(string)
      implicit none
#include "errquit.fh"
      character*(*) string !< [In/Output] On input the string to be 
                           !< converted, and on output the same string
                           !< in lower case characters
      intrinsic ichar, len
      integer i, length, uca, ucz, lca, shift, test
c
      uca = ichar('A')          ! MUST be uppercase A
      ucz = ichar('Z')          ! MUST be uppercase Z
      lca = ichar('a')          ! MUST be lowercase a
      shift = lca - uca
      if (shift .eq. 0) 
     $     call errquit('inp_lcase: check case of program source', 0,
     &       INPUT_ERR)
c
      length = len(string)
      do i = 1, length
         test = ichar(string(i:i))
         if (test.ge.uca .and. test.le.ucz) 
     $        string(i:i) = char(test+shift)
      enddo
c
      end
c
C> \brief Remove any leading spaces from a string
C>
C> If a string has any white space characters at the start, remove
C> those characters and return the string starting with non-white space
C> characters.
c
      subroutine inp_adjustl(a)
      implicit none
#include "errquit.fh"
#include "inpP.fh"
      character*(*) a !< [In/Output] On input the string to be adjusted,
                      !< on output the same string but without any 
                      !< leading white space characters.
      intrinsic  len
      logical ois_ws
      character*1 xtest
      integer i, length,s
c
      ois_ws(xtest) = (xtest.eq.xblnk .or. xtest.eq.xtab)
c
      length = len(a)
      do i = 1, length
        if (.not. ois_ws(a(i:i))) goto 10
      end do
10    continue
      s=i-1
      do i=s+1,length
       a(i-s:i-s) = a(i:i) 
       a(i:i)=""
      enddo
c
      end
c
C> \brief Convert a string to upper case
C>
C> Replaces the contents of a string with the same value but in upper
C> case characters only.
c
      subroutine inp_ucase(string)
      implicit none
#include "errquit.fh"
      character*(*) string !< [In/Output] On input the string to be 
                           !< converted, and on output the same string
                           !< in upper case characters
      intrinsic ichar, len
      integer i, length, lca, lcz, uca, shift, test
c
      lca = ichar('a')          ! MUST be lowercase A
      lcz = ichar('z')          ! MUST be lowercase Z
      uca = ichar('A')          ! MUST be uppercase a
      shift = uca - lca
      if (shift .eq. 0) 
     $     call errquit('inp_ucase: check case of program source', 0,
     &       INPUT_ERR)
c
      length = len(string)
      do i = 1, length
         test = ichar(string(i:i))
         if (test.ge.lca .and. test.le.lcz) 
     $        string(i:i) = char(test+shift)
      enddo
c
      end
c
C> \brief Check if any of the strings specified match any of the fields
C> in the input file
C>
C> Searches the full input file for the presence of any of the provided
C> strings. The search stops at the first field to match any of the 
C> strings. The cursor is moved to just before the matching field.
C>
C> \return Return .true. if a string matches a field, and .false. 
C> otherwise
c
      logical function inp_search(ocase, z, nz)
      implicit none
#include "errquit.fh"
      integer nz          !< [Input] The length of z
      character*(*) z(nz) !< [Input] The array of strings
      logical ocase       !< [Input] If .true. do case sensitive
                          !< matching
      character*1024 tmp
      integer inp_strlen
      logical inp_read, inp_a, inp_compare
      external inp_read, inp_a, inp_compare, inp_strlen
c
      integer i
      integer maxz
      parameter (maxz = 100)
      integer length(maxz)
c
      if (maxz .lt. nz)
     $     call errquit('inp_search: hard dim fail',nz, INPUT_ERR)
      do i = 1, nz
         length(i) = inp_strlen(z(i))
      enddo
c
 10   if (inp_read()) then
         if (inp_a(tmp)) then
            do i = 1, nz
               if (inp_compare(ocase, z(i)(1:length(i)), tmp)) then
                  call inp_prev_field()
                  inp_search = .true.
                  return
               endif
            enddo
         endif
         goto 10
      endif
c
      inp_search = .false.
c
      end
c
C> \brief Search the input for a string using a fast algorithm
C>
C> Searches the input file for a given string. The given string must be
C> at least 3 characters long. To make the algorithm fast the search is
C> simplified in a number of ways:
C>
C> - All searches are case sensitive
C>
C> - Only attempt to match against the beginning of a line
C>
C> - Ignore continuation line, comments, quotes, and such like
C>
C> On return the cursor is positioned at the beginning of the line
C> where the match was found. If no match was found the cursor is
C> positioned at the end of the file.
C>
C> \return Return .true. if the string was found in the input, and
C> .false. otherwise.
C
      logical function inp_search_fast(z)
      implicit none
#include "errquit.fh"
#include "inpP.fh"
      character*(*) z !< [Input] The string
      integer length
      integer inp_strlen
      logical inp_read, inp_a, inp_compare
      external inp_read, inp_a, inp_compare, inp_strlen
c
c     Quicker search that 
c
c     1) matches case 
c     2) assumes the token being searched for is at the beginning of the line
c     3) ignores continutation lines, comments, quotes etc. 
c     4) Still attempts to track line numbers EOF.
c
c     Only called from inside the basis set input routine?
c
      inp_search_fast = .false.
      length = max(3,inp_strlen(z)) ! 3 for EOF/eof detection
c
 10   read(iread,'(a)',end=300) ja(1:length)
      input_line = input_line + 1
      if (ja(1:3).eq.'EOF' .or. ja(1:3).eq.'eof') goto 300
      if (z(1:length) .eq. ja(1:length)) then
         backspace(iread)       ! Re-read line with full input routine
         input_line = input_line - 1
         if (.not. inp_read()) call errquit('inp_search_fast: inp?',0,
     &       INPUT_ERR)
         inp_search_fast = .true.
         return
      endif
      goto 10
c
 300  oswit = .true.            ! EOF code copied from inp_read
      ierrpos = -1
      errmsg = 'unexpected end of data file'
      jump=0
      jrec=0
c
      end
c
C> \brief Push the state for the current input file on stack
C>
C> Save the state of the current input file on a stack so that the code
C> can start reading an include file. Currently the stack can hold
C> only 3 entries.
C
      subroutine inp_save_state
#include "inpP.fh"
c
      character*1 cdata(1)
      integer idata(1)
      integer fdata(2)
      integer i
c
      equivalence (cdata(1),ia)
      equivalence (idata(1),jrec)
      equivalence (fdata(1),iread)
c
#ifdef SOLARIS
c     Needed with -stackvar compilation option in pre WS5 compilers
c      save cdata, idata, fdata
#endif
c
      if (save_level.lt.0 .or. save_level.gt.3) call errquit
     $     ('inp_save_sate: invalid level ', save_level, INPUT_ERR)
      save_level = save_level + 1
c
      do i = 1, 3*max_width+80+6
         csave(i,save_level) = cdata(i)
      enddo
      do i = 1, 9+4*max_field
         isave(i,save_level) = idata(i)
      enddo
      do i = 1, 2
         fsave(i,save_level) = fdata(i)
      enddo
c
      end
c
C> \brief Pop the state of an input file from the stack
C>
C> Restore the state of an input file from the stack so that the code
C> can start continue reading a file after finishing with an include 
C> file.
C
      subroutine inp_restore_state
#include "inpP.fh"
c
      character*1 cdata(1)
      integer idata(1)
      integer fdata(2)
      integer i
c
      equivalence (cdata(1),ia)
      equivalence (idata(1),jrec)
      equivalence (fdata(1),iread)
#ifdef SOLARIS
c     Needed with -stackvar compilation option in pre WS5 compilers
c      save cdata, idata, fdata
#endif
c
      if (save_level.lt.1 .or. save_level.gt.3) call errquit
     $     ('inp_restore_sate: invalid level ', save_level, INPUT_ERR)
c
      do i = 1, 3*max_width+80+6
         cdata(i) = csave(i,save_level)
      enddo
      do i = 1, 9+4*max_field
         idata(i) = isave(i,save_level)
      enddo
      do i = 1, 2
         fdata(i) = fsave(i,save_level)
      enddo
c
      save_level = save_level - 1
c
      end
c
C> \brief Copy as much of the current input line as possible into a
C> string
C>
C> Copies as much of the current input line as will fit into the
C> string Z. The number of characters copied is returned in LEN, and
C> the success or failure of this operation is reported in SUCCESS.
C>
      subroutine inp_cline(z, len, success)
      implicit none
#include "inpP.fh"
      character*(*) z !< [Output] The current input line
      logical success !< [Output] If .true. to copy operation succeeded
      integer len     !< [Output] The number of characters copied to Z
      integer flen
c     
c     set the variable z to be as much of the current input line
c     that it can hold
c     
      flen = max_width
      if (jump .gt. 0) then
         call c_cnvt(len, flen, z, ia)
         success = .true.
         ierrpos = -1
         errmsg = ' '
      else
         errmsg = 'no input line available'
         ierrpos = -1
         success = .false.
      endif
      return
      end
C
C> \brief Find the next token on a string
C>
C> Returns the number of the start and end character of the
C> next token in the character string. Tokens are separated
C> by one of the characters in SEP. Note that all characters
C> in SEP are used including any trailing blanks.
C>
C> Before the first call initialize ISTART to zero, and leave
C> ISTART and IEND <b>UNCHANGED</b> for subsequent calls.
C> Repeated calls return the next token and .true., or .false. if
C> there are no more tokens. The separators may be changed
C> between calls.
C>
C> No internal state is maintained (which is ISTART and IEND
C> must not be modified between calls) so multiple strings
C> may be parsed simultaneously.
C>
C> E.g., to split LIST = 'robert:rick:jeff' into tokens separated 
C> by ':'. You execute
C>
C> ~~~~
C> istart = 0
C> do while (inp_strtok(list, ':', istart, iend))
C>    write(6,*) list(istart:iend)
C> enddo
C> ~~~~
C>
C> \return Returns .true. if a next token was found, and .false.
C> otherwise
c
      logical function inp_strtok(z, sep, istart, iend)
      implicit none
      character*(*) z   !< [Input] The string to parse
      character*(*) sep !< [Input] The token separators
      integer istart    !< [Output] The start of next token
      integer iend      !< [Output] The end of next token
c
c     Returns the number of the start and end character of the
c     next token in the character string.  Tokens are separated
c     by one of the characters in sep.  Note that all characters
c     in sep are used including any trailing blanks.
c
c     Before the first call initialize istart to zero, and leave
c     istart and iend UNCHANGED for subsequent calls.
c     Repeated calls return the next token and true, or false if
c     there are no more tokens.  The separators may be changed
c     between calls.
c
c     No internal state is maintained (which is istart and iend
c     must not be modified between calls) so multiple strings
c     may be parsed simultaneously.
c
c     E.g., to split list = 'robert:rick:jeff' into tokens separated 
c     by ':'. You execute
c
c     istart = 0
c  10 if (inp_strtok(list, ':', istart, iend)) then
c     write(6,*) list(istart:iend)
c     goto 10
c     endif
c
      integer i, k, length, nsep
c
      if (istart .eq. 0) then
         istart = 1
      else
         istart = iend + 1
      endif
c
c     Scan start forward to next non-separator
c
      length = len(z)
      nsep   = len(sep)
c
      do i = istart, length
         do k = 1, nsep
            if (z(i:i) .eq. sep(k:k)) goto 10
         enddo
         goto 20
 10      continue
      enddo
      inp_strtok = .false.      ! No more tokens
      return
c
 20   istart = i                ! Beginning of next token
c
c     Scan end forward to one-before next separator
c
      do i = istart+1, length
         do k = 1, nsep
            if (z(i:i) .eq. sep(k:k)) goto 30
         enddo
      enddo
 30   iend = i - 1
c
      inp_strtok = .true.
c
      end
C> @}

      
