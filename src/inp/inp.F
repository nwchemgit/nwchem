      block data inp_data
C$Id: inp.F,v 1.50 2008-10-23 20:19:43 d3p307 Exp $
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
#if defined(DECOSF) || defined(SGITFP) || defined(SGI_N32) || defined(CRAY) || defined(HPUX) || defined(WIN32) || defined(PSCALE) ||defined(GCC4) 
      data xtab  /'	'/      ! Tab ... no backslash necessary
#elif (defined(LINUX) || defined(MACX)) && !defined(PGLINUX) && !defined(XLFLINUX) &&!defined(GCC4)
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
      subroutine inp_init(ir,iw)
      implicit none
#include "inpP.fh"
      integer ir, iw
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
      integer function inp_n_field()
      implicit none
#include "inpP.fh"
c     
c     return no. of fields in the input line ... 0 = EOF
c     
      inp_n_field = jump
c     
      end
      integer function inp_cur_field()
      implicit none
#include "inpP.fh"
c     
c     return no. of fields processed so far (0,...,inp_n_field())
c     
      inp_cur_field = jrec
c     
      end
      subroutine inp_set_field(ivalue)
      implicit none
#include "errquit.fh"
#include "inpP.fh"
      integer ivalue
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
      logical function inp_line(z)
      implicit none
#include "inpP.fh"
      character*(*) z
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
      logical function inp_eof()
      implicit none
#include "inpP.fh"
c      
      inp_eof = oswit
c
      end
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
      logical function inp_a(a)
      implicit none
#include "inpP.fh"
      integer i1, i2, length
      character*(*) a
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
      logical function inp_a_trunc(a)
      implicit none
#include "inpP.fh"
      integer i1, i2, length
      character*(*) a
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
      logical function inp_f (buf)
      implicit none
#include "inpP.fh"
      double precision ten, buf, dtmp
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
      subroutine inp_mark_err(message)
      implicit none
#include "inpP.fh"
      character*(*) message
c
c     Mark an input error at the beginning of the current input field
c
      ierrpos = istrt(min(max_field,max(1,jrec)))
      errmsg  = message
c
      end
      logical function inp_i(jbuf)
      implicit none
#include "inpP.fh"
      character*1 xchar(12)
      integer n, ifact, ist, nstrt, i, j
      character*1 xtemp
      integer jbuf, jtmp
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
      logical function inp_compare(ocase, a, b)
      implicit none
      logical ocase
      character*(*) a, b
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
      logical function inp_match(nrec, ocase, test, array, ind)
      implicit none
      integer nrec, ind
      logical ocase, inp_compare
      character*(*) test, array(*)
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
      logical function inp_contains(ocase, a, b,ipos)
c
c  check if string a is contained in b. return starting
c  location of string a in b.
c
      implicit none
      logical ocase
      character*(*) a, b
      integer la, lb, i, j, ipos
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
            return
00011    continue
         enddo
         return
      endif
c
      end
      subroutine inp_prev_field()
      implicit none
#include "inp.fh"
c
      call inp_set_field(max(0,inp_cur_field()-1))
c
      end
      integer function inp_strlen(a)
      implicit none
#include "inpP.fh"
      character*(*) a
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
      subroutine inp_lcase(string)
      implicit none
#include "errquit.fh"
      character*(*) string
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
      subroutine inp_ucase(string)
      implicit none
#include "errquit.fh"
      character*(*) string
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
      logical function inp_search(ocase, z, nz)
      implicit none
#include "errquit.fh"
      integer nz
      character*(*) z(nz)
      logical ocase
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
      logical function inp_search_fast(z)
      implicit none
#include "errquit.fh"
#include "inpP.fh"
      character*(*) z
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
      subroutine inp_save_state
#include "inpP.fh"
c
      character*1 cdata(1)
      integer idata(1)
      integer fdata(1)
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
      subroutine inp_restore_state
#include "inpP.fh"
c
      character*1 cdata(1)
      integer idata(1)
      integer fdata(1)
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

      subroutine inp_cline(z, len, success)
      implicit none
#include "inpP.fh"
      character*(*) z
      logical success
      integer len, flen
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
      logical function inp_strtok(z, sep, istart, iend)
      implicit none
      character*(*) z           ! [input] string to parse
      character*(*) sep         ! [input] token separators
      integer istart, iend      ! [output] start/end of next token
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

      
