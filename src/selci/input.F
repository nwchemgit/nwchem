C> \ingroup selci
C> @{
      subroutine setinp(ir,iw)
*
* $Id$
*
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
#include "work.fh"
      common/iofile/iread,iwrite
c
      iread = ir
      iwrite = iw
      jwidth = 80
      jrec = -1
      oswit = .false. 
      nerr = 999
      nline = 0
      noline = 0
c
      end
      integer function nfield()
#include "work.fh"
c
c     return no. of fields in the input
c
      nfield = jump
c
      end
      integer function ifield()
#include "work.fh"
c
c     return no. of fields processed so far (0,...,nfield())
c
      ifield = jrec
c
      end
      subroutine setfld(ivalue)
#include "work.fh"
c
c     set field to be read next (ivalue=0,...,nfield)
c
      if (ivalue.lt.0 .or. ivalue.gt.nfield())
     $     call errquit('setfld: stupid field value',ivalue, INPUT_ERR)
      jrec = ivalue
c
      end
      subroutine getlin(z)
      character*(*) z
#include "workc.fh"
c
c     set the variable z to be as much of the input line
c     that it can hold
c
      z = ia
c
      end
      subroutine input
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
c.... this routine reads a data card and scans it for non - space fields
c    the number of fields is stored in jump, the starting point of a
c    field in istrt(i) and the number of characters in that field
c    in inumb(i).
#include "workc.fh"
#include "work.fh"
      common/iofile/iread,iwrite
      dimension xchar(4),ncol(40),xcol(2)
      data xblnk,xstp,xcol/' ','.',':','/'/
      data nchar,xchar/4,'*','?','>','<'/
c
      nline=nline+1
      if(nline.le.noline)go to 150
      if(oswit)
     * call errquit('unexpected end of data file',nline, INPUT_ERR)
 100  read(iread,'(a80)',end=300)ja
      do 70 i=1,nchar
      if(ja(1:1).eq.xchar(i))go to 110
 70   continue
      go to 80
 110  ia(1:79)=ja(2:80)
      write(iwrite,90)ia
 90   format(/
     *' comment :-',1x,a79)
      go to 100
 80   k=jwidth
      mark=0
      do 130 i=1,jwidth
      if(ja(i:i).ne.xcol(1).and.ja(i:i).ne.xcol(2)) go to 130
      mark=mark+1
      ncol(mark)=i
 130  continue
      noline=1
      if(mark.ne.0)go to 140
      nstart(noline)=1
      nend(noline)=jwidth
      go to 200
 140  i=ncol(mark)+1
      if(i.gt.jwidth)go to 180
      do 160 j=i,jwidth
      if(ja(j:j).ne.xblnk) go to 170
 160  continue
 180  k=ncol(mark)-1
      mark=mark-1
c
 170  noline=mark+1
      nstart(1)=1
      do 190 i=1,mark
      j=ncol(i)
      nend(i)=j-1
 190  nstart(i+1)=j+1
      nend(noline)=k
 200  nline=1
c
 150  jump=0
      jrec=0
      isw=0
      nbegin=nstart(nline)
      nfini=nend  (nline)
      iwidth=nfini-nbegin+1
      ia(1:iwidth)=ja(nbegin:nfini)
c
      j=iwidth+1
      if(j.ge.80)go to 210
      do 220 i=j,80
 220  ia(i:i)=xblnk
c
 210  do 40 i = 1,iwidth
      if(ia(i:i).eq.xblnk)go to 30
      if(isw)10,10,20
10    jump = jump +1
      istrt(jump) = i
      inumb(jump) = 0
      isw=1
20    inumb(jump) = inumb(jump) + 1
      go to 40
30    isw = 0
40    continue
      return
 300  oswit=.true.
      jump=0
      jrec=0
      return
      end
      subroutine errout(n)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "workc.fh"
#include "work.fh"
      common/iofile/iread,iwrite
      data xpt,xstp/    '*','.'/
      jrec=-1
      write(iwrite,50)ia
50    format(1x,a80)
      do 60 i=1,iwidth
60    ia(i:i)=xstp
      ia(n:n)=xpt
      write(iwrite,50)ia
      return
      end
      subroutine outrec
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*80 iaa
#include "workc.fh"
      common/iofile/iread,iwrite
#include "work.fh"
c
      iaa=ia(1:iwidth)
      write(iwrite,50)iaa
50    format(1x,a80)
      return
      end
      subroutine inpa8(zguf)
c.... this routine examines the contents of ia  and extracts a
c    character string of 8 chars. this string is stored in iguf .
c    characters beyond the eighth in any field are ignored
c      dimension ibuf(8)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "workc.fh"
#include "work.fh"
      data xblnk/' '/
      jrec = jrec + 1
      zguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.8)n=8
      zguf = ia(istrt(jrec):istrt(jrec)+n-1)
      return
      end
      subroutine inpa(zguf)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character*(*) zguf
#include "workc.fh"
#include "work.fh"
      data xblnk/' '/
      jrec = jrec + 1
      nguf=len(zguf)
      zguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.nguf)n=nguf
      zguf = ia(istrt(jrec):istrt(jrec)+n-1)
      return
      end
      subroutine inpa4(yguf)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "workc.fh"
#include "work.fh"
      data xblnk/' '/
      jrec = jrec + 1
      yguf = xblnk
      if(jrec .gt. jump) return
      n = inumb(jrec)
      if(n.gt.4)n=4
      yguf = ia(istrt(jrec):istrt(jrec)+n-1)
      return
      end
      subroutine inpf (buf)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
      dimension xchar(15)
#include "workc.fh"
#include "work.fh"
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1,'+','-','.','e','d'/
      data ten/10.0d0/
      buf=0.0d0
      jrec=jrec+1
      if (jrec.gt.jump) return
      i1=istrt(jrec)
      i2=i1+inumb(jrec)-1
      ie2=i2
c...  sign
      isign=1
      if (ia(i1:i1).eq.xchar(12))isign=-1
      if (ia(i1:i1).eq.xchar(12).or.ia(i1:i1).eq.xchar(11)) i1=i1+1
c...  exponent
      do 10 ie=i1,i2
      if (ia(ie:ie).eq.xchar(14) .or. ia(ie:ie).eq.xchar(15)) go to 20
10    continue
      iexp=0
      go to 50
20    i2=ie-1
      iexp=1
      ie1=ie+1
      if (ia(ie1:ie1).eq.xchar(12))iexp=-1
      if (ia(ie1:ie1).eq.xchar(12).or.ia(ie1:ie1).eq.xchar(11))
     * ie1=ie1+1
      ibuff=0
      do 40 i=ie1,ie2
      do 30 j=1,10
      if (ia(i:i).eq.xchar(j)) go to 40
30    continue
      goto 100
40    ibuff=ibuff*10+j-1
      iexp=iexp*ibuff
c.... the number itself
 50   orep=.false.
      do 90 i=i1,i2
      if(ia(i:i).eq.xchar(13)) go to 80
      do 60 j=1,10
      if (ia(i:i).eq.xchar(j)) go to 70
60    continue
      goto 100
70    buf=buf*ten+ dble(j-1)
      go to 90
 80   if(orep)go to 100
      iexp=iexp+i-i2
      orep=.true.
90    continue
      buf = buf* dble(isign) * ten**iexp
      return
100   call errout(i)
      call errquit('illegal character reading floating point number',
     $     i, INPUT_ERR)
      return
      end
      subroutine inpi(junke)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
      dimension xchar(12)
#include "workc.fh"
#include "work.fh"
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1,'+','-'/
c.... subroutine for reading integers from the array ia,
c     starting at ia(istrt(jrec)) and going on for inumb(jrec))
c     elements. plus signs are ignored, the answer is accumulated
c     in jbuf and transferred to junke
      jbuf = 0
      jrec = jrec + 1
      if(jrec.gt.jump)go to 160
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      do 150 i = 1,n
      xtemp = ia(nstrt:nstrt)
      do 110 j=1,12
      if(xchar(j).eq.xtemp)go to 130
110   continue
120   call errout(nstrt)
      call errquit('illegal character when reading integer',nstrt,
     &       INPUT_ERR)
130   if(j.lt.11)go to 140
      if(nstrt.ne.ist)go to 120
      if(j.ge.12)jbuf=-jbuf
      go to 160
140   jbuf=jbuf+(j-1)*ifact
      ifact = ifact * 10
150   nstrt=nstrt-1
160   junke=jbuf
      return
      end
      subroutine inpwid(iwid)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
#include "errquit.fh"
#include "workc.fh"
#include "work.fh"
      common/iofile/iread,iwrite
      if((iwid.lt.1).or.(iwid.ge.81))call errquit(
     1 'illegal line width specified',iwid, INPUT_ERR)
      iwidth=iwid
      write(iwrite,10)iwidth
10    format(///' input line width set to',i4,' characters')
      return
      end
      character*4 function ytrunc(ztext)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      ytrunc(1:4)=ztext(1:4)
      return
      end
      integer function locatc(label,nf,itext)
      implicit real*8    (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      character * (*) label,itext
      dimension label(*)
      do 1 i=1,nf
      if(label(i).eq.itext)go to 2
 1    continue
      locatc=0
      return
 2    locatc=i
      return
      end
C> @}
