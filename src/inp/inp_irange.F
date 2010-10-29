      logical function inp_irange(JFirst, JLast, JStride)
*
* $Id$
*
      implicit none
#include "inpP.fh"
      integer n, ifact, ist, nstrt, i, j
      character*1 xtemp
      integer jFirst, JLast, JStride, jtmp(3), term
      Logical Expect_Sep, Expect_Dig
      character*1 xchar(13)
      data xchar /'0','1','2','3','4','5','6','7','8','9'
     1     ,'+','-',':'/
c     
c     subroutine for reading an integer range specification (in Fortran90-
C     style triplet notation, stride optional) from the array ia,
c     starting at ia(istrt(jrec)) and going on for inumb(jrec))
c     elements. plus signs are ignored, the answer is accumulated
c     in jtmp.
C
C     F90 triple notation:  <first>[:<last>[:<stride>]]
C
C     Note that a simple integer "<first>" will be accepted and interpreted
C     as "<first>:<first>:<stride>"
c     
      ierrpos = -1
      errmsg = ' '
C
      term = 1
      jtmp(1) = 0
      jtmp(2) = 0
      jtmp(3) = 0
C
      if(jrec.ge.jump) then
         inp_irange = .false.
         ierrpos = -1
         errmsg = 'at end of line looking for integer range'
         return
      endif
      jrec = jrec + 1
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      Expect_Sep = .FALSE.
      Expect_Dig = .TRUE.   ! Should be a digit before anything else appears
      do i = 1,n
         xtemp = ia(nstrt:nstrt)
         do j=1,13
            if(xchar(j).eq.xtemp)go to 130
         enddo
         goto 120
c
 130     Continue
         If ( J .eq. 13) then                    ! ":" separating terms
            If ( Expect_Dig ) Goto 120
            Term = Term + 1
            IFact = 1
            Expect_Dig = .TRUE.   ! Must have a digit next (or EOR)
            Expect_Sep = .FALSE.
C
         ElseIf(j.eq.11 .OR. J .eq. 12) then     ! +/- signs
            If ( Expect_Sep .OR. Expect_Dig) Goto 120
            if(j.eq.12)jtmp(term) =-jtmp(term) 
            Expect_sep = .TRUE. ! Need separator (or EOR) next
C
         Else                                    ! Add new digit
            If ( Expect_Sep) Goto 120
            jtmp(term) =jtmp(term) +(j-1)*ifact
            ifact = ifact * 10
            Expect_Dig = .FALSE.
         EndIf
         nstrt=nstrt-1
      enddo
C
C     Finished the loop over the character in this field.  See if
C     we are missing anything (like ending with a : at the front).
C
      If ( Expect_Dig ) Goto 120
C
C     All done.  Now we must sort out how much of the triplet we
C     were given.  Remember, the terms are back to front.
C
 160  continue
      inp_irange = .true.
C
      If (term .eq. 1) then       ! Just an integer
         JFirst = JTmp(1)
         JLast  = JTmp(1)
         JStride = 1
      ElseIf (Term .eq. 2) then   ! Start and end, stride defaults to 1
         If ( JTmp(2) .gt. JTmp(1) ) Goto 125 ! Stride 1 will not work here
         JFirst = JTmp(2)
         JLast  = JTmp(1)
         JStride = 1
      ElseIf (Term .eq. 3) then   ! Full triplet
         If ( JTmp(1) .eq. 0) Goto 126 ! Stride 0 not allowed
         If ( JTmp(2) .ne. JTmp(3) ) then
            If ( (JTmp(2) - JTmp(3)) * JTmp(1) .lt. 0) Goto 127
         EndIf
         JFirst = JTmp(3)
         JLast  = JTmp(2)
         JStride = JTmp(1)
      EndIf
c
      return
c
 120  ierrpos = nstrt
      errmsg  = 'illegal character when reading integer range'
      inp_irange = .false.
      jrec = jrec-1
      return
c     
 125  ierrpos = ist
      errmsg = 'invalid range spec. -- first > last'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return
C
 126  ierrpos = ist
      errmsg = 'invalid range spec. -- stride can not be 0'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return
C
 127  ierrpos = ist
      errmsg = 'invalid range spec. -- stride incompatable with range'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return
C
      end
