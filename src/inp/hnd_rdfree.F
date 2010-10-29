c$$$      program test
c$$$      implicit none
c$$$#include "inp.fh"
c$$$c
c$$$      integer ir
c$$$      character*255 string
c$$$      integer ierr, i
c$$$c
c$$$      CHARACTER*80  PRSWRD      ! Returns tokens as ASCII
c$$$      INTEGER       FLGWRD      ! Returns type (1=float,2=int,3=string)
c$$$      INTEGER       NUMWRD      ! Returns number of words
c$$$      INTEGER       NUMCHR      ! Returns length of character string
c$$$      COMMON/HND_FREERD/PRSWRD(40),NUMCHR(40),FLGWRD(40),NUMWRD
c$$$c
c$$$      character*8 types(3)
c$$$      data types/'float','integer','string'/
c$$$c
c$$$      call inp_init(5,6)
c$$$c
c$$$ 10   call hnd_rdfree(ir, string, ierr)
c$$$      if (ierr .eq. 0) then
c$$$         write(6,*) ' STRING = "',string(1:inp_strlen(string)),'"'
c$$$         do i = 1, numwrd
c$$$            write(6,*) i, ' ', types(flgwrd(i)), ' ', numchr(i), '=',
c$$$     $           prswrd(i)(1:numchr(i))
c$$$         enddo
c$$$         goto 10
c$$$      endif
c$$$c
c$$$      end
      subroutine hnd_rdfree(ir, string, ierr)
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "inp.fh"
      integer ir                ! [input] IGNORED
      character*(*) string      ! [output] input line with quoted char strings
      integer ierr              ! [output] Error (0=OK, 1=EOF)
c     
c     Routine to mimic hnd_rdfree
c     
      CHARACTER*80  PRSWRD      ! Returns tokens as ASCII
      INTEGER       FLGWRD      ! Returns type (1=float,2=int,3=string)
      INTEGER       NUMWRD      ! Returns number of words
      INTEGER       NUMCHR      ! Returns length of character string
      COMMON/HND_FREERD/PRSWRD(40),NUMCHR(40),FLGWRD(40),NUMWRD
c     
      character*255 aval
      double precision dval
      integer ival
      integer next, i
c$$$  character*8 types(3)
c$$$  data types/'float','integer','string'/
c     
      ierr = 0
      if (.not. inp_read()) then
         if (inp_eof()) then
            ierr = 1
            return              ! Return EOF condition
         else
            call errquit('hnd_rdfree: unknown error reading input',0,
     &       INPUT_ERR)
         endif
      endif
c
      string = ' '
      next = 1
c     
      numwrd = inp_n_field()
      do i = 1, numwrd
c     
c     Figure out the type
c     
         if (inp_i(ival)) then
            flgwrd(i) = 2
         else if (inp_f(dval)) then
            flgwrd(i) = 1
         else if (inp_a(aval)) then
            flgwrd(i) = 3
         else
            call errquit('hnd_rdfree: failed reading token',i,
     &       INPUT_ERR)
         endif
c     
c     Stuff the ASCII value and length into prswrd
c     
         call inp_prev_field()
         if (.not. inp_a(aval)) call errquit
     $        ('hnd_rdfree: failed reading ascii value',i, INPUT_ERR)
         if (flgwrd(i) .eq. 3) then
            prswrd(i)(1:1) = ''''
            prswrd(i)(2:)  = aval
            numchr(i) = inp_strlen(prswrd(i)) + 1
            prswrd(i)(numchr(i):numchr(i)) = ''''
         else
            prswrd(i) = aval
            numchr(i) = inp_strlen(prswrd(i))
         endif
c
c     Concatenate it onto the end of string
c
         if (next+numchr(i) .le. len(string)) then
            string(next:) = prswrd(i)
            next = next + numchr(i) + 1
         else
            call errquit('hnd_rdfree: string too short',len(string),
     &       INPUT_ERR)
         endif
c
      enddo
c
c$$$      write(6,*) ' STRING "',string(1:inp_strlen(string)),'"'
c$$$      do i = 1, numwrd
c$$$         write(6,*) i, ' ', types(flgwrd(i)), ' ', numchr(i), '=',
c$$$     $        prswrd(i)(1:numchr(i))
c$$$      enddo
c     
      end
