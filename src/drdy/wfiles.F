C $Id$
      Subroutine drdy_wfiles (swf,iflg)
C
      Implicit none
#include "errquit.fh"
      Integer iflg
      Double Precision swf
C
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "stdio.fh"
C
      Integer i,nstrs,ilen
      Character*1 chsval(5),schar(6),char2(255)
      Character*32 theory
      Character*255 cmd,cmd1
      Character*255 dummy_file
      Character*255 filein,fileout
      Double Precision absswf
C
      Integer util_system
      External util_system
C
      Equivalence (char2,cmd)
C
      logical status_rtdb, status_ignore
      status_rtdb   = rtdb_parallel(.false.)

      nstrs = 6
      if (swf.eq.sreact) then
         schar(1) = 's'
         schar(2) = 'r'
         schar(3) = 'e'
         schar(4) = 'a'
         schar(5) = 'c'
         schar(6) = 't'
      else if (swf.eq.sprod) then
         schar(1) = 's'
         schar(2) = 'p'
         schar(3) = 'r'
         schar(4) = 'o'
         schar(5) = 'd'
         schar(6) = 't'
      else
         if(swf.lt.0.0) then
            schar(1) = 'M'
         else
            schar(1) = 'P'
         endif
         absswf = abs(swf)
         call drdy_getstr(absswf,chsval,nstrs,2,2)
         do i = 1,nstrs
            schar(i+1) = chsval(i)
         enddo
         nstrs = nstrs + 1
      endif
      if (nstrs.ne.6) then
         write (fu6,*) ' expected string of length 6 but got',nstrs
         write (fu6,*) ' could not save NWChem files'
      else
         if (iflg.eq.0) then
            if(lgs(38).eq.1.or.lgs(38).ge.3) then
               if(.not.rtdb_cget(my_rtdb,'drdy:theory',1,theory))
     *           call errquit('drdy_wfiles:cannot get theory',555,
     &       RTDB_ERR)
               call drdy_get_movecname(my_rtdb,theory,dummy_file)
               ilen = inp_strlen(dummy_file)
               filein = dummy_file(1:ilen)
               fileout = dummy_file(1:ilen)//"."
               ilen = inp_strlen(fileout)
               do i=1,nstrs
                 fileout(i+ilen:i+ilen) = schar(i)
               enddo
               write(luout,*)'copy '//filein(1:inp_strlen(filein))//
     +                       ' --> '//fileout(1:inp_strlen(fileout))
               call util_file_copy(filein,fileout)
            endif
            if(lgs(38).ge.2) then
               call util_file_name('hess',  .false., .false.,dummy_file)
               ilen = inp_strlen(dummy_file)
               filein = dummy_file(1:ilen)
               fileout = dummy_file(1:ilen)//"."
               ilen = inp_strlen(fileout)
               do i=1,nstrs
                 fileout(i+ilen:i+ilen) = schar(i)
               enddo
               write(luout,*)'copy '//filein(1:inp_strlen(filein))//
     +                       ' --> '//fileout(1:inp_strlen(fileout))
               call util_file_copy(filein,fileout)
            endif
         else
            if(.not.rtdb_cget(my_rtdb,'drdy:sptheory',1,theory))
     *        call errquit('drdy_wfiles:cannot get sptheory',555,
     &       RTDB_ERR)
            call drdy_get_movecname(my_rtdb,theory,dummy_file)
            ilen = inp_strlen(dummy_file)
            filein = dummy_file(1:ilen)
            fileout = dummy_file(1:ilen)//".sp."
            ilen = inp_strlen(fileout)
            do i=1,nstrs
              fileout(i+ilen:i+ilen) = schar(i)
            enddo
            write(luout,*)'copy '//filein(1:inp_strlen(filein))//
     +                    ' --> '//fileout(1:inp_strlen(fileout))
            call util_file_copy(filein,fileout)
         endif
      endif
C
      status_ignore = rtdb_parallel(status_rtdb)
      return
600   format(10x,'save file, system command=  ',A255)
601   format(10x,'error saving file with system command=  ',A255)
      end
c
      subroutine drdy_get_movecname(rtdb,theory,movecs)
c
      Implicit none
#include "errquit.fh"
c
#include "inp.fh"
#include "rtdb.fh"
c
c Get the full name of the movecs file given the theory and rtdb
c
      integer rtdb
      character*32 theory, key
      character*255 movecs
c
      if (inp_compare(.false.,'scf',theory)) then
        key = 'scf:output vectors'
      elseif (inp_compare(.false.,'dft',theory)) then
        key = 'dft:output vectors'
      elseif (inp_compare(.false.,'mcscf',theory)) then
        key = 'mcscf:output vectors'
      elseif (inp_compare(.false.,'selci',theory)) then
        key = 'mcscf:output vectors'
      endif
c
c Actually get the movecs name
c
      if (.not. rtdb_cget(rtdb,key,1,movecs)) call errquit
     *   ('drdy_wfiles:failed getting movecs name',555, RTDB_ERR)
c
      end
