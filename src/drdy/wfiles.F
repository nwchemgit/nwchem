C $Id: wfiles.F,v 1.1 2001-01-19 17:43:15 windus Exp $
      Subroutine drdy_wfiles (swf,iflg)
C
      Implicit none
      Integer iflg
      Double Precision swf
C
#include "drdyP.fh"
C
      Integer i,n,nstrs
      Character*1 chsval(5),schar(6),char1(30),char2(30),char3(30)
      Character*30 cmd1,cmd2,cmd3
      Double Precision absswf
C
      Integer util_system
      External util_system
C
      Equivalence (char1,cmd1),(char2,cmd2),(char3,cmd3)
C
      Data cmd1/'/bin/mv gauss.out s      .gss '/
      Data cmd2/'/bin/cp Test.FChk s      .chk '/
      Data cmd3/'/bin/mv gausp.out s      .spc '/
C
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
         n = 19
         do i = 1,nstrs
            char1(i+n) = schar(i)
            char2(i+n) = schar(i)
            char3(i+n) = schar(i)
         enddo
C
         if (iflg.eq.0) then
            if(lgs(38).eq.1.or.lgs(38).ge.3) then
               if (util_system (cmd1) .eq. 0) then
                  Write (fu6,600) cmd1
               else
                  Write (fu6,601) cmd1
                  call errquit('drdy_wfiles: problem with system call',
     *                          555)
               endif
            endif
            if(lgs(38).ge.2) then
               if (util_system (cmd2) .eq. 0) then
                  Write (fu6,600) cmd2
               else
                  Write (fu6,601) cmd2
                  call errquit('drdy_wfiles: problem with system call',
     *                          555)
               endif
            endif
         else
            if (util_system (cmd3) .eq. 0) then
               Write (fu6,600) cmd3
            else
               Write (fu6,601) cmd3
               call errquit('drdy_wfiles: problem with system call',
     *                       555)
            endif
         endif
      endif
C
      return
600   format(10x,'save file, system command=  ',A30)
601   format(10x,'error saving file with system command=  ',A30)
      end
