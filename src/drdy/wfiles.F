C $Id: wfiles.F,v 1.2 2001-02-01 17:42:48 windus Exp $
      Subroutine drdy_wfiles (swf,iflg)
C
      Implicit none
      Integer iflg
      Double Precision swf
C
#include "drdyP.fh"
#include "nwc_drdyP.fh"
#include "inp.fh"
#include "rtdb.fh"
C
      Integer i,n,nstrs,ilen
      Character*1 chsval(5),schar(6),char1(255),char2(255),char3(255)
      Character*32 theory
      Character*255 cmd,cmd1,cmd2,cmd3
      Character*255 dummy_file
      Double Precision absswf
C
      Integer util_system
      External util_system
C
      Equivalence (char1,cmd1),(char2,cmd),(char3,cmd3)
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
         if (iflg.eq.0) then
            if(lgs(38).eq.1.or.lgs(38).ge.3) then
               if(.not.rtdb_cget(my_rtdb,'drdy:theory',1,theory))
     *           call errquit('drdy_wfiles:cannot get theory',555)
               call drdy_get_movecname(my_rtdb,theory,dummy_file)
               ilen = inp_strlen(dummy_file)
               cmd = 'cp '
               write(fu6,*) 'first cmd'
               write(fu6,600) cmd
               write(cmd1,'(a,'' '',a,'' '',a,''.'')') 
     *           cmd(1:inp_strlen(cmd)),dummy_file(1:ilen),
     *           dummy_file(1:ilen)
               write(fu6,*) 'second cmd'
               write(fu6,600) cmd1
               cmd = cmd1
               ilen = inp_strlen(cmd)
               do i=1,nstrs
                 char2(i+ilen) = schar(i)
               enddo
               write(fu6,*) 'fifth cmd'
               write(fu6,600) cmd
c              inquire(file=dummy_file,exist=oexists)
c              if (.not.oexists) then
c                write(fu6,*) 'The movecs file does not exist!'
c                call errquit(
c    *             'drdy_wfiles: problem saving movecs file',555)
c              endif
               if (util_system (cmd) .eq. 0) then
                  Write (fu6,600) cmd
               else
                  Write (fu6,601) cmd
                  call errquit('drdy_wfiles: problem with system call',
     *                          555)
               endif
            endif
            if(lgs(38).ge.2) then
               call util_file_name('hess',  .false., .false.,dummy_file)
               ilen = inp_strlen(dummy_file)
               cmd = 'cp '
               write(fu6,*) 'first cmd'
               write(fu6,600) cmd
               write(cmd1,'(a,'' '',a,'' '',a,''.'')')
     *           cmd(1:inp_strlen(cmd)),dummy_file(1:ilen),
     *           dummy_file(1:ilen)
               write(fu6,*) 'second cmd'
               write(fu6,600) cmd1
               cmd = cmd1
               ilen = inp_strlen(cmd)
               do i=1,nstrs
                 char2(i+ilen) = schar(i)
               enddo
               write(fu6,*) 'fifth cmd'
               write(fu6,600) cmd 
               if (util_system (cmd) .eq. 0) then
                  Write (fu6,600) cmd
               else
                  Write (fu6,601) cmd
                  call errquit('drdy_wfiles: problem with system call',
     *                          555)
               endif
            endif
         else
            if(.not.rtdb_cget(my_rtdb,'drdy:sptheory',1,theory))
     *        call errquit('drdy_wfiles:cannot get sptheory',555)
            call drdy_get_movecname(my_rtdb,theory,dummy_file)
            ilen = inp_strlen(dummy_file)
            cmd = 'cp '
               write(fu6,*) 'first cmd'
               write(fu6,600) cmd
               write(cmd1,'(a,'' '',a,'' '',a,''.'')')
     *           cmd(1:inp_strlen(cmd)),dummy_file(1:ilen),
     *           dummy_file(1:ilen)
               write(fu6,*) 'second cmd'
               write(fu6,600) cmd1
               cmd = cmd1
               ilen = inp_strlen(cmd)
               do i=1,nstrs
                 char2(i+ilen) = schar(i)
               enddo
               write(fu6,*) 'fifth cmd'
               write(fu6,600) cmd
            if (util_system (cmd) .eq. 0) then
               Write (fu6,600) cmd
            else
               Write (fu6,601) cmd
               call errquit('drdy_wfiles: problem with system call',
     *                       555)
            endif
         endif
      endif
C
      return
600   format(10x,'save file, system command=  ',A255)
601   format(10x,'error saving file with system command=  ',A255)
      end
c
      subroutine drdy_get_movecname(rtdb,theory,movecs)
c
      Implicit none
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
     *   ('drdy_wfiles:failed getting movecs name',555)
c
      end
