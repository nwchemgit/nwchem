      block data scrperm
      implicit none
#include "cfileprefix.fh"
c
      data scratch_dir /'null'/
      data permanent_dir /'null'/
c
      end
      subroutine util_file_info_rtdb(rtdb)
      implicit none
#include "errquit.fh"
#include "cfileprefix.fh"
#include "rtdb.fh"
      integer rtdb
c
c     Store the file prefix in the database.  If the scratch_dir
c     and permanent_dir have been set by the user then store them
c     otherwise restore them from the database.  They are stored
c     as 'scratch_dir' and 'permanent_dir' so that unset can
c     be used to delete them so that defaults may be restored.
c
      if (.not. rtdb_cput(rtdb, 'file_prefix', 1, file_prefix))
     $     call errquit('rtdb_put of file_prefix failed', 0, RTDB_ERR)
c
      if (scratch_dir .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'scratch_dir', 1, scratch_dir))
     $        call errquit('rtdb_put of scratch_dir failed', 0,
     &       RTDB_ERR)
      else
         if (.not. rtdb_cget(rtdb, 'scratch_dir', 1, scratch_dir))
     $        scratch_dir = ' '
      endif
c
      if (permanent_dir .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'permanent_dir', 1, permanent_dir))
     $        call errquit('rtdb_put of permanent_dir failed', 0,
     &       RTDB_ERR)
      else
         if (.not. rtdb_cget(rtdb, 'permanent_dir', 1, permanent_dir))
     $        permanent_dir = ' '
      endif
c
      end
      subroutine util_file_prefix_get(fullname)
      implicit none
#include "cfileprefix.fh"
      character*(*) fullname

      fullname = file_prefix
c
      end
      subroutine util_file_prefix_set(fullname)
      implicit none
#include "cfileprefix.fh"
      character*(*) fullname

      file_prefix = fullname
c
      end
      subroutine util_file_prefix(name, fullname)
      implicit none
#include "errquit.fh"
#include "cfileprefix.fh"
#include "inp.fh"
#include "stdio.fh"
      character*(*) name, fullname
c
c     prepend the file_prefix onto name as <file_prefix>.name
c     returning the result in fullname.
c
      if ((inp_strlen(name)+inp_strlen(file_prefix)+1) .gt.
     $     len(fullname)) then
         write(LuOut,*) ' file_prefix: name = ', name
         write(LuOut,*) ' file_prefix: prfx = ', file_prefix
         call util_flush(LuOut)
         call errquit('file_prefix: insufficient space ', len(fullname),
     &       DISK_ERR)
      endif
c
      fullname = ' '
      write(fullname,'(a,''.'',a)')
     $     file_prefix(1:inp_strlen(file_prefix)),
     $     name(1:inp_strlen(name))
c
      end
      subroutine util_file_name(stub, oscratch, oparallel, name)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"      
#include "cfileprefix.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) stub      ! [input] stub name for file
      logical oscratch        ! [input] true=scratch, false=permanent
      logical oparallel       ! [input] true=append .nodeid
      character*(*) name      ! [output] full filename
c

      call util_file_name0(stub, oscratch, oparallel, name, -1)
c
      end
      subroutine util_file_name0(stub, oscratch, oparallel, name, nodgs)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"      
#include "cfileprefix.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) stub      ! [input] stub name for file
      logical oscratch        ! [input] true=scratch, false=permanent
      logical oparallel       ! [input] true=append .nodeid
      character*(*) name      ! [output] full filename
      integer nodgs           ! [input] no. of digits for prefix
c
      character*(nw_max_path_len) dir, tmp
      integer ltmp, ldir, me
      logical util_file_parse_dir
      external util_file_parse_dir
c
      me = ga_nodeid()
c
      call util_directory_name(dir, oscratch, me)
c
*     write(LuOut,*) 'a stub= |',stub(1:inp_strlen(stub)),'|'
*     write(LuOut,*) 'a dir = |', dir(1:inp_strlen(dir)),'|'
c
      call util_file_prefix(stub, tmp)
      ltmp = inp_strlen(tmp)
      ldir  = inp_strlen(dir)
      if (ltmp+ldir+1 .gt. len(name)) then
         write(LuOut,*) ' util_file_name: stub = ', stub
         write(LuOut,*) ' util_file_name: ltmp, ldir, lname',
     $        ltmp, ldir, len(name)
         call util_flush(LuOut)
         call errquit('util_file_name: name too small', ltmp+ldir+1,
     &       INPUT_ERR)
      endif
      if (dir .ne. ' ') then
         name = dir
         name(ldir+1:ldir+1) = '/'
         name(ldir+2:) = tmp
      else
         name = tmp
      endif
c
      if (oparallel) then
         if (inp_strlen(name) .gt. len(tmp)) then
            write(LuOut,*) ' util_file_name: name = ', name
            call util_flush(LuOut)
            call errquit('util_file_name: tmp too small',
     $           inp_strlen(name), INPUT_ERR)
         endif
         tmp = name
         if(nodgs.eq.-1) then
            call util_pname(tmp, name)
         else
            call util_pname0(tmp, name,ga_nodeid(),10**nodgs)
         endif
            
      endif
c
      end
      subroutine util_file_print_dirs()
      implicit none
#include "util.fh"
#include "inp.fh"      
#include "cfileprefix.fh"
#include "stdio.fh"
#include "global.fh"
c
c     Print a summary of the permanent and scratch file directories
c
      character*(nw_max_path_len) sdir, pdir, prevsdir, prevpdir
      integer node
      logical util_file_parse_dir
      external util_file_parse_dir
c
      prevpdir = ' '
      prevsdir = ' '
      if (ga_nodeid().eq.0) then
         do node = 0, ga_nnodes()-1
c
            call util_directory_name(sdir, .true., node)
            call util_directory_name(pdir, .false., node)
c
            if ((pdir.ne.prevpdir .or. sdir.ne.prevsdir)) then
               write(LuOut,1) node, pdir(1:inp_strlen(pdir)), node,
     $              sdir(1:inp_strlen(sdir))
 1             format(i3,' permanent = ', a/
     $              i3,' scratch   = ', a)
            endif
            prevpdir = pdir
            prevsdir = sdir
         enddo
         write(LuOut,*)
         call util_flush(LuOut)
      endif
c
      end
      logical function util_file_parse_dir(dirlist, dir, nodeid)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) dirlist     ! [input] List of dirs (by host/proc)
      character*(*) dir         ! [output] Returns matching dir
c
c     Dirlist is the input line from the scratch_dir/permanent_dir
c     directives ... attempt to find a match for process nodeid.
c
c     Return true if a match was found, or false (with dir=' ')
c
      character*256 hostname
      character*1024 default, specific
      character*1 numbers(10)
      integer nodeid
c
      integer istart, iend, i, ind, p, ihostend, def1, spe1, nspe, ndef
      data numbers /'0','1','2','3','4','5','6','7','8','9'/
c
*      write(LuOut,*) nodeid, ' dirlist in ufpd |',dirlist,'|'
c
#if defined(CYGNUS) || defined(WIN32)
      call fix_windows_path(dirlist)
#endif
      util_file_parse_dir = .true.
      ndef = 0                  ! No. of default dirs found
      nspe = 0                  ! No. of host specific dirs found
      def1 = 1                  ! Pointer to end of default list 
      spe1 = 1                  ! Pointer to end of specific list 
      default = ' '
      specific = ' '
      hostname = ' '
c
c     Go thru and assemble a space separated list of directories that
c     are either specific to this host or are defaults for all processes.
c     If we encounter a process specific directory immediately return.
c
      istart = 0
 10   if (inp_strtok(dirlist, ' ', istart, iend)) then ! While loop
c
c     Check if there is a host/id present ... a colon in the token
c     On Windows platforms this could be part of a legitimate local
c     path name, so we let it through - note that this precludes
c     specification of hosts this way here
#if !defined(CYGNUS) && !defined(WIN32)
         do i = istart, iend
            if (dirlist(i:i) .eq. ':') goto 20
         enddo
#endif
         ndef = ndef + 1
         default(def1:) = dirlist(istart:iend) ! No colon=add to default list
         def1 = def1 + iend - istart + 2
         goto 10
 20      ihostend = i-1
         if (ihostend .lt. istart) call errquit
     $        ('util_dir_parse: colon at start of dirname?',0,
     &       INPUT_ERR)
c
c     Found host/process ID in dirlist(istart:ihostend).  If the
c     first character is a number, then assume it's all a process id
c     
         if (inp_match(10, .true., dirlist(istart:istart),
     $        numbers, ind)) then
            p = 0
            do i = istart, ihostend
               if (.not. inp_match(10, .true., dirlist(i:i),
     $              numbers, ind)) call errquit
     $              ('util_dir_parse: bad character in process id',0,
     &       INPUT_ERR)
               p = p*10 + ind - 1
            enddo
*            write(LuOut,*) ' p ',p
            if (p .eq. nodeid) then ! Dir for me and only me
               dir = dirlist(ihostend+2:iend)
               return
            endif
         else
            if (hostname .eq. ' ') call util_hostname(hostname)
*            write(LuOut,*) ' hostname ', hostname
*            write(LuOut,*) ' ........ ', dirlist(istart:ihostend)
            if (inp_compare(.false., hostname,
     $           dirlist(istart:ihostend))) then
               specific(spe1:) = dirlist(ihostend+2:iend)
               nspe = nspe + 1
               spe1 = spe1 + iend - ihostend-2 + 2
*               write(LuOut,*) ' set spe to |',specific,'|'
            endif
         endif
         goto 10
      endif                     ! End of while
c
c     Round robin allocation from either host specific or default lists.
c     Exploit sequential number of processes on a given host.
c
      if (nspe .gt. 0) then
         istart = 0
         do i = 0, mod(nodeid,nspe)
            if (.not. inp_strtok(specific, ' ', istart, iend))
     $           call errquit('util_file_parse_dir: internal err?',0,
     &       INPUT_ERR)
         enddo
         dir = specific(istart:iend)
         return
      endif
      if (ndef .gt. 0) then
         istart = 0
         do i = 0, mod(nodeid,ndef)
            if (.not. inp_strtok(default, ' ', istart, iend))
     $           call errquit('util_file_parse_dir: internal err?',1,
     &       INPUT_ERR)
         enddo
         dir = default(istart:iend)
         return
      endif
c
c     Nothing matched
c
      dir = ' '
      util_file_parse_dir = .false.
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine util_set_default_scratch_dir(scratch_dir)
      implicit none
c
c     Sets the default scratch directory name (this may yet be 
c     overriden by the name specified in the input). 
c
c     The order of precedence is as follows:
c     1) the value of NWCHEM_SCRATCH_DIR environment variable
c     2) the value of scratch_dir key in nwchemrc file(s)
c     3) the value of the compiled in name
c
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) scratch_dir ! output
c
      logical from_environment
      logical from_nwchemrc
      logical from_compile
      logical debug
c
      debug = .false.
      from_nwchemrc = .false.
c
c     1: check for NWCHEM_SCRATCH_DIR environment variable
c
      call util_getenv('NWCHEM_SCRATCH_DIR',scratch_dir)
      if (debug) then
        write(luout,*)
     &      'env return value of NWCHEM_SCRATCH_DIR <',
     &      scratch_dir(1:inp_strlen(scratch_dir)),'>'
      endif
      from_environment = (inp_strlen(scratch_dir).gt.0)
c
c     2: check for scratch_dir defined in nwchemrc config file(s)
c
      if (.not.from_environment) then
        if (.not.util_nwchemrc_get('scratch_dir',scratch_dir)) then
          if (debug) then
            write(luout,*)'util_nwchemrc_get failed for scratch_dir'
          endif
        else
          from_nwchemrc = .true.
          if (debug) then
            write(luout,*)
     &          'nwchemrc return value of scratch_dir <',
     &          scratch_dir(1:inp_strlen(scratch_dir)),'>'
          endif
        endif
      endif
c
c     3: use compiled in default setting
c
      if (.not.from_environment.and..not.from_nwchemrc) then
        scratch_dir = ' '
        from_compile = .true.
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine util_set_default_permanent_dir(permanent_dir)
      implicit none
c
c     Sets the default permanent directory name (this may yet be 
c     overriden by the name specified in the input). 
c
c     The order of precedence is as follows:
c     1) the value of NWCHEM_PERMANENT_DIR environment variable
c     2) the value of permanent_dir key in nwchemrc file(s)
c     3) the value of the compiled in name
c
#include "stdio.fh"
#include "inp.fh"
#include "util.fh"
c
      character*(*) permanent_dir ! output
c
      logical from_environment
      logical from_nwchemrc
      logical from_compile
      logical debug
c
      debug = .false.
      from_nwchemrc = .false.
      from_environment = .false.
c
c     1: check for NWCHEM_PERMANENT_DIR environment variable
c
      call util_getenv('NWCHEM_PERMANENT_DIR',permanent_dir)
      if (debug) then
        write(luout,*)
     &      'env return value of NWCHEM_PERMANENT_DIR <',
     &      permanent_dir(1:inp_strlen(permanent_dir)),'>'
      endif
      from_environment = (inp_strlen(permanent_dir).gt.0)
c
c     2: check for scratch_dir defined in nwchemrc config file(s)
c
      if (.not.from_environment) then
        if (.not.util_nwchemrc_get('permanent_dir',permanent_dir)) then
          if (debug) then
            write(luout,*)'util_nwchemrc_get failed for permanent_dir'
          endif
        else
          from_nwchemrc = .true.
          if (debug) then
            write(luout,*)
     &          'nwchemrc return value of permanent_dir <',
     &          permanent_dir(1:inp_strlen(permanent_dir)),'>'
          endif
        endif
      endif
c
c     3: use compiled in default setting
c
      if (.not.from_environment.and..not.from_nwchemrc) then
        permanent_dir = ' '
        from_compile = .true.
      endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine input_file_info(input_filename, 
     $     rtdb_name, ostartup, ocontinue)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "cfileprefix.fh"
#include "util.fh"
#include "stdio.fh"
      character*(*) input_filename ! [input]
      character*(*) rtdb_name   ! [output]
      logical ostartup           ! [output]
      logical ocontinue          ! [output]
c
      character*(nw_max_path_len) ecce_file_name
      character*(nw_max_path_len) a_temporary_file
      logical status, odirective, echo
      logical bad_permanent_dir, bad_scratch_dir
      logical already
      integer nkeys, istart, iend
      parameter (nkeys = 7)
      logical iocheckk
      integer mitob1
      character*16 keys(nkeys), field
      data keys/'start','restart','continue',
     &    'scratch_dir','permanent_dir', 'ecce_print',
     $     'echo'/
c
c     Scan the input for start/restart directives and attempt
c     to figure out the name of the desired data base, if the
c     job is a startup or a restart, what the file_prefix is.
c
c     While we're doing this also scan for scratch_dir and permanent_dir
c
c     (start || restart) [<file_prefix> = 'from input file base'] \
c                        [rtdb <rtdb_file_name>]
c
c     scratch_dir <read rest of line as character string>
c     permanent_dir <read rest of line as character string>
c
c     Only process 0 reads ... everyone else jumps to the broadcast
c
      mitob1=MA_sizeof(MT_INT,1,MT_BYTE)
      scratch_dir = ' '
      permanent_dir = ' '
      call util_set_default_scratch_dir(scratch_dir)
      call util_set_default_permanent_dir(permanent_dir)
      iocheckk=.true.
#if defined(NOFSCHECK) || defined(CRAYXT) || defined(BGP) || defined(BGQ)
c     on catamount all fs are parallel, therefore we need only node0
c     on BGP/BGQ all fs are parallel, therefore we need only node0
      iocheckk=ga_nodeid().eq.0
#endif     
c
      rtdb_name  = ' '
      
      if (ga_nodeid() .gt. 0) goto 10000
c
c     default is a startup with name extracted from that of the input
c     file unless a database of that name is present in which case
c     you get a restart.  Overriden by presenting start/restart.
c
      odirective = .false.      ! True if find a start/restart/continue
      ostartup   = .true.
      ocontinue  = .false.
      echo       = .false.
      call input_default_file_prefix(input_filename,file_prefix)
c
      rewind LuIn
      call inp_init(LuIn,LuOut)
 10   if (inp_search(.false., keys, nkeys)) then ! While
         if (.not. inp_a(field)) call errquit('input_start: inp?',0,
     &       INPUT_ERR)
         if (inp_compare(.false.,'start',field)) then
            odirective = .true.
            ostartup = .true.
            ocontinue = .false.
            call util_read_start_dir(file_prefix, rtdb_name)
         else if (inp_compare(.false.,'continue',field)) then
            call errquit('continue directive is no longer supported',
     *        555, INPUT_ERR)
            odirective = .true.
            ostartup = .false.
            ocontinue = .true.
            call util_read_start_dir(file_prefix, rtdb_name)
         else if (inp_compare(.false.,'restart',field)) then
            odirective = .true.
            ostartup = .false.
            ocontinue = .false.
            call util_read_start_dir(file_prefix, rtdb_name)
         else if (inp_compare(.false.,'scratch_dir',field)) then
            status = inp_line(scratch_dir)
            istart = 0
            status = inp_strtok(scratch_dir,' ',istart, iend)
            scratch_dir(istart:iend) = ' '
         else if (inp_compare(.false.,'permanent_dir',field)) then
            status = inp_line(permanent_dir)
            istart = 0
            status = inp_strtok(permanent_dir,' ',istart, iend)
            permanent_dir(istart:iend) = ' '
         else if (inp_compare(.false.,'echo',field)) then
            echo = .true.
         else if (inp_compare(.false.,'ecce_print', field)) then
            if (inp_a(ecce_file_name)) then
               call ecce_print_file_open(ecce_file_name)
               call ecce_print_echo_input(input_filename)
            endif
         else
            call errquit('input_start_opt: wierd error',0, INPUT_ERR)
         endif
         goto 10                ! End while
      endif
      rewind LuIn
      call inp_init(LuIn,LuOut)
c
      if (echo) call input_echo(LuIn,LuOut)
c
      if (rtdb_name .eq. ' ')
     $     call util_file_name('db',.false.,.false.,rtdb_name)
c
      if (.not. odirective) then
c
c     No start/restart directive presented.  See if we can find a database,
c     if so, assume a restart.  Otherwise it must be a startup.
c
         inquire(file=rtdb_name,exist=status)
         ostartup = .not. status
      endif
c
      rewind LuIn
      call inp_init(LuIn,LuOut)
c
c     Broadcast start options to everyone else
c
10000 call ga_brdcst(Msg_StartUp, ostartup, mitob1, 0)
      call ga_brdcst(Msg_StartUp, ocontinue, mitob1, 0)
      call util_char_ga_brdcst(Msg_startup, file_prefix, 0)
      call util_char_ga_brdcst(Msg_startup, scratch_dir, 0)
      call util_char_ga_brdcst(Msg_startup, permanent_dir, 0)
*
* now confirm (on each node) that scratch_dir and permanent_dir
* can have files
*
      call ga_sync()
*check permanent directory
      if(iocheckk) then
        call util_file_name('dir_check_p',.false.,.true.,
     &                      a_temporary_file)
        inquire(file=a_temporary_file, exist=already)
        bad_permanent_dir = .true.
        if(already) then
          iend = inp_strlen(a_temporary_file)
          write(luout,*)' Warning: test file already existed: ',
     &                  a_temporary_file(1:iend)
C         Assume all is well in the world, since file might be 
C         chmod 000
          bad_permanent_dir = .false.
        else
          open(UNIT=42,FILE=a_temporary_file,STATUS="new",ERR=91111)
          close(UNIT=42,STATUS="delete")
          call util_file_unlink(a_temporary_file)
          bad_permanent_dir = .false.
        endif
91111   continue
*check scratch directory
        call util_file_name('dir_check_s',.true.,.true.,
     &                      a_temporary_file)
        inquire(file=a_temporary_file, exist=already)
        bad_scratch_dir = .true.
        if(already) then
          iend = inp_strlen(a_temporary_file)
          write(luout,*)' Warning: test file already existed: ',
     &                  a_temporary_file(1:iend)
C         Assume all is well in the world, since file might be
C         chmod 000
          bad_scratch_dir = .false.
        else
          open(UNIT=43,FILE=a_temporary_file,STATUS="new",ERR=91122)
          close(UNIT=43,STATUS="delete")
          call util_file_unlink(a_temporary_file)
          bad_scratch_dir = .false.
        endif
91122   continue
      else
        bad_permanent_dir=.false.
        bad_scratch_dir=.false.
      endif
*
      if (bad_permanent_dir) then
        iend = inp_strlen(permanent_dir)
        write(luout,*)' could not open a file in permanent directory: ',
     &      permanent_dir(1:iend)
      endif
      if (bad_scratch_dir) then
        iend = inp_strlen(scratch_dir)
        write(luout,*)' could not open a file in scratch directory: ',
     &      scratch_dir(1:iend)
      endif
      if (bad_permanent_dir.and.bad_scratch_dir) then
        write(luout,*)
     &      ' Both permanent and scratch directory not accessible'
        call errquit('******** Fatal Error ********',911, INPUT_ERR)
      else if (bad_permanent_dir) then
        call errquit
     &      ('Fatal Error: permanent directory not accessible',911,
     &       INPUT_ERR)
      else if (bad_scratch_dir) then
        call errquit
     &      ('Fatal Error: scratch directory not accessible',911,
     &       INPUT_ERR)
      else
        return
      endif
c
      end
      subroutine input_default_file_prefix(input_file_name,file_prefix)
      implicit none
#include "inp.fh"
      character*(*) input_file_name, file_prefix
c
      integer i, start, end
c
      end = inp_strlen(input_file_name)
      do start = end,1,-1           ! Ignore any directories in the path
         if (input_file_name(start:start) .eq. '/') goto 10
      enddo
 10   start = start + 1
c
      do i = end,start,-1       ! Remove last trailing .*
         if (input_file_name(i:i) .eq. '.') then
            end = i - 1
            goto 20
         endif
      enddo
c
 20   if (end .lt. start) then
         file_prefix = 'calc'   ! Confused ... just punt
      else
         file_prefix = input_file_name(start:end)
      endif
c
      end
      subroutine util_read_start_dir(file_prefix, rtdb_name)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "util.fh"
      character*(*) file_prefix, rtdb_name 
c
c     (start|restart|continue) [<file_prefix>] [rtdb <rtdb_name>]
c
      character*(nw_max_path_len) test
c
 10   if (inp_a(test)) then
         if (inp_compare(.false.,test,'rtdb')) then
            if (.not. inp_a(rtdb_name)) call errquit
     $           ('util_read_start_directive: missing rtdb name',0,
     &       INPUT_ERR)
         else
            file_prefix = test
         endif
         goto 10
      endif
c
      end
      subroutine util_directory_name(dir, oscratch, node)
      implicit none
#include "cfileprefix.fh"
#include "inp.fh"
#include "global.fh"
      character*(*) dir         ! [output]
      logical oscratch          ! [input]
      integer node              ! [input]
      logical util_file_parse_dir
      external util_file_parse_dir
      character*1024 envscr
      integer util_getblnk
      external util_getblnk
      logical envok,util_getenvf
      external util_getenvf
c
c     return the name of the scratch/permanent directory for the
c     specified process
c
      dir='. '
      envok=.false.
      if (oscratch) then
         if(util_getenvf('SCRATCH_DIR',envscr)) then
            envok =
     U           util_file_parse_dir(envscr, dir, ga_nodeid())
         endif
         if(.not.envok)   then
            if(scratch_dir.eq.'null') call
     E           errquit('scratch_dir not set ',0,0)
            if (.not. util_file_parse_dir(scratch_dir, 
     d           dir, node)) then
               dir = '. '       ! Final default is blank
            endif
         endif
      else
         if(util_getenvf('PERMANENT_DIR',envscr)) then
            envok = 
     U           util_file_parse_dir(envscr, dir, ga_nodeid())
         endif
         if(.not.envok) then
            if(permanent_dir.eq.'null') call
     E           errquit('permanent_dir not set ',0,0)
            if (.not. util_file_parse_dir(permanent_dir, 
     D           dir, node)) then
               dir = '. '       ! Final default is blank
            endif
         endif
      endif
c
      end
      logical function util_getenvf(env_name,env_value)
      implicit none
      character*(*) env_name  ! [input] the name of the environment variable to check for
      character*(*) env_value ! [output] the value of the named environement variable
c
      integer len, status
c
      call get_environment_variable(env_name, env_value, len, status)
cdbg      write(6,*) ' value is !', trim(env_value),'!'
cdbg      write(6,*) ' len= ',len,' ; status= ',status
      util_getenvf=status.eq.0
      return
      end
      double precision function util_scratch_dir_avail_for_me()
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "util.fh"
#include "eaf.fh"
#include "cfileprefix.fh"
#include "inp.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "msgids.fh"
c
c     Return the amount of space in Kb available in the scratch
c     directory for this process.
c
c     Eventually this will be hooked up to input control.
c     Presently, it does the following.
c
c     For the IBM SP if the directory is /scratch or the same
c     as the default scratch directory, then  it assumes that a 
c     local (non-shared) scratch directory is being used.  
c     FOR ALL OTHER MACHINES it tries to determine the number of
c     processes sharing the directory by looping thru all nodes
c     and seeing if they map to the same directory.  This is valid
c     for machines with shared filesystems.
c
      integer me, nproc,  ierr
      character*(nw_max_path_len) mine
      integer nuse
      integer avail0,avail1
      integer fd
      character*8 fstype
      integer l1megabyte,i_k,l_k,nuse_fail,nattpt,
     ,     availmax
      character*255 dirscr
      parameter(l1megabyte=2*1000000)
      integer fiftytb! 50TB
      parameter(fiftytb=50*1024*1024)
      logical util_file_parse_dir,oprint,io_node
      integer ppn,toolarge
      double precision my8mega
      external util_file_parse_dir
c
      oprint=util_print('available disk',print_high)
c
c     Construct a name in the scratch directory of the current 
c     process and inquire how much space is available
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      nattpt=0
      avail0=0
      avail1=0
      toolarge=0
#ifdef NOIO
      avail0=10**8
      nuse=ga_nnodes()
      avail1=-16d0*nuse+avail0
#else
      call util_getppn(ppn)
      io_node=mod(me,ppn).eq.0
      if(io_node) then
         call util_file_name('junk',.true.,.true.,mine)
         ierr=eaf_delete(mine)
         ierr=eaf_open(mine, eaf_rw, fd)
         ierr = eaf_stat(mine, avail0, fstype)
         if(avail0.gt.fiftytb) toolarge=1
      else
         ierr=0
      endif
      call ga_igop(msg_utscr,toolarge,1,'+')
      if(io_node.and.oprint) 
     c     write(6,*) ga_nodeid(),': toolarge ',toolarge
      if(toolarge.gt.0) then
         nuse=ga_nnodes()
      else
 1025 call ga_sync()
      if(ierr.ne.0) call errquit('utilscratchavail: eaf_open',ierr,
     &     DISK_ERR)
      if (io_node) then
c
c     now write 16M of doubles
c
      if (.not.ma_push_get(MT_Dbl,l1megabyte,'cc',l_k,i_k))
     &   call errquit('utilfname: cannot allocate ',0, MA_ERR)
      call dcopy(l1megabyte,0d0,0,dbl_mb(i_k),1)
      ierr = eaf_write(fd, 0d0, dbl_mb(i_k),8*l1megabyte)
      my8mega=8*l1megabyte/1024d0/1024d0
      if (.not.ma_pop_stack(l_k))
     &   call errquit('utilfname: cannot deallocate ',0, MA_ERR)
      ierr=eaf_close(fd)
      endif
      call ga_sync()
      nuse_fail=0
      if(io_node) then
         ierr = eaf_stat(mine, avail1, fstype)
         if (ierr .ne. 0) call errquit(
     U        'util_scratch_avail: eaf_stat',ierr,DISK_ERR)
         ierr=eaf_delete(mine)
         nuse=nint((avail0-avail1)/my8mega)
         if (nuse .le. 0) then
            write(6,*) ' avail0,avail1,avail0-avail1 ',
     R           avail0,avail1,avail0-avail1
            nuse_fail= 1
         endif
      else
         nuse=0
      endif
      call ga_igop(msg_utscr,nuse_fail,1,'+')
      if(nuse_fail.ne.0) then
        if(oprint.and.ga_nodeid().eq.0) write(luout,*) ' negative nuse'
        nattpt=nattpt+1
        if(nattpt.gt.3) then
           if(oprint.and.ga_nodeid().eq.0) then
              write(luout,*)'utilscravail: too many attempts',nattpt
              write(luout,*)'utilscravail:  set nuse = ga_nnodes'
           endif
c
c    something wrong with filesystem: go for upper bound, ie nuse=nproc
c
           nuse=ga_nnodes()
        else
          goto 1025
        endif
      endif
      call ga_igop(msg_utscr,nuse,1,'+')
      nuse=nuse*ppn
      endif
c     io_node.eq.f => avail0=0
      call ga_igop(msg_utscr2,avail0,1,'max')
#endif
      if(oprint) then
         call util_directory_name(dirscr, .true., ga_nodeid())
         write(luout,11) ga_nodeid(),nuse,
     .     avail0/nuse,
     .     dirscr(1:inp_strlen(dirscr))
 11      format(i4,': nuse=',i4,' avail=',i9,
     .   'Mb on ',A)
      endif
c
c     now eaf_stat returns Mb instead of kb
c
!      avail=avail0*1024

c
      availmax=avail0/nuse
      util_scratch_dir_avail_for_me = 1024d0*dble(availmax)
c
      end
      subroutine util_file_name_resolve(filename, oscratch)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
      character*(*) filename    ! [input/output] 
      logical oscratch          ! [input] 
c
c     If the given filename is not a full path (begins with /)
c     or explicitly relative to the current directory (./ or ../)
c     then resolve it to the scratch/permanent directory for the
c     current process according to oscratch.
c
c     Note that this resolution cannot happen at input time since only
c     process 0 reads the input and the directories are process specific.
cc
      character*(nw_max_path_len) dir
      integer flen, dlen
c
      if ((filename(1:1).ne.'/')  .and. (filename(1:2).ne.'./')
     $     .and. (filename(1:3).ne.'../')) then
         call util_directory_name(dir, oscratch, ga_nodeid())
c
*        write(LuOut,*) 'b fnm = |',filename(1:inp_strlen(filename)),'|'
*        write(LuOut,*) 'b dir = |',dir(1:inp_strlen(dir)),'|'
c
         dlen = inp_strlen(dir)
         if (dlen .gt. 0) then
            flen = inp_strlen(filename)
            if ((flen+dlen+1).gt.len(filename)) call errquit
     $           ('util_file_name_resolve: filename too small',
     $           flen+dlen+1, INPUT_ERR)
            dir(dlen+1:dlen+1) = '/'
            dir(dlen+2:) = filename
c     
*            write(6,*) ' RESOLVED ', filename(1:flen), ' TO ', 
*     $           dir(1:inp_strlen(dir))
c     
            filename = dir
         endif
      endif
c
      end
      subroutine fix_windows_path(path)
      implicit none
#include "inp.fh"
      character*(*) path        ! [input/output] 
#if defined(CYGNUS) || defined(WIN32)
      integer i, l
c     On Windows platforms a backslash is a valid directory separator.
c     Replace backslashes with forward slashes so these pathnames are
c     accepted.
      l = inp_strlen(path)
      do i = 1, l
#if defined(CYGNUS)
         if (path(i:i) .eq. '\\') then
#elif defined(WIN32)
         if (path(i:i) .eq. '\') then
#endif
            path(i:i) = '/'
         endif
      enddo
#endif
      return
      end

c     **** added by EJB 11/7/00 ****
      subroutine util_file_name_noprefix(stub, oscratch, 
     >                                         oparallel, 
     >                                         name)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"      
#include "cfileprefix.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) stub      ! [input] stub name for file
      logical oscratch        ! [input] true=scratch, false=permanent
      logical oparallel       ! [input] true=append .nodeid
      character*(*) name      ! [output] full filename
c
      character*(nw_max_path_len) dir, tmp
      integer ltmp, ldir, me
      logical util_file_parse_dir
      external util_file_parse_dir
c
      me = ga_nodeid()
c
      call util_directory_name(dir, oscratch, me)
c
*     write(LuOut,*) 'c stub= |',stub(1:inp_strlen(stub)),'|'
*     write(LuOut,*) 'c dir = |', dir(1:inp_strlen(dir)),'|'
c

      tmp   = stub
      ltmp  = inp_strlen(tmp)
      ldir  = inp_strlen(dir)
      if (ltmp+ldir+1 .gt. len(name)) then
         write(LuOut,*) ' util_file_name: stub = ', stub
         write(LuOut,*) ' util_file_name: ltmp, ldir, lname',
     $        ltmp, ldir, len(name)
         call util_flush(LuOut)
         call errquit('util_file_name: name too small', ltmp+ldir+1,
     &       INPUT_ERR)
      endif
      if (dir .ne. ' ') then
         name = dir
         name(ldir+1:ldir+1) = '/'
         name(ldir+2:) = tmp
      else
         name = tmp
      endif
c
      if (oparallel) then
         if (inp_strlen(name) .gt. len(tmp)) then
            write(LuOut,*) ' util_file_name: name = ', name
            call util_flush(LuOut)
            call errquit('util_file_name: tmp too small',
     $           inp_strlen(name), INPUT_ERR)
         endif
         tmp = name
         call util_pname(tmp, name)
      endif
c
      end
      logical function util_find_dir(dname)
      implicit none
c
c     hack that returns true if directory exists
c
#include "eaf.fh"
#include "inp.fh"
      character*(*) dname
c
      integer availkb
      character*20 fstype
c
      util_find_dir=eaf_stat(dname(1:inp_strlen(dname)),
     .     availkb, fstype).eq.0
      return
      end
c
      subroutine util_full_file_name(filename, oscratch, name)
      implicit none
#include "errquit.fh"
#include "util.fh"
#include "inp.fh"
#include "cfileprefix.fh"
#include "global.fh"
#include "stdio.fh"
c
      character*(*) filename          ! [input] raw file name
      logical oscratch                ! [input] true=scratch, false=permanent
      character*(*) name              ! [output] full filename including scratch path
c
      character*(nw_max_path_len) dir, tmp
      integer ltmp, ldir, me
      logical util_file_parse_dir
      external util_file_parse_dir
c
      me = ga_nodeid()
c
      call util_directory_name(dir, oscratch, me)
c
*     write(LuOut,*) 'd fnm = |',filename(1:inp_strlen(filename)),'|'
*     write(LuOut,*) 'd dir = |', dir(1:inp_strlen(dir)),'|'
c
      tmp   = filename
      ltmp  = inp_strlen(tmp)
      ldir  = inp_strlen(dir)
      if (ltmp+ldir+1 .gt. len(name)) then
        write(LuOut,*) ' util_full_file_name: filename = ', filename
        write(LuOut,*) ' util_full_file_name: ltmp, ldir, lname',
     $        ltmp, ldir, len(name)
        call util_flush(LuOut)
        call errquit('util_full_file_name: name too small', ltmp+ldir+1,
     &       INPUT_ERR)
      endif
      if (dir .ne. ' ') then
         name = dir
         name(ldir+1:ldir+1) = '/'
         name(ldir+2:) = tmp
      else
         name = tmp
      endif
c
      end
      subroutine cphf_fname(cphf_str1,cphf_str2)
      implicit none
#include "inp.fh"
#include "global.fh"
      character*(*) cphf_str1 ! [in]
      character*(*) cphf_str2 ! [out]
c     file is managed serially (node 0) on perm_dir
c
      integer str1len
c     
      str1len=inp_strlen(cphf_str1)
      call util_file_name0(cphf_str1(1:str1len),
     L     .false.,.false.,
     S     cphf_str2,1)
      return
      end
      subroutine cphf_fname_parallel(cphf_str1,cphf_str2)
      implicit none
#include "inp.fh"
#include "global.fh"
      character*(*) cphf_str1 ! [in]
      character*(*) cphf_str2 ! [out]
c     file is managed in parallel on scratch_dir
c
      integer str1len
c     
      str1len=inp_strlen(cphf_str1)
      call util_file_name(cphf_str1(1:str1len),
c      scratch_dir, parallel
     L     .true.,.true.,
     S     cphf_str2)
      write(6,*) ga_nodeid(),' fparal ',cphf_str1(1:str1len),
     A     cphf_str2(1:inp_strlen(cphf_str2))
      return
      end
