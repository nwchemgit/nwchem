      block data bfileprefix
      implicit none
#include "cfileprefix.fh"
      data file_prefix /'calc'/
      data scratch_dir/' '/, permanent_dir/' '/
c
      end
      subroutine util_file_info_rtdb(rtdb)
      implicit none
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
     $     call errquit('rtdb_put of file_prefix failed', 0)
c
      if (scratch_dir .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'scratch_dir', 1, scratch_dir))
     $        call errquit('rtdb_put of scratch_dir failed', 0)
      else
         if (.not. rtdb_cget(rtdb, 'scratch_dir', 1, scratch_dir))
     $        scratch_dir = ' '
      endif
c
      if (permanent_dir .ne. ' ') then
         if (.not. rtdb_cput(rtdb, 'permanent_dir', 1, permanent_dir))
     $        call errquit('rtdb_put of permanent_dir failed', 0)
      else
         if (.not. rtdb_cget(rtdb, 'permanent_dir', 1, permanent_dir))
     $        permanent_dir = ' '
      endif
c
      end
      subroutine util_file_prefix(name, fullname)
      implicit none
#include "cfileprefix.fh"
#include "inp.fh"
      character*(*) name, fullname
      external bfileprefix  ! This for T3D linker
c
c     prepend the file_prefix onto name as <file_prefix>.name
c     returning the result in fullname.
c
      if ((inp_strlen(name)+inp_strlen(file_prefix)+1) .gt.
     $     len(fullname)) then
         write(6,*) ' file_prefix: name = ', name
         write(6,*) ' file_prefix: prfx = ', file_prefix
         call util_flush(6)
         call errquit('file_prefix: insufficient space ', len(fullname))
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
#include "util.fh"
#include "inp.fh"      
#include "cfileprefix.fh"
#include "global.fh"
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
*      write(6,*) ' dir = |', dir,'|'
c
      call util_file_prefix(stub, tmp)
      ltmp = inp_strlen(tmp)
      ldir  = inp_strlen(dir)
      if (ltmp+ldir+1 .gt. len(name)) then
         write(6,*) ' util_file_name: stub = ', stub
         write(6,*) ' util_file_name: ltmp, ldir, lname',
     $        ltmp, ldir, len(name)
         call util_flush(6)
         call errquit('util_file_name: name too small', ltmp+ldir+1)
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
            write(6,*) ' util_file_name: name = ', name
            call util_flush(6)
            call errquit('util_file_name: tmp too small',
     $           inp_strlen(name))
         endif
         tmp = name
         call util_pname(tmp, name)
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
               write(6,1) node, pdir(1:inp_strlen(pdir)), node,
     $              sdir(1:inp_strlen(sdir))
 1             format(i3,' permanent = ', a/
     $              i3,' scratch   = ', a)
            endif
            prevpdir = pdir
            prevsdir = sdir
         enddo
         write(6,*)
         call util_flush(6)
      endif
c
      end
      logical function util_file_parse_dir(dirlist, dir, nodeid)
      implicit none
#include "inp.fh"
#include "util.fh"
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
*      write(6,*) nodeid, ' dirlist in ufpd |',dirlist,'|'
c
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
c     
         do i = istart, iend
            if (dirlist(i:i) .eq. ':') goto 20
         enddo
         ndef = ndef + 1
         default(def1:) = dirlist(istart:iend) ! No colon=add to default list
         def1 = def1 + iend - istart + 2
         goto 10
 20      ihostend = i-1
         if (ihostend .lt. istart) call errquit
     $        ('util_dir_parse: colon at start of dirname?',0)
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
     $              ('util_dir_parse: bad character in process id',0)
               p = p*10 + ind - 1
            enddo
*            write(6,*) ' p ',p
            if (p .eq. nodeid) then ! Dir for me and only me
               dir = dirlist(ihostend+2:iend)
               return
            endif
         else
            if (hostname .eq. ' ') call util_hostname(hostname)
*            write(6,*) ' hostname ', hostname
*            write(6,*) ' ........ ', dirlist(istart:ihostend)
            if (inp_compare(.false., hostname,
     $           dirlist(istart:ihostend))) then
               specific(spe1:) = dirlist(ihostend+2:iend)
               nspe = nspe + 1
               spe1 = spe1 + iend - ihostend-2 + 2
*               write(6,*) ' set spe to |',specific,'|'
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
     $           call errquit('util_file_parse_dir: internal err?',0)
         enddo
         dir = specific(istart:iend)
         return
      endif
      if (ndef .gt. 0) then
         istart = 0
         do i = 0, mod(nodeid,ndef)
            if (.not. inp_strtok(default, ' ', istart, iend))
     $           call errquit('util_file_parse_dir: internal err?',1)
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
      subroutine input_file_info(input_filename, 
     $     rtdb_name, ostartup, ocontinue)
C$Id: util_file_name.F,v 1.13 1999-06-30 22:48:36 d3e129 Exp $
      implicit none
#include "inp.fh"
#include "global.fh"
#include "tcgmsg.fh"
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
      integer nkeys, istart, iend
      parameter (nkeys = 7)
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
      scratch_dir = ' '
      permanent_dir = ' '
c
      if (nodeid() .gt. 0) goto 10000
c
c     default is a startup with name extracted from that of the input
c     file unless a database of that name is present in which case
c     you get a restart.  Overriden by presenting start/restart.
c
      odirective = .false.      ! True if find a start/restart/continue
      ostartup   = .true.
      ocontinue  = .false.
      echo       = .false.
      rtdb_name  = ' '
      call input_default_file_prefix(input_filename,file_prefix)
c
      rewind 5
      call inp_init(5,6)
 10   if (inp_search(.false., keys, nkeys)) then ! While
         if (.not. inp_a(field)) call errquit('input_start: inp?',0)
         if (inp_compare(.false.,'start',field)) then
            odirective = .true.
            ostartup = .true.
            ocontinue = .false.
            call util_read_start_dir(file_prefix, rtdb_name)
         else if (inp_compare(.false.,'continue',field)) then
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
            call errquit('input_start_opt: wierd error',0)
         endif
         goto 10                ! End while
      endif
      rewind 5
      call inp_init(5,6)
c
      if (echo) call input_echo(5,6)
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
      rewind 5
      call inp_init(5,6)
c
c     Broadcast start options to everyone else
c
10000 call ga_brdcst(Msg_StartUp, ostartup, mitob(1), 0)
      call ga_brdcst(Msg_StartUp, ocontinue, mitob(1), 0)
      call util_char_ga_brdcst(Msg_startup, file_prefix, 0)
      call util_char_ga_brdcst(Msg_startup, scratch_dir, 0)
      call util_char_ga_brdcst(Msg_startup, permanent_dir, 0)
*
* now confirm (on each node) that scratch_dir and permanent_dir can have files
*
      call ga_sync()
*check permanent directory
      call util_file_name('dir_check_p',.false.,.true.,a_temporary_file)
      open(unit=42,file=a_temporary_file,status='new',err=91101)
      close(unit=42,status='delete')
      call util_file_unlink(a_temporary_file)
      bad_permanent_dir = .false.
      goto 91111
91101 bad_permanent_dir = .true.
91111 continue
*check scratch directory
      call util_file_name('dir_check_s',.true.,.true.,a_temporary_file)
      open(unit=43,file=a_temporary_file,status='new',err=91102)
      close(unit=43,status='delete')
      call util_file_unlink(a_temporary_file)
      bad_scratch_dir = .false.
      goto 91122
91102 bad_scratch_dir = .true.
91122 continue
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
        call errquit('******** Fatal Error ********',911)
      else if (bad_permanent_dir) then
        call errquit
     &      ('Fatal Error: permanent directory not accessible',911)
      else if (bad_scratch_dir) then
        call errquit
     &      ('Fatal Error: scratch directory not accessible',911)
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
     $           ('util_read_start_directive: missing rtdb name',0)
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
      character*(*) dir         ! [output]
      logical oscratch          ! [input]
      integer node              ! [input]
      logical util_file_parse_dir
      external util_file_parse_dir
c
c     return the name of the scratch/permanent directory for the
c     specified process
c
      if (oscratch) then
         if (.not. util_file_parse_dir(scratch_dir, dir, node)) then
            if (.not. util_file_parse_dir(
     $ SCRATCH_DEF_DIR
     $           , dir, node)) dir = ' ' ! Final default is blank
         endif
      else
         if (.not. util_file_parse_dir(permanent_dir, dir, node)) then
            if (.not. util_file_parse_dir(
     $ PERM_DEF_DIR
     $           , dir, node)) dir = ' ' ! Final default is blank
         endif
      endif
c
      end
      integer function util_scratch_dir_avail_for_me()
      implicit none
#include "global.fh"
#include "util.fh"
#include "eaf.fh"
#include "cfileprefix.fh"
#include "inp.fh"
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
      integer me, nproc, you, ierr
      character*(nw_max_path_len) mine, yours
      integer avail, nuse
      character*8 fstype
      logical util_file_parse_dir
      external util_file_parse_dir
c
c     Construct a name in the scratch directory of the current 
c     process and inquire how much space is available
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      call util_file_name('junk',.true.,.true.,mine)
      ierr = eaf_stat(mine, avail, fstype)
      if (ierr .ne. 0) call errquit('util_scratch_avail: eaf_stat',ierr)
c
c     Now get the name of my scratch directory
c
      call util_directory_name(mine, .true., me)
c
      nuse = 0
c
#if defined(SP) || defined(SP1) || defined(SP2)
      if (inp_compare(.true.,mine,'/scratch') .or.
     $     inp_compare(.true.,mine,
     $ SCRATCH_DEF_DIR
     $     )) then
         nuse = 1
      endif
#endif
c
      if (nuse .eq. 0) then
c     
c     Count other processes using this directory
c     
         do you = 0, nproc-1
            call util_directory_name(yours, .true., you)
            if (yours .eq. mine) nuse = nuse + 1
         enddo
      endif
c     
      util_scratch_dir_avail_for_me = avail / nuse
c
      end
      subroutine util_file_name_resolve(filename, oscratch)
      implicit none
#include "inp.fh"
#include "util.fh"
#include "global.fh"
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
         dlen = inp_strlen(dir)
         if (dlen .gt. 0) then
            flen = inp_strlen(filename)
            if ((flen+dlen+1).gt.len(filename)) call errquit
     $           ('util_file_name_resolve: filename too small',
     $           flen+dlen+1)
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
