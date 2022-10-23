      subroutine bummer( text, ierr, errtyp )
c
c  process a program error.
c  input:
c  text  = character string to be printed.
c  ierr  = internal program error to be printed.
c  errtyp = 0 for warning.  traceback may be generated. execution
c             continues.
c         = 1 for nonfatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to allow
c             subsequent program steps to continue if possible.
c         = 2 for fatal error.  traceback may be generated.
c             execution is stopped. jcl condition code is set to abort
c             subsequent program steps if possible.
c         = 3 only writing the bummer file (normal termination)
c
c  entry ibummr must be called prior to bummer() to set the output
c  unit and to perform any additional initialization.
c
c  version log:
c  13-oct-94 remove 'call exit' for ibm rs6000. (ahhc)
c  24-apr-92 %val() blocks added for ibm rs6000. -rls
c  11-sep-91 parerr() calls added. -rjh
c  13-mar-91 posix code added. -rls
c  01-mar-89 write(stderr,*) and call exit() for unix machines (rls).
c  05-jul-88 unicos version (rls).
c  03-nov-87 ibm version (dcc).
c  10-sep-87 written by ron shepard.
c
      implicit integer(a-z)
c
      character*(*) text
      integer iunit, ierr, errtyp
c
c     # bummer error types.
c     wrnerr = warning message
c     nfterr = non fatal error
c     faterr = fatal error
c     termerr = writes "normally terminated" onto file bummer
      integer   wrnerr,  nfterr,  faterr, termerr
      parameter(wrnerr=0,nfterr=1,faterr=2, termerr=3)
c
      integer pxferr
      integer nlist,        stderr
      save    nlist,        stderr
      data    nlist / 6 /,  stderr / 0 /
c
c*******************************************************************
      if ( errtyp .eq. wrnerr ) then
c*******************************************************************
c
c        # print a warning message and continue execution.
c
         open(unit=99,file='bummer')
         write(99,*) 'warning error encountered'
         close(99)

         write(nlist,6010) 'bummer (warning):', text, ierr
#if defined(POSIX)
         write(stderr,6010) 'bummer (warning):', text, ierr
#elif defined(UNICOS)
         write(stderr,6010) 'bummer (warning):', text, ierr
c        call tracebk( nlist )
#elif defined(UNIX)
         write(stderr,6010) 'bummer (warning):', text, ierr
#endif 
         return
c
c*******************************************************************
      elseif ( errtyp .eq. nfterr ) then
c*******************************************************************
c
c        # print a warning message, stop execution.
c
         open(unit=99,file='bummer')
         write(99,*) 'nonfatal error encountered'
         close(99)

         write(nlist,6010) 'bummer (nonfatal):', text, ierr
#if defined(POSIX)
         write(stderr,6010) 'bummer (nonfatal):', text, ierr
c        call parerr( ierr )
         call pxfexit( 0 )
#elif defined(UNICOS)
         write(stderr,6010) 'bummer (nonfatal):', text, ierr
         call tracebk( nlist )
c        call parerr( ierr )
         stop 'program error'
#elif defined(HP) || defined(RS6000) || defined(UNIX)
         write(stderr,6010) 'bummer (nonfatal):', text, ierr
c        call parerr( ierr )
         stop 'program error'
#elif defined(VAX)
         call sys$exit(%val(42))
#elif defined(IBM)
         stop 901
#else 
c        call parerr( ierr )
         stop 'program error'
#endif 
c
c**********************************************************************
      elseif ( errtyp .eq. faterr ) then
c**********************************************************************
c
c        # print an error message, stop execution, and abort job
c        # sequence.
c
         write(nlist,6010) 'bummer (fatal):', text, ierr
         open(unit=99,file='bummer')
         write(99,*) 'fatal error encountered'
         close(99)

#if defined(POSIX)
         write(stderr,6010) 'bummer (fatal):', text, ierr
c        call parerr( ierr )
         call pxfexit( 1 )
#elif defined(UNICOS)
         write(stderr,6010) 'bummer (fatal):', text, ierr
c        call parerr( ierr )
c        # abort() generates tracebacks automatically.
         call abort
#elif defined(HP) || defined(RS6000)
         write(stderr,6010) 'bummer (fatal):', text, ierr
c        call parerr( ierr )
c        # abort() generates tracebacks automatically.
         stop 'program error'
#elif defined(UNIX)
         write(stderr,6010) 'bummer (fatal):', text, ierr
c        call parerr( ierr )
         call exit( 1 )
         stop 'bummer'
#elif defined(VAX)
         call sys$exit( %val(44) )
#elif defined(IBM)
         stop 902
#else 
         stop 'program error'
#endif 
c
c*******************************************************************
      elseif (errtyp .eq. termerr) then
c*******************************************************************
          open(unit=99,file='bummer')
          write(99,*) 'normally terminated'
          close (99)
          return
c*******************************************************************
      else
c*******************************************************************
c
c        # unknown error level.  treat as a fatal error.
c
         write(nlist,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
         open(unit=99,file='bummer')
         write(99,*) 'unknown error encountered'
         close(99)

#if defined(POSIX)
         write(stderr,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
c        call parerr( ierr )
         call pxfexit( 1 )
#elif defined(HP) || defined(RS6000) || defined(UNIX)
         write(stderr,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
c        call parerr( ierr )
c        # abort() generates tracebacks automatically.
         stop
#elif defined(UNICOS)
         write(stderr,6020) 'bummer (unknown): errtyp=',
     &    errtyp, text, ierr
c        call parerr( ierr )
c        # abort() generates tracebacks automatically.
         call abort
#elif defined(VAX)
         call sys$exit( %val(44) )
#elif defined(IBM)
         stop 903
#else 
         stop 'program error'
#endif 
c*******************************************************************
      endif
c*******************************************************************
c
c     # this statement is not executed, it is included
c     # just to avoid compiler warnings. -rls
      stop 'bummer error'
c
c     # initialization...
c
      entry ibummr( iunit )
c
c     # save the listing unit for use later.
c
      nlist = iunit
c
#if defined(POSIX)
c     # set stderr to the correct value.
      call pxfconst( 'stderr', stderr, pxferr )
      if ( pxferr .ne. 0 ) then
         write(*,6010) 'ibummr pxfconst() error=', pxferr
         call pxfexit( 1 )
      endif
#endif 
c
      return
6010  format(1x,a,a,i10)
6020  format(1x,a,i10,a,i10)
      end
