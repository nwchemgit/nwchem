      subroutine vectors_input(rtdb)
C$Id: vectors_inp.F,v 1.1 1996-07-26 16:24:32 d3g681 Exp $
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "context.fh"
#include "mafdecls.fh"
      integer rtdb
c
c     Parse the vectors directive which specifies input/output
c     of MO vectors.  Insert options into the database using
c     the current context so that can handle multiple tasks
c     and wavefunctions in the same input deck
c
c     Input vectors may be specified as
c
c     vectors [input] atomic
c     vectors [input] hcore
c     vectors [input] filename
c
c     with an application specific default
c
c     output vectors default to either the same filename as the
c     input vectors if they were read from a file or to the
c     filename movecs (which is searched for by the application)
c
c     vectors output filename
c
      integer num_dirs, ind
      parameter (num_dirs = 7)
      character*8 dirs(num_dirs)
      integer max_swap
      parameter (max_swap = 20)
      integer swap(max_swap), nswap, itest
      character*255 test, vectors_in, vectors_out
      character*80 project(2)
      character*20 swapname
      logical olock
      data dirs /'input', 'output', 'atomic', 'hcore', 'swap',
     $     'project', 'lock'/
c
      vectors_in = ' '
      vectors_out = ' '
      olock = .false.
c
 10   if (.not. inp_a(test)) goto 999
c
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c
c     Does not match a keyword ... assume it is the name of the
c     file with the input vectors
c     
         if (vectors_in .ne. ' ') goto 10000
         vectors_in = test
         goto 10
      endif
c
      goto (100, 200, 300, 400, 500, 600, 700, 999) ind
      call errquit('scf_input: unimplemented directive', ind)
c
c     input ... quietly ignore this
c
 100  goto 10
c
c     output filename
c
 200  if (.not. inp_a(test)) call errquit
     $     ('vectors_in: vectors output filename', 0)
      if (vectors_out .ne. ' ')  goto 10000
      vectors_out = test
      goto 10
c
c     atomic
c     
 300  if (vectors_in .ne. ' ')  goto 10000
      vectors_in = 'atomic'
      goto 10
c
c     hcore
c     
 400  if (vectors_in .ne. ' ')  goto 10000
      vectors_in = 'hcore'
      goto 10
c
c     swap [alpha|beta] <integer list of pairs>
c
 500  nswap = 0
      swapname = ' '
      if (.not. inp_i(itest)) then
         if (.not. inp_a(swapname)) goto 10000
         call inp_lcase(swapname)
         if (swapname.ne.'alpha' .and. swapname.ne.'beta') goto 10000
      else
         call inp_prev_field()
      endif
 510  if (.not. inp_i(itest)) then
c
c     Have hit end of pair list ... EOL or another directive
c
         if (mod(nswap,2).ne.0 .or. nswap.le.0) call errquit
     $        ('vectors_input: invalid no. of orbitals to swap', nswap)
         if (swapname .eq. ' ') then
            swapname = 'scf:swap'
         else if (swapname .eq. 'alpha') then
            swapname = 'scf:swap alpha'
         else
            swapname = 'scf:swap beta'
         endif
         if (.not. rtdb_put(rtdb, swapname, mt_int, nswap,swap))
     $        call errquit('vectors_input: put of swap?', nswap)
         goto 10
      else
         if ((nswap+2) .gt. max_swap) call errquit
     $        ('vectors_input: swapping too many orbitls', nswap)
         swap(nswap+1) = itest
         if (.not. inp_i(swap(nswap+2))) call errquit
     $        ('vectors input: it takes two to swap!', 0)
         nswap = nswap + 2
         goto 510
      endif
c
c     Project basisname filename
c
 600  if (.not. inp_a(project(1))) goto 10000
      if (.not. inp_a(project(2))) goto 10000
      vectors_in = 'project'
      if (.not. rtdb_cput(rtdb, 'vectors:project', 2, project))
     $     call errquit('vectors_input: failed to store project',0)
      goto 10
c
c     Force MO locking so output matches ordering of input
c
 700  olock = .true.
      goto 10
c
 999  continue
      if (vectors_in .ne. ' ') then
         if (.not. context_prefix('input vectors', test)) call errquit
     $     ('vectors_input: context_prefix failed', 0)
         if (.not. rtdb_cput(rtdb, test, 1, vectors_in)) call errquit
     $     ('vectors_input: failed to store vectors_input', 0)
      endif
      if (vectors_out .ne. ' ') then
         if (.not. context_prefix('output vectors', test)) call errquit
     $     ('vectors_out: context_prefix failed', 0)
         if (.not. rtdb_cput(rtdb, test, 1, vectors_out)) call errquit
     $     ('vectors_out: failed to store vectors_out', 0)
      endif
      if (olock) then
        if (.not. context_prefix('movecs lock', test)) call errquit
     $    ('movecs_lock: context_prefix failed', 0)
        if (.not. rtdb_put(rtdb, test, MT_LOG, 1, olock))
     $     call errquit('vectors: rtdb_put of olock failed',0)
      endif
c
      return
c
10000 write(6,10001)
10001 format(/' vectors [[input] filename|atomic|hcore] '
     $       /'         [output  filename] [lock]'
     $       /'         [swap [alpha|beta] pair_list]'/)
      call util_flush(6)
      call errquit('vectors_input: invalid format', 0)
c
      end

