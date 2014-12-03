
cC     $Id$
c

      subroutine nwpw_vectors_input(rtdb, module)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "bafdecls.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the vectors directive which specifies input/output
c     of MO vectors.  Insert options into the database prefixed
c     with module:...
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
      parameter (num_dirs = 10)
      character*12 dirs(num_dirs)
      integer max_swap
      parameter (max_swap = 20)
      integer swap(max_swap), nswap, itest, mlen
      character*255 test, vectors_in, vectors_out
      character*255 vvectors_in, vvectors_out
      integer maxfrag
      parameter (maxfrag = 60)
      character*255 fragments(maxfrag)
      integer nfrag
      character*80 project(2)
      character*50 swapname
      logical olock
      data dirs /'input', 'output', 'atomic', 'hcore', 'swap',
     $     'project', 'lock', 'fragment','vinput','voutput'/
c     
      nfrag = 0  
      vectors_in = ' '
      vectors_out = ' '

      vvectors_in = ' '
      vvectors_out = ' '

      olock = .false.
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 9999
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
      goto (100,200,300,400,500,600,700,800,900,1000,
     >      9999) ind
      call errquit('nwpw_input: unimplemented directive', ind,
     &       INPUT_ERR)
c     
c     input ... quietly ignore this
c     
 100  goto 10
c     
c     output filename
c     
 200  if (.not. inp_a(test)) call errquit
     $     ('nwpw_vectors_input: vectors output filename', 0, INPUT_ERR)
      if (vectors_out .ne. ' ')  goto 10000
      vectors_out = test
      goto 10
c     
c     atomic
c     
 300  if (vectors_in .ne. ' ')  goto 10000
c     vectors_in = 'atomic'
      call util_file_prefix('movecs',vectors_in)

      goto 10
c     
c     hcore
c     
 400  if (vectors_in .ne. ' ')  goto 10000
c     vectors_in = 'hcore'
      call util_file_prefix('movecs',vectors_in)
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
     $        ('vectors_input: invalid no. of orbitals to swap', nswap,
     &       INPUT_ERR)
         if (swapname .eq. ' ') then
            swapname = ' '
            write(swapname,511) module(1:mlen), 'swap'
         else if (swapname .eq. 'alpha') then
            swapname = ' '
            write(swapname,511) module(1:mlen), 'swap alpha'
         else
            swapname = ' '
            write(swapname,511) module(1:mlen), 'swap beta'
         endif
 511     format(a,':',a)
         if (.not. rtdb_put(rtdb, swapname, mt_int, nswap,swap))
     $        call errquit('vectors_input: put of swap?', nswap,
     &       RTDB_ERR)
         goto 10
      else
         if ((nswap+2) .gt. max_swap) call errquit
     $        ('vectors_input: swapping too many orbitls', nswap,
     &       INPUT_ERR)
         swap(nswap+1) = itest
         if (.not. inp_i(swap(nswap+2))) call errquit
     $        ('vectors input: it takes two to swap!', 0, INPUT_ERR)
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
     $     call errquit('vectors_input: failed to store project',0,
     &       RTDB_ERR)
      goto 10
c     
c     Force MO locking so output matches ordering of input
c     
 700  olock = .true.
      goto 10
c
c
c     fragment frag1 frag2 ...
c
 800  continue
 801  if (inp_a(test)) then
         if (inp_match(num_dirs, .false., test, dirs, ind)) then
            call inp_prev_field()
         else
            nfrag = nfrag + 1
            if(nfrag.gt.maxfrag) call errquit(
     '         ' vec_inp: too many fragment ',maxfrag, INPUT_ERR)
            fragments(nfrag) = test
            goto 801
         endif
      endif
      if (nfrag .gt. 0) vectors_in = 'fragment'
      goto 10
c     
c     vinput filename
c     
 900  if (.not. inp_a(test)) call errquit
     $ ('nwpw_vectors_output:vvectors input filename',0,INPUT_ERR)
      if (vvectors_in .ne. ' ')  goto 10000
      vvectors_in = test
      goto 10
c
c     voutput filename
c
 1000 if (.not. inp_a(test)) call errquit
     $ ('nwpw_vectors_output:vvectors output filename',0,INPUT_ERR)
      if (vvectors_out .ne. ' ')  goto 10000
      vvectors_out = test
      goto 10

c     
 9999 continue
      if (vectors_in .ne. ' ') then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen), 'input vectors'
         if (.not. rtdb_cput(rtdb, test, 1, vectors_in)) call errquit
     $        ('vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
         write(test,'(a,'':'',a)') 'cgsd',
     $                             'input_wavefunction_filename'
         if (.not. rtdb_cput(rtdb, test, 1, vectors_in)) call errquit
     $        ('vectors_input: failed to store vectors_input', 0,
     &       RTDB_ERR)
      endif
      if (vectors_in .eq. 'fragment') then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen), 'fragments'
         if (.not. rtdb_cput(rtdb, test, nfrag, fragments)) call errquit
     $        ('vectors_input: failed to store fragment info', 0,
     &       RTDB_ERR)
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen), 'nfrag'
         if (.not. rtdb_put(rtdb, test, mt_int, 1, nfrag)) call errquit
     $        ('vectors_input: failed to store fragment info', 0,
     &       RTDB_ERR)
      endif
      if (vectors_out .ne. ' ') then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen), 'output vectors'
         if (.not. rtdb_cput(rtdb, test, 1, vectors_out)) call errquit
     $        ('vectors_out: failed to store vectors_out', 0,
     &       RTDB_ERR)
         write(test,'(a,'':'',a)') 'cgsd',
     $                             'output_wavefunction_filename'
         if (.not. rtdb_cput(rtdb, test, 1, vectors_out)) call errquit
     $        ('vectors_out: failed to store vectors_out', 0,
     &       RTDB_ERR)
      endif
      if (olock) then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen), 'movecs lock'
         if (.not. rtdb_put(rtdb, test, MT_LOG, 1, olock))
     $        call errquit('vectors: rtdb_put of olock failed',0,
     &       RTDB_ERR)
      endif

      if (vvectors_in .ne. ' ') then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen),
     >                             'input vvectors'
         if (.not. rtdb_cput(rtdb, test, 1, vvectors_in))
     >    call errquit
     >     ('vvectors_in: failed to store vvectors_in',0,
     >     RTDB_ERR)
      endif
      if (vvectors_out .ne. ' ') then
         test = ' '
         write(test,'(a,'':'',a)') module(1:mlen),
     >                             'output vvectors'
         if (.not. rtdb_cput(rtdb, test, 1, vvectors_out))
     >    call errquit
     >     ('vvectors_out: failed to store vvectors_out',0,
     >     RTDB_ERR)
      endif

c     
      return
c     
10000 write(6,10001)
10001 format(/' vectors [[input] filename|atomic|hcore] '
     $     /'         [output  filename] [lock]'
     $     /'         [swap [alpha|beta] pair_list]'/)
      call util_flush(6)
      call errquit('nwpw_vectors_input: invalid format', 0, INPUT_ERR)
c     
      end

