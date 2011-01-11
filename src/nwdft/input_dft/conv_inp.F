      subroutine conv_input(rtdb, module)
C     $Id: 
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the convergence directive which specifies how to control
c     the iterative SCF process.  Insert options into the database 
c     prefixed with module:...
c     
c     Possible variables are:
c     
c     energy
c     density
c     gradient
c     dampon/dampoff
c     diison/diisoff
c     levlon/levloff
c     ncydp
c     ncyds
c     ncysh
c     damp/nodamping
c     diis/nodiis/nfock
c     lshift/nolevelshifting
c     pcglevl initial/tolerance/final levelshift
c     nrlevl initial/tolerance/final levelshift
c     hl_tol
c     rabuck
c     nr tolerance
c     
      integer num_dirs, ind, mlen
      parameter (num_dirs = 25)
      character*15 dirs(num_dirs)
      character*255 test
      double precision e_conv, d_conv, g_conv, dampon, dampoff, diison,
     &                 diisoff, levlon, levloff, rlshift, hl_tol,
     &                 scfshift(6),tolnr
      integer ncydp, ncyds, ncysh, ndamp, nfock,n_rabuck
      logical damp, nodamping, levelshift, nolevelshifting, diis, nodiis
      data dirs /'energy', 'density', 'gradient', 'dampon', 'dampoff',
     &           'diison', 'diisoff', 'levlon', 'levloff', 'ncydp',
     &           'ncyds', 'ncysh', 'damp', 'nodamping', 'diis',
     &           'nodiis', 'lshift', 'nolevelshifting', 'hl_tol',
     .       'fast','stable','rabuck','pcglevl','nrlevl','nr'/
c
      call dfill(6,-1.0d0,scfshift,1)
c     
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 1999
c     
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c     
c        Does not match a keyword ... 
c     
         goto 10000
      endif
c
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 
     &      1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 
     .      2000, 2010,2020,2030,2040,2050,
     &      1999) ind
      call errquit('conv_inp: unimplemented directive', ind, INPUT_ERR)
c     
c     energy; threshold for energy convergence
c
  100 if (.not. inp_f(e_conv)) then
         write(LuOut,*)' conv_inp: energy convergence value ',
     &             'not found; defaulting to 1d-6 Hartree.'
         e_conv = 1.d-6
      endif
      if (.not. rtdb_put(rtdb, 'dft:e_conv', mt_dbl, 1, e_conv))
     &   call errquit('conv_inp: rtdb_put failed', 100, RTDB_ERR)
c
      goto 10
c     
c     density; threshold for density convergence
c
  200 if (.not. inp_f(d_conv)) then
         write(LuOut,*)' conv_inp: density convergence value ',
     &             'not found; defaulting to 1d-5.'
         d_conv = 1.d-5
      endif 
      if (.not. rtdb_put(rtdb, 'dft:d_conv', mt_dbl, 1, d_conv))
     &   call errquit('conv_inp: rtdb_put failed', 200, RTDB_ERR)
c
      goto 10
c     
c     gradient; threshold for gradient convergence
c
  300 if (.not. inp_f(g_conv)) then
         write(LuOut,*)' conv_inp: gradient convergence value ',
     &             'not found; defaulting to 5d-4.'
         g_conv = 5.d-4
      endif 
      if (.not. rtdb_put(rtdb, 'dft:g_conv', mt_dbl, 1, g_conv))
     &   call errquit('conv_inp: rtdb_put failed', 300, RTDB_ERR)
c
      goto 10
c     
c     dampon; turn damping on when change in energy 
c             per iteration is less than this value
c
  400 if (.not. inp_f(dampon)) then
         write(LuOut,*)' conv_inp: dampon value not found; defaulting',
     &             ' to 0.d0 Hartree'
         dampon = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:dampon', mt_dbl, 1, dampon))
     &   call errquit('conv_inp: rtdb_put failed', 400, RTDB_ERR)
c
      goto 10
c     
c     dampoff; turn damping off when change in energy 
c              per iteration is less than this value
c
  500 if (.not. inp_f(dampoff)) then
         write(LuOut,*)' conv_inp: dampoff value not found; defaulting',
     &             ' to 0.0d0 Hartree'
         dampoff = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:dampoff', mt_dbl, 1, dampoff))
     &   call errquit('conv_inp: rtdb_put failed', 500, RTDB_ERR)
c
      goto 10
c     
c     diison; turn diis on when change in energy 
c             per iteration is less than this value
c
  600 if (.not. inp_f(diison)) then
         write(LuOut,*)' conv_inp: diison value not found; defaulting',
     &             ' to 0.d0 Hartree'
         diison = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:diison', mt_dbl, 1, diison))
     &   call errquit('conv_inp: rtdb_put failed', 600, RTDB_ERR)
c
      goto 10
c     
c     diisoff; turn diis off when change in energy 
c              per iteration is less than this value
c
  700 if (.not. inp_f(diisoff)) then
         write(LuOut,*)' conv_inp: diisoff value not found; defaulting',
     &             ' to 0.0d0 Hartree'
         diisoff = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:diisoff', mt_dbl, 1, diisoff))
     &   call errquit('conv_inp: rtdb_put failed', 700, RTDB_ERR)
c
      goto 10
c     
c     levlon; turn levelshifting on when change in energy 
c             per iteration is less than this value
c
  800 if (.not. inp_f(levlon)) then
         write(LuOut,*)' conv_inp: levlon value not found; defaulting',
     &             ' to 0.d0 Hartree'
         levlon = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:levlon', mt_dbl, 1, levlon))
     &   call errquit('conv_inp: rtdb_put failed', 800, RTDB_ERR)
c
      goto 10
c     
c     levloff; turn levelshifting off when change in energy 
c              per iteration is less than this value
c
  900 if (.not. inp_f(levloff)) then
         write(LuOut,*)' conv_inp: levloff value not found; defaulting',
     &             ' to 0.0d0 Hartree'
         levloff = 0.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:levloff', mt_dbl, 1, levloff))
     &   call errquit('conv_inp: rtdb_put failed', 900, RTDB_ERR)
c
      goto 10
c     
c     ncydp; use damping for this number of iterations
c
 1000 if (.not. inp_i(ncydp)) then
         write(LuOut,*)' conv_inp: ncydp value not found; defaulting',
     &             ' to 2 iterations'
         ncydp = 2
      endif            
      if (.not. rtdb_put(rtdb, 'dft:ncydp', mt_int, 1, ncydp))
     &   call errquit('conv_inp: rtdb_put failed', 1000, RTDB_ERR)
c
      goto 10
c     
c     ncyds; use diis for this number of iterations
c
 1100 if (.not. inp_i(ncyds)) then
         write(LuOut,*)' conv_inp: ncyds value not found; defaulting',
     &             ' to 9999 iterations'
         ncyds = 9999
      endif            
      if (.not. rtdb_put(rtdb, 'dft:ncyds', mt_int, 1, ncyds))
     &   call errquit('conv_inp: rtdb_put failed', 1100, RTDB_ERR)
c
      goto 10
c     
c     ncysh; use levelshifting for this number of iterations
c
 1200 if (.not. inp_i(ncysh)) then
         write(LuOut,*)' conv_inp: ncysh value not found; defaulting',
     &             ' to 9999 iterations'
         ncysh = 9999
      endif            
      if (.not. rtdb_put(rtdb, 'dft:ncysh', mt_int, 1, ncysh))
     &   call errquit('conv_inp: rtdb_put failed', 1200, RTDB_ERR)
c
      goto 10
c     
c     damp/ndamp; use damping with the new density = ndamp% 
c                 of old density + (1-ndamp)% of new density 
c
 1300 damp = .true.             
      if (.not. rtdb_put(rtdb, 'dft:damp', mt_log, 1, damp))
     &   call errquit('conv_inp: rtdb_put failed', 1300, RTDB_ERR)
      if (.not. inp_i(ndamp)) then
         write(LuOut,*)' conv_inp: ndamp value not found; defaulting',
     &             ' to 40%.'
         ndamp = 40
      endif            
      if (.not. rtdb_put(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     &   call errquit('conv_inp: rtdb_put failed', 1300, RTDB_ERR)
c
      goto 10
c     
c     nodamping; set no damping to be used for the entire calculation
c
 1400 nodamping = .true.
      if (.not. rtdb_put(rtdb, 'dft:nodamping', mt_log, 1, nodamping))
     &   call errquit('conv_inp: rtdb_put failed', 1400, RTDB_ERR)
c
      damp = .false.
      if (.not. rtdb_put(rtdb, 'dft:damp', mt_log, 1, damp))
     &   call errquit('dft_input: rtdb_put failed', 1400, RTDB_ERR)
      goto 10
c     
c     diis/nfock; use diis procedure with nfock Fock matrices 
c                 saved and used in the extrapolation
c
 1500 diis = .true.
      if (.not. rtdb_put(rtdb, 'dft:diis', mt_log, 1, diis))
     &   call errquit('conv_inp: rtdb_put failed', 1500, RTDB_ERR)
      if (.not. inp_i(nfock)) then
         write(LuOut,*) ' conv_inp: nfock value not found; defaulting',
     &             ' to 10.'
         nfock = 10
      endif            
c
c     nfock=0 same as nodiis
c
      if(nfock.eq.0) goto 1600
      if (.not. rtdb_put(rtdb, 'dft:nfock', mt_int, 1, nfock))
     &   call errquit('conv_inp: rtdb_put failed', 1500, RTDB_ERR)
c
      goto 10
c     
c     nodiis; set no diis to be used for the entire calculation
c
 1600 nodiis = .true.
      if (.not. rtdb_put(rtdb, 'dft:nodiis', mt_log, 1, nodiis))
     &   call errquit('conv_inp: rtdb_put failed', 1600, RTDB_ERR)
      diis = .false.
      if (.not. rtdb_put(rtdb, 'dft:diis', mt_log, 1, diis))
     &   call errquit('conv_inp: rtdb_put failed', 160, RTDB_ERR)
c
      goto 10
c     
c     lshift/rlshift; use levelshifting with a shift of rlshift applied
c                     to the diagonal elements of the unoccupied block
c                     of the Fock matrix.
c
 1700 levelshift = .true.
      if (.not. rtdb_put(rtdb, 'dft:levelshift', mt_log, 1, levelshift))
     &   call errquit('conv_inp: rtdb_put failed', 1700, RTDB_ERR)
      if (.not. inp_f(rlshift)) then
         write(LuOut,*)' conv_inp: lshift value not found; defaulting',
     &             ' to  0.5 Hartree'
         rlshift = 0.5d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:rlshift', mt_dbl, 1,rlshift))
     &   call errquit('conv_inp: rtdb_put failed', 1700, RTDB_ERR)
c
      goto 10
c     
c     nolevelshifting; set no levelshifting to be used for the 
c                      entire calculation
c
 1800 nolevelshifting = .true.
      if (.not. rtdb_put(rtdb, 'dft:nolevelshifting', mt_log, 1, 
     &   nolevelshifting))
     &   call errquit('conv_inp: rtdb_put failed', 1800, RTDB_ERR)
      levelshift = .false.
      if (.not. rtdb_put(rtdb, 'dft:levelshift', mt_log, 1, levelshift))
     &   call errquit('conv_inp: rtdb_put failed', 1800, RTDB_ERR)
c
      goto 10
c     
c     hl_tol; a Homo - Lumo gap tolerance used to invoke "automatic"
c             levelshifting and also monitoring of misuse of maximum
c             orbital overlaps
c
 1900 if (.not. inp_f(hl_tol)) then
         write(LuOut,*)' conv_inp: hl_tol value not found; defaulting',
     &             ' to 0.1 Hartree'
         hl_tol = 0.01d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:hl_tol', mt_dbl, 1, hl_tol))
     &   call errquit('conv_inp: rtdb_put failed', 1900, RTDB_ERR)
c
      goto 10
c     
c     fast (usually faster, but less safe)
c                      
 2000 levelshift = .false.
      if (.not. rtdb_put(rtdb, 'dft:levelshift', mt_log, 1, levelshift))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      damp = .true.             
      if (.not. rtdb_put(rtdb, 'dft:damp', mt_log, 1, damp))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      ndamp = 40
      if (.not. rtdb_put(rtdb, 'dft:ndamp', mt_int, 1, ndamp))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      ncydp = 0
      if (.not. rtdb_put(rtdb, 'dft:ncydp', mt_int, 1, ncydp))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      dampon=0d0
      if (.not. rtdb_put(rtdb, 'dft:dampon', mt_dbl, 1, dampon))
     &   call errquit('conv_inp: rtdb_put failed', 400, RTDB_ERR)
      dampoff=1d-4
      if (.not. rtdb_put(rtdb, 'dft:dampoff', mt_dbl, 1, dampoff))
     &   call errquit('conv_inp: rtdb_put failed', 400, RTDB_ERR)
      rlshift=0.0d0
      if (.not. rtdb_put(rtdb, 'dft:rlshift', mt_dbl, 1,rlshift))
     &   call errquit('conv_inp: rtdb_put failed', 1700, RTDB_ERR)
c
      goto 10
c     
c     stable (defaults)
c                      
 2010 levelshift = .true.
      if (.not. rtdb_put(rtdb, 'dft:levelshift', mt_log, 1, levelshift))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'dft:ncydp', mt_int, 1, ncydp))
     &   call errquit('conv_inp: rtdb_put failed', 2000, RTDB_ERR)
      rlshift=0.5d0
      if (.not. rtdb_put(rtdb, 'dft:rlshift', mt_dbl, 1,rlshift))
     &   call errquit('conv_inp: rtdb_put failed', 1700, RTDB_ERR)
      hl_tol = 0.01d0
      if (.not. rtdb_put(rtdb, 'dft:hl_tol', mt_dbl, 1, hl_tol))
     &   call errquit('conv_inp: rtdb_put failed', 2010, RTDB_ERR)
c
      goto 10
c     
c     Rabuck convergence method by annealed smearing
c                      
 2020 continue
      if (.not. inp_i(n_rabuck)) then
         write(LuOut,*)' conv_inp: energy convergence value not found;',
     &             ' defaulting to 1d-6 Hartree.'
         n_rabuck = 25
      endif
      if (.not. rtdb_put(rtdb, 'dft:smear_sigma', mt_dbl, 1, 1d-3))
     &   call errquit('conv_inp: rtdb_put failed', 2021, RTDB_ERR)
      if (.not.rtdb_put(rtdb, 'dft:n_rabuck', mt_int, 1, n_rabuck))
     &   call errquit('conv_inp: rtdb_put failed', 2022, RTDB_ERR)
c
      goto 10
c
c     Pcglevl (preconditioned conjugate gradient) level shifters
c
 2030 continue
      if (inp_f(scfshift(1))) then
        if (inp_f(scfshift(2))) then
          if (inp_f(scfshift(3))) then
            continue
          endif
        endif
      endif
      if (.not.rtdb_put(rtdb,'dft:level shift info',MT_DBL,6,scfshift))
     &  call errquit('conv_inp: rtdb_put of pcg scfshift failed',0,
     &               RTDB_ERR)
c
      goto 10
c
c     Nrlevl (Newton-Raphson) level shifters
c
 2040 continue
      if (inp_f(scfshift(4))) then
        if (inp_f(scfshift(5))) then
          if (inp_f(scfshift(6))) then
            continue
          endif
        endif
      endif
      if (.not.rtdb_put(rtdb,'dft:level shift info',MT_DBL,6,scfshift))
     &  call errquit('conv_inp: rtdb_put of nr scfshift failed',0,
     &               RTDB_ERR)
c
      goto 10
c
c     Nr <switch>
c
 2050 if (.not. inp_f(tolnr))
     &    call errquit('dft_input: nr switch',0,INPUT_ERR)
      if (.not. rtdb_put(rtdb,'dft:full hessian switch',MT_DBL,1,tolnr))
     &    call errquit('dft_input: rtdb put tolnr failed',0,RTDB_ERR)
c
      goto 10
 1999 continue
c     
      return
c     
10000 write(LuOut,10001)

10001 format(/,' convergence [[energy rval], [density rval], ',
     &       /,' [gradient rval], [dampon rval], [dampoff rval], ',
     &       /,' [diison rval], [diisoff rval], [levlon rval], ',
     &       /,' [levloff rval], [ncydp ival], [ncyds ival], ',
     &       /,' [ncysh ival], [damp ival], [nodamping], ',
     &       /,' [diis nfock ival], [nodiis], [lshift rval], ',
     &       /,' [nolevelshifting], [hl_tol rval]], ',
     &       /,' [pcg [rinit [rtol [rfinal]]]], ',
     &       /,' [nr [rinit [rtol [rfinal]]]] ')
      call util_flush(LuOut)
      call errquit('conv_input: invalid format', 0, INPUT_ERR)
c     
      end
