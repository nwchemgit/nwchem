      subroutine grid_input(rtdb, module)
C     $Id: 
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      character*(*) module
c     
c     Parse the grid directive which specifies how to construct the
c     quadrature for the numerical integrations needed in the DFT.  
c     Insert options into the database prefixed with module:...
c     
c     Possible variables are:
c
c     delley
c     becke
c     gausleg
c     lebedev
c     coarse
c     medium
c     fine
c     xfine
c     user
c     store_wght
c     nquad_task
c     rm
c     
      integer num_dirs, ind, mlen
      parameter (num_dirs = 12)
      character*10 dirs(num_dirs), gridtype
      character*255 test
      logical ldelley, lbecke, store_wght
      integer iangquad, nquad_task, nagrid, nradpts
      double precision rm
      data dirs / 'delley', 'becke', 'gausleg', 'lebedev', 'coarse', 
     &            'medium', 'fine', 'xfine', 'user', 'store_wght', 
     &            'nquad_task', 'rm'/
c     
      mlen = inp_strlen(module)
c     
 10   if (.not. inp_a(test)) goto 1999
c     
      if (.not. inp_match(num_dirs, .false., test, dirs, ind)) then
c     
c        Does not match a keyword ... 
c     
         goto 10
      endif
c     
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 
     &      1100, 1200, 1999) ind
      call errquit('grid_inp: unimplemented directive', ind)
c     
c     delley; use delley spatial weights
c             
  100 ldelley = .true.
      if (.not. rtdb_put(rtdb, 'dft:delley', mt_log, 1, ldelley))
     &   call errquit('grid_inp: rtdb_put failed', 100)
c
      goto 10
c     
c     becke; use becke spatial weights
c             
  200 lbecke = .true.
      if (.not. rtdb_put(rtdb, 'dft:becke', mt_log, 1, lbecke))
     &   call errquit('grid_inp: rtdb_put failed', 200)
c
      goto 10
c     
c     gausleg; use Gauss-Legendre for angular quadrature (ala Handy) 
c     
  300 iangquad=0
      if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     &   call errquit('grid_inp: rtdb_put failed', 300)
c
      goto 10
c     
c     lebedev; use Lebedev angular quadrature
c     
  400 if (.not. inp_i(iangquad)) then
         write(6,*)' grid_inp: lebedev value not found; defaulting ',
     &             ' to 4.'
         iangquad=4
      endif
      if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     &   call errquit('grid_inp: rtdb_put failed', 400)
c
      goto 10
c     
c     coarse; use coarse grid definition
c             
  500 gridtype='coarse'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 500) 
c
      goto 10
c     
c     medium; use medium grid definition
c             
  600 gridtype='medium'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 600) 
c
      goto 10
c     
c     fine; use fine grid definition
c             
  700 gridtype='fine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 700) 
c
      goto 10
c     
c     xfine; use xfine grid definition
c             
  800 gridtype='xfine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 800) 
c
      goto 10
c     
c     user; use user defined grid
c             
  900 gridtype='user'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 900) 
      if (.not. inp_i(nagrid)) then
         write(6,*)' grid_inp: nagrid value not found; ',
     &             ' defaulting to 8 '
         nagrid = 8
      endif
      if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &   call errquit('grid_inp: rtdb_put failed', 900)
      if (.not. inp_i(nradpts)) then
         write(6,*)' grid_inp: nradpts value not found; ',
     &             ' defaulting to 40 '
         nradpts = 40
      endif
      if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &   call errquit('grid_inp: rtdb_put failed', 900)
c
      goto 10
c     
c     store_wght; store quadrature weights
c     
 1000 store_wght = .true.
      if (.not. rtdb_put(rtdb, 'dft:store_wght', mt_log, 1, store_wght))
     &   call errquit('grid_inp: rtdb_put failed', 1000)
c
      goto 10
c     
c     nquad_task; redefine granularity of tasks in parallelism of grid
c     
 1100 if (.not. inp_i(nquad_task)) then
         write(6,*)' grid_inp: nquad_task value not found; ',
     &             ' defaulting to 1 '
         nquad_task = 1
      endif            
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, nquad_task))
     &   call errquit('grid_inp: rtdb_put failed', 1100)
c
      goto 10
c     
c     rm; the exponent in defining the Euler-Maclaurin radial quadrature
c     
 1200 if (.not. inp_f(rm)) then
         write(6,*)' grid_inp: rm value not found; defaulting ',
     &             ' to 2.0d0 '
         rm = 2.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:rm', mt_dbl, 1, rm))
     &   call errquit('grid_inp: rtdb_put failed', 1200)
c
      goto 10
c     
 1999 continue
c     
      return
c     
10000 write(6,10001)
10001 format(/,' grid [[delley], [becke], [gausleg], [lebedev ival], ',
     &       /,' [coarse], [medium], [fine], [xfine], ',
     &       /,' [user ival ival], [store_wght], [nquad_task ival], ',
     &       /,' [rm rval] ')
      call util_flush(6)
      call errquit('grid_input: invalid format', 0)
c     
      end


