      subroutine grid_input(rtdb, module)
C     $Id: 
      implicit none
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "cdft.fh"
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
c     xcoarse = energy accuracy 1.d-4
c     coarse  = energy accuracy 1.d-5
c     medium  = energy accuracy 1.d-6
c     fine    = energy accuracy 1.d-7
c     xfine   = energy accuracy 1.d-8
c     user
c     store_wght
c     nquad_task
c     rm
c     
      integer num_dirs, ind, mlen
      integer atn
      character*2 symbol
      parameter (num_dirs = 13)
      character*10 dirs(num_dirs)
      character*16 test,element
      logical lbecke
      integer accqrad
      data dirs / 'delley', 'becke', 'gausleg', 'lebedev', 'xcoarse', 
     &            'coarse', 'medium', 'fine', 'xfine', 'user', 
     &            'store_wght', 'nquad_task', 'rm'/
c     
      ngridin=0
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
     &      1200, 1300, 1999) ind
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
  300 gridtype='gausleg'
c
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 300) 
c
      if (.not. inp_i(nradpts)) then
         write(6,*)' grid_inp: nradpts value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 900)
      endif
c
      if (.not. inp_i(nagrid)) then
         write(6,*)' grid_inp: nagrid value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &      call errquit('grid_inp: rtdb_put failed', 900)
      endif
c
      goto 10
c     
c     lebedev; use Lebedev angular quadrature
c     
  400 gridtype='lebedev'
      leb=.true.
c
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 400) 
c
      if (.not. inp_i(nradpts)) then
         write(6,*)' grid_inp: nradpts value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 400)
      endif
c
      if (.not. inp_i(iangquad)) then
         write(6,*)' grid_inp: iangquad value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     &      call errquit('grid_inp: rtdb_put failed', 400)
      endif
c
      goto 10
c     
c     coarse; use coarse grid definition
c             
  500 gridtype='xcoarse'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 500) 
c
      goto 10
c     
c     coarse; use coarse grid definition
c             
  600 gridtype='coarse'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 600) 
c
      goto 10
c     
c     medium; use medium grid definition
c             
  700 gridtype='medium'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 700) 
c
      goto 10
c     
c     fine; use fine grid definition
c             
  800 gridtype='fine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 800) 
c
      goto 10
c     
c     xfine; use xfine grid definition
c             
  900 gridtype='xfine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 900) 
c
      goto 10
c     
c     user; use user defined grid
c             
 1000 gridtype='user'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 1000) 
      if (.not. inp_i(nradpts)) then
         write(6,*)' grid_inp: nradpts value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 1000)
      endif
      if (.not. inp_i(nagrid)) then
         write(6,*)' grid_inp: nagrid value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &      call errquit('grid_inp: rtdb_put failed', 1000)
      endif
c
      goto 10
c     
c     store_wght; store quadrature weights
c     
 1100 store_wght = .true.
      if (.not. rtdb_put(rtdb, 'dft:store_wght', mt_log, 1, store_wght))
     &   call errquit('grid_inp: rtdb_put failed', 1100)
c
      goto 10
c     
c     nquad_task; redefine granularity of tasks in parallelism of grid
c     
 1200 if (.not. inp_i(nquad_task)) then
         write(6,*)' grid_inp: nquad_task value not found; ',
     &             ' defaulting to 1 '
         nquad_task = 1
      endif            
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, nquad_task))
     &   call errquit('grid_inp: rtdb_put failed', 1200)
c
      goto 10
c     
c     rm; the exponent in defining the Euler-Maclaurin radial quadrature
c     
 1300 if (.not. inp_f(rm)) then
         write(6,*)' grid_inp: rm value not found; defaulting ',
     &             ' to 2.0d0 '
         rm = 2.0d0
      endif            
      if (.not. rtdb_put(rtdb, 'dft:rm', mt_dbl, 1, rm))
     &   call errquit('grid_inp: rtdb_put failed', 1300)
c
      goto 10
c     
 1999 continue
c     
      return
c     
10000 continue
c
c     see if we have tag + rad + ang + accqrad
c
c            if (.not. inp_a(tag)) goto 1688
      if (.not. geom_tag_to_element(test, symbol, element, atn))
     .     goto 1688
!     $     call errquit
!     &    ('grid_inp: tag does not refer to an atom', 0)
      ngridin=ngridin+1
      znucin(ngridin)=atn
      if (.not. inp_i(nradpts)) then
         write(6,*)' grid_inp: nradpts value not found. '
      else
         nradin(ngridin)=nradpts
      endif
c
      if (.not. inp_i(iangquad)) then
         write(6,*)' grid_inp: iangquad value not found. '
      else
         nangin(ngridin)=iangquad
      endif
      if (.not. inp_i(accqrad)) then
         write(6,*)' grid_inp: iangquad value not found. '
      else
         accqrin(ngridin)=accqrad
      endif
      goto 10
 1688       write(6,10001)
10001 format(/,' grid [[delley], [becke], [gausleg], 
     &       [lebedev ival ival], ',
     &       /,' [xcoarse], [coarse], [medium], [fine], [xfine], ',
     &       /,' [user ival ival], [store_wght], [nquad_task ival], ',
     &       /,' [rm rval] ')
      call util_flush(6)
      call errquit('grid_input: invalid format', 0)
c     
      end


