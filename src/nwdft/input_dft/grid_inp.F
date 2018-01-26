      subroutine grid_input(rtdb, module)
C     $Id: 
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "cdft.fh"
#include "stdio.fh"
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
c     huge    = energy accuracy 1.d-10
c     user
c     store_wght
c     nquad_task
c     rm
c     
      integer num_dirs, ind, mlen
      double precision ddumm
      integer atn
      character*2 symbol
      parameter (num_dirs = 26)
      character*10 dirs(num_dirs)
      character*16 test,element
      logical aintgot_rad,aintgot_ang
      logical lbecke
      data dirs / 'delley', 'becke', 'gausleg', 'lebedev', 'xcoarse', 
     &     'coarse', 'medium', 'fine', 'xfine', 'user',
     H     'huge',
C     .            'euler','mura','treutler','lindh',
     .            'euler','mura','treutler','xxxx',
     .            'ssf','erf1','erf2',
     &            'store_wght', 'nquad_task', 'rm',
     .     'old','new','nosymm','symm','nodisk'/
c     
      ngridin=0
c
c     Lebedev on by default
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
      goto (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
     H     1111,
     .     2001,2100,2200,2300,1700, 1800,1900,1100, 
     &      1200, 1300, 1400,1500, 1600, 1650,1660,1999) ind
      call errquit('grid_inp: unimplemented directive', ind,
     &       INPUT_ERR)
c     
c     delley; use delley spatial weights
c             
  100 continue
c      ldelley = .true.
       call errquit(' Delley weights no longer working', 100,
     &       INPUT_ERR)
      if (.not. rtdb_put(rtdb, 'dft:delley', mt_log, 1, ldelley))
     &   call errquit('grid_inp: rtdb_put failed', 100, RTDB_ERR)
c
      goto 10
c     
c     becke; use becke spatial weights
c             
  200 lbecke = .true.
      if (.not. rtdb_put(rtdb, 'dft:becke', mt_log, 1, lbecke))
     &   call errquit('grid_inp: rtdb_put failed', 200, RTDB_ERR)
      lssw=.false.
      if (.not. rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw))
     &   call errquit('grid_inp: rtdb_put failed', 1700, RTDB_ERR)
c
      goto 10

c     
c     ssf; use ssf spatial weights
c             
 1700 lssw=.true.
      if (.not. rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw))
     &   call errquit('grid_inp: rtdb_put failed', 1700, RTDB_ERR)
      if (.not.rtdb_cput(rtdb, 'dft:ssweights',  1, 'ssf '))
     &   call errquit('grid_inp: rtdb_cput failed', 1700, RTDB_ERR)
      whichssw='ssf '
c
      goto 10
c     
c     erf1; use erf1 spatial weights
c             
 1800 lssw=.true.
      if (.not. rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw))
     &   call errquit('grid_inp: rtdb_put failed', 1800, RTDB_ERR)
      if (.not.rtdb_cput(rtdb, 'dft:ssweights',  1, 'erf1'))
     &   call errquit('grid_inp: rtdb_cput failed', 1800, RTDB_ERR)
      whichssw='erf1'
c
      goto 10
c     
c     erf2; use erf2 spatial weights
c             
 1900 lssw=.true.
      if (.not. rtdb_put(rtdb, 'dft:ssw', mt_log, 1, lssw))
     &   call errquit('grid_inp: rtdb_put failed', 1900, RTDB_ERR)
      if (.not.rtdb_cput(rtdb, 'dft:ssweights',  1, 'erf2'))
     &   call errquit('grid_inp: rtdb_cput failed', 1900, RTDB_ERR)
      whichssw='erf2'
c
      goto 10
c     
c     gausleg; use Gauss-Legendre for angular quadrature (ala Handy) 
c     
  300 gridtype='gausleg'
      leb=.false.
      if (.not. rtdb_put(rtdb, 'dft:lebedev', mt_log, 1, leb))
     &   call errquit('grid_inp: rtdb_put failed', 590, RTDB_ERR)
c
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 300, RTDB_ERR)
c
      if (.not. inp_i(nradpts)) then
         write(LuOut,*)' grid_inp: nradpts value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 900, RTDB_ERR)
      endif
c
      if (.not. inp_i(nagrid)) then
         write(LuOut,*)' grid_inp: nagrid value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &      call errquit('grid_inp: rtdb_put failed', 900, RTDB_ERR)
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
     &   call errquit('grid_inp: rtdb_put failed', 400, RTDB_ERR)
c
      if (inp_i(nradpts)) then
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 400, RTDB_ERR)
      endif
c
      if (inp_i(iangquad)) then
         if (.not. rtdb_put(rtdb, 'dft:iangquad', mt_int, 1, iangquad))
     &      call errquit('grid_inp: rtdb_put failed', 400, RTDB_ERR)
      endif
c
      goto 10
c     
c     coarse; use coarse grid definition
c             
  500 gridtype='xcoarse'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 500, RTDB_ERR)
c
      goto 10
c     
c     coarse; use coarse grid definition
c             
  600 gridtype='coarse'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 600, RTDB_ERR)
c
      goto 10
c     
c     medium; use medium grid definition
c             
  700 gridtype='medium'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 700, RTDB_ERR)
c
      goto 10
c     
c     fine; use fine grid definition
c             
  800 gridtype='fine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 800, RTDB_ERR)
c
      goto 10
c     
c     xfine; use xfine grid definition
c             
  900 gridtype='xfine'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 900, RTDB_ERR)
c
      goto 10
c     
c     
c     huge; use huge grid definition
c             
 1111 gridtype='huge'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 900, RTDB_ERR)
c
      goto 10
c     
c     user; use user defined grid
c             
 1000 gridtype='user'
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
      if (.not. inp_i(nradpts)) then
         write(LuOut,*)' grid_inp: nradpts value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nradpts', mt_int, 1, nradpts))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
      endif
      if (.not. inp_i(nagrid)) then
         write(LuOut,*)' grid_inp: nagrid value not found. '
      else
         if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
      endif
c
      goto 10
c     
c     store_wght; store quadrature weights
c     
 1100 store_wght = .true.
      if (.not. rtdb_put(rtdb, 'dft:store_wght', mt_log, 1, store_wght))
     &   call errquit('grid_inp: rtdb_put failed', 1100, RTDB_ERR)
c
      goto 10
c     
c     nquad_task; redefine granularity of tasks in parallelism of grid
c     
 1200 if (.not. inp_i(nquad_task)) then
         write(LuOut,*)' grid_inp: nquad_task value not found; ',
     &             ' defaulting to 1 '
         nquad_task = 1
      endif            
      if (.not. rtdb_put(rtdb, 'dft:nquad_task', mt_int, 1, nquad_task))
     &   call errquit('grid_inp: rtdb_put failed', 1200, RTDB_ERR)
c
      goto 10
c     
c     rm; the exponent in defining the Euler-Maclaurin radial quadrature
c     
 1300 if (.not. inp_f(ddumm)) then
      endif            
      write(LuOut,*)'WARNING: rm keyword NOT recognized'
c
      goto 10
c
c     Euler_McLaurin
c
 2001 continue
      radgrid='eumac'
      goto 2301
c
c     Mura-Knowles radial grid
c
 2100 continue
      radgrid='murakn'
      goto 2301
c
c     Treutler-Ahlrichs radial grid
c
 2200 continue
      radgrid='tarad'
      goto 2301
c
c     Roland Lindh radial grid
c
 2300 continue
      call errquit(' wrong keyword ',999, INPUT_ERR)
      radgrid='rlindh'

 2301 if (.not. rtdb_cput(rtdb,'dft:radgrid', 1,radgrid)) 
     &   call errquit('grid_inp: rtdb_put failed', 2300, RTDB_ERR)
      goto 10

 1400 continue
c
c     old grid schem
c
      call errquit(' old grid scheme no longer available',0,0)
c
      goto 10
 1500 continue
c
c     new grid schem
c
c
      goto 10
 1600 continue
c
c     symm off (lebedev)
c
      if(.not.rtdb_put(rtdb,'dft:gridsym',mt_log,1,.false.))
     &   call errquit('grid_inp: rtdb_put failed', 1600, RTDB_ERR)
c
      goto 10
 1650 continue
c
c     symm on (lebedev)
c
      if(.not.rtdb_put(rtdb,'dft:gridsym',mt_log,1,.true.))
     &   call errquit('grid_inp: rtdb_put failed', 1700, RTDB_ERR)
c
      goto 10
 1660 continue
c
c     dont use disk for gridpts
c
      if(.not.rtdb_put(rtdb,'dft:largenode',mt_log,1,.true.))
     &   call errquit('grid_inp: rtdb_put failed', 1660, RTDB_ERR)
c
      goto 10
 1999 continue
c
      if(ngridin.ne.0) then
         if (.not. rtdb_put(rtdb, 'dft:ngridin', mt_int, 1, ngridin))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
         if(.not. rtdb_cput(rtdb,'dft:grtagin',ngridin,grtagin))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
         if(.not.rtdb_put(rtdb,'dft:nradin',mt_int,ngridin,nradin))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
         if(.not.rtdb_put(rtdb,'dft:nangin',mt_int,ngridin,nangin))
     &      call errquit('grid_inp: rtdb_put failed', 1000, RTDB_ERR)
      endif
      return
c     
10000 continue
c
c     see if we have tag + rad + ang 
c
c            if (.not. inp_a(tag)) goto 1688
      if (.not. geom_tag_to_element(test, symbol, element, atn))
     .     goto 1688
!     $     call errquit
!     &    ('grid_inp: tag does not refer to an atom', 0)
      ngridin=ngridin+1
c     lowercase
cdbg      write(6,*) ' symbol ',symbol,'@@'
      aintgot_ang=.true.
      aintgot_rad=.true.
      call inp_lcase(test)
      grtagin(ngridin)=test
      if (.not. inp_i(nradpts)) then
         write(LuOut,*)' grid_inp: nradpts value not found. '
      else
         aintgot_rad=.false.
         nradin(ngridin)=nradpts
      endif
c
      if (.not. inp_i(iangquad)) then
         write(LuOut,*)' grid_inp: iangquad value not found. '
      else
         aintgot_ang=.false.
         nangin(ngridin)=iangquad
      endif
      if(aintgot_ang.and.aintgot_rad) goto 1688
      goto 10
 1688       write(LuOut,10001)
10001 format(/,' grid [ [becke], [gausleg], 
     &       [lebedev ival ival], ',
     &       /,' [xcoarse],[coarse],[medium],[fine],[xfine],[huge], ',
     &       /,' [user ival ival], [store_wght], [nquad_task ival], ')
      call util_flush(LuOut)
      call errquit('grid_input: invalid format', 0, INPUT_ERR)
c     
      end


