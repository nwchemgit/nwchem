      logical function task_pspw(rtdb)
c
c $Id: task_pspw.F,v 1.1 2001-08-30 00:22:02 edo Exp $
c
      implicit none
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
      logical value,dummy
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 11)
      character*30 dirs(num_dirs)
      character*255 test
      data dirs / 'psp_generator',
     >            'psp_formatter',
     >            'wavefunction_initializer',
     >            'wavefunction_converter',
     >            'steepest_descent',
     >            'conjugate_gradient',
     >            'Car-Parrinello',
     >            'v_wavefunction_initializer',
     >            'wavefunction_expander',
     >            'cell_optimize',
     >            'pspw_dplot' /

c     *** define external functions ****
      logical  psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5
      logical  v_wvfnc_init,wvfnc_expander,cell_optimize
      logical  pspw_dplot
      external psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5
      external v_wvfnc_init,wvfnc_expander,cell_optimize
      external pspw_dplot
c
c
      value = .false.
      task_pspw = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task pspw: unknown directive', 0)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900,
     >       1000, 1100 ) ind
      call errquit('task pspw: unimplemented directive', ind)


c     
c  psp_generator
c     
 100  value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> PSPW Serial Module - psp_generator <<<<"
c       value = psp_generator(rtdb)
      endif
      call ga_sync()
      goto 9999


c
c  psp_formatter
c
 200  value = .true.
      
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> PSPW parallel Module - psp_formatter <<<<"  
      end if
      value = rtdb_parallel(.true.)
      value = psp_formatter(rtdb)
      call ga_sync()
      goto 9999


c
c  wavefunction_initializer
c 
 300  value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then 
         write(*,*) 
     >      ">>>> PSPW Serial Module - wavefunction_initializer <<<<"
         value = wvfnc_init(rtdb)
      end if
      call ga_sync()
      goto 9999


c
c  wavefunction_converter
c
 400  value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)  
     >     ">>>> PSPW Serial Module - wavefunction_converter <<<<"
        write(*,*) "not finished"
c       value = wavefunction_converter(rtdb)
      endif
      call ga_sync()
      goto 9999


c
c  steepest_descent
c
 500  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - steepest_descent <<<<"
      end if

      value = .true.
      value = rtdb_parallel(.true.)
      value = cpsdv5(rtdb)
      call ga_sync()
      goto 9999


c
c  conjugate gradient
c
 600  if (ga_nodeid().eq.MASTER) then
        write(*,*) ">>>> PSPW Parallel Module - conjugate_gradient <<<<"
      end if
      value =  cgsdv5(rtdb,0)
      call ga_sync()
      goto 9999


c
c  Car-Parrinello
c
 700  if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - Car-Parrinello <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = cpmdv5(rtdb)
      call ga_sync()
      goto 9999

c
c  v_wavefunction_initializer
c 
 800  value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then 
         write(*,*) 
     >      ">>>> PSPW Serial Module - v_wavefunction_initializer <<<<"
         value = v_wvfnc_init(rtdb)
      end if
      call ga_sync()
      goto 9999

c
c  wavefunction_expander
c
 900  value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)  
     >     ">>>> PSPW Serial Module - wavefunction_expander <<<<"
        write(*,*) "not finished"
        value = wvfnc_expander(rtdb)
      endif
      call ga_sync()
      goto 9999

c
c  cell_optimize
c
 1000 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - cell_optimize <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = cell_optimize(rtdb)
      call ga_sync()
      goto 9999

c
c  pspw_dplot
c
 1100 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - pspw_dplot <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = pspw_dplot(rtdb)
      call ga_sync()
      goto 9999


 9999 continue
      dummy = rtdb_parallel(.true.)
      task_pspw = value

      return
      end
