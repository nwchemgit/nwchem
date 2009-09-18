*
* $Id: task_pspw.F,v 1.15 2008-11-03 06:06:10 bylaska Exp $
*

      logical function task_pspw(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
      logical value,dummy,md_only
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 19)
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
     >            'pspw_dplot',
     >            'stress',
     >            'noit_energy',
     > 		  'wannier',
     >            'expand_cell',
     >            'ionize',
     >            'convert_restricted',
     >            'Born-Oppenheimer',
     >            'rdf'/

c     *** define external functions ****
      logical  psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5,cgmdv5
      logical  v_wvfnc_init,wvfnc_expander,cell_optimize,mmsdv1,mmmdv1
      logical  pspw_dplot,pspw_stress,pspw_noit_energy,cpmdv6
      logical  pspw_wannier,expand_cell,wvfnc_ionize
      logical  wvfnc_convert_restricted,nwpw_rdf
      external psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5,cgmdv5
      external v_wvfnc_init,wvfnc_expander,cell_optimize,mmsdv1,mmmdv1
      external pspw_dplot,pspw_stress,pspw_noit_energy,cpmdv6
      external pspw_wannier,expand_cell,wvfnc_ionize
      external wvfnc_convert_restricted,nwpw_rdf
c
c
      value = .false.
      task_pspw = .false. 
 10   value = rtdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task pspw: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900,
     >       1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     >       1900 ) ind
      call errquit('task pspw: unimplemented directive', ind, INPUT_ERR)


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
c  psp_formatter - This is outdated code, but kept for legacy purposes.
c
 200  value = .true.
      
      if (ga_nodeid().eq.MASTER) then 
        write(*,*)  ">>>> PSPW parallel Module - psp_formatter <<<<"  
      end if
      value = rtdb_parallel(.true.)
      call Parallel_Init()        !** Parallel calls need since psp_formatter is also called directly from codes
      value = psp_formatter(rtdb)
      call Parallel_Finalize()    !** Parallel calls need since psp_formatter is also called directly from codes
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

      if (.not.rtdb_get(rtdb,'nwpw:md_only',mt_log,1,md_only))
     >   md_only = .false.

      if (md_only) then
         value = mmsdv1(rtdb)
      else
         value = cpsdv5(rtdb)
      end if

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

c
c  stress - debug
c
 1200 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - stress <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = pspw_stress(rtdb)
      call ga_sync()
      goto 9999

c
c  noit_energy - debug
c
 1300 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - noit_energy <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = pspw_noit_energy(rtdb)
      call ga_sync()
      goto 9999

c
c  wannier
c
 1400 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - Wannier <<<<"
      end if
      value = .true.
      value = rtdb_parallel(.true.)
      value = pspw_wannier(rtdb)
      call ga_sync()
      goto 9999

c
c  expand_cell
c
 1500 value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)
     >     ">>>> PSPW Serial Module - expand_cell <<<<"
        write(*,*) "not finished"
        value = expand_cell(rtdb)
      endif
      call ga_sync()
      goto 9999

c
c  ionize
c
 1600 value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)
     >     ">>>> PSPW Serial Module - ionize <<<<"
        write(*,*) "not finished"
        value = wvfnc_ionize(rtdb)
      endif
      call ga_sync()
      goto 9999


c
c  convert_restricted
c
 1700 value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)
     >     ">>>> PSPW Serial Module - convert_restricted <<<<"
        write(*,*) "not finished"
        value = wvfnc_convert_restricted(rtdb)
      endif
      call ga_sync()
      goto 9999

c
c  Born-Oppenheimer
c
 1800 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - Born-Oppenheimer <<<<"
      end if

      value = .true.
      value = rtdb_parallel(.true.)
      if (.not.rtdb_get(rtdb,'nwpw:md_only',mt_log,1,md_only))
     >   md_only = .false.

      if (md_only) then
         value = mmmdv1(rtdb,0)
      else
         value = cgmdv5(rtdb,0)
      end if

      call ga_sync()
      goto 9999

c
c  rdf
c
 1900 value = .true.
      value = rtdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)
     >     ">>>> PSPW Serial Module - rdf <<<<"
        write(*,*) "not finished"
        value = nwpw_rdf(rtdb)
      endif
      call ga_sync()
      goto 9999


 9999 continue
      dummy = rtdb_parallel(.true.)
      task_pspw = value

      return
      end
