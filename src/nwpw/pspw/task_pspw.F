*
* $Id$
*
C> \ingroup task
C> @{
      logical function task_pspw(rtdb)
      implicit none
#include "errquit.fh"
#include "inp.fh"
#include "global.fh"
#include "bafdecls.fh"
#include "btdb.fh"
c
      integer rtdb
      logical value,dummy,md_only,wannier2
      integer MASTER
      parameter (MASTER=0)
c
      integer ind               ! Index of matched directive
      integer num_dirs          ! No. of known directives
      parameter (num_dirs = 27)
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
     >            'wannier',
     >            'expand_cell',
     >            'ionize',
     >            'convert_restricted',
     >            'Born-Oppenheimer',
     >            'rdf',
     >            'lcao',
     >            'pspw_et',
     >            'Metropolis',
     >            'nwpw-optimize',
     >            'aimd_properties',
     >            'exafs',
     >            'pauls_matrix_diag',
     >            'translate'/

c     *** define external functions ****
      logical  psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5,cgmdv5
      logical  v_wvfnc_init,wvfnc_expander,cell_optimize,mmsdv1,mmmdv1
      logical  pspw_dplot,pspw_stress,pspw_noit_energy,cpmdv6
      logical  pspw_wannier,pspw_wannier_old,expand_cell,wvfnc_ionize
      logical  wvfnc_convert_restricted,nwpw_rdf,lcaov1,pspw_et
      logical  cgmontecarlov1,cgoptimize1,pspw_translate
      logical  task_aimd_properties,task_md_xs,task_paulsmatrix
      logical  pspw_wannier2
      external psp_formatter,wvfnc_init,cpsdv5,cpmdv5,cgsdv5,cgmdv5
      external v_wvfnc_init,wvfnc_expander,cell_optimize,mmsdv1,mmmdv1
      external pspw_dplot,pspw_stress,pspw_noit_energy,cpmdv6
      external pspw_wannier,pspw_wannier_old,expand_cell,wvfnc_ionize
      external wvfnc_convert_restricted,nwpw_rdf,lcaov1,pspw_et
      external cgmontecarlov1,cgoptimize1,pspw_translate
      external task_aimd_properties,task_md_xs,task_paulsmatrix
      external pspw_wannier2
c
c
      value = .false.
      task_pspw = .false. 
 10   value = btdb_cget(rtdb,'task:operation',1,test)

      if (.not. inp_match(num_dirs, .false., test, dirs, ind))
     >     call errquit('task pspw: unknown directive', 0, INPUT_ERR)
c
c
c
      goto ( 100,  200, 300, 400, 500, 600, 700, 800, 900,
     >       1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
     >       1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700 ) ind
      call errquit('task pspw: unimplemented directive', ind, INPUT_ERR)


c     
c  psp_generator
c     
 100  value = .true.
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.true.)
      call Parallel_Init()        !** Parallel calls need since psp_formatter is also called directly from codes
      value = psp_formatter(rtdb)
      call Parallel_Finalize()    !** Parallel calls need since psp_formatter is also called directly from codes
      call ga_sync()
      goto 9999


c
c  wavefunction_initializer
c 
 300  value = .true.
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.true.)

      if (.not.btdb_get(rtdb,'nwpw:md_only',mt_log,1,md_only))
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
      value = btdb_parallel(.true.)
      value = cpmdv5(rtdb)
      call ga_sync()
      goto 9999

c
c  v_wavefunction_initializer
c 
 800  value = .true.
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
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
      value = btdb_parallel(.true.)
      value = pspw_noit_energy(rtdb)
      call ga_sync()
      goto 9999

c
c  wannier
c
 1400 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - Wannier <<<<"
      end if
      if (.not.btdb_get(rtdb,'nwpw:wannier_old',mt_log,1,dummy))
     >  dummy = .false.
      if (.not.btdb_get(rtdb,'nwpw:wannier2',mt_log,1,wannier2))
     >  wannier2 = .false.
      value = .true.
      value = btdb_parallel(.true.)
      if (dummy) then
         value = pspw_wannier_old(rtdb)
      else
         if (wannier2) then
            value = pspw_wannier2(rtdb)
         else
            value = pspw_wannier(rtdb)
         end if
      end if
      call ga_sync()
      goto 9999

c
c  expand_cell
c
 1500 value = .true.
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.false.)
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
      value = btdb_parallel(.true.)
      if (.not.btdb_get(rtdb,'nwpw:md_only',mt_log,1,md_only))
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
      value = btdb_parallel(.false.)
      if (ga_nodeid().eq.MASTER) then
        write(*,*)
     >     ">>>> PSPW Serial Module - rdf <<<<"
        write(*,*) "not finished"
        value = nwpw_rdf(rtdb)
      endif
      call ga_sync()
      goto 9999

c
c  lcao
c
 2000 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - LCAO <<<<"
      end if
      value = lcaov1(rtdb)
      call ga_sync()
      goto 9999
c
c  et
c
 2100 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - ET <<<<"
      end if
      value = .true.
      value = btdb_parallel(.true.)
      value = pspw_et(rtdb)
      call ga_sync()
      goto 9999

c
c  Metropolis
c
 2200 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - Metropolis <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)
      if (.not.btdb_get(rtdb,'nwpw:md_only',mt_log,1,md_only))
     >   md_only = .false.

      if (md_only) then
         !value = mmmdv1(rtdb,0)
         value = .false.
      else
         value = cgmontecarlov1(rtdb,0)
      end if

      call ga_sync()
      goto 9999


c
c  nwpw-optimize
c
 2300 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - nwpw-optimize <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)
      value = cgoptimize1(rtdb,0)

      call ga_sync()
      goto 9999

c
c  aimd_properties
c
 2400 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - aimd_properties <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)

      value = task_aimd_properties(rtdb)

      call ga_sync()
      goto 9999

c
c  exafs          
c
 2500 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - exafs <<<<"
      end if

      value = .true.
      value = btdb_parallel(.true.)

      value = task_md_xs(rtdb)

      call ga_sync()
      goto 9999
c
c  Pauls_matrix_diag
c
 2600 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - pauls_matrix_diag <<<<"
      end if

      value = task_paulsmatrix(rtdb)
      call ga_sync()
      goto 9999

c
c  pspw translate        
c
 2700 if (ga_nodeid().eq.MASTER) then
        write(*,*)  ">>>> PSPW Parallel Module - translate  <<<<"
      end if

      value = pspw_translate(rtdb)
      call ga_sync()


 9999 continue
      dummy = btdb_parallel(.true.)
      task_pspw = value

      return
      end
C> @}
