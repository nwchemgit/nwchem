c     $Id$
      logical function band_energy(rtdb)
      implicit none
      integer rtdb
      logical  band_minimizer
      external band_minimizer

      call band_cutoff_boot_psi(rtdb)

      band_energy = band_minimizer(rtdb,0)
      return
      end

*     *********************************************
*     *                                           *
*     *             band_cutoff_boot_psi          *
*     *                                           *
*     *********************************************

      subroutine band_cutoff_boot_psi(rtdb)
      implicit none
      integer rtdb

#include "rtdb.fh"
#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical value,wcutset,ecutset,runme,use_sym
      integer i,steps
      real*8  wcut0,ecut0,wcut,ecut

*     **** external functions ****
      logical  band_minimizer,control_read,psi_filefind
      external band_minimizer,control_read,psi_filefind
      real*8   lattice_wcut,lattice_ecut
      external lattice_wcut,lattice_ecut

               
      if (.not.rtdb_get(rtdb,'nwpw:rho_use_symmetry',mt_log,1,use_sym))
     >   use_sym = .false.

      if (.not.rtdb_get(rtdb,'nwpw:cutoff_boot_psi',mt_log,1,
     >                  runme)) 
     >   runme = .true.
         
      !**** wavefunction is new then run a sequence of cutoffs ***
      if (runme) then
         call Parallel_Init()
         if (.not.control_read(3,rtdb))
     >      call errquit('band_cutoff_boot_psi:error reading control',0,
     >                   RTDB_ERR)
         call Parallel_Finalize()

         if (.not.psi_filefind()) then

            if (use_sym) then
               if (.not.rtdb_put(rtdb,'nwpw:rho_use_symmetry',
     >            mt_log,1,.false.))
     >            call errquit(
     >         'band_cutoff_boot_psi:error writing rho_use_symmetry',2,
     >         RTDB_ERR)
            end if

            call lattice_init()
            wcut0 = lattice_wcut()
            ecut0 = lattice_ecut()

            if (rtdb_get(rtdb,'band:wcut',mt_dbl,1,wcut)) 
     >         wcutset = .true.
            if (rtdb_get(rtdb,'band:ecut',mt_dbl,1,ecut))
     >         ecutset = .true.
            if (wcut0.lt.5.0d0) then
               steps = 0
            else if (wcut0.le.20.0d0) then
               steps = 2
            else if (wcut0.le.30.0d0) then
               steps = 3
            else 
               steps = 4
            end if

            do i=1,steps-1
               wcut = (i/dble(steps))*wcut0
               ecut = (i/dble(steps))*ecut0
               value = rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut)
     >            .and.rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut)
               if (.not.value) call errquit(
     >           'band_cutoff_boot_psi:error writing wcut,ecut',1,
     >            RTDB_ERR)
               value = band_minimizer(rtdb,0)
            end do

            if (wcutset) then
               if (.not.rtdb_put(rtdb,'band:wcut',mt_dbl,1,wcut0))
     >         call errquit('band_cutoff_boot_psi:error writing wcut',2,
     >                      RTDB_ERR)
            else
               if (.not.rtdb_delete(rtdb,'band:wcut'))
     >         call errquit('band_cutoff_boot_psi:error delete wcut',3,
     >                      RTDB_ERR)
            end if
            if (ecutset) then
               if (.not.rtdb_put(rtdb,'band:ecut',mt_dbl,1,ecut0))
     >         call errquit('band_cutoff_boot_psi:error writing ecut',4,
     >                      RTDB_ERR)
            else
               if (.not.rtdb_delete(rtdb,'band:ecut'))
     >         call errquit('band_cutoff_boot_psi:error delete ecut',5,
     >                      RTDB_ERR)
            end if

            if (use_sym) then
               if (.not.rtdb_put(rtdb,'nwpw:rho_use_symmetry',
     >            mt_log,1,.true.))
     >            call errquit(
     >         'band_cutoff_boot_psi:error writing rho_use_symmetry',3,
     >         RTDB_ERR)
            end if


         end if

         !*** set the runme variable to false ****
         if (.not.rtdb_put(rtdb,'nwpw:cutoff_boot_psi',mt_log,1,
     >               .false.)) 
     >   call errquit('band_cutoff_boot_psi:error setting runme',6,
     >                RTDB_ERR)

      end if !*runme*

      return
      end

