      logical function grad_dft(rtdb)
*
* $Id: grad_dft.F,v 1.26 1998-03-03 02:13:24 d3h449 Exp $
*
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "schwarz.fh"
#include "stdio.fh"
#include "cdft.fh"
c
c     dftgrad module.
c
c     Context is '...:dftgrad' --> changed to '...:dft'   JAN
c
c     Assumes DFT has been completed, MO vectors stored
c     and all information is still in the RTDB
c
      logical int_normalize
      external int_normalize
      logical int_norm_2c
      external int_norm_2c
c
      integer rtdb              ! [input] database handle
      integer nbases
      logical converged, status
      integer bases(3), plevel
c
c----------------------------------------------------------------------- 
c
c     Push context down to DFT
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push()
      call util_print_rtdb_load(rtdb,'dft')
      call ecce_print_module_entry ('dft')
c
      status = rtdb_get(rtdb, 'dft:converged', MT_LOG, 1, converged) 
      if (.not.( status .and. converged ))then
         if ((ga_nodeid() .eq. 0) 
     &      .and. util_print('information', print_none) )then
            write (luout,*)'status: ', status, '   converged: ',
     &                     converged
        endif
        call errquit(
     &       'dft grad_dft: no converged DFT wavefunction available',
     &     0)
      endif 
c
c     Extract high level info from the data-base setting defaults
c
c     load DFT parameters into common; turn off DFT printing
c
      call util_print_get_level(plevel)
      call util_print_set_level(print_none)
      call dft_rdinput(rtdb)
      call util_print_set_level(plevel)
c
      if (XCFIT.and.CDFIT)then
         nbases = 3
         bases(1) = AO_bas_han
         bases(2) = CD_bas_han
         bases(3) = XC_bas_han
      elseif((.not.XCFIT).and.CDFIT)then
         nbases = 2
         bases(1) = AO_bas_han
         bases(2) = CD_bas_han
      elseif((.not.CDFIT).and.XCFIT)then
         nbases = 2
         bases(1) = AO_bas_han
         bases(2) = XC_bas_han
      else
         nbases = 1
         bases(1) = AO_bas_han
      endif
c
c     initialize for schwarz screening
c
      call int_app_set_no_texas(rtdb)
      call int_init(rtdb, 1, ao_bas_han)
      call schwarz_init (geom, ao_bas_han)
      call int_terminate()
c
c     initialize for derivative integrals
c
      call intd_init(rtdb,nbases,bases)
c
      if (ga_nodeid() .eq. 0)then
         if (util_print('information',print_default) )then
            write(luout,*)
            write(luout,*)
            call util_print_centered(luout,
     &          'NWChem DFT Gradient Module',40,.true.)
            write(luout,*)
            write(luout,*)
            if (title .ne. ' ')then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
               write(luout,*)
            endif

            if (ipol .eq. 1)then
               write(luout,1) rcharge, 'closed shell'
            else 
               write(luout,1) rcharge, 'open shell'
            endif 
            call util_flush(luout)

         endif

 1       format(/
     &        '  charge          = ', f6.2/
     &        '  wavefunction    = ', a/)
         if (util_print('debug',print_debug) )then
            if (.not. geom_print(geom)) 
     &         call errquit('grad_dft: geom_print ?',0)
            if (.not. bas_print(ao_bas_han)) 
     &         call errquit('grad_dft: bas_print ?',0)
         endif
      endif
c
      call dftg_force(rtdb)
c
      call schwarz_tidy ()
c
c     grad_dft is done destroy basis and geometry handles  
c     (e.g., preserve the memory available to other modules!!)
c
      if(CDFIT)then
        if (.not. bas_destroy(cd_bas_han))
     &     call errquit('grad_dft:not able to destroy CD_bas:',86)
      endif
      if(XCFIT)then
        if (.not.bas_destroy(XC_bas_han))
     &     call errquit('grad_dft:not able to destroy XC_bas:',86)
      endif
      if (.not.((bas_destroy(ao_bas_han)) .and. (geom_destroy(geom)) ) )
     &   call errquit
     &   ('grad_dft:error destroying geom and ao_bas_han handles',911)
c
c     terminate integral scope
c
      call intd_terminate()
      call int_app_unset_no_texas(rtdb)
c
      call ecce_print_module_exit ('dft','ok')
      call util_print_pop
      grad_dft = .true.
c
      return
      end
