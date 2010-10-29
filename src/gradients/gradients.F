      logical function gradients(rtdb)
c$Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "util.fh"
#include "stdio.fh"
c
c     gradients module.
c
c     Assumes SCF has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c
      logical int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical status
      character*255 title
c      double precision charge
c
c     cosmo flag
c
      logical odbug
      logical ocosmo
      logical osome
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call ecce_print_module_entry ('gradients')
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0, GEOM_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('gradients: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('gradients: no ao basis ', 0, BASIS_ERR)
      if(.not. int_normalize(rtdb,basis))
     $     call errquit('gradients: normalization failed',911, INT_ERR)
c
c     Figure out the numer of electrons from the required total
c     charge and the sum of nuclear charges
c
c      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
c     $     charge = 0.0d0
c
      if (ga_nodeid().eq.0) then
         if (util_print('information',print_low)) then
            call util_print_centered(luout,'NWChem Gradients Module',
     $           40,.true.)
            write(luout,*)
            call util_flush(luout)
          endif
          if (util_print('information',print_medium)) then
            write(luout,*)
            if (title .ne. ' ') then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
            endif
            call util_flush(luout)
         endif
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom)) 
     $           call errquit('gradients: geom_print ?',0, GEOM_ERR)
            call util_flush(luout)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('gradients: bas_print ?',0, BASIS_ERR)
            call util_flush(luout)
         endif
      endif
c
c     ----- cosmo check -----
c
      odbug=.false.
      odbug=odbug.and.ga_nodeid().eq.0
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
         if(odbug) then
            write(luout,*) '-cosmo- ... found in -gradients-',
     $                     ocosmo,ga_nodeid()
         endif
         if(ocosmo) then
            if(odbug) then
               osome=.true.
            else
               osome=.false.
            endif
            osome=osome.and.ga_nodeid().eq.0
            if(odbug) then
               write(luout,*) '-cosmo- ... found and .true. ',
     $                        ocosmo,ga_nodeid()
            endif
         else
            if(odbug) then
               write(luout,*) '-cosmo- ... found but .false.',
     $                        ocosmo,ga_nodeid()
            endif
         endif
      else
         if(odbug) then
            write(luout,*) '-cosmo- not found in -gradients-'
         endif
      endif
      call ga_sync()
c
c     go for it ... finally ...
c
      call grad_force ( rtdb, basis, geom )
c
c gradients is done destroy basis and geometry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('gradients:error destroying geom and basis handles',911,
     &       GEOM_ERR)
c
      call ecce_print_module_exit ('gradients','ok')
      gradients = .true.
c
      return
      end
