      subroutine mp2_grad(rtdb)
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
      integer rtdb 
      integer geom, basis
      character*255 title
      character*255 movecs_in
      logical status,converged
      status=rtdb_parallel(.true.)
      if(.not.context_push('mp2_grad'))
     $     call errquit('mp2_grad: context push failed',0)
      call util_print_push
      call util_print_rtdb_load(rtdb,'mp2_grad')
      status=rtdb_get(rtdb,'scf:converged',mt_log,1,converged)
      if(.not.(status .and. converged)) call errquit
     $     ('mp2_grad: no converged scf wavefunction available',0)
      if(.not.rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if(.not.geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0)
      if(.not.geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0)
      if(.not.bas_create(basis, 'ao basis'))
     $     call errquit('gradients: bas_create?', 0)
      if(.not.bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('gradients: no ao basis ', 0)
      if (ga_nodeid().eq.0) then
         if (util_print('information',print_low)) then
            call util_print_centered(6,
     $           'NWCHEM MP2 semi-direct energy/gradient module',
     $           40,.true.)
            write(6,*)
            write(6,*)
            if (title .ne. ' ') then
               call util_print_centered(6, title, 40, .false.)
               write(6,*)
               write(6,*)
            endif
            call util_flush(6)
         endif
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom)) 
     $           call errquit('mp2_grad: geom_print ?',0)
            call util_flush(6)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('mp2_grad: bas_print ?',0)
            call util_flush(6)
         endif
      endif
      call util_file_prefix('movecs', movecs_in)
      call int_init(rtdb,1,basis)
      call int_normalize(basis)
      call schwarz_init(geom,basis)
      call mp2_energy(rtdb,geom,basis,movecs_in)
      call schwarz_tidy()
      if(.not.((bas_destroy(basis)).and.(geom_destroy(geom))))
     &    call errquit
     &    ('mp2_grad:error destroying geom and basis handles',0)
      call int_terminate()
      call util_print_pop
      if(.not. context_pop('mp2_grad')) 
     $     call errquit('context_pop failed',0)
      end

















