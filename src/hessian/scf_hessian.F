      logical function scf_hessian(rtdb)
* $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "stdio.fh"
#include "cfock.fh"
c
c     hessian module.
c
c     Assumes SCF has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c
      logical  int_normalize, scf
      logical  hess_anal
      external hess_anal
      external int_normalize, scf

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status
      character*255 title
c
c temp
c
      status = hess_anal(rtdb)
      if (.not. status) call errquit
     $  ('analytic hessian did not work',555, CALC_ERR)
      scf_hessian = status
      return
c
      if (.not. scf(rtdb))
     $     call errquit('scf_hessian: scf failed', 0, CALC_ERR)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push
c     call util_print_rtdb_load(rtdb,'hessian')

c
      status = rtdb_get(rtdb, 'scf:converged', MT_LOG, 1, converged) 
      if ( .not.( status .and. converged )) call errquit
     $     ('hessian: no converged SCF wavefunction available', 0,
     &       CALC_ERR)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('hessian: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('hessian: no geometry ', 0, RTDB_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('hessian: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('hessian: no ao basis ', 0, RTDB_ERR)
      if(.not. int_normalize(rtdb,basis))
     $     call errquit('hessian: normalization failed',911, INT_ERR)
c
      if (ga_nodeid().eq.0) then
         if (util_print('information',print_low)) then
            call util_print_centered(luout,'NWChem Hessian Module',
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
     $           call errquit('hessian: geom_print ?',0, GEOM_ERR)
            call util_flush(luout)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis))

     $           call errquit('hessian: bas_print ?',0, BASIS_ERR)
            call util_flush(luout)
         endif
      endif
c
c modify the code to pass both basis and geom
c
      call scf_hnd_dder(rtdb,geom,basis)
c
c hessian is done destroy basis and geometry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('hesian:error destroying geom and basis handles',911,
     &       GEOM_ERR)
c
      call util_print_pop
      scf_hessian = .true.
c
      return
      end
