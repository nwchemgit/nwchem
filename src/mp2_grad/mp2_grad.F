      logical function semi_dir_mp2(rtdb)
*
* $Id: mp2_grad.F,v 1.24 2000-02-04 20:07:42 mg201 Exp $
*
      implicit none
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "schwarz.fh"
#include "util.fh"
#include "cmp2ps.fh"
      integer rtdb 
      integer geom, basis
      character*255 title
      character*255 movecs_in
      logical status
      character*8 task
      double precision threshscf
c
      logical cphf, scf, int_normalize
      external cphf, scf, int_normalize
c
      status=rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb,'mp2')
      call ecce_print_module_entry('mp2')
c
c     Make sure the SCF is being converged accurately enough
c
      if (.not. rtdb_get(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $     threshscf = 1.0d-4
      if (threshscf .gt. 1.0d-6) then
         threshscf = 1.0d-6
         if (.not. rtdb_put(rtdb,'scf:thresh',mt_dbl,1,threshscf))
     $        call errquit('mp2_grad: failed setting scf:thresh',0)
         status = rtdb_delete(rtdb,'scf:converged')
      endif
c
      if (.not. scf(rtdb)) then
         semi_dir_mp2 = .false.
         call ecce_print_module_exit('mp2', 'failed')
         Call Util_Print_pop
         return
      endif
c
      call mp2_pstat_init()
c
      if(.not.rtdb_cget(rtdb,'sdmp2:task',1,task)) call errquit
     $     ('could not determine task for semi-direct mp2',1)
      status = rtdb_delete(rtdb,'sdmp2:task')
c      
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
c
      if(ga_nodeid().eq.0)then
         if(util_print('information',print_low))then
            call util_print_centered(6,
     $           'NWChem MP2 Semi-direct Energy/Gradient Module',
     $           40,.true.)
            write(6,*)
            write(6,*)
            if(title.ne.' ')then
               call util_print_centered(6, title, 40, .false.)
               write(6,*)
               write(6,*)
            endif
            call util_flush(6)
         endif
         if(util_print('geometry',print_high))then
            if(.not.geom_print(geom)) 
     $           call errquit('mp2_grad: geom_print ?',0)
            call util_flush(6)
         endif
         if(util_print('basis',print_high))then
            if(.not.bas_print(basis)) 
     $           call errquit('mp2_grad: bas_print ?',0)
            call util_flush(6)
         endif
      endif
      if (.not.int_normalize(rtdb,basis))
     &    call errquit('mp2_grad: int_normalize failed',0)
c
c     Make sure that we use the same vectors as the SCF
c
      if (.not. rtdb_cget(rtdb,'scf:input vectors',1,movecs_in)) 
     $     call errquit('mp2_grad: no SCF MO vectors',0)
c
      call mp2_energy(rtdb, geom, basis, movecs_in, task)
c
      if(.not.((bas_destroy(basis)).and.(geom_destroy(geom))))
     $    call errquit
     $    ('mp2_grad:error destroying geom or basis handles',0)
c
      call mp2_pstat_tidy()
c
      call util_print_pop
c
      semi_dir_mp2 = .true.
      call ecce_print_module_exit('mp2', 'ok')
c
      end
