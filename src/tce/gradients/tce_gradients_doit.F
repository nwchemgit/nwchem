      logical function tce_gradients_doit(rtdb,model,
     1        d_1pdm, k_1pdm_offset, size_1pdm,
     1        d_1pwdm,k_1pwdm_offset,size_1pwdm,
     1        d_2pdm, k_2pdm_offset, size_2pdm,
     1        noab,nvab,k_range,k_offset,k_spin,
     1        k_movecs_sorted,restricted)
c
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
c
      integer rtdb         ! [input] database handle
      character*10 model
      integer d_1pdm, k_1pdm_offset, size_1pdm
      integer d_1pwdm,k_1pwdm_offset,size_1pwdm
      integer d_2pdm, k_2pdm_offset, size_2pdm
      integer noab,nvab,k_range,k_offset,k_spin
      integer k_movecs_sorted
      logical restricted
c
      logical int_normalize
      external int_normalize

      integer geom, basis       ! handles
      logical status
      character*255 title
      double precision charge
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
     $     call errquit('tce_gradients_doit: geom_create?', 0, GEOM_ERR)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('tce_gradients_doit: no geometry ', 0, GEOM_ERR)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('tce_gradients_doit: bas_create?', 0, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('tce_gradients_doit: no ao basis ',0,BASIS_ERR)
      if(.not. int_normalize(rtdb,basis))
     $     call errquit('tce_gradients_doit: normalization failed',911, 
     1     INT_ERR)
c
c     Figure out the numer of electrons from the required total
c     charge and the sum of nuclear charges
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
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
     $           call errquit('tce_gradients_doit: geom_print ?',0, 
     1           GEOM_ERR)
            call util_flush(luout)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis)) 
c            if (.not. bas_print(ao_bas_han))
     $           call errquit('tce_gradients_doit: bas_print ?',0,
     1           BASIS_ERR)
            call util_flush(luout)
         endif
      endif
c
c     go for it ... finally ...
c
       call tce_grad_force(rtdb,model,basis,geom,
     1      d_1pdm, k_1pdm_offset,size_1pdm,
     1      d_1pwdm,k_1pwdm_offset,size_1pwdm,
     1      d_2pdm, k_2pdm_offset,size_2pdm,
     1      noab,nvab,k_range,k_offset,k_spin,
     1      k_movecs_sorted,restricted)

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
     &    ('tce_gradients_doit:error destroying geom and basis handles',
     1    911,GEOM_ERR)
c
      call ecce_print_module_exit ('gradients','ok')
      tce_gradients_doit = .true.
c
      return
      end
