      logical function gradients(rtdb)
c$Id: gradients.F,v 1.12 1995-03-31 01:39:31 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
#include "util.fh"
c
c     gradients module.
c
c     Context is '...:gradients'
c
c     Assumes SCF (currently only RHF) has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c
      logical int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status
      character*255 title
      double precision charge
c
      character*255 movecs_in, movecs_out
c
c     Push context down to GRADIENTS
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      if (.not. context_push('gradients'))
     $     call errquit('gradients: context_push failed', 0)
      call util_print_push
      call util_print_rtdb_load(rtdb,'gradients')
c
      status = rtdb_get(rtdb, 'scf:converged', MT_LOG, 1, converged) 
      if ( .not.( status .and. converged )) call errquit
     $     ('gradients: no converged SCF wavefunction available', 0)
c
c     Extract high level info from the data-base setting defaults
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0)
      if (.not. bas_create(basis, 'ao basis'))
     $     call errquit('gradients: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis'))
     $     call errquit('gradients: no ao basis ', 0)
      if(.not. int_normalize(basis))
     $     call errquit('gradients: normalization failed',911)

c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'scf:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
            movecs_out = 'movecs'
         else
            movecs_out = movecs_in
         endif
      endif
c
c     Figure out the numer of electrons from the required total
c     charge and the sum of nuclear charges
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
c
      if (ga_nodeid().eq.0) then
         if (util_print('information',print_low)) then
            call util_print_centered(6,'NWCHEM GRADIENTS Module',
     $           40,.true.)
            write(6,*)
            write(6,*)
            if (title .ne. ' ') then
               call util_print_centered(6, title, 40, .false.)
               write(6,*)
               write(6,*)
            endif
            write(6,1) charge, 'rhf'
 1          format(/
     $           '  charge          = ', f6.2/
     $           '  wavefunction    = ', a/)
            call util_flush(6)
         endif
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom)) 
     $           call errquit('gradients: geom_print ?',0)
            call util_flush(6)
         endif
         if (util_print('basis',print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('gradients: bas_print ?',0)
            call util_flush(6)
         endif
      endif
c
c initialize for schwarz screening
      call intd_init(rtdb,1,basis)
      call schwarz_init (geom, basis)
c
      call force(rtdb,basis,geom,movecs_out)
c
      call schwarz_tidy ()
c
c gradients is done destroy basis and geomtry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if (.not.(
     &    (bas_destroy(basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('gradients:error destroying geom and basis handles',911)
c
c terminate integral scope
c
      call int_terminate()
c
      call util_print_pop
      if (.not. context_pop('gradients')) 
     $     call errquit('context_pop failed',0)
      gradients = .true.
c
      return
      end
