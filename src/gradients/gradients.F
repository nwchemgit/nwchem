      logical function gradients(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "schwarz.fh"
c
c     gradients module.
c
c     Context is '...:gradients'
c
c     Assumes SCF (currently only RHF) has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c     Parameters passed in the rtdb
c
c       'geom' -> handle to geometry to use for the calculation
c                     (no default)
c
c       'basis' -> handle to basis set to use for the calculation
c                     (defaults to 'mo basis')
c
c       'scftype' -> name of function to use (rhf,rohf,uhf)
c                    (defaults to rhf currently)
c                    (more detailed info will be needed for rohf and uhf)
c                    (open shells, coupling, multiplicity, ...)
c
c       'charge' -> net charge on the system (double)
c                   (defaults to zero)
c
c       'title' -> title for the calculation
c                  (defaults to ' ')
c
c       'mo vector file' -> generic source/sink of mo vectors
c                           (defaults to <prefix>.movec)
c                           at the moment just 'movecs'
c
c       'converged' -> true/false if  SCF calculation converged OK
c
c       'integral accuracy' -> Screening threshold for integrals.
c                              Defaults to being determined from 
c                              'wavefunction accuracy'
c
c       not yet implemented for gradients:
c       'print level' -> text description of required output level.
c                        This is a parameter that is understood by all modules.
c                        None   -> no output except to report problems
c                        Low    -> Module name + 
c                                  brief summary of results
c                        Medium -> Module name +
c                                  title + 
c                                  brief description of the calculation + 
c                                  brief convergence information +
c                                  detailed summary of results + 
c                                  brief timing information
c                        High   -> Module name + 
c                                  title + 
c                                  detailed description of calculation + 
c                                  detailed convergence information + 
c                                  detailed listing of results +
c                                  detailed timing information
c                        Debug  -> Print anything appropriate to finding bugs
c
c
      logical int_normalize
      external int_normalize

      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status, print_debug
      character*255 title
      double precision charge
c
      character*255 movecs_in, movecs_out
      data print_debug /.false./
c
c     Push context down to GRADIENTS
c
      if (.not. context_push('gradients'))
     $     call errquit('gradients: context_push failed', 0)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
      status = rtdb_get(rtdb, 'scf:converged', MT_LOG, 1, converged) 
      if ( .not.( status .and. converged )) call errquit
     $     ('gradients: no converged SCF wavefunction available', 0)

c     Extract high level info from the data-base setting defaults
c
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('gradients: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0)
      if (.not. bas_create(basis, 'mo basis'))
     $     call errquit('gradients: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $     call errquit('gradients: no mo basis ', 0)
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
c     Eventually print level will be under input control. For now
c     print verbosely
c 
      if (ga_nodeid() .eq. 0) then
         call util_print_centered(6,'NWCHEM GRADIENTS Module',40,.true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
         write(6,1) charge, 'rhf'
 1       format(/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/)
         if (print_debug) then
           if (.not. geom_print(geom)) 
     $         call errquit('gradients: geom_print ?',0)
           if (.not. bas_print(basis)) 
     $         call errquit('gradients: bas_print ?',0)
         endif
         call util_flush(6)
      endif
c
c initialize for schwarz screening
      call intd_init(1,basis)
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
      if (.not. context_pop('gradients')) 
     $     call errquit('context_pop failed',0)
      gradients = .true.
c
      return
      end
