      logical function gradients(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
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
      logical status
      character*255 title
      double precision charge
c
      character*255 movecs
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
     $     ('grad: no converged SCF wavefunction available', 0)

c     Extract high level info from the data-base setting defaults
c
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      status = rtdb_get(rtdb, 'basis', MT_INT, 1, basis)
      if (.not.status) call errquit('scf: error storing basis handle',0)
      status = rtdb_get(rtdb, 'geom', MT_INT, 1, geom)
      if (.not.status) call errquit('scf: error storing geom handle', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('gradients: no geometry ', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $     call errquit('gradients: no mo basis ', 0)
      if(.not. int_normalize(basis))
     $     call errquit('gradients: normalization failed',911)

      movecs = 'movecs'
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
         if (.not. geom_print(geom)) 
     $        call errquit('gradients: geom_print ?',0)
         write(6,1) charge, 'rhf'
 1       format(/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/)
         if (.not. bas_print(basis)) 
     $        call errquit('gradients: bas_print ?',0)
         call util_flush(6)
      endif
c
      call force(rtdb,basis,geom,movecs)
c
      if (.not. context_pop('gradients')) 
     $     call errquit('context_pop failed',0)
      gradients = .true.
c
      return
      end
