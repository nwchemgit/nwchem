      logical function scf(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "inp.fh"
c
c     RHF module.
c
c     Context is '...:scf'
c
c     Currently performs closed shell ground state RHF calculation
c     but will at some point be extended to include ROHF and UHF, locking
c     of orbitals to converge excited states ..
c
c     Parameters passed in the rtdb
c
c       'geometry' -> name of geometry to use for the calculation
c                     (defaults to 'geometry')
c
c       'mo basis' -> name of basis set to use for the calculation
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
c       'convergence' -> required convergence threshold ... largest element
c                        of gradient vector (defaults to 10^-4).
c
c       'input mo vector file' -> name of file to look for input mo vectors
c                             (defaults to value of 'mo vector file')
c
c       'output mo vector file' -> name of file to look for output mo vectors
c                               (defaults to value of 'input mo vector file')
c
c       'mo vector file' -> generic source/sink of mo vectors
c                           (defaults to <prefix>.movec)
c
c       'energy' -> final energy (stores final converged energy)
c
c       'converged' -> true/false if calculation converged OK
c
c       'integral accuracy' -> Screening threshold for integrals.
c                              Defaults to being determined from 
c                              'wavefunction accuracy'
c
c       'wavefunction convergence' -> Convergence thre
c
c       'wavefunction accuracy' -> text description of required accuracy
c                                  for the wavefunction.  This is a
c                                  parameter that is understood by all
c                                  modules.
c                               Low -> result is suitable for starting
c                                      guess for subsequent calculation
c                               Medium -> Energy converged to at least
c                                         10^-6 au
c                               High -> Energy accurate to accuracy of the
c                                       integrals and wavefunction accurate 
c                                       to at least 1 part in 10^-6
c                   
c
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
c       Other options (fock-build method, ...)
c       are managed in the lower subroutines
c
      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status
      character*255 name, title
      double precision energy
      double precision charge
c
      character*255 movecs_in, movecs_out
      logical rhf
      external rhf
c
c     Push context down to SCF
c
      if (.not. context_push('scf'))
     $     call errquit('scf: context_push failed', 0)
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c     Extract high level info from the data-base setting defaults
c
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf: no geometry ', 0)
      if (.not. bas_create(basis, 'mo basis'))
     $     call errquit('scf: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $     call errquit('scf: no mo basis ', 0)
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
c     !! NEED to change context_rtdb_match -> context_rtdb_translate
c     !! and also provide context_rtdb_locate
c
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
c
c     Eventually print level will be under input control. For now
c     print verbosely
c 
      if (ga_nodeid() .eq. 0) then
         call util_print_centered(6, 'NWCHEM SCF Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
*     if (.not. geom_print(geom)) call errquit('scf: geom_print ?',0)
*     if (.not. bas_print(basis)) call errquit('scf: bas_print ?',0)
         write(6,1) charge, 'rhf',
     $        movecs_in(1:inp_strlen(movecs_in)),
     $        movecs_out(1:inp_strlen(movecs_out))
 1       format(/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/
     $        '  input vectors   = ', a/
     $        '  output vectors  = ', a/)
         call util_flush(6)
      endif
c     
c     When do other than RHF need to check that WFN type is valid
c     in case have inherited the wavefunction description from
c     an upper context that has not been correctly overridden
c
c     Here branch to the appropriate wavefunction type ... returns
c     true if happily converged, false otherwise (e.g., if restart
c     is needed or problems encountered).
c
      converged = rhf(rtdb, geom, basis, charge, movecs_in, movecs_out,
     $     energy)
c
c     Store energy and convergence status ... not sure if better to
c     do this here or have each wavefunction type do it
c
      if (.not. context_prefix('energy', name)) call errquit
     $     ('scf: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_DBL, 1, energy)) call errquit
     $     ('scf: failed to store energy in rtdb', 0)
      if (.not. context_prefix('converged', name)) call errquit
     $     ('scf: context_prefix failed', 0)
      if (.not. rtdb_put(rtdb, name, MT_LOG, 1, converged)) call errquit
     $     ('scf: failed to store converged in rtdb', 0)
c
c     Overwrite the input movecs option with where the final movecs
c     ended up so that a restart automatically picks them up
c
      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, movecs_out))
     $     call errquit('scf:rtdb_cput failed scf:input vectors',0)
c
c     SCF is done destroy basis and geomtry handles  
c
      if (.not.(bas_destroy(basis) .and. geom_destroy(geom)))
     $     call errquit ('scf:destroying geom and basis handles',0)
c
c terminate integrals 
c
      call int_terminate
c
c     Nothing else for now ... return convergence status
c
      if (.not. context_pop('scf')) call errquit('context_pop failed',0)
      scf = converged
c
      end
