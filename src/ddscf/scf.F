      logical function scf(rtdb)
C$Id: scf.F,v 1.22 1995-02-02 23:12:07 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
#include "inp.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "cskeleton.fh"
#include "sym.fh"
#include "util.fh"
c
c     RHF module.
c
c     Context is '...:scf'
c
c     Currently performs closed shell ground state RHF calculation
c     but will at some point be extended to include ROHF and UHF, locking
c     of orbitals to converge excited states ..
c
      integer rtdb              ! [input] database handle
      integer geom, basis       ! handles
      logical converged
      logical status
      character*255 name, title, trans
      double precision energy
      double precision charge
      logical oprint
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
c     Initialize performance statistics
c
      call scf_pstat_init(rtdb)
c
c     Extract high level info from the data-base setting defaults
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'scf')
      oprint = util_print('information',print_low)
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
      if (.not. bas_name(basis, name, trans))
     $     call errquit('scf: bas_name?', 0)
      if (.not.rtdb_get(rtdb, 'scf:skeleton',
     $     MT_LOG, 1, oskel)) then
        oskel = sym_number_ops(geom) .gt. 0
      endif
      if (name .eq. ' ') name = 'mo basis'
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'scf:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project') then
            movecs_out = 'movecs'
         else
            movecs_out = movecs_in
         endif
      endif
c
c     Figure out the number of electrons from the required total
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
      
      if (ga_nodeid() .eq. 0 .and. oprint) then
         call util_print_centered(6, 'NWCHEM SCF Module', 40, .true.)
         write(6,*)
         write(6,*)
         if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(6,*)
            write(6,*)
         endif
c
         write(6,1) trans(1:inp_strlen(trans)), charge, 'rhf',
     $        movecs_in(1:inp_strlen(movecs_in)),
     $        movecs_out(1:inp_strlen(movecs_out)), oskel
 1       format(/
     $        '  mo basis        = "',a,'"'/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/
     $        '  input vectors   = ', a/
     $        '  output vectors  = ', a/
     $        '  use symmetry    = ', l/)
         call util_flush(6)
      endif
c
c     Print out other junk
c
      if (ga_nodeid() .eq. 0) then
         if (util_print('geometry',print_high)) then
            if (.not. geom_print(geom))
     $        call errquit('scf: geom_print', 0)
            write(6,*)
         endif
         if (util_print('symmetry', print_debug)) then
            call sym_print_all(geom, .true., .true., .true., .true.)
            write(6,*)
         endif
         if (util_print('basis', print_high)) then
            if (.not. bas_print(basis)) 
     $           call errquit('scf: bas_print', 0)
            write(6,*)
         endif
         if (util_print('geombas', print_debug)) then
            if (.not. gbs_map_print(basis)) 
     $           call errquit('scf:gmp', 0)
            write(6,*)
         endif
         call util_flush(6)
      endif
c
c     Overwrite the input movecs option with where the final movecs
c     will end up so that a restart automatically picks them up
c
      if (.not. rtdb_cput(rtdb, 'scf:input vectors', 1, movecs_out))
     $     call errquit('scf:rtdb_cput failed scf:input vectors',0)
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
c     SCF is done destroy basis and geometry handles  
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
      call util_print_pop
      if (.not. context_pop('scf')) call errquit('context_pop failed',0)
      scf = converged
c
      call scf_pstat_print
c
      end
