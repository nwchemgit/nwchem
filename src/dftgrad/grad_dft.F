      logical function grad_dft(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "util.fh"
#include "schwarz.fh"
#include "stdio.fh"
c****
c**** dft gradients module hacked from the scf gradients module
c****
c
c     dftgrad module.
c
c     Context is '...:dftgrad'
c
c     Assumes SCF (currently only RHF) has been completed, MO vectors stored
C     and all information is still in the RTDB
c
c     Parameters passed in the rtdb
c
c       'geom' -> handle to geometry to use for the calculation
c                     (no default)
c
c       'ao_basis' -> handle to basis set to use for the calculation
c                     (defaults to 'ao basis')
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
      Logical int_normalize
      External int_normalize
      Logical int_norm_2c
      External int_norm_2c

      Integer rtdb              ! [input] database handle
      Integer geom, ao_basis,cd_basis       ! handles
      Integer ipol,nbases
      Logical converged, CDFIT, status
      Character*255 title
      Double precision charge
      integer bases(3)
c
      Character*255 movecs_in, movecs_out
c----------------------------------------------------------------------- 
c
c     Push context down to DFTGRAD
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      call util_print_push('dftgrad')
      call util_print_rtdb_load(rtdb,'dft_grad')
      call ecce_print_module_entry ('dft')

c**** debug
C      status = rtdb_print(rtdb,.true.)
c
c
      status = rtdb_get(rtdb, 'dft:converged', MT_LOG, 1, converged) 
      if ( .not.( status .and. converged )) then
        if ( (ga_nodeid() .eq. 0) 
     $       .and. util_print('information', print_none) ) then
          write (luout,*) 'status: ', status, '   converged: ',converged
        endif
        call errquit(
     &       'dft grad_dft: no converged DFT wavefunction available',
     &     0)
      endif 

c
c     Extract high level info from the data-base setting defaults
c
c     load DFT parameters into common
c     
c      call dft_rdinput(rtdb)
c
      nbases=1
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('grad_dft: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('grad_dft: no geometry ', 0)
      if (.not. bas_create(ao_basis, 'ao basis'))
     $     call errquit('grad_dft: bas_create?', 0)
        bases(1)=ao_basis
      if (.not. bas_rtdb_load(rtdb, geom, ao_basis, 'ao basis'))
     $     call errquit('grad_dft: no ao basis ', 0)
      if(.not. int_normalize(ao_basis))
     $     call errquit('grad_dft: normalization failed',911)
C
C     initialize CD basis (if exists)

      if (.not. bas_create(cd_basis, 'cd basis'))
     &     call errquit('grad_dft: bas_create?',102)
      CDFIT=bas_rtdb_load(rtdb, geom, cd_basis, 'cd basis')
      if(CDFIT) then
        nbases=nbases+1
        bases(nbases)=cd_basis
        if (.not. int_norm_2c(cd_basis))
     &       call errquit('grad_dft: int_norm_2c failed', 208)
        else
        if (.not. bas_destroy(cd_basis))
     &       call errquit('dft_grad:not able to destroy CD_bas:',86)
        cd_basis=-99999
      endif
c
c     Figure output input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c
      if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'dft:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore') then
            call util_file_name('movecs',.false.,.false.,movecs_out)
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
      if (.not. rtdb_get(rtdb, 'dft:ipol', MT_INT, 1, ipol ))
     $     call errquit('dft_grad: dft:ipol get failed', 0)
      if (ga_nodeid() .eq. 0) then
         if ( util_print('information',print_default) ) then
            call util_print_centered(luout,
     $          'NWCHEM DFTGRAD Module',40,.true.)
            write(luout,*)
            write(luout,*)
            if (title .ne. ' ') then
               call util_print_centered(luout, title, 40, .false.)
               write(luout,*)
               write(luout,*)
            endif

            if (ipol .eq. 1) then
               write(luout,1) charge, 'closed shell'
            else 
               write(luout,1) charge, 'open shell'
            endif 
            call util_flush(luout)

         endif

 1       format(/
     $        '  charge          = ', f6.2/
     $        '  wavefunction    = ', a/)
         if ( util_print('debug',print_debug) ) then
           if (.not. geom_print(geom)) 
     $         call errquit('grad_dft: geom_print ?',0)
           if (.not. bas_print(ao_basis)) 
     $         call errquit('grad_dft: bas_print ?',0)
         endif
      endif
c
c initialize NWints
c
      call intd_init(rtdb,nbases,bases)
c
c initialize for schwarz screening
c
      call schwarz_init (geom, ao_basis)
c
      call dftg_force(rtdb,ao_basis,cd_basis,geom, 
     &     CDFIT,movecs_out)
c
      call schwarz_tidy ()
c
c grad_dft is done destroy basis and geometry handles  
c (e.g., preserve the memory available to other modules!!)
c
      if(CDFIT) then
        if (.not. bas_destroy(cd_basis))
     &       call errquit('dft_grad:not able to destroy CD_bas:',86)
      endif
      if (.not.(
     &    (bas_destroy(ao_basis))
     &    .and.
     &    (geom_destroy(geom))
     &    ))
     &    call errquit
     &    ('grad_dft:error destroying geom and ao_basis handles',911)
c
c terminate integral scope
c
      call intd_terminate()
c
      call ecce_print_module_exit ('dft','ok')
      call util_print_pop
      grad_dft = .true.
c
      return
      end
