      block data cscfdata
      implicit none
#include "cscf.fh"
c
      data oinitialized /.false./
      data g_movecs /0/
      data k_occ /10000000/
      data k_eval /10000000/
c
      end
      subroutine scf_init(rtdb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "sym.fh"
#include "util.fh"
#include "cscf.fh"
c     
      integer rtdb              ! database handle
c     
      double precision nuclear_charge
      character*255 name
c     
      if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $     title = ' '
c
c     load geometry and symmetry info
c
      if (.not. geom_create(geom, 'geometry'))
     $     call errquit('scf_init: geom_create?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, 'geometry'))
     $     call errquit('scf_init: no geometry ', 0)
      if (.not.rtdb_get(rtdb, 'scf:skeleton',MT_LOG, 1, oskel)) then
         oskel = sym_number_ops(geom) .gt. 0
      endif
c     
c     load the basis set and get info about it
c
      if (.not. bas_create(basis, 'mo basis'))
     $     call errquit('scf_init: bas_create?', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, 'ao basis')) then
         if (.not. bas_rtdb_load(rtdb, geom, basis, 'mo basis'))
     $        call errquit('scf_init: no mo or ao basis set', 0)
      end if
      if (.not. bas_name(basis, name, trans))
     $     call errquit('scf_init: bas_name?', 0)
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('scf_init: basis info',0)
c
c     For now just set NMO=NBF but need to detect linear dependence
c
      nmo = nbf
c
c     Figure input/output MO vectors ... all this should be
c     context sensitive but that design is incomplete.
c     
      if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $     movecs_in = 'atomic'
      if (.not. rtdb_cget(rtdb, 'scf:output vectors', 1, movecs_out))
     $     movecs_out = ' '
      if (movecs_out.eq.' ') then
         if (movecs_in.eq.'atomic' .or. movecs_in.eq.'hcore' .or. 
     $        movecs_in.eq.'project') then
            call util_file_prefix('movecs', movecs_out)
         else
            movecs_out = movecs_in
         endif
      endif
c     
c     Figure out the number of electrons from the required total
c     charge and the sum of nuclear charges
c     
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      if (.not. geom_nuc_charge(geom, nuclear_charge))
     $     call errquit('scf: geom_nuc_charge failed', 0)
      nelec = nint(nuclear_charge - charge)
      if (nelec .le. 0) call errquit
     $     ('scf: negative no. of electrons ?', nelec)
      if (abs(nuclear_charge - charge - dble(nelec)) .gt. 1d-8)
     $     call errquit('scf: non-integral no. of electrons ?', 0)
c     
c     Determine no. of open and closed shells ... default is to run closed
c     shell unless told otherwise
c     
      if (.not. rtdb_get(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     nopen = 0
      if (nopen .gt. nelec) call errquit
     $     ('scf: nopen > nelec ', nopen)
      if (mod(nelec-nopen,2) .ne. 0) call errquit
     $     ('scf: no. of closed-shell electrons is not even!',nopen)
      nclosed = (nelec-nopen) / 2
c     
      if (.not. rtdb_cget(rtdb, 'scf:scftype', 1, scftype)) then
         if (nopen .eq. 0) then
            scftype = 'RHF'
         else
            scftype = 'ROHF'
         endif
      endif
      if (scftype.ne.'ROHF' .and. scftype.ne.'RHF') call errquit
     $     ('scf: only ROHF and RHF currently supported', 0)
c
c     Store the derived info in the database for other wavefunction
c     modules and/or restart to access
c
      if (.not. rtdb_cput(rtdb, 'scf:scftype', 1, scftype)) 
     $     call errquit('scf_init: put of scftyp failed',0)
      if (.not. rtdb_put(rtdb, 'scf:nopen', MT_INT, 1, nopen))
     $     call errquit('scf_init: put of nopen failed',0)
      if (.not. rtdb_put(rtdb, 'scf:nclosed', MT_INT, 1, nclosed))
     $     call errquit('scf_init: put of nclosed failed',0)
      if (.not. rtdb_put(rtdb, 'scf:nelec', MT_INT, 1, nelec))
     $     call errquit('scf_init: put of nelec failed',0)
c     
c     Allocate persistent local and global arrays
c
      if (.not. ga_create(MT_DBL, nbf, nbf, 'scf_init: MOs',
     $     1, 1, g_movecs)) call errquit('scf_init: MOs', 0)
c     
      if (.not. ma_push_get(mt_dbl, nbf, 'scf_init: mo evals',
     $     l_eval, k_eval)) call errquit
     $     ('rihf_init: insufficient memory?',nbf)
c     
      if (.not. ma_push_get(mt_dbl, nbf, 'scf_init: mo occ',
     $     l_occ, k_occ)) call errquit
     $     ('rihf_init: insufficient memory?',nbf)
c
c     Fill in the SCF convergence info
c     
      call scf_get_conv_info(rtdb)
      call scf_get_fock_param(rtdb)
c
      oinitialized = .true.
c     
      end
      subroutine scf_get_conv_info(rtdb)
C$Id: scf_get_info.F,v 1.2 1995-04-03 23:37:49 d3g681 Exp $
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "cscf.fh"
c     
      integer rtdb
c
      double precision mone
      parameter (mone = -1.0d0)
c     
      if (.not. rtdb_get(rtdb, 'scf:maxiter',
     $     MT_INT, 1, maxiter)) maxiter = 8
      if (.not. rtdb_get(rtdb, 'scf:tol2e',
     $     MT_DBL, 1, tol2e)) tol2e = 1.0d-7
      if (.not.rtdb_get(rtdb, 'scf:thresh',
     $     MT_DBL, 1, gnorm_tol)) gnorm_tol = 1.d-5
      if (.not.rtdb_get(rtdb, 'scf:level shift info',
     $     MT_DBL, 6, shifts)) call dfill(6, mone, shifts, 1)
      if (.not.rtdb_get(rtdb, 'scf:full hessian switch',
     $     MT_DBL, 1, nr_gswitch)) nr_gswitch = 0.1d0
      if (.not.rtdb_get(rtdb, 'scf:disable conjugacy',
     $     MT_LOG, 1, odisable_conjugacy)) odisable_conjugacy = .false.
c
c     Apply defaults
c
      if (shifts(1) .eq. -1.0d0) shifts(1) = 5.0d0
      if (shifts(2) .eq. -1.0d0) shifts(2) = 0.5d0
      if (shifts(3) .eq. -1.0d0) shifts(3) = 0.0d0
      if (shifts(4) .eq. -1.0d0) shifts(4) = 0.0d0
      if (shifts(5) .eq. -1.0d0) shifts(5) = 0.0d0
      if (shifts(6) .eq. -1.0d0) shifts(6) = 0.0d0
c
      end
      subroutine scf_get_fock_param(rtdb)
      implicit none
#include "mafdecls.fh"
#include "rtdb.fh"
#include "cfock.fh"
#include "global.fh"
#include "util.fh"
      integer rtdb
      integer chunktmp
c
      if (rtdb_get(rtdb, 'scf:task_chunks',
     $     MT_INT, 1, chunktmp)) then
         task_chunks = chunktmp         
         if (ga_nodeid().eq.0 .and.
     $        util_print('information',print_default)) then
            write(6,2) task_chunks
 2          format(/'Setting fock-build task_chunks to', i4)
            call util_flush(6)
         endif
      endif
c
      end
      subroutine scf_tidy(rtdb)
      implicit none
#include "cscf.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
      integer rtdb
c
      logical status
c
      if (oinitialized) then
         if (.not. geom_destroy(geom)) call errquit
     $        ('scf_tidy: geom destroy failed', 0)
         if (.not. bas_destroy(basis)) call errquit
     $        ('scf_tidy: basis destroy failed',0)
         status = ma_pop_stack(l_occ)
         status = ma_pop_stack(l_eval) .and. status
         if (.not. status) call errquit
     $        ('scf_tidy: failed to free occupation/evals',0)
         if (.not. ga_destroy(g_movecs)) call errquit
     $        ('scf_tidy: failed to free movecs',0)
c
         oinitialized = .false.
c
      endif
c
      end
