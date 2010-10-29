      block data Dscfps
#include "cscfps.fh"
C$Id$
      data oscfps /.false./     ! So that other codes do not break
c
      end
      subroutine scf_pstat_init(rtdb)
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "cscf.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
      integer rtdb
      integer junk
      external dscfps ! For T3D linker
c
      if (.not. rtdb_get(rtdb, 'scf:pstat',
     $     MT_LOG, 1, oscfps)) oscfps = .false.
c
      if (.not. oscfps) return
c
      if (.not. pstat_allocate('scf: solve', pstat_qstat, 0,junk,
     $     ps_solve)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: search', pstat_qstat, 0, junk,
     $     ps_search)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: fock-2e', pstat_qstat, 0, junk,
     $     ps_fock_2e)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: fock-add', pstat_qstat, 0, junk,
     $     ps_fock_add)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: fock-io', pstat_qstat, 0, junk,
     $     ps_fock_io)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: int-1e', pstat_qstat, 0, junk,
     $     ps_int_1e)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: int-2e', pstat_qstat, 0, junk,
     $     ps_int_2e)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: ga ops', pstat_qstat, 0, junk,
     $     ps_gaops)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: matexp', pstat_qstat, 0, junk,
     $     ps_matexp)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: orthog', pstat_qstat, 0, junk,
     $     ps_orthog)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: diag', pstat_qstat, 0, junk,
     $     ps_diag)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c     
      if (.not. pstat_allocate('scf: vecio', pstat_qstat, 0, junk,
     $     ps_vecio)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: intio', pstat_qstat, 0, junk,
     $     ps_intio)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: dipole', pstat_qstat, 0, junk,
     $     ps_int_dip)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: sym ad', pstat_qstat, 0, junk,
     $     ps_sym_ad)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: sym   ', pstat_qstat, 0, junk,
     $     ps_sym_sym)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: k2cf  ', pstat_qstat, 0, junk,
     $     ps_k2cf)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: mullik', pstat_qstat, 0, junk,
     $     ps_mull)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: multip', pstat_qstat, 0, junk,
     $     ps_multip)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: moanal', pstat_qstat, 0, junk,
     $     ps_moanal)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
      if (.not. pstat_allocate('scf: total', pstat_qstat, 0, junk,
     $     ps_total)) call errquit('scf: scf_pstat_init', 0, MEM_ERR)
c
C     additional for RI-SCF
      if ( nriscf.ne.0 ) then
c       
        if (.not. pstat_allocate('riscf: Vhalf', pstat_qstat, 0,junk,
     $       ps_vhalf)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: index', pstat_qstat, 0, junk,
     $       ps_index)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: transform', pstat_qstat,0,junk,
     $       ps_transform)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: extract', pstat_qstat, 0, junk,
     $       ps_extract)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: exchange', pstat_qstat, 0,junk,
     $       ps_exchange)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: coulomb', pstat_qstat, 0, junk,
     $       ps_coulomb)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: transp', pstat_qstat, 0, junk,
     $       ps_transpose)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
c       
        if (.not. pstat_allocate('riscf: copy', pstat_qstat, 0, junk,
     $       ps_copy)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)

        if (.not. pstat_allocate('riscf: read', pstat_qstat, 0, junk,
     $       ps_riscf_read)) call errquit('riscf: ri_scf_pstat_init', 0,
     &       MEM_ERR)
      endif
c       
      call pstat_on(ps_total)
      oscfps = .true.
c
      end
      subroutine scf_pstat_print()
      implicit none
#include "errquit.fh"
#include "pstat.fh"
#include "cscfps.fh"
#include "cscf.fh"
c     
      if (oscfps) then
         call pstat_off(ps_total)
         call pstat_print_all
         oscfps = .false.
c
c     delete the SCF timers
c
         if (.not. pstat_free(ps_solve))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_search))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_fock_2e))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_int_1e))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_int_2e))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_gaops))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_fock_add))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_matexp))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_orthog))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_diag))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_vecio))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_int_dip))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_sym_ad))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_sym_sym))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_k2cf))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_fock_io))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_moanal))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_mull))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_multip))call errquit('scf_pstat?', 0,
     &       MEM_ERR)
         if (.not. pstat_free(ps_total))call errquit('scf_pstat?', 0,
     &       MEM_ERR)

C        additional for RI-SCF
         if ( nriscf.ne.0 ) then
           if (.not. pstat_free(ps_vhalf))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_index))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_transform))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_extract))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_exchange))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_coulomb))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_copy))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_transpose))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
           if (.not. pstat_free(ps_riscf_read))
     $          call errquit('scf_pstat?', 0, MEM_ERR)
         endif
c
      endif
c
      end
