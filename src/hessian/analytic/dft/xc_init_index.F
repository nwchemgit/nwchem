      subroutine xc_init_index(geom,ipol,nbf_ao,ao_bas_han,
     .     g_dens, l_cntoce, k_cntoce, l_cntobfr,
     &                         k_cntobfr, l_cetobfr, k_cetobfr,
     &                         l_rdens_atom, k_rdens_atom)
c
c     $Id$
c
c     Sets up indexing and reduced density arrays used by several
c     XC routines
c
c     BGJ - 8/98
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "bas.fh"
      integer geom ,ipol ! [in]
      integer nbf_ao,ao_bas_han     ! [in]
      integer g_dens(*)         ! [input] GA handles for density matrices
      integer l_cntoce, k_cntoce, l_cntobfr, k_cntobfr,
     &        l_cetobfr, k_cetobfr, l_rdens_atom, k_rdens_atom
      integer nbf_ao_mxnbf_ce,nshells_ao
c
c
c     Local variables
c
      integer natoms, nscr, l_scr, k_scr
c
      if (.not. bas_nbf_ce_max(AO_bas_han,nbf_ao_mxnbf_ce))
     .     call errquit('xcinitindx: cemax failed ',0, BASIS_ERR)
      if (.not. bas_numcont(AO_bas_han,nshells_ao))
     .     call errquit('xcinitindx: numcont failed ',1, BASIS_ERR)
c
c     The following code was taken from dft_scf.F
c
c     Pre-compute mapping vectors
c 
      if (.not. geom_ncent(geom, natoms))
     &     call errquit('xc_init_index: geom_ncent failed',1, GEOM_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao,'cntoce map',l_cntoce,k_cntoce))
     &   call errquit('xc_init_index:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,nshells_ao*2,'cntoce map',l_cntobfr,k_cntobfr))
     &   call errquit('xc_init_index:push_get failed', 13, MA_ERR)
      if (.not.ma_push_get
     &   (mt_int,natoms*2,'cntoce map',l_cetobfr,k_cetobfr))
     &   call errquit('xc_init_index:push_get failed', 13, MA_ERR)
c     
      call build_maps(ao_bas_han, int_mb(k_cntoce), int_mb(k_cntobfr), 
     &                int_mb(k_cetobfr), natoms, nshells_ao)
c
c     Set aside some memory for reduced density matrix
c
      if (.not.MA_Push_Get(MT_Dbl,ipol*natoms*natoms,'rdens_atom',
     &   l_rdens_atom,k_rdens_atom))
     &   call errquit('xc_init_index: cannot allocate rdens_atom',0,
     &       MA_ERR)
c
c     Pre-compute reduced total density matrices over atoms
c 
      call dfill(ipol*natoms*natoms, 0.0d0, dbl_mb(k_rdens_atom), 1)
      nscr = nbf_ao_mxnbf_ce*nbf_ao_mxnbf_ce
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr',l_scr,k_scr))
     &   call errquit('dft_scf: cannot allocate scr',0, MA_ERR)
      call util_ga_mat_reduce(nbf_ao, natoms, int_mb(k_cetobfr), g_dens, 
     &                        ipol, dbl_mb(k_rdens_atom), 'rms', 
     &                        dbl_mb(k_scr), nbf_ao_mxnbf_ce,.true.)
      if (.not.ma_pop_stack(l_scr))
     &   call errquit('xc_init_index: cannot pop stack',0, MA_ERR)
c      write(*,*)'rdens_atom'
c      call output(dbl_mb(k_rdens_atom),1,natoms,1,natoms,natoms,natoms,
c     &     1)
c!!!  Until this is going properly, put values in which won't trigger
c!!!  cutoffs
      call dfill(ipol*natoms*natoms, 10.0d0, dbl_mb(k_rdens_atom), 1)
c
      return
      end
