c 
c $Id$ 
c 
      subroutine get_density(g_dens,g_vecs,theory,scftype,nalpha,
     *                       nbeta,nbfs,rtdb)
      implicit none
#include "errquit.fh"
c
c This routine gets alpha and beta densities and passes them back in g_dens.
c 
#include "mafdecls.fh"
#include "rtdb.fh"
      integer g_dens(*)    ! [output] density handles
      integer g_vecs(*)    ! [output] vector handles
      character*3 theory   ! [input]  HF or DFT
      character*4 scftype  ! [input]  RHF, UHF, or ROHF
      integer nalpha       ! [input]  number of alpha occupieds
      integer nbeta        ! [input]  number of beta occupieds
      integer nbfs         ! [input]  number of basis functions
      integer rtdb         ! [input]  RTDB handle
c
      integer l_evals, k_evals, l_occ, k_occ
      double precision rhffact
      character*255 movecs
c
      logical movecs_read
      external movecs_read
c
c Temporarily disable ROHF
c
      if (scftype.eq.'ROHF')
     *  call errquit('get_density: ROHF is not supported yet',555,
     &       CAPMIS_ERR)
c
      if (theory.eq.'HF') then
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs))
     *    call errquit('get_density: SCF MO vectors not defined',555,
     &       RTDB_ERR)
      else if (theory.eq.'DFT') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs))
     *   call errquit('get_density: DFT MO vectors not defined',555,
     &       RTDB_ERR)
      else
        call errquit('get_density:theory is unrecognized',555,
     &       INPUT_ERR)
      endif
c
c get additional space to store the eigenvalues and occupation numbers
c
      if (.not. ma_push_get(mt_dbl, nbfs,'MO evals', l_evals, k_evals))
     *  call errquit('get_density: could not allocate l_evals',nbfs,
     &       MA_ERR)
      if (.not. ma_push_get(mt_dbl, nbfs,'occ. numbers', l_occ, k_occ))
     *  call errquit('get_density: could not allocate l_occ',nbfs,
     &       MA_ERR)
c
c read vectors
c
      if (.not. movecs_read (movecs, 1, dbl_mb(k_occ), 
     *  dbl_mb(k_evals), g_vecs(1)))
     *  call errquit('get_density: could not read mo vectors', 555,
     &       DISK_ERR)
      if (scftype.eq.'UHF') then
        if (.not. movecs_read (movecs, 2, dbl_mb(k_occ), 
     *    dbl_mb(k_evals), g_vecs(2))) 
     *    call errquit('get_density: could not read beta vectors', 555,
     &       DISK_ERR)
      endif
c
c free additional space
c
      if (.not.ma_pop_stack(l_occ))
     *  call errquit('get_density:ma free occ',555, MA_ERR)
      if (.not.ma_pop_stack(l_evals))
     *  call errquit('get_density:ma free eval',555, MA_ERR)
c
c form the densities
c
      rhffact = 2.0
      if (scftype.eq.'UHF') rhffact = 1.0
      call ga_dgemm('n', 't', nbfs, nbfs, nalpha, rhffact,
     *   g_vecs(1), g_vecs(1), 0.0d0, g_dens(1))
      call ga_symmetrize(g_dens(1))
      if (scftype.eq.'UHF') then
        call ga_dgemm('n', 't', nbfs, nbfs, nbeta, rhffact,
     *    g_vecs(2), g_vecs(2), 0.0d0, g_dens(2))
        call ga_symmetrize(g_dens(2))
      endif
c
      end
