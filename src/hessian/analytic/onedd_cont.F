C
C $Id: onedd_cont.F,v 1.8 2000-07-18 15:00:57 windus Exp $
C
      subroutine onedd_cont(rtdb)
c
c  Add the one electron derivative (first and seconds) contributions to 
c  the hessian matrix and the RHS
c  Start with RHF and work our way up.
c
      implicit none
c
#include "bas.fh"
#include "global.fh"
#include "hess_info.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb  ! [input] RTDB handle
c
      integer i
      integer g_wdens
      integer g_eigen_diag
      integer max1e, max2e, mscratch_1e, mscratch_2e
      integer l_buf, k_buf, l_scr, k_scr
      integer matb, lsqatom, l_dens, k_dens, l_wdens, k_wdens
      integer l_evals, k_evals, l_occ, k_occ
      integer k_pdm2, k_coeff
      integer nbf2, nocc
      logical doS, doT, doV, frac_occ, omp2
c
      character*255 movecs_in                                                   
c
      integer ga_create_atom_blocked
      external ga_create_atom_blocked
c
c Get memory for the integral buffers
c
      call int_mem(max1e, max2e, mscratch_1e, mscratch_2e)
      if (.not. ma_push_get(mt_dbl,max1e,'intdd buffer',l_buf,k_buf))
     *  call errquit('onedd_cont:could not get intdd buffer',max1e)
      if (.not. ma_push_get(mt_dbl,mscratch_1e,'intdd scratch',
     *  l_scr, k_scr))call errquit('onedd_cont: scratch alloc failed',
     *  mscratch_1e)
c
c Need to form the weighted density matrix to multiply by the S 
c second derivative integrals.
c
      g_wdens = ga_create_atom_blocked(geom, basis, 'weighted density')
      nbf2 = nbf
      if (scftype .eq. 'UHF') nbf2 = nbf + nbf
C     eigenvalue diagonal matrix for forming energy weighted density
      nocc = nopen + nclosed                                                    
      if (.not. ga_create(mt_dbl, nocc, nocc,'eigen_diag', 0, nocc,
     $  g_eigen_diag)) call errquit('onedd_cont: ga diags?',nocc*nocc)
      call ga_zero(g_eigen_diag)                                                
C     eigenvalues
      if (.not. ma_push_get(mt_dbl, nbf2,'MO evals', l_evals, k_evals))
     $     call errquit('onedd_cont: could not allocate l_evals',nbf2)
C     occupation numbers
      if (.not. ma_push_get(mt_dbl, nbf2,'occ. numbers', l_occ, k_occ))
     $     call errquit('onedd_cont: could not allocate l_occ',nbf2)
c
c     check if smearing is on
c
      if (.not.rtdb_get(rtdb, 'dft:fractional_occup', mt_log, 1,
     & frac_occ  )) frac_occ=.false.
      omp2 = .false.
      k_pdm2 = 1
      k_coeff = 1
c
c get MO vectors
c
      if (theory.eq.'dft') then
        if (.not. rtdb_cget(rtdb, 'dft:input vectors', 1, movecs_in))
     *    call errquit('onedd_cont: DFT MO vectors not defined',555)
      else
        if (.not. rtdb_cget(rtdb, 'scf:input vectors', 1, movecs_in))
     $    call errquit('onedd_cont: SCF MO vectors not defined',0)           
      endif
c
c get density and weighted density
c
      call grad_dens(geom, basis, g_dens, g_wdens, g_eigen_diag,
     $     dbl_mb(k_evals), dbl_mb(k_occ), ndens, nbf, nbf2,
     $     nopen, nclosed, scftype, movecs_in, omp2,
     $     dbl_mb(k_pdm2), dbl_mb(k_coeff),frac_occ)                            
c
C     free temporary arrays
      if (.not. ga_destroy(g_eigen_diag))
     $  call errquit('onedd_cont: could not destroy g_eigen_diag',1)
 
      if (.not.ma_pop_stack(l_occ))
     $     call errquit('onedd_cont:ma free occ',1)
      if (.not.ma_pop_stack(l_evals))
     $     call errquit('onedd_cont:ma free eval',1)                                  
c
      write(6,*) 'The weighted density is'
      call ga_print(g_wdens)
      write(6,*) 'The density is'
      call ga_print(g_dens(1))
c
      if (.not. bas_nbf_ce_max(basis,matb))
     *  call errquit('onedd_cont: could not get matb',555)
      lsqatom = matb * matb
      if (.not. ma_push_get(mt_dbl,lsqatom,'piece_density',
     *  l_dens,k_dens))
     *  call errquit('onedd_cont:could not allocate l_dens',lsqatom)
      if (.not. ma_push_get(mt_dbl,lsqatom,'piece_wdensity',l_wdens,
     *     k_wdens)) call errquit
     *    ('onedd_cont:could not allocate l_wdens',lsqatom) 
c
c  Call the actual routine to add the S, V and T contributions
c

c  first print out the rhs values (should be all zeroes)!
      write(6,*) 'g_rhs before going into onedd_calc'
      do i = 1, ncent*3
        call ga_print(g_rhs(i))
      enddo
c
c     write(6,*) 'before onedd_calc'
c     call flush(6)
c
      doS = .true. 
      doV = .true.
      doT = .true.
      call onedd_calc(dbl_mb(k_exy),ncent*3,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,g_wdens,g_dens,geom,basis,ncent,dbl_mb(k_dens),
     *  dbl_mb(k_wdens),matb)
c
c should put the next line in when done with some testing
c     call sym_hess_symmetrize(geom,dbl_mb(k_exy))
      write(6,*) 'The Hessian after derivative S contributions:'
      call flush(6)
      call hnd_hssout(dbl_mb(k_exy),ncent*3,ncent*3,ncent*3)
c
c Get rid of the weighted density and the other temporary arrays
c
      if (.not. ma_pop_stack(l_wdens)) 
     *  call errquit('onedd_cont: failed to pop l_wdens',555)
      if (.not. ma_pop_stack(l_dens)) 
     *  call errquit('onedd_cont: failed to pop l_dens',555)
      if (.not. ga_destroy(g_wdens)) 
     *  call errquit('onedd_cont: failed to free g_wdens',555)
c
c Now handle the Hx and Sx terms for the rhs
c First handle Sx term

      doS = .true.
      doV = .false.
      doT = .false.
      call oned_calc(ncent*3,g_rhs,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,geom,basis,ncent)
c
      write(6,*) 'g_rhs for Sx term'
      do i = 1, ncent*3
        call ga_print(g_rhs(i))
      enddo
c
c write out Sx terms to disk
c
      call sx_write(g_rhs,ncent*3,nbf,sx_dra_handle)
c
c Clear g_rhs again
c
      do i = 1, ncent*3
        call ga_zero(g_rhs(i))
      enddo
c
      doS = .false.
      doV = .true.
      doT = .true.
c  first print out the rhs values 
      write(6,*) 'g_rhs before oned_calc'
      do i = 1, ncent*3
        call ga_print(g_rhs(i))
      enddo
      call oned_calc(ncent*3,g_rhs,dbl_mb(k_buf),
     *  max1e,dbl_mb(k_scr),mscratch_1e,skelsym,
     *  doS,doT,doV,geom,basis,ncent)
c  first print out the rhs values 
      write(6,*) 'g_rhs after oned_calc'
      do i = 1, ncent*3
        call ga_print(g_rhs(i))
      enddo
c
c  Get rid of the temporary buffers for the integrals
c
      if (.not. ma_pop_stack(l_scr)) 
     *  call errquit('onedd_cont: failed to pop l_scr',555)
      if (.not. ma_pop_stack(l_buf)) 
     *  call errquit('onedd_cont: failed to pop l_buf',555)
c
      return
      end
