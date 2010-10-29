      subroutine tddft_transfm(ivector,g_vector,g_movecs,
     1  nbf_ao,nocc,nmo,nfc,nfv,ipol,g_density,itype)
c
c $Id$
c
c itype=1
c
c Get a (trial) vector from GA and transform it into a (trial)
c density matrix in an AO basis and return it as an array.
c
c itype=2
c
c Receive a (product) density matrix in an AO basis as an input
c and transform it into a (product) vector in MO basis and put
c it in a GA.
c
c ! Caution ! Never destroy ivector
c
c Written by So Hirata, Feb 2002. (c) Battelle, PNNL, 2002.
c
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
c
      integer itype            ! Subroutine function type =1 or =2
      integer ivector          ! ivector-th vector to be processed
      integer nbf_ao           ! Number of AO basis functions
      integer nocc(2)          ! Number of occupied orbitals
      integer nmo(2)           ! Number of orbitals
      integer nfc(2),nfv(2)    ! Number of orbitals
      integer ipol             ! =1 (restricted), =2 (unrestricted)
      integer g_vector(2)      ! GA handle for trial/product vectors
      integer g_movecs(2)      ! GA handle for MO vectors
      integer g_density(2)     ! GA handle for density matrix
c
      integer l_w1,k_w1        ! MA handle for a work array
      integer l_w2,k_w2        ! MA handle for a work array
      integer l_w3,k_w3        ! MA handle for MO vectors
      integer i,j,k,l
      integer g_dd,g_vpatch,l0,l1,l1l0,ilo,ihi
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c --------------------
c Allocate work arrays
c --------------------
c
c To avoid confusion and mistake in case when nmo(i) < nbf_ao,
c I choose to carry out full matrix multiplication in the space
c of nbf_ao*nbf_ao, while leaving the meaningless part of the 
c space zeroed.
c
      if (.not.ma_push_get(mt_dbl,nbf_ao**2,'Work array',
     1  l_w1,k_w1)) call errquit
     2  ('tddft_transfm: failed to allocate w1',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao**2,'Work array',
     1  l_w2,k_w2)) call errquit
     2  ('tddft_transfm: failed to allocate w2',0, MA_ERR)
      if (.not.ma_push_get(mt_dbl,nbf_ao**2,'Work MO vectors',
     1  l_w3,k_w3)) call errquit
     2  ('tddft_transfm: failed to allocate w2',0, MA_ERR)
c
c ------
c Type 1
c ------
c
      if (itype.eq.1) then
c
c -------------------------
c Retrieve a vector from GA
c -------------------------
c
        do i=1,ipol
#if 0
          if (nodezero.and.util_print('transform',print_debug)) 
     1      write(LuOut,9000) ivector,i
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w1),1)
          l=0
          do j=nfc(i)+1,nocc(i)
            do k=nocc(i)+1,nmo(i)-nfv(i)
              l=l+1
              call ga_get(g_vector(i),l,l,ivector,ivector,
     1          dbl_mb(k_w1+(j-1)*nbf_ao+k-1),1)
            enddo
          enddo
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w3),1)
          call ga_get(g_movecs(i),1,nbf_ao,1,nmo(i),dbl_mb(k_w3),
     1      nbf_ao)
          if (nodezero.and.util_print('transform',print_debug))
     1      call ma_print(dbl_mb(k_w3),nbf_ao,nbf_ao,'MO vectors')
          if (nodezero.and.util_print('transform',print_debug))
     1      call ma_print(dbl_mb(k_w1),nbf_ao,nbf_ao,'Trial vector')
c
c ------------------
c Transform MO to AO
c ------------------
c
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w2),1)
          call dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,
     1      1.0d0,dbl_mb(k_w3),nbf_ao,dbl_mb(k_w1),nbf_ao,
     2      0.0d0,dbl_mb(k_w2),nbf_ao)
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w1),1)
          call dgemm('N','T',nbf_ao,nbf_ao,nbf_ao,
     1      1.0d0,dbl_mb(k_w2),nbf_ao,dbl_mb(k_w3),nbf_ao,
     2      0.0d0,dbl_mb(k_w1),nbf_ao)
          call ga_put(g_density(i),1,nbf_ao,1,nbf_ao,
     1      dbl_mb(k_w1),nbf_ao)
#else
          if(.not.ga_duplicate(g_density(i), g_dd,' tddft '))
     A         call errquit('tddfttransfm: gadupl failed ',0,0)
          call ga_zero(g_dd)
          if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,
     1      'trial vectors1',-1,-1,g_vpatch))
     A         call errquit('tddfttransfm: gadupl failed ',0,0)
          call ga_zero(g_vpatch)
          l0=nocc(i)+1
          l1=nmo(i)-nfv(i)
          l1l0=l1-l0
          ilo=1
          ihi=ilo+l1l0
          l1l0=l1l0+1
          call ga_sync()
          do j=nfc(i)+1,nocc(i)
             call ga_mask_sync(.false.,.false.)
             call ga_copy_patch('N', g_vector(i), 
     A            ilo, ihi, ivector, ivector,
     B            g_vpatch, 
     C            l0, l1, j, j)
             ilo=ilo+l1l0
             ihi=ihi+l1l0
          enddo
          call ga_sync()
C     MO*V
          call ga_dgemm('N','N',nbf_ao,nocc(i),l1,1d0,
     G         g_movecs(i),g_vpatch,
     ,     0d0,g_dd)
C     (MO*V)*MO(transp)
          call ga_dgemm('N','T',nbf_ao,nbf_ao,nocc(i),1d0,
     G         g_dd,g_movecs(i),
     ,     0d0,g_density(i))
          if(.not.ga_destroy(g_dd))
     A         call errquit('tddfttransfm: gadestr failed ',0,0)
          if(.not.ga_destroy(g_vpatch))
     A         call errquit('tddfttransfm: gadestr failed ',0,0)
#endif
          if (util_print('transform',print_debug)) 
     1      call ga_print(g_density(i))
        enddo
c
c ------
c Type 2
c ------
c
      else if (itype.eq.2) then
c
c ------------------
c Transform AO to MO
c ------------------
c
        do i=1,ipol
#if 0
          if (nodezero.and.util_print('transform',print_debug)) 
     1      write(LuOut,9000) ivector,i
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w3),1)
          call ga_get(g_movecs(i),1,nbf_ao,1,nmo(i),
     1      dbl_mb(k_w3),nbf_ao)
          call ga_get(g_density(i),1,nbf_ao,1,nbf_ao,
     1      dbl_mb(k_w1),nbf_ao)
          if (nodezero.and.util_print('transform',print_debug))
     1      call ma_print(dbl_mb(k_w3),nbf_ao,nbf_ao,'MO vectors')
          if (nodezero.and.util_print('transform',print_debug))
     1      call ma_print(dbl_mb(k_w1),nbf_ao,nbf_ao,
     2      'Product density matrix')
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w2),1)
          call dgemm('T','N',nbf_ao,nbf_ao,nbf_ao,
     1      1.0d0,dbl_mb(k_w3),nbf_ao,dbl_mb(k_w1),nbf_ao,
     2      0.0d0,dbl_mb(k_w2),nbf_ao)
          call dfill(nbf_ao**2,0.0d0,dbl_mb(k_w1),1)
          call dgemm('N','N',nbf_ao,nbf_ao,nbf_ao,
     1      1.0d0,dbl_mb(k_w2),nbf_ao,dbl_mb(k_w3),nbf_ao,
     2      0.0d0,dbl_mb(k_w1),nbf_ao)
          if (nodezero.and.util_print('transform',print_debug))
     1      call ma_print(dbl_mb(k_w1),nbf_ao,nbf_ao,
     2      'Product vector')
c
c -------------------
c Copy an array to GA
c -------------------
c
          l=0
          do j=nfc(i)+1,nocc(i)
            do k=nocc(i)+1,nmo(i)-nfv(i)
              l=l+1
              call ga_put(g_vector(i),l,l,ivector,ivector,
     1          dbl_mb(k_w1+(j-1)*nbf_ao+k-1),1)
            enddo
          enddo
#else
          if(.not.ga_duplicate(g_density(i), g_dd,' tddft '))
     A         call errquit('tddfttransfm: gadupl failed ',0,0)
          if (.not.ga_create(mt_dbl,nbf_ao,nbf_ao,
     1      'trial vectors1',-1,-1,g_vpatch))
     A         call errquit('tddfttransfm: gadupl failed ',0,0)
          l0=nocc(i)+1
          l1=nmo(i)-nfv(i)
          call ga_zero(g_vpatch)
C     MO(tranps)*DM
          call ga_dgemm('T','N',l1,nbf_ao,nbf_ao,1d0,
     G         g_movecs(i),g_density(i),
     ,     0d0,g_dd)
C     (MO(transp)*DM)*MO
          call ga_dgemm('N','N',l1,nocc(i),nbf_ao,1d0,
     G         g_dd,g_movecs(i),
     ,     0d0,g_vpatch)
          if(.not.ga_destroy(g_dd))
     A         call errquit('tddfttransfm: gadestr failed ',0,0)
          l1l0=l1-l0
          ilo=1
          ihi=ilo+l1l0
          l1l0=l1l0+1
          call ga_sync()
          do j=nfc(i)+1,nocc(i)
             call ga_mask_sync(.false.,.false.)
             call ga_copy_patch('N', 
     B            g_vpatch,  l0, l1, j, j,
     G            g_vector(i), ilo, ihi, ivector, ivector)
             ilo=ilo+l1l0
             ihi=ihi+l1l0
          enddo
          call ga_sync()
          if(.not.ga_destroy(g_vpatch))
     A         call errquit('tddfttransfm: gadestr failed ',0,0)
#endif
        enddo
c
c ------------
c Unknown type
c ------------
c
      else
        call errquit('tddft_transfm: illegal itype value',0, INPUT_ERR)
      endif
c
c ------
c Return
c ------
c
      if (.not.ma_pop_stack(l_w3)) call errquit
     1  ('tddft_transfm: failed to deallocate w3',0, MA_ERR)
      if (.not.ma_pop_stack(l_w2)) call errquit
     1  ('tddft_transfm: failed to deallocate w2',0, MA_ERR)
      if (.not.ma_pop_stack(l_w1)) call errquit
     1  ('tddft_transfm: failed to deallocate w1',0, MA_ERR)
      return
 9000 format(' Vector: ',i5,' Spin: ',i5)
      end
