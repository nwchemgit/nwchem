C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Print the overlap between components of the next step
C> wavefunction
C>
C> Print the overlap among the alpha and beta occupation functions
C> and the alpha and beta natural orbitals.
C>
      subroutine wfn1_nxt_print_overlap(nbf,state,ov)
      implicit none
C
#include "errquit.fh"
#include "mafdecls.fh"
#include "wfn1_nxt.fh"
C
      integer nbf !< [Input] The number of basis functions
      double precision state(8*nbf*nbf) !< [Input] The wavefunction
      double precision ov(nbf,nbf) !< [Input] The overlap matrix
c
      integer l_m1, k_m1
      integer l_m2, k_m2
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"scratch m1",l_m1,k_m1))
     &  call errquit("wfn1_nxt_print_overlap: could not alloc m1",
     &               nbf*nbf,MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"scratch m2",l_m2,k_m2))
     &  call errquit("wfn1_nxt_print_overlap: could not alloc m2",
     &               nbf*nbf,MA_ERR)
c
      write(*,*)
      write(*,*)'alpha occupation functions'
      write(*,*)'--------------------------'
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,
     &           state(wfn1_nxt_aocc()),nbf,state(wfn1_nxt_aocc()),nbf,
     &           0.0d0,dbl_mb(k_m2),nbf)
      call hess_hssout(dbl_mb(k_m2),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta occupation functions'
      write(*,*)'-------------------------'
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,
     &           state(wfn1_nxt_bocc()),nbf,state(wfn1_nxt_bocc()),nbf,
     &           0.0d0,dbl_mb(k_m2),nbf)
      call hess_hssout(dbl_mb(k_m2),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'alpha natural orbitals'
      write(*,*)'----------------------'
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &           ov,nbf,state(wfn1_nxt_anat()),nbf,
     &           0.0d0,dbl_mb(k_m1),nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,
     &           state(wfn1_nxt_anat()),nbf,dbl_mb(k_m1),nbf,
     &           0.0d0,dbl_mb(k_m2),nbf)
      call hess_hssout(dbl_mb(k_m2),nbf,nbf,nbf)
      write(*,*)
      write(*,*)'beta natural orbitals'
      write(*,*)'---------------------'
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,
     &           ov,nbf,state(wfn1_nxt_bnat()),nbf,
     &           0.0d0,dbl_mb(k_m1),nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,
     &           state(wfn1_nxt_bnat()),nbf,dbl_mb(k_m1),nbf,
     &           0.0d0,dbl_mb(k_m2),nbf)
      call hess_hssout(dbl_mb(k_m2),nbf,nbf,nbf)
c
      if (.not.ma_pop_stack(l_m2))
     &  call errquit("wfn1_nxt_print_overlap: could not dealloc m2",
     &               nbf*nbf,MA_ERR)
      if (.not.ma_pop_stack(l_m1))
     &  call errquit("wfn1_nxt_print_overlap: could not dealloc m1",
     &               nbf*nbf,MA_ERR)
C
      end
C> @}
