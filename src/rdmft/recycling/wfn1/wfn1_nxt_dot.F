C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the dot product of two state vectors
C> 
C> Compute the dot product of two state vectors or similar vectors. In
C> most ways this is the same as a regular dot-product except that for
C> the orbital part we need to compute
C> \f{eqnarray*}{
C>    dot &=& \sum_\sigma\mathrm{tr}((O^\sigma_{ip})^T(N^\sigma_{ai})^T,
C>            S_{ab}N^\sigma_{bj}O^\sigma_{jq})
C> \f}
C> where repeated indices are summed. For the Lagrangian multiplier
C> parts we use a simple dot-product as it is not clear at this point
C> what else would make sense.
C>
      double precision function wfn1_nxt_dot(nbf,ov,state1,state2)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "wfn1_nxt.fh"
c
      integer nbf !< [Input] The number of basis functions
c
      double precision ov(nbf,nbf) !< [Input] The overlap matrix
      double precision state1(8*nbf*nbf) !< [Input] State vector 1
      double precision state2(8*nbf*nbf) !< [Input] State vector 2
c
      integer l_t1, k_t1 !< Memory handle for temporary matrix
      integer l_t2, k_t2 !< Memory handle for temporary matrix
c
      double precision result !< Accumulator for the result
c
      double precision ddot
      external         ddot
c
      result = 0.0d0
c
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"T1",l_t1,k_t1))
     &  call errquit("wfn1_nxt_dot: could not allocate T1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_push_get(MT_DBL,nbf*nbf,"T2",l_t2,k_t2))
     &  call errquit("wfn1_nxt_dot: could not allocate T2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,state2(wfn1_nxt_anat()),nbf,
     &           state2(wfn1_nxt_aocc()),nbf,0.0d0,dbl_mb(k_t1),nbf)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     &           dbl_mb(k_t1),nbf,0.0d0,dbl_mb(k_t2),nbf)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,state1(wfn1_nxt_anat()),nbf,
     &           state1(wfn1_nxt_aocc()),nbf,0.0d0,dbl_mb(k_t1),nbf)
      result = result + ddot(nbf*nbf,dbl_mb(k_t1),1,dbl_mb(k_t2),1)
c
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,state2(wfn1_nxt_bnat()),nbf,
     &           state2(wfn1_nxt_bocc()),nbf,0.0d0,dbl_mb(k_t1),nbf)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     &           dbl_mb(k_t1),nbf,0.0d0,dbl_mb(k_t2),nbf)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,state1(wfn1_nxt_bnat()),nbf,
     &           state1(wfn1_nxt_bocc()),nbf,0.0d0,dbl_mb(k_t1),nbf)
      result = result + ddot(nbf*nbf,dbl_mb(k_t1),1,dbl_mb(k_t2),1)
c
      result = result + ddot(nbf*nbf,state1(wfn1_nxt_anatl()),1,
     &                               state2(wfn1_nxt_anatl()),1)
      result = result + ddot(nbf*nbf,state1(wfn1_nxt_bnatl()),1,
     &                               state2(wfn1_nxt_bnatl()),1)
c
      result = result + ddot(nbf*nbf,state1(wfn1_nxt_aoccl()),1,
     &                               state2(wfn1_nxt_aoccl()),1)
      result = result + ddot(nbf*nbf,state1(wfn1_nxt_boccl()),1,
     &                               state2(wfn1_nxt_boccl()),1)
c
      if (.not.ma_pop_stack(l_t2))
     &  call errquit("wfn1_nxt_dot: could not deallocate T2",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
      if (.not.ma_pop_stack(l_t1))
     &  call errquit("wfn1_nxt_dot: could not deallocate T1",
     &  ma_sizeof(MT_DBL,nbf*nbf,MT_BYTE),MA_ERR)
c
      wfn1_nxt_dot = result
c
      return
      end
C>
C> @}
