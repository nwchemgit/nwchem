C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Prepare the gradient to update the state
C>
C> In the linesearch algorithm updating the penalty coefficients
C> often drives the energy up so that the optimization returns a
C> zero step length. Hence we will try what a linesearch with 
C> constant penalty coefficients does. To do this we must zero all
C> gradient components for the penalty coefficients.
C>
      subroutine wfn1_nxt_prep_grad(nbf,gradin,gradout)
      implicit none
C
#include "wfn1_nxt.fh"
C
      integer nbf !< [Input] The number of basis functions
C
      double precision gradin(8*nbf*nbf) !< [Input] The gradient
      double precision gradout(8*nbf*nbf) !< [Input] The prepped gradient
C
      call dcopy(wfn1_nxt_size1(),gradin,1,gradout,1)
      call dfill(nbf*nbf,0.0d0,gradout(wfn1_nxt_aoccl()),1)
      call dfill(nbf*nbf,0.0d0,gradout(wfn1_nxt_boccl()),1)
      call dfill(nbf*nbf,0.0d0,gradout(wfn1_nxt_anatl()),1)
      call dfill(nbf*nbf,0.0d0,gradout(wfn1_nxt_bnatl()),1)
C
      end
C> @}
      
