C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the first Lagrangian (L1) energy contribution and
C> optionally the gradient within the wfn1 formalism
C> 
C> This routine adds the first Lagrangian (L1) energy and optionally the
C> corresponding gradient terms within the wfn1 formalism to its
C> arguments. The current formulation assumes that all quantities are
C> kept in memory and the implementation supports only serial execution.
C>
C> Within the wfn1 formalism the L1 energy is expressed in terms
C> the natural orbitals as:
C> \f{eqnarray*}{
C>  E^{L1} &=& \sum_{\sigma}\sum_{ij}
C>  \lambda^{1\sigma}_{ij}\left(\delta_{ij}-\sum_{ab}
C>  C^{D\sigma}_{ai}S_{ab}C^{D\sigma}_{bj}\right)
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> where \f$i\f$ and \f$j\f$ run over occupied natural orbitals only.
C> The argument list contains special tables that indicate whether
C> a given natural orbital is occupied or not.
C>
C> The gradient requires the derivatives with respect to
C> the Lagrangian multipliers as well as the natural orbitals. I.e.
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L1}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& 2\sum_{\sigma}\sum_{ij}\delta_{\sigma\sigma"}\delta_{im}
C>    \lambda^{1\sigma}_{ij}\left(\delta_{ij}-\sum_{ab}
C>    \delta_{ea}S_{ab}C^{D\sigma}_{bj}\right)
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>    &=& 2\sum_{j}\sum_{b}
C>    -\lambda^{1\sigma"}_{mj}
C>    S_{eb}C^{D\sigma"}_{bj}
C>    \;\;\;\;\;\;\;\;\;\;(3) \\\\
C> \f}
C> and
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L1}}{\mathrm{d}\lambda^{1\sigma"}_{mn}}
C>    &=& \sum_{\sigma}\sum_{ij}\delta_{\sigma\sigma"}
C>    \delta_{im}\delta_{jn}
C>    \left(\delta_{ij}-\sum_{ab}
C>    C^{D\sigma}_{ai}S_{ab}C^{D\sigma}_{bj}\right)
C>    \;\;\;\;\;\;\;\;\;\;(4) \\\\
C>    &=& \left(\delta_{mn}-\sum_{ab}
C>    C^{D\sigma"}_{am}S_{ab}C^{D\sigma"}_{bn}\right)
C>    \;\;\;\;\;\;\;\;\;\;(5) \\\\
C> \f}
C>
      subroutine wfn1_l1(do_grad,nbf,nst,l_occa,l_occb,r_noa,r_nob,
     +                   r_l1a,r_l1b,ov,
     +                   En,r_dnoa,r_dnob,r_dl1a,r_dl1b,
     +                   dc1,dc2,dc3)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      logical do_grad  !< [Input] .True. if the energy gradient should
                       !< be calculated, otherwise only the energy is
                       !< evaluated.
      integer nbf      !< [Input] The number of basis functions
      integer nst      !< [Input] The number of states
c
      logical l_occa(nbf,nst) !< [Input] .True. if an 
                              !< \f$\alpha\f$-orbital is occupied
      logical l_occb(nbf,nst) !< [Input] .True. if a 
                              !< \f$\beta\f$-orbital is occupied
c
      double precision r_noa(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron natural orbitals
      double precision r_nob(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron natural orbitals
      double precision r_l1a(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron first Lagrangian
                                          !< multipliers
      double precision r_l1b(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron first Lagrangian
                                          !< multipliers
c
      double precision ov(nbf,nbf) !< [Input] The actual overlap
                                   !< integrals
c
      double precision En(nst) !< [Output] The 1-electron
                               !< energy for each state
c
      double precision r_dnoa(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< natural orbitals.
      double precision r_dnob(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< natural orbitals.
      double precision r_dl1a(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< first Lagrangian
                                           !< multipliers.
      double precision r_dl1b(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< first Lagrangian
                                           !< multipliers.
      double precision dc1(nbf,nbf)    !< [Scratch] Partially
                                       !< transformed matrix
      double precision dc2(nbf,nbf)    !< [Scratch] Partially
                                       !< transformed matrix
      double precision dc3(nbf,nbf)    !< [Scratch] Copy of natural
                                       !< orbital matrix
c
c     Local variables
c
      integer ia, ib ! Counters over basis functions
      integer ii, jj ! Counter over natural orbitals
      integer ir     ! Counter over occupation functions
      integer ist    ! Counter over states
      integer ione   ! The number 1
      parameter (ione = 1)
c
c     Functions
c
      double precision ddot
      external         ddot
c
c     Code
c
      do ist = 1, nst ! loop over all states
c
c       Compute the overlap matrices for this state
c
        call dcopy(nbf*nbf,r_noa(1,1,ist),1,dc3,1)
        do ii = 1, nbf
          if (.not.l_occa(ii,ist)) then
            call dfill(nbf,0.0d0,dc3(1,ii),1)
          endif
        enddo
        call dgemm('T','N',nbf,nbf,nbf,-1.0d0,dc3,nbf,
     +             ov,nbf,0.0d0,dc1,nbf)
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,dc3,nbf,
     +             0.0d0,dc2,nbf)
cDEBUG
c      write(*,*)'nat orb lagran'
c      call hess_hssout(dc2,nbf,nbf,nbf)
cDEBUG
        do ii = 1, nbf
          if (l_occa(ii,ist)) then
            dc2(ii,ii) = dc2(ii,ii) + 1.0d0
          endif
        enddo
c
c       Contract the 1-electron density matrix with the
c       1-electron integrals and add the energy contribution
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l1a(1,1,ist),ione,dc2,ione)
c
        if (do_grad) then
          call dcopy(nbf*nbf,dc2,ione,r_dl1a(1,1,ist),ione)
c
c         Compute the gradient wrt. the natural orbital
c         coefficients, see Eq.(3).
c
c         = A-A: Compute F^A_am = sum_b S_ab C^DA_bm
c
          call dgemm('N','N',nbf,nbf,nbf,-2.0d0,ov,nbf,
     +               dc3,nbf,0.0d0,dc2,nbf)
c
c         = A-A: Multiply with the Lagrangian multipliers and
c           add to the gradient (assuming the Lagrangian matrix is
c           symmetric)
c
          call dcopy(nbf*nbf,r_l1a(1,1,ist),ione,dc3,ione)
          do jj = 1, nbf
            if (l_occa(jj,ist)) then
              do ii = 1, nbf
                if (.not.l_occa(ii,ist)) then
                  dc3(ii,jj) = 0.0d0
                endif
              enddo
            else
              call dfill(nbf,0.0d0,dc3(1,jj),ione)
            endif
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc2,nbf,
     +               dc3,nbf,1.0d0,r_dnoa(1,1,ist),nbf)
        endif
c
c       Compute the overlap matrices for this state
c
        call dcopy(nbf*nbf,r_nob(1,1,ist),1,dc3,1)
        do ii = 1, nbf
          if (.not.l_occb(ii,ist)) then
            call dfill(nbf,0.0d0,dc3(1,ii),1)
          endif
        enddo
        call dgemm('T','N',nbf,nbf,nbf,-1.0d0,dc3,nbf,
     +             ov,nbf,0.0d0,dc1,nbf)
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,dc3,nbf,
     +             0.0d0,dc2,nbf)
        do ii = 1, nbf
          if (l_occb(ii,ist)) then
            dc2(ii,ii) = dc2(ii,ii) + 1.0d0
          endif
        enddo
c
c       Contract the 1-electron density matrix with the
c       1-electron integrals and add the energy contribution
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l1b(1,1,ist),ione,dc2,ione)
c
        if (do_grad) then
c
          call dcopy(nbf*nbf,dc2,ione,r_dl1b(1,1,ist),ione)
c
c         Compute the gradient wrt. the natural orbital
c         coefficients, see Eq.(3).
c
c         = B-B: Compute F^B_am = sum_b S_ab C^DB_bm
c
          call dgemm('N','N',nbf,nbf,nbf,-2.0d0,ov,nbf,
     +               dc3,nbf,0.0d0,dc2,nbf)
c
c         = A-A: Multiply with the Lagrangian multipliers and
c           add to the gradient (assuming the Lagrangian matrix is
c           symmetric)
c
          call dcopy(nbf*nbf,r_l1b(1,1,ist),ione,dc3,ione)
          do jj = 1, nbf
            if (l_occb(jj,ist)) then
              do ii = 1, nbf
                if (.not.l_occb(ii,ist)) then
                  dc3(ii,jj) = 0.0d0
                endif
              enddo
            else
              call dfill(nbf,0.0d0,dc3(1,jj),ione)
            endif
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc2,nbf,
     +               dc3,nbf,1.0d0,r_dnob(1,1,ist),nbf)
        endif
c
      enddo
c
      end
C> @}
