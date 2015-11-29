C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the second Lagrangian (L2) energy contribution and
C> optionally the gradient within the wfn1 formalism
C> 
C> This routine adds the second Lagrangian (L2) energy and optionally
C> the corresponding gradient terms within the wfn1 formalism to its
C> arguments. The current formulation assumes that all quantities are
C> kept in memory and the implementation supports only serial execution.
C>
C> Within the wfn1 formalism the L2 energy is expressed in terms of
C> the occupation functions as:
C> \f{eqnarray*}{
C>  E^{L2} &=& \sum_{\sigma}\sum_{ij}
C>  \lambda^{2\sigma}_{ij}\left(1-\delta_{ij}\right)\sum_{abcd}
C>  C^{D\sigma}_{ai}S_{ab}D^{\sigma}_{bc}S_{cd}C^{D\sigma}_{dj},
C>  \;\;\;\;\;\;\;\;\;\;(1) \\\\
C>  D^{\sigma}_{bc}
C>  &=& \sum_{k}\sum_{p\in\{occ\}}
C>  C^{D\sigma}_{bk}C^{D\sigma}_{ck}C^{O\sigma}_{kp}C^{O\sigma}_{kp}
C> \f}
C> The gradient requires the derivatives with respect to
C> the Lagrangian multipliers as well as the natural orbitals and the
C> occupation functions. I.e.
C> \f{eqnarray*}{
C>  \frac{\mathrm{d}E^{L2}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& 2\sum_{\sigma}\sum_{ij}
C>    \delta_{\sigma\sigma"}\delta_{ae}\delta_{im}
C>    \lambda^{2\sigma}_{ij}\left(1-\delta_{ij}\right)\sum_{abcd}
C>    S_{ab}D^{\sigma}_{bc}S_{cd}C^{D\sigma}_{dj} \\\\
C>    &+& \sum_{\sigma}\sum_{ij}
C>    \lambda^{2\sigma}_{ij}\left(1-\delta_{ij}\right)\sum_{abcd}
C>    C^{D\sigma}_{ai}S_{ab}\frac{D^{\sigma}_{bc}}{C^{D\sigma"}_{em}}
C>    S_{cd}C^{D\sigma}_{dj},
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>  \frac{D^{\sigma}_{bc}}{C^{D\sigma"}_{em}}
C>    &=& 2\sum_{k}\sum_{p\in\{occ\}}
C>    \delta_{\sigma\sigma"}\delta_{be}\delta_{km}
C>    C^{D\sigma}_{ck}C^{O\sigma}_{kp}C^{O\sigma}_{kp} \\\\
C>  \frac{\mathrm{d}E^{L2}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& 2\sum_{j}
C>    \lambda^{2\sigma"}_{mj}\left(1-\delta_{mj}\right)\sum_{bcd}
C>    S_{eb}D^{\sigma"}_{bc}S_{cd}C^{D\sigma"}_{dj} \\\\
C>    &+& 2\sum_{ij}
C>    \lambda^{2\sigma"}_{ij}\left(1-\delta_{ij}\right)\sum_{acd}
C>    C^{D\sigma"}_{ai}S_{ae}
C>    \left(\sum_{p\in\{occ\}}
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}\right)
C>    S_{cd}C^{D\sigma"}_{dj}
C>    \;\;\;\;\;\;\;\;\;\;(3) \\\\
C>  \frac{\mathrm{d}E^{L2}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& 2\sum_{jbcd}S_{eb}D^{\sigma"}_{bc}S_{cd}C^{D\sigma"}_{dj}
C>    \lambda^{2\sigma"}_{jm}\left(1-\delta_{jm}\right) \\\\
C>    &+& 2\sum_{ijacd}S_{ae}C^{D\sigma"}_{ai}
C>    \lambda^{2\sigma"}_{ij}\left(1-\delta_{ij}\right)
C>    C^{D\sigma"}_{dj}S_{dc}
C>    \left(\sum_{p\in\{occ\}}
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}\right)
C>    \;\;\;\;\;\;\;\;\;\;(4) \\\\
C> \f}
C> Considering the derivative wrt. the occupation functions we have
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L2}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>      &=& \sum_{\sigma}\sum_{ij}
C>      \lambda^{2\sigma}_{ij}\left(1-\delta_{ij}\right)\sum_{abcd}
C>      C^{D\sigma}_{ai}S_{ab}
C>      \frac{\mathrm{d}D^{\sigma}_{bc}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>      S_{cd}C^{D\sigma}_{dj},
C>      \;\;\;\;\;\;\;\;\;\;(5) \\\\
C>    \frac{\mathrm{d}D^{\sigma}_{bc}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>      &=& 2\sum_k\sum_{p\in\{occ\}}
C>      \delta_{\sigma\sigma"}\delta_{km}\delta_{pu}
C>      C^{D\sigma}_{bk}C^{D\sigma}_{ck}C^{O\sigma}_{kp} \\\\
C>    \frac{\mathrm{d}D^{\sigma"}_{bc}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>      &=& 2C^{D\sigma"}_{bm}C^{D\sigma"}_{cm}C^{O\sigma"}_{mu}
C>    \;\;\;\;\;\;\;\;\;\;(6) \\\\
C>    \frac{\mathrm{d}E^{L2}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>      &=& 2\sum_{ij}\sum_{abcd}
C>      C^{D\sigma"}_{bm}S_{ab}C^{D\sigma"}_{ai}
C>      \{\lambda^{2\sigma"}_{ij}(1-\delta_{ij})\}
C>      C^{D\sigma"}_{dj}S_{cd}C^{D\sigma"}_{cm}C^{O\sigma"}_{mu}
C>      \;\;\;\;\;\;\;\;\;\;(7) \\\\
C> \f}
C> Finally, the derivatives wrt. the second Lagrangians are
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L2}}{\mathrm{d}\lambda^{2\sigma"}_{kl}}
C>    &=& \left(1-\delta_{kl}\right)\sum_{abcd}
C>    C^{D\sigma"}_{ak}S_{ab}D^{\sigma"}_{bc}S_{cd}C^{D\sigma"}_{dl},
C>    \;\;\;\;\;\;\;\;\;\;(8) \\\\
C> \f}
C>
      subroutine wfn1_l2(do_grad,nbf,nst,nea,neb,r_noa,r_nob,
     +                   r_ofa,r_ofb,r_l2a,r_l2b,ov,
     +                   En,r_dnoa,r_dnob,r_dofa,r_dofb,r_dl2a,r_dl2b,
     +                   dc1,dc2,dc3,docc)
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
      integer nea(nst) !< [Input] The number of alpha electrons
      integer neb(nst) !< [Input] The number of beta electrons
c
      double precision r_noa(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron natural orbitals
      double precision r_nob(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron natural orbitals
      double precision r_ofa(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron occupation
                                          !< functions
      double precision r_ofb(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron occupation
                                          !< functions
      double precision r_l2a(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron second Lagrangian
                                          !< multipliers
      double precision r_l2b(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron second Lagrangian
                                          !< multipliers
      double precision ov(nbf,nbf)        !< [Input] The overlap
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
      double precision r_dofa(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< occupation functions.
      double precision r_dofb(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< occupation functions.
      double precision r_dl2a(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< second Lagrangian
                                           !< multipliers.
      double precision r_dl2b(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< second Lagrangian
                                           !< multipliers.
      double precision dc1(nbf,nbf) !< [Scratch] Temporary matrix
      double precision dc2(nbf,nbf) !< [Scratch] Temporary matrix
      double precision dc3(nbf,nbf) !< [Scratch] Temporary matrix
      double precision docc(nbf)    !< [Scratch] Temporary vector
c
c     Local variables
c
      integer ip             ! Counter over occupation functions
      integer ii, ij         ! Counters over natural orbitals
      integer ia, ib, ic, id ! Counters over natural orbitals
      integer ist            ! Counter over states
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
c       == Alpha electrons first ==
c
c       Calculate occupation numbers
c
        call dfill(nbf,0.0d0,docc,1)
        do ip = 1, nea(ist)
          do ii = 1, nbf
            docc(ii) = docc(ii) + r_ofa(ii,ip,ist)*r_ofa(ii,ip,ist)
          enddo
        enddo
c
c       Calculate A_{cj} = \sum_d S_{cd}C^{Da}_{dj} first
c
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +             r_noa(1,1,ist),nbf,0.0d0,dc1,nbf)
c
c       Calculate B_{kj} = \sum_c (C^{Da}_{ck})^T A_{cj} next
c
        call dgemm('T','N',nbf,nbf,nbf,1.0d0,
     +             r_noa(1,1,ist),nbf,dc1,nbf,0.0d0,dc2,nbf)
c
c       Calculate C_{jk} = P_k * B_{kj}
c
        do ii = 1, nbf
          do ij = 1, nbf
            dc1(ij,ii) = docc(ii)*dc2(ii,ij)
          enddo
        enddo
c
c       Calculate the Lagrangian factor
c
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,
     +             dc2,nbf,0.0d0,dc3,nbf)
c
c       Apply the factor (1-\delta_{ij})
c
        do ii = 1, nbf
          dc3(ii,ii) = 0.0d0
        enddo
c
c       Calculate energy contribution 
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l2a(1,1,ist),1,dc3,1)
c
        if (do_grad) then
          call daxpy(nbf*nbf,1.0d0,dc3,1,r_dl2a(1,1,ist),1)
        endif
c
        if (do_grad) then
c
c         evaluate Eq.(4), first term
c
          call dfill(nbf*nbf,0.0d0,dc3,1)
          do ii = 1, nbf
            do ib = 1, nbf
              do ia = 1, nbf
                dc3(ia,ib) = dc3(ia,ib)
     +                     + r_noa(ia,ii,ist)*r_noa(ib,ii,ist)*docc(ii)
              enddo
            enddo
          enddo
          do ii = 1, nbf
            r_l2a(ii,ii,ist) = 0.0d0
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,
     +               r_l2a(1,1,ist),nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc3,nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,2.0d0,ov,nbf,
     +               dc1,nbf,1.0d0,r_dnoa(1,1,ist),nbf)
c
c         evaluate Eq.(4), second term
c
          do ii = 1, nbf
            do ia = 1, nbf
              dc3(ia,ii) = r_noa(ia,ii,ist)*docc(ii)
            enddo
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               dc3,nbf,0.0d0,dc1,nbf)
          call dgemm('T','N',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_l2a(1,1,ist),nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,2.0d0,ov,nbf,
     +               dc2,nbf,1.0d0,r_dnoa(1,1,ist),nbf)
c
        endif
c
c       == Beta electrons next ==
c
c       Calculate occupation numbers
c
        call dfill(nbf,0.0d0,docc,1)
        do ip = 1, neb(ist)
          do ii = 1, nbf
            docc(ii) = docc(ii) + r_ofb(ii,ip,ist)*r_ofb(ii,ip,ist)
          enddo
        enddo
c
c       Calculate A_{cj} = \sum_d S_{cd}C^{Da}_{dj} first
c
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +             r_nob(1,1,ist),nbf,0.0d0,dc1,nbf)
c
c       Calculate B_{kj} = \sum_c (C^{Da}_{ck})^T A_{cj} next
c
        call dgemm('T','N',nbf,nbf,nbf,1.0d0,
     +             r_nob(1,1,ist),nbf,dc1,nbf,0.0d0,dc2,nbf)
c
c       Calculate C_{jk} = P_k * B_{kj}
c
        do ii = 1, nbf
          do ij = 1, nbf
            dc1(ij,ii) = docc(ii)*dc2(ii,ij)
          enddo
        enddo
c
c       Calculate the Lagrangian factor
c
        call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,
     +             dc2,nbf,0.0d0,dc3,nbf)
c
c       Apply the factor (1-\delta_{ij})
c
        do ii = 1, nbf
          dc3(ii,ii) = 0.0d0
        enddo
c
c       Calculate energy contribution 
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l2b(1,1,ist),1,dc3,1)
c
        if (do_grad) then
          call daxpy(nbf*nbf,1.0d0,dc3,1,r_dl2b(1,1,ist),1)
        endif
c
        if (do_grad) then
c
c         evaluate Eq.(4), first term
c
          call dfill(nbf*nbf,0.0d0,dc3,1)
          do ii = 1, nbf
            do ib = 1, nbf
              do ia = 1, nbf
                dc3(ia,ib) = dc3(ia,ib)
     +                     + r_nob(ia,ii,ist)*r_nob(ib,ii,ist)*docc(ii)
              enddo
            enddo
          enddo
          do ii = 1, nbf
            r_l2b(ii,ii,ist) = 0.0d0
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,
     +               r_l2b(1,1,ist),nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc3,nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,2.0d0,ov,nbf,
     +               dc1,nbf,1.0d0,r_dnob(1,1,ist),nbf)
c
c         evaluate Eq.(4), second term
c
          do ii = 1, nbf
            do ia = 1, nbf
              dc3(ia,ii) = r_nob(ia,ii,ist)*docc(ii)
            enddo
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               dc3,nbf,0.0d0,dc1,nbf)
          call dgemm('T','N',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_l2b(1,1,ist),nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,2.0d0,ov,nbf,
     +               dc2,nbf,1.0d0,r_dnob(1,1,ist),nbf)
c
        endif
c
c       Calculate derivatives wrt the occupation functions
c
        if (do_grad) then
c
c         Alpha-electrons
c
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               r_noa(1,1,ist),nbf,0.0d0,dc2,nbf)
          call dgemm('T','N',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          do ii = 1, nbf
            r_l2a(ii,ii,ist) = 0.0d0
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_l2a(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,
     +               dc2,nbf,0.0d0,dc3,nbf)
          do ip = 1, nea(ist)
            do ii = 1, nbf
              r_dofa(ii,ip,ist) = r_dofa(ii,ip,ist)
     +                          + dc3(ii,ii)*r_ofa(ii,ip,ist)
            enddo
          enddo
c
c         Beta-electrons
c
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,
     +               r_nob(1,1,ist),nbf,0.0d0,dc2,nbf)
          call dgemm('T','N',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,
     +               dc2,nbf,0.0d0,dc1,nbf)
          do ii = 1, nbf
            r_l2b(ii,ii,ist) = 0.0d0
          enddo
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,r_l2b(1,1,ist),nbf,
     +               dc1,nbf,0.0d0,dc2,nbf)
          call dgemm('N','N',nbf,nbf,nbf,1.0d0,dc1,nbf,
     +               dc2,nbf,0.0d0,dc3,nbf)
          do ip = 1, neb(ist)
            do ii = 1, nbf
              r_dofb(ii,ip,ist) = r_dofb(ii,ip,ist)
     +                          + dc3(ii,ii)*r_ofb(ii,ip,ist)
            enddo
          enddo
        endif
c
      enddo
c
      end
C> @}
