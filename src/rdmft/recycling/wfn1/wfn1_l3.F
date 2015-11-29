C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the third Lagrangian (L3) energy contribution and
C> optionally the gradient within the wfn1 formalism
C> 
C> This routine adds the third Lagrangian (L3) energy and optionally the
C> corresponding gradient terms within the wfn1 formalism to its
C> arguments. The current formulation assumes that all quantities are
C> kept in memory and the implementation supports only serial execution.
C>
C> Within the wfn1 formalism the L3 energy is expressed in terms of
C> the occupation functions as:
C> \f{eqnarray*}{
C>  E^{L3} &=& \sum_{\sigma}\sum_{pq}
C>  \lambda^{3\sigma}_{pq}\left(\delta_{pq}-\sum_{i}
C>  C^{O\sigma}_{ip}C^{O\sigma}_{iq}\right)
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> where \f$p\f$ and \f$q\f$ are filled occupation functions.
C>
C> The gradient requires the derivatives with respect to
C> the Lagrangian multipliers as well as the occupation functions. I.e.
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L3}}{\mathrm{d}C^{O\sigma"}_{mu}}
C>    &=& -2\sum_{\sigma}\sum_{pq}\sum_i
C>    \delta_{\sigma\sigma"}\delta_{im}\delta_{qu}
C>    \lambda^{3\sigma}_{pq}
C>    C^{O\sigma}_{ip}
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>    &=& -2\sum_{p}
C>    C^{O\sigma"}_{mp}\lambda^{3\sigma"}_{pu}
C>    \;\;\;\;\;\;\;\;\;\;(3) \\\\
C> \f}
C> and
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{L3}}{\mathrm{d}\lambda^{3\sigma"}_{uv}}
C>    &=& \sum_{\sigma}\sum_{pq}
C>    \delta_{\sigma\sigma"}\delta_{pu}\delta_{qv}
C>    \left(\delta_{pq}-\sum_{i}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{iq}\right)
C>    \;\;\;\;\;\;\;\;\;\;(4) \\\\
C>    &=& \left(\delta_{uv}-\sum_{i}
C>    C^{O\sigma"}_{iu}C^{O\sigma"}_{iv}\right)
C>    \;\;\;\;\;\;\;\;\;\;(5) \\\\
C> \f}
C>
      subroutine wfn1_l3(do_grad,nbf,nst,nea,neb,r_ofa,r_ofb,
     +                   r_l3a,r_l3b,
     +                   En,r_dofa,r_dofb,r_dl3a,r_dl3b,dcd)
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
      integer nea(nst) !< [Input] The number of \f$\alpha\f$-electrons
                       !< for each state
      integer neb(nst) !< [Input] The number of \f$\beta\f$-electrons
                       !< for each state
c
      double precision r_ofa(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron occupation
                                          !< functions
      double precision r_ofb(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron occupation
                                          !< functions
      double precision r_l3a(nbf,nbf,nst) !< [Input] The \f$\alpha\f$
                                          !< electron third Lagrangian
                                          !< multipliers
      double precision r_l3b(nbf,nbf,nst) !< [Input] The \f$\beta\f$
                                          !< electron third Lagrangian
                                          !< multipliers
c
      double precision En(nst) !< [Output] The 1-electron
                               !< energy for each state
c
      double precision r_dofa(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< occupation functions.
      double precision r_dofb(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< occupation functions.
      double precision r_dl3a(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\alpha\f$-electron
                                           !< third Lagrangian
                                           !< multipliers.
      double precision r_dl3b(nbf,nbf,nst) !< [Output] The gradient of
                                           !< the energy with respect to
                                           !< the \f$\beta\f$-electron
                                           !< third Lagrangian
                                           !< multipliers.
      double precision dcd(nbf,nbf) !< [Scratch] Temporary matrix
c
c     Local variables
c
      integer ip, iq ! Counters over occupation functions
      integer ii     ! Counter over natural orbitals
      integer ist    ! Counter over states
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
c       Do the matrix-matrix multiplication of the filled occupation
c       functions first
c
        call dfill(nbf*nbf,0.0d0,dcd,1)
        call dgemm('T','N',nea(ist),nea(ist),nbf,-1.0d0,
     +             r_ofa(1,1,ist),nbf,r_ofa(1,1,ist),nbf,0.0d0,dcd,nbf)
c
c       Add 1 to the diagonal
c
        do ip = 1, nea(ist)
          dcd(ip,ip) = dcd(ip,ip) + 1.0d0
        enddo
c
c       Contract the scratch matrix with the third Lagrangian
c       and add the energy contribution
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l3a(1,1,ist),1,dcd,1)
c
        if (do_grad) then
          call dcopy(nbf*nbf,dcd,1,r_dl3a(1,1,ist),1)
        endif
c
c       Do the matrix-matrix multiplication of the occupation functions
c       first
c
        call dfill(nbf*nbf,0.0d0,dcd,1)
        call dgemm('T','N',neb(ist),neb(ist),nbf,-1.0d0,
     +             r_ofb(1,1,ist),nbf,r_ofb(1,1,ist),nbf,0.0d0,dcd,nbf)
c
c       Add 1 to the diagonal
c
        do ip = 1, neb(ist)
          dcd(ip,ip) = dcd(ip,ip) + 1.0d0
        enddo
c
c       Contract the scratch matrix with the third Lagrangian
c       and add the energy contribution
c
        En(ist) = En(ist) + ddot(nbf*nbf,r_l3b(1,1,ist),1,dcd,1)
c
        if (do_grad) then
c
          call dcopy(nbf*nbf,dcd,1,r_dl3b(1,1,ist),1)
c
c         First compute the gradient wrt. the occupation function
c         coefficients, see Eq.(3).
c
          call dgemm('N','N',nbf,nea(ist),nea(ist),-2.0d0,
     +               r_ofa(1,1,ist),nbf,r_l3a(1,1,ist),nbf,1.0d0,
     +               r_dofa(1,1,ist),nbf)
c
          call dgemm('N','N',nbf,neb(ist),neb(ist),-2.0d0,
     +               r_ofb(1,1,ist),nbf,r_l3b(1,1,ist),nbf,1.0d0,
     +               r_dofb(1,1,ist),nbf)
c
        endif
c
      enddo
c
      end
C> @}
