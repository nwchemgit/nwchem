C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the 1-electron energy and optionally the gradient
C> within the wfn1 formalism
C> 
C> This routine adds the 1-electron energy and optionally the
C> corresponding gradient terms within the wfn1 formalism to its
C> arguments. The current formulation assumes that all quantities are
C> kept in memory and the implementation supports only serial execution.
C>
C> Within the wfn1 formalism the 1-electron energy is expressed in terms
C> of the occupation functions and the natural orbitals as:
C> \f{eqnarray*}{
C>  E^{1} &=& \sum_{\sigma}\sum_{i}
C>  \sum_{p\in\{occ\}}\sum_{ab}
C>  C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>  h_{ab}
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> Obviously this energy expression can be summarized in terms of the
C> total 1-electron density matrix as
C> \f{eqnarray*}{
C>    E^{1} &=& \sum_{\sigma}\sum_{ab}
C>    D^{T\sigma}_{ab}h_{ab},
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>    D^{T\sigma}_{ab} &=& \sum_i \sum_{p}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C> \f}
C> Although the energy can be summarized this way the gradient of the
C> energy cannot. The gradient requires the derivatives with respect to
C> the occupation functions as well as the natural orbitals. I.e.
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{1}}{\mathrm{d}C^{O\sigma"}_{mt}}
C>    &=& \sum_{\sigma}\sum_{ab}
C>    \sum_{i}\sum_{p} \delta_{\sigma\sigma"}\delta_{im}\delta_{pt}
C>    2C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}
C>    h_{ab}
C>    \;\;\;\;\;\;\;\;\;\;(3)\\\\
C>    &=& \sum_{ab}
C>    \sum_{j}
C>    2C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    h_{ab}
C>    \;\;\;\;\;\;\;\;\;\;(4)\\\\
C> \f}
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{1}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& \sum_{\sigma}\sum_{ab}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{ae}\delta_{im}
C>    C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    h_{ab} \\\\
C>    &&+ \sum_{\sigma}\sum_{ab}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{be}\delta_{im}
C>    C^{D\sigma}_{ai}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    h_{ab}
C>    \;\;\;\;\;\;\;\;\;\;(5)\\\\
C>    &=& \sum_{b}\sum_{p}
C>    C^{D\sigma"}_{bm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    h_{eb} \\\\
C>    &&+ \sum_{a}\sum_{p}
C>    C^{D\sigma"}_{am}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    h_{ae} 
C>    \;\;\;\;\;\;\;\;\;\;(6) \\\\
C> \f}
      subroutine wfn1_1el(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,h1,
     +                    En,da,db,docca,doccb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$ electrons
                  !< for each state
      integer neb !< [Input] The number of \f$\beta\f$ electrons
                  !< for each state
c
      double precision r_noa(nbf,nbf) !< [Input] The \f$\alpha\f$
                                      !< electron natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] The \f$\beta\f$
                                      !< electron natural orbitals
      double precision r_ofa(nbf,nbf) !< [Input] The \f$\alpha\f$
                                      !< electron occupation
                                      !< functions
      double precision r_ofb(nbf,nbf) !< [Input] The \f$\beta\f$
                                      !< electron occupation
                                      !< functions
c
      double precision h1(nbf,nbf) !< [Input] The actual
                                   !< 1-electron integrals
c
      double precision En !< [Output] The 1-electron
                          !< energy for each state
c
      double precision da(nbf,nbf)     !< [Scratch] Alpha-electron
                                       !< density matrix.
      double precision db(nbf,nbf)     !< [Scratch] Beta-electron
                                       !< density matrix.
      double precision docca(nbf)      !< [Scratch] Alpha-electron
                                       !< occupation numbers
      double precision doccb(nbf)      !< [Scratch] Beta-electron
                                       !< occupation numbers
c
c     Local variables
c
      integer ia, ib ! Counters over basis functions
      integer ii     ! Counter over natural orbitals
      integer ir     ! Counter over occupation functions
      integer ist    ! Counter over states
c
c     Functions
c
      double precision ddot
      external         ddot
c
c     Code
c
c     Compute the total 1-electron density matrices for this state
c
      do ii = 1, nbf
        docca(ii) = 0.0d0
        doccb(ii) = 0.0d0
      enddo
      do ir = 1, nea
        do ii = 1, nbf
          docca(ii) = docca(ii) + r_ofa(ii,ir)*r_ofa(ii,ir)
        enddo
      enddo
      do ir = 1, neb
        do ii = 1, nbf
          doccb(ii) = doccb(ii) + r_ofb(ii,ir)*r_ofb(ii,ir)
        enddo
      enddo
      call dfill(nbf*nbf,0.0d0,da,1)
      call dfill(nbf*nbf,0.0d0,db,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = da(ia,ib)
     +                + r_noa(ia,ii)*r_noa(ib,ii)*docca(ii)
            db(ia,ib) = db(ia,ib)
     +                + r_nob(ia,ii)*r_nob(ib,ii)*doccb(ii)
          enddo
        enddo
      enddo
c
c     Contract the 1-electron density matrix with the
c     1-electron integrals and add the energy contribution
c
      En = En + ddot(nbf*nbf,da,1,h1,1)
     +        + ddot(nbf*nbf,db,1,h1,1)
c
      end
C> @}
