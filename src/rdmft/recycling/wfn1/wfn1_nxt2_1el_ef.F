C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Calculate the 1-electron energy and the Fork matrix
C> within the wfn1 formalism
C> 
C> This routine adds the 1-electron energy and the corresponding Fock
C> matrix terms within the wfn1 formalism to its arguments. The current
C> formulation assumes that all quantities are kept in memory and the
C> implementation supports only serial execution.
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
      subroutine wfn1_nxt2_1el_ef(nbf,nea,neb,h1,occa,occb,da,db,
     +           oa,na,ca,ob,nb,cb,fa,fb,E1el)
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
      double precision h1(nbf,nbf) !< [Input] The actual
                                   !< 1-electron integrals
c
      double precision occa(nbf)   !< [Scratch] The alpha occupation
                                   !< numbers
      double precision occb(nbf)   !< [Scratch] The beta occupation
                                   !< numbers
      double precision da(nbf,nbf) !< [Scratch] The alpha density
                                   !< matrix
      double precision db(nbf,nbf) !< [Scratch] The beta density
                                   !< matrix
c
      double precision oa(nbf,nbf) !< [Input] The alpha occupation
                                   !< functions
      double precision ob(nbf,nbf) !< [Input] The beta occupation
                                   !< functions
      double precision na(nbf,nbf) !< [Input] The alpha natural
                                   !< orbitals
      double precision nb(nbf,nbf) !< [Input] The beta natural
                                   !< orbitals
      double precision ca(nbf)     !< [Input] The alpha electron
                                   !< configuration
      double precision cb(nbf)     !< [Input] The beta electron
                                   !< configuration
c
      double precision fa(nbf,nbf) !< [In/Output] The alpha Fock
                                   !< matrix
      double precision fb(nbf,nbf) !< [In/Output] The beta Fock
                                   !< matrix
c
      double precision E1el !< [Output] The 1-electron
                            !< energy for each state
c
c     Local variables
c
      integer ia, ib ! Counters over basis functions
      integer ii     ! Counter over natural orbitals
      integer ir     ! Counter over occupation functions
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
        occa(ii) = 0.0d0
        occb(ii) = 0.0d0
      enddo
      do ir = 1, nbf
        do ii = 1, nbf
          occa(ii) = occa(ii) + oa(ii,ir)*oa(ii,ir)*ca(ir)
        enddo
      enddo
      do ir = 1, nbf
        do ii = 1, nbf
          occb(ii) = occb(ii) + ob(ii,ir)*ob(ii,ir)*cb(ir)
        enddo
      enddo
      call dfill(nbf*nbf,0.0d0,da,1)
      call dfill(nbf*nbf,0.0d0,db,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = da(ia,ib)
     +                + na(ia,ii)*na(ib,ii)*occa(ii)
            db(ia,ib) = db(ia,ib)
     +                + nb(ia,ii)*nb(ib,ii)*occb(ii)
          enddo
        enddo
      enddo
c
c     Contract the 1-electron density matrix with the
c     1-electron integrals and add the energy contribution
c
      E1el = ddot(nbf*nbf,da,1,h1,1)
     +     + ddot(nbf*nbf,db,1,h1,1)
c
c     Add the contributions to the Fock matrices
c
      call daxpy(nbf*nbf,1.0d0,h1,1,fa,1)
      call daxpy(nbf*nbf,1.0d0,h1,1,fb,1)
c
      end
C> @}
