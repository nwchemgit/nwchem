C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the Coulomb orbital self repulsion energy within
C> the wfn1 formalism
C> 
C> This routine adds the Coulomb orbital self repulsion energies within
C> the wfn1 formalism to its arguments. The current formulation assumes
C> that all quantities are kept in memory and the implementation
C> supports only serial execution.
C>
C> Within the wfn1 formalism the Coulomb orbital energy is expressed in
C> terms of the occupation functions and the natural orbitals as:
C> \f{eqnarray*}{
C>  E^{2C} &=& \sum_{\sigma,\sigma'}\sum_{ij}
C>  \sum_{pq\in\{occ\}}\frac{1}{2}
C>  C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>  (ab|cd)
C>  C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}C^{O\sigma'}_{jq}C^{O\sigma'}_{jq}
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> Obviously this energy expression can be summarized in terms of the
C> total 1-electron density matrix as
C> \f{eqnarray*}{
C>    E^{2C} &=& \sum_{\sigma\sigma'}\sum_{abcd}\frac{1}{2}
C>    D^{T\sigma}_{ab}(ab|cd)D^{T\sigma'}_{cd}
C>    \;\;\;\;\;\;\;\;\;\;(2)
C> \f}
      subroutine wfn1_2el_Cslf_orb(nbf,nea,neb,r_noa,r_nob,
     +                             r_ofa,r_ofb,eri,
     +                             Ea,Eb,da,db,d_fa,d_fb,docca,doccb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$ electrons
      integer neb !< [Input] The number of \f$\beta\f$ electrons
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
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The actual
                                            !< 2-electron integrals
c
      double precision Ea(nbf) !< [Output] The 2-electron Coulomb 
                               !< energy for each \f$\alpha\f$ orbital
      double precision Eb(nbf) !< [Output] The 2-electron Coulomb 
                               !< energy for each \f$\beta\f$ orbital
c
      double precision da(nbf,nbf)     !< [Scratch] Alpha-electron
                                       !< density matrix.
      double precision db(nbf,nbf)     !< [Scratch] Beta-electron
                                       !< density matrix.
      double precision d_fa(nbf,nbf)   !< [Scratch] Alpha-electron
                                       !< Fock like matrix.
      double precision d_fb(nbf,nbf)   !< [Scratch] Beta-electron
                                       !< Fock like matrix.
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
c
c     Functions
c
      double precision ddot
      external         ddot
c
c     Code
c
c     Calculate the alpha single-electron density matrices
c
      do ir = 1, nea
        do ii = 1, nbf
          docca(ii) = r_ofa(ii,ir)*r_ofa(ii,ir)
        enddo
        call dfill(nbf*nbf,0.0d0,da,1)
        do ib = 1, nbf
          do ia = 1, nbf
            do ii = 1, nbf
              da(ia,ib) = da(ia,ib)
     +                  + r_noa(ia,ii)*r_noa(ib,ii)*docca(ii)
            enddo
          enddo
        enddo
c
c       Contract the density matrices with the 2-electron integrals
c       to generate Fock-like quantities
c
        call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +             0.0d0,d_fa,1)
c
c       Calculate and subtract the self-interaction energy
c
        Ea(ir) = Ea(ir) - ddot(nbf*nbf,da,1,d_fa,1)
c
      enddo
c
c     Calculate the beta single-electron density matrices
c
      do ir = 1, neb
        do ii = 1, nbf
          doccb(ii) = r_ofb(ii,ir)*r_ofb(ii,ir)
        enddo
        call dfill(nbf*nbf,0.0d0,db,1)
        do ib = 1, nbf
          do ia = 1, nbf
            do ii = 1, nbf
              db(ia,ib) = db(ia,ib)
     +                  + r_nob(ia,ii)*r_nob(ib,ii)*doccb(ii)
            enddo
          enddo
        enddo
c
c       Contract the density matrices with the 2-electron integrals
c       to generate Fock-like quantities
c
        call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +             0.0d0,d_fb,1)
c
c       Calculate and subtract the self-interaction energy
c
        Eb(ir) = Eb(ir) - ddot(nbf*nbf,db,1,d_fb,1)
c
      enddo
c
      end
C>
C> \brief Calculate the self exchange orbital energy as in Hartree-Fock
C>
C> The Hartree-Fock exchange orbital energy is different in that it is
C> computed from the total alpha and beta density matrices. This
C> formulation is similar to the Coulomb term apart from the fact that
C> alpha electrons interact only with alpha electrons, likewise for the
C> beta electrons, and finally the 2-electron integrals are ordered 
C> differently, i.e. the j and l labels are swapped.
C>
      subroutine wfn1_2el_hfXslf_orb(nbf,nea,neb,r_noa,r_nob,
     +                               r_ofa,r_ofb,eri,Ea,Eb,da,db,
     +                               d_fa,d_fb,docca,doccb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$ electrons
      integer neb !< [Input] The number of \f$\beta\f$ electrons
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
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The actual
                                            !< 2-electron integrals
c
      double precision Ea(nbf) !< [Output] The 2-electron exchange 
                               !< energy for each alpha orbital
      double precision Eb(nbf) !< [Output] The 2-electron exchange 
                               !< energy for each beta orbital
c
      double precision da(nbf,nbf)     !< [Scratch] Alpha-electron
                                       !< density matrix.
      double precision db(nbf,nbf)     !< [Scratch] Beta-electron
                                       !< density matrix.
      double precision d_fa(nbf,nbf)   !< [Scratch] Alpha-electron
                                       !< Fock like matrix.
      double precision d_fb(nbf,nbf)   !< [Scratch] Beta-electron
                                       !< Fock like matrix.
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
c
c     Functions
c
      double precision ddot
      external         ddot
c
c     Code
c
c     Calculate the alpha single-electron density matrices
c
      do ir = 1, nea
        do ii = 1, nbf
          docca(ii) = r_ofa(ii,ir)*r_ofa(ii,ir)
        enddo
        call dfill(nbf*nbf,0.0d0,da,1)
        do ib = 1, nbf
          do ia = 1, nbf
            do ii = 1, nbf
              da(ia,ib) = da(ia,ib)
     +                  + r_noa(ia,ii)*r_noa(ib,ii)*docca(ii)
            enddo
          enddo
        enddo
c
c       Contract the density matrices with the 2-electron integrals
c       to generate Fock-like quantities
c
        call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +             0.0d0,d_fa,1)
c
        Ea(ir) = Ea(ir) + ddot(nbf*nbf,da,1,d_fa,1)
c
      enddo
c
c     Calculate the beta single-electron density matrices
c
      do ir = 1, neb
        do ii = 1, nbf
          doccb(ii) = r_ofb(ii,ir)*r_ofb(ii,ir)
        enddo
        call dfill(nbf*nbf,0.0d0,db,1)
        do ib = 1, nbf
          do ia = 1, nbf
            do ii = 1, nbf
              db(ia,ib) = db(ia,ib)
     +                  + r_nob(ia,ii)*r_nob(ib,ii)*doccb(ii)
            enddo
          enddo
        enddo
c
c       Contract the density matrices with the 2-electron integrals
c       to generate Fock-like quantities
c
        call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +             0.0d0,d_fb,1)
c
        Eb(ir) = Eb(ir) + ddot(nbf*nbf,db,1,d_fb,1)
c
      enddo
c
      end
C> @}
