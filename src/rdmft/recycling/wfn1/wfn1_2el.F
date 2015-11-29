C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the Coulomb energy and optionally the gradient
C> within the wfn1 formalism
C> 
C> This routine adds the Coulomb energy and optionally the corresponding
C> gradient terms within the wfn1 formalism to its arguments. The
C> current formulation assumes that all quantities are kept in memory
C> and the implementation supports only serial execution.
C>
C> Within the wfn1 formalism the Coulomb energy is expressed in terms
C> of the occupation functions and the natural orbitals as:
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
C> Although the energy can be summarized this way the gradient of the
C> energy cannot. The gradient requires the derivatives with respect to
C> the occupation functions as well as the natural orbitals. I.e.
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{2C}}{\mathrm{d}C^{O\sigma"}_{mt}}
C>    &=& \sum_{\sigma\sigma'}\sum_{abcd}
C>    \sum_{ij}\sum_{pq} \delta_{\sigma\sigma"}\delta_{im}\delta_{pt}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq}\\\\
C>    &&+ 
C>    \sum_{\sigma\sigma'}\sum_{abcd}
C>    \sum_{ij}\sum_{pq} \delta_{\sigma'\sigma"}\delta_{jm}\delta_{qt}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}
C>    \;\;\;\;\;\;\;\;\;\;(3)\\\\
C>    &=& \sum_{\sigma'}\sum_{abcd}
C>    \sum_{j}\sum_{q}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq}\\\\
C>    &&+ 
C>    \sum_{\sigma}\sum_{abcd}
C>    \sum_{i}\sum_{p}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt}
C>    \;\;\;\;\;\;\;\;\;\;(4)\\\\
C>    &=& \sum_{\sigma'}\sum_{abcd}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    (ab|cd)D^{T\sigma'}_{cd} \\\\
C>    &&+ 
C>    \sum_{\sigma}\sum_{abcd}
C>    D^{T\sigma}_{ab}
C>    (ab|cd)
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt}
C>    \;\;\;\;\;\;\;\;\;\;(5)\\\\
C>    &=& \sum_{\sigma'}\sum_{ab}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    F^{\sigma'}_{ab} \\\\
C>    &&+ 
C>    \sum_{\sigma}\sum_{cd}
C>    F^{\sigma}_{cd}
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt},
C>    \;\;\;\;\;\;\;\;\;\;(6)\\\\
C>    F^{\sigma'}_{ab} &=& \sum_{cd}(ab|cd)D^{T\sigma'}_{cd}, \\\\
C>    F^{\sigma}_{cd} &=& \sum_{ab}D^{T\sigma}_{ab}(ab|cd) \\\\
C> \f}
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{2C}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& \frac{1}{2}\sum_{\sigma\sigma'}\sum_{abcd}\sum_{ij}\sum_{pq}
C>    \delta_{\sigma\sigma"}\delta_{ae}\delta_{im}
C>    C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma\sigma'}\sum_{abcd}\sum_{ij}\sum_{pq}
C>    \delta_{\sigma\sigma"}\delta_{be}\delta_{im}
C>    C^{D\sigma}_{ai}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma\sigma'}\sum_{abcd}\sum_{ij}\sum_{pq}
C>    \delta_{\sigma\sigma"}\delta_{ce}\delta_{jm}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{dj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma\sigma'}\sum_{abcd}\sum_{ij}\sum_{pq}
C>    \delta_{\sigma\sigma"}\delta_{de}\delta_{jm}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma'}_{cj}
C>    C^{O\sigma'}_{jq}C^{O\sigma'}_{jq}
C>    \;\;\;\;\;\;\;\;\;\;(7) \\\\
C>    &=& \frac{1}{2}\sum_{\sigma'}\sum_{bcd}\sum_{p}
C>    C^{D\sigma"}_{bm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    (eb|cd)D^{T\sigma'}_{cd} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma'}\sum_{acd}\sum_{p}
C>    C^{D\sigma"}_{am}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    (ae|cd)D^{T\sigma'}_{cd} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma}\sum_{abd}\sum_{q}
C>    D^{T\sigma}_{ab}(ab|ed)
C>    C^{D\sigma"}_{dm}C^{O\sigma"}_{mq}C^{O\sigma"}_{mq} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma}\sum_{abc}\sum_{q}
C>    D^{T\sigma}_{ab}(ab|ce)
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mq}C^{O\sigma"}_{mq}
C>    \;\;\;\;\;\;\;\;\;\;(8) \\\\
C>    &=& \frac{1}{2}\sum_{\sigma'}\sum_{b}\sum_{p}
C>    C^{D\sigma"}_{bm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    F^{\sigma'}_{ab}\delta_{ae} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma'}\sum_{a}\sum_{p}
C>    C^{D\sigma"}_{am}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    F^{\sigma'}_{ab}\delta_{be} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma}\sum_{d}\sum_{q}
C>    F^{\sigma}_{cd}\delta_{ce}
C>    C^{D\sigma"}_{dm}C^{O\sigma"}_{mq}C^{O\sigma"}_{mq} \\\\
C>    &&+ \frac{1}{2}\sum_{\sigma}\sum_{c}\sum_{q}
C>    F^{\sigma}_{cd}\delta_{de}
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mq}C^{O\sigma"}_{mq}
C>    \;\;\;\;\;\;\;\;\;\;(9) \\\\
C>    F^{\sigma'}_{ab} &=& \sum_{cd}(ab|cd)D^{T\sigma'}_{cd}, \\\\
C>    F^{\sigma}_{cd} &=& \sum_{ab}D^{T\sigma}_{ab}(ab|cd) \\\\
C> \f}
      subroutine wfn1_2el_C(nbf,nea,neb,r_noa,r_nob,
     +                      r_ofa,r_ofb,eri,
     +                      En,da,db,d_fa,d_fb,docca,doccb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nst !< [Input] The number of states
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
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The actual
                                            !< 2-electron integrals
c
      double precision En !< [Output] The 2-electron Coulomb 
                          !< energy for each state
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
c     Contract the density matrices with the 2-electron integrals
c     to generate Fock-like quantities
c
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +           0.0d0,d_fa,1)
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +           0.0d0,d_fb,1)
      En = En + ddot(nbf*nbf,da,1,d_fa,1)/2.0d0
      En = En + ddot(nbf*nbf,da,1,d_fb,1)/2.0d0
      En = En + ddot(nbf*nbf,db,1,d_fa,1)/2.0d0
      En = En + ddot(nbf*nbf,db,1,d_fb,1)/2.0d0
c
      end
C>
C> \brief Calculate the Coulomb self repulsion energy within the wfn1
C> formalism
C> 
C> This routine adds the Coulomb self repulsion energy within the wfn1
C> formalism to its arguments. The current formulation assumes that all
C> quantities are kept in memory and the implementation supports only
C> serial execution.
C>
C> Within the wfn1 formalism the exchange energy is expressed in terms
C> of the occupation functions and the natural orbitals as:
C> \f{eqnarray*}{
C>  E^{2X} &=& \frac{1}{2}\sum_{\sigma}\sum_{i}
C>  \sum_{pq\in\{occ\}}
C>  C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{iq}
C>  (ab|cd)
C>  C^{D\sigma}_{ci}C^{D\sigma}_{di}C^{O\sigma}_{iq}C^{O\sigma}_{ip}
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> Obviously this energy expression can be summarized in terms of the
C> single 1-electron density matrix and using the orthonormality of the
C> occupation functions as
C> \f{eqnarray*}{
C>    E^{2X} &=& \frac{1}{2}\sum_{\sigma}\sum_{p\in\{occ\}}\sum_{abcd}
C>    D^{p\sigma}_{ab}(ab|cd)D^{p\sigma}_{cd},
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>    D^{p\sigma}_{ab} &=& \sum_i 
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C> \f}
C> Although the energy can be summarized this way the gradient of the
C> energy cannot. The gradient requires the derivatives with respect to
C> the occupation functions as well as the natural orbitals. I.e.
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{2X}}{\mathrm{d}C^{O\sigma"}_{mt}}
C>    &=& \sum_{\sigma}\sum_{abcd}
C>    \sum_{i}\sum_{p} \delta_{\sigma\sigma"}\delta_{im}\delta_{pt}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{ci}C^{D\sigma}_{di}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip}\\\\
C>    &&+ 
C>    \sum_{\sigma}\sum_{abcd}
C>    \sum_{i}\sum_{p} \delta_{\sigma\sigma"}\delta_{jm}\delta_{pt}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{ci}C^{D\sigma}_{di}
C>    C^{O\sigma}_{ip}
C>    \;\;\;\;\;\;\;\;\;\;(3)\\\\
C>    &=& \sum_{abcd}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    (ab|cd)
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt}C^{O\sigma"}_{mt}\\\\
C>    &&+ 
C>    \sum_{abcd}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}
C>    C^{O\sigma"}_{mt}C^{O\sigma"}_{mt}
C>    (ab|cd)
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt}
C>    \;\;\;\;\;\;\;\;\;\;(4)\\\\
C>    &=& \sum_{abcd}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    (ab|cd)D^{t\sigma"}_{cd} \\\\
C>    &&+ 
C>    \sum_{abcd}
C>    D^{t\sigma"}_{ab}
C>    (ab|cd)
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt}
C>    \;\;\;\;\;\;\;\;\;\;(5)\\\\
C>    &=& \sum_{ab}
C>    C^{D\sigma"}_{am}C^{D\sigma"}_{bm}C^{O\sigma"}_{mt}
C>    F^{\sigma"}_{ab} \\\\
C>    &&+ 
C>    \sum_{cd}
C>    F^{\sigma"}_{cd}
C>    C^{D\sigma"}_{cm}C^{D\sigma"}_{dm}
C>    C^{O\sigma"}_{mt},
C>    \;\;\;\;\;\;\;\;\;\;(6)\\\\
C>    F^{\sigma"}_{ab} &=& \sum_{cd}(ab|cd)D^{t\sigma"}_{cd}, \\\\
C>    F^{\sigma"}_{cd} &=& \sum_{ab}D^{t\sigma"}_{ab}(ab|cd) \\\\
C> \f}
C> \f{eqnarray*}{
C>    \frac{\mathrm{d}E^{2X}}{\mathrm{d}C^{D\sigma"}_{em}}
C>    &=& \sum_{\sigma}\sum_{abcd}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{ae}\delta_{im}
C>    C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{ci}C^{D\sigma}_{di}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip} \\\\
C>    &&+ \sum_{\sigma}\sum_{abcd}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{be}\delta_{im}
C>    C^{D\sigma}_{ai}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{ci}C^{D\sigma}_{di}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip} \\\\
C>    &&+ \sum_{\sigma}\sum_{abcd}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{ce}\delta_{im}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{di}
C>    C^{O\sigma}_{ip}C^{O\sigma}_{ip} \\\\
C>    &&+ \sum_{\sigma}\sum_{abcd}\sum_{i}\sum_{p}
C>    \delta_{\sigma\sigma"}\delta_{de}\delta_{im}
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>    (ab|cd)
C>    C^{D\sigma}_{ci}
C>    C^{O\sigma}_{iq}C^{O\sigma}_{ip}
C>    \;\;\;\;\;\;\;\;\;\;(7) \\\\
C>    &=& \sum_{bcd}\sum_{p}
C>    C^{D\sigma"}_{bm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    (eb|cd)D^{p\sigma"}_{cd} \\\\
C>    &&+ \sum_{acd}\sum_{p}
C>    C^{D\sigma"}_{am}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    (ae|cd)D^{p\sigma"}_{cd} \\\\
C>    &&+ \sum_{abd}\sum_{p}
C>    D^{p\sigma"}_{ab}(ab|ed)
C>    C^{D\sigma"}_{dm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp} \\\\
C>    &&+ \sum_{abc}\sum_{p}
C>    D^{p\sigma"}_{ab}(ab|ce)
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    \;\;\;\;\;\;\;\;\;\;(8) \\\\
C>    &=& \sum_{b}\sum_{p}
C>    C^{D\sigma"}_{bm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    F^{\sigma"}_{ab}\delta_{ae} \\\\
C>    &&+ \sum_{a}\sum_{p}
C>    C^{D\sigma"}_{am}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp}
C>    F^{\sigma"}_{ab}\delta_{be} \\\\
C>    &&+ \sum_{d}\sum_{p}
C>    F^{\sigma"}_{cd}\delta_{ce}
C>    C^{D\sigma"}_{dm}C^{O\sigma"}_{mp}C^{O\sigma"}_{mp} \\\\
C>    &&+ \sum_{c}\sum_{p}
C>    F^{\sigma"}_{cd}\delta_{de}
C>    C^{D\sigma"}_{cm}C^{O\sigma"}_{mq}C^{O\sigma"}_{mq}
C>    \;\;\;\;\;\;\;\;\;\;(9) \\\\
C>    F^{\sigma"}_{ab} &=& \sum_{cd}(ab|cd)D^{p\sigma"}_{cd}, \\\\
C>    F^{\sigma"}_{cd} &=& \sum_{ab}D^{p\sigma"}_{ab}(ab|cd) \\\\
C> \f}
      subroutine wfn1_2el_Csr(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,eri,
     +                        En,da,db,d_fa,d_fb,docca,doccb)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer nst !< [Input] The number of states
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
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The actual
                                            !< 2-electron integrals
c
      double precision En !< [Output] The 2-electron exchange 
                          !< energy for each state
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
      double precision Et
      double precision Eij
c
c     Local variables
c
      integer ia, ib ! Counters over basis functions
      integer ii     ! Counter over natural orbitals
      integer ip     ! Counter over occupation functions
      integer iq     ! Counter over occupation functions
      integer ist    ! Counter over states
c
c     Functions
c
      double precision ddot
      external         ddot
c
c     Code
c
      Eij = 0.0d0
c
c     First do the alpha-electron part
c
      do ip = 1, nea
c       do iq = 1, nea
        do iq = ip, ip
          call dfill(nbf*nbf,0.0d0,da,1)
          do ib = 1, nbf
            do ia = 1, nbf
              do ii = 1, nbf
                da(ia,ib) = da(ia,ib)
     +                    + r_noa(ia,ii)*r_noa(ib,ii)*
     +                      r_ofa(ii,ip)*r_ofa(ii,iq)
              enddo
            enddo
          enddo
c
c         Contract the single electron density matrix with the
c         2-electron integrals to generate Fock-like quantities
c
          call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +               0.0d0,d_fa,1)
          Et = ddot(nbf*nbf,da,1,d_fa,1)/8.0d0
          En = En + Et
          if (ip.ne.iq) Eij = Eij + Et
c
        enddo ! iq
      enddo ! ip
c
c     Next do the beta-electron part
c
      do ip = 1, neb
c       do iq = 1, neb
        do iq = ip, ip
          call dfill(nbf*nbf,0.0d0,db,1)
          do ib = 1, nbf
            do ia = 1, nbf
              do ii = 1, nbf
                db(ia,ib) = db(ia,ib)
     +                    + r_nob(ia,ii)*r_nob(ib,ii)*
     +                      r_ofb(ii,ip)*r_ofb(ii,iq)
              enddo
            enddo
          enddo
c
c         Contract the single electron density matrix with the
c         2-electron integrals to generate Fock-like quantities
c
          call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +               0.0d0,d_fb,1)
          Et = ddot(nbf*nbf,db,1,d_fb,1)/8.0d0
          En = En + Et
          if (ip.ne.iq) Eij = Eij + Et
c
        enddo ! iq
      enddo ! ip
c
      end
C>
C> \brief Calculate the exchange energy as in Hartree-Fock
C>
C> The Hartree-Fock exchange energy is different in that it is computed
C> from the total alpha and beta density matrices. This formulation
C> is similar to the Coulomb term apart from the fact that alpha
C> electrons interact only with alpha electrons, likewise for the beta
C> electrons, and finally the 2-electron integrals are ordered 
C> differently, i.e. the j and l labels are swapped.
C>
      subroutine wfn1_2el_hfX(nbf,nea,neb,r_noa,r_nob,r_ofa,r_ofb,ov,
     +                        eri,En,da,db,d_fa,d_fb,docca,doccb,power)
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
      double precision ov(nbf,nbf)    !< [Input] The overlap matrix
c
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The actual
                                            !< 2-electron integrals
c
      double precision En !< [Output] The 2-electron Coulomb 
                          !< energy for each state
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
      double precision power           !< [Input] Matrix power
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
      do ii = 1, nbf
        docca(ii) = docca(ii)**power
        doccb(ii) = doccb(ii)**power
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
c     If power <> 1.0d0 then compute the power of the density
c     matrices
c
c     if (power.ne.1.0d0) then
c       call wfn1_mat_power(nbf,da,nbf,ov,nbf,power)
c       call wfn1_mat_power(nbf,db,nbf,ov,nbf,power)
c     endif
c
c     Contract the density matrices with the 2-electron integrals
c     to generate Fock-like quantities
c
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +           0.0d0,d_fa,1)
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +           0.0d0,d_fb,1)
      En = En - ddot(nbf*nbf,da,1,d_fa,1)/2.0d0
      En = En - ddot(nbf*nbf,db,1,d_fb,1)/2.0d0
c
      end
C> @}
