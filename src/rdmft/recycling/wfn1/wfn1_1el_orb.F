C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the 1-electron energy for all extended orbitals
C> 
C> This routine adds the 1-electron energy for every extended orbital
C> within the wfn1 formalism to its arguments. The current formulation
C> assumes that all quantities are kept in memory and the implementation
C> supports only serial execution.
C>
C> Within the wfn1 formalism the 1-electron energy for a given orbital
C> is expressed in terms of the occupation functions and the natural
C> orbitals as:
C> \f{eqnarray*}{
C>  E^{1\sigma}(p) &=& \sum_{i}
C>  \sum_{ab}
C>  C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C>  h_{ab}
C>  \;\;\;\;\;\;\;\;\;\;(1)
C> \f}
C> Obviously this energy expression can be summarized in terms of the
C> total 1-electron density matrix as
C> \f{eqnarray*}{
C>    E^{1\sigma}(p) &=& \sum_{\sigma}\sum_{ab}
C>    D^{T\sigma}_{ab}(p)h_{ab},
C>    \;\;\;\;\;\;\;\;\;\;(2) \\\\
C>    D^{T\sigma}_{ab}(p) &=& \sum_i 
C>    C^{D\sigma}_{ai}C^{D\sigma}_{bi}C^{O\sigma}_{ip}C^{O\sigma}_{ip}
C> \f}
      subroutine wfn1_1el_orb(nbf,ne,r_no,r_of,h1,
     +                    E1,da,docc)
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
c
      integer nbf !< [Input] The number of basis functions
      integer ne  !< [Input] The number of occupied functions in a spin
                  !< channel
c
      double precision r_no(nbf,nbf) !< [Input] The spin channel
                                     !< electron natural orbitals
      double precision r_of(nbf,nbf) !< [Input] The spin channel
                                     !< electron occupation
                                     !< functions
c
      double precision h1(nbf,nbf) !< [Input] The actual
                                   !< 1-electron integrals
c
      double precision E1(nbf) !< [Output] The 1-electron
                               !< energy for extended orbitals
c
      double precision da(nbf,nbf)     !< [Scratch] single-electron
                                       !< density matrix.
      double precision docc(nbf)       !< [Scratch] Alpha-electron
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
      do ir = 1, nbf
c
c       Compute the total 1-electron density matrices for this state
c
        do ii = 1, nbf
          docc(ii) = 0.0d0
        enddo
        do ii = 1, nbf
          docc(ii) = docc(ii) + r_of(ii,ir)*r_of(ii,ir)
        enddo
        call dfill(nbf*nbf,0.0d0,da,1)
        do ib = 1, nbf
          do ia = 1, nbf
            do ii = 1, nbf
              da(ia,ib) = da(ia,ib)
     +                  + r_no(ia,ii)*r_no(ib,ii)*docc(ii)
            enddo
          enddo
        enddo
c
c       Contract the 1-electron density matrix with the
c       1-electron integrals and add the energy contribution
c
        E1(ir) = E1(ir) + ddot(nbf*nbf,da,1,h1,1)
c
      enddo
c
      end
C> @}
