C> \ingroup wfn1
C> @{
C>
C> \brief Calculate \f$ S^2 \f$ for the current wavefunction
C>
C> According to Szabo and Ostlund [1] the expectation value of \f$S^2\f$
C> of an unrestricted wavefunction can be calculated as:
C> \f{eqnarray*}{
C>   \left\langle S^2 \right\rangle &=& S^2_{\mathrm{Exact}} + N^{\beta} +
C>   \sum_{i=1}^{N^\alpha}\sum_{j=1}^{N^\beta}
C>   \left\langle \psi_i^\alpha\right|\left.\psi_j^\beta\right\rangle^2
C> \f}
C> Here we express this equation in terms of 1-electron density matrices
C> as 
C> \f{eqnarray*}{
C>   \left\langle S^2 \right\rangle &=& S^2_{\mathrm{Exact}} + N^{\beta}
C>   \sum_{a,b,c,d,e}
C>   S_{ab}D^\alpha_{bc}S_{cd}D^\beta_{de}
C> \f}
C>
      double precision function wfn1_s2(nbf,nea,neb,ofa,ofb,noa,nob,ov,
     +                                  m1,m2,dd,occ)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$ \alpha \f$-electrons
      integer neb !< [Input] The number of \f$ \beta \f$-electrons
c
      double precision ofa(nbf,nbf) !< [Input] The \f$\alpha\f$
      !< occupation functions
      double precision ofb(nbf,nbf) !< [Input] The \f$\beta\f$
      !< occupation functions
      double precision noa(nbf,nbf) !< [Input] The \f$\alpha\f$
      !< natural orbitals
      double precision nob(nbf,nbf) !< [Input] The \f$\beta\f$
      !< natural orbitals
      double precision ov(nbf,nbf)  !< [Input] The overlap matrix
c
      double precision m1(nbf,nbf)  !< [Scratch] A matrix
      double precision m2(nbf,nbf)  !< [Scratch] A matrix
      double precision dd(nbf,nbf)  !< [Scratch] density matrix
      double precision occ(nbf)     !< [Scratch] occupation vector
c
      double precision sz !< Sz spin quantum number
      double precision t2 !< Temporary variable for S2
c
      integer ia, ib !< Counters
c
      sz = 0.5d0*(nea-neb)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,nob,nbf,0.0d0,m1,nbf)
      call dgemm('T','N',nbf,nbf,nbf,1.0d0,noa,nbf,m1,nbf,0.0d0,m2,nbf)
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,m2,nbf,ofb,nbf,0.0d0,m1,nbf)
      call dgemm('T','N',nbf,nbf,nbf,1.0d0,ofa,nbf,m1,nbf,0.0d0,m2,nbf)
c
      t2 = 0.0d0
      do ib = 1, neb
        do ia = 1, nea
          t2 = t2 + m2(ia,ib)**2
        enddo
      enddo
c
      wfn1_s2 = sz*(sz+1.0d0) + neb - t2
c
      return
      end
C>
C> @}

