C> \ingroup wfn1
C> @{
C>
C> \brief Add \f$\alpha\f$ and \f$\beta\f$ density matrices
C>
C> This routine calculates the \f$\alpha\f$ and \f$\beta\f$ density
C> matrix from a specified WFN1 wavefunction. The density matrices
C> are multiplied by a provided weight and added onto the matrices
C> in the argument list.
C>
C> In equations this routine calculates
C> \f{eqnarray*}{
C>    D^\sigma_{ab} &=& w2*D^\sigma_{ab} + w1*D^{\sigma+}_{ab}
C> \f}
C> where \f$D^\sigma\f$ is a density matrix given in the argument and
C> \f$D^{\sigma+}\f$ is the term calculated from the specified 
C> wavefunction.
C>
      subroutine wfn1_add_dmat(nbf,nea,neb,w1,r_noa,r_nob,r_ofa,r_ofb,
     +                         w2,da,db,docca,doccb)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$-electrons
      integer neb !< [Input] The number of \f$\beta\f$-electrons
c
      double precision w1 !< [Input] The weight of the new contributions
      double precision w2 !< [Input] The weight of the old density 
                          !< matrices
c
      double precision r_noa(nbf,nbf) !< [Input] The \f$\alpha\f$
      !< natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] The \f$\beta\f$
      !< natural orbitals
      double precision r_ofa(nbf,nbf) !< [Input] The \f$\alpha\f$
      !< occupation functions
      double precision r_ofb(nbf,nbf) !< [Input] The \f$\beta\f$
      !< occupation functions
c
      double precision da(nbf,nbf) !< [In/Output] The \f$\alpha\f$
      !< density matrix
      double precision db(nbf,nbf) !< [In/Output] The \f$\beta\f$
      !< density matrix
c
      double precision docca(nbf) !< [Scratch] The \f$\alpha\f$
      !< occupation numbers
      double precision doccb(nbf) !< [Scratch] The \f$\beta\f$
      !< occupation numbers
c
c     Local variables
c
      integer ii !< Counter over natural orbitals
      integer ir !< Counter over occupation functions
      integer ia !< Counter over basis functions
      integer ib !< Counter over basis functions
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
        docca(ii) = w1*docca(ii)
        doccb(ii) = w1*doccb(ii)
      enddo
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = w2*da(ia,ib)
     +                + r_noa(ia,ii)*r_noa(ib,ii)*docca(ii)
            db(ia,ib) = w2*db(ia,ib)
     +                + r_nob(ia,ii)*r_nob(ib,ii)*doccb(ii)
          enddo
        enddo
      enddo
c
      end
C>
C> \brief Add \f$\alpha\f$ and \f$\beta\f$ transition density matrices
C>
C> This routine takes a left and a right WFN1 wavefunction transformed
C> to a common basis. The overlap of the orbitals in the common basis
C> is given in the variables `va` and `vb`. The corresponding 
C> \f$\alpha\f$ and \f$\beta\f$ transition density matrices are 
C> calculated, multiplied by the weight `w1` and added to the old
C> density matrix multiplied with weight `w2`. It is assumed that
C> the vectors `va` and `vb` contain the overlaps in decreasing order.
C> This assumption has consequences for the numerical stability.
C>
C> Note that `w1` should include the parity of the orbital permutations
C> to ensure the correct sign overall.
C>
C> This routine exploits the fact that for all occupation functions
C> in the same spin channel the natural orbitals are the same.
C>
C> In equations this routine calculates
C> \f{eqnarray*}{
C>    D^\sigma_{ab} &=& w2*D^\sigma_{ab} + w1*D^{\sigma+}_{ab}
C> \f}
C> where \f$D^\sigma\f$ is a density matrix given in the argument and
C> \f$D^{\sigma+}\f$ is the term calculated from the specified 
C> wavefunction.
C>
      subroutine wfn1_add_tdmat(nbf,nea,neb,w1,va,vb,
     +                          l_noa,l_nob,l_ofa,l_ofb,
     +                          r_noa,r_nob,r_ofa,r_ofb,
     +                          w2,da,db,docca,doccb)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$-electrons
      integer neb !< [Input] The number of \f$\beta\f$-electrons
c
      double precision w1 !< [Input] The weight of the new contributions
      double precision w2 !< [Input] The weight of the old density 
                          !< matrices
c
      double precision va(nea) !< [Input] The overlap between the 
      !< left and right \f$\alpha\f$ orbitals in the common basis for
      !< all occupied functions.
      double precision vb(neb) !< [Input] The overlap between the 
      !< left and right \f$\beta\f$ orbitals in the common basis for
      !< all occupied functions.
c
      double precision l_noa(nbf,nbf) !< [Input] The left \f$\alpha\f$
      !< natural orbitals
      double precision l_nob(nbf,nbf) !< [Input] The left \f$\beta\f$
      !< natural orbitals
      double precision l_ofa(nbf,nbf) !< [Input] The left \f$\alpha\f$
      !< occupation functions
      double precision l_ofb(nbf,nbf) !< [Input] The left \f$\beta\f$
      !< occupation functions
      double precision r_noa(nbf,nbf) !< [Input] The right \f$\alpha\f$
      !< natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] The right \f$\beta\f$
      !< natural orbitals
      double precision r_ofa(nbf,nbf) !< [Input] The right \f$\alpha\f$
      !< occupation functions
      double precision r_ofb(nbf,nbf) !< [Input] The right \f$\beta\f$
      !< occupation functions
c
      double precision da(nbf,nbf) !< [In/Output] The \f$\alpha\f$
      !< density matrix
      double precision db(nbf,nbf) !< [In/Output] The \f$\beta\f$
      !< density matrix
c
      double precision docca(nbf) !< [Scratch] The \f$\alpha\f$
      !< occupation numbers
      double precision doccb(nbf) !< [Scratch] The \f$\beta\f$
      !< occupation numbers
c
c     Local variables
c
      integer ii !< Counter over natural orbitals
      integer ir !< Counter over occupation functions
      integer is !< Counter over occupation functions
      integer ia !< Counter over basis functions
      integer ib !< Counter over basis functions
c
      double precision sa !< The product of all `va` entries
      double precision sb !< The product of all `vb` entries
      double precision ta !< The product of all `va` entries but one
      double precision tb !< The product of all `vb` entries but one
c
c     Compute the total 1-electron density matrices for this state
c
      sa = va(1)
      sb = va(2)
      do ir = 2, nea
        sa = sa * va(ir)
      enddo
      do ir = 2, neb
        sb = sb * vb(ir)
      enddo
c
      do ii = 1, nbf
        docca(ii) = 0.0d0
        doccb(ii) = 0.0d0
      enddo
      do ir = 1, nea
        ta = 1.0d0
        do is = 1, nea
          if (is.ne.ir) ta = ta * va(is)
        enddo
        do ii = 1, nbf
          docca(ii) = docca(ii) + l_ofa(ii,ir)*r_ofa(ii,ir)*ta*sb
        enddo
      enddo
      do ir = 1, neb
        tb = 1.0d0
        do is = 1, neb
          if (is.ne.ir) tb = tb * vb(is)
        enddo
        do ii = 1, nbf
          doccb(ii) = doccb(ii) + l_ofb(ii,ir)*r_ofb(ii,ir)*tb*sa
        enddo
      enddo
      do ii = 1, nbf
        docca(ii) = w1*docca(ii)
        doccb(ii) = w1*doccb(ii)
      enddo
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = w2*da(ia,ib)
     +                + l_noa(ia,ii)*r_noa(ib,ii)*docca(ii)
            db(ia,ib) = w2*db(ia,ib)
     +                + l_nob(ia,ii)*r_nob(ib,ii)*doccb(ii)
          enddo
        enddo
      enddo
c
      end
C>
C> @}
