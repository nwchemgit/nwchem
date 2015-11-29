C> \ingroup wfn1
C> @{
C>
C> \brief Select the active natural orbitals by occupation number
C>
C> In order to ensure that the full energy expression (the energy
C> expression including all the Lagrangian terms) is bound from below
C> the Lagrangians can only contain functions that are present in the
C> physical part of the energy expression. This means that we need to
C> select only the occupied natural orbitals to be part of the 
C> Lagrangians. Here all natural orbitals with an occupation number
C> exceeding the specified threshold are selected to be included in
C> the Lagrangian.
C>
      subroutine wfn1_occ_natorb(nbf,nst,nea,neb,r_ofa,r_ofb,
     +                           l_occa,l_occb,tol_occ)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nst !< [Input] The number of states
      integer nea(nst) !< [Input] The number of \f$\alpha\f$-electrons
                       !< per state
      integer neb(nst) !< [Input] The number of \f$\beta\f$-electrons
                       !< per state
c
      double precision r_ofa(nbf,nbf,nst) !< [Input] The 
      !< \f$\alpha\f$-occupation functions
      double precision r_ofb(nbf,nbf,nst) !< [Input] The
      !< \f$\beta\f$-occupation functions
      double precision tol_occ !< [Input] Threshold for active 
      !< natural orbitals
c
      logical l_occa(nbf,nst) !< [Output] .True. for every occupied
      !< \f$\alpha\f$-electron natural orbital
      logical l_occb(nbf,nst) !< [Output] .True. for every occupied
      !< \f$\beta\f$-electron natural orbital
c
      integer ist !< Counter over states
      integer ii  !< Counter over natural orbitals
      integer ir  !< Counter over occupation functions
c
      double precision occ !< The occupation number
c
      do ist = 1, nst
        do ii = 1, nbf
c
          occ = 0.0d0
          do ir = 1, nea(ist)
            occ = occ + r_ofa(ii,ir,ist)*r_ofa(ii,ir,ist)
          enddo
          l_occa(ii,ist) = occ.gt.tol_occ
c
          occ = 0.0d0
          do ir = 1, neb(ist)
            occ = occ + r_ofb(ii,ir,ist)*r_ofb(ii,ir,ist)
          enddo
          l_occb(ii,ist) = occ.gt.tol_occ
c
        enddo
      enddo
c
      end
C> @}
