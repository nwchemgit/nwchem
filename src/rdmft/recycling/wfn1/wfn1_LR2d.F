C> \ingroup wfn1
C> @{
C>
C> \brief Transform the left and right wavefunction to the common
C> representation
C>
C> In order to compute transition density matrices the left and 
C> right wavefunctions must be transformed to a common representation.
C> This routine assumes that the left and right transformations 
C> \f$\mathrm{L}\f$ and \f$\mathrm{R}\f$ have been computed already.
C> Also the permutations \f$P_L\f$ and \f$P_R\f$ due to the pivotting
C> have been provided.
C> Finally, all left occupation functions are expressed in terms
C> of the same left natural orbitals. Likewise all right occupation
C> functions are expressed in terms of the same right natural orbitals.
C> Hence we only need to transform the occupation functions.
C> 
C> So we compute for the left wavefunction
C> \f{eqnarray*}{
C>    \Psi_L^{common} &=& L P_L \Psi_L 
C> \f}
C> and for the right wavefunction
C> \f{eqnarray*}{
C>    \Psi_R^{common} &=& \Psi_R P_R R
C> \f}
C>
      subroutine wfn1_LR2d(nbf,ne,lp,rp,lm,rm,lof,rof,lofc,rofc,m1)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer ne  !< [Input] The number of electrons in the current
                  !< spin channel
c
      double precision lp(ne,ne) !< [Input] The left permutation
      double precision rp(ne,ne) !< [Input] The right permutation
      double precision lm(ne,ne) !< [Input] The left matrix
      !< \f$\mathrm{L}\f$
      double precision rm(ne,ne) !< [Input] The right matrix
      !< \f$\mathrm{R}\f$
      double precision lof(nbf,ne) !< [Input] The left occupation
      !< functions
      double precision rof(nbf,ne) !< [Input] The right occupation
      !< functions
c
      double precision lofc(nbf,ne) !< [Output] The left occupation
      !< functions in the common representation
      double precision rofc(nbf,ne) !< [Output] The right occupation
      !< functions in the common representation
c
      double precision m1(ne,ne) !< [Scratch] Matrix
c
c     Do the left transformation
c
      call dgemm('n','n',ne,ne,ne,1.0d0,lm,ne,lp,ne,0.0d0,m1,ne)
      call dgemm('n','t',ne,nbf,ne,1.0d0,lof,nbf,m1,ne,0.0d0,lofc,nbf)
c
c     Do the right transformation
c
      call dgemm('n','n',ne,ne,ne,1.0d0,rp,ne,rm,ne,0.0d0,m1,ne)
      call dgemm('n','n',ne,nbf,ne,1.0d0,rof,nbf,m1,ne,0.0d0,rofc,nbf)
c
      end
C>
C> @}
