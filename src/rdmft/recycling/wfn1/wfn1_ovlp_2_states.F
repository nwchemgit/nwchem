C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the overlap matrix of the extended orbitals of 
C> two states
C>
C> This routine assumes there are two states and for each one a wavefunction
C> is given. In this context this subroutine calculates the overlap matrix
C> between the occupied orbitals. This information may then be used to 
C> calculate the overlap between the two states (see wfn1_overlap).
C>
      subroutine wfn1_ovlp_2_states(ne,nbf,of1,no1,of2,no2,ov,m1,m2,ovl)
      implicit none
c
      integer ne  !< [Input] The number of electrons
      integer nbf !< [Input] The number of basis functions
c
      double precision of1(nbf,ne)  !< [Input] The occupation functions of 
      !< wavefunction 1
      double precision no1(nbf,nbf) !< [Input] The natural orbitals of 
      !< wavefunction 1
c
      double precision of2(nbf,ne)  !< [Input] The occupation functions of 
      !< wavefunction 2
      double precision no2(nbf,nbf) !< [Input] The natural orbitals of 
      !< wavefunction 2
c
      double precision ov(nbf,nbf) !< [Input] The overlap
      !< matrix in terms of atomic orbitals
c
      double precision m1(nbf,nbf) !< [Scratch] A matrix
      double precision m2(nbf,nbf) !< [Scratch] A matrix
c
      double precision ovl(ne,ne) !< [Output] The overlap matrix between
      !< state 1 and 2
c
      if (ne.eq.0) return
      call dgemm('N','N',nbf,nbf,nbf,1.0d0,ov,nbf,no2,nbf,0.0d0,m1,nbf)
      call dgemm('T','N',nbf,nbf,nbf,1.0d0,no1,nbf,m1,nbf,0.0d0,m2,nbf)
      call dgemm('N','N',nbf,ne,nbf,1.0d0,m2,nbf,of2,nbf,0.0d0,m1,nbf)
      call dgemm('T','N',ne,ne,nbf,1.0d0,of1,nbf,m1,nbf,0.0d0,ovl,ne)
c
c     call dgemm('T','N',ne,ne,nbf,1.0d0,of1,nbf,of2,nbf,0.0d0,ovl,ne)
c
      end
C>
C> @}
