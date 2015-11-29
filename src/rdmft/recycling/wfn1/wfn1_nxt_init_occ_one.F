C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Initialize the occupation functions
C>
C> The occupation functions are initialized so that the electrons
C> are distributed as in the Hartree-Fock wavefunction. I.e. the
C> occupation function vectors form the identity matrix.
C> 
      subroutine wfn1_nxt_init_occ_one(nbf,oa,ob)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision oa(nbf,nbf) !< [Output] The alpha occupation
                                   !< functions
      double precision ob(nbf,nbf) !< [Output] The beta occupation
                                   !< functions
c
      integer k              !< Counter over frequencies
      integer i              !< Counter over coefficients
      integer p              !< Counter over occupation functions
      integer q              !< Counter over occupation functions
      double precision delta !< The width of an interval
      double precision x     !< The position associated with an
                             !< occupation function coefficient
      double precision dd    !< The overlap
      double precision pi    !< Pi
c
      double precision ddot
      external         ddot
c
c     Fill the occupation functions with unit coefficients
c
      call dfill(nbf*nbf,0.0d0,oa,1)
      call dfill(nbf*nbf,0.0d0,ob,1)
c
c     Set the diagonal elements to one
c
      call dfill(nbf,1.0d0,oa,nbf+1)
      call dfill(nbf,1.0d0,ob,nbf+1)
c
      end
C>
C> @}
