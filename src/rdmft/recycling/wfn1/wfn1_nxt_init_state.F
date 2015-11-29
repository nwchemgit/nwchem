C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Initialize the state vector given the SCF orbitals
C>
C> Initializes the state vector given the SCF orbitals, which means 
C> setting the following
C>
C> 1. Set the occupation functions
C> 2. Set the natural orbitals to the SCF vectors
C> 3. Set the occupation function Lagrangian to zero
C> 4. Set the natural orbital Lagrangian to zero
C>
C> Initializing the occupation functions is a slightly tricky matter.
C> Because of the quadratic form by which the occupation function
C> coefficients contribute to the density matrix the coefficient
C> itself is a factor in the gradient. This means that zero coefficients
C> imply zero gradients for those coefficients, and therefore those
C> coefficients are trapped. Hence a good starting wavefunction should
C> have only non-zero occupation function coefficients. An extreme
C> approach would be to distribute the electrons evenly across all 
C> natural orbitals, but this seems very far removed from a reasonable
C> state. A better alternative is to use Fermi smearing with a 
C> sufficiently high temperature. The procedure is then
C>
C> 1. choose initial temperature T
C> 2. do while (oa(1,1).lt.occtol.or.oa(nbf,1).gt.1.0d0-occtol)
C>    1. increase the temperature
C>    2. distribute `nea` electrons over the natural orbitals
C> 3. copy the first occupation function `nea` times to fill
C>    `oa(:,1:nea)`
C> 4. fill the virtual occupation functions `oa(:,nea+1:nbf)`
C>    with `1-oa(:,1)`
C> 5. change the signs of occupation function coefficients to generate
C>    a linearly independent set
C> 6. orthonormalize the occupation functions
C>
C> Although the Fermi smearing approach can be done it is rather
C> involved. So as a first test it seems better to do something simple
C> and therefore we'll distribute the electrons evenly across the
C> natural orbitals for now.
C>
      subroutine wfn1_nxt_init_state(g_movecs,nbf,ehfa,ehfb,
     &                               oa,na,loa,lna,
     &                               ob,nb,lob,lnb)
      implicit none
c
cinclude "global.fh"
c
      integer g_movecs(2)           !< [Input] GA handles with the
                                    !< SCF vectors
      integer nbf                   !< [Input] The number of basis
                                    !< functions
      double precision ehfa(nbf)    !< [Input] Hartree-Fock alpha 
                                    !< orbital energies
      double precision ehfb(nbf)    !< [Input] Hartree-Fock beta 
                                    !< orbital energies
      double precision oa(nbf,nbf)  !< [Output] The alpha occupation
                                    !< functions
      double precision na(nbf,nbf)  !< [Output] The alpha natural
                                    !< orbitals
      double precision loa(nbf,nbf) !< [Output] The alpha occupation
                                    !< function Lagrangian
      double precision lna(nbf,nbf) !< [Output] The alpha natural
                                    !< orbital Lagrangian
      double precision ob(nbf,nbf)  !< [Output] The beta occupation
                                    !< functions
      double precision nb(nbf,nbf)  !< [Output] The beta natural
                                    !< orbitals
      double precision lob(nbf,nbf) !< [Output] The beta occupation
                                    !< function Lagrangian
      double precision lnb(nbf,nbf) !< [Output] The beta natural
                                    !< orbital Lagrangian
c
      integer ii !< Counter
c
c     Code
c
c     1. Set the occupation functions to unit vectors
c
      call wfn1_nxt_init_occ(nbf,oa,ob)
c
c     2. Set the natural orbitals to SCF vectors
c
      call ga_get(g_movecs(1),1,nbf,1,nbf,na,nbf)
      call ga_get(g_movecs(2),1,nbf,1,nbf,nb,nbf)
c
c     3. Set the occupation function Lagrangian to zero
c
c     call dfill(nbf*nbf,0.0d0,loa,1)
c     call dfill(nbf*nbf,0.0d0,lob,1)
c
c     3. Set the occupation function penalty coefficients to one
c
      call dfill(nbf*nbf,1.0d0,loa,1)
      call dfill(nbf*nbf,1.0d0,lob,1)
c
c     4. Set the natural orbital Lagrangian to zero
c
c     call dfill(nbf*nbf,0.0d0,lna,1)
c     call dfill(nbf*nbf,0.0d0,lnb,1)
c
c     4. Set the natural orbital penalty coefficients to one
c
      call dfill(nbf*nbf,1.0d0,lna,1)
      call dfill(nbf*nbf,1.0d0,lnb,1)
c
      end
 
C>
C> @}
