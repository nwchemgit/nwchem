C> \ingroup wfn1
C> @{
C>
C> \brief Updates the current wavefunction with a scaled term
C>
C> Update the current wavefunction by adding a scaled term. Typically
C> the term being added is somehow related to the gradient.
C>
      subroutine wfn1_update(nbf,fac,r_noa,r_nob,r_ofa,r_ofb,
     +                       r_l1a,r_l1b,r_l2a,r_l2b,r_l3a,r_l3b,
     +                       r_dnoa,r_dnob,r_dofa,r_dofb,r_dl1a,r_dl1b,
     +                       r_dl2a,r_dl2b,r_dl3a,r_dl3b)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
c
      integer nbf      !< [Input] The number of basis functions
c
      double precision fac            !< [Input] The scale factor
c
      double precision r_noa(nbf,nbf) !< [Input] The 
      !< \f$\alpha\f$-spin natural orbitals. This is a 2D array with the
      !< dimensions:
      !< - 1. The basis functions
      !< - 2. The natural orbitals
      double precision r_nob(nbf,nbf) !< [Input] The 
      !< \f$\beta\f$-spin natural orbitals. This is a 2D array with the
      !< same dimensions as `r_noa`.
      double precision r_ofa(nbf,nbf) !< [Input] The
      !< \f$\alpha\f$-spin occupation functions. This is a 2D array with
      !< the dimension:
      !< - 1. The natural orbitals
      !< - 2. The occupation functions
      double precision r_ofb(nbf,nbf) !< [Input] The 
      !< \f$\beta\f$-spin occupation functions. This is a 2D array with
      !< the same dimensions as `r_ofa`.
      double precision r_l1a(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the orthonormality of the \f$\alpha\f$-spin
      !< natural orbitals. This is a 2D array with dimensions:
      !< - 1. The natural orbitals
      !< - 2. The natural orbitals
      double precision r_l1b(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the orthonormality of the \f$\beta\f$-spin
      !< natural orbitals. This is a 2D array with the same dimensions
      !< as `r_l1a`.
      double precision r_l2a(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the \f$\alpha\f$-spin natural orbitals being
      !< eigenfunctions of the density matrix. This is
      !< a 2D GA with the dimensions:
      !< - 1. The natural orbitals
      !< - 2. The natural orbitals
      double precision r_l2b(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the \f$\beta\f$-spin natural orbitals being
      !< eigenfunctions of the density matrix. This is a 2D array with
      !< the same dimensions as `r_l2a`.
      double precision r_l3a(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the orthonormality of the \f$\alpha\f$-spin
      !< occupation functions. This is a 2D array with dimensions:
      !< - 1. The occupation functions
      !< - 2. The occupation functions
      double precision r_l3b(nbf,nbf) !< [Input] The Lagrangian
      !< multipliers for the orthonormality of the \f$\beta\f$-spin
      !< natural orbitals. This is a 2D array with the same dimensions
      !< as `r_l1a`.
c
      double precision r_dnoa(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_noa`. The array dimensions are defined
      !< analogously to `r_noa`.
      double precision r_dnob(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_nob`. The array dimensions are defined
      !< analogously to `r_nob`.
      double precision r_dofa(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_ofa`. The array dimensions are defined
      !< analogously to `r_ofa`.
      double precision r_dofb(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_ofb`. The array dimensions are defined
      !< analogously to `r_ofb`.
      double precision r_dl1a(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l1a`. The array dimensions are defined
      !< analogously to `r_l1a`.
      double precision r_dl1b(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l1b`. The array dimensions are defined
      !< analogously to `r_l1b`.
      double precision r_dl2a(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l2a`. The array dimensions are defined
      !< analogously to `r_l2a`.
      double precision r_dl2b(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l2b`. The array dimensions are defined
      !< analogously to `r_l2b`.
      double precision r_dl3a(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l3a`. The array dimensions are defined
      !< analogously to `r_l3a`.
      double precision r_dl3b(nbf,nbf) !< [Output] The energy
      !< gradient wrt `r_l3b`. The array dimensions are defined
      !< analogously to `r_l3b`.
C
C     Local variables
C
      integer len
c
c     Code
c
      len = nbf*nbf
      call daxpy(len,fac,r_dnoa,1,r_noa,1)
      call daxpy(len,fac,r_dnob,1,r_nob,1)
      call daxpy(len,fac,r_dofa,1,r_ofa,1)
      call daxpy(len,fac,r_dofb,1,r_ofb,1)
      call daxpy(len,fac,r_dl1a,1,r_l1a,1)
      call daxpy(len,fac,r_dl1b,1,r_l1b,1)
      call daxpy(len,fac,r_dl2a,1,r_l2a,1)
      call daxpy(len,fac,r_dl2b,1,r_l2b,1)
      call daxpy(len,fac,r_dl3a,1,r_l3a,1)
      call daxpy(len,fac,r_dl3b,1,r_l3b,1)
c
      end
C> @}

