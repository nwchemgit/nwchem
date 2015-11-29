C> \ingroup wfn1
C> @{
C>
C> \brief Go and compute the energy in wfn1 formalism
C>
C> Go and actually compute the total energy within the wfn1
C> formalism. It is assumed that all the memory has been arranged in
C> the calling routine.
C>
      subroutine wfn1_energy_doit(rtdb,geom,basis,tol,nbf,nst,nea,neb,
     +                            h1,eri,erix,ov,
     +                            r_noa,r_nob,r_ofa,r_ofb,
     +                            r_l1a,r_l1b,r_l2a,r_l2b,r_l3a,r_l3b,
     +                            En,r_dnoa,r_dnob,r_dofa,r_dofb,
     +                            r_dl1a,r_dl1b,r_dl2a,r_dl2b,
     +                            r_dl3a,r_dl3b,l_acta,l_actb,tol_act)
      implicit none
c
#include "stdio.fh"
c
      integer rtdb     !< [Input] The runtime database handle
      integer geom     !< [Input] The geometry handle
      integer basis    !< [Input] The basis set handle
      integer nbf      !< [Input] The number of basis functions
      integer nst      !< [Input] The maximum number of states
      integer nea(nst) !< [Input] The number of alpha electrons
      integer neb(nst) !< [Input] The number of beta electrons
c
      double precision h1(nbf,nbf) !< [Scratch] The 1-electron 
                                   !< Hamiltonian
      double precision ov(nbf,nbf) !< [Scratch] The overlap integrals
      double precision eri(nbf,nbf,nbf,nbf) !< [Scratch] The 2-electron 
                                            !< integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Scratch] The 2-electron
                                             !< integrals
c
      double precision tol                !< Convergence tolerance
      double precision r_noa(nbf,nbf,nst) !< [In/Output] The alpha
      !< electron natural orbitals
      double precision r_nob(nbf,nbf,nst) !< [In/Output] The beta
      !< electron natural orbitals
      double precision r_ofa(nbf,nbf,nst) !< [In/Output] The alpha
      !< electron occupation functions
      double precision r_ofb(nbf,nbf,nst) !< [In/Output] The beta
      !< electron occupation functions
      double precision r_l1a(nbf,nbf,nst) !< [In/Output] The alpha
      !< electron 1st Lagrangian multipliers
      double precision r_l1b(nbf,nbf,nst) !< [In/Output] The beta
      !< electron 1st Lagrangian multipliers
      double precision r_l2a(nbf,nbf,nst) !< [In/Output] The alpha
      !< electron 2nd Lagrangian multipliers
      double precision r_l2b(nbf,nbf,nst) !< [In/Output] The beta
      !< electron 2nd Lagrangian multipliers
      double precision r_l3a(nbf,nbf,nst) !< [In/Output] The alpha
      !< electron 3rd Lagrangian multipliers
      double precision r_l3b(nbf,nbf,nst) !< [In/Output] The beta
      !< electron 3rd Lagrangian multipliers
c
      double precision En(nst) !< [Output] The total energies
c
      double precision r_dnoa(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron natural orbitals derivatives
      double precision r_dnob(nbf,nbf,nst) !< [Scratch] The beta
      !< electron natural orbitals derivatives
      double precision r_dofa(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron occupation functions derivatives
      double precision r_dofb(nbf,nbf,nst) !< [Scratch] The beta
      !< electron occupation functions derivatives
      double precision r_dl1a(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron 1st Lagrangian multipliers derivatives
      double precision r_dl1b(nbf,nbf,nst) !< [Scratch] The beta
      !< electron 1st Lagrangian multipliers derivatives
      double precision r_dl2a(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron 2nd Lagrangian multipliers derivatives
      double precision r_dl2b(nbf,nbf,nst) !< [Scratch] The beta
      !< electron 2nd Lagrangian multipliers derivatives
      double precision r_dl3a(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron 3rd Lagrangian multipliers derivatives
      double precision r_dl3b(nbf,nbf,nst) !< [Scratch] The beta
      !< electron 3rd Lagrangian multipliers derivatives
      logical l_acta(nbf,nst) !< [Scratch] Which \f$\alpha\f$
      !< natural orbitals are active
      logical l_actb(nbf,nst) !< [Scratch] Which \f$\beta\f$
      !< natural orbitals are active
c
      double precision tol_act !< [Input] The threshold for active
      !< natural orbitals
c
c     Local variables
c
      integer iteration      !< The iteration counter
      double precision dnorm !< The gradient norm
      double precision damp  !< The damping factor, has to be negative
                             !< for an energy minimization
c
c     Functions
c
      double precision wfn1_norm
      external         wfn1_norm
c
c     For now we use a steepest decent implementation
c
      call dfill(nbf*nbf,0.0d0,h1,1)
      call dfill(nbf*nbf,0.0d0,ov,1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,eri,1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,erix,1)
cDEBUG
c     write(*,*)'1-electron Hamiltonian'
c     call hess_hssout(h1,nbf,nbf,nbf)
      call errquit("wfn1_energy_doit: code disabled",0,0)
cDEBUG
c
      call int_init(rtdb,1,basis)
      call wfn1_1e_tv(basis,nbf,h1)
      call wfn1_1e_s(basis,nbf,ov)
      call wfn1_2e_rep(basis,nbf,eri,erix)
      call int_terminate
cDEBUG
      write(*,*)'1-electron Hamiltonian'
      call hess_hssout(h1,nbf,nbf,nbf)
      write(*,*)'overlap matrix'
      call hess_hssout(ov,nbf,nbf,nbf)
cDEBUG
c
      iteration = 0
      damp = -0.005d0
      write(luout,'(1x,"iter",3x,"energy",6x,"norm",6x,"damping")')
      call wfn1_engrad(.true.,nbf,1,nea,neb,h1,eri,erix,ov,r_noa,r_nob,
     +                 r_ofa,r_ofb,r_l1a,r_l1b,r_l2a,r_l2b,r_l3a,r_l3b,
     +                 En,r_dnoa,r_dnob,r_dofa,r_dofb,r_dl1a,r_dl1b,
     +                 r_dl2a,r_dl2b,r_dl3a,r_dl3b,l_acta,l_actb,
     +                 tol_act)
      dnorm = wfn1_norm(nbf,r_dnoa,r_dnob,r_dofa,r_dofb,r_dl1a,r_dl1b,
     +                  r_dl2a,r_dl2b,r_dl3a,r_dl3b)
      do while (dnorm.gt.tol)
cDEBUG
        write(*,*)
        write(*,*)'natural orbitals'
        call hess_hssout(r_noa,nbf,nbf,nbf)
        write(*,*)'first Lagrangian'
        call hess_hssout(r_l1a,nbf,nbf,nbf)
        write(*,*)'occupation functions'
        call hess_hssout(r_ofa,nbf,nbf,nbf)
        write(*,*)'third Lagrangian'
        call hess_hssout(r_l3a,nbf,nbf,nbf)
        write(*,*)
cDEBUG
        iteration = iteration + 1
        write(luout,'(1x,i5,f22.8,f22.8,f8.3)')iteration,en(1),dnorm,
     +                                         damp
        call wfn1_update(nbf,damp,r_noa,r_nob,r_ofa,r_ofb,r_l1a,r_l1b,
     +                   r_l2a,r_l2b,r_l3a,r_l3b,r_dnoa,r_dnob,
     +                   r_dofa,r_dofb,r_dl1a,r_dl1b,r_dl2a,r_dl2b,
     +                   r_dl3a,r_dl3b)
        call wfn1_engrad(.true.,nbf,1,nea,neb,h1,eri,erix,ov,
     +                   r_noa,r_nob,
     +                   r_ofa,r_ofb,r_l1a,r_l1b,r_l2a,r_l2b,
     +                   r_l3a,r_l3b,En,r_dnoa,r_dnob,r_dofa,r_dofb,
     +                   r_dl1a,r_dl1b,r_dl2a,r_dl2b,r_dl3a,r_dl3b,
     +                   l_acta,l_actb,tol_act)
        dnorm = wfn1_norm(nbf,r_dnoa,r_dnob,r_dofa,r_dofb,r_dl1a,r_dl1b,
     +                    r_dl2a,r_dl2b,r_dl3a,r_dl3b)
      enddo
      write(luout,'(1x,i5,f22.8,f22.8,f8.3)')iteration,en(1),dnorm,
     +                                       damp
      
      end
C> @}
