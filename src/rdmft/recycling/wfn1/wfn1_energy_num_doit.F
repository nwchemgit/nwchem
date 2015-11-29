C> \ingroup wfn1
C> @{
C>
C> \brief Go and numerically optimize the energy in wfn1 formalism
C>
C> Go and actually compute the total energy within the wfn1
C> formalism. It is assumed that all the memory has been arranged in
C> the calling routine.
C>
C> Because the gradient of the energy is somehow non-trivial to
C> formulate correctly this routine uses numerical differentation
C> to optimize the energy. Obviously this is rather inefficient but
C> it will
C> - provide an opportunity to assess the quality of the energy
C>   expression,
C> - provide a reference implementation that can be used to test
C>   analytical derivatives.
C> The routine uses the 3-point rule for the gradient:
C> \f{eqnarray*}{
C>   \frac{\partial f}{\partial x} 
C>   &=& \frac{f(x+h)-f(x-h)}{2h} + O(h^2)
C> \f}
C> as well as the 3-point rule for the second derivative:
C> \f{eqnarray*}{
C>   \frac{\partial^2 f}{\partial x^2} 
C>   &=& \frac{f(x+h)-2f(x)+f(x-h)}{h^2} + O(h)
C> \f}
C> This information allows a reasonable estimate of the displacement
C> \f$d\f$ such that
C> \f{eqnarray*}{
C>   \left.\frac{\partial f(y)}{\partial y}\right|_{y=x+d} &=& 0
C> \f}
C> The question of picking the optimal \f$h\f$ involves the trade-off
C> between two kinds of errors:
C> - the truncation error is smaller for smaller \f$h\f$
C> - the numerical round-off error is smaller for larger \f$h\f$
C> On balance \f$h_{\mathrm{opt}} = max(|x|,1)\sqrt{\epsilon}\f$,
C> where \f$\epsilon\f$ is the machine precision.
C>
C> The basic algorithm behind this routine is to go through the four
C> sets of rotations:
C> - \f$\alpha\f$-occupation functions
C> - \f$\beta\f$-occupation functions
C> - \f$\alpha\f$-natural orbitals functions
C> - \f$\beta\f$-natural orbitals functions
C> and for each set loop over all function pairs evaluating for each
C> what a rotation by a small angle would do to the energy. These 
C> iterations are repeated until the gradient gets really small.
C>
C> If we do the \f$\alpha\f$ rotations first and then the \f$\beta\f$
C> rotations spin symmetry breaking occurs. The reason is that the 
C> optimization of the beta orbitals takes place in a different 
C> environment than the optimization of the alpha orbitals. Hence we
C> need a different algorithm. The new algorithms is as follows:
C> - For the occupation functions
C>   - for all orbital pairs work out the rotations of alpha functions
C>     (but do not apply the rotation)
C>   - for all orbital pairs work out the rotations of beta functions
C>     (but do not apply the rotation)
C>   - having constructed the alpha and beta rotation matrices
C>     apply both rotations to the wavefunction
C> - For the natural orbitals
C>   - for all orbital pairs work out the rotations of alpha functions
C>     (but do not apply the rotation)
C>   - for all orbital pairs work out the rotations of beta functions
C>     (but do not apply the rotation)
C>   - having constructed the alpha and beta rotation matrices
C>     apply both rotations to the wavefunction
C> In addition we need to start the calculation off in a smart way
C> as the Hartree-Fock solution is typically meta-stable. We do this
C> by starting with small rotations and increase the rotation angle
C> until the entropy has a favorable effect. In this phase we do not
C> calculate the gradient, but pick the lowest energy from the three
C> points we calculate and use the corresponding rotation. Using the
C> gradient could once again make the step too small and we get no 
C> energy lowering. Finally to get a sensible distribution of the
C> electrons over the natural orbitals we want to reduce the angle
C> of rotation between orbitals with very different energies. For this
C> purpose we can use the expression
C> \f{eqnarray*}{
C>   A &=& \frac{A0}{1+abs(\epsilon_s-\epsilon_t)}
C> \f}
C> where \f$A0\f$ is the base rotation, \f$A\f$ is the actual rotation,
C> and \f$\epsilon_s\f$ and \f$\epsilon_t\f$ are orbital energies.
C>
C> To deal with ground states and excited states we need to enforce
C> some constraints. I.e. an excited state is a minimal energy state
C> subject to the constraint that it is orthogonal to the ground state.
C> For this purpose we can define the total energy with a Lagrange
C> multiplier as:
C> \f{eqnarray*}{
C>   E &=& H + \lambda S
C> \f}
C> from which we have by differentation
C> \f{eqnarray*}{
C>   \frac{\partial E}{\partial C_{ij}} &=&
C>      \frac{\partial H}{\partial C_{ij}} + 
C>      \lambda \frac{\partial S}{\partial C_{ij}}  \\\\
C>   \frac{\partial E}{\partial \lambda} &=&
C>      S \\\\
C> \f}
C> The 2-norm of this expression amounts to
C> \f{eqnarray*}{
C>   ||dE|| &=&
C>      \sqrt{\sum_{ij} \left[\left(\frac{\partial H}{\partial C_{ij}}
C>                            \right)^2
C>      +2\lambda\frac{\partial H}{\partial C_{ij}}
C>               \frac{\partial S}{\partial C_{ij}}
C>      +\lambda^2\left(\frac{\partial S}{\partial C_{ij}}\right)^2
C>       \right]
C>      + S^2}
C> \f}
C> Next we choose \f$\lambda\f$ so as to minimize this norm
C> \f{eqnarray*}{
C>   \begin{array}{cc}
C>     \textrm{min}&||dE(\lambda)||^2 \\\\
C>     \lambda\in R&
C>   \end{array}
C> \f}
C> where, for convenience, we chose to minimize the square of the norm
C> as that leads to the same results but eliminates the square root.
C> Differentiation gives
C> \f{eqnarray*}{
C>    0 &=& 2\sum_{ij}\frac{\partial H}{\partial C_{ij}}
C>                    \frac{\partial S}{\partial C_{ij}} +
C>          2\lambda\sum_{ij}\left(
C>                    \frac{\partial H}{\partial C_{ij}}\right)^2 \\\\
C>    \lambda &=& -\frac{\sum_{ij}\frac{\partial H}{\partial C_{ij}}
C>                          \frac{\partial S}{\partial C_{ij}}}
C>                 {\sum_{ij}\left(
C>                    \frac{\partial S}{\partial C_{ij}}\right)^2 }
C> \f}
C> With this the general algorithm becomes for a given state:
C> * Compute \f$\frac{\partial H}{\partial C_{ij}}\f$ and
C>   \f$\frac{\partial S}{\partial C_{ij}}\f$
C> * Compute \f$\lambda\f$
C> * Compute \f$\frac{\partial E}{\partial C_{ij}}\f$
C> * Compute \f$||dE||\f$
C> * Compute \f$\frac{\partial^2 E}{\partial C^2_{ij}}\f$
C> * Precondition the gradient using the Hessian
C> * Compute the rotation matrices
C> * Apply the rotations
C> This algorithm is repeated until all states are converged.
C> 
      subroutine wfn1_energy_num_doit(rtdb,geom,basis,tol,uwfn1,nbf,nst,
     +         nea,neb,h1,eri,erix,ov,
     +         r_noa,r_nob,r_ofa,r_ofb,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +         d_noa,d_nob,d_ofa,d_ofb,d_noa2,d_nob2,d_ofa2,d_ofb2,
     +         rnoa,rnob,rofa,rofb,sm1,
     +         En0,On0,r_u,
     +         r_ohpa,r_ohpb,r_ohma,r_ohmb,r_olpa,r_olpb,r_olma,r_olmb,
     +         r_nhpa,r_nhpb,r_nhma,r_nhmb,r_nlpa,r_nlpb,r_nlma,r_nlmb,
     +         r_doha,r_dohb,r_dola,r_dolb,r_dnha,r_dnhb,r_dnla,r_dnlb,
     +         ovnoa,ovnob,ovofa,ovofb,m1,m2,m3,
     +         l_acta,l_actb,Ehfa,Ehfb,Ewfa,Ewfb,enta,entb,
     +         temperature,power,
     +         nperma,npermb,perma,permb,signa,signb,ovla,ovlb,
     +         tol_act)
      implicit none
c
#include "stdio.fh"
#include "mafdecls.fh"
c
      integer rtdb     !< [Input] The runtime database handle
      integer geom     !< [Input] The geometry handle
      integer basis    !< [Input] The basis set handle
      logical uwfn1    !< [Input] If .true. do unrestricted WFN1
      integer nbf      !< [Input] The number of basis functions
      integer nst      !< [Input] The number of states
      integer nea      !< [Input] The number of alpha electrons
      integer neb      !< [Input] The number of beta electrons
      integer nperma   !< [Input] The number of alpha electron
                       !< permutations
      integer npermb   !< [Input] The number of beta electron
                       !< permutations
c
      double precision temperature !< [Input] The electron temperature
      double precision power       !< [Input] The matrix power
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
      double precision r_noa2(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron natural orbitals
      double precision r_nob2(nbf,nbf,nst) !< [Scratch] The beta
      !< electron natural orbitals
      double precision r_ofa2(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron occupation functions
      double precision r_ofb2(nbf,nbf,nst) !< [Scratch] The beta
      !< electron occupation functions
c
      double precision d_noa(nbf,nbf) !< [Scratch] The alpha
      !< electron natural orbitals gradient
      double precision d_nob(nbf,nbf) !< [Scratch] The beta
      !< electron natural orbitals gradient
      double precision d_ofa(nbf,nbf) !< [Scratch] The alpha
      !< electron occupation functions gradient
      double precision d_ofb(nbf,nbf) !< [Scratch] The beta
      !< electron occupation functions gradient
      double precision d_noa2(nbf,nbf) !< [Scratch] The alpha
      !< electron natural orbitals gradient
      double precision d_nob2(nbf,nbf) !< [Scratch] The beta
      !< electron natural orbitals gradient
      double precision d_ofa2(nbf,nbf) !< [Scratch] The alpha
      !< electron occupation functions gradient
      double precision d_ofb2(nbf,nbf) !< [Scratch] The beta
      !< electron occupation functions gradient
c
      double precision rnoa(nbf,nbf) !< [Scratch] The alpha natural
      !< orbital rotation matrix
      double precision rnob(nbf,nbf) !< [Scratch] The beta natural
      !< orbital rotation matrix
      double precision rofa(nbf,nbf) !< [Scratch] The alpha occupation
      !< function rotation matrix
      double precision rofb(nbf,nbf) !< [Scratch] The beta occupation
      !< function rotation matrix
      double precision sm1(nbf,nbf) !< [Scratch] Matrix
c
      double precision En0(nst) !< [Output] The total energies
      double precision On0(nst) !< [Output] The overlaps
      double precision On1(nst) !< The old overlaps
      double precision On2(nst) !< The old overlaps
      double precision En2(nst) !< The old total energies
c
      double precision r_u(nbf,nbf) !< [Scratch] The rotation matrix
      double precision r_ohpa(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_ohpb(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_ohma(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_ohmb(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_olpa(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_olpb(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_olma(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_olmb(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_nhpa(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_nhpb(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_nhma(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_nhmb(nbf,nbf) !< [Scratch] The perturbed
                                       !< Hamiltonian matrix
      double precision r_nlpa(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_nlpb(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_nlma(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
      double precision r_nlmb(nbf,nbf,nst) !< [Scratch] The perturbed
                                           !< overlap matrix
c
      double precision r_doha(nbf,nbf) !< [Scratch] The gradient of the
      !< Hamiltonian wrt occupation functions
      double precision r_dohb(nbf,nbf) !< [Scratch] The gradient of the
      !< Hamiltonian wrt occupation functions
      double precision r_dola(nbf,nbf,nst) !< [Scratch] The gradient of
      !< the overlap wrt occupation functions
      double precision r_dolb(nbf,nbf,nst) !< [Scratch] The gradient of
      !< the overlap wrt occupation functions
      double precision r_dnha(nbf,nbf) !< [Scratch] The gradient of the
      !< Hamiltonian wrt natural orbitals
      double precision r_dnhb(nbf,nbf) !< [Scratch] The gradient of the
      !< Hamiltonian wrt natural orbitals
      double precision r_dnla(nbf,nbf,nst) !< [Scratch] The gradient of
      !< the overlap wrt natural orbitals
      double precision r_dnlb(nbf,nbf,nst) !< [Scratch] The gradient of
      !< the overlap wrt natural orbitals
c
      double precision ovnoa(nbf,nbf) !< [Scratch] The alpha natural
      !< orbital overlap matrix
      double precision ovnob(nbf,nbf) !< [Scratch] The beta natural
      !< orbital overlap matrix
      double precision ovofa(nbf,nbf) !< [Scratch] The alpha occupation
      !< function overlap matrix
      double precision ovofb(nbf,nbf) !< [Scratch] The beta occupation
      !< function overlap matrix
c
      double precision m1(nbf,nbf) !< [Scratch] 
      double precision m2(nbf,nbf) !< [Scratch] 
      double precision m3(nbf,nbf) !< [Scratch] 
c
      logical l_acta(nbf,nst) !< [Scratch] Which \f$\alpha\f$
      !< natural orbitals are active
      logical l_actb(nbf,nst) !< [Scratch] Which \f$\beta\f$
      !< natural orbitals are active
c
      double precision Ehfa(nbf) !< [Input] The Hartree-Fock alpha 
      !< orbital energies
      double precision Ehfb(nbf) !< [Input] The Hartree-Fock beta
      !< orbital energies
      double precision Ewfa(nbf) !< [Scratch] The wfn1 alpha 
      !< orbital energies
      double precision Ewfb(nbf) !< [Scratch] The wfn1 beta
      !< orbital energies
      double precision enta(nbf) !< [Scratch] The alpha orbital
      !< entropies
      double precision entb(nbf) !< [Scratch] The beta orbital entropies
c
      double precision dnh, dnl
c
      integer perma(nea,nperma) !< [Input] The alpha electron
                                !< permutations
      integer permb(neb,npermb) !< [Input] The beta electron
                                !< permutations
c
      double precision signa(nperma) !< [Input] The sign of the alpha
                                     !< permutations
      double precision signb(npermb) !< [Input] The sign of the beta
                                     !< permutations
      double precision ovla(nea,nea) !< [Scratch] Overlap of alpha
                                     !< orbitals
      double precision ovlb(neb,neb) !< [Scratch] Overlap of beta
                                     !< orbitals
c
      double precision tol_act !< [Input] The threshold for active
      !< natural orbitals
c
c     Local variables
c
      integer iteration           !< The iteration counter
      integer it_ent              !< The iteration counter for entropy
      double precision dnorm      !< The gradient norm
      double precision dnrmg      !< The gradient norm
      double precision dnrmm      !< The gradient,hessian overlap
      double precision dnrm       !< The norm
      double precision d2nrm      !< The norm of the Hessian
      double precision dnrmh      !< The norm
      double precision dnrml      !< The norm
      double precision dnorm1(nst)!< The gradient norm for each state
      double precision dlambda(nst) !< Aggregate Lagrange multiplier
                                    !< for each state
      double precision dlambda_i(nst) !< Lagrange multiplier
                                      !< for each state
      double precision dlambda_k(nst) !< Lagrange multiplier
                                      !< for line search
      double precision damp(nst)  !< The damping factor
                                  !< for an energy minimization
      double precision dampl      !< The damping factor
                                  !< for an overlap minimization
      double precision h          !< The step size
      double precision ho         !< The step size for
                                  !< occupation functions
      double precision sum        !< An accumulator
      double precision d          !< The displacement
      double precision rotup(2,2) !< Rotation over positive angle
      double precision rotdn(2,2) !< Rotation over negative angle
      double precision rotupo(2,2) !< Rotation over positive
                                   !< angle occupation
      double precision rotdno(2,2) !< Rotation over negative
                                   !< angle occupation
      double precision Em1(nst)   !< The total energies
      double precision Ep1(nst)   !< The total energies
      double precision Om1(nst)   !< The overlaps
      double precision Op1(nst)   !< The overlaps
      double precision dS         !< The overlap gradient
      double precision dE         !< The energy gradient
      double precision d2E        !< The energy Hessian
      double precision serr       !< The error in S^2
c     double precision oerr(nst)  !< The error in the orbital ordering
      double precision sfac       !< The scale factor for the spin
                                  !< penality function
      double precision ofac       !< The scale factor for the orbital
                                  !< order penality function
      double precision pi         !< The number Pi.
      double precision step
      integer is                  !< Counter over occupation functions
      integer it                  !< Counter over occupation functions
      integer iu                  !< Counter over occupation functions
      integer im                  !< Counter over natural orbitals
      integer in                  !< Counter over natural orbitals
      integer ist                 !< Counter over states
      integer jst                 !< Counter over states
      integer mst                 !< The current state being optimized
      logical oconv(nst)          !< .True. if a state is converged
      logical oprint              !< .True. if more details should be
                                  !< printed
cDEBUG
      logical stat
cDEBUG
c
c     Parameters
c
      integer max_it  !< The maximum number of iterations in pre-opt
                      !< loop.
      parameter (max_it = 20)
c
c     Functions
c
      double precision wfn1_norm
      external         wfn1_norm
c
c     Optimal h is the square root of the machine precision
c
      pi = acos(-1.0d0)
      sfac = 1.0d0
      ofac = 1.0d0
      tol_act = -1.0d0
      h = 1.0d-8
cDEBUG
      write(*,*)'convergence tolerance = ',tol
      h = 1.0d-6
      h = 1.0d-4
      write(*,*)"*** HVD A: ma_verify"
      call util_flush(6)
      stat = ma_verify_allocator_stuff()
cDEBUG
      rotup(1,1) =  cos(h)
      rotup(2,1) = -sin(h)
      rotup(1,2) =  sin(h)
      rotup(2,2) =  cos(h)
      rotdn(1,1) =  cos(-h)
      rotdn(2,1) = -sin(-h)
      rotdn(1,2) =  sin(-h)
      rotdn(2,2) =  cos(-h)
c
c     For now we use a steepest decent implementation
c
      call dfill(nbf*nbf,0.0d0,h1,1)
      call dfill(nbf*nbf,0.0d0,ov,1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,eri,1)
      call dfill(nbf*nbf*nbf*nbf,0.0d0,erix,1)
c
      call int_init(rtdb,1,basis)
      call wfn1_1e_tv(basis,nbf,h1)
      call wfn1_1e_s(basis,nbf,ov)
      call wfn1_2e_rep(basis,nbf,eri,erix)
      call int_terminate
cDEBUG
      write(*,*)"*** HVD B: ma_verify"
      call util_flush(6)
      stat = ma_verify_allocator_stuff()
cDEBUG
c
      dnorm = 10.0d0
      serr  = 10.0d0
      iteration = 0
      dampl = 1.0000d0
      mst = 1
      oprint = .false.
      do ist = 1, nst
        damp(ist) = 0.250d0*(0.25d0**(ist-1))
        dnorm1(ist) = 10.0d0
        dlambda(ist) = 0.0d0
        dlambda_i(ist) = 0.0d0
        oconv(ist)  = .false.
        On1(ist) = 0.0d0
        En0(ist) = 0.0d0
        En2(ist) = 0.0d0
        On0(ist) = 0.0d0
      enddo
cDEBUG
      write(*,*)'*** uwfn1 = ',uwfn1
cDEBUG
      write(luout,'(3x,"tolerance = ",f22.8)')tol
      write(luout,'(3x,"iter",3x,"energy",6x,"norm",6x,"damping")')
cDEBUG
c     write(*,*)"*** HVD C: ma_verify"
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
      do while (dnorm.gt.tol)
c
c       Now move on to actually optimize the wavefunctions:
c       - Minimize energy of states in increasing state number order,
c       i.e. state 1 first, then state 2, etc.
c       - Higher states have penalty function added for overlap with
c       lower states.
c
        iteration = iteration + 1
cDEBUG
        oprint = .true.
        if (iteration.eq.1.or.mod(iteration,1000).eq.0.or.oprint) then
          oprint = .true.
          do ist = 1, nst
            write(*,*)
            write(*,*)'alpha natural orbitals',ist
            call hess_hssout(r_noa(1,1,ist),nbf,nbf,nbf)
            write(*,*)
            write(*,*)'beta natural orbitals',ist
            call hess_hssout(r_nob(1,1,ist),nbf,nbf,nbf)
          enddo
          do ist = 1, nst
            write(*,*)
            write(*,*)'alpha occupation functions',ist
            call hess_hssout(r_ofa(1,1,ist),nbf,nbf,nbf)
            write(*,*)
            write(*,*)'beta occupation functions',ist
            call hess_hssout(r_ofb(1,1,ist),nbf,nbf,nbf)
            write(*,*)
          enddo
          call dfill(nbf*nbf,0.0d0,m2,1)
          do ist = 1, nea
            m2(ist,1) = 1.0d0
          enddo
          do ist = 1, neb
            m2(ist,2) = 1.0d0
          enddo
          do ist = 1, nst
            write(*,*)
            write(*,*)'occupation numbers',ist
            call wfn1_print_occ(nbf,nea,neb,
     +           r_ofa(1,1,ist),r_ofb(1,1,ist),m2,m2(1,2),dnrm,m1)
          enddo
        endif
cDEBUG
        if (iteration.eq.1.or.mod(iteration,100).eq.0.or.oprint) then
          do ist = 1, nst
            write(luout,'(1x,i7," iter ",i3,f22.8,f22.8,f22.16,f22.8,
     +                    2f10.5,l4)')
     +            iteration,ist,En0(ist),dnorm1(ist),On0(ist),
     +            dlambda(ist),dlambda_i(ist),dampl,oconv(ist)
c           dnorm1(ist) = 0.0d0
            On0(ist) = 0.0d0
            oconv(ist)  = .true.
          enddo
          write(luout,*)
          dnorm = 0.0d0
        else
          do ist = 1, nst
c           dnorm1(ist) = 0.0d0
            on0(ist) = 0.0d0
            oconv(ist)  = .true.
          enddo
          dnorm = 0.0d0
        endif
c
c***
c       do ist = 1, nst
c
c         Initialize Lagrangian or the gradient of the overlap
c
          call dfill(nbf*nbf,0.0d0,r_ohpa,1)
          call dfill(nbf*nbf,0.0d0,r_ohpb,1)
          call dfill(nbf*nbf,0.0d0,r_ohma,1)
          call dfill(nbf*nbf,0.0d0,r_ohmb,1)
          call dfill(nbf*nbf,0.0d0,r_olpa,1)
          call dfill(nbf*nbf,0.0d0,r_olpb,1)
          call dfill(nbf*nbf,0.0d0,r_olma,1)
          call dfill(nbf*nbf,0.0d0,r_olmb,1)
c
          call dfill(nbf*nbf,0.0d0,r_nhpa,1)
          call dfill(nbf*nbf,0.0d0,r_nhpb,1)
          call dfill(nbf*nbf,0.0d0,r_nhma,1)
          call dfill(nbf*nbf,0.0d0,r_nhmb,1)
          call dfill(nbf*nbf,0.0d0,r_nlpa,1)
          call dfill(nbf*nbf,0.0d0,r_nlpb,1)
          call dfill(nbf*nbf,0.0d0,r_nlma,1)
          call dfill(nbf*nbf,0.0d0,r_nlmb,1)
c
c         Rotate occupation functions
c
          if (mst.eq.1) then
          call dcopy(nbf*nbf*nst,r_ofb,1,r_ofb2,1)
          call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
          call dcopy(nbf*nbf*nst,r_nob,1,r_nob2,1)
          do it = 1, nbf
            do is = 1, it-1
c
c             Rotate alpha occupation functions
c
              call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
cDEBUG
c     write(*,*)"*** HVD D: ma_verify"
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         En0(mst),serr,On2,
     +                         (is.eq.1.and.it.eq.2).and.oprint)
cDEBUG
c     write(*,*)"*** HVD E: ma_verify"
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
              if (is.eq.1.and.it.eq.2) then
                dlambda(mst) = 0.0d0
                do jst = 1, mst-1
                  dlambda_i(1) = dlambda_i(1) + dampl*On2(jst)
                enddo
                do jst = 1, mst
                  dlambda_i(jst) = dlambda_i(1)
                enddo
                do jst = 1, mst-1
                  dlambda(mst) = dlambda(mst) + On2(jst)
                enddo
              endif
              oprint = .false.
              On0(mst) = 0.0d0
              do jst = 1, mst-1
                On0(mst) = On0(mst) + abs(On2(jst))
              enddo
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(is,is) = rotup(1,1)
              r_u(it,it) = rotup(2,2)
              r_u(it,is) = rotup(2,1)
              r_u(is,it) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_ofa2(1,1,mst),nbf)
              if (.not.uwfn1) then
                call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_ohpa(is,it),serr,On1,.false.)
              do jst = 1, mst
                r_olpa(is,it,jst) = On1(jst)
              enddo
              r_u(is,is) = rotdn(1,1)
              r_u(it,it) = rotdn(2,2)
              r_u(it,is) = rotdn(2,1)
              r_u(is,it) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_ofa2(1,1,mst),nbf)
              if (.not.uwfn1) then
                call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_ohma(is,it),serr,On1,.false.)
              do jst = 1, mst
                r_olma(is,it,jst) = On1(jst)
              enddo
            enddo
          enddo
c
c         Rotate beta occupation functions
c
          if (uwfn1) then
          call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
          call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
          call dcopy(nbf*nbf*nst,r_nob,1,r_nob2,1)
          do it = 1, nbf
            do is = 1, it-1
c
c             Rotate beta occupation functions
c
              call dcopy(nbf*nbf*nst,r_ofb,1,r_ofb2,1)
cDEBUG
c     write(*,*)"*** HVD F: ma_verify"
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         En0(mst),serr,On2,.false.)
cDEBUG
c     write(*,*)"*** HVD G: ma_verify"
c     call util_flush(6)
c     stat = ma_verify_allocator_stuff()
cDEBUG
              On0(mst) = 0.0d0
              do jst = 1, mst-1
                On0(mst) = On0(mst) + abs(On2(jst))
              enddo
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(is,is) = rotup(1,1)
              r_u(it,it) = rotup(2,2)
              r_u(it,is) = rotup(2,1)
              r_u(is,it) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_ofb2(1,1,mst),nbf)
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_ohpb(is,it),serr,On1,.false.)
              do jst = 1, mst
                r_olpb(is,it,jst) = On1(jst)
              enddo
              r_u(is,is) = rotdn(1,1)
              r_u(it,it) = rotdn(2,2)
              r_u(it,is) = rotdn(2,1)
              r_u(is,it) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_ofb2(1,1,mst),nbf)
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_ohmb(is,it),serr,On1,.false.)
              do jst = 1, mst
                r_olmb(is,it,jst) = On1(jst)
              enddo
            enddo
          enddo
          else
            call dcopy(nbf*nbf,r_ohpa,1,r_ohpb,1)
            call dcopy(nbf*nbf*mst,r_olpa,1,r_olpb,1)
            call dcopy(nbf*nbf,r_ohma,1,r_ohmb,1)
            call dcopy(nbf*nbf*mst,r_olma,1,r_olmb,1)
          endif
          endif ! mst.eq.1
c
c         Rotate natural orbitals
c
          call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
          call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
          call dcopy(nbf*nbf*nst,r_ofb,1,r_ofb2,1)
          call dcopy(nbf*nbf*nst,r_nob,1,r_nob2,1)
          do in = 1, nbf
            do im = 1, in-1
c
c             Rotate alpha natural orbitals
c
              call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
              call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         En0(mst),serr,On2,.false.)
              if (mst.gt.1.and.im.eq.1.and.in.eq.2) then
                dlambda(mst) = 0.0d0
                do jst = 1, mst-1
                  dlambda_i(1) = dlambda_i(1) + dampl*On2(jst)
                enddo
                do jst = 1, mst
                  dlambda_i(jst) = dlambda_i(1)
                enddo
                do jst = 1, mst-1
                  dlambda(mst) = dlambda(mst) + On2(jst)
                enddo
              endif
              oprint = .false.
              On0(mst) = 0.0d0
              do jst = 1, mst-1
                On0(mst) = On0(mst) + abs(On2(jst))
              enddo
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(im,im) = rotup(1,1)
              r_u(in,in) = rotup(2,2)
              r_u(in,im) = rotup(2,1)
              r_u(im,in) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_noa2(1,1,mst),nbf)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_u,nbf,
     +                     r_ofa(1,1,mst),nbf,0.0d0,r_ofa2(1,1,mst),nbf)
              else
                call dcopy(nbf*nbf,r_ofa(1,1,mst),1,r_ofa2(1,1,mst),1)
              endif
              if (.not.uwfn1) then
                call dcopy(nbf*nbf,r_noa2(1,1,mst),1,r_nob2(1,1,mst),1)
                call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_nhpa(im,in),serr,On1,.false.)
              do jst = 1, mst
                r_nlpa(im,in,jst) = On1(jst)
              enddo
              r_u(im,im) = rotdn(1,1)
              r_u(in,in) = rotdn(2,2)
              r_u(in,im) = rotdn(2,1)
              r_u(im,in) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_noa2(1,1,mst),nbf)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_u,nbf,
     +                     r_ofa(1,1,mst),nbf,0.0d0,r_ofa2(1,1,mst),nbf)
              else
                call dcopy(nbf*nbf,r_ofa(1,1,mst),1,r_ofa2(1,1,mst),1)
              endif
              if (.not.uwfn1) then
                call dcopy(nbf*nbf,r_noa2(1,1,mst),1,r_nob2(1,1,mst),1)
                call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_nhma(im,in),serr,On1,.false.)
              do jst = 1, mst
                r_nlma(im,in,jst) = On1(jst)
              enddo
            enddo
          enddo
c
c         Rotate beta natural orbitals
c
          if (uwfn1) then
          call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
          call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
          do in = 1, nbf
            do im = 1, in-1
c
c             Rotate beta natural orbitals
c
              call dcopy(nbf*nbf*nst,r_nob,1,r_nob2,1)
              call dcopy(nbf*nbf*nst,r_ofb,1,r_ofb2,1)
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         En0(mst),serr,On2,.false.)
              On0(mst) = 0.0d0
              do jst = 1, mst-1
                On0(mst) = On0(mst) + abs(On2(jst))
              enddo
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(im,im) = rotup(1,1)
              r_u(in,in) = rotup(2,2)
              r_u(in,im) = rotup(2,1)
              r_u(im,in) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_nob2(1,1,mst),nbf)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_u,nbf,
     +                     r_ofb(1,1,mst),nbf,0.0d0,r_ofb2(1,1,mst),nbf)
              else
                call dcopy(nbf*nbf,r_ofb(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_nhpb(im,in),serr,On1,.false.)
              do jst = 1, mst
                r_nlpb(im,in,jst) = On1(jst)
              enddo
              r_u(im,im) = rotdn(1,1)
              r_u(in,in) = rotdn(2,2)
              r_u(in,im) = rotdn(2,1)
              r_u(im,in) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,mst),nbf,
     +                   r_u,nbf,0.0d0,r_nob2(1,1,mst),nbf)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_u,nbf,
     +                     r_ofb(1,1,mst),nbf,0.0d0,r_ofb2(1,1,mst),nbf)
              else
                call dcopy(nbf*nbf,r_ofb(1,1,mst),1,r_ofb2(1,1,mst),1)
              endif
              call wfn1_engrad(rtdb,geom,mst,nbf,nea,neb,h1,eri,erix,ov,
     +                         sfac,ofac,
     +                         r_noa2,r_nob2,
     +                         r_ofa2,r_ofb2,
     +                         ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                         temperature,power,
     +                         nperma,npermb,perma,permb,signa,signb,
     +                         ovla,ovlb,
     +                         r_nhmb(im,in),serr,On1,.false.)
              do jst = 1, mst
                r_nlmb(im,in,jst) = On1(jst)
              enddo
            enddo
          enddo
          else
            call dcopy(nbf*nbf,r_nhpa,1,r_nhpb,1)
            call dcopy(nbf*nbf*mst,r_nlpa,1,r_nlpb,1)
            call dcopy(nbf*nbf,r_nhma,1,r_nhmb,1)
            call dcopy(nbf*nbf*mst,r_nlma,1,r_nlmb,1)
          endif
c
c         Calculate Lambda
c
          call dfill(nbf*nbf,0.0d0,r_doha,1)
          call dfill(nbf*nbf,0.0d0,r_dohb,1)
          call dfill(nbf*nbf*nst,0.0d0,r_dola,1)
          call dfill(nbf*nbf*nst,0.0d0,r_dolb,1)
          call dfill(nbf*nbf,0.0d0,r_dnha,1)
          call dfill(nbf*nbf,0.0d0,r_dnhb,1)
          call dfill(nbf*nbf*nst,0.0d0,r_dnla,1)
          call dfill(nbf*nbf*nst,0.0d0,r_dnlb,1)
cNEW
          En2(mst) = En0(mst)
          do jst = 1, mst
            En0(mst) = En0(mst) + dlambda_i(jst)*abs(On2(jst))
          enddo
          if (mst.eq.1) then
            do it = 1, nbf
              do is = 1, it-1
                do jst = 1, mst
                  r_ohpa(is,it) = r_ohpa(is,it)
     +                          + dlambda_i(jst)*abs(r_olpa(is,it,jst))
                  r_ohpb(is,it) = r_ohpb(is,it)
     +                          + dlambda_i(jst)*abs(r_olpb(is,it,jst))
                  r_ohma(is,it) = r_ohma(is,it)
     +                          + dlambda_i(jst)*abs(r_olma(is,it,jst))
                  r_ohmb(is,it) = r_ohmb(is,it)
     +                          + dlambda_i(jst)*abs(r_olmb(is,it,jst))
                enddo
              enddo
            enddo
          else
            do it = 1, nbf
              do is = 1, it-1
                r_ohpa(is,it) = En0(mst)+0.1d0
                r_ohpb(is,it) = En0(mst)+0.1d0
                r_ohma(is,it) = En0(mst)+0.1d0
                r_ohmb(is,it) = En0(mst)+0.1d0
              enddo
            enddo
          endif
          do in = 1, nbf
            do im = 1, in-1
              do jst = 1, mst
                r_nhpa(im,in) = r_nhpa(im,in)
     +                        + dlambda_i(jst)*abs(r_nlpa(im,in,jst))
                r_nhpb(im,in) = r_nhpb(im,in)
     +                        + dlambda_i(jst)*abs(r_nlpb(im,in,jst))
                r_nhma(im,in) = r_nhma(im,in)
     +                        + dlambda_i(jst)*abs(r_nlma(im,in,jst))
                r_nhmb(im,in) = r_nhmb(im,in)
     +                        + dlambda_i(jst)*abs(r_nlmb(im,in,jst))
              enddo
            enddo
          enddo
cNEW
          do it = 1, nbf
            do is = 1, it-1
c
              r_doha(is,it) = (r_ohpa(is,it)-r_ohma(is,it))/(2*h)
c             do jst = 1, mst
c               r_dola(is,it,jst)=(r_olpa(is,it,jst)-r_olma(is,it,jst))
c    +                           /(2*h)
c             enddo
              if (abs(r_doha(is,it)).lt.tol.and.
     +            r_ohma(is,it).lt.En0(mst)) then
                r_doha(is,it) = (En0(mst)-r_ohma(is,it))/(h)
c               do jst = 1, mst
c                 r_dola(is,it,jst) = (On2(jst)-r_olma(is,it,jst))/(h)
c               enddo
              endif
c
              r_dohb(is,it) = (r_ohpb(is,it)-r_ohmb(is,it))/(2*h)
c             do jst = 1, mst
c               r_dolb(is,it,jst)=(r_olpb(is,it,jst)-r_olmb(is,it,jst))
c    +                           /(2*h)
c             enddo
              if (abs(r_dohb(is,it)).lt.tol.and.
     +            r_ohmb(is,it).lt.En0(mst)) then
                r_dohb(is,it) = (En0(mst)-r_ohmb(is,it))/(h)
c               do jst = 1, mst
c                 r_dolb(is,it,jst) = (On2(jst)-r_olmb(is,it,jst))/(h)
c               enddo
              endif
c
              r_doha(it,is) = -r_doha(is,it)
              r_dohb(it,is) = -r_dohb(is,it)
c             do jst = 1, mst
c               r_dola(it,is,jst) = -r_dola(is,it,jst)
c               r_dolb(it,is,jst) = -r_dolb(is,it,jst)
c             enddo
c
            enddo
          enddo
c         do it = 1, nbf
c           do is = 1, nbf
c             dnrmh = dnrmh + r_doha(is,it)*r_dola(is,it)
c    +                      + r_dohb(is,it)*r_dolb(is,it)
c             dnrml = dnrml + r_dola(is,it)**2 + r_dolb(is,it)**2
c           enddo
c         enddo
c         d2nrm = 0.0d0
cDEBUG
c         write(*,*)'*** d2nrm A = ',d2nrm
cDEBUG
c         d2nrm = max(1.0d0,d2nrm)
          do in = 1, nbf
            do im = 1, in-1
c
              r_dnha(im,in) = (r_nhpa(im,in)-r_nhma(im,in))/(2*h)
c             do jst = 1, mst
c               r_dnla(im,in,jst)=(r_nlpa(im,in,jst)-r_nlma(im,in,jst))
c    +                           /(2*h)
c             enddo
              if (abs(r_dnha(im,in)).lt.tol.and.
     +            r_nhma(im,in).lt.En0(mst)) then
                r_dnha(im,in) = (En0(mst)-r_nhma(im,in))/(h)
c               do jst = 1, mst
c                 r_dnla(im,in,jst) = (On2(jst)-r_nlma(im,in,jst))/(h)
c               enddo
              endif
c
              r_dnhb(im,in) = (r_nhpb(im,in)-r_nhmb(im,in))/(2*h)
c             do jst = 1, mst
c               r_dnlb(im,in,jst)=(r_nlpb(im,in,jst)-r_nlmb(im,in,jst))
c    +                           /(2*h)
c             enddo
              if (abs(r_dnhb(im,in)).lt.tol.and.
     +            r_nhmb(im,in).lt.En0(mst)) then
                r_dnhb(im,in) = (En0(mst)-r_nhmb(im,in))/(h)
c               do jst = 1, mst
c                 r_dnlb(im,in,jst) = (On2(jst)-r_nlmb(im,in,jst))/(h)
c               enddo
              endif
c
              r_dnha(in,im) = -r_dnha(im,in)
              r_dnhb(in,im) = -r_dnhb(im,in)
c             do jst = 1, mst
c               r_dnla(in,im,jst) = -r_dnla(im,in,jst)
c               r_dnlb(in,im,jst) = -r_dnlb(im,in,jst)
c             enddo
c
            enddo
          enddo
c
          call dcopy(nbf*nbf,r_doha,1,d_ofa,1)
          call dcopy(nbf*nbf,r_dohb,1,d_ofb,1)
          call dcopy(nbf*nbf,r_dnha,1,d_noa,1)
          call dcopy(nbf*nbf,r_dnhb,1,d_nob,1)
c
c         dlambda(mst) = 0.0d0
c         d2nrm = 1.0d7
c         do jst = 1, mst-1
c           dnrmh = 0.0d0
c           dnrml = 0.0d0
c           do it = 1, nbf
c             do is = it+1, nbf
c               dnrmh = dnrmh + d_ofa(is,it)*r_dola(is,it,jst)
c    +                        + d_ofb(is,it)*r_dolb(is,it,jst)
c               dnrml = dnrml + r_dola(is,it,jst)**2
c    +                        + r_dolb(is,it,jst)**2
c             enddo
c           enddo
c           do in = 1, nbf
c             do im = in+1, nbf
c               dnrmh = dnrmh + d_noa(im,in)*r_dnla(im,in,jst)
c    +                        + d_nob(im,in)*r_dnlb(im,in,jst)
c               dnrml = dnrml + r_dnla(im,in,jst)**2
c    +                        + r_dnlb(im,in,jst)**2
c             enddo
c           enddo
c           dlambda_i(jst) = 0.0d0
c           if (dnrml.gt.1.0d-7*abs(dnrmh)) then
c             dlambda_i(jst) = -dnrmh/dnrml
c           else
c             dlambda_i(jst) =  0.0d0
c           endif
c           d2nrm = min(d2nrm,abs(dlambda_i(jst)))
c           dlambda(mst) = dlambda(mst) + abs(dlambda_i(jst))
c           do it = 1, nbf
c             do is = 1, nbf
c               d_ofa(is,it) = d_ofa(is,it)
c    +                       + dlambda_i(jst)*r_dola(is,it,jst)
c               d_ofb(is,it) = d_ofb(is,it)
c    +                       + dlambda_i(jst)*r_dolb(is,it,jst)
c             enddo
c           enddo
c           do in = 1, nbf
c             do im = 1, nbf
c               d_noa(im,in) = d_noa(im,in)
c    +                       + dlambda_i(jst)*r_dnla(im,in,jst)
c               d_nob(im,in) = d_nob(im,in)
c    +                       + dlambda_i(jst)*r_dnlb(im,in,jst)
c             enddo
c           enddo
c         enddo ! jst
c
          dnrm = 0.0d0
          do it = 1, nbf
            do is = 1, nbf
              dnrm = dnrm + d_ofa(is,it)**2 + d_ofb(is,it)**2
            enddo
          enddo
cDEBUG
c         write(*,*)'dnrm(ofa) = ',dnrm
c         write(*,*)'dnrmh,dnrml = ',dnrmh,dnrml
c         write(*,*)'grad H: nat orb ',iteration
c         write(*,*)'alpha'
c         call hess_hssout(r_dnha,nbf,nbf,nbf)
c         write(*,*)'beta'
c         call hess_hssout(r_dnhb,nbf,nbf,nbf)
c         do jst=1,mst-1
c           write(*,*)'grad S: nat orb'
c           write(*,*)'alpha: ',jst
c           call hess_hssout(r_dnla(1,1,jst),nbf,nbf,nbf)
c           write(*,*)'beta : ',jst
c           call hess_hssout(r_dnlb(1,1,jst),nbf,nbf,nbf)
c         enddo
c         write(*,*)
c         write(*,*)'grad H: occ fun ',iteration
c         write(*,*)'alpha'
c         call hess_hssout(r_doha,nbf,nbf,nbf)
c         write(*,*)'beta'
c         call hess_hssout(r_dohb,nbf,nbf,nbf)
c         do jst=1,mst-1
c           write(*,*)'grad S: occ fun'
c           write(*,*)'alpha: ',jst
c           call hess_hssout(r_dola(1,1,jst),nbf,nbf,nbf)
c           write(*,*)'beta : ',jst
c           call hess_hssout(r_dolb(1,1,jst),nbf,nbf,nbf)
c         enddo
c         write(*,*)
cDEBUG
          do in = 1, nbf
            do im = 1, nbf
              dnrm = dnrm + d_noa(im,in)**2 + d_nob(im,in)**2
            enddo
          enddo
cDEBUG
c         write(*,*)'dnrm(noa) = ',dnrm
c         write(*,*)'*** hvd: dnrm,dnh,dnl = ',dsqrt(dnrm),dsqrt(dnh),
c    +                                         dsqrt(dnl)
cDEBUG
          dnrmg = dnrm
          do jst = 1, mst-1
            dnrm = dnrm + On2(jst)**2
          enddo
cDEBUG
c         write(*,*)'dnrm(ovl) = ',dnrm
cDEBUG
          dnrm = sqrt(dnrm)
cDEBUG
c         write(*,*)'*** hvd: dnrmh, dnrml = ',dnrmh,dnrml
c         write(*,*)'*** hvd: lambda, norm = ',ist,dlambda,dnrm,On0(mst)
cDEBUG
c
          dnorm1(mst) = dnrm
c         if (iteration.gt.1) then
c           dnrmm = 0.0d0
c           do it = 1, nbf
c             do is = 1, it-1
c               d2nrm = (((r_ohpa(is,it)-2*En0(mst)+r_ohma(is,it))
c    +                + dlambda(mst)
c    +                * (r_olpa(is,it)-2*On0(mst)+r_olma(is,it)))
c    +                / (h*h))**1
c               dnrmm = dnrmm + (r_ofa2(is,it,1)*d2nrm)**2/dnrmg
c               d2nrm = max(1.0d0,abs(d2nrm))
c               if (.not.(abs(r_ofa2(is,it,1)).lt.h.and.
c    +                    abs(d2nrm).lt.h)) then
c                 r_ofa2(is,it,1) = r_ofa2(is,it,1)/d2nrm
c                 r_ofa2(it,is,1) = r_ofa2(it,is,1)/d2nrm
c               endif
c               d2nrm = (((r_ohpb(is,it)-2*En0(mst)+r_ohmb(is,it))
c    +                + dlambda(mst)
c    +                * (r_olpb(is,it)-2*On0(mst)+r_olmb(is,it)))
c    +                / (h*h))**1
c               dnrmm = dnrmm + (r_ofb2(is,it,1)*d2nrm)**2/dnrmg
c               d2nrm = max(1.0d0,abs(d2nrm))
c               if (.not.(abs(r_ofb2(is,it,1)).lt.h.and.
c    +                    abs(d2nrm).lt.h)) then
c                 r_ofb2(is,it,1) = r_ofb2(is,it,1)/d2nrm
c                 r_ofb2(it,is,1) = r_ofb2(it,is,1)/d2nrm
c               endif
c             enddo
c           enddo
c           do im = 1, nbf
c             do in = 1, im-1
c               d2nrm = (((r_nhpa(in,im)-2*En0(mst)+r_nhma(in,im))
c    +                + dlambda(mst)
c    +                * (r_nlpa(in,im)-2*On0(mst)+r_nlma(in,im)))
c    +                / (h*h))**1
c               dnrmm = dnrmm + (r_noa2(is,it,1)*d2nrm)**2/dnrmg
c               d2nrm = max(1.0d0,abs(d2nrm))
c               if (.not.(abs(r_noa2(in,im,1)).lt.h.and.
c    +                    abs(d2nrm).lt.h)) then
c                 r_noa2(in,im,1) = r_noa2(in,im,1)/d2nrm
c                 r_noa2(im,in,1) = r_noa2(im,in,1)/d2nrm
c               endif
c               d2nrm = (((r_nhpb(in,im)-2*En0(mst)+r_nhmb(in,im))
c    +                + dlambda(mst)
c    +                * (r_nlpb(in,im)-2*On0(mst)+r_nlmb(in,im)))
c    +                / (h*h))**1
c               dnrmm = dnrmm + (r_nob2(is,it,1)*d2nrm)**2/dnrmg
c               d2nrm = max(1.0d0,abs(d2nrm))
c               if (.not.(abs(r_nob2(in,im,1)).lt.h.and.
c    +                    abs(d2nrm).lt.h)) then
c                 r_nob2(in,im,1) = r_nob2(in,im,1)/d2nrm
c                 r_nob2(im,in,1) = r_nob2(im,in,1)/d2nrm
c               endif
c             enddo
c           enddo
cDEBUG
c           write(*,*)'*** hvd: dnrmm = ',mst,sqrt(dnrmm)
cDEBUG
c           dnrmm = min(max(0.50d0,sqrt(dnrmm)),25.0d0)
c           do it = 1, nbf
c             do is = 1, it-1
c               r_ofa2(is,it,1) = r_ofa2(is,it,1)/dnrmm
c               r_ofa2(it,is,1) = r_ofa2(it,is,1)/dnrmm
c               r_ofb2(is,it,1) = r_ofb2(is,it,1)/dnrmm
c               r_ofb2(it,is,1) = r_ofb2(it,is,1)/dnrmm
c             enddo
c           enddo
c           do im = 1, nbf
c             do in = 1, im-1
c               r_noa2(in,im,1) = r_noa2(in,im,1)/dnrmm
c               r_noa2(im,in,1) = r_noa2(im,in,1)/dnrmm
c               r_nob2(in,im,1) = r_nob2(in,im,1)/dnrmm
c               r_nob2(im,in,1) = r_nob2(im,in,1)/dnrmm
c             enddo
c           enddo
c         endif
c         d2nrm = 2*d2nrm
c         d2nrm = sqrt(d2nrm)
cDEBUG
c         write(*,*)'*** d2nrm B = ',d2nrm
cDEBUG
          do jst = 1, mst-1
c           en0(mst) = en0(mst)+dlambda_i(jst)*abs(on2(jst))
            dlambda_k(jst) = 1.0d8
          enddo
          d2nrm = 1.0d0
c         if (d2nrm.gt.1.0d-10) then
            call wfn1_linesearch(rtdb,geom,uwfn1,mst,nbf,nea,neb,
     +      h1,eri,erix,ov,
     +      En0(mst),dlambda_i,temperature,
     +      power,tol,nperma,npermb,perma,permb,signa,signb,
     +      ovla,ovlb,r_noa,r_nob,r_ofa,r_ofb,
     +      r_noa2,r_nob2,r_ofa2,r_ofb2,
     +      d_noa,d_nob,d_ofa,d_ofb,d_noa2,d_nob2,d_ofa2,d_ofb2,
     +      rnoa,rnob,rofa,rofb,m1,m2,m3,step)
c         else
c           step = 1.0d-7
c         endif
c
c         Apply the rotations. To maintain orthonormality of all 
c         states the rotations have to be applied to the current
c         state and all higher excited states. For the higher excited
c         states this requires some extra work as they are represented
c         in a different basis. So we need to transform the rotation.
c         To limit the amount of workspace needed we first do the
c         alpha-electron transformations and then the beta electron
c         ones.
cDEBUG
c         write(*,*)"*** mem: r_noa  = ",loc(r_noa(1,1,1)),
c    +                                   loc(r_noa(nbf,nbf,nst))
c         write(*,*)"*** mem: r_nob  = ",loc(r_nob(1,1,1)),
c    +                                   loc(r_nob(nbf,nbf,nst))
c         write(*,*)"*** mem: r_ofa  = ",loc(r_ofa(1,1,1)),
c    +                                   loc(r_ofa(nbf,nbf,nst))
c         write(*,*)"*** mem: r_ofb  = ",loc(r_ofb(1,1,1)),
c    +                                   loc(r_ofb(nbf,nbf,nst))
c         write(*,*)"*** mem: r_noa2 = ",loc(r_noa2(1,1,1)),
c    +                                   loc(r_noa2(nbf,nbf,nst))
c         write(*,*)"*** mem: r_nob2 = ",loc(r_nob2(1,1,1)),
c    +                                   loc(r_nob2(nbf,nbf,nst))
c         write(*,*)"*** mem: r_ofa2 = ",loc(r_ofa2(1,1,1)),
c    +                                   loc(r_ofa2(nbf,nbf,nst))
c         write(*,*)"*** mem: r_ofb2 = ",loc(r_ofb2(1,1,1)),
c    +                                   loc(r_ofb2(nbf,nbf,nst))
c         write(*,*)"*** mem: d_noa2 = ",loc(d_noa2(1,1)),
c    +                                   loc(d_noa2(nbf,nbf))
c         write(*,*)"*** mem: d_nob2 = ",loc(d_nob2(1,1)),
c    +                                   loc(d_nob2(nbf,nbf))
c         write(*,*)"*** mem: d_ofa2 = ",loc(d_ofa2(1,1)),
c    +                                   loc(d_ofa2(nbf,nbf))
c         write(*,*)"*** mem: d_ofb2 = ",loc(d_ofb2(1,1)),
c    +                                   loc(d_ofb2(nbf,nbf))
c         write(*,*)"*** mem: m1     = ",loc(m1(1,1)),
c    +                                   loc(m1(nbf,nbf))
c         write(*,*)"*** mem: m2     = ",loc(m2(1,1)),
c    +                                   loc(m2(nbf,nbf))
c         write(*,*)"*** mem: m3     = ",loc(m3(1,1)),
c    +                                   loc(m3(nbf,nbf))
c         write(*,*)"*** mem: rnoa   = ",loc(rnoa(1,1)),
c    +                                   loc(rnoa(nbf,nbf))
c         write(*,*)"*** mem: rnob   = ",loc(rnob(1,1)),
c    +                                   loc(rnob(nbf,nbf))
c         write(*,*)"*** mem: rofa   = ",loc(rofa(1,1)),
c    +                                   loc(rofa(nbf,nbf))
c         write(*,*)"*** mem: rofb   = ",loc(rofb(1,1)),
c    +                                   loc(rofb(nbf,nbf))
c         write(*,*)"*** mem: ovnoa  = ",loc(ovnoa(1,1)),
c    +                                   loc(ovnoa(nbf,nbf))
c         write(*,*)"*** mem: ovnob  = ",loc(ovnob(1,1)),
c    +                                   loc(ovnob(nbf,nbf))
c         write(*,*)"*** mem: ovofa  = ",loc(ovofa(1,1)),
c    +                                   loc(ovofa(nbf,nbf))
c         write(*,*)"*** mem: ovofb  = ",loc(ovofb(1,1)),
c    +                                   loc(ovofb(nbf,nbf))
c         write(*,*)"*** mem: d_noa  = ",loc(d_noa(1,1)),
c    +                                   loc(d_noa(nbf,nbf))
c         write(*,*)"*** mem: d_nob  = ",loc(d_nob(1,1)),
c    +                                   loc(d_nob(nbf,nbf))
c         write(*,*)"*** mem: d_ofa  = ",loc(d_ofa(1,1)),
c    +                                   loc(d_ofa(nbf,nbf))
c         write(*,*)"*** mem: d_ofb  = ",loc(d_ofb(1,1)),
c    +                                   loc(d_ofb(nbf,nbf))
c         write(*,*)
cDEBUG
c
c         Calculate the natural orbital rotation matrix
c
c         d2nrm = 1.0d0
          call dscal(nbf*nbf,-step,d_noa,1)
          call dscal(nbf*nbf,-step,d_nob,1)
cDEBUG
c         if (mst.gt.1) then
c           call dfill(nbf*nbf,0.0d0,r_noa2,1)
c           call dfill(nbf*nbf,0.0d0,r_nob2,1)
c         endif
cDEBUG
          call wfn1_exp(nbf,d_noa,m1,m2,rnoa)
          call wfn1_exp(nbf,d_nob,m1,m2,rnob)
cDEBUG
c         if (mst.eq.1) then
c             call dgemm('n','t',nbf,nbf,nbf,1.0d0,rnoa,nbf,
c    +                   rnoa,nbf,0.0d0,m1,nbf)
c             write(*,*)'*** rnoa unitary?'
c             call hess_hssout(m1,nbf,nbf,nbf)
c             call dgemm('n','t',nbf,nbf,nbf,1.0d0,rnob,nbf,
c    +                   rnob,nbf,0.0d0,m1,nbf)
c             write(*,*)'*** rnob unitary?'
c             call hess_hssout(m1,nbf,nbf,nbf)
c         endif
cDEBUG
c
c         Calculate the occupation function rotation matrix
c
c         d2nrm = 1.0d0
          call dscal(nbf*nbf,-step,d_ofa,1)
          call dscal(nbf*nbf,-step,d_ofb,1)
cDEBUG
          if (mst.gt.1) then
            call dfill(nbf*nbf,0.0d0,d_ofa,1)
            call dfill(nbf*nbf,0.0d0,d_ofb,1)
          endif
cDEBUG
          call wfn1_exp(nbf,d_ofa,m1,m2,rofa)
          call wfn1_exp(nbf,d_ofb,m1,m2,rofb)
cDEBUG
c         write(*,*)'*** rnoa:'
c         call hess_hssout(rnoa,nbf,nbf,nbf)
c         write(*,*)'*** rnob:'
c         call hess_hssout(rnob,nbf,nbf,nbf)
c         write(*,*)'*** rofa:'
c         call hess_hssout(rofa,nbf,nbf,nbf)
c         write(*,*)'*** rofb:'
c         call hess_hssout(rofb,nbf,nbf,nbf)
c         write(*,*)
cDEBUG
c
c         Now actually rotate the functions
c
          if (dnorm1(mst).gt.tol)
     +    then
            oprint = .false.
c
c           Rotate higher excited states
c
            do jst = mst+1, nst
              oconv(jst) = .false.
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,jst),nbf,
     +                   ov,nbf,0.0d0,m1,nbf)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   r_noa(1,1,mst),nbf,0.0d0,ovnoa,nbf)
c
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,jst),nbf,
     +                   ov,nbf,0.0d0,m1,nbf)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   r_nob(1,1,mst),nbf,0.0d0,ovnob,nbf)
c
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,jst),nbf,
     +                   ovnoa,nbf,0.0d0,m1,nbf)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   r_ofa(1,1,mst),nbf,0.0d0,ovofa,nbf)
c
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,jst),nbf,
     +                   ovnob,nbf,0.0d0,m1,nbf)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   r_ofb(1,1,mst),nbf,0.0d0,ovofb,nbf)
c
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,ovnoa,nbf,
     +                   rnoa,nbf,0.0d0,m1,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   ovnoa,nbf,0.0d0,m2,nbf)
cDEBUG
c             write(50,*)'wfn1_u a: ist,jst=',ist,jst
cDEBUG
              call wfn1_unitarize(nbf,m2,ov,m1,m3,.false.)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,jst),nbf,
     +                   m2,nbf,0.0d0,r_noa2(1,1,jst),nbf)
              call dcopy(nbf*nbf,r_noa2(1,1,jst),1,r_noa(1,1,jst),1)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,m2,nbf,
     +                     r_ofa(1,1,jst),nbf,0.0d0,r_ofa2(1,1,jst),nbf)
                call dcopy(nbf*nbf,r_ofa2(1,1,jst),1,r_ofa(1,1,jst),1)
              endif
c
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,ovnob,nbf,
     +                   rnob,nbf,0.0d0,m1,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   ovnob,nbf,0.0d0,m2,nbf)
cDEBUG
c             write(50,*)'wfn1_u b: ist,jst=',ist,jst
cDEBUG
              call wfn1_unitarize(nbf,m2,ov,m1,m3,.false.)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,jst),nbf,
     +                   m2,nbf,0.0d0,r_nob2(1,1,jst),nbf)
              call dcopy(nbf*nbf,r_nob2(1,1,jst),1,r_nob(1,1,jst),1)
              if (mst.eq.1) then
                call dgemm('t','n',nbf,nbf,nbf,1.0d0,m2,nbf,
     +                     r_ofb(1,1,jst),nbf,0.0d0,r_ofb2(1,1,jst),nbf)
                call dcopy(nbf*nbf,r_ofb2(1,1,jst),1,r_ofb(1,1,jst),1)
              endif
c
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,ovofa,nbf,
     +                   rofa,nbf,0.0d0,m1,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   ovofa,nbf,0.0d0,m2,nbf)
cDEBUG
c             write(50,*)'wfn1_u c: ist,jst=',ist,jst
cDEBUG
              call wfn1_unitarize(nbf,m2,ov,m1,m3,.false.)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,jst),nbf,
     +                   m2,nbf,0.0d0,r_ofa2(1,1,jst),nbf)
              call dcopy(nbf*nbf,r_ofa2(1,1,jst),1,r_ofa(1,1,jst),1)
c
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,ovofb,nbf,
     +                   rofb,nbf,0.0d0,m1,nbf)
              call dgemm('n','t',nbf,nbf,nbf,1.0d0,m1,nbf,
     +                   ovofb,nbf,0.0d0,m2,nbf)
cDEBUG
c             write(50,*)'wfn1_u d: ist,jst=',ist,jst
cDEBUG
              call wfn1_unitarize(nbf,m2,ov,m1,m3,.false.)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,jst),nbf,
     +                   m2,nbf,0.0d0,r_ofb2(1,1,jst),nbf)
              call dcopy(nbf*nbf,r_ofb2(1,1,jst),1,r_ofb(1,1,jst),1)
c
              if (.not.uwfn1) then
                call dcopy(nbf*nbf,r_noa(1,1,jst),1,r_nob(1,1,jst),1)
                call dcopy(nbf*nbf,r_ofa(1,1,jst),1,r_ofb(1,1,jst),1)
              endif
c
            enddo
            oconv(mst) = .false.
c
c           Rotate natural orbitals
c
c           - First rotate the natural orbitals among each other
c
            call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,mst),nbf,
     +                 rnoa,nbf,0.0d0,r_noa2(1,1,mst),nbf)
            call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,mst),nbf,
     +                 rnob,nbf,0.0d0,r_nob2(1,1,mst),nbf)
            call dcopy(nbf*nbf,r_noa2(1,1,mst),1,r_noa(1,1,mst),1)
            call dcopy(nbf*nbf,r_nob2(1,1,mst),1,r_nob(1,1,mst),1)
c
c           - Second rotate the basis of the occupation functions
c
            if (mst.eq.1) then
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnoa,nbf,
     +                   r_ofa(1,1,mst),nbf,0.0d0,r_ofa2(1,1,mst),nbf)
              call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnob,nbf,
     +                   r_ofb(1,1,mst),nbf,0.0d0,r_ofb2(1,1,mst),nbf)
              call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofa(1,1,mst),1)
              call dcopy(nbf*nbf,r_ofb2(1,1,mst),1,r_ofb(1,1,mst),1)
            endif
c
c           Rotate occupation functions among each other
c
            call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,mst),nbf,
     +                 rofa,nbf,0.0d0,r_ofa2(1,1,mst),nbf)
            call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,mst),nbf,
     +                 rofb,nbf,0.0d0,r_ofb2(1,1,mst),nbf)
            call dcopy(nbf*nbf,r_ofa2(1,1,mst),1,r_ofa(1,1,mst),1)
            call dcopy(nbf*nbf,r_ofb2(1,1,mst),1,r_ofb(1,1,mst),1)
          else
            oprint = .true.
            do jst = mst+1, nst
              oconv(jst) = .false.
            enddo
          endif
          call dcopy(nbf*nbf*nst,r_ofa,1,r_ofa2,1)
          call dcopy(nbf*nbf*nst,r_ofb,1,r_ofb2,1)
          call dcopy(nbf*nbf*nst,r_noa,1,r_noa2,1)
          call dcopy(nbf*nbf*nst,r_nob,1,r_nob2,1)
c
          if (.not.uwfn1) then
            call dcopy(nbf*nbf,r_noa(1,1,mst),1,r_nob(1,1,mst),1)
            call dcopy(nbf*nbf,r_ofa(1,1,mst),1,r_ofb(1,1,mst),1)
          endif
c
c         call wfn1_engrad(rtdb,geom,ist,nbf,nea,neb,h1,eri,erix,ov,
c    +                     sfac,ofac,
c    +                     r_noa2,r_nob2,
c    +                     r_ofa2,r_ofb2,
c    +                     ehfa,ehfb,ewfa,ewfb,enta,entb,
c    +                     temperature,power,
c    +                     nperma,npermb,perma,permb,signa,signb,
c    +                     ovla,ovlb,
c    +                     En0(mst),serr,On0(mst),.true.)
cDEBUG
c         write(*,*)'*** Rot done: ',ist
cDEBUG
c         sfac = sfac * (1.0d0 + max(0.0d0,serr))
c         ofac = ofac * (1.0d0 + max(0.0d0,oerr(mst)))
cDEBUG
c         call dfill(nbf*nbf,0.0d0,r_noa2,1)
c         Em1 = 0.0d0
c         call wfn1_print_occ(nbf,nea,neb,r_ofa,r_ofb,r_noa2,r_noa2,
c    +                        Em1,r_nob2)
c         do is = 1, nbf
c           write(*,'(" ehf a,b: ",3f20.8,3x,3f20.8)')
c    +            enta(is),ehfa(is),ewfa(is),
c    +            entb(is),ehfb(is),ewfb(is)
c         enddo
          write(*,*)'*** dnorm,dnorm1=',dnorm,dnorm1(mst)
cDEBUG
c         dnorm = max(dnorm,dnorm1(mst))
          dnorm = dnorm1(mst)
          if (.not.(dnorm.gt.tol).and.mst.lt.nst) then
            dnorm = 1.0d0
            mst = mst + 1
            do jst = 1, nst
              dlambda_i(jst) = max(dlambda_i(jst),abs(En0(mst-1)))
            enddo
          endif
c***    enddo
c
c       Evaluate the current energy for every state
c
        do ist = 1, nst
          call wfn1_engrad(rtdb,geom,ist,nbf,nea,neb,h1,eri,erix,ov,
     +                     sfac,ofac,
     +                     r_noa,r_nob,
     +                     r_ofa,r_ofb,
     +                     ehfa,ehfb,ewfa,ewfb,enta,entb,
     +                     temperature,power,
     +                     nperma,npermb,perma,permb,signa,signb,
     +                     ovla,ovlb,
     +                     En0(ist),serr,On2,
     +                     (is.eq.1.and.it.eq.2).and..false.)
          On0(ist) = 0.0d0
          dlambda(ist) = 0.0d0
          do jst = 1, ist-1
            On0(ist) = On0(ist) + abs(On2(jst))
            dlambda(ist) = dlambda(ist) + abs(On2(jst))
          enddo
        enddo
      enddo
      write(*,*)
      write(*,*)"Final results"
      write(*,*)"============="
      write(*,*)
      do ist = 1, nst
        write(luout,'(1x,i5,f22.8,2f22.8,f8.3)')iteration,en0(ist),
     +              dnorm1(ist),On0(ist),damp(ist)
      enddo
c
      write(*,*)
      write(*,*)"Wavefunction analysis"
      write(*,*)"====================="
      write(*,*)
      do ist = 1, nst
        call dfill(nbf*nbf,0.0d0,r_noa2,1)
        call dfill(nbf*nbf,0.0d0,r_nob2,1)
        write(*,*)'State: ',ist
        call wfn1_print_occ(nbf,nea,neb,r_ofa(1,1,ist),r_ofb(1,1,ist),
     +                      r_noa2,r_noa2(1,2,1),en0,r_nob2)
        write(*,*)'Natural Orbital State: ',ist
        call hess_hssout(r_noa(1,1,ist),nbf,nbf,nbf)
        call hess_hssout(r_nob(1,1,ist),nbf,nbf,nbf)
        write(*,*)'Occupation Functions State: ',ist
        call hess_hssout(r_ofa(1,1,ist),nbf,nbf,nbf)
        call hess_hssout(r_ofb(1,1,ist),nbf,nbf,nbf)
      enddo
c
      end
C> @}
