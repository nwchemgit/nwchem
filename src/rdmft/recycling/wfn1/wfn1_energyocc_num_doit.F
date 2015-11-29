C> \ingroup wfn1
C> @{
C>
C> \brief Go and numerically optimize the energy in wfn1 formalism
C> under the constraint that the occupation numbers match the input
C>
C> Go and actually compute the total energy within the wfn1
C> formalism under the constraint that the occupation numbers.
C> It is assumed that all the memory has been arranged in
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
      subroutine wfn1_energyocc_num_doit(rtdb,geom,basis,tol,nbf,nst,
     +           nea,neb,h1,eri,erix,ov,r_noa,r_nob,r_ofa,r_ofb,
     +           r_noa2,r_nob2,r_ofa2,r_ofb2,rota,rotb,sm1,En0,r_u,r_l,
     +           ehfa,ehfb,ewfa,ewfb,enta,entb,power,occa,occb,
     +           tol_act)
      implicit none
c
#include "stdio.fh"
c
      integer rtdb     !< [Input] The runtime database handle
      integer geom     !< [Input] The geometry handle
      integer basis    !< [Input] The basis set handle
      integer nbf      !< [Input] The number of basis functions
      integer nst      !< [Input] The number of states
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
      double precision r_noa2(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron natural orbitals
      double precision r_nob2(nbf,nbf,nst) !< [Scratch] The beta
      !< electron natural orbitals
      double precision r_ofa2(nbf,nbf,nst) !< [Scratch] The alpha
      !< electron occupation functions
      double precision r_ofb2(nbf,nbf,nst) !< [Scratch] The beta
      !< electron occupation functions
      double precision rota(nbf,nbf) !< [Scratch] The alpha
      !< orbital rotation matrix
      double precision rotb(nbf,nbf) !< [Scratch] The beta
      !< orbital rotation matrix
      double precision sm1(nbf,nbf) !< [Scratch] Matrix
c
      double precision En0(nst) !< [Output] The total energies
c
      double precision r_u(nbf,nbf) !< [Scratch] The rotation matrix
      double precision r_l(nbf,nbf) !< [Scratch] The lagrangian matrix
c
      double precision ehfa(nbf) !< [Input] The Hartree-Fock alpha
      !< orbital energies.
      double precision ehfb(nbf) !< [Input] The Hartree-Fock beta
      !< orbital energies.
      double precision ewfa(nbf) !< [Output] The WFN1 alpha
      !< orbital energies.
      double precision ewfb(nbf) !< [Output] The WFN1 beta
      !< orbital energies.
      double precision enta(nbf) !< [Output] The alpha orbital entropies
      double precision entb(nbf) !< [Output] The beta orbital entropies
      double precision power     !< [Input] The matrix power
c
      double precision occa(nbf) !< [Input] The \f$\alpha\f$ electron
      !< occupation numbers
      double precision occb(nbf) !< [Input] The \f$\beta\f$ electron
      !< occupation numbers
c
      double precision tol_act !< [Input] The threshold for active
      !< natural orbitals
c
c     Local variables
c
      double precision oerr       !< The orbital order error
      double precision occ_diff   !< The deviation of the occupation
                                  !< numbers
      double precision occdmax    !< The maximum deviation
      double precision factor     !< The scale factor for the occupation
                                  !< penalty function
      double precision ofac       !< The scale factor for the orbital
                                  !< order penalty function
      integer iteration           !< The iteration counter
      double precision dnorm      !< The gradient norm
      double precision damp       !< The damping factor
                                  !< for an energy minimization
      double precision h          !< The step size
      double precision h0         !< The step size
      double precision d          !< The displacement
      double precision rotup(2,2) !< Rotation over positive angle
      double precision rotdn(2,2) !< Rotation over negative angle
      double precision Em1(nst)   !< The total energies
      double precision Ep1(nst)   !< The total energies
      double precision dE         !< The energy gradient
      double precision d2E        !< The energy Hessian
      double precision serr       !< The error in S^2
      double precision sfac       !< The scale factor for the spin
                                  !< penality function
      double precision dv0        !< The deviation at the reference
                                  !< point
      double precision dvm        !< The deviation negative rotation
      double precision dvp        !< The deviation positive rotation
      integer is                  !< Counter over occupation functions
      integer it                  !< Counter over occupation functions
      integer iu                  !< Counter over occupation functions
      integer im                  !< Counter over natural orbitals
      integer in                  !< Counter over natural orbitals
      integer ist                 !< Counter over states
c
c     Functions
c
      double precision wfn1_norm, ddot
      external         wfn1_norm, ddot
c
c     Optimal h is the square root of the machine precision
c
      factor = 1.0d0
      ofac   = 3.0d0
      tol_act = -1.0d0
      h = 1.0d-8
cDEBUG
      write(*,*)'convergence tolerance = ',tol
      h = 1.0d-6
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
      call dfill(nbf*nbf,0.0d0,r_l,1)
c
      call int_init(rtdb,1,basis)
      call wfn1_1e_tv(basis,nbf,h1)
      call wfn1_1e_s(basis,nbf,ov)
      call wfn1_2e_rep(basis,nbf,eri,erix)
      call int_terminate
cDEBUG
      write(*,*)'2-electron ints:',ddot(nbf**4,eri,1,eri,1),
     +                             ddot(nbf**4,erix,1,erix,1)

c     write(*,*)'1-electron Hamiltonian'
c     call hess_hssout(h1,nbf,nbf,nbf)
c     write(*,*)'overlap matrix'
c     call hess_hssout(ov,nbf,nbf,nbf)
cDEBUG
c
      do ist = 1, nst
        iteration = 0
        damp = 1.000d0
        dnorm = 10.0d0
        occdmax = 0.0d0
        occ_diff = 0.0d0
        call wfn1_print_occ(nbf,nea,neb,r_ofa(1,1,ist),r_ofb(1,1,ist),
     +                      occa,occb,En0(ist),r_ofa2)
        call dcopy(nbf*nbf,r_ofa(1,1,ist),1,r_ofa2,1)
        call dcopy(nbf*nbf,r_ofb(1,1,ist),1,r_ofb2,1)
        call dcopy(nbf*nbf,r_noa(1,1,ist),1,r_noa2,1)
        call dcopy(nbf*nbf,r_nob(1,1,ist),1,r_nob2,1)
        call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,factor,ofac,
     +                      r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                      ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                      occa,occb,En0(ist),dv0,oerr,.true.)
        occdmax = dv0
cDEBUG
        do is = 1, nbf
          write(*,'(" ent a,b: ",3f20.8,3x,3f20.8)')
     +          enta(is),ehfa(is),ewfa(is),
     +          entb(is),ehfb(is),ewfb(is)
        enddo
cDEBUG
c
        write(luout,'(2x,"iter",17x,"energy",6x,"grad",6x,"occ diff",
     +                6x,"order diff",6x,"tolerance",4x,"damping",
     +                4x,"occ factor",4x,"order factor")')
        do while (dnorm.gt.tol.or.occdmax.gt.tol.or.iteration.lt.3)
cDEBUG
c         call wfn1_print_occ(nbf,nea,neb,r_ofa,r_ofb,occa,occb,En0,
c    +                        r_ofa2)
c         write(*,*)
c         write(*,*)'alpha natural orbitals'
c         call hess_hssout(r_noa,nbf,nbf,nbf)
c         write(*,*)
c         write(*,*)'beta natural orbitals'
c         call hess_hssout(r_nob,nbf,nbf,nbf)
c         write(*,*)
c         write(*,*)'alpha occupation functions'
c         call hess_hssout(r_ofa,nbf,nbf,nbf)
c         write(*,*)
c         write(*,*)'beta occupation functions'
c         call hess_hssout(r_ofb,nbf,nbf,nbf)
c         write(*,*)
cDEBUG
          h0 = h/sqrt(factor*ofac)
          rotup(1,1) =  cos(h0)
          rotup(2,1) = -sin(h0)
          rotup(1,2) =  sin(h0)
          rotup(2,2) =  cos(h0)
          rotdn(1,1) =  cos(-h0)
          rotdn(2,1) = -sin(-h0)
          rotdn(1,2) =  sin(-h0)
          rotdn(2,2) =  cos(-h0)
          iteration = iteration + 1
          write(luout,'(1x,i5,i3," iter ",5f20.8,3f10.3)')
     +          iteration,ist,en0(ist),dnorm,occdmax,oerr,tol,
     +          damp,factor,ofac
          dnorm = 0.0d0
          occdmax = 0.0d0
c
c         Rotate alpha occupation functions
c
          call dfill(nbf*nbf,0.0d0,rota,1)
          call dfill(nbf*nbf,0.0d0,rotb,1)
          do iu = 1, nbf
             rota(iu,iu) = 1.0d0
             rotb(iu,iu) = 1.0d0
          enddo
          call dcopy(nbf*nbf,r_ofb(1,1,ist),1,r_ofb2,1)
          call dcopy(nbf*nbf,r_noa(1,1,ist),1,r_noa2,1)
          call dcopy(nbf*nbf,r_nob(1,1,ist),1,r_nob2,1)
          do it = 1, nbf
            do is = 1, it-1
              call dcopy(nbf*nbf,r_ofa(1,1,ist),1,r_ofa2,1)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,En0(ist),dv0,oerr,.false.)
              occdmax = max(occdmax,dv0)
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(is,is) = rotup(1,1)
              r_u(it,it) = rotup(2,2)
              r_u(it,is) = rotup(2,1)
              r_u(is,it) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,ist),nbf,
     +                   r_u,nbf,0.0d0,r_ofa2,nbf)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,Ep1(ist),dvp,oerr,.false.)
c             occdmax = max(occdmax,occ_diff)
              r_u(is,is) = rotdn(1,1)
              r_u(it,it) = rotdn(2,2)
              r_u(it,is) = rotdn(2,1)
              r_u(is,it) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,ist),nbf,
     +                   r_u,nbf,0.0d0,r_ofa2,nbf)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,Em1(ist),dvm,oerr,.false.)
c             occdmax = max(occdmax,occ_diff)
              dE = (dvp-dvm)/(2*h0)
              d2E = (dvp-2*dv0+dvm)/(h0*h0)
              d = -dE/max(1.0d0,abs(d2E))
cDEBUG
c             write(*,*)'*** ofa ',it,is,Em1(1),En0(1),Ep1(1)
c             write(*,*)'*** ofa ',it,is,dE,d2E,d
cDEBUG
c             damp = 1.0d0/(1.0d0 + abs(ehfa(it)-ehfa(is)))
              damp = exp(-abs(ehfa(it)-ehfa(is)))/sqrt(factor*ofac)
              if (abs(dE).lt.tol.and.dvm.lt.dv0) then
c               This is a meta-stable point, pick any rotation to get
c               out of it.
                d = -h0/damp
              endif
              d = damp*d
              dnorm = max(abs(dE),dnorm)
              r_u(is,is) =  cos(d)
              r_u(it,it) =  cos(d)
              r_u(it,is) = -sin(d)
              r_u(is,it) =  sin(d)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,rota,nbf,r_u,nbf,
     +                   0.0d0,sm1,nbf)
              call dcopy(nbf*nbf,sm1,1,rota,1)
c             call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
c    +                            r_noa2,r_nob2,r_ofa2,r_ofb2,
c    +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
c    +                            occa,occb,En0,.false.)
c             occdmax = max(occdmax,occ_diff)
c             call dcopy(nbf*nbf,r_ofa2,1,r_ofa,1)
            enddo
          enddo
cDEBUG
c         write(*,*)
cDEBUG
c
c         Rotate beta occupation functions
c
          call dcopy(nbf*nbf,r_ofa(1,1,ist),1,r_ofa2,1)
          call dcopy(nbf*nbf,r_noa(1,1,ist),1,r_noa2,1)
          call dcopy(nbf*nbf,r_nob(1,1,ist),1,r_nob2,1)
          do it = 1, nbf
            do is = 1, it-1
              call dcopy(nbf*nbf,r_ofb(1,1,ist),1,r_ofb2,1)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,En0(ist),dv0,oerr,.false.)
              occdmax = max(occdmax,dv0)
              call dfill(nbf*nbf,0.0d0,r_u,1)
              do iu = 1, nbf
                r_u(iu,iu) = 1.0d0
              enddo
              r_u(is,is) = rotup(1,1)
              r_u(it,it) = rotup(2,2)
              r_u(it,is) = rotup(2,1)
              r_u(is,it) = rotup(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,ist),nbf,
     +                   r_u,nbf,0.0d0,r_ofb2,nbf)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,Ep1(ist),dvp,oerr,.false.)
c             occdmax = max(occdmax,occ_diff)
              r_u(is,is) = rotdn(1,1)
              r_u(it,it) = rotdn(2,2)
              r_u(it,is) = rotdn(2,1)
              r_u(is,it) = rotdn(1,2)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,ist),nbf,
     +                   r_u,nbf,0.0d0,r_ofb2,nbf)
              call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
     +                            factor,
     +                            ofac,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                            occa,occb,Em1(ist),dvm,oerr,.false.)
c             occdmax = max(occdmax,occ_diff)
              dE = (dvp-dvm)/(2*h0)
              d2E = (dvp-2*dv0+dvm)/(h0*h0)
              d = -dE/max(1.0d0,abs(d2E))
cDEBUG
c             write(*,*)'*** ofb ',it,is,Em1(1),En0(1),Ep1(1)
c             write(*,*)'*** ofb ',it,is,dE,d2E,d
cDEBUG
c             damp = 1.0d0/(1.0d0 + abs(ehfb(it)-ehfb(is)))
              damp = exp(-abs(ehfb(it)-ehfb(is)))/sqrt(factor*ofac)
              if (abs(dE).lt.tol.and.dvm.lt.dv0) then
c               This is a meta-stable point, pick any rotation to get
c               out of it.
                d = -h0/damp
              endif
              d = damp*d
              dnorm = max(abs(dE),dnorm)
              r_u(is,is) =  cos(d)
              r_u(it,it) =  cos(d)
              r_u(it,is) = -sin(d)
              r_u(is,it) =  sin(d)
              call dgemm('n','n',nbf,nbf,nbf,1.0d0,rotb,nbf,r_u,nbf,
     +                   0.0d0,sm1,nbf)
              call dcopy(nbf*nbf,sm1,1,rotb,1)
c             call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb,nbf,r_u,nbf,
c    +                   0.0d0,r_ofb2,nbf)
c             call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,
c    +                            r_noa2,r_nob2,r_ofa2,r_ofb2,
c    +                            ehfa,ehfb,ewfa,ewfb,enta,entb,power,
c    +                            occa,occb,En0,.false.)
c             occdmax = max(occdmax,occ_diff)
c             call dcopy(nbf*nbf,r_ofb2,1,r_ofb,1)
            enddo
          enddo
          call dcopy(nbf*nbf,r_noa,1,r_noa2,1)
          call dcopy(nbf*nbf,r_nob,1,r_nob2,1)
          call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofa(1,1,ist),nbf,rota,
     +               nbf,0.0d0,r_ofa2,nbf)
          call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_ofb(1,1,ist),nbf,rotb,
     +               nbf,0.0d0,r_ofb2,nbf)
          call wfn1_engradocc(rtdb,nbf,nea,neb,h1,eri,erix,ov,factor,
     +                        ofac,
     +                        r_noa2,r_nob2,r_ofa2,r_ofb2,
     +                        ehfa,ehfb,ewfa,ewfb,enta,entb,power,
     +                        occa,occb,En0(ist),dv0,oerr,.true.)
          occdmax = max(occdmax,dv0)
          call dcopy(nbf*nbf,r_ofa2,1,r_ofa(1,1,ist),1)
          call dcopy(nbf*nbf,r_ofb2,1,r_ofb(1,1,ist),1)
c         factor = factor*(1.0d0 + max(0.00d0,occdmax)/factor)
c         ofac   = ofac  *(1.0d0 + max(0.00d0,oerr)/ofac)
          factor = min(factor,100.0d0*nbf)
          ofac   = min(ofac,10.0d0*nbf)
cDEBUG
c         write(*,*)
          do is = 1, nbf
            write(*,'(" ent a,b: ",3f20.8,3x,3f20.8)')
     +            enta(is),ehfa(is),ewfa(is),
     +            entb(is),ehfb(is),ewfb(is)
          enddo
          call wfn1_print_occ(nbf,nea,neb,r_ofa(1,1,ist),r_ofb(1,1,ist),
     +                        occa,occb,Em1(ist),r_ofa2)
cDEBUG
        enddo
        write(luout,'(1x,i5,i3," iter ",5f20.8,3f10.3)')
     +        iteration,ist,en0(ist),dnorm,occdmax,oerr,tol,
     +        damp,factor,ofac
        write(luout,*)"========================================="
        write(luout,*)"Occupation number optimization converged."
        write(luout,*)"========================================="
        call wfn1_print_occ(nbf,nea,neb,r_ofa(1,1,ist),r_ofb(1,1,ist),
     +                      occa,occb,En0(ist),r_ofa2)
      enddo ! ist
      
      end
C>
C>

C> @}
