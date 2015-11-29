C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Calculate the 2-electron exchange energy and Fock matrix
C> contributions
C>
C> Compute the 2-electron exchange energy contributions as
C> \f{eqnarray*}{
C>    E_{2X} &=& \sum_{\sigma}\sum_{ab}F\left(D^\sigma\right)_{ab}
C>               \left(\sum_{cd}(ad|cb)F\left(D^\sigma\right)_{cd}
C>               \right)
C> \f}
C> I.e. after storing the integrals in a particular way the summation
C> over \f$c, d\f$ can be thought of as a matrix-vector product. The 
C> final summation over \f$a, b\f$ is then a dot-product.
C>
C> The density matrices in this expression are subject to functions.
C> Typically these functions involve taking the input to a low
C> power, e.g. \f$F(X) = \sqrt{X}\f$, or other powers close to 
C> \f$1/2\f$.
C>
C> The gradients of this quantity can be obtained by directly
C> differentiating the density matrix functions wrt. the occupation
C> function and the natural orbital coefficients (see 
C> wfn1_next_step.dox for details).
C>
C> Here we want to calculate the Fock matrix. The Fock matrix can be
C> defined in multiple ways:
C> \f{eqnarray*}{
C>    F_{ab} &=& \frac{\partial E_{2X}}{\partial D_{ab}} 
C> \f}
C> where we differentiate the energy wrt. a density matrix element. We
C> repeat this for all atomic orbitals \f$a\f$ and \f$b\f$. 
C> Alternatively we can use
C> \f{eqnarray*}{
C>   E &=& tr(F,D) 
C> \f}
C> as the definition of the Fock matrix given that the energy expression
C> is known and so is the density matrix.
C>
C> Interestingly the first definition is problematic in the case where
C> the exchange energy is considered in terms of functions of the
C> density matrix. E.g. for Hartree-Fock I could consider the density
C> as
C> \f{eqnarray*}{
C>   D^X &=& D^1
C> \f}
C> In Hartree-Fock theory we use that
C> \f{eqnarray*}{
C>   \frac{\partial D_{ab}}{\partial D_{ab}} &=& 1
C> \f}
C> which is trivial if we consider simply two matrix elements. If I 
C> consider this equation to be equivalent to 
C> \f{eqnarray*}{
C>   \frac{\partial D^1_{ab}}{\partial D_{ab}} &=& 1
C> \f}
C> one could hope that this means
C> \f{eqnarray*}{
C>   \partial D^0_{ab} &=& 1
C> \f}
C> but given that \f$D\f$ is defined as
C> \f{eqnarray*}{
C>   D_{ab} &=& \sum_{i=1}^{N_{basis}} N_{ai}d_{i}N_{bi}
C> \f}
C> this is clearly not true as
C> \f{eqnarray*}{
C>   D^0_{ab} &=& \sum_{i=1}^{N_{basis}} N_{ai}(d_{i})^0N_{bi} \\\\
C>            &=& \sum_{i=1}^{N_{basis}} N_{ai}N_{bi} \\\\
C> \f}
C> I.e. to have any hope of success one would at least have to define
C> \f{eqnarray*}{
C>   \frac{\partial F(D)_{ab}}{\partial D_{ab}}
C>   &=& \frac{F'(D)_{ab}}{D^0_{ab}}
C> \f}
C>
      subroutine wfn1_nxt2_2el_X_ef(nbf,nea,neb,pow,erix,ov,ov12,ovm12,
     +           occa,occb,da,db,oa,na,ca,ob,nb,cb,fa,fb,E2elX,
     +           ta,tb,t4)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision pow !< [Input] The power of the density matrix
                           !< function
c
      double precision ov(nbf,nbf) !< [Input] The overlap matrix
      double precision ov12(nbf,nbf) !< [Input] The sqrt(overlap)
      double precision ovm12(nbf,nbf) !< [Input] The 1/sqrt(overlap)
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                             !< exchange integrals
c
      double precision occa(nbf)   !< [Scratch] The alpha occupation
                                   !< numbers
      double precision occb(nbf)   !< [Scratch] The beta occupation
                                   !< numbers
      double precision da(nbf,nbf) !< [Scratch] The alpha density
                                   !< matrix
      double precision db(nbf,nbf) !< [Scratch] The beta density
                                   !< matrix
c
      double precision oa(nbf,nbf) !< [Input] The alpha occupation
                                   !< functions
      double precision ob(nbf,nbf) !< [Input] The beta occupation
                                   !< functions
      double precision na(nbf,nbf) !< [Input] The alpha natural
                                   !< orbitals
      double precision nb(nbf,nbf) !< [Input] The beta natural
                                   !< orbitals
      double precision ca(nbf)     !< [Input] The alpha electron
                                   !< configuration
      double precision cb(nbf)     !< [Input] The beta electron
                                   !< configuration
c
      double precision fa(nbf,nbf) !< [In/Output] The alpha occupation
                                   !< functions gradient
      double precision fb(nbf,nbf) !< [In/Output] The beta occupation
                                   !< functions gradient
c
      double precision E2elX !< [Output] The Coulomb energy
c
      double precision ta(nbf,nbf) !< [Scratch] The alpha Fock like
                                   !< matrix
      double precision tb(nbf,nbf) !< [Scratch] The beta Fock like
                                   !< matrix
      double precision t4(nbf,nbf,nbf,nbf) !< [Scratch] Derivative
                                   !< tensor
c
      double precision contract
c
      integer ir, ip !< Counters over occupation functions
      integer ii, ij !< Counters over natural orbitals
      integer ia, ib !< Counters over basis functions
      integer ic, id !< Counters over basis functions
c
      double precision ddot
      external         ddot
c
c     Code
c
      call dfill(nbf,0.0d0,occa,1)
      call dfill(nbf,0.0d0,occb,1)
      do ir = 1, nbf
        do ii = 1, nbf
          occa(ii) = occa(ii) + oa(ii,ir)*oa(ii,ir)*ca(ir)
        enddo
      enddo
      do ir = 1, nbf
        do ii = 1, nbf
          occb(ii) = occb(ii) + ob(ii,ir)*ob(ii,ir)*cb(ir)
        enddo
      enddo
      do ii = 1, nbf
        occa(ii) = occa(ii)**pow
        occb(ii) = occb(ii)**pow
      enddo
      call dfill(nbf*nbf,0.0d0,da,1)
      call dfill(nbf*nbf,0.0d0,db,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = da(ia,ib)
     +                + na(ia,ii)*na(ib,ii)*occa(ii)
            db(ia,ib) = db(ia,ib)
     +                + nb(ia,ii)*nb(ib,ii)*occb(ii)
          enddo
        enddo
      enddo
c
c     Contract the density matrix functions with the 2-electron
c     integrals to generate Fock-like quatities
c
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,erix,nbf*nbf,da,1,
     +           0.0d0,ta,1)
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,erix,nbf*nbf,db,1,
     +           0.0d0,tb,1)
c
c     Calculate the energy contributions
c
      E2elX = 0.0d0
      E2elX = E2elX - ddot(nbf*nbf,da,1,ta,1)/2.0d0
      E2elX = E2elX - ddot(nbf*nbf,db,1,tb,1)/2.0d0
c
c     Attempt 2 use automatic differentiation to generate matrix
c     function gradient
c
      call dfill(nbf,0.0d0,occa,1)
      call dfill(nbf,0.0d0,occb,1)
      do ir = 1, nbf
        do ii = 1, nbf
          occa(ii) = occa(ii) + oa(ii,ir)*oa(ii,ir)*ca(ir)
        enddo
      enddo
      do ir = 1, nbf
        do ii = 1, nbf
          occb(ii) = occb(ii) + ob(ii,ir)*ob(ii,ir)*cb(ir)
        enddo
      enddo
      call dfill(nbf*nbf,0.0d0,da,1)
      call dfill(nbf*nbf,0.0d0,db,1)
      do ib = 1, nbf
        do ia = 1, nbf
          do ii = 1, nbf
            da(ia,ib) = da(ia,ib)
     +                + na(ia,ii)*na(ib,ii)*occa(ii)
            db(ia,ib) = db(ia,ib)
     +                + nb(ia,ii)*nb(ib,ii)*occb(ii)
          enddo
        enddo
      enddo
      call wfn1_d_matfunc(nbf,da,ov,ov12,ovm12,1,t4)
      do id = 1, nbf
        do ic = 1, nbf
          contract = 0.0d0
          do ib = 1, nbf
            do ia = 1, nbf
              contract = contract + ta(ia,ib)*t4(ia,ib,ic,id)
            enddo
          enddo
          fa(ic,id) = fa(ic,id) - contract
        enddo
      enddo
      call wfn1_d_matfunc(nbf,db,ov,ov12,ovm12,1,t4)
      do id = 1, nbf
        do ic = 1, nbf
          contract = 0.0d0
          do ib = 1, nbf
            do ia = 1, nbf
              contract = contract + tb(ia,ib)*t4(ia,ib,ic,id)
            enddo
          enddo
          fb(ic,id) = fb(ic,id) - contract
        enddo
      enddo
c
      end
C>
C> @}
