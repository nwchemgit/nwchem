C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the 2-electron exchange energy and gradient
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
      subroutine wfn1_nxt_2el_X_eg(nbf,nea,neb,pow,erix,occa,occb,da,db,
     +           oa,na,ob,nb,doa,dna,dob,dnb,E2elX,fa,fb)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision pow !< [Input] The power of the density matrix
                           !< function
c
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
c
      double precision doa(nbf,nbf) !< [In/Output] The alpha occupation
                                    !< functions gradient
      double precision dob(nbf,nbf) !< [In/Output] The beta occupation
                                    !< functions gradient
      double precision dna(nbf,nbf) !< [In/Output] The alpha natural
                                    !< orbitals gradient
      double precision dnb(nbf,nbf) !< [In/Output] The beta natural
                                    !< orbitals gradient
c
      double precision E2elX !< [Output] The Coulomb energy
c
      double precision fa(nbf,nbf) !< [Scratch] The alpha Fock like
                                   !< matrix
      double precision fb(nbf,nbf) !< [Scratch] The beta Fock like
                                   !< matrix
      double precision fda(nbf)    !< [Scratch] The alpha Fock diagonal
      double precision fdb(nbf)    !< [Scratch] The beta Fock diagonal
c
      integer ir, ip !< Counters over occupation functions
      integer ii, ij !< Counters over natural orbitals
      integer ia, ib !< Counters over basis functions
c
      double precision ddot
      external         ddot
c
c     Code
c
      call dfill(nbf,0.0d0,occa,1)
      call dfill(nbf,0.0d0,occb,1)
      do ir = 1, nea
        do ii = 1, nbf
          occa(ii) = occa(ii) + oa(ii,ir)*oa(ii,ir)
        enddo
      enddo
      do ir = 1, neb
        do ii = 1, nbf
          occb(ii) = occb(ii) + ob(ii,ir)*ob(ii,ir)
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
     +           0.0d0,fa,1)
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,erix,nbf*nbf,db,1,
     +           0.0d0,fb,1)
c
c     Calculate the energy contributions
c
      E2elX = 0.0d0
      E2elX = E2elX - ddot(nbf*nbf,da,1,fa,1)/2.0d0
      E2elX = E2elX - ddot(nbf*nbf,db,1,fb,1)/2.0d0
c
c     Calculate the gradient contributions
c
c     - First calculate the gradients wrt. the natural orbitals
c
      do ii = 1, nbf
        do ia = 1, nbf
          do ib = 1, nbf
            dna(ia,ii) = dna(ia,ii) - 2*fa(ia,ib)*na(ib,ii)*occa(ii)
            dnb(ia,ii) = dnb(ia,ii) - 2*fb(ia,ib)*nb(ib,ii)*occb(ii)
          enddo
        enddo
      enddo
c
c     - Second calculate the gradients wrt. the occupation functions
c
      call dfill(nbf,0.0d0,fda,1)
      call dfill(nbf,0.0d0,fdb,1)
      do ii = 1, nbf
        do ia = 1, nbf
          do ib = 1, nbf
            fda(ii) = fda(ii) + na(ia,ii)*(2*fa(ia,ib))*na(ib,ii)
            fdb(ii) = fdb(ii) + nb(ia,ii)*(2*fb(ia,ib))*nb(ib,ii)
          enddo
        enddo
      enddo
c
      call dfill(nbf,0.0d0,occa,1)
      call dfill(nbf,0.0d0,occb,1)
      do ir = 1, nea
        do ii = 1, nbf
          occa(ii) = occa(ii) + oa(ii,ir)*oa(ii,ir)
        enddo
      enddo
      do ir = 1, neb
        do ii = 1, nbf
          occb(ii) = occb(ii) + ob(ii,ir)*ob(ii,ir)
        enddo
      enddo
      do ii = 1, nbf
        occa(ii) = pow*(occa(ii)**(pow-1.0d0))
        occb(ii) = pow*(occb(ii)**(pow-1.0d0))
      enddo
      do ir = 1, nea
        do ii = 1, nbf
          doa(ii,ir) = doa(ii,ir) - 2*fda(ii)*occa(ii)*oa(ii,ir)
        enddo
      enddo
      do ir = 1, neb
        do ii = 1, nbf
          dob(ii,ir) = dob(ii,ir) - 2*fdb(ii)*occb(ii)*ob(ii,ir)
        enddo
      enddo
c
      end
C>
C> @}
