C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Calculate the 2-electron Coulomb energy and Fock matrix
C> contributions
C>
C> Compute the 2-electron Coulomb energy contribution as
C> \f{eqnarray*}{
C>    E_{2C} &=& \sum_{\sigma\sigma'}\sum_{ab}D^\sigma_{ab}
C>               \left(\sum_{cd}(ab|cd)D^{\sigma'}_{cd}\right)
C> \f}
C> I.e. the summation over \f$c, d\f$ can be thought of as a
C> matrix-vector product. The final summation over \f$a, b\f$ is then a
C> dot-product.
C>
C> The gradient of this quantity can be obtained by directly
C> differentiating the density matrix expressions wrt. the occupation
C> function and natural orbital coefficients (see wfn1_next_step.dox
C> for details).
C>
C> The Fock matrix contributions are calculated in the natural 
C> orbital basis.
C>
      subroutine wfn1_nxt2_2el_n_C_ef(nbf,nea,neb,eri,occa,occb,da,db,
     +           oa,na,ca,ob,nb,cb,fa,fb,E2elC,Ef2elC,ta,tb,t1,t2)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The 2-electron
                                            !< Coulomb integrals
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
      double precision E2elC  !< [Output] The Coulomb energy
      double precision Ef2elC !< [Output] The Coulomb energy from the
                              !< Fock matrix
c
      double precision ta(nbf,nbf) !< [Scratch] The alpha Fock like
                                   !< matrix
      double precision tb(nbf,nbf) !< [Scratch] The beta Fock like
                                   !< matrix
      double precision t1(nbf,nbf) !< [Scratch] 
      double precision t2(nbf,nbf) !< [Scratch]
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
c
c     Contract the density matrices with the 2-electron integrals
c     to generate Fock-like quantities
c
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,da,1,
     +           0.0d0,ta,1)
      call dgemv('N',nbf*nbf,nbf*nbf,1.0d0,eri,nbf*nbf,db,1,
     +           0.0d0,tb,1)
c
c     Calculate the energy contributions
c
      E2elC  = 0.0d0
      Ef2elC = 0.0d0
      E2elC  = E2elC + ddot(nbf*nbf,da,1,ta,1)/2.0d0
      E2elC  = E2elC + ddot(nbf*nbf,da,1,tb,1)/2.0d0
      E2elC  = E2elC + ddot(nbf*nbf,db,1,ta,1)/2.0d0
      E2elC  = E2elC + ddot(nbf*nbf,db,1,tb,1)/2.0d0
c
c     Add contributions to the Fock matrices
c
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,na,nbf,ta,nbf,0.0d0,
     +           t1,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,t1,nbf,na,nbf,0.0d0,
     +           t2,nbf)
      call daxpy(nbf*nbf,1.0d0,t2,1,fa,1)
      do ii = 1, nbf
        Ef2elC = Ef2elC + t2(ii,ii)*occa(ii)
      enddo
c
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,nb,nbf,ta,nbf,0.0d0,
     +           t1,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,t1,nbf,nb,nbf,0.0d0,
     +           t2,nbf)
      call daxpy(nbf*nbf,1.0d0,t2,1,fb,1)
      do ii = 1, nbf
        Ef2elC = Ef2elC + t2(ii,ii)*occb(ii)
      enddo
c
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,nb,nbf,tb,nbf,0.0d0,
     +           t1,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,t1,nbf,nb,nbf,0.0d0,
     +           t2,nbf)
      call daxpy(nbf*nbf,1.0d0,t2,1,fb,1)
      do ii = 1, nbf
        Ef2elC = Ef2elC + t2(ii,ii)*occb(ii)
      enddo
c
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,na,nbf,tb,nbf,0.0d0,
     +           t1,nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,t1,nbf,na,nbf,0.0d0,
     +           t2,nbf)
      call daxpy(nbf*nbf,1.0d0,t2,1,fa,1)
      do ii = 1, nbf
        Ef2elC = Ef2elC + t2(ii,ii)*occa(ii)
      enddo
c
      end
C>
C> @}
