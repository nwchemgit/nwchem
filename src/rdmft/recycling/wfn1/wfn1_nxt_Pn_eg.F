C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the natural orbital penalty energy and
C> gradient contributions
C>
C> Compute the natural orbital penalty energy contribution as
C> \f{eqnarray*}{
C>   L_2 &=& \sum_\sigma\sum_{ij}\lambda^{N\sigma}_{ij}\left(
C>           \sum_{ab} N^\sigma_{ai}S_{ab}N^\sigma_{bj} - \delta_{ij}
C>           \right)^2
C> \f}
C> and the corresponding gradient contributions to the natural
C> orbital and penalty coefficient components. See wfn1_next_step.dox
C> for details.
C>
      subroutine wfn1_nxt_Pn_eg(nbf,ov,na,nb,pna,pnb,
     +           dna,dnb,dpna,dpnb,EPn,tv,ta,tb)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision ov(nbf,nbf)  !< [Input] The overlap matrix
c
      double precision na(nbf,nbf)  !< [Input] The alpha natural
                                    !< orbitals
      double precision nb(nbf,nbf)  !< [Input] The beta natural
                                    !< orbitals
      double precision pna(nbf,nbf) !< [Input] The alpha natural
                                    !< orbitals penalty coefficients
      double precision pnb(nbf,nbf) !< [Input] The beta natural
                                    !< orbitals penalty coefficients
c
      double precision dna(nbf,nbf)  !< [In/Output] The alpha natural
                                     !< orbitals gradient
      double precision dnb(nbf,nbf)  !< [In/Output] The beta natural
                                     !< orbitals gradient
      double precision dpna(nbf,nbf) !< [In/Output] The alpha natural
                                     !< orbitals penalty coefficient
                                     !< gradient
      double precision dpnb(nbf,nbf) !< [In/Output] The beta natural
                                     !< orbitals penalty coefficient
                                     !< gradient
c
      double precision EPn !< [Output] The natural orbitals
                           !< penalty energy
c
      double precision tv(nbf) !< [Scratch] A temporary vector
      double precision ta(nbf,nbf) !< [Scratch]
      double precision tb(nbf,nbf) !< [Scratch]
c
      double precision tt
c
      double precision ddot
      external         ddot
c
      integer ii, ij
      integer ia, ib
c
c     Compute the Lagrangian contributions
c
      EPn = 0.0d0
      do ij = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,na(1,ij),1,0.0d0,tv,1)
        do ii = 1, nbf
          tt = ddot(nbf,na(1,ii),1,tv,1)
          if (ii.eq.ij) tt = tt - 1.0d0
          ta(ii,ij) = tt
          EPn = EPn + pna(ii,ij)*(tt**2)
          dpna(ii,ij) = dpna(ii,ij) + tt**2
        enddo
      enddo
      do ij = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,nb(1,ij),1,0.0d0,tv,1)
        do ii = 1, nbf
          tt = ddot(nbf,nb(1,ii),1,tv,1)
          if (ii.eq.ij) tt = tt - 1.0d0
          tb(ii,ij) = tt
          EPn = EPn + pnb(ii,ij)*(tt**2)
          dpnb(ii,ij) = dpnb(ii,ij) + tt**2
        enddo
      enddo
c
c     Compute derivatives wrt occupation function coefficients
c
      do ij = 1, nbf
        do ii = 1, nbf
          ta(ii,ij) = 2*ta(ii,ij)*pna(ii,ij)
          tb(ii,ij) = 2*tb(ii,ij)*pnb(ii,ij)
        enddo
      enddo
      do ii = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,na(1,ii),1,0.0d0,tv,1)
        do ij = 1, nbf
          do ia = 1, nbf
            dna(ia,ij) = dna(ia,ij) + tv(ia)*(ta(ii,ij)+ta(ij,ii))
          enddo
        enddo
      enddo
      do ii = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,nb(1,ii),1,0.0d0,tv,1)
        do ij = 1, nbf
          do ia = 1, nbf
            dnb(ia,ij) = dnb(ia,ij) + tv(ia)*(tb(ii,ij)+tb(ij,ii))
          enddo
        enddo
      enddo
c
      end
C>
C> @}
