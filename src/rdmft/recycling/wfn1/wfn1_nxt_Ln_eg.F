C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the natural orbital Lagrangian energy and
C> gradient contributions
C>
C> Compute the natural orbital Lagrangian energy contribution as
C> \f{eqnarray*}{
C>   L_2 &=& \sum_\sigma\sum_{ij}\lambda^{N\sigma}_{ij}\left(
C>           \sum_{ab} N^\sigma_{ai}S_{ab}N^\sigma_{bj} - \delta_{ij}
C>           \right)
C> \f}
C> and the corresponding gradient contributions to the natural
C> orbital and Lagrangian components. See wfn1_next_step.dox for
C> details.
C>
      subroutine wfn1_nxt_Ln_eg(nbf,ov,na,nb,lna,lnb,
     +           dna,dnb,dlna,dlnb,ELn,tv)
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
      double precision lna(nbf,nbf) !< [Input] The alpha natural
                                    !< orbitals Lagrangian multipliers
      double precision lnb(nbf,nbf) !< [Input] The beta natural
                                    !< orbitals Lagrangian multipliers
c
      double precision dna(nbf,nbf)  !< [In/Output] The alpha natural
                                     !< orbitals gradient
      double precision dnb(nbf,nbf)  !< [In/Output] The beta natural
                                     !< orbitals gradient
      double precision dlna(nbf,nbf) !< [In/Output] The alpha natural
                                     !< orbitals Lagrangian gradient
      double precision dlnb(nbf,nbf) !< [In/Output] The beta natural
                                     !< orbitals Lagrangian gradient
c
      double precision ELn !< [Output] The natural orbitals
                           !< Lagrangian energy
c
      double precision tv(nbf) !< [Scratch] A temporary vector
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
      ELn = 0.0d0
      do ij = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,na(1,ij),1,0.0d0,tv,1)
        do ii = 1, nbf
          tt = ddot(nbf,na(1,ii),1,tv,1)
          if (ii.eq.ij) tt = tt - 1.0d0
          ELn = ELn + lna(ii,ij)*tt
          dlna(ii,ij) = dlna(ii,ij) + tt
        enddo
      enddo
      do ij = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,nb(1,ij),1,0.0d0,tv,1)
        do ii = 1, nbf
          tt = ddot(nbf,nb(1,ii),1,tv,1)
          if (ii.eq.ij) tt = tt - 1.0d0
          ELn = ELn + lnb(ii,ij)*tt
          dlnb(ii,ij) = dlnb(ii,ij) + tt
        enddo
      enddo
c
c     Compute derivatives wrt occupation function coefficients
c
      do ii = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,na(1,ii),1,0.0d0,tv,1)
        do ij = 1, nbf
          do ia = 1, nbf
            dna(ia,ij) = dna(ia,ij) + tv(ia)*(lna(ii,ij)+lna(ij,ii))
          enddo
        enddo
      enddo
      do ii = 1, nbf
        call dgemv('N',nbf,nbf,1.0d0,ov,nbf,nb(1,ii),1,0.0d0,tv,1)
        do ij = 1, nbf
          do ia = 1, nbf
            dnb(ia,ij) = dnb(ia,ij) + tv(ia)*(lnb(ii,ij)+lnb(ij,ii))
          enddo
        enddo
      enddo
c
      end
C>
C> @}
