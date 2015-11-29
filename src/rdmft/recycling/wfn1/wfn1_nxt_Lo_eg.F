C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the occupation function Lagrangian energy and
C> gradient contributions
C>
C> Compute the occupation function Lagrangian energy contribution as
C> \f{eqnarray*}{
C>   L_1 &=& \sum_\sigma\sum_{pq}\lambda^{O\sigma}_{pq}\left(
C>           \sum_i O^\sigma_{ip}O^\sigma_{iq} - \delta_{pq}\right)
C> \f}
C> and the corresponding gradient contributions to the occupation
C> function and Lagrangian components. See wfn1_next_step.dox for
C> details.
C>
      subroutine wfn1_nxt_Lo_eg(nbf,oa,ob,loa,lob,
     +           doa,dob,dloa,dlob,ELo)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision oa(nbf,nbf)  !< [Input] The alpha occupation
                                    !< functions
      double precision ob(nbf,nbf)  !< [Input] The beta occupation
                                    !< functions
      double precision loa(nbf,nbf) !< [Input] The alpha occupation
                                    !< functions Lagrangian multipliers
      double precision lob(nbf,nbf) !< [Input] The beta occupation
                                    !< functions Lagrangian multipliers
c
      double precision doa(nbf,nbf)  !< [In/Output] The alpha occupation
                                     !< function gradient
      double precision dob(nbf,nbf)  !< [In/Output] The beta occupation
                                     !< function gradient
      double precision dloa(nbf,nbf) !< [In/Output] The alpha occupation
                                     !< function Lagrangian gradient
      double precision dlob(nbf,nbf) !< [In/Output] The beta occupation
                                     !< function Lagrangian gradient
c
      double precision ELo !< [Output] The occupation function
                           !< Lagrangian energy
c
      double precision tt
c
      double precision ddot
      external         ddot
c
      integer ip, iq
      integer ii
c
c     Compute the Lagrangian contributions
c
      ELo = 0.0d0
      do ip = 1, nbf
        do iq = 1, nbf
          tt = ddot(nbf,oa(1,ip),1,oa(1,iq),1)
          if (ip.eq.iq) tt = tt - 1.0d0
          ELo = ELo + loa(ip,iq)*tt
          dloa(ip,iq) = dloa(ip,iq) + tt
        enddo
      enddo
      do ip = 1, nbf
        do iq = 1, nbf
          tt = ddot(nbf,ob(1,ip),1,ob(1,iq),1)
          if (ip.eq.iq) tt = tt - 1.0d0
          ELo = ELo + lob(ip,iq)*tt
          dlob(ip,iq) = dlob(ip,iq) + tt
        enddo
      enddo
c
c     Compute derivatives wrt occupation function coefficients
c
      do ip = 1, nbf
        do ii = 1, nbf
          doa(ii,ip) = doa(ii,ip) + ddot(nbf,oa(ii,1),nbf,loa(1,ip),1)
          doa(ii,ip) = doa(ii,ip) + ddot(nbf,oa(ii,1),nbf,loa(ip,1),nbf)
        enddo
      enddo
      do ip = 1, nbf
        do ii = 1, nbf
          dob(ii,ip) = dob(ii,ip) + ddot(nbf,ob(ii,1),nbf,lob(1,ip),1)
          dob(ii,ip) = dob(ii,ip) + ddot(nbf,ob(ii,1),nbf,lob(ip,1),nbf)
        enddo
      enddo
c
      end
C>
C> @}
