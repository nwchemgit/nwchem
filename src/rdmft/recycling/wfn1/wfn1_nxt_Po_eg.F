C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the occupation function penalty energy and
C> gradient contributions
C>
C> Compute the occupation function penalty energy contribution as
C> \f{eqnarray*}{
C>   L_1 &=& \sum_\sigma\sum_{pq}\lambda^{O\sigma}_{pq}\left(
C>           \sum_i O^\sigma_{ip}O^\sigma_{iq} - \delta_{pq}\right)^2
C> \f}
C> and the corresponding gradient contributions to the occupation
C> function and penalty coefficient components. See wfn1_next_step.dox
C> for details.
C>
      subroutine wfn1_nxt_Po_eg(nbf,oa,ob,poa,pob,
     +           doa,dob,dpoa,dpob,EPo,ta,tb)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
c
      double precision oa(nbf,nbf)  !< [Input] The alpha occupation
                                    !< functions
      double precision ob(nbf,nbf)  !< [Input] The beta occupation
                                    !< functions
      double precision poa(nbf,nbf) !< [Input] The alpha occupation
                                    !< functions penalty coefficients
      double precision pob(nbf,nbf) !< [Input] The beta occupation
                                    !< functions penalty coefficients
c
      double precision doa(nbf,nbf)  !< [In/Output] The alpha occupation
                                     !< function gradient
      double precision dob(nbf,nbf)  !< [In/Output] The beta occupation
                                     !< function gradient
      double precision dpoa(nbf,nbf) !< [In/Output] The alpha occupation
                                     !< function penalty coefficient
                                     !< gradient
      double precision dpob(nbf,nbf) !< [In/Output] The beta occupation
                                     !< function penalty coefficient
                                     !< gradient
c
      double precision EPo !< [Output] The occupation function
                           !< penalty energy
c
      double precision ta(nbf,nbf) !< [Scratch] 
      double precision tb(nbf,nbf) !< [Scratch] 
c
      double precision tt
c
      double precision ddot
      external         ddot
c
      integer ip, iq
      integer ii
c
c     Compute the penalty function contributions
c
      EPo = 0.0d0
      do ip = 1, nbf
        do iq = 1, nbf
          tt = ddot(nbf,oa(1,ip),1,oa(1,iq),1)
          if (ip.eq.iq) tt = tt - 1.0d0
          ta(ip,iq) = tt
          EPo = EPo + poa(ip,iq)*(tt**2)
          dpoa(ip,iq) = dpoa(ip,iq) + tt**2
        enddo
      enddo
      do ip = 1, nbf
        do iq = 1, nbf
          tt = ddot(nbf,ob(1,ip),1,ob(1,iq),1)
          if (ip.eq.iq) tt = tt - 1.0d0
          tb(ip,iq) = tt
          EPo = EPo + pob(ip,iq)*(tt**2)
          dpob(ip,iq) = dpob(ip,iq) + tt**2
        enddo
      enddo
c
c     Compute derivatives wrt occupation function coefficients
c
      do iq = 1, nbf
        do ip = 1, nbf
          ta(ip,iq) = 2*ta(ip,iq)*poa(ip,iq)
          tb(ip,iq) = 2*tb(ip,iq)*pob(ip,iq)
        enddo
      enddo
      do ip = 1, nbf
        do ii = 1, nbf
          doa(ii,ip) = doa(ii,ip) + ddot(nbf,oa(ii,1),nbf,ta(1,ip),1)
          doa(ii,ip) = doa(ii,ip) + ddot(nbf,oa(ii,1),nbf,ta(ip,1),nbf)
        enddo
      enddo
      do ip = 1, nbf
        do ii = 1, nbf
          dob(ii,ip) = dob(ii,ip) + ddot(nbf,ob(ii,1),nbf,tb(1,ip),1)
          dob(ii,ip) = dob(ii,ip) + ddot(nbf,ob(ii,1),nbf,tb(ip,1),nbf)
        enddo
      enddo
c
      end
C>
C> @}
