C> \ingroup wfn1
C> @{
C>
C> \brief Compute the deviation from the true occupation numbers
C>
C> Formulate a penalty function of the occupation that is non-negative
C> and zero only when the total current occupation numbers match 
C> the given occupation numbers.
C>
      subroutine wfn1_occ(nbf,nea,neb,r_ofa,r_ofb,occa,occb,
     +                    factor,En,deviation,v1)
      implicit none
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision r_ofa(nbf,nbf) !< [Output] The alpha electron
      !< occupation functions
      double precision r_ofb(nbf,nbf) !< [Output] The beta electron
      !< occupation functions
      double precision occa(nbf) !< [Input] The actual occupation
      !< numbers
      double precision occb(nbf) !< [Input] The actual occupation
      !< numbers
      double precision factor !< [Input] Scale factor for the penalty
      !< function
c
      double precision En !< [Output] The value of the functional
      double precision deviation !< [Output] Distance to the desired
      !< outcome
c
      double precision v1(nbf) !< [Scratch]
c
c     Local
c
      integer ip !< Counter over extended orbitals
      integer ii !< Counter over natural orbitals
c
      En = 0.0d0
      do ii = 1, nbf
        v1(ii) = 0.0d0
        do ip = 1, nea
          v1(ii) = v1(ii) + r_ofa(ii,ip)*r_ofa(ii,ip)
        enddo
      enddo
      do ii = 1, nbf
        En = En + (v1(ii)-occa(ii))**2
      enddo
c
      do ii = 1, nbf
        v1(ii) = 0.0d0
        do ip = 1, neb
          v1(ii) = v1(ii) + r_ofb(ii,ip)*r_ofb(ii,ip)
        enddo
      enddo
      do ii = 1, nbf
        En = En + (v1(ii)-occb(ii))**2
      enddo
      deviation = sqrt(En)
      En = En * factor
c
      end
C>
C> \brief Print the deviation from the true occupation numbers
C>
C> Print the actual occupation numbers as well as the occupation
C> numbers from the wavefunction.
C>
      subroutine wfn1_print_occ(nbf,nea,neb,r_ofa,r_ofb,occa,occb,En,v1)
      implicit none
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
c
      double precision r_ofa(nbf,nbf) !< [Output] The alpha electron
      !< occupation functions
      double precision r_ofb(nbf,nbf) !< [Output] The beta electron
      !< occupation functions
      double precision occa(nbf) !< [Input] The actual occupation
      !< numbers
      double precision occb(nbf) !< [Input] The actual occupation
      !< numbers
c
      double precision En !< [Output] The value of the functional
c
      double precision v1(nbf) !< [Scratch]
c
c     Local
c
      integer ip !< Counter over extended orbitals
      integer ii !< Counter over natural orbitals
c
c     En = 0.0d0
      do ii = 1, nbf
        v1(ii) = 0.0d0
        do ip = 1, nea
          v1(ii) = v1(ii) + r_ofa(ii,ip)*r_ofa(ii,ip)
        enddo
      enddo
      write(*,*)'Alpha electron occupation numbers:'
      do ii = 1, nbf
        write(*,'(1x,i5,2(2x,f12.4))')ii,occa(ii),v1(ii)
c       write(*,'(1x,i5,1(2x,f12.4))')ii,v1(ii)
        En = En + (v1(ii)-occa(ii))**2
      enddo
c
      do ii = 1, nbf
        v1(ii) = 0.0d0
        do ip = 1, neb
          v1(ii) = v1(ii) + r_ofb(ii,ip)*r_ofb(ii,ip)
        enddo
      enddo
      write(*,*)'Beta electron occupation numbers:'
      do ii = 1, nbf
        write(*,'(1x,i5,2(2x,f12.4))')ii,occb(ii),v1(ii)
c       write(*,'(1x,i5,1(2x,f12.4))')ii,v1(ii)
        En = En + (v1(ii)-occb(ii))**2
      enddo
c
      end
C> @}
