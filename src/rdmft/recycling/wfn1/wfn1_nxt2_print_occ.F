C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Print the occupation numbers from a next step format
C> wavefunction 
C>
C> Print the occupation number of the wavefunction
C>
      subroutine wfn1_nxt2_print_occ(nea,neb,state,nbf,occa,occb)
      implicit none
C
#include "wfn1_nxt.fh"
C
      integer nea !< [Input] The number of alpha electrons
      integer neb !< [Input] The number of beta electrons
      integer nbf !< [Input] The number of basis functions
      double precision state(8*nbf*nbf+2*nbf) !< [Input] The wavefunction
      double precision occa(nbf) !< [Scratch] The alpha occupation
                                 !< numbers
      double precision occb(nbf) !< [Scratch] The beta occupation
                                 !< numbers
c
      integer ii, ir
c
      call dfill(nbf,0.0d0,occa,1)
      call dfill(nbf,0.0d0,occb,1)
      do ii = 1, nbf
        do ir = 1, nbf
          occa(ii) = occa(ii)
     &             + state(wfn1_nxt_aconf()-1+ir)
     &             * state(wfn1_nxt_aocc()-1+ii+(ir-1)*nbf)**2
        enddo
        do ir = 1, nbf
          occb(ii) = occb(ii)
     &             + state(wfn1_nxt_bconf()-1+ir)
     &             * state(wfn1_nxt_bocc()-1+ii+(ir-1)*nbf)**2
        enddo
      enddo
c
      write(*,*)
      write(*,*)'alpha, beta occupation numbers'
      write(*,*)'------------------------------'
      do ii = 1, nbf
        write(*,'(i4,2f14.8)')ii,occa(ii),occb(ii)
      enddo
      write(*,*)
C
      end
C> @}
