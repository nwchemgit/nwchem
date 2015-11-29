C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Take a Fock matrix and construct the gradient matrices
C>
      subroutine wfn1_nxt2_fg(nbf,nex,ox,nx,cx,fx,fo,fn,go,gn,occx,
     +                        ex,tx)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nex !< [Input] The number of electrons in this spin
                  !< channel
c
      double precision ox(nbf,nbf) !< [Input] The occupation functions
                                   !< for the current spin channel
      double precision nx(nbf,nbf) !< [Input] The natural orbitals
                                   !< for the current spin channel
      double precision cx(nbf)     !< [Input] The electron configuration
                                   !< for the current spin channel
c
      double precision fx(nbf,nbf) !< [Input] The Fock matrix for the
                                   !< current spin channel
c
      double precision fo(nbf,nbf) !< [Output] The Fock matrix for the
                                   !< occupation functions
      double precision fn(nbf,nbf) !< [Output] The Fock matrix for the
                                   !< natural orbitals
c
      double precision go(nbf,nbf) !< [Output] The gradient for the
                                   !< occupation functions
      double precision gn(nbf,nbf) !< [Output] The gradient for the
                                   !< natural orbitals
c
      double precision occx(nbf)   !< [Scratch] The occupation numbers
      double precision ex(nbf)     !< [Scratch] Expectation values
      double precision tx(nbf,nbf) !< [Scratch] Temporary matrix
c
      integer ii !< Counter
      integer jj !< Counter
      integer iq !< Counter over occupation functions
      integer ir !< Counter over occupation functions
c
      double precision tt !< Temporary
c
      double precision wfn1_ddot3 !< Ddot function with weights
      external         wfn1_ddot3
c
c     First calculate the occupation numbers
c
      call dfill(nbf,0.0d0,occx,1)
      do ir = 1, nbf
        do ii = 1, nbf
          occx(ii) = occx(ii) + ox(ii,ir)*ox(ii,ir)*cx(ir)
        enddo
      enddo
c
c     Transform the Fock matrix from AO basis into the natural orbital
c     basis
c
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,fx,nbf,nx,nbf,0.0d0,tx,nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,nx,nbf,tx,nbf,0.0d0,fn,nbf)
c
c     Now generate the Fock matrix in occupation function basis
c
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,fn,nbf,ox,nbf,0.0d0,tx,nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,ox,nbf,tx,nbf,0.0d0,fo,nbf)
c
c     Extract the expectation values of the natural orbital energies
c
c     do ii = 1, nbf
c       ex(ii) = fn(ii,ii)
c     enddo
c
c     Scale the columns with the occupation numbers
c
      do ii = 1, nbf
        call dscal(nbf,occx(ii),fn(1,ii),1)
      enddo
c
c     Symmetrize the matrix
c
      do jj = 1, nbf
        do ii = jj+1, nbf
          tt = 0.5d0*(fn(ii,jj)+fn(jj,ii))
          fn(ii,jj) = tt
          fn(jj,ii) = tt
        enddo
      enddo
c
c     Copy the result to the gradient
c
      call dcopy(nbf*nbf,fn,1,gn,1)
c
c     Compute the gradient for the occupation functions
c
c     do iq = 1, nbf
c       do ir = 1, nbf
c         fo(ir,iq) = wfn1_ddot3(nbf,ex,1,ox(1,ir),1,ox(1,iq),1)
c       enddo
c     enddo
c
c     Scale the columns with the occupation numbers
c
      do ir = 1, nbf
        call dscal(nbf,cx(ir),fo(1,ir),1)
      enddo
c
c     Symmetrize the matrix
c
      do iq = 1, nbf
        do ir = iq+1, nbf
          tt = 0.5d0*(fo(ir,iq)+fo(iq,ir))
          fo(ir,iq) = tt
          fo(iq,ir) = tt
        enddo
      enddo
c
c     Copy the result to the gradient and zero the virtual-virtual
c     block (the virtual orbitals do not contribute to the energy)
c
      call dcopy(nbf*nbf,fo,1,go,1)
c     do iq = 1, nbf
c       if (cx(iq).le.0.0d0) then
c         do ir = 1, nbf
c           if (cx(ir).le.0.0d0) then
c             go(ir,iq) = 0.0d0
c           endif
c         enddo
c       endif
c     enddo
c
      end
C> @}
