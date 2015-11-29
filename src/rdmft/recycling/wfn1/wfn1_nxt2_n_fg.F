C> \ingroup wfn1_nxt2
C> @{
C>
C> \brief Take a Fock matrix and construct the gradient matrices
C>
      subroutine wfn1_nxt2_n_fg(nbf,nex,ox,nx,cx,fx,fo,fn,go,gn,occx,
     +                          ex,tx,Ef,levelshift)
      implicit none
c
#include "wfn1_nxt.fh"
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
                                   !< current spin channel in natural
                                   !< orbital basis
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
      double precision Ef          !< [Output] Electronic energy as
                                   !< calculated from trace(F,D)
      double precision levelshift  !< [In/Output] The levelshifter
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
c     call dgemm('n','n',nbf,nbf,nbf,1.0d0,fx,nbf,nx,nbf,0.0d0,tx,nbf)
c     call dgemm('t','n',nbf,nbf,nbf,1.0d0,nx,nbf,tx,nbf,0.0d0,fn,nbf)
      call dcopy(nbf*nbf,fx,1,fn,1)
      if (levelshift.ge.0.0d0) then
        do ii = 2, nbf
cDEBUG
c         write(*,'("lvlshft occx:",i4,2f14.6,3x,3f14.6)')ii,
c    +          occx(ii-1),occx(ii),
c    +          fn(ii-1,ii-1),fn(ii,ii),fn(ii-1,ii-1)-fn(ii,ii)
cDEBUG
          if (occx(ii-1).ge.occx(ii)) then
            levelshift = max(levelshift,fn(ii-1,ii-1)-fn(ii,ii))
          endif
        enddo
      endif
c
c     Now generate the Fock matrix in occupation function basis
c
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,fn,nbf,ox,nbf,0.0d0,tx,nbf)
      call dgemm('t','n',nbf,nbf,nbf,1.0d0,ox,nbf,tx,nbf,0.0d0,fo,nbf)
      if (levelshift.ge.0.0d0) then
        do ir = 2, nbf
cDEBUG
c         write(*,'("lvlshft   cx:",i4,2f14.6,3x,3f14.6)')ir,
c    +          cx(ir-1),cx(ir),
c    +          fo(ir-1,ir-1),fo(ir,ir),fo(ir-1,ir-1)-fo(ir,ir)
cDEBUG
          if (cx(ir-1).ge.cx(ir)) then
            levelshift = max(levelshift,fo(ir-1,ir-1)-fo(ir,ir))
          endif
        enddo
      endif
cDEBUG
c     write(*,*)'*** Fock diagonal occ, nat'
c     do ii = 1, nbf
c       write(*,*)ii,fo(ii,ii),fn(ii,ii)
c     enddo
cDEBUG
c
c     Extract the expectation values of the natural orbital energies
c
      do ii = 1, nbf
        ex(ii) = fn(ii,ii)
      enddo
c
c     Scale the columns with the occupation numbers
c
      if (wfn1_nxt_sc) then
        do ii = 1, nbf
          call dscal(nbf,occx(ii),fn(1,ii),1)
        enddo
      endif
c
c     Symmetrize the matrix
c
c     do jj = 1, nbf
c       do ii = jj+1, nbf
c         tt = 0.5d0*(fn(ii,jj)+fn(jj,ii))
c         fn(ii,jj) = tt
c         fn(jj,ii) = tt
c       enddo
c       Ef = Ef + fn(jj,jj)
c     enddo
c
c     Copy the result to the gradient
c
c     call dcopy(nbf*nbf,fn,1,gn,1)
      if (wfn1_nxt_sc) then
        do jj = 1, nbf
          do ii = 1, nbf
            gn(ii,jj) = fn(ii,jj)-fn(jj,ii)
            if (ii.lt.jj) gn(ii,jj)=-gn(ii,jj) ! wfn1_nxt2_sym_2_skw
          enddo
        enddo
      else
        call dcopy(nbf*nbf,fn,1,gn,1)
      endif
c
c     Compute the gradient for the occupation functions
c
      do iq = 1, nbf
        do ir = 1, nbf
          fo(ir,iq) = wfn1_ddot3(nbf,ex,1,ox(1,ir),1,ox(1,iq),1)
        enddo
      enddo
c
c     Scale the columns with the occupation numbers
c
      if (wfn1_nxt_sc) then
        do ir = 1, nbf
cDEBUG
c         write(*,*)'fg occ: ',ir,cx(ir)
cDEBUG
          call dscal(nbf,cx(ir),fo(1,ir),1)
        enddo
      endif
c
c     Symmetrize the matrix
c
c     do iq = 1, nbf
c       do ir = iq+1, nbf
c         tt = 0.5d0*(fo(ir,iq)+fo(iq,ir))
c         fo(ir,iq) = tt
c         fo(iq,ir) = tt
c       enddo
c     enddo
c
c     Copy the result to the gradient and zero the virtual-virtual
c     block (the virtual orbitals do not contribute to the energy)
c
c     call dcopy(nbf*nbf,fo,1,go,1)
c     do iq = 1, nbf
c       if (cx(iq).le.0.0d0) then
c         do ir = 1, nbf
c           if (cx(ir).le.0.0d0) then
c             go(ir,iq) = 0.0d0
c           endif
c         enddo
c       else if (cx(iq).ge.1.0d0) then
c         do ir = 1, nbf
c           if (cx(ir).ge.1.0d0) then
c             go(ir,iq) = 0.0d0
c           endif
c         enddo
c       endif
c     enddo
      if (wfn1_nxt_sc) then
        do iq = 1, nbf
          do ir = 1, nbf
            go(ir,iq) = fo(ir,iq)-fo(iq,ir)
            if (ir.lt.iq) go(ir,iq)=-go(ir,iq) ! wfn1_nxt2_sym_2_skw
          enddo
        enddo
      else
        call dcopy(nbf*nbf,fo,1,go,1)
      endif
cDEBUG
c     write(*,*)'*** occ-scaled Fock diagonal occ, nat'
c     do ii = 1, nbf
c       write(*,*)ii,fo(ii,ii),fn(ii,ii)
c     enddo
cDEBUG
c
      end
C> @}
