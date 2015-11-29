C> \ingroup wfn1
C> @{
C>
C> \brief Evaluate the penalty function that ensures the proper
C> orbital ordering is maintained
C>
C> When the optimization is given full freedom there are several ways
C> in which the "free energy" contribution can be significantly 
C> inflated. This can happen by swapping occupied orbitals (the orbital
C> with the negative temperature will correspond to exactly 1 natural
C> orbital so the entropy is 0, the other orbital will acquire a
C> massively high temperature, combining to give an unphysically
C> large "free energy"), or by swapping occupied and virtual orbitals
C> (again the occupied orbital will acquire a spurious massively high
C> temperature). 
C>
C> To block the calculation from pursuing these spurious solutions 
C> a penalty function is introduced. This function ensures that 
C> - all occupied orbital have non-negative temperatures
C> - the orbital energy ordering is maintained 
C> (we cannot enforce that all orbital temperatures be positive as
C> moving away from the Hartree-Fock solution compresses the spectrum
C> of the Fock matrix, i.e. some virtual orbitals will come in energy).
C>
C> To implement the penalty function we use
C> \f{eqnarray*}{
C>    p(x;n) &=& (\mathrm{abs}(x)-x)*(-x)^n
C> \f}
C> which has the properties that \f$ p(x;n) = 0 \f$ if \f$ x \ge 0 \f$,
C> and \f$ p(x;n) = 2\mathrm{abs}(x^{n+1}) \f$ if \f$ x < 0 \f$.
C>
      subroutine wfn1_order(nbf,nea,neb,ehfa,ehfb,ewfa,ewfb,factor,
     +                      En,oerr)
      implicit none
c
      integer nbf !< [Input] The number of basis functions
      integer nea !< [Input] The number of \f$\alpha\f$-electrons
      integer neb !< [Input] The number of \f$\beta\f$-electrons
c
      double precision ehfa(nbf) !< [Input] The Hartree-Fock
      !< \f$\alpha\f$-electron orbital energies
      double precision ehfb(nbf) !< [Input] The Hartree-Fock
      !< \f$\beta\f$-electron orbital energies
      double precision ewfa(nbf) !< [Input] The WFN1
      !< \f$\alpha\f$-electron orbital energies
      double precision ewfb(nbf) !< [Input] The WFN1
      !< \f$\beta\f$-electron orbital energies
      double precision factor    !< [Input] The scale factor for the
      !< penalty function
c
      double precision En !< [Output] The penalty energy term
      double precision oerr !< [Output] The order error
c
c     Local variables
c
      integer nta !< Temporary num of alpha orbitals
      integer ntb !< Temporary num of beta orbitals
      integer ii !< Counter
      integer jj !< Counter
      integer n  !< Order of function p(x;n)
      double precision x !< Argument of statement function
c
c     Statement function
c
      double precision p
      p(x,n) = (abs(x)-x)*(-x)**n
c
c     Code
c
      En = 0.0d0
      oerr = 0.0d0
      do ii = 2, nbf
        oerr = oerr + p(ewfa(ii)-ewfa(ii-1),1)
      enddo
      do ii = 2, nbf
        oerr = oerr + p(ewfb(ii)-ewfb(ii-1),1)
      enddo
      En = oerr * factor
      return



c
c     Establish any degenerate orbitals
c
      nta = max(nea,nbf/2)
      ntb = max(neb,nbf/2)
      jj  = nta
      do ii = jj+1, nbf
        if (abs(ehfa(ii)-ehfa(jj)).lt.1.0d-4) nta = nta + 1
      enddo
      jj  = ntb
      do ii = jj+1, nbf
        if (abs(ehfb(ii)-ehfb(jj)).lt.1.0d-4) ntb = ntb + 1
      enddo
c
      En = 0.0d0
      do ii = 1, nta
        En = En + p(ewfa(ii)-ehfa(ii),0)
      enddo
cDEBUG
c     if (En.gt.0.0d0) write(*,*)'wfn1_order A: ',En
cDEBUG
      do ii = 1, ntb
        En = En + p(ewfb(ii)-ehfb(ii),0)
      enddo
cDEBUG
c     if (En.gt.0.0d0) write(*,*)'wfn1_order B: ',En
cDEBUG
      do ii = 1, nta
        do jj = nta+1, nbf
          En = En + p(ewfa(jj)-ehfa(ii),0)
        enddo
      enddo
cDEBUG
c     if (En.gt.0.0d0) write(*,*)'wfn1_order C: ',En
cDEBUG
      do ii = 1, ntb
        do jj = ntb+1, nbf
          En = En + p(ewfb(jj)-ehfb(ii),0)
        enddo
      enddo
cDEBUG
c     if (En.gt.0.0d0) write(*,*)'wfn1_order D: ',En
cDEBUG
      En = En * factor / 10.0d0
cDEBUG
c     if (En.gt.0.0d0) write(*,*)'wfn1_order E: ',En,factor
cDEBUG
c
      end
C>
C> @}
