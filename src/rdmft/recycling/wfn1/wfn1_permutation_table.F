C> \ingroup wfn1
C> @{
C>
C> \brief Compute the table of all permutation of a number of things
C>
C> This routine takes as input the number of things to permute \f$nnum\f$
C> and then starts from the first permutation \f$1, 2, \ldots, nnum\f$
C> and then generates all other permutations.
C>
C> The algorithm uses \f$nnum\f$ counters all starting at \f$0\f$. 
C> The value for a particular position is given by the counter value
C> for that position \f$Cn\f$ and the position number as
C> \f$mod(Cn+n,nnum)+1\f$. Each time we check that the resulting number
C> does not occur for any lower \f$n\f$. When we have picked a value
C> for position \f$n\f$ we increase \f$n\f$. If \f$n\f$ reaches
C> \f$nnum\f$ we have found a new permutation. All permutations have
C> encountered when all counters reach \f$nnum-1\f$.
C> 
      subroutine wfn1_permutation_table(ne,nprm,tbl)
      implicit none
c
      integer ne !< [Input] The number of electrons
      integer nprm !< [Input] The number of permutations
c
      integer tbl(ne,nprm) !< [Output] All permutations of the available
                           !< electrons. This table is a 2D data
                           !< structure with the dimensions:
                           !< - 1. The electron labels for a given
                           !<      permutation
                           !< - 2. The different permutations
c
c     Local
c
      logical oused(ne) !< .True. if already in a given permutation
      integer Cn(ne+1)  !< Counters for each electron
      integer ii        !< Counter
      integer jj        !< Counter
      integer iperm     !< The number of the permutation
c
c     Statement functions
c
      integer l, icnt, ipos, n
      l(icnt,ipos,n) = mod(icnt+ipos,n)+1
c
c     Code
c
      nprm = 0
      if (ne.eq.0) return
      iperm = 1
      do ii = 1, ne
        oused(ii) = .true.
        Cn(ii) = 0
        tbl(ii,iperm) = l(Cn(ii),ii,ne)
      enddo
      ii = ne+1
      do while (iperm.lt.nprm)
        ii = ii - 1
        do while (ii.le.ne)
          oused(l(Cn(ii),ii,ne)) = .false.
          Cn(ii) = Cn(ii) + 1
          do while (oused(l(Cn(ii),ii,ne)).and.Cn(ii).lt.ne)
            Cn(ii) = Cn(ii) + 1
          enddo ! Cn(ii).lt.ne
          if (Cn(ii).ge.ne) then
c
c           Found no opening so go down with the index
c
            ii = ii - 1
c
          else if (.not.oused(l(Cn(ii),ii,ne))) then
c
c           Found opening so claim it and fix the higher indices
c
            oused(l(Cn(ii),ii,ne)) = .true.
            ii = ii + 1
            Cn(ii) = 0
c
          endif
        enddo ! ii.le.ne
        iperm = iperm + 1
        do jj = 1, ne
          tbl(jj,iperm) = l(Cn(jj),jj,ne)
        enddo
      enddo ! iperm.le.nprm
c
      end
C>
C> @}
