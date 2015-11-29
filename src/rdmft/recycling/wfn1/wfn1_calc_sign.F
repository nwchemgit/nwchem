C> \ingroup wfn1
C> @{
C>
C> \brief Calculate the sign for each permutation in the table
C>
C> Each permutation has a sign depending on the number of pair
C> permutations needed to get all entries in increasing order.
C> This routine sorts the entries in each permutation to establish
C> the sign.
C>
      subroutine wfn1_calc_sign(ne,nprm,tbl,sgn)
      implicit none
c
      integer ne           !< [Input] The number of particles in each
                           !< permutation
      integer nprm         !< [Input] The number of permutations
      integer tbl(ne,nprm) !< [Input] The table of permutations
c
      double precision sgn(nprm) !< [Output] The sign of each
                                 !< permutation
c
      integer iprm    !< The current permutation
      integer tmp(ne) !< Temporary memory for permutation
      integer ii      !< Counter
      integer jj      !< Counter
      integer small   !< The smallest number
      integer iismall !< The index of the smallest in permutation
c
      do iprm = 1, nprm
        sgn(iprm) = 1.0d0
      enddo
c
      do iprm = 1, nprm
        do ii = 1, ne
          tmp(ii) = tbl(ii,iprm)
        enddo ! ii
        do ii = 1, ne-1
          small = tmp(ii)
          iismall = ii
          do jj = ii + 1, ne
            if (tmp(jj).lt.small) then
              small = tmp(jj)
              iismall = jj
            endif
          enddo ! jj
          if (iismall.ne.ii) then
            tmp(iismall) = tmp(ii)
            tmp(ii) = small
            sgn(iprm) = -sgn(iprm)
          endif
        enddo ! ii
      enddo ! iprm
c
      end
C>
C> @}
