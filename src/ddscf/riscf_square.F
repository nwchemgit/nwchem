      subroutine riscf_square (ao_basis, g_three, g_sthree, three, 
     $     sthree, g_indx, indx, t, in, nbf, nsh, mxbf, nproc)

C$Id: riscf_square.F,v 1.2 1995-12-29 02:04:14 vg038 Exp $

C     get integrals from compressed array g_three and blow up to square matrix
C     g_sthree (both GAs)

      implicit none
      
#include "tcgmsg.fh"
#include "bas.fh"

      integer ao_basis, g_three, g_sthree, g_indx, t, first, indx, in, 
     $     nbf, nsh, mxbf, nproc

      double precision three, sthree

      dimension three(nbf*mxbf), sthree(nbf,mxbf), indx(nsh)

      integer ish, jsh, next, i, j, ind, last, ifirst, ilast, jfirst, 
     $     jlast, ioffs, ni, nint

      double precision zero

      logical status

      parameter (zero=0.0D0)

C     offset in multiple-t array
      ioffs = (in-1) * nbf

      next = nxtval(nproc) + 1
      do ish = 1, nsh

        if ( ish .eq. next ) then
          status = bas_cn2bfr(ao_basis, ish, ifirst, ilast)
          call ga_get (g_indx, 1, ish, ish, ish, indx, ish)

C         first count the integrals
          nint = 0
          first = 0
          ni = ilast - ifirst + 1
          do jsh = 1, ish - 1
            if (indx(jsh).ne.0) then
              if (first.eq.0) then
                first = indx(jsh)
              endif
              status = bas_cn2bfr(ao_basis, jsh, jfirst, jlast)
              nint = nint + (jlast-jfirst+1) * ni
            endif
          enddo
          if (indx(ish).ne.0) then
            if (first.eq.0) then
              first = indx(ish)
            endif
            nint = nint + (ni*(ni+1)/2)
          endif

C         get the integrals
          last = first + nint - 1
          call ga_get (g_three, first, last, t, t, three, nint)

C         place them in sthree
          call dfill (nbf*mxbf, zero, sthree, 1 )
          ind = 1
          do jsh = 1, ish - 1
            if (indx(jsh).ne.0) then
              status = bas_cn2bfr(ao_basis, jsh, jfirst, jlast)
              do i = ifirst, ilast
                do j = jfirst, jlast
                  sthree(j,i-ifirst+1) = three(ind)
                  ind = ind + 1
                enddo
              enddo
            endif
          enddo

          if (indx(ish).ne.0) then
            do i = ifirst, ilast
              do j = ifirst, i - 1
                sthree(i,j-ifirst+1) = three(ind)
                sthree(j,i-ifirst+1) = three(ind)
                ind = ind + 1
              enddo
              sthree(i,i-ifirst+1) = three(ind)
              ind = ind + 1
            enddo
          endif
          next = nxtval(nproc) + 1

C         place patch in global square array
          call ga_put (g_sthree, ioffs+1, ioffs+ilast, ifirst, ilast, 
     $         sthree, nbf)
C         put the transpose of it in the other triangle
          jlast = ilast - ni
          if ( jlast .gt. 0 ) then
C           use three for transpose of sthree
            call transpose (sthree, three, jlast, ni, nbf, ni)
            call ga_put (g_sthree, ioffs+ifirst, ioffs+ilast, 1, jlast, 
     $           three, ni)
          endif 
        endif
      enddo
      next = nxtval(-nproc)

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine transpose (a, b, idim1, idim2, ld1, ld2)
      implicit none
      integer idim1, idim2, ld1, ld2, i, j
      double precision a, b, x
      dimension a(ld1,idim2), b(ld2,idim1)

      do i=1, idim1
        do j=1, idim2
          x = a(i,j)
          a(i,j) = b(j,i)
          b(j,i) = x
        enddo
      enddo

      return
      end


