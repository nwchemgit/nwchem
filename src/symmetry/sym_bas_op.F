      subroutine sym_bas_op(geom, op, r, maxf, ang_max)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c     
      integer geom
      integer op                ! Desired operator
      integer maxf              ! Leading dimension of r
      integer ang_max           ! Max. ang. momentum of shell
c     
c     Return the transformation matrices for basis functions
c     under the specified group operation.
c     
c     NOTE that the identity operation is NOT included
c     
c     NOTE that only cartesian shells are supported
c     
c     Let X(I,L) be the I'th function in the shell with angular
c     momentum L.  The application of a symmetry operator will
c     map shell X into an equivalent shell on a possibly different
c     center and will also mix up the components of the shell
c     according to
c     
c     R op X(I,L) = sum(J) X(J,L)*R(J,I,L)
c     
c     A function X(J,L) is just a product of the cartesian
c     coordinates raised to some power.  This product may be
c     expressed as
c     
c     X(I,L) = X(IM,L-1) * X(M,1)
c     
c     X(M,1) is of course just a component of a P function and
c     the indices IM and M are determined so that the product
c     yields X(I,L).
c     
c     Thus, the action of R may be expressed in terms of lower
c     angular momenta shells
c     
c     R op X(I,L) = R op X(IM,L-1) * X(M,1)
c     
c     = sum(JK) sum(K) X(JK,L-1) * X(K,1) * 
c     R(JK,IM,L-1) * R(K,M,1)
c     
c     = sum(J) X(J,L)*R(J,I,L) (by definition)
c     
c     We also have from the transformation of coordinates
c     (note that for consistency coordinates are written as
c     a 3 component column vector and basis functions in a
c     shell as a row vector)
c     
c     (R op r)(i) = sum(j) O(i,j) r(j)
c     
c     or
c     
c     R op r = O * r
c     
c     where O is the transformation matrix for coordinates
c     under operation R.
c     
c     And we also have the following result (obtained by
c     considering transformations of the coordinates instead
c     of the functions)
c     
c     R op f(r) = f(R^-1 op r)
c     
c     Considering how the P (L=1) functions transform
c     (note a p-shell is the row vector X(I,1) = (rT)(I))
c     
c     R op X(r;I,1) = X(R^-1 op r; I, 1)
c     = (O(R^-1)*r)T
c     = rT * O(R^-1)T
c     = rT * O(R)
c     = X(r;J,1) * R(J,I,1)
c     
c     Hence, R(J,I,1) = O(J,I) (but note that they are
c     applied from opposite sides).  Also, note that any
c     translational components of R are irrelevant.
c     
c     Finally, R(1,1,1) = 1 (s functions).
c     
      double precision R(1:maxf,1:maxf,0:ang_max)
      double precision matrix(3,4)
      integer II, JJ, KK, LL, MM, IIMM, JJKK, i, j,  nLL, nLLm1
      integer ijk(3), pqr(3), lmn(3)
c     
c     Ensure cartesian basis function polynomials are defined
c     
      call defnxyz(ang_max)
c     
*     call util_print_centered(6,
*     $     'Transformation of basis shells under symmetry operators', 
*     $     30, .true.)
c     
c     Loop through operators in the group
c     
*     write(6,1) op
*     1       format(/' Symmetry operator ', i3)
      call dfill(maxf*maxf*(ang_max+1), 0.0d0, R, 1)
c     
c     Loop thru angular momentua
c     
      do LL = 0, ang_max
         nLL = ((LL+1)*(LL+2))/2
         if (LL .eq. 0) then
c     
c     S functions are easy
c     
            R(1,1,0) = 1.0d0
c     
         else if (LL .eq. 1) then
c     
c     P functions transform using the same matrix that
c     transform coordinates without the translational
c     components ... except of course that the matrix
c     is applied from the other side.
c     
            call sym_get_cart_op(geom, op, matrix)
            do i = 1, 3
               do j = 1, 3
                  R(j,i,1) = matrix(i,j)
*                  R(j,i,1) = matrix(j,i)
               enddo
            enddo
c     
         else
c     
c     Higher functions are built recursively
c     
            nLLm1 = ((LL)*(LL+1))/2
c     
c     Loop thru cartesian functions with sum of powers LL
c     and determine the powers of x, y, z
c     
            do II = 1, nLL
               call getNxyz(LL, II, ijk)
*     *                  write(6,*) ' ijk ', ijk
c     
c     Express this function as a product of functions in (LL-1)*(1) and
c     determine the indices of these functions (IIMM and MM)
c     
               do MM = 1, 3
                  if (ijk(MM) .gt. 0) then
                     ijk(MM) = ijk(MM)-1
*     *                        write(6,*) ' MM ', MM
*     *                        write(6,*) ' reduced ijk ', ijk
                     call matchNxyz(LL-1, IIMM, ijk)
*     *                        write(6,*) ' IIMM ', IIMM
                     goto 10
                  endif
               enddo
               call errquit('sym_bas_op: should never get here ', 0)
 10            continue
c     
c     Loop thru LL-1 functions in the result index
c     
               do JJKK = 1, nLLm1
                  call getNxyz(LL-1, JJKK, pqr)
*     *                     write(6,*) ' pqr ', pqr
c     
c     Loop thru P functions and add onto the LL-1 functions
c     to get a new LL function
c     
                  do KK = 1, 3
                     lmn(1) = pqr(1)
                     lmn(2) = pqr(2)
                     lmn(3) = pqr(3)
                     lmn(KK) = lmn(KK) + 1
*     *                        write(6,*) ' KK ', KK
*     *                        write(6,*) ' lnm ', lmn
                     call matchNxyz(LL, JJ, lmn)
*     *                        write(6,*) ' JJ ', JJ
c     
                     R(JJ,II,LL) = R(JJ,II,LL) + 
     $                    R(JJKK, IIMM, LL-1) * R(KK, MM, 1)
c     
                  enddo
               enddo
            enddo
         endif
c     
*         write(6,2) LL
* 2       format(/'      Angular momentum ', i2)
*         call output(R(1,1,LL), 1, nLL, 1, nLL, maxf, maxf, 1)
      enddo
c     
      end
