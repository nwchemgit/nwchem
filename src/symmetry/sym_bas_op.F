      subroutine sym_bas_op(geom, op, r, maxf, ang_max,
     &    otransform)
C$Id$
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c     
      integer geom         ! [input]  geometry handle
      integer op           ! [input]  Desired operator
      integer maxf         ! [input]  Leading dimension of r
      integer ang_max      ! [input]  Max. ang. momentum of shell
      logical otransform   ! [input]  transform to spherical basis
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
      double precision  ! [output] basis function transformation under operation op
     &    R(1:maxf,1:maxf,0:ang_max) 
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
*      write(6,1) op
*    1    format(/' Symmetry operator ', i3)
      call dfill(maxf*maxf*(ang_max+1), 0.0d0, R, 1)
c     
c     Loop thru angular momenta
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
c     transforms coordinates without the translational
c     components ... except of course that the matrix
c     is applied from the other side.
c     
            call sym_get_cart_op(geom, op, matrix)
            do i = 1, 3
               do j = 1, 3
                  R(j,i,1) = matrix(i,j)
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
c     
c     Express this function as a product of functions in (LL-1)*(1) and
c     determine the indices of these functions (IIMM and MM)
c     
               do MM = 1, 3
                  if (ijk(MM) .gt. 0) then
                     ijk(MM) = ijk(MM)-1
                     call matchNxyz(LL-1, IIMM, ijk)
                     goto 10
                  endif
               enddo
               call errquit('sym_bas_op: should never get here ', 0,
     &       BASIS_ERR)
 10            continue
c     
c     Loop thru LL-1 functions in the result index
c     
               do JJKK = 1, nLLm1
                  call getNxyz(LL-1, JJKK, pqr)
c     
c     Loop thru P functions and add onto the LL-1 functions
c     to get a new LL function
c     
                  do KK = 1, 3
                     lmn(1) = pqr(1)
                     lmn(2) = pqr(2)
                     lmn(3) = pqr(3)
                     lmn(KK) = lmn(KK) + 1
                     call matchNxyz(LL, JJ, lmn)
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
      if (.not.otransform) return
c
      call sym_bas_op_spherical(maxf,ang_max,r)
c
      end
      subroutine sym_bas_op_spherical(maxf,ang_max,r)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "spcartP.fh"
      integer maxf    ! [input] leading dimension(s) of R
      integer ang_max ! [input] trailing dimension of R
      double precision R(1:maxf,1:maxf,0:ang_max)  ! [input/output] cartesian/spherical
c::local
      double precision sum(36,36)  ! intermediate matrix (handles j functions)
      double precision val
*     double precision scale
      integer lval
      integer l2s
      integer s1, s2, c1, c2, r1, r2
      double precision sover(6,6,0:2)
c::functions
      logical spcart_init
      external spcart_init
c::statement function ----- start
      integer iic,iis,iil,lval_max
      double precision Dtrans, Dinvtrans, Dscale
      integer sindex
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
      Dinvtrans(iis,iic,iil) =
     &    dbl_mb((int_mb(k_invsp2c_lindx+iil))+
     &    (iic-1)*(2*iil+1)+iis+iil)
      sindex(iis,lval) = lval+iis+1
c--(2*lm+1)*l + l+i + (lm-l):m=2,a=5
c--(2*lm+1)*l + i + lm      :m=2,a=3
c--lm*(2*l+1) + l + i       :m=2,a=3
      Dscale(iis,lval,lval_max) =
     &    dbl_mb(k_cart_norm_scale + lval_max*(2*lval+1)+lval+iis)
c::statement function ----- end
c
      if (.not.spcart_init(ang_max,.true.,.false.))
     &    call errquit
     &    ('sym_bas_op_spherical: could not initialize spcart',911,
     &       BASIS_ERR)
c
      call dfill ((6*6*3),0.0d00,sover,1)
      sover(1,1,0) = 1.0d00
      sover(1,1,1) = 1.0d00
      sover(2,2,1) = 1.0d00
      sover(3,3,1) = 1.0d00
      sover(1,1,2) = 1.0d00
      sover(2,2,2) = 1.0d00/3.0d00
      sover(3,3,2) = 1.0d00/3.0d00
      sover(4,4,2) = 1.0d00
      sover(5,5,2) = 1.0d00/3.0d00
      sover(6,6,2) = 1.0d00
      sover(1,4,2) = 1.0d00/3.0d00
      sover(1,6,2) = 1.0d00/3.0d00
      sover(4,6,2) = 1.0d00/3.0d00
      sover(4,1,2) = 1.0d00/3.0d00
      sover(6,1,2) = 1.0d00/3.0d00
      sover(6,4,2) = 1.0d00/3.0d00
c
      do lval = 2,ang_max         ! skip s and p transformations
        l2s = (lval+1)*(lval+2)/2
*        write(6,*)' cartesian matrix for angular momentum ',lval
*        call output(R(1,1,lval),1,l2s,1,l2s,maxf,maxf,1)
        call dfill((36*36),0.0d0,sum,1)
        do c2 = 1,l2s
          do s1 = -lval,lval
            val = 0.0d00
            do c1 = 1,l2s
              val = val +
     &            dinvtrans(s1,c1,lval)*R(c1,c2,lval)
            enddo
            sum(sindex(s1,lval),c2) = val
          enddo
        enddo
* now sum is sum(spherical,cartesian,implied lval)
*        write(6,*)' intermediate sum(sph,car)'
*        call output(sum,1,(2*lval+1),1,l2s,36,36,1)
        do s1 = -lval,lval
          do s2 = -lval,lval
            val = 0.0d00
            do c2 = 1,l2s
              val = val +
     &            dtrans(c2,s2,lval)*sum(sindex(s1,lval),c2)
            enddo
            r1 = sindex(s1,lval)
            r2 = sindex(s2,lval)
            R(r1,r2,lval) = val
          enddo
        enddo
*        write(6,*)' spherical matrix for angular momentum ',lval
*        call output(R(1,1,lval),1,(2*lval+1),1,(2*lval+1),maxf,maxf,1)
*        do s1 = 2*lval+2,maxf
*          do s2 = 2*lval+2,maxf
*            R(s1,s2,lval) = 1234.0d00
*          enddo
*        enddo
      enddo
c      
      end
