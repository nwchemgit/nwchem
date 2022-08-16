C> \ingroup nwdft_xc
C> @{
C>
C> \file xc_eval_basis.F
C> Evaluate basis functions and derivatives
C>
C> \brief Evaluate the basis functions and their derivatives on a grid
C>
C> This routine calculates the values of contracted Gaussians on a grid. 
C> Upon request the subroutine can evaluate up to 3rd order derivatives
C> in addition to the basis function values. In addition the basis 
C> functions may be Cartesian or spherical harmonic basis functions.
C>
C> The derivatives of the basis functions are simple derivatives with
C> respect to the Cartesian coordinates.
C>
      subroutine xc_eval_basis(basis_hand, maxder, chi, d1chi, d2chi,
     &     d3chi, rq, qxyz, ccoord, nq, natoms, iniz, docset,
     &     bas_cent_info, bas_cset_info)
c
c     $Id$
c
c     Evaluate basis functions and their derivatives on a grid
c
c     BGJ - 7/00
c
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
c!!! Will we still need this?
#include "bas.fh"
c
      integer basis_hand        !< [Input] basis set handle
      integer maxder            !< [Input] max derivative order to compute
      integer nq                !< [Input] number of quadrature points
      integer natoms            !< [Input] number of atoms
      double precision chi(nq,*) !< [output] basis function values
      double precision d1chi(nq,3,*) !< [output] basis function 1st derivs
      double precision d2chi(nq,6,*) !< [output] basis function 2nd derivs
      double precision d3chi(nq,10,*) !< [output] basis function 3rd derivs
      double precision rq(nq,natoms) !< [Input] squared pt-ctr distances
      double precision qxyz(3,nq) !< [Input] grid point coordinates
      double precision ccoord(3,natoms) !< [Input] nuclear coordinates
      integer iniz(*)           !< [Input] iniz(i).eq.0 means the basis 
                                !< functions on atom \f$i\f$ are not evaluated,
                                !< only used in derivatives where atoms
                                !< positions may be frozen
      logical docset(*)         !< [Input] whether cset set is significant
      integer bas_cent_info(3,natoms) !< [Input] Basis set info for each
                                      !< center
                                      !< - 1. The maximum angular momentum
                                      !< - 2. The lowest rank shell on the
                                      !<   atom
                                      !< - 3. The highest rank shell on the
                                      !<   atom
      integer bas_cset_info(6,*) !< [Input] Basis set info for each shell
                                 !< - 1. The lowest rank basis function in the
                                 !<   shell
                                 !< - 2. The highest rank basis function in the
                                 !<   shell
                                 !< - 3. The angular momentum of the shell
                                 !< - 4. The number of primitives in the shell
                                 !< - 5. The number of contractions in the shell
                                 !< - 6. Whether the shell is Cartesian or
                                 !<   spherical harmonic
c
      external nang_cart, nang_cart_tot
      integer nang_cart, nang_cart_tot
c
c     Local declarations
c
      integer gen_l
      integer MXprim,MXang,MXcont
      integer izprim,iccoef,iccoef1
      integer max_elem_ang_scr,max_der_component,max_ang_for_der
      integer iang,id1ang,id2ang,id3ang,il2k,idrd,irad,itop
      integer iscr,lscr
      integer ictr,icset,ichi
      integer ifirst,ilast,ltype,nprimo,ncontr,nshbf(2)
      integer l,nang,iang_beg,iang_end
      integer icontr,l0(2)
      integer mx_l2k,k_l2k,l_l2k
      logical spherical,LResult
#define Ax 2
#define Ay 3
#define Az 4
#define Axx 5
#define Axy 6
#define Axz 7
#define Ayy 8
#define Ayz 9
#define Azz 10
#define Axxx 11
#define Axxy 12
#define Axxz 13
#define Axyy 14
#define Axyz 15
#define Axzz 16
#define Ayyy 17
#define Ayyz 18
#define Ayzz 19
#define Azzz 20
c
#define Ad1 5
#define Ad2 6
#define Ad3 7
#define Ad4 8
#define Ad5 9
#define Af1 10
#define Af2 11
#define Af3 12
#define Af4 13
#define Af5 14
#define Af6 15
#define Af7 16
c
#define dzero 1
#define dx 2
#define dy 3
#define dz 4
#define dxx 5
#define dxy 6
#define dxz 7
#define dyy 8
#define dyz 9
#define dzz 10
#define dxxx 11
#define dxxy 12
#define dxxz 13
#define dxyy 14
#define dxyz 15
#define dxzz 16
#define dyyy 17
#define dyyz 18
#define dyzz 19
#define dzzz 20
c
c     Determine when to use general code
c
      spherical = bas_is_spherical(basis_hand)
      if (spherical) then
         if (maxder .eq. 0) then
            gen_l = 4
         elseif (maxder .eq. 1) then
            gen_l = 3
         else
            gen_l = 2
         endif
      else
         if (maxder .le. 1) then
            gen_l = 4
         else
            gen_l = 3
         endif
      endif
c
c     Get the max number of contractions, number of primitives
c     and angular momentum for any given shell of this basis set.
c
      if( .not. bas_ncontr_cn_max(basis_hand,MXcont) )then
         call errquit('Exiting in xc_eval_basis.',1, BASIS_ERR)
      end if
      if( .not. bas_nprim_cn_max(basis_hand,MXprim) )then
         call errquit('Exiting in xc_eval_basis.',2, BASIS_ERR)
      end if
      if( .not. bas_high_angular(basis_hand,MXang) )then
         call errquit('Exiting in xc_eval_basis.',3, BASIS_ERR)
      end if
c
c     Allocate scratch space for primitives for a single shell
c
      izprim = 0
      iccoef = izprim + MXprim
      itop = iccoef + MXprim*MXcont
c
c     Allocate scratch space for angular and radial functions
c     Allocate cartesian no matter what since in the general case
c     the cartesian functions are evaluated first and then transformed
c
      max_elem_ang_scr = nang_cart_tot(MXang)
      max_der_component = nang_cart_tot(maxder)
c
c     If basis function derivatives are needed, always evaluate x, y and z
c     angular functions for derivative purposes no matter what the angular
c     momentum is
c
      max_ang_for_der = min(maxder,1)
      max_elem_ang_scr = max(max_elem_ang_scr,
     &                       nang_cart_tot(max_ang_for_der))
c
      iang = itop
      itop = iang + nq*max_elem_ang_scr
      id1ang = iang
      id2ang = iang
      id3ang = iang
      il2k = iang
      if (MXang .ge. gen_l) then
c!!! Could subtract off everything below gen_l here and change the
c!!! dimensioning of the d*ang arrays to save the space if needed
         if (maxder .ge. 1) then
            id1ang = itop
            il2k = id1ang + nq*3*max_elem_ang_scr
            itop = il2k + (MXang+1)*(MXang+1)*(MXang+1)
         endif
         if (maxder .ge. 2) then
            id2ang = itop
            itop = id2ang + nq*6*max_elem_ang_scr
         endif
         if (maxder .ge. 3) then
            id3ang = itop
            itop = id3ang + nq*10*max_elem_ang_scr
         endif
      endif
      idrd = itop
      irad = idrd + nq*max_der_component*MXcont
      itop = irad + nq*(maxder+1)*MXcont
      if (.not. MA_Push_Get(mt_dbl, itop, 'xc_eval_basis', lscr, iscr))
     &     call errquit('Could not allocate scratch in xc_eval_basis',
     &     itop, MA_ERR)
      mx_l2k=max_ang_for_der
      do ictr=1,natoms
         mx_l2k=max(bas_cent_info(1,ictr),mx_l2k)
      enddo
      if (.not. MA_Push_Get(mt_int,(mx_l2k+1)**3,' l2k',
     Z     l_l2k, k_l2k))
     &     call errquit('Could not allocate l2k in xc_eval_basis',
     &     mx_l2k, MA_ERR)
      izprim = iscr
      iccoef = iccoef + iscr
      iang = iang + iscr
      id1ang = id1ang + iscr
      id2ang = id2ang + iscr
      id3ang = id3ang + iscr
      idrd = idrd + iscr
      irad = irad + iscr
c
c     Offset for basis functions in chi arrays
c
      ichi = 1
c
      do 100 ictr = 1,natoms
         if (iniz(ictr) .eq. 0) goto 100
c
c     Evaluate the angular functions up through the max angular momentum
c     on the current center, so that the angular function values can be
c     for all functions on this center - also make sure that the angular
c     functions needed for derivatives get evaluated
c
         call eval_angular(dbl_mb(iang), dbl_mb(id1ang), dbl_mb(id2ang),
     &        dbl_mb(id3ang), qxyz, ccoord(1,ictr), nq, int_mb(k_l2k),
     &        max(bas_cent_info(1,ictr),max_ang_for_der), maxder, gen_l,
     &        spherical)
c
c     Loop over all contraction sets on this center
c
         do 200 icset = bas_cent_info(2,ictr),bas_cent_info(3,ictr)
c
c     Check if we have neglected this contraction
c
            if (.not. docset(icset)) goto 200
c
c     Get info about current contraction set - first and last basis function,
c     angular momentum, number of primitives, number of contractions and
c     whether spherical harmomic
c
            ifirst = bas_cset_info(1,icset)
            ilast  = bas_cset_info(2,icset)
            ltype  = bas_cset_info(3,icset)
            nprimo = bas_cset_info(4,icset)
            ncontr = bas_cset_info(5,icset)
            spherical = bas_cset_info(6,icset) .ne. 0
c
c     Number of basis functions counting all contractions
c
            nshbf(1) = ilast-ifirst+1
c
c     Angular momentum
c
            if (ltype .lt. 0) then
               l0(1)=0
               l0(2)=1
               nshbf(1)=1
               nshbf(2)=3
            else
               l0(1) = ltype
cold               l = ltype
            endif
c
c     Get exponents and contraction coefficients for this contraction set
c
c!!! These calls could be causing too much overhead - check
            LResult = bas_get_exponent(basis_hand,icset,dbl_mb(izprim))
            LResult = bas_get_coeff(basis_hand,icset,dbl_mb(iccoef))
            do icontr=1,ncontr
               l = l0(icontr)
c
c     Cartesian/spherical harmonic
c
            if (spherical) then
               nang = 2*l+1
               iang_end = (l+1)*(l+1)
               iang_beg = iang_end - nang + 1
            else
               nang = nang_cart(l)
               iang_end = nang_cart_tot(l)
               iang_beg = iang_end - nang + 1
            endif
c
c     Evaluate the radial components of the contracted Gaussian functions
c
            iccoef1=iccoef+(icontr-1)*nprimo
            if(maxder.eq.0) then
               call eval_radial0(dbl_mb(irad), rq(1,ictr), nq, 
     &              dbl_mb(izprim), dbl_mb(iccoef1), nprimo, ncontr)
            else
               call eval_radial(dbl_mb(irad), rq(1,ictr), nq, maxder,
     &              dbl_mb(izprim), dbl_mb(iccoef1), nprimo, ncontr)
            endif
c
c     Compute and store the derivatives of the radial function
c
            call eval_der_radial(dbl_mb(idrd), dbl_mb(irad),
     &           dbl_mb(iang), nq, maxder, max_der_component)

c     Evaluate basis functions and begin evaluation of basis function
c     derivative by multiplying angular function by the radial function
c     and its derivatives
            call eval_ang_x_drad(chi(1,ichi), d1chi(1,1,ichi),
     &           d2chi(1,1,ichi), d3chi(1,1,ichi), dbl_mb(idrd),
     &           dbl_mb(iang), nq, iang_beg, iang_end,
     &           max_der_component)
c
c     Finish off with terms involving angular function derivatives
c
            if (l .gt. 0 .and. maxder .gt. 0) then
               call eval_dang_x_rad(chi(1,ichi), d1chi(1,1,ichi),
     &              d2chi(1,1,ichi), d3chi(1,1,ichi), dbl_mb(idrd),
     &              dbl_mb(iang), dbl_mb(id1ang), dbl_mb(id2ang),
     &              dbl_mb(id3ang), nq, iang_beg, iang_end,
     &              max_der_component, gen_l)
            endif
c
c     Advance offset into chi arrays by number of basis functions
c     just processed
c
            ichi = ichi + nshbf(icontr)
            enddo
c
 200     continue
 100  continue
c
      if (.not. MA_chop_Stack(lscr))
     &     call errquit('xc_eval_basis: could not pop stack',0, MA_ERR)
c
      return
      end
c
C> \brief Evaluate the radial factor of a contracted Gaussian
C>
C> This routine evaluates the radial factor of a contracted Gaussian
C> basis function given the square of the atom-grid point distances.
C> In addition the derivatives of the basis function are calculated
C> also. All derivatives from 0-th order up to and including 3-rd
C> order can be calculated.
C>
C> To calculate the derivatives this routine uses the fact that
C> the derivative of a Gaussian function is again a Gaussian function,
C> i.e.
C> \f{eqnarray*}{
C>   f(r) &=& \sum_i c_i e^{-\alpha_i r^2} \\\\
C>   \frac{\partial f(r)}{\partial r} &=& 
C>       \sum_i -2\alpha_i c_i r e^{-\alpha_i r^2} \\\\
C>   \frac{\partial^2 f(r)}{\partial r^2} 
C>       &=& \sum_i \left\{ -2\alpha_i c_i e^{-\alpha_i r^2}+
C>                  (2\alpha_i r)^2 c_i e^{-\alpha_i r^2}\right\} \\\\
C>   &\ldots&
C> \f}
C> where \f$\alpha_i\f$ are the exponents and \f$c_i\f$ are the 
C> contraction coefficients.
C> This routine however does not evaluate the complete derivatives but
C> only the basic terms, i.e.
C> \f{eqnarray*}{
C>   \left[\frac{\partial^n f(r)}{\partial r^n}\right]_{\mathrm{basic}}
C>   &=& \sum_i (-2\alpha_i)^n c_i e^{-\alpha_i r^2}
C> \f}
C> These quantities are to be combined with the angular quantities 
C> elsewhere to construct the correct total derivatives.
C>
      subroutine eval_radial(crad, r2, nq, maxd, zprim, ccoef,
     *     nprimo, ncontr)
      implicit none
#include "errquit.fh"
      integer nq     !< [Input] The number of grid points
      integer maxd   !< [Input] The maximum order of derivative
      integer nprimo !< [Input] The number of primitives
      integer ncontr !< [Input] The number of contractions
c
      double precision r2(nq) !< [Input] The square of the distance
                              !< between the atom and the grid points
      double precision zprim(nprimo) !< [Input] The exponents
      double precision ccoef(nprimo,ncontr) !< [Input] The contraction
                                            !< coefficients
c
      double precision crad(0:maxd,nq) !< [Output] The radial function
                                       !< and its derivatives, where
                                       !< the derivatives run fastest,
                                       !< and the grid points slowest
c
      integer nd
      integer n, iprimo
      double precision zeta, t(0:3), tmp
#ifdef CHKUNDFLW
      double precision rmax,undovl,zz
c 20 digits ln(10)=2.3025
      parameter(undovl=-20d0*2.3025d0)
      rmax=r2(1)
      do n=2,nq
         rmax=max(rmax,r2(n))
      enddo
#endif
c
c
c     Store results for first primitive directly and then sum the rest
c
      if (maxd .gt. 3) then
        call errquit('Unexpected in eval_radial', maxd, INPUT_ERR)
      endif
      zeta = -zprim(1)
      t(0) = ccoef(1,1)
      do n = 1, maxd
         t(n) = zeta*t(n-1)*2.0d0
      enddo
#ifdef CHKUNDFLW
      if(rmax*zeta.lt.undovl) then
         zz=undovl/zeta
         do n = 1,nq
            if(r2(n).lt.zz) then
               tmp = exp(zeta*r2(n))
            else
               tmp =0d0
            endif
            do nd = 0,maxd
               crad(nd,n) = t(nd)*tmp
            enddo
         enddo
      else
#endif
      do n = 1,nq
         tmp = exp(zeta*r2(n))
         do nd = 0,maxd
            crad(nd,n) = t(nd)*tmp
         enddo
      enddo
#ifdef CHKUNDFLW
      endif
#endif
c
      do iprimo = 2, nprimo
         zeta = -zprim(iprimo)
         t(0) = ccoef(iprimo,1)
         do n = 1, maxd
            t(n) = 2d0*zeta*t(n-1)
         enddo
         do n = 1,nq
            tmp = exp(zeta*r2(n))
            do nd = 0,maxd
               crad(nd,n) = crad(nd,n) + t(nd)*tmp
            enddo
         enddo
      enddo
c
      return
      end
c
C> \brief Evaluate the angular factor of a contracted Gaussian
C>
C> This routine evaluates the angular factor of a contracted Gaussian
C> basis function. If desired it will also calculate up-to 3-rd order
C> derivatives of the angular function.
C>
C> To evaluate the angular components the subroutine is given the
C> coordinates of the atom on which the basis functions are centered
C> as well as the coordinates of the grid points. The angular momentum
C> functions can be evaluated as spherical harmonics or Cartesian 
C> functions.
C>
C> Examples:
C> - S-functions: Angular factor is 1, all derivatives are 0.
C> - P-functions: Angular factors are \f$x\f$, \f$y\f$, and \f$z\f$,
C>   the derivatives for \f$x\f$ are 1, 0, and 0, for \f$y\f$ they are
C>   0, 1, and 0, and for \f$z\f$ they are 0, 0, and 1, the 2nd 
C>   derivatives and higher are all 0.
C>
      subroutine eval_angular(ang, d1ang, d2ang, d3ang, q, ctr, nq,
     &     l2k, maxl, maxder, gen_l, spherical)
      implicit none
#include "errquit.fh"
c
#include "mafdecls.fh"
c
      integer nq        !< [Input] The number of grid points
      integer maxl      !< [Input] The maximum angular momentum in basis
                        !< functions
      integer maxder    !< [Input] The maximum order of derivative
      integer gen_l     !< [Input] The momentum above which the angular
                        !< momentum part should be evaluated with the
                        !< general expression
      logical spherical !< [Input] If .true. the angular momentum 
                        !< functions will be evaluated in spherical
                        !< harmonics, otherwise Cartesian functions 
                        !< are evaluated
c
      double precision q(3,nq) !< [Input] The grid point coordinates
      double precision ctr(3)  !< [Input] The atom coordinates
c
      integer l2k(0:maxl,0:maxl,0:maxl) !< [Scratch] Maps the power of
                                        !< x, y, and z to the rank of
                                        !< a particular basis function
      double precision ang(nq,*)             !< [Output] The angular
                                             !< momentum functions
      double precision d1ang(nq,dx:dz,*)     !< [Output] The 1st
                                             !< derivative of the
                                             !< angular momentum
                                             !< functions
      double precision d2ang(nq,dxx:dzz,*)   !< [Output] The 2nd
                                             !< derivative of the
                                             !< angular momentum
                                             !< functions
      double precision d3ang(nq,dxxx:dzzz,*) !< [Output] The 3rd
                                             !< derivative of the
                                             !< angular momentum
                                             !< functions
c
      double precision cd(2), cf(7)
c
      integer nl2
      integer n, l, nl, nl1, i1s, i1, ix, iy, i, nx, nxy, lx, ly, lz,
     &     k, k1(dx:dz), k2(dxx:dzz), k3(dxxx:dzzz), iscr, lscr, nscr
      double precision x, y, z,
     &     fac1(dx:dz), fac2(dxx:dzz), fac3(dxxx:dzzz)
c
      external nang_cart, nang_cart_tot
      integer nang_cart, nang_cart_tot,nxyz(3)
c
      if (maxl .eq. 0) then
         return
      elseif (maxl .eq. 1 .or. gen_l .eq. 2) then
         do n = 1, nq
            ang(n,Ax) = q(1,n) - ctr(1)
            ang(n,Ay) = q(2,n) - ctr(2)
            ang(n,Az) = q(3,n) - ctr(3)
         enddo
      elseif (maxl .eq. 2 .or. gen_l .eq. 3) then
         if (spherical .and.
     &        (maxl .lt. gen_l .or. gen_l - maxder .gt. 2)) then
            cd(1) = sqrt(3.d0)
            cd(2) = cd(1)*0.5d0
            do n = 1, nq
               x = q(1,n) - ctr(1)
               y = q(2,n) - ctr(2)
               z = q(3,n) - ctr(3)
c
               ang(n,Ax) = x
               ang(n,Ay) = y
               ang(n,Az) = z
c
               ang(n,Ad1) = cd(1)*x*y
               ang(n,Ad2) = cd(1)*y*z
               ang(n,Ad3) = 0.5d0*(2d0*z*z-x*x-y*y)
               ang(n,Ad4) =-cd(1)*x*z
               ang(n,Ad5) = cd(2)*(x*x-y*y)
            enddo
         else
            do n = 1, nq
               x = q(1,n) - ctr(1)
               y = q(2,n) - ctr(2)
               z = q(3,n) - ctr(3)
c
               ang(n,Ax) = x
               ang(n,Ay) = y
               ang(n,Az) = z
c
               ang(n,Axx) = x*x
               ang(n,Axy) = x*y
               ang(n,Axz) = x*z
               ang(n,Ayy) = y*y
               ang(n,Ayz) = y*z
               ang(n,Azz) = z*z
            enddo
         endif
      elseif (maxl .eq. 3 .or. gen_l .eq. 4) then
         if (spherical .and.
     &        (maxl .lt. gen_l .or. gen_l - maxder .gt. 3)) then
            cd(1) = sqrt(3.d0)
            cd(2) = cd(1)*0.5d0
            cf(1) = sqrt(2.5d0)*0.5d0
            cf(2) = cf(1)*3d0
            cf(3) = sqrt(15.d0)
            cf(4) = sqrt(1.5d0)*0.5d0
            cf(5) = sqrt(6.d0)
            cf(6) = 1.5d0
            cf(7) = cf(3)*0.5d0
            do n = 1, nq
               x = q(1,n) - ctr(1)
               y = q(2,n) - ctr(2)
               z = q(3,n) - ctr(3)
c
               ang(n,Ax) = x
               ang(n,Ay) = y
               ang(n,Az) = z
c
               ang(n,Ad1) = cd(1)*x*y
               ang(n,Ad2) = cd(1)*y*z
               ang(n,Ad3) = 0.5d0*(2d0*z*z-x*x-y*y)
               ang(n,Ad4) =-cd(1)*x*z
               ang(n,Ad5) = cd(2)*(x*x-y*y)
c
               ang(n,Af1) = y*(cf(2)*x*x-cf(1)*y*y)
               ang(n,Af2) = cf(3)*x*y*z
               ang(n,Af3) = y*(cf(5)*z*z-cf(4)*(x*x+y*y))
               ang(n,Af4) = z*(z*z-cf(6)*(x*x+y*y))
               ang(n,Af5) =-x*(cf(5)*z*z-cf(4)*(x*x+y*y))
               ang(n,Af6) = z*cf(7)*(x*x-y*y)
               ang(n,Af7) = x*(cf(2)*y*y-cf(1)*x*x)
            enddo
         else
            do n = 1, nq
               x = q(1,n) - ctr(1)
               y = q(2,n) - ctr(2)
               z = q(3,n) - ctr(3)
c
               ang(n,Ax) = x
               ang(n,Ay) = y
               ang(n,Az) = z
c
               ang(n,Axx) = x*x
               ang(n,Axy) = x*y
               ang(n,Axz) = x*z
               ang(n,Ayy) = y*y
               ang(n,Ayz) = y*z
               ang(n,Azz) = z*z
c
               ang(n,Axxx) = x*ang(n,Axx)
               ang(n,Axxy) = x*ang(n,Axy)
               ang(n,Axxz) = x*ang(n,Axz)
               ang(n,Axyy) = x*ang(n,Ayy)
               ang(n,Axyz) = x*ang(n,Ayz)
               ang(n,Axzz) = x*ang(n,Azz)
               ang(n,Ayyy) = y*ang(n,Ayy)
               ang(n,Ayyz) = y*ang(n,Ayz)
               ang(n,Ayzz) = y*ang(n,Azz)
               ang(n,Azzz) = z*ang(n,Azz)
            enddo
         endif
      else
         call errquit('Unexpected in eval_angular', maxl, INPUT_ERR)
      endif
c
c     Finish off maxl >= gen_l with general code
c
      if (maxl .ge. gen_l) then

c     Fill the l2k indexing array and initialize derivative values if needed

         if (maxder .gt. 0) then
            k = 0
            do l = 0, maxl
               do lx = l, 0, -1
                  do ly = l - lx, 0, -1
                     lz = l - lx - ly
                     k = k + 1
                     l2k(lx,ly,lz) = k
                  enddo
               enddo
            enddo
c     Just to be on the safe side if we need it in the general code,
c     set the 0th angular function to 1
            call dfill(nq, 1d0, ang(1,1), 1)
            nl1 = nang_cart_tot(gen_l-1)
            nl = nang_cart_tot(maxl) - nl1
            call dfill(nq*3*nl, 0d0, d1ang(1,dx,nl1+1), 1)
            if (maxder .ge. 2) then
               call dfill(nq*6*nl, 0d0, d2ang(1,dxx,nl1+1), 1)
            endif
            if (maxder .ge. 3) then
               call dfill(nq*10*nl, 0d0, d3ang(1,dxxx,nl1+1), 1)
            endif
         endif

         do l = gen_l, maxl

            nl  = nang_cart(l)
            nl1 = nang_cart(l-1)
            nl2 = nang_cart_tot(l-2)
            i1s = nang_cart_tot(l-1)
            ix  = nl1
            iy  = nl - 1
            nx  = nang_cart(l-2)
            nxy = nl1 - nx
            if(.not.spherical .or. (l.gt.gen_l)) then
c     Evaluate the remaining cartesian functions

              do i=1,nl1
                ang(:,i1s+i) = ang(:,Ax)*ang(:,nl2+i)
              enddo

              do i=1,nxy
                ang(:,i1s+nl1+i) = ang(:,Ay)*ang(:,nl2+nx+i)
              enddo

              ang(:,i1s+nl1+nxy+1) = ang(:,Az)*ang(:,i1s)
         else
            nl1 = nang_cart_tot(l-1)
            do i=1,((l+1)*(l+2))/2
               call getNxyz(l,i,Nxyz)

               if (nxyz(1).gt.1) then
                 ang(:,nl1+i) = ang(:,Ax)**nxyz(1)
               elseif(nxyz(1).eq.1) then
                 ang(:,nl1+i) = ang(:,Ax)
               else
                 ang(:,nl1+i) = 1d0
               endif

               if (nxyz(2).gt.1) then
                 ang(:,nl1+i) = ang(:,nl1+i)*ang(:,Ay)**nxyz(2)
               elseif(nxyz(2).eq.1) then
                 ang(:,nl1+i) = ang(:,nl1+i)*ang(:,Ay)
               endif

               if (nxyz(3).gt.1) then
                 ang(:,nl1+i) = ang(:,nl1+i)*ang(:,Az)**nxyz(3)
               elseif(nxyz(3).eq.1) then
                 ang(:,nl1+i) = ang(:,nl1+i)*ang(:,Az)
               endif

            enddo
         endif
#ifdef DEBUG
         write(0,*) 'l l1l*l ',l,nl1,l*l
         do i=1,((l+1)*(l+2))/2
            do n=1,nq
               write(0,*) ' cart l q a ',i,n,ang(n,nl1+i)
            enddo
         enddo
#endif

c     Construct cartesian angular derivatives if needed
c     Only need to evaluate the nonzero terms - everything has been
c     zeroed out initially

            if (maxder .ge. 1) then
               do lx = l-1, 0, -1
                  do ly = (l-1) - lx, 0, -1
                     lz = (l-1) - lx - ly
                     k = l2k(lx,ly,lz)
                     k1(dx) = l2k(lx+1,ly,lz)
                     k1(dy) = l2k(lx,ly+1,lz)
                     k1(dz) = l2k(lx,ly,lz+1)
                     fac1(dx) = lx + 1
                     fac1(dy) = ly + 1
                     fac1(dz) = lz + 1
                     do n = 1, nq
                        d1ang(n,dx,k1(dx)) = fac1(dx) * ang(n,k)
                        d1ang(n,dy,k1(dy)) = fac1(dy) * ang(n,k)
                        d1ang(n,dz,k1(dz)) = fac1(dz) * ang(n,k)
                     enddo
                  enddo
               enddo
            endif
            if (maxder .ge. 2) then
               do lx = l-2, 0, -1
                  do ly = (l-2) - lx, 0, -1
                     lz = (l-2) - lx - ly
                     k = l2k(lx,ly,lz)
                     k2(dxx) = l2k(lx+2,ly,lz)
                     k2(dxy) = l2k(lx+1,ly+1,lz)
                     k2(dxz) = l2k(lx+1,ly,lz+1)
                     k2(dyy) = l2k(lx,ly+2,lz)
                     k2(dyz) = l2k(lx,ly+1,lz+1)
                     k2(dzz) = l2k(lx,ly,lz+2)
                     fac2(dxx) = (lx + 1)*(lx + 2)
                     fac2(dxy) = (lx + 1)*(ly + 1)
                     fac2(dxz) = (lx + 1)*(lz + 1)
                     fac2(dyy) = (ly + 1)*(ly + 2)
                     fac2(dyz) = (ly + 1)*(lz + 1)
                     fac2(dzz) = (lz + 1)*(lz + 2)
                     do n = 1, nq
                        d2ang(n,dxx,k2(dxx)) = fac2(dxx) * ang(n,k)
                        d2ang(n,dxy,k2(dxy)) = fac2(dxy) * ang(n,k)
                        d2ang(n,dxz,k2(dxz)) = fac2(dxz) * ang(n,k)
                        d2ang(n,dyy,k2(dyy)) = fac2(dyy) * ang(n,k)
                        d2ang(n,dyz,k2(dyz)) = fac2(dyz) * ang(n,k)
                        d2ang(n,dzz,k2(dzz)) = fac2(dzz) * ang(n,k)
                     enddo
                  enddo
               enddo
            endif
            if (maxder .ge. 3) then
               do lx = l-3, 0, -1
                  do ly = (l-3) - lx, 0, -1
                     lz = (l-3) - lx - ly
                     k = l2k(lx,ly,lz)
                     k3(dxxx) = l2k(lx+3,ly,lz)
                     k3(dxxy) = l2k(lx+2,ly+1,lz)
                     k3(dxxz) = l2k(lx+2,ly,lz+1)
                     k3(dxyy) = l2k(lx+1,ly+2,lz)
                     k3(dxyz) = l2k(lx+1,ly+1,lz+1)
                     k3(dxzz) = l2k(lx+1,ly,lz+2)
                     k3(dyyy) = l2k(lx,ly+3,lz)
                     k3(dyyz) = l2k(lx,ly+2,lz+1)
                     k3(dyzz) = l2k(lx,ly+1,lz+2)
                     k3(dzzz) = l2k(lx,ly,lz+3)
                     fac3(dxxx) = (lx + 1)*(lx + 2)*(lx + 3)
                     fac3(dxxy) = (lx + 1)*(lx + 2)*(ly + 1)
                     fac3(dxxz) = (lx + 1)*(lx + 2)*(lz + 1)
                     fac3(dxyy) = (lx + 1)*(ly + 1)*(ly + 2)
                     fac3(dxyz) = (lx + 1)*(ly + 1)*(lz + 1)
                     fac3(dxzz) = (lx + 1)*(lz + 1)*(lz + 2)
                     fac3(dyyy) = (ly + 1)*(ly + 2)*(ly + 3)
                     fac3(dyyz) = (ly + 1)*(ly + 2)*(lz + 1)
                     fac3(dyzz) = (ly + 1)*(lz + 1)*(lz + 2)
                     fac3(dzzz) = (lz + 1)*(lz + 2)*(lz + 3)
                     do n = 1, nq
                        d3ang(n,dxxx,k3(dxxx)) = fac3(dxxx) * ang(n,k)
                        d3ang(n,dxxy,k3(dxxy)) = fac3(dxxy) * ang(n,k)
                        d3ang(n,dxxz,k3(dxxz)) = fac3(dxxz) * ang(n,k)
                        d3ang(n,dxyy,k3(dxyy)) = fac3(dxyy) * ang(n,k)
                        d3ang(n,dxyz,k3(dxyz)) = fac3(dxyz) * ang(n,k)
                        d3ang(n,dxzz,k3(dxzz)) = fac3(dxzz) * ang(n,k)
                        d3ang(n,dyyy,k3(dyyy)) = fac3(dyyy) * ang(n,k)
                        d3ang(n,dyyz,k3(dyyz)) = fac3(dyyz) * ang(n,k)
                        d3ang(n,dyzz,k3(dyzz)) = fac3(dyzz) * ang(n,k)
                        d3ang(n,dzzz,k3(dzzz)) = fac3(dzzz) * ang(n,k)
                     enddo
                  enddo
               enddo
            endif
         enddo

c     Transform to spherical if needed

         if (spherical) then

c     Allocate scratch space for transformation

            nscr = nq*(2*maxl+1)
            if (maxder .eq. 1) then
               nscr = nscr * 3
            elseif (maxder .eq. 2) then
               nscr = nscr * 6
            elseif (maxder .eq. 3) then
               nscr = nscr * 10
            endif
            if (.not. MA_Push_Get(mt_dbl, nscr, 'eval_angular',
     &           lscr, iscr)) call errquit(
     &           'Could not allocate scratch in eval_angular',nscr,
     &       MA_ERR)

            do 100 l = max(2,gen_l-maxder), maxl
               nl1 = nang_cart_tot(l-1) + 1
               nl = l*l + 1
               call spcart_a_s(ang(1,nl1),dbl_mb(iscr),nq,l,
     &              1,.false.,.false.)
               call dcopy(nq*(2*l+1),dbl_mb(iscr),1,ang(1,nl),1)
#ifdef DEBUG
               do i=0,2*l
                  do n=1,nq
                     write(0,*) ' sphe l q a ',i,n,ang(i,nl+n)
                  enddo
               enddo
#endif
               if (l .lt. gen_l) goto 100

               if (maxder .ge. 1) then
                  call spcart_a_s(d1ang(1,dx,nl1),dbl_mb(iscr),nq*3,l,
     &                 1,.false.,.false.)
                  call dcopy(nq*3*(2*l+1),dbl_mb(iscr),1,
     &                 d1ang(1,dx,nl),1)
               endif
               if (maxder .ge. 2) then
                  call spcart_a_s(d2ang(1,dxx,nl1),dbl_mb(iscr),nq*6,l,
     &                 1,.false.,.false.)
                  call dcopy(nq*6*(2*l+1),dbl_mb(iscr),1,
     &                 d2ang(1,dxx,nl),1)
               endif
               if (maxder .ge. 3) then
                  call spcart_a_s(d3ang(1,dxxx,nl1),dbl_mb(iscr),nq*10,
     &                 l,1,.false.,.false.)
                  call dcopy(nq*10*(2*l+1),dbl_mb(iscr),1,
     &                 d3ang(1,dxxx,nl),1)
               endif
 100        continue

            if (.not. MA_Pop_Stack(lscr))
     &           call errquit('eval_angular: could not pop stack',0,
     &       MA_ERR)
         endif
      endif
c
      return
      end
c
c!!! There has to be an api call for this !!!
      function nang_cart_tot(l)
      implicit none
      integer nang_cart_tot, l
c
      nang_cart_tot = (l+1)*(l+2)*(l+3)/6
c
      return
      end
c
c!!! There has to be an api call for this !!!
      function nang_cart(l)
      implicit none
      integer nang_cart, l
      external nang_cart_tot
      integer nang_cart_tot
c
      nang_cart = nang_cart_tot(l) - nang_cart_tot(l-1)
c
      return
      end
c
      subroutine eval_der_radial(drad, rad, ang, nq, maxd, maxcomp)
      implicit none
#include "errquit.fh"
      integer nq      !< [Input] The number of grid points
      integer maxd    !< [Input] The maximum order of derivative
      integer maxcomp !< [Input] The leading dimension of `drad`
c
      double precision rad(0:maxd,nq) !< [Input] The radial derivatives
      double precision ang(nq,*)      !< [Input] The angular momentum
                                      !< functions
c
      double precision drad(maxcomp,nq) !< [Output] The radial
                                        !< derivatives times the angular
                                        !< momentum functions
c
      integer n
      double precision r1, r2, r3, x, y, z
c
      if (maxd .eq. 0) then
         call dcopy(nq, rad, 1, drad, 1)
      elseif (maxd .eq. 1) then
         do n = 1, nq
            drad(dzero,n) = rad(0,n)
            r1 = rad(1,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(dx,n) = x*r1
            drad(dy,n) = y*r1
            drad(dz,n) = z*r1
         enddo
      elseif (maxd .eq. 2) then
         do n = 1, nq
            drad(dzero,n) = rad(0,n)
            r1 = rad(1,n)
            r2 = rad(2,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(dx,n) = x*r1
            drad(dy,n) = y*r1
            drad(dz,n) = z*r1
            drad(dxx,n) = x*x*r2 + r1
            drad(dxy,n) = x*y*r2
            drad(dxz,n) = x*z*r2
            drad(dyy,n) = y*y*r2 + r1
            drad(dyz,n) = y*z*r2
            drad(dzz,n)= z*z*r2 + r1
         enddo
      elseif (maxd .eq. 3) then
         do n = 1, nq
            drad(dzero,n) = rad(0,n)
            r1 = rad(1,n)
            r2 = rad(2,n)
            r3 = rad(3,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(dx,n) = x*r1
            drad(dy,n) = y*r1
            drad(dz,n) = z*r1
            drad(dxx,n) = x*x*r2 + r1
            drad(dxy,n) = x*y*r2
            drad(dxz,n) = x*z*r2
            drad(dyy,n) = y*y*r2 + r1
            drad(dyz,n) = y*z*r2
            drad(dzz,n)= z*z*r2 + r1
            drad(dxxx,n)= x*(3d0*r2 + x*x*r3)
            drad(dxxy,n)= y*(r2 + x*x*r3)
            drad(dxxz,n)= z*(r2 + x*x*r3)
            drad(dxyy,n)= x*(r2 + y*y*r3)
            drad(dxyz,n)= x*y*z*r3
            drad(dxzz,n)= x*(r2 + z*z*r3)
            drad(dyyy,n)= y*(3d0*r2 + y*y*r3)
            drad(dyyz,n)= z*(r2 + y*y*r3)
            drad(dyzz,n)= y*(r2 + z*z*r3)
            drad(dzzz,n)= z*(3d0*r2 + z*z*r3)
         enddo
      else
         call errquit('Bad value for maxd in eval_der_radial',maxd,
     &       INPUT_ERR)
      endif
c
      return
      end
c
      subroutine eval_ang_x_drad(chi, d1chi, d2chi, d3chi, drad, ang,
     &     nq, iabeg, iaend, max_der_component)
c
      implicit none
#include "mafdecls.fh"
      integer nq, iabeg, iaend, max_der_component
      double precision chi(nq,iabeg:iaend),
     &     d1chi(nq,dx:dz,iabeg:iaend), d2chi(nq,dxx:dzz,iabeg:iaend),
     &     d3chi(nq,dxxx:dzzz,iabeg:iaend), drad(max_der_component,nq),
     &     ang(nq,*)
c
      integer n, iang, irem
c
      if (iabeg .eq. 1) then

c     Special case: s functions - just copy

         if (max_der_component .eq. dzero) then
           call dcopy(nq, drad, 1, chi, 1)
         elseif (max_der_component .eq. dz) then
           call dcopy(nq, drad, max_der_component, chi, 1)
           call dcopy(
     .           nq,drad(dx,1),max_der_component,d1chi(1,dx,1),1)
           call dcopy(
     .           nq,drad(dy,1),max_der_component,d1chi(1,dy,1),1)
           call dcopy(
     .           nq,drad(dz,1),max_der_component,d1chi(1,dz,1),1)
!           do n = 1, nq
!               chi(n,1) = drad(dzero,n)
!              d1chi(n,dx,1) = drad(dx,n)
!              d1chi(n,dy,1) = drad(dy,n)
!              d1chi(n,dz,1) = drad(dz,n)
!           enddo
         elseif (max_der_component .eq. dzz) then
            do n = 1, nq
               chi(n,1) = drad(dzero,n)
               d1chi(n,dx,1) = drad(dx,n)
               d1chi(n,dy,1) = drad(dy,n)
               d1chi(n,dz,1) = drad(dz,n)
               d2chi(n,dxx,1) = drad(dxx,n)
               d2chi(n,dxy,1) = drad(dxy,n)
               d2chi(n,dxz,1) = drad(dxz,n)
               d2chi(n,dyy,1) = drad(dyy,n)
               d2chi(n,dyz,1) = drad(dyz,n)
               d2chi(n,dzz,1) = drad(dzz,n)
            enddo
         elseif (max_der_component .eq. dzzz) then
            do n = 1, nq
               chi(n,1) = drad(dzero,n)
               d1chi(n,dx,1) = drad(dx,n)
               d1chi(n,dy,1) = drad(dy,n)
               d1chi(n,dz,1) = drad(dz,n)
               d2chi(n,dxx,1) = drad(dxx,n)
               d2chi(n,dxy,1) = drad(dxy,n)
               d2chi(n,dxz,1) = drad(dxz,n)
               d2chi(n,dyy,1) = drad(dyy,n)
               d2chi(n,dyz,1) = drad(dyz,n)
               d2chi(n,dzz,1) = drad(dzz,n)
               d3chi(n,dxxx,1) = drad(dxxx,n)
               d3chi(n,dxxy,1) = drad(dxxy,n)
               d3chi(n,dxxz,1) = drad(dxxz,n)
               d3chi(n,dxyy,1) = drad(dxyy,n)
               d3chi(n,dxyz,1) = drad(dxyz,n)
               d3chi(n,dxzz,1) = drad(dxzz,n)
               d3chi(n,dyyy,1) = drad(dyyy,n)
               d3chi(n,dyyz,1) = drad(dyyz,n)
               d3chi(n,dyzz,1) = drad(dyzz,n)
               d3chi(n,dzzz,1) = drad(dzzz,n)
            enddo
         endif
      else

c     General case: multiply derivatives of radial functions by angular
c     functions

c     Do remainder first

         irem = (iaend-iabeg+1)-((iaend-iabeg+1)/3)*3
         if (irem .eq. 1) then
            if (max_der_component .eq. dzero) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
               enddo
            elseif (max_der_component .eq. dz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
c
                  d2chi(n,dxx,iabeg  ) = ang(n,iabeg  )*drad(dxx,n)
                  d2chi(n,dxy,iabeg  ) = ang(n,iabeg  )*drad(dxy,n)
                  d2chi(n,dxz,iabeg  ) = ang(n,iabeg  )*drad(dxz,n)
                  d2chi(n,dyy,iabeg  ) = ang(n,iabeg  )*drad(dyy,n)
                  d2chi(n,dyz,iabeg  ) = ang(n,iabeg  )*drad(dyz,n)
                  d2chi(n,dzz,iabeg  ) = ang(n,iabeg  )*drad(dzz,n)
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
c
                  d2chi(n,dxx,iabeg  ) = ang(n,iabeg  )*drad(dxx,n)
                  d2chi(n,dxy,iabeg  ) = ang(n,iabeg  )*drad(dxy,n)
                  d2chi(n,dxz,iabeg  ) = ang(n,iabeg  )*drad(dxz,n)
                  d2chi(n,dyy,iabeg  ) = ang(n,iabeg  )*drad(dyy,n)
                  d2chi(n,dyz,iabeg  ) = ang(n,iabeg  )*drad(dyz,n)
                  d2chi(n,dzz,iabeg  ) = ang(n,iabeg  )*drad(dzz,n)
c
                  d3chi(n,dxxx,iabeg  ) = ang(n,iabeg  )*drad(dxxx,n)
                  d3chi(n,dxxy,iabeg  ) = ang(n,iabeg  )*drad(dxxy,n)
                  d3chi(n,dxxz,iabeg  ) = ang(n,iabeg  )*drad(dxxz,n)
                  d3chi(n,dxyy,iabeg  ) = ang(n,iabeg  )*drad(dxyy,n)
                  d3chi(n,dxyz,iabeg  ) = ang(n,iabeg  )*drad(dxyz,n)
                  d3chi(n,dxzz,iabeg  ) = ang(n,iabeg  )*drad(dxzz,n)
                  d3chi(n,dyyy,iabeg  ) = ang(n,iabeg  )*drad(dyyy,n)
                  d3chi(n,dyyz,iabeg  ) = ang(n,iabeg  )*drad(dyyz,n)
                  d3chi(n,dyzz,iabeg  ) = ang(n,iabeg  )*drad(dyzz,n)
                  d3chi(n,dzzz,iabeg  ) = ang(n,iabeg  )*drad(dzzz,n)
               enddo
            endif
         elseif (irem .eq. 2) then
            if (max_der_component .eq. dzero) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(dzero,n)
               enddo
            elseif (max_der_component .eq. dz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
                  d1chi(n,dx,iabeg+1) = ang(n,iabeg+1)*drad(dx,n)
                  d1chi(n,dy,iabeg+1) = ang(n,iabeg+1)*drad(dy,n)
                  d1chi(n,dz,iabeg+1) = ang(n,iabeg+1)*drad(dz,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
                  d1chi(n,dx,iabeg+1) = ang(n,iabeg+1)*drad(dx,n)
                  d1chi(n,dy,iabeg+1) = ang(n,iabeg+1)*drad(dy,n)
                  d1chi(n,dz,iabeg+1) = ang(n,iabeg+1)*drad(dz,n)
c
                  d2chi(n,dxx,iabeg  ) = ang(n,iabeg  )*drad(dxx,n)
                  d2chi(n,dxy,iabeg  ) = ang(n,iabeg  )*drad(dxy,n)
                  d2chi(n,dxz,iabeg  ) = ang(n,iabeg  )*drad(dxz,n)
                  d2chi(n,dyy,iabeg  ) = ang(n,iabeg  )*drad(dyy,n)
                  d2chi(n,dyz,iabeg  ) = ang(n,iabeg  )*drad(dyz,n)
                  d2chi(n,dzz,iabeg  ) = ang(n,iabeg  )*drad(dzz,n)
                  d2chi(n,dxx,iabeg+1) = ang(n,iabeg+1)*drad(dxx,n)
                  d2chi(n,dxy,iabeg+1) = ang(n,iabeg+1)*drad(dxy,n)
                  d2chi(n,dxz,iabeg+1) = ang(n,iabeg+1)*drad(dxz,n)
                  d2chi(n,dyy,iabeg+1) = ang(n,iabeg+1)*drad(dyy,n)
                  d2chi(n,dyz,iabeg+1) = ang(n,iabeg+1)*drad(dyz,n)
                  d2chi(n,dzz,iabeg+1) = ang(n,iabeg+1)*drad(dzz,n)
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(dzero,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(dzero,n)
c
                  d1chi(n,dx,iabeg  ) = ang(n,iabeg  )*drad(dx,n)
                  d1chi(n,dy,iabeg  ) = ang(n,iabeg  )*drad(dy,n)
                  d1chi(n,dz,iabeg  ) = ang(n,iabeg  )*drad(dz,n)
                  d1chi(n,dx,iabeg+1) = ang(n,iabeg+1)*drad(dx,n)
                  d1chi(n,dy,iabeg+1) = ang(n,iabeg+1)*drad(dy,n)
                  d1chi(n,dz,iabeg+1) = ang(n,iabeg+1)*drad(dz,n)
c
                  d2chi(n,dxx,iabeg  ) = ang(n,iabeg  )*drad(dxx,n)
                  d2chi(n,dxy,iabeg  ) = ang(n,iabeg  )*drad(dxy,n)
                  d2chi(n,dxz,iabeg  ) = ang(n,iabeg  )*drad(dxz,n)
                  d2chi(n,dyy,iabeg  ) = ang(n,iabeg  )*drad(dyy,n)
                  d2chi(n,dyz,iabeg  ) = ang(n,iabeg  )*drad(dyz,n)
                  d2chi(n,dzz,iabeg  ) = ang(n,iabeg  )*drad(dzz,n)
                  d2chi(n,dxx,iabeg+1) = ang(n,iabeg+1)*drad(dxx,n)
                  d2chi(n,dxy,iabeg+1) = ang(n,iabeg+1)*drad(dxy,n)
                  d2chi(n,dxz,iabeg+1) = ang(n,iabeg+1)*drad(dxz,n)
                  d2chi(n,dyy,iabeg+1) = ang(n,iabeg+1)*drad(dyy,n)
                  d2chi(n,dyz,iabeg+1) = ang(n,iabeg+1)*drad(dyz,n)
                  d2chi(n,dzz,iabeg+1) = ang(n,iabeg+1)*drad(dzz,n)
c
                  d3chi(n,dxxx,iabeg  ) = ang(n,iabeg  )*drad(dxxx,n)
                  d3chi(n,dxxy,iabeg  ) = ang(n,iabeg  )*drad(dxxy,n)
                  d3chi(n,dxxz,iabeg  ) = ang(n,iabeg  )*drad(dxxz,n)
                  d3chi(n,dxyy,iabeg  ) = ang(n,iabeg  )*drad(dxyy,n)
                  d3chi(n,dxyz,iabeg  ) = ang(n,iabeg  )*drad(dxyz,n)
                  d3chi(n,dxzz,iabeg  ) = ang(n,iabeg  )*drad(dxzz,n)
                  d3chi(n,dyyy,iabeg  ) = ang(n,iabeg  )*drad(dyyy,n)
                  d3chi(n,dyyz,iabeg  ) = ang(n,iabeg  )*drad(dyyz,n)
                  d3chi(n,dyzz,iabeg  ) = ang(n,iabeg  )*drad(dyzz,n)
                  d3chi(n,dzzz,iabeg  ) = ang(n,iabeg  )*drad(dzzz,n)
                  d3chi(n,dxxx,iabeg+1) = ang(n,iabeg+1)*drad(dxxx,n)
                  d3chi(n,dxxy,iabeg+1) = ang(n,iabeg+1)*drad(dxxy,n)
                  d3chi(n,dxxz,iabeg+1) = ang(n,iabeg+1)*drad(dxxz,n)
                  d3chi(n,dxyy,iabeg+1) = ang(n,iabeg+1)*drad(dxyy,n)
                  d3chi(n,dxyz,iabeg+1) = ang(n,iabeg+1)*drad(dxyz,n)
                  d3chi(n,dxzz,iabeg+1) = ang(n,iabeg+1)*drad(dxzz,n)
                  d3chi(n,dyyy,iabeg+1) = ang(n,iabeg+1)*drad(dyyy,n)
                  d3chi(n,dyyz,iabeg+1) = ang(n,iabeg+1)*drad(dyyz,n)
                  d3chi(n,dyzz,iabeg+1) = ang(n,iabeg+1)*drad(dyzz,n)
                  d3chi(n,dzzz,iabeg+1) = ang(n,iabeg+1)*drad(dzzz,n)
               enddo
            endif
         endif

c     Do the bulk

         do iang = iabeg+irem, iaend, 3
            if (max_der_component .eq. dzero) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(dzero,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(dzero,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(dzero,n)
cbadhere
               enddo
            elseif (max_der_component .eq. dz) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(dzero,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(dzero,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(dzero,n)
c
                  d1chi(n,dx,iang  ) = ang(n,iang  )*drad(dx,n)
                  d1chi(n,dy,iang  ) = ang(n,iang  )*drad(dy,n)
                  d1chi(n,dz,iang  ) = ang(n,iang  )*drad(dz,n)
                  d1chi(n,dx,iang+1) = ang(n,iang+1)*drad(dx,n)
                  d1chi(n,dy,iang+1) = ang(n,iang+1)*drad(dy,n)
                  d1chi(n,dz,iang+1) = ang(n,iang+1)*drad(dz,n)
                  d1chi(n,dx,iang+2) = ang(n,iang+2)*drad(dx,n)
                  d1chi(n,dy,iang+2) = ang(n,iang+2)*drad(dy,n)
                  d1chi(n,dz,iang+2) = ang(n,iang+2)*drad(dz,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(dzero,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(dzero,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(dzero,n)
c
                  d1chi(n,dx,iang  ) = ang(n,iang  )*drad(dx,n)
                  d1chi(n,dy,iang  ) = ang(n,iang  )*drad(dy,n)
                  d1chi(n,dz,iang  ) = ang(n,iang  )*drad(dz,n)
                  d1chi(n,dx,iang+1) = ang(n,iang+1)*drad(dx,n)
                  d1chi(n,dy,iang+1) = ang(n,iang+1)*drad(dy,n)
                  d1chi(n,dz,iang+1) = ang(n,iang+1)*drad(dz,n)
                  d1chi(n,dx,iang+2) = ang(n,iang+2)*drad(dx,n)
                  d1chi(n,dy,iang+2) = ang(n,iang+2)*drad(dy,n)
                  d1chi(n,dz,iang+2) = ang(n,iang+2)*drad(dz,n)
c
                  d2chi(n,dxx,iang  ) = ang(n,iang  )*drad(dxx,n)
                  d2chi(n,dxy,iang  ) = ang(n,iang  )*drad(dxy,n)
                  d2chi(n,dxz,iang  ) = ang(n,iang  )*drad(dxz,n)
                  d2chi(n,dyy,iang  ) = ang(n,iang  )*drad(dyy,n)
                  d2chi(n,dyz,iang  ) = ang(n,iang  )*drad(dyz,n)
                  d2chi(n,dzz,iang  ) = ang(n,iang  )*drad(dzz,n)
                  d2chi(n,dxx,iang+1) = ang(n,iang+1)*drad(dxx,n)
                  d2chi(n,dxy,iang+1) = ang(n,iang+1)*drad(dxy,n)
                  d2chi(n,dxz,iang+1) = ang(n,iang+1)*drad(dxz,n)
                  d2chi(n,dyy,iang+1) = ang(n,iang+1)*drad(dyy,n)
                  d2chi(n,dyz,iang+1) = ang(n,iang+1)*drad(dyz,n)
                  d2chi(n,dzz,iang+1) = ang(n,iang+1)*drad(dzz,n)
                  d2chi(n,dxx,iang+2) = ang(n,iang+2)*drad(dxx,n)
                  d2chi(n,dxy,iang+2) = ang(n,iang+2)*drad(dxy,n)
                  d2chi(n,dxz,iang+2) = ang(n,iang+2)*drad(dxz,n)
                  d2chi(n,dyy,iang+2) = ang(n,iang+2)*drad(dyy,n)
                  d2chi(n,dyz,iang+2) = ang(n,iang+2)*drad(dyz,n)
                  d2chi(n,dzz,iang+2) = ang(n,iang+2)*drad(dzz,n)
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(dzero,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(dzero,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(dzero,n)
c
                  d1chi(n,dx,iang  ) = ang(n,iang  )*drad(dx,n)
                  d1chi(n,dy,iang  ) = ang(n,iang  )*drad(dy,n)
                  d1chi(n,dz,iang  ) = ang(n,iang  )*drad(dz,n)
                  d1chi(n,dx,iang+1) = ang(n,iang+1)*drad(dx,n)
                  d1chi(n,dy,iang+1) = ang(n,iang+1)*drad(dy,n)
                  d1chi(n,dz,iang+1) = ang(n,iang+1)*drad(dz,n)
                  d1chi(n,dx,iang+2) = ang(n,iang+2)*drad(dx,n)
                  d1chi(n,dy,iang+2) = ang(n,iang+2)*drad(dy,n)
                  d1chi(n,dz,iang+2) = ang(n,iang+2)*drad(dz,n)
c
                  d2chi(n,dxx,iang  ) = ang(n,iang  )*drad(dxx,n)
                  d2chi(n,dxy,iang  ) = ang(n,iang  )*drad(dxy,n)
                  d2chi(n,dxz,iang  ) = ang(n,iang  )*drad(dxz,n)
                  d2chi(n,dyy,iang  ) = ang(n,iang  )*drad(dyy,n)
                  d2chi(n,dyz,iang  ) = ang(n,iang  )*drad(dyz,n)
                  d2chi(n,dzz,iang  ) = ang(n,iang  )*drad(dzz,n)
                  d2chi(n,dxx,iang+1) = ang(n,iang+1)*drad(dxx,n)
                  d2chi(n,dxy,iang+1) = ang(n,iang+1)*drad(dxy,n)
                  d2chi(n,dxz,iang+1) = ang(n,iang+1)*drad(dxz,n)
                  d2chi(n,dyy,iang+1) = ang(n,iang+1)*drad(dyy,n)
                  d2chi(n,dyz,iang+1) = ang(n,iang+1)*drad(dyz,n)
                  d2chi(n,dzz,iang+1) = ang(n,iang+1)*drad(dzz,n)
                  d2chi(n,dxx,iang+2) = ang(n,iang+2)*drad(dxx,n)
                  d2chi(n,dxy,iang+2) = ang(n,iang+2)*drad(dxy,n)
                  d2chi(n,dxz,iang+2) = ang(n,iang+2)*drad(dxz,n)
                  d2chi(n,dyy,iang+2) = ang(n,iang+2)*drad(dyy,n)
                  d2chi(n,dyz,iang+2) = ang(n,iang+2)*drad(dyz,n)
                  d2chi(n,dzz,iang+2) = ang(n,iang+2)*drad(dzz,n)
c
                  d3chi(n,dxxx,iang  ) = ang(n,iang  )*drad(dxxx,n)
                  d3chi(n,dxxy,iang  ) = ang(n,iang  )*drad(dxxy,n)
                  d3chi(n,dxxz,iang  ) = ang(n,iang  )*drad(dxxz,n)
                  d3chi(n,dxyy,iang  ) = ang(n,iang  )*drad(dxyy,n)
                  d3chi(n,dxyz,iang  ) = ang(n,iang  )*drad(dxyz,n)
                  d3chi(n,dxzz,iang  ) = ang(n,iang  )*drad(dxzz,n)
                  d3chi(n,dyyy,iang  ) = ang(n,iang  )*drad(dyyy,n)
                  d3chi(n,dyyz,iang  ) = ang(n,iang  )*drad(dyyz,n)
                  d3chi(n,dyzz,iang  ) = ang(n,iang  )*drad(dyzz,n)
                  d3chi(n,dzzz,iang  ) = ang(n,iang  )*drad(dzzz,n)
                  d3chi(n,dxxx,iang+1) = ang(n,iang+1)*drad(dxxx,n)
                  d3chi(n,dxxy,iang+1) = ang(n,iang+1)*drad(dxxy,n)
                  d3chi(n,dxxz,iang+1) = ang(n,iang+1)*drad(dxxz,n)
                  d3chi(n,dxyy,iang+1) = ang(n,iang+1)*drad(dxyy,n)
                  d3chi(n,dxyz,iang+1) = ang(n,iang+1)*drad(dxyz,n)
                  d3chi(n,dxzz,iang+1) = ang(n,iang+1)*drad(dxzz,n)
                  d3chi(n,dyyy,iang+1) = ang(n,iang+1)*drad(dyyy,n)
                  d3chi(n,dyyz,iang+1) = ang(n,iang+1)*drad(dyyz,n)
                  d3chi(n,dyzz,iang+1) = ang(n,iang+1)*drad(dyzz,n)
                  d3chi(n,dzzz,iang+1) = ang(n,iang+1)*drad(dzzz,n)
                  d3chi(n,dxxx,iang+2) = ang(n,iang+2)*drad(dxxx,n)
                  d3chi(n,dxxy,iang+2) = ang(n,iang+2)*drad(dxxy,n)
                  d3chi(n,dxxz,iang+2) = ang(n,iang+2)*drad(dxxz,n)
                  d3chi(n,dxyy,iang+2) = ang(n,iang+2)*drad(dxyy,n)
                  d3chi(n,dxyz,iang+2) = ang(n,iang+2)*drad(dxyz,n)
                  d3chi(n,dxzz,iang+2) = ang(n,iang+2)*drad(dxzz,n)
                  d3chi(n,dyyy,iang+2) = ang(n,iang+2)*drad(dyyy,n)
                  d3chi(n,dyyz,iang+2) = ang(n,iang+2)*drad(dyyz,n)
                  d3chi(n,dyzz,iang+2) = ang(n,iang+2)*drad(dyzz,n)
                  d3chi(n,dzzz,iang+2) = ang(n,iang+2)*drad(dzzz,n)
               enddo
            endif
         enddo
      endif
c
      return
      end
c
      subroutine eval_dang_x_rad(chi, d1chi, d2chi, d3chi, drad, ang,
     &     d1ang, d2ang, d3ang, nq, iabeg, iaend, max_der_component,
     &     gen_l)
c
      implicit none
#include "errquit.fh"
      integer nq, iabeg, iaend, max_der_component, gen_l
      double precision chi(nq,iabeg:iaend),
     &     d1chi(nq,dx:dz,iabeg:iaend),
     &     d2chi(nq,dxx:dzz,iabeg:iaend),
     &     d3chi(nq,dxxx:dzzz,iabeg:iaend),
     &     drad(max_der_component,nq), ang(nq,*),
     &     d1ang(nq,dx:dz,*), d2ang(nq,dxx:dzz,*), d3ang(nq,dxxx:dzzz,*)
c
      integer n, irem, iang
      double precision t(Ax:Azz), u(dx:dz,Ax:Azz), v(dxx:dzz,Ax:Az),
     &     s(2)
c
c     Nothing to do if s functions or no derivatives
c
      if (iabeg .eq. 1 .or. max_der_component .eq. dzero) return
c
      if (iabeg .eq. 2) then

c     p functions

         if (max_der_component .eq. dz) then
            do n = 1, nq
               d1chi(n,dx,Ax) = d1chi(n,dx,Ax) + drad(dzero,n)
               d1chi(n,dy,Ay) = d1chi(n,dy,Ay) + drad(dzero,n)
               d1chi(n,dz,Az) = d1chi(n,dz,Az) + drad(dzero,n)
            enddo
         elseif (max_der_component .eq. dzz) then
            do n = 1, nq
               d1chi(n,dx,Ax) = d1chi(n,dx,Ax) + drad(dzero,n)
               d1chi(n,dy,Ay) = d1chi(n,dy,Ay) + drad(dzero,n)
               d1chi(n,dz,Az) = d1chi(n,dz,Az) + drad(dzero,n)
c
               d2chi(n,dxx,Ax) = d2chi(n,dxx,Ax) + drad(dx,n)*2d0
               d2chi(n,dxy,Ax) = d2chi(n,dxy,Ax) + drad(dy,n)
               d2chi(n,dxz,Ax) = d2chi(n,dxz,Ax) + drad(dz,n)
               d2chi(n,dxy,Ay) = d2chi(n,dxy,Ay) + drad(dx,n)
               d2chi(n,dyy,Ay) = d2chi(n,dyy,Ay) + drad(dy,n)*2d0
               d2chi(n,dyz,Ay) = d2chi(n,dyz,Ay) + drad(dz,n)
               d2chi(n,dxz,Az) = d2chi(n,dxz,Az) + drad(dx,n)
               d2chi(n,dyz,Az) = d2chi(n,dyz,Az) + drad(dy,n)
               d2chi(n,dzz,Az) = d2chi(n,dzz,Az) + drad(dz,n)*2d0
            enddo
         elseif (max_der_component .eq. dzzz) then
            do n = 1, nq
               d1chi(n,dx,Ax) = d1chi(n,dx,Ax) + drad(dzero,n)
               d1chi(n,dy,Ay) = d1chi(n,dy,Ay) + drad(dzero,n)
               d1chi(n,dz,Az) = d1chi(n,dz,Az) + drad(dzero,n)
c
               d2chi(n,dxx,Ax) = d2chi(n,dxx,Ax) + drad(dx,n)*2d0
               d2chi(n,dxy,Ax) = d2chi(n,dxy,Ax) + drad(dy,n)
               d2chi(n,dxz,Ax) = d2chi(n,dxz,Ax) + drad(dz,n)
               d2chi(n,dxy,Ay) = d2chi(n,dxy,Ay) + drad(dx,n)
               d2chi(n,dyy,Ay) = d2chi(n,dyy,Ay) + drad(dy,n)*2d0
               d2chi(n,dyz,Ay) = d2chi(n,dyz,Ay) + drad(dz,n)
               d2chi(n,dxz,Az) = d2chi(n,dxz,Az) + drad(dx,n)
               d2chi(n,dyz,Az) = d2chi(n,dyz,Az) + drad(dy,n)
               d2chi(n,dzz,Az) = d2chi(n,dzz,Az) + drad(dz,n)*2d0
c
               d3chi(n,dxxx,Ax) = d3chi(n,dxxx,Ax) + drad(dxx,n)*3d0
               d3chi(n,dxxy,Ax) = d3chi(n,dxxy,Ax) + drad(dxy,n)*2d0
               d3chi(n,dxxz,Ax) = d3chi(n,dxxz,Ax) + drad(dxz,n)*2d0
               d3chi(n,dxyy,Ax) = d3chi(n,dxyy,Ax) + drad(dyy,n)
               d3chi(n,dxyz,Ax) = d3chi(n,dxyz,Ax) + drad(dyz,n)
               d3chi(n,dxzz,Ax) = d3chi(n,dxzz,Ax) + drad(dzz,n)
               d3chi(n,dxxy,Ay) = d3chi(n,dxxy,Ay) + drad(dxx,n)
               d3chi(n,dxyy,Ay) = d3chi(n,dxyy,Ay) + drad(dxy,n)*2d0
               d3chi(n,dxyz,Ay) = d3chi(n,dxyz,Ay) + drad(dxz,n)
               d3chi(n,dyyy,Ay) = d3chi(n,dyyy,Ay) + drad(dyy,n)*3d0
               d3chi(n,dyyz,Ay) = d3chi(n,dyyz,Ay) + drad(dyz,n)*2d0
               d3chi(n,dyzz,Ay) = d3chi(n,dyzz,Ay) + drad(dzz,n)
               d3chi(n,dxxz,Az) = d3chi(n,dxxz,Az) + drad(dxx,n)
               d3chi(n,dxyz,Az) = d3chi(n,dxyz,Az) + drad(dxy,n)
               d3chi(n,dxzz,Az) = d3chi(n,dxzz,Az) + drad(dxz,n)*2d0
               d3chi(n,dyyz,Az) = d3chi(n,dyyz,Az) + drad(dyy,n)
               d3chi(n,dyzz,Az) = d3chi(n,dyzz,Az) + drad(dyz,n)*2d0
               d3chi(n,dzzz,Az) = d3chi(n,dzzz,Az) + drad(dzz,n)*3d0
            enddo
         endif
      elseif (iaend .eq. 9 .and. gen_l .gt. 2) then

c     Spherical d functions

         s(1) = sqrt(3.d0)
         s(2) = 1d0/s(1)
         if (max_der_component .eq. dz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(dzero,n)*s(1)
               t(Ay) = ang(n,Ay)*drad(dzero,n)*s(1)
               t(Az) = ang(n,Az)*drad(dzero,n)*s(1)
c
               d1chi(n,dx,Ad1) = d1chi(n,dx,Ad1) +     t(Ay)
               d1chi(n,dy,Ad1) = d1chi(n,dy,Ad1) +     t(Ax)
               d1chi(n,dy,Ad2) = d1chi(n,dy,Ad2) +     t(Az)
               d1chi(n,dz,Ad2) = d1chi(n,dz,Ad2) +     t(Ay)
               d1chi(n,dx,Ad3) = d1chi(n,dx,Ad3) -     t(Ax)*s(2)
               d1chi(n,dy,Ad3) = d1chi(n,dy,Ad3) -     t(Ay)*s(2)
               d1chi(n,dz,Ad3) = d1chi(n,dz,Ad3) + 2d0*t(Az)*s(2)
               d1chi(n,dx,Ad4) = d1chi(n,dx,Ad4) -     t(Az)
               d1chi(n,dz,Ad4) = d1chi(n,dz,Ad4) -     t(Ax)
               d1chi(n,dx,Ad5) = d1chi(n,dx,Ad5) +     t(Ax)
               d1chi(n,dy,Ad5) = d1chi(n,dy,Ad5) -     t(Ay)
            enddo
         else
            call errquit('d fn gen_l error in eval_dang_x_rad',1,
     &       UNKNOWN_ERR)
         endif
      elseif (iaend .eq. 10) then

c     Cartesian d functions

         if (max_der_component .eq. dz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(dzero,n)
               t(Ay) = ang(n,Ay)*drad(dzero,n)
               t(Az) = ang(n,Az)*drad(dzero,n)
c
               d1chi(n,dx,Axx) = d1chi(n,dx,Axx) + 2d0*t(Ax)
               d1chi(n,dx,Axy) = d1chi(n,dx,Axy) +     t(Ay)
               d1chi(n,dy,Axy) = d1chi(n,dy,Axy) +     t(Ax)
               d1chi(n,dx,Axz) = d1chi(n,dx,Axz) +     t(Az)
               d1chi(n,dz,Axz) = d1chi(n,dz,Axz) +     t(Ax)
               d1chi(n,dy,Ayy) = d1chi(n,dy,Ayy) + 2d0*t(Ay)
               d1chi(n,dy,Ayz) = d1chi(n,dy,Ayz) +     t(Az)
               d1chi(n,dz,Ayz) = d1chi(n,dz,Ayz) +     t(Ay)
               d1chi(n,dz,Azz) = d1chi(n,dz,Azz) + 2d0*t(Az)
            enddo
         elseif (max_der_component .eq. dzz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(dzero,n)
               t(Ay) = ang(n,Ay)*drad(dzero,n)
               t(Az) = ang(n,Az)*drad(dzero,n)
               u(dx,Ax) = ang(n,Ax)*drad(dx,n)
               u(dy,Ax) = ang(n,Ax)*drad(dy,n)
               u(dz,Ax) = ang(n,Ax)*drad(dz,n)
               u(dx,Ay) = ang(n,Ay)*drad(dx,n)
               u(dy,Ay) = ang(n,Ay)*drad(dy,n)
               u(dz,Ay) = ang(n,Ay)*drad(dz,n)
               u(dx,Az) = ang(n,Az)*drad(dx,n)
               u(dy,Az) = ang(n,Az)*drad(dy,n)
               u(dz,Az) = ang(n,Az)*drad(dz,n)
c
               d1chi(n,dx,Axx) = d1chi(n,dx,Axx) + 2d0*t(Ax)
               d1chi(n,dx,Axy) = d1chi(n,dx,Axy) +     t(Ay)
               d1chi(n,dy,Axy) = d1chi(n,dy,Axy) +     t(Ax)
               d1chi(n,dx,Axz) = d1chi(n,dx,Axz) +     t(Az)
               d1chi(n,dz,Axz) = d1chi(n,dz,Axz) +     t(Ax)
               d1chi(n,dy,Ayy) = d1chi(n,dy,Ayy) + 2d0*t(Ay)
               d1chi(n,dy,Ayz) = d1chi(n,dy,Ayz) +     t(Az)
               d1chi(n,dz,Ayz) = d1chi(n,dz,Ayz) +     t(Ay)
               d1chi(n,dz,Azz) = d1chi(n,dz,Azz) + 2d0*t(Az)
c
               d2chi(n,dxx,Axx) = d2chi(n,dxx,Axx) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dx,Ax)
               d2chi(n,dxy,Axx) = d2chi(n,dxy,Axx) + 2d0*u(dy,Ax)
               d2chi(n,dxz,Axx) = d2chi(n,dxz,Axx) + 2d0*u(dz,Ax)

               d2chi(n,dxx,Axy) = d2chi(n,dxx,Axy) + 2d0*u(dx,Ay)
               d2chi(n,dxy,Axy) = d2chi(n,dxy,Axy) + drad(dzero,n)
     &                                             + u(dy,Ay) + u(dx,Ax)
               d2chi(n,dxz,Axy) = d2chi(n,dxz,Axy) + u(dz,Ay)
               d2chi(n,dyy,Axy) = d2chi(n,dyy,Axy) + 2d0*u(dy,Ax)
               d2chi(n,dyz,Axy) = d2chi(n,dyz,Axy) + u(dz,Ax)

               d2chi(n,dxx,Axz) = d2chi(n,dxx,Axz) + 2d0*u(dx,Az)
               d2chi(n,dxy,Axz) = d2chi(n,dxy,Axz) + u(dy,Az)
               d2chi(n,dxz,Axz) = d2chi(n,dxz,Axz) + drad(dzero,n)
     &                                             + u(dz,Az) + u(dx,Ax)
               d2chi(n,dyz,Axz) = d2chi(n,dyz,Axz) + u(dy,Ax)
               d2chi(n,dzz,Axz) = d2chi(n,dzz,Axz) + 2d0*u(dz,Ax)

               d2chi(n,dxy,Ayy) = d2chi(n,dxy,Ayy) + 2d0*u(dx,Ay)
               d2chi(n,dyy,Ayy) = d2chi(n,dyy,Ayy) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dy,Ay)
               d2chi(n,dyz,Ayy) = d2chi(n,dyz,Ayy) + 2d0*u(dz,Ay)

               d2chi(n,dxy,Ayz) = d2chi(n,dxy,Ayz) + u(dx,Az)
               d2chi(n,dxz,Ayz) = d2chi(n,dxz,Ayz) + u(dx,Ay)
               d2chi(n,dyy,Ayz) = d2chi(n,dyy,Ayz) + 2d0*u(dy,Az)
               d2chi(n,dyz,Ayz) = d2chi(n,dyz,Ayz) + drad(dzero,n)
     &                                             + u(dy,Ay) + u(dz,Az)
               d2chi(n,dzz,Ayz) = d2chi(n,dzz,Ayz) + 2d0*u(dz,Ay)

               d2chi(n,dxz,Azz) = d2chi(n,dxz,Azz) + 2d0*u(dx,Az)
               d2chi(n,dyz,Azz) = d2chi(n,dyz,Azz) + 2d0*u(dy,Az)
               d2chi(n,dzz,Azz) = d2chi(n,dzz,Azz) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dz,Az)
            enddo
         elseif (max_der_component .eq. dzzz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(dzero,n)
               t(Ay) = ang(n,Ay)*drad(dzero,n)
               t(Az) = ang(n,Az)*drad(dzero,n)
               u(dx,Ax) = ang(n,Ax)*drad(dx,n)
               u(dy,Ax) = ang(n,Ax)*drad(dy,n)
               u(dz,Ax) = ang(n,Ax)*drad(dz,n)
               u(dx,Ay) = ang(n,Ay)*drad(dx,n)
               u(dy,Ay) = ang(n,Ay)*drad(dy,n)
               u(dz,Ay) = ang(n,Ay)*drad(dz,n)
               u(dx,Az) = ang(n,Az)*drad(dx,n)
               u(dy,Az) = ang(n,Az)*drad(dy,n)
               u(dz,Az) = ang(n,Az)*drad(dz,n)
               v(dxx,Ax) = ang(n,Ax)*drad(dxx,n)
               v(dxy,Ax) = ang(n,Ax)*drad(dxy,n)
               v(dxz,Ax) = ang(n,Ax)*drad(dxz,n)
               v(dyy,Ax) = ang(n,Ax)*drad(dyy,n)
               v(dyz,Ax) = ang(n,Ax)*drad(dyz,n)
               v(dzz,Ax) = ang(n,Ax)*drad(dzz,n)
               v(dxx,Ay) = ang(n,Ay)*drad(dxx,n)
               v(dxy,Ay) = ang(n,Ay)*drad(dxy,n)
               v(dxz,Ay) = ang(n,Ay)*drad(dxz,n)
               v(dyy,Ay) = ang(n,Ay)*drad(dyy,n)
               v(dyz,Ay) = ang(n,Ay)*drad(dyz,n)
               v(dzz,Ay) = ang(n,Ay)*drad(dzz,n)
               v(dxx,Az) = ang(n,Az)*drad(dxx,n)
               v(dxy,Az) = ang(n,Az)*drad(dxy,n)
               v(dxz,Az) = ang(n,Az)*drad(dxz,n)
               v(dyy,Az) = ang(n,Az)*drad(dyy,n)
               v(dyz,Az) = ang(n,Az)*drad(dyz,n)
               v(dzz,Az) = ang(n,Az)*drad(dzz,n)
c
               d1chi(n,dx,Axx) = d1chi(n,dx,Axx) + 2d0*t(Ax)
               d1chi(n,dx,Axy) = d1chi(n,dx,Axy) +     t(Ay)
               d1chi(n,dy,Axy) = d1chi(n,dy,Axy) +     t(Ax)
               d1chi(n,dx,Axz) = d1chi(n,dx,Axz) +     t(Az)
               d1chi(n,dz,Axz) = d1chi(n,dz,Axz) +     t(Ax)
               d1chi(n,dy,Ayy) = d1chi(n,dy,Ayy) + 2d0*t(Ay)
               d1chi(n,dy,Ayz) = d1chi(n,dy,Ayz) +     t(Az)
               d1chi(n,dz,Ayz) = d1chi(n,dz,Ayz) +     t(Ay)
               d1chi(n,dz,Azz) = d1chi(n,dz,Azz) + 2d0*t(Az)
c
               d2chi(n,dxx,Axx) = d2chi(n,dxx,Axx) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dx,Ax)
               d2chi(n,dxy,Axx) = d2chi(n,dxy,Axx) + 2d0*u(dy,Ax)
               d2chi(n,dxz,Axx) = d2chi(n,dxz,Axx) + 2d0*u(dz,Ax)

               d2chi(n,dxx,Axy) = d2chi(n,dxx,Axy) + 2d0*u(dx,Ay)
               d2chi(n,dxy,Axy) = d2chi(n,dxy,Axy) + drad(dzero,n)
     &                                             + u(dy,Ay) + u(dx,Ax)
               d2chi(n,dxz,Axy) = d2chi(n,dxz,Axy) + u(dz,Ay)
               d2chi(n,dyy,Axy) = d2chi(n,dyy,Axy) + 2d0*u(dy,Ax)
               d2chi(n,dyz,Axy) = d2chi(n,dyz,Axy) + u(dz,Ax)

               d2chi(n,dxx,Axz) = d2chi(n,dxx,Axz) + 2d0*u(dx,Az)
               d2chi(n,dxy,Axz) = d2chi(n,dxy,Axz) + u(dy,Az)
               d2chi(n,dxz,Axz) = d2chi(n,dxz,Axz) + drad(dzero,n)
     &                                             + u(dz,Az) + u(dx,Ax)
               d2chi(n,dyz,Axz) = d2chi(n,dyz,Axz) + u(dy,Ax)
               d2chi(n,dzz,Axz) = d2chi(n,dzz,Axz) + 2d0*u(dz,Ax)

               d2chi(n,dxy,Ayy) = d2chi(n,dxy,Ayy) + 2d0*u(dx,Ay)
               d2chi(n,dyy,Ayy) = d2chi(n,dyy,Ayy) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dy,Ay)
               d2chi(n,dyz,Ayy) = d2chi(n,dyz,Ayy) + 2d0*u(dz,Ay)

               d2chi(n,dxy,Ayz) = d2chi(n,dxy,Ayz) + u(dx,Az)
               d2chi(n,dxz,Ayz) = d2chi(n,dxz,Ayz) + u(dx,Ay)
               d2chi(n,dyy,Ayz) = d2chi(n,dyy,Ayz) + 2d0*u(dy,Az)
               d2chi(n,dyz,Ayz) = d2chi(n,dyz,Ayz) + drad(dzero,n)
     &                                             + u(dy,Ay) + u(dz,Az)
               d2chi(n,dzz,Ayz) = d2chi(n,dzz,Ayz) + 2d0*u(dz,Ay)

               d2chi(n,dxz,Azz) = d2chi(n,dxz,Azz) + 2d0*u(dx,Az)
               d2chi(n,dyz,Azz) = d2chi(n,dyz,Azz) + 2d0*u(dy,Az)
               d2chi(n,dzz,Azz) = d2chi(n,dzz,Azz) + 2d0*drad(dzero,n)
     &                                             + 4d0*u(dz,Az)
c
               d3chi(n,dxxx,Axx) = d3chi(n,dxxx,Axx) + 6d0*drad(dx,n)
     &                                               + 6d0*v(dxx,Ax)
               d3chi(n,dxxy,Axx) = d3chi(n,dxxy,Axx) + 2d0*drad(dy,n)
     &                                               + 4d0*v(dxy,Ax)
               d3chi(n,dxxz,Axx) = d3chi(n,dxxz,Axx) + 2d0*drad(dz,n)
     &                                               + 4d0*v(dxz,Ax)
               d3chi(n,dxyy,Axx) = d3chi(n,dxyy,Axx) + 2d0*v(dyy,Ax)
               d3chi(n,dxyz,Axx) = d3chi(n,dxyz,Axx) + 2d0*v(dyz,Ax)
               d3chi(n,dxzz,Axx) = d3chi(n,dxzz,Axx) + 2d0*v(dzz,Ax)

               d3chi(n,dxxx,Axy) = d3chi(n,dxxx,Axy) + 3d0*v(dxx,Ay)
               d3chi(n,dxxy,Axy) = d3chi(n,dxxy,Axy) + 2d0*drad(dx,n)
     &                               + 2d0*v(dxy,Ay) + v(dxx,Ax)
               d3chi(n,dxxz,Axy) = d3chi(n,dxxz,Axy) + 2d0*v(dxz,Ay)
               d3chi(n,dxyy,Axy) = d3chi(n,dxyy,Axy) + 2d0*drad(dy,n)
     &                               + 2d0*v(dxy,Ax) + v(dyy,Ay)
               d3chi(n,dxyz,Axy) = d3chi(n,dxyz,Axy) + drad(dz,n)
     &                                   + v(dyz,Ay) + v(dxz,Ax)
               d3chi(n,dxzz,Axy) = d3chi(n,dxzz,Axy) + v(dzz,Ay)
               d3chi(n,dyyy,Axy) = d3chi(n,dyyy,Axy) + 3d0*v(dyy,Ax)
               d3chi(n,dyyz,Axy) = d3chi(n,dyyz,Axy) + 2d0*v(dyz,Ax)
               d3chi(n,dyzz,Axy) = d3chi(n,dyzz,Axy) + v(dzz,Ax)

               d3chi(n,dxxx,Axz) = d3chi(n,dxxx,Axz) + 3d0*v(dxx,Az)
               d3chi(n,dxxy,Axz) = d3chi(n,dxxy,Axz) + 2d0*v(dxy,Az)
               d3chi(n,dxxz,Axz) = d3chi(n,dxxz,Axz) + 2d0*drad(dx,n)
     &                               + 2d0*v(dxz,Az) + v(dxx,Ax)
               d3chi(n,dxyy,Axz) = d3chi(n,dxyy,Axz) + v(dyy,Az)
               d3chi(n,dxyz,Axz) = d3chi(n,dxyz,Axz) + drad(dy,n)
     &                                   + v(dyz,Az) + v(dxy,Ax)
               d3chi(n,dxzz,Axz) = d3chi(n,dxzz,Axz) + 2d0*drad(dz,n)
     &                               + 2d0*v(dxz,Ax) + v(dzz,Az)
               d3chi(n,dyyz,Axz) = d3chi(n,dyyz,Axz) + v(dyy,Ax)
               d3chi(n,dyzz,Axz) = d3chi(n,dyzz,Axz) + 2d0*v(dyz,Ax)
               d3chi(n,dzzz,Axz) = d3chi(n,dzzz,Axz) + 3d0*v(dzz,Ax)

               d3chi(n,dxxy,Ayy) = d3chi(n,dxxy,Ayy) + 2d0*v(dxx,Ay)
               d3chi(n,dxyy,Ayy) = d3chi(n,dxyy,Ayy) + 2d0*drad(dx,n)
     &                                               + 4d0*v(dxy,Ay)
               d3chi(n,dxyz,Ayy) = d3chi(n,dxyz,Ayy) + 2d0*v(dxz,Ay)
               d3chi(n,dyyy,Ayy) = d3chi(n,dyyy,Ayy) + 6d0*drad(dy,n)
     &                                               + 6d0*v(dyy,Ay)
               d3chi(n,dyyz,Ayy) = d3chi(n,dyyz,Ayy) + 2d0*drad(dz,n)
     &                                               + 4d0*v(dyz,Ay)
               d3chi(n,dyzz,Ayy) = d3chi(n,dyzz,Ayy) + 2d0*v(dzz,Ay)

               d3chi(n,dxxy,Ayz) = d3chi(n,dxxy,Ayz) + v(dxx,Az)
               d3chi(n,dxxz,Ayz) = d3chi(n,dxxz,Ayz) + v(dxx,Ay)
               d3chi(n,dxyy,Ayz) = d3chi(n,dxyy,Ayz) + 2d0*v(dxy,Az)
               d3chi(n,dxyz,Ayz) = d3chi(n,dxyz,Ayz) + drad(dx,n)
     &                                   + v(dxz,Az) + v(dxy,Ay)
               d3chi(n,dxzz,Ayz) = d3chi(n,dxzz,Ayz) + 2d0*v(dxz,Ay)
               d3chi(n,dyyy,Ayz) = d3chi(n,dyyy,Ayz) + 3d0*v(dyy,Az)
               d3chi(n,dyyz,Ayz) = d3chi(n,dyyz,Ayz) + 2d0*drad(dy,n)
     &                               + 2d0*v(dyz,Az) + v(dyy,Ay)
               d3chi(n,dyzz,Ayz) = d3chi(n,dyzz,Ayz) + 2d0*drad(dz,n)
     &                               + 2d0*v(dyz,Ay) + v(dzz,Az)
               d3chi(n,dzzz,Ayz) = d3chi(n,dzzz,Ayz) + 3d0*v(dzz,Ay)

               d3chi(n,dxxz,Azz) = d3chi(n,dxxz,Azz) + 2d0*v(dxx,Az)
               d3chi(n,dxyz,Azz) = d3chi(n,dxyz,Azz) + 2d0*v(dxy,Az)
               d3chi(n,dxzz,Azz) = d3chi(n,dxzz,Azz) + 2d0*drad(dx,n)
     &                                               + 4d0*v(dxz,Az)
               d3chi(n,dyyz,Azz) = d3chi(n,dyyz,Azz) + 2d0*v(dyy,Az)
               d3chi(n,dyzz,Azz) = d3chi(n,dyzz,Azz) + 2d0*drad(dy,n)
     &                                               + 4d0*v(dyz,Az)
               d3chi(n,dzzz,Azz) = d3chi(n,dzzz,Azz) + 6d0*drad(dz,n)
     &                                               + 6d0*v(dzz,Az)
            enddo
         endif
      elseif (iaend .eq. 16 .and. gen_l .gt. 3) then

c     Spherical f functions

         call errquit('f fn gen_l error in eval_dang_x_rad',1,
     &       UNKNOWN_ERR)
      elseif (iaend .eq. 20 .and. gen_l .gt. 3) then

c     Cartesian f functions

         if (max_der_component .eq. dz) then
            do n = 1, nq
               t(Axx) = ang(n,Axx)*drad(dzero,n)
               t(Axy) = ang(n,Axy)*drad(dzero,n)
               t(Axz) = ang(n,Axz)*drad(dzero,n)
               t(Ayy) = ang(n,Ayy)*drad(dzero,n)
               t(Ayz) = ang(n,Ayz)*drad(dzero,n)
               t(Azz) = ang(n,Azz)*drad(dzero,n)
c
               d1chi(n,dx,Axxx) = d1chi(n,dx,Axxx) + 3d0*t(Axx)
               d1chi(n,dx,Axxy) = d1chi(n,dx,Axxy) + 2d0*t(Axy)
               d1chi(n,dy,Axxy) = d1chi(n,dy,Axxy) +     t(Axx)
               d1chi(n,dx,Axxz) = d1chi(n,dx,Axxz) + 2d0*t(Axz)
               d1chi(n,dz,Axxz) = d1chi(n,dz,Axxz) +     t(Axx)
               d1chi(n,dx,Axyy) = d1chi(n,dx,Axyy) +     t(Ayy)
               d1chi(n,dy,Axyy) = d1chi(n,dy,Axyy) + 2d0*t(Axy)
               d1chi(n,dx,Axyz) = d1chi(n,dx,Axyz) +     t(Ayz)
               d1chi(n,dy,Axyz) = d1chi(n,dy,Axyz) +     t(Axz)
               d1chi(n,dz,Axyz) = d1chi(n,dz,Axyz) +     t(Axy)
               d1chi(n,dx,Axzz) = d1chi(n,dx,Axzz) +     t(Azz)
               d1chi(n,dz,Axzz) = d1chi(n,dz,Axzz) + 2d0*t(Axz)
               d1chi(n,dy,Ayyy) = d1chi(n,dy,Ayyy) + 3d0*t(Ayy)
               d1chi(n,dy,Ayyz) = d1chi(n,dy,Ayyz) + 2d0*t(Ayz)
               d1chi(n,dz,Ayyz) = d1chi(n,dz,Ayyz) +     t(Ayy)
               d1chi(n,dy,Ayzz) = d1chi(n,dy,Ayzz) +     t(Azz)
               d1chi(n,dz,Ayzz) = d1chi(n,dz,Ayzz) + 2d0*t(Ayz)
               d1chi(n,dz,Azzz) = d1chi(n,dz,Azzz) + 3d0*t(Azz)
            enddo
         else
            call errquit('f fn gen_l error in eval_dang_x_rad',1,
     &       UNKNOWN_ERR)
         endif
      else

c     General case

c     Do remainder first

         irem = (iaend-iabeg+1)-((iaend-iabeg+1)/3)*3
         if (irem .eq. 1) then
            if (max_der_component .eq. dz) then
               do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
c
                  d2chi(n,dxx,iabeg  ) = d2chi(n,dxx,iabeg  )
     &                 + d2ang(n,dxx,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg  ) = d2chi(n,dxy,iabeg  )
     &                 + d2ang(n,dxy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dx,n)
                  d2chi(n,dxz,iabeg  ) = d2chi(n,dxz,iabeg  )
     &                 + d2ang(n,dxz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dx,n)
                  d2chi(n,dyy,iabeg  ) = d2chi(n,dyy,iabeg  )
     &                 + d2ang(n,dyy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg  ) = d2chi(n,dyz,iabeg  )
     &                 + d2ang(n,dyz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dy,n)
                  d2chi(n,dzz,iabeg  ) = d2chi(n,dzz,iabeg  )
     &                 + d2ang(n,dzz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dz,n)*2d0
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
c
                  d2chi(n,dxx,iabeg  ) = d2chi(n,dxx,iabeg  )
     &                 + d2ang(n,dxx,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg  ) = d2chi(n,dxy,iabeg  )
     &                 + d2ang(n,dxy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dx,n)
                  d2chi(n,dxz,iabeg  ) = d2chi(n,dxz,iabeg  )
     &                 + d2ang(n,dxz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dx,n)
                  d2chi(n,dyy,iabeg  ) = d2chi(n,dyy,iabeg  )
     &                 + d2ang(n,dyy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg  ) = d2chi(n,dyz,iabeg  )
     &                 + d2ang(n,dyz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dy,n)
                  d2chi(n,dzz,iabeg  ) = d2chi(n,dzz,iabeg  )
     &                 + d2ang(n,dzz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dz,n)*2d0
c
                  d3chi(n,dxxx,iabeg  ) = d3chi(n,dxxx,iabeg  )
     &                 + d3ang(n,dxxx,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iabeg  ) = d3chi(n,dxxy,iabeg  )
     &                 + d3ang(n,dxxy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dy,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iabeg  )*drad(dxx,n)
                  d3chi(n,dxxz,iabeg  ) = d3chi(n,dxxz,iabeg  )
     &                 + d3ang(n,dxxz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iabeg  )*drad(dxx,n)
                  d3chi(n,dxyy,iabeg  ) = d3chi(n,dxyy,iabeg  )
     &                 + d3ang(n,dxyy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dyy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iabeg  ) = d3chi(n,dxyz,iabeg  )
     &                 + d3ang(n,dxyz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dy,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dyz,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dxz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dxy,n)
                  d3chi(n,dxzz,iabeg  ) = d3chi(n,dxzz,iabeg  )
     &                 + d3ang(n,dxzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iabeg  ) = d3chi(n,dyyy,iabeg  )
     &                 + d3ang(n,dyyy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg  )*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iabeg  )*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iabeg  ) = d3chi(n,dyyz,iabeg  )
     &                 + d3ang(n,dyyz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iabeg  )*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iabeg  )*drad(dyy,n)
                  d3chi(n,dyzz,iabeg  ) = d3chi(n,dyzz,iabeg  )
     &                 + d3ang(n,dyzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iabeg  ) = d3chi(n,dzzz,iabeg  )
     &                 + d3ang(n,dzzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dzz,iabeg  )*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iabeg  )*drad(dzz,n)*3d0
               enddo
            endif
         elseif (irem .eq. 2) then
            if (max_der_component .eq. dz) then
               do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
                  d1chi(n,dx,iabeg+1) = d1chi(n,dx,iabeg+1)
     &                 + d1ang(n,dx,iabeg+1)*drad(dzero,n)
                  d1chi(n,dy,iabeg+1) = d1chi(n,dy,iabeg+1)
     &                 + d1ang(n,dy,iabeg+1)*drad(dzero,n)
                  d1chi(n,dz,iabeg+1) = d1chi(n,dz,iabeg+1)
     &                 + d1ang(n,dz,iabeg+1)*drad(dzero,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
                  d1chi(n,dx,iabeg+1) = d1chi(n,dx,iabeg+1)
     &                 + d1ang(n,dx,iabeg+1)*drad(dzero,n)
                  d1chi(n,dy,iabeg+1) = d1chi(n,dy,iabeg+1)
     &                 + d1ang(n,dy,iabeg+1)*drad(dzero,n)
                  d1chi(n,dz,iabeg+1) = d1chi(n,dz,iabeg+1)
     &                 + d1ang(n,dz,iabeg+1)*drad(dzero,n)
c
                  d2chi(n,dxx,iabeg  ) = d2chi(n,dxx,iabeg  )
     &                 + d2ang(n,dxx,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg  ) = d2chi(n,dxy,iabeg  )
     &                 + d2ang(n,dxy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dx,n)
                  d2chi(n,dxz,iabeg  ) = d2chi(n,dxz,iabeg  )
     &                 + d2ang(n,dxz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dx,n)
                  d2chi(n,dyy,iabeg  ) = d2chi(n,dyy,iabeg  )
     &                 + d2ang(n,dyy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg  ) = d2chi(n,dyz,iabeg  )
     &                 + d2ang(n,dyz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dy,n)
                  d2chi(n,dzz,iabeg  ) = d2chi(n,dzz,iabeg  )
     &                 + d2ang(n,dzz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dz,n)*2d0
                  d2chi(n,dxx,iabeg+1) = d2chi(n,dxx,iabeg+1)
     &                 + d2ang(n,dxx,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg+1) = d2chi(n,dxy,iabeg+1)
     &                 + d2ang(n,dxy,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dy,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dx,n)
                  d2chi(n,dxz,iabeg+1) = d2chi(n,dxz,iabeg+1)
     &                 + d2ang(n,dxz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dx,n)
                  d2chi(n,dyy,iabeg+1) = d2chi(n,dyy,iabeg+1)
     &                 + d2ang(n,dyy,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg+1) = d2chi(n,dyz,iabeg+1)
     &                 + d2ang(n,dyz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dy,n)
                  d2chi(n,dzz,iabeg+1) = d2chi(n,dzz,iabeg+1)
     &                 + d2ang(n,dzz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dz,n)*2d0
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  d1chi(n,dx,iabeg  ) = d1chi(n,dx,iabeg  )
     &                 + d1ang(n,dx,iabeg  )*drad(dzero,n)
                  d1chi(n,dy,iabeg  ) = d1chi(n,dy,iabeg  )
     &                 + d1ang(n,dy,iabeg  )*drad(dzero,n)
                  d1chi(n,dz,iabeg  ) = d1chi(n,dz,iabeg  )
     &                 + d1ang(n,dz,iabeg  )*drad(dzero,n)
                  d1chi(n,dx,iabeg+1) = d1chi(n,dx,iabeg+1)
     &                 + d1ang(n,dx,iabeg+1)*drad(dzero,n)
                  d1chi(n,dy,iabeg+1) = d1chi(n,dy,iabeg+1)
     &                 + d1ang(n,dy,iabeg+1)*drad(dzero,n)
                  d1chi(n,dz,iabeg+1) = d1chi(n,dz,iabeg+1)
     &                 + d1ang(n,dz,iabeg+1)*drad(dzero,n)
c
                  d2chi(n,dxx,iabeg  ) = d2chi(n,dxx,iabeg  )
     &                 + d2ang(n,dxx,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg  ) = d2chi(n,dxy,iabeg  )
     &                 + d2ang(n,dxy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dx,n)
                  d2chi(n,dxz,iabeg  ) = d2chi(n,dxz,iabeg  )
     &                 + d2ang(n,dxz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dx,n)
                  d2chi(n,dyy,iabeg  ) = d2chi(n,dyy,iabeg  )
     &                 + d2ang(n,dyy,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg  ) = d2chi(n,dyz,iabeg  )
     &                 + d2ang(n,dyz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dy,n)
                  d2chi(n,dzz,iabeg  ) = d2chi(n,dzz,iabeg  )
     &                 + d2ang(n,dzz,iabeg  )*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dz,n)*2d0
                  d2chi(n,dxx,iabeg+1) = d2chi(n,dxx,iabeg+1)
     &                 + d2ang(n,dxx,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dx,n)*2d0
                  d2chi(n,dxy,iabeg+1) = d2chi(n,dxy,iabeg+1)
     &                 + d2ang(n,dxy,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dy,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dx,n)
                  d2chi(n,dxz,iabeg+1) = d2chi(n,dxz,iabeg+1)
     &                 + d2ang(n,dxz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dx,n)
                  d2chi(n,dyy,iabeg+1) = d2chi(n,dyy,iabeg+1)
     &                 + d2ang(n,dyy,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dy,n)*2d0
                  d2chi(n,dyz,iabeg+1) = d2chi(n,dyz,iabeg+1)
     &                 + d2ang(n,dyz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dy,n)
                  d2chi(n,dzz,iabeg+1) = d2chi(n,dzz,iabeg+1)
     &                 + d2ang(n,dzz,iabeg+1)*drad(dzero,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dz,n)*2d0
c
                  d3chi(n,dxxx,iabeg  ) = d3chi(n,dxxx,iabeg  )
     &                 + d3ang(n,dxxx,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iabeg  ) = d3chi(n,dxxy,iabeg  )
     &                 + d3ang(n,dxxy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dy,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iabeg  )*drad(dxx,n)
                  d3chi(n,dxxz,iabeg  ) = d3chi(n,dxxz,iabeg  )
     &                 + d3ang(n,dxxz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg  )*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iabeg  )*drad(dxx,n)
                  d3chi(n,dxyy,iabeg  ) = d3chi(n,dxyy,iabeg  )
     &                 + d3ang(n,dxyy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dyy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iabeg  ) = d3chi(n,dxyz,iabeg  )
     &                 + d3ang(n,dxyz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dy,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dyz,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dxz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dxy,n)
                  d3chi(n,dxzz,iabeg  ) = d3chi(n,dxzz,iabeg  )
     &                 + d3ang(n,dxzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dxz,iabeg  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg  )*drad(dx,n)
     &                 + d1ang(n,dx,iabeg  )*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iabeg  ) = d3chi(n,dyyy,iabeg  )
     &                 + d3ang(n,dyyy,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg  )*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iabeg  )*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iabeg  ) = d3chi(n,dyyz,iabeg  )
     &                 + d3ang(n,dyyz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg  )*drad(dz,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iabeg  )*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iabeg  )*drad(dyy,n)
                  d3chi(n,dyzz,iabeg  ) = d3chi(n,dyzz,iabeg  )
     &                 + d3ang(n,dyzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dyz,iabeg  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg  )*drad(dy,n)
     &                 + d1ang(n,dy,iabeg  )*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg  )*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iabeg  ) = d3chi(n,dzzz,iabeg  )
     &                 + d3ang(n,dzzz,iabeg  )*drad(dzero,n)
     &                 + d2ang(n,dzz,iabeg  )*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iabeg  )*drad(dzz,n)*3d0
                  d3chi(n,dxxx,iabeg+1) = d3chi(n,dxxx,iabeg+1)
     &                 + d3ang(n,dxxx,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg+1)*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iabeg+1)*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iabeg+1) = d3chi(n,dxxy,iabeg+1)
     &                 + d3ang(n,dxxy,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg+1)*drad(dy,n)
     &                 + d2ang(n,dxy,iabeg+1)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg+1)*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iabeg+1)*drad(dxx,n)
                  d3chi(n,dxxz,iabeg+1) = d3chi(n,dxxz,iabeg+1)
     &                 + d3ang(n,dxxz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iabeg+1)*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg+1)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iabeg+1)*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iabeg+1)*drad(dxx,n)
                  d3chi(n,dxyy,iabeg+1) = d3chi(n,dxyy,iabeg+1)
     &                 + d3ang(n,dxyy,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg+1)*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iabeg+1)*drad(dx,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dyy,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iabeg+1) = d3chi(n,dxyz,iabeg+1)
     &                 + d3ang(n,dxyz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxy,iabeg+1)*drad(dz,n)
     &                 + d2ang(n,dxz,iabeg+1)*drad(dy,n)
     &                 + d2ang(n,dyz,iabeg+1)*drad(dx,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dyz,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dxz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dxy,n)
                  d3chi(n,dxzz,iabeg+1) = d3chi(n,dxzz,iabeg+1)
     &                 + d3ang(n,dxzz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dxz,iabeg+1)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg+1)*drad(dx,n)
     &                 + d1ang(n,dx,iabeg+1)*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iabeg+1) = d3chi(n,dyyy,iabeg+1)
     &                 + d3ang(n,dyyy,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg+1)*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iabeg+1)*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iabeg+1) = d3chi(n,dyyz,iabeg+1)
     &                 + d3ang(n,dyyz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dyy,iabeg+1)*drad(dz,n)
     &                 + d2ang(n,dyz,iabeg+1)*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iabeg+1)*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iabeg+1)*drad(dyy,n)
                  d3chi(n,dyzz,iabeg+1) = d3chi(n,dyzz,iabeg+1)
     &                 + d3ang(n,dyzz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dyz,iabeg+1)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iabeg+1)*drad(dy,n)
     &                 + d1ang(n,dy,iabeg+1)*drad(dzz,n)
     &                 + d1ang(n,dz,iabeg+1)*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iabeg+1) = d3chi(n,dzzz,iabeg+1)
     &                 + d3ang(n,dzzz,iabeg+1)*drad(dzero,n)
     &                 + d2ang(n,dzz,iabeg+1)*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iabeg+1)*drad(dzz,n)*3d0
               enddo
            endif
         endif

c     Do the bulk

         do iang = iabeg+irem, iaend, 3
            if (max_der_component .eq. dz) then
               do n = 1, nq
                  d1chi(n,dx,iang  ) = d1chi(n,dx,iang  )
     &                 + d1ang(n,dx,iang  )*drad(dzero,n)
                  d1chi(n,dy,iang  ) = d1chi(n,dy,iang  )
     &                 + d1ang(n,dy,iang  )*drad(dzero,n)
                  d1chi(n,dz,iang  ) = d1chi(n,dz,iang  )
     &                 + d1ang(n,dz,iang  )*drad(dzero,n)
                  d1chi(n,dx,iang+1) = d1chi(n,dx,iang+1)
     &                 + d1ang(n,dx,iang+1)*drad(dzero,n)
                  d1chi(n,dy,iang+1) = d1chi(n,dy,iang+1)
     &                 + d1ang(n,dy,iang+1)*drad(dzero,n)
                  d1chi(n,dz,iang+1) = d1chi(n,dz,iang+1)
     &                 + d1ang(n,dz,iang+1)*drad(dzero,n)
                  d1chi(n,dx,iang+2) = d1chi(n,dx,iang+2)
     &                 + d1ang(n,dx,iang+2)*drad(dzero,n)
                  d1chi(n,dy,iang+2) = d1chi(n,dy,iang+2)
     &                 + d1ang(n,dy,iang+2)*drad(dzero,n)
                  d1chi(n,dz,iang+2) = d1chi(n,dz,iang+2)
     &                 + d1ang(n,dz,iang+2)*drad(dzero,n)
               enddo
            elseif (max_der_component .eq. dzz) then
                do n = 1, nq
                  d1chi(n,dx,iang  ) = d1chi(n,dx,iang  )
     &                 + d1ang(n,dx,iang  )*drad(dzero,n)
                  d1chi(n,dy,iang  ) = d1chi(n,dy,iang  )
     &                 + d1ang(n,dy,iang  )*drad(dzero,n)
                  d1chi(n,dz,iang  ) = d1chi(n,dz,iang  )
     &                 + d1ang(n,dz,iang  )*drad(dzero,n)
                  d1chi(n,dx,iang+1) = d1chi(n,dx,iang+1)
     &                 + d1ang(n,dx,iang+1)*drad(dzero,n)
                  d1chi(n,dy,iang+1) = d1chi(n,dy,iang+1)
     &                 + d1ang(n,dy,iang+1)*drad(dzero,n)
                  d1chi(n,dz,iang+1) = d1chi(n,dz,iang+1)
     &                 + d1ang(n,dz,iang+1)*drad(dzero,n)
                  d1chi(n,dx,iang+2) = d1chi(n,dx,iang+2)
     &                 + d1ang(n,dx,iang+2)*drad(dzero,n)
                  d1chi(n,dy,iang+2) = d1chi(n,dy,iang+2)
     &                 + d1ang(n,dy,iang+2)*drad(dzero,n)
                  d1chi(n,dz,iang+2) = d1chi(n,dz,iang+2)
     &                 + d1ang(n,dz,iang+2)*drad(dzero,n)
c
                  d2chi(n,dxx,iang  ) = d2chi(n,dxx,iang  )
     &                 + d2ang(n,dxx,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iang  ) = d2chi(n,dxy,iang  )
     &                 + d2ang(n,dxy,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dy,n)
     &                 + d1ang(n,dy,iang  )*drad(dx,n)
                  d2chi(n,dxz,iang  ) = d2chi(n,dxz,iang  )
     &                 + d2ang(n,dxz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dz,n)
     &                 + d1ang(n,dz,iang  )*drad(dx,n)
                  d2chi(n,dyy,iang  ) = d2chi(n,dyy,iang  )
     &                 + d2ang(n,dyy,iang  )*drad(dzero,n)
     &                 + d1ang(n,dy,iang  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iang  ) = d2chi(n,dyz,iang  )
     &                 + d2ang(n,dyz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dy,iang  )*drad(dz,n)
     &                 + d1ang(n,dz,iang  )*drad(dy,n)
                  d2chi(n,dzz,iang  ) = d2chi(n,dzz,iang  )
     &                 + d2ang(n,dzz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dz,iang  )*drad(dz,n)*2d0
                  d2chi(n,dxx,iang+1) = d2chi(n,dxx,iang+1)
     &                 + d2ang(n,dxx,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dx,n)*2d0
                  d2chi(n,dxy,iang+1) = d2chi(n,dxy,iang+1)
     &                 + d2ang(n,dxy,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dy,n)
     &                 + d1ang(n,dy,iang+1)*drad(dx,n)
                  d2chi(n,dxz,iang+1) = d2chi(n,dxz,iang+1)
     &                 + d2ang(n,dxz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dx,n)
                  d2chi(n,dyy,iang+1) = d2chi(n,dyy,iang+1)
     &                 + d2ang(n,dyy,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+1)*drad(dy,n)*2d0
                  d2chi(n,dyz,iang+1) = d2chi(n,dyz,iang+1)
     &                 + d2ang(n,dyz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+1)*drad(dz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dy,n)
                  d2chi(n,dzz,iang+1) = d2chi(n,dzz,iang+1)
     &                 + d2ang(n,dzz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dz,iang+1)*drad(dz,n)*2d0
                  d2chi(n,dxx,iang+2) = d2chi(n,dxx,iang+2)
     &                 + d2ang(n,dxx,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dx,n)*2d0
                  d2chi(n,dxy,iang+2) = d2chi(n,dxy,iang+2)
     &                 + d2ang(n,dxy,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dy,n)
     &                 + d1ang(n,dy,iang+2)*drad(dx,n)
                  d2chi(n,dxz,iang+2) = d2chi(n,dxz,iang+2)
     &                 + d2ang(n,dxz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dx,n)
                  d2chi(n,dyy,iang+2) = d2chi(n,dyy,iang+2)
     &                 + d2ang(n,dyy,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+2)*drad(dy,n)*2d0
                  d2chi(n,dyz,iang+2) = d2chi(n,dyz,iang+2)
     &                 + d2ang(n,dyz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+2)*drad(dz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dy,n)
                  d2chi(n,dzz,iang+2) = d2chi(n,dzz,iang+2)
     &                 + d2ang(n,dzz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dz,iang+2)*drad(dz,n)*2d0
               enddo
            elseif (max_der_component .eq. dzzz) then
               do n = 1, nq
                  d1chi(n,dx,iang  ) = d1chi(n,dx,iang  )
     &                 + d1ang(n,dx,iang  )*drad(dzero,n)
                  d1chi(n,dy,iang  ) = d1chi(n,dy,iang  )
     &                 + d1ang(n,dy,iang  )*drad(dzero,n)
                  d1chi(n,dz,iang  ) = d1chi(n,dz,iang  )
     &                 + d1ang(n,dz,iang  )*drad(dzero,n)
                  d1chi(n,dx,iang+1) = d1chi(n,dx,iang+1)
     &                 + d1ang(n,dx,iang+1)*drad(dzero,n)
                  d1chi(n,dy,iang+1) = d1chi(n,dy,iang+1)
     &                 + d1ang(n,dy,iang+1)*drad(dzero,n)
                  d1chi(n,dz,iang+1) = d1chi(n,dz,iang+1)
     &                 + d1ang(n,dz,iang+1)*drad(dzero,n)
                  d1chi(n,dx,iang+2) = d1chi(n,dx,iang+2)
     &                 + d1ang(n,dx,iang+2)*drad(dzero,n)
                  d1chi(n,dy,iang+2) = d1chi(n,dy,iang+2)
     &                 + d1ang(n,dy,iang+2)*drad(dzero,n)
                  d1chi(n,dz,iang+2) = d1chi(n,dz,iang+2)
     &                 + d1ang(n,dz,iang+2)*drad(dzero,n)
c
                  d2chi(n,dxx,iang  ) = d2chi(n,dxx,iang  )
     &                 + d2ang(n,dxx,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dx,n)*2d0
                  d2chi(n,dxy,iang  ) = d2chi(n,dxy,iang  )
     &                 + d2ang(n,dxy,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dy,n)
     &                 + d1ang(n,dy,iang  )*drad(dx,n)
                  d2chi(n,dxz,iang  ) = d2chi(n,dxz,iang  )
     &                 + d2ang(n,dxz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dx,iang  )*drad(dz,n)
     &                 + d1ang(n,dz,iang  )*drad(dx,n)
                  d2chi(n,dyy,iang  ) = d2chi(n,dyy,iang  )
     &                 + d2ang(n,dyy,iang  )*drad(dzero,n)
     &                 + d1ang(n,dy,iang  )*drad(dy,n)*2d0
                  d2chi(n,dyz,iang  ) = d2chi(n,dyz,iang  )
     &                 + d2ang(n,dyz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dy,iang  )*drad(dz,n)
     &                 + d1ang(n,dz,iang  )*drad(dy,n)
                  d2chi(n,dzz,iang  ) = d2chi(n,dzz,iang  )
     &                 + d2ang(n,dzz,iang  )*drad(dzero,n)
     &                 + d1ang(n,dz,iang  )*drad(dz,n)*2d0
                  d2chi(n,dxx,iang+1) = d2chi(n,dxx,iang+1)
     &                 + d2ang(n,dxx,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dx,n)*2d0
                  d2chi(n,dxy,iang+1) = d2chi(n,dxy,iang+1)
     &                 + d2ang(n,dxy,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dy,n)
     &                 + d1ang(n,dy,iang+1)*drad(dx,n)
                  d2chi(n,dxz,iang+1) = d2chi(n,dxz,iang+1)
     &                 + d2ang(n,dxz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+1)*drad(dz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dx,n)
                  d2chi(n,dyy,iang+1) = d2chi(n,dyy,iang+1)
     &                 + d2ang(n,dyy,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+1)*drad(dy,n)*2d0
                  d2chi(n,dyz,iang+1) = d2chi(n,dyz,iang+1)
     &                 + d2ang(n,dyz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+1)*drad(dz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dy,n)
                  d2chi(n,dzz,iang+1) = d2chi(n,dzz,iang+1)
     &                 + d2ang(n,dzz,iang+1)*drad(dzero,n)
     &                 + d1ang(n,dz,iang+1)*drad(dz,n)*2d0
                  d2chi(n,dxx,iang+2) = d2chi(n,dxx,iang+2)
     &                 + d2ang(n,dxx,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dx,n)*2d0
                  d2chi(n,dxy,iang+2) = d2chi(n,dxy,iang+2)
     &                 + d2ang(n,dxy,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dy,n)
     &                 + d1ang(n,dy,iang+2)*drad(dx,n)
                  d2chi(n,dxz,iang+2) = d2chi(n,dxz,iang+2)
     &                 + d2ang(n,dxz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dx,iang+2)*drad(dz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dx,n)
                  d2chi(n,dyy,iang+2) = d2chi(n,dyy,iang+2)
     &                 + d2ang(n,dyy,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+2)*drad(dy,n)*2d0
                  d2chi(n,dyz,iang+2) = d2chi(n,dyz,iang+2)
     &                 + d2ang(n,dyz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dy,iang+2)*drad(dz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dy,n)
                  d2chi(n,dzz,iang+2) = d2chi(n,dzz,iang+2)
     &                 + d2ang(n,dzz,iang+2)*drad(dzero,n)
     &                 + d1ang(n,dz,iang+2)*drad(dz,n)*2d0
c
                  d3chi(n,dxxx,iang  ) = d3chi(n,dxxx,iang  )
     &                 + d3ang(n,dxxx,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iang  )*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iang  )*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iang  ) = d3chi(n,dxxy,iang  )
     &                 + d3ang(n,dxxy,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iang  )*drad(dy,n)
     &                 + d2ang(n,dxy,iang  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang  )*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iang  )*drad(dxx,n)
                  d3chi(n,dxxz,iang  ) = d3chi(n,dxxz,iang  )
     &                 + d3ang(n,dxxz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxx,iang  )*drad(dz,n)
     &                 + d2ang(n,dxz,iang  )*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang  )*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iang  )*drad(dxx,n)
                  d3chi(n,dxyy,iang  ) = d3chi(n,dxyy,iang  )
     &                 + d3ang(n,dxyy,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iang  )*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iang  )*drad(dx,n)
     &                 + d1ang(n,dx,iang  )*drad(dyy,n)
     &                 + d1ang(n,dy,iang  )*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iang  ) = d3chi(n,dxyz,iang  )
     &                 + d3ang(n,dxyz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxy,iang  )*drad(dz,n)
     &                 + d2ang(n,dxz,iang  )*drad(dy,n)
     &                 + d2ang(n,dyz,iang  )*drad(dx,n)
     &                 + d1ang(n,dx,iang  )*drad(dyz,n)
     &                 + d1ang(n,dy,iang  )*drad(dxz,n)
     &                 + d1ang(n,dz,iang  )*drad(dxy,n)
                  d3chi(n,dxzz,iang  ) = d3chi(n,dxzz,iang  )
     &                 + d3ang(n,dxzz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dxz,iang  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang  )*drad(dx,n)
     &                 + d1ang(n,dx,iang  )*drad(dzz,n)
     &                 + d1ang(n,dz,iang  )*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iang  ) = d3chi(n,dyyy,iang  )
     &                 + d3ang(n,dyyy,iang  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iang  )*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iang  )*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iang  ) = d3chi(n,dyyz,iang  )
     &                 + d3ang(n,dyyz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dyy,iang  )*drad(dz,n)
     &                 + d2ang(n,dyz,iang  )*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iang  )*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iang  )*drad(dyy,n)
                  d3chi(n,dyzz,iang  ) = d3chi(n,dyzz,iang  )
     &                 + d3ang(n,dyzz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dyz,iang  )*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang  )*drad(dy,n)
     &                 + d1ang(n,dy,iang  )*drad(dzz,n)
     &                 + d1ang(n,dz,iang  )*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iang  ) = d3chi(n,dzzz,iang  )
     &                 + d3ang(n,dzzz,iang  )*drad(dzero,n)
     &                 + d2ang(n,dzz,iang  )*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iang  )*drad(dzz,n)*3d0
                  d3chi(n,dxxx,iang+1) = d3chi(n,dxxx,iang+1)
     &                 + d3ang(n,dxxx,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+1)*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iang+1)*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iang+1) = d3chi(n,dxxy,iang+1)
     &                 + d3ang(n,dxxy,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+1)*drad(dy,n)
     &                 + d2ang(n,dxy,iang+1)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang+1)*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iang+1)*drad(dxx,n)
                  d3chi(n,dxxz,iang+1) = d3chi(n,dxxz,iang+1)
     &                 + d3ang(n,dxxz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+1)*drad(dz,n)
     &                 + d2ang(n,dxz,iang+1)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang+1)*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iang+1)*drad(dxx,n)
                  d3chi(n,dxyy,iang+1) = d3chi(n,dxyy,iang+1)
     &                 + d3ang(n,dxyy,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxy,iang+1)*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iang+1)*drad(dx,n)
     &                 + d1ang(n,dx,iang+1)*drad(dyy,n)
     &                 + d1ang(n,dy,iang+1)*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iang+1) = d3chi(n,dxyz,iang+1)
     &                 + d3ang(n,dxyz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxy,iang+1)*drad(dz,n)
     &                 + d2ang(n,dxz,iang+1)*drad(dy,n)
     &                 + d2ang(n,dyz,iang+1)*drad(dx,n)
     &                 + d1ang(n,dx,iang+1)*drad(dyz,n)
     &                 + d1ang(n,dy,iang+1)*drad(dxz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dxy,n)
                  d3chi(n,dxzz,iang+1) = d3chi(n,dxzz,iang+1)
     &                 + d3ang(n,dxzz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dxz,iang+1)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang+1)*drad(dx,n)
     &                 + d1ang(n,dx,iang+1)*drad(dzz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iang+1) = d3chi(n,dyyy,iang+1)
     &                 + d3ang(n,dyyy,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dyy,iang+1)*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iang+1)*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iang+1) = d3chi(n,dyyz,iang+1)
     &                 + d3ang(n,dyyz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dyy,iang+1)*drad(dz,n)
     &                 + d2ang(n,dyz,iang+1)*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iang+1)*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iang+1)*drad(dyy,n)
                  d3chi(n,dyzz,iang+1) = d3chi(n,dyzz,iang+1)
     &                 + d3ang(n,dyzz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dyz,iang+1)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang+1)*drad(dy,n)
     &                 + d1ang(n,dy,iang+1)*drad(dzz,n)
     &                 + d1ang(n,dz,iang+1)*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iang+1) = d3chi(n,dzzz,iang+1)
     &                 + d3ang(n,dzzz,iang+1)*drad(dzero,n)
     &                 + d2ang(n,dzz,iang+1)*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iang+1)*drad(dzz,n)*3d0
                  d3chi(n,dxxx,iang+2) = d3chi(n,dxxx,iang+2)
     &                 + d3ang(n,dxxx,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+2)*drad(dx,n)*3d0
     &                 + d1ang(n,dx,iang+2)*drad(dxx,n)*3d0
                  d3chi(n,dxxy,iang+2) = d3chi(n,dxxy,iang+2)
     &                 + d3ang(n,dxxy,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+2)*drad(dy,n)
     &                 + d2ang(n,dxy,iang+2)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang+2)*drad(dxy,n)*2d0
     &                 + d1ang(n,dy,iang+2)*drad(dxx,n)
                  d3chi(n,dxxz,iang+2) = d3chi(n,dxxz,iang+2)
     &                 + d3ang(n,dxxz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxx,iang+2)*drad(dz,n)
     &                 + d2ang(n,dxz,iang+2)*drad(dx,n)*2d0
     &                 + d1ang(n,dx,iang+2)*drad(dxz,n)*2d0
     &                 + d1ang(n,dz,iang+2)*drad(dxx,n)
                  d3chi(n,dxyy,iang+2) = d3chi(n,dxyy,iang+2)
     &                 + d3ang(n,dxyy,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxy,iang+2)*drad(dy,n)*2d0
     &                 + d2ang(n,dyy,iang+2)*drad(dx,n)
     &                 + d1ang(n,dx,iang+2)*drad(dyy,n)
     &                 + d1ang(n,dy,iang+2)*drad(dxy,n)*2d0
                  d3chi(n,dxyz,iang+2) = d3chi(n,dxyz,iang+2)
     &                 + d3ang(n,dxyz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxy,iang+2)*drad(dz,n)
     &                 + d2ang(n,dxz,iang+2)*drad(dy,n)
     &                 + d2ang(n,dyz,iang+2)*drad(dx,n)
     &                 + d1ang(n,dx,iang+2)*drad(dyz,n)
     &                 + d1ang(n,dy,iang+2)*drad(dxz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dxy,n)
                  d3chi(n,dxzz,iang+2) = d3chi(n,dxzz,iang+2)
     &                 + d3ang(n,dxzz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dxz,iang+2)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang+2)*drad(dx,n)
     &                 + d1ang(n,dx,iang+2)*drad(dzz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dxz,n)*2d0
                  d3chi(n,dyyy,iang+2) = d3chi(n,dyyy,iang+2)
     &                 + d3ang(n,dyyy,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dyy,iang+2)*drad(dy,n)*3d0
     &                 + d1ang(n,dy,iang+2)*drad(dyy,n)*3d0
                  d3chi(n,dyyz,iang+2) = d3chi(n,dyyz,iang+2)
     &                 + d3ang(n,dyyz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dyy,iang+2)*drad(dz,n)
     &                 + d2ang(n,dyz,iang+2)*drad(dy,n)*2d0
     &                 + d1ang(n,dy,iang+2)*drad(dyz,n)*2d0
     &                 + d1ang(n,dz,iang+2)*drad(dyy,n)
                  d3chi(n,dyzz,iang+2) = d3chi(n,dyzz,iang+2)
     &                 + d3ang(n,dyzz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dyz,iang+2)*drad(dz,n)*2d0
     &                 + d2ang(n,dzz,iang+2)*drad(dy,n)
     &                 + d1ang(n,dy,iang+2)*drad(dzz,n)
     &                 + d1ang(n,dz,iang+2)*drad(dyz,n)*2d0
                  d3chi(n,dzzz,iang+2) = d3chi(n,dzzz,iang+2)
     &                 + d3ang(n,dzzz,iang+2)*drad(dzero,n)
     &                 + d2ang(n,dzz,iang+2)*drad(dz,n)*3d0
     &                 + d1ang(n,dz,iang+2)*drad(dzz,n)*3d0
               enddo
            endif
         enddo
      endif
c
      return
      end
      subroutine eval_radial0(crad, r2, nq, zprim, ccoef,
     *     nprimo, ncontr)
      implicit none
      integer nq, nprimo, ncontr
      double precision crad(nq), r2(nq),
     *     zprim(nprimo), ccoef(nprimo,ncontr),t0
c
      integer n, iprimo
      double precision zeta
#ifdef CHKUNDFLW
      double precision rmax,undovl,zz
cedo      parameter(undovl=-700d0)
!      parameter(undovl=-160d0)
c 20 digits ln(10)=2.3025
      parameter(undovl=-20d0*2.3025d0)
#endif
c
c
c     Store results for first primitive directly and then sum the rest
c
      call dcopy(nq,0d0,0,crad,1) 
#ifdef CHKUNDFLW
      rmax=r2(1)
      do n=2,nq
         rmax=max(rmax,r2(n))
      enddo
#endif
      if(nprimo.eq.1) then
         zeta = -zprim(1)
         t0 = ccoef(1,1)
#ifdef CHKUNDFLW
         if(rmax*zeta.lt.undovl) then
            zz=undovl/zeta
            do n = 1,nq
               if(r2(n).lt.zz)
     .              crad(n) = exp(zeta*r2(n))*t0
            enddo
         else
#endif
            do n = 1,nq
               crad(n) = exp(zeta*r2(n))*t0
            enddo
#ifdef CHKUNDFLW
         endif
#endif
      else
         do iprimo = 1, nprimo
            zeta = -zprim(iprimo)
            t0 = ccoef(iprimo,1)
#ifdef CHKUNDFLW
            if(rmax*zeta.lt.undovl) then
               zz=undovl/zeta
               do n = 1,nq
                  if(r2(n).lt.zz)
     .                 crad(n) = crad(n)+exp(zeta*r2(n))*t0
               enddo
            else
#endif
               do n = 1,nq
                  crad(n) = crad(n) + exp(zeta*r2(n))*t0
               enddo
#ifdef CHKUNDFLW
            endif
#endif
         enddo
      endif
      return
      end
C>
C> @}
