      subroutine xc_eval_basis(basis_hand, maxder, chi, d1chi, d2chi,
     &     d3chi, rq, qxyz, ccoord, nq, natoms, iniz, docset,
     &     bas_cent_info, bas_cset_info)
c
c     $Id: xc_eval_basis.F,v 1.3 2001-06-02 00:42:22 edo Exp $
c
c     Evaluate basis functions and their derivatives on a grid
c
c     BGJ - 7/00
c
      implicit none
c
#include "mafdecls.fh"
c!!! Will we still need this?
#include "bas.fh"
c
      integer basis_hand        ! [input] basis set handle
      integer maxder            ! [input] max derivative order to compute
      integer nq                ! [input] number of quadrature points
      integer natoms            ! [input] number of atoms
      double precision chi(nq,*) ! [output] basis function values
      double precision d1chi(nq,3,*) ! [output] basis function 1st derivs
      double precision d2chi(nq,6,*) ! [output] basis function 2nd derivs
      double precision d3chi(nq,10,*) ! [output] basis function 3rd derivs
      double precision rq(nq,natoms) ! [input] squared pt-ctr distances
      double precision qxyz(3,nq) ! [input] grid point coordinates
      double precision ccoord(3,natoms) ! [input] nuclear coordinates
      integer iniz(*)           ! [input] ***???
      logical docset(*)         ! [input] whether cset set is significant
      integer bas_cent_info(3,natoms) ! [input] ***???
      integer bas_cset_info(6,*) ! [input] ***???
c
      external nang_cart, nang_cart_tot
      integer nang_cart, nang_cart_tot
c
c     Local declarations
c
      integer gen_l
      integer MXprim,MXang,MXcont
      integer izprim,iccoef
      integer max_elem_ang_scr,max_der_component,max_ang_for_der
      integer iang,id1ang,id2ang,id3ang,il2k,idrd,irad,itop
      integer iscr,lscr
      integer ictr,icset,ichi
      integer ifirst,ilast,ltype,nprimo,ncontr,nshbf
      integer l,nang,iang_beg,iang_end
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
#define D0 1
#define Dx 2
#define Dy 3
#define Dz 4
#define Dxx 5
#define Dxy 6
#define Dxz 7
#define Dyy 8
#define Dyz 9
#define Dzz 10
#define Dxxx 11
#define Dxxy 12
#define Dxxz 13
#define Dxyy 14
#define Dxyz 15
#define Dxzz 16
#define Dyyy 17
#define Dyyz 18
#define Dyzz 19
#define Dzzz 20
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
         call errquit('Exiting in xc_eval_basis.',1)
      end if
      if( .not. bas_nprim_cn_max(basis_hand,MXprim) )then
         call errquit('Exiting in xc_eval_basis.',2)
      end if
      if( .not. bas_high_angular(basis_hand,MXang) )then
         call errquit('Exiting in xc_eval_basis.',3)
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
     &     itop)
      izprim = iscr
      iccoef = iccoef + iscr
      iang = iang + iscr
      id1ang = id1ang + iscr
      id2ang = id2ang + iscr
      id3ang = id3ang + iscr
      il2k = il2k + iscr
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
     &        dbl_mb(id3ang), qxyz, ccoord(1,ictr), nq, dbl_mb(il2k),
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
            nshbf = ilast-ifirst+1
c
c     Angular momentum
c
            if (ltype .lt. 0) then
               call errquit(
     &              'xc_eval_basis: sp-type orbital not coded',1)
            else
               l = ltype
            endif
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
c     Number of contractions
c
            if (ncontr.ne.1) then
               call errquit(
     &              'xc_eval_basis: general contraction not ready',1)
            endif
c
c     Get exponents and contraction coefficients for this contraction set
c
c!!! These calls could be causing too much overhead - check
            LResult = bas_get_exponent(basis_hand,icset,dbl_mb(izprim))
            LResult = bas_get_coeff(basis_hand,icset,dbl_mb(iccoef))
c
c     Evaluate the radial components of the contracted Gaussian functions
c
            call eval_radial(dbl_mb(irad), rq(1,ictr), nq, maxder,
     &           dbl_mb(izprim), dbl_mb(iccoef), nprimo, ncontr)
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
            ichi = ichi + nshbf
c
 200     continue
 100  continue
c
      if (.not. MA_Pop_Stack(lscr))
     &     call errquit('xc_eval_basis: could not pop stack',0)
c
      return
      end
c
      subroutine eval_radial(crad, r2, nq, maxd, zprim, ccoef,
     *     nprimo, ncontr)
      implicit none
      integer nq, maxd, nprimo, ncontr
      integer nd
      double precision crad(0:maxd,nq), r2(nq),
     *     zprim(nprimo), ccoef(nprimo,ncontr)
c
      integer n, iprimo
      double precision zeta, t(0:3), tmp
c
*      write(6,*)' eval_radial new used '
      if (ncontr .ne. 1) then
         call errquit('Implement ncontr .ne. 1 in eval_radial', ncontr)
      endif
c
c     Store results for first primitive directly and then sum the rest
c
      if (maxd .gt. 3) then
        call errquit('Unexpected in eval_radial', maxd)
      endif
      zeta = -zprim(1)
      t(0) = ccoef(1,1)
      do n = 1, maxd
         t(n) = 2d0*zeta*t(n-1)
      enddo
      do n = 1,nq
        tmp = exp(zeta*r2(n))
        do nd = 0,maxd
          crad(nd,n) = t(nd)*tmp
        enddo
      enddo
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
      subroutine eval_angular(ang, d1ang, d2ang, d3ang, q, ctr, nq,
     &     l2k, maxl, maxder, gen_l, spherical)
      implicit none
c
#include "mafdecls.fh"
c
      integer nq, maxl, maxder, l2k(0:maxl,0:maxl,0:maxl), gen_l
      double precision ang(nq,*), d1ang(nq,Dx:Dz,*),
     &     d2ang(nq,Dxx:Dzz,*), d3ang(nq,Dxxx:Dzzz,*),
     &     q(3,nq), ctr(3), cd(2), cf(7)
      logical spherical
c
      integer n, l, nl, nl1, i1s, i1, ix, iy, i, nx, nxy, lx, ly, lz,
     &     k, k1(Dx:Dz), k2(Dxx:Dzz), k3(Dxxx:Dzzz), iscr, lscr, nscr
      double precision x, y, z,
     &     fac1(Dx:Dz), fac2(Dxx:Dzz), fac3(Dxxx:Dzzz)
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
            cd(1) = sqrt(3d0)
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
            cd(1) = sqrt(3d0)
            cd(2) = cd(1)*0.5d0
            cf(1) = sqrt(2.5d0)*0.5d0
            cf(2) = cf(1)*3d0
            cf(3) = sqrt(15d0)
            cf(4) = sqrt(1.5d0)*0.5d0
            cf(5) = sqrt(6d0)
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
         call errquit('Unexpected in eval_angular', maxl)
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
            call dfill(nq*3*nl, 0d0, d1ang(1,Dx,nl1+1), 1)
            if (maxder .ge. 2) then
               call dfill(nq*6*nl, 0d0, d2ang(1,Dxx,nl1+1), 1)
            endif
            if (maxder .ge. 3) then
               call dfill(nq*10*nl, 0d0, d3ang(1,Dxxx,nl1+1), 1)
            endif
         endif

         do l = gen_l, maxl

            nl  = nang_cart(l)
            nl1 = nang_cart(l-1)
            i1s = nang_cart_tot(l-1)
            ix  = nl1
            iy  = nl - 1
            nx  = nang_cart(l-2)
            nxy = nl1 - nx - 1
            if(.not.spherical) then
c     Evaluate the remaining cartesian functions

            do n = 1, nq
               i1 = i1s
c
               x = ang(n,Ax)
               y = ang(n,Ay)
               z = ang(n,Az)
c
               ang(n,i1+ix) = x*ang(n,i1)
               ang(n,i1+iy) = y*ang(n,i1)
               ang(n,i1+nl) = z*ang(n,i1)
c
               do i = 1, nxy
                  i1 = i1 - 1
                  ang(n,i1+ix) = x*ang(n,i1)
                  ang(n,i1+iy) = y*ang(n,i1)
               enddo
c
               do i = 1, nx
                  i1 = i1 - 1
                  ang(n,i1+ix) = x*ang(n,i1)
               enddo
            enddo
         else
            nl1 = nang_cart_tot(l-1)
            do i=1,((l+1)*(l+2))/2
               call getNxyz(l,i,Nxyz)
               do n = 1,nq

               x = q(1,n) - ctr(1)
               y = q(2,n) - ctr(2)
               z = q(3,n) - ctr(3)

               ang(n,nl1+i) = x**nxyz(1)*y**nxyz(2)*z**nxyz(3)
              enddo
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
                     k1(Dx) = l2k(lx+1,ly,lz)
                     k1(Dy) = l2k(lx,ly+1,lz)
                     k1(Dz) = l2k(lx,ly,lz+1)
                     fac1(Dx) = lx + 1
                     fac1(Dy) = ly + 1
                     fac1(Dz) = lz + 1
                     do n = 1, nq
                        d1ang(n,Dx,k1(Dx)) = fac1(Dx) * ang(n,k)
                        d1ang(n,Dy,k1(Dy)) = fac1(Dy) * ang(n,k)
                        d1ang(n,Dz,k1(Dz)) = fac1(Dz) * ang(n,k)
                     enddo
                  enddo
               enddo
            endif
            if (maxder .ge. 2) then
               do lx = l-2, 0, -1
                  do ly = (l-2) - lx, 0, -1
                     lz = (l-2) - lx - ly
                     k = l2k(lx,ly,lz)
                     k2(Dxx) = l2k(lx+2,ly,lz)
                     k2(Dxy) = l2k(lx+1,ly+1,lz)
                     k2(Dxz) = l2k(lx+1,ly,lz+1)
                     k2(Dyy) = l2k(lx,ly+2,lz)
                     k2(Dyz) = l2k(lx,ly+1,lz+1)
                     k2(Dzz) = l2k(lx,ly,lz+2)
                     fac2(Dxx) = (lx + 1)*(lx + 2)
                     fac2(Dxy) = (lx + 1)*(ly + 1)
                     fac2(Dxz) = (lx + 1)*(lz + 1)
                     fac2(Dyy) = (ly + 1)*(ly + 2)
                     fac2(Dyz) = (ly + 1)*(lz + 1)
                     fac2(Dzz) = (lz + 1)*(lz + 2)
                     do n = 1, nq
                        d2ang(n,Dxx,k2(Dxx)) = fac2(Dxx) * ang(n,k)
                        d2ang(n,Dxy,k2(Dxy)) = fac2(Dxy) * ang(n,k)
                        d2ang(n,Dxz,k2(Dxz)) = fac2(Dxz) * ang(n,k)
                        d2ang(n,Dyy,k2(Dyy)) = fac2(Dyy) * ang(n,k)
                        d2ang(n,Dyz,k2(Dyz)) = fac2(Dyz) * ang(n,k)
                        d2ang(n,Dzz,k2(Dzz)) = fac2(Dzz) * ang(n,k)
                     enddo
                  enddo
               enddo
            endif
            if (maxder .ge. 3) then
               do lx = l-3, 0, -1
                  do ly = (l-3) - lx, 0, -1
                     lz = (l-3) - lx - ly
                     k = l2k(lx,ly,lz)
                     k3(Dxxx) = l2k(lx+3,ly,lz)
                     k3(Dxxy) = l2k(lx+2,ly+1,lz)
                     k3(Dxxz) = l2k(lx+2,ly,lz+1)
                     k3(Dxyy) = l2k(lx+1,ly+2,lz)
                     k3(Dxyz) = l2k(lx+1,ly+1,lz+1)
                     k3(Dxzz) = l2k(lx+1,ly,lz+2)
                     k3(Dyyy) = l2k(lx,ly+3,lz)
                     k3(Dyyz) = l2k(lx,ly+2,lz+1)
                     k3(Dyzz) = l2k(lx,ly+1,lz+2)
                     k3(Dzzz) = l2k(lx,ly,lz+3)
                     fac3(Dxxx) = (lx + 1)*(lx + 2)*(lx + 3)
                     fac3(Dxxy) = (lx + 1)*(lx + 2)*(ly + 1)
                     fac3(Dxxz) = (lx + 1)*(lx + 2)*(lz + 1)
                     fac3(Dxyy) = (lx + 1)*(ly + 1)*(ly + 2)
                     fac3(Dxyz) = (lx + 1)*(ly + 1)*(lz + 1)
                     fac3(Dxzz) = (lx + 1)*(lz + 1)*(lz + 2)
                     fac3(Dyyy) = (ly + 1)*(ly + 2)*(ly + 3)
                     fac3(Dyyz) = (ly + 1)*(ly + 2)*(lz + 1)
                     fac3(Dyzz) = (ly + 1)*(lz + 1)*(lz + 2)
                     fac3(Dzzz) = (lz + 1)*(lz + 2)*(lz + 3)
                     do n = 1, nq
                        d3ang(n,Dxxx,k3(Dxxx)) = fac3(Dxxx) * ang(n,k)
                        d3ang(n,Dxxy,k3(Dxxy)) = fac3(Dxxy) * ang(n,k)
                        d3ang(n,Dxxz,k3(Dxxz)) = fac3(Dxxz) * ang(n,k)
                        d3ang(n,Dxyy,k3(Dxyy)) = fac3(Dxyy) * ang(n,k)
                        d3ang(n,Dxyz,k3(Dxyz)) = fac3(Dxyz) * ang(n,k)
                        d3ang(n,Dxzz,k3(Dxzz)) = fac3(Dxzz) * ang(n,k)
                        d3ang(n,Dyyy,k3(Dyyy)) = fac3(Dyyy) * ang(n,k)
                        d3ang(n,Dyyz,k3(Dyyz)) = fac3(Dyyz) * ang(n,k)
                        d3ang(n,Dyzz,k3(Dyzz)) = fac3(Dyzz) * ang(n,k)
                        d3ang(n,Dzzz,k3(Dzzz)) = fac3(Dzzz) * ang(n,k)
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
     &           'Could not allocate scratch in eval_angular',nscr)

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
                  call spcart_a_s(d1ang(1,Dx,nl1),dbl_mb(iscr),nq*3,l,
     &                 1,.false.,.false.)
                  call dcopy(nq*3*(2*l+1),dbl_mb(iscr),1,
     &                 d1ang(1,Dx,nl),1)
               endif
               if (maxder .ge. 2) then
                  call spcart_a_s(d2ang(1,Dxx,nl1),dbl_mb(iscr),nq*6,l,
     &                 1,.false.,.false.)
                  call dcopy(nq*6*(2*l+1),dbl_mb(iscr),1,
     &                 d2ang(1,Dxx,nl),1)
               endif
               if (maxder .ge. 3) then
                  call spcart_a_s(d3ang(1,Dxxx,nl1),dbl_mb(iscr),nq*10,
     &                 l,1,.false.,.false.)
                  call dcopy(nq*10*(2*l+1),dbl_mb(iscr),1,
     &                 d3ang(1,Dxxx,nl),1)
               endif
 100        continue

            if (.not. MA_Pop_Stack(lscr))
     &           call errquit('eval_angular: could not pop stack',0)
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
      integer nq, maxd, maxcomp
      double precision drad(maxcomp,nq), rad(0:maxd,nq), ang(nq,*)
c
      integer n
      double precision r1, r2, r3, x, y, z
c
      if (maxd .eq. 0) then
         call dcopy(nq, rad, 1, drad, 1)
      elseif (maxd .eq. 1) then
         do n = 1, nq
            drad(D0,n) = rad(0,n)
            r1 = rad(1,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(Dx,n) = x*r1
            drad(Dy,n) = y*r1
            drad(Dz,n) = z*r1
         enddo
      elseif (maxd .eq. 2) then
         do n = 1, nq
            drad(D0,n) = rad(0,n)
            r1 = rad(1,n)
            r2 = rad(2,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(Dx,n) = x*r1
            drad(Dy,n) = y*r1
            drad(Dz,n) = z*r1
            drad(Dxx,n) = x*x*r2 + r1
            drad(Dxy,n) = x*y*r2
            drad(Dxz,n) = x*z*r2
            drad(Dyy,n) = y*y*r2 + r1
            drad(Dyz,n) = y*z*r2
            drad(Dzz,n)= z*z*r2 + r1
         enddo
      elseif (maxd .eq. 3) then
         do n = 1, nq
            drad(D0,n) = rad(0,n)
            r1 = rad(1,n)
            r2 = rad(2,n)
            r3 = rad(3,n)
            x = ang(n,Ax)
            y = ang(n,Ay)
            z = ang(n,Az)
            drad(Dx,n) = x*r1
            drad(Dy,n) = y*r1
            drad(Dz,n) = z*r1
            drad(Dxx,n) = x*x*r2 + r1
            drad(Dxy,n) = x*y*r2
            drad(Dxz,n) = x*z*r2
            drad(Dyy,n) = y*y*r2 + r1
            drad(Dyz,n) = y*z*r2
            drad(Dzz,n)= z*z*r2 + r1
            drad(Dxxx,n)= x*(3d0*r2 + x*x*r3)
            drad(Dxxy,n)= y*(r2 + x*x*r3)
            drad(Dxxz,n)= z*(r2 + x*x*r3)
            drad(Dxyy,n)= x*(r2 + y*y*r3)
            drad(Dxyz,n)= x*y*z*r3
            drad(Dxzz,n)= x*(r2 + z*z*r3)
            drad(Dyyy,n)= y*(3d0*r2 + y*y*r3)
            drad(Dyyz,n)= z*(r2 + y*y*r3)
            drad(Dyzz,n)= y*(r2 + z*z*r3)
            drad(Dzzz,n)= z*(3d0*r2 + z*z*r3)
         enddo
      else
         call errquit('Bad value for maxd in eval_der_radial',maxd)
      endif
c
      return
      end
c
      subroutine eval_ang_x_drad(chi, d1chi, d2chi, d3chi, drad, ang,
     &     nq, iabeg, iaend, max_der_component)
c
      implicit none
      integer nq, iabeg, iaend, max_der_component
      double precision chi(nq,iabeg:iaend),
     &     d1chi(nq,Dx:Dz,iabeg:iaend), d2chi(nq,Dxx:Dzz,iabeg:iaend),
     &     d3chi(nq,Dxxx:Dzzz,iabeg:iaend), drad(max_der_component,nq),
     &     ang(nq,*)
c
      integer n, iang, irem
c
      if (iabeg .eq. 1) then

c     Special case: s functions - just copy

         if (max_der_component .eq. D0) then
            call dcopy(nq, drad, 1, chi, 1)
         elseif (max_der_component .eq. Dz) then
            do n = 1, nq
               chi(n,1) = drad(D0,n)
               d1chi(n,Dx,1) = drad(Dx,n)
               d1chi(n,Dy,1) = drad(Dy,n)
               d1chi(n,Dz,1) = drad(Dz,n)
            enddo
         elseif (max_der_component .eq. Dzz) then
            do n = 1, nq
               chi(n,1) = drad(D0,n)
               d1chi(n,Dx,1) = drad(Dx,n)
               d1chi(n,Dy,1) = drad(Dy,n)
               d1chi(n,Dz,1) = drad(Dz,n)
               d2chi(n,Dxx,1) = drad(Dxx,n)
               d2chi(n,Dxy,1) = drad(Dxy,n)
               d2chi(n,Dxz,1) = drad(Dxz,n)
               d2chi(n,Dyy,1) = drad(Dyy,n)
               d2chi(n,Dyz,1) = drad(Dyz,n)
               d2chi(n,Dzz,1) = drad(Dzz,n)
            enddo
         elseif (max_der_component .eq. Dzzz) then
            do n = 1, nq
               chi(n,1) = drad(D0,n)
               d1chi(n,Dx,1) = drad(Dx,n)
               d1chi(n,Dy,1) = drad(Dy,n)
               d1chi(n,Dz,1) = drad(Dz,n)
               d2chi(n,Dxx,1) = drad(Dxx,n)
               d2chi(n,Dxy,1) = drad(Dxy,n)
               d2chi(n,Dxz,1) = drad(Dxz,n)
               d2chi(n,Dyy,1) = drad(Dyy,n)
               d2chi(n,Dyz,1) = drad(Dyz,n)
               d2chi(n,Dzz,1) = drad(Dzz,n)
               d3chi(n,Dxxx,1) = drad(Dxxx,n)
               d3chi(n,Dxxy,1) = drad(Dxxy,n)
               d3chi(n,Dxxz,1) = drad(Dxxz,n)
               d3chi(n,Dxyy,1) = drad(Dxyy,n)
               d3chi(n,Dxyz,1) = drad(Dxyz,n)
               d3chi(n,Dxzz,1) = drad(Dxzz,n)
               d3chi(n,Dyyy,1) = drad(Dyyy,n)
               d3chi(n,Dyyz,1) = drad(Dyyz,n)
               d3chi(n,Dyzz,1) = drad(Dyzz,n)
               d3chi(n,Dzzz,1) = drad(Dzzz,n)
            enddo
         endif
      else

c     General case: multiply derivatives of radial functions by angular
c     functions

c     Do remainder first

         irem = (iaend-iabeg+1)-((iaend-iabeg+1)/3)*3
         if (irem .eq. 1) then
            if (max_der_component .eq. D0) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
c
                  d2chi(n,Dxx,iabeg  ) = ang(n,iabeg  )*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg  ) = ang(n,iabeg  )*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg  ) = ang(n,iabeg  )*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg  ) = ang(n,iabeg  )*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg  ) = ang(n,iabeg  )*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg  ) = ang(n,iabeg  )*drad(Dzz,n)
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
c
                  d2chi(n,Dxx,iabeg  ) = ang(n,iabeg  )*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg  ) = ang(n,iabeg  )*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg  ) = ang(n,iabeg  )*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg  ) = ang(n,iabeg  )*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg  ) = ang(n,iabeg  )*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg  ) = ang(n,iabeg  )*drad(Dzz,n)
c
                  d3chi(n,Dxxx,iabeg  ) = ang(n,iabeg  )*drad(Dxxx,n)
                  d3chi(n,Dxxy,iabeg  ) = ang(n,iabeg  )*drad(Dxxy,n)
                  d3chi(n,Dxxz,iabeg  ) = ang(n,iabeg  )*drad(Dxxz,n)
                  d3chi(n,Dxyy,iabeg  ) = ang(n,iabeg  )*drad(Dxyy,n)
                  d3chi(n,Dxyz,iabeg  ) = ang(n,iabeg  )*drad(Dxyz,n)
                  d3chi(n,Dxzz,iabeg  ) = ang(n,iabeg  )*drad(Dxzz,n)
                  d3chi(n,Dyyy,iabeg  ) = ang(n,iabeg  )*drad(Dyyy,n)
                  d3chi(n,Dyyz,iabeg  ) = ang(n,iabeg  )*drad(Dyyz,n)
                  d3chi(n,Dyzz,iabeg  ) = ang(n,iabeg  )*drad(Dyzz,n)
                  d3chi(n,Dzzz,iabeg  ) = ang(n,iabeg  )*drad(Dzzz,n)
               enddo
            endif
         elseif (irem .eq. 2) then
            if (max_der_component .eq. D0) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
                  d1chi(n,Dx,iabeg+1) = ang(n,iabeg+1)*drad(Dx,n)
                  d1chi(n,Dy,iabeg+1) = ang(n,iabeg+1)*drad(Dy,n)
                  d1chi(n,Dz,iabeg+1) = ang(n,iabeg+1)*drad(Dz,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
                  d1chi(n,Dx,iabeg+1) = ang(n,iabeg+1)*drad(Dx,n)
                  d1chi(n,Dy,iabeg+1) = ang(n,iabeg+1)*drad(Dy,n)
                  d1chi(n,Dz,iabeg+1) = ang(n,iabeg+1)*drad(Dz,n)
c
                  d2chi(n,Dxx,iabeg  ) = ang(n,iabeg  )*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg  ) = ang(n,iabeg  )*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg  ) = ang(n,iabeg  )*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg  ) = ang(n,iabeg  )*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg  ) = ang(n,iabeg  )*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg  ) = ang(n,iabeg  )*drad(Dzz,n)
                  d2chi(n,Dxx,iabeg+1) = ang(n,iabeg+1)*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg+1) = ang(n,iabeg+1)*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg+1) = ang(n,iabeg+1)*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg+1) = ang(n,iabeg+1)*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg+1) = ang(n,iabeg+1)*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg+1) = ang(n,iabeg+1)*drad(Dzz,n)
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  chi(n,iabeg  ) = ang(n,iabeg  )*drad(D0,n)
                  chi(n,iabeg+1) = ang(n,iabeg+1)*drad(D0,n)
c
                  d1chi(n,Dx,iabeg  ) = ang(n,iabeg  )*drad(Dx,n)
                  d1chi(n,Dy,iabeg  ) = ang(n,iabeg  )*drad(Dy,n)
                  d1chi(n,Dz,iabeg  ) = ang(n,iabeg  )*drad(Dz,n)
                  d1chi(n,Dx,iabeg+1) = ang(n,iabeg+1)*drad(Dx,n)
                  d1chi(n,Dy,iabeg+1) = ang(n,iabeg+1)*drad(Dy,n)
                  d1chi(n,Dz,iabeg+1) = ang(n,iabeg+1)*drad(Dz,n)
c
                  d2chi(n,Dxx,iabeg  ) = ang(n,iabeg  )*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg  ) = ang(n,iabeg  )*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg  ) = ang(n,iabeg  )*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg  ) = ang(n,iabeg  )*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg  ) = ang(n,iabeg  )*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg  ) = ang(n,iabeg  )*drad(Dzz,n)
                  d2chi(n,Dxx,iabeg+1) = ang(n,iabeg+1)*drad(Dxx,n)
                  d2chi(n,Dxy,iabeg+1) = ang(n,iabeg+1)*drad(Dxy,n)
                  d2chi(n,Dxz,iabeg+1) = ang(n,iabeg+1)*drad(Dxz,n)
                  d2chi(n,Dyy,iabeg+1) = ang(n,iabeg+1)*drad(Dyy,n)
                  d2chi(n,Dyz,iabeg+1) = ang(n,iabeg+1)*drad(Dyz,n)
                  d2chi(n,Dzz,iabeg+1) = ang(n,iabeg+1)*drad(Dzz,n)
c
                  d3chi(n,Dxxx,iabeg  ) = ang(n,iabeg  )*drad(Dxxx,n)
                  d3chi(n,Dxxy,iabeg  ) = ang(n,iabeg  )*drad(Dxxy,n)
                  d3chi(n,Dxxz,iabeg  ) = ang(n,iabeg  )*drad(Dxxz,n)
                  d3chi(n,Dxyy,iabeg  ) = ang(n,iabeg  )*drad(Dxyy,n)
                  d3chi(n,Dxyz,iabeg  ) = ang(n,iabeg  )*drad(Dxyz,n)
                  d3chi(n,Dxzz,iabeg  ) = ang(n,iabeg  )*drad(Dxzz,n)
                  d3chi(n,Dyyy,iabeg  ) = ang(n,iabeg  )*drad(Dyyy,n)
                  d3chi(n,Dyyz,iabeg  ) = ang(n,iabeg  )*drad(Dyyz,n)
                  d3chi(n,Dyzz,iabeg  ) = ang(n,iabeg  )*drad(Dyzz,n)
                  d3chi(n,Dzzz,iabeg  ) = ang(n,iabeg  )*drad(Dzzz,n)
                  d3chi(n,Dxxx,iabeg+1) = ang(n,iabeg+1)*drad(Dxxx,n)
                  d3chi(n,Dxxy,iabeg+1) = ang(n,iabeg+1)*drad(Dxxy,n)
                  d3chi(n,Dxxz,iabeg+1) = ang(n,iabeg+1)*drad(Dxxz,n)
                  d3chi(n,Dxyy,iabeg+1) = ang(n,iabeg+1)*drad(Dxyy,n)
                  d3chi(n,Dxyz,iabeg+1) = ang(n,iabeg+1)*drad(Dxyz,n)
                  d3chi(n,Dxzz,iabeg+1) = ang(n,iabeg+1)*drad(Dxzz,n)
                  d3chi(n,Dyyy,iabeg+1) = ang(n,iabeg+1)*drad(Dyyy,n)
                  d3chi(n,Dyyz,iabeg+1) = ang(n,iabeg+1)*drad(Dyyz,n)
                  d3chi(n,Dyzz,iabeg+1) = ang(n,iabeg+1)*drad(Dyzz,n)
                  d3chi(n,Dzzz,iabeg+1) = ang(n,iabeg+1)*drad(Dzzz,n)
               enddo
            endif
         endif

c     Do the bulk

         do iang = iabeg+irem, iaend, 3
            if (max_der_component .eq. D0) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(D0,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(D0,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dz) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(D0,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(D0,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(D0,n)
c
                  d1chi(n,Dx,iang  ) = ang(n,iang  )*drad(Dx,n)
                  d1chi(n,Dy,iang  ) = ang(n,iang  )*drad(Dy,n)
                  d1chi(n,Dz,iang  ) = ang(n,iang  )*drad(Dz,n)
                  d1chi(n,Dx,iang+1) = ang(n,iang+1)*drad(Dx,n)
                  d1chi(n,Dy,iang+1) = ang(n,iang+1)*drad(Dy,n)
                  d1chi(n,Dz,iang+1) = ang(n,iang+1)*drad(Dz,n)
                  d1chi(n,Dx,iang+2) = ang(n,iang+2)*drad(Dx,n)
                  d1chi(n,Dy,iang+2) = ang(n,iang+2)*drad(Dy,n)
                  d1chi(n,Dz,iang+2) = ang(n,iang+2)*drad(Dz,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(D0,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(D0,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(D0,n)
c
                  d1chi(n,Dx,iang  ) = ang(n,iang  )*drad(Dx,n)
                  d1chi(n,Dy,iang  ) = ang(n,iang  )*drad(Dy,n)
                  d1chi(n,Dz,iang  ) = ang(n,iang  )*drad(Dz,n)
                  d1chi(n,Dx,iang+1) = ang(n,iang+1)*drad(Dx,n)
                  d1chi(n,Dy,iang+1) = ang(n,iang+1)*drad(Dy,n)
                  d1chi(n,Dz,iang+1) = ang(n,iang+1)*drad(Dz,n)
                  d1chi(n,Dx,iang+2) = ang(n,iang+2)*drad(Dx,n)
                  d1chi(n,Dy,iang+2) = ang(n,iang+2)*drad(Dy,n)
                  d1chi(n,Dz,iang+2) = ang(n,iang+2)*drad(Dz,n)
c
                  d2chi(n,Dxx,iang  ) = ang(n,iang  )*drad(Dxx,n)
                  d2chi(n,Dxy,iang  ) = ang(n,iang  )*drad(Dxy,n)
                  d2chi(n,Dxz,iang  ) = ang(n,iang  )*drad(Dxz,n)
                  d2chi(n,Dyy,iang  ) = ang(n,iang  )*drad(Dyy,n)
                  d2chi(n,Dyz,iang  ) = ang(n,iang  )*drad(Dyz,n)
                  d2chi(n,Dzz,iang  ) = ang(n,iang  )*drad(Dzz,n)
                  d2chi(n,Dxx,iang+1) = ang(n,iang+1)*drad(Dxx,n)
                  d2chi(n,Dxy,iang+1) = ang(n,iang+1)*drad(Dxy,n)
                  d2chi(n,Dxz,iang+1) = ang(n,iang+1)*drad(Dxz,n)
                  d2chi(n,Dyy,iang+1) = ang(n,iang+1)*drad(Dyy,n)
                  d2chi(n,Dyz,iang+1) = ang(n,iang+1)*drad(Dyz,n)
                  d2chi(n,Dzz,iang+1) = ang(n,iang+1)*drad(Dzz,n)
                  d2chi(n,Dxx,iang+2) = ang(n,iang+2)*drad(Dxx,n)
                  d2chi(n,Dxy,iang+2) = ang(n,iang+2)*drad(Dxy,n)
                  d2chi(n,Dxz,iang+2) = ang(n,iang+2)*drad(Dxz,n)
                  d2chi(n,Dyy,iang+2) = ang(n,iang+2)*drad(Dyy,n)
                  d2chi(n,Dyz,iang+2) = ang(n,iang+2)*drad(Dyz,n)
                  d2chi(n,Dzz,iang+2) = ang(n,iang+2)*drad(Dzz,n)
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  chi(n,iang  ) = ang(n,iang  )*drad(D0,n)
                  chi(n,iang+1) = ang(n,iang+1)*drad(D0,n)
                  chi(n,iang+2) = ang(n,iang+2)*drad(D0,n)
c
                  d1chi(n,Dx,iang  ) = ang(n,iang  )*drad(Dx,n)
                  d1chi(n,Dy,iang  ) = ang(n,iang  )*drad(Dy,n)
                  d1chi(n,Dz,iang  ) = ang(n,iang  )*drad(Dz,n)
                  d1chi(n,Dx,iang+1) = ang(n,iang+1)*drad(Dx,n)
                  d1chi(n,Dy,iang+1) = ang(n,iang+1)*drad(Dy,n)
                  d1chi(n,Dz,iang+1) = ang(n,iang+1)*drad(Dz,n)
                  d1chi(n,Dx,iang+2) = ang(n,iang+2)*drad(Dx,n)
                  d1chi(n,Dy,iang+2) = ang(n,iang+2)*drad(Dy,n)
                  d1chi(n,Dz,iang+2) = ang(n,iang+2)*drad(Dz,n)
c
                  d2chi(n,Dxx,iang  ) = ang(n,iang  )*drad(Dxx,n)
                  d2chi(n,Dxy,iang  ) = ang(n,iang  )*drad(Dxy,n)
                  d2chi(n,Dxz,iang  ) = ang(n,iang  )*drad(Dxz,n)
                  d2chi(n,Dyy,iang  ) = ang(n,iang  )*drad(Dyy,n)
                  d2chi(n,Dyz,iang  ) = ang(n,iang  )*drad(Dyz,n)
                  d2chi(n,Dzz,iang  ) = ang(n,iang  )*drad(Dzz,n)
                  d2chi(n,Dxx,iang+1) = ang(n,iang+1)*drad(Dxx,n)
                  d2chi(n,Dxy,iang+1) = ang(n,iang+1)*drad(Dxy,n)
                  d2chi(n,Dxz,iang+1) = ang(n,iang+1)*drad(Dxz,n)
                  d2chi(n,Dyy,iang+1) = ang(n,iang+1)*drad(Dyy,n)
                  d2chi(n,Dyz,iang+1) = ang(n,iang+1)*drad(Dyz,n)
                  d2chi(n,Dzz,iang+1) = ang(n,iang+1)*drad(Dzz,n)
                  d2chi(n,Dxx,iang+2) = ang(n,iang+2)*drad(Dxx,n)
                  d2chi(n,Dxy,iang+2) = ang(n,iang+2)*drad(Dxy,n)
                  d2chi(n,Dxz,iang+2) = ang(n,iang+2)*drad(Dxz,n)
                  d2chi(n,Dyy,iang+2) = ang(n,iang+2)*drad(Dyy,n)
                  d2chi(n,Dyz,iang+2) = ang(n,iang+2)*drad(Dyz,n)
                  d2chi(n,Dzz,iang+2) = ang(n,iang+2)*drad(Dzz,n)
c
                  d3chi(n,Dxxx,iang  ) = ang(n,iang  )*drad(Dxxx,n)
                  d3chi(n,Dxxy,iang  ) = ang(n,iang  )*drad(Dxxy,n)
                  d3chi(n,Dxxz,iang  ) = ang(n,iang  )*drad(Dxxz,n)
                  d3chi(n,Dxyy,iang  ) = ang(n,iang  )*drad(Dxyy,n)
                  d3chi(n,Dxyz,iang  ) = ang(n,iang  )*drad(Dxyz,n)
                  d3chi(n,Dxzz,iang  ) = ang(n,iang  )*drad(Dxzz,n)
                  d3chi(n,Dyyy,iang  ) = ang(n,iang  )*drad(Dyyy,n)
                  d3chi(n,Dyyz,iang  ) = ang(n,iang  )*drad(Dyyz,n)
                  d3chi(n,Dyzz,iang  ) = ang(n,iang  )*drad(Dyzz,n)
                  d3chi(n,Dzzz,iang  ) = ang(n,iang  )*drad(Dzzz,n)
                  d3chi(n,Dxxx,iang+1) = ang(n,iang+1)*drad(Dxxx,n)
                  d3chi(n,Dxxy,iang+1) = ang(n,iang+1)*drad(Dxxy,n)
                  d3chi(n,Dxxz,iang+1) = ang(n,iang+1)*drad(Dxxz,n)
                  d3chi(n,Dxyy,iang+1) = ang(n,iang+1)*drad(Dxyy,n)
                  d3chi(n,Dxyz,iang+1) = ang(n,iang+1)*drad(Dxyz,n)
                  d3chi(n,Dxzz,iang+1) = ang(n,iang+1)*drad(Dxzz,n)
                  d3chi(n,Dyyy,iang+1) = ang(n,iang+1)*drad(Dyyy,n)
                  d3chi(n,Dyyz,iang+1) = ang(n,iang+1)*drad(Dyyz,n)
                  d3chi(n,Dyzz,iang+1) = ang(n,iang+1)*drad(Dyzz,n)
                  d3chi(n,Dzzz,iang+1) = ang(n,iang+1)*drad(Dzzz,n)
                  d3chi(n,Dxxx,iang+2) = ang(n,iang+2)*drad(Dxxx,n)
                  d3chi(n,Dxxy,iang+2) = ang(n,iang+2)*drad(Dxxy,n)
                  d3chi(n,Dxxz,iang+2) = ang(n,iang+2)*drad(Dxxz,n)
                  d3chi(n,Dxyy,iang+2) = ang(n,iang+2)*drad(Dxyy,n)
                  d3chi(n,Dxyz,iang+2) = ang(n,iang+2)*drad(Dxyz,n)
                  d3chi(n,Dxzz,iang+2) = ang(n,iang+2)*drad(Dxzz,n)
                  d3chi(n,Dyyy,iang+2) = ang(n,iang+2)*drad(Dyyy,n)
                  d3chi(n,Dyyz,iang+2) = ang(n,iang+2)*drad(Dyyz,n)
                  d3chi(n,Dyzz,iang+2) = ang(n,iang+2)*drad(Dyzz,n)
                  d3chi(n,Dzzz,iang+2) = ang(n,iang+2)*drad(Dzzz,n)
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
      integer nq, iabeg, iaend, max_der_component, gen_l
      double precision chi(nq,iabeg:iaend),
     &     d1chi(nq,Dx:Dz,iabeg:iaend),
     &     d2chi(nq,Dxx:Dzz,iabeg:iaend),
     &     d3chi(nq,Dxxx:Dzzz,iabeg:iaend),
     &     drad(max_der_component,nq), ang(nq,*),
     &     d1ang(nq,Dx:Dz,*), d2ang(nq,Dxx:Dzz,*), d3ang(nq,Dxxx:Dzzz,*)
c
      integer n, irem, iang
      double precision t(Ax:Azz), u(Dx:Dz,Ax:Azz), v(Dxx:Dzz,Ax:Az),
     &     s(2)
c
c     Nothing to do if s functions or no derivatives
c
      if (iabeg .eq. 1 .or. max_der_component .eq. D0) then
         return
      endif
c
      if (iabeg .eq. 2) then

c     p functions

         if (max_der_component .eq. Dz) then
            do n = 1, nq
               d1chi(n,Dx,Ax) = d1chi(n,Dx,Ax) + drad(D0,n)
               d1chi(n,Dy,Ay) = d1chi(n,Dy,Ay) + drad(D0,n)
               d1chi(n,Dz,Az) = d1chi(n,Dz,Az) + drad(D0,n)
            enddo
         elseif (max_der_component .eq. Dzz) then
            do n = 1, nq
               d1chi(n,Dx,Ax) = d1chi(n,Dx,Ax) + drad(D0,n)
               d1chi(n,Dy,Ay) = d1chi(n,Dy,Ay) + drad(D0,n)
               d1chi(n,Dz,Az) = d1chi(n,Dz,Az) + drad(D0,n)
c
               d2chi(n,Dxx,Ax) = d2chi(n,Dxx,Ax) + drad(Dx,n)*2d0
               d2chi(n,Dxy,Ax) = d2chi(n,Dxy,Ax) + drad(Dy,n)
               d2chi(n,Dxz,Ax) = d2chi(n,Dxz,Ax) + drad(Dz,n)
               d2chi(n,Dxy,Ay) = d2chi(n,Dxy,Ay) + drad(Dx,n)
               d2chi(n,Dyy,Ay) = d2chi(n,Dyy,Ay) + drad(Dy,n)*2d0
               d2chi(n,Dyz,Ay) = d2chi(n,Dyz,Ay) + drad(Dz,n)
               d2chi(n,Dxz,Az) = d2chi(n,Dxz,Az) + drad(Dx,n)
               d2chi(n,Dyz,Az) = d2chi(n,Dyz,Az) + drad(Dy,n)
               d2chi(n,Dzz,Az) = d2chi(n,Dzz,Az) + drad(Dz,n)*2d0
            enddo
         elseif (max_der_component .eq. Dzzz) then
            do n = 1, nq
               d1chi(n,Dx,Ax) = d1chi(n,Dx,Ax) + drad(D0,n)
               d1chi(n,Dy,Ay) = d1chi(n,Dy,Ay) + drad(D0,n)
               d1chi(n,Dz,Az) = d1chi(n,Dz,Az) + drad(D0,n)
c
               d2chi(n,Dxx,Ax) = d2chi(n,Dxx,Ax) + drad(Dx,n)*2d0
               d2chi(n,Dxy,Ax) = d2chi(n,Dxy,Ax) + drad(Dy,n)
               d2chi(n,Dxz,Ax) = d2chi(n,Dxz,Ax) + drad(Dz,n)
               d2chi(n,Dxy,Ay) = d2chi(n,Dxy,Ay) + drad(Dx,n)
               d2chi(n,Dyy,Ay) = d2chi(n,Dyy,Ay) + drad(Dy,n)*2d0
               d2chi(n,Dyz,Ay) = d2chi(n,Dyz,Ay) + drad(Dz,n)
               d2chi(n,Dxz,Az) = d2chi(n,Dxz,Az) + drad(Dx,n)
               d2chi(n,Dyz,Az) = d2chi(n,Dyz,Az) + drad(Dy,n)
               d2chi(n,Dzz,Az) = d2chi(n,Dzz,Az) + drad(Dz,n)*2d0
c
               d3chi(n,Dxxx,Ax) = d3chi(n,Dxxx,Ax) + drad(Dxx,n)*3d0
               d3chi(n,Dxxy,Ax) = d3chi(n,Dxxy,Ax) + drad(Dxy,n)*2d0
               d3chi(n,Dxxz,Ax) = d3chi(n,Dxxz,Ax) + drad(Dxz,n)*2d0
               d3chi(n,Dxyy,Ax) = d3chi(n,Dxyy,Ax) + drad(Dyy,n)
               d3chi(n,Dxyz,Ax) = d3chi(n,Dxyz,Ax) + drad(Dyz,n)
               d3chi(n,Dxzz,Ax) = d3chi(n,Dxzz,Ax) + drad(Dzz,n)
               d3chi(n,Dxxy,Ay) = d3chi(n,Dxxy,Ay) + drad(Dxx,n)
               d3chi(n,Dxyy,Ay) = d3chi(n,Dxyy,Ay) + drad(Dxy,n)*2d0
               d3chi(n,Dxyz,Ay) = d3chi(n,Dxyz,Ay) + drad(Dxz,n)
               d3chi(n,Dyyy,Ay) = d3chi(n,Dyyy,Ay) + drad(Dyy,n)*3d0
               d3chi(n,Dyyz,Ay) = d3chi(n,Dyyz,Ay) + drad(Dyz,n)*2d0
               d3chi(n,Dyzz,Ay) = d3chi(n,Dyzz,Ay) + drad(Dzz,n)
               d3chi(n,Dxxz,Az) = d3chi(n,Dxxz,Az) + drad(Dxx,n)
               d3chi(n,Dxyz,Az) = d3chi(n,Dxyz,Az) + drad(Dxy,n)
               d3chi(n,Dxzz,Az) = d3chi(n,Dxzz,Az) + drad(Dxz,n)*2d0
               d3chi(n,Dyyz,Az) = d3chi(n,Dyyz,Az) + drad(Dyy,n)
               d3chi(n,Dyzz,Az) = d3chi(n,Dyzz,Az) + drad(Dyz,n)*2d0
               d3chi(n,Dzzz,Az) = d3chi(n,Dzzz,Az) + drad(Dzz,n)*3d0
            enddo
         endif
      elseif (iaend .eq. 9 .and. gen_l .gt. 2) then

c     Spherical d functions

         s(1) = sqrt(3d0)
         s(2) = 1d0/s(1)
         if (max_der_component .eq. Dz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(D0,n)*s(1)
               t(Ay) = ang(n,Ay)*drad(D0,n)*s(1)
               t(Az) = ang(n,Az)*drad(D0,n)*s(1)
c
               d1chi(n,Dx,Ad1) = d1chi(n,Dx,Ad1) +     t(Ay)
               d1chi(n,Dy,Ad1) = d1chi(n,Dy,Ad1) +     t(Ax)
               d1chi(n,Dy,Ad2) = d1chi(n,Dy,Ad2) +     t(Az)
               d1chi(n,Dz,Ad2) = d1chi(n,Dz,Ad2) +     t(Ay)
               d1chi(n,Dx,Ad3) = d1chi(n,Dx,Ad3) -     t(Ax)*s(2)
               d1chi(n,Dy,Ad3) = d1chi(n,Dy,Ad3) -     t(Ay)*s(2)
               d1chi(n,Dz,Ad3) = d1chi(n,Dz,Ad3) + 2d0*t(Az)*s(2)
               d1chi(n,Dx,Ad4) = d1chi(n,Dx,Ad4) -     t(Az)
               d1chi(n,Dz,Ad4) = d1chi(n,Dz,Ad4) -     t(Ax)
               d1chi(n,Dx,Ad5) = d1chi(n,Dx,Ad5) +     t(Ax)
               d1chi(n,Dy,Ad5) = d1chi(n,Dy,Ad5) -     t(Ay)
            enddo
         else
            call errquit('d fn gen_l error in eval_dang_x_rad',1)
         endif
      elseif (iaend .eq. 10) then

c     Cartesian d functions

         if (max_der_component .eq. Dz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(D0,n)
               t(Ay) = ang(n,Ay)*drad(D0,n)
               t(Az) = ang(n,Az)*drad(D0,n)
c
               d1chi(n,Dx,Axx) = d1chi(n,Dx,Axx) + 2d0*t(Ax)
               d1chi(n,Dx,Axy) = d1chi(n,Dx,Axy) +     t(Ay)
               d1chi(n,Dy,Axy) = d1chi(n,Dy,Axy) +     t(Ax)
               d1chi(n,Dx,Axz) = d1chi(n,Dx,Axz) +     t(Az)
               d1chi(n,Dz,Axz) = d1chi(n,Dz,Axz) +     t(Ax)
               d1chi(n,Dy,Ayy) = d1chi(n,Dy,Ayy) + 2d0*t(Ay)
               d1chi(n,Dy,Ayz) = d1chi(n,Dy,Ayz) +     t(Az)
               d1chi(n,Dz,Ayz) = d1chi(n,Dz,Ayz) +     t(Ay)
               d1chi(n,Dz,Azz) = d1chi(n,Dz,Azz) + 2d0*t(Az)
            enddo
         elseif (max_der_component .eq. Dzz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(D0,n)
               t(Ay) = ang(n,Ay)*drad(D0,n)
               t(Az) = ang(n,Az)*drad(D0,n)
               u(Dx,Ax) = ang(n,Ax)*drad(Dx,n)
               u(Dy,Ax) = ang(n,Ax)*drad(Dy,n)
               u(Dz,Ax) = ang(n,Ax)*drad(Dz,n)
               u(Dx,Ay) = ang(n,Ay)*drad(Dx,n)
               u(Dy,Ay) = ang(n,Ay)*drad(Dy,n)
               u(Dz,Ay) = ang(n,Ay)*drad(Dz,n)
               u(Dx,Az) = ang(n,Az)*drad(Dx,n)
               u(Dy,Az) = ang(n,Az)*drad(Dy,n)
               u(Dz,Az) = ang(n,Az)*drad(Dz,n)
c
               d1chi(n,Dx,Axx) = d1chi(n,Dx,Axx) + 2d0*t(Ax)
               d1chi(n,Dx,Axy) = d1chi(n,Dx,Axy) +     t(Ay)
               d1chi(n,Dy,Axy) = d1chi(n,Dy,Axy) +     t(Ax)
               d1chi(n,Dx,Axz) = d1chi(n,Dx,Axz) +     t(Az)
               d1chi(n,Dz,Axz) = d1chi(n,Dz,Axz) +     t(Ax)
               d1chi(n,Dy,Ayy) = d1chi(n,Dy,Ayy) + 2d0*t(Ay)
               d1chi(n,Dy,Ayz) = d1chi(n,Dy,Ayz) +     t(Az)
               d1chi(n,Dz,Ayz) = d1chi(n,Dz,Ayz) +     t(Ay)
               d1chi(n,Dz,Azz) = d1chi(n,Dz,Azz) + 2d0*t(Az)
c
               d2chi(n,Dxx,Axx) = d2chi(n,Dxx,Axx) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dx,Ax)
               d2chi(n,Dxy,Axx) = d2chi(n,Dxy,Axx) + 2d0*u(Dy,Ax)
               d2chi(n,Dxz,Axx) = d2chi(n,Dxz,Axx) + 2d0*u(Dz,Ax)

               d2chi(n,Dxx,Axy) = d2chi(n,Dxx,Axy) + 2d0*u(Dx,Ay)
               d2chi(n,Dxy,Axy) = d2chi(n,Dxy,Axy) + drad(D0,n)
     &                                             + u(Dy,Ay) + u(Dx,Ax)
               d2chi(n,Dxz,Axy) = d2chi(n,Dxz,Axy) + u(Dz,Ay)
               d2chi(n,Dyy,Axy) = d2chi(n,Dyy,Axy) + 2d0*u(Dy,Ax)
               d2chi(n,Dyz,Axy) = d2chi(n,Dyz,Axy) + u(Dz,Ax)

               d2chi(n,Dxx,Axz) = d2chi(n,Dxx,Axz) + 2d0*u(Dx,Az)
               d2chi(n,Dxy,Axz) = d2chi(n,Dxy,Axz) + u(Dy,Az)
               d2chi(n,Dxz,Axz) = d2chi(n,Dxz,Axz) + drad(D0,n)
     &                                             + u(Dz,Az) + u(Dx,Ax)
               d2chi(n,Dyz,Axz) = d2chi(n,Dyz,Axz) + u(Dy,Ax)
               d2chi(n,Dzz,Axz) = d2chi(n,Dzz,Axz) + 2d0*u(Dz,Ax)

               d2chi(n,Dxy,Ayy) = d2chi(n,Dxy,Ayy) + 2d0*u(Dx,Ay)
               d2chi(n,Dyy,Ayy) = d2chi(n,Dyy,Ayy) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dy,Ay)
               d2chi(n,Dyz,Ayy) = d2chi(n,Dyz,Ayy) + 2d0*u(Dz,Ay)

               d2chi(n,Dxy,Ayz) = d2chi(n,Dxy,Ayz) + u(Dx,Az)
               d2chi(n,Dxz,Ayz) = d2chi(n,Dxz,Ayz) + u(Dx,Ay)
               d2chi(n,Dyy,Ayz) = d2chi(n,Dyy,Ayz) + 2d0*u(Dy,Az)
               d2chi(n,Dyz,Ayz) = d2chi(n,Dyz,Ayz) + drad(D0,n)
     &                                             + u(Dy,Ay) + u(Dz,Az)
               d2chi(n,Dzz,Ayz) = d2chi(n,Dzz,Ayz) + 2d0*u(Dz,Ay)

               d2chi(n,Dxz,Azz) = d2chi(n,Dxz,Azz) + 2d0*u(Dx,Az)
               d2chi(n,Dyz,Azz) = d2chi(n,Dyz,Azz) + 2d0*u(Dy,Az)
               d2chi(n,Dzz,Azz) = d2chi(n,Dzz,Azz) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dz,Az)
            enddo
         elseif (max_der_component .eq. Dzzz) then
            do n = 1, nq
               t(Ax) = ang(n,Ax)*drad(D0,n)
               t(Ay) = ang(n,Ay)*drad(D0,n)
               t(Az) = ang(n,Az)*drad(D0,n)
               u(Dx,Ax) = ang(n,Ax)*drad(Dx,n)
               u(Dy,Ax) = ang(n,Ax)*drad(Dy,n)
               u(Dz,Ax) = ang(n,Ax)*drad(Dz,n)
               u(Dx,Ay) = ang(n,Ay)*drad(Dx,n)
               u(Dy,Ay) = ang(n,Ay)*drad(Dy,n)
               u(Dz,Ay) = ang(n,Ay)*drad(Dz,n)
               u(Dx,Az) = ang(n,Az)*drad(Dx,n)
               u(Dy,Az) = ang(n,Az)*drad(Dy,n)
               u(Dz,Az) = ang(n,Az)*drad(Dz,n)
               v(Dxx,Ax) = ang(n,Ax)*drad(Dxx,n)
               v(Dxy,Ax) = ang(n,Ax)*drad(Dxy,n)
               v(Dxz,Ax) = ang(n,Ax)*drad(Dxz,n)
               v(Dyy,Ax) = ang(n,Ax)*drad(Dyy,n)
               v(Dyz,Ax) = ang(n,Ax)*drad(Dyz,n)
               v(Dzz,Ax) = ang(n,Ax)*drad(Dzz,n)
               v(Dxx,Ay) = ang(n,Ay)*drad(Dxx,n)
               v(Dxy,Ay) = ang(n,Ay)*drad(Dxy,n)
               v(Dxz,Ay) = ang(n,Ay)*drad(Dxz,n)
               v(Dyy,Ay) = ang(n,Ay)*drad(Dyy,n)
               v(Dyz,Ay) = ang(n,Ay)*drad(Dyz,n)
               v(Dzz,Ay) = ang(n,Ay)*drad(Dzz,n)
               v(Dxx,Az) = ang(n,Az)*drad(Dxx,n)
               v(Dxy,Az) = ang(n,Az)*drad(Dxy,n)
               v(Dxz,Az) = ang(n,Az)*drad(Dxz,n)
               v(Dyy,Az) = ang(n,Az)*drad(Dyy,n)
               v(Dyz,Az) = ang(n,Az)*drad(Dyz,n)
               v(Dzz,Az) = ang(n,Az)*drad(Dzz,n)
c
               d1chi(n,Dx,Axx) = d1chi(n,Dx,Axx) + 2d0*t(Ax)
               d1chi(n,Dx,Axy) = d1chi(n,Dx,Axy) +     t(Ay)
               d1chi(n,Dy,Axy) = d1chi(n,Dy,Axy) +     t(Ax)
               d1chi(n,Dx,Axz) = d1chi(n,Dx,Axz) +     t(Az)
               d1chi(n,Dz,Axz) = d1chi(n,Dz,Axz) +     t(Ax)
               d1chi(n,Dy,Ayy) = d1chi(n,Dy,Ayy) + 2d0*t(Ay)
               d1chi(n,Dy,Ayz) = d1chi(n,Dy,Ayz) +     t(Az)
               d1chi(n,Dz,Ayz) = d1chi(n,Dz,Ayz) +     t(Ay)
               d1chi(n,Dz,Azz) = d1chi(n,Dz,Azz) + 2d0*t(Az)
c
               d2chi(n,Dxx,Axx) = d2chi(n,Dxx,Axx) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dx,Ax)
               d2chi(n,Dxy,Axx) = d2chi(n,Dxy,Axx) + 2d0*u(Dy,Ax)
               d2chi(n,Dxz,Axx) = d2chi(n,Dxz,Axx) + 2d0*u(Dz,Ax)

               d2chi(n,Dxx,Axy) = d2chi(n,Dxx,Axy) + 2d0*u(Dx,Ay)
               d2chi(n,Dxy,Axy) = d2chi(n,Dxy,Axy) + drad(D0,n)
     &                                             + u(Dy,Ay) + u(Dx,Ax)
               d2chi(n,Dxz,Axy) = d2chi(n,Dxz,Axy) + u(Dz,Ay)
               d2chi(n,Dyy,Axy) = d2chi(n,Dyy,Axy) + 2d0*u(Dy,Ax)
               d2chi(n,Dyz,Axy) = d2chi(n,Dyz,Axy) + u(Dz,Ax)

               d2chi(n,Dxx,Axz) = d2chi(n,Dxx,Axz) + 2d0*u(Dx,Az)
               d2chi(n,Dxy,Axz) = d2chi(n,Dxy,Axz) + u(Dy,Az)
               d2chi(n,Dxz,Axz) = d2chi(n,Dxz,Axz) + drad(D0,n)
     &                                             + u(Dz,Az) + u(Dx,Ax)
               d2chi(n,Dyz,Axz) = d2chi(n,Dyz,Axz) + u(Dy,Ax)
               d2chi(n,Dzz,Axz) = d2chi(n,Dzz,Axz) + 2d0*u(Dz,Ax)

               d2chi(n,Dxy,Ayy) = d2chi(n,Dxy,Ayy) + 2d0*u(Dx,Ay)
               d2chi(n,Dyy,Ayy) = d2chi(n,Dyy,Ayy) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dy,Ay)
               d2chi(n,Dyz,Ayy) = d2chi(n,Dyz,Ayy) + 2d0*u(Dz,Ay)

               d2chi(n,Dxy,Ayz) = d2chi(n,Dxy,Ayz) + u(Dx,Az)
               d2chi(n,Dxz,Ayz) = d2chi(n,Dxz,Ayz) + u(Dx,Ay)
               d2chi(n,Dyy,Ayz) = d2chi(n,Dyy,Ayz) + 2d0*u(Dy,Az)
               d2chi(n,Dyz,Ayz) = d2chi(n,Dyz,Ayz) + drad(D0,n)
     &                                             + u(Dy,Ay) + u(Dz,Az)
               d2chi(n,Dzz,Ayz) = d2chi(n,Dzz,Ayz) + 2d0*u(Dz,Ay)

               d2chi(n,Dxz,Azz) = d2chi(n,Dxz,Azz) + 2d0*u(Dx,Az)
               d2chi(n,Dyz,Azz) = d2chi(n,Dyz,Azz) + 2d0*u(Dy,Az)
               d2chi(n,Dzz,Azz) = d2chi(n,Dzz,Azz) + 2d0*drad(D0,n)
     &                                             + 4d0*u(Dz,Az)
c
               d3chi(n,Dxxx,Axx) = d3chi(n,Dxxx,Axx) + 6d0*drad(Dx,n)
     &                                               + 6d0*v(Dxx,Ax)
               d3chi(n,Dxxy,Axx) = d3chi(n,Dxxy,Axx) + 2d0*drad(Dy,n)
     &                                               + 4d0*v(Dxy,Ax)
               d3chi(n,Dxxz,Axx) = d3chi(n,Dxxz,Axx) + 2d0*drad(Dz,n)
     &                                               + 4d0*v(Dxz,Ax)
               d3chi(n,Dxyy,Axx) = d3chi(n,Dxyy,Axx) + 2d0*v(Dyy,Ax)
               d3chi(n,Dxyz,Axx) = d3chi(n,Dxyz,Axx) + 2d0*v(Dyz,Ax)
               d3chi(n,Dxzz,Axx) = d3chi(n,Dxzz,Axx) + 2d0*v(Dzz,Ax)

               d3chi(n,Dxxx,Axy) = d3chi(n,Dxxx,Axy) + 3d0*v(Dxx,Ay)
               d3chi(n,Dxxy,Axy) = d3chi(n,Dxxy,Axy) + 2d0*drad(Dx,n)
     &                               + 2d0*v(Dxy,Ay) + v(Dxx,Ax)
               d3chi(n,Dxxz,Axy) = d3chi(n,Dxxz,Axy) + 2d0*v(Dxz,Ay)
               d3chi(n,Dxyy,Axy) = d3chi(n,Dxyy,Axy) + 2d0*drad(Dy,n)
     &                               + 2d0*v(Dxy,Ax) + v(Dyy,Ay)
               d3chi(n,Dxyz,Axy) = d3chi(n,Dxyz,Axy) + drad(Dz,n)
     &                                   + v(Dyz,Ay) + v(Dxz,Ax)
               d3chi(n,Dxzz,Axy) = d3chi(n,Dxzz,Axy) + v(Dzz,Ay)
               d3chi(n,Dyyy,Axy) = d3chi(n,Dyyy,Axy) + 3d0*v(Dyy,Ax)
               d3chi(n,Dyyz,Axy) = d3chi(n,Dyyz,Axy) + 2d0*v(Dyz,Ax)
               d3chi(n,Dyzz,Axy) = d3chi(n,Dyzz,Axy) + v(Dzz,Ax)

               d3chi(n,Dxxx,Axz) = d3chi(n,Dxxx,Axz) + 3d0*v(Dxx,Az)
               d3chi(n,Dxxy,Axz) = d3chi(n,Dxxy,Axz) + 2d0*v(Dxy,Az)
               d3chi(n,Dxxz,Axz) = d3chi(n,Dxxz,Axz) + 2d0*drad(Dx,n)
     &                               + 2d0*v(Dxz,Az) + v(Dxx,Ax)
               d3chi(n,Dxyy,Axz) = d3chi(n,Dxyy,Axz) + v(Dyy,Az)
               d3chi(n,Dxyz,Axz) = d3chi(n,Dxyz,Axz) + drad(Dy,n)
     &                                   + v(Dyz,Az) + v(Dxy,Ax)
               d3chi(n,Dxzz,Axz) = d3chi(n,Dxzz,Axz) + 2d0*drad(Dz,n)
     &                               + 2d0*v(Dxz,Ax) + v(Dzz,Az)
               d3chi(n,Dyyz,Axz) = d3chi(n,Dyyz,Axz) + v(Dyy,Ax)
               d3chi(n,Dyzz,Axz) = d3chi(n,Dyzz,Axz) + 2d0*v(Dyz,Ax)
               d3chi(n,Dzzz,Axz) = d3chi(n,Dzzz,Axz) + 3d0*v(Dzz,Ax)

               d3chi(n,Dxxy,Ayy) = d3chi(n,Dxxy,Ayy) + 2d0*v(Dxx,Ay)
               d3chi(n,Dxyy,Ayy) = d3chi(n,Dxyy,Ayy) + 2d0*drad(Dx,n)
     &                                               + 4d0*v(Dxy,Ay)
               d3chi(n,Dxyz,Ayy) = d3chi(n,Dxyz,Ayy) + 2d0*v(Dxz,Ay)
               d3chi(n,Dyyy,Ayy) = d3chi(n,Dyyy,Ayy) + 6d0*drad(Dy,n)
     &                                               + 6d0*v(Dyy,Ay)
               d3chi(n,Dyyz,Ayy) = d3chi(n,Dyyz,Ayy) + 2d0*drad(Dz,n)
     &                                               + 4d0*v(Dyz,Ay)
               d3chi(n,Dyzz,Ayy) = d3chi(n,Dyzz,Ayy) + 2d0*v(Dzz,Ay)

               d3chi(n,Dxxy,Ayz) = d3chi(n,Dxxy,Ayz) + v(Dxx,Az)
               d3chi(n,Dxxz,Ayz) = d3chi(n,Dxxz,Ayz) + v(Dxx,Ay)
               d3chi(n,Dxyy,Ayz) = d3chi(n,Dxyy,Ayz) + 2d0*v(Dxy,Az)
               d3chi(n,Dxyz,Ayz) = d3chi(n,Dxyz,Ayz) + drad(Dx,n)
     &                                   + v(Dxz,Az) + v(Dxy,Ay)
               d3chi(n,Dxzz,Ayz) = d3chi(n,Dxzz,Ayz) + 2d0*v(Dxz,Ay)
               d3chi(n,Dyyy,Ayz) = d3chi(n,Dyyy,Ayz) + 3d0*v(Dyy,Az)
               d3chi(n,Dyyz,Ayz) = d3chi(n,Dyyz,Ayz) + 2d0*drad(Dy,n)
     &                               + 2d0*v(Dyz,Az) + v(Dyy,Ay)
               d3chi(n,Dyzz,Ayz) = d3chi(n,Dyzz,Ayz) + 2d0*drad(Dz,n)
     &                               + 2d0*v(Dyz,Ay) + v(Dzz,Az)
               d3chi(n,Dzzz,Ayz) = d3chi(n,Dzzz,Ayz) + 3d0*v(Dzz,Ay)

               d3chi(n,Dxxz,Azz) = d3chi(n,Dxxz,Azz) + 2d0*v(Dxx,Az)
               d3chi(n,Dxyz,Azz) = d3chi(n,Dxyz,Azz) + 2d0*v(Dxy,Az)
               d3chi(n,Dxzz,Azz) = d3chi(n,Dxzz,Azz) + 2d0*drad(Dx,n)
     &                                               + 4d0*v(Dxz,Az)
               d3chi(n,Dyyz,Azz) = d3chi(n,Dyyz,Azz) + 2d0*v(Dyy,Az)
               d3chi(n,Dyzz,Azz) = d3chi(n,Dyzz,Azz) + 2d0*drad(Dy,n)
     &                                               + 4d0*v(Dyz,Az)
               d3chi(n,Dzzz,Azz) = d3chi(n,Dzzz,Azz) + 6d0*drad(Dz,n)
     &                                               + 6d0*v(Dzz,Az)
            enddo
         endif
      elseif (iaend .eq. 16 .and. gen_l .gt. 3) then

c     Spherical f functions

         call errquit('f fn gen_l error in eval_dang_x_rad',1)
      elseif (iaend .eq. 20 .and. gen_l .gt. 3) then

c     Cartesian f functions

         if (max_der_component .eq. Dz) then
            do n = 1, nq
               t(Axx) = ang(n,Axx)*drad(D0,n)
               t(Axy) = ang(n,Axy)*drad(D0,n)
               t(Axz) = ang(n,Axz)*drad(D0,n)
               t(Ayy) = ang(n,Ayy)*drad(D0,n)
               t(Ayz) = ang(n,Ayz)*drad(D0,n)
               t(Azz) = ang(n,Azz)*drad(D0,n)
c
               d1chi(n,Dx,Axxx) = d1chi(n,Dx,Axxx) + 3d0*t(Axx)
               d1chi(n,Dx,Axxy) = d1chi(n,Dx,Axxy) + 2d0*t(Axy)
               d1chi(n,Dy,Axxy) = d1chi(n,Dy,Axxy) +     t(Axx)
               d1chi(n,Dx,Axxz) = d1chi(n,Dx,Axxz) + 2d0*t(Axz)
               d1chi(n,Dz,Axxz) = d1chi(n,Dz,Axxz) +     t(Axx)
               d1chi(n,Dx,Axyy) = d1chi(n,Dx,Axyy) +     t(Ayy)
               d1chi(n,Dy,Axyy) = d1chi(n,Dy,Axyy) + 2d0*t(Axy)
               d1chi(n,Dx,Axyz) = d1chi(n,Dx,Axyz) +     t(Ayz)
               d1chi(n,Dy,Axyz) = d1chi(n,Dy,Axyz) +     t(Axz)
               d1chi(n,Dz,Axyz) = d1chi(n,Dz,Axyz) +     t(Axy)
               d1chi(n,Dx,Axzz) = d1chi(n,Dx,Axzz) +     t(Azz)
               d1chi(n,Dz,Axzz) = d1chi(n,Dz,Axzz) + 2d0*t(Axz)
               d1chi(n,Dy,Ayyy) = d1chi(n,Dy,Ayyy) + 3d0*t(Ayy)
               d1chi(n,Dy,Ayyz) = d1chi(n,Dy,Ayyz) + 2d0*t(Ayz)
               d1chi(n,Dz,Ayyz) = d1chi(n,Dz,Ayyz) +     t(Ayy)
               d1chi(n,Dy,Ayzz) = d1chi(n,Dy,Ayzz) +     t(Azz)
               d1chi(n,Dz,Ayzz) = d1chi(n,Dz,Ayzz) + 2d0*t(Ayz)
               d1chi(n,Dz,Azzz) = d1chi(n,Dz,Azzz) + 3d0*t(Azz)
            enddo
         else
            call errquit('f fn gen_l error in eval_dang_x_rad',1)
         endif
      else

c     General case

c     Do remainder first

         irem = (iaend-iabeg+1)-((iaend-iabeg+1)/3)*3
         if (irem .eq. 1) then
            if (max_der_component .eq. Dz) then
               do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
c
                  d2chi(n,Dxx,iabeg  ) = d2chi(n,Dxx,iabeg  )
     &                 + d2ang(n,Dxx,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg  ) = d2chi(n,Dxy,iabeg  )
     &                 + d2ang(n,Dxy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dx,n)
                  d2chi(n,Dxz,iabeg  ) = d2chi(n,Dxz,iabeg  )
     &                 + d2ang(n,Dxz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dx,n)
                  d2chi(n,Dyy,iabeg  ) = d2chi(n,Dyy,iabeg  )
     &                 + d2ang(n,Dyy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg  ) = d2chi(n,Dyz,iabeg  )
     &                 + d2ang(n,Dyz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dy,n)
                  d2chi(n,Dzz,iabeg  ) = d2chi(n,Dzz,iabeg  )
     &                 + d2ang(n,Dzz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dz,n)*2d0
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
c
                  d2chi(n,Dxx,iabeg  ) = d2chi(n,Dxx,iabeg  )
     &                 + d2ang(n,Dxx,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg  ) = d2chi(n,Dxy,iabeg  )
     &                 + d2ang(n,Dxy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dx,n)
                  d2chi(n,Dxz,iabeg  ) = d2chi(n,Dxz,iabeg  )
     &                 + d2ang(n,Dxz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dx,n)
                  d2chi(n,Dyy,iabeg  ) = d2chi(n,Dyy,iabeg  )
     &                 + d2ang(n,Dyy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg  ) = d2chi(n,Dyz,iabeg  )
     &                 + d2ang(n,Dyz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dy,n)
                  d2chi(n,Dzz,iabeg  ) = d2chi(n,Dzz,iabeg  )
     &                 + d2ang(n,Dzz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dz,n)*2d0
c
                  d3chi(n,Dxxx,iabeg  ) = d3chi(n,Dxxx,iabeg  )
     &                 + d3ang(n,Dxxx,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iabeg  ) = d3chi(n,Dxxy,iabeg  )
     &                 + d3ang(n,Dxxy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dy,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxx,n)
                  d3chi(n,Dxxz,iabeg  ) = d3chi(n,Dxxz,iabeg  )
     &                 + d3ang(n,Dxxz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxx,n)
                  d3chi(n,Dxyy,iabeg  ) = d3chi(n,Dxyy,iabeg  )
     &                 + d3ang(n,Dxyy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dyy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iabeg  ) = d3chi(n,Dxyz,iabeg  )
     &                 + d3ang(n,Dxyz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dy,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dyz,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxy,n)
                  d3chi(n,Dxzz,iabeg  ) = d3chi(n,Dxzz,iabeg  )
     &                 + d3ang(n,Dxzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iabeg  ) = d3chi(n,Dyyy,iabeg  )
     &                 + d3ang(n,Dyyy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iabeg  ) = d3chi(n,Dyyz,iabeg  )
     &                 + d3ang(n,Dyyz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dyy,n)
                  d3chi(n,Dyzz,iabeg  ) = d3chi(n,Dyzz,iabeg  )
     &                 + d3ang(n,Dyzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iabeg  ) = d3chi(n,Dzzz,iabeg  )
     &                 + d3ang(n,Dzzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dzz,n)*3d0
               enddo
            endif
         elseif (irem .eq. 2) then
            if (max_der_component .eq. Dz) then
               do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
                  d1chi(n,Dx,iabeg+1) = d1chi(n,Dx,iabeg+1)
     &                 + d1ang(n,Dx,iabeg+1)*drad(D0,n)
                  d1chi(n,Dy,iabeg+1) = d1chi(n,Dy,iabeg+1)
     &                 + d1ang(n,Dy,iabeg+1)*drad(D0,n)
                  d1chi(n,Dz,iabeg+1) = d1chi(n,Dz,iabeg+1)
     &                 + d1ang(n,Dz,iabeg+1)*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
                  d1chi(n,Dx,iabeg+1) = d1chi(n,Dx,iabeg+1)
     &                 + d1ang(n,Dx,iabeg+1)*drad(D0,n)
                  d1chi(n,Dy,iabeg+1) = d1chi(n,Dy,iabeg+1)
     &                 + d1ang(n,Dy,iabeg+1)*drad(D0,n)
                  d1chi(n,Dz,iabeg+1) = d1chi(n,Dz,iabeg+1)
     &                 + d1ang(n,Dz,iabeg+1)*drad(D0,n)
c
                  d2chi(n,Dxx,iabeg  ) = d2chi(n,Dxx,iabeg  )
     &                 + d2ang(n,Dxx,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg  ) = d2chi(n,Dxy,iabeg  )
     &                 + d2ang(n,Dxy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dx,n)
                  d2chi(n,Dxz,iabeg  ) = d2chi(n,Dxz,iabeg  )
     &                 + d2ang(n,Dxz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dx,n)
                  d2chi(n,Dyy,iabeg  ) = d2chi(n,Dyy,iabeg  )
     &                 + d2ang(n,Dyy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg  ) = d2chi(n,Dyz,iabeg  )
     &                 + d2ang(n,Dyz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dy,n)
                  d2chi(n,Dzz,iabeg  ) = d2chi(n,Dzz,iabeg  )
     &                 + d2ang(n,Dzz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iabeg+1) = d2chi(n,Dxx,iabeg+1)
     &                 + d2ang(n,Dxx,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg+1) = d2chi(n,Dxy,iabeg+1)
     &                 + d2ang(n,Dxy,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dx,n)
                  d2chi(n,Dxz,iabeg+1) = d2chi(n,Dxz,iabeg+1)
     &                 + d2ang(n,Dxz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dx,n)
                  d2chi(n,Dyy,iabeg+1) = d2chi(n,Dyy,iabeg+1)
     &                 + d2ang(n,Dyy,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg+1) = d2chi(n,Dyz,iabeg+1)
     &                 + d2ang(n,Dyz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dy,n)
                  d2chi(n,Dzz,iabeg+1) = d2chi(n,Dzz,iabeg+1)
     &                 + d2ang(n,Dzz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dz,n)*2d0
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  d1chi(n,Dx,iabeg  ) = d1chi(n,Dx,iabeg  )
     &                 + d1ang(n,Dx,iabeg  )*drad(D0,n)
                  d1chi(n,Dy,iabeg  ) = d1chi(n,Dy,iabeg  )
     &                 + d1ang(n,Dy,iabeg  )*drad(D0,n)
                  d1chi(n,Dz,iabeg  ) = d1chi(n,Dz,iabeg  )
     &                 + d1ang(n,Dz,iabeg  )*drad(D0,n)
                  d1chi(n,Dx,iabeg+1) = d1chi(n,Dx,iabeg+1)
     &                 + d1ang(n,Dx,iabeg+1)*drad(D0,n)
                  d1chi(n,Dy,iabeg+1) = d1chi(n,Dy,iabeg+1)
     &                 + d1ang(n,Dy,iabeg+1)*drad(D0,n)
                  d1chi(n,Dz,iabeg+1) = d1chi(n,Dz,iabeg+1)
     &                 + d1ang(n,Dz,iabeg+1)*drad(D0,n)
c
                  d2chi(n,Dxx,iabeg  ) = d2chi(n,Dxx,iabeg  )
     &                 + d2ang(n,Dxx,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg  ) = d2chi(n,Dxy,iabeg  )
     &                 + d2ang(n,Dxy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dx,n)
                  d2chi(n,Dxz,iabeg  ) = d2chi(n,Dxz,iabeg  )
     &                 + d2ang(n,Dxz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dx,n)
                  d2chi(n,Dyy,iabeg  ) = d2chi(n,Dyy,iabeg  )
     &                 + d2ang(n,Dyy,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg  ) = d2chi(n,Dyz,iabeg  )
     &                 + d2ang(n,Dyz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dy,n)
                  d2chi(n,Dzz,iabeg  ) = d2chi(n,Dzz,iabeg  )
     &                 + d2ang(n,Dzz,iabeg  )*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iabeg+1) = d2chi(n,Dxx,iabeg+1)
     &                 + d2ang(n,Dxx,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iabeg+1) = d2chi(n,Dxy,iabeg+1)
     &                 + d2ang(n,Dxy,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dx,n)
                  d2chi(n,Dxz,iabeg+1) = d2chi(n,Dxz,iabeg+1)
     &                 + d2ang(n,Dxz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dx,n)
                  d2chi(n,Dyy,iabeg+1) = d2chi(n,Dyy,iabeg+1)
     &                 + d2ang(n,Dyy,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iabeg+1) = d2chi(n,Dyz,iabeg+1)
     &                 + d2ang(n,Dyz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dy,n)
                  d2chi(n,Dzz,iabeg+1) = d2chi(n,Dzz,iabeg+1)
     &                 + d2ang(n,Dzz,iabeg+1)*drad(D0,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dz,n)*2d0
c
                  d3chi(n,Dxxx,iabeg  ) = d3chi(n,Dxxx,iabeg  )
     &                 + d3ang(n,Dxxx,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iabeg  ) = d3chi(n,Dxxy,iabeg  )
     &                 + d3ang(n,Dxxy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dy,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxx,n)
                  d3chi(n,Dxxz,iabeg  ) = d3chi(n,Dxxz,iabeg  )
     &                 + d3ang(n,Dxxz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg  )*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxx,n)
                  d3chi(n,Dxyy,iabeg  ) = d3chi(n,Dxyy,iabeg  )
     &                 + d3ang(n,Dxyy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dyy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iabeg  ) = d3chi(n,Dxyz,iabeg  )
     &                 + d3ang(n,Dxyz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dy,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dyz,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dxz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxy,n)
                  d3chi(n,Dxzz,iabeg  ) = d3chi(n,Dxzz,iabeg  )
     &                 + d3ang(n,Dxzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dxz,iabeg  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iabeg  ) = d3chi(n,Dyyy,iabeg  )
     &                 + d3ang(n,Dyyy,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iabeg  ) = d3chi(n,Dyyz,iabeg  )
     &                 + d3ang(n,Dyyz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg  )*drad(Dz,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iabeg  )*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dyy,n)
                  d3chi(n,Dyzz,iabeg  ) = d3chi(n,Dyzz,iabeg  )
     &                 + d3ang(n,Dyzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dyz,iabeg  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg  )*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iabeg  ) = d3chi(n,Dzzz,iabeg  )
     &                 + d3ang(n,Dzzz,iabeg  )*drad(D0,n)
     &                 + d2ang(n,Dzz,iabeg  )*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iabeg  )*drad(Dzz,n)*3d0
                  d3chi(n,Dxxx,iabeg+1) = d3chi(n,Dxxx,iabeg+1)
     &                 + d3ang(n,Dxxx,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg+1)*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iabeg+1) = d3chi(n,Dxxy,iabeg+1)
     &                 + d3ang(n,Dxxy,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg+1)*drad(Dy,n)
     &                 + d2ang(n,Dxy,iabeg+1)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dxx,n)
                  d3chi(n,Dxxz,iabeg+1) = d3chi(n,Dxxz,iabeg+1)
     &                 + d3ang(n,Dxxz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iabeg+1)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg+1)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dxx,n)
                  d3chi(n,Dxyy,iabeg+1) = d3chi(n,Dxyy,iabeg+1)
     &                 + d3ang(n,Dxyy,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg+1)*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iabeg+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dyy,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iabeg+1) = d3chi(n,Dxyz,iabeg+1)
     &                 + d3ang(n,Dxyz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxy,iabeg+1)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iabeg+1)*drad(Dy,n)
     &                 + d2ang(n,Dyz,iabeg+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dyz,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dxz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dxy,n)
                  d3chi(n,Dxzz,iabeg+1) = d3chi(n,Dxzz,iabeg+1)
     &                 + d3ang(n,Dxzz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dxz,iabeg+1)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iabeg+1)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iabeg+1) = d3chi(n,Dyyy,iabeg+1)
     &                 + d3ang(n,Dyyy,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg+1)*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iabeg+1) = d3chi(n,Dyyz,iabeg+1)
     &                 + d3ang(n,Dyyz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dyy,iabeg+1)*drad(Dz,n)
     &                 + d2ang(n,Dyz,iabeg+1)*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dyy,n)
                  d3chi(n,Dyzz,iabeg+1) = d3chi(n,Dyzz,iabeg+1)
     &                 + d3ang(n,Dyzz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dyz,iabeg+1)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iabeg+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iabeg+1)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iabeg+1) = d3chi(n,Dzzz,iabeg+1)
     &                 + d3ang(n,Dzzz,iabeg+1)*drad(D0,n)
     &                 + d2ang(n,Dzz,iabeg+1)*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iabeg+1)*drad(Dzz,n)*3d0
               enddo
            endif
         endif

c     Do the bulk

         do iang = iabeg+irem, iaend, 3
            if (max_der_component .eq. Dz) then
               do n = 1, nq
                  d1chi(n,Dx,iang  ) = d1chi(n,Dx,iang  )
     &                 + d1ang(n,Dx,iang  )*drad(D0,n)
                  d1chi(n,Dy,iang  ) = d1chi(n,Dy,iang  )
     &                 + d1ang(n,Dy,iang  )*drad(D0,n)
                  d1chi(n,Dz,iang  ) = d1chi(n,Dz,iang  )
     &                 + d1ang(n,Dz,iang  )*drad(D0,n)
                  d1chi(n,Dx,iang+1) = d1chi(n,Dx,iang+1)
     &                 + d1ang(n,Dx,iang+1)*drad(D0,n)
                  d1chi(n,Dy,iang+1) = d1chi(n,Dy,iang+1)
     &                 + d1ang(n,Dy,iang+1)*drad(D0,n)
                  d1chi(n,Dz,iang+1) = d1chi(n,Dz,iang+1)
     &                 + d1ang(n,Dz,iang+1)*drad(D0,n)
                  d1chi(n,Dx,iang+2) = d1chi(n,Dx,iang+2)
     &                 + d1ang(n,Dx,iang+2)*drad(D0,n)
                  d1chi(n,Dy,iang+2) = d1chi(n,Dy,iang+2)
     &                 + d1ang(n,Dy,iang+2)*drad(D0,n)
                  d1chi(n,Dz,iang+2) = d1chi(n,Dz,iang+2)
     &                 + d1ang(n,Dz,iang+2)*drad(D0,n)
               enddo
            elseif (max_der_component .eq. Dzz) then
                do n = 1, nq
                  d1chi(n,Dx,iang  ) = d1chi(n,Dx,iang  )
     &                 + d1ang(n,Dx,iang  )*drad(D0,n)
                  d1chi(n,Dy,iang  ) = d1chi(n,Dy,iang  )
     &                 + d1ang(n,Dy,iang  )*drad(D0,n)
                  d1chi(n,Dz,iang  ) = d1chi(n,Dz,iang  )
     &                 + d1ang(n,Dz,iang  )*drad(D0,n)
                  d1chi(n,Dx,iang+1) = d1chi(n,Dx,iang+1)
     &                 + d1ang(n,Dx,iang+1)*drad(D0,n)
                  d1chi(n,Dy,iang+1) = d1chi(n,Dy,iang+1)
     &                 + d1ang(n,Dy,iang+1)*drad(D0,n)
                  d1chi(n,Dz,iang+1) = d1chi(n,Dz,iang+1)
     &                 + d1ang(n,Dz,iang+1)*drad(D0,n)
                  d1chi(n,Dx,iang+2) = d1chi(n,Dx,iang+2)
     &                 + d1ang(n,Dx,iang+2)*drad(D0,n)
                  d1chi(n,Dy,iang+2) = d1chi(n,Dy,iang+2)
     &                 + d1ang(n,Dy,iang+2)*drad(D0,n)
                  d1chi(n,Dz,iang+2) = d1chi(n,Dz,iang+2)
     &                 + d1ang(n,Dz,iang+2)*drad(D0,n)
c
                  d2chi(n,Dxx,iang  ) = d2chi(n,Dxx,iang  )
     &                 + d2ang(n,Dxx,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang  ) = d2chi(n,Dxy,iang  )
     &                 + d2ang(n,Dxy,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dx,n)
                  d2chi(n,Dxz,iang  ) = d2chi(n,Dxz,iang  )
     &                 + d2ang(n,Dxz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dx,n)
                  d2chi(n,Dyy,iang  ) = d2chi(n,Dyy,iang  )
     &                 + d2ang(n,Dyy,iang  )*drad(D0,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang  ) = d2chi(n,Dyz,iang  )
     &                 + d2ang(n,Dyz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dy,n)
                  d2chi(n,Dzz,iang  ) = d2chi(n,Dzz,iang  )
     &                 + d2ang(n,Dzz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iang+1) = d2chi(n,Dxx,iang+1)
     &                 + d2ang(n,Dxx,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang+1) = d2chi(n,Dxy,iang+1)
     &                 + d2ang(n,Dxy,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dx,n)
                  d2chi(n,Dxz,iang+1) = d2chi(n,Dxz,iang+1)
     &                 + d2ang(n,Dxz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dx,n)
                  d2chi(n,Dyy,iang+1) = d2chi(n,Dyy,iang+1)
     &                 + d2ang(n,Dyy,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang+1) = d2chi(n,Dyz,iang+1)
     &                 + d2ang(n,Dyz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dy,n)
                  d2chi(n,Dzz,iang+1) = d2chi(n,Dzz,iang+1)
     &                 + d2ang(n,Dzz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iang+2) = d2chi(n,Dxx,iang+2)
     &                 + d2ang(n,Dxx,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang+2) = d2chi(n,Dxy,iang+2)
     &                 + d2ang(n,Dxy,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dx,n)
                  d2chi(n,Dxz,iang+2) = d2chi(n,Dxz,iang+2)
     &                 + d2ang(n,Dxz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dx,n)
                  d2chi(n,Dyy,iang+2) = d2chi(n,Dyy,iang+2)
     &                 + d2ang(n,Dyy,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang+2) = d2chi(n,Dyz,iang+2)
     &                 + d2ang(n,Dyz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dy,n)
                  d2chi(n,Dzz,iang+2) = d2chi(n,Dzz,iang+2)
     &                 + d2ang(n,Dzz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dz,n)*2d0
               enddo
            elseif (max_der_component .eq. Dzzz) then
               do n = 1, nq
                  d1chi(n,Dx,iang  ) = d1chi(n,Dx,iang  )
     &                 + d1ang(n,Dx,iang  )*drad(D0,n)
                  d1chi(n,Dy,iang  ) = d1chi(n,Dy,iang  )
     &                 + d1ang(n,Dy,iang  )*drad(D0,n)
                  d1chi(n,Dz,iang  ) = d1chi(n,Dz,iang  )
     &                 + d1ang(n,Dz,iang  )*drad(D0,n)
                  d1chi(n,Dx,iang+1) = d1chi(n,Dx,iang+1)
     &                 + d1ang(n,Dx,iang+1)*drad(D0,n)
                  d1chi(n,Dy,iang+1) = d1chi(n,Dy,iang+1)
     &                 + d1ang(n,Dy,iang+1)*drad(D0,n)
                  d1chi(n,Dz,iang+1) = d1chi(n,Dz,iang+1)
     &                 + d1ang(n,Dz,iang+1)*drad(D0,n)
                  d1chi(n,Dx,iang+2) = d1chi(n,Dx,iang+2)
     &                 + d1ang(n,Dx,iang+2)*drad(D0,n)
                  d1chi(n,Dy,iang+2) = d1chi(n,Dy,iang+2)
     &                 + d1ang(n,Dy,iang+2)*drad(D0,n)
                  d1chi(n,Dz,iang+2) = d1chi(n,Dz,iang+2)
     &                 + d1ang(n,Dz,iang+2)*drad(D0,n)
c
                  d2chi(n,Dxx,iang  ) = d2chi(n,Dxx,iang  )
     &                 + d2ang(n,Dxx,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang  ) = d2chi(n,Dxy,iang  )
     &                 + d2ang(n,Dxy,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dx,n)
                  d2chi(n,Dxz,iang  ) = d2chi(n,Dxz,iang  )
     &                 + d2ang(n,Dxz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dx,n)
                  d2chi(n,Dyy,iang  ) = d2chi(n,Dyy,iang  )
     &                 + d2ang(n,Dyy,iang  )*drad(D0,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang  ) = d2chi(n,Dyz,iang  )
     &                 + d2ang(n,Dyz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dy,n)
                  d2chi(n,Dzz,iang  ) = d2chi(n,Dzz,iang  )
     &                 + d2ang(n,Dzz,iang  )*drad(D0,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iang+1) = d2chi(n,Dxx,iang+1)
     &                 + d2ang(n,Dxx,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang+1) = d2chi(n,Dxy,iang+1)
     &                 + d2ang(n,Dxy,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dx,n)
                  d2chi(n,Dxz,iang+1) = d2chi(n,Dxz,iang+1)
     &                 + d2ang(n,Dxz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dx,n)
                  d2chi(n,Dyy,iang+1) = d2chi(n,Dyy,iang+1)
     &                 + d2ang(n,Dyy,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang+1) = d2chi(n,Dyz,iang+1)
     &                 + d2ang(n,Dyz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dy,n)
                  d2chi(n,Dzz,iang+1) = d2chi(n,Dzz,iang+1)
     &                 + d2ang(n,Dzz,iang+1)*drad(D0,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dz,n)*2d0
                  d2chi(n,Dxx,iang+2) = d2chi(n,Dxx,iang+2)
     &                 + d2ang(n,Dxx,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dx,n)*2d0
                  d2chi(n,Dxy,iang+2) = d2chi(n,Dxy,iang+2)
     &                 + d2ang(n,Dxy,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dx,n)
                  d2chi(n,Dxz,iang+2) = d2chi(n,Dxz,iang+2)
     &                 + d2ang(n,Dxz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dx,n)
                  d2chi(n,Dyy,iang+2) = d2chi(n,Dyy,iang+2)
     &                 + d2ang(n,Dyy,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dy,n)*2d0
                  d2chi(n,Dyz,iang+2) = d2chi(n,Dyz,iang+2)
     &                 + d2ang(n,Dyz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dy,n)
                  d2chi(n,Dzz,iang+2) = d2chi(n,Dzz,iang+2)
     &                 + d2ang(n,Dzz,iang+2)*drad(D0,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dz,n)*2d0
c
                  d3chi(n,Dxxx,iang  ) = d3chi(n,Dxxx,iang  )
     &                 + d3ang(n,Dxxx,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iang  )*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iang  )*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iang  ) = d3chi(n,Dxxy,iang  )
     &                 + d3ang(n,Dxxy,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iang  )*drad(Dy,n)
     &                 + d2ang(n,Dxy,iang  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang  )*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iang  )*drad(Dxx,n)
                  d3chi(n,Dxxz,iang  ) = d3chi(n,Dxxz,iang  )
     &                 + d3ang(n,Dxxz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxx,iang  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang  )*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang  )*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iang  )*drad(Dxx,n)
                  d3chi(n,Dxyy,iang  ) = d3chi(n,Dxyy,iang  )
     &                 + d3ang(n,Dxyy,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iang  )*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iang  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dyy,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iang  ) = d3chi(n,Dxyz,iang  )
     &                 + d3ang(n,Dxyz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxy,iang  )*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang  )*drad(Dy,n)
     &                 + d2ang(n,Dyz,iang  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dyz,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dxz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dxy,n)
                  d3chi(n,Dxzz,iang  ) = d3chi(n,Dxzz,iang  )
     &                 + d3ang(n,Dxzz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dxz,iang  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang  )*drad(Dx,n)
     &                 + d1ang(n,Dx,iang  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iang  ) = d3chi(n,Dyyy,iang  )
     &                 + d3ang(n,Dyyy,iang  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iang  )*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iang  )*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iang  ) = d3chi(n,Dyyz,iang  )
     &                 + d3ang(n,Dyyz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dyy,iang  )*drad(Dz,n)
     &                 + d2ang(n,Dyz,iang  )*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iang  )*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iang  )*drad(Dyy,n)
                  d3chi(n,Dyzz,iang  ) = d3chi(n,Dyzz,iang  )
     &                 + d3ang(n,Dyzz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dyz,iang  )*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang  )*drad(Dy,n)
     &                 + d1ang(n,Dy,iang  )*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang  )*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iang  ) = d3chi(n,Dzzz,iang  )
     &                 + d3ang(n,Dzzz,iang  )*drad(D0,n)
     &                 + d2ang(n,Dzz,iang  )*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iang  )*drad(Dzz,n)*3d0
                  d3chi(n,Dxxx,iang+1) = d3chi(n,Dxxx,iang+1)
     &                 + d3ang(n,Dxxx,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+1)*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iang+1)*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iang+1) = d3chi(n,Dxxy,iang+1)
     &                 + d3ang(n,Dxxy,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+1)*drad(Dy,n)
     &                 + d2ang(n,Dxy,iang+1)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang+1)*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iang+1)*drad(Dxx,n)
                  d3chi(n,Dxxz,iang+1) = d3chi(n,Dxxz,iang+1)
     &                 + d3ang(n,Dxxz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+1)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang+1)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang+1)*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iang+1)*drad(Dxx,n)
                  d3chi(n,Dxyy,iang+1) = d3chi(n,Dxyy,iang+1)
     &                 + d3ang(n,Dxyy,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxy,iang+1)*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iang+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dyy,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iang+1) = d3chi(n,Dxyz,iang+1)
     &                 + d3ang(n,Dxyz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxy,iang+1)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang+1)*drad(Dy,n)
     &                 + d2ang(n,Dyz,iang+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dyz,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dxz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dxy,n)
                  d3chi(n,Dxzz,iang+1) = d3chi(n,Dxzz,iang+1)
     &                 + d3ang(n,Dxzz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dxz,iang+1)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang+1)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+1)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iang+1) = d3chi(n,Dyyy,iang+1)
     &                 + d3ang(n,Dyyy,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dyy,iang+1)*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iang+1)*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iang+1) = d3chi(n,Dyyz,iang+1)
     &                 + d3ang(n,Dyyz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dyy,iang+1)*drad(Dz,n)
     &                 + d2ang(n,Dyz,iang+1)*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iang+1)*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iang+1)*drad(Dyy,n)
                  d3chi(n,Dyzz,iang+1) = d3chi(n,Dyzz,iang+1)
     &                 + d3ang(n,Dyzz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dyz,iang+1)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang+1)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+1)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang+1)*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iang+1) = d3chi(n,Dzzz,iang+1)
     &                 + d3ang(n,Dzzz,iang+1)*drad(D0,n)
     &                 + d2ang(n,Dzz,iang+1)*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iang+1)*drad(Dzz,n)*3d0
                  d3chi(n,Dxxx,iang+2) = d3chi(n,Dxxx,iang+2)
     &                 + d3ang(n,Dxxx,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+2)*drad(Dx,n)*3d0
     &                 + d1ang(n,Dx,iang+2)*drad(Dxx,n)*3d0
                  d3chi(n,Dxxy,iang+2) = d3chi(n,Dxxy,iang+2)
     &                 + d3ang(n,Dxxy,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+2)*drad(Dy,n)
     &                 + d2ang(n,Dxy,iang+2)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang+2)*drad(Dxy,n)*2d0
     &                 + d1ang(n,Dy,iang+2)*drad(Dxx,n)
                  d3chi(n,Dxxz,iang+2) = d3chi(n,Dxxz,iang+2)
     &                 + d3ang(n,Dxxz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxx,iang+2)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang+2)*drad(Dx,n)*2d0
     &                 + d1ang(n,Dx,iang+2)*drad(Dxz,n)*2d0
     &                 + d1ang(n,Dz,iang+2)*drad(Dxx,n)
                  d3chi(n,Dxyy,iang+2) = d3chi(n,Dxyy,iang+2)
     &                 + d3ang(n,Dxyy,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxy,iang+2)*drad(Dy,n)*2d0
     &                 + d2ang(n,Dyy,iang+2)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dyy,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dxy,n)*2d0
                  d3chi(n,Dxyz,iang+2) = d3chi(n,Dxyz,iang+2)
     &                 + d3ang(n,Dxyz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxy,iang+2)*drad(Dz,n)
     &                 + d2ang(n,Dxz,iang+2)*drad(Dy,n)
     &                 + d2ang(n,Dyz,iang+2)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dyz,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dxz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dxy,n)
                  d3chi(n,Dxzz,iang+2) = d3chi(n,Dxzz,iang+2)
     &                 + d3ang(n,Dxzz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dxz,iang+2)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang+2)*drad(Dx,n)
     &                 + d1ang(n,Dx,iang+2)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dxz,n)*2d0
                  d3chi(n,Dyyy,iang+2) = d3chi(n,Dyyy,iang+2)
     &                 + d3ang(n,Dyyy,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dyy,iang+2)*drad(Dy,n)*3d0
     &                 + d1ang(n,Dy,iang+2)*drad(Dyy,n)*3d0
                  d3chi(n,Dyyz,iang+2) = d3chi(n,Dyyz,iang+2)
     &                 + d3ang(n,Dyyz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dyy,iang+2)*drad(Dz,n)
     &                 + d2ang(n,Dyz,iang+2)*drad(Dy,n)*2d0
     &                 + d1ang(n,Dy,iang+2)*drad(Dyz,n)*2d0
     &                 + d1ang(n,Dz,iang+2)*drad(Dyy,n)
                  d3chi(n,Dyzz,iang+2) = d3chi(n,Dyzz,iang+2)
     &                 + d3ang(n,Dyzz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dyz,iang+2)*drad(Dz,n)*2d0
     &                 + d2ang(n,Dzz,iang+2)*drad(Dy,n)
     &                 + d1ang(n,Dy,iang+2)*drad(Dzz,n)
     &                 + d1ang(n,Dz,iang+2)*drad(Dyz,n)*2d0
                  d3chi(n,Dzzz,iang+2) = d3chi(n,Dzzz,iang+2)
     &                 + d3ang(n,Dzzz,iang+2)*drad(D0,n)
     &                 + d2ang(n,Dzz,iang+2)*drad(Dz,n)*3d0
     &                 + d1ang(n,Dz,iang+2)*drad(Dzz,n)*3d0
               enddo
            endif
         enddo
      endif
c
      return
      end
