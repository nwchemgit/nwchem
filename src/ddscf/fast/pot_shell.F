      subroutine potential_from_shell_pair(basis, ish, jsh,
     $     den, ld, tol, d, dinv, lmaxd, npts, xyz, pot, fact)
      implicit none
#include "bas.fh"
#include "geom.fh"
      integer basis, ish, jsh   ! [input] Basis and shell info
      integer ld
      double precision den(ld,*) ! [input] Density matrix
      integer npts              ! [input] No. of points
      double precision tol      ! [input] Used to screen primitive pairs
      integer lmaxd             ! [input]
      double precision d(((lmaxd+1)*(lmaxd+2))/2,-lmaxd:lmaxd,0:lmaxd), 
     $     dinv(*)              ! [input] xlm_coeff_inv(lmaxd,d,dinv)
      double precision xyz(3,npts) ! [input] Coordinates
      double precision pot(npts) ! [output] Potential at points
      double precision fact
c
c     Compute the potential due to the charge distribution
c
c     P(r) = Xi(r) Xj(r) Dij 
c
c     where Xi(r) are basis functions in shell ish, ditto for j,
c     and Dij is the density matrix.
c
c     Products of primitives are neglected if 
c     .  (GPT radial prefactor) * density * (contraction coeffs) < tol
c
c     The efficiency of the routine increases a lot for large npts.
c
c     THE RESULT IS JUST ADDED INTO POT WHICH MUST THEREFORE BE 
c     ZEROED BY THE CALLING ROUTINE IF DESIRED.
c
      integer maxprim           ! Max. no. primitives in a contraction
      parameter (maxprim=25)
      integer lmax, lmax2, lencart, lenprod, maxn
      parameter (lmax = 3)      ! Max. angular momentum of shell
      parameter (lmax2 = 2*lmax)
      parameter (lencart = ((lmax2+1)*(lmax2+2))/2)
      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax2+3))/6)
      parameter (maxn = 1024)   ! Max. grid dimension
c
      double precision q(-lmax2:lmax2,0:lmax2)
      double precision work(lenprod**2), ndens(lenprod**2),
     $     dens(lencart**2)
      double precision expi(maxprim), coeffi(maxprim)
      double precision expj(maxprim), coeffj(maxprim)
      double precision a, b, c
c
      integer itype,inprim,ingen,numi,ibflo,ibfhi,i,iprim,iat
      integer jtype,jnprim,jngen,numj,jbflo,jbfhi,j,jprim,jat
      integer ijtype,nlm, ipt, n, l, m, geom, ind
      logical iosph, josph
      double precision ri(3), rj(3), x, y, z, alpha, prefac, sqrta, ar,
     $     dmax, radial, factor, sum, ppot, r, charge, ralpha
      double precision anls((lmax2+2)*(lmax2+3)/4)
      character*16 tag
      double precision anl_fit
      external anl_fit
c
c     Get info on the shells
c
      if (.not. bas_geom(basis, geom)) call errquit('fastJ: basis?',0)
c
      if (.not.bas_continfo(basis, jsh, jtype, jnprim,jngen, 
     $     josph)) call errquit('fastJ: bas_cont_info', jsh)
      if (jtype .lt. 0 .or. jngen.ne.1) call errquit
     $     ('fastJ: no general/SP functions yet',0)
      if (jtype .gt. lmax) call errquit
     $     ('fastJ: hard dimension lmax exceeded ', lmax)
      if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $     call errquit('fastJ dens: bas_cn2bfr', jsh)
      if (.not.bas_get_exponent(basis, jsh, expj))
     $     call errquit('fastJ: bas_get_exponent?',jsh)
      if (.not.bas_get_coeff(basis, jsh, coeffj))
     $     call errquit('fastJ: bas_get_coeff?',jsh)
      if (.not. bas_cn2ce(basis, jsh, jat))
     $     call errquit('fastJ: bas_cn2ce?',jsh)
      if (.not. geom_cent_get(geom, jat, tag, rj, charge))
     $     call errquit('fastJ: geom_cart_get ', jat)
      numj = jbfhi - jbflo + 1
c
      if (.not.bas_continfo(basis,ish,itype,inprim,ingen, 
     $     iosph)) call errquit('fastJ: bas_cont_info',ish)
      if (itype .lt. 0 .or. ingen.ne.1) call errquit
     $     ('fastJ: no general/SP functions yet',0)
      if (itype .gt. lmax) call errquit
     $     ('fastJ: hard dimension lmax exceeded ', lmax)
      if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $     call errquit('fastJ dens: bas_cn2bfr', ish)
      if (.not.bas_get_exponent(basis, ish, expi))
     $     call errquit('fastJ: bas_get_exponent?',ish)
      if (.not.bas_get_coeff(basis, ish, coeffi))
     $     call errquit('fastJ: bas_get_coeff?',ish)
      if (.not. bas_cn2ce(basis, ish, iat))
     $     call errquit('fastJ: bas_cn2ce?',ish)
      if (.not. geom_cent_get(geom, iat, tag, ri, charge))
     $     call errquit('fastJ: geom_cart_get ', iat)
      numi = ibfhi - ibflo + 1
c
      ijtype = itype + jtype
c
c     Copy the density into a local array so that all elements
c     are contiguous, transforming from spherical as necessary.
c
      if (iosph .and. josph) then
         numi = ((itype+1)*(itype+2))/2
         numj = ((jtype+1)*(jtype+2))/2
         call dens_sph_to_cart(itype, jtype, den, ld, dens, 
     $        work, numi, numj)
      else if (.not. (iosph .and. josph)) then
         ind = 1
         do j = 1, numj
            do i = 1, numi
               dens(ind) = den(i,j)
               ind = ind + 1
            enddo
         enddo
      else
         call errquit('fastJ: I am confused about spherical',0)
      endif
c
      dmax = 0.0d0
      ind = 1
      do j = 1, numj
         do i = 1, numi
            dmax = max(dmax,abs(dens(ind)))
            ind = ind + 1
         enddo
      enddo
c
c     Loop thru primitive distributions
c
      do jprim = 1, jnprim
         do iprim = 1, inprim
            call gaussian_product(
     $           expi(iprim), ri(1), ri(2), ri(3),
     $           expj(jprim), rj(1), rj(2), rj(3),
     $           alpha, a, b, c, prefac)
c
            prefac = prefac*coeffi(iprim)*coeffj(jprim)*
     $           4.0d0*3.1415926535897932d0 / alpha
c
c     Screen
c
            if (abs(prefac*dmax) .gt. tol) then
c
      prefac=fact*prefac
c
c     Translate density to product center and reform as spherical
c
               if (ijtype .eq. 0) then
                  ndens(1) = dens(1)
               else
                  call cart_dens_trans_prod_sph(
     $                 itype, ri(1), ri(2), ri(3), 
     $                 jtype, rj(1), rj(2), rj(3), 
     $                 a, b, c, work, dinv, lmax2, 
     $                 dens, ndens)
               endif
c
c     Compute the potential at each point
c
               sqrta = sqrt(alpha)
c
               if (ijtype .eq. 0) then
                  factor = ndens(1)*prefac
                  do ipt = 1, npts
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     pot(ipt) = pot(ipt) - factor*anl_fit(0,0,ar)
                  enddo
               else if (ijtype .eq. 1) then
                  do ipt = 1, npts
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     call anl_fit_group(1, ar, anls)
                     pot(ipt) = pot(ipt) - 
     $                    prefac*ndens(1)*anls(1)
                     pot(ipt) = pot(ipt) - prefac*anls(2)*
     $                    (y*ndens(2) + z*ndens(3) + x*ndens(4))
                  enddo
               else
                  ralpha = 1.0d0 / alpha
                  do ipt = 1, npts
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     call xlm(ijtype, x, y, z, q, lmax2)
                     call anl_fit_group(ijtype, ar, anls)
c     
                     nlm = 0
                     ppot = 0.0d0
                     ind = 1
                     do n = 0, ijtype
                        radial = prefac
                        do l = n, 0, -2
                           factor = radial * anls(ind) ! anl_fit(n,l,ar)
                           ind = ind + 1
                           sum = 0.0d0
                           do m = -l, l
                              nlm = nlm + 1
                              sum = sum + q(m,l)*ndens(nlm)
                           enddo
                           ppot = ppot + factor*sum
                           radial = radial * ralpha
                        enddo
                     enddo
                     pot(ipt) = pot(ipt) - ppot
                  enddo
               endif
            endif
         enddo
      enddo
c
      end
      subroutine dens_sph_to_cart(itype, jtype, dens, ld, denc,
     $     work, numi, numj)
      implicit none
#include "mafdecls.fh"
#include "spcartP.fh"
      integer itype, jtype, ld, numi, numj
      double precision dens(ld,*)
      double precision work(numi,*)
      double precision denc(numi,numj)
c
c     Given a density block over a pair of spherical shells 
c     transform it to cartesians
c
c     numi and numj should be input as the dimension of the cartesian 
c     shells.
c
c     The integral routines MUST be initialized so that the
c     common block in spcartP is setup.
c
      integer i, j, numj_sph, icart, jcart
c
      integer iic,iis,iil
      double precision Dtrans
      Dtrans(iic,iis,iil) =
     &    dbl_mb((int_mb(k_sp2c_lindx+iil))+
     &           ((iis+iil)*(iil+1)*(iil+2)/2)
     &           + iic - 1)
c
c     First transform i
c
      numj_sph = jtype+jtype+1
c
      if (itype .le. 1) then
         call dcopy(numi*numj_sph, dens, 1, work, 1)
      else
         call dfill(numi*numj_sph, 0.0d0, work, 1)
         do j = 1, numj_sph
            do i = 1, itype+itype+1
               do icart = 1, numi
                  work(icart,j) = work(icart,j) + 
     $                 Dtrans(icart,i-itype-1,itype)*dens(i,j)
               enddo
            enddo
         enddo
      endif
c
      if (jtype .le. 1) then
         call dcopy(numi*numj, work, 1, denc, 1)
      else
         call dfill(numi*numj, 0.0d0, denc, 1)
         do j = 1, jtype+jtype+1
            do jcart = 1, numj
               do icart = 1, numi
                  denc(icart,jcart) = denc(icart,jcart) + 
     $                 work(icart,j)*Dtrans(jcart,j-jtype-1,jtype)
               enddo
            enddo
         enddo
      endif
c
      end

