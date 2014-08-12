      block data block_pottimes
#include "pottimes.fh"
      data otime /.false./
      end
      subroutine potential_from_shell_pair(basis, ish, jsh,
     $     den, ld, tol, d, dinv, lmaxd, npts, xyz, pot, fact,
     $     omask, opresent, ldp)
*
* $Id$
*
      implicit none
#include "errquit.fh"
c
#include "bas.fh"
#include "geom.fh"
#include "canlfit.fh"
#include "pottimes.fh"
#include "util.fh"
#include "mafdecls.fh"
c
#define NOT_DONE(i) int_mb(k_not_done - 1 + i)
c
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
      double precision fact     ! [input]
      integer ldp
      logical omask             ! [input] If true mask primtive pairs
      logical opresent(ldp,*)   ! [input] Optional mask for prim. pairs
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
c     If (omask) primtive pairs are only included if (opresent(i,j))
c
      integer maxprim           ! Max. no. primitives in a contraction
      parameter (maxprim=25)
      integer lmax, lmax2, lencart, lenprod
      parameter (lmax = 5)      ! Max. angular momentum of shell
      parameter (lmax2 = 2*lmax)
      parameter (lencart = ((lmax2+1)*(lmax2+2))/2)
      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax+3))/6)
c
      integer maxl              ! Max order of multipole
      parameter (maxl = 2*lmax2)
      double precision mpoles(-maxl:maxl,0:maxl)
      double precision    tmp(-maxl:maxl,0:maxl)
      double precision     qn(-maxl:maxl,0:maxl)
c
      integer num_not_done, k_not_done, l_not_done
c
      double precision q(-lmax2:lmax2,0:lmax2)
      double precision work(lenprod**2), ndens(lenprod**2),
     $     dens(lencart**2)
      integer k_expi,l_expi,k_coeffi,l_coeffi
      integer k_expj,l_expj,k_coeffj,l_coeffj
      double precision a, b, c, center(3), value
c
      logical odoit
c
      integer itype,inprim,ingen,numi,ibflo,ibfhi,i,iprim,iat
      integer jtype,jnprim,jngen,numj,jbflo,jbfhi,j,jprim,jat
      integer ijtype,nlm, ipt, iipt, n, l, m, geom, ind, luse, leff,
     $     luse2, luse4
      logical iosph, josph
      double precision ri(3), rj(3), x, y, z, alpha, prefac, sqrta, ar,
     $     dmax, radial, factor, sum, ppot, r, rsq, charge, ralpha,
     $     radius2_sq, radius4_sq
      double precision qm2, qm1, q0, qp1, qp2, start
      double precision anls((lmax2+2)*(lmax2+3)/4), a00, a11, 
     $     dd, dd2, dd3, rar, radius, radius_sq, pi
      character*16 tag
      integer nprim_max
      double precision anl_fit, fastj_gaussian_range
      external anl_fit, double_factorial, fastj_gaussian_range
c
      if (npts .le. 0) return
      if (.not. ma_push_get(mt_int, npts, 'pot_shell', 
     $     l_not_done, k_not_done)) call errquit('pot_shell: ma', npts,
     &       MA_ERR)

      if (otime) start = util_cpusec()
c
      pi = 4.0d0*atan(1.0d0)
      if(.not.(bas_nprim_cn_max(basis, nprim_max))) call errquit(
     E     ' pot_shell: nprim_cn_max failed ',basis,BASIS_ERR)
      if (.not. ma_push_get(mt_dbl, nprim_max, 'expi', 
     $     l_expi, k_expi)) call errquit(
     E     'pot_shell: expi alloc failed', nprim_max,    MA_ERR)
      if (.not. ma_push_get(mt_dbl, nprim_max, 'expj', 
     $     l_expj, k_expj)) call errquit(
     E     'pot_shell: expj alloc failed', nprim_max,    MA_ERR)
      if (.not. ma_push_get(mt_dbl, nprim_max, 'coeffi', 
     $     l_coeffi, k_coeffi)) call errquit(
     E     'pot_shell: coeffi alloc failed', nprim_max,    MA_ERR)
      if (.not. ma_push_get(mt_dbl, nprim_max, 'coeffj', 
     $     l_coeffj, k_coeffj)) call errquit(
     E     'pot_shell: coeffj alloc failed', nprim_max,    MA_ERR)

c
c     Get info on the shells
c
      if (.not. bas_geom(basis, geom)) call errquit('fastJ: basis?',0,
     &       BASIS_ERR)
c
      if (.not.bas_continfo(basis, jsh, jtype, jnprim,jngen, 
     $     josph)) call errquit('fastJ: bas_cont_info', jsh,
     &       BASIS_ERR)
      if(jnprim.gt.nprim_max) call errquit(
     A     ' pot_shell: increas maxprim from ',nprim_max,0)
      if (jtype .lt. 0 .or. jngen.ne.1) call errquit
     $     ('fastJ: no general/SP functions yet',0, BASIS_ERR)
      if (jtype .gt. lmax) call errquit
     $     ('fastJ: hard dimension lmax exceeded ', lmax, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, jsh, jbflo, jbfhi))
     $     call errquit('fastJ dens: bas_cn2bfr', jsh, BASIS_ERR)
      if (.not.bas_get_exponent(basis, jsh, dbl_mb(k_expj)))
     $     call errquit('fastJ: bas_get_exponent?',jsh, BASIS_ERR)
      if (.not.bas_get_coeff(basis, jsh, dbl_mb(k_coeffj)))
     $     call errquit('fastJ: bas_get_coeff?',jsh,
     &       BASIS_ERR)
      if (.not. bas_cn2ce(basis, jsh, jat))
     $     call errquit('fastJ: bas_cn2ce?',jsh, BASIS_ERR)
      if (.not. geom_cent_get(geom, jat, tag, rj, charge))
     $     call errquit('fastJ: geom_cart_get ', jat, GEOM_ERR)
      numj = jbfhi - jbflo + 1
c
      if (.not.bas_continfo(basis,ish,itype,inprim,ingen, 
     $     iosph)) call errquit('fastJ: bas_cont_info',ish, BASIS_ERR)
      if(inprim.gt.nprim_max) call errquit(
     A     ' pot_shell: increas maxprim from ',nprim_max,0)
      if (itype .lt. 0 .or. ingen.ne.1) call errquit
     $     ('fastJ: no general/SP functions yet',0, BASIS_ERR)
      if (itype .gt. lmax) call errquit
     $     ('fastJ: hard dimension lmax exceeded ', lmax, BASIS_ERR)
      if (.not. bas_cn2bfr(basis, ish, ibflo, ibfhi))
     $     call errquit('fastJ dens: bas_cn2bfr', ish, BASIS_ERR)
      if (.not.bas_get_exponent(basis, ish, dbl_mb(k_expi)))
     $     call errquit('fastJ: bas_get_exponent?',ish, BASIS_ERR)
      if (.not.bas_get_coeff(basis, ish, dbl_mb(k_coeffi)))
     $     call errquit('fastJ: bas_get_coeff?',ish, BASIS_ERR)
      if (.not. bas_cn2ce(basis, ish, iat))
     $     call errquit('fastJ: bas_cn2ce?',ish, BASIS_ERR)
      if (.not. geom_cent_get(geom, iat, tag, ri, charge))
     $     call errquit('fastJ: geom_cart_get ', iat, GEOM_ERR)
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
            end do
         end do
      else
         call errquit('fastJ: I am confused about spherical',0,
     &       BASIS_ERR)
      end if
c
c     Scale the density by fact and find maximum value
c
      dmax = 0.0d0
      ind = 1
      do j = 1, numj
         do i = 1, numi
            dens(ind) = dens(ind)*fact
            dmax = max(dmax,abs(dens(ind)))
            ind = ind + 1
         end do
      end do
c
      if (otime) then
         time_init = time_init + util_cpusec()-start
         start = util_cpusec()
      endif
c
      call pot_shell_pair_multipoles(
     $     dens,
     $     itype, inprim, ingen, dbl_mb(k_expi), dbl_mb(k_coeffi), 
     R     ri, numi,
     $     jtype, jnprim, jngen, dbl_mb(k_expj), dbl_mb(k_coeffj), 
     R     rj, numj,
     $     dinv, lmaxd,
     $     tol, 
     $     qn, tmp,
     $     mpoles, maxl, center, radius, luse, luse2, luse4,
     $     omask, opresent, ldp)
c
      if (radius .eq. 0.0d0) goto 99999
      radius_sq = radius*radius
      radius2_sq = 4.0d0*radius_sq
      radius4_sq = 16.0d0*radius_sq
c
c     Apply field terms to the potential to simply potential eval
c
      if (luse .ge. 0) then
         call xlm_norm(luse, qn, maxl)
         factor = 4.0d0*pi
         do l = 0, luse
            do m = -l, l
               mpoles(m,l) = mpoles(m,l)*factor*qn(m,l)*qn(m,l) / 
     $              dble(2*l+1)
            enddo
            factor = -factor
         end do
      end if
c
      if (otime) then
         time_make_poles = time_make_poles + util_cpusec() - start
         start = util_cpusec()
      endif
c
c     Now that we have the multipoles make a first pass thru
c     the points doing as many as possible with the MPOLE expansion.
c
      num_not_done = 0
      if (luse .eq. 0) then
         do ipt = 1, npts
            x = center(1) - xyz(1,ipt)
            y = center(2) - xyz(2,ipt)
            z = center(3) - xyz(3,ipt)
            rsq = x*x + y*y + z*z
            if (rsq .gt. radius_sq) then
               pot(ipt) = pot(ipt) - mpoles(0,0)/sqrt(rsq)
            else
               num_not_done = num_not_done + 1
               NOT_DONE(num_not_done) = ipt
            end if
         end do
      else if (luse .eq. 1) then
         do ipt = 1, npts
            x = center(1) - xyz(1,ipt)
            y = center(2) - xyz(2,ipt)
            z = center(3) - xyz(3,ipt)
            rsq = x*x + y*y + z*z
            if (rsq .gt. radius_sq) then
               rsq = 1.0d0/rsq
               r = sqrt(rsq)
               pot(ipt) = pot(ipt) - mpoles(0,0)*r -
     $              (y*mpoles(-1,1)+z*mpoles(0,1)+x*mpoles(1,1)) *
     $              (r*rsq)
            else
               num_not_done = num_not_done + 1
               NOT_DONE(num_not_done) = ipt
            end if
         end do
      else if (luse .eq. 2) then
         do ipt = 1, npts
            x = center(1) - xyz(1,ipt)
            y = center(2) - xyz(2,ipt)
            z = center(3) - xyz(3,ipt)
            rsq = x*x + y*y + z*z
            if (rsq .gt. radius_sq) then
               rsq = 1.0d0/rsq
               r = sqrt(rsq)
               qm2 = 6.0d0*x*y
               qm1 = 3.0d0*y*z
               q0  = z*z - 0.5d0*(x*x+y*y)
               qp1 = 3.0d0*x*z
               qp2 = 3.0d0*(x*x-y*y)
               pot(ipt) = pot(ipt) - 
     $              mpoles(0,0)*r -
     $              (y*mpoles(-1,1)+z*mpoles(0,1)+x*mpoles(1,1)) *
     $              (r*rsq) -
     $              (qm2*mpoles(-2,2)+qm1*mpoles(-1,2)+q0*mpoles(0,2)+
     $              qp1*mpoles(1,2)+qp2*mpoles(2,2)) * (r*rsq*rsq)
            else
               num_not_done = num_not_done + 1
               NOT_DONE(num_not_done) = ipt
            end if
         end do
      else
         do ipt = 1, npts
            x = center(1) - xyz(1,ipt)
            y = center(2) - xyz(2,ipt)
            z = center(3) - xyz(3,ipt)
            rsq = x*x + y*y + z*z
            if (rsq .gt. radius_sq) then
               if (rsq .gt. radius4_sq) then
                  leff = luse4
               else if (rsq .gt. radius2_sq) then
                  leff = luse2
               else
                  leff = luse
               endif
               if (leff .eq. 0) then
                  pot(ipt) = pot(ipt) - mpoles(0,0)/sqrt(rsq)
               else if (leff .eq. 1) then
                  rsq = 1.0d0/rsq
                  r = sqrt(rsq)
                  pot(ipt) = pot(ipt) - mpoles(0,0)*r -
     $                 (y*mpoles(-1,1)+z*mpoles(0,1)+x*mpoles(1,1)) *
     $                 (r*rsq)
               else if (leff .eq. 2) then
                  rsq = 1.0d0/rsq
                  r = sqrt(rsq)
                  qm2 = 6.0d0*x*y
                  qm1 = 3.0d0*y*z
                  q0  = z*z - 0.5d0*(x*x+y*y)
                  qp1 = 3.0d0*x*z
                  qp2 = 3.0d0*(x*x-y*y)
                  pot(ipt) = pot(ipt) - 
     $                 mpoles(0,0)*r -
     $                 (y*mpoles(-1,1)+z*mpoles(0,1)+x*mpoles(1,1)) *
     $                 (r*rsq) -
     $                 (qm2*mpoles(-2,2)+qm1*mpoles(-1,2)+
     $                 q0*mpoles(0,2)+
     $                 qp1*mpoles(1,2)+qp2*mpoles(2,2)) * (r*rsq*rsq)
               else if (leff .ge. 3) then
                  value = 0.0d0
                  rsq = 1.0d0/rsq
                  r = sqrt(rsq)
                  call xlm(leff, x, y, z, qn, maxl)
                  value = 0.0d0
                  do l = 0, leff
                     sum = 0.0d0
                     do m = -l, l
                        sum = sum + mpoles(m,l)*qn(m,l)
                     enddo
                     value = value + r * sum
                     r = r * rsq
                  enddo
                  pot(ipt) = pot(ipt) - value
               end if
            else
               num_not_done = num_not_done + 1
               NOT_DONE(num_not_done) = ipt
            end if
         end do
      endif
c
*      write(6,131)100.0d0*dble(npts-num_not_done)/dble(npts)
* 131  format(' Done ', f8.2,'%')
c

      if (otime) then
         time_poles = time_poles + util_cpusec() - start
         start = util_cpusec()
      endif

      if (num_not_done .eq. 0) goto 99999
c
c     Loop thru primitive distributions
c
      odoit = .true.
      do jprim = 1, jnprim
         do iprim = 1, inprim
            if (omask) odoit = opresent(iprim,jprim)
            
            if (odoit) then
               call gaussian_product(
     $              dbl_mb(k_expi+iprim-1), ri(1), ri(2), ri(3),
     $              dbl_mb(k_expj+jprim-1), rj(1), rj(2), rj(3),
     $              alpha, a, b, c, prefac)
c     
               prefac = prefac*dbl_mb(k_coeffi+iprim-1)*
     C              dbl_mb(k_coeffj+jprim-1)*
     $              4.0d0*pi / alpha
            else
               prefac = 0.0d0
            endif
c
c     Screen
c
            if (odoit .and. (abs(prefac*dmax).gt.tol*1d-2)) then
c
c     Translate density to product center and reform as spherical
c
               if (ijtype .eq. 0) then
                  ndens(1) = dens(1)
               else
                  call cart_dens_trans_prod_sph(
     $                 itype, ri(1), ri(2), ri(3), 
     $                 jtype, rj(1), rj(2), rj(3), 
     $                 a, b, c, work, dinv, lmaxd, 
     $                 dens, ndens)
               end if
c
c     Compute the potential at each point
c
               sqrta = sqrt(alpha)
c
               if (ijtype .eq. 0) then
                  factor = ndens(1)*prefac
                  do iipt = 1, num_not_done
                     ipt = NOT_DONE(iipt)
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     if (ar .lt. anl_xcut(0)) then
                        i = int(anl_scale*ar)
                        dd = ar - anl_table(0,i,1)
                        dd2 = dd*dd
                        dd3 = dd2*dd
                        a00 = anl_table(1,i,1)+dd*anl_table(2,i,1)+
     $                       dd2*anl_table(3,i,1)+dd3*anl_table(4,i,1)
                     else
                        a00 = anl_fac(0,0)/ar
                     end if
                     pot(ipt) = pot(ipt) - factor*a00
                  end do
               else if (ijtype .eq. 1) then
                  do iipt = 1, num_not_done
                     ipt = NOT_DONE(iipt)
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     if (ar .lt. anl_xcut(1)) then
                        i = int(anl_scale*ar)
                        dd = ar - anl_table(0,i,1)
                        dd2 = dd*dd
                        dd3 = dd2*dd
                        a00 = anl_table(1,i,1) + dd*anl_table(2,i,1) +
     $                       dd2*anl_table(3,i,1) + dd3*anl_table(4,i,1)
                        a11 = anl_table(1,i,2) + dd*anl_table(2,i,2) +
     $                       dd2*anl_table(3,i,2) + dd3*anl_table(4,i,2)
                     else
                        rar = 1.0d0/ar
                        a00 = anl_fac(0,0)*rar
                        a11 = anl_fac(1,1)*rar*rar*rar
                     end if
                     pot(ipt) = pot(ipt) - prefac*(
     $                    a00*ndens(1) + 
     $                    a11*(y*ndens(2)+z*ndens(3)+x*ndens(4)))
                  end do
               else if (ijtype .eq. 2) then
                  ralpha = 1.0d0 / alpha
                  do iipt = 1, num_not_done
                     ipt = NOT_DONE(iipt)
                     x = xyz(1,ipt)-a
                     y = xyz(2,ipt)-b
                     z = xyz(3,ipt)-c
                     r = sqrt(x*x + y*y + z*z)
                     ar = sqrta * r
                     qm2 = 6.0d0*x*y
                     qm1 = 3.0d0*y*z
                     q0  = z*z - 0.5d0*(x*x+y*y)
                     qp1 = 3.0d0*x*z
                     qp2 = 3.0d0*(x*x-y*y)
                     call anl_fit_group_2(ijtype, ar, anls)
                     pot(ipt) = pot(ipt) - prefac*(
     $                    anls(1)*ndens(1) + 
     $                    anls(2)*(y*ndens(2)+z*ndens(3)+x*ndens(4)) +
     $                    anls(3)*(qm2*ndens(5)+qm1*ndens(6)+
     $                    q0*ndens(7)+qp1*ndens(8) + qp2*ndens(9)) +
     $                    anls(4)*ndens(10)*ralpha)
                  end do
               else 
                  ralpha = 1.0d0 / alpha
                  do iipt = 1, num_not_done
                     ipt = NOT_DONE(iipt)
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
                           end do
                           ppot = ppot + factor*sum
                           radial = radial * ralpha
                        end do
                     end do
                     pot(ipt) = pot(ipt) - ppot
                  end do
               end if
            end if              ! screen and odoit
         end do                 ! iprim
      end do                    ! jprim
c
      if (otime) time_direct = time_direct + util_cpusec() - start
c
99999 if (.not. ma_chop_stack(l_not_done)) call errquit
     $     ('pot_shell: ma corrupt?',0, MA_ERR)
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
               end do
            end do
         end do
      end if
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
               end do
            end do
         end do
      end if
c
      end
      subroutine pot_shell_pair_multipoles(
     $     dens,
     $     itype, inprim, ingen, expi, coeffi, ri, numi,
     $     jtype, jnprim, jngen, expj, coeffj, rj, numj,
     $     dinv, lmaxd,
     $     tol, 
     $     qn, tmp,
     $     mpoles, maxl, center, radius, luse, luse2, luse4,
     $     omask, opresent, ldp)
      implicit none
#include "errquit.fh"
c
c     Given a block of the density matrix over a pair of CARTESIAN
c     shells and a detailed description of those shells, return
c
c     .   mpoles() = multipole moments about center WITHOUT field terms
c     .   center() = natural center of the shell pair
c     .   radius   = radius beyond which mpole expansion is valid
c     .   luse/2/4 = order of mpole required to obtain desired precision
c     .              at distances of 1*, 2* and 4*radius.
c     .              These may be -1, indicating potential<tol
c
c     If radius is returned as 0.0d0 then the is no significant overlap
c     between the shells and the entire pair may be neglected.
c     
c     If spherical basis functions are used, the density matrix block
c     must first be transformed using dens_cart_to_sph().
c
c     Input shell parameters (numi/j = dimension of cartesian shells)
c
      integer itype, inprim, ingen, numi
      integer jtype, jnprim, jngen, numj
      double precision expi(inprim), coeffi(inprim), ri(3)
      double precision expj(jnprim), coeffj(jnprim), rj(3)
c
      double precision dens(numi,*) ! Cartesian density block
      integer lmaxd
      double precision dinv(*)  ! [input] dinv from xlm_coeff_inv
      double precision tol      ! [input] required accuracy
      integer maxl
      double precision qn(-maxl:maxl,0:maxl),tmp(-maxl:maxl,0:maxl) ! [scratch]
      double precision mpoles(-maxl:maxl,0:maxl) ! [output] mpoles about center
      double precision center(3) ! [output] natural center
      double precision radius   ! [output] Radius beyond which mpole valid
      integer luse, luse2, luse4 ! [output] Order of mpole for r=radius*{1,2,4}
      integer ldp
      logical omask             ! [input]
      logical opresent(ldp,*)   ! [input]
c
c     Local variables
c
      integer lmax, lmax2, lenprod,dfmax
      parameter (lmax = 5)      ! Max. angular momentum of shell
      parameter (lmax2 = 2*lmax)
      parameter (lenprod = ((lmax2+1)*(lmax2+2)*(lmax2+3))/6)
      parameter(dfmax=2*lmax2+1)
      double precision df(dfmax)
      double precision work(lenprod**2), ndens(lenprod**2) 
c
      logical odoit
      integer ijtype
      integer i, j, iprim, jprim
      integer l, m, n, nlm, leff
      double precision expi_min, expj_min, a, b, c, alpha, prefac, dmax
      double precision factor, tolp, range, pi, rootpi, scale, x, y, z
      double precision distance, q0, radial
c
c     Externals
c
      double precision double_factorial, fastj_gaussian_range
      external double_factorial, fastj_gaussian_range
c
c     Sanity checks
c
      if (itype.gt.lmax .or. jtype.gt.lmax) call errquit
     $     ('pot_shell_pair_multipole: lmax?',itype*1000+jtype,
     &       BASIS_ERR)
      if (ingen.ne.1 .or. jngen.ne.1) call errquit
     $     ('p_s_p_m: no generally conracted functions ', 0,
     &       BASIS_ERR)
c     
      ijtype = itype + jtype
      if((2*ijtype + 1).gt.dfmax) call errquit(
     .     ' potshellmultip: dfmax must be ',2*ijtype + 1,
     &       BASIS_ERR)
      pi = 4.0d0*atan(1.0d0)
      rootpi = sqrt(pi)
      do i = 1, 2*ijtype + 1
         df(i) = double_factorial(i)
      enddo
c
      dmax = 0.0d0
      do j = 1, numj
         do i = 1, numi
            dmax = max(dmax,abs(dens(i,j)))
         end do
      end do
c
*      write(6,*) ' INPUT DENSITY ', numi, numj
*      call output(dens, 1, numi, 1, numj, numi, numj, 1)
c
c     Compute the multipoles, center and radius of convergence
c
c     First determine the center as being the GPT center of 
c     the most diffuse exponents. 
c
      expj_min = 1d300
      expi_min = 1d300
      if (omask) then
         do jprim = 1, jnprim
            do iprim = 1, inprim
               if (opresent(iprim,jprim)) then
                  expj_min = min(expj_min, expj(jprim))
                  expi_min = min(expi_min, expi(iprim))
               endif
            enddo
         enddo
      else
         do jprim = 1, jnprim
            expj_min = min(expj_min, expj(jprim))
         end do
         do iprim = 1, inprim
            expi_min = min(expi_min, expi(iprim))
         end do
      endif
      call gaussian_product(
     $     expi_min, ri(1), ri(2), ri(3),
     $     expj_min, rj(1), rj(2), rj(3),
     $     alpha, center(1), center(2), center(3), prefac)
c
c     Don't compute too high an order mpole.  Check for precision later.
c
c     THIS NEEDS TO BE MADE TOL DEPENDENT.
c
      leff = maxl ! The old stuff is not accurate enuf min(maxl, max(4,2*(ijtype+1)))
c
      call xlm_norm(leff, qn, maxl)
      do l = 0, leff
         do m = -l,l
            qn(m,l) = 1.0d0/(qn(m,l)*qn(m,l))
         enddo
      enddo
c
      call dfill((2*maxl+1)*(leff+1), 0.0d0, mpoles, 1)
      radius = 0.0d0
c
c     Next, form mpoles for primitives and update the radius.
c
      odoit = .true.
      do jprim = 1, jnprim
         do iprim = 1, inprim
            if (omask) odoit = opresent(iprim,jprim)
            if (odoit) then
               call gaussian_product(
     $              expi(iprim), ri(1), ri(2), ri(3),
     $              expj(jprim), rj(1), rj(2), rj(3),
     $              alpha, a, b, c, prefac)
c     
               prefac = prefac*coeffi(iprim)*coeffj(jprim)
            else
               prefac = 0.0d0
            endif
c
c     Screen with the same test as used below.
c
            if (abs(prefac*dmax*4d0*pi/alpha) .gt. tol*1d-2 
     $           .and. odoit) then 
               prefac = prefac*rootpi/(4.0d0*alpha*sqrt(alpha))
c
c     Translate density to product center and reform as spherical
c
               if (ijtype .eq. 0) then
                  ndens(1) = dens(1,1)
               else
                  call cart_dens_trans_prod_sph(
     $                 itype, ri(1), ri(2), ri(3), 
     $                 jtype, rj(1), rj(2), rj(3), 
     $                 a, b, c, work, dinv, lmaxd, 
     $                 dens, ndens)
               end if
c
c     Form multipoles at the GPT center (a,b,c) for this shell
c
               do l = 0, leff
                  do m = -l, l
                     tmp(m,l) = 0.0d0
                  end do
               end do
               nlm = 0
               scale = prefac
               do n = 0, ijtype
                  radial = scale
                  do l = n, 0, -2
                     factor = radial*df(n+l+1)
                     do m = -l, l
                        nlm = nlm + 1
                        tmp(m,l) = tmp(m,l)+ndens(nlm)*factor
                     end do
                     radial = radial * (2.0d0*alpha)
                  end do
                  scale = scale / (2.0d0*alpha)
               end do
               do l = 0, ijtype
                  do m = -l, l
                     tmp(m,l) = tmp(m,l)*qn(m,l)
                  enddo
               enddo
c
c     Translate to the desired center and add on to the total
c
               x = center(1)-a
               y = center(2)-b
               z = center(3)-c
               distance = sqrt(x*x + y*y + z*z)
c
               if (distance .gt. 1d-12) then 
                  if (ijtype .eq. 0) then
                     q0 = tmp(0,0)
                     call xlm(leff, -x, -y, -z, tmp, maxl)
                     do l = 0, leff
                        do m = -l, l
                           mpoles(m,l) = mpoles(m,l) + tmp(m,l)*q0
                        end do
                     end do
                  else
                     call xlm_translate2(leff, x, y, z, tmp, maxl)
                     do l = 0, leff
                        do m = -l, l
                           mpoles(m,l) = mpoles(m,l) + tmp(m,l)
                        end do
                     end do
                  endif
               else
                  do l = 0, ijtype
                     do m = -l, l
                        mpoles(m,l) = mpoles(m,l) + tmp(m,l)
                     end do
                  end do
               end if
c
c     Finally, adjust the radius of convergence to include this function
c
               tolp =  min(tol*1d-4/abs(dmax*prefac),1d-6)
               range = fastj_gaussian_range(ijtype, alpha, tolp) 
               range = range + distance
               radius = max(radius,range)
c
            end if
         end do
      end do
c
      if (radius .eq. 0.0d0) then
         luse = -1
         luse2 = -1
         luse4 = -1
      else
c
c     Adjust the maximum angular momentum so that we barely obtain 
c     the desired precision beyond radius, radius*2 and radius*4
c
         call xlm_thresh(leff, mpoles, maxl, tol, radius, luse)
         call xlm_thresh(luse, mpoles, maxl, tol, radius*2d0, luse2)
         call xlm_thresh(luse, mpoles, maxl, tol, radius*4d0, luse4)
      endif
c
      end
