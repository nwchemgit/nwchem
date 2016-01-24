      subroutine dimqm_screening(rtdb, geom, basis, dimxyz, screen)
c
c   Subroutine to calculate the electric field due to the QM electrons
c   
c   Called from dimqm_main.F

c      use constants
      implicit none
#include "nwc_const.fh"
#include "errquit.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "bas_ibs_dec.fh"
#include "int_nbf.fh"
#include "stdio.fh"
#include "apiP.fh"
#include "util.fh"
#include "bas.fh"
#include "dimqm.fh"

c    Input Variables
      integer rtdb
      integer geom
      integer basis
      double precision dimxyz(3, nDIM)
      double precision screen(nDIM)
c
      integer nshell, nbf_max
      integer ishell, ilo, ihi, jatom
      integer ibas, ucont, itype, inp, igen, iexp, icf, iatom, igeom
      integer lens
      double precision S(isz_1e)
      double precision scr(isz_1e)
      double precision overlap(nDIM)
      integer i_nbf_x, i_nbf_s, j_nbf_x, j_nbf_s
#include "bas_exndcf_sfn.fh"
#include "bas_ibs_sfn.fh"

      if (.not. bas_geom(basis, geom)) call errquit
     $   ('screening: bad basis', 555, BASIS_ERR)
      if (.not. bas_numcont(basis, nshell)) call errquit
     $   ('screening: bas_numcont failed for basis', basis, BASIS_ERR)
      if (.not. bas_nbf_cn_max(basis,nbf_max)) call errquit
     &   ('screening: bas_nbf_cn_max failed',555, BASIS_ERR)
      lens = sqrt(real(isz_1e))
      S = 0.0d0
c      write(luout,*) "len?", isz_1e
      overlap = 0.0d0
      do ishell = 1, nshell
        if (.not. bas_cn2bfr(basis, ishell, ilo, ihi)) call errquit
     &      ('screening: bas_cn2bfr failed for i',basis,BASIS_ERR)
        call int_nogencont_check(basis,'screening:basis')
        ibas = basis + basis_handle_offset
        ucont = (sf_ibs_cn2ucn(ishell,ibas))
        itype = infbs_cont(CONT_TYPE ,ucont,ibas)
        inp   = infbs_cont(CONT_NPRIM,ucont,ibas)
        igen  = infbs_cont(CONT_NGEN ,ucont,ibas)
        iexp  = infbs_cont(CONT_IEXP ,ucont,ibas)
        icf   = infbs_cont(CONT_ICFP ,ucont,ibas)
        iatom = (sf_ibs_cn2ce(ishell,ibas))
        igeom = ibs_geom(ibas)
        do jatom = 1, nDIM
c          write(luout,*) "DIM: ", jatom
c          write(luout,*) "Shell: ", ishell
          call dimqm_screen2(coords(1,iatom,igeom),
     $                     dbl_mb(mb_exndcf(iexp,ibas)),
     $                     dbl_mb(mb_exndcf(icf, ibas)),
     $                     inp, igen, itype, dimxyz(:, jatom),
     $                     S, lens)
c          i_nbf_x = int_nbf_x(itype)
c          i_nbf_s = int_nbf_s(itype)
c          j_nbf_x = int_nbf_x(0)
c          j_nbf_s = j_nbf_x
c          call spcart_tran1e(S, scr,
c     $           j_nbf_x, i_nbf_x, 0, 1,
c     $           j_nbf_s, i_nbf_s, itype, igen, .false.)
c          write(luout,*) "Transform"
c          write(luout,*) S(1:lens)
          overlap(jatom) = overlap(jatom) + sum(abs(S))
        end do

      end do
c
      do jatom = 1, nDIM
          screen(jatom) = erfc(overlap(jatom))
      end do
c      write(luout,*) "Total Overlap:"
c      write(luout,*) overlap
c      write(luout,*) "Screen factors:"
c      write(luout,*) erfc(overlap)
      end subroutine dimqm_screening

      subroutine dimqm_screen2(xyzi, expi, coefi, i_nprim, i_ngen, Li,
     $                         dimxyz, s, lens)
      implicit none
#include "nwc_const.fh"
#include "hnd_rys.fh"
#include "stdio.fh"
#include "hnd_tol.fh"
#include "dimqm.fh"
#include "dimqm_constants.fh"

      double precision xyzi(3)
      double precision expi(i_nprim)
      double precision coefi(i_nprim, i_ngen)
      integer i_nprim
      integer i_ngen
      integer Li
      double precision dimxyz(3)
      double precision s(lens)
      integer lens
c
      double precision xi, yi, zi
      integer lit, maxi
      double precision xj, yj, zj
      integer ljt, maxj, ljtmod
      double precision rr
      double precision tol
c
      integer ig
      double precision ai, arri, axi, ayi, azi, csi, csj
      double precision aj, aa, aa1, dum, ax, ay, az, fac
c
      double precision dij, t
      integer i, j, ij
      double precision xint, yint, zint
      double precision xs(Li+1, 3), ys(Li+1, 3), zs(Li+1, 3) ! Lj = 0 since all DIM atoms are spherical gaussians
      integer Nxyz(3), ix, iy, iz, jx, jy, jz

      tol = 2.30258d+00*itol
c
c     ---- i shell ----
c
      xi = xyzi(1)
      yi = xyzi(2)
      zi = xyzi(3)
      lit = Li + 1
      maxi = lit*(lit+1)/2
c      write(luout,*) "ngen:", i_ngen
c
c     ---- DIM atom ----
c
      xj = dimxyz(1)
      yj = dimxyz(2)
      zj = dimxyz(3)
      ljt = 1 ! All DIM atoms are spherical (Lj = 0, so Lj+1 = 1)
      ljtmod = ljt + 2
      maxj = ljt*(ljt+1)/2  ! Always 1
      rr = (xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2
      nroots = (lit+ljt - 2)/2 + 1
      s = zero
c
c    --- i primitive ---
c
      do 7000 ig = 1, i_nprim
        ai = expi(ig)
        arri = ai*rr
        axi  = ai*xi
        ayi  = ai*yi
        azi  = ai*zi
        csi  = coefi(ig, i_ngen)
c        write(luout,*) csi
c
c    --- DIM primitive ---
c
        ! No loop cause DIM atoms only have 1 gaussian
        aj = 1.0d0  ! Width of DIM gaussian
        aa = ai+aj
        aa1 = one/aa
        dum = aj*arri*aa1
        if(dum .gt. tol) then
c            write(luout,*) "value", dum
c            write(luout,*) "greater than tolerance", tol
c            write(luout,*) "Overlap is zero"
            go to 7000
        end if
        fac = exp(-dum)
        csj = one ! Contraction coefficent to DIM gaussian
        ax = (axi + aj*xj)*aa1
        ay = (ayi + aj*yj)*aa1
        az = (azi + aj*zj)*aa1
c
c    ---- density factor ----
c
        dij = fac * csi * csj
c        write(luout,*) "dij:", dij
c
c     --- overlap ---
c
        t = sqrt(aa1)
        do j = 1, ljtmod
          do i = 1, lit
            call dimqm_overlap(ax, ay, az, xi, yi, zi,
     $                         xj, yj, zj, i, j, t,
     $                         xint, yint, zint)
            xs(i,j) = xint*t
            ys(i,j) = yint*t
            zs(i,j) = zint*t
c            write(luout,*) xint, yint, zint
          end do ! i = 1, lit
        end do ! j = 1, ljtmod
c
        ij = 0
        do i = 1, maxi
          call getNxyz(Li, i, Nxyz)
          ix = Nxyz(1) + 1
          iy = Nxyz(2) + 1
          iz = Nxyz(3) + 1
          do j = 1, maxj
            call getNxyz(0, j, Nxyz)
            jx = Nxyz(1) + 1 ! Should always be 1
            jy = Nxyz(2) + 1 ! Should always be 1
            jz = Nxyz(3) + 1 ! Should always be 1
            ij = ij + 1
            dum = xs(ix,jx)*ys(iy,jy)*zs(iz,jz)
            s(ij) = s(ij) + dum*dij
          end do ! j = 1, maxj
        end do ! i = 1, maxi
7000  end do ! ig = 1, i_nprim
c      write(luout,*) s

      end subroutine dimqm_screen2

      subroutine dimqm_overlap(x0, y0, z0, xi, yi, zi, xj, yj, zj,
     $                         ni, nj, t, xint, yint, zint)
      implicit none
#include "hnd_whermt.fh"
#include "dimqm_constants.fh"
#include "stdio.fh"
      double precision x0, y0, z0
      double precision xi, yi, zi
      double precision xj, yj, zj
      integer ni, nj 
      double precision t, xint, yint, zint
c
      integer i, imin, imax, ii, jj, npts
      double precision dum
      double precision px, py, pz
      double precision ax, ay, az
      double precision bx, by, bz
      double precision ptx, pty, ptz
      
      xint = zero
      yint = zero
      zint = zero
      npts = (ni + nj - 2)/2 + 1
c      write(luout,*) "npts:", npts
      imin = hermin(npts)
      imax = hermax(npts)
c      write(luout,*) "min, max", imin, imax
      do i = imin, imax
        dum = w(i)
c        write(luout,*) "w", dum
        px = dum
        py = dum
        pz = dum
        dum = h(i) * t
        ptx = dum+x0
        pty = dum+y0
        ptz = dum+z0
        ax =  ptx-xi
        ay =  pty-yi
        az =  ptz-zi
        bx =  ptx-xj
        by =  pty-yj
        bz =  ptz-zj
        do ii = 1, ni-1
          px = px*ax
          py = py*ay
          pz = pz*az
        end do
        do jj = 1, nj-1
          px = px*bx
          py = py*by
          pz = pz*bz
        end do
        xint = xint + px
        yint = yint + py
        zint = zint + pz
      end do
      return

      end subroutine dimqm_overlap
