      subroutine sym_movecs_apply_op_1(basis, op, v, t)
C$Id$
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "sym.fh"
#include "bas.fh"
#include "mafdecls.fh"
      integer basis             ! [input]
      integer  op               ! [input]
      double precision v(*)     ! [input]
      double precision t(*)     ! [output]
c
c     Apply the group operation op to the vector of basis function
c     coefficients in v(*), returning the result in t(*)
c
      integer nbf, nbf_per_sh
      integer geom
      integer maxf, maxr, ang_max
      integer l_r, k_r
      integer i
      integer ncent, icent, jcent, ish, jsh, ilo, ihi, jlo, jhi
      integer ishlo, ishhi, jshlo, jshhi
      integer type, ngen, nprim, sphcart
      integer opinv
      logical do_spherical
c
c     Get basis and geom info
c     
      if (.not. bas_geom(basis, geom)) call errquit
     $     ('sym_movecs_apply_op: no geometry from basis?', basis,
     &       BASIS_ERR)
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('sym_movecs_apply_op: bas_numbf?', basis, BASIS_ERR)
      if (.not. bas_high_angular(basis, ang_max))
     $     call errquit('sym_movecs_apply_op: bas_high_angular', basis,
     &       BASIS_ERR)
      if (.not. geom_ncent(geom, ncent)) call errquit
     $     ('sym_bas_irreps: geom_ncent ?', geom, GEOM_ERR)
c     
c     Get the transformation matrices for shells
c
      do_spherical = bas_is_spherical(basis)
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
      if (.not. ma_push_get(MT_DBL, maxr, 'r', l_r, k_r))
     $     call errquit('sym_movecs_apply_op: ma failed', maxr, MA_ERR)
c
      call sym_bas_op
     &    (geom, op, dbl_mb(k_r), maxf, ang_max, do_spherical)
      call sym_inv_op(geom, op, opinv)
c
**      call dfill(nbf, 0.0d0, t, 1) ! Now defined to add in to result
c
      do icent = 1, ncent
         jcent = sym_center_map(geom, icent, opinv)
         if (.not. bas_ce2cnr(basis, icent, ishlo, ishhi))
     $        call errquit('sym_movecs_apply_op: bas_ce2cnr',0,
     &       BASIS_ERR)
         if (.not. bas_ce2cnr(basis, jcent, jshlo, jshhi))
     $        call errquit('sym_movecs_apply_op: bas_ce2cnr',0,
     &       BASIS_ERR)
c
         do ish = ishlo, ishhi
            jsh = jshlo + ish - ishlo ! OP ish -> jsh
c
            if (.not. bas_continfo(basis, ish, type, nprim, ngen,
     $           sphcart)) call errquit
     $           ('sym_movecs_apply_op: bas_continfo', ish,
     &       BASIS_ERR)
            if (.not. bas_cn2bfr(basis, ish, ilo, ihi))
     $           call errquit('sym_movecs_apply_op: bas_cn2bfr', ish,
     &       BASIS_ERR)
            if (.not. bas_cn2bfr(basis, jsh, jlo, jhi))
     $           call errquit('sym_movecs_apply_op: bas_cn2bfr', jsh,
     &       BASIS_ERR)
c
            if (type .eq. -1) then
*               call sym_movecs_app_op2(dbl_mb(k_r), maxf, ang_max,
*     $              0, v(ilo), 1, 1, t(jlo))
               t(jlo) = t(jlo) + v(ilo)
               call sym_movecs_app_op2(dbl_mb(k_r), maxf, ang_max,
     $              1, v(ilo+1), 3, 1, t(jlo+1))
            else if (type .eq. 0) then
               do i = 0, ngen-1
                  t(jlo+i) = t(jlo+i) + v(ilo+i)  ! S functions are easy
               enddo
            else
               nbf_per_sh = (ihi-ilo+1)/ngen
               call sym_movecs_app_op2(dbl_mb(k_r), maxf, ang_max,
     $              type, v(ilo), nbf_per_sh, ngen, t(jlo))
            endif
         enddo
      enddo
c
c      write(6,*) ' sym_movecs_apply_op: ', op
c      do i = 1, nbf
c         write(6,11) v(i), t(i)
c 11      format(2f12.6)
c      enddo
c
      if (.not. ma_pop_stack(l_r)) call errquit
     $     ('sym_movecs_apply_op: ma?', 0, MA_ERR)
c
      end
      subroutine sym_movecs_app_op2_so(r, maxf, ang_max,
     $              type, v, nbf_per_sh, ngen, t)
      implicit none
      integer maxf, ang_max
      double precision r(maxf, maxf, 0:ang_max)
      integer type
      integer nbf_per_sh
      integer ngen
      double precision v(nbf_per_sh, ngen)
      double precision t(nbf_per_sh, ngen)
c
      integer i, j, k
      double precision s
c
c     Apply symmetry operation to MO coefficients of a generally
c     contracted shell
c
      do k = 1, ngen
         do i = 1, nbf_per_sh
            s = v(i,k)
            if (s .ne. 0.0d0) then
               do j = 1, nbf_per_sh
                  t(j,k) = t(j,k) + r(j,i,type)*s
               enddo
            endif
         enddo
      enddo
c
      end
