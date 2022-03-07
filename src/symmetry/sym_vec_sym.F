      subroutine sym_vec_symmetrize(geom, basis, vec)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
      integer geom, basis
      double precision vec(*)   ! [input/output]
c
c     Symmetrize a skeleton vector in the given basis set.
c
c     V <- (1/h) * sum(R) V * R
c
c     where h = the order of the group and R = operators of the
c     group (including the identity)
c
c     Vec is assumed to be a Hamiltonian-like vector NOT a 
c     vector of coefficients (must use RT instead).
c
c     This operation executes in the local processor with no
c     communicaiton or syncrhonization.
c
      integer l_tmp, k_tmp, l_r, k_r
      integer ncent, ncont
      integer maxf, maxr, ang_max, nbf
      integer l_ce_cnr, l_ce_bfr, l_cn_info, l_cn_bfr
      integer k_ce_cnr, k_ce_bfr, k_cn_info, k_cn_bfr
      logical status
c
c     Get all basis set information for efficiency (for small
c     basis C60 test case bas routine dominate timings, but for
c     all other cases bas times are small).
c
      if (.not. geom_ncent(geom, ncent)) call errquit
     $     ('sym_symmetrize: geom_ncent ?', geom, GEOM_ERR)
      if (.not. bas_numcont(basis, ncont)) call errquit
     $     ('sym_symmetrize: basis handle?', 0, BASIS_ERR)
      status = .true.
      status = status .and. ma_push_get(MT_INT, 2*ncent, 'cecnr',
     $     l_ce_cnr, k_ce_cnr)
      status = status .and. ma_push_get(MT_INT, 2*ncent, 'cebfr',
     $     l_ce_bfr, k_ce_bfr)
      status = status .and. ma_push_get(MT_INT, 4*ncont, 'cninf',
     $     l_cn_info, k_cn_info)
      status = status .and. ma_push_get(MT_INT, 2*ncont, 'cnbfr',
     $     l_cn_bfr, k_cn_bfr)
      call bas_all_ce_info(basis, int_mb(k_ce_cnr), int_mb(k_ce_bfr))
      call bas_all_cn_info(basis, int_mb(k_cn_info),int_mb(k_cn_bfr))
c
      if (.not. bas_numbf(basis, nbf)) call errquit
     $     ('sym_vec_symmetrize: bas_numbf failed',0, BASIS_ERR)
c
      if (.not. bas_high_angular(basis, ang_max))
     $     call errquit('sym_symmetrize: bas_high_angular', basis,
     &       BASIS_ERR)
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
c
      if (.not. ma_push_get(MT_DBL, nbf, 'tmp', l_tmp, k_tmp))
     $     call errquit('sym_vec_symmetrize: ma failed', nbf, MA_ERR)
      if (.not. ma_push_get(MT_DBL, maxr, 'r', l_r, k_r))
     $     call errquit('sym_vec_symmetrize: ma failed', maxr, MA_ERR)
c
      call sym_vec_symm2(geom, basis, vec, dbl_mb(k_tmp), 
     $     dbl_mb(k_r), maxf, ang_max,
     $     int_mb(k_ce_cnr), int_mb(k_ce_bfr),
     $     int_mb(k_cn_info), int_mb(k_cn_bfr), nbf, ncent)
c
      status = status .and. ma_pop_stack(l_r)
      status = status .and. ma_pop_stack(l_tmp)
      status = status .and. ma_pop_stack(l_cn_bfr)
      status = status .and. ma_pop_stack(l_cn_info)
      status = status .and. ma_pop_stack(l_ce_bfr)
      status = status .and. ma_pop_stack(l_ce_cnr)
      if (.not. status) call errquit('sym_symmetrize: ma pop?', 0,
     &       MA_ERR)
c
      end
      subroutine sym_vec_symm2(geom, basis, vec, 
     $     tmp, r, maxf, ang_max, ce_cnr,
     $     ce_bfr, cn_info, cn_bfr, nbf, ncent)
      implicit none
#include "errquit.fh"
#include "bas.fh"
      integer geom, basis
      integer icent, icent_new, maxf, ang_max, nbf, ncent
      double precision vec(nbf), tmp(nbf)
      double precision r(maxf, maxf, 0:ang_max)
      integer ce_cnr(2,*), ce_bfr(2,*), cn_info(4,*), cn_bfr(2,*)
c
      integer i, j, ilo, ihi, ni, ioff, joff
      integer icn, icnlo, icnhi
      integer icebflo, icebfhi, jcebflo
      integer ish, ish_type
      integer ingen, itype, inprim
      integer isphcart, i_nbf_per_sh
      integer op
      logical do_spherical
c
      double precision sum, scale
      integer sym_center_map, sym_number_ops
      external sym_center_map, sym_number_ops
c
c     Identity
c
      do i = 1, nbf
         tmp(i) = vec(i)
      end do
      do_spherical = bas_is_spherical(basis)
c
c     Other operators
c
      do op = 1, sym_number_ops(geom)
         call sym_bas_op(geom, op, r, maxf, ang_max, do_spherical)
c
c     Centers
c
         do icent = 1, ncent
            icent_new = sym_center_map(geom, icent, op)
c
c     Figure out range of basis functions on this atom
c
            icebflo = ce_bfr(1,icent)
            jcebflo = ce_bfr(1,icent_new)
c
            icebfhi = ce_bfr(2,icent)
            icnlo   = ce_cnr(1,icent)
            icnhi   = ce_cnr(2,icent)
c
c     loop thru general contractions on the center
c     
            do icn = icnlo, icnhi
               itype    = cn_info(1,icn)
               inprim   = cn_info(2,icn)
               ingen    = cn_info(3,icn)
               isphcart = cn_info(4,icn)
               ilo      = cn_bfr(1,icn)
               ihi      = cn_bfr(2,icn)
               if (itype .ge. 0) then
                  i_nbf_per_sh = (ihi - ilo + 1) / ingen
               else if (itype .eq. -1) then
                  i_nbf_per_sh = 1
               else
                  call errquit('sym_sym2: invalid shell type', 
     .                 itype, 0)
               end if
c
c     Loop thru shells in each general contraction and construct
c     indexes and types of a single shell block
c
               do ish = 1, ingen
                  ioff = ilo + i_nbf_per_sh*(ish-1) - 1
                  joff = ioff + jcebflo - icebflo
                  if (itype .ge. 0) then
                     ish_type = itype
                  else if (ish .eq. 1) then
                     ish_type = 0 ! S component of sp shell
                  else 
                     ish_type = 1 ! P component of sp shell
                  end if
c
c     Finally have all of the necessary indices ... apply the
c     operator ... number of functions here assumes that we
c     have cartesian functions
c
                  if (do_spherical) then
                    ni = 2*ish_type+1
                  else
                    ni = ((ish_type+1)*(ish_type+2))/2
                  endif
c
c     tmp(j) <- sum(i) vec(i) R(i,j) for j on icent_new, i on icent
c
                  do j = 1, ni
                     sum = 0.0d0
                     do i = 1, ni
                        sum = sum + vec(i+ioff)*R(i,j,ish_type)
                     end do
                     tmp(j+joff) = tmp(j+joff) + sum
                  end do
               end do
            end do 
         end do
      end do
c
      scale = 1.0d0 / (sym_number_ops(geom)+1)
      do i = 1, nbf
         vec(i) = tmp(i) * scale
      end do
c
      end
