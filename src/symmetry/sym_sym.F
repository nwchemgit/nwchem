      subroutine sym_symmetrize(geom, basis, odensity, g_a)
C$Id$
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
      integer geom, basis, g_a
      logical odensity
c
c     Symmetrize a skeleton matrix (in a global array) in the
c     given basis set.
c
c     B = (1/2h) * sum(R) [RT * (A + AT) * R]
c
c     where h = the order of the group and R = operators of the
c     group (including the identity)
c
c     Note that density matrices tranform according to slightly
c     different rules to Hamiltonian matrices if components
c     of a shell (e.g., cartesian d's) are not orthonormal.
c     (see Dupuis and King, IJQC 11, 613-625, 1977)
c
      integer g_b
      integer icent, jcent, idim, jdim
      integer l_aij, k_aij, l_aji, k_aji, l_tmp, k_tmp, l_r, k_r
      integer max_at_bf, maxd,  ncent, ilo, ihi, ncont
      integer maxf, maxr, ang_max
      integer l_ce_cnr, l_ce_bfr, l_cn_info, l_cn_bfr
      integer k_ce_cnr, k_ce_bfr, k_cn_info, k_cn_bfr
      logical status
c
      integer next, ijcent, nproc
      integer nxtask
      external nxtask
c
*: ang_max used below
      if (.not. bas_high_angular(basis, ang_max))
     $     call errquit('sym_symmetrize: bas_high_angular', basis,
     &       BASIS_ERR)
*rak:      status = bas_is_spherical(basis)
*rak:      if (status.and.(ang_max.gt.1)) then
*rak:        if (ga_nodeid().eq.0) then
*rak:          write(6,*)' spherical basis set and symmetry ',
*rak:     &        'are not allowed yet'
*rak:          write(6,*)
*rak:     &        ' you can use no symmetry and spherical basis sets '
*rak:          write(6,*)'or'
*rak:          write(6,*)' cartesian basis sets with symmetry'
*rak:        endif
*rak:        call errquit('sym_symmetrize: fatal error',911)
*rak:      endif
*
      call ga_sync()
      if (.not. ga_duplicate(g_a, g_b, 'sym_sym b')) call errquit
     $     ('sym_symmetrize: failed to allocate temp GA', 0, GA_ERR)
      call ga_zero(g_b)
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
c     Figure out size of largest atom block and highest angular
c     momentum and then allocate memory
c
      max_at_bf = 0
      do icent = 1, ncent
         if (.not. bas_ce2bfr(basis, icent, ilo, ihi))
     $        call errquit('rhf_fock_2e: bas_ce2bfr failed', icent,
     &       BASIS_ERR)
         max_at_bf = max(max_at_bf, ihi-ilo+1)
      end do
      maxd = max_at_bf**2       ! Size of largest atom block
c
      maxf = ((ang_max+1)*(ang_max+2))/2 ! Size of largest cart. shell
      maxr = (ang_max+1)*maxf*maxf ! Size of transformation matrices
c
      status = .true.
      status = status .and. ma_push_get(MT_DBL, maxd, 'aij',
     $     l_aij, k_aij)
      status = status .and. ma_push_get(MT_DBL, maxd, 'aji',
     $     l_aji, k_aji)
      status = status .and. ma_push_get(MT_DBL, maxd, 'tmp',
     $     l_tmp, k_tmp)
      status = status .and. ma_push_get(MT_DBL, maxr, 'r',
     $     l_r, k_r)
      if (.not. status) call errquit
     $     ('sym_symmetrize: ma failed', 2*maxd + maxr, MA_ERR)
c
      nproc  = ga_nnodes()
      ijcent = 0
      next   = nxtask(nproc, 1)
      do icent = 1, ncent
         do jcent = 1, ncent
            if (next .eq. ijcent) then
               call get_atom_block(g_a, basis, icent, jcent,
     $              dbl_mb(k_aij), idim, jdim)
               call get_atom_block(g_a, basis, jcent, icent,
     $              dbl_mb(k_aji), jdim, idim)
               
               call sym_symm2(geom, basis, g_b, icent, jcent,
     $              dbl_mb(k_aij), dbl_mb(k_aji), dbl_mb(k_tmp),
     $              idim, jdim, dbl_mb(k_r), maxf, ang_max,
     $              odensity, int_mb(k_ce_cnr), int_mb(k_ce_bfr),
     $              int_mb(k_cn_info), int_mb(k_cn_bfr))
c
               next = nxtask(nproc, 1)
            end if
            ijcent = ijcent + 1
         end do
      end do
c
      next = nxtask(-nproc, 1)
c
      status = status .and. ma_pop_stack(l_r)
      status = status .and. ma_pop_stack(l_tmp)
      status = status .and. ma_pop_stack(l_aji)
      status = status .and. ma_pop_stack(l_aij)
      status = status .and. ma_pop_stack(l_cn_bfr)
      status = status .and. ma_pop_stack(l_cn_info)
      status = status .and. ma_pop_stack(l_ce_bfr)
      status = status .and. ma_pop_stack(l_ce_cnr)
      if (.not. status) call errquit('sym_symmetrize: ma pop?', 0,
     &       MA_ERR)
c
      call ga_copy(g_b, g_a)
      if (.not. ga_destroy(g_b)) call errquit
     $     ('sym_symmetrize: ga_destroy of b', 0, GA_ERR)
      call ga_sync()
c
c*$*optimize(1)
c*$*scalar optimize(1)
c
      end
      subroutine sym_symm2(geom, basis, g_b, icent, jcent, aij, aji,
     $     tmp, idim, jdim, r, maxf, ang_max, odensity, ce_cnr,
     $     ce_bfr, cn_info, cn_bfr)
      implicit none
#include "errquit.fh"
#include "bas.fh"
      integer geom, basis, g_b
      integer icent, jcent, idim, jdim, maxf, ang_max
      double precision aij(idim,jdim), aji(jdim,idim), tmp(idim, jdim),
     $     r(maxf, maxf, 0:ang_max)
      logical odensity
      integer ce_cnr(2,*), ce_bfr(2,*), cn_info(4,*), cn_bfr(2,*)
c
c     Arrive in here with atomic blocks Aij and Aji ... fold them
c     together and loop thru the operations in the group to symmetrize
c
      integer i, j, k, l, ilo, ihi, jlo, jhi, ni, nj, ioff, joff
      integer icn, jcn, icnlo, icnhi, jcnlo, jcnhi
      integer icebflo, icebfhi, jcebflo, jcebfhi
      integer ish, jsh, ish_type, jsh_type
      integer ingen, jngen, itype, jtype, inprim, jnprim
      integer isphcart, jsphcart, i_nbf_per_sh, j_nbf_per_sh
      integer op, icent_new, jcent_new, opinv
      logical do_spherical
c
      double precision sum(28,28), scale, factor, test
      integer sym_center_map, sym_number_ops
      external sym_center_map, sym_number_ops
c
      scale = 1.0d0 / (2.0d0 * (sym_number_ops(geom)+1))
c
      test = 0.0d0
      do i = 1, idim
         do j = 1, jdim
            aij(i,j) = scale*(aij(i,j) + aji(j,i))
            test = test + aij(i,j)*aij(i,j)
         end do
      end do
      if (test .eq. 0.0d0) return
c
      do_spherical = bas_is_spherical(basis)
c
c     Figure out range of basis functions on this atom
c
      icebflo = ce_bfr(1,icent)
      icebfhi = ce_bfr(2,icent)
      jcebflo = ce_bfr(1,jcent)
      jcebfhi = ce_bfr(2,jcent)
      icnlo   = ce_cnr(1,icent)
      icnhi   = ce_cnr(2,icent)
      jcnlo   = ce_cnr(1,jcent)
      jcnhi   = ce_cnr(2,jcent)
c     
c     handle the identity
c
      call upd_atom_block(g_b, basis, icent, jcent, aij)
c
c     loop thru other operators
c
      do op = 1, sym_number_ops(geom)
         call dfill(idim*jdim, 0.0d0, tmp, 1)
c
         if (.not. odensity) then
            icent_new = sym_center_map(geom, icent, op)
            jcent_new = sym_center_map(geom, jcent, op)
         else
            call sym_inv_op(geom, op, opinv)
            icent_new = sym_center_map(geom, icent, opinv)
            jcent_new = sym_center_map(geom, jcent, opinv)
         end if
c
         call sym_bas_op(geom, op, r, maxf, ang_max, do_spherical)
c     
c     loop thru general contractions on each center
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
               i_nbf_per_sh = 0   ! to get rid of compiler warning
               call errquit('sym_sym2: invalid shell type', itype,
     &       INPUT_ERR)
            end if
            do jcn = jcnlo, jcnhi
               jtype    = cn_info(1,jcn)
               jnprim   = cn_info(2,jcn)
               jngen    = cn_info(3,jcn)
               jsphcart = cn_info(4,jcn)
               jlo      = cn_bfr(1,jcn)
               jhi      = cn_bfr(2,jcn)
               if (jtype .ge. 0) then
                  j_nbf_per_sh = (jhi - jlo + 1) / jngen
               else if (jtype .eq. -1) then
                  j_nbf_per_sh = 1
               else
                  j_nbf_per_sh = 0  ! to get rid of compiler warning
                  call errquit('sym_sym2: invalid shell type', jtype,
     &       INPUT_ERR)
               end if
c
               if (i_nbf_per_sh.gt.28 .or. j_nbf_per_sh.gt.28)
     $              call errquit('sym_symmetrize: sum not big enuf',0,
     &       INPUT_ERR)
c
c     Loop thru shells in each general contraction and construct
c     indexes and types of a single shell block
c
               do ish = 1, ingen
                  ioff = (ilo - icebflo) + i_nbf_per_sh*(ish-1)
                  if (itype .ge. 0) then
                     ish_type = itype
                  else if (ish .eq. 1) then
                     ish_type = 0 ! S component of sp shell
                  else 
                     ish_type = 1 ! P component of sp shell
                  end if
                  do jsh = 1, jngen
                     joff = (jlo - jcebflo) + j_nbf_per_sh*(jsh-1)
                     if (jtype .ge. 0) then
                        jsh_type = jtype
                     else if (jsh .eq. 1) then
                        jsh_type = 0 ! S component of sp shell
                     else 
                        jsh_type = 1 ! P component of sp shell
                     end if
c
c     Finally have all of the necessary indices ... apply the
c     operator ... number of functions here assumes that we
c     have cartesian functions
c
                     if (do_spherical) then
                       ni = 2*ish_type+1
                       nj = 2*jsh_type+1
                     else
                       ni = ((ish_type+1)*(ish_type+2))/2
                       nj = ((jsh_type+1)*(jsh_type+2))/2
                     endif
c
c     tmp(i,j) <- sum(k,l) R(k,i)*A(k,l)*R(l,j)
c
c
c     In the following it is important to use the sparsity of
c     the operators
c
                     if (odensity) then
c
c     Density matrix ... R*D*RT
c
                        do j = 1, nj
                           do k = 1, ni
                              sum(k,j)=0.0d0
                           end do
                           do l = 1, nj
                              factor = R(j,l,jsh_type)
                              if (factor .ne. 0.0d0) then
                                 do k = 1, ni
                                    sum(k,j) = sum(k,j) + 
     $                                   aij(ioff+k,joff+l)*factor
                                 end do
                              end if
                           end do
                        end do
                        do k = 1, ni
                           do i = 1, ni
                              factor = R(i,k,ish_type)
                              if (factor .ne. 0.0d0) then
                                 do j = 1, nj
                                    tmp(ioff+i,joff+j) =
     $                                   tmp(ioff+i,joff+j) +
     $                                   factor*sum(k,j)
                                 end do
                              end if
                           end do
                        end do
                     else
c
c     Hamiltonian matrix ... RT*H*R
c
                        do j = 1, nj
                           do k = 1, ni
                              sum(k,j)=0.0d0
                           end do
                           do l = 1, nj
                              factor = R(l,j,jsh_type)
                              if (factor .ne. 0.0d0) then
                                 do k = 1, ni
                                    sum(k,j) = sum(k,j) + 
     $                                   aij(ioff+k,joff+l)*factor
                                 end do
                              end if
                           end do
                        end do
                        do k = 1, ni
                           do i = 1, ni
                              factor = R(k,i,ish_type)
                              if (factor .ne. 0.0d0) then
                                 do j = 1, nj
                                    tmp(ioff+i,joff+j) =
     $                                   tmp(ioff+i,joff+j) +
     $                                   factor*sum(k,j)
                                 end do
                              end if
                           end do
                        end do
                     end if
                  end do         ! jsh  
               end do            ! ish  
            end do               ! jcn
         end do                  ! icn
         call upd_atom_block(g_b, basis, icent_new, jcent_new, tmp)
      end do                     ! op
c
      end
