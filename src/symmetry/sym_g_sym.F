      subroutine sym_grad_symmetrize(geom, grad)
C     $Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "geom.fh"
      integer geom
      double precision grad(3,*)
      logical geom_isystype_set
      external geom_isystype_set
c     
c     Apply a projection operator to the gradient
c     so that it posesses the symmetry of the group
c     to machine precision.
c     
      integer iat, jat, k, op
      integer nops, nat, l_tmp, k_tmp
      double precision rnew(3)
      integer save_sys
c     
      nops = sym_number_ops(geom)
      if (nops .eq. 0) return
c     
      if (.not. geom_ncent(geom, nat)) call errquit('s_g_s: geom?',0,
     &       GEOM_ERR)
      if (.not. ma_push_get(mt_dbl, 3*nat, 'gradtmp', l_tmp, k_tmp))
     $     call errquit('sym_grad_sym: scratch space', 3*nat, MA_ERR)
c     
      if (util_print('sym_grad_sym', print_debug)) then
         write(6,*) ' sym_grad_sym: Input gradient'
         call output(grad, 1, 3, 1, nat, 3, nat, 1)
      endif
c
      call dcopy(3*nat, grad, 1, dbl_mb(k_tmp), 1)
c
! sym_apply_cart_op() wraps back into the box, which is wrong for gradients
      if (.not. geom_systype_get(geom, save_sys))
     $     call errquit('sym_grad_symmetrize:systype?',0, GEOM_ERR)
      if (.not.geom_isystype_set(geom,0))
     >     call errquit('sym_grad_symmetrize:failed to set isystype',
     >                  0,GEOM_ERR)
! Call with minu op to tell sym_apply_op() to not apply translatsion
! which is also wrong for gradients
      do iat = 1, nat
         do op = 1, nops
            jat = sym_center_map(geom, iat, op)
            call sym_apply_cart_op(geom,-op, grad(1,iat), rnew)
            do k = 1, 3
               dbl_mb(k_tmp-1 + k + (jat-1)*3) = 
     $              dbl_mb(k_tmp-1 + k + (jat-1)*3) + rnew(k)
            enddo
         end do
      end do
      if (.not.geom_isystype_set(geom,save_sys))
     >     call errquit('sym_grad_symmetrize:failed to set isystype',
     >                  0,GEOM_ERR)
c     
      call dscal(3*nat, 1.0d0/dble(nops+1), dbl_mb(k_tmp), 1)
c
      call dcopy(3*nat, dbl_mb(k_tmp), 1, grad, 1)
c
      if (util_print('sym_grad_sym', print_debug)) then
         write(6,*) ' sym_grad_sym: Output gradient'
         call output(grad, 1, 3, 1, nat, 3, nat, 1)
      endif
c
      if (.not.ma_pop_stack(l_tmp)) call errquit('s_g_s: ma?',0, MA_ERR)
c
      end
      subroutine sym_hess_symmetrize(geom, h)
C     $Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "geom.fh"
      integer geom
      double precision h(*)     ! Assumed (3*nat,3*nat)
c     
c     Apply a projection operator to the hessian
c     so that it posesses the symmetry of the group
c     to machine precision.
c     
      integer iat, jat, op, iatnew, jatnew
      integer nops, nat, l_tmp, k_tmp
      double precision hold(3,3), hnew(3,3)
      integer nat3
c
      nops = sym_number_ops(geom)
      if (nops .eq. 0) return
c     
      if (.not. geom_ncent(geom, nat)) call errquit('s_g_s: geom?',0,
     &       GEOM_ERR)
      nat3 = nat*3
      if (.not. ma_push_get(mt_dbl, (nat3)**2, 'htmp', l_tmp, k_tmp))
     $     call errquit('sym_hess_sym: scratch space', (nat3)**2,
     &       MA_ERR)
c     
      if (util_print('sym_hess_sym', print_debug)) then
         write(6,*) ' sym_hess_sym: Input hessian'
         call output(h, 1, nat3, 1, nat3, nat3, nat3, 1)
      endif
c
      call dcopy(nat3**2, h, 1, dbl_mb(k_tmp), 1)
c
      do iat = 1, nat
         do jat = 1, nat
            call sym_copy_h_block(iat,jat,nat,h,hold)
            do op = 1, nops
               iatnew = sym_center_map(geom, iat, op)
               jatnew = sym_center_map(geom, jat, op)
               call sym_apply_hess_op(geom, op, hold, hnew)
               call sym_add_h_block(iatnew, jatnew, nat, 
     &             dbl_mb(k_tmp), hnew)
            enddo
         end do
      end do
c     
      call dscal(nat3**2, 1.0d0/dble(nops+1), dbl_mb(k_tmp), 1)
c
      call dcopy(nat3**2, dbl_mb(k_tmp), 1, h, 1)
c
      if (util_print('sym_hess_sym', print_debug)) then
         write(6,*) ' sym_hess_sym: Output hessian'
         call output(h, 1, nat3, 1, nat3, nat3, nat3, 1)
      endif
c
      if (.not.ma_pop_stack(l_tmp)) call errquit('s_h_s: ma?',0, MA_ERR)
c
      end
      subroutine sym_der_dip_symmetrize(geom, d)
C     $Id$
      implicit none
#include "errquit.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "geom.fh"
      integer geom
      double precision d(*)     ! Assumed (3,3,nat)
c     
c     Apply a projection operator to the derivative dipole moment
c     so that it posesses the symmetry of the group
c     to machine precision.
c
      integer iat, op, iatnew
      integer nops, nat, l_tmp, k_tmp
      double precision dold(3,3), dnew(3,3)
      integer nat3
c
      nops = sym_number_ops(geom)
c     
      if (.not. geom_ncent(geom, nat)) call errquit('s_d_s: geom?',0,
     &       GEOM_ERR)
      nat3 = nat*3
      if (.not. ma_push_get(mt_dbl, 3*nat3, 'dtmp', l_tmp, k_tmp))
     $     call errquit('sym_dip_sym: scratch space', (nat3)**2, MA_ERR)
c     
      if (util_print('sym_dip_sym', print_debug)) then
         write(6,*) ' sym_dip_sym: Input derivative dipole'
         call output(d, 1, 3, 1, nat3, 3, nat3, 1)
      endif
c
      call dcopy(3*nat3, d, 1, dbl_mb(k_tmp), 1)
c
      do iat = 1, nat
         call sym_copy_h_block(1,iat,1,d,dold)
         do op = 1, nops
            iatnew = sym_center_map(geom, iat, op)
            call sym_apply_hess_op(geom, op, dold, dnew)
            call sym_add_h_block(1,iatnew, 1, 
     &           dbl_mb(k_tmp), dnew)
         enddo
      end do
c     
      call dscal(3*nat3, 1.0d0/dble(nops+1), dbl_mb(k_tmp), 1)
c
      call dcopy(3*nat3, dbl_mb(k_tmp), 1, d, 1)
c
      if (util_print('sym_dip_sym', print_debug)) then
         write(6,*) ' sym_dip_sym: Output deriviative dipole'
         call output(d, 1, 3, 1, nat3, 3, nat3, 1)
      endif
c
      if (.not.ma_pop_stack(l_tmp)) call errquit('s_d_s: ma?',0, MA_ERR)
c
      end
      subroutine sym_apply_hess_op(geom, opnum, h, h1)
      implicit none
c
#include "nwc_const.fh"
#include "geomP.fh"
c
      integer geom, opnum
      double precision h(3,3), h1(3,3)
      integer i, j, k, l
c
c     Apply operators to atom block of Hessian ... R*H*RT
c
c     Hik =  sum(jl) Rij hjl  Rkl
c
c     THIS HAS ONLY BE DESIGNED WITH MOLECULAR SYSTEM IN MIND
c     .... TRANSLATION TERMS NEED ADDING FOR PERIODIC
c
c     stupid 3**4 version for easy debugging
c
      do i = 1, 3
         do k = 1, 3
            h1(i,k) = 0.0d0
            do j = 1, 3
               do l = 1, 3
                  h1(i,k) = h1(i,k) + 
     $                 sym_ops((opnum-1)*3+i,j,geom)*
     $                 sym_ops((opnum-1)*3+k,l,geom)*h(j,l)
               enddo
            enddo
         enddo
      enddo
c
      end
      subroutine sym_copy_h_block(iat,jat,nat,h,hold)
      implicit none
      integer iat, jat, nat
      double precision h(3,nat,3,*)
      double precision hold(3,3)
c
      integer i, j
c
      do i = 1, 3
         do j = 1, 3
            hold(i,j) = h(i,iat,j,jat)
         enddo
      enddo
c
      end
      subroutine sym_add_h_block(iat,jat,nat,h,hnew)
      implicit none
      integer iat, jat, nat
      double precision h(3,nat,3,*)
      double precision hnew(3,3)
c
      integer i, j
c
      do i = 1, 3
         do j = 1, 3
             h(i,iat,j,jat) =  h(i,iat,j,jat) + hnew(i,j)
         enddo
      enddo
c
      end

