      subroutine sym_grad_symmetrize(geom, grad)
C     $Id: sym_g_sym.F,v 1.1 1995-12-13 01:28:13 d3g681 Exp $
      implicit none
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "geom.fh"
      integer geom
      double precision grad(3,*)
c     
c     Apply a projection operator to the gradient
c     so that it posesses the symmetry of the group
c     to machine precision.
c     
      integer iat, jat, k, op
      integer nops, nat, l_tmp, k_tmp
      double precision rnew(3)
c     
      nops = sym_number_ops(geom)
c     
      if (.not. geom_ncent(geom, nat)) call errquit('s_g_s: geom?',0)
      if (.not. ma_push_get(mt_dbl, 3*nat, 'gradtmp', l_tmp, k_tmp))
     $     call errquit('sym_grad_sym: scratch space', 3*nat)
c     
      if (util_print('sym_grad_sym', print_debug)) then
         write(6,*) ' sym_grad_sym: Input gradient'
         call output(grad, 1, 3, 1, nat, 3, nat, 1)
      endif
c
      call dcopy(3*nat, grad, 1, dbl_mb(k_tmp), 1)
c
      do iat = 1, nat
         do op = 1, nops
            jat = sym_center_map(geom, iat, op)
            call sym_apply_cart_op(geom, op, grad(1,iat), rnew)
            do k = 1, 3
               dbl_mb(k_tmp-1 + k + (jat-1)*3) = 
     $              dbl_mb(k_tmp-1 + k + (jat-1)*3) + rnew(k)
            enddo
         end do
      end do
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
      if (.not.ma_pop_stack(l_tmp)) call errquit('s_g_s: ma?',0)
c
      end
