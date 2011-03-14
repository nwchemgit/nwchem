      subroutine sym_geom_project(geom, tol)
C     $Id$
      implicit none
#include "errquit.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "sym.fh"
#include "mafdecls.fh"
#include "util.fh"
      integer geom
      double precision tol
c     
c     Apply a projection operator to the geometry
c     so that it posesses the symmetry of the group
c     to machine precision.  An atom and the image of
c     that atom under the operations of the group are
c     considered to be identical iff they are less than
c     tol distant from each other.  
c     
      integer iat, jat, k, op
      integer nops, nat, l_tmp, k_tmp
      double precision rnew(3), r, dx, dy, dz
      logical geom_check_handle, geom_print
      external geom_check_handle, geom_print
c     
      if (.not. geom_check_handle(geom, 'sym_geom_project'))
     $     call errquit('sym_geom_project: bad geom',0, GEOM_ERR)
c
      nops = sym_number_ops(geom)
c     
      nat = ncenter(geom)
      if (.not. ma_push_get(mt_dbl, 3*nat, 'coordtmp', l_tmp, k_tmp))
     $     call errquit('sym_geom_project: scratch space', 3*nat,
     &       MA_ERR)
c     
      if (util_print('sym_geom_project', print_debug)) then
         write(6,*) ' sym_geom_project: Input geometry'
         if (.not. geom_print(geom)) call errquit
     $        ('sym_geom_project: geom print failed',0, GEOM_ERR)
      endif
c
c     Map periodic system back to main unit cell with identity symmetry
c
      if(isystype(geom).gt.0) then
         do iat = 1, nat
            op = -9999 ! Secret identity flag - maps back to unit cell only
            call sym_apply_cart_op(geom, op, coords(1,iat,geom), rnew)
            coords(1,iat,geom)=rnew(1)
            coords(2,iat,geom)=rnew(2)
            coords(3,iat,geom)=rnew(3)
         end do
      end if
c
      call dcopy(3*nat, coords(1,1,geom), 1, dbl_mb(k_tmp), 1)
c
      do iat = 1, nat
         do op = 1, nops
            jat = sym_center_map(geom, iat, op)
C           This call maps back to unit cell in periodic cells
            call sym_apply_cart_op(geom, op, coords(1,iat,geom), rnew)
            dx = rnew(1)-coords(1,jat,geom)
            dy = rnew(2)-coords(2,jat,geom)
            dz = rnew(3)-coords(3,jat,geom)
            r = sqrt(dx**2 + dy**2 + dz**2)
            if (r.le.tol) then
               do k = 1, 3
                  dbl_mb(k_tmp-1 + k + (jat-1)*3) = 
     $                 dbl_mb(k_tmp-1 + k + (jat-1)*3) + rnew(k)
               enddo
            else
               if (.not. geom_print(geom)) call errquit
     $            ('sym_geom_project: geom print failed',0, GEOM_ERR)
               call errquit('sym_geom_project: sym_center_map'//
     $             ' is inconsistent with requested accuracy',
     $              op, INPUT_ERR)
            end if
         end do
      end do
c     
      call dscal(3*nat, 1.0d0/dble(nops+1), dbl_mb(k_tmp), 1)
c
      call dcopy(3*nat, dbl_mb(k_tmp), 1, coords(1,1,geom), 1)
c     
c     Should also munge any applied electric field which must be
c     totally symmetric and therefore invariant to all group
c     operations.
c     
c     
      if (util_print('sym_geom_project', print_debug)) then
         write(6,*) ' sym_geom_project: Output geometry'
         if (.not. geom_print(geom)) call errquit
     $        ('sym_geom_project: geom print failed',0, GEOM_ERR)
      endif
c
      if (.not.ma_pop_stack(l_tmp)) call errquit('s_g_p: ma?',0, MA_ERR)
c
c     May have changed the coords ... recompute derived info
c
        call geom_compute_values(geom)

      end

