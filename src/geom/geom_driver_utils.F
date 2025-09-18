#ifdef NEW_DRIVER_UTILS
      subroutine geom_update_cart_from_int_new(geom, ds, alpha, err)
#else      
      subroutine geom_update_cart_from_int(geom, ds, alpha, err)
#endif      
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "util_params.fh"
      integer geom
      double precision alpha, ds(*)
      double precision err      ! [output] Returns the error
c
c     Adjust the cartesian coordinates in geom to reflect a
c     displacement of the internal coordinates alpha*ds().  
c     ds() is in units of bohr and radians.
c
c     Also enforce symmetry.
c
      integer max_cent, max_nzvar
      parameter (max_cent = nw_max_atom)
      parameter (max_nzvar= 3*max_cent)
      double precision  p(max_nzvar) ! Current internals
      double precision  q(max_nzvar) ! Target internals
      double precision dq(max_nzvar) ! Step in internals
      double precision  x(max_nzvar) ! Cartesian coordinates
      double precision dx(max_nzvar) ! Step in Cartesian coordinates
      double precision xsave(max_nzvar) ! Save first order step for recovery
c
      integer l_bi, k_bi, i, ipass,nat, ncart, nzvar
      double precision err1, xmax, xmax1
      logical odebug, oprint
      logical geom_zmt_get_nzvar, geom_compute_zmatrix
      double precision bohr, deg
      external geom_impose_constraints_on_q
c
      odebug = (ga_nodeid().eq.0).and.util_print('optitoc',print_never)
      oprint = (ga_nodeid().eq.0).and.util_print('xupdate',print_high)
      oprint = oprint .or. odebug
c
      if (.not. geom_zmt_get_nzvar(geom, nzvar))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      if (.not. geom_ncent(geom, nat))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      if (.not. geom_compute_zmatrix(geom, q))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      ncart = nat * 3
c
      do i = 1, nzvar
         dq(i) = alpha*ds(i)
      enddo
c
c     Form the target internals in bohr & degrees
c
      bohr = cau2ang
      deg  = bohr*180d0/(4d0*atan(1d0))
      if (.not. geom_compute_zmatrix(geom, q))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      call geom_zmat_ico_scale(geom, dq, bohr, deg)
*      if (.not. geom_print_zmatrix(geom, dq, 'Step  ', .true.))
*     $     call errquit('fjl',0)
      do i = 1, nzvar
         q(i) = q(i) + dq(i)
      enddo
      call geom_zmat_ico_scale(geom, dq, 1d0/bohr, 1d0/deg)
      err = 0.0d0
      do i = 1, nzvar
         err = max(err, abs(dq(i)))
      enddo
      if (odebug) then
         write(LuOut,*) ' Target internals '
         call output(q,1,nzvar,1,1,nzvar,1,1)
         call util_flush(LuOut)
      endif
c
      if (.not. ma_push_get(mt_dbl,ncart*nzvar,'mem bi',l_bi,k_bi))
     $     call errquit('opt_int_to_cart: ma', ncart*nzvar, MA_ERR)
c
c     In a non-redundant set of coordinates the iteration will be 
c     quadratically convergent once very large steps have been damped 
c     out ... well this would be true except some torsions and
c     angles can be very strongly coupled.
c
#ifdef NEW_DRIVER_UTILS      
      call geom_bandbi(geom)
      call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
      if (.not. geom_cart_coords_get(geom, x))
     $     call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)

#endif      
      do ipass = 1, 5
c
c     dq -> dx -> x -> new q
c
#ifndef NEW_DRIVER_UTILS         
         call geom_bandbi(geom)
         call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
#endif         
         call dgemv('n', ncart, nzvar, 1d0, dbl_mb(k_bi), ncart,
     $        dq, 1, 0.0d0, dx, 1)
         call sym_grad_symmetrize(geom, dx)
#ifndef NEW_DRIVER_UTILS         
         if (.not. geom_cart_coords_get(geom, x))
     $        call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)
#endif         
         xmax = 0d0
         do i = 1, ncart
            xmax = max(xmax, abs(dx(i)))
            x(i) = x(i) + dx(i)
         enddo
         if (ipass .eq. 1) then
            call dcopy(ncart, x, 1, xsave, 1)
            xmax1 = xmax
         endif
c
         if (oprint) write(LuOut,77) ipass, xmax, err
 77      format(' Cartesian to internals: iter =',i2,'  dxmax =',
     $        1p,d8.1,'  dqmax =',d8.1)
c
         if (.not. geom_cart_coords_set(geom, x))
     $        call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)
         call sym_geom_project(geom, max(1d-6,xmax*0.1d0))
         if (.not. geom_compute_zmatrix(geom, p))
     $        call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
         if (odebug) then 
            write(LuOut,*) ' Current internals '
            call output(p,1,nzvar,1,1,nzvar,1,1)
         endif
c
c     Get zmat returns bonds in angstrom and angles in degrees
c     (-180..180).  Take care of angles that have changed sign.
c
         do i = 1, nzvar
            dq(i) = q(i) - p(i)
         enddo
         if (odebug) then 
            write(LuOut,*) ' Step in internals before angle restrict'
            call doutput(dq,1,nzvar,1,1,nzvar,1,1)
         endif
         call geom_zmat_sane_step(geom,dq)
         if (odebug) then 
            write(LuOut,*) ' Step in internals after angle restrict'
            call doutput(dq,1,nzvar,1,1,nzvar,1,1)
         endif
c
c     Convert change in internals into hondo internal units
c     (Bohr for bonds, and radians for angles)
c
         call geom_zmat_ico_scale(geom, dq, 1d0/bohr, 1d0/deg)
c
c     Compute the max error in any element
c
         err = 0.0d0
         do i = 1, nzvar
            err = max(err, abs(dq(i)))
         enddo
         if (odebug) write(LuOut,*) ' ERR in internals ', err
c
         if (ipass .eq. 1) err1 = err
c
         if (err.lt.1d-6 .and. xmax.lt.1d-6) goto 100
*         if (xmax.lt.1d-6) goto 100
c
      enddo
c
c     Take first order step here if things are not good
c
      if (err .gt. err1) then
         if (.not. geom_cart_coords_set(geom, xsave))
     $        call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)
         call sym_geom_project(geom, max(0.1d0*xmax1,1d-6))
         err = err1
         if (oprint) write(LuOut,78)
 78      format(' Cartesian to internals: did not converge. Taking',
     $        ' first order step')
      endif
c
c     If the iteration did not converge, then repimpose constraints
c
#ifdef NEW_DRIVER_UTILS
      if (err .gt. 1d-6) then
          call geom_impose_constraints_new(geom,
     $     geom_impose_constraints_on_q)
      else
        call geom_bandbi(geom)
      endif
#else      
      if (err .gt. 1d-6) 
     $     call geom_impose_constraints(geom,
     $     geom_impose_constraints_on_q)
#endif      
c
 100  if (.not. ma_chop_stack(l_bi))
     $     call errquit('opt_int_to_cart: ma chop', 0, MA_ERR)
#ifndef NEW_DRIVER_UTILS      
c
      call geom_bandbi(geom)
#endif      
c
      end
#ifdef NEW_DRIVER_UTILS
      subroutine geom_impose_constraints_new(geom, impose_constraints)
#else      
      subroutine geom_impose_constraints(geom, impose_constraints)
#endif      
      implicit none
#include "errquit.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "nwc_const.fh"
#include "stdio.fh"
#include "util_params.fh"
      integer geom
      external impose_constraints
c     
c     Impose initial values or constraints in internal coordinates
c     upon the cartesian coordinates
c     
c     impose_constraints() is one of 
c     .   geom_impose_constraints_on_q(geom,q)
c     .   geom_impose_initial_values_on_q(geom,q)
c     which given a set of internals munge them to the desired values.
c     
c     Also enforce symmetry.
c     
      integer max_cent, max_nzvar
      parameter (max_cent = nw_max_atom)
      parameter (max_nzvar= 3*max_cent)
      double precision  p(max_nzvar) ! Current internals
      double precision  q(max_nzvar) ! Target internals
      double precision dq(max_nzvar) ! Step in internals
      double precision  x(max_nzvar) ! Cartesian coordinates
      double precision dx(max_nzvar) ! Step in Cartesian coordinates
c     
      integer l_bi, k_bi, i, ipass,nat, ncart, nzvar
      logical odebug, oprint
      logical geom_zmt_get_nzvar, geom_compute_zmatrix
      double precision bohr, deg, err, xmax
c     
      odebug = (ga_nodeid().eq.0).and.util_print('optitoc',print_never)
      oprint = (ga_nodeid().eq.0).and.util_print('xupdate',print_high)
      oprint = oprint .or. odebug
c     
      if (.not. geom_zmt_get_nzvar(geom, nzvar))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      if (.not. geom_ncent(geom, nat))
     $     call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
      ncart = nat * 3
c     
      bohr = cau2ang
      deg  = bohr*180d0/(4d0*atan(1d0))
c     
      if (.not. ma_push_get(mt_dbl,ncart*nzvar,'mem bi',l_bi,k_bi))
     $     call errquit('opt_int_to_cart: ma', ncart*nzvar, MA_ERR)
c     
c     In a non-redundant set of coordinates the iteration will be 
c     quadratically convergent once very large steps have been damped 
c     out ... well this would be true except some torsions and
c     angles can be very strongly coupled.
c
      do ipass = 1, 50
c     
c     Get current q, impose constraints. Compute dq
c     
         if (.not. geom_compute_zmatrix(geom, q))
     $        call errquit('driver_u_c_f_i: geom?',0, GEOM_ERR)
*         write(LuOut,*) ' CURRENT q '
*         call output(q,1,nvar,1,1,nvar,1,1)
         call dcopy(nzvar, q, 1, p, 1)
         call impose_constraints(geom,q) ! q is the target
         do i = 1, nzvar
            dq(i) = q(i) - p(i)
         enddo
         call geom_zmat_sane_step(geom,dq)
         call geom_zmat_ico_scale(geom, dq, 1d0/bohr, 1d0/deg)
*         write(LuOut,*) ' CURRENT dq '
*         call output(dq,1,nvar,1,1,nvar,1,1)
         err = 0.0d0
         do i = 1, nzvar
            err = max(err, abs(dq(i)))
         enddo
#ifdef NEW_DRIVER_UTILS         
         if (err .lt. 1d-6) goto 100
#endif
         if (err .gt. 0.25d0) then ! Small steps move coupled angles together?
            if (oprint) write(LuOut,*) '   restricting step', err
            call dscal(nzvar, 0.5d0/err, dq, 1)
         endif
c     
c     dq -> dx -> x -> new q
c     
         call geom_bandbi(geom)
         call geom_hnd_get_data('b^-1', dbl_mb(k_bi), ncart*nzvar)
         call dgemv('n', ncart, nzvar, 1d0, dbl_mb(k_bi), ncart,
     $        dq, 1, 0.0d0, dx, 1)
         call sym_grad_symmetrize(geom, dx)
         if (.not. geom_cart_coords_get(geom, x))
     $        call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)
         xmax = 0d0
         do i = 1, ncart
            xmax = max(xmax,abs(dx(i)))
            x(i) = x(i) + dx(i)
         enddo
         if (.not. geom_cart_coords_set(geom, x))
     $        call errquit('opt_int_to_cart: geom get/set',0, GEOM_ERR)
         call sym_geom_project(geom, max(1d-6,xmax*0.1d0))
c     
         if (oprint) then
            write(LuOut,33) ipass, err
 33         format(' Imposing constraints:   iter =',i2,'  dqmax =',
     $           1p,d8.1)
         endif
         if (err .lt. 1d-6) goto 100
      enddo
      if (ga_nodeid().eq.0 .and. util_print('warnings',print_low))then
         write(LuOut,44) err
 44      format(/,'!! warning, imposition of constraints did not',
     $        ' converge.  dqmax=',1p,d8.1)
      endif
      call errquit('geom_impose_constraints: failure',0,GEOM_ERR)
c     
 100  if (.not. ma_chop_stack(l_bi))
     $     call errquit('opt_int_to_cart: ma chop', 0, MA_ERR)
c     
      call geom_bandbi(geom)
c     
      end
#ifndef NEW_DRIVER_UTILS
#define NEW_DRIVER_UTILS
#include "geom_driver_utils.F"
#endif
