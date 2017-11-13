C $Id: ecp_int2.F 19696 2010-10-29 16:53:42Z d3y133 $
      subroutine ecp_int2_k(
     R     r_bc,r_ac,pre,
     G     gam,zeta_c,zeta_b,zeta_a,
     B     bet,alp,bi,ai,
     N     n_co,n_prim_b,n_prim_a,i_z)
      implicit none
      integer n_co,n_prim_b,n_prim_a
      integer i_z
      double precision r_bc,r_ac,pre(*)
      double precision gam(*),zeta_c(*),zeta_b(*),zeta_a(*)
      double precision bet(*),alp(*),bi(*),ai(*)
c     
      integer i,i_c,i_b,i_a
      double precision log_prefactor
      double precision rz_b,r2z_b
      i = 1
!DEC$ LOOP COUNT MAX=30, MIN=2
      do i_c = 0,n_co-1
!DEC$ LOOP COUNT MAX=100, MIN=10
         do i_b = 1,n_prim_b
            rz_b=R_bc*zeta_b(i_b)
            r2z_b=rz_b*R_bc
!DEC$ LOOP COUNT MAX=100, MIN=10
            do i_a = 1,n_prim_a
               gam(i) = 1d0/sqrt(zeta_c(i_z+i_c)+zeta_b(i_b)
     &              +zeta_a(i_a))
               bet(i) = rz_b*gam(i)
               alp(i) = R_ac*zeta_a(i_a)*gam(i)
               bi(i) = 0.5d0/rz_b
               ai(i) = 1d0/(2d0*R_ac*zeta_a(i_a))
               pre(i) = exp(
     E              (alp(i)+bet(i))**2
     &              -zeta_a(i_a)*R_ac**2
     &              -r2z_b)
               i = i+1
            enddo
         enddo
      enddo
      return
      end
