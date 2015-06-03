      subroutine ecp_local2_k(
     A     ab,xpc,ypc,zpc,zeta_a,zeta_b,zeta_c,
     b     prod,gam,alp,ai,rpc,pre,
     X     x_ac,x_bc,y_ac,y_bc,z_ac,z_bc,
     I     n_prim_a,n_prim_b,n_coef_c,n_ab
     S     )
      implicit none
#include "ecp_consts.fh"
      double precision ab(*),xpc(*),ypc(*),zpc(*),
     Z     zeta_a(*),zeta_b(*),zeta_c(*),
     b     prod(*),gam(*),alp(*),ai(*),rpc(*),pre(*),
     X     x_ac,x_bc,y_ac,y_bc,z_ac,z_bc
      integer n_prim_a,n_prim_b,n_coef_c,n_ab

c
      double precision rsq,ab_in,log_prefac
      integer i,i_a,i_b,i_c,j
c
      Rsq = (X_AC-X_BC)**2+(Y_AC-Y_BC)**2+(Z_AC-Z_BC)**2
        i = 1
        do i_b = 1,n_prim_b
          do i_a = 1,n_prim_a
            ab_in = zeta_b(i_b)+zeta_a(i_a)
            ab(i) = ab_in
            xpc(i) = (X_AC*zeta_a(i_a)+X_BC*zeta_b(i_b))/ab_in
            ypc(i) = (Y_AC*zeta_a(i_a)+Y_BC*zeta_b(i_b))/ab_in
            zpc(i) = (Z_AC*zeta_a(i_a)+Z_BC*zeta_b(i_b))/ab_in
            rpc(i) = sqrt(xpc(i)**2+ypc(i)**2
     &          +zpc(i)**2)
            prod(i) = Rsq*zeta_b(i_b)*zeta_a(i_a)/ab_in
            i = i+1
          end do
        end do
*
*     Set up argument values for radial integrals
*
        i = 1
        do i_c = 1,n_coef_c
          do j = 1,n_ab
            gam(i) = 1d0/sqrt(zeta_c(i_c)+ab(j))
            alp(i) = rpc(j)*ab(j)*gam(i)
            if (rpc(j) .eq. 0d0) then
              ai(i) = 0d0
            else
              ai(i) = 1d0/(2d0*rpc(j)*ab(j))
            end if
            log_prefac = alp(i)**2-ab(j)*rpc(j)**2
     &          -prod(j)
            pre(i) = 0d0
            if (log_prefac.gt.ln_thr_ecp) 
     &         pre(i) = exp(log_prefac)
            i = i+1
          end do
        end do
        return
        end
