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
      i = 1
      do i_c = 0,n_co-1
         do i_b = 1,n_prim_b
            do i_a = 1,n_prim_a
               gam(i) = 1d0/sqrt(zeta_c(i_z+i_c)+zeta_b(i_b)
     &              +zeta_a(i_a))
               bet(i) = R_bc*zeta_b(i_b)*gam(i)
               alp(i) = R_ac*zeta_a(i_a)*gam(i)
               bi(i) = 1d0/(2d0*R_bc*zeta_b(i_b))
               ai(i) = 1d0/(2d0*R_ac*zeta_a(i_a))
               log_prefactor = (alp(i)+bet(i))**2
     &              -zeta_a(i_a)*R_ac**2
     &              -zeta_b(i_b)*R_bc**2
               pre(i) = exp(log_prefactor)
#ifdef DEBUG
               if (debug_arrays) then
                  write (LuOut,'(3i5,1p3e20.10)') i_a,i_b,i_c,
     &                 zeta_a(i_a),zeta_b(i_b),zeta_c(i_z+i_c)
                  write (LuOut,'(10x,i5,1p3e20.10)') i,tmp(i_alp+i),
     &                 tmp(i_bet+i),tmp(i_gam+i)
                  write (LuOut,'(15x,1p2e20.10)') tmp(i_ai+i),
     &                  tmp(i_bi+i)
                  write (LuOut,'(15x,1p2e20.10)') log_prefactor,
     &                 tmp(i_pre+i)
               end if
#endif
               i = i+1
            enddo
         enddo
      enddo
      return
      end
