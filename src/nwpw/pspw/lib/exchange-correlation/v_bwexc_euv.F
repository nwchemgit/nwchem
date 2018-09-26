*
* $Id$
*

      

*    ************************************
*    *					*
*    *		v_bwexc_euv	        *
*    *					*
*    ************************************
*
      
      subroutine v_bwexc_euv(gga,n2ft3d,ispin,dn,
     >                       x_parameter,c_parameter,
     >                       stress)
      implicit none
      integer gga
      integer n2ft3d
      integer  ispin
      real*8  dn(n2ft3d,2)
      real*8  x_parameter,c_parameter
      real*8  stress(3,3)


#include "bafdecls.fh"
#include "nwpwxc.fh"
#include "util.fh"


c     **** local variables ****      
      logical value, use_nwpwxc
      integer nx,ny,nz,u,v,s
      real*8  scal1,pi,scal,sum
      integer rho(2),grx(2),gry(2),grz(2)
      integer agr(3),fn(2),fdn(3),tmp(2),rhog(2),xce(2)

      integer rhoup(2),grupx(2),grupy(2),grupz(2)
      integer rhodn(2),grdnx(2),grdny(2),grdnz(2)
      integer          grallx(2),grally(2),grallz(2)
      integer          grad(2),gtmp(2)
      integer xagr(2),xfn(2),xfdn(2)
      real*8 hm(3,3),W(3,3),omega
      real*8 dncut,dumtau
      parameter(dncut = 1.0d-30)

*     **** external functions ****
      integer  G_indx
      real*8   lattice_unitg,lattice_omega
      external G_indx
      external lattice_unitg,lattice_omega

      call nwpw_timing_start(4)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      omega = lattice_omega()


*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

      use_nwpwxc = .false.
      use_nwpwxc = nwpwxc_is_on()
                

*     **********************************
*     ***** restricted calculation *****
*     **********************************
      if (ispin.eq.1) then

c        ***** tempory variables needed rho,grx,gry,grz *****
c        *****                          agr,fn,fdn      *****
        value = BA_push_get(mt_dbl,n2ft3d,'rho', rho(2), rho(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grx',grx(2),grx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'gry',gry(2),gry(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grz',grz(2),grz(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'agr',agr(2),agr(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'fn',fn(2),fn(1))
        value = value.and.
     >        BA_push_get(mt_dbl, 2*n2ft3d,'fdn',fdn(2),fdn(1))
        tmp(1) = fdn(1)
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'rhog',rhog(2),rhog(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'xce',xce(2),xce(1))
      if (.not. value) call errquit('out of stack memory',0,0)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(rho(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(agr(1)),1)



c        ***** calculate rho tmp1=rho(g) ****
         call D3dB_rr_Sum(1,dn(1,1),dn(1,1),dbl_mb(rho(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho(1)))
         call D3dB_r_SMul(1,scal1,dbl_mb(rho(1)),dbl_mb(rhog(1)))
         call D3dB_rc_fft3f(1,dbl_mb(rhog(1)))
         call mask_C(0,dbl_mb(rhog(1)))

c        ***** calculated  grup= grad n ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(grx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(gry(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhog(1)),
     >                      dbl_mb(grz(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(gry(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grz(1)))

c        ***** calculate agr = |grad n| ****
         call D3dB_rr_Sqr(1,dbl_mb(grx(1)),
     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sqr(1,dbl_mb(gry(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         call D3dB_rr_Sqr(1,dbl_mb(grz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         if (use_nwpwxc) then
c
c          Evaluate the functional
c
           call nwpwxc_eval_df(1,n2ft3d,dbl_mb(rho(1)),dbl_mb(agr(1)),
     >                       dumtau,xce,
     >                       dbl_mb(fn(1)),dbl_mb(fdn(1)),dumtau)
cc
cc          Combine (df/d|grad a|) with (df/d(grad a|grad b))
cc
           call D3dB_rr_daxpy(1,0.5d0,dbl_mb(fdn(1)+n2ft3d),
     >                                dbl_mb(fdn(1)))
cc
cc          Calculate the energy density from the functional
c
           do u = 1, n2ft3d
             xce(u) = xce(u)/(dbl_mb(rho(1)-1+u)+dncut)
           enddo
c
c          Tranform the gradient terms from (df/d|grad n|^2) to 
c          (df/d|grad n|)
c
           call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))
           call D3dB_rr_Mul2(1,dbl_mb(agr(1)),dbl_mb(fdn(1)))
c
         else
c         call D3dB_rr_Sqrt(1,dbl_mb(agr(1)),
c     >                       dbl_mb(agr(1)))
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))


    
c        ***** calculate fdn=df/d|grad n|  ****
         if (gga.eq.10) then
         call gen_PBE96_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))

         else if (gga.eq.11) then
         call gen_BLYP_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.12) then
         call gen_revPBE_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.13) then
         call gen_PBEsol_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >                                x_parameter,c_parameter,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.14) then
         call gen_HSE_BW_restricted(n2ft3d,
     >                              dbl_mb(rho(1)),
     >                              dbl_mb(agr(1)),
     >                              x_parameter,c_parameter,
     >                              dbl_mb(xce(1)), !*** not used ??***
     >                              dbl_mb(fn(1)),  !*** not used ??***
     >                              dbl_mb(fdn(1)))
         else if (gga.eq.15) then
         call gen_B3LYP_BW_restricted(n2ft3d,
     >                              dbl_mb(rho(1)),
     >                              dbl_mb(agr(1)),
     >                              x_parameter,c_parameter,
     >                              dbl_mb(xce(1)), !*** not used ??***
     >                              dbl_mb(fn(1)),  !*** not used ??***
     >                              dbl_mb(fdn(1)))
         else if (gga.eq.16) then
         call gen_BEEF_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >        x_parameter,c_parameter,0.6001664769d0,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))
         else if (gga.eq.17) then
         call gen_BEEF_BW_restricted(n2ft3d,
     >                                dbl_mb(rho(1)),
     >                                dbl_mb(agr(1)),
     >        x_parameter,c_parameter,0.0d0,
     >                                dbl_mb(xce(1)), !*** not used ??***
     >                                dbl_mb(fn(1)),  !*** not used ??***
     >                                dbl_mb(fdn(1)))
         else
          call errquit('bad gga',0,0)
         end if
         endif


     
c        ***** calculate df/d|grad n| *(grad n)/|grad n| ****
c         call D3dB_rr_Divide(1,dbl_mb(grx(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(gry(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(gry(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grz(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(gry(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grz(1)))

c         call D3dB_rr_Mul(1,dbl_mb(grx(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(gry(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(gry(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(gry(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grz(1)))

c         call D3dB_r_SMul(1,scal1,dbl_mb(grx(1)),
c     >                            dbl_mb(grx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(gry(1)),
c     >                            dbl_mb(gry(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grz(1)),
c     >                            dbl_mb(grz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(gry(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grz(1)))

         call D3dB_r_Zero_Ends(1,dbl_mb(grx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(gry(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grz(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(gry(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grz(1)))


*        **** W(u,s) = Sum(G) [i*G(u)*dcongj(rhog)*gr(s)] ****
*        ****         where gr(1)=grx,gr(2)=gry,gr(3)=grz ****
          
         !call Pack_c_pack(0,dbl_mb(rhog(1)))
         call mask_C(0,dbl_mb(rhog(1)))
         do u=1,3

*          **** agr = i*G(u)*grx ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grx(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           W(u,1) = sum*omega

*          **** agr = i*G(u)*gry ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(gry(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           W(u,2) = sum*omega

*          **** agr = i*G(u)*grz ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grz(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhog(1)),dbl_mb(agr(1)),sum)
           W(u,3) = sum*omega
         end do




*        **** deallocate temporary memory ****
         value = BA_pop_stack(xce(2))
         value = value.and.BA_pop_stack(rhog(2))
         value = value.and.BA_pop_stack(fdn(2))
         value = value.and.BA_pop_stack(fn(2))
         value = value.and.BA_pop_stack(agr(2))
         value = value.and.BA_pop_stack(grz(2))
         value = value.and.BA_pop_stack(gry(2))
         value = value.and.BA_pop_stack(grx(2))
         value = value.and.BA_pop_stack(rho(2))
         if (.not. value) call errquit('cannot pop stack memory',0,0)



*     *******************************************************
*     ***** unrestricted calculation                    *****
*     *******************************************************
      else

c        ***** tempory variables needed rho,grx,gry,grz *****
c        *****                          agr,fn,fdn      *****
        value = BA_push_get(mt_dbl,n2ft3d,'rhoup', rhoup(2), rhoup(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupx',grupx(2),grupx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupy',grupy(2),grupy(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grupz',grupz(2),grupz(1))

        value = value.and.
     >        BA_push_get(mt_dbl,n2ft3d,'rhodn', rhodn(2), rhodn(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdnx',grdnx(2),grdnx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdny',grdny(2),grdny(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grdnz',grdnz(2),grdnz(1))

        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grallx',grallx(2),grallx(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grally',grally(2),grally(1))
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'grallz',grallz(2),grallz(1))

        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'xagr',xagr(2),xagr(1))
        agr(1) = xagr(1)
        agr(2) = xagr(1) +   n2ft3d
        agr(3) = xagr(1) + 2*n2ft3d
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'grad',grad(2),grad(1))
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'gtmp',gtmp(2),gtmp(1))
        value = value.and.
     >        BA_push_get(mt_dbl, 2*n2ft3d,'xfn',xfn(2),xfn(1))
        fn(1) = xfn(1)
        fn(2) = xfn(1)+n2ft3d
        value = value.and.
     >        BA_push_get(mt_dbl, 3*n2ft3d,'xfdn',xfdn(2),xfdn(1))
        fdn(1) = xfdn(1)
        fdn(2) = xfdn(1) +   n2ft3d
        fdn(3) = xfdn(1) + 2*n2ft3d
        tmp(1) = xfdn(1)
        value = value.and.
     >        BA_push_get(mt_dbl, n2ft3d,'xce',xce(2),xce(1))
        if (.not. value) call errquit('out of stack memory',0,0)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhoup(1)),1)
      call dcopy(n2ft3d,0.0d0,0,dbl_mb(rhodn(1)),1)
      call dcopy(3*n2ft3d,0.0d0,0,dbl_mb(xagr(1)),1)


c        ***** calculate rhoup  ****
         call D3dB_r_Copy(1,dn(1,1),dbl_mb(rhoup(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rhoup(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(rhoup(1)),
c     >                            dbl_mb(rhoup(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(rhoup(1)))
         call D3dB_rc_fft3f(1,dbl_mb(rhoup(1)))
         call mask_C(0,dbl_mb(rhoup(1)))

c        ***** calculate   grup= grad nup ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupy(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhoup(1)),
     >                      dbl_mb(grupz(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grupx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grupy(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grupz(1)))

c        ***** calculate agrup = |grad nup| ****
         call D3dB_rr_Sqr(1,dbl_mb(grupx(1)),
     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sqr(1,dbl_mb(grupy(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

         call D3dB_rr_Sqr(1,dbl_mb(grupz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(1)),
c     >                      dbl_mb(agr(1)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(1)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(1)),
c     >                       dbl_mb(agr(1)))     
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(1)))

c        ***** calculate rhodn  ****
         call D3dB_r_Copy(1,dn(1,2),dbl_mb(rhodn(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rhodn(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(rhodn(1)),
c     >                            dbl_mb(rhodn(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(rhodn(1)))
         call D3dB_rc_fft3f(1,dbl_mb(rhodn(1)))
         call mask_C(0,dbl_mb(rhodn(1)))


c        ***** calculate   grdn = grad ndn ****
         call D3dB_ic_Mul(1,dbl_mb(G_indx(1)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdnx(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(2)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdny(1)))
         call D3dB_ic_Mul(1,dbl_mb(G_indx(3)),
     >                      dbl_mb(rhodn(1)),
     >                      dbl_mb(grdnz(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grdnx(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grdny(1)))
         call D3dB_cr_fft3b(1,dbl_mb(grdnz(1)))

c        ***** calculate agrdn = |grad ndn| ****
         call D3dB_rr_Sqr(1,dbl_mb(grdnx(1)),
     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sqr(1,dbl_mb(grdny(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(2)),
c     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(2)))

         call D3dB_rr_Sqr(1,dbl_mb(grdnz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(2)),
c     >                      dbl_mb(agr(2)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(2)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(2)),
c     >                       dbl_mb(agr(2)))     
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(2)))


c        ***** calculate agr = |grad nup +grad ndn| ****
         call D3dB_rr_Sum(1,dbl_mb(grupx(1)),
     >                      dbl_mb(grdnx(1)),
     >                      dbl_mb(grallx(1)))
         call D3dB_rr_Sum(1,dbl_mb(grupy(1)),
     >                      dbl_mb(grdny(1)),
     >                      dbl_mb(grally(1)))
         call D3dB_rr_Sum(1,dbl_mb(grupz(1)),
     >                      dbl_mb(grdnz(1)),
     >                      dbl_mb(grallz(1)))

         call D3dB_rr_Sqr(1,dbl_mb(grallx(1)),
     >                      dbl_mb(agr(3)))
         call D3dB_rr_Sqr(1,dbl_mb(grally(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(3)),
c     >                      dbl_mb(agr(3)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(3)))

         call D3dB_rr_Sqr(1,dbl_mb(grallz(1)),
     >                      dbl_mb(tmp(1)))

c         call D3dB_rr_Sum(1,dbl_mb(tmp(1)),
c     >                      dbl_mb(agr(3)),
c     >                      dbl_mb(agr(3)))
         call D3dB_rr_Sum2(1,dbl_mb(tmp(1)),dbl_mb(agr(3)))

c         call D3dB_rr_Sqrt(1,dbl_mb(agr(3)),
c     >                       dbl_mb(agr(3)))
         call D3dB_rr_Sqrt1(1,dbl_mb(agr(3)))
         if (use_nwpwxc) then
c
c          Copy |grad nup|->grad(1), |grad n|->grad(2), and |grad ndn|->grad(3)
c
           call D3dB_r_Copy(1,dbl_mb(agr(1)),dbl_mb(grad(1)))
           call D3dB_r_Copy(1,dbl_mb(agr(2)),dbl_mb(grad(1)+2*n2ft3d))
           call D3dB_r_Copy(1,dbl_mb(agr(3)),dbl_mb(grad(1)+n2ft3d))
c
c          Replace |grad x| with |grad x|^2
c
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)))
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)+2*n2ft3d))
           call D3dB_rr_Sqr1(1,dbl_mb(grad(1)+n2ft3d))
c
c          Replace |grad n|^2 with |(grad nup|grad ndn)|^2
c
           call D3dB_rr_Sub2(1,dbl_mb(grad(1)),dbl_mb(grad(1)+n2ft3d))
           call D3dB_rr_Sub2(1,dbl_mb(grad(1)+2*n2ft3d),
     +                         dbl_mb(grad(1)+n2ft3d))
           call D3dB_r_SMul1(1,0.5d0,dbl_mb(grad(1)+n2ft3d))
c
c          Evaluate the functional
c
           call nwpwxc_eval_df(2,n2ft3d,dn,dbl_mb(grad(1)),
     >                       dumtau,xce,
     >                       dbl_mb(fn(1)),dbl_mb(fdn(1)),dumtau)
c
c          Replace f with the energy density f/n
c
           do u = 1, n2ft3d
             xce(u) = xce(u)/(dn(u,1)+dn(u,2)+dncut)
           enddo
c          
c          Replace (df/d|grad nup|^2) with (df/d|grad nup|)
c
           call D3dB_rr_daxpy(1,(-0.5d0),dbl_mb(fdn(2)),dbl_mb(fdn(1)))
           call D3dB_rr_Mul2(1,dbl_mb(agr(1)),dbl_mb(fdn(1)))
           call D3dB_r_SMul1(1,2.0d0,dbl_mb(fdn(1)))
c
c          Replace (df/d|grad ndn|^2) with (df/d|grad ndn|)
c
           call D3dB_rr_daxpy(1,(-0.5d0),dbl_mb(fdn(2)),dbl_mb(fdn(3)))
           call D3dB_rr_Mul(1,dbl_mb(agr(2)),dbl_mb(fdn(3)),
     +                      dbl_mb(gtmp(1)+n2ft3d))
           call D3dB_r_SMul1(1,2.0d0,dbl_mb(gtmp(1)+n2ft3d))
c
c          Replace (df/d|(grad nup|grad ndn)|^2) with (df/d|grad n|)
c
           call D3dB_rr_Mul(1,dbl_mb(agr(3)),dbl_mb(fdn(2)),
     +                      dbl_mb(gtmp(1)+2*n2ft3d))
c
c          Put the results back into fdn
c
           call dcopy(n2ft3d,dbl_mb(gtmp(1)+n2ft3d),1,
     +                       dbl_mb(fdn(2)),1)
           call dcopy(n2ft3d,dbl_mb(gtmp(1)+2*n2ft3d),1,
     +                       dbl_mb(fdn(3)),1)
c
         else
 
c        ***** calculate    ****
c        *****  fdn=(dfx/d|grad nup|,dfx/d|grad ndn|,dfc/d|grad n|) ****
         if (gga.eq.10) then
         call gen_PBE96_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.11) then
         call gen_BLYP_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.12) then
         call gen_revPBE_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.13) then
         call gen_PBEsol_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.14) then
         call gen_HSE_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.15) then
         call gen_B3LYP_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >                                  x_parameter,c_parameter,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.16) then
         call gen_BEEF_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >          x_parameter,c_parameter,0.6001664769d0,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else if (gga.eq.17) then
         call gen_BEEF_BW_unrestricted(n2ft3d,dn,
     >                                  dbl_mb(agr(1)),
     >          x_parameter,c_parameter,0.0d0,
     >                                  dbl_mb(xce(1)), !*** not used ***
     >                                  dbl_mb(fn(1)), !*** not used ***
     >                                  dbl_mb(fdn(1)))
         else 
          call errquit('bad gga',0,0)
         end if
         end if ! use_nwpwxc

         
*        **** calculate df/d|grad nup|* (grad nup)/|grad nup|  ****
*        **** calculate df/d|grad ndn|* (grad ndn)/|grad ndn|  ****
*        **** calculate df/d|grad n|  * (grad n)/|grad n|  ****
c         call D3dB_rr_Divide(1,dbl_mb(grupx(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grupy(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupy(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grupz(1)),
c     >                         dbl_mb(agr(1)),
c     >                         dbl_mb(grupz(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdnx(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdnx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdny(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdny(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grdnz(1)),
c     >                         dbl_mb(agr(2)),
c     >                         dbl_mb(grdnz(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grallx(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grallx(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grally(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grally(1)))
c         call D3dB_rr_Divide(1,dbl_mb(grallz(1)),
c     >                         dbl_mb(agr(3)),
c     >                         dbl_mb(grallz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdnx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdny(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(2)),dbl_mb(grdnz(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grallx(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grally(1)))
         call D3dB_rr_Divide2(1,dbl_mb(agr(3)),dbl_mb(grallz(1)))


c         call D3dB_rr_Mul(1,dbl_mb(grupx(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grupy(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupy(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grupz(1)),
c     >                      dbl_mb(fdn(1)),
c     >                      dbl_mb(grupz(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdnx(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdnx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdny(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdny(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grdnz(1)),
c     >                      dbl_mb(fdn(2)),
c     >                      dbl_mb(grdnz(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grallx(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grallx(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grally(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grally(1)))
c         call D3dB_rr_Mul(1,dbl_mb(grallz(1)),
c     >                      dbl_mb(fdn(3)),
c     >                      dbl_mb(grallz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdnx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdny(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(2)),dbl_mb(grdnz(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grallx(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grally(1)))
         call D3dB_rr_Mul2(1,dbl_mb(fdn(3)),dbl_mb(grallz(1)))

*        **** calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
*        ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****

*        **** calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
*        ****         + (df/d|grad n|  * (grad n)/|grad n|)      ****
c         call D3dB_rr_Sum(1,dbl_mb(grupx(1)),
c     >                      dbl_mb(grallx(1)),
c     >                      dbl_mb(grupx(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grupy(1)),
c     >                      dbl_mb(grally(1)),
c     >                      dbl_mb(grupy(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grupz(1)),
c     >                      dbl_mb(grallz(1)),
c     >                      dbl_mb(grupz(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdnx(1)),
c     >                      dbl_mb(grallx(1)),
c     >                      dbl_mb(grdnx(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdny(1)),
c     >                      dbl_mb(grally(1)),
c     >                      dbl_mb(grdny(1)))
c         call D3dB_rr_Sum(1,dbl_mb(grdnz(1)),
c     >                      dbl_mb(grallz(1)),
c     >                      dbl_mb(grdnz(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallx(1)),dbl_mb(grupx(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grally(1)),dbl_mb(grupy(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallz(1)),dbl_mb(grupz(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallx(1)),dbl_mb(grdnx(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grally(1)),dbl_mb(grdny(1)))
         call D3dB_rr_Sum2(1,dbl_mb(grallz(1)),dbl_mb(grdnz(1)))

c         call D3dB_r_SMul(1,scal1,dbl_mb(grupx(1)),
c     >                            dbl_mb(grupx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grupy(1)),
c     >                            dbl_mb(grupy(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grupz(1)),
c     >                            dbl_mb(grupz(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdnx(1)),
c     >                            dbl_mb(grdnx(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdny(1)),
c     >                            dbl_mb(grdny(1)))
c         call D3dB_r_SMul(1,scal1,dbl_mb(grdnz(1)),
c     >                            dbl_mb(grdnz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupy(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grupz(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdnx(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdny(1)))
         call D3dB_r_SMul1(1,scal1,dbl_mb(grdnz(1)))

         call D3dB_r_Zero_Ends(1,dbl_mb(grupx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grupy(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grupz(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdnx(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdny(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(grdnz(1)))

*        **** put sums in k-space ***
         call D3dB_rc_fft3f(1,dbl_mb(grupx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grupy(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grupz(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdnx(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdny(1)))
         call D3dB_rc_fft3f(1,dbl_mb(grdnz(1)))


*        **** W(u,s) = Sum(G) [i*G(u)*dcongj(rhoup)*grup(s)] ****
*        ****         where grup(1)=grupx,grup(2)=grupy,grup(3)=grupz ****
         !call Pack_c_pack(0,dbl_mb(rhoup(1)))
         call mask_C(0,dbl_mb(rhoup(1)))
         do u=1,3

*          **** agr = i*G(u)*grupx ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grupx(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           W(u,1) = sum*omega

*          **** agr = i*G(u)*grupy ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grupy(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           W(u,2) = sum*omega

*          **** agr = i*G(u)*grupz ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grupz(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhoup(1)),dbl_mb(agr(1)),sum)
           W(u,3) = sum*omega
         end do


*        **** W(u,s) = Sum(G) [i*G(u)*dcongj(rhodn)*grup(s)] ****
*        ****         where grdn(1)=grdnx,grup(2)=grdny,grup(3)=grdnz ****
         !call Pack_c_pack(0,dbl_mb(rhodn(1)))
         call mask_C(0,dbl_mb(rhodn(1)))
         do u=1,3

*          **** agr = i*G(u)*grdnx ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grdnx(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           W(u,1) = W(u,1) + sum*omega

*          **** agr = i*G(u)*grdny ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grdny(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           W(u,2) = W(u,2) + sum*omega

*          **** agr = i*G(u)*grdnz ****
           call D3dB_ic_Mul(1,dbl_mb(G_indx(u)),
     >                        dbl_mb(grdnz(1)),
     >                        dbl_mb(agr(1)))
           call mask_C(0,dbl_mb(agr(1)))
           call D3dB_cc_dot(1,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           !call Pack_c_pack(0,dbl_mb(agr(1)))
           !call Pack_cc_dot(0,dbl_mb(rhodn(1)),dbl_mb(agr(1)),sum)
           W(u,3) = W(u,3) + sum*omega
         end do



*        **** deallocate temporary memory ****
         value = BA_pop_stack(xce(2))
         value = value.and.BA_pop_stack(xfdn(2))
         value = value.and.BA_pop_stack(xfn(2))
         value = value.and.BA_pop_stack(gtmp(2))
         value = value.and.BA_pop_stack(grad(2))
         value = value.and.BA_pop_stack(xagr(2))

         value = value.and.BA_pop_stack(grallz(2))
         value = value.and.BA_pop_stack(grally(2))
         value = value.and.BA_pop_stack(grallx(2))
         value = value.and.BA_pop_stack(grdnz(2))
         value = value.and.BA_pop_stack(grdny(2))
         value = value.and.BA_pop_stack(grdnx(2))
         value = value.and.BA_pop_stack(rhodn(2))
         value = value.and.BA_pop_stack(grupz(2))
         value = value.and.BA_pop_stack(grupy(2))
         value = value.and.BA_pop_stack(grupx(2))
         value = value.and.BA_pop_stack(rhoup(2))
         if (.not. value) call errquit('cannot pop stack memory',0,0)


      end if


*     **** stress(u,v) =  Sum(s){W(u,s)*hm(s,v) }  ****
      do v=1,3
      do u=1,3
        stress(u,v) = 0.0d0
        do s=1,3
           stress(u,v) = stress(u,v) + W(u,s)*hm(s,v)
        end do
      end do
      end do

      
      call nwpw_timing_end(4)
        
      return
      end



 
