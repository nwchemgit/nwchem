*
* $Id: paw_mult.F,v 1.11 2003-02-22 03:10:44 bylaska Exp $
*
*
* $Log: not supported by cvs2svn $
* Revision 1.10  2003/02/21 22:37:26  bylaska
* find_comp_coeff subroutine has been added....EJB
*
* Revision 1.9  2003/02/21 19:44:22  bylaska
* Routines for computing the mult_energy_coeff have been added to paw_mult
* ...EJB
*

!**************************************************
!
!       Name: paw_mult_init
!
!       Purpose: initializes  paw_mult
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_init()
      implicit none

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer i,j,k,nfft3d,npack0,nion,nkatm
      integer lm_size,lmax,v_mult_size
      integer ia,ii,jj
      integer l,m
      integer Gx(2),Gy(2),Gz(2),Ylm(2)
      real*8  scal,gg,fourpi,omega
      complex*16 itol,cscal
    
      !**** external functions ****
      integer  control_ncut,G_indx
      real*8   control_rcut,double_factorial,lattice_omega
      real*8   lattice_unita
      external control_ncut,G_indx
      external control_rcut,double_factorial,lattice_omega
      external lattice_unita


      fourpi = 16.0d0*datan(1.0d0)
      omega = lattice_omega()

      !*** allocate paw mult memory from heap ***
      call D3dB_nfft3d(1,nfft3d)
      nion = ion_nion()
      nkatm = ion_nkatm()
      call Pack_npack(0,npack0)
      lmax    = paw_basis_max_mult_l()
      lm_size = (paw_basis_max_mult_l()+1)**2

      !*** allocate gk_smooth,gk,and glm ***
      ok =        my_alloc(mt_dbl,npack0,"gk_smooth",gk_smooth)
      ok = ok.and.my_alloc(mt_dbl,npack0*nkatm,"gk",gk)
      ok = ok.and.my_alloc(mt_dcpl,npack0*lm_size,"g_lm",g_lm)

      !*** allocate paw mult arrays ***
      ok = ok.and.my_alloc(mt_int,nion,"i_v_mult",i_v_mult) !** same as i_paw_qlm **
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",0)

      v_mult_size = 0
      do ii=1,nion
        ia = ion_katm(ii)

        int_mb(i_v_mult(1) + ii - 1) = v_mult_size
        v_mult_size = v_mult_size 
     >              + (paw_basis_mult_l(ia)+1)**2
      end do
      ok = my_alloc(mt_dcpl,v_mult_size,"v_mult",v_mult)
      ok = ok.and.
     >     my_alloc(mt_dcpl,v_mult_size,"comp_coeff",comp_coeff)
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",1)


      !*** allocate self_energy_coeff and mult_energy_coeff arrays ***
      ok = my_alloc(mt_dbl,(lmax+1)*nkatm,
     >              'self_energy_coeff',self_energy_coeff)

      lm_size  = (lmax+1)**2
      ok = ok.and.
     >     my_alloc(mt_dcpl,nion*lm_size*nion*lm_size,
     >              'mult_energy_coeff',mult_energy_coeff)
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",2)






      !**** initialize sigma_smooth and ncut***
      sigma_smooth = control_rcut()
      ncut         = control_ncut()



      !**** initialize gk_smooth, gk, and g_lm ****

      !**** allocate stack memory ****
      ok = MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >     MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'Ylm',Ylm(2),Ylm(1))
      if (.not.ok) 
     > call errquit("paw_mult_init: out of stack memory",2)

      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))

      !**** g_lm ****
      jj = 0
      itol = dcmplx(1.0d0,0.0d0)
      do l=0,paw_basis_max_mult_l()
        k = 2*l+1
        cscal = itol/double_factorial(k)
        itol = itol*dcmplx(0.0d0,-1.0d0)

        do m = -l,l

          !*** generate Ylm ***
          call spher_harmonics_generate(l,m,npack0,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(Ylm(1)))

          !**** (-i)**l*Ylm(k)*|k|**l/(2l+1)!! ****
          do k=1,npack0
            gg = dbl_mb(Gx(1)+k-1)**2
     >         + dbl_mb(Gy(1)+k-1)**2
     >         + dbl_mb(Gz(1)+k-1)**2
            dcpl_mb(g_lm(1)+k-1+(jj)*npack0)
     >            =cscal
     >            *dsqrt(gg)**l
     >            *dcpl_mb(Ylm(1)+k-1)
     >
          end do !*k*
          jj = jj + 1

        end do !*m*
      end do !*l*

      !**** gk_smooth and gk ****
      do k=1,npack0
         gg = dbl_mb(Gx(1)+k-1)**2
     >      + dbl_mb(Gy(1)+k-1)**2
     >      + dbl_mb(Gz(1)+k-1)**2
         scal = 0.25d0 * sigma_smooth**2
         dbl_mb(gk_smooth(1)+k-1) = fourpi*dexp(-scal*gg)
     >                             /omega
         do ia=1,nkatm
           scal = 0.25d0 * paw_basis_sigma(ia)**2
           dbl_mb(gk(1)+k-1+(ia-1)*npack0) = fourpi*dexp(-scal*gg)
     >                                      /omega
         end do !*ia*
      end do !*k*

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(Ylm(2))
      ok = ok.and.MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      if (.not.ok) 
     > call errquit('paw_mult_init:error popping stack',2)



*     **** allocate rcell memory ****
      nshl3d=(2*ncut+1)**3
      ok = my_alloc(mt_dbl,(3*nshl3d),'rcell',rcell)
      if (.not. ok) call errquit('out of heap memory',0)

*     **** get lattice vectors in real space ****
      l=0
      do k=-ncut,ncut
        do j=-ncut,ncut
          do i=-ncut,ncut
             l = l+1
             dbl_mb(rcell(1)+3*(l-1) )
     >                = i*lattice_unita(1,1)
     >                + j*lattice_unita(1,2)
     >                + k*lattice_unita(1,3)
             dbl_mb(rcell(1)+3*(l-1)+1)
     >                = i*lattice_unita(2,1)
     >                + j*lattice_unita(2,2)
     >                + k*lattice_unita(2,3)
             dbl_mb(rcell(1)+3*(l-1)+2)
     >                = i*lattice_unita(3,1)
     >                + j*lattice_unita(3,2)
     >                + k*lattice_unita(3,3)
          end do
        end do
      end do


      !*** intitalize self_energy_coeff and mult_energy_coeff ****
      call find_self_energy_coeff(lmax,nkatm,self_energy_coeff)
      call paw_set_mult_energy_coeff() !*needs to be recalled when geometry changes*


      return
      end


*************************************************
!
!       Name: paw_mult_end
!
!       Purpose: deallocates heap memory
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_end()
      implicit none

#include "paw_mult_data.fh"
#include "paw_ma.fh"

      !**** local varables ****
      logical ok

      ok =        my_dealloc(i_v_mult)
      ok = ok.and.my_dealloc(v_mult)
      ok = ok.and.my_dealloc(comp_coeff)
      ok = ok.and.my_dealloc(g_lm)
      ok = ok.and.my_dealloc(gk)
      ok = ok.and.my_dealloc(gk_smooth)
      ok = ok.and.my_dealloc(self_energy_coeff)
      ok = ok.and.my_dealloc(mult_energy_coeff)
      ok = ok.and.my_dealloc(rcell)
      if (.not.ok)
     > call errquit("paw_mult_end: error freeing heap",0)

      return
      end 




!**************************************************
!
!       Name: paw_mult_dn_cmp_get
!
!       Purpose: returns dn_cmp and dn_cmp_smooth 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_dn_cmp_get(dn_cmp,
     >                               dn_cmp_smooth)
      implicit none
      complex*16 dn_cmp(*)
      complex*16 dn_cmp_smooth(*)

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer ia,ii,jj,kk,l,m,nion,npack0,mult_l
      integer exi(2),tmp(2),QYlm(2)
      real*8  sum

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'QYlm',QYlm(2),QYlm(1))
      tmp(1) = exi(1)
      tmp(2) = exi(2)
      if (.not.ok)
     > call errquit(
     >  'paw_mult_dn_cmp_get: out of stack memory',0)


      call dcopy(2*npack0,0.0d0,0,dn_cmp,1)
      call dcopy(2*npack0,0.0d0,0,dn_cmp_smooth,1)
      sum = 0.0d0
      nion = ion_nion()
      do ii=1,nion
         ia = ion_katm(ii)

         !**** Define QYlm ****
         mult_l = paw_basis_mult_l(ia)
         jj     = int_mb(i_paw_qlm(1)+ii-1)
         write(*,*) "Qlm(0):",ii,dcpl_mb(paw_qlm(1)+jj)
         sum = sum + dble(dcpl_mb(paw_qlm(1)+jj))
         kk     = 0
         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QYlm(1)),1)
         do l=0,mult_l
         do m=-l,l
           call Pack_cc_zaxpy(0,
     >               dcpl_mb(paw_qlm(1)+jj),
     >               dcpl_mb(g_lm(1)+(kk)*npack0),
     >               dcpl_mb(QYlm(1)))

           jj = jj + 1
           kk = kk + 1
         end do
         end do

         !**** Multiply by Structure Factor ****
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         call Pack_cc_Mul(0,
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(QYlm(1)))


         !**** add up ncmp_smooth^ii  ****
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
         call Pack_cc_Sum(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp_smooth,
     >                    dn_cmp_smooth)

         !**** add up ncmp^ii ***
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
         call Pack_cc_Sum(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp,
     >                    dn_cmp)


      end do !*ii*

      write(*,*) "Qtotal:",sum*dsqrt(16.0d0*datan(1.0d0))
      sum = 2*dbl_mb(paw_q00(1)) + dbl_mb(paw_q00(1)+1)
      write(*,*) "Qcore :",sum*dsqrt(16.0d0*datan(1.0d0))

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(QYlm(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_dn_cmp_get: error popping stack',1)

      return
      end

!**************************************************
!
!       Name: paw_mult_coeff_set
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_coeff_set(vh,vcmp)
      implicit none
      complex*16 vh(*)
      complex*16 vcmp(*)

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"


      !**** local variables ***
      logical ok
      integer ia,ii,ja,jj,kk,l,m,nion,npack0,mult_l,mabs
      integer il,li,mi,mult_li,ill,jll
      integer jl,lj,mj,mult_lj
      integer  isgn,lmax,lmax2,indx
      integer exi(2),tmp1(2),gls(2),gl(2),t_mult(2)
      real*8  omega
      complex*16 csum1,csum2

      !**** external functions ****
      real*8   lattice_omega
      external lattice_omega

      omega = lattice_omega()

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = MA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'tmp1',tmp1(2),tmp1(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'gl',gl(2),gl(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,npack0,'gls',gls(2),gls(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,v_mult(3),
     >                 't_mult',t_mult(2),t_mult(1))
      if (.not.ok)
     > call errquit(
     >  'paw_mult_coeff_set: out of stack memory',0)

      nion = ion_nion()
      do ii=1,nion
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
         ia = ion_katm(ii)

         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(gl(1)))
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(gls(1)))

         mult_l   = paw_basis_mult_l(ia)
         jj       = int_mb(i_v_mult(1) + ii - 1) 
         kk       = 0
         do l=0,mult_l
           do m=-l,l
             call Pack_cc_Mul(0,
     >                    dcpl_mb(gl(1)),
     >                    dcpl_mb(g_lm(1)+(kk)*npack0),
     >                    dcpl_mb(tmp1(1)))
             call Pack_cc_izdot(0,
     >                        vh,
     >                        dcpl_mb(tmp1(1)),
     >                        csum1)
             call Pack_cc_Mul(0,
     >                    dcpl_mb(gls(1)),
     >                    dcpl_mb(g_lm(1)+(kk)*npack0),
     >                    dcpl_mb(tmp1(1)))
             call Pack_cc_izdot(0,
     >                        vcmp,
     >                        dcpl_mb(tmp1(1)),
     >                        csum2)

             !**** v_mult(l,m,ii) = <g_lm^a|vh> + <tg_lm^a|vcmp> ****
             dcpl_mb(t_mult(1)+jj) = (csum1 + csum2)*omega

             jj = jj + 1
             kk = kk + 1
           end do !*m*
         end do !*l*
      end do !*ii*


      !**** unscramble multipole coefficients ****
      write(*,*)
      write(*,*)
      do ii=1,nion
         ia = ion_katm(ii)
         mult_l   = paw_basis_mult_l(ia)
         jj       = int_mb(i_v_mult(1) + ii - 1) 
         do l=0,mult_l
         do m=-l,l
            mabs = abs(m)
            if (mod(mabs,2).eq.0) then
              isgn = 1
            else
              isgn = -1
            end if
           kk = jj - 2*m
           dcpl_mb(v_mult(1)+jj) = 0.5d0*(dcpl_mb(t_mult(1)+jj) 
     >                           + isgn*dconjg(dcpl_mb(t_mult(1)+kk)))

           write(*,*) "v_mult 1:",ii,l,m,dcpl_mb(v_mult(1)+jj)
           jj = jj + 1
         end do
         end do
      end do
      call D3dB_Vector_SumAll(2*v_mult(3),dcpl_mb(v_mult(1)))
      write(*,*)
      write(*,*)

      write(*,*)
      write(*,*)
      lmax  = paw_basis_max_mult_l()
      lmax2 = (lmax+1)**2
      do ii=1,nion
         ia = ion_katm(ii)
         mult_li  = paw_basis_mult_l(ia)
         il       = int_mb(i_v_mult(1) + ii - 1) 
         ill      = 0
         do li=0,mult_li
         do mi=-li,li

           do jj=1,nion
             ja = ion_katm(jj)
             mult_lj  = paw_basis_mult_l(ja)
             jl       = int_mb(i_paw_qlm(1)+jj-1)
             jll       = 0
             do lj=0,mult_lj
             do mj=-lj,lj

               indx = (ii-1) 
     >              + ill*nion
     >              + (jj-1)*nion*lmax2
     >              + jll*nion*lmax2*nion

               dcpl_mb(v_mult(1)+il) 
     >         = dcpl_mb(v_mult(1)+il) 
     >         + dcpl_mb(paw_qlm(1)+jl)
     >           *dcpl_mb(mult_energy_coeff(1)+indx)

       write(*,*) "v_mult 1.5:",ii,ill+1,jj,jll+1,
     >            dcpl_mb(paw_qlm(1)+jl),
     >            dcpl_mb(mult_energy_coeff(1)+indx)

               jl  = jl+1
               jll = jll+1
             end do !*mj*
             end do !*lj*
           end do !*jj*

           il = il + 1
           ill = ill + 1
         end do !*mi*
         end do !*li*
      end do !*ii*
      write(*,*)
      write(*,*)

      !*** DEBUG ***
      write(*,*)
      write(*,*)
      do ii=1,nion
         ia = ion_katm(ii)
         mult_li  = paw_basis_mult_l(ia)
         il       = int_mb(i_v_mult(1) + ii - 1) 
         do li=0,mult_li
         do mi=-li,li
           write(*,*) "v_mult 2:",ii,li,mi,dcpl_mb(v_mult(1)+il)
           il = il + 1
         end do !*mi*
         end do !*li*
      end do !*ii*
      write(*,*)
      write(*,*)
      !*** DEBUG ***


      call find_comp_coeff()

      write(*,*)
      write(*,*)
      do ii=1,nion
         ia = ion_katm(ii)
         mult_li  = paw_basis_mult_l(ia)
         il       = int_mb(i_v_mult(1)+ii-1) 
         do li=0,mult_li
           indx = li + (ia-1)*(lmax+1)
           do mi=-li,li
             dcpl_mb(v_mult(1)+il) 
     >         = dcpl_mb(v_mult(1)+il) 
     >         - dconjg(dcpl_mb(paw_qlm(1)+il))
     >           *dbl_mb(self_energy_coeff(1)+indx)
     >         + dcpl_mb(comp_coeff(1)+il)

       write(*,*) "v_mult 2.5:",ii,li,mi,dcpl_mb(paw_qlm(1)+il),
     >     dbl_mb(self_energy_coeff(1)+indx),
     >     dcpl_mb(comp_coeff(1)+il)

             il = il + 1
           end do !*mi*
         end do !*li*
      end do !*ii*
      write(*,*)
      write(*,*)

      !*** DEBUG ***
      write(*,*)
      write(*,*)
      do ii=1,nion
         ia = ion_katm(ii)
         mult_li  = paw_basis_mult_l(ia)
         il       = int_mb(i_v_mult(1) + ii - 1) 
         do li=0,mult_li
         do mi=-li,li
           write(*,*) "v_mult 3:",ii,li,mi,dcpl_mb(v_mult(1)+il)
           il = il + 1
         end do !*mi*
         end do !*li*
      end do !*ii*
      write(*,*)
      write(*,*)
      !*** DEBUG ***

      !**** deallocate stack memory ****
      ok =        MA_pop_stack(t_mult(2))
      ok = ok.and.MA_pop_stack(gls(2))
      ok = ok.and.MA_pop_stack(gl(2))
      ok = ok.and.MA_pop_stack(tmp1(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_coeff_set: error popping stack',1)

      return
      end

!**************************************************
!
!       Name: paw_mult_vzero
!
!       Purpose: 
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_vzero(vzero)
      implicit none
      real*8 vzero

#include "mafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ****
      integer ia,ii,jj,nion
      real*8  fourpi
 
      !**** external functions ****
      real*8   lattice_omega
      external lattice_omega

      fourpi = 16.0d0*datan(1.0d0)
      nion = ion_nion()
      vzero = 0.0d0
      do ii=1,nion
         ia = ion_katm(ii)
         jj = int_mb(i_paw_qlm(1)+ii-1)
         vzero = vzero 
     >         + dble(dcpl_mb(paw_qlm(1)+jj))
     >          *(sigma_smooth**2-paw_basis_sigma(ia)**2)
      end do
      vzero = vzero*fourpi*dsqrt(fourpi)/lattice_omega()/4.0d0

      return
      end


  !*************************************************
  !
  !   Name    : find_self_energy_coeff
  !
  !   Purpose :
  !
  !   Created :
  !*************************************************
      subroutine find_self_energy_coeff(lmax,nkatm,coeff)
      implicit none
      integer lmax,nkatm
      real*8 coeff(lmax+1,nkatm)

#include "paw_basis.fh"

      !*** local variables ***
      integer ia,l,mult_l
      real*8 sigma_tmp,twopi

      !*** external functions ***
      integer  paw_double_factorial
      external paw_double_factorial

      twopi = 8.0d0*datan(1.0d0)
      call dcopy(lmax*nkatm,0.0d0,0,coeff,1)

      do ia=1,nkatm
        sigma_tmp = paw_basis_sigma(ia)
        mult_l    = paw_basis_mult_l(ia)
        do l=0,mult_l

          coeff(l+1,ia) =  4.0d0*dsqrt(twopi)
     >     /(dble((2*l+1)*paw_double_factorial(2*l+1))
     >      *sigma_tmp**(2*l+1) )
        end do
      end do

      return
      end 

  !*************************************************
  !
  !   Name    : find_comp_coeff
  !
  !   Purpose :
  !
  !   Created :
  !*************************************************
      subroutine find_comp_coeff()
      implicit none

#include "mafdecls.fh"
#include "paw_gaunt.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_mult_data.fh"
#include "paw_matrix_comp_pot.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"


      !*** local variables **
      integer i,j,ii,ia,jj,il,jl
      integer l,m,li,mi,lj,mj,mult_l
      integer indx,indx_2,mtrx_ptr,mtrx_ptr2
      integer nion,nbasis,basis_nbasis


      call dcopy(2*comp_coeff(3),0.0d0,0,
     >           dcpl_mb(comp_coeff(1)),1)

      nion = ion_nion()
      do ii=1,nion
        ia = ion_katm(ii)
        mult_l       = paw_basis_mult_l(ia)
        basis_nbasis = paw_basis_nbasis(ia)
        nbasis       = paw_proj_nbasis(ia)
       
        jj        = int_mb(i_v_mult(1)+ii-1)
        mtrx_ptr  = int_mb(i_paw_comp_pot_matrix(1)+ia-1)
        mtrx_ptr2 = int_mb(i_paw_ovlp_w(1)+ii-1)
        do l=0,mult_l
        do m=-l,l

          jl = 0
          do j=1,basis_nbasis
            lj=paw_basis_orb_l(j,ia)
            do mj=-lj,lj

            il = 0
            do i=1,basis_nbasis
              li=paw_basis_orb_l(i,ia)
              do mi=-li,li

              !*** check for non-zero gaunt coefficient ***
              if ( (l.le.(li+lj))   .and. 
     >             (l.ge.abs(li-lj)).and.
     >             (m.eq.(mi-mj)) ) then
                indx = mtrx_ptr + (i-1) 
     >                          + (j-1)*basis_nbasis
     >                          +     l*basis_nbasis**2
                indx_2 = mtrx_ptr2 + (il-1) 
     >                             + (jl-1)*nbasis 

                dcpl_mb(comp_coeff(1)+jj)
     >            = dcpl_mb(comp_coeff(1)+jj)
     >            + gen_gaunt_coeff(l,m,li,mi,lj,mj)
     >             *dcpl_mb(paw_comp_pot_matrix(1)+indx)
     >             *dcpl_mb(paw_ovlp_w(1)+indx_2)
              end if !*non-zero gaunt*

              il = il + 1
              end do !*mi*
            end do !*i*

            jl = jl + 1
            end do !*mj*
          end do !*j*

          jj = jj + 1
        end do !*m*
        end do !*l*

      end do !*ii*

      return
      end 
