*
* $Id$
*
*
* $Log: not supported by cvs2svn $
* Revision 1.27  2006/01/12 00:54:01  bylaska
* Added charge component analysis wavefunctions, i.e. the following is now outputed in the results section:
*
*  number of electrons: spin up=    8.00000  down=    6.00000 (real space)
*      plane-wave part:             7.96734           6.02652 (real space)
*       augmented part:             0.03266          -0.02652 (real space)
*
* ...EJB
*
* Revision 1.26  2005/02/23 21:40:18  edo
* fixed x1 fpe caused by 0**0
*
* Revision 1.25.4.1  2005/02/23 21:39:50  edo
* fixed x1 fpe caused by 0**0
*
* Revision 1.25  2003/10/28 19:50:51  edo
* errquizzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
*
* Revision 1.24  2003/10/21 02:05:17  marat
* switched to new errquit by running global replace operation
* see the script below (note it will not work on multiline errquit calls)
* *********************************************************
* #!/bin/sh
*
* e=`find . -name "*F" -print`
*
* for f in $e
* do
* cp $f $f.bak
* sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
* #rm $f.bak
* done
* **********************************************************
*
* Revision 1.23  2003/03/21 23:41:13  bylaska
*
* paw updates ...EJB
*
* Revision 1.22  2003/03/11 17:56:03  bylaska
* rcut algorithm has been reinstated....EJB
*
* Revision 1.21  2003/03/05 23:16:32  bylaska
* Commented out write statements and other minor fixes.....
* self-consistent loop looks like it is working.....
* ....EJB
*
* Revision 1.20  2003/03/05 20:35:04  bylaska
* bug fix.....(m.eq.(mi-mj)) changed to (m.eq.(mj-mi)) in the  paw_pot_mult loop...
* Eigevalues now agree!
*
* ....EJB
*
* Revision 1.19  2003/03/04 00:04:04  marat
* added printouts for atomic potenitials
* for debug purposes
* MV
*
* Revision 1.18  2003/02/27 02:19:25  bylaska
* The electron gradient has been added.
* The eigenvalues are not yet agreeing with the F90 code....EJB
*
* Revision 1.17  2003/02/26 20:34:06  marat
* fixed bug related to calculation of
* comp_coeff, the i and j loop were incorrectly switched
* MV
*
* Revision 1.16  2003/02/25 00:21:00  bylaska
* debug write statements commented out...EJB
*
* Revision 1.15  2003/02/24 20:59:58  bylaska
* paw_mult_rcut and paw_mult_ncut functions have been added....EJB
*
* Revision 1.14  2003/02/24 20:52:52  bylaska
* fixed initialization of rcut and ncut.....EJB
*
* Revision 1.13  2003/02/23 21:37:07  bylaska
* routines for calculating atomic multipole energies have been added....EJB
*
* Revision 1.12  2003/02/23 20:53:08  bylaska
* bug fix in find_comp_coeff ....EJB
*
* Revision 1.11  2003/02/22 03:10:44  bylaska
* debugging multipole coefficients...There is currently a bug in
* find_comp_coeff...EJB
*
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

#include "bafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer i,j,k,nfft3d,npack0,nion,nkatm
      integer lm_size,lmax,v_mult_size
      integer ia,ii,jj
      integer l,m
      integer Gx(2),Gy(2),Gz(2),Ylm(2)
      integer paw_pot_mult_size
      real*8  scal,gg,fourpi,omega,rs,w
      complex*16 itol,cscalf
    
      !**** external functions ****
      integer  control_ncut,G_indx,control_version
      real*8   control_rcut,double_factorial,lattice_omega
      real*8   lattice_unita
      external control_ncut,G_indx,control_version
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
     > call errquit("paw_mult_init: out of heap memory",0,0)

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
     > call errquit("paw_mult_init: out of heap memory",0,1)


      !*** allocate self_energy_coeff and mult_energy_coeff arrays ***
      ok = my_alloc(mt_dbl,(lmax+1)*nkatm,
     >              'self_energy_coeff',self_energy_coeff)

      lm_size  = (lmax+1)**2
      ok = ok.and.
     >     my_alloc(mt_dcpl,nion*lm_size*nion*lm_size,
     >              'mult_energy_coeff',mult_energy_coeff)
      if (.not.ok) 
     > call errquit("paw_mult_init: out of heap memory",0,2)

      !*** allocate multiple potential nonlocal matrix ***
      ok = my_alloc(mt_int,nion,"i_paw_pot_mult",
     >                                     i_paw_pot_mult)
      if (.not.ok)
     > call errquit('init_paw_pot_mult:out of heap memory',0,3)

      paw_pot_mult_size = 0
      do ii=1,nion
         int_mb(i_paw_pot_mult(1)+ii-1) = paw_pot_mult_size
         ia = ion_katm(ii)
         paw_pot_mult_size = paw_pot_mult_size
     >                     + paw_proj_nbasis(ia)**2
      end do
      ok = my_alloc(mt_dcpl,paw_pot_mult_size,
     >               "paw_pot_mult",paw_pot_mult)
      if (.not.ok)
     > call errquit("paw_mult_init:out of heap memory",0,4)





      !**** initialize sigma_smooth and ncut***
      sigma_smooth = control_rcut()
      if (control_version().eq.3) then
         ncut         = control_ncut()
      else
         ncut         = 0
         if (sigma_smooth.le.0.0d0) sigma_smooth=1.0d0
      end if

      if (ncut.lt.0)     ncut=0
      if (sigma_smooth.le.0.0d0) then
         rs = lattice_unita(1,1)**2
     >      + lattice_unita(2,1)**2
     >      + lattice_unita(3,1)**2
         rs = dsqrt(rs)
         sigma_smooth=4.0d0*rs/fourpi

         rs = lattice_unita(1,2)**2
     >      + lattice_unita(2,2)**2
     >      + lattice_unita(3,2)**2
         rs = dsqrt(rs)
         w=4.0d0*rs/fourpi
         if (w.lt.sigma_smooth) sigma_smooth = w

         rs = lattice_unita(1,3)**2
     >      + lattice_unita(2,3)**2
     >      + lattice_unita(3,3)**2
         rs = dsqrt(rs)
         w=4.0d0*rs/fourpi
         if (w.lt.sigma_smooth) sigma_smooth = w
      end if



      !**** initialize gk_smooth, gk, and g_lm ****

      !**** allocate stack memory ****
      ok = BA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >     BA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >     BA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'Ylm',Ylm(2),Ylm(1))
      if (.not.ok) 
     > call errquit("paw_mult_init: out of stack memory",0,2)

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
        cscalf = itol/double_factorial(k)
        itol = itol*dcmplx(0.0d0,-1.0d0)

        do m = -l,l

          !*** generate Ylm ***
          call spher_harmonics_generate(l,m,npack0,
     >                      dbl_mb(Gx(1)),
     >                      dbl_mb(Gy(1)),
     >                      dbl_mb(Gz(1)),
     >                      dcpl_mb(Ylm(1)))

          !**** (-i)**l*Ylm(k)*|k|**l/(2l+1)!! ****
          if(l.eq.0) then
          do k=1,npack0
            dcpl_mb(g_lm(1)+k-1+(jj)*npack0)
     >            =cscalf
     >            *dcpl_mb(Ylm(1)+k-1)
          enddo
          else
          do k=1,npack0
            gg = dbl_mb(Gx(1)+k-1)**2
     >         + dbl_mb(Gy(1)+k-1)**2
     >         + dbl_mb(Gz(1)+k-1)**2
            dcpl_mb(g_lm(1)+k-1+(jj)*npack0)
     >            =cscalf
     >            *dsqrt(gg)**l
     >            *dcpl_mb(Ylm(1)+k-1)
     >
          end do !*k*
          endif
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
      ok =        BA_pop_stack(Ylm(2))
      ok = ok.and.BA_pop_stack(Gz(2))
      ok = ok.and.BA_pop_stack(Gy(2))
      ok = ok.and.BA_pop_stack(Gx(2))
      if (.not.ok) 
     > call errquit('paw_mult_init:error popping stack',0,2)



*     **** allocate rcell memory ****
      nshl3d=(2*ncut+1)**3
      ok = my_alloc(mt_dbl,(3*nshl3d),'rcell',rcell)
      if (.not. ok) call errquit('out of heap memory',0,0)

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
      call find_self_energy_coeff(lmax,nkatm,
     >                            dbl_mb(self_energy_coeff(1)))
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
      ok = ok.and.my_dealloc(paw_pot_mult)
      ok = ok.and.my_dealloc(i_paw_pot_mult)
      if (.not.ok)
     > call errquit("paw_mult_end: error freeing heap",0,0)

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

#include "bafdecls.fh"
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
      ok = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'QYlm',QYlm(2),QYlm(1))
      tmp(1) = exi(1)
      tmp(2) = exi(2)
      if (.not.ok)
     > call errquit(
     >  'paw_mult_dn_cmp_get: out of stack memory',0,0)


      call dcopy(2*npack0,0.0d0,0,dn_cmp,1)
      call dcopy(2*npack0,0.0d0,0,dn_cmp_smooth,1)
      sum = 0.0d0
      nion = ion_nion()
      do ii=1,nion
         ia = ion_katm(ii)

         !**** Define QYlm ****
         mult_l = paw_basis_mult_l(ia)
         jj     = int_mb(i_paw_qlm(1)+ii-1)
         sum    = sum + dble(dcpl_mb(paw_qlm(1)+jj))
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
c         call Pack_cc_Mul(0,
c     >                    dcpl_mb(exi(1)),
c     >                    dcpl_mb(QYlm(1)),
c     >                    dcpl_mb(QYlm(1)))
         call Pack_cc_Mul2(0,dcpl_mb(exi(1)),
     >                       dcpl_mb(QYlm(1)))


         !**** add up ncmp_smooth^ii  ****
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
c         call Pack_cc_Sum(0,
c     >                    dcpl_mb(tmp(1)),
c     >                    dn_cmp_smooth,
c     >                    dn_cmp_smooth)
         call Pack_cc_Sum2(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp_smooth)

         !**** add up ncmp^ii ***
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
c         call Pack_cc_Sum(0,
c     >                    dcpl_mb(tmp(1)),
c     >                    dn_cmp,
c     >                    dn_cmp)
         call Pack_cc_Sum2(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp)


      end do !*ii*


      !**** deallocate stack memory ****
      ok =        BA_pop_stack(QYlm(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_dn_cmp_get: error popping stack',0,1)

      return
      end


!**************************************************
!
!       Name: paw_mult_pw_force
!
!       Purpose: returns 
!           Sum(G) (i*G*vh(G)*ncmp^a(G)+i*G*vcmp(G)*n_cmp_smooth^a(G))
!
!       Created:        2/16/2003
!**************************************************
      subroutine paw_mult_pw_force(vh,vcmp,fion)
      implicit none
      complex*16 vh(*)
      complex*16 vcmp(*)
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer i,ia,ii,jj,kk,l,m,nion,npack0,nfft3d,mult_l
      integer exi(2),QYlm(2),ncmp(2),ncmp_smooth(2)
      integer xtmp(2),Gx(2),Gy(2),Gz(2),G(3)
      real*8  fx,fy,fz,omega

      !*** external functions
      integer  G_indx
      real*8   lattice_omega
      external G_indx
      external lattice_omega


      omega = lattice_omega()

      !**** allocate stack memory ****
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      ok = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'QYlm',QYlm(2),QYlm(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'ncmp',ncmp(2),ncmp(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,
     >     'ncmp_smooth',ncmp_smooth(2),ncmp_smooth(1))
      ok = ok.and.
     >       BA_push_get(mt_dbl, npack0,'xtmp',xtmp(2),xtmp(1))
      ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not.ok)
     > call errquit(
     >  'paw_mult_pw_force: out of stack memory',0,0)

      !**** define Gx,Gy and Gz in packed space ****
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)
      call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))


      nion = ion_nion()
      do ii=1,nion
         ia = ion_katm(ii)

         !**** Define QYlm ****
         mult_l = paw_basis_mult_l(ia)
         jj     = int_mb(i_paw_qlm(1)+ii-1)
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
c         call Pack_cc_Mul(0,
c     >                    dcpl_mb(exi(1)),
c     >                    dcpl_mb(QYlm(1)),
c     >                    dcpl_mb(QYlm(1)))
         call Pack_cc_Mul2(0,
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(QYlm(1)))


         !**** add up ncmp_smooth^ii  ****
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(ncmp_smooth(1)))

         !**** add up ncmp^ii ***
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk(1)+(ia-1)*npack0),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(ncmp(1)))
          do i=1,npack0
             dbl_mb(xtmp(1)+i-1)
     >          = dimag(vh(i))* dble(dcpl_mb(ncmp(1)+i-1))
     >          - dble(vh(i))*dimag(dcpl_mb(ncmp(1)+i-1))
     >          + dimag(vcmp(i))* dble(dcpl_mb(ncmp_smooth(1)+i-1))
     >          - dble(vcmp(i))*dimag(dcpl_mb(ncmp_smooth(1)+i-1))
          end do

         call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fx)
         call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fy)
         call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fz)
         fion(1,ii) = fion(1,ii) + fx*omega
         fion(2,ii) = fion(2,ii) + fy*omega
         fion(3,ii) = fion(3,ii) + fz*omega

      end do !*ii*

      !**** deallocate stack memory ****
      ok =        BA_pop_stack(Gz(2))
      ok = ok.and.BA_pop_stack(Gy(2))
      ok = ok.and.BA_pop_stack(Gx(2))
      ok = ok.and.BA_pop_stack(xtmp(2))
      ok = ok.and.BA_pop_stack(ncmp_smooth(2))
      ok = ok.and.BA_pop_stack(ncmp(2))
      ok = ok.and.BA_pop_stack(QYlm(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_pw_force: error popping stack',0,1)

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

#include "bafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_comp_charge_matrix.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"


      !**** local variables ***
      logical ok
      integer ia,ii,ja,jj,kk,l,m,nion,npack0,mult_l,mabs
      integer i,il,li,mi,mult_li,ill,jll
      integer j,jl,lj,mj,mult_lj
      integer  isgn,lmax,lmax2,indx
      integer exi(2),tmp1(2),gls(2),gl(2),t_mult(2)
      real*8  omega
      complex*16 csum1,csum2

      integer i_mtr,i_mtr0
      integer nb,nb2
      integer ilm
      integer i_cp,i_cp0
      integer nilm,njlm
      integer i_qlm,i_qlm0
      complex*16 tmp_mult_pot

      !**** external functions ****
      real*8   lattice_omega,gen_gaunt_coeff
      external lattice_omega,gen_gaunt_coeff

      omega = lattice_omega()

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'tmp1',tmp1(2),tmp1(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'gl',gl(2),gl(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'gls',gls(2),gls(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,v_mult(3),
     >                 't_mult',t_mult(2),t_mult(1))
      if (.not.ok)
     > call errquit(
     >  'paw_mult_coeff_set: out of stack memory',0,0)

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
      !write(*,*)
      !write(*,*)
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

           !write(*,*) "v_mult 1:",ii,l,m,dcpl_mb(v_mult(1)+jj)
           jj = jj + 1
         end do
         end do
      end do
      call D3dB_Vector_SumAll(2*v_mult(3),dcpl_mb(v_mult(1)))



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





      call find_comp_coeff()

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


             il = il + 1
           end do !*mi*
         end do !*li*
      end do !*ii*






      !**** multipole potential nonlocal operator ****
      call dcopy(2*paw_pot_mult(3),
     >           0.0d0,0,
     >           dcpl_mb(paw_pot_mult(1)),1)

      do ii=1,nion
        ia     = ion_katm(ii)
        nb     = paw_basis_nbasis(ia)
        nb2    = nb*nb
        mult_l = paw_basis_mult_l(ia)

        i_mtr0 = int_mb(i_comp_charge_matrix(1) + ia - 1)
        i_qlm0 = int_mb(i_v_mult(1) + ii - 1)
        i_cp0  = int_mb(i_paw_pot_mult(1) + ii - 1)

        i_qlm = i_qlm0
        do l=0,mult_l
        do m=-l,l

          i_cp = i_cp0
          do j=1,nb
          lj = paw_basis_orb_l(j,ia)
          do mj=-lj,lj

          do i=1,nb
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
 
            if ( (l.le.(li+lj)   )    .and.
     >           (l.ge.abs(li-lj))    .and.
     >           (m.eq.(mj-mi)   )  ) then


             i_mtr = i_mtr0 + (i-1) + (j-1)*nb + l*nb2

             tmp_mult_pot = dbl_mb(comp_charge_matrix(1)+i_mtr)*
     >               gen_gaunt_coeff(l,m,lj,mj,li,mi)*
     >               dcpl_mb(v_mult(1) + i_qlm)

                  dcpl_mb(paw_pot_mult(1)+i_cp) =
     >                  dcpl_mb(paw_pot_mult(1)+i_cp)+
     >                  tmp_mult_pot
            end if

            i_cp = i_cp + 1
          end do !*mj*
          end do !*j*

          end do !*mi*
          end do !*i*

          i_qlm = i_qlm + 1
        end do !*m*
        end do !*l*

      end do !*ii*



      !**** deallocate stack memory ****
      ok =        BA_pop_stack(t_mult(2))
      ok = ok.and.BA_pop_stack(gls(2))
      ok = ok.and.BA_pop_stack(gl(2))
      ok = ok.and.BA_pop_stack(tmp1(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit('paw_mult_coeff_set: error popping stack',0,1)

      return
      end


!      subroutine paw_pot_mult_print()
!      implicit none
!
!#include "bafdecls.fh"
!#include "paw_mult_data.fh"
!#include "paw_comp_charge_data.fh"
!#include "paw_comp_charge_matrix.fh"
!#include "paw_geom.fh"
!#include "paw_basis.fh"
!#include "paw_proj.fh"
!
!
!      !**** local variables ***
!      logical ok
!      integer ia,ii,ja,jj,kk,l,m,nion,npack0,mult_l,mabs
!      integer i,il,li,mi,mult_li,ill,jll
!      integer j,jl,lj,mj,mult_lj
!      integer  isgn,lmax,lmax2,indx
!      integer exi(2),tmp1(2),gls(2),gl(2),t_mult(2)
!      real*8  omega
!      complex*16 csum1,csum2
!
!      integer i_mtr,i_mtr0
!      integer nb,nb2
!      integer ilm
!      integer i_cp,i_cp0
!      integer nilm,njlm
!      integer i_qlm,i_qlm0
!      complex*16 tmp_mult_pot
!
!
!
!      write(48,*) paw_pot_mult(3)
!      do ii=1,ion_nion()
!        ia = ion_katm(ii)
!        nb = paw_basis_nbasis(ia)
!        nb2 = nb*nb
!        mult_l = paw_basis_mult_l(ia)
!        i_cp0 = int_mb(i_paw_pot_mult(1) + ii - 1)
!        nilm = 0
!        do i=1,nb
!        li = paw_basis_orb_l(i,ia)
!        njlm = 0
!        do j=1,nb
!        lj = paw_basis_orb_l(j,ia)
!          do mi=-li,li
!          do mj=-lj,lj
!             i_cp  = i_cp0-1+(njlm+lj+mj+1)+
!     >                (nilm+li+mi)*paw_proj_nbasis(ia)
!
!
!                  write(48,*) i,mi,j,mj,ii,dcpl_mb(paw_pot_mult(1)+i_cp)
!
!
!          end do !mi
!          end do !mj
!
!       njlm = njlm + 2*lj+1
!        end do !j
!        nilm = nilm + 2*li+1
!        end do !i
!      end do !ii
!
!
!      return
!      end



!**************************************************
!
!       Name: paw_mult_rcut
!
!       Purpose:
!
!       Created:        2/16/2003
!**************************************************
      function paw_mult_rcut()
      implicit none
      real*8 paw_mult_rcut !*RESULT*

#include "paw_mult_data.fh"

      paw_mult_rcut = sigma_smooth
      return
      end


!**************************************************
!
!       Name: paw_mult_ncut
!
!       Purpose:
!
!       Created:        2/16/2003
!**************************************************
      function paw_mult_ncut()
      implicit none
      integer paw_mult_ncut !*RESULT*

#include "paw_mult_data.fh"

      paw_mult_ncut = ncut
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

#include "bafdecls.fh"
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
      call dcopy((lmax+1)*nkatm,0.0d0,0,coeff,1)

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

#include "bafdecls.fh"
#include "paw_gaunt.fh"
#include "paw_geom.fh"
#include "paw_ovlp_data.fh"
#include "paw_mult_data.fh"
#include "paw_matrix_comp_pot.fh"
#include "paw_proj.fh"
#include "paw_basis.fh"


      !*** local variables **
      integer i,j,ii,ia,jj
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

          indx_2 = mtrx_ptr2
          do i=1,basis_nbasis
            li=paw_basis_orb_l(i,ia)
            do mi=-li,li

            do j=1,basis_nbasis
              lj=paw_basis_orb_l(j,ia)
              do mj=-lj,lj

              !*** check for non-zero gaunt coefficient ***
              if ( (l.le.(li+lj))   .and. 
     >             (l.ge.abs(li-lj)).and.
     >             (m.eq.(mi-mj)) ) then
                indx = mtrx_ptr + (j-1) 
     >                          + (i-1)*basis_nbasis
     >                          +     l*basis_nbasis**2
                dcpl_mb(comp_coeff(1)+jj)
     >            = dcpl_mb(comp_coeff(1)+jj)
     >            - gen_gaunt_coeff(l,m,li,mi,lj,mj)
     >             *dbl_mb(paw_comp_pot_matrix(1)+indx)
     >             *dcpl_mb(paw_ovlp_w(1)+indx_2)
              end if !*non-zero gaunt*

              indx_2 = indx_2 + 1
              end do !*mi*
            end do !*i*

            end do !*mj*
          end do !*j*

          jj = jj + 1
        end do !*m*
        end do !*l*

      end do !*ii*

      return
      end 



      subroutine paw_mult_pot_ptr(ptr)
      implicit none
      integer ptr

#include "paw_mult_data.fh"

      ptr = paw_pot_mult(1)
      return
      end





!**************************************************
!
!       Name: paw_mult_dn_cmp_smooth_spin_get
!
!       Purpose: returns  dn_cmp_smooth for each spin w/o core charge
!
!       Created:        1/11/2006
!**************************************************
      subroutine paw_mult_dn_cmp_smooth_spin_get(ms,dn_cmp_smooth)
      implicit none
      integer    ms
      complex*16 dn_cmp_smooth(*)

#include "bafdecls.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"

      !**** local variables ***
      logical ok
      integer ia,ii,jj,kk,l,m,nion,npack0,mult_l
      integer exi(2),tmp(2),QYlm(2)

      !**** allocate stack memory ****
      call Pack_npack(0,npack0)
      ok = BA_push_get(mt_dcpl,npack0,'exi',exi(2),exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,npack0,'QYlm',QYlm(2),QYlm(1))
      tmp(1) = exi(1)
      tmp(2) = exi(2)
      if (.not.ok)
     > call errquit(
     >  'paw_mult_dn_cmp_smooth_get2: out of stack memory',0,0)


      call dcopy(2*npack0,0.0d0,0,dn_cmp_smooth,1)
      nion = ion_nion()
      do ii=1,nion
         ia = ion_katm(ii)

         !**** Define QYlm ****
         mult_l = paw_basis_mult_l(ia)
         jj     = int_mb(i_paw_qlm(1)+ii-1)
         
         kk     = 0
         call dcopy(2*npack0,0.0d0,0,dcpl_mb(QYlm(1)),1)
         do l=0,mult_l
         do m=-l,l
           call Pack_cc_zaxpy(0,
     >               dcpl_mb(paw_qlm_spin(1,ms)+jj),
     >               dcpl_mb(g_lm(1)+(kk)*npack0),
     >               dcpl_mb(QYlm(1)))

           !write(*,*) " qlm=",l,m,dcpl_mb(paw_qlm_spin(1,ms)+jj)
           jj = jj + 1
           kk = kk + 1
         end do
         end do

         !**** Multiply by Structure Factor ****
         call strfac_pack(0,ii,dcpl_mb(exi(1)))
c         call Pack_cc_Mul(0,
c     >                    dcpl_mb(exi(1)),
c     >                    dcpl_mb(QYlm(1)),
c     >                    dcpl_mb(QYlm(1)))
         call Pack_cc_Mul2(0,
     >                    dcpl_mb(exi(1)),
     >                    dcpl_mb(QYlm(1)))


         !**** add up ncmp_smooth^ii  ****
         call Pack_tc_Mul(0,
     >                    dbl_mb(gk_smooth(1)),
     >                    dcpl_mb(QYlm(1)),
     >                    dcpl_mb(tmp(1)))
c         call Pack_cc_Sum(0,
c     >                    dcpl_mb(tmp(1)),
c     >                    dn_cmp_smooth,
c     >                    dn_cmp_smooth)
         call Pack_cc_Sum2(0,
     >                    dcpl_mb(tmp(1)),
     >                    dn_cmp_smooth)


      end do !*ii*


      !**** deallocate stack memory ****
      ok =        BA_pop_stack(QYlm(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not.ok)
     >  call errquit(
     >  'paw_mult_dn_cmp_smooth_get2:error popping stack',0,1)

      return
      end

