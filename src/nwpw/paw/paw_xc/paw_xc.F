!$Id$
!
!$Log: not supported by cvs2svn $
!Revision 1.20  2006/02/11 02:50:49  bylaska
!GGA's using 1st derivative formulas have been added in core part of PAW....EJB
!
!Revision 1.19  2005/01/30 07:16:18  bylaska
!Parallelized radial integrations...EJB
!
!Revision 1.18  2004/07/29 15:55:37  bylaska
!Temporary array for Gaunt coefficients added.  Speeds up the program considerably, but it is extrememely memory intensive.  Also added timing routines to multipole calculations and fixed the initial total charge calculation.
!
! ...EJB
!
!Revision 1.17  2003/10/28 01:52:30  edo
!errquizzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
!Revision 1.16  2003/10/21 02:05:18  marat
!switched to new errquit by running global replace operation
!see the script below (note it will not work on multiline errquit calls)
!*********************************************************
!#!/bin/sh
!
!e=`find . -name "*F" -print`
!
!for f in $e
!do
!cp $f $f.bak
!sed  's|\(^[ ].*call[ ]*errquit([^,]*\)\(,[^,]*\)\()\)|\1,0\2\3|' $f.bak > $f
!#rm $f.bak
!done
!**********************************************************
!
!Revision 1.15  2003/03/07 20:51:10  bylaska
!Code cleanup...0.0 changed to 0.0d0 in paw_xc.F
!Tangent vector now used for SD with Gram-schmidt.
!....EJB
!
!Revision 1.14  2003/03/05 23:16:32  bylaska
!Commented out write statements and other minor fixes.....
!self-consistent loop looks like it is working.....
!....EJB
!
!Revision 1.13  2003/03/04 00:04:04  marat
!added printouts for atomic potenitials
!for debug purposes
!MV
!
!Revision 1.12  2003/02/27 02:19:26  bylaska
!The electron gradient has been added.
!The eigenvalues are not yet agreeing with the F90 code....EJB
!
!Revision 1.11  2003/02/24 23:35:04  marat
!removed write statements
!MV
!
!Revision 1.10  2003/02/24 01:13:55  bylaska
!Memory deallocation completed...EJB
!
!Revision 1.9  2003/02/21 21:09:38  marat
!...
!MV
!
!Revision 1.8  2003/02/20 18:19:26  edo
!removed f90-isms
!
!Revision 1.7  2003/02/19 01:10:46  marat
!fixed vosko
!
!MV
!
!Revision 1.6  2003/02/19 01:02:04  marat
!finished generation of exch-corr
!tested only with l=0 with exchange only
!
!MV
!
!Revision 1.5  2003/02/18 01:57:31  marat
!generated atomic exc energies.
!there is a difference between paw and nwchem
!due different results from vosko exch-corr routine
!
!MV
!
!Revision 1.4  2003/02/17 22:03:36  marat
!almost done
!MV
!
!Revision 1.3  2003/02/17 05:54:59  marat
!added i_paw_vxc to common block
!MV
!
!Revision 1.2  2003/02/17 05:43:37  marat
!added memory deallocation
!MV
!
!Revision 1.1  2003/02/17 05:28:19  marat
!initial revision
!MV
!

**********************************
*
*
*
*
*
**********************************
      subroutine init_paw_xc(ispin,lmax)
      implicit none

      integer ispin
      integer lmax

#include "paw_xc_work_data.fh"
#include "paw_xc_data.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "bafdecls.fh"      
#include "paw_ma.fh"
#include "paw_geom.fh"
#include "paw_spher_func.fh"

c     **** local variables ****
      logical ok
      integer in,ia
      integer nion
      integer ic
      integer l,m
      integer i_p,i_t
      integer i_ylm,i_plm,i_rlm
      integer mtr_size
      integer nb
      double precision tmp_theta,cs_theta
      double precision angle_phi
      integer paw_vxc_size

c     **** external functions ****
      integer control_gga
      external control_gga
c    !**** external functions ****
      double precision theta_lm_div
      external         theta_lm_div



      call nwpw_timing_start(4)
      ok = .true.

      nion = ion_nion()
*     *** spin treatment flag (LDA=1,LSDA=2) ***
      paw_xc_ispin = ispin      

*     *** gga calculation ***
      paw_xc_gga = control_gga()

*     *** lmax for multipole expansion of the density***      
      paw_xc_lmax = lmax
      
*     ***paw_xc energies ***
      ok = ok .and. my_alloc(mt_dbl,nion,"paw_xc_e",
     >                                paw_xc_e)

     
      
*     *** xc matrix arrays ***
      mtr_size = 0
      do in=1,nion
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        if (mtr_size.lt.(nb**2)) mtr_size = (nb**2)
      end do
      mtr_size = mtr_size*(paw_xc_lmax+1)**2

      ok = ok .and. my_alloc(mt_dcpl,ispin*mtr_size,"paw_xc_matr_u",
     >                               paw_xc_matr)

      if (paw_xc_gga.ge.10) then
      ok = ok .and. my_alloc(mt_dcpl,3*ispin*mtr_size,"paw_xc_dmatr_u",
     >                               paw_xc_dmatr)
      end if
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_matr",
     >   mtr_size,0)
 
*     *** xc potential arrays ***
      ok = ok .and. my_alloc(mt_int,nion,"i_paw_xc_pot",
     >                                i_paw_xc_pot)
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating i_paw_xc_pot",0,0)

      mtr_size = 0
      do in=1,nion
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        int_mb(i_paw_xc_pot(1)+in-1)=mtr_size
        mtr_size=mtr_size+paw_proj_nbasis(ia)**2
      end do
      paw_xc_pot_size = mtr_size
      ok = ok .and. my_alloc(mt_dcpl,ispin*mtr_size,"paw_xc_pot_u",
     >                               paw_xc_pot) 
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_pot",
     >   mtr_size,0)



*     *** spherical grid arrays ***
      if(paw_xc_lmax .eq. 0 ) then
        !paw_xc_nphi   = 2
        paw_xc_nphi   = 1
        paw_xc_ntheta = 1
      else
        paw_xc_nphi   = 3*paw_xc_lmax
        paw_xc_ntheta = 3*paw_xc_lmax                      
      end if
      

      ok = ok .and. my_alloc(mt_dbl,paw_xc_nphi,"paw_xc_angle_phi",
     >                               paw_xc_angle_phi)

      ok = ok .and. my_alloc(mt_dbl,paw_xc_ntheta,"paw_xc_cos_theta",
     >                               paw_xc_cos_theta)

      ok = ok .and. my_alloc(mt_dbl,paw_xc_nphi,"paw_xc_w_phi",
     >                               paw_xc_w_phi)

      ok = ok .and. my_alloc(mt_dbl,paw_xc_ntheta,"paw_xc_w_theta",
     >                               paw_xc_w_theta)

      ok = ok .and. my_alloc(mt_dcpl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_ylm",
     >              paw_xc_ylm)


*     **** used for generating derivatives of ylm's ****
      if (paw_xc_gga.ge.10) then

c     **** derivatives wrt to theta ****
      ok = ok .and. my_alloc(mt_dcpl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_dylm_theta",
     >              paw_xc_dylm_theta)

c     **** derivatives wrt to phi ****
      ok = ok .and. my_alloc(mt_dcpl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_dylm_phi",
     >              paw_xc_dylm_phi)

      end if


      call get_spher_grid(paw_xc_ntheta,paw_xc_nphi,
     >                    dbl_mb(paw_xc_angle_phi(1)),
     >                    dbl_mb(paw_xc_cos_theta(1)),
     >                    dbl_mb(paw_xc_w_theta(1)),
     >                    dbl_mb(paw_xc_w_phi(1)))

c     **** define ylm's ****
      i_ylm = 0
      do i_t=1,paw_xc_ntheta
      do i_p=1,paw_xc_nphi
         do l=0,paw_xc_lmax
         do m=-l,l
           tmp_theta = theta_lm(l,m,
     >              dbl_mb(paw_xc_cos_theta(1)+i_t-1))
           

           angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
           dcpl_mb(paw_xc_ylm(1)+i_ylm) = 
     >             tmp_theta*exp(dcmplx(0.0d0,m*angle_phi))


           i_ylm = i_ylm + 1                          
                                        
         end do
         end do
      end do
      end do

      if (paw_xc_gga.ge.10) then

c       **** define derivative wrt to theta ****
        i_ylm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,paw_xc_lmax
           do m=-l,l
              cs_theta = dbl_mb(paw_xc_cos_theta(1)+i_t-1)
              angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
              dcpl_mb(paw_xc_dylm_theta(1)+i_ylm) 
     >        = dtheta_lm(l,m,cs_theta)
     >         *exp(dcmplx(0.0d0,m*angle_phi))
              i_ylm = i_ylm + 1 
           end do
           end do
        end do
        end do

c       **** define derivative wrt to phi ****
        i_ylm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,paw_xc_lmax
           do m=-l,l
              cs_theta = dbl_mb(paw_xc_cos_theta(1)+i_t-1)
              angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
              if (m.eq.0) then 
              dcpl_mb(paw_xc_dylm_phi(1)+ i_ylm) = dcmplx(0.0d0,0.0d0)
              else
              dcpl_mb(paw_xc_dylm_phi(1)+ i_ylm) 
     >        = theta_lm_div(l,m,cs_theta)
     >         * dcmplx(-m*dsin(m*angle_phi),
     >                         m*dcos(m*angle_phi))
              end if
              i_ylm = i_ylm + 1
           end do
           end do
        end do
        end do

      end if

*     *** temp arrays ***
      paw_xc_work_size = 0
      do in=1,nion
        ia = ion_katm(in)
c        ic = paw_basis_i_cut(ia)
        ic = paw_basis_npoints(ia)
        if(ic.gt.paw_xc_work_size) then
           paw_xc_work_size = ic
        end if
      end do
 
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_rho_ae",
     >                               rho_ae)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_rho_ps",
     >                               rho_ps)
      
*     **** allocate gradient's and agr's ****
      if (paw_xc_gga.ge.10) then
       ok = ok.and.my_alloc(mt_dbl,ispin*3*paw_xc_work_size,
     >                               "paw_tmp_drho_ae",
     >                               rho_ae_prime)
       ok = ok.and.my_alloc(mt_dbl,ispin*3*paw_xc_work_size,
     >                               "paw_tmp_drho_ps",
     >                               rho_ps_prime)
       ok = ok.and.my_alloc(mt_dbl,(2*ispin-1)*paw_xc_work_size,
     >                               "paw_tmp_agr_ae",
     >                               agr_ae)
       ok = ok.and.my_alloc(mt_dbl,(2*ispin-1)*paw_xc_work_size,
     >                               "paw_tmp_agr_ps",
     >                               agr_ps)
       ok = ok.and.my_alloc(mt_dbl,(2*ispin-1)*paw_xc_work_size,
     >                               "paw_tmp_fdn_ae",
     >                               fdn_ae)
       ok = ok.and.my_alloc(mt_dbl,(2*ispin-1)*paw_xc_work_size,
     >                               "paw_tmp_fdn_ps",
     >                               fdn_ps)
      end if

      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_vxc_ae",
     >                               vxc_ae)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_vxc_ps",
     >                               vxc_ps)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_exc_ae",
     >                               exc_ae)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_exc_ps",
     >                               exc_ps)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size,"paw_tmp_x",
     >                               xc_temp)

      ok = ok .and. my_alloc(mt_dcpl,paw_xc_work_size,
     >                      "paw_tmp_x_cmplx",
     >                       xc_c_temp)


*     *** allocate vxclm multipole expansion  arrays ****
      paw_vxc_size = paw_xc_work_size*(paw_xc_lmax+1)**2

      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,"paw_vxc_u_ae",
     >                               paw_vxc_ae)
      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,"paw_vxc_u_ps",
     >                               paw_vxc_ps)

*     *** allocate dvxclm multipole expansion  arrays ****
      if (paw_xc_gga.ge.10) then
      ok = ok .and. my_alloc(mt_dcpl,3*ispin*paw_vxc_size,
     >                               "paw_dvxc_u_ae",
     >                               paw_dvxc_ae)
      ok = ok .and. my_alloc(mt_dcpl,3*ispin*paw_vxc_size,
     >                               "paw_dvxc_u_ps",
     >                               paw_dvxc_ps)
      end if

*     *** allocate rholm multipole expansion  arrays ****
      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,"paw_rho2_ae",
     >                               paw_rho2_ae)
      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,"paw_rho2_ps",
     >                               paw_rho2_ps)

*     *** allocate rholm_prime multipole expansion arrays - need for GGA's****
      if (paw_xc_gga.ge.10) then
      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,
     >                               "paw_rho2_ae_prime",
     >                               paw_rho2_ae_prime)
      ok = ok .and. my_alloc(mt_dcpl,ispin*paw_vxc_size,
     >                               "paw_rho2_ps_prime",
     >                               paw_rho2_ps_prime)
      end if

      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating work arrays",0,0)
      

      call paw_init_gntxc(paw_xc_lmax)
      if ((paw_xc_gga.ge.10).and.(paw_xc_lmax.ge.1)) then
        call paw_init_gntxc2(paw_xc_lmax)
        call paw_init_gntxc3(paw_xc_lmax)
      end if

      call nwpw_timing_end(4)
      return
      end 

**********************************
*
*
*
*
*
**********************************
      subroutine paw_xc_solve()
      implicit none

#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"

#include "paw_xc_work_data.fh"
#include "paw_xc_data.fh"
#include "integrate.fh"
#include "bafdecls.fh"

c      double precision gen_gaunt_coeff,gen_gaunt2_coeff,gen_gaunt3_coeff
c      external         gen_gaunt_coeff,gen_gaunt2_coeff,gen_gaunt3_coeff
      integer i,j,np
      integer in,ia
      integer i_t,i_p
      integer nion 
      integer l,m,ms
      integer ng,ic,ig,igg,istart
      integer i_rho0,i_rho
      integer n_s
      integer i_ylm,i_ylm1
      integer nlmi
      integer i_core_ps,i_core_ae
      integer i_core_ps_prime,i_core_ae_prime
      integer i_vxc0,i_vxc
      integer i_pot0,i_pot
      integer i_r
      integer nb,nb2
      integer i_phi_ae0,i_phi_ps0
      integer i_phi_ae0_prime,i_phi_ps0_prime
      integer i_mtr0,i_mtr
      integer li,lj,mi,mj
      integer nilm,njlm
      integer lmax2
      double precision exc_tmp,gnt1
      double precision log_amesh
      double precision tmp_ae,tmp_ps
      complex*16 tmp
      complex*16 ylm

      call nwpw_timing_start(4)
      call nwpw_timing_start(22)
      nion = ion_nion() 
      call Parallel_np(np)


*     *** zero out vxc arrays ***
      call dcopy(2*paw_xc_pot(3),0.0d0,0,dcpl_mb(paw_xc_pot(1)),1)

c      do i=1,paw_vxc_u_ae(3)
c       dcpl_mb(paw_vxc_u_ae(1)+i-1)=dcmplx(0.0d0,0.0d0)  
c       dcpl_mb(paw_vxc_d_ae(1)+i-1)=dcmplx(0.0d0,0.0d0)
c       dcpl_mb(paw_vxc_u_ps(1)+i-1)=dcmplx(0.0d0,0.0d0) 
c       dcpl_mb(paw_vxc_d_ps(1)+i-1)=dcmplx(0.0d0,0.0d0)
c      end do
c      do i=1,paw_xc_pot_u(3)
c       dcpl_mb(paw_xc_pot_u(1)+i-1)=dcmplx(0.0d0,0.0d0)  
c       dcpl_mb(paw_xc_pot_d(1)+i-1)=dcmplx(0.0d0,0.0d0)  
c      end do


*     *** index for spin down in temp arrays ***      
      n_s = paw_xc_work_size 

      lmax2  = (paw_xc_lmax+1)**2

     
      do in=1,nion

        ia = ion_katm(in)
        ng = paw_basis_ngrid(ia)
        ic = paw_basis_npoints(ia)
        istart = paw_basis_i_start(ia)
c        i_rho0 = int_mb(i_paw_rho(1) + in - 1)
        i_core_ae = paw_basis_i_core_density(ia)
        i_core_ps = paw_basis_i_ps_core_density(ia)
        i_core_ae_prime = paw_basis_i_core_density_prime(ia)
        i_core_ps_prime = paw_basis_i_ps_core_density_prime(ia)
        i_r = paw_basis_i_rgrid(ia)
        log_amesh = paw_basis_log_amesh(ia)

        dbl_mb(paw_xc_e(1)+in-1)=0.0d0
        i_ylm = 0
        i_ylm1 = 0
c

*     *** zero out vxclm multipole arrays ***
      call dcopy(2*paw_vxc_ae(3),0.0d0,0,dcpl_mb(paw_vxc_ae(1)),1)
      call dcopy(2*paw_vxc_ps(3),0.0d0,0,dcpl_mb(paw_vxc_ps(1)),1)

*     *** zero out dvxclm multipole arrays ***
      if (paw_xc_gga.ge.10) then
      call dcopy(2*paw_dvxc_ae(3),0.0d0,0,dcpl_mb(paw_dvxc_ae(1)),1)
      call dcopy(2*paw_dvxc_ps(3),0.0d0,0,dcpl_mb(paw_dvxc_ps(1)),1)
      end if

      do i_t = 1,paw_xc_ntheta
      do i_p = 1,paw_xc_nphi
      
*       *** zero out temp arrays ***
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(rho_ae(1)),1)
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(rho_ps(1)),1)
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(vxc_ae(1)),1)
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(vxc_ps(1)),1)
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(exc_ae(1)),1)
        call dcopy(2*paw_xc_work_size,0.0d0,0,dbl_mb(exc_ps(1)),1)

       
*       *** find rholm multipole expansion on spherical grid ***          
        call paw_density_solve2(in,paw_xc_lmax,
     >                          ic,lmax2,paw_xc_ispin,
     >                          dcpl_mb(paw_rho2_ae(1)),
     >                          dcpl_mb(paw_rho2_ps(1)))

*       *** generate atomic densities on spherical grid ***          
        call paw_xc_gen_atomic_densities(ic,lmax2,paw_xc_ispin,
     >                                    dcpl_mb(paw_rho2_ae(1)),
     >                                    dcpl_mb(paw_xc_ylm(1)+i_ylm),
     >                                    dbl_mb(i_core_ae+istart-1),
     >                                    dbl_mb(rho_ae(1)))
        call paw_xc_gen_atomic_densities(ic,lmax2,paw_xc_ispin,
     >                                    dcpl_mb(paw_rho2_ps(1)),
     >                                    dcpl_mb(paw_xc_ylm(1)+i_ylm),
     >                                    dbl_mb(i_core_ps+istart-1),
     >                                    dbl_mb(rho_ps(1)))


*       ***************************************************
*       *** find exchange-correlation on spherical grid ***
*       ***************************************************


*       **** LDA functionals on spherical grid ****
        if (paw_xc_gga.lt.10) then

*          **** LDA functionals of ae and ps atomic densities ****
           call paw_vosko(ic,ic,paw_xc_ispin,
     >                    dbl_mb(rho_ae(1)),        
     >                    dbl_mb(exc_ae(1)),        
     >                    dbl_mb(vxc_ae(1)),        
     >                    dbl_mb(xc_temp(1)) )       

           call paw_vosko(ic,ic,paw_xc_ispin,
     >                    dbl_mb(rho_ps(1)),        
     >                    dbl_mb(exc_ps(1)),        
     >                    dbl_mb(vxc_ps(1)),        
     >                    dbl_mb(xc_temp(1)) )       



*       **** GGA functionals on spherical grid ****
        else

*         *** find rholm_prime multipole expansion on spherical grid ***          
          call paw_density_prime_solve2(in,paw_xc_lmax,
     >                          ic,lmax2,paw_xc_ispin,
     >                          dcpl_mb(paw_rho2_ae_prime(1)),
     >                          dcpl_mb(paw_rho2_ps_prime(1)))

*       *** generate the gradient of atomic densities in spherical ****
*       **** coordinates  on spherical grid                        ****          
          call paw_xc_gen_atomic_gradients(ic,lmax2,paw_xc_ispin,
     >                 dcpl_mb(paw_rho2_ae(1)),
     >                 dcpl_mb(paw_rho2_ae_prime(1)),
     >                 dcpl_mb(paw_xc_ylm(1)+i_ylm),
     >                 dcpl_mb(paw_xc_dylm_theta(1)+i_ylm),
     >                 dcpl_mb(paw_xc_dylm_phi(1)+i_ylm),
     >                 dbl_mb(i_r+istart-1),
     >                 dbl_mb(i_core_ae+istart-1),
     >                 dbl_mb(i_core_ae_prime+istart-1),
     >                 dbl_mb(rho_ae_prime(1)))
          call paw_xc_gen_atomic_gradients(ic,lmax2,paw_xc_ispin,
     >                 dcpl_mb(paw_rho2_ps(1)),
     >                 dcpl_mb(paw_rho2_ps_prime(1)),
     >                 dcpl_mb(paw_xc_ylm(1)+i_ylm),
     >                 dcpl_mb(paw_xc_dylm_theta(1)+i_ylm),
     >                 dcpl_mb(paw_xc_dylm_phi(1)+i_ylm),
     >                 dbl_mb(i_r+istart-1),
     >                 dbl_mb(i_core_ps+istart-1),
     >                 dbl_mb(i_core_ps_prime+istart-1),
     >                 dbl_mb(rho_ps_prime(1)))


*         **** generate agr ****
          call paw_xc_gen_atomic_agr(ic,lmax2,paw_xc_ispin,
     >                               dbl_mb(rho_ae_prime(1)),
     >                               dbl_mb(agr_ae(1)))
          call paw_xc_gen_atomic_agr(ic,lmax2,paw_xc_ispin,
     >                               dbl_mb(rho_ps_prime(1)),
     >                               dbl_mb(agr_ps(1)))

*          **** GGA functionals of ae and ps atomic densities ****
           call paw_gga(paw_xc_gga,ic,paw_xc_ispin,
     >                    dbl_mb(rho_ae(1)),        
     >                    dbl_mb(agr_ae(1)),        
     >                    dbl_mb(exc_ae(1)),
     >                    dbl_mb(vxc_ae(1)),    !* df/dnup, df/dndn
     >                    dbl_mb(fdn_ae(1)),    !* df/d|grad nup|,df/d|grad ndn|, df/d|grad n|
     >                    dbl_mb(xc_temp(1)))       

           call paw_gga(paw_xc_gga,ic,paw_xc_ispin,
     >                    dbl_mb(rho_ps(1)),        
     >                    dbl_mb(agr_ps(1)),        
     >                    dbl_mb(exc_ps(1)),
     >                    dbl_mb(vxc_ps(1)),   !* df/dnup, df/dndn
     >                    dbl_mb(fdn_ps(1)),   !* df/d|grad nup|,df/d|grad ndn|, df/d|grad n|
     >                    dbl_mb(xc_temp(1)))       


c          **** if restricted calculate df/d|grad n| *(grad n)/|grad n|  ****
*          *****   and put it in rho_prime                               ****

*          **** if unrestricted calculate (df/d|grad nup|* (grad nup)/|grad nup|)  ****
*          ****                         + (df/d|grad n|  * (grad n)/|grad n|)      ****
*          ****             and calculate (df/d|grad ndn|* (grad ndn)/|grad ndn|)  ****
*          ****                         + (df/d|grad n|  * (grad n)/|grad n|)      ****
*          *****   and put it in rho_prime                                         ****
           call paw_xc_gen_dvxc(ic,lmax2,paw_xc_ispin,
     >                          dbl_mb(fdn_ae(1)),
     >                          dbl_mb(agr_ae(1)),
     >                          dbl_mb(rho_ae_prime(1)))
           call paw_xc_gen_dvxc(ic,lmax2,paw_xc_ispin,
     >                          dbl_mb(fdn_ps(1)),
     >                          dbl_mb(agr_ps(1)),
     >                          dbl_mb(rho_ps_prime(1)))

c          **** find dvxclm multipole expansion ****
           call paw_xc_addto_dvxclm(ic,lmax2,paw_xc_ispin,
     >                          (dbl_mb(paw_xc_w_theta(1)+i_t-1)
     >                          *dbl_mb(paw_xc_w_phi(1)  +i_p-1)),
     >                          dcpl_mb(paw_xc_ylm(1)+i_ylm),
     >                          dbl_mb(rho_ae_prime(1)),
     >                          dbl_mb(rho_ps_prime(1)),
     >                          dcpl_mb(paw_dvxc_ae(1)),
     >                          dcpl_mb(paw_dvxc_ps(1)))

        end if
        i_ylm = i_ylm + lmax2

*       ************************************************************************************
*       *** find vxclm the multipole expansion of atomic vxc=df/dn  or (df/dnup,df/dndn) ***          
*       ************************************************************************************
        call paw_xc_addto_vxclm(ic,lmax2,paw_xc_ispin,
     >                          (dbl_mb(paw_xc_w_theta(1)+i_t-1)
     >                          *dbl_mb(paw_xc_w_phi(1)  +i_p-1)),
     >                          dcpl_mb(paw_xc_ylm(1)+i_ylm1),
     >                          dbl_mb(vxc_ae(1)),
     >                          dbl_mb(vxc_ps(1)),
     >                          dcpl_mb(paw_vxc_ae(1)),
     >                          dcpl_mb(paw_vxc_ps(1)))
        i_ylm1 = i_ylm1 + lmax2





*        ******************************************************
*        *** compute the atomic exchange correlation energy ***
*        ******************************************************
         do ig=1,ic
          dbl_mb(xc_temp(1)+ig-1)=
     >             (
     >             dbl_mb(rho_ae(1)+ig-1)  + 
     >             dbl_mb(rho_ae(1)+(paw_xc_ispin-1)*ic+ig-1)
     >             )*dbl_mb(exc_ae(1)+ig-1)
     >             -
     >             (
     >             dbl_mb(rho_ps(1)+ig-1)+ 
     >             dbl_mb(rho_ps(1)+(paw_xc_ispin-1)*ic+ig-1)
     >             )*
     >             dbl_mb(exc_ps(1)+ig-1)
        end do
        exc_tmp = def_integr1(0,dbl_mb(xc_temp(1)),2,
     >                       dbl_mb(i_r+istart-1),
     >                       log_amesh,ic,istart)

        dbl_mb(paw_xc_e(1)+in-1)=dbl_mb(paw_xc_e(1)+in-1)+
     >                           exc_tmp*
     >                           dbl_mb(paw_xc_w_theta(1)+i_t-1)*
     >                           dbl_mb(paw_xc_w_phi(1)+i_p-1)



      end do !i_phi
      end do !i_theta


*       *********************************************************
*       ****    non-local operator computation - LDA part    ****
*       *********************************************************
        nb  = paw_basis_nbasis(ia)
        nb2 = nb*nb
        i_phi_ae0 = paw_basis_i_phi_ae(ia)
        i_phi_ps0 = paw_basis_i_phi_ps(ia)
 
*       **** compute (vxc^a)_nln'l'^lm radial integrals *****
        call paw_xc_gen_matr(ng,nb,ic,istart,lmax2,paw_xc_ispin,
     >                       log_amesh,
     >                       dbl_mb(i_phi_ae0),
     >                       dbl_mb(i_phi_ps0),
     >                       dbl_mb(i_r+istart-1),
     >                       dcpl_mb(paw_vxc_ae(1)),
     >                       dcpl_mb(paw_vxc_ps(1)),
     >                       dcpl_mb(xc_c_temp(1)),
     >                       dcpl_mb(paw_xc_matr(1)))


*       **** xc potential  non-local matrix elements ****
        i_pot0 = int_mb(i_paw_xc_pot(1) + in - 1)
        do ms=1,paw_xc_ispin
           call paw_addto_potential_gntxc(ia,
     >           dcpl_mb(paw_xc_matr(1)+(ms-1)*nb2*lmax2),
     >           dcpl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
        end do

c        nlmi = 0
c        do l=0,paw_xc_lmax
c        nilm = 0
c        do i=1,nb
c        li = paw_basis_orb_l(i,ia)
c        njlm = 0
c        do j=1,nb
c        lj = paw_basis_orb_l(j,ia)
c        if( ((l.le.(li+lj)) .and. (l.ge.ABS(li-lj))).or.
c     >      (mod(l+li+lj,2).eq.0)                        ) then
c
c        do m=-l,l
c          do mi=-li,li
c          do mj=-lj,lj
c          if(m.eq.(mi-mj)) then
c            gnt1 = gen_gaunt_coeff(l,m,li,mi,lj,mj)       
c
c            do ms=1,paw_xc_ispin
cc            i_mtr = i_mtr0-1+j+(i-1)*nb+(nlmi+l+m)*nb2+(ms-1)*nb2*lmax2
c            i_mtr = (j-1)+(i-1)*nb+(nlmi+l+m)*nb2+(ms-1)*nb2*lmax2
c            i_pot = i_pot0
c     >            + (njlm+lj+mj)
c     >            + (nilm+li+mi)*paw_proj_nbasis(ia)
c     >            + (ms-1)*paw_xc_pot_size
c              dcpl_mb(paw_xc_pot(1)+i_pot)
c     >           = dcpl_mb(paw_xc_pot(1)+i_pot) 
c     >           +  gnt1*dcpl_mb(paw_xc_matr(1)+i_mtr)          
c            end do
c
c          end if
c          end do !mi
c          end do !mj
c        end do !m
c
c        end if
c
c        njlm = njlm + 2*lj+1
c        end do !j
c        nilm = nilm + 2*li+1
c        end do !i  
c        nlmi=nlmi+2*l+1
c        end do !l


*       *********************************************************
*       ****    non-local operator computation - GGA part    ****
*       *********************************************************
        if (paw_xc_gga.ge.10) then

          i_phi_ae0_prime = paw_basis_i_phi_ae_prime(ia)
          i_phi_ps0_prime = paw_basis_i_phi_ps_prime(ia)

*         **** compute (dvxc^a)_nln'l'^lm radial integrals *****
          call paw_xc_gen_dmatr(ng,nb,ic,istart,lmax2,paw_xc_ispin,
     >                       log_amesh,
     >                       dbl_mb(i_phi_ae0),
     >                       dbl_mb(i_phi_ps0),
     >                       dbl_mb(i_phi_ae0_prime),
     >                       dbl_mb(i_phi_ps0_prime),
     >                       dbl_mb(i_r+istart-1),
     >                       dcpl_mb(paw_dvxc_ae(1)),
     >                       dcpl_mb(paw_dvxc_ps(1)),
     >                       dcpl_mb(xc_c_temp(1)),
     >                       dcpl_mb(paw_xc_dmatr(1)))


*         **** xc potential  non-local matrix elements ****
          !i_pot0 = int_mb(i_paw_xc_pot(1) + in - 1)
          do ms=1,paw_xc_ispin
           call paw_addto_potential_gntxc(ia,
     >           dcpl_mb(paw_xc_dmatr(1)+(ms-1)*nb2*lmax2*3),
     >           dcpl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
          end do

c         **** theta and phi parts ****
          if (lmax2.gt.1) then

c           **** theta derivatives ****
            do ms=1,paw_xc_ispin
             call paw_addto_potential_gntxc2(ia,
     >         dcpl_mb(paw_xc_dmatr(1)+nb2*lmax2+(ms-1)*nb2*lmax2*3),
     >         dcpl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
            end do

c           **** phi derivatives ****
            do ms=1,paw_xc_ispin
             call paw_addto_potential_gntxc3(ia,
     >         dcpl_mb(paw_xc_dmatr(1)+nb2*lmax2*2+(ms-1)*nb2*lmax2*3),
     >         dcpl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
            end do


         end if !**lmax2>1**

        
        end if !**gga**
 
      end do !in

      if (np.gt.1)  then

c       **** xc non-local matrix elements ****
        call D3dB_Vector_Sumall(2*paw_xc_pot(3),
     >                          dcpl_mb(paw_xc_pot(1)))

c       **** atomic exchange-correlation energies ****
        call D3dB_Vector_Sumall(paw_xc_e(3),dbl_mb(paw_xc_e(1)))
      end if
     

      call nwpw_timing_end(4)
      call nwpw_timing_end(22)
      return
      end 


**********************************
*
*
*
*
*
**********************************
      subroutine paw_xc_end()
      implicit none

#include "paw_xc_work_data.fh"
#include "paw_xc_data.fh"
#include "paw_ma.fh"

      logical ok
      integer ms,i
      
      call nwpw_timing_start(4)
      ok = .true.

      ok = ok .and. my_dealloc(exc_ps)
      ok = ok .and. my_dealloc(exc_ae)
      ok = ok .and. my_dealloc(vxc_ps)
      ok = ok .and. my_dealloc(vxc_ae)
      ok = ok .and. my_dealloc(rho_ps)
      ok = ok .and. my_dealloc(rho_ae)
      ok = ok .and. my_dealloc(xc_temp)

      ok = ok .and. my_dealloc(i_paw_xc_pot)
c      ok = ok .and. my_dealloc(i_paw_xc_matr)
      ok = ok .and. my_dealloc(paw_xc_e)
      ok = ok. and. my_dealloc(paw_xc_ylm) 
      ok = ok. and. my_dealloc(paw_xc_w_theta) 
      ok = ok. and. my_dealloc(paw_xc_w_phi) 
      ok = ok. and. my_dealloc(paw_xc_cos_theta) 
      ok = ok. and. my_dealloc(paw_xc_angle_phi) 
c      ok = ok .and. my_dealloc(i_paw_vxc)

      ok = ok .and. my_dealloc(xc_c_temp)
      ok = ok .and. my_dealloc(paw_xc_pot)
      ok = ok .and. my_dealloc(paw_xc_matr)
      ok = ok .and. my_dealloc(paw_vxc_ae)
      ok = ok. and. my_dealloc(paw_vxc_ps)
      ok = ok. and. my_dealloc(paw_rho2_ae)
      ok = ok. and. my_dealloc(paw_rho2_ps)

      if (paw_xc_gga.ge.10) then
         ok = ok. and. my_dealloc(paw_rho2_ae_prime)
         ok = ok. and. my_dealloc(paw_rho2_ps_prime)
         ok = ok .and. my_dealloc(rho_ps_prime)
         ok = ok .and. my_dealloc(rho_ae_prime)
         ok = ok .and. my_dealloc(agr_ae)
         ok = ok .and. my_dealloc(agr_ps)
         ok = ok .and. my_dealloc(fdn_ae)
         ok = ok .and. my_dealloc(fdn_ps)
         ok = ok .and. my_dealloc(paw_xc_dylm_theta)
         ok = ok .and. my_dealloc(paw_xc_dylm_phi)
         ok = ok .and. my_dealloc(paw_dvxc_ae)
         ok = ok. and. my_dealloc(paw_dvxc_ps)
         ok = ok .and. my_dealloc(paw_xc_dmatr)
      end if
      if (.not.ok)
     > call errquit("paw_vxc_end: error freeing heap",0,0)


      call paw_end_gntxc()
      if ((paw_xc_gga.ge.10).and.(paw_xc_lmax.ge.1)) then
        call paw_end_gntxc2()
        call paw_end_gntxc3()
      end if

      call nwpw_timing_end(4)
      return
      end 


**********************************
*
*
*
*
*
**********************************
      function paw_energy_xc_atom()
      implicit none

      double precision paw_energy_xc_atom

#include "paw_xc_data.fh"
#include "bafdecls.fh"

      integer in

      paw_energy_xc_atom = 0
      do in=1,paw_xc_e(3)
       paw_energy_xc_atom = paw_energy_xc_atom + 
     >                      dbl_mb(paw_xc_e(1)+in-1)  
      end do

      return
      end



      subroutine paw_xc_pot_up_ptr(ptr)
      implicit none
      integer ptr

#include "paw_xc_data.fh"

      ptr = paw_xc_pot(1)
      return
      end

      subroutine paw_xc_pot_down_ptr(ptr)
      implicit none
      integer ptr

#include "paw_xc_data.fh"

      ptr = paw_xc_pot(1) + paw_xc_pot_size
      return
      end


      subroutine paw_xc_gen_sphere_rho(lmax2,ic,rholm,Ylm,rho)
      implicit none
      integer lmax2,ic
      real*8     rholm(*)
      real*8     Ylm(*)
      real*8     rho(*)

      integer    l,indx1,igg1,ig

      indx1 = 1
      igg1  = 1
      do l=1,lmax2
      do ig=1,ic
        rho(ig) = rho(ig) + rholm(igg1)*Ylm(indx1) 
     >                    - rholm(igg1+1)*Ylm(indx1+1) 
        igg1 = igg1 + 2
      end do
      indx1 = indx1 + 2
      end do
      
      return
      end




      subroutine paw_xc_rho_div_r(ic,r,rho)
      implicit none
      integer ic
      real*8     r(*)
      real*8     rho(*)

      integer    ig

      do ig=1,ic
        rho(ig) = rho(ig)/r(ig)
      end do

      return
      end



cc     *************************************************
cc     *                                               *
cc     *            derivative_ng_ic                  *
cc     *                                               *
cc     *************************************************
cc
cc  This routine computes the seven point derivative of f.
cc  where f and df are stored on a logarithmic grid. The
cc  dimensions of f and df are, f(1:ng), and df(1:ic) and the 
cc  indexing is such that df  distributed in parallel while
cc  the complete f is stored on each processor,
cc  i.e. for the operation df <-- f
cc    df(i) = f(i+istart-1) or df(i-istart+1) = f(i)
cc
cc     Note- this routine still needs to be vectorized.
c
c      subroutine derivative_ng_ic(ng,istart,ic,log_amesh,r,f,df)
c      implicit none
c      integer           ng,istart,ic
c      double precision  log_amesh
c      double precision  r(ng)
c      double precision  f(ng)
c      double precision df(ic)
c
c      double precision one_over_60
c      parameter (one_over_60 = 1.0d0/60.0d0)
c
c      integer i,shift,n1,n2,m1,m2
c      double precision aa
c
c      aa = one_over_60/log_amesh
c      n1 = istart-1+1
c      n2 = istart-1+ic
c      m1 = n1
c      m2 = n2
c      shift = -(istart-1)
c
c      if (n1.le.3) then
c        if ((n1.eq.1).and.(n1.ge.m1).and.(n1.le.m2)) then
c          df(1) = aa*(-147.0d0*f(1)
c     >               + 360.0d0*f(2)
c     >               - 450.0d0*f(3)
c     >               + 400.0d0*f(4)
c     >               - 225.0d0*f(5)
c     >               +  72.0d0*f(6)
c     >               -  10.0d0*f(7))/r(1)
c          n1 = n1+1
c        end if 
c        if ((n1.eq.2).and.(n1.ge.m1).and.(n1.le.m2)) then
c          df(2) = aa*( -10.0d0*f(1)
c     >               -  77.0d0*f(2)
c     >               + 150.0d0*f(3)
c     >               - 100.0d0*f(4)
c     >               +  50.0d0*f(5)
c     >               -  15.0d0*f(6)
c     >               +   2.0d0*f(7))/r(2)
c          n1 = n1+1
c        end if 
c        if ((n1.eq.3.and.(n1.ge.m1).and.(n1.le.m2))) then
c          df(3) = aa*(  +2.0d0*f(1)
c     >               -  24.0d0*f(2)
c     >               -  35.0d0*f(3)
c     >               +  80.0d0*f(4)
c     >               -  30.0d0*f(5)
c     >               +   8.0d0*f(6)
c     >               -   1.0d0*f(7))/r(3)
c          n1 = n1+1
c        end if 
c      end if
c
c      if (n2.ge.(ng-2)) then
c        if ((n2.eq.ng).and.(n2.ge.m1).and.(n2.le.m2)) then
c          df(ng) = aa*( +147.0d0*f(ng)
c     >                - 360.0d0*f(ng-1)
c     >                + 450.0d0*f(ng-2)
c     >                - 400.0d0*f(ng-3)
c     >                + 225.0d0*f(ng-4)
c     >                -  72.0d0*f(ng-5)
c     >                +  10.0d0*f(ng-6))/r(ng)
c          n2 = n2-1
c        end if
c        if ((n2.eq.(ng-1).and.(n2.ge.m1).and.(n2.le.m2))) then
c          df(ng-1) = aa*( +10.0d0*f(ng)
c     >                  +  77.0d0*f(ng-1)
c     >                  - 150.0d0*f(ng-2)
c     >                  + 100.0d0*f(ng-3)
c     >                  -  50.0d0*f(ng-4)
c     >                  +  15.0d0*f(ng-5)
c     >                  -   2.0d0*f(ng-6))/r(ng-1)
c          n2 = n2-1
c        end if
c        if ((n2.eq.(ng-2).and.(n2.ge.m1).and.(n2.le.m2))) then
c          df(ng-2) = aa*(  -2.0d0*f(ng)
c     >                  +  24.0d0*f(ng-1)
c     >                  +  35.0d0*f(ng-2)
c     >                  -  80.0d0*f(ng-3)
c     >                  +  30.0d0*f(ng-4)
c     >                  -   8.0d0*f(ng-5)
c     >                  +   1.0d0*f(ng-6))/r(ng-2)
c          n2 = n2-1
c        end if
c      end if
c
c      do i=n1,n2
c        df(shift+i) = aa*(  -1.0d0*f(i-3)
c     >                   +   9.0d0*f(i-2)
c     >                   -  45.0d0*f(i-1)
c     >                   +  45.0d0*f(i+1)
c     >                   -   9.0d0*f(i+2)
c     >                   +   1.0d0*f(i+3))/r(i)
c      end do
c
c      return
c      end
      



      subroutine paw_xc_addto_vxclm(ic,lmax2,ispin,
     >                              alpha,
     >                              ylm,
     >                              vxc_ae,vxc_ps,
     >                              vxclm_ae,vxclm_ps)
      implicit none
      integer ic,lmax2,ispin
      double precision alpha
      double precision ylm(2,*)
      double precision vxc_ae(ic,ispin)
      double precision vxc_ps(ic,ispin)

      double precision vxclm_ae(2,ic,lmax2,ispin)
      double precision vxclm_ps(2,ic,lmax2,ispin)

      integer    i,lm,ms

      do ms = 1,ispin
      do lm = 1,lmax2
      do i  = 1,ic
        vxclm_ae(1,i,lm,ms) = vxclm_ae(1,i,lm,ms) 
     >                      + vxc_ae(i,ms)*(ylm(1,lm)*alpha)
        vxclm_ae(2,i,lm,ms) = vxclm_ae(2,i,lm,ms) 
     >                      - vxc_ae(i,ms)*(ylm(2,lm)*alpha)
      end do
      end do
      end do

      do ms = 1,ispin
      do lm = 1,lmax2
      do i  = 1,ic
        vxclm_ps(1,i,lm,ms) = vxclm_ps(1,i,lm,ms) 
     >                      + vxc_ps(i,ms)*(ylm(1,lm)*alpha)
        vxclm_ps(2,i,lm,ms) = vxclm_ps(2,i,lm,ms) 
     >                      - vxc_ps(i,ms)*(ylm(2,lm)*alpha)
      end do
      end do
      end do

      return
      end

 
      subroutine paw_xc_gen_matr(ng,nb,ic,istart,lmax2,ispin,log_amesh,
     >                           phi_ae,phi_ps,r,
     >                           vxclm_ae,vxclm_ps,
     >                           tmpC,
     >                           matr) 
      implicit none
      integer ng,nb,ic,istart,lmax2,ispin
      double precision log_amesh

      double precision  phi_ae(ng,nb)
      double precision  phi_ps(ng,nb)
      double precision  r(ic)
      double precision  vxclm_ae(2,ic,lmax2,ispin)
      double precision  vxclm_ps(2,ic,lmax2,ispin)
      double precision  tmpC(2,ic)

      double complex   matr(nb,nb,lmax2,ispin)

#include "integrate.fh"

*     **** local varialbles ****
      integer ig,i,j,lm,ms
      double precision tmp_ae,tmp_ps


      do ms=1,ispin
      do lm=1,lmax2


        do i=1,nb
        do j=i,nb

          do ig=1,ic
             tmp_ae = phi_ae(ig+istart-1,i)
     >               *phi_ae(ig+istart-1,j)
     >               /(r(ig)**2)
             tmp_ps = phi_ps(ig+istart-1,i)
     >               *phi_ps(ig+istart-1,j)
     >               /(r(ig)**2)

             tmpC(1,ig) = vxclm_ae(1,ig,lm,ms)*tmp_ae
     >                     - vxclm_ps(1,ig,lm,ms)*tmp_ps
             tmpC(2,ig) = vxclm_ae(2,ig,lm,ms)*tmp_ae
     >                     - vxclm_ps(2,ig,lm,ms)*tmp_ps
          end do
          matr(i,j,lm,ms)= c_def_integr1(0,tmpC,2,r,log_amesh,ic,istart)
          if (i.ne.j) matr(j,i,lm,ms) = matr(i,j,lm,ms)
        end do
        end do



      end do
      end do

      return
      end



      subroutine paw_xc_gen_atomic_densities(ic,lmax2,ispin,
     >                                 rholm,
     >                                 ylm,
     >                                 rhocore,
     >                                 rho)
      implicit none
      integer ic,lmax2,ispin
      double complex rholm(ic,lmax2,ispin)
      double complex ylm(lmax2)
      double precision rhocore(ic)
      double precision rho(ic,ispin)

      integer ms 

      do ms=1,ispin
         call paw_xc_gen_sphere_rho(lmax2,ic,
     >                              rholm(1,1,ms),
     >                              ylm,
     >                              rho(1,ms))
         call daxpy(ic,0.5d0,rhocore,1,rho(1,ms),1)
      end do

      return
      end





      subroutine paw_xc_gen_atomic_gradients(ic,lmax2,ispin,
     >                                       rholm,rholm_prime,
     >                                       ylm,dylm_theta,dylm_phi,
     >                                       r,rhocore,rhocore_prime,
     >                                       rho_prime)
      implicit none
      integer ic,lmax2,ispin
      double complex rholm(ic,lmax2,ispin)
      double complex rholm_prime(ic,lmax2,ispin)
      double complex ylm(lmax2)
      double complex dylm_theta(lmax2)
      double complex dylm_phi(lmax2)
      double precision r(ic)
      double precision rhocore(ic)
      double precision rhocore_prime(ic)

      double precision rho_prime(ic,3,ispin)


      integer ms

      call dcopy(3*ic*ispin,0.0d0,0,rho_prime,1)

      do ms=1,ispin

*        *** find [d/dr paw_rho_ae] on spherical grid ***
         call paw_xc_gen_sphere_rho(lmax2,ic,
     >                              rholm_prime(1,1,ms),
     >                              ylm,
     >                              rho_prime(1,1,ms))

*        *** add core d/dr densities***
         call daxpy(ic,0.5d0,rhocore_prime,1,rho_prime(1,1,ms),1)


         !*** only computate radial derivatives if lmax2==1 ****
         if (lmax2.gt.1) then

*          *** find  (1/r) * d/dtheta paw_rho_ae on spherical grid ***
           call paw_xc_gen_sphere_rho(lmax2,ic,
     >               rholm(1,1,ms),
     >               dylm_theta,
     >               rho_prime(1,2,ms))
           call paw_xc_rho_div_r(ic,r,rho_prime(1,2,ms))

*          *** find  [(r*sin(theta)) * d/dphi paw_rho_ae] on spherical grid ***
           call paw_xc_gen_sphere_rho(lmax2,ic,
     >               rholm(1,1,ms),
     >               dylm_phi,
     >               rho_prime(1,3,ms))
           call paw_xc_rho_div_r(ic,r,rho_prime(1,3,ms))

         end if


      end do

      return
      end


*    ************************************
*    *                                  *
*    *      paw_xc_gen_atomic_agr       *
*    *                                  *
*    ************************************
*
*   This function returns  the absolute values of the gradient.
*
*   Entry - ic     : number of grid points
*           ispin  : restricted/unrestricted
*           rho_prime(ic,3,ispin) : gradient in spherical coordinates
*                                   of atomic spin densites nup and ndn
*
*   Exit - agr_in(*,1): |grad n| if restricted
*   Exit - agr_in(*,3): |grad nup|, |grad ndn|, and |grad n| if unrestricted

      subroutine paw_xc_gen_atomic_agr(ic,lmax2,ispin,
     >                                 rho_prime,
     >                                 agr)
      implicit none
      integer ic,lmax2,ispin
      double precision rho_prime(ic,3,ispin)
      double precision agr(ic,*)                !*(ic,2*ispin-1)

      integer ig,ms


      !*** only computate radial derivatives if lmax2==1 ****
      if (lmax2.eq.1) then

c        **** compute |grad n| ****
         do ig=1,ic
            agr(ig,2*ispin-1)
     >       = dsqrt( (rho_prime(ig,1,1)+rho_prime(ig,1,ispin))**2)
         end do

c        **** compute |grad nup| and |grad ndn| ****
         if (ispin.eq.2) then
           do ms=1,ispin
           do ig=1,ic
              agr(ig,ms) = dsqrt(rho_prime(ig,1,ms)**2)
           end do
           end do
         end if



      !***  computate theta and phi derivatives if lmax2>1 ****
      else
c        **** compute |grad n| ****
         do ig=1,ic   
            agr(ig,2*ispin-1) 
     >       = dsqrt( (rho_prime(ig,1,1)+rho_prime(ig,1,ispin))**2
     >              + (rho_prime(ig,2,1)+rho_prime(ig,2,ispin))**2
     >              + (rho_prime(ig,3,1)+rho_prime(ig,3,ispin))**2)
         end do

c        **** compute |grad nup| and |grad ndn| ****
         if (ispin.eq.2) then
           do ms=1,ispin
           do ig=1,ic
              agr(ig,ms) = dsqrt( rho_prime(ig,1,ms)**2
     >                          + rho_prime(ig,2,ms)**2
     >                          + rho_prime(ig,3,ms)**2)
           end do
           end do
         end if

      end if
      return
      end
                                       




      subroutine paw_xc_gen_dvxc(ic,lmax2,ispin,
     >                           fdn,agr,rho_prime)
      implicit none
      integer ic,lmax2,ispin
      double precision fdn(ic,*)
      double precision agr(ic,*)
      double precision rho_prime(ic,3,ispin)

      integer    i,j,lm,ms,jmax
      double precision drho1,drho2,drhoa


      !*** only computate radial derivatives if lmax2=1 ****
      if (lmax2.eq.1) then
         jmax = 1
      else
         jmax = 3
      end if

*     *** restricted ****
      if (ispin.eq.1) then
         do j=1,jmax
          do i=1,ic
           rho_prime(i,j,1) = (rho_prime(i,j,1)+rho_prime(i,j,1))
     >                       *(fdn(i,1)/agr(i,1))
          end do
         end do

*     *** unrestricted ****
      else
         do j=1,jmax
          do i=1,ic
            drho1 = rho_prime(i,j,1)
            drho2 = rho_prime(i,j,2)
            drhoa = drho1+drho2
            rho_prime(i,j,1) = (drho1/agr(i,1))*fdn(i,1)
     >                       + (drhoa/agr(i,3))*fdn(i,3)
            rho_prime(i,j,2) = (drho2/agr(i,2))*fdn(i,2)
     >                       + (drhoa/agr(i,3))*fdn(i,3)
          end do
         end do
        

      end if

      return
      end





      subroutine paw_xc_addto_dvxclm(ic,lmax2,ispin,
     >                              alpha,
     >                              ylm,
     >                              dvxc_ae,dvxc_ps,
     >                              dvxclm_ae,dvxclm_ps)
      implicit none
      integer ic,lmax2,ispin
      double precision alpha
      double precision ylm(2,*)
      double precision dvxc_ae(ic,3,ispin)
      double precision dvxc_ps(ic,3,ispin)

      double precision dvxclm_ae(2,ic,lmax2,3,ispin)
      double precision dvxclm_ps(2,ic,lmax2,3,ispin)

      integer    i,j,lm,ms,jmax

      !*** only computate radial derivatives if lmax2=1 ****
      if (lmax2.eq.1) then
         jmax = 1
      else
         jmax = 3
      end if

      do ms = 1,ispin
      do j  = 1,jmax
      do lm = 1,lmax2
      do i  = 1,ic
        dvxclm_ae(1,i,lm,j,ms) = dvxclm_ae(1,i,lm,j,ms)
     >                         + dvxc_ae(i,j,ms)*(ylm(1,lm)*alpha)
        dvxclm_ae(2,i,lm,j,ms) = dvxclm_ae(2,i,lm,j,ms)
     >                         - dvxc_ae(i,j,ms)*(ylm(2,lm)*alpha)
      end do
      end do
      end do
      end do

      do ms = 1,ispin
      do j  = 1,jmax
      do lm = 1,lmax2
      do i  = 1,ic
        dvxclm_ps(1,i,lm,j,ms) = dvxclm_ps(1,i,lm,j,ms)
     >                         + dvxc_ps(i,j,ms)*(ylm(1,lm)*alpha)
        dvxclm_ps(2,i,lm,j,ms) = dvxclm_ps(2,i,lm,j,ms)
     >                         - dvxc_ps(i,j,ms)*(ylm(2,lm)*alpha)
      end do
      end do
      end do
      end do

      return
      end





      subroutine paw_xc_gen_dmatr(ng,nb,ic,istart,lmax2,ispin,log_amesh,
     >                           phi_ae,phi_ps,
     >                           phi_ae_prime,phi_ps_prime,
     >                           r,
     >                           dvxclm_ae,dvxclm_ps,
     >                           tmpC,
     >                           dmatr)
      implicit none
      integer ng,nb,ic,istart,lmax2,ispin
      double precision log_amesh

      double precision  phi_ae(ng,nb)
      double precision  phi_ps(ng,nb)
      double precision  phi_ae_prime(ng,nb)
      double precision  phi_ps_prime(ng,nb)
      double precision  r(ic)
      double precision  dvxclm_ae(2,ic,lmax2,3,ispin)
      double precision  dvxclm_ps(2,ic,lmax2,3,ispin)
      double precision  tmpC(2,ic)

      double complex   dmatr(nb,nb,lmax2,3,ispin)

#include "integrate.fh"

*     **** local varialbles ****
      integer ig,i,j,lm,ms
      double precision tmp_ae,tmp_ps


*     **** radial derivative integral ****
      do ms=1,ispin
      do lm=1,lmax2
        do i=1,nb
        do j=i,nb

          do ig=1,ic
            tmp_ae 
     >        = ( phi_ae_prime(ig+istart-1,i)*phi_ae(ig+istart-1,j)
     >          + phi_ae(ig+istart-1,i)*phi_ae_prime(ig+istart-1,j))
     >          /(r(ig)**2)
     >        - 2.0d0*phi_ae(ig+istart-1,i)*phi_ae(ig+istart-1,j)
     >          /(r(ig)**3)
            tmp_ps 
     >        = ( phi_ps_prime(ig+istart-1,i)*phi_ps(ig+istart-1,j)
     >          + phi_ps(ig+istart-1,i)*phi_ps_prime(ig+istart-1,j))
     >          /(r(ig)**2)
     >        - 2.0d0*phi_ps(ig+istart-1,i)*phi_ps(ig+istart-1,j)
     >          /(r(ig)**3)


             tmpC(1,ig) = dvxclm_ae(1,ig,lm,1,ms)*tmp_ae
     >                  - dvxclm_ps(1,ig,lm,1,ms)*tmp_ps
             tmpC(2,ig) = dvxclm_ae(2,ig,lm,1,ms)*tmp_ae
     >                  - dvxclm_ps(2,ig,lm,1,ms)*tmp_ps
          end do
          dmatr(i,j,lm,1,ms)
     >       = c_def_integr1(0,tmpC,2,r,log_amesh,ic,istart)
          if (i.ne.j) dmatr(j,i,lm,1,ms) = dmatr(i,j,lm,1,ms)
        end do
        end do
      end do
      end do


*     *** only comnpute theta and phi integrals if lmax2>1 ****
      if (lmax2.gt.1) then

*     **** theta derivative integral ****
      do ms=1,ispin
      do lm=1,lmax2
        do i=1,nb
        do j=i,nb

          do ig=1,ic
             tmp_ae = phi_ae(ig+istart-1,i)
     >               *phi_ae(ig+istart-1,j)
     >               /(r(ig)**3)
             tmp_ps = phi_ps(ig+istart-1,i)
     >               *phi_ps(ig+istart-1,j)
     >               /(r(ig)**3)

             tmpC(1,ig) = dvxclm_ae(1,ig,lm,2,ms)*tmp_ae
     >                  - dvxclm_ps(1,ig,lm,2,ms)*tmp_ps
             tmpC(2,ig) = dvxclm_ae(2,ig,lm,2,ms)*tmp_ae
     >                  - dvxclm_ps(2,ig,lm,2,ms)*tmp_ps
          end do
          dmatr(i,j,lm,2,ms)
     >      = c_def_integr1(0,tmpC,2,r,log_amesh,ic,istart)
          if (i.ne.j) dmatr(j,i,lm,2,ms) = dmatr(i,j,lm,2,ms)
        end do
        end do
      end do
      end do

*     **** phi   derivative integral ****
      do ms=1,ispin
      do lm=1,lmax2
        do i=1,nb
        do j=i,nb

          do ig=1,ic
             tmp_ae = phi_ae(ig+istart-1,i)
     >               *phi_ae(ig+istart-1,j)
     >               /(r(ig)**3)
             tmp_ps = phi_ps(ig+istart-1,i)
     >               *phi_ps(ig+istart-1,j)
     >               /(r(ig)**3)

             tmpC(1,ig) = dvxclm_ae(1,ig,lm,3,ms)*tmp_ae
     >                  - dvxclm_ps(1,ig,lm,3,ms)*tmp_ps
             tmpC(2,ig) = dvxclm_ae(2,ig,lm,3,ms)*tmp_ae
     >                  - dvxclm_ps(2,ig,lm,3,ms)*tmp_ps
          end do
          dmatr(i,j,lm,3,ms)
     >      = c_def_integr1(0,tmpC,2,r,log_amesh,ic,istart)
          if (i.ne.j) dmatr(j,i,lm,3,ms) = dmatr(i,j,lm,3,ms)
        end do
        end do
      end do
      end do

      end if


      return
      end

