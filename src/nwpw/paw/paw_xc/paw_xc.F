!$Id: paw_xc.F,v 1.5 2003-02-18 01:57:31 marat Exp $
!
!$Log: not supported by cvs2svn $
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
#include "mafdecls.fh"      
#include "paw_ma.fh"
#include "paw_geom.fh"
#include "paw_spher_func.fh"

      logical ok
      integer in,ia
      integer nion
      integer ic
      integer l,m
      integer i_p,i_t
      integer i_ylm
      double precision tmp_theta
      double precision angle_phi
      integer paw_vxc_size

      ok = .true.

      nion = ion_nion()
*     *** spin treatment flag (LDA=1,LSDA=2) ***
      paw_xc_ispin = ispin      

*     *** lmax for multipole expansion of the density***      
      paw_xc_lmax = lmax
      
*     ***paw_xc energies ***
      ok = ok .and. my_alloc(mt_dbl,nion,"paw_xc_e",
     >                                paw_xc_e)

*     *** size for multipole density arrays ***
      ok = ok .and. my_alloc(mt_int,nion,"i_paw_vxc",
     >                                i_paw_vxc)
      if (.not.ok)
     > call errquit("init_paw_vxc_end: error allocating heap",0)

      paw_vxc_size = 0     
      do in=1,nion
        ia = ion_katm(in)
        int_mb(i_paw_vxc(1) + in - 1) = paw_vxc_size
        paw_vxc_size = paw_vxc_size + paw_basis_i_cut(ia)*
     >                                (paw_xc_lmax+1)**2
      end do

*     *** allocate multipole density arrays ****
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_u_ae",
     >                               paw_vxc_u_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_d_ae",
     >                               paw_vxc_d_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_u_ps",
     >                               paw_vxc_u_ps) 
      ok = ok .and. my_alloc(mt_dcpl,paw_vxc_size,"paw_vxc_d_ps",
     >                               paw_vxc_d_ps) 


*     *** spherical grid arrays ***
      if(paw_xc_lmax .eq. 0 ) then
        paw_xc_nphi   = 2
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


      call get_spher_grid(paw_xc_ntheta,paw_xc_nphi,
     >                    dbl_mb(paw_xc_angle_phi(1)),
     >                    dbl_mb(paw_xc_cos_theta(1)),
     >                    dbl_mb(paw_xc_w_theta(1)),
     >                    dbl_mb(paw_xc_w_phi(1)))

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

*     *** temp arrays ***
      paw_xc_work_size = 0
      do in=1,nion
        ia = ion_katm(in)
        ic = paw_basis_i_cut(ia)
        if(ic > paw_xc_work_size) then
           paw_xc_work_size = ic
        end if
      end do
 
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_rho_ae",
     >                               rho_ae)
      ok = ok .and. my_alloc(mt_dbl,paw_xc_work_size*2,"paw_tmp_rho_ps",
     >                               rho_ps)
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
      write(*,*) "paw_xc initialized"
      
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
#include "paw_density_data.fh"
#include "paw_xc_work_data.fh"
#include "paw_xc_data.fh"
#include "integrate.fh"
#include "mafdecls.fh"

      integer in,ia
      integer i_t,i_p
      integer nion 
      integer l,m
      integer ng,ic,ig
      integer i_rho0,i_rho
      integer n_s
      integer i_ylm
      integer nlmi
      integer i_core_ps,i_core_ae
      double precision exc_tmp
      integer i_r
      double precision log_amesh
      complex*16 ylm

      nion = ion_nion() 


*     *** index for spin down in temp arrays ***      
      n_s = paw_xc_work_size 

     
      do in=1,nion
      ia = ion_katm(in)
      ng = paw_basis_ngrid(ia)
      ic = paw_basis_i_cut(ia)
      i_rho0 = int_mb(i_paw_rho(1) + in - 1)
      i_core_ae = paw_basis_i_core_density(ia)
      i_core_ps = paw_basis_i_ps_core_density(ia)
      i_r = paw_basis_i_rgrid(ia)
      log_amesh = paw_basis_log_amesh(ia)

      dbl_mb(paw_xc_e(1)+in-1)=0.0d0
      i_ylm = 0
      do i_t = 1,paw_xc_ntheta
      do i_p = 1,paw_xc_nphi
      
*     *** zero out temp arrays ***
        do ig=1,paw_xc_work_size*2
         dbl_mb(rho_ae(1)+ig-1) = 0.d0
         dbl_mb(rho_ps(1)+ig-1) = 0.d0
         dbl_mb(vxc_ae(1)+ig-1) = 0.d0
         dbl_mb(vxc_ps(1)+ig-1) = 0.d0
         dbl_mb(exc_ae(1)+ig-1) = 0.d0
         dbl_mb(exc_ps(1)+ig-1) = 0.d0
        end do
       
*       *** find densities on spher grid ***          
        nlmi = 0
        do l=0,paw_xc_lmax
        do m=-l,l
          ylm=dcpl_mb(paw_xc_ylm(1)+i_ylm)
          i_ylm = i_ylm + 1
            do ig=1,ic
              i_rho = i_rho0-1+ig+
     >            (nlmi+l+m)*ic

              dbl_mb(rho_ae(1)+ig-1)=
     >             dbl_mb(rho_ae(1)+ig-1)  + 
     >             dble(dcpl_mb(paw_rho_u_ae(1) + i_rho)*ylm)          

              dbl_mb(rho_ae(1)+n_s+ig-1)=
     >             dbl_mb(rho_ae(1)+n_s+ig-1) + 
     >             dble(dcpl_mb(paw_rho_d_ae(1) + i_rho)*ylm)          

              dbl_mb(rho_ps(1)+ig-1)=
     >             dbl_mb(rho_ps(1)+ig-1) + 
     >             dble(dcpl_mb(paw_rho_u_ps(1) + i_rho)*ylm)          

              dbl_mb(rho_ps(1)+n_s+ig-1)=
     >             dbl_mb(rho_ps(1)+n_s+ig-1) + 
     >             dble(dcpl_mb(paw_rho_d_ps(1) + i_rho)*ylm)          

            end do !ig
        end do !m
        nlmi = nlmi+2*l+1
        end do !l

        write(*,*) "******************************"
        write(*,*) "na,i_phi,i_theta",in,i_p,i_t
        write(*,*) "******************************"
        write(*,*) "paw_exc: rho_up_ae(1)",dbl_mb(rho_ae(1))
        write(*,*) "paw_exc: rho_do_ae(1)",dbl_mb(rho_ae(1)+n_s)
        write(*,*) "paw_exc: rho_up_ps(1)",dbl_mb(rho_ps(1))
        write(*,*) "paw_exc: rho_do_ps(1)",dbl_mb(rho_ps(1)+n_s)
        
*       *** add core densities***        
        do ig=1,ic
          dbl_mb(rho_ae(1)+ig-1)=
     >             dbl_mb(rho_ae(1)+ig-1)  + 
     >             0.5d0*dbl_mb(i_core_ae+ig-1)

          dbl_mb(rho_ae(1)+n_s+ig-1)=
     >             dbl_mb(rho_ae(1)+n_s+ig-1) + 
     >             0.5d0*dbl_mb(i_core_ae+ig-1)

              dbl_mb(rho_ps(1)+ig-1)=
     >             dbl_mb(rho_ps(1)+ig-1) + 
     >             0.5d0*dbl_mb(i_core_ps+ig-1)

              dbl_mb(rho_ps(1)+n_s+ig-1)=
     >             dbl_mb(rho_ps(1)+n_s+ig-1) + 
     >             0.5d0*dbl_mb(i_core_ps+ig-1)

        end do !ig

        write(*,*) "paw_exc: rho_up_ae(1)+core",dbl_mb(rho_ae(1))
       
*       *** find exchange-correlation on spher grid ***
        call paw_vosko(n_s,ic,paw_xc_ispin,
     >                 dbl_mb(rho_ae(1)),        
     >                 dbl_mb(vxc_ae(1)),        
     >                 dbl_mb(exc_ae(1)),        
     >                 dbl_mb(xc_temp(1)) )       

        call paw_vosko(n_s,ic,paw_xc_ispin,
     >                 dbl_mb(rho_ps(1)),        
     >                 dbl_mb(vxc_ps(1)),        
     >                 dbl_mb(exc_ps(1)),        
     >                 dbl_mb(xc_temp(1)) )       

        write(*,*) "paw_exc: vxc_up_ae(1)",dbl_mb(vxc_ae(1))
        write(*,*) "paw_exc: vxc_do_ae(1)",dbl_mb(vxc_ae(1)+n_s)
        write(*,*) "paw_exc: vxc_up_ps(1)",dbl_mb(vxc_ps(1))
        write(*,*) "paw_exc: vxc_do_ps(1)",dbl_mb(vxc_ps(1)+n_s)
        
        write(*,*) "paw_exc: exc_up_ae(1)",dbl_mb(exc_ae(1))
        write(*,*) "paw_exc: exc_up_ps(1)",dbl_mb(exc_ps(1))

         do ig=1,ic
          dbl_mb(xc_temp(1)+ig-1)=
     >             (
     >             dbl_mb(rho_ae(1)+ig-1)  + 
     >             dbl_mb(rho_ae(1)+n_s+ig-1)
     >             )*dbl_mb(exc_ae(1)+ig-1)
     >             -
     >             (
     >             dbl_mb(rho_ps(1)+ig-1)+ 
     >             dbl_mb(rho_ps(1)+n_s+ig-1)
     >             )*
     >             dbl_mb(exc_ps(1)+ig-1)

        end do

        write(*,*) "paw_exc: xc_temp(1)",dbl_mb(xc_temp(1))

        exc_tmp = def_integr(0,dbl_mb(xc_temp(1)),2,dbl_mb(i_r),
     >                       log_amesh,ic)

        write(*,*) "paw_exc: exc_tmp",exc_tmp
        dbl_mb(paw_xc_e(1)+in-1)=dbl_mb(paw_xc_e(1)+in-1)+
     >                           exc_tmp*
     >                           dbl_mb(paw_xc_w_theta(1)+i_t-1)*
     >                           dbl_mb(paw_xc_w_phi(1)+i_p-1)


      end do !i_phi
      end do !i_theta
      write(*,*) "paw_xc_e",in,dbl_mb(paw_xc_e(1)+in-1)
      end do !in
      
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
      
      ok = .true.

      ok = ok .and. my_dealloc(exc_ps)
      ok = ok .and. my_dealloc(exc_ae)
      ok = ok .and. my_dealloc(vxc_ps)
      ok = ok .and. my_dealloc(vxc_ae)
      ok = ok .and. my_dealloc(rho_ps)
      ok = ok .and. my_dealloc(rho_ae)

      ok = ok. and. my_dealloc(paw_xc_ylm) 
      if (.not.ok)
     > call errquit("paw_vxc_end1: error freeing heap",0)
      ok = ok. and. my_dealloc(paw_xc_w_theta) 
      if (.not.ok)
     > call errquit("paw_vxc_end2: error freeing heap",0)
      ok = ok. and. my_dealloc(paw_xc_w_phi) 
      if (.not.ok)
     > call errquit("paw_vxc_end3: error freeing heap",0)
      ok = ok. and. my_dealloc(paw_xc_cos_theta) 
      if (.not.ok)
     > call errquit("paw_vxc_end4: error freeing heap",0)
      ok = ok. and. my_dealloc(paw_xc_angle_phi) 
      if (.not.ok)
     > call errquit("paw_vxc_end5: error freeing heap",0)
      ok = ok. and. my_dealloc(paw_vxc_d_ps) 
      if (.not.ok)
     > call errquit("paw_vxc_end6: error freeing heap",0)
      ok = ok .and. my_dealloc(paw_vxc_u_ps) 
      if (.not.ok)
     > call errquit("paw_vxc_end7: error freeing heap",0)
      ok = ok .and. my_dealloc(paw_vxc_d_ae) 
      if (.not.ok)
     > call errquit("paw_vxc_end8: error freeing heap",0)
      ok = ok .and. my_dealloc(paw_vxc_u_ae) 
      if (.not.ok)
     > call errquit("paw_vxc_end9: error freeing heap",0)
      ok = ok .and. my_dealloc(i_paw_vxc)

      if (.not.ok)
     > call errquit("paw_vxc_end: error freeing heap",0)


      end 
