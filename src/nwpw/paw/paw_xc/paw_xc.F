!$Id: paw_xc.F,v 1.11 2003-02-24 23:35:04 marat Exp $
!
!$Log: not supported by cvs2svn $
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
      integer mtr_size
      integer nb
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

     
      
*     *** xc matrix arrays ***
      ok = ok .and. my_alloc(mt_int,nion,"i_paw_xc_matr",
     >                                i_paw_xc_matr)
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating i_paw_xc_matr",0)

      mtr_size = 0
      do in=1,nion
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        int_mb(i_paw_xc_matr(1)+in-1)=mtr_size
        mtr_size=mtr_size+(nb**2)*
     >         ((paw_xc_lmax+1)**2)
      end do

      ok = ok .and. my_alloc(mt_dcpl,mtr_size,"paw_xc_matr_u",
     >                               paw_xc_matr_u) 
      ok = ok .and. my_alloc(mt_dcpl,mtr_size,"paw_xc_matr_d",
     >                               paw_xc_matr_d) 
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_matr",
     >   mtr_size)
 
*     *** xc potential arrays ***
      ok = ok .and. my_alloc(mt_int,nion,"i_paw_xc_pot",
     >                                i_paw_xc_pot)
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating i_paw_xc_pot",0)

      mtr_size = 0
      do in=1,nion
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        int_mb(i_paw_xc_pot(1)+in-1)=mtr_size
        mtr_size=mtr_size+paw_proj_nbasis(ia)**2
      end do

      ok = ok .and. my_alloc(mt_dcpl,mtr_size,"paw_xc_pot_u",
     >                               paw_xc_pot_u) 
      ok = ok .and. my_alloc(mt_dcpl,mtr_size,"paw_xc_pot_d",
     >                               paw_xc_pot_d) 
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_pot",
     >   mtr_size)
 
*     *** multipole density arrays ***
      ok = ok .and. my_alloc(mt_int,nion,"i_paw_vxc",
     >                                i_paw_vxc)
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating i_paw_vxc",0)

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
        if(ic.gt.paw_xc_work_size) then
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
      ok = ok .and. my_alloc(mt_dcpl,paw_xc_work_size,"paw_tmp_x_cmplx",
     >                               xc_c_temp)
      ok = ok .and. my_alloc(mt_dcpl,paw_xc_work_size,
     >                                "paw_tmp_x_cmplx1",
     >                               xc_c_temp1)
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating work arrays",0)
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
#include "paw_proj.fh"
#include "paw_density_data.fh"
#include "paw_xc_work_data.fh"
#include "paw_xc_data.fh"
#include "integrate.fh"
#include "mafdecls.fh"

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      integer i,j
      integer in,ia
      integer i_t,i_p
      integer nion 
      integer l,m
      integer ng,ic,ig
      integer i_rho0,i_rho
      integer n_s
      integer i_ylm,i_ylm1
      integer nlmi
      integer i_core_ps,i_core_ae
      integer i_vxc0,i_vxc
      integer i_pot0,i_pot
      integer i_r
      integer nb,nb2
      integer i_phi_ae0,i_phi_ps0
      integer i_mtr0,i_mtr
      integer li,lj,mi,mj
      integer nilm,njlm
      double precision exc_tmp
      double precision log_amesh
      double precision tmp_ae,tmp_ps
      complex*16 tmp
      complex*16 ylm

      nion = ion_nion() 

*     *** zero out vxc arrays ***
      do i=1,paw_vxc_u_ae(3)
       dcpl_mb(paw_vxc_u_ae(1)+i-1)=dcmplx(0.0,0.0)  
       dcpl_mb(paw_vxc_d_ae(1)+i-1)=dcmplx(0.0,0.0)
       dcpl_mb(paw_vxc_u_ps(1)+i-1)=dcmplx(0.0,0.0) 
       dcpl_mb(paw_vxc_d_ps(1)+i-1)=dcmplx(0.0,0.0)
      end do

      do i=1,paw_xc_pot_u(3)
       dcpl_mb(paw_xc_pot_u(1)+i-1)=dcmplx(0.0,0.0)  
       dcpl_mb(paw_xc_pot_d(1)+i-1)=dcmplx(0.0,0.0)  
      end do


*     *** index for spin down in temp arrays ***      
      n_s = paw_xc_work_size 

     
      do in=1,nion
      ia = ion_katm(in)
      ng = paw_basis_ngrid(ia)
      ic = paw_basis_i_cut(ia)
      i_rho0 = int_mb(i_paw_rho(1) + in - 1)
      i_vxc0 = int_mb(i_paw_vxc(1) + in - 1)
      i_core_ae = paw_basis_i_core_density(ia)
      i_core_ps = paw_basis_i_ps_core_density(ia)
      i_r = paw_basis_i_rgrid(ia)
      log_amesh = paw_basis_log_amesh(ia)

      dbl_mb(paw_xc_e(1)+in-1)=0.0d0
      i_ylm = 0
      i_ylm1 = 0
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

c        write(*,*) "******************************"
c        write(*,*) "na,i_phi,i_theta",in,i_p,i_t
c        write(*,*) "******************************"
c        write(*,*) "paw_exc: rho_up_ae(1)",dbl_mb(rho_ae(1))
c        write(*,*) "paw_exc: rho_do_ae(1)",dbl_mb(rho_ae(1)+n_s)
c        write(*,*) "paw_exc: rho_up_ps(1)",dbl_mb(rho_ps(1))
c        write(*,*) "paw_exc: rho_do_ps(1)",dbl_mb(rho_ps(1)+n_s)
        
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

c        write(*,*) "paw_exc: rho_up_ae(1)+core",dbl_mb(rho_ae(1))
       
*       *** find exchange-correlation on spher grid ***
c        call paw_dirac_x(n_s,ic,paw_xc_ispin,
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

c        write(*,*) "paw_exc: vxc_up_ae(1)",dbl_mb(vxc_ae(1))
c        write(*,*) "paw_exc: vxc_do_ae(1)",dbl_mb(vxc_ae(1)+n_s)
c        write(*,*) "paw_exc: vxc_up_ps(1)",dbl_mb(vxc_ps(1))
c        write(*,*) "paw_exc: vxc_do_ps(1)",dbl_mb(vxc_ps(1)+n_s)
c        
c        write(*,*) "paw_exc: exc_up_ae(1)",dbl_mb(exc_ae(1))
c        write(*,*) "paw_exc: exc_up_ps(1)",dbl_mb(exc_ps(1))

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

c        write(*,*) "paw_exc: xc_temp(1)",dbl_mb(xc_temp(1))

        exc_tmp = def_integr(0,dbl_mb(xc_temp(1)),2,dbl_mb(i_r),
     >                       log_amesh,ic)

        dbl_mb(paw_xc_e(1)+in-1)=dbl_mb(paw_xc_e(1)+in-1)+
     >                           exc_tmp*
     >                           dbl_mb(paw_xc_w_theta(1)+i_t-1)*
     >                           dbl_mb(paw_xc_w_phi(1)+i_p-1)




*       *** find xc mult. moments ***          
        nlmi = 0
        do l=0,paw_xc_lmax
        do m=-l,l
          ylm=dcpl_mb(paw_xc_ylm(1)+i_ylm1)
          ylm=conjg(ylm)*dbl_mb(paw_xc_w_theta(1)+i_t-1)*
     >                   dbl_mb(paw_xc_w_phi(1)+i_p-1)
          i_ylm1 = i_ylm1 + 1

            do ig=1,ic

              i_vxc = i_vxc0-1+ig+
     >            (nlmi+l+m)*ic

              dcpl_mb(paw_vxc_u_ae(1) + i_vxc)=
     >             dcpl_mb(paw_vxc_u_ae(1) + i_vxc) +
     >             dbl_mb(vxc_ae(1)+ig-1)*ylm
               
               dcpl_mb(paw_vxc_d_ae(1) + i_vxc)=
     >             dcpl_mb(paw_vxc_d_ae(1) + i_vxc) +
     >             dbl_mb(vxc_ae(1)+n_s+ig-1)*ylm
             
              dcpl_mb(paw_vxc_u_ps(1) + i_vxc)=
     >             dcpl_mb(paw_vxc_u_ps(1) + i_vxc) +
     >             dbl_mb(vxc_ps(1)+ig-1)*ylm

               dcpl_mb(paw_vxc_d_ps(1) + i_vxc)=
     >             dcpl_mb(paw_vxc_d_ps(1) + i_vxc) +
     >             dbl_mb(vxc_ps(1)+n_s+ig-1)*ylm

             end do

 
        end do
        nlmi = nlmi+2*l+1
        end do

      end do !i_phi
      end do !i_theta

c        write(*,*) "******************************"
c        write(*,*) " paw_xc mult moments
c        write(*,*) "******************************"
c        nlmi = 0
c        do l=0,paw_xc_lmax
c        do m=-l,l
c
c            ig=1
c
c              i_vxc = i_vxc0-1+ig+
c     >            (nlmi+l+m)*ic
c
c              write(*,*) "vxc_u_ae",dcpl_mb(paw_vxc_u_ae(1) + i_vxc),
c     >             in,l,m
c              write(*,*) "vxc_d_ae",dcpl_mb(paw_vxc_d_ae(1) + i_vxc),
c     >             in,l,m
c              write(*,*) "vxc_u_ps",dcpl_mb(paw_vxc_u_ps(1) + i_vxc),
c     >             in,l,m
c              write(*,*) "vxc_d_ps",dcpl_mb(paw_vxc_d_ps(1) + i_vxc),
c     >             in,l,m
c               
c
c 
c        end do
c        nlmi = nlmi+2*l+1
c        end do



      end do !in
      

      do in=1,nion
        ia = ion_katm(in)
        ng = paw_basis_ngrid(ia)
        ic = paw_basis_i_cut(ia)
        i_r = paw_basis_i_rgrid(ia)
        log_amesh = paw_basis_log_amesh(ia)
        nb = paw_basis_nbasis(ia)
        nb2 = nb*nb
        i_phi_ae0 = paw_basis_i_phi_ae(ia)
        i_phi_ps0 = paw_basis_i_phi_ps(ia)
        i_vxc0 = int_mb(i_paw_vxc(1) + in - 1)
        i_mtr0 = int_mb(i_paw_xc_matr(1) + in - 1)
 
        nlmi = 0
        do l=0,paw_xc_ lmax
        do m=-l,l
        do i=1,nb
        li = paw_basis_orb_l(i,ia)
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)
c          write(*,*) "paw_xc loop start",i,j,nb,in
          do ig=1,ic
          
            i_vxc = i_vxc0-1+ig+
     >            (nlmi+l+m)*ic

            tmp_ae = dbl_mb(i_phi_ae0+(i-1)*ng+ig-1)*       
     >               dbl_mb(i_phi_ae0+(j-1)*ng+ig-1)/
     >               dbl_mb(i_r+ig-1)**2 

            tmp_ps = dbl_mb(i_phi_ps0+(i-1)*ng+ig-1)*       
     >               dbl_mb(i_phi_ps0+(j-1)*ng+ig-1)/
     >               dbl_mb(i_r+ig-1)**2 


          dcpl_mb(xc_c_temp(1)+ig-1)= 
     >        dcpl_mb(paw_vxc_u_ae(1) + i_vxc)*tmp_ae -
     >        dcpl_mb(paw_vxc_u_ps(1) + i_vxc)*tmp_ps
 
          dcpl_mb(xc_c_temp1(1)+ig-1)= 
     >        dcpl_mb(paw_vxc_d_ae(1) + i_vxc)*tmp_ae -
     >        dcpl_mb(paw_vxc_d_ps(1) + i_vxc)*tmp_ps
 

          end do !ig

c          write(*,*) "xc_c_temp(1)",i,j,l,m,dcpl_mb(xc_c_temp(1))
          i_mtr = i_mtr0-1+j+(i-1)*nb+(nlmi+l+m)*nb2

          tmp = c_def_integr(0,dcpl_mb(xc_c_temp(1)),2,dbl_mb(i_r),
     >                       log_amesh,ic)

          dcpl_mb(paw_xc_matr_u(1)+i_mtr)=tmp

          tmp = c_def_integr(0,dcpl_mb(xc_c_temp1(1)),2,dbl_mb(i_r),
     >                       log_amesh,ic)

          dcpl_mb(paw_xc_matr_d(1)+i_mtr)=tmp
cc          write(*,*) "exc_matrix_up",tmp

        end do !j
        end do !i
        end do !m
        nlmi=nlmi+2*l+1
        end do !l



        i_pot0 = int_mb(i_paw_xc_pot(1) + in - 1)
        nlmi = 0
        do l=0,paw_xc_lmax
        nilm = 0
        do i=1,nb
        li = paw_basis_orb_l(i,ia)
        njlm = 0
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)
        if( (l.le.(li+lj)) .and. (l.ge.ABS(li-lj))) then
        do m=-l,l
          do mi=-li,li
          do mj=-lj,lj
          if(m.eq.(mi-mj)) then
            i_mtr = i_mtr0-1+j+(i-1)*nb+(nlmi+l+m)*nb2
            i_pot = i_pot0-1+(njlm+lj+mj+1)+
     >                (nilm+li+mi)*paw_proj_nbasis(ia)


c            write(*,*) "xc_pot",i_pot
c            write(*,*) "xc_mat",i_mtr
            dcpl_mb(paw_xc_pot_u(1)+i_pot)=
     >            dcpl_mb(paw_xc_pot_u(1)+i_pot) +
     >            gen_gaunt_coeff(l,m,li,mi,lj,mj)*       
     >            dcpl_mb(paw_xc_matr_u(1)+i_mtr)          
             
c            write(*,*) "i,mi,j,mj,l,m",i,mi,j,mj,l,m
c            write(*,*) "xc_matrix",dcpl_mb(paw_xc_matr_u(1)+i_mtr)
c            write(*,*) "xc_pot",dcpl_mb(paw_xc_pot_u(1)+i_pot)
            
            dcpl_mb(paw_xc_pot_d(1)+i_pot)=
     >            dcpl_mb(paw_xc_pot_d(1)+i_pot) +
     >            gen_gaunt_coeff(l,m,li,mi,lj,mj)*       
     >            dcpl_mb(paw_xc_matr_d(1)+i_mtr)          

          end if
            end do !mi
          end do !mj
        end do !m

        end if
        njlm = njlm + 2*lj+1
        end do !j
        nilm = nilm + 2*li+1
        end do !i  
        nlmi=nlmi+2*l+1
        end do !l
        
*       *** print out exc potenetial ***
        i_pot0 = int_mb(i_paw_xc_pot(1) + in - 1)
        nilm = 0
        do i=1,nb
        li = paw_basis_orb_l(i,ia)
        njlm = 0
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)
          do mi=-li,li
          do mj=-lj,lj
            i_pot = i_pot0-1+(njlm+lj+mj+1)+
     >                (nilm+li+mi)*paw_proj_nbasis(ia)



c            write(*,*) "i,mi,j,mj,na",i,mi,j,mj,in
c            write(*,*) "paw_xc_pot_up",dcpl_mb(paw_xc_pot_u(1)+i_pot)
c            write(*,*) "paw_xc_pot_down",dcpl_mb(paw_xc_pot_u(1)+i_pot)
             
            

            end do !mi
          end do !mj
        njlm = njlm + 2*lj+1
        end do !j
        nilm = nilm + 2*li+1
        end do !i  
 
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
     
      ok = ok .and. my_dealloc(xc_temp)
      ok = ok .and. my_dealloc(xc_c_temp)
      ok = ok .and. my_dealloc(xc_c_temp1)

      ok = ok .and. my_dealloc(paw_xc_pot_d)
      ok = ok .and. my_dealloc(paw_xc_pot_u)
      ok = ok .and. my_dealloc(i_paw_xc_pot)
      ok = ok .and. my_dealloc(paw_xc_matr_d)
      ok = ok .and. my_dealloc(paw_xc_matr_u)
      ok = ok .and. my_dealloc(i_paw_xc_matr)
      ok = ok .and. my_dealloc(paw_xc_e)


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
#include "mafdecls.fh"

      integer in

      paw_energy_xc_atom = 0
      do in=1,paw_xc_e(3)
       paw_energy_xc_atom = paw_energy_xc_atom + 
     >                      dbl_mb(paw_xc_e(1)+in-1)  
      end do

      end
