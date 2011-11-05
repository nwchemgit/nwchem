c $Id$

*     *********************************
*     *                               *
*     *         nwpw_xc_init          *
*     *                               *
*     *********************************
      subroutine nwpw_xc_init(ispin,lmax,nion,nkatm,
     >                        nprj,nbasis,n1dgrid,psp_type,
     >                        nprj_max0,l_prj,m_prj,b_prj)
      implicit none
      integer ispin,lmax,nion,nkatm
      integer nprj(*),nbasis(*),n1dgrid(*),psp_type(*)
      integer nprj_max0
      integer l_prj(nprj_max0,*),m_prj(nprj_max0,*),b_prj(nprj_max0,*)


#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

c     **** local variables ****
      logical ok
      integer ii,ia,nsize
      integer ic
      integer l,m
      integer i_p,i_t
      integer i_tlm,i_plm,i_rlm
      integer mtr_size
      integer nb
      real*8  tmp_theta,cs_theta,tmp_phi
      real*8  angle_phi
      integer paw_vxc_size

c     **** external functions ****
      integer  control_gga
      external control_gga
      real*8   rtheta_lm,drtheta_lm
      external rtheta_lm,drtheta_lm


      call nwpw_timing_start(4)
      ok = .true.

      nprj_max    = 0
      nbasis_max  = 0
      n1dgrid_max = 0
      do ia=1,nkatm
         if (nprj(ia)   .gt.nprj_max)    nprj_max = nprj(ia)
         if (nbasis(ia) .gt.nbasis_max)  nbasis_max = nbasis(ia)
         if (n1dgrid(ia).gt.n1dgrid_max) n1dgrid_max = n1dgrid(ia)
      end do

*     *** spin treatment flag (LDA=1,LSDA=2) ***
      paw_xc_ispin = ispin      

*     *** gga calculation ***
      paw_xc_gga = control_gga()

*     *** lmax for multipole expansion of the density***      
      paw_xc_lmax = lmax
      
*     ***paw_xc energies ***
      ok = MA_alloc_get(mt_dbl,nion,"paw_xc_e",paw_xc_e(2),paw_xc_e(1))

      
*     *** xc matrix arrays - nbasis_max**2 * (lmax+1)**2 ***
      mtr_size = (nbasis_max**2)*((paw_xc_lmax+1)**2)
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*mtr_size,"paw_xc_matr_u",
     >                  paw_xc_matr(2),paw_xc_matr(1))

      if (paw_xc_gga.ge.10) then
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,3*ispin*mtr_size,"paw_xc_dmatr_u",
     >                  paw_xc_dmatr(2),paw_xc_dmatr(1))
      end if
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_matr",
     >   mtr_size,0)
 

      mtr_size        = nprj_max*nprj_max
      paw_xc_pot_size = mtr_size
      ok = MA_alloc_get(mt_dbl,ispin*mtr_size,"paw_xc_pot_u",
     >                  paw_xc_pot(2),paw_xc_pot(1)) 
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_pot",
     >   mtr_size,0)



*     *** spherical grid arrays ***
      if(paw_xc_lmax .eq. 0 ) then
        paw_xc_nphi   = 1
        paw_xc_ntheta = 1
      else
        paw_xc_nphi   = 3*paw_xc_lmax
        paw_xc_ntheta = 3*paw_xc_lmax                      
      end if
      

      ok = MA_alloc_get(mt_dbl,paw_xc_nphi,"paw_xc_angle_phi",
     >                  paw_xc_angle_phi(2),paw_xc_angle_phi(1))

      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,paw_xc_ntheta,"paw_xc_cos_theta",
     >                  paw_xc_cos_theta(2),paw_xc_cos_theta(1))

      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,paw_xc_nphi,"paw_xc_w_phi",
     >                  paw_xc_w_phi(2),paw_xc_w_phi(1))

      ok = ok.and.
     >     MA_alloc_get(mt_dbl,paw_xc_ntheta,"paw_xc_w_theta",
     >                  paw_xc_w_theta(2),paw_xc_w_theta(1))

      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_tlm",
     >              paw_xc_tlm(2),paw_xc_tlm(1))


*     **** used for generating derivatives of ylm's ****
      if (paw_xc_gga.ge.10) then

c     **** derivatives wrt to theta ****
      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_dtlm_theta",
     >              paw_xc_dtlm_theta(2),paw_xc_dtlm_theta(1))

c     **** derivatives wrt to phi ****
      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(paw_xc_lmax+1)**2,
     >              "paw_xc_dtlm_phi",
     >              paw_xc_dtlm_phi(2),paw_xc_dtlm_phi(2))

      end if

      call nwpw_get_spher_grid(paw_xc_ntheta,paw_xc_nphi,
     >                         dbl_mb(paw_xc_angle_phi(1)),
     >                         dbl_mb(paw_xc_cos_theta(1)),
     >                         dbl_mb(paw_xc_w_theta(1)),
     >                         dbl_mb(paw_xc_w_phi(1)))

c     **** define tlm's ****
      i_tlm = 0
      do i_t=1,paw_xc_ntheta
      do i_p=1,paw_xc_nphi
         do l=0,paw_xc_lmax
         do m=-l,l
           tmp_theta = rtheta_lm(l,m,
     >              dbl_mb(paw_xc_cos_theta(1)+i_t-1))
           angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
           if (m.lt.0) then
              tmp_phi = dsin(abs(m)*angle_phi)
           else if (m.gt.0) then
              tmp_phi = dcos(abs(m)*angle_phi)
           else
              tmp_phi = 1.0d0
           end if
           dcpl_mb(paw_xc_tlm(1)+i_tlm) = tmp_theta*tmp_phi
           i_tlm = i_tlm + 1                          
         end do
         end do
      end do
      end do

      if (paw_xc_gga.ge.10) then

c       **** define derivative wrt to theta ****
        i_tlm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,paw_xc_lmax
           do m=-l,l
              tmp_theta = drtheta_lm(l,m,
     >                               dbl_mb(paw_xc_cos_theta(1)+i_t-1))
              angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
              if (m.lt.0) then
                 tmp_phi = dsin(abs(m)*angle_phi)
              else if (m.gt.0) then
                 tmp_phi = dcos(abs(m)*angle_phi)
              else
                 tmp_phi = 1.0d0
              end if
              dcpl_mb(paw_xc_dtlm_theta(1)+i_tlm) = tmp_theta*tmp_phi
              i_tlm = i_tlm + 1 
           end do
           end do
        end do
        end do

c       **** define derivative wrt to phi ****
        i_tlm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,paw_xc_lmax
           do m=-l,l
              tmp_theta = rtheta_lm(l,m,
     >                              dbl_mb(paw_xc_cos_theta(1)+i_t-1))

              angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
              if (m.lt.0) then
                 tmp_phi = -abs(m)*dcos(abs(m)*angle_phi)
              else if (m.gt.0) then
                 tmp_phi = abs(m)*dsin(abs(m)*angle_phi)
              else
                 tmp_phi = 0.0d0
              end if
              dbl_mb(paw_xc_dtlm_phi(1)+ i_tlm) = tmp_theta*tmp_phi
              i_tlm = i_tlm + 1
           end do
           end do
        end do
        end do

      end if

*     *** temp arrays ***
      nsize = 2*n1dgrid_max
      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_rho_ae",
     >                  rho_ae(2),rho_ae(1))
      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_rho_ps",
     >                  rho_ps(2),rho_ps(1))
      
*     **** allocate gradient's and agr's ****
      if (paw_xc_gga.ge.10) then
       nsize = ispin*3*n1dgrid_max
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >              "paw_tmp_drho_ae",
     >              rho_ae_prime(2),rho_ae_prime(1))
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_drho_ps",
     >                   rho_ps_prime(2),rho_ps_prime(1))
       nsize = (2*ispin-1)*n1dgrid_max
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_agr_ae",
     >                   agr_ae(2),agr_ae(1))
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_agr_ps",
     >                    agr_ps(2),agr_ps(1))
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_fdn_ae",
     >                   fdn_ae(2),fdn_ae(1))
       ok = ok.and.
     >      MA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_fdn_ps",
     >                   fdn_ps(2),fdn_ps(1))
      end if

      nsize = 2*n1dgrid_max
      ok = ok.and. 
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_vxc_ae",
     >                  vxc_ae(2),vxc_ae(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_vxc_ps",
     >                  vxc_ps(2),vxc_ps(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_exc_ae",
     >                  exc_ae(2),exc_ae(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_exc_ps",
     >                  exc_ps(2),exc_ps(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,nsize,"paw_tmp_x",
     >                   xc_temp(2),xc_temp(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,nsize,
     >                  "paw_tmp_x_cmplx",
     >                  xc_c_temp(2),xc_c_temp(1))


*     *** allocate vxclm multipole expansion  arrays ****
      nsize = n1dgrid_max*(paw_xc_lmax+1)**2

      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,"paw_vxc_u_ae",
     >                  paw_vxc_ae(2),paw_vxc_ae(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,"paw_vxc_u_ps",
     >                  paw_vxc_ps(2),paw_vxc_ps(1))

*     *** allocate dvxclm multipole expansion  arrays ****
      if (paw_xc_gga.ge.10) then
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,3*ispin*nsize,
     >                  "paw_dvxc_u_ae",
     >                  paw_dvxc_ae(2),paw_dvxc_ae(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,3*ispin*nsize,
     >                  "paw_dvxc_u_ps",
     >                  paw_dvxc_ps(2),paw_dvxc_ps(1))
      end if

*     *** allocate rholm multipole expansion  arrays ****
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,"paw_rho2_ae",
     >                  paw_rho2_ae(2),paw_rho2_ae(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,"paw_rho2_ps",
     >                  paw_rho2_ps(2),paw_rho2_ps(1))

*     *** allocate rholm_prime multipole expansion arrays - need for GGA's****
      if (paw_xc_gga.ge.10) then
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,
     >                  "paw_rho2_ae_prime",
     >                   paw_rho2_ae_prime(2),paw_rho2_ae_prime(1))
      ok = ok.and.
     >     MA_alloc_get(mt_dbl,ispin*nsize,
     >                  "paw_rho2_ps_prime",
     >                  paw_rho2_ps_prime(2),paw_rho2_ps_prime(1))
      end if
      if (.not.ok)
     > call errquit("nwpw_xc_init: error allocating work arrays",0,0)
      
      call nwpw_init_gntxc(paw_xc_lmax)
      if ((paw_xc_gga.ge.10).and.(paw_xc_lmax.ge.1)) then
        call nwpw_init_gntxc2(paw_xc_lmax)
        call nwpw_init_gntxc3(paw_xc_lmax)
      end if

      call nwpw_timing_end(4)
      return
      end 

*     ********************************************
*     *                                          *
*     *             nwpw_xc_solve                *
*     *                                          *
*     ********************************************
      subroutine nwpw_xc_solve(ii,ia,n1dgrid,nbasis,
     >               phi_ae,phi_ps,phi_ae_prime,phi_ps_prime,
     >               core_ae,core_ps,core_ae_prime,core_ps_prime,
     >               regrid,log_amesh,
     >               ispin,ne,nprj,sw1,sw2)
      implicit none
      integer ii,ia
      integer n1dgrid,nbasis
      real*8  phi_ae(n1dgrid,nbasis)
      real*8  phi_ps(n1dgrid,nbasis)
      real*8  phi_ae_prime(n1dgrid,nbasis)
      real*8  phi_ps_prime(n1dgrid,nbasis)
      real*8  core_ae(n1dgrid)
      real*8  core_ps(n1dgrid)
      real*8  core_ae_prime(n1dgrid)
      real*8  core_ps_prime(n1dgrid)
      real*8  rgrid(n1dgrid)
      real*8  log_amesh

      integer ispin,ne(2)
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      integer i,j,np
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
      integer i_mtr0,i_mtr
      integer li,lj,mi,mj
      integer nilm,njlm
      integer lmax2
      double precision exc_tmp,gnt1
      double precision tmp_ae,tmp_ps
      real*8 tmp
      real*8 tlm

      call nwpw_timing_start(4)
      call nwpw_timing_start(22)
      call Parallel_np(np)

      i_tlm = 0
      i_tlm1 = 0


*     *** index for spin down in temp arrays ***      
      lmax2  = (paw_xc_lmax+1)**2
      dbl_mb(paw_xc_e(1)+ii-1)=0.0d0
     
*     *** zero out vxclm multipole arrays ***
      nsize = ispin*n1dgrid_max*lmax2
      call dcopy(nsize,0.0d0,0,dbl_mb(paw_vxc_ae(1)),1)
      call dcopy(nsize,0.0d0,0,dbl_mb(paw_vxc_ps(1)),1)

*     *** zero out dvxclm multipole arrays ***
      if (paw_xc_gga.ge.10) then
      call dcopy(3*nsize,0.0d0,0,dbl_mb(paw_dvxc_ae(1)),1)
      call dcopy(3*nsize,0.0d0,0,dbl_mb(paw_dvxc_ps(1)),1)
      end if

      do i_t = 1,paw_xc_ntheta
      do i_p = 1,paw_xc_nphi
      
*       *** zero out temp arrays ***
        nsize = ispin*n1dgrid_max
        call dcopy(nsize,0.0d0,0,dbl_mb(rho_ae(1)),1)
        call dcopy(nsize,0.0d0,0,dbl_mb(rho_ps(1)),1)
        call dcopy(nsize,0.0d0,0,dbl_mb(vxc_ae(1)),1)
        call dcopy(nsize,0.0d0,0,dbl_mb(vxc_ps(1)),1)
        call dcopy(nsize,0.0d0,0,dbl_mb(exc_ae(1)),1)
        call dcopy(nsize,0.0d0,0,dbl_mb(exc_ps(1)),1)

       
*       *** find rholm multipole expansion on spherical grid ***          
        call paw_density_solve2(in,paw_xc_lmax,
     >                          ic,lmax2,paw_xc_ispin,
     >                          dbl_mb(paw_rho2_ae(1)),
     >                          dbl_mb(paw_rho2_ps(1)))

*       *** generate atomic densities on spherical grid ***          
        call paw_xc_gen_atomic_densities(ic,lmax2,paw_xc_ispin,
     >                                    dcpl_mb(paw_rho2_ae(1)),
     >                                    dcpl_mb(paw_xc_tlm(1)+i_ylm),
     >                                    dbl_mb(i_core_ae+istart-1),
     >                                    dbl_mb(rho_ae(1)))
        call paw_xc_gen_atomic_densities(ic,lmax2,paw_xc_ispin,
     >                                    dcpl_mb(paw_rho2_ps(1)),
     >                                    dcpl_mb(paw_xc_tlm(1)+i_ylm),
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
     >                 dbl_mb(paw_rho2_ae(1)),
     >                 dbl_mb(paw_rho2_ae_prime(1)),
     >                 dbl_mb(paw_xc_tlm(1)+i_ylm),
     >                 dbl_mb(paw_xc_dtlm_theta(1)+i_ylm),
     >                 dbl_mb(paw_xc_dtlm_phi(1)+i_ylm),
     >                 dbl_mb(i_r+istart-1),
     >                 dbl_mb(i_core_ae+istart-1),
     >                 dbl_mb(i_core_ae_prime+istart-1),
     >                 dbl_mb(rho_ae_prime(1)))
          call paw_xc_gen_atomic_gradients(ic,lmax2,paw_xc_ispin,
     >                 dbl_mb(paw_rho2_ps(1)),
     >                 dbl_mb(paw_rho2_ps_prime(1)),
     >                 dbl_mb(paw_xc_tlm(1)+i_ylm),
     >                 dbl_mb(paw_xc_dtlm_theta(1)+i_ylm),
     >                 dbl_mb(paw_xc_dtlm_phi(1)+i_ylm),
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
        i_tlm = i_tlm + lmax2

*       ************************************************************************************
*       *** find vxclm the multipole expansion of atomic vxc=df/dn  or (df/dnup,df/dndn) ***          
*       ************************************************************************************
        call paw_xc_addto_vxclm(ic,lmax2,paw_xc_ispin,
     >                          (dbl_mb(paw_xc_w_theta(1)+i_t-1)
     >                          *dbl_mb(paw_xc_w_phi(1)  +i_p-1)),
     >                          dcpl_mb(paw_xc_tlm(1)+i_ylm1),
     >                          dbl_mb(vxc_ae(1)),
     >                          dbl_mb(vxc_ps(1)),
     >                          dcpl_mb(paw_vxc_ae(1)),
     >                          dcpl_mb(paw_vxc_ps(1)))
        i_tlm1 = i_tlm1 + lmax2





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
     >                       dbl_mb(paw_vxc_ae(1)),
     >                       dbl_mb(paw_vxc_ps(1)),
     >                       dbl_mb(xc_c_temp(1)),
     >                       dbl_mb(paw_xc_matr(1)))


*       **** xc potential  non-local matrix elements ****
        i_pot0 = int_mb(i_paw_xc_pot(1) + in - 1)
        do ms=1,paw_xc_ispin
           call paw_addto_potential_gntxc(ia,
     >           dcpl_mb(paw_xc_matr(1)+(ms-1)*nb2*lmax2),
     >           dcpl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
        end do


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
     >           dbl_mb(paw_xc_dmatr(1)+(ms-1)*nb2*lmax2*3),
     >           dbl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
          end do

c         **** theta and phi parts ****
          if (lmax2.gt.1) then

c           **** theta derivatives ****
            do ms=1,paw_xc_ispin
             call paw_addto_potential_gntxc2(ia,
     >         dbl_mb(paw_xc_dmatr(1)+nb2*lmax2+(ms-1)*nb2*lmax2*3),
     >         dbl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
            end do

c           **** phi derivatives ****
            do ms=1,paw_xc_ispin
             call paw_addto_potential_gntxc3(ia,
     >         dbl_mb(paw_xc_dmatr(1)+nb2*lmax2*2+(ms-1)*nb2*lmax2*3),
     >         dbl_mb(paw_xc_pot(1)+i_pot0+(ms-1)*paw_xc_pot_size))
            end do


         end if !**lmax2>1**

        
        end if !**gga**
 
      end do !in

      if (np.gt.1)  then

c       **** xc non-local matrix elements ****
        call D3dB_Vector_Sumall(2*paw_xc_pot(3),
     >                          dbl_mb(paw_xc_pot(1)))

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
      subroutine nwpw_xc_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      logical ok
      integer ms,i
      
      call nwpw_timing_start(4)
      ok = .true.
      ok = ok.and.MA_free_heap(exc_ps(2))
      ok = ok.and.MA_free_heap(exc_ae(2))
      ok = ok.and.MA_free_heap(vxc_ps(2))
      ok = ok.and.MA_free_heap(vxc_ae(2))
      ok = ok.and.MA_free_heap(rho_ps(2))
      ok = ok.and.MA_free_heap(rho_ae(2))
      ok = ok.and.MA_free_heap(xc_temp(2))

      ok = ok.and.MA_free_heap(i_paw_xc_pot(2))
      ok = ok.and.MA_free_heap(paw_xc_e(2))
      ok = ok.and.MA_free_heap(paw_xc_ylm(2)) 
      ok = ok.and.MA_free_heap(paw_xc_w_theta)(2) 
      ok = ok.and.MA_free_heap(paw_xc_w_phi(2))
      ok = ok.and.MA_free_heap(paw_xc_cos_theta(2)) 
      ok = ok.and.MA_free_heap(paw_xc_angle_phi(2))

      ok = ok.and.MA_free_heap(xc_c_temp(2))
      ok = ok.and.MA_free_heap(paw_xc_pot(2))
      ok = ok.and.MA_free_heap(paw_xc_matr(2))
      ok = ok.and.MA_free_heap(paw_vxc_ae(2))
      ok = ok.and.MA_free_heap(paw_vxc_ps(2))
      ok = ok.and.MA_free_heap(paw_rho2_ae(2))
      ok = ok.and.MA_free_heap(paw_rho2_ps(2))

      if (paw_xc_gga.ge.10) then
         ok = ok.and.MA_free_heap(paw_rho2_ae_prime(2))
         ok = ok.and.MA_free_heap(paw_rho2_ps_prime(2))
         ok = ok.and.MA_free_heap(rho_ps_prime(2))
         ok = ok.and.MA_free_heap(rho_ae_prime(2))
         ok = ok.and.MA_free_heap(agr_ae(2))
         ok = ok.and.MA_free_heap(agr_ps(2))
         ok = ok.and.MA_free_heap(fdn_ae(2))
         ok = ok.and.MA_free_heap(fdn_ps(2))
         ok = ok.and.MA_free_heap(paw_xc_dylm_theta(2))
         ok = ok.and.MA_free_heap(paw_xc_dylm_phi(2))
         ok = ok.and.MA_free_heap(paw_dvxc_ae(2))
         ok = ok.and.MA_free_heap(paw_dvxc_ps(2))
         ok = ok.and.MA_free_heap(paw_xc_dmatr(2))
      end if
      if (.not.ok)
     > call errquit("nwpw_xc_end: error freeing heap",0,0)


      call nwpw_end_gntxc()
      if ((paw_xc_gga.ge.10).and.(paw_xc_lmax.ge.1)) then
        call nwpw_end_gntxc2()
        call nwpw_end_gntxc3()
      end if

      call nwpw_timing_end(4)
      return
      end 


*      ******************************************
*      *                                        *
*      *          nwpw_energy_xc_atom           *
*      *                                        *
*      ******************************************
      real*8 function nwpw_energy_xc_atom()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      integer ii

      nwpw_energy_xc_atom = 0
      do ii=1,paw_xc_e(3)
       nwpw_energy_xc_atom = nwpw_energy_xc_atom + 
     >                      dbl_mb(paw_xc_e(1)+ii-1)  
      end do
      return
      end



      subroutine paw_xc_pot_up_ptr(ptr)
      implicit none
      integer ptr

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"


      ptr = paw_xc_pot(1)
      return
      end

      subroutine paw_xc_pot_down_ptr(ptr)
      implicit none
      integer ptr

#include "mafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      ptr = paw_xc_pot(1) + paw_xc_pot_size
      return
      end


      subroutine nwpw_xc_gen_sphere_rho(lmax2,ic,rholm,Ylm,rho)
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




      subroutine nwpw_xc_rho_div_r(ic,r,rho)
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



*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_addto_vxclm                      *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_addto_vxclm(ic,lmax2,ispin,
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

*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_gen_matr                         *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_gen_matr(ng,nb,ic,istart,lmax2,ispin,log_amesh,
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
      double precision  vxclm_ae(ic,lmax2,ispin)
      double precision  vxclm_ps(ic,lmax2,ispin)
      double precision  tmpC(ic)

      double complex   matr(nb,nb,lmax2,ispin)


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
             tmpC(ig) = vxclm_ae(ig,lm,ms)*tmp_ae
     >                - vxclm_ps(ig,lm,ms)*tmp_ps
          end do
          matr(i,j,lm,ms)=log_integrate_def(0,tmpC,2,r,log_amesh,ic)
          if (i.ne.j) matr(j,i,lm,ms) = matr(i,j,lm,ms)
        end do
        end do
      end do
      end do

      return
      end


*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_gen_atomic_densities             *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_gen_atomic_densities(ic,lmax2,ispin,
     >                                 rholm,
     >                                 tlm,
     >                                 rhocore,
     >                                 rho)
      implicit none
      integer ic,lmax2,ispin
      real*8 rholm(ic,lmax2,ispin)
      real*8 tlm(lmax2)
      real*8 rhocore(ic)
      real*8 rho(ic,ispin)

      integer ms 

      do ms=1,ispin
         call nwpw_xc_gen_sphere_rho(lmax2,ic,
     >                              rholm(1,1,ms),
     >                              ylm,
     >                              rho(1,ms))
         call daxpy(ic,0.5d0,rhocore,1,rho(1,ms),1)
      end do

      return
      end



*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_gen_atomic_gradients             *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_gen_atomic_gradients(ic,lmax2,ispin,
     >                                       rholm,rholm_prime,
     >                                       ylm,dylm_theta,dylm_phi,
     >                                       r,rhocore,rhocore_prime,
     >                                       rho_prime)
      implicit none
      integer ic,lmax2,ispin
      real*8 rholm(ic,lmax2,ispin)
      real*8 rholm_prime(ic,lmax2,ispin)
      real*8 tlm(lmax2)
      real*8 dtlm_theta(lmax2)
      real*8 dtlm_phi(lmax2)
      real*8 r(ic)
      real*8 rhocore(ic)
      real*8 rhocore_prime(ic)

      real*8 rho_prime(ic,3,ispin)


      integer ms

      call dcopy(3*ic*ispin,0.0d0,0,rho_prime,1)

      do ms=1,ispin

*        *** find [d/dr paw_rho_ae] on spherical grid ***
         call nwpw_xc_gen_sphere_rho(lmax2,ic,
     >                              rholm_prime(1,1,ms),
     >                              tlm,
     >                              rho_prime(1,1,ms))

*        *** add core d/dr densities***
         call daxpy(ic,0.5d0,rhocore_prime,1,rho_prime(1,1,ms),1)


         !*** only computate radial derivatives if lmax2==1 ****
         if (lmax2.gt.1) then

*          *** find  (1/r) * d/dtheta paw_rho_ae on spherical grid ***
           call nwpw_xc_gen_sphere_rho(lmax2,ic,
     >               rholm(1,1,ms),
     >               dtlm_theta,
     >               rho_prime(1,2,ms))
           call nwpw_xc_rho_div_r(ic,r,rho_prime(1,2,ms))

*          *** find  [(r*sin(theta)) * d/dphi paw_rho_ae] on spherical grid ***
           call nwpw_xc_gen_sphere_rho(lmax2,ic,
     >               rholm(1,1,ms),
     >               dtlm_phi,
     >               rho_prime(1,3,ms))
           call nwpw_xc_rho_div_r(ic,r,rho_prime(1,3,ms))

         end if


      end do

      return
      end


*    ************************************
*    *                                  *
*    *      nwpw_xc_gen_atomic_agr      *
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

      subroutine nwpw_xc_gen_atomic_agr(ic,lmax2,ispin,
     >                                 rho_prime,
     >                                 agr)
      implicit none
      integer ic,lmax2,ispin
      real*8  rho_prime(ic,3,ispin)
      real*8  agr(ic,*)                !*(ic,2*ispin-1)

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
                                       


*    ************************************
*    *                                  *
*    *      nwpw_xc_gen_dvxc            *
*    *                                  *
*    ************************************
      subroutine nwpw_xc_gen_dvxc(ic,lmax2,ispin,
     >                           fdn,agr,rho_prime)
      implicit none
      integer ic,lmax2,ispin
      real*8  fdn(ic,*)
      real*8  agr(ic,*)
      real*8  rho_prime(ic,3,ispin)

      integer i,j,lm,ms,jmax
      real*8  drho1,drho2,drhoa

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

*    ************************************
*    *                                  *
*    *      nwpw_xc_addto_dvxclm        *
*    *                                  *
*    ************************************
      subroutine nwpw_xc_addto_dvxclm(ic,lmax2,ispin,
     >                              alpha,
     >                              tlm,
     >                              dvxc_ae,dvxc_ps,
     >                              dvxclm_ae,dvxclm_ps)
      implicit none
      integer ic,lmax2,ispin
      real*8  alpha
      real*8  tlm(*)
      real*8  dvxc_ae(ic,3,ispin)
      real*8  dvxc_ps(ic,3,ispin)

      real*8  dvxclm_ae(ic,lmax2,3,ispin)
      real*8  dvxclm_ps(ic,lmax2,3,ispin)

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
        dvxclm_ae(i,lm,j,ms) = dvxclm_ae(i,lm,j,ms)
     >                         + dvxc_ae(i,j,ms)*(tlm(lm)*alpha)
      end do
      end do
      end do
      end do

      do ms = 1,ispin
      do j  = 1,jmax
      do lm = 1,lmax2
      do i  = 1,ic
        dvxclm_ps(i,lm,j,ms) = dvxclm_ps(i,lm,j,ms)
     >                         + dvxc_ps(i,j,ms)*(tlm(lm)*alpha)
      end do
      end do
      end do
      end do

      return
      end




*    ************************************
*    *                                  *
*    *      nwpw_xc_gen_dmatr           *
*    *                                  *
*    ************************************

      subroutine nwpw_xc_gen_dmatr(ng,nb,ic,istart,lmax2,
     >                           ispin,log_amesh,
     >                           phi_ae,phi_ps,
     >                           phi_ae_prime,phi_ps_prime,
     >                           r,
     >                           dvxclm_ae,dvxclm_ps,
     >                           tmpC,
     >                           dmatr)
      implicit none
      integer ng,nb,ic,istart,lmax2,ispin
      real*8  log_amesh

      real*8  phi_ae(ng,nb)
      real*8  phi_ps(ng,nb)
      real*8  phi_ae_prime(ng,nb)
      real*8  phi_ps_prime(ng,nb)
      real*8  r(ic)
      real*8  dvxclm_ae(ic,lmax2,3,ispin)
      real*8  dvxclm_ps(ic,lmax2,3,ispin)
      real*8  tmpC(ic)

      real*8   dmatr(nb,nb,lmax2,3,ispin)

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


             tmpC(ig) = dvxclm_ae(ig,lm,1,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,1,ms)*tmp_ps

          end do
          dmatr(i,j,lm,1,ms)=log_integrate_def(0,tmpC,2,r,log_amesh,ic)
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
             tmpC(ig) = dvxclm_ae(ig,lm,2,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,2,ms)*tmp_ps
          end do
          dmatr(i,j,lm,2,ms) =
     >         log_integrate_def(0,tmpC,2,r,log_amesh,ic,istart)
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
             tmpC(ig) = dvxclm_ae(ig,lm,3,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,3,ms)*tmp_ps
          end do
          dmatr(i,j,lm,3,ms)
     >      = log_integrate_def(0,tmpC,2,r,log_amesh,ic)
          if (i.ne.j) dmatr(j,i,lm,3,ms) = dmatr(i,j,lm,3,ms)
        end do
        end do
      end do
      end do

      end if

      return
      end


*     *****************************************************
*     *                                                   *
*     *          nwpw_get_spher_grid                      *
*     *                                                   *
*     *****************************************************
      subroutine nwpw_get_spher_grid(ntheta,nphi,angle_phi,
     >                          cos_theta,w_theta,w_phi)
      implicit none

      integer ntheta
      integer nphi 
      double precision cos_theta(ntheta)
      double precision angle_phi(nphi)
      double precision w_theta(ntheta)
      double precision w_phi(nphi)

*     *** local variables ***
      integer i
      real*8 pi

      pi = 4.0d0*datan(1.0d0)

*     *** gaussian quadrature angular grid for cos_theta ***
      call nwpw_gauss_weights(-1.0d0,1.0d0,cos_theta,w_theta,ntheta)

      if (nphi.gt.1) then
*       *** linear angular grid for angle_phi***
        do i=1,nphi
         angle_phi(i) = 2.0d0*pi*(i-1)/dble(nphi-1)
         w_phi(i) = 2.0d0*pi/dble(nphi-1)
        end do
        w_phi(1)    = 0.5d0*w_phi(1)
        w_phi(nphi) = w_phi(1)
      else
        angle_phi(1) = 0.0d0
        w_phi(1)     = 2.0d0*pi
      end if
      return
      end  

*     *****************************************************
*     *                                                   *
*     *          nwpw_gauss_weights                       *
*     *                                                   *
*     *****************************************************
      subroutine nwpw_gauss_weights(x1, x2, x, w, n)
      implicit none
      integer n
      double precision x1, x2
      double precision x(*), w(*)

!    *** local variables ***
      integer i, j, m, niter
      double precision eps
      parameter (eps = 3.0d-14)
      double precision p1, p2, p3, pp, xl, xm, z, z1,pi

      pi = 4.0d0*datan(1.0d0)
      m = (n + 1)/2  
      xm = 0.5d0*(x2 + x1)
      xl = 0.5d0*(x2 - x1)

      do i = 1, m
         z = cos(pi*(i - 0.25d0)/(n + 0.5d0))
         niter = 0
    1    continue

         niter = niter + 1
         if (niter .ge. 1000000)
     >     call errquit('cannot converge in gauss_weights',0,1)

         p1 = 1.0d0
         p2 = 0.0d0
         do j = 1, n
            p3 = p2
            p2 = p1
            p1 = ((2.0d0*j - 1.0d0)*z*p2 - (j - 1.0d0)*p3)/j
         end do

         pp = n*(z*p1 - p2)/(z*z - 1.0d0)
         z1 = z
         z = z1 - p1/pp

         if (abs(z - z1) .gt. eps) go to 1
         x(i) = xm - xl*z
         x(n+1-i) = xm + xl*z
         w(i) = 2.0d0*xl/((1.0d0 - z*z)*pp*pp)
         w(n+1-i) = w(i)
      end do
      return
      end

