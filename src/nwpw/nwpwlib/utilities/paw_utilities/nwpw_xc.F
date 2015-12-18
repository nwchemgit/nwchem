c $Id$

*     *********************************
*     *                               *
*     *         nwpw_xc_init          *
*     *                               *
*     *********************************
      subroutine nwpw_xc_init(nion0,nkatm0,
     >                        nprj,nbasis,n1dgrid,psp_type,lmax0,
     >                        nprj_max0,l_prj,m_prj,b_prj)
      implicit none
      integer nion0,nkatm0
      integer nprj(*),nbasis(*),n1dgrid(*),psp_type(*),lmax0(*)
      integer nprj_max0
      integer l_prj(nprj_max0,*),m_prj(nprj_max0,*),b_prj(nprj_max0,*)


#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

c     **** local variables ****
      logical ok
      integer ii,ia,nsize
      integer ic
      integer l,m,li,mi,lj,mj,bi,bj,iprj,jprj,lm
      integer i_p,i_t
      integer i_tlm,i_plm,i_rlm
      integer mtr_size
      integer nb
      real*8  tmp_theta,cs_theta,tmp_phi,tmp_gaunt,tmp_gaunt2,tmp_gaunt3
      real*8  angle_phi
      integer paw_vxc_size

c     **** external functions ****
      integer  control_gga,control_lmax_multipole
      external control_gga,control_lmax_multipole
      real*8   rtheta_lm,drtheta_lm,nwpw_gaunt,nwpw_gaunt2,nwpw_gaunt3
      real*8   rtheta_lm_div
      external rtheta_lm,drtheta_lm,nwpw_gaunt,nwpw_gaunt2,nwpw_gaunt3
      external rtheta_lm_div


      call nwpw_timing_start(4)
      ok = .true.

      nion  = nion0
      nkatm = nkatm0
      gga   = control_gga()

      nprj_max    = 0
      nbasis_max  = 0
      n1dgrid_max = 0
      do ia=1,nkatm
         if (nprj(ia)   .gt.nprj_max)    nprj_max = nprj(ia)
         if (nbasis(ia) .gt.nbasis_max)  nbasis_max = nbasis(ia)
         if (n1dgrid(ia).gt.n1dgrid_max) n1dgrid_max = n1dgrid(ia)
      end do

      mult_l_max = control_lmax_multipole()
      if (mult_l_max.lt.0) then
         mult_l_max = 0
         do ia=1,nkatm
            if (psp_type(ia).eq.4) then
               if (mult_l_max.lt.(2*lmax0(ia))) mult_l_max = 2*lmax0(ia)
            end if
         end do
      end if
      !call nwpw_gaunt2_init(.false.,mult_l_max)
      !call nwpw_gaunt3_init(.false.,mult_l_max)
      
*     ***paw_xc energies ***
      ok = BA_alloc_get(mt_dbl,nion,"paw_xc_e",paw_xc_e(2),paw_xc_e(1))
      call dcopy(nion,0.0d0,0,dbl_mb(paw_xc_e(1)),1)

      
*     *** xc matrix arrays - nbasis_max**2 * (mult_l_max+1)**2 ***
      mtr_size = 2*(nbasis_max**2)*((mult_l_max+1)**2)
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,mtr_size,"paw_xc_matr",
     >                  paw_xc_matr(2),paw_xc_matr(1))

      if (gga.ge.10) then
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,3*mtr_size,"paw_xc_dmatr_u",
     >                  paw_xc_dmatr(2),paw_xc_dmatr(1))
      end if
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_matr",
     >   mtr_size,0)
 

      mtr_size        = 2*nprj_max*nprj_max
      ok = BA_alloc_get(mt_dbl,mtr_size,"paw_xc_pot",
     >                  paw_xc_pot(2),paw_xc_pot(1)) 
      if (.not.ok)
     > call errquit("init_paw_vxc: error allocating paw_xc_pot",
     >   mtr_size,0)



*     *** spherical grid arrays ***
      if(mult_l_max .eq. 0 ) then
        paw_xc_nphi   = 1
        paw_xc_ntheta = 1
      else
        paw_xc_nphi   = 3*mult_l_max
        paw_xc_ntheta = 3*mult_l_max                      
      end if
      

      ok = BA_alloc_get(mt_dbl,paw_xc_nphi,"paw_xc_angle_phi",
     >                  paw_xc_angle_phi(2),paw_xc_angle_phi(1))

      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,paw_xc_ntheta,"paw_xc_cos_theta",
     >                  paw_xc_cos_theta(2),paw_xc_cos_theta(1))

      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,paw_xc_nphi,"paw_xc_w_phi",
     >                  paw_xc_w_phi(2),paw_xc_w_phi(1))

      ok = ok.and.
     >     BA_alloc_get(mt_dbl,paw_xc_ntheta,"paw_xc_w_theta",
     >                  paw_xc_w_theta(2),paw_xc_w_theta(1))

      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(mult_l_max+1)**2,
     >              "paw_xc_tlm",
     >              paw_xc_tlm(2),paw_xc_tlm(1))


*     **** used for generating derivatives of tlm's ****
      if (gga.ge.10) then

c     **** derivatives wrt to theta ****
      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(mult_l_max+1)**2,
     >              "paw_xc_dtlm_theta",
     >              paw_xc_dtlm_theta(2),paw_xc_dtlm_theta(1))

c     **** derivatives wrt to phi ****
      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,
     >              paw_xc_ntheta*paw_xc_nphi*(mult_l_max+1)**2,
     >              "paw_xc_dtlm_phi",
     >              paw_xc_dtlm_phi(2),paw_xc_dtlm_phi(1))

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
         do l=0,mult_l_max
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
           dbl_mb(paw_xc_tlm(1)+i_tlm) = tmp_theta*tmp_phi
           i_tlm = i_tlm + 1                          
         end do
         end do
      end do
      end do

      if (gga.ge.10) then

c       **** define derivative wrt to theta ****
        i_tlm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,mult_l_max
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
              dbl_mb(paw_xc_dtlm_theta(1)+i_tlm) = tmp_theta*tmp_phi
              i_tlm = i_tlm + 1 
           end do
           end do
        end do
        end do

c       **** define derivative wrt to phi ****
        i_tlm = 0
        do i_t=1,paw_xc_ntheta
        do i_p=1,paw_xc_nphi
           do l=0,mult_l_max
           do m=-l,l
              if (m.eq.0) then
                 dbl_mb(paw_xc_dtlm_phi(1)+ i_tlm) = 0.0d0
              else
                 tmp_theta = rtheta_lm_div(l,m,
     >                              dbl_mb(paw_xc_cos_theta(1)+i_t-1))

                 angle_phi=dbl_mb(paw_xc_angle_phi(1)+i_p-1)
                 if (m.lt.0) then
                    tmp_phi = abs(m)*dcos(abs(m)*angle_phi)
                 else if (m.gt.0) then
                    tmp_phi = -abs(m)*dsin(abs(m)*angle_phi)
                 else
                    tmp_phi = 0.0d0
                 end if
                 dbl_mb(paw_xc_dtlm_phi(1)+ i_tlm) = tmp_theta*tmp_phi
              end if
              i_tlm = i_tlm + 1
           end do
           end do
        end do
        end do

      end if

*     *** temp arrays ***
      nsize = 2*n1dgrid_max
      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,nsize,"rho_ae",rho_ae(2),rho_ae(1))
      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,nsize,"rho_ps",rho_ps(2),rho_ps(1))
      
*     **** allocate gradient's and agr's ****
      if (gga.ge.10) then
       nsize = 2*3*n1dgrid_max
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,
     >              "paw_tmp_drho_ae",
     >              rho_ae_prime(2),rho_ae_prime(1))
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,
     >                   "paw_tmp_drho_ps",
     >                   rho_ps_prime(2),rho_ps_prime(1))
       nsize = (2*2-1)*n1dgrid_max
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,"agr_ae",agr_ae(2),agr_ae(1))
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,"agr_ps",agr_ps(2),agr_ps(1))
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,"fdn_ae",fdn_ae(2),fdn_ae(1))
       ok = ok.and.
     >      BA_alloc_get(mt_dbl,nsize,"fdn_ps",fdn_ps(2),fdn_ps(1))
      end if

      nsize = 2*n1dgrid_max
      ok = ok.and. 
     >     BA_alloc_get(mt_dbl,nsize,"vxc_ae",vxc_ae(2),vxc_ae(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"vxc_ps",vxc_ps(2),vxc_ps(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"exc_ae",exc_ae(2),exc_ae(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"exc_ps",exc_ps(2),exc_ps(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"xc_temp",
     >                   xc_temp(2),xc_temp(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"xc_temp2",
     >                  xc_temp2(2),xc_temp2(1))


*     *** allocate vxclm multipole expansion  arrays ****
      nsize = 2*n1dgrid_max*(mult_l_max+1)**2

      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"paw_vxc_ae",
     >                  paw_vxc_ae(2),paw_vxc_ae(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"paw_vxc_ps",
     >                  paw_vxc_ps(2),paw_vxc_ps(1))

*     *** allocate dvxclm multipole expansion  arrays ****
      if (gga.ge.10) then
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,3*nsize,
     >                  "paw_dvxc_u_ae",
     >                  paw_dvxc_ae(2),paw_dvxc_ae(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,3*nsize,
     >                  "paw_dvxc_u_ps",
     >                  paw_dvxc_ps(2),paw_dvxc_ps(1))
      end if

*     *** allocate rholm multipole expansion  arrays ****
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"paw_rho2_ae",
     >                  paw_rho2_ae(2),paw_rho2_ae(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"paw_rho2_ps",
     >                  paw_rho2_ps(2),paw_rho2_ps(1))

*     *** allocate rholm_prime multipole expansion arrays - need for GGA's****
      if (gga.ge.10) then
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,
     >                  "paw_rho2_ae_prime",
     >                   paw_rho2_ae_prime(2),paw_rho2_ae_prime(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,
     >                  "paw_rho2_ps_prime",
     >                  paw_rho2_ps_prime(2),paw_rho2_ps_prime(1))
      end if
      if (.not.ok)
     > call errquit("nwpw_xc_init: error allocating work arrays",
     >               1,MA_ERR)
      
c      call nwpw_init_gntxc(mult_l_max)
c      if ((gga.ge.10).and.(mult_l_max.ge.1)) then
c        call nwpw_init_gntxc2(mult_l_max)
c        call nwpw_init_gntxc3(mult_l_max)
c      end if


*      ***** set up index arrays for nwpw_xc_density_solve2 *****
      ok = BA_alloc_get(mt_int,nkatm,"nindx_rholm",
     >                  nindx_rholm(2),nindx_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nkatm,"shift_rholm",
     >                  shift_rholm(2),shift_rholm(1))
      if (.not.ok)
     > call errquit("nwpw_xc_init: error allocating work arrays",2,0)

      nsize = 0
      do ia=1,nkatm
         int_mb(shift_rholm(1)+ia-1) = nsize
         do jprj = 1,nprj(ia)
            lj = l_prj(jprj,ia)
            mj = m_prj(jprj,ia)
            do iprj = 1,nprj(ia)
               li = l_prj(iprj,ia)
               mi = m_prj(iprj,ia)
               do l=0,mult_l_max
                  do m=-l,l
                      if ((l.le.(li+lj)) .and. (l.ge.abs(li-lj))) then
                        tmp_gaunt = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
                        if (dabs(tmp_gaunt).gt.1.0d-11) then
                           nsize = nsize + 1
                        end if
                      end if
                  end do
               end do

            end do
         end do
        int_mb(nindx_rholm(1)+ia-1) = nsize-int_mb(shift_rholm(1)+ia-1)
      end do

      ok = BA_alloc_get(mt_int,nsize,"bi_rholm",bi_rholm(2),bi_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"bj_rholm",bj_rholm(2),bj_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"lm_rholm",lm_rholm(2),lm_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"iprj_rholm",
     >                 iprj_rholm(2),iprj_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_int,nsize,"jprj_rholm",
     >                 jprj_rholm(2),jprj_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"coeff_rholm",
     >                 coeff_rholm(2),coeff_rholm(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"coeff_rholm2",
     >                 coeff_rholm2(2),coeff_rholm2(1))
      ok = ok.and.
     >     BA_alloc_get(mt_dbl,nsize,"coeff_rholm3",
     >                 coeff_rholm3(2),coeff_rholm3(1))
      if (.not.ok)
     > call errquit("nwpw_xc_init: error allocating work arrays",3,0)

      nsize = 0
      do ia=1,nkatm
         do jprj = 1,nprj(ia)
            lj = l_prj(jprj,ia)
            mj = m_prj(jprj,ia)
            bj = b_prj(jprj,ia)
            do iprj = 1,nprj(ia)
               li = l_prj(iprj,ia)
               mi = m_prj(iprj,ia)
               bi = b_prj(iprj,ia)
               lm = 0
               do l=0,mult_l_max
                  do m=-l,l
                     if ((l.le.(li+lj)) .and. (l.ge.abs(li-lj))) then
                        tmp_gaunt  = nwpw_gaunt(.false.,l,m,li,mi,lj,mj)
                        if ((gga.ge.10).and.(mult_l_max.gt.0)) then
                         tmp_gaunt2=nwpw_gaunt2(.false.,l,m,li,mi,lj,mj)
                         tmp_gaunt3=nwpw_gaunt3(.false.,l,m,li,mi,lj,mj)
                        else
                         tmp_gaunt2 = 0.0d0
                         tmp_gaunt3 = 0.0d0
                        end if
                        if ((dabs(tmp_gaunt) .gt.1.0d-11).or.
     >                      (dabs(tmp_gaunt2).gt.1.0d-11).or.
     >                      (dabs(tmp_gaunt3).gt.1.0d-11)) then
                           dbl_mb(coeff_rholm(1) +nsize) = tmp_gaunt
                           dbl_mb(coeff_rholm2(1)+nsize) = tmp_gaunt2
                           dbl_mb(coeff_rholm3(1)+nsize) = tmp_gaunt3
                           int_mb(lm_rholm(1)+nsize)   = lm
                           int_mb(iprj_rholm(1)+nsize) = iprj
                           int_mb(jprj_rholm(1)+nsize) = jprj
                           int_mb(bi_rholm(1)+nsize)   = bi
                           int_mb(bj_rholm(1)+nsize)   = bj
                           nsize = nsize + 1
                        end if
                     end if
                     lm = lm + 1
                  end do
               end do
            end do
         end do
      end do


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
     >               rgrid,log_amesh,
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

      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      logical value
      integer i,j,np,i_tlm,i_tlm1,lmax2,nsize,shift,ig
      integer i_t,i_p,lm
      real*8 exc_tmp

*     **** external functions ****
      real*8   log_integrate_def_destroyf
      external log_integrate_def_destroyf

      call nwpw_timing_start(36)

*     *** index for spin down in temp arrays ***      
      lmax2  = (mult_l_max+1)**2
!$OMP MASTER
      dbl_mb(paw_xc_e(1)+ii-1)=0.0d0
!$OMP END MASTER
     
*     *** zero out vxclm multipole arrays ***
      nsize = ispin*n1dgrid_max*lmax2
      !call dcopy(nsize,0.0d0,0,dbl_mb(paw_vxc_ae(1)),1)
      !call dcopy(nsize,0.0d0,0,dbl_mb(paw_vxc_ps(1)),1)
      call Parallel_shared_vector_zero(.false.,nsize,
     >                                 dbl_mb(paw_vxc_ae(1)))
      call Parallel_shared_vector_zero(.true.,nsize,
     >                                 dbl_mb(paw_vxc_ps(1)))

*     *** zero out dvxclm multipole arrays ***
      if (gga.ge.10) then
         !call dcopy(3*nsize,0.0d0,0,dbl_mb(paw_dvxc_ae(1)),1)
         !call dcopy(3*nsize,0.0d0,0,dbl_mb(paw_dvxc_ps(1)),1)
         call Parallel_shared_vector_zero(.false.,3*nsize,
     >                                    dbl_mb(paw_dvxc_ae(1)))
         call Parallel_shared_vector_zero(.true.,3*nsize,
     >                                    dbl_mb(paw_dvxc_ps(1)))
      end if

*     *** find rholm multipole expansion on spherical grid ***
      shift = int_mb(shift_rholm(1)+ia-1)
      call nwpw_xc_density_solve2(ispin,ne,nprj,sw1,
     >                            n1dgrid,nbasis,lmax2,
     >                            int_mb(nindx_rholm(1)+ia-1),
     >                            int_mb(lm_rholm(1)+shift),
     >                            int_mb(iprj_rholm(1)+shift),
     >                            int_mb(jprj_rholm(1)+shift),
     >                            int_mb(bi_rholm(1)+shift),
     >                            int_mb(bj_rholm(1)+shift),
     >                            dbl_mb(coeff_rholm(1)+shift),
     >                            phi_ae,phi_ps,rgrid,
     >                            dbl_mb(paw_rho2_ae(1)),
     >                            dbl_mb(paw_rho2_ps(1)))

*     *** find rholm_prime multipole expansion on spherical grid ***
      if (gga.ge.10) then
          shift = int_mb(shift_rholm(1)+ia-1)
          call nwpw_xc_density_prime_solve2(ispin,ne,nprj,sw1,
     >                                  n1dgrid,nbasis,lmax2,
     >                                  int_mb(nindx_rholm(1)+ia-1),
     >                                  int_mb(lm_rholm(1)+shift),
     >                                  int_mb(iprj_rholm(1)+shift),
     >                                  int_mb(jprj_rholm(1)+shift),
     >                                  int_mb(bi_rholm(1)+shift),
     >                                  int_mb(bj_rholm(1)+shift),
     >                                  dbl_mb(coeff_rholm(1)+shift),
     >                                  phi_ae,phi_ps,
     >                                  phi_ae_prime,phi_ps_prime,rgrid,
     >                                  dbl_mb(paw_rho2_ae_prime(1)),
     >                                  dbl_mb(paw_rho2_ps_prime(1)))
      end if



      i_tlm = 0
      i_tlm1 = 0
      do i_t = 1,paw_xc_ntheta
      do i_p = 1,paw_xc_nphi
      
*       *** zero out temp arrays ***
        nsize = ispin*n1dgrid_max
        !call dcopy(nsize,0.0d0,0,dbl_mb(rho_ae(1)),1)
        !call dcopy(nsize,0.0d0,0,dbl_mb(rho_ps(1)),1)
        !call dcopy(nsize,0.0d0,0,dbl_mb(vxc_ae(1)),1)
        !call dcopy(nsize,0.0d0,0,dbl_mb(vxc_ps(1)),1)
        !call dcopy(nsize,0.0d0,0,dbl_mb(exc_ae(1)),1)
        !call dcopy(nsize,0.0d0,0,dbl_mb(exc_ps(1)),1)
        call Parallel_shared_vector_zero(.false.,nsize,
     >                                   dbl_mb(rho_ae(1)))
        call Parallel_shared_vector_zero(.false.,nsize,
     >                                   dbl_mb(rho_ps(1)))
        call Parallel_shared_vector_zero(.false.,nsize,
     >                                   dbl_mb(vxc_ae(1)))
        call Parallel_shared_vector_zero(.false.,nsize,
     >                                   dbl_mb(vxc_ps(1)))
        call Parallel_shared_vector_zero(.false.,nsize,
     >                                   dbl_mb(exc_ae(1)))
        call Parallel_shared_vector_zero(.true., nsize,
     >                                   dbl_mb(exc_ps(1)))

*       *** generate atomic densities on spherical grid ***          
        call nwpw_xc_gen_atomic_densities(n1dgrid,lmax2,ispin,
     >                                    dbl_mb(paw_rho2_ae(1)),
     >                                    dbl_mb(paw_xc_tlm(1)+i_tlm),
     >                                    core_ae,
     >                                    dbl_mb(rho_ae(1)))
        call nwpw_xc_gen_atomic_densities(n1dgrid,lmax2,ispin,
     >                                    dbl_mb(paw_rho2_ps(1)),
     >                                    dbl_mb(paw_xc_tlm(1)+i_tlm),
     >                                    core_ps,
     >                                    dbl_mb(rho_ps(1)))

*       ***************************************************
*       *** find exchange-correlation on spherical grid ***
*       ***************************************************

*       **** LDA functionals on spherical grid ****
        if (gga.lt.10) then
           
*          **** LDA functionals of ae and ps atomic densities ****
           call nwpw_vosko(n1dgrid,n1dgrid,ispin,
     >              dbl_mb(rho_ae(1)),        
     >              dbl_mb(exc_ae(1)),        
     >              dbl_mb(vxc_ae(1)),        
     >              dbl_mb(xc_temp(1)))       
           call nwpw_vosko(n1dgrid,n1dgrid,ispin,
     >              dbl_mb(rho_ps(1)),        
     >              dbl_mb(exc_ps(1)),        
     >              dbl_mb(vxc_ps(1)),        
     >              dbl_mb(xc_temp(1)))       


*       **** GGA functionals on spherical grid ****
        else

*       *** generate the gradient of atomic densities in spherical ****
*       **** coordinates  on spherical grid                        ****          
          call nwpw_xc_gen_atomic_gradients(n1dgrid,lmax2,ispin,
     >                 dbl_mb(paw_rho2_ae(1)),
     >                 dbl_mb(paw_rho2_ae_prime(1)),
     >                 dbl_mb(paw_xc_tlm(1)+i_tlm),
     >                 dbl_mb(paw_xc_dtlm_theta(1)+i_tlm),
     >                 dbl_mb(paw_xc_dtlm_phi(1)+i_tlm),
     >                 rgrid,core_ae,core_ae_prime,
     >                 dbl_mb(rho_ae_prime(1)))

          call nwpw_xc_gen_atomic_gradients(n1dgrid,lmax2,ispin,
     >                 dbl_mb(paw_rho2_ps(1)),
     >                 dbl_mb(paw_rho2_ps_prime(1)),
     >                 dbl_mb(paw_xc_tlm(1)+i_tlm),
     >                 dbl_mb(paw_xc_dtlm_theta(1)+i_tlm),
     >                 dbl_mb(paw_xc_dtlm_phi(1)+i_tlm),
     >                 rgrid,core_ps,core_ps_prime,
     >                 dbl_mb(rho_ps_prime(1)))

*         **** generate agr ****
          call nwpw_xc_gen_atomic_agr(n1dgrid,lmax2,ispin,
     >                               dbl_mb(rho_ae_prime(1)),
     >                               dbl_mb(agr_ae(1)))
          call nwpw_xc_gen_atomic_agr(n1dgrid,lmax2,ispin,
     >                               dbl_mb(rho_ps_prime(1)),
     >                               dbl_mb(agr_ps(1)))

*          **** GGA functionals of ae and ps atomic densities ****
           call nwpw_gga(gga,n1dgrid,ispin,
     >                    dbl_mb(rho_ae(1)),        
     >                    dbl_mb(agr_ae(1)),        
     >                    dbl_mb(exc_ae(1)),
     >                    dbl_mb(vxc_ae(1)),    !* df/dnup, df/dndn
     >                    dbl_mb(fdn_ae(1)),    !* df/d|grad nup|,df/d|grad ndn|, df/d|grad n|
     >                    dbl_mb(xc_temp(1)))       
           call nwpw_gga(gga,n1dgrid,ispin,
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
           call nwpw_xc_gen_dvxc(n1dgrid,lmax2,ispin,
     >                          dbl_mb(fdn_ae(1)),
     >                          dbl_mb(agr_ae(1)),
     >                          dbl_mb(rho_ae_prime(1)))
           call nwpw_xc_gen_dvxc(n1dgrid,lmax2,ispin,
     >                          dbl_mb(fdn_ps(1)),
     >                          dbl_mb(agr_ps(1)),
     >                          dbl_mb(rho_ps_prime(1)))

c          **** find dvxclm multipole expansion ****
           call nwpw_xc_addto_dvxclm(n1dgrid,lmax2,ispin,
     >                          (dbl_mb(paw_xc_w_theta(1)+i_t-1)
     >                          *dbl_mb(paw_xc_w_phi(1)  +i_p-1)),
     >                          dbl_mb(paw_xc_tlm(1)+i_tlm),
     >                          dbl_mb(rho_ae_prime(1)),
     >                          dbl_mb(rho_ps_prime(1)),
     >                          dbl_mb(paw_dvxc_ae(1)),
     >                          dbl_mb(paw_dvxc_ps(1)))

        end if
        i_tlm = i_tlm + lmax2


*       ************************************************************************************
*       *** find vxclm the multipole expansion of atomic vxc=df/dn  or (df/dnup,df/dndn) ***          
*       ************************************************************************************
         call nwpw_xc_addto_vxclm(n1dgrid,lmax2,ispin,
     >                          (dbl_mb(paw_xc_w_theta(1)+i_t-1)
     >                          *dbl_mb(paw_xc_w_phi(1)  +i_p-1)),
     >                          dbl_mb(paw_xc_tlm(1)+i_tlm1),
     >                          dbl_mb(vxc_ae(1)),
     >                          dbl_mb(vxc_ps(1)),
     >                          dbl_mb(paw_vxc_ae(1)),
     >                          dbl_mb(paw_vxc_ps(1)))
         i_tlm1 = i_tlm1 + lmax2

*        ******************************************************
*        *** compute the atomic exchange correlation energy ***
*        ******************************************************
!$OMP DO
         do ig=1,n1dgrid
            dbl_mb(xc_temp(1)+ig-1)=
     >             (
     >             dbl_mb(rho_ae(1)+ig-1)  + 
     >             dbl_mb(rho_ae(1)+(ispin-1)*n1dgrid+ig-1)
     >             )*dbl_mb(exc_ae(1)+ig-1)
     >             -
     >             (
     >             dbl_mb(rho_ps(1)+ig-1)+ 
     >             dbl_mb(rho_ps(1)+(ispin-1)*n1dgrid+ig-1)
     >             )*
     >             dbl_mb(exc_ps(1)+ig-1)
         end do
!$OMP END DO
!$OMP MASTER
         exc_tmp = log_integrate_def_destroyf(0,dbl_mb(xc_temp(1)),
     >                               2,rgrid,log_amesh,n1dgrid)
         dbl_mb(paw_xc_e(1)+ii-1)=dbl_mb(paw_xc_e(1)+ii-1)+
     >                           exc_tmp*
     >                           dbl_mb(paw_xc_w_theta(1)+i_t-1)*
     >                           dbl_mb(paw_xc_w_phi(1)+i_p-1)
!$OMP END MASTER
!$OMP BARRIER

      end do !i_p
      end do !i_t


*       *********************************************************
*       ****    non-local operator computation - LDA part    ****
*       *********************************************************

*       **** compute (vxc^a)_nln'l'^lm radial integrals *****
        call nwpw_xc_gen_matr(n1dgrid,nbasis,lmax2,ispin,
     >                        log_amesh,
     >                        phi_ae,phi_ps,rgrid,
     >                        dbl_mb(paw_vxc_ae(1)),
     >                        dbl_mb(paw_vxc_ps(1)),
     >                        dbl_mb(xc_temp(1)),
     >                        dbl_mb(paw_xc_matr(1)))


*       **** xc potential  non-local matrix elements ****
        shift = int_mb(shift_rholm(1)+ia-1)
        call nwpw_xc_sw1tosw2(ispin,ne,nprj,sw1,sw2,
     >                        n1dgrid,nbasis,lmax2,
     >                        int_mb(nindx_rholm(1)+ia-1),
     >                        int_mb(lm_rholm(1)+shift),
     >                        int_mb(iprj_rholm(1)+shift),
     >                        int_mb(jprj_rholm(1)+shift),
     >                        int_mb(bi_rholm(1)+shift),
     >                        int_mb(bj_rholm(1)+shift),
     >                         dbl_mb(coeff_rholm(1)+shift),
     >                         dbl_mb(paw_xc_matr(1)))


*       *********************************************************
*       ****    non-local operator computation - GGA part    ****
*       *********************************************************
        if (gga.ge.10) then

*         **** compute (dvxc^a)_nln'l'^lm radial integrals *****
          call nwpw_xc_gen_dmatr(n1dgrid,nbasis,lmax2,ispin,
     >                       log_amesh,
     >                       phi_ae,phi_ps,
     >                       phi_ae_prime,phi_ps_prime,
     >                       rgrid,
     >                       dbl_mb(paw_dvxc_ae(1)),
     >                       dbl_mb(paw_dvxc_ps(1)),
     >                       dbl_mb(xc_temp(1)),
     >                       dbl_mb(paw_xc_dmatr(1)))

          shift = int_mb(shift_rholm(1)+ia-1)
          call nwpw_xc_sw1tosw2_dmatr(ispin,ne,nprj,sw1,sw2,
     >                        n1dgrid,nbasis,lmax2,
     >                        int_mb(nindx_rholm(1)+ia-1),
     >                        int_mb(lm_rholm(1)+shift),
     >                        int_mb(iprj_rholm(1)+shift),
     >                        int_mb(jprj_rholm(1)+shift),
     >                        int_mb(bi_rholm(1)+shift),
     >                        int_mb(bj_rholm(1)+shift),
     >                        dbl_mb(coeff_rholm(1)+shift),
     >                        dbl_mb(coeff_rholm2(1)+shift),
     >                        dbl_mb(coeff_rholm3(1)+shift),
     >                        dbl_mb(paw_xc_dmatr(1)))
        end if !**gga**

c      if (np.gt.1)  then
c
cc       **** xc non-local matrix elements ****
c        call D3dB_Vector_Sumall(2*paw_xc_pot(3),
c     >                          dbl_mb(paw_xc_pot(1)))
c
cc       **** atomic exchange-correlation energies ****
c        call D3dB_Vector_Sumall(paw_xc_e(3),dbl_mb(paw_xc_e(1)))
c      end if
c     
c
      call nwpw_timing_end(36)

      return
      end 


*     *****************************************
*     *                                       *
*     *             nwpw_xc_end               *
*     *                                       *
*     *****************************************
      subroutine nwpw_xc_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      logical ok
      integer ms,i
      
      call nwpw_timing_start(4)

      !call nwpw_gaunt2_end()
      !call nwpw_gaunt3_end()

      ok = .true.
      ok = ok.and.BA_free_heap(coeff_rholm(2))
      ok = ok.and.BA_free_heap(coeff_rholm2(2))
      ok = ok.and.BA_free_heap(coeff_rholm3(2))
      ok = ok.and.BA_free_heap(jprj_rholm(2))
      ok = ok.and.BA_free_heap(iprj_rholm(2))
      ok = ok.and.BA_free_heap(lm_rholm(2))
      ok = ok.and.BA_free_heap(bj_rholm(2))
      ok = ok.and.BA_free_heap(bi_rholm(2))
      ok = ok.and.BA_free_heap(shift_rholm(2))
      ok = ok.and.BA_free_heap(nindx_rholm(2))


      ok = ok.and.BA_free_heap(exc_ps(2))
      ok = ok.and.BA_free_heap(exc_ae(2))
      ok = ok.and.BA_free_heap(vxc_ps(2))
      ok = ok.and.BA_free_heap(vxc_ae(2))
      ok = ok.and.BA_free_heap(rho_ps(2))
      ok = ok.and.BA_free_heap(rho_ae(2))
      ok = ok.and.BA_free_heap(xc_temp(2))
      ok = ok.and.BA_free_heap(xc_temp2(2))

      ok = ok.and.BA_free_heap(paw_xc_e(2))
      ok = ok.and.BA_free_heap(paw_xc_tlm(2)) 
      ok = ok.and.BA_free_heap(paw_xc_w_theta(2) )
      ok = ok.and.BA_free_heap(paw_xc_w_phi(2))
      ok = ok.and.BA_free_heap(paw_xc_cos_theta(2)) 
      ok = ok.and.BA_free_heap(paw_xc_angle_phi(2))

      ok = ok.and.BA_free_heap(paw_xc_pot(2))
      ok = ok.and.BA_free_heap(paw_xc_matr(2))
      ok = ok.and.BA_free_heap(paw_vxc_ae(2))
      ok = ok.and.BA_free_heap(paw_vxc_ps(2))
      ok = ok.and.BA_free_heap(paw_rho2_ae(2))
      ok = ok.and.BA_free_heap(paw_rho2_ps(2))

      if (gga.ge.10) then
         ok = ok.and.BA_free_heap(paw_rho2_ae_prime(2))
         ok = ok.and.BA_free_heap(paw_rho2_ps_prime(2))
         ok = ok.and.BA_free_heap(rho_ps_prime(2))
         ok = ok.and.BA_free_heap(rho_ae_prime(2))
         ok = ok.and.BA_free_heap(agr_ae(2))
         ok = ok.and.BA_free_heap(agr_ps(2))
         ok = ok.and.BA_free_heap(fdn_ae(2))
         ok = ok.and.BA_free_heap(fdn_ps(2))
         ok = ok.and.BA_free_heap(paw_xc_dtlm_phi(2))
         ok = ok.and.BA_free_heap(paw_xc_dtlm_theta(2))
         ok = ok.and.BA_free_heap(paw_dvxc_ae(2))
         ok = ok.and.BA_free_heap(paw_dvxc_ps(2))
         ok = ok.and.BA_free_heap(paw_xc_dmatr(2))
      end if


      if (.not.ok)
     > call errquit("nwpw_xc_end: error freeing heap",0,MA_ERR)


      call nwpw_timing_end(4)
      return
      end 

*      ******************************************
*      *                                        *
*      *          nwpw_xc_mult_l_max            *
*      *                                        *
*      ******************************************
      integer function nwpw_xc_mult_l_max()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      nwpw_xc_mult_l_max = mult_l_max
      return
      end



*      ******************************************
*      *                                        *
*      *          nwpw_xc_energy_atom           *
*      *                                        *
*      ******************************************
      real*8 function nwpw_xc_energy_atom()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_xc.fh"

      integer ii
      real*8 e

      e = 0.0d0
      do ii=1,nion
         e = e + dbl_mb(paw_xc_e(1)+ii-1)
      end do
      nwpw_xc_energy_atom = e
      return
      end



*     ******************************************
*     *                                        *
*     *          nwpw_xc_gen_sphere_rho        *
*     *                                        *
*     ******************************************
      subroutine nwpw_xc_gen_sphere_rho(ng,lmax2,rholm,Tlm,rho)
      implicit none
      integer ng,lmax2
      real*8  rholm(ng,lmax2)
      real*8  Tlm(lmax2)
      real*8  rho(ng)

      integer    lm,ig

!$OMP DO
      do ig=1,ng
         do lm=1,lmax2
           rho(ig) = rho(ig) + rholm(ig,lm)*Tlm(lm)
         end do
      end do
!$OMP END DO
      
      return
      end


*     ******************************************
*     *                                        *
*     *          nwpw_xc_rho_div_r             *
*     *                                        *
*     ******************************************
      subroutine nwpw_xc_rho_div_r(ic,r,rho)
      implicit none
      integer ic
      real*8     r(*)
      real*8     rho(*)

      integer    ig
!$OMP DO
      do ig=1,ic
        rho(ig) = rho(ig)/r(ig)
      end do
!$OMP END DO
      return
      end



*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_addto_vxclm                      *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_addto_vxclm(ic,lmax2,ispin,
     >                              alpha,
     >                              tlm,
     >                              vxc_ae,vxc_ps,
     >                              vxclm_ae,vxclm_ps)
      implicit none
      integer ic,lmax2,ispin
      real*8 alpha
      real*8 tlm(*)
      real*8 vxc_ae(ic,ispin)
      real*8 vxc_ps(ic,ispin)
      real*8 vxclm_ae(ic,lmax2,ispin)
      real*8 vxclm_ps(ic,lmax2,ispin)

*     **** local variables ****
      integer    i,lm,ms

      do ms = 1,ispin
         do lm = 1,lmax2
!$OMP DO
            do i  = 1,ic
               vxclm_ae(i,lm,ms) = vxclm_ae(i,lm,ms) 
     >                           + vxc_ae(i,ms)*(tlm(lm)*alpha)
               vxclm_ps(i,lm,ms) = vxclm_ps(i,lm,ms) 
     >                          + vxc_ps(i,ms)*(tlm(lm)*alpha)
            end do
!$OMP END DO
         end do
      end do
      return
      end

*     ********************************************************
*     *                                                      *
*     *             nwpw_xc_gen_matr                         *
*     *                                                      *
*     ********************************************************
      subroutine nwpw_xc_gen_matr(ng,nb,lmax2,ispin,log_amesh,
     >                           phi_ae,phi_ps,r,
     >                           vxclm_ae,vxclm_ps,
     >                           tmpC,
     >                           matr) 
      implicit none
      integer ng,nb,lmax2,ispin
      real*8  log_amesh
      real*8  phi_ae(ng,nb)
      real*8  phi_ps(ng,nb)
      real*8  r(ng)
      real*8  vxclm_ae(ng,lmax2,ispin)
      real*8  vxclm_ps(ng,lmax2,ispin)
      real*8  tmpC(ng)
      real*8 matr(nb,nb,lmax2,ispin)

*     **** local varialbles ****
      integer ig,i,j,lm,ms
      real*8  tmp_ae,tmp_ps

*     **** external functions ****
      real*8   log_integrate_def_destroyf
      external log_integrate_def_destroyf

      do ms=1,ispin
         do lm=1,lmax2
           do i=1,nb
           do j=i,nb
!$OMP DO
              do ig=1,ng
                 tmp_ae = phi_ae(ig,i)
     >                   *phi_ae(ig,j)/(r(ig)**2)
                 tmp_ps = phi_ps(ig,i)
     >                   *phi_ps(ig,j)/(r(ig)**2)
                 tmpC(ig) = vxclm_ae(ig,lm,ms)*tmp_ae
     >                    - vxclm_ps(ig,lm,ms)*tmp_ps
              end do
!$OMP END DO
!$OMP MASTER
              matr(i,j,lm,ms)
     >         = log_integrate_def_destroyf(0,tmpC,2,r,log_amesh,ng)
              if (i.ne.j) matr(j,i,lm,ms) = matr(i,j,lm,ms)
!$OMP END MASTER
!$OMP BARRIER
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
      subroutine nwpw_xc_gen_atomic_densities(ng,lmax2,ispin,
     >                                 rholm,
     >                                 tlm,
     >                                 rhocore,
     >                                 rho)
      implicit none
      integer ng,lmax2,ispin
      real*8 rholm(ng,lmax2,ispin)
      real*8 tlm(lmax2)
      real*8 rhocore(ng)
      real*8 rho(ng,ispin)

      integer ms,i

      do ms=1,ispin
         call nwpw_xc_gen_sphere_rho(ng,lmax2,
     >                              rholm(1,1,ms),
     >                              tlm,
     >                              rho(1,ms))
         !call daxpy(ng,0.5d0,rhocore,1,rho(1,ms),1)
!$OMP DO
         do i=1,ng
            rho(i,ms) = dabs(rho(i,ms)+0.5d0*rhocore(i))
         end do
!$OMP END DO
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
     >                                       tlm,dtlm_theta,dtlm_phi,
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


      integer ms,i

      !call dcopy(3*ic*ispin,0.0d0,0,rho_prime,1)
      call Parallel_shared_vector_zero(.true.,3*ic*ispin,rho_prime)

      do ms=1,ispin

*        *** find [d/dr paw_rho_ae] on spherical grid ***
         call nwpw_xc_gen_sphere_rho(ic,lmax2,
     >                              rholm_prime(1,1,ms),
     >                              tlm,
     >                              rho_prime(1,1,ms))

*        *** add core d/dr densities***
         !call daxpy(ic,0.5d0,rhocore_prime,1,rho_prime(1,1,ms),1)
!$OMP DO
         do i=1,ic
            rho_prime(i,1,ms) = rho_prime(i,1,ms) 
     >                        + 0.5d0*rhocore_prime(i)
         end do
!$OMP END DO


         !*** only computate radial derivatives if lmax2==1 ****
         if (lmax2.gt.1) then

*          *** find  (1/r) * d/dtheta paw_rho_ae on spherical grid ***
           call nwpw_xc_gen_sphere_rho(ic,lmax2,
     >               rholm(1,1,ms),
     >               dtlm_theta,
     >               rho_prime(1,2,ms))
           call nwpw_xc_rho_div_r(ic,r,rho_prime(1,2,ms))

*          *** find  [(r*sin(theta)) * d/dphi paw_rho_ae] on spherical grid ***
           call nwpw_xc_gen_sphere_rho(ic,lmax2,
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
!$OMP DO
         do ig=1,ic
            agr(ig,2*ispin-1)
     >       = dsqrt( (rho_prime(ig,1,1)+rho_prime(ig,1,ispin))**2)
         end do
!$OMP END DO

c        **** compute |grad nup| and |grad ndn| ****
         if (ispin.eq.2) then
           do ms=1,ispin
!$OMP DO
           do ig=1,ic
              agr(ig,ms) = dsqrt(rho_prime(ig,1,ms)**2)
           end do
!$OMP END DO
           end do
         end if



      !***  computate theta and phi derivatives if lmax2>1 ****
      else
c        **** compute |grad n| ****
!$OMP DO
         do ig=1,ic   
            agr(ig,2*ispin-1) 
     >       = dsqrt( (rho_prime(ig,1,1)+rho_prime(ig,1,ispin))**2
     >              + (rho_prime(ig,2,1)+rho_prime(ig,2,ispin))**2
     >              + (rho_prime(ig,3,1)+rho_prime(ig,3,ispin))**2)
         end do
!$OMP END DO

c        **** compute |grad nup| and |grad ndn| ****
         if (ispin.eq.2) then
           do ms=1,ispin
!$OMP DO
           do ig=1,ic
              agr(ig,ms) = dsqrt( rho_prime(ig,1,ms)**2
     >                          + rho_prime(ig,2,ms)**2
     >                          + rho_prime(ig,3,ms)**2)
           end do
!$OMP END DO
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
      real*8  eta
      parameter (eta=1.0d-20)

      !*** only computate radial derivatives if lmax2=1 ****
      if (lmax2.eq.1) then
         jmax = 1
      else
         jmax = 3
      end if

*     *** restricted ****
      if (ispin.eq.1) then
         do j=1,jmax
!$OMP DO
          do i=1,ic
           rho_prime(i,j,1) = (rho_prime(i,j,1)+rho_prime(i,j,1))
     >                       *(fdn(i,1)/(agr(i,1)+eta))
          end do
!$OMP END DO
         end do

*     *** unrestricted ****
      else
         do j=1,jmax
!$OMP DO
          do i=1,ic
            drho1 = rho_prime(i,j,1)
            drho2 = rho_prime(i,j,2)
            drhoa = drho1+drho2
            rho_prime(i,j,1) = (drho1/(agr(i,1)+eta))*fdn(i,1)
     >                       + (drhoa/(agr(i,3)+eta))*fdn(i,3)
            rho_prime(i,j,2) = (drho2/(agr(i,2)+eta))*fdn(i,2)
     >                       + (drhoa/(agr(i,3)+eta))*fdn(i,3)
          end do
!$OMP END DO
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
!$OMP DO
               do i  = 1,ic
                 dvxclm_ae(i,lm,j,ms) = dvxclm_ae(i,lm,j,ms)
     >                                + dvxc_ae(i,j,ms)*(tlm(lm)*alpha)
                 dvxclm_ps(i,lm,j,ms) = dvxclm_ps(i,lm,j,ms)
     >                                + dvxc_ps(i,j,ms)*(tlm(lm)*alpha)
               end do
!$OMP END DO
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

      subroutine nwpw_xc_gen_dmatr(ng,nb,lmax2,
     >                           ispin,log_amesh,
     >                           phi_ae,phi_ps,
     >                           phi_ae_prime,phi_ps_prime,
     >                           r,
     >                           dvxclm_ae,dvxclm_ps,
     >                           tmpC,
     >                           dmatr)
      implicit none
      integer ng,nb,lmax2,ispin
      real*8  log_amesh

      real*8  phi_ae(ng,nb)
      real*8  phi_ps(ng,nb)
      real*8  phi_ae_prime(ng,nb)
      real*8  phi_ps_prime(ng,nb)
      real*8  r(ng)
      real*8  dvxclm_ae(ng,lmax2,3,ispin)
      real*8  dvxclm_ps(ng,lmax2,3,ispin)
      real*8  tmpC(ng)

      real*8   dmatr(nb,nb,lmax2,3,ispin)

*     **** local varialbles ****
      integer ig,i,j,lm,ms
      real*8  tmp_ae,tmp_ps

*     **** external functions ****
      real*8   log_integrate_def_destroyf
      external log_integrate_def_destroyf

      !call dcopy(nb*nb*lmax2*3*ispin,0.0d0,0,dmatr,1)
      call Parallel_shared_vector_zero(.true.,nb*nb*lmax2*3*ispin,dmatr)

*     **** radial derivative integral ****
      do ms=1,ispin
      do lm=1,lmax2
        do i=1,nb
        do j=i,nb
!$OMP DO
          do ig=1,ng
            tmp_ae 
     >        = ( phi_ae_prime(ig,i)*phi_ae(ig,j)
     >          + phi_ae(ig,i)*phi_ae_prime(ig,j))
     >          /(r(ig)**2)
     >        - 2.0d0*phi_ae(ig,i)*phi_ae(ig,j)
     >          /(r(ig)**3)
            tmp_ps 
     >        = ( phi_ps_prime(ig,i)*phi_ps(ig,j)
     >          + phi_ps(ig,i)*phi_ps_prime(ig,j))
     >          /(r(ig)**2)
     >        - 2.0d0*phi_ps(ig,i)*phi_ps(ig,j)
     >          /(r(ig)**3)

             tmpC(ig) = dvxclm_ae(ig,lm,1,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,1,ms)*tmp_ps
          end do
!$OMP END DO
!$OMP MASTER
          dmatr(i,j,lm,1,ms)
     >     = log_integrate_def_destroyf(0,tmpC,2,r,log_amesh,ng)
          if (i.ne.j) dmatr(j,i,lm,1,ms) = dmatr(i,j,lm,1,ms)
!$OMP END MASTER
!$OMP BARRIER
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
!$OMP DO
          do ig=1,ng
             tmp_ae = phi_ae(ig,i)
     >               *phi_ae(ig,j)
     >               /(r(ig)**3)
             tmp_ps = phi_ps(ig,i)
     >               *phi_ps(ig-1,j)
     >               /(r(ig)**3)
             tmpC(ig) = dvxclm_ae(ig,lm,2,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,2,ms)*tmp_ps
          end do
!$OMP END DO
!$OMP MASTER
          dmatr(i,j,lm,2,ms)
     >     = log_integrate_def_destroyf(0,tmpC,2,r,log_amesh,ng)
          if (i.ne.j) dmatr(j,i,lm,2,ms) = dmatr(i,j,lm,2,ms)
!$OMP END MASTER
!$OMP BARRIER
        end do
        end do
      end do
      end do

*     **** phi derivative integral ****
      do ms=1,ispin
      do lm=1,lmax2
        do i=1,nb
        do j=i,nb

!$OMP DO
          do ig=1,ng
             tmp_ae = phi_ae(ig,i)
     >               *phi_ae(ig,j)
     >               /(r(ig)**3)
             tmp_ps = phi_ps(ig,i)
     >               *phi_ps(ig,j)
     >               /(r(ig)**3)
             tmpC(ig) = dvxclm_ae(ig,lm,3,ms)*tmp_ae
     >                - dvxclm_ps(ig,lm,3,ms)*tmp_ps
          end do
!$OMP END DO
!$OMP MASTER
          dmatr(i,j,lm,3,ms)
     >     = log_integrate_def_destroyf(0,tmpC,2,r,log_amesh,ng)
          if (i.ne.j) dmatr(j,i,lm,3,ms) = dmatr(i,j,lm,3,ms)
!$OMP END MASTER
!$OMP BARRIER
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



c     *********************************************
c     *                                           *
c     *           nwpw_xc_density_solve2          *
c     *                                           *
c     *********************************************
c
c   Calculates the atomic density lm expansions  from the
c   overlap coefficients.

      subroutine nwpw_xc_density_solve2(ispin,ne,nprj,sw1,
     >                                  n1dgrid,nbasis,lmax2,
     >                                  nindx,lm_rholm,
     >                                  iprj_rholm,jprj_rholm,
     >                                  bi_rholm,bj_rholm,
     >                                  coeff_rholm,
     >                                  phi_ae,phi_ps,
     >                                  rgrid,
     >                                  rholm_ae,rholm_ps)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      integer n1dgrid,nbasis,lmax2
      integer nindx,lm_rholm(*),iprj_rholm(*),jprj_rholm(*)
      integer bi_rholm(*),bj_rholm(*)
      real*8  coeff_rholm(*)
      real*8  phi_ae(n1dgrid,nbasis),phi_ps(n1dgrid,nbasis)
      real*8  rgrid(n1dgrid)
      real*8  rholm_ae(n1dgrid,lmax2,ispin)
      real*8  rholm_ps(n1dgrid,lmax2,ispin)

      integer ms,n,i,lm,iprj,jprj,bi,bj,n1(2),n2(2)
      real*8  tmp_ms,coeff,w,scal

*     **** parallel_wshared common block - used for simple OMP reductions ***
      real*8 wshared(100),wshared1
      common /parallel_wshared/ wshared,wshared1

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega
      
      call nwpw_timing_start(21)
      n1(1) = 1 
      n1(2) = ne(1)+1 
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      scal = 1.0d0/lattice_omega()

*     ***init to zero***
      !call dcopy(ispin*lmax2*n1dgrid,0.0d0,0,rholm_ae,1)
      !call dcopy(ispin*lmax2*n1dgrid,0.0d0,0,rholm_ps,1)
      call Parallel_shared_vector_zero(.false.,ispin*lmax2*n1dgrid,
     >                                 rholm_ae)
      call Parallel_shared_vector_zero(.true.,ispin*lmax2*n1dgrid,
     >                                 rholm_ps)

      do i=1,nindx
         lm    = lm_rholm(i)
         iprj  = iprj_rholm(i)
         jprj  = jprj_rholm(i)
         bi    = bi_rholm(i)
         bj    = bj_rholm(i)
         coeff = coeff_rholm(i)*scal
         do ms=1,ispin
            if (ne(ms).gt.0) then
!$OMP MASTER
               wshared1 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:wshared1)
               do n=n1(ms),n2(ms)
                  wshared1 = wshared1 + sw1(n,iprj)*sw1(n,jprj)
               end do
!$OMP END DO
               call D1dB_SumAll(wshared1)
               tmp_ms = wshared1*coeff
               call nwpw_xc_density_gen_rho(n1dgrid,tmp_ms,
     >                  rgrid,phi_ae(1,bi),phi_ae(1,bj),
     >                  rholm_ae(1,lm+1,ms))
               call nwpw_xc_density_gen_rho(n1dgrid,tmp_ms,
     >                  rgrid,phi_ps(1,bi),phi_ps(1,bj),
     >                  rholm_ps(1,lm+1,ms))
            end if
         end do
      end do


      call nwpw_timing_end(21)
      return
      end


      subroutine nwpw_xc_density_gen_rho(ic,alpha,r,phi1,phi2,rho)
      implicit none
      integer ic
      double precision alpha
      double precision r(ic)
      double precision phi1(ic)
      double precision phi2(ic)
      double precision rho(ic)

*     **** local variables ****
      integer i
      double precision tmp

!$OMP DO
      do i=1,ic
         tmp = (phi1(i)*phi2(i))/(r(i)**2)
         rho(i) = rho(i) + alpha*tmp
      end do
!$OMP END DO

      return
      end


c     *********************************************
c     *                                           *
c     *        nwpw_xc_density_prime_solve2       *
c     *                                           *
c     *********************************************
c
c   Calculates the atomic density lm expansions  from the
c   overlap coefficients.
      subroutine nwpw_xc_density_prime_solve2(ispin,ne,nprj,sw1,
     >                                  n1dgrid,nbasis,lmax2,
     >                                  nindx,lm_rholm,
     >                                  iprj_rholm,jprj_rholm,
     >                                  bi_rholm,bj_rholm,
     >                                  coeff_rholm,
     >                                  phi_ae,phi_ps,
     >                                  dphi_ae,dphi_ps,rgrid,
     >                                  rholm_ae_prime,rholm_ps_prime)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      integer n1dgrid,nbasis,lmax2
      integer nindx,lm_rholm(*),iprj_rholm(*),jprj_rholm(*)
      integer bi_rholm(*),bj_rholm(*)
      real*8  coeff_rholm(*)
      real*8  phi_ae(n1dgrid,nbasis),phi_ps(n1dgrid,nbasis)
      real*8  dphi_ae(n1dgrid,nbasis),dphi_ps(n1dgrid,nbasis)
      real*8  rgrid(n1dgrid)
      real*8  rholm_ae_prime(n1dgrid,lmax2,ispin)
      real*8  rholm_ps_prime(n1dgrid,lmax2,ispin)

      integer ms,n,i,lm,iprj,jprj,bi,bj,n1(2),n2(2)
      real*8  tmp_ms,coeff,w,scal

*     **** parallel_wshared common block - used for simple OMP reductions ***
      real*8 wshared(100),wshared1
      common /parallel_wshared/ wshared,wshared1


*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      scal = 1.0d0/lattice_omega()

*     ***init to zero***
      !call dcopy(n1dgrid*lmax2*ispin,0.0d0,0,rholm_ae_prime,1)
      !call dcopy(n1dgrid*lmax2*ispin,0.0d0,0,rholm_ps_prime,1)
      call Parallel_shared_vector_zero(.false.,n1dgrid*lmax2*ispin,
     >                                 rholm_ae_prime)
      call Parallel_shared_vector_zero(.true.,n1dgrid*lmax2*ispin,
     >                                 rholm_ps_prime)

      do i=1,nindx
         lm    = lm_rholm(i)
         iprj  = iprj_rholm(i)
         jprj  = jprj_rholm(i)
         bi    = bi_rholm(i)
         bj    = bj_rholm(i)
         coeff = coeff_rholm(i)*scal
         do ms=1,ispin
            if (ne(ms).gt.0) then
!$OMP MASTER
               wshared1 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:wshared1)
               do n=n1(ms),n2(ms)
                  wshared1 = wshared1 + sw1(n,iprj)*sw1(n,jprj)
               end do
!$OMP END DO
               call D1dB_SumAll(wshared1)
               tmp_ms = wshared1*coeff
               call nwpw_xc_density_gen_drho(n1dgrid,tmp_ms,
     >                  rgrid,
     >                  phi_ae(1,bi),dphi_ae(1,bi),
     >                  phi_ae(1,bj),dphi_ae(1,bj),
     >                  rholm_ae_prime(1,lm+1,ms))
               call nwpw_xc_density_gen_drho(n1dgrid,tmp_ms,
     >                  rgrid,
     >                  phi_ps(1,bi),dphi_ps(1,bi),
     >                  phi_ps(1,bj),dphi_ps(1,bj),
     >                  rholm_ps_prime(1,lm+1,ms))
            end if
         end do
      end do

      call nwpw_timing_end(21)
      return
      end

      subroutine nwpw_xc_density_gen_drho(ic,alpha,r,
     >                                phi1,dphi1,
     >                                phi2,dphi2,
     >                                drho)
      implicit none
      integer ic
      double precision alpha
      double precision r(ic)
      double precision phi1(ic),dphi1(ic)
      double precision phi2(ic),dphi2(ic)
      double precision drho(ic)

*     **** local variables ****
      integer i
      double precision tmp

!$OMP DO
      do i=1,ic
         tmp = (dphi1(i)*phi2(i)+phi1(i)*dphi2(i))/(r(i)**2)
     >       - 2.0d0*(phi1(i)*phi2(i))/(r(i)**3)
         drho(i) = drho(i) + alpha*tmp
      end do
!$OMP END DO
      return
      end


c     *********************************************
c     *                                           *
c     *           nwpw_xc_sw1tosw2                *
c     *                                           *
c     *********************************************
c
      subroutine nwpw_xc_sw1tosw2(ispin,ne,nprj,sw1,sw2,
     >                            n1dgrid,nbasis,lmax2,
     >                            nindx,lm_rholm,
     >                            iprj_rholm,jprj_rholm,
     >                            bi_rholm,bj_rholm,
     >                            coeff_rholm,
     >                            matr)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer n1dgrid,nbasis,lmax2
      integer nindx,lm_rholm(*),iprj_rholm(*),jprj_rholm(*)
      integer bi_rholm(*),bj_rholm(*)
      real*8  coeff_rholm(*)
      real*8  matr(nbasis,nbasis,lmax2,ispin)

      integer ms,n,i,lm,iprj,jprj,bi,bj,n1(2),n2(2)
      real*8  tmp_ms,coeff,w,scal

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      n1(1) = 1 
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      !scal = 1.0d0/dsqrt(lattice_omega())

*     ***init to zero***
      do i=1,nindx
         lm    = lm_rholm(i)
         iprj  = iprj_rholm(i)
         jprj  = jprj_rholm(i)
         bi    = bi_rholm(i)
         bj    = bj_rholm(i)
         !coeff = coeff_rholm(i)*scal
         coeff = coeff_rholm(i)
         do ms=1,ispin
            if (ne(ms).gt.0) then
!$OMP DO
               do n=n1(ms),n2(ms)
                  sw2(n,iprj) = sw2(n,iprj) 
     >                        + coeff*matr(bi,bj,lm+1,ms)*sw1(n,jprj)
               end do 
!$OMP END DO
            end if
         end do
      end do
      call nwpw_timing_end(21)
      return
      end 



c     *********************************************
c     *                                           *
c     *           nwpw_xc_sw1tosw2_dmatr          *
c     *                                           *
c     *********************************************
c
      subroutine nwpw_xc_sw1tosw2_dmatr(ispin,ne,nprj,sw1,sw2,
     >                            n1dgrid,nbasis,lmax2,
     >                            nindx,lm_rholm,
     >                            iprj_rholm,jprj_rholm,
     >                            bi_rholm,bj_rholm,
     >                            coeff_rholm,
     >                            coeff_rholm2,
     >                            coeff_rholm3,
     >                            dmatr)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      integer n1dgrid,nbasis,lmax2
      integer nindx,lm_rholm(*),iprj_rholm(*),jprj_rholm(*)
      integer bi_rholm(*),bj_rholm(*)
      real*8  coeff_rholm(*)
      real*8  coeff_rholm2(*)
      real*8  coeff_rholm3(*)
      real*8  dmatr(nbasis,nbasis,lmax2,3,ispin)

      integer ms,n,i,lm,iprj,jprj,bi,bj,n1(2),n2(2)
      real*8  tmp_ms,coeff,coeff2,coeff3,w,scal

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(21)
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      !scal = 1.0d0/dsqrt(lattice_omega())

*     ***init to zero***
      do i=1,nindx
         lm    = lm_rholm(i)
         iprj  = iprj_rholm(i)
         jprj  = jprj_rholm(i)
         bi    = bi_rholm(i)
         bj    = bj_rholm(i)
         !coeff = coeff_rholm(i)*scal
         coeff  = coeff_rholm(i)
         coeff2 = coeff_rholm2(i)
         coeff3 = coeff_rholm3(i)
         do ms=1,ispin
            if (ne(ms).gt.0) then
!$OMP DO
               do n=n1(ms),n2(ms)
                  sw2(n,iprj) = sw2(n,iprj)
     >                       + coeff *dmatr(bi,bj,lm+1,1,ms)*sw1(n,jprj)
     >                       + coeff2*dmatr(bi,bj,lm+1,2,ms)*sw1(n,jprj)
     >                       + coeff3*dmatr(bi,bj,lm+1,3,ms)*sw1(n,jprj)
               end do
!$OMP END DO
            end if
         end do
      end do
      call nwpw_timing_end(21)
      return
      end


