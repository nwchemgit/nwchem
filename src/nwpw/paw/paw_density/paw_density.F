!$Id: paw_density.F,v 1.5 2003-02-18 01:58:33 marat Exp $
!
!$Log: not supported by cvs2svn $
!Revision 1.4  2003/02/17 05:58:08  marat
!added lmax as an argument to init_paw_density
!MV
!
!Revision 1.3  2003/02/17 05:53:40  marat
!fixed bug in deallocation
!MV
!
!Revision 1.2  2003/02/17 00:52:31  marat
!generated multipole expansion of the atomic densities
!freed the memory
!
!MV
!
!Revision 1.1  2003/02/16 21:40:00  marat
!initial version
!MV
!
**********************************
*
*
*
*
*
**********************************
      subroutine init_paw_density(lmax)
      implicit none

      integer lmax
#include "paw_density_data.fh"
#include "paw_basis.fh"
#include "mafdecls.fh"      
#include "paw_ma.fh"
#include "paw_geom.fh"

      logical ok
      integer in,ia
      integer nion
      integer paw_rho_size

      nion = ion_nion()
*     *** lmax for multipole expansion of the density***      
      paw_rho_lmax = lmax
      
*     *** size for multipole density arrays ***
      ok = my_alloc(mt_int,nion,"i_paw_rho",
     >                                i_paw_rho)

      paw_rho_size = 0     
      do in=1,nion
        ia = ion_katm(in)
        int_mb(i_paw_rho(1) + in - 1) = paw_rho_size
        paw_rho_size = paw_rho_size + paw_basis_i_cut(ia)*
     >                                (paw_rho_lmax+1)**2
      end do

*     *** allocate multipole density arrays ****
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ae",
     >                               paw_rho_u_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ae",
     >                               paw_rho_d_ae) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ps",
     >                               paw_rho_u_ps) 
      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ps",
     >                               paw_rho_d_ps) 




      end 
**********************************
*
*
*
*
*
**********************************
      subroutine paw_density_solve()
      implicit none

#include "mafdecls.fh"
#include "paw_density_data.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
      

      integer i,li,mi,nilm
      integer j,lj,mj,njlm
      integer l,m,nlmi
      integer in,ia
      integer nion
      integer nb,nb2
      integer ig,ng,ic
      integer i_w,i_w0
      integer i_rgrid
      integer i_rho,i_rho0
      integer i_phi_ae0
      integer i_phi_ps0
      double precision tmp
      double precision tmp_gaunt
      double precision tmp_u,tmp_d

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      
      nion = ion_nion()

*     ***init to zero***
      do i=1,paw_rho_u_ae(3)
        dcpl_mb(paw_rho_u_ae(1)+i-1) = dcmplx(0,0) 
        dcpl_mb(paw_rho_d_ae(1)+i-1) = dcmplx(0,0) 
        dcpl_mb(paw_rho_u_ps(1)+i-1) = dcmplx(0,0) 
        dcpl_mb(paw_rho_d_ps(1)+i-1) = dcmplx(0,0) 
      end do


      do in=1,nion
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        nb2 = nb*nb
        ng = paw_basis_ngrid(ia)
        ic = paw_basis_i_cut(ia)
        i_rgrid = paw_basis_i_rgrid(ia)
        i_phi_ae0 = paw_basis_i_phi_ae(ia)
        i_phi_ps0 = paw_basis_i_phi_ps(ia)
        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
        i_rho0 = int_mb(i_paw_rho(1) + in - 1)
        nilm = 0
        do i=1,nb
        li = paw_basis_orb_l(i,ia)
        njlm = 0
        do j=1,nb
        lj = paw_basis_orb_l(j,ia)
        nlmi = 0
        do l=0,paw_rho_lmax
        if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
        do m=-l,l
        do mi=-li,li
        do mj=-lj,lj
        if(m.eq.(mj-mi)) then
          i_w   = i_w0-1+(njlm+lj+mj+1)+
     >            (nilm+li+mi)*paw_proj_nbasis(ia)
          tmp_gaunt = gen_gaunt_coeff(l,m,lj,mj,li,mi)
          tmp_u = tmp_gaunt* 
     >             dcpl_mb(paw_ovlp_w_u(1) + i_w)       
          tmp_d = tmp_gaunt*
     >             dcpl_mb(paw_ovlp_w_d(1) + i_w)       
          
          do ig=1,ic

            i_rho   = i_rho0-1+ig+
     >            (nlmi+l+m)*ic
          
            tmp = dbl_mb(i_phi_ae0+(i-1)*ng+ig-1)*       
     >                 dbl_mb(i_phi_ae0+(j-1)*ng+ig-1)/
     >                 dbl_mb(i_rgrid+ig-1)**2 

            dcpl_mb(paw_rho_u_ae(1)+i_rho)=
     >                 dcpl_mb(paw_rho_u_ae(1)+i_rho)+
     >                 tmp_u*tmp
                      
            dcpl_mb(paw_rho_d_ae(1)+i_rho)=
     >                 dcpl_mb(paw_rho_d_ae(1)+i_rho)+
     >                 tmp_d*tmp

            tmp = dbl_mb(i_phi_ps0+(i-1)*ng+ig-1)*       
     >                 dbl_mb(i_phi_ps0+(j-1)*ng+ig-1)/
     >                 dbl_mb(i_rgrid+ig-1)**2 

            dcpl_mb(paw_rho_u_ps(1)+i_rho)=
     >                 dcpl_mb(paw_rho_u_ps(1)+i_rho)+
     >                 tmp_u*tmp
                      
            dcpl_mb(paw_rho_d_ps(1)+i_rho)=
     >                 dcpl_mb(paw_rho_d_ps(1)+i_rho)+
     >                 tmp_d*tmp
            end do

         end if
          end do
          end do
          end do !m
          end if
          nlmi = nlmi + 2*l+1
          end do !l
          njlm = njlm + 2*lj+1
          end do !j
          nilm = nilm + 2*li+1
          end do !i
 
      end do!in

c      write(*,*) "final in paw density"
c      do in=1,nion
c        ia = ion_katm(in)
c        ng = paw_basis_ngrid(ia)
c        ic = paw_basis_i_cut(ia)
c        i_rho0 = int_mb(i_paw_rho(1) + in - 1)
c        nlmi = 0
c        do l=0,paw_rho_lmax
c        do m=-l,l
c          write(*,*) "na,l,m",in,l,m
c          do ig=1,ic
c
c          i_rho   = i_rho0-1+ig+
c     >            (nlmi+l+m)*ic
c                     
c            write(*,*) ig,dcpl_mb(paw_rho_u_ps(1)+i_rho),i_rho
c          end do
c
c          end do !m
c          nlmi = nlmi + 2*l+1
c          end do !l
c      end do!in

      end

      subroutine paw_density_end()
      implicit none

#include "paw_density_data.fh"
#include "mafdecls.fh"      
#include "paw_ma.fh"

      logical ok

      
      ok = my_dealloc(paw_rho_d_ps) 
      ok = ok .and. my_dealloc(paw_rho_u_ps) 
      ok = ok .and. my_dealloc(paw_rho_d_ae) 
      ok = ok .and. my_dealloc(paw_rho_u_ae) 
      ok = ok .and. my_dealloc(i_paw_rho)

      if (.not.ok)
     > call errquit("paw_density_end: error freeing heap",0)


      end 
