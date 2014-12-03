!$Id$
!
!$Log: not supported by cvs2svn $
!Revision 1.7  2005/01/30 07:16:17  bylaska
!Parallelized radial integrations...EJB
!
!Revision 1.6  2003/10/21 02:05:16  marat
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
!Revision 1.5  2003/02/18 01:58:33  marat
!removed write statement
!MV
!
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


Cc     *************************************************
Cc     *                                               *
Cc     *             init_paw_density                  *
Cc     *                                               *
Cc     *************************************************
C
C      subroutine init_paw_density(ispin,lmax)
C      implicit none
C
C      integer ispin,lmax
C#include "paw_density_data.fh"
C#include "paw_basis.fh"
C#include "bafdecls.fh"      
C#include "paw_ma.fh"
C#include "paw_geom.fh"
C
C      logical ok
C      integer in,ia
C      integer nion
C      integer paw_rho_size
C
C*     **** external functions ****
C      integer  control_gga
C      external control_gga
C
C      nion = ion_nion()
C
C*     *** lmax for multipole expansion of the density***      
C      paw_rho_lmax  = lmax
C      paw_rho_ispin = ispin
C
C      paw_rho_gga = control_gga()
C      
C*     *** size for multipole density arrays ***
C      ok = my_alloc(mt_int,nion,"i_paw_rho",
C     >                                i_paw_rho)
C
C      paw_rho_size = 0     
C      do in=1,nion
C        ia = ion_katm(in)
C        int_mb(i_paw_rho(1) + in - 1) = paw_rho_size
Cc        paw_rho_size = paw_rho_size + paw_basis_i_cut(ia)*
Cc     >                                (paw_rho_lmax+1)**2
C        paw_rho_size = paw_rho_size + paw_basis_npoints(ia)*
C     >                                (paw_rho_lmax+1)**2
C      end do
C
C*     *** allocate multipole density arrays ****
C      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ae",
C     >                               paw_rho_ae(1,1)) 
C      ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ps",
C     >                               paw_rho_ps(1,1)) 
C
C      if (paw_rho_gga.ge.10) then
C        ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ae",
C     >                                 paw_rho_ae_prime(1,1)) 
C        ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_u_ps",
C     >                                 paw_rho_ps_prime(1,1)) 
C      end if
C
C      if (paw_rho_ispin.eq.2) then
C        ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ae",
C     >                                 paw_rho_ae(1,2)) 
C        ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ps",
C     >                               paw_rho_ps(1,2)) 
C        if (paw_rho_gga.ge.10) then
C          ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ae",
C     >                                 paw_rho_ae_prime(1,2)) 
C          ok = ok .and. my_alloc(mt_dcpl,paw_rho_size,"paw_rho_d_ps",
C     >                               paw_rho_ps_prime(1,2)) 
C        end if
C
C      end if
C
C      end 
C
Cc     *********************************************
Cc     *                                           *
Cc     *           paw_density_solve               *
Cc     *                                           *
Cc     *********************************************
Cc
Cc   Calculates the atomic densities using from the
Cc   overlap coefficients.
C
C      subroutine paw_density_solve()
C      implicit none
C
C#include "bafdecls.fh"
C#include "paw_density_data.fh"
C#include "paw_geom.fh"
C#include "paw_basis.fh"
C#include "paw_proj.fh"
C#include "paw_ovlp_data.fh"
C      
C
C      integer i,li,mi,nilm
C      integer j,lj,mj,njlm
C      integer l,m,nlmi
C      integer in,ia,ms
C      integer nion
C      integer nb,nb2
C      integer ig,ng,ic,igg,istart
C      integer i_w,i_w0
C      integer i_rgrid
C      integer i_rho,i_rho0
C      integer i_phi_ae0
C      integer i_phi_ps0
C      double precision tmp
C      double precision tmp_gaunt
C      double complex    tmp_ms
C
C      double precision gen_gaunt_coeff
C      external gen_gaunt_coeff
C      
C      call nwpw_timing_start(21)
C      nion = ion_nion()
C
C*     ***init to zero***
C      do ms=1,paw_rho_ispin
C      call dcopy(2*paw_rho_ae(3,1),0.0d0,0,dcpl_mb(paw_rho_ae(1,ms)),1)
C      call dcopy(2*paw_rho_ae(3,1),0.0d0,0,dcpl_mb(paw_rho_ps(1,ms)),1)
C      end do
C
Cc      do i=1,paw_rho_u_ae(3)
Cc        dcpl_mb(paw_rho_ae(1,1)+i-1) = dcmplx(0,0) 
Cc        dcpl_mb(paw_rho_ae(1,2)+i-1) = dcmplx(0,0) 
Cc        dcpl_mb(paw_rho_ps(1,1)+i-1) = dcmplx(0,0) 
Cc        dcpl_mb(paw_rho_ps(1,2)+i-1) = dcmplx(0,0) 
Cc      end do
C
C
C      do in=1,nion
C        ia = ion_katm(in)
C        nb = paw_basis_nbasis(ia)
C        nb2 = nb*nb
C        ng = paw_basis_ngrid(ia)
Cc        ic = paw_basis_i_cut(ia)
C        ic = paw_basis_npoints(ia)
C        istart = paw_basis_i_start(ia)
C        i_rgrid = paw_basis_i_rgrid(ia)
C        i_phi_ae0 = paw_basis_i_phi_ae(ia)
C        i_phi_ps0 = paw_basis_i_phi_ps(ia)
C        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
C        i_rho0 = int_mb(i_paw_rho(1) + in - 1)
C        nilm = 0
C        do i=1,nb
C        li = paw_basis_orb_l(i,ia)
C        njlm = 0
C        do j=1,nb
C        lj = paw_basis_orb_l(j,ia)
C        nlmi = 0
C        do l=0,paw_rho_lmax
C        if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
C        do m=-l,l
C        do mi=-li,li
C        do mj=-lj,lj
C        if(m.eq.(mj-mi)) then
C          i_w   = i_w0-1+(njlm+lj+mj+1)+
C     >            (nilm+li+mi)*paw_proj_nbasis(ia)
C          tmp_gaunt = gen_gaunt_coeff(l,m,lj,mj,li,mi)
C
C          do ms=1,paw_rho_ispin
C            tmp_ms = tmp_gaunt* 
C     >             dcpl_mb(paw_ovlp_w_spin(1,ms) + i_w)       
C          
C            do ig=1,ic
C            igg = ig+istart-1
C
C            i_rho   = i_rho0-1+ig+
C     >            (nlmi+l+m)*ic
C          
C            tmp = dbl_mb(i_phi_ae0+(i-1)*ng+igg-1)*       
C     >                 dbl_mb(i_phi_ae0+(j-1)*ng+igg-1)/
C     >                 dbl_mb(i_rgrid+igg-1)**2 
C
C            dcpl_mb(paw_rho_ae(1,ms)+i_rho)=
C     >                 dcpl_mb(paw_rho_ae(1,ms)+i_rho)+
C     >                 tmp_ms*tmp
C                      
C            tmp = dbl_mb(i_phi_ps0+(i-1)*ng+igg-1)*       
C     >                 dbl_mb(i_phi_ps0+(j-1)*ng+igg-1)/
C     >                 dbl_mb(i_rgrid+igg-1)**2 
C
C            dcpl_mb(paw_rho_ps(1,ms)+i_rho)=
C     >                 dcpl_mb(paw_rho_ps(1,ms)+i_rho)+
C     >                 tmp_ms*tmp
C                      
C            end do
C           end do
C
C        end if
C        end do
C        end do
C        end do !m
C        end if
C        nlmi = nlmi + 2*l+1
C        end do !l
C        njlm = njlm + 2*lj+1
C        end do !j
C        nilm = nilm + 2*li+1
C        end do !i
C 
C      end do!in
C
C      call nwpw_timing_end(21)
C      return
C      end
C
Cc     **********************************************
Cc     *                                            *
Cc     *             paw_density_end                *
Cc     *                                            *
Cc     **********************************************
C
C      subroutine paw_density_end()
C      implicit none
C
C#include "paw_density_data.fh"
C#include "bafdecls.fh"      
C#include "paw_ma.fh"
C
C      logical ok
C      integer ms
C
C      ok =          my_dealloc(i_paw_rho)
C      do ms=1,paw_rho_ispin
C        ok = ok .and. my_dealloc(paw_rho_ae(1,ms)) 
C        ok = ok .and. my_dealloc(paw_rho_ps(1,ms))
C        if (paw_rho_gga.ge.10) then
C          ok = ok .and. my_dealloc(paw_rho_ae_prime(1,ms))
C          ok = ok .and. my_dealloc(paw_rho_ps_prime(1,ms))
C        end if
C      end do 
C      
C      if (.not.ok)
C     > call errquit("paw_density_end: error freeing heap",0,0)
C
C
C      end 
C
C
C
Cc     ********************************************
Cc     *                                          *
Cc     *        paw_density_prime_solve           *
Cc     *                                          *
Cc     ********************************************
C
Cc  This routine computes the radial gradient of the 
Cc  atomic densities from the overlap coefficients.
Cc
C
C      subroutine paw_density_prime_solve()
C      implicit none
C
C#include "bafdecls.fh"
C#include "paw_density_data.fh"
C#include "paw_geom.fh"
C#include "paw_basis.fh"
C#include "paw_proj.fh"
C#include "paw_ovlp_data.fh"
C      
C
C      integer i,li,mi,nilm
C      integer j,lj,mj,njlm
C      integer l,m,nlmi
C      integer in,ia,ms
C      integer nion
C      integer nb,nb2
C      integer ig,ng,ic,igg,istart
C      integer i_w,i_w0
C      integer i_rgrid
C      integer i_rho,i_rho0
C      integer i_phi_ae0,i_phi_ae0_prime
C      integer i_phi_ps0,i_phi_ps0_prime
C      double precision tmp
C      double precision tmp_gaunt
C      double complex   tmp_ms
C
C      double precision gen_gaunt_coeff
C      external gen_gaunt_coeff
C      
C      call nwpw_timing_start(21)
C      nion = ion_nion()
C
C*     ***init to zero***
C      do ms=1,paw_rho_ispin
C      call dcopy(2*paw_rho_ae_prime(3,1),0.0d0,0,
C     >           dcpl_mb(paw_rho_ae_prime(1,ms)),1)
C      call dcopy(2*paw_rho_ps_prime(3,1),0.0d0,0,
C     >           dcpl_mb(paw_rho_ps_prime(1,ms)),1)
C      end do
C
C
C
C      do in=1,nion
C        ia = ion_katm(in)
C        nb = paw_basis_nbasis(ia)
C        nb2 = nb*nb
C        ng = paw_basis_ngrid(ia)
Cc        ic = paw_basis_i_cut(ia)
C        ic = paw_basis_npoints(ia)
C        istart = paw_basis_i_start(ia)
C        i_rgrid = paw_basis_i_rgrid(ia)
C        i_phi_ae0       = paw_basis_i_phi_ae(ia)
C        i_phi_ps0       = paw_basis_i_phi_ps(ia)
C        i_phi_ae0_prime = paw_basis_i_phi_ae_prime(ia)
C        i_phi_ps0_prime = paw_basis_i_phi_ps_prime(ia)
C        i_w0 = int_mb(i_paw_ovlp_w(1) + in - 1)
C        i_rho0 = int_mb(i_paw_rho(1) + in - 1)
C        nilm = 0
C        do i=1,nb
C        li = paw_basis_orb_l(i,ia)
C        njlm = 0
C        do j=1,nb
C        lj = paw_basis_orb_l(j,ia)
C        nlmi = 0
C        do l=0,paw_rho_lmax
C        if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
C        do m=-l,l
C        do mi=-li,li
C        do mj=-lj,lj
C        if(m.eq.(mj-mi)) then
C          i_w   = i_w0-1+(njlm+lj+mj+1)+
C     >            (nilm+li+mi)*paw_proj_nbasis(ia)
C          tmp_gaunt = gen_gaunt_coeff(l,m,lj,mj,li,mi)
C
C          do ms=1,paw_rho_ispin
C            tmp_ms = tmp_gaunt* 
C     >             dcpl_mb(paw_ovlp_w_spin(1,ms) + i_w)       
C          
Cc            call paw_density_sub2(ic,tmp_ms,
Cc     >           dbl_mb(i_rgrid+istart-1),
Cc     >           dbl_mb(i_phi_ae0       +istart-1 +(i-1)*ng),
Cc     >           dbl_mb(i_phi_ae0_prime +istart-1 +(i-1)*ng),
Cc     >           dbl_mb(i_phi_ae0       +istart-1 +(j-1)*ng),
Cc     >           dbl_mb(i_phi_ae0_prime +istart-1 +(j-1)*ng),
Cc     >           dcpl_mb(paw_rho_ae_prime(1,ms)+irho0+(nlmi+l+m)*ic))
C
C            do ig=1,ic
C            igg = ig+istart-1
C
C            i_rho   = i_rho0-1+ig+
C     >            (nlmi+l+m)*ic
C          
C            tmp = (dbl_mb(i_phi_ae0_prime+(i-1)*ng+igg-1)       
C     >            *dbl_mb(i_phi_ae0+(j-1)*ng+igg-1)
C     >          +  dbl_mb(i_phi_ae0+(i-1)*ng+igg-1)
C     >            *dbl_mb(i_phi_ae0_prime+(j-1)*ng+igg-1))
C     >           /dbl_mb(i_rgrid+igg-1)**2 
C     >          - 2.0d0*(dbl_mb(i_phi_ae0+(i-1)*ng+igg-1)
C     >                 *dbl_mb(i_phi_ae0+(j-1)*ng+igg-1)
C     >                 /dbl_mb(i_rgrid+igg-1)**3)
C
C            dcpl_mb(paw_rho_ae_prime(1,ms)+i_rho)=
C     >                 dcpl_mb(paw_rho_ae_prime(1,ms)+i_rho)+
C     >                 tmp_ms*tmp
C
C            tmp =  (dbl_mb(i_phi_ps0_prime+(i-1)*ng+igg-1)       
C     >             *dbl_mb(i_phi_ps0+(j-1)*ng+igg-1)
C     >          +   dbl_mb(i_phi_ps0+(i-1)*ng+igg-1)
C     >             *dbl_mb(i_phi_ps0_prime+(j-1)*ng+igg-1))
C     >            /dbl_mb(i_rgrid+igg-1)**2 
C     >          - 2.0d0*(dbl_mb(i_phi_ps0+(i-1)*ng+igg-1)
C     >                 *dbl_mb(i_phi_ps0+(j-1)*ng+igg-1)
C     >                 /dbl_mb(i_rgrid+igg-1)**3)
C                      
C
C            dcpl_mb(paw_rho_ps_prime(1,ms)+i_rho)=
C     >                 dcpl_mb(paw_rho_ps_prime(1,ms)+i_rho)+
C     >                 tmp_ms*tmp
C                      
C            end do
C           end do
C
C        end if
C        end do
C        end do
C        end do !m
C        end if
C        nlmi = nlmi + 2*l+1
C        end do !l
C        njlm = njlm + 2*lj+1
C        end do !j
C        nilm = nilm + 2*li+1
C        end do !i
C 
C      end do!in
C
C      call nwpw_timing_end(21)
C      return
C      end
C
Cc      subroutine paw_density_sub2(ic,alpha,
Cc     >                            r,
Cc     >                            phi_i,dphi_i,
Cc     >                            phi_j,dphi_j,
Cc     >                            rho_prime)
Cc      implicit none
Cc      integer ic
Cc      double precision alpha
Cc      double precision r(*)
Cc      double precision phi_i(*),dphi_i(*)
Cc      double precision phi_j(*),dphi_j(*)
Cc      double precision rho_prime(*)
Cc
Cc*     **** local variables ****
Cc      integer          ig
Cc      double precision tmp
Cc
Cc      do ig=1,ic
Cc         tmp = dphi_i(ig)*phi_j(ig)/r(ig)**2
Cc     >       + phi_i(ig)*dphi_j(ig)/r(ig)**2
Cc     >       - 2.0d0*phi_i(ig)*phi_j(ig)/r(ig)**3
Cc
Cc         rho_prime(ig)=rho_prime(ig) + alpha*tmp
Cc      end do
Cc
Cc      return
Cc      end



c     *********************************************
c     *                                           *
c     *           paw_density_solve2              *
c     *                                           *
c     *********************************************
c
c   Calculates the atomic density lm expansions  from the
c   overlap coefficients.

      subroutine paw_density_solve2(ii,lmax,ic,lmax2,ispin0,
     >                              paw_rho_ae,paw_rho_ps)
      implicit none
      integer ii
      integer lmax
      integer ic,lmax2,ispin0
      double complex paw_rho_ae(ic,lmax2,ispin0)
      double complex paw_rho_ps(ic,lmax2,ispin0)

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
      

      integer i,li,mi,nilm
      integer j,lj,mj,njlm
      integer l,m,nlmi
      integer ia,ms,l2
      integer nion
      integer nb,nb2
      integer ig,ng,igg,istart
      integer i_w,i_w0
      integer i_rgrid
      integer i_rho,i_rho0
      integer i_phi_ae0
      integer i_phi_ps0
      double precision tmp_ps,tmp_ae
      double precision tmp_gaunt
      double complex    tmp_ms

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      
      call nwpw_timing_start(21)


*     ***init to zero***
      call dcopy(2*ic*lmax2*ispin0,0.0d0,0,paw_rho_ae,1)
      call dcopy(2*ic*lmax2*ispin0,0.0d0,0,paw_rho_ps,1)


      ia = ion_katm(ii)
      nb = paw_basis_nbasis(ia)
      nb2 = nb*nb
      ng = paw_basis_ngrid(ia)
c      ic = paw_basis_npoints(ia)
      istart = paw_basis_i_start(ia)
      i_rgrid = paw_basis_i_rgrid(ia)
      i_phi_ae0 = paw_basis_i_phi_ae(ia)
      i_phi_ps0 = paw_basis_i_phi_ps(ia)
      i_w0 = int_mb(i_paw_ovlp_w(1) + ii - 1)


      nilm = 0
      do i=1,nb
      li = paw_basis_orb_l(i,ia)
      njlm = 0
      do j=1,nb
      lj = paw_basis_orb_l(j,ia)
      nlmi = 0
      do l=0,lmax
        if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
        do m=-l,l
        do mi=-li,li
        do mj=-lj,lj
        if(m.eq.(mj-mi)) then

          i_w   = i_w0-1+(njlm+lj+mj+1)+
     >            (nilm+li+mi)*paw_proj_nbasis(ia)
          tmp_gaunt = gen_gaunt_coeff(l,m,lj,mj,li,mi)

          l2 = (nlmi+l+m)+1

          do ms=1,ispin0
            tmp_ms = tmp_gaunt* 
     >             dcpl_mb(paw_ovlp_w_spin(1,ms) + i_w)       
          
            call paw_density_gen_rho(ic,tmp_ms,
     >                              dbl_mb(i_rgrid+istart-1),
     >                              dbl_mb(i_phi_ae0+(i-1)*ng+istart-1),
     >                              dbl_mb(i_phi_ae0+(j-1)*ng+istart-1),
     >                              paw_rho_ae(1,l2,ms))
            call paw_density_gen_rho(ic,tmp_ms,
     >                              dbl_mb(i_rgrid+istart-1),
     >                              dbl_mb(i_phi_ps0+(i-1)*ng+istart-1),
     >                              dbl_mb(i_phi_ps0+(j-1)*ng+istart-1),
     >                              paw_rho_ps(1,l2,ms))
c            do ig=1,ic
c            igg = ig+istart-1
c
c            tmp_ae = dbl_mb(i_phi_ae0+(i-1)*ng+igg-1)*       
c     >                 dbl_mb(i_phi_ae0+(j-1)*ng+igg-1)/
c     >                 dbl_mb(i_rgrid+igg-1)**2 
c            tmp_ps = dbl_mb(i_phi_ps0+(i-1)*ng+igg-1)*       
c     >                 dbl_mb(i_phi_ps0+(j-1)*ng+igg-1)/
c     >                 dbl_mb(i_rgrid+igg-1)**2 
c
c            paw_rho_ae(ig,l2,ms) = paw_rho_ae(ig,l2,ms) + tmp_ms*tmp_ae
c            paw_rho_ps(ig,l2,ms) = paw_rho_ps(ig,l2,ms) + tmp_ms*tmp_ps
c            end do

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
 
      call nwpw_timing_end(21)
      return
      end


      subroutine paw_density_gen_rho(ic,alpha,r,phi1,phi2,rho)
      implicit none
      integer ic
      double precision alpha(2)
      double precision r(ic)
      double precision phi1(ic)
      double precision phi2(ic)
      double precision rho(2,ic)

*     **** local variables ****
      integer i
      double precision tmp

      do i=1,ic
         tmp = (phi1(i)*phi2(i))/(r(i)**2)
         rho(1,i) = rho(1,i) + alpha(1)*tmp
         rho(2,i) = rho(2,i) + alpha(2)*tmp
      end do

      return
      end





c     *********************************************
c     *                                           *
c     *           paw_density_prime_solve2        *
c     *                                           *
c     *********************************************
c
c   Calculates the atomic density lm expansions  from the
c   overlap coefficients.

      subroutine paw_density_prime_solve2(ii,lmax,ic,lmax2,ispin0,
     >                              paw_rho_ae_prime,paw_rho_ps_prime)
      implicit none
      integer ii
      integer lmax
      integer ic,lmax2,ispin0
      double complex paw_rho_ae_prime(ic,lmax2,ispin0)
      double complex paw_rho_ps_prime(ic,lmax2,ispin0)

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"
      

      integer i,li,mi,nilm
      integer j,lj,mj,njlm
      integer l,m,nlmi
      integer ia,ms,l2
      integer nion
      integer nb,nb2
      integer ig,ng,igg,istart
      integer i_w,i_w0
      integer i_rgrid
      integer i_rho,i_rho0
      integer i_phi_ae0,i_phi_ae0_prime
      integer i_phi_ps0,i_phi_ps0_prime
      double precision tmp_ps,tmp_ae
      double precision tmp_gaunt
      double complex    tmp_ms

      double precision gen_gaunt_coeff
      external gen_gaunt_coeff
      
      call nwpw_timing_start(21)


*     ***init to zero***
      call dcopy(2*ic*lmax2*ispin0,0.0d0,0,paw_rho_ae_prime,1)
      call dcopy(2*ic*lmax2*ispin0,0.0d0,0,paw_rho_ps_prime,1)


      ia = ion_katm(ii)
      nb = paw_basis_nbasis(ia)
      nb2 = nb*nb
      ng = paw_basis_ngrid(ia)
c      ic = paw_basis_npoints(ia)
      istart  = paw_basis_i_start(ia)
      i_rgrid = paw_basis_i_rgrid(ia)

      i_phi_ae0       = paw_basis_i_phi_ae(ia)
      i_phi_ps0       = paw_basis_i_phi_ps(ia)
      i_phi_ae0_prime = paw_basis_i_phi_ae_prime(ia)
      i_phi_ps0_prime = paw_basis_i_phi_ps_prime(ia)

      i_w0 = int_mb(i_paw_ovlp_w(1) + ii - 1)


      nilm = 0
      do i=1,nb
      li = paw_basis_orb_l(i,ia)
      njlm = 0
      do j=1,nb
      lj = paw_basis_orb_l(j,ia)
      nlmi = 0
      do l=0,lmax
        if(l.le.(li+lj) .and. l.ge.abs(li-lj)) then
        do m=-l,l
        do mi=-li,li
        do mj=-lj,lj
        if(m.eq.(mj-mi)) then

          i_w   = i_w0-1+(njlm+lj+mj+1)+
     >            (nilm+li+mi)*paw_proj_nbasis(ia)
          tmp_gaunt = gen_gaunt_coeff(l,m,lj,mj,li,mi)

          l2 = (nlmi+l+m)+1

          do ms=1,ispin0
            tmp_ms = tmp_gaunt* 
     >             dcpl_mb(paw_ovlp_w_spin(1,ms) + i_w)       
          
            call paw_density_gen_drho(ic,tmp_ms,
     >              dbl_mb(i_rgrid                 +istart-1),
     >              dbl_mb(i_phi_ae0      +(i-1)*ng+istart-1),
     >              dbl_mb(i_phi_ae0_prime+(i-1)*ng+istart-1),
     >              dbl_mb(i_phi_ae0      +(j-1)*ng+istart-1),
     >              dbl_mb(i_phi_ae0_prime+(j-1)*ng+istart-1),
     >              paw_rho_ae_prime(1,l2,ms))
            call paw_density_gen_drho(ic,tmp_ms,
     >              dbl_mb(i_rgrid                 +istart-1),
     >              dbl_mb(i_phi_ps0      +(i-1)*ng+istart-1),
     >              dbl_mb(i_phi_ps0_prime+(i-1)*ng+istart-1),
     >              dbl_mb(i_phi_ps0      +(j-1)*ng+istart-1),
     >              dbl_mb(i_phi_ps0_prime+(j-1)*ng+istart-1),
     >              paw_rho_ps_prime(1,l2,ms))

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
 
      call nwpw_timing_end(21)
      return
      end






      subroutine paw_density_gen_drho(ic,alpha,r,
     >                                phi1,dphi1,
     >                                phi2,dphi2,
     >                                drho)
      implicit none
      integer ic
      double precision alpha(2)
      double precision r(ic)
      double precision phi1(ic),dphi1(ic)
      double precision phi2(ic),dphi2(ic)
      double precision drho(2,ic)

*     **** local variables ****
      integer i
      double precision tmp

      do i=1,ic
         tmp = (dphi1(i)*phi2(i)+phi1(i)*dphi2(i))/(r(i)**2)
     >       - 2.0d0*(phi1(i)*phi2(i))/(r(i)**3)
         drho(1,i) = drho(1,i) + alpha(1)*tmp
         drho(2,i) = drho(2,i) + alpha(2)*tmp
      end do

      return
      end


