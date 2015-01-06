*
* $Id$
*


!**************************************************
!
!       Name: paw_core_atom_energy
!
!       Purpose: calculates atom portion of the coreetic energy
!
!       Created: Mon Feb 10 13:46:23 PST 2003       
!**************************************************
      function paw_energy_core_atom()
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_core_data.fh"
#include "paw_basis.fh"
#include "paw_ovlp_data.fh"

      double precision  paw_energy_core_atom
      integer in,ia
      integer i_w,i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer nb
      double precision ecore

      ecore = 0.D0
      i_w = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        i_mtr0 = int_mb(i_paw_core_matrix(1) + ia - 1)
        do i=1,nb
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
            do j=1,nb
              lj = paw_basis_orb_l(j,ia)
                do mj=-lj,lj
                i_mtr = i_mtr0 +i+(j-1)*nb
                if(mi.eq.mj) then
                  ecore = ecore + dcpl_mb(paw_ovlp_w(1)+i_w)*
     >                          dbl_mb(paw_core_matrix(1) + i_mtr - 1) 
                end if
                i_w = i_w +1
                end do
             end do
           end do   
         end do
      end do
 
      paw_energy_core_atom = ecore

      end

!**************************************************
!
!       Name: paw_energy_core_ion_atom
!
!       Purpose: calculates core ion interac. energy
!
!       Created: Mon Feb 10 13:46:23 PST 2003       
!**************************************************
      function paw_energy_core_ion_atom()
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_core_data.fh"
#include "paw_basis.fh"
#include "integrate.fh"
#include "paw_params.fh"

      double precision  paw_energy_core_ion_atom
      integer k,igg,np,istart
      integer in,ia
      integer i_w,i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer nb
      integer i_r
      integer ic,max_ic
      integer i_core_ae
      integer  i_f1,l_f1
      logical ok
      double precision log_amesh
      double precision ecore

      call Parallel_np(np)
      !*** temporary array for calculations ***
      max_ic = paw_basis_max_i_cut()
      ok = ma_push_get(mt_dbl,max_ic,"tmp_array_core",l_f1,i_f1)
      if (.not.ok)
     > call errquit(
     >  'paw_energy_core_ion_atom: out of stack memory',0,0)

      ecore = 0.D0
      i_w = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
c        ic = paw_basis_i_cut(ia)
        ic = paw_basis_npoints(ia)
        istart = paw_basis_i_start(ia)
        i_core_ae = paw_basis_i_core_density(ia)
        i_r = paw_basis_i_rgrid(ia)
        log_amesh = paw_basis_log_amesh(ia)
        do k=1,ic
          igg=k+istart-1

          dbl_mb(i_f1+k-1)=-paw_basis_ion_charge(ia)*
     >                      dbl_mb(i_core_ae+igg-1)/
     >                      dbl_mb(i_r+igg-1)       

        end do
c        ecore = ecore + def_integr(-1,dbl_mb(i_f1),2,dbl_mb(i_r),
c     >                       log_amesh,ic)
        ecore = ecore + def_integr1(-1,dbl_mb(i_f1),2,
     >                       dbl_mb(i_r+istart-1),
     >                       log_amesh,ic,istart)

      end do
      if (np.gt.1) call D3dB_SumAll(ecore)
 
      ok = ma_pop_stack(l_f1)
      if (.not.ok)
     > call errquit(
     >  'paw_energy_core_ion_atom:error popping stack memory',0,1)
      paw_energy_core_ion_atom = ecore*fourpi
 
      end

!**************************************************
!
!       Name: paw_energy_core_self_atom
!
!       Purpose: calculates core self energy
!
!       Created: Mon Feb 10 13:46:23 PST 2003       
!**************************************************
      function paw_energy_core_self_atom()
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_core_data.fh"
#include "paw_basis.fh"
#include "integrate.fh"
#include "paw_params.fh"
#include "paw_utils.fh"

      double precision  paw_energy_core_self_atom
      integer k,igg,np
      integer in,ia
      integer i_w,i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer nb
      integer i_r
      integer ic,max_ic
      integer i_core_ae
      integer  i_f1,l_f1
      logical ok
      double precision log_amesh
      double precision ecore

      !*** temporary array for calculations ***
      max_ic = paw_basis_max_i_cut()
      ok = ma_push_get(mt_dbl,max_ic,"tmp_array",l_f1,i_f1)
      if (.not.ok)
     > call errquit(
     >  'paw_energy_core_self_atom: out of stack memory',0,0)

      ecore = 0.0d0
      i_w = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        ic = paw_basis_i_cut(ia)
        i_core_ae = paw_basis_i_core_density(ia)
        i_r = paw_basis_i_rgrid(ia)
        log_amesh = paw_basis_log_amesh(ia)
        do k=1,ic

          dbl_mb(i_f1+k-1)= dbl_mb(i_core_ae+k-1)*
     >                      dbl_mb(i_r+k-1)**2       

        end do
        ecore = ecore + find_multipole_energy(0,ic,dbl_mb(i_r),2,
     >                                  dbl_mb(i_f1),2,
     >                                  dbl_mb(i_f1),log_amesh)

      end do
 
      ok = ma_pop_stack(l_f1)
      if (.not.ok)
     > call errquit(
     >  'paw_energy_core_self_atom: error popping stack',0,1)

      paw_energy_core_self_atom = 0.5d0*ecore*fourpi
 
      end      
