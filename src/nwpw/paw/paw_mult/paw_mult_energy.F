*
* $Id: paw_mult_energy.F,v 1.7 2005-01-31 02:09:38 bylaska Exp $
*
*
* $Log: not supported by cvs2svn $
* Revision 1.6  2004/07/29 15:55:37  bylaska
* Temporary array for Gaunt coefficients added.  Speeds up the program considerably, but it is extrememely memory intensive.  Also added timing routines to multipole calculations and fixed the initial total charge calculation.
*
*  ...EJB
*
* Revision 1.5  2003/03/05 23:16:32  bylaska
* Commented out write statements and other minor fixes.....
* self-consistent loop looks like it is working.....
* ....EJB
*
* Revision 1.4  2003/02/26 20:34:06  marat
* fixed bug related to calculation of
* comp_coeff, the i and j loop were incorrectly switched
* MV
*
* Revision 1.3  2003/02/23 23:06:48  bylaska
* energy subroutines changed to functions...EJB
*
* Revision 1.2  2003/02/23 22:45:52  bylaska
* bug fixes - energies agree with fortran90 code...EJB
*
* Revision 1.1  2003/02/23 21:37:07  bylaska
* routines for calculating atomic multipole energies have been added....EJB
*


  !*************************************************
  !
  !   Name    : paw_mult_energy_atom_comp
  !
  !   Purpose :
  !
  !   Created :
  !
  !*************************************************

      function paw_mult_energy_atom_comp()
      implicit none
      real*8 paw_mult_energy_atom_comp !*RESULT*

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer ii,ia,li,mi,il,mult_li
      integer nion
      real*8 e_comp

      nion=ion_nion()
      e_comp = 0.0d0
      do ii=1,nion
        ia      = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il      = int_mb(i_v_mult(1)+ii-1)
        do li=0,mult_li
        do mi=-li,li

          e_comp = e_comp
     >           + dble(dcpl_mb(paw_qlm(1)+il)
     >                 *dcpl_mb(comp_coeff(1)+il))

!         write(33,*) li,mi,ia,e_comp,dcpl_mb(paw_qlm(1)+il),
!    >                dcpl_mb(comp_coeff(1)+il)     

          il = il + 1
        end do !*mi*
        end do !*li*
      end do !*ii*

      paw_mult_energy_atom_comp = e_comp
      !write(*,*) " paw_mult_energy_atom_comp",e_comp
      return
      end


  !*************************************************
  !
  !   Name    : paw_mult_energy_atom_self
  !
  !   Purpose :
  !
  !   Created :
  !
  !*************************************************

      function paw_mult_energy_atom_self()
      implicit none
      real*8 paw_mult_energy_atom_self !*RESULT*

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer ii,ia,li,mi,il,mult_li
      integer indx,lmax,nion
      real*8 e_self

      nion=ion_nion()
      lmax=paw_basis_max_mult_l()

      e_self = 0.0d0
      do ii=1,nion
        ia     = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il = int_mb(i_v_mult(1)+ii-1)
        do li=0,mult_li
          indx = li + (ia-1)*(lmax+1)
          do mi=-li,li

          e_self
     >     = e_self
     >     - 0.5d0*dble(
     >             dconjg(dcpl_mb(paw_qlm(1)+il))
     >            *dcpl_mb(paw_qlm(1)+il))
     >            *dbl_mb(self_energy_coeff(1)+indx)
          il = il +1
          end do !*mi*
        end do !*li*
      end do !*ii*

      paw_mult_energy_atom_self = e_self
      return
      end


  !*************************************************
  !
  !   Name    : paw_mult_energy_atom_mult
  !
  !   Purpose :
  !
  !   Created :
  !
  !*************************************************

      function paw_mult_energy_atom_mult()
      implicit none
      real*8 paw_mult_energy_atom_mult !*RESULT*

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer taskid,np,tcount
      integer ii,ia,li,mi,il,ill,mult_li
      integer jj,ja,lj,mj,jl,jll,mult_lj
      integer nion,indx,lmax,lmax2
      real*8 e_mult

      call nwpw_timing_start(14)
      call Parallel_np(np)
      call Parallel_taskid(taskid)

      nion  = ion_nion()
      lmax  = paw_basis_max_mult_l()
      lmax2 = (lmax+1)**2

      tcount = 0
      e_mult = 0.0d0
      do ii=1,nion
        ia     = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il     = int_mb(i_v_mult(1)+ii-1)
        ill    = 0
        do li=0,mult_li
        do mi=-li,li

          do jj=1,nion
          if (tcount.eq.taskid) then
            ja     = ion_katm(jj)
            mult_lj = paw_basis_mult_l(ja)
            jl     = int_mb(i_v_mult(1)+jj-1)
            jll    = 0
            do lj=0,mult_lj
            do mj=-lj,lj

               indx = (ii-1)
     >              + ill*nion
     >              + (jj-1)*nion*lmax2
     >              + jll*nion*lmax2*nion

               e_mult
     >          = e_mult
     >          + 0.5d0*dble(
     >             dcpl_mb(paw_qlm(1)+jl)
     >            *dcpl_mb(paw_qlm(1)+il)
     >            *dcpl_mb(mult_energy_coeff(1)+indx))

              jl  = jl+1
              jll = jll+1
            end do !*mj*
            end do !*lj*
          end if
          tcount=mod(tcount+1,np)
          end do !*jj*

          il  = il +1
          ill = ill +1
        end do !*mi*
        end do !*li*
      end do !*ii*
      if (np.gt.1) call D3dB_SumAll(e_mult)
      call nwpw_timing_end(14)

      paw_mult_energy_atom_mult = e_mult
      return
      end

