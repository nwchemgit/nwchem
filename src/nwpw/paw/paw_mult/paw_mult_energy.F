*
* $Id: paw_mult_energy.F,v 1.2 2003-02-23 22:45:52 bylaska Exp $
*
*
* $Log: not supported by cvs2svn $
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

      subroutine paw_mult_energy_atom_comp(e_comp)
      implicit none
      real*8 e_comp(*)

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer ii,ia,li,mi,il,mult_li
      integer nion

      nion=ion_nion()
      call dcopy(nion,0.0d0,0,e_comp,1)

      do ii=1,nion
        ia      = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il      = int_mb(i_v_mult(1)+ii-1)
        do li=0,mult_li
        do mi=-li,li

          e_comp(ii) = e_comp(ii) 
     >               + dble(dcpl_mb(paw_qlm(1)+il)
     >                      *dcpl_mb(comp_coeff(1)+il))
          il = il + 1
        end do !*mi*
        end do !*li*
      end do !*ii*

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

      subroutine paw_mult_energy_atom_self(e_self)
      implicit none
      real*8 e_self(*)

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer ii,ia,li,mi,il,mult_li
      integer indx,lmax,nion

      nion=ion_nion()
      lmax=paw_basis_max_mult_l()

      call dcopy(nion,0.0d0,0,e_self,1)
      do ii=1,nion
        ia     = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il = int_mb(i_v_mult(1)+ii-1)
        do li=0,mult_li
          indx = li + (ia-1)*(lmax+1)
          do mi=-li,li

          e_self(ii) 
     >     = e_self(ii)
     >     - 0.5d0*dble(
     >             dconjg(dcpl_mb(paw_qlm(1)+il))
     >            *dcpl_mb(paw_qlm(1)+il))
     >            *dbl_mb(self_energy_coeff(1)+indx)
          il = il +1
          end do !*mi*
        end do !*li*
      end do !*ii*

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

      subroutine paw_mult_energy_atom_mult(e_mult)
      implicit none
      real*8 e_mult

#include "mafdecls.fh"
#include "paw_geom.fh"
#include "paw_mult_data.fh"
#include "paw_comp_charge_data.fh"
#include "paw_basis.fh"

      !*** local variables ***
      integer ii,ia,li,mi,il,ill,mult_li
      integer jj,ja,lj,mj,jl,jll,mult_lj
      integer nion,indx,lmax,lmax2

      nion  = ion_nion()
      lmax  = paw_basis_max_mult_l()
      lmax2 = (lmax+1)**2

      e_mult = 0.0d0
      do ii=1,nion
        ia     = ion_katm(ii)
        mult_li = paw_basis_mult_l(ia)
        il     = int_mb(i_v_mult(1)+ii-1)
        ill    = 0
        do li=0,mult_li
        do mi=-li,li

          do jj=1,nion
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
          end do !*jj*

          il  = il +1
          ill = ill +1
        end do !*mi*
        end do !*li*
      end do !*ii*

      return
      end

