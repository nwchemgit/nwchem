*
* $Id$
*


!**************************************************
!
!       Name: paw_kin_atom_energy
!
!       Purpose: calculates atom portion of the kinetic energy
!
!       Created: Mon Feb 10 13:46:23 PST 2003       
!**************************************************
      function paw_energy_kin_atom()
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_matrix_kin_data.fh"
#include "paw_basis.fh"
#include "paw_ovlp_data.fh"

      double precision  paw_energy_kin_atom
      integer in,ia
      integer i_w,i_mtr,i_mtr0
      integer i,j,li,lj,mi,mj
      integer nb
      double precision ekin

      ekin = 0.D0
      i_w = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        nb = paw_basis_nbasis(ia)
        i_mtr0 = int_mb(i_paw_kin_matrix(1) + ia - 1)
        ekin = ekin + paw_basis_core_kin_energy(ia)
        do i=1,nb
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
            do j=1,nb
              lj = paw_basis_orb_l(j,ia)
                do mj=-lj,lj
                i_mtr = i_mtr0 +i+(j-1)*nb
                if(mi.eq.mj) then
                  ekin = ekin + dcpl_mb(paw_ovlp_w(1)+i_w)*
     >                          dbl_mb(paw_kin_matrix(1) + i_mtr - 1) 
                end if
                i_w = i_w +1
                end do
             end do
           end do   
         end do
      end do
 
      paw_energy_kin_atom = ekin

      end

 
