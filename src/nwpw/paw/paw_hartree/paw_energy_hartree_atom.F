*
* $Id$
*


!**************************************************
!
!       Name: paw_hartree_atom_energy
!
!       Purpose: calculates atom portion of the hartree energy
!
!       Created: Mon Feb 10 13:46:23 PST 2003       
!**************************************************
      function paw_energy_hartree_atom()
      implicit none

#include "bafdecls.fh"
#include "paw_geom.fh"
#include "paw_ma.fh"
#include "paw_pot_hartree_data.fh"
#include "paw_basis.fh"
#include "paw_proj.fh"
#include "paw_ovlp_data.fh"

      double precision  paw_energy_hartree_atom
      integer in,ia
      integer i_w
      integer i_hp0,i_hp
      integer i,j,li,lj,mi,mj
      integer nilm,njlm
      integer nb
      double precision ehartree

      
      ehartree = 0.d0
      i_w = 0
      do in=1,ion_nion()
        ia = ion_katm(in)
        i_hp0 = int_mb(i_paw_pot_hartree(1) + in - 1)
        nb = paw_basis_nbasis(ia)
        nilm = 0
        do i=1,nb
          li = paw_basis_orb_l(i,ia)
          do mi=-li,li
            njlm = 0
            do j=1,nb
              lj = paw_basis_orb_l(j,ia)
                do mj=-lj,lj
                 i_hp  = i_hp0-1+(njlm+lj+mj+1)+
     >                   (nilm+li+mi)*paw_proj_nbasis(ia)

                  ehartree = ehartree + dble(dcpl_mb(paw_ovlp_w(1)+i_w)*
     >                         dcpl_mb(paw_pot_hartree(1)+i_hp) )


                  i_w = i_w + 1
                end do !mj
                njlm=njlm + 2*lj+1
             end do !j
           end do !mi   
           nilm=nilm + 2*li+1
         end do !j
      end do
  
      paw_energy_hartree_atom = 0.5d0*ehartree

      end

 
