*
* $Id$
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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
      real*8 xxx

c      complex*16 q1,q2,c12
c      integer indx1,indx2,indx3,indx4

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

c          if (ii.eq.jj) then
c               etest2
c     >          = etest2
c     >          + 0.5d0*dble(
c     >             dcpl_mb(paw_qlm(1)+jl)
c     >            *dcpl_mb(paw_qlm(1)+il)
c     >            *dcpl_mb(mult_energy_coeff(1)+indx))
c          end if
c
c          if (mi.lt.0) then
c             q1 = dcpl_mb(paw_qlm(1)+il+2*abs(mi))
c             q1 =  (-1)**abs(mi)*q1
c             q1 = q1 - dcpl_mb(paw_qlm(1)+il)
c             q1 = q1*dcmplx(0.0d0,1.0d0)/dsqrt(2.0d0)
c          else if (mi.gt.0) then
c             q1 = dcpl_mb(paw_qlm(1)+il)
c             q1 =  (-1)**abs(mi)*q1
c             q1 = dcpl_mb(paw_qlm(1)+il-2*abs(mi)) + q1
c             q1 = q1/dsqrt(2.0d0)
c          else
c             q1 = dble(dcpl_mb(paw_qlm(1)+il))
c          end if
c
c          if (mj.lt.0) then
c             q2 = dcpl_mb(paw_qlm(1)+jl+2*abs(mj))
c             q2 =  (-1)**abs(mj)*q2
c             q2 = q2 - dcpl_mb(paw_qlm(1)+jl)
c             q2 = q2*dcmplx(0.0d0,1.0d0)/dsqrt(2.0d0)
c          else if (mj.gt.0) then
c             q2 = dcpl_mb(paw_qlm(1)+jl)
c             q2 =  (-1)**abs(mj)*q2
c             q2 = dcpl_mb(paw_qlm(1)+jl-2*abs(mj)) + q2
c             q2 = q2/dsqrt(2.0d0)
c          else
c             q2 = dble(dcpl_mb(paw_qlm(1)+jl))
c          end if
c
c          indx1 = indx 
c          indx2 = indx - 2*mi*nion
c          indx3 = indx - 2*mj*nion*lmax2*nion
c          indx4 = indx - 2*mi*nion - 2*mj*nion*lmax2*nion
c          if ((mi.gt.0).and.(mj.gt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)*(-1)**abs(mi+mj)
c             c12 = c12 
c     >           + dcpl_mb(mult_energy_coeff(1)+indx2)*(-1)**abs(mj)
c             c12 = c12 
c     >           + dcpl_mb(mult_energy_coeff(1)+indx3)*(-1)**abs(mi)
c             c12 = c12
c     >           + dcpl_mb(mult_energy_coeff(1)+indx4)
c              c12 = 0.5d0*c12
c          else if ((mi.eq.0).and.(mj.gt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)*(-1)**abs(mj)
c             c12 = c12 
c     >           + dcpl_mb(mult_energy_coeff(1)+indx3)
c             c12 = c12/dsqrt(2.0d0)
c          else if ((mi.lt.0).and.(mj.gt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)*(-1)**abs(mj)
c             c12 = c12 
c     >           - dcpl_mb(mult_energy_coeff(1)+indx2)*(-1)**abs(mi+mj)
c             c12 = c12 
c     >           + dcpl_mb(mult_energy_coeff(1)+indx3)
c             c12 = c12 
c     >           - dcpl_mb(mult_energy_coeff(1)+indx4)*(-1)**abs(mi)
c             c12 = dcmplx(0.0d0,0.5d0)*c12
c
c          else if ((mi.lt.0).and.(mj.eq.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx2)*(-1)**abs(mi)
c             c12 = c12*dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
c          else if ((mi.eq.0).and.(mj.eq.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)
c          else if ((mi.lt.0).and.(mj.eq.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx2)*(-1)**abs(mi)
c             c12 = c12*dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
c
c          else if ((mi.gt.0).and.(mj.lt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)*(-1)**abs(mi)
c             c12 = c12
c     >           + dcpl_mb(mult_energy_coeff(1)+indx2)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx3)*(-1)**abs(mi+mj)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx4)*(-1)**abs(mj)
c             c12 = c12*dcmplx(0.0d0,0.5d0)
c          else if ((mi.eq.0).and.(mj.lt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx3)*(-1)**abs(mj)
c             c12 = c12*dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
c          else if ((mi.lt.0).and.(mj.lt.0)) then
c             c12 = dcpl_mb(mult_energy_coeff(1)+indx1)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx2)*(-1)**abs(mi)
c             c12 = c12
c     >           - dcpl_mb(mult_energy_coeff(1)+indx3)*(-1)**abs(mj)
c             c12 = c12
c     >           + dcpl_mb(mult_energy_coeff(1)+indx4)*(-1)**abs(mi+mj)
c             c12 = -0.5d0*c12
c          end if
c
c          write(23,'(A,6I3,E15.6,6E14.6)') 
c     >           "paw emult=",ii,li,mi,jj,lj,mj,
c     >            0.5d0*dble(
c     >             dcpl_mb(paw_qlm(1)+jl)
c     >            *dcpl_mb(paw_qlm(1)+il)
c     >            *dcpl_mb(mult_energy_coeff(1)+indx)),
c     >             q1,q2,c12
c
c          if (ii.eq.jj) then
c          write(33,'(A,6I3,E15.6,6E14.6)') 
c     >           "paw emult=",ii,li,mi,jj,lj,mj,
c     >            0.5d0*dble(
c     >             dcpl_mb(paw_qlm(1)+jl)
c     >            *dcpl_mb(paw_qlm(1)+il)
c     >            *dcpl_mb(mult_energy_coeff(1)+indx)),
c     >             q1,q2,c12*0.5d0
c          end if

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

