C
C     rt_tddft_calc.F
C
C     Contains various simple calculation routines which are universal
C     to open shell, closed shell, etc.
C

C====================================================================
C
C     Sums the total energy of the system.
C     
      subroutine rt_tddft_calc_Etot (params, E)
      implicit none

#include "rt_tddft.fh"

C     == In/out ==
      type(rt_params_t) params
      type(rt_energies_t) E

      E%tot = E%core + E%coul + E%xc(1) + E%xc(2) + params%Enuc
      
      end subroutine
      


C====================================================================
C
C     Calculates the dipole moment of system from the real part of the
C     density matrix in AO basis.
C
C     (parts ripped from rt_tddft_utils.F)
C
C     If supplied the total dens mat, will return total dip moment.  Can
C     also pass alpha, beta matrices separately and it will return the
C     corresponding dipole moments for these spins.
C     
C     
      subroutine rt_tddft_calc_dip (params, g_densre, dipmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  ! real part of complex dens mat in AO basis

      
C     == Outputs ==
      type(rt_vector_t), intent(out) :: dipmom

      integer geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments)
      double precision elecmoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      integer i
C      data center/3*0.0d0/
      character(*), parameter :: pname = "rt_tddft_calc_dip: "


C     == Checks ==
C      if ( (abs(params%origin(1)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(2)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(3)) .ge. params%tol_zero) )
C     $     call errquit (pname//"system origin not at 0,0,0",0,0)

      if (.not. bas_cando_mpoles(params%ao_bas_han))
     &     call errquit(pname//'cant do mpoles',0, BASIS_ERR)

      if (.not. bas_geom(params%ao_bas_han, geom)) 
     &     call errquit(pname//'bas_geom failed',0, BASIS_ERR)
 
C     (used to be hardcoded 0,0,0)
      center(1) = params%origin(1)
      center(2) = params%origin(2)
      center(3) = params%origin(3)

c     Initialisation
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, elecmoments, 1)
      
c     Calculate the multipole components of the density and nucleus
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(params%ao_bas_han, center, lmax, g_densre, 
     &     elecmoments, len_moments)
      
      do k = 1, len_moments
         totalmoments(k) = elecmoments(k) + nuclmoments(k)
      enddo

C     return x,y,z components of 2nd moment, i.e. dipole moment
      dipmom%x = totalmoments(2)
      dipmom%y = totalmoments(3)
      dipmom%z = totalmoments(4)

      end subroutine



C====================================================================
C
C     Computes electric dipole moment of system on geometry fragment
C     "ig" via density matrix in AO basis--P(t).
C
C    \mu_e (t) = Tr [ D P(t) ]
C
C     where we first project P(t) onto the fragment via a mask.
C     
      subroutine rt_tddft_calc_dip_geom (params, g_densre, ig, dipmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: ig         ! index of geometry to calculate on
      integer, intent(in)           :: g_densre   ! real part of complex dens mat in AO basis

      
C     == Outputs ==
      type(rt_vector_t), intent(out) :: dipmom


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_dip_geom: "


C     == Variables ==
      integer g_dip_elem    !, g_densre_mask
      integer g_dip_maskd
      

      if (.not. ga_duplicate (params%g_dip%x, g_dip_elem,
     $     "dip_elem"))
     $     call errquit (pname//"failed to alloc g_dip_elem",0,GA_ERR)

      if (.not. ga_duplicate (g_densre, g_dip_maskd,
     $     "dip_maskd"))
     $     call errquit(pname//"failed to alloc g_dip_maskd",0,GA_ERR)


      
      if ( (ig .lt. 1).or.(ig .gt. params%ngeoms) )
     $     call errquit (pname//"bad ig: ", ig, 0)


C     (x dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%x, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%x)
      dipmom%x = params%geom(ig)%ndip%x - dipmom%x

      
C     (y dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%y, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%y)
      dipmom%y = params%geom(ig)%ndip%y - dipmom%y

      
C     (z dipole moment)
      call ga_zero (g_dip_elem)
      call rt_tddft_geom_mask (params, params%g_dip%z, ig, g_dip_maskd)
      call ga_dgemm ("N", "N", params%ns_ao, params%ns_ao, params%ns_ao,
     $     1d0, g_dip_maskd, g_densre, 0d0, g_dip_elem)
      call mat_dtrace (g_dip_elem, dipmom%z)
      dipmom%z = params%geom(ig)%ndip%z - dipmom%z


      if (.not. ga_destroy (g_dip_elem))
     $     call errquit (pname//"failed to destroy g_dip_elem",0,GA_ERR)

      if (.not. ga_destroy (g_dip_maskd))
     $     call errquit (pname//"failed to destroy g_dip_maskd",
     $     0,GA_ERR)

      end subroutine





C====================================================================
C
C     Calculates the quadrupole moment of system from the real part of
C     the density matrix in AO basis.
C
C     If supplied the total dens mat, will return total dip moment.  Can
C     also pass alpha, beta matrices separately and it will return the
C     corresponding dipole moments for these spins.
C     
C     
      subroutine rt_tddft_calc_quad (params, g_densre, qmom)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "geom.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  ! real part of complex dens mat in AO basis

C     == Outputs ==
      type(rt_quad_t), intent(out) :: qmom

      integer  geom
      integer len_moments, lmax
      parameter (lmax = 2)
      parameter (len_moments = (lmax+1)*(lmax+2)*(lmax+3)/6)
      double precision totalmoments(len_moments)
      double precision elecmoments (len_moments)
      double precision nuclmoments (len_moments)
      double precision center(3)
      integer k, l, m, ltotal, ind
      integer i
C      data center/3*0.0d0/
      character(*), parameter :: pname = "rt_tddft_calc_quad: "

C     == Checks ==
C      if ( (abs(params%origin(1)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(2)) .ge. params%tol_zero) .or.
C     $     (abs(params%origin(3)) .ge. params%tol_zero) )
C     $     call errquit (pname//"system origin not at 0,0,0",0,0)

      if (.not. bas_cando_mpoles(params%ao_bas_han))
     &     call errquit(pname//'cant do mpoles',0, BASIS_ERR)

      if (.not. bas_geom(params%ao_bas_han, geom)) 
     &     call errquit(pname//'bas_geom failed',0, BASIS_ERR)

C     (used to be hardcoded 0,0,0)
      center(1) = params%origin(1)
      center(2) = params%origin(2)
      center(3) = params%origin(3)
 
c     Initialisation
      call dfill(len_moments, 0.0d0, totalmoments, 1)
      call dfill(len_moments, 0.0d0, nuclmoments, 1)
      call dfill(len_moments, 0.0d0, elecmoments, 1)
      
c     Calculate the multipole components of the density and nucleus
      call geom_nuc_mpole(geom, center, lmax, nuclmoments, len_moments)
      call multipole_density(params%ao_bas_han, center, lmax, g_densre, 
     &     elecmoments, len_moments)
      
      do k = 1, len_moments
         totalmoments(k) = elecmoments(k) + nuclmoments(k)
      enddo

C
C     return xx,xy,xz,yy,yz,zz components 
C
C     recall:
C     totalmoments(1) = monopole,
C     totalmoments(2:4) = dipole
C     totalmoments(5:10) = quadrupoles
C
C
      qmom%xx = totalmoments(5)
      qmom%xy = totalmoments(6)
      qmom%xz = totalmoments(7)
      qmom%yy = totalmoments(8)
      qmom%yz = totalmoments(9)
      qmom%zz = totalmoments(10)
      
      end subroutine



C====================================================================
C
C     Computes electronic charge on geometry fragment "i":
C
C     q_i = - e_i = Tr [ Re[P] S ]_i,
C     
C     where the trace is projected onto fragment via the mask.  It does
C     not subtract the nuclear charge.
C
C
      subroutine rt_tddft_calc_charge_geom(params, g_densre, ig, charge)
      implicit none

#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre  !re part of dens mat AO bas, nbf_ao x nbf_ao
      integer, intent(in)           :: ig

      
C     == Outputs ==
      double precision, intent(out) :: charge   !electronic charge on this fragment


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_charge_geom: "


C     == Variables ==
      integer n
      integer g_dscr1, g_dscr2

      
      if (.not. ga_duplicate (g_densre, g_dscr1,
     $     "dscr1"))
     $     call errquit (pname//"failed to alloc g_dscr1",0,GA_ERR)

      if (.not. ga_duplicate (g_densre, g_dscr2,
     $     "dscr2"))
     $     call errquit (pname//"failed to alloc g_dscr2",0,GA_ERR)

      
      call ga_zero (g_dscr1)
      call ga_zero (g_dscr2)
      
      n = params%ns_ao

      
C     (note negative to make electronic charge -ve)
      call ga_dgemm ("N", "N", n, n, n, -1d0, g_densre, params%g_s,
     $     0d0, g_dscr1)
      
      call rt_tddft_geom_mask (params, g_dscr1, ig, g_dscr2)
 
      call mat_dtrace (g_dscr2, charge)

      if (.not. ga_destroy (g_dscr1))
     $     call errquit (pname//"failed to destroy g_dscr1",0,GA_ERR)
      
      if (.not. ga_destroy (g_dscr2))
     $     call errquit (pname//"failed to destroy g_dscr2",0,GA_ERR)

      end subroutine



C====================================================================
C
C     Computes <S^2> value (for entire active geometry).
C     
C     <S^2 (t)> = [(Na - Nb)/2]^2 + (Na + Nb)/2 - Tr [P'a(t).P'b(t)]
C
C     (P' is real part of dens mat in MO basis)h
C
      subroutine rt_tddft_calc_s2 (params, g_densre_mo, s2val)
      implicit none
      
#include "errquit.fh"
#include "util.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      integer, intent(in)           :: g_densre_mo(2)! real part dens mat; alpha,beta

      
C     == Outputs ==
      double precision, intent(out) :: s2val


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_calc_s2: "


C     == Variables ==
      integer          :: g_pp
      double precision :: trpp


      if (.not. ga_duplicate (g_densre_mo(1), g_pp,
     $     "dscr"))
     $     call errquit (pname//"failed to alloc g_pp",0,GA_ERR)

      call ga_zero (g_pp)

C
C     Note params%nbf_mo size since if we are doing spin-orbit
C     calculation, we will send the "alpha" and "beta" parts separately.
C
      call ga_dgemm ("N", "N",
     $     params%nbf_mo, params%nbf_mo, params%nbf_mo,
     $     1d0, g_densre_mo(1), g_densre_mo(2), 0d0, g_pp)

      call mat_dtrace (g_pp, trpp)

      s2val = 0.25d0*(params%nalpha - params%nbeta)**2
     $     + 0.5d0*(params%nalpha + params%nbeta)
     $     - trpp

      if (.not. ga_destroy (g_pp))
     $     call errquit (pname//"destroy failed", 0, 0)

      end subroutine
c $Id$
