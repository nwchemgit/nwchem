C
C     Initializes dipole and quadrupole integrals.
C
      subroutine rt_tddft_init_dip_quad_ints (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_dip_quad_ints: "

      
C     == Variables ==
      integer n
      logical ok


C
C     All matricies are nbf_ao x nbf_ao *not* ns_ao.
C
      n = params%nbf_ao
      

C
C     Dipole integrals.
C
      ok = .true.

      ok=ok.and.ga_create(mt_dbl,n,n,"D_x", 0,0,
     $     params%g_dip%x)

      ok=ok.and.ga_create(mt_dbl,n,n,"D_y", 0,0,
     $     params%g_dip%y)

      ok=ok.and.ga_create(mt_dbl,n,n,"D_z", 0,0,
     $     params%g_dip%z)

      if (.not.ok) call errquit(pname//"alloc dip ints failed",0,GA_ERR)

      
C     (for now make excite and measure dipole mats both full system)
      call int_dip_ga (ao_bas_han, ao_bas_han, params%g_dip%x,
     $     params%g_dip%y, params%g_dip%z)



C     XXX
c$$$      call rt_tddft_ao_mask (params, params%g_dip%x,
c$$$     $     1, params%g_dip%x)
c$$$      call rt_tddft_ao_mask (params, params%g_dip%y,
c$$$     $     1, params%g_dip%y)
c$$$      call rt_tddft_ao_mask (params, params%g_dip%z,
c$$$     $     1, params%g_dip%z)
c$$$      if (ga_nodeid().eq.0) then
c$$$         write (luout, *) "XXX MASKED G_DIP XXX"
c$$$      endif

      
C
C     Quadrupole integrals.
C
      ok = .true.
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xx", 0,0, params%g_quad%xx)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xy", 0,0, params%g_quad%xy)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xz", 0,0, params%g_quad%xz)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_yy", 0,0, params%g_quad%yy)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_yz", 0,0, params%g_quad%yz)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_zz", 0,0, params%g_quad%zz)

      if (.not.ok) call errquit(pname//"alloc quad ints fail",0,GA_ERR)

      call int_qdr_ga (ao_bas_han, ao_bas_han,
     $     params%g_quad%xx, params%g_quad%xy, params%g_quad%xz,
     $     params%g_quad%yy, params%g_quad%yz, params%g_quad%zz)


C     (force quad integrals to be symmetric)
      call ga_symmetrize (params%g_quad%xx)
      call ga_symmetrize (params%g_quad%xy)
      call ga_symmetrize (params%g_quad%xz)
      call ga_symmetrize (params%g_quad%yy)
      call ga_symmetrize (params%g_quad%yz)
      call ga_symmetrize (params%g_quad%zz)

      end subroutine rt_tddft_init_dip_quad_ints

      




#if 0
C
C     Initializes dipole and quadrupole integrals for closed shell or
C     open shell system.
C
      subroutine rt_tddft_init_dip_quad_ints_new (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_dip_quad_ints: "

      
C     == Variables ==
      integer n
      logical ok


      n = params%ns_ao            !alias for clarity

C
C     Dipole integrals.
C
      ok = .true.


C      XXXX MAKE ALLOCATION ROUTINE


      ok=ok.and.ga_create(mt_dbl,n,n,"D_x", 0,0,
     $     params%g_dip%x)

      ok=ok.and.ga_create(mt_dbl,n,n,"D_y", 0,0,
     $     params%g_dip%y)

      ok=ok.and.ga_create(mt_dbl,n,n,"D_z", 0,0,
     $     params%g_dip%z)

      if (.not.ok) call errquit(pname//"alloc dip ints failed",0,GA_ERR)

      
C     (for now make excite and measure dipole mats both full system)
      call int_dip_ga (ao_bas_han, ao_bas_han, params%g_dip%x,
     $     params%g_dip%y, params%g_dip%z)


      
C
C     Quadrupole integrals.
C
      ok = .true.
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xx", 0,0, params%g_quad%xx)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xy", 0,0, params%g_quad%xy)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_xz", 0,0, params%g_quad%xz)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_yy", 0,0, params%g_quad%yy)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_yz", 0,0, params%g_quad%yz)
      ok=ok.and.ga_create(mt_dbl,n,n,"Q_zz", 0,0, params%g_quad%zz)

      if (.not.ok) call errquit(pname//"alloc quad ints fail",0,GA_ERR)

      call int_qdr_ga (ao_bas_han, ao_bas_han,
     $     params%g_quad%xx, params%g_quad%xy, params%g_quad%xz,
     $     params%g_quad%yy, params%g_quad%yz, params%g_quad%zz)


C     (force quad integrals to be symmetric)
      call ga_symmetrize (params%g_quad%xx)
      call ga_symmetrize (params%g_quad%xy)
      call ga_symmetrize (params%g_quad%xz)
      call ga_symmetrize (params%g_quad%yy)
      call ga_symmetrize (params%g_quad%yz)
      call ga_symmetrize (params%g_quad%zz)

      end subroutine rt_tddft_init_dip_quad_ints_new
#endif
    
c $Id$
