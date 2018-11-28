C     
C     Check that rt-tddft parameters are legal.
C     
      subroutine rt_tddft_init_check_params (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t) params

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_init_check_params: "
      double precision, parameter :: max_tol = 1d-2
      integer, parameter :: max_terms = 100


C     == Variables ==
      double precision tol_schwartz


C     (now hardcoded to be -5*dt)
C      if (params%tmin.lt.0d0)
C     $     call errquit(pname//"tmin must be > 0", 0, 0)
      
      if (params%tmax.lt.0d0)
     $     call errquit(pname//"tmax must be > 0", 0, 0)
      
      if (params%dt.lt.0d0)
     $     call errquit(pname//"dt must be > 0", 0, 0)
      
      if (params%tmax.lt.params%tmin)
     $     call errquit(pname//"tmax must be > tmin", 0, 0)
      
      call rt_tddft_init_check_fields (params)

      if   ((params%prop_method.ne.1.).and.
     $     (params%prop_method.ne.2.).and.
     $     (params%prop_method.ne.3.).and.
     $     (params%prop_method.ne.4.))
     $     call errquit (pname//"bad prop_method", 0, 0)


      if ( (params%checklvl .ne. 1).and.
     $     (params%checklvl .ne. 2).and.
     $     (params%checklvl .ne. 3)) call errquit(pname//
     $     "checklvl must be 1, 2, 3", 0, 0)

      

C     
C     Check tolerances, etc
C     

      if ( (params%tol_zero .lt. 0d0).or.(params%tol_zero .ge. max_tol))
     $     call errquit (pname//"bad tol_zero", 0, 0)

      if ( (params%tol_series .lt. 0d0).or.
     $     (params%tol_series .gt. max_tol))
     $     call errquit (pname//"bad tol_series", 0, 0)

      if ( (params%tol_interpol .lt. 0d0).or.
     $     (params%tol_interpol .gt. max_tol))
     $     call errquit (pname//"bad tol_interpol", 0, 0)

      if ( (params%terms_series .lt. 1).or.
     $     (params%terms_series .gt. max_terms))
     $     call errquit (pname//"bad terms_series", 0, 0)

      if ( (params%terms_interpol .lt. 1).or.
     $     (params%terms_interpol .gt. max_terms))
     $     call errquit (pname//"bad terms_interpol", 0, 0)


      tol_schwartz = 10d0**(-itol2e) 

      if (tol_schwartz .gt. params%tol_interpol) then
         call rt_tddft_print_warning ("Schwartz screening tolerance "//
     $        "is looser than Magnus interpolation tolerance--"//
     $        "convergence will be difficult.  You should loosen "//
     $        "tol_interpol.")
         
      elseif (10d0*tol_schwartz .gt. params%tol_interpol) then
         call rt_tddft_print_warning ("Schwartz screening tolerance "//
     $        "is similar to the Magnus interpolation tolerance--"//
     $        "this may cause problems with convergence.  If so, "//
     $        "try loosening tol_interpol.")

      endif

      end subroutine



C====================================================================
C
C     Check that all fields have valid parameters
C
      subroutine rt_tddft_init_check_fields (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t) params


C     == Parameters ==
      character(len=*), parameter :: pname =
     $     "rt_tddft_init_check_fields: "


C     == Variables ==
      integer ifield
      type(rt_field_t) field

      do ifield = 1, params%nfields
         call rt_tddft_field_rtdb_get (params%rtdb, ifield, field)

         if ( (field%polarization .ne. "x").and.
     $        (field%polarization .ne. "y").and.
     $        (field%polarization .ne. "z").and.
     $        (field%polarization .ne. "angle"))
     $        call errquit (pname//"invalid polarization: "//
     $        field%polarization, 0, 0)

C     (actually its ok to have a negative "max")
C         if (field%max .lt. 0d0) call errquit (pname//
C     $        "cannot have negative field max", 0, 0)
 
C     (its ok to have a negative freq, too)
C         if ( (field%type .eq. "cw").or.(field%type.eq."pulse") ) then
C            if (field%frequency .lt. 0d0) call errquit (pname//
C     $           "cannot have negative frequency",0,0)
C         endif

         if (field%type.eq."pulse") then
            if (field%center .lt. params%tmin) call errquit (pname//
     $           "pulse center must be after tmin",0,0)
            if (field%width .lt. 0d0) call errquit (pname//
     $           "pulse width must be positive",0,0)
         endif

         if (params%theory.eq."CS") then
            if (field%spin .ne. "t") call errquit (pname//
     $           "fields can only act on total spin " //
     $           "for a closed shell calculation: "//field%spin, 0, 0)
         endif
         
      enddo

CXXX  [KAL]: check that field center is before tmax      

      end subroutine
c $Id$
