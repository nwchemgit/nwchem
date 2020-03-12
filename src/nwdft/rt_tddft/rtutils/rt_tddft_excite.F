C
C     rt_tddft_excite.F
C
C     Routines for computing fields, dipole interaction potentials, etc.
C     


C====================================================================
C
C    Compute the applied field potential acting on spin "spin", and also
C    updates the total applied field on each geometry.
C    
C
      subroutine rt_tddft_calc_excite (params, tt, g_zpot)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      double precision, intent(in)  :: tt
C      character, intent(in)         :: spin !"a"-alpha, "b"-beta


C     == Outputs ==
      integer, intent(in) :: g_zpot(2) !total complex pot due to ALL applied fields; alph,beta

      
C     == Parameters ==
      character(*), parameter     :: pname = "rt_tddft_excite: "


C     == Variables ==
      integer iexcite, ifield, igeom
      integer g_scr, g_dip_maskd, g_zscr
      double precision elapsed
      integer is
      integer iao, jao
      integer imo, jmo


      if (params%prof) call prof_start (elapsed)
      

C      if ( (spin.ne."a").and.(spin.ne."b") )
C     $     call errquit (pname//"invalid spin: "//spin, 0, 0)


C
C     Initialze GA's for masked potential matrix.
C
      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "scr", 0, 0, g_scr))
     $     call errquit ("couldnt create scr", 0, GA_ERR)

      if (.not. ga_create(mt_dbl, params%ns_ao, params%ns_ao,
     $     "dip_maskd", 0, 0, g_dip_maskd))
     $     call errquit ("couldnt create dip_maskd", 0, GA_ERR)

      if (.not. ga_create(mt_dcpl, params%ns_ao, params%ns_ao,
     $     "zscr", 0, 0, g_zscr))
     $     call errquit ("couldnt create zscr", 0, GA_ERR)

      call ga_zero (g_scr)
      call ga_zero (g_zscr)


C
C     Update the values of all the fields.  Doesnt really matter how
C     many times we call this...
C
CXXX  [KAL]: move to propagation routines?
      call rt_tddft_update_fields (params, tt)

      
C
C     Loop over the excitation rules (will do nothing if no exciation
C     rules set up), and compute the total applied field (from all
C     fields) on each geom.  Note, the values of the fields at the
C     current time should already be stored in the corresponding
C     structs.
C     
C     If this field acts on this spin, or if it acts on total (both),
C     then add it.
C
      do igeom = 1, params%ngeoms
         call rt_tddft_vec_zero (params%geom(igeom)%vec(1)) !alpha
         call rt_tddft_vec_zero (params%geom(igeom)%vec(2)) !beta
         call rt_tddft_quad_zero (params%geom(igeom)%quad)  !xxx make alpha beta quads
      enddo

      
      do iexcite = 1, params%nexcites
         ifield = params%excite(iexcite)%field_indx
         igeom = params%excite(iexcite)%geom_indx
         
C         if ( (params%field(ifield)%spin .eq. spin) .or.
C     $        (params%field(ifield)%spin .eq. "t") ) then
C            call rt_tddft_vec_add (params%geom(igeom)%vec,
C     $           params%field(ifield)%vec, params%geom(igeom)%vec)

         call rt_tddft_vec_add (params%geom(igeom)%vec(1),
     $        params%field(ifield)%vec(1), params%geom(igeom)%vec(1))
         call rt_tddft_vec_add (params%geom(igeom)%vec(2),
     $        params%field(ifield)%vec(2), params%geom(igeom)%vec(2))
         
      enddo


C
C     Now compute D.Etot for each geom, and accumulate in g_zpot.  We do
C     this by masking the D matrix for each geometry, and computing
C     Dmask.Etot on each geom.  This is valid because:
C
C     Vappl_x = - D.E_1 - D.E_2 - ... - D.E_N, 
C           = - D.(E_1 + E_2 + ... + E_N)
C
C     We add both the alpha and beta spin parts.
C
      call ga_zero (g_zpot(1))
      call ga_zero (g_zpot(2))
         
      iao = 1
      jao = params%nbf_ao

      imo = 1
      jmo = 1 !params%nbf_mo
      
      do igeom = 1, params%ngeoms
         do is = 1, 2
            call ga_zero (g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%x,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           3, 3, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%x, g_dip_maskd,
     $           1d0, g_scr, g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%y,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           3, 3, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%y, g_dip_maskd,
     $           1d0, g_scr, g_scr)

            call rt_tddft_geom_mask (params, params%g_dip%z,
     $           igeom, g_dip_maskd)
c$$$            call rt_tddft_aomo_mask (params, g_dip_maskd,
c$$$     $           imo, jmo, g_dip_maskd)
c$$$            call rt_tddft_ao_mask (params, g_dip_maskd,
c$$$     $           1, 1, g_dip_maskd)
            call ga_dadd (params%geom(igeom)%vec(is)%z, g_dip_maskd,
     $           1d0, g_scr, g_scr)

         
C     
C     Cast resulting potential to complex so we can add to Fock mat, and
C     accumulate in g_zpot.  Note negative sign!
C     
            call ga_zero (g_zscr)
            call ga_scale (g_scr, -1d0)
            call convert_d2z (1d0, g_scr, 0d0, g_scr, g_zscr)
            call ga_zadd (z1, g_zscr, z1, g_zpot(is), g_zpot(is))
         enddo
      enddo
      
         
C
C     Clean up
C     
      if (.not. ga_destroy (g_scr))
     $     call errquit ("failed to destroy g_scr", 0, GA_ERR)
      
      if (.not. ga_destroy (g_dip_maskd))
     $     call errquit ("failed to destroy g_dip_maskd", 0, GA_ERR)
      
      if (.not. ga_destroy (g_zscr))
     $     call errquit ("failed to destroy g_zscr", 0, GA_ERR)

      if (params%prof) call prof_end (elapsed,
     $     "Fields and excitation potential")

      end subroutine 



C====================================================================
C
C     Computes Hann function for windowing pulse.
C
      double precision function rt_tddft_hann_window (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_hann_window: "

C     == Variables ==
      double precision tstart, tend, tshift, hout, arg

      tstart = field%center - field%width / 2d0
      tend = tstart + field%width
      
      tshift = tt - tstart

      if ((tt .gt. tend).or.(tt .lt. tstart)) then
         hout = 0d0
      else
         arg = dpi / field%width * tshift
         hout = sin(arg) * sin(arg)
      endif

      rt_tddft_hann_window = hout

      end function



C====================================================================
C
C     Computes Gaussian function for windowing pulse.
C
      double precision function rt_tddft_gaussian_window (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_gaussian_window: "

C     == Variables ==
      double precision arg, gout
      

      arg = -1d0 * (tt - field%center)*(tt - field%center) /
     $     (2d0*field%width*field%width)

c$$$            c = field%max / sqrt (2d0 * dpi) / s  !normalization

      gout = exp (arg)
      rt_tddft_gaussian_window = gout

      end function


C====================================================================
C
C     Computes sin^2 ramp up function for windowing applied field.
C
      double precision function rt_tddft_sin2ramp_window (field, tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"
      
C     == Inputs ==
      type(rt_field_t), intent(in)  :: field
      double precision, intent(in)  :: tt
      
C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_tddft_sin2ramp_window: "

C     == Variables ==
      double precision ww, wid, t0, t1, arg
      double precision ramp, ss

      ww = field%frequency
      t0 = field%start      !start of ramp up 
      wid = field%width     !field%ramp_wid
      t1 = t0 + wid

      if (tt.le.t0) then
         ramp = 0d0
      elseif ((tt.gt.t0) .and. (tt.lt.t1)) then
         ss = sin(dpi * (tt-t0) / (2d0*wid))
         ramp = ss*ss
      else
         ramp = 1d0
      endif
      
      rt_tddft_sin2ramp_window = ramp

      end function

      

C====================================================================
C
C     Compute the values of all the fields at this time and store values
C     in the corresponding structs.
C     
      subroutine rt_tddft_update_fields (params, tt)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      type(rt_params_t), intent(inout) :: params  !XXXINTENT
      double precision, intent(in)     :: tt

C     == Parameters ==
      character(len=*), parameter :: pname = "rt_tddft_update_fields: "

C     == External ==
      double precision, external :: rt_tddft_hann_window
      double precision, external :: rt_tddft_gaussian_window
      double precision, external :: rt_tddft_sin2ramp_window

C     == Variables ==
      integer ifield
      type(rt_field_t) field
      type(rt_vector_t) fvec
      double precision val
      double precision delta_time
      double precision :: theta, phi

      do ifield = 1, params%nfields
         field = params%field(ifield)

         if ( (field%polarization .ne. "x").and.
     $        (field%polarization .ne. "y").and.
     $        (field%polarization .ne. "z").and.
     $        (field%polarization .ne. "angle"))
     $        call errquit (pname//
     $        "only works for field values (no quads yet)", 0, 0)
         
         val = 0d0

         
C     (frequency-containing pulses, potentially windowed)
         if ( (field%type .eq. "cw") .or.
     $        (field%type .eq. "gaussian") .or.
     $        (field%type .eq. "hann").or.
     $        (field%type .eq. "sin2ramp")) then
            
            val = field%max * sin (field%frequency*tt + field%phase)
            
            if (field%type .eq. "hann") then
               val = val * rt_tddft_hann_window (field, tt)
            elseif (field%type .eq. "gaussian") then
               val = val * rt_tddft_gaussian_window (field, tt)
            elseif (field%type .eq. "sin2ramp") then
               val = val * rt_tddft_sin2ramp_window (field, tt)
            endif

         elseif (field%type .eq. "delta") then
            
            delta_time = field%center

            if ( abs(tt - delta_time) .lt. 0.5d0*params%dt) then
               val = field%max
            else
               val = 0d0
            endif

         else
            call errquit (pname//"bad excite type: "//field%type,0,0)
         endif


C
C     negative times never have excitations
C
         if (tt .lt. -0.1d0*params%dt) then
            val = 0d0
         endif

         
C     
C     Store the applied field value (vector) in the field struct, note
C     alpha, beta spin parts.
C
C     XXX TODO: Instead of cases, just use the theta,phi from the input
C         
         call rt_tddft_vec_zero (fvec)
         if (field%polarization .eq. "x") then
            fvec%x = val
         elseif (field%polarization .eq. "y") then
            fvec%y = val
         elseif (field%polarization .eq. "z") then
            fvec%z = val
         elseif (field%polarization .eq. "angle") then
            theta = field%theta
            phi = field%phi
c$$$            if (.not.rtdb_get(params%rtdb, "rt_tddft:efield:theta",
c$$$     $           mt_dbl, 1, theta))
c$$$     $           call errquit (pname//
c$$$     $           "couldnt read rt_tddft:efield:theta",0,0)
c$$$
c$$$            if (.not.rtdb_get(params%rtdb, "rt_tddft:efield:phi",
c$$$     $           mt_dbl, 1, phi))
c$$$     $           call errquit (pname//
c$$$     $           "couldnt read rt_tddft:efield:phi",0,0)

            if (ga_nodeid().eq.0) then
               write(luout, *) "Efield angles (theta, phi) = ",
     $              theta, phi
            endif

            theta = theta / deg2rad
            phi = phi / deg2rad
            
            fvec%x = val * sin(theta) * cos(phi) 
            fvec%y = val * sin(theta) * sin(phi)
            fvec%z = val * cos(theta)

            if (abs(fvec%x) > abs(val)) then
               call errquit("Ex should not be larger than Emax",0,0)
            endif

            if (abs(fvec%y) > abs(val)) then
               call errquit("Ey should not be larger than Emax",0,0)
            endif

            if (abs(fvec%z) > abs(val)) then
               call errquit("Ez should not be larger than Emax",0,0)
            endif
            
            
         else
            call errquit (pname//"bad polarization: "//
     $           field%polarization, 0, 0)
         endif

         call rt_tddft_vec_zero (field%vec(1)) !alpha part
         call rt_tddft_vec_zero (field%vec(2)) !beta part

         if ((field%spin .eq. "a") .or. (field%spin .eq. "t")) then
            field%vec(1) = fvec
         endif
         if ((field%spin .eq. "b") .or. (field%spin .eq. "t")) then
            field%vec(2) = fvec
         endif


C
C     Save field info back into params struct
C
         params%field(ifield) = field
      enddo

      
      end subroutine




c $Id$
