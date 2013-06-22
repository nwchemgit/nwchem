C
C     rt_tddft_output.F
C
C     Contains various routines for quantities to stdout which are
C     universal for open shell, closed shell, etc.
C


C====================================================================
      subroutine rt_tddft_print_val (params, tt, val, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      double precision, intent(in)  :: val
      character(len=*), intent(in)  ::  msg


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_val: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      if (me.eq.0) then 
C         write (luout, "(1x, a, i3, a, 2x, 1f12.4, 1e22.12e3, 5x, a)")
C     $        params%tag, it*100/params%nt, "% ", tt, val, "# "//trim (msg)

c$$$         write (luout, "(a, 2x, i3, a, 2x, 1f12.4, 1es22.12e3, 5x, a)")
c$$$     $        trim(params%tag), (it-1)*100/params%nt, "% ",
c$$$     $        tt, val, "# "//trim (msg)
         
         write (luout, "(a, 2x, 1f11.5, 1es22.12e3, 5x, a)")
     $        trim(params%tag),  tt, val, "# "//trim (msg)


         call util_flush (luout)
      endif
      
      end subroutine



C====================================================================
      subroutine rt_tddft_print_ival (params, tt, ival, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      integer, intent(in)           :: ival
      character(len=*), intent(in)  ::  msg


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_ival: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      if (me.eq.0) then 

c$$$         write (luout, "(a, 2x, i3, a, 2x, 1f12.4, 5x, i17, 5x, a)")
c$$$     $        trim(params%tag), (it-1)*100/params%nt, "% ",
c$$$     $        tt, ival,"# "//trim (msg)

         write (luout, "(a, 2x, 1f11.5, 5x, i17, 5x, a)")
     $        trim(params%tag), tt, ival,"# "//trim (msg)
         
         call util_flush (luout)
      endif
      
      end subroutine

      

C====================================================================
      subroutine rt_tddft_print_quad (params, tt, quad, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      type(rt_quad_t), intent(in)   :: quad
      character(len=*), intent(in)  ::  msg


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_quad: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      if (me.eq.0) then 
         write (luout, "(a, 2x, 1f11.5, 6es22.12e3, 5x, a)")
     $        trim(params%tag), 
     $        tt, quad%xx, quad%xy, quad%xz, quad%yy, quad%yz, quad%zz,
     $        "# "//trim (msg)
         call util_flush (luout)
      endif
      
      end subroutine



C====================================================================
      subroutine rt_tddft_print_vector (params, tt, vec, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      type(rt_vector_t), intent(in) :: vec
      character(len=*), intent(in)  ::  msg


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_vector: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

C         write (luout, "(a, 2x, 1f11.5, 1es22.12e3, 5x, a)")
C     $        trim(params%tag),  tt, val, "# "//trim (msg)


      if (me.eq.0) then 
         write (luout, "(a, 2x, 1f11.5, 3es22.12e3, 5x, a)")
     $        trim(params%tag), tt,
     $        vec%x, vec%y, vec%z, "# "//trim (msg)
         call util_flush (luout)
      endif

      end subroutine




      
C====================================================================
C
C     Prints tagged message with simulation time.
C
      subroutine rt_tddft_print_msg (params, tt, msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      character(len=*), intent(in)  :: msg


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_msg: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      if (me.eq.0) then 
c$$$         write (luout, "(a, 2x, i3, a, 2x, 1f12.4, 8x, a)")
c$$$     $        trim(params%tag), (it-1)*100/params%nt, "% ", tt,
c$$$     $        "### "//trim(msg)//" ###"

         write (luout, "(a, 2x, 1f11.5, 8x, a)")
     $        trim(params%tag), tt, "### "//trim(msg)//" ###"

         call util_flush (luout)
      endif
      
      end subroutine


C====================================================================
C
C     Prints a tagged notice (without the simulation time but with CPU
C     time)
C
      subroutine rt_tddft_print_timestamp (msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      character(len=*), intent(in) :: msg

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_timestamp: "

C     == Variables ==      
      integer me
      double precision cpusec
      integer isec

      me = ga_nodeid ()
      cpusec = util_cpusec ()

      isec = floor (cpusec)

      if (me.eq.0) then
         write(luout, "(1x,a,i0,a)")
C     $        params%tag//"### "//msg//" (CPU time = ", isec, " s) ###"
     $        "### "//msg//" (CPU time = ", isec, " s) ###"
         call util_flush (luout)
      endif

      
      end subroutine


C====================================================================
C
C     Prints a tagged notice (without the simulation time or CPU time)
C
      subroutine rt_tddft_print_notice (msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      character(len=*), intent(in) :: msg

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_notice: "

C     == Variables ==      
      integer me

      me = ga_nodeid ()

      if (me.eq.0) then
C         write(luout, "(1x,a)") params%tag//"=== "//msg//" ==="
         write(luout, "(1x,a)") "=== "//msg//" ==="
         call util_flush (luout)
      endif

      end subroutine


C====================================================================
C
C     Prints a tagged warning notice (without the simulation time or CPU time)
C
      subroutine rt_tddft_print_warning (msg)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      character(len=*), intent(in) :: msg

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_notice: "
      integer, parameter      :: linewidth = 60

      
C     == Variables ==      
      integer me
      integer i

      
      me = ga_nodeid ()

      if (me.eq.0) then
         write (luout, *) ""
         
         do i = 1, linewidth / 2 
            write (luout, "(a)", advance="no") "*"
         enddo
         
         write (luout, "(a)", advance="no") " WARNING "

         do i = 1, linewidth / 2
            write (luout, "(a)", advance="no") "*"
         enddo
         
         write (luout, *) ""
         
      endif
      
      call rt_tddft_print_justified (msg, linewidth-5)


      if (me.eq.0) then
         write (luout, *) ""         
         do i = 1, linewidth + 9
            write (luout, "(a)", advance="no") "*"
         enddo
         write (luout, *) ""
         write (luout, *) ""
      endif


      end subroutine


C====================================================================
      subroutine rt_tddft_print_justified (msg, linewidth)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "util.fh"
#include "rt_tddft.fh"


C     == Inputs ==
      character(len=*), intent(in) :: msg
      integer, intent(in)          :: linewidth
      

C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_justified: "


C     == Variables ==      
      integer me
      integer i
      integer j
      character outchar

      me = ga_nodeid ()

      if (me.eq.0) then
         j = 0
         write (luout, "(a)", advance="no") "  "
         do i = 1, len_trim(msg)
            j = j + 1
            outchar = msg(i:i)
            write (luout, "(A)", advance="no") msg(i:i)
            if ( (j .gt. linewidth).and.(outchar .eq. " ") ) then
               write (luout, *) ""
               write (luout, "(a)", advance="no") "  "
               j = 0
            endif
         enddo
         
         call util_flush (luout)
      endif
      end subroutine




C====================================================================
      subroutine rt_tddft_print_energies (params, tt, energies)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in)   :: params
C      integer, intent(in)             :: it
      double precision, intent(in)    :: tt
      type(rt_energies_t), intent(in) :: energies


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_energies: "

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      call rt_tddft_print_val (params, tt, params%Enuc, "Enuc")
      call rt_tddft_print_val (params, tt, energies%core, "Ecore")
      call rt_tddft_print_val (params, tt, energies%coul, "Ecoul")
      call rt_tddft_print_val (params, tt, energies%xc(1), "Exc(1)")
      call rt_tddft_print_val (params, tt, energies%xc(2), "Exc(2)")
      call rt_tddft_print_val (params, tt, energies%tot, "Etot")
      call rt_tddft_print_val (params, tt,
     $     energies%tot - energies%gs, "Eadded")
      
      end subroutine



C====================================================================

      subroutine rt_tddft_print_fields (params, tt)
      implicit none
      
#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt

      
C     == Variables ==
      integer igeom
      character*60 outstr
      type(rt_vector_t) eval
      
      call rt_tddft_update_fields (params, tt) !not really needed, but doesnt hurt

C
C     Loop over all geom and print the total applied field on each.
C
      if (params%nexcites .eq. 0) then
         call rt_tddft_print_msg (params, tt, "No applied E-field")
         
      else
         
         do igeom = 1, params%ngeoms
            
            outstr = "Applied E-field [" //
     $           trim (params%geom(igeom)%name) // "] (alpha spin)"
            eval = params%geom(igeom)%vec(1)
            call rt_tddft_print_vector (params, tt,
     $           eval, trim(outstr))

            outstr = "Applied E-field [" //
     $           trim (params%geom(igeom)%name) // "] (beta spin)"
            eval = params%geom(igeom)%vec(2)
            call rt_tddft_print_vector (params, tt,
     $           eval, trim(outstr))

         enddo
      endif
         
      end subroutine




C====================================================================
      subroutine rt_tddft_print_canorg (params, tt)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_print_canorg: "


C     == Variables ==      
      integer me

      me = ga_nodeid()

      if (params%nt < 1) call errquit (pname//"nt must be > 0", 0, 0)

      if (me.eq.0) then 
         write (luout,
     $        "(1x, a, 2x, 1f11.5, 5x, i0, 5x, i0, 5x, a)")
     $        tt, params%ns_ao, params%ns_mo,
     $        "# Canorg: ns_ao, ns_mo"

C     $        params%tag, (it-1)*100/params%nt, "% ",


         call util_flush (luout)
      endif
      
      end subroutine




C====================================================================

      subroutine rt_tddft_print_initialization (params)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(len=*), parameter :: pname =
     $     "rt_tddft_print_initialization: "

      
C     == Variables ==
      integer me
      integer active_geom_indx


      me = ga_nodeid ()

C
C     Print brief summary to stdout.
C
      active_geom_indx = params%geom_active

      if (me.eq.0) then
         write (luout, *)
         write (luout, "(a)")
     $        "=============== Geometry Initialization ==============="

         write (luout, "(2x,a)") 'Active geometry: "'//
     $        trim (params%geom(active_geom_indx)%name)//'"'

         if (params%nbf_ao .eq. params%nbf_mo) then
            write (luout, "(2x,a)") "No linear dependencies"
         elseif (params%nbf_ao .gt. params%nbf_mo) then
            write (luout, "(2x,a,i0,a)") "Found ",
     $           params%nbf_ao - params%nbf_mo, " linear dependencies."
         else
            call errquit (pname//"found nbf_mo > nbf_ao ?!", 0, 0)
         endif
         
         write (luout, "(2x,a,i0)") "Number of atomic orbitals    = ",
     $        params%nbf_ao
         write (luout, "(2x,a,i0)") "Number of molecular orbitals = ",
     $        params%nbf_mo
         
         write (luout, "(a)")
     $        "======================================================="
         write (luout, *)
         call util_flush (luout)
      endif
      
      
      end subroutine





C====================================================================

      subroutine rt_tddft_print_charge(params, tt, g_densre_ao, str)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
      type(rt_params_t), intent(in) :: params
C      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      integer, intent(in)           :: g_densre_ao !re part of dens mat ao bas (tot, a or b)
      character(len=*), intent(in)  :: str !tag to be output

      
C     == Parameters ==
      character(*), parameter :: pname  = "rt_tddft_print_charge: "


C     == Variables ==
      integer ig
      double precision charge_elec, charge_nucl, charge_tot
      character*60 outstr      


      do ig = 1, params%ngeoms
         if (params%geom(ig)%measure) then

            outstr = "Charge [" // trim (params%geom(ig)%name) //
     $           "]" // " " // trim(str)
                           
            charge_nucl = params%geom(ig)%charge_nuc
               
            call rt_tddft_calc_charge_geom (params,
     $           g_densre_ao, ig, charge_elec)
               
C            charge_tot = charge_nucl - charge_elec
CXXX  [KAL]: dont subtract nuclear charge for now
            charge_tot = charge_elec
            
            call rt_tddft_print_val (params, tt,
     $           charge_tot, trim(outstr))
         endif
      enddo

      
      end subroutine
      
      
c $Id$
