C
C     rt_tddft_os_print.F

      subroutine rt_tddft_os_print_charge (params, it, tt, g_densre_ao)
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
      integer, intent(in)           :: it
      double precision, intent(in)  :: tt
      integer, intent(in)           :: g_densre_ao(2) !alpha beta re part of dens mat ao bas

      
C     == Parameters ==
      character(*), parameter :: pname  = "rt_tddft_os_print_charge: "


C     == Variables ==
      integer ig
      double precision charge, charge_elec, charge_nuc, charge_tot
      character*60 outstr      



C     Total charge on each geom, note that we divide nuclear charge by
C     two since we have two spins--this way charge is 0 for an overall
C     neutral system.
      if (params%print_charge) then
         
         do ig = 1, params%ngeoms
            if (params%geom(ig)%measure) then
               
               charge_nuc = params%geom(ig)%charge_nuc
               
               charge_tot = 0d0 !sum of alpha+bet
               
C     (alpha spin)
               outstr = "Charge [" //
     $              trim (params%geom(ig)%name) // "]"//
     $              " (alpha spin)"
               
               call rt_tddft_calc_charge_geom (params,
     $              g_densre_ao(1), ig, charge_elec)
               
               charge = charge_nuc/2d0 - charge_elec
               charge_tot = charge_tot + charge
               
               call rt_tddft_print_val(params, tt, charge, trim(outstr))
               
               
C     (beta spin)
               outstr = "Charge [" //
     $              trim (params%geom(ig)%name) // "]"//
     $              " (beta spin)"
               
               call rt_tddft_calc_charge_geom (params,
     $              g_densre_ao(2), ig, charge_elec)
               
               charge = charge_nuc/2d0 - charge_elec
               charge_tot = charge_tot + charge
               
               call rt_tddft_print_val(params, tt, charge, trim(outstr))
               
               
C     (alpha+beta total spin)
               outstr = "Charge [" //
     $              trim (params%geom(ig)%name) // "]"//
     $              " (total spin)"
               call rt_tddft_print_val (params, tt,
     $              charge_tot, trim(outstr))
               
            endif
         enddo
      endif
      
      end subroutine
      
c $Id$
