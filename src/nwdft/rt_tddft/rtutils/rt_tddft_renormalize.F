C     
C     Re-normalize the time-dependent density matrix such that the trace
C     sums up to a desired number of electrons.  Mainly used in
C     ionization calculations.
C     
      subroutine rt_tddft_renormalize (params, g_zdens_mo, type)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      type(rt_params_t), intent(in) :: params
      character(len=*), intent(in) :: type       ! = "alpha", "beta", "closedshell"

C     == Outputs ==
      integer, intent(in) :: g_zdens_mo         ! MO density matrix

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_renormalize: "

      
C     == Variables ==
      integer :: me
      double precision :: nelec_alpha, nelec_beta, nelec_total
      double complex :: ztrace, zscale
      double precision :: dtrace


      me = ga_nodeid()

C
C     trace before normalization
C
      call mat_ztrace (g_zdens_mo, ztrace)
      dtrace = dble (ztrace)

      if (type.eq."alpha") then
         if (.not.rtdb_get(params%rtdb,
     $        "rt_tddft:renormalize_nelec_alpha",
     $        mt_dbl, 1, nelec_alpha))
     $        call errquit (pname//"failed to read nelec_alpha",0,0)

         if (nelec_alpha > 1e-6) then
            zscale = nelec_alpha / dtrace 
            call ga_scale (g_zdens_mo, zscale)
            if (me.eq.0) then 
               write(luout,*) pname//"Scaled alpha", 
     $              " density matrix by ", zscale
            endif
         else
            if (me.eq.0) then 
               write(luout,*) pname//"Skipping ", trim(type),
     $              " normalization due to small number of electrons"
            endif
         endif


      elseif (type.eq."beta") then
         if (.not.rtdb_get(params%rtdb,
     $        "rt_tddft:renormalize_nelec_beta",
     $        mt_dbl, 1, nelec_beta))
     $        call errquit (pname//"failed to read nelec_beta",0,0)

         if (nelec_beta > 1e-6) then
            zscale = nelec_beta / dtrace 
            call ga_scale (g_zdens_mo, zscale)
            if (me.eq.0) then 
               write(luout,*) pname//"Scaled beta",
     $              " density matrix by ", zscale
            endif
         else
            if (me.eq.0) then 
               write(luout,*) pname//"Skipping ", trim(type),
     $              " normalization due to small number of electrons"
            endif
         endif

         
      elseif (type.eq."closedshell") then
         if (.not.rtdb_get(params%rtdb,
     $        "rt_tddft:renormalize_nelec_total",
     $        mt_dbl, 1, nelec_total))
     $        call errquit (pname//"failed to read nelec_total",0,0)

         if (nelec_total > 1e-6) then
            zscale = nelec_total / dtrace 
            call ga_scale (g_zdens_mo, zscale)
            if (me.eq.0) then 
               write(luout,*) pname//"Scaled total",
     $              " density matrix by ", zscale
            endif
         else
            if (me.eq.0) then 
               write(luout,*) pname//"Skipping ", trim(type),
     $              " normalization due to small number of electrons"
            endif
         endif

      else
         call errquit("invalid type", 0, 0)
      endif



      
c$$$      if (me.eq.0) then
c$$$         write (luout, *) "trace P' (", trim(type), ") = ", ztrace
c$$$      endif

      end subroutine
