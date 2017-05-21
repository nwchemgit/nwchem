C> @file rt_tddft.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C
C  Overall RT-TDDFT documentation:
C
C> @defgroup rt_tddft Real-time TDDFT
C>    @ingroup nwdft
C--------------------------------------------------------------------

C--------------------------------------------------------------------
C> Driver routine for real-time TDDFT.
C>       
C> Call this from main DFT routine after convergence.
C> It will call the closed-shell, open-shell, spin-orbit, etc versions of RT-TDDFT.
C--------------------------------------------------------------------
      logical function rt_tddft (rtdb)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "case.fh"
#include "rtdb.fh"
#include "headers/rt_tddft.fh"

C ^^^ strange that I had to do this


C     == Inputs ==
      integer, intent(in) :: rtdb !< Runtime database

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft: "


C     == Variables ==
      type(rt_params_t) :: params
      character*80 theory
      integer me


      rt_tddft = .false.

      me = ga_nodeid ()
c
c     == Coulomb Attenuation Method (CAM/LC) parameters ==
      if (.not.rtdb_get(rtdb, 'dft:cam_exch', mt_log, 1,
     &   cam_exch)) cam_exch=.false.
      if (.not.rtdb_get(rtdb, 'dft:cam_srhf', mt_log, 1,
     &   cam_srhf)) cam_srhf=.false.
      if (.not.rtdb_get(rtdb, 'dft:cam_omega', mt_dbl, 1,
     &   cam_omega)) cam_omega=0.d0
      if (.not.rtdb_get(rtdb, 'dft:cam_alpha', mt_dbl, 1,
     &   cam_alpha)) cam_alpha=0.d0
      if (.not.rtdb_get(rtdb, 'dft:cam_beta', mt_dbl, 1,
     &   cam_beta)) cam_beta=0.d0
c
C     Store rtdb in params struct, then determine the theory we are using.
C
      params%rtdb = rtdb

      if (.not. rtdb_cget(rtdb, "dft:theory", 1, theory))
     $     call errquit(pname//"theory not specified",0, RTDB_ERR)

      if (theory .eq. "dft") then
         
         if (ipol.eq.1) then             !closed shell DFT
            if (me.eq.0) then
               write (luout, *) " " 
               write (luout, *) "*************************"
               write (luout, *) "* Closed-shell RT-TDDFT *"
               write (luout, *) "*************************"
               write (luout, *) " "
               call util_flush (luout)
            endif
            call rt_tddft_cs_init (params)
            call rt_tddft_header (params)
            call rt_tddft_cs_prop (params)
            
         elseif (ipol.eq.2) then         !open shell DFT
            if (me.eq.0) then
               write (luout, *) " "
               write (luout, *) "***********************"
               write (luout, *) "* Open-shell RT-TDDFT *"
               write (luout, *) "***********************"
               write (luout, *) " "
               call util_flush (luout)
            endif
            call rt_tddft_os_init (params)
            call rt_tddft_header (params)
            call rt_tddft_os_prop (params)
         else
            call errquit (pname//"bad ipol value", 0, 0)
         endif
         
      elseif (theory .eq. "sodft") then  !spin orbit DFT
         call errquit (pname//"spin-orbit not implemented yet",0,0)
      else
         call errquit (pname//"Invalid theory: "//trim(theory), 0, 0)
      endif


      
C
C     Clean up.
C
      call rt_tddft_clean (params)

      

C
C     Finally post-process dens mats to cube files.  Note we explicitly pass subgs option...
C
      if (params%viz_dplot) call rt_tddft_snapshot_dplot (params,
     $     params%viz_subgs)


      rt_tddft = .true.
      end function
c $Id$
