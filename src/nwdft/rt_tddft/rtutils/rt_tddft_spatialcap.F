C
C     rt_tddft_spatialcap.F
C
C     Outputs a nao x nao complex absorbing potential (CAP) build in
C     real space via projection on the AO basis.
C
      subroutine rt_tddft_spatialcap (params, nao, g_zcap)
      implicit none
      
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "util.fh"
#include "cdft.fh"
#include "geomP.fh"
#include "geom.fh"
#include "bas.fh"
#include "rtdb.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==
C      integer, intent(in)            :: rtdb
      type(rt_params_t), intent(in) :: params
      integer, intent(in)            :: nao


C     == Outputs ==
      integer, intent(in)           :: g_zcap


C     == Parameters ==
      character(len=*), parameter :: pname="rt_tddft_spatialcap: "

      
C     == Variables ==
      integer  :: rtdb
      integer  ::  g_pot
      character(255) :: fname
      logical  :: lreload

      rtdb = params%rtdb

      if (.not. ga_create(mt_dbl, nao, nao, "potential",
     $     0, 0, g_pot))
     $     call errquit(pname//"alloc pot failed", 0, 0)


C      
C     Note, the following routine reads in CAP parameters from the rtdb.
C      
      call rt_tddft_spatial_potential (params, nao, g_pot)


C
C     Load into complex GA for shipping out
C
C
      call convert_d2z (1d0, g_pot, 0d0, g_pot, g_zcap)

      
C
C     Scale potential by -i => imaginary absorbing pot
C     
      call ga_scale (g_zcap, zni)


      if (.not. ga_destroy (g_pot))
     $     call errquit(pname//"destroy g_pot", 0, 0)
      
      end subroutine
