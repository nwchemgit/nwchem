c===============================================================c
c                                                               c
c       NCC - a New Coupled-Cluster code for NWChem             c
c                                                               c
c       Developed by:                                           c
c                                                               c
c               Jeff R. Hammond                                 c
c               Leadership Computing Facility                   c
c               Argonne National Laboratory                     c
c               jhammond@mcs.anl.gov                            c
c                                                               c
c               Karol Kowalski                                  c
c               Environmental Molecular Sciences Laboratory     c
c               Pacific Northwest National Laboratory           c
c               karol.kowalski@pnl.gov                          c
c                                                               c
c               Marta Włoch                                     c
c               Department of Chemistry                         c
c               Michigan Technological University               c
c               wloch@mtu.edu                                   c
c                                                               c
c===============================================================c
      subroutine ncc_doubles_create(x2info, ! output array
     1                              label,  ! character name for this array
     2                              nocc,   ! number of occupied orbitals
     3                              nvir)   ! number of virtual orbitals
c
c $Id: ncc_doubles.F,v 1.0 2009/20/10 23:48:58 jhammond Exp $
c
      implicit none
#include "mafdecls.fh"
c#include "global.fh"
#include "errquit.fh"
c#include "stdio.fh"
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
c                       x2info(1) = MA size   ~ s_x2
c                       x2info(2) = MA handle ~ k_x2
c                       x2info(3) = MA handle ~ l_x2
c                       x2info(4) = first dimension ~ nocc
c                       x2info(5) = second dimension ~ nvir
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
      character*(*)      ! label
      integer nocc       ! number of occupied orbitals
      integer nvir       ! number of virtual orbitals
c
c     internal variables
c
      integer k_x2 ! MA handle
      integer l_x2 ! MA handle
      integer s_x2 ! MA size
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_doubles_create'
#endif
c
      s_x2 = nocc*nvir
c
      if (.not.ma_push_get(mt_dbl,s_x2,label,l_x2,k_x2)) then
          call errquit ('ncc_doubles_create: ma_push_get',s_x2,MA_ERR)
      endif
c
      x2info(1) = s_x2
      x2info(2) = k_x2
      x2info(3) = l_x2
      x2info(4) = nocc
      x2info(5) = nvir
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_create'
#endif
c
      return
      end




      subroutine ncc_doubles_destroy(x2info)
c
      implicit none
#include "mafdecls.fh"
c#include "global.fh"
#include "errquit.fh"
c#include "stdio.fh"
c
c     interface variables
c
      integer x2info(*)  ! array containing doubles descriptor:
c
c                       x2info(1) = MA size   ~ s_x2
c                       x2info(2) = MA handle ~ k_x2
c                       x2info(3) = MA handle ~ l_x2
c                       x2info(4) = first dimension ~ nocc
c                       x2info(5) = second dimension ~ nvir
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
c     internal variables
c
      integer l_x2 ! MA handle
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_doubles_destroy'
#endif
c
      l_x2 = x2info(3)

      if (.not.ma_pop_stack(l_x2)) then
          call errquit('ncc_doubles_destroy: ma_pop_stack',l_x2,MA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_doubles_destroy'
#endif
c
      return
      end
