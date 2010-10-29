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
      subroutine ncc_singles_create(x1info, ! output array
     1                              label,  ! character name for this array
     2                              nocc,   ! number of occupied orbitals
     3                              nvir)   ! number of virtual orbitals
c
c $Id$
c
      implicit none
#include "mafdecls.fh"
c#include "global.fh"
#include "errquit.fh"
c#include "stdio.fh"
c
c     interface variables
c
      integer x1info(*)  ! array containing singles descriptor:
c
c                       x1info(1) = MA handle ~ k_x1
c                       x1info(2) = MA handle ~ l_x1
c                       x1info(3) = MA size   ~ s_x1
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
      character*(*) label  ! label
      integer nocc         ! number of occupied orbitals
      integer nvir         ! number of virtual orbitals
c
c     internal variables
c
      integer s_x1 ! MA size
      integer k_x1 ! MA handle
      integer l_x1 ! MA handle
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_singles_create'
#endif
c
      s_x1 = nocc*nvir
c
      if (.not.ma_push_get(mt_dbl,s_x1,label,l_x1,k_x1)) then
          call errquit ('ncc_singles_create: ma_push_get',s_x1,MA_ERR)
      endif
c
      x1info(1) = k_x1
      x1info(2) = l_x1
      x1info(3) = s_x1
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_singles_create'
#endif
c
      return
      end




      subroutine ncc_singles_destroy(x1info)
c
      implicit none
#include "mafdecls.fh"
c#include "global.fh"
#include "errquit.fh"
c#include "stdio.fh"
c
c     interface variables
c
      integer x1info(*)  ! array containing singles descriptor:
c
c                       x1info(1) = MA size   ~ k_x1
c                       x1info(2) = MA handle ~ l_x1
c                       x1info(3) = MA handle ~ s_x1
c
c     other handles, such as for check-pointing, should be 
c     added to this array, so be careful to not hard-core
c     the length too many places
c
c     internal variables
c
      integer l_x1 ! MA handle
c
#ifdef DEBUG_PRINT
      print*,'top of ncc_singles_destroy'
#endif
c
      l_x1 = x1info(2)
c
      if (.not.ma_pop_stack(l_x1)) then
          call errquit('ncc_singles_destroy: ma_pop_stack',l_x1,MA_ERR)
      endif
c
#ifdef DEBUG_PRINT
      print*,'end of ncc_singles_destroy'
#endif
c
      return
      end
