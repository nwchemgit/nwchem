c
c $Id: pspw_qmmm.F,v 1.22 2007-03-22 20:46:22 bylaska Exp $
c

***********************************************************
*                                                         *
*                    PSPW-QMMM  module                    *
*                                                         *
* Top-level interface to Model_potential and Model_charge *
* modules                                                 *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*    -- developed by Eric J. Bylaska on October 18,2001   *
*                                                         *
***********************************************************
*
*
*  

*     **********************************
*     *	                               *
*     *          pspw_qmmm_init        *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"


*     **** local variables ****
      logical value
      integer i,ii,ia,ni,frag_size
      character*80 rtdb_name

*     **** external functions ****
      integer     ion_nion,ion_nion_qm,ion_nion_mm,control_version
      integer     ewald_nshl3d,ewald_rcell_ptr
      character*7 c_index_name
      external    ion_nion,ion_nion_qm,ion_nion_mm,control_version
      external    ewald_nshl3d,ewald_rcell_ptr
      external    c_index_name


      qmmm_found =  (ion_nion_mm().gt.0)

*     **** don't do anything if there are no MM atoms ****
      if (qmmm_found) then

*        **** set the nshl3d and rcell ****
         if (control_version().eq.3) then
            nshl3d   = ewald_nshl3d()
            rcell(1) = ewald_rcell_ptr()
         else
            nshl3d = 1
            if (.not.MA_alloc_get(mt_dbl,3,'rcell',rcell(2),rcell(1))) 
     >      call errquit(
     >      'pspw_qmmm_init:failed allocating rcell from heap',0,MA_ERR)
            dbl_mb(rcell(1))   = 0.0d0
            dbl_mb(rcell(1)+1) = 0.0d0
            dbl_mb(rcell(1)+2) = 0.0d0
         end if

*        **** reading number of frags ****
         rtdb_name = 'pspw_qmmm_nfrag'
         value = rtdb_get(rtdb,rtdb_name,mt_int,1,nfrag)
         if (.not. value)
     >   call errquit(
     >   'pspw_qmmm_init:failed reading nfrag',0,RTDB_ERR)

*        **** reading number of kinds of frags ****
         rtdb_name = 'pspw_qmmm_nkfrag'
         value = rtdb_get(rtdb,rtdb_name,mt_int,1,nkfrag)
         if (.not. value)
     >   call errquit(
     >   'pspw_qmmm_init:failed reading nkfrag',0,RTDB_ERR)


         value = MA_alloc_get(mt_int,nfrag,'index_frag_start',
     >                        indx_frag_start(2),indx_frag_start(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nfrag,'size_frag',
     >                        size_frag(2),size_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nfrag,'kfrag',kfrag(2),kfrag(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,nkfrag,'switch_Rin',
     >                        switch_Rin(2),switch_Rin(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,nkfrag,'switch_Rout',
     >                        switch_Rout(2),switch_Rout(1))
         if (.not.value) call errquit('cannot free heap memory',0,
     >       MA_ERR)



*        **** set index_start and size_frag****
         ii = 0
         do ia=1,nkfrag

            rtdb_name = 'pspw_qmmm_frag_size:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,frag_size)
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading frag_size',0,RTDB_ERR)

            rtdb_name = 'pspw_qmmm_frag_nindex_start:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,ni)
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading nindx_start',0,RTDB_ERR)
            
            do i=1,ni
               int_mb(size_frag(1)+ii+i-1) = frag_size
               int_mb(kfrag(1)+ii+i-1)     = ia
            end do

            rtdb_name = 'pspw_qmmm_frag_index_start:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,
     >                       ni,int_mb(indx_frag_start(1)+ii))
            if (.not. value)
     >        call errquit(
     >       'pspw_qmmm_init:failed reading indx_start',0,RTDB_ERR)
            ii = ii + ni


            rtdb_name = 'pspw_qmmm_frag_swtich_Rin:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                       dbl_mb(switch_Rin(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rin(1)+ia-1) =(2.0160d0/0.529177d0)

            rtdb_name = 'pspw_qmmm_frag_swtich_Rout:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_int,1,
     >                       dbl_mb(switch_Rout(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rout(1)+ia-1) =(3.1287d0/0.529177d0)

         end do


*       ***** initialize Model, LJ, Pol, VIB, and CAV ****
        call pspw_LJ_init(rtdb)
      end if

      return
      end


*     **********************************
*     *	                               *
*     *         pspw_qmmm_end 	       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

      logical value
      integer  control_version
      external control_version

      if (qmmm_found) then
        value =           MA_free_heap(indx_frag_start(2))
        value = value.and.MA_free_heap(size_frag(2))
        value = value.and.MA_free_heap(kfrag(2))
        value = value.and.MA_free_heap(switch_Rin(2))
        value = value.and.MA_free_heap(switch_Rout(2))
        if (control_version().ne.3) then
        value = value.and.MA_free_heap(rcell(2))
        end if
        if (.not.value) call errquit('cannot free heap memory',0,
     &       MA_ERR)

        call pspw_LJ_end()
      end if

      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_found        *
*     *                                *
*     **********************************
      logical function pspw_qmmm_found()
      implicit none

#include "pspw_qmmm.fh"

      pspw_qmmm_found = qmmm_found
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_LJ_E 	       *
*     *                                *
*     **********************************

      real*8 function pspw_qmmm_LJ_E()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      real*8   pspw_LJ_E
      external pspw_LJ_E

      pspw_qmmm_LJ_E = pspw_LJ_E(ion_nion(),
     >                           ion_nion_qm(),
     >                           int_mb(ion_katm_ptr()),
     >                           nfrag,
     >                           int_mb(indx_frag_start(1)),
     >                           int_mb(size_frag(1)),
     >                           nshl3d,dbl_mb(rcell(1)),
     >                           dbl_mb(ion_rion_ptr()))
      return
      end



*     **********************************
*     *                                *
*     *         pspw_qmmm_Q_E          *
*     *                                *
*     **********************************

      real*8 function pspw_qmmm_Q_E()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  psp_zv_ptr,ion_amass_ptr
      external psp_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      pspw_qmmm_Q_E = pspw_Q_E(ion_nion(),
     >                         ion_nion_qm(),
     >                         int_mb(ion_katm_ptr()),
     >                         dbl_mb(psp_zv_ptr()),
     >                         dbl_mb(ion_amass_ptr()),
     >                         nfrag,
     >                         int_mb(indx_frag_start(1)),
     >                         int_mb(size_frag(1)),
     >                         int_mb(kfrag(1)),
     >                         dbl_mb(switch_Rin(1)),
     >                         dbl_mb(switch_Rout(1)),
     >                         nshl3d,dbl_mb(rcell(1)),
     >                         dbl_mb(ion_rion_ptr()))
      return
      end





*     **********************************
*     *                                *
*     *         pspw_qmmm_LJ_fion      *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_LJ_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr

      call  pspw_LJ_fion(ion_nion(),
     >                   ion_nion_qm(),
     >                   int_mb(ion_katm_ptr()),
     >                   nfrag,
     >                   int_mb(indx_frag_start(1)),
     >                   int_mb(size_frag(1)),
     >                   nshl3d,dbl_mb(rcell(1)),
     >                   dbl_mb(ion_rion_ptr()),fion)
      return
      end



*     **********************************
*     *                                *
*     *         pspw_qmmm_Q_fion       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_Q_fion(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

c     **** external functions ***
      integer  ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      external ion_nion,ion_nion_qm,ion_katm_ptr,ion_rion_ptr
      integer  psp_zv_ptr,ion_amass_ptr
      external psp_zv_ptr,ion_amass_ptr
      real*8   pspw_Q_E
      external pspw_Q_E

      call pspw_Q_fion(ion_nion(),
     >                 ion_nion_qm(),
     >                 int_mb(ion_katm_ptr()),
     >                 dbl_mb(psp_zv_ptr()),
     >                 dbl_mb(ion_amass_ptr()),
     >                 nfrag,
     >                 int_mb(indx_frag_start(1)),
     >                 int_mb(size_frag(1)),
     >                 int_mb(kfrag(1)),
     >                 dbl_mb(switch_Rin(1)),
     >                 dbl_mb(switch_Rout(1)),
     >                 nshl3d,dbl_mb(rcell(1)),
     >                 dbl_mb(ion_rion_ptr()),fion)
      return
      end




*     **********************************
*     *                                *
*     *         pspw_qmmm_fion         *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_fion(fion)
      implicit none
      real*8 fion(3,*)

      call pspw_qmmm_LJ_fion(fion)
      call pspw_qmmm_Q_fion(fion)

      return
      end


