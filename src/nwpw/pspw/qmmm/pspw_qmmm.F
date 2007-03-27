c
c $Id: pspw_qmmm.F,v 1.23 2007-03-27 02:02:50 bylaska Exp $
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
      integer taskid,MASTER
      parameter (MASTER=0)
      logical value
      integer i,j,ii,jj,ia,ni,nj,frag_size,shake_ptr
      real*8  volume
      character*80 rtdb_name

*     **** external functions ****
      integer     ion_nion,ion_nion_qm,ion_nion_mm,control_version
      integer     ewald_nshl3d,ewald_rcell_ptr,ion_rion_ptr
      character*7 c_index_name
      real*8      lattice_unita
      external    ion_nion,ion_nion_qm,ion_nion_mm,control_version
      external    ewald_nshl3d,ewald_rcell_ptr,ion_rion_ptr
      external    c_index_name
      external    lattice_unita


      qmmm_found  =  (ion_nion_mm().gt.0)
      shake_found = .false.

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
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'nb_shake',
     >                        nb_shake(2),nb_shake(1))
         value = value.and.
     >           MA_alloc_get(mt_int,nkfrag,'indx_shake_start',
     >                        indx_shake_start(2),indx_shake_start(1))
         if (.not.value) call errquit('cannot free heap memory',0,
     >       MA_ERR)



*        **** set index_start and size_frag****
         ii = 0
         jj = 0
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


            rtdb_name = 'pspw_qmmm_frag_switch_Rin:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rin(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rin(1)+ia-1) =(2.0160d0/0.529177d0)

            rtdb_name = 'pspw_qmmm_frag_switch_Rout:'//c_index_name(ia)
            value = rtdb_get(rtdb,rtdb_name,mt_dbl,1,
     >                       dbl_mb(switch_Rout(1)+ia-1))
            if (.not.value) 
     >         dbl_mb(switch_Rout(1)+ia-1) =(3.1287d0/0.529177d0)

*           **** shake stuff ****
            rtdb_name = 'pspw_qmmm_frag_n_sdist:'//c_index_name(ia)
            if (.not.rtdb_get(rtdb,rtdb_name,mt_int,1,nj)) nj=0
            int_mb(indx_shake_start(1)+ia-1) = jj+1
            int_mb(nb_shake(1)+ia-1)         = nj
            jj = jj + nj

         end do
*        *** more shake ****
         shake_found = (jj.gt.0)
         if (shake_found) then
            value = MA_alloc_get(mt_dbl,jj,'dsq_shake',
     >                        dsq_shake(2),dsq_shake(1))
            if (.not.value) call errquit(
     >       'cannot allocate heap memory',0,MA_ERR)

            do ia=1,nkfrag
               rtdb_name = 'pspw_qmmm_frag_sdist:'//c_index_name(ia)
               value = rtdb_get(rtdb,rtdb_name,mt_dbl,
     >          int_mb(nb_shake(1)+ia-1),
     >          dbl_mb(dsq_shake(1)+int_mb(indx_shake_start(1)+ia-1)-1))
               if (.not. value) call errquit(
     >          'pspw_qmmm_init:failed reading frag_sdist',0,RTDB_ERR)
            end do
         end if


*       ***** initialize Model, LJ, Pol, VIB, and CAV ****
        call pspw_LJ_init(rtdb)

*       ***** make sure mm fragments are in cell ****
        call pspw_qmmm_incell1(dbl_mb(ion_rion_ptr()))


*       ***** write out header info ****
        call Parallel_taskid(taskid)
        if (taskid.eq.MASTER) then
           do ia=1,nkfrag
             write(*,*) "fragment ",ia
             write(*,*) "  - Switching Parameters = ",
     >                  dbl_mb(switch_Rin(1)+ia-1),
     >                  dbl_mb(switch_Rout(1)+ia-1)

             nj = int_mb(nb_shake(1)+ia-1)
             if (nj.gt.0) then
             shake_ptr = dsq_shake(1)+int_mb(indx_shake_start(1)+ia-1)-1
             write(*,*) "  - shake                = ",
     >                  (dsqrt(dbl_mb(shake_ptr+j-1)),j=1,nj)
             end if
           end do
        end if


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
        value = value.and.MA_free_heap(nb_shake(2))
        value = value.and.MA_free_heap(indx_shake_start(2))
        if (shake_found) value = value.and.MA_free_heap(dsq_shake(2))
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
*     *                                *
*     *    pspw_qmmm_shake_found       *
*     *                                *
*     **********************************
      logical function pspw_qmmm_shake_found()
      implicit none

#include "pspw_qmmm.fh"

      pspw_qmmm_shake_found = shake_found
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


*     ******************************
*     *                            *
*     *       pspw_qmmm_frag_cm    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_frag_cm(size_frag,kfrag_s,rion,rcm,mtotal)
      implicit none
      integer size_frag,kfrag_s
      real*8 rion(3,*),rcm(3),mtotal

#include "mafdecls.fh"

*     **** local variables ****
      integer k,kk,ptr
      real*8  m

*     **** external functions ****
      integer  ion_amass_ptr
      external ion_amass_ptr

      ptr = ion_amass_ptr()

      mtotal = 0.0d0
      rcm(1) = 0.0d0
      rcm(2) = 0.0d0
      rcm(3) = 0.0d0
      kk = kfrag_s
      do k=1,size_frag
         m = dbl_mb(ptr+kk-1)
         rcm(1) = rcm(1) + m*rion(1,kk)
         rcm(2) = rcm(2) + m*rion(2,kk)
         rcm(3) = rcm(3) + m*rion(3,kk)
         mtotal = mtotal + m
         kk = kk + 1
      end do
      rcm(1) = rcm(1)/mtotal
      rcm(2) = rcm(2)/mtotal
      rcm(3) = rcm(3)/mtotal

      return
      end



*     ******************************
*     *                            *
*     *       pspw_qmmm_incell1    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell1(r1)
      implicit none
      real*8 r1(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1
      real*8  rcm(3),mtotal

      do w1=1,nfrag
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell1_frag(rcm,n1,r1(1,ks1))
      end do
      return
      end 

*     ******************************
*     *                            *
*     *       pspw_qmmm_incell2    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell2(r1,r2)
      implicit none
      real*8 r1(3,*),r2(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1
      real*8  rcm(3),mtotal

      do w1=1,nfrag
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell2_frag(rcm,n1,r1(1,ks1),r2(1,ks1))
      end do
      return
      end


*     ******************************
*     *                            *
*     *       pspw_qmmm_incell3    *
*     *                            *
*     ******************************
*   Computes the center of mass and total mass of a fragment.
*
*    Entry - size_frag - fragment size
*            kfrag_s   - index to start of fragment
*            rion      - total atom list
*    Exit - rcm    - center of mass of fragment
*           mtotal - total mass of fragment
*
      subroutine pspw_qmmm_incell3(r1,r2,r3)
      implicit none
      real*8 r1(3,*),r2(3,*),r3(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer w1,n1,ks1
      real*8  rcm(3),mtotal

      do w1=1,nfrag
        n1  = int_mb(size_frag(1)+w1-1)
        ks1 = int_mb(indx_frag_start(1)+w1-1)
        call pspw_qmmm_frag_cm(n1,ks1,r1,rcm,mtotal)
        call lattice_incell3_frag(rcm,n1,r1(1,ks1),r2(1,ks1),r3(1,ks1))
      end do
      return
      end



*     ******************************
*     *                            *
*     *       pspw_qmmm_shake      *
*     *                            *
*     ******************************
*   shakes the fragments 
*
      subroutine pspw_qmmm_shake(r2,r1)
      implicit none
      real*8 r2(3,*),r1(3,*)

#include "mafdecls.fh"
#include "errquit.fh"
#include "pspw_qmmm.fh"

*     **** local variables ****
      integer maxit
      double precision tol
      parameter (maxit=1000,tol=1.0d-4)

      integer n1,nb,w1,wk,wks,ks1

*     **** external functions ****
      integer  ion_amass_ptr
      external ion_amass_ptr

      do w1=1,nfrag
         wk = int_mb(kfrag(1)+w1-1)
         nb = int_mb(nb_shake(1)+wk-1)
         if (nb.gt.0) then
            n1  = int_mb(size_frag(1)+w1-1)
            ks1 = int_mb(indx_frag_start(1)+w1-1)
            wks = int_mb(indx_shake_start(1)+wk-1)
            call shake_chain2(n1,nb,
     >                        tol,maxit,
     >                        dbl_mb(dsq_shake(1)+wks-1),
     >                        dbl_mb(ion_amass_ptr()+ks1-1),
     >                        r2(1,ks1),r1(1,ks1))
         end if
      end do

      return
      end
