c
c $Id: pspw_qmmm.F,v 1.15 2002-02-28 00:56:26 bylaska Exp $
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
#include "geom.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"


*     **** local variables ****
      logical value
      integer i
      integer h_nfrag,h_indx,h_hndl,ma_type
      double precision q
      character*16 t

      logical     ion_mmexist
      character*2 pspw_qmmm_aname
      integer     control_version
      external    ion_mmexist
      external    pspw_qmmm_aname
      external    control_version


*     *************************
*     **** read in MM data ****
*     *************************
      if (ion_mmexist()) then
    
*        **** model potential is the default, qmmm_type = 1 ****
         if (.not.
     >      rtdb_get(rtdb,'pspw:qmmm_type',mt_int,1,qmmm_type))
     >      qmmm_type = 1

         value = geom_create(geom,'qmmmgeometry')
         value = value.and.geom_rtdb_load(rtdb,geom,'qmmmgeometry')
         value = value.and.geom_ncent(geom,nion)
         if (.not. value) call errquit('opening qmmgeometry',0)
      else
         nion = 0
      end if



*     **** don't do anything if there are no MM atoms ****
      if (nion.gt.0) then

c*        ***** Boundary condition check ****
c         if (control_version().ne.4) then
c          value = geom_destroy(geom)
c          call errquit(
c     >    'Aperiodic boundary conditions must be used with QM/MM',0)
c         end if


*        ***** Allocate pspw_qmmm geometries, velocities, and masses *****
         value = MA_alloc_get(mt_dbl,(3*nion),
     >                        'mm_r2',r2(2),r2(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,(3*nion),
     >                        'mm_r1',r1(2),r1(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,(3*nion),
     >                        'mm_r0',r0(2),r0(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,(nion),
     >                        'mm_mass',mass(2),mass(1))
         if (.not. value) call errquit('out of heap memory',0)


*        **** read in qmmm geometries, velocities, and masses ****
         do i=1,nion
             value = value.and.
     >            geom_cent_get(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
         end do
         call dcopy(3*nion,dbl_mb(r1(1)),1,dbl_mb(r2(1)),1)

         call dcopy(3*nion,0.0d0,0,dbl_mb(r0(1)),1)
         value = value.and.geom_vel_get(geom,dbl_mb(r0(1)))
         value = value.and.geom_masses_get(geom,nion,dbl_mb(mass(1)))
         if (.not. value) call errquit('error reading qmmmgeometry', 0)

*        **** read in fragment list ****
         if (rtdb_ma_get(rtdb, 'pspw:qmmm_frag', ma_type,
     >        h_nfrag, h_hndl)) then


           if (.not.MA_get_index(h_hndl,h_indx))
     >      call errquit(
     >      'pspw_qmmm_init: ma_get_index failed for qmmm_frag',911)
         end if

*        **** find maximum fragment length and nfrag  ****
         call set_max_frag_length_and_nfrag(h_nfrag,int_mb(h_indx),
     >                                      max_frag_length,nfrag)

*        **** Allocate the rest of pspw_qmmm: frag_from_indx(2,nion), ****
*        **** indx_from_frag(max_frag_length,nfrag), size_frag(nfrag) ****
*        **** dsqm_frag(max_frag_length,nfrag), shake_frag(nfrag)     ****
         value = value.and.
     >           MA_alloc_get(mt_int,(2*nion),
     >           'frag_from_indx',frag_from_indx(2),frag_from_indx(1))
         value = value.and.
     >           MA_alloc_get(mt_int,(max_frag_length*nfrag),
     >           'indx_from_frag',indx_from_frag(2),indx_from_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_int,(nfrag),
     >           'size_frag',size_frag(2),size_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_dbl,(max_frag_length*nfrag),
     >           'dsqm_frag',dsqm_frag(2),dsqm_frag(1))
         value = value.and.
     >           MA_alloc_get(mt_log,(nfrag),
     >           'shake_frag',shake_frag(2),shake_frag(1))
         if (.not. value) call errquit('out of heap memory',0)

*        ***** set indx_from_frag,size_frag, and shake_frag ****
         call set_indx_from_frag(h_nfrag,int_mb(h_indx),
     >                           max_frag_length,nfrag,
     >                           int_mb(size_frag(1)),
     >                           int_mb(indx_from_frag(1)),
     >                           int_mb(shake_frag(1)))

*        **** set frag_from_indx ****
         call set_frag_from_indx(max_frag_length,nfrag,
     >                           int_mb(size_frag(1)),
     >                           int_mb(indx_from_frag(1)),
     >                           int_mb(frag_from_indx(1)))

*        **** set dsqm_frag ****
         call set_dsqm_frag(max_frag_length,nfrag,
     >                      int_mb(size_frag(1)),
     >                      int_mb(indx_from_frag(1)),
     >                      dbl_mb(r1(1)),
     >                      dbl_mb(dsqm_frag(1)))


*       **** free h_array ****
        value = MA_free_heap(h_hndl)
        if (.not.value) call errquit('cannot free heap memory',0)

*       **** set katm_list: nkatm, katm, natm, atom ****
        call set_katm_list(2,nion,pspw_qmmm_aname,
     >                     nkatm,katm,natm,atom)


     

*       ***** initialize Model, LJ, Pol, VIB, and CAV ****
        call pspw_LJ_init(rtdb)
        call pspw_Q_init(rtdb)
        if (qmmm_type.eq.1) call Model_Potential_init(rtdb)
        if (qmmm_type.eq.2) call Switching_Charge_init(rtdb)
      end if

      return
      end
*     **************************************
*     *                                    *
*     *    set_max_frag_length_and_nfrag   *
*     *                                    *
*     **************************************
      subroutine set_max_frag_length_and_nfrag(nh,h_array,
     >                                         max_frag_length,nfrag)      
      implicit none
      integer nh
      integer h_array(*)
      integer max_frag_length,nfrag

      integer h,last

      max_frag_length = 0
      nfrag = 0
      last  = 0
      do h=1,nh
         if (h_array(h).lt.0) then
           nfrag = nfrag + 1
           if ((h-last-1).gt.max_frag_length) 
     >            max_frag_length = (h-last-1)
           last = h
         end if
      end do
   
      return
      end 
*     ****************************
*     *                          *
*     *    set_indx_from_frag    *
*     *                          *
*     ****************************
      subroutine set_indx_from_frag(nh,h_array,
     >                              max_frag_length,nfrag,
     >                              size_frag,
     >                              indx_from_frag,
     >                              shake_frag)
      implicit none
      integer nh
      integer h_array(nh)
      integer max_frag_length,nfrag
      integer size_frag(nfrag)
      integer indx_from_frag(max_frag_length,nfrag)
      logical shake_frag(nfrag)
 
      integer h,j,k

      j = 1
      k = 1
      do h=1,nh 
        if (h_array(h).gt.0) then
             indx_from_frag(j,k) = h_array(h)
             j=j+1
           else
              size_frag(k) = j-1
              if (h_array(h).eq. (-1)) then
                shake_frag(k) = .true.
              else
                shake_frag(k) = .false. 
              end if
              k=k+1
              j=1 
           end if
      end do

      return
      end
*     ****************************
*     *                          *
*     *    set_frag_from_indx    *
*     *                          *
*     ****************************
      subroutine set_frag_from_indx(max_frag_length,nfrag,
     >                              size_frag,
     >                              indx_from_frag,
     >                              frag_from_indx)
      implicit none
      integer max_frag_length,nfrag
      integer size_frag(nfrag)
      integer indx_from_frag(max_frag_length,nfrag)
      integer frag_from_indx(2,*)

      integer j,k
      do k=1,nfrag
      do j=1,size_frag(k)
         frag_from_indx(1,indx_from_frag(j,k)) = k
         frag_from_indx(2,indx_from_frag(j,k)) = j
      end do
      end do

      return
      end 

*     ***********************
*     *                     *
*     *    set_dsqm_frag    *
*     *                     *
*     ***********************
      subroutine set_dsqm_frag(max_frag_length,nfrag,
     >                         size_frag,
     >                         indx_from_frag,
     >                         r1,
     >                         dsqm_frag)
      implicit none
      integer max_frag_length,nfrag
      integer size_frag(nfrag)
      integer indx_from_frag(max_frag_length,nfrag)
      real*8  r1(3,*)
      real*8  dsqm_frag(max_frag_length,nfrag)

      integer i,j,k,ii,jj

      do k=1,nfrag
      do j=1,size_frag(k)
         i = j+1
         if (i.gt.size_frag(k)) i = 1
         ii = indx_from_frag(i,k)
         jj = indx_from_frag(j,k)
         dsqm_frag(j,k) = (r1(1,ii)-r1(1,jj))**2
     >                  + (r1(2,ii)-r1(2,jj))**2
     >                  + (r1(3,ii)-r1(3,jj))**2
      end do
      end do
      return
      end


*     **********************************
*     *	                               *
*     *        pspw_qmmm_destroy       *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_destroy()
      implicit none

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** local variables ****
      logical value

      if (nion.gt.0) then
        value = geom_destroy(geom)
        if (.not. value) 
     >     call errquit('error destroying qmmmgeometry', 0)
      end if

      return
      end

*     **********************************
*     *	                               *
*     *        pspw_qmmm_write         *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_write(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** local variables ****
      logical value
      integer i
      double precision rxyz(3),q
      character*16 t

      integer  control_code
      external control_code

*     ***************************
*     **** write out MM data ****
*     ***************************
      if (nion.gt.0) then
        value = .true.
        do i=1,nion
           value = value.and.geom_cent_get(geom,i,t,rxyz,q)
           value = value.and.
     >             geom_cent_set(geom,i,t,dbl_mb(r1(1)+(i-1)*3),q)
        end do
        value = value.and.geom_vel_set(geom,dbl_mb(r0(1)))

        value = value.and.geom_rtdb_delete(rtdb,'qmmmgeometry')
        if ((control_code().eq.1).or.
     >      (control_code().eq.2)) then
        value = value.and.geom_rtdb_store(rtdb,geom,'qmmmgeometry')
        end if
        value = value.and.geom_destroy(geom)
        if (.not. value) call errquit('error writing qmmmgeometry', 0)
      end if

      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_nkatm        *
*     *                                *
*     **********************************
      integer function pspw_qmmm_nkatm()
      implicit none

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"
 
      pspw_qmmm_nkatm = nkatm
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_katm        *
*     *                                *
*     **********************************
      integer function pspw_qmmm_katm(i)
      implicit none
      integer i

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"
 
      pspw_qmmm_katm = int_mb(katm(1)+i-1)
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_natm        *
*     *                                *
*     **********************************
      integer function pspw_qmmm_natm(i)
      implicit none
      integer i

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"
 
      pspw_qmmm_natm = int_mb(natm(1)+i-1)
      return
      end

*     **********************************
*     *	                               *
*     *         pspw_qmmm_atom        *
*     *                                *
*     **********************************
      character*2 function pspw_qmmm_atom(i)
      implicit none
      integer i

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      character*2 aname

      call pspw_copy_atom(1,byte_mb(atom(1)+2*(i-1)),aname)
      pspw_qmmm_atom = aname
      return
      end


*     **********************************
*     *	                               *
*     *         pspw_qmmm_aname        *
*     *                                *
*     **********************************
      character*2 function pspw_qmmm_aname(i)
      implicit none
      integer i

#include "stdio.fh"
#include "geom.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** local variables ****
      character*2  symbol
      character*16 t,name
      real*8 q

      if (.not. geom_cent_tag(geom,i,t))
     >        call errquit(' pspw_qmmm_aname  failed ',i)

      if (.not. geom_tag_to_element(t,symbol,name,q)) then
        symbol = 'bq'
      end if

      pspw_qmmm_aname = symbol
      return
      end


*     ***************************
*     *                         *
*     *      pspw_qmmm_amass    *
*     *                         *
*     ***************************
      real*8 function pspw_qmmm_amass(i)
      implicit none
      integer i

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"


      pspw_qmmm_amass = dbl_mb(mass(1)+i-1)*1822.89d0
      return
      end


*     **********************************
*     *	                               *
*     *         pspw_qmmm_frag         *
*     *	                               *
*     **********************************
      integer function pspw_qmmm_frag(i)
      implicit none
      integer i

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_frag = int_mb(frag_from_indx(1)+2*(i-1))
      return
      end

*     **********************************
*     *	                               *
*     *   pspw_qmmm_different_frag     *
*     *	                               *
*     **********************************
      logical function pspw_qmmm_different_frag(i,j)
      implicit none
      integer i,j

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      integer mi,mj
      mi = int_mb(frag_from_indx(1)+2*(i-1))
      mj = int_mb(frag_from_indx(1)+2*(j-1))
      pspw_qmmm_different_frag = ((mi.eq.0).or.(mj.eq.0).or.(mi.ne.mj))
      return
      end



*     **********************************
*     *	                               *
*     *         pspw_qmmm_found        *
*     *                                *
*     **********************************
      logical function pspw_qmmm_found()
      implicit none

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"
 
      logical value

      value = .false.
      if (nion.gt.0) value = .true.

      pspw_qmmm_found = value
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
****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      logical value

      if (nion.gt.0) then
        value = MA_free_heap(r2(2))
        value = value.and.MA_free_heap(r1(2))
        value = value.and.MA_free_heap(r0(2))
        value = value.and.MA_free_heap(mass(2))
        value = value.and.MA_free_heap(frag_from_indx(2))
        value = value.and.MA_free_heap(indx_from_frag(2))
        value = value.and.MA_free_heap(dsqm_frag(2))
        value = value.and.MA_free_heap(size_frag(2))
        value = value.and.MA_free_heap(shake_frag(2))
        value = value.and.MA_free_heap(atom(2))
        value = value.and.MA_free_heap(natm(2))
        value = value.and.MA_free_heap(katm(2))
        if (.not.value) call errquit('cannot free heap memory',0)

        if (qmmm_type.eq.1) call Model_Potential_end()
        if (qmmm_type.eq.2) call Switching_Charge_end()
        call pspw_LJ_end()
        call pspw_Q_end()
      end if


      return
      end

*     **********************************
*     *	                               *
*     *        pspw_qmmm_nion          *
*     *                                *
*     **********************************
      integer function pspw_qmmm_nion()
      implicit none

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_nion = nion
      return
      end

*     **********************************
*     *	                               *
*     *        pspw_qmmm_ke	           *
*     *                                *
*     **********************************
      real*8 function pspw_qmmm_ke()
      implicit none

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_ke = ekw
      return
      end

*     **********************************
*     *	                               *
*     *       pspw_qmmm_Temperature    *
*     *                                *
*     **********************************
      real*8 function pspw_qmmm_Temperature()
      implicit none

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      real*8 kb
      parameter (kb=3.16679d-6)

      pspw_qmmm_Temperature = 2.0d0*ekw_total/dble(ekw_count)
     >                                      /(3.0d0*nion-6.0d0)
     >                                      /kb
      return
      end


*     **********************************
*     *	                               *
*     *       pspw_qmmm_Print          *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_Print(unit)
      implicit none
      integer unit

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_qmmm_aname
      real*8      pspw_qmmm_amass
      external    pspw_qmmm_aname
      external    pspw_qmmm_amass

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        write(unit,1180)   
        do ii=1,nion
          if (int_mb(frag_from_indx(1)+2*(ii-1)).eq.0) then
            write(unit,1190) ii,pspw_qmmm_aname(ii),
     >                       (dbl_mb(r1(1)+3*(ii-1)+k-1),k=1,3),
     >                       pspw_qmmm_amass(ii)/1822.89d0
          else
            write(unit,1191) ii,pspw_qmmm_aname(ii),
     >                       (dbl_mb(r1(1)+3*(ii-1)+k-1),k=1,3),
     >                       pspw_qmmm_amass(ii)/1822.89d0,
     >                       int_mb(frag_from_indx(1)+2*(ii-1)),
     >                       int_mb(frag_from_indx(1)+2*(ii-1)+1)
          end if 
        end do
      end if

      end if
      return

 1180 FORMAT(/' position of MM ions (au): ')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3)
 1191 FORMAT(5X, I4, A3  ,' (',3F11.5,' ) - atomic mass= ',F6.3,
     >      ' fragment:',2I6)
      end

*     **********************************
*     *	                               *
*     *        pspw_qmmm_Print2        *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_Print2(unit)
      implicit none
      integer unit

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_qmmm_aname
      external    pspw_qmmm_aname

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then
        write(unit,1180)   
        do ii=1,nion
          write(unit,1190) ii,pspw_qmmm_aname(ii),
     >        (dbl_mb(r0(1)+ 3*(ii-1)+ k-1),k=1,3)    
        end do
      end if

      end if
      return

 1180 FORMAT(/' velocity of pseudopotential water (au):')
 1190 FORMAT(5X, I4, A3  ,' (',3F11.5,' )')
      end




*     **********************************
*     *	                               *
*     *       pspw_qmmm_PrintXYZ       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_PrintXYZ(unit)
      implicit none
      integer unit

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     ***** local variables ****
      integer taskid,MASTER
      parameter (MASTER=0)
      integer ii,k

      character*2 pspw_qmmm_aname
      external    pspw_qmmm_aname

      if (nion.gt.0) then
      call Parallel_taskid(taskid)

      if (taskid.eq.MASTER) then  
       do ii=1,nion
        write(unit,*) pspw_qmmm_aname(ii),'      ',
     >               (dbl_mb(r1(1)+3*(ii-1)+k-1)*0.529177d0,k=1,3)
       end do
      end if

      end if

      return
      end


*     **********************************
*     *	                               *
*     *        pspw_qmmm_rion          *
*     *                                *
*     **********************************
      real*8 function pspw_qmmm_rion(i,ii)
      implicit none
      integer i,ii

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_rion = dbl_mb(r1(1)+3*(ii-1)+i-1)
      return
      end

*     **********************************
*     *	                               *
*     *        pspw_qmmm_vion          *
*     *                                *
*     **********************************
      real*8 function pspw_qmmm_vion(i,ii)
      implicit none
      integer i,ii

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_vion = dbl_mb(r0(1)+3*(ii-1)+i-1)
      return
      end


*     **********************************
*     *	                               *
*     *        pspw_qmmm_shift         *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_shift()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      call dcopy((3*nion),dbl_mb(r1(1)),1,dbl_mb(r0(1)),1)
      call dcopy((3*nion),dbl_mb(r2(1)),1,dbl_mb(r1(1)),1)

      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Generate_V       *
*     *                                *
*     **********************************

      subroutine pspw_qmmm_Generate_V(n2ft3d,rgrid,Vqm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 Vqm(*)

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      if (nion.gt.0) then

*         **** Model Charge ****
c          if (qmmm_type.eq.0)
c     >      call Model_Charge_Generate_V(nion,dbl_mb(r1(1)),
c     >                                   n2ft3d,rgrid,Vqm)           

*         **** Model Potential ****
          if (qmmm_type.eq.1)
     >      call Model_Potential_Generate_V(n2ft3d,rgrid,Vqm)

*         **** Switching Charge Potential ****
          if (qmmm_type.eq.2)
     >      call Switching_Charge_Generate_V(n2ft3d,rgrid,Vqm)
 
      end if

      return
      end

*     **********************************
*     *                                *
*     *        pspw_qmmm_Update        *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_Update(algorithm,
     >                           n2ft3d,rgrid,rho,dv,
     >                           dt,fion,deltar)
      implicit none
      integer algorithm
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 dt
      real*8 fion(3,*)
      real*8 deltar

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** local variables ****
      logical value
      integer fmm(2)


*     **** allocate fmm force ***
      value = MA_push_get(mt_dbl,(3*nion),
     >                     'fmm',fmm(2),fmm(1))
      if (.not.value)
     > call errquit('pspw_qmmm_Update:error push stack',0)

*     **** Get Forces ****
      call pspw_qmmm_Fion(fion)
      call pspw_qmmm_Fmm(n2ft3d,rgrid,rho,dv,dbl_mb(fmm(1)))
      deltar = 0.0d0


*     **** steepest descent update of water ****
      if (algorithm.eq.0) then
      call pspw_qmmm_SD_subUpdate(nion,dbl_mb(r2(1)),
     >                                  dbl_mb(r1(1)),
     >                                  dbl_mb(fmm(1)),
     >                                  dbl_mb(mass(1)),
     >                                  dt,
     >                                  deltar)

*     **** Newton step update of water ****
      else if (algorithm.eq.1) then
        call pspw_qmmm_Newton_subUpdate(nion,dbl_mb(r2(1)),
     >                                  dbl_mb(r1(1)),
     >                                  dbl_mb(r0(1)),
     >                                  dbl_mb(fmm(1)),
     >                                  dbl_mb(mass(1)),
     >                                  dt,ekw)
        ekw_total = ekw_total + ekw
        ekw_count = ekw_count + 1

*     **** Verlet step update of water ****
      else if (algorithm.eq.2) then
        call pspw_qmmm_Verlet_subUpdate(nion,dbl_mb(r2(1)),
     >                                  dbl_mb(r1(1)),
     >                                  dbl_mb(r0(1)),
     >                                  dbl_mb(fmm(1)),
     >                                  dbl_mb(mass(1)),
     >                                  dt,ekw)
        ekw_total = ekw_total + ekw
        ekw_count = ekw_count + 1
      end if


*     **** deallocate fmm force ***
      value = MA_pop_stack(fmm(2))
      if (.not.value)
     > call errquit('pspw_qmmm_Update:error pop stack',0)
      return
      end


*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_ionMM     *
*     *                                *
*     **********************************

*  This function returns the electronstatic
*  interaction energy between the QM subsystem's
* ions and the MM subsystem.

      real*8 function pspw_qmmm_Energy_ionMM()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"
*     **** external functions ****
      real*8   pspw_Q_Energy_IonMM
      external pspw_Q_Energy_IonMM


      pspw_qmmm_Energy_ionMM = pspw_Q_Energy_IonMM()
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_QMMM      *
*     *                                *
*     **********************************

*  This function returns the Lenard-Jones (LJ)
* dispersion and repulsion energy between the
* QM subsystem and the MM subsystem.

      real*8 function pspw_qmmm_Energy_QMMM()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** external functions ****
      real*8   pspw_LJ_Energy_IonMM
      external pspw_LJ_Energy_IonMM

      pspw_qmmm_Energy_QMMM = pspw_LJ_Energy_IonMM()
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_MMMM      *
*     *                                *
*     **********************************

*  This function returns the electrostatic
* and Lenard-Jones (LJ) dispersion and repulsion 
* energy between the QM subsystem and the MM subsystem.

      real*8 function pspw_qmmm_Energy_MMMM()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

*     **** external functions ****
      real*8   pspw_LJ_Energy_MM,pspw_Q_Energy_MM
      external pspw_LJ_Energy_MM,pspw_Q_Energy_MM

      pspw_qmmm_Energy_MMMM = pspw_LJ_Energy_MM()
     >                      + pspw_Q_Energy_MM()
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_pol       *
*     *                                *
*     **********************************

*  This function returns the polarization energy
* of the MM subsystem.

      real*8 function pspw_qmmm_Energy_pol()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_Energy_pol = 0.0d0
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_vib       *
*     *                                *
*     **********************************

*  This function returns the vibration energy
* of the MM subsystem.

      real*8 function pspw_qmmm_Energy_vib()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_Energy_vib = 0.0d0
      return
      end

*     **********************************
*     *	                               *
*     *     pspw_qmmm_Energy_cav       *
*     *                                *
*     **********************************

*  This function returns the polarization energy
* of the MM subsystem.

      real*8 function pspw_qmmm_Energy_cav()
      implicit none

#include "mafdecls.fh"

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      pspw_qmmm_Energy_cav = 0.0d0
      return
      end


*     **********************************
*     *                                *
*     *        pspw_qmmm_Fion          *
*     *                                *
*     **********************************

*     This routine calculates the force on the QM ions.
*
*     Exit - fion: force on the QM ions
*
*     Uses - ion_nion
*
*     Author - Eric Bylaska

      subroutine pspw_qmmm_Fion(fion)
      implicit none
      real*8 fion(3,*)

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      if (nion.gt.0) then

*       **** get LJ forces ****
        call pspw_LJ_Fion(fion)

*       **** get Q forces ****
        call pspw_Q_Fion(fion)
      end if

      return
      end
     
*     **********************************
*     *                                *
*     *        pspw_qmmm_Fmm          *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_Fmm(n2ft3d,rgrid,rho,dv,f_mm)
      implicit none
      integer n2ft3d
      real*8 rgrid(3,*)
      real*8 rho(*)
      real*8 dv
      real*8 f_mm(3,*)

****** pspw_qmmm common block ****
#include "pspw_qmmm.fh"

      if (nion.gt.0) then

*        **** get forces from model charge ****

*        **** get forces from model potential ****
         if (qmmm_type.eq.1)
     >      call Model_Potential_Fmm(n2ft3d,rgrid,rho,dv,f_mm)
 
         if (qmmm_type.eq.2)
     >    call Switching_Charge_Potential_Fmm(n2ft3d,rgrid,rho,dv,f_mm)

*        **** get LJ forces ****
         call pspw_LJ_Fmm(f_mm)

*        **** get Q forces ****
         call pspw_Q_Fmm(f_mm)

      end if
      return
      end

*     **********************************
*     *                                *
*     *      pspw_qmmm_SD_subUpdate    *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_SD_subUpdate(nion,r2,r1,fmm,mass,
     >                                  dt,deltar)
      implicit none
      integer nion
      real*8 r2(3,*)
      real*8 r1(3,*)
      real*8 fmm(3,*)
      real*8 mass(*)
      real*8 dt
      real*8 deltar

*     **** local variables ***
      integer ii
      real*8 dti,sum


*     **** steepest descent update of MM subsystem ****
      do ii=1,nion
        dti = dt/dsqrt(mass(ii))
        r2(1,ii) = r1(1,ii) + dti*fmm(1,ii)
        r2(2,ii) = r1(2,ii) + dti*fmm(2,ii)
        r2(3,ii) = r1(3,ii) + dti*fmm(3,ii)

      end do

*    **** find maximum force ****
      deltar = 0.0d0
      do ii=1,nion
        sum = dsqrt(fmm(1,ii)**2
     >             +fmm(2,ii)**2
     >             +fmm(3,ii)**2)
        if (sum.gt.deltar) deltar = sum
      end do

      return
      end

*     **********************************
*     *                                *
*     *  pspw_qmmm_Newton_subUpdate    *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_Newton_subUpdate(nion,
     >                                      r2,r1,v,fmm,mass,
     >                                      dt,ekw)
      implicit none
      integer nion
      real*8 r2(3,*)
      real*8 r1(3,*)
      real*8 v(3,*)
      real*8 fmm(3,*)
      real*8 mass(*)
      real*8 dt,ekw

*     **** local variables ****
      integer ii
      real*8  dti


*     **** steepest descent update of water ****
      do ii=1,nion
        dti = 0.5d0*dt*dt/(mass(ii))
        r2(1,ii) = r1(1,ii) + dt *v(1,ii) + dti*fmm(1,ii)
        r2(2,ii) = r1(2,ii) + dt *v(2,ii) + dti*fmm(2,ii)
        r2(3,ii) = r1(3,ii) + dt *v(3,ii) + dti*fmm(3,ii)
      end do

*     ***** find kinetic energy ****
      ekw = 0.0d0
      do ii=1,nion
         ekw = ekw + mass(ii)*(v(1,ii)**2+v(2,ii)**2+v(3,ii)**2)
      end do
      ekw = 0.5d0*ekw

      return
      end



*     **********************************
*     *                                *
*     *  pspw_qmmm_Verlet_subUpdate    *
*     *                                *
*     **********************************
      subroutine pspw_qmmm_Verlet_subUpdate(nion,
     >                                 r2,r1,r0,fmm,mass,
     >                                 dt,ekw)
      implicit none
      integer nion
      real*8 r2(3,*)
      real*8 r1(3,*)
      real*8 r0(3,*)
      real*8 fmm(3,*)
      real*8 mass(*)
      real*8 dt,ekw

*     **** local variables ****
      integer i,ii
      real*8 h,dti

*     **** verlet update of water ****
      do ii=1,nion
        dti = dt*dt/(mass(ii))
        r2(1,ii) = 2*r1(1,ii) - r0(1,ii) + dti*fmm(1,ii)
        r2(2,ii) = 2*r1(2,ii) - r0(2,ii) + dti*fmm(2,ii)
        r2(3,ii) = 2*r1(3,ii) - r0(3,ii) + dti*fmm(3,ii)
      end do

*     **** make rwater0 the velocity - note that velocity is deleted after ****
*     **** pspw_qmmm_shift call                                            ****
      h = 1.0d0/(2.0d0*dt)
      do ii=1,nion
      do i=1,3
         r0(i,ii) = h*(r2(i,ii) - r0(i,ii))
      end do
      end do

*     ***** find kinetic energy ****
      ekw = 0.0d0
      do ii=1,nion
         ekw = ekw + mass(ii)*(r0(1,ii)**2+r0(2,ii)**2+r0(3,ii)**2)
      end do
      ekw = 0.5d0*ekw

      return
      end



