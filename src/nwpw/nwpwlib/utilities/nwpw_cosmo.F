
*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_init              *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_init(rtdb0,geom)
      implicit none
      integer rtdb0,geom

#include "bafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      logical  omaster,oprint,value,restart
      integer  taskid,np,i,j,pcount,ia,ne(2),ispin,cosmo_sccor
      integer  MASTER
      parameter(MASTER=0)

*     ***** external functions *****
      logical  control_print
      external control_print
      integer  Parallel_threadid
      external Parallel_threadid
      integer  ion_natm,ion_nkatm0
      external ion_natm,ion_nkatm0
      real*8   psp_zv
      external psp_zv

      rtdb = rtdb0
c
c     == Cosmo initialization ==
      call Parallel_taskid(taskid)
      omaster =  (taskid.eq.MASTER).and.(Parallel_threadid().eq.MASTER)
      oprint = omaster.and.control_print(print_medium)

      do_cosmo_smd = .false.
      nefc = 0
      gcds = 0.0d0
      dielec = 0.0d0
      screen = 0.0d0
      cosmo1_on = .false.
      cosmo2_on = .false.

      if (.not.btdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) 
     >   cosmo_on = .false.

      if(cosmo_on) then

         if(.not.btdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc)) nefc = 0
         if (nefc.eq.0) then
            call cosmo_initialize(rtdb,geom,-1,oprint)
            firstrestart = .false.
            firstrun = .true.
         else
            firstrestart = .true.
            firstrun = .false.
         end if

c        **** set the cosmo type ****
         if (.not.btdb_get(rtdb,'nwpw:cosmo2',mt_log,1,cosmo2_on))
     >      cosmo2_on = .false.
         cosmo1_on = .not.cosmo2_on


c
c        Turn cosmo on, we want to run the calculation
c        Start with gas_phase run unless told otherwise
c


         cosmo_on = .true.
         if (.not.btdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     >      cosmo_phase = 1

         if (.not.btdb_get(rtdb,'cosmo:do_cosmo_smd',
     >                     mt_log,1,do_cosmo_smd))
     >      do_cosmo_smd = .false.
     
         gcds = 0.0d0
         if (do_cosmo_smd) then
            if (.not.btdb_get(rtdb,'cosmo:gcds',mt_dbl,1,gcds)) gcds = 0.0d0
         end if
         if (.not.btdb_get(rtdb,'cosmo:dielec',mt_dbl,1,dielec)) 
     >      dielec = 0.0d0
         if (.not.btdb_get(rtdb,'cosmo:sccor',mt_int,1,cosmo_sccor)) 
     >      cosmo_sccor = 11

         screen = (dielec - 1.0d0)/(dielec+0.5d0)

         if(.not.btdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc)) nefc = 0
         if (oprint) then
            !write(*,*) "nefc,restart=",nefc,firstrestart
            if (firstrestart) then
               write(luout,*) 
               write(luout,*) "restarting -cosmo-"
               write(luout,'(1x,A,I5)') 
     >              "number of -cosmo- surface points = ",nefc
            else
               write(luout,*) "starting -cosmo-"
            end if
         end if

         value = BA_alloc_get(mt_dbl,nefc,'qc0',qc0(2),qc0(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'qc1',qc1(2),qc1(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'BQc',BQc(2),BQc(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'qs',  qs(2), qs(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,
     >                                  'qold',qold(2),qold(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nefc*nefc,'AA',AA(2),AA(1))

         if (cosmo2_on) 
     >      value = value.and.BA_alloc_get(mt_dbl,3*nefc,
     >                                    'rqion',rqion(2),rqion(1))
         if (.not. value) 
     >     call errquit('nwpw_cosmo_init:out of heap memory',0,MA_ERR)


         if(.not.btdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,dbl_mb(qc0(1))))
     >     call errquit('nwpw_cosmo_qc0: failed',0,rtdb_err)

         if(.not.btdb_get(rtdb,'cosmo:efcs',mt_dbl,nefc,dbl_mb(qs(1))))
     >     call errquit('nwpw_cosmo_qsurfaces: failed',0,rtdb_err)


        
*        **** set qthreshhold ****
         if (.not.btdb_get(rtdb,'nwpw_cosmo:qthresh',
     >                   mt_dbl,1,qthresh)) then
            if (cosmo1_on) then
               qthresh = 0.5d-8
            else
               qthresh = 0.5d-11
            end if
         end if
  

         Qt=.true.
         Qion = 0.0d0
         Qelc = 0.0d0
         call Parallel_shared_vector_zero(.true.,nefc,dbl_mb(qold(1)))
         if (firstrestart) then
            if(.not.btdb_get(rtdb,'nwpw_cosmo:qc1',mt_dbl,
     >                       nefc,dbl_mb(qc1(1))))
     >        call errquit('nwpw_cosmo_qc: failed',0,rtdb_err)

            if(.not.btdb_get(rtdb,'nwpw_cosmo:qold',mt_dbl,
     >                       nefc,dbl_mb(qold(1))))
     >        call errquit('nwpw_cosmo_qold: failed',0,rtdb_err)

            if(.not.btdb_get(rtdb,'nwpw_cosmo:AA',mt_dbl,
     >                       nefc*nefc,dbl_mb(AA(1))))
     >        call errquit('nwpw_cosmo_AA: failed',0,rtdb_err)

            if(.not.btdb_get(rtdb,'nwpw_cosmo:BQc',mt_dbl,
     >                       nefc,dbl_mb(BQc(1))))
     >        call errquit('nwpw_cosmo_BQc: failed',0,rtdb_err)

            if(.not.btdb_get(rtdb,'nwpw_cosmo:Qt',mt_log,1,Qt)) 
     >        Qt=.true.
            if(.not.btdb_get(rtdb,'nwpw_cosmo:Qion',mt_dbl,1,Qion))
     >        Qion = 0.0d0
            if(.not.btdb_get(rtdb,'nwpw_cosmo:Qelc',mt_dbl,1,Qelc))
     >        Qelc = 0.0d0
         end if
         Qt0 = Qt

         if (oprint) then
            write(luout,'(1x,A,E12.6)') "-cosmo- q thresh = ",qthresh
            if (.not.Qt) write(luout,*) "-cosmo- q fixed"
         end if

         Qset = (cosmo_sccor.gt.10)
         Qlmd = (cosmo_sccor.gt.11)

         call Parallel_taskid(taskid)
         call Parallel_np(np)
         nindx  = 0
         pcount = 0
         do i=1,nefc
            do j=1,i
               if (pcount.eq.taskid) nindx = nindx + 1
               pcount = mod(pcount+1,np)
            end do
         end do

         value = BA_alloc_get(mt_int,nindx,'i_indx',i_indx(2),i_indx(1))
         value = value.and.
     >           BA_alloc_get(mt_int,nindx,'j_indx',j_indx(2),j_indx(1))
         if (.not. value) 
     >     call errquit('nwpw_cosmo_init:out of heap memory',1,MA_ERR)

         nindx  = 0
         pcount = 0
         do i=1,nefc
            do j=1,i
               if (pcount.eq.taskid) then
                  int_mb(i_indx(1)+nindx) = i
                  int_mb(j_indx(1)+nindx) = j
                  nindx = nindx + 1
               end if
               pcount = mod(pcount+1,np)
            end do
         end do

         if (cosmo2_on) then
            call nwpw_cosmo_qcoords(dbl_mb(rqion(1)))
            call nwpw_cosmo_generate_Ainv(dbl_mb(rqion(1)))
         end if
         
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_end               *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_end()
      implicit none
      
#include "bafdecls.fh"
#include "btdb.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      logical value


c     == Clean up cosmo stuff, if necessary ==
      if (cosmo_on) then
         value = btdb_put(rtdb,'nwpw_cosmo:Qt',mt_log,1,(Qt.and.Qt0))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:Qion',mt_dbl,1,Qion)
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:Qelc',mt_dbl,1,Qelc)
         value = value.and.
     >           btdb_put(rtdb,'cosmo:efcz',mt_dbl,nefc,dbl_mb(qc0(1)))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:qc1',mt_dbl,
     >                    nefc,dbl_mb(qc1(1)))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:qold',mt_dbl,
     >                    nefc,dbl_mb(qold(1)))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:BQc',mt_dbl,
     >                    nefc,dbl_mb(BQc(1)))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:AA',mt_dbl,
     >                    nefc*nefc,dbl_mb(AA(1)))
         value = value.and.
     >           btdb_put(rtdb,'nwpw_cosmo:invAsum',mt_dbl,1,invAsum)
         value = value.and.BA_free_heap(i_indx(2))
         value = value.and.BA_free_heap(j_indx(2))
         value = value.and.BA_free_heap(qc0(2))
         value = value.and.BA_free_heap(qc1(2))
         value = value.and.BA_free_heap(BQc(2))
         value = value.and.BA_free_heap(qs(2))
         value = value.and.BA_free_heap(qold(2))
         value = value.and.BA_free_heap(AA(2))
         if (cosmo2_on) value = value.and.BA_free_heap(rqion(2))
         if (.not. value) 
     >      call errquit('nwpw_cosmo_end:heap memory',0,MA_ERR)

         if (firstrun) call cosmo_tidy(rtdb)
         cosmo_on = .false.
         cosmo_phase = 1
         nefc = 0
         gcds = 0.0d0
      endif
      return
      end

      real*8 function nwpw_cosmo_screen()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

      nwpw_cosmo_screen = screen
      return
      end


      real*8 function nwpw_cosmo_screen_qc0(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

      nwpw_cosmo_screen_qc0 = screen*dbl_mb(qc0(1)+i-1)
      return
      end


      real*8 function nwpw_cosmo_qc0(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

      nwpw_cosmo_qc0 = dbl_mb(qc0(1)+i-1)
      return
      end

      real*8 function nwpw_cosmo_qc1(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

      nwpw_cosmo_qc1 = dbl_mb(qc1(1)+i-1)
      return
      end


      integer function nwpw_cosmo_BQ_ptr()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

      nwpw_cosmo_BQ_ptr = BQc(1)
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_nq                *
*     *                                         *
*     *******************************************
      integer function nwpw_cosmo_nq()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo_nq = 0
      if (cosmo_on) then
         nwpw_cosmo_nq = nefc
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_phase             *
*     *                                         *
*     *******************************************
      integer function nwpw_cosmo_phase()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo_phase = 0
      if (cosmo_on) then
         nwpw_cosmo_phase = cosmo_phase
      endif
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_firsttime         *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_firsttime()
      implicit none

#include "nwpw_cosmo.fh"

      logical firsttime

      firsttime = .false.

      if (cosmo_on) then
         firsttime = (.not.firstrestart)
      end if

      nwpw_cosmo_firsttime = firsttime
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_smd               *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_smd()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo_smd = .false.
      if (cosmo_on) then
         nwpw_cosmo_smd = do_cosmo_smd
      end if
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_off               *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_off()
      implicit none

#include "nwpw_cosmo.fh"

      cosmo_on = .false.
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_on                *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_on()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo_on = cosmo_on
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo1_on               *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo1_on()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo1_on = cosmo1_on
      return
      end 

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo2_on               *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo2_on()
      implicit none

#include "nwpw_cosmo.fh"

      nwpw_cosmo2_on = cosmo2_on
      return
      end 

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_reset_Qt          *
*     *                                         *
*     *******************************************
      subroutine  nwpw_cosmo_reset_Qt()
      implicit none

#include "bafdecls.fh"
#include "nwpw_cosmo.fh"

      call Parallel_shared_vector_zero(.true.,nefc,dbl_mb(qold(1)))
      Qt = .true.
      return
      end 


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qcoords           *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qcoords(qcoord)
      implicit none
      real*8 qcoord(*)

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "nwpw_cosmo.fh"

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,qcoord))
     >     call errquit('nwpw_cosmo_qcoord: failed',0,rtdb_err)
      endif
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qcoord_ptr        *
*     *                                         *
*     *******************************************
      integer function nwpw_cosmo_qcoord_ptr()
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "nwpw_cosmo.fh"

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion

      if (cosmo1_on) then
         nwpw_cosmo_qcoord_ptr = ion_rion_ptr()+3*ion_nion()
      else
         nwpw_cosmo_qcoord_ptr = rqion(1)
      endif
      return
      end





*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_qa_indexing       *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_qa_indexing(indx)
      implicit none
      integer indx

#include "errquit.fh"
#include "bafdecls.fh"
#include "btdb.fh"
#include "nwpw_cosmo.fh"

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efciat',mt_int,nefc,indx))
     >     call errquit('nwpw_cosmo_qa_indexing: failed',0,rtdb_err)
      endif
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_Eqq               *
*     *                                         *
*     *******************************************
      real*8 function nwpw_cosmo_Eqq()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     *** local variables ****
      integer i,j,k
      real*8  e,x,y,z,r
      real*8  sai,saj
      integer qcoord_ptr

      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion
      external ion_rion_ptr,ion_nion

      if (cosmo1_on) then
         qcoord_ptr = ion_rion_ptr()+3*ion_nion()
      else
         qcoord_ptr = rqion(1)
      end if
      
!$OMP MASTER
            sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sb)
      do k=1,nindx
         i = int_mb(i_indx(1)+k-1)
         j = int_mb(j_indx(1)+k-1)
         !sai = dbl_mb(qs(1)+i-1)
         !saj = dbl_mb(qs(1)+j-1)
         if (i.eq.j) then
            sb = sb + 0.5d0*3.8d0
     >                /dsqrt(dbl_mb(qs(1)+i-1))
     >                *dbl_mb(qc0(1)+i-1)
     >                *dbl_mb(qc1(1)+i-1)
         else
            x = dbl_mb(qcoord_ptr+3*(i-1))
     >        - dbl_mb(qcoord_ptr+3*(j-1))
            y = dbl_mb(qcoord_ptr+3*(i-1)+1)
     >        - dbl_mb(qcoord_ptr+3*(j-1)+1)
            z = dbl_mb(qcoord_ptr+3*(i-1)+2)
     >        - dbl_mb(qcoord_ptr+3*(j-1)+2)
            r = dsqrt(x*x + y*y + z*z)
            sb = sb + dbl_mb(qc1(1)+i-1)
     >               *dbl_mb(qc0(1)+j-1)/r
         end if
      end do
!$OMP END DO
      call Parallel_SumAll(sb)

      nwpw_cosmo_Eqq = sb*screen
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_EQionq            *
*     *                                         *
*     *******************************************
      real*8 function nwpw_cosmo_EQionq()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     *** local variables ***
      integer i,j,nion
      real*8 x,y,z,r,qi,qj
      integer qcoord_ptr,rcoord_ptr

      !common /eiond_ejtmp/ e

      real*8 e1,e2
      common /pack_conjg_tcc_indot_sasb/ e1,e2
     

*     *** external functions ***
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm
      integer  ion_rion_ptr
      external ion_rion_ptr
      real*8   psp_zv,nwpw_cosmo_qc0
      external psp_zv,nwpw_cosmo_qc0
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      nion = ion_nion()
      if (cosmo1_on) then
         qcoord_ptr = ion_rion_ptr() + 3*nion
      else
         qcoord_ptr = rqion(1)
      end if
      rcoord_ptr = ion_rion_ptr()

!$OMP MASTER
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO reduction(+:e2)
      do i=1,nefc
         !qi = nwpw_cosmo_qc0(i-nion)
         !qi = nwpw_cosmo_screen_qc0(i-nion)
         qi = nwpw_cosmo_screen_qc0(i)
         do j=1,nion
           qj = psp_zv(ion_katm(j))
           x=dbl_mb(qcoord_ptr+3*(i-1))  -dbl_mb(rcoord_ptr+3*(j-1))
           y=dbl_mb(qcoord_ptr+3*(i-1)+1)-dbl_mb(rcoord_ptr+3*(j-1)+1)
           z=dbl_mb(qcoord_ptr+3*(i-1)+2)-dbl_mb(rcoord_ptr+3*(j-1)+2)
           r= dsqrt(x*x + y*y + z*z)
           if (r.gt.1.0d-6) then
              e2 = e2 + qi*qj/r
           end if
        end do
      end do
!$OMP END DO

      nwpw_cosmo_EQionq = e2
      return
      end 


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_FQionq            *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_FQionq(fion)
      implicit none
      real*8 fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     *** local variables ***
      integer tid,nthr
      integer i,j,nion
      real*8 x,y,z,r,qi,qj,fx,fy,fz,v
      integer qcoord_ptr,rcoord_ptr

      !common /eiond_ejtmp/ e

      real*8 e1,e2
      common /pack_conjg_tcc_indot_sasb/ e1,e2


*     *** external functions ***
      integer  ion_nion,ion_katm,Parallel_threadid,Parallel_nthreads
      external ion_nion,ion_katm,Parallel_threadid,Parallel_nthreads
      integer  ion_rion_ptr
      external ion_rion_ptr
      real*8   psp_zv,nwpw_cosmo_qc0
      external psp_zv,nwpw_cosmo_qc0
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nion = ion_nion()
      if (cosmo1_on) then
         qcoord_ptr = ion_rion_ptr() + 3*nion
      else
         qcoord_ptr = rqion(1)
      end if
      rcoord_ptr = ion_rion_ptr()

      do j=tid+1,nion,nthr
         qj = psp_zv(ion_katm(j))
         fx = 0.0d0
         fy = 0.0d0
         fz = 0.0d0
         do i=1,nefc
            qi = nwpw_cosmo_screen_qc0(i)
            x=dbl_mb(qcoord_ptr+3*(i-1))  -dbl_mb(rcoord_ptr+3*(j-1))
            y=dbl_mb(qcoord_ptr+3*(i-1)+1)-dbl_mb(rcoord_ptr+3*(j-1)+1)
            z=dbl_mb(qcoord_ptr+3*(i-1)+2)-dbl_mb(rcoord_ptr+3*(j-1)+2)
            r= dsqrt(x*x + y*y + z*z)
            if (r.gt.1.0d-6) then
               v =  qi*qj/r**3
               fx = fx + (x*v)
               fy = fy + (y*v)
               fz = fz + (z*v)
            end if
         end do
         fion(1,j) = fion(1,j) - fx
         fion(2,j) = fion(2,j) - fy
         fion(3,j) = fion(3,j) - fz
      end do

      return
      end





*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_EQelcq            *
*     *                                         *
*     *******************************************
      real*8 function nwpw_cosmo_EQelcq()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     *** local variables ***
      integer i,j,nion
      real*8 x,y,z,r,qi,qj
      integer qcoord_ptr,rcoord_ptr

      !common /eiond_ejtmp/ e

      real*8 e1,e2
      common /pack_conjg_tcc_indot_sasb/ e1,e2


*     *** external functions ***
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm
      integer  ion_rion_ptr
      external ion_rion_ptr
      real*8   pspw_getQtot_APC,nwpw_cosmo_qc0
      external pspw_getQtot_APC,nwpw_cosmo_qc0
      real*8   nwpw_cosmo_screen_qc0
      external nwpw_cosmo_screen_qc0

      nion = ion_nion()
      if (cosmo1_on) then
         qcoord_ptr = ion_rion_ptr() + 3*nion
      else
         qcoord_ptr = rqion(1)
      end if
      rcoord_ptr = ion_rion_ptr()

!$OMP MASTER
      e2 = 0.0d0
!$OMP END MASTER
!$OMP BARRIER

!$OMP DO reduction(+:e2)
      do i=1,nefc
         !qi = nwpw_cosmo_qc0(i-nion)
         !qi = nwpw_cosmo_screen_qc0(i-nion)
         qi = nwpw_cosmo_screen_qc0(i)
         do j=1,nion
           qj = pspw_getQtot_APC(j)
           x=dbl_mb(qcoord_ptr+3*(i-1))  -dbl_mb(rcoord_ptr+3*(j-1))
           y=dbl_mb(qcoord_ptr+3*(i-1)+1)-dbl_mb(rcoord_ptr+3*(j-1)+1)
           z=dbl_mb(qcoord_ptr+3*(i-1)+2)-dbl_mb(rcoord_ptr+3*(j-1)+2)
           r= dsqrt(x*x + y*y + z*z)
           if (r.gt.1.0d-6) then
              e2 = e2 + qi*qj/r
           end if
        end do
      end do
!$OMP END DO


      nwpw_cosmo_EQelcq = e2
      return
      end




*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_calc_resdual      *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_calc_residual(rion,SA,BQ,q,rq)
      implicit none
      real*8 rion(3,*),SA(*),BQ(*),q(*)
      real*8 rq(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      integer i,j,k
      integer tid,nthreads
      real*8  x,y,z,r,alpha

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      alpha = 0.005d0
      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()

      call Parallel_shared_vector_zero(.true.,nefc,rq)
      do k=tid+1,nindx,nthreads
         i = int_mb(i_indx(1)+k-1)
         j = int_mb(j_indx(1)+k-1)
         if (i.eq.j) then
            rq(i) = rq(i) + alpha*(-BQ(i)-3.8d0/dsqrt(SA(i))*q(i))+q(i)
            !rq(i) = rq(i) + BQ(i) + q(i)
         else
            x = rion(1,i)-rion(1,j)
            y = rion(2,i)-rion(2,j)
            z = rion(3,i)-rion(3,j)
            r = dsqrt(x*x + y*y + z*z)
            rq(i) = rq(i) - alpha*q(j)/r
            rq(j) = rq(j) - alpha*q(i)/r
         end if
      end do
      call Parallel_Vector_SumAll(nefc,rq)

      !do i=1,nefc
      !   write(*,*) "i,rq=",i,rq(i),q(i)
      !end do

      return
      end



*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_calc_Aq           *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_calc_Aq(rion,SA,q,Aq)
      implicit none
      real*8 rion(3,*),SA(*),q(*)
      real*8 Aq(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      integer i,j,k
      integer tid,nthreads
      real*8  x,y,z,r

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()

      call Parallel_shared_vector_zero(.true.,nefc,Aq)
      do k=tid+1,nindx,nthreads
         i = int_mb(i_indx(1)+k-1)
         j = int_mb(j_indx(1)+k-1)
         if (i.eq.j) then
            Aq(i) = Aq(i) + 3.8d0/dsqrt(SA(i))*q(i)
            !rq(i) = rq(i) + BQ(i) + q(i)
         else
            x = rion(1,i)-rion(1,j)
            y = rion(2,i)-rion(2,j)
            z = rion(3,i)-rion(3,j)
            r = dsqrt(x*x + y*y + z*z)
            Aq(i) = Aq(i) + q(j)/r
            Aq(j) = Aq(j) + q(i)/r
         end if
      end do
      call Parallel_Vector_SumAll(nefc,Aq)

      return
      end



*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_solve_q           *
*     *                                         *
*     *******************************************
      subroutine nwpw_cosmo_solve_q()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      integer it,qcoord_ptr,ia,ne(2),ispin
      real*8  esum,esum1,esum2,qsum1,qsum2,scal,sa,error

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion,ion_nkatm0,ion_natm
      external ion_rion_ptr,ion_nion,ion_nkatm0,ion_natm
      real*8   psp_zv
      external psp_zv


      !write(*,*) "nwpw_cosmo_solve_q Qt,Qlmd,Qset=",Qt,Qlmd,Qset
      if (Qt.and.Qt0) then
         if ((dabs(Qion).lt.1.0d-6).or.(dabs(Qelc).lt.1.0d-6)) then
            call Dneall_ispin(ispin)
            call Dneall_ne(ne)
            Qelc = -(ne(1)+ne(ispin))
            Qion = 0.0d0
            do ia=1,ion_nkatm0()
              Qion = Qion + ion_natm(ia)*psp_zv(ia)
            end do
         end if
      
         !**** really simple solver to start ****
         qcoord_ptr = ion_rion_ptr()+3*ion_nion()
         call DGEMM_OMP('N','N',nefc,1,nefc,
     >                 (-1.0d0),
     >                 dbl_mb(AA(1)),nefc,
     >                 dbl_mb(BQc(1)),nefc,
     >                 (0.0d0),
     >                 dbl_mb(qc0(1)),nefc)

!$OMP MASTER
         qsum2 = 0.0d0
         do it=1,nefc
            qsum2 = qsum2+(dbl_mb(qc0(1)+it-1)-dbl_mb(qold(1)+it-1))**2
         end do
         error = dsqrt(qsum2)/dble(nefc)
         if ((qthresh.gt.0.0).and.(error.lt.qthresh)) Qt = .false.
!$OMP END MASTER


         call Parallel_shared_vector_copy(.true.,nefc,
     >                       dbl_mb(qc0(1)),
     >                       dbl_mb(qold(1)))

         call Parallel_shared_vector_copy(.true.,nefc,
     >                    dbl_mb(qc0(1)),
     >                    dbl_mb(qc1(1)))


*        !**** rescale the total surface charge ****
         if (Qset) then
!$OMP MASTER
            if (.not.Qlmd) then
               qsum1 = 0.0d0
               sa = 0.0d0
               do it=1,nefc
                  qsum1 = qsum1 + dbl_mb(qc0(1)+it-1)
                  sa = sa + dbl_mb(qs(1)+it-1)
               end do
               !write(*,*) "qsum1=",qsum1,-(Qelc+Qion)
               qsum2 = -(Qelc+Qion) - qsum1
               qsum2 = qsum2/sa
               do it=1,nefc
                  dbl_mb(qc0(1)+it-1) = dbl_mb(qc0(1)+it-1) 
     >                                + qsum2*dbl_mb(qs(1)+it-1)
               end do

            !*** lagrange ***
            else
               call nwpw_cosmo_lagrange(-(Qelc+Qion),invAsum,
     >                               nefc,dbl_mb(AA(1)),dbl_mb(qc0(1)))
            end if
        

c  E = 0.5QCQ + QBq + 0.5qAq + Lambda(sum(i) qi - charge)
c    = 0.5*Sum(I,J) Q(I) C(I,J) Q(J)
c    +  sum(i,J) q(i) B(i,J) Q(J)
c    +  0.5*Sum(i,j) q(i) A(i,j) q(j)
c    +  Lambda*(Sum(i) q(i) - charge)
c  dE/dqk = sum(J) B(k,J)*Q(J)
c         + sum(j) A(k,j) q(j)
c         + Lambda
c
c   0 = Q*B + A*q + Lambda*vec(1)
c    q = -1*Am*B*Q - Am*vec(Lambda)
c
c    qt = -Am*B*Q
c
c    sum(i) q(i) = Charge = sum(i) (qt(i) - Lambda*sum(j) Am(j,i))
c
c    lambda = (-Charge + sum(i) qt(i)) / (sum(i,j) Am(i,j) ) 
c  Bq + Aq = 0
c         qsum1 = 0.0d0
c         qsum2 = 0.0d0
c         do it=1,nefc
c            if (dbl_mb(qc0(1)+it-1).gt.0.0d0) then
c               qsum1 = qsum1 + dbl_mb(qc0(1)+it-1)
c            else
c               qsum2 = qsum2 + dbl_mb(qc0(1)+it-1)
c            end if
c         end do
c         if (dabs(Qion+Qelc).gt.1.0d-6) then
c            scal = -(Qelc+Qion)/(qsum1+qsum2)
c            call dscal(nefc,scal,dbl_mb(qc0(1)),1)
c         else
c            if (dabs(qsum1).gt.dabs(qsum2)) then
c               scal = -qsum1/qsum2
c               do it=1,nefc
c                  if (dbl_mb(qc0(1)+it-1).lt.0.0d0) then
c                     dbl_mb(qc0(1)+it-1) = scal*dbl_mb(qc0(1)+it-1)
c                  end if
c               end do
c            else
c               scal = -qsum2/qsum1
c               do it=1,nefc
c                  if (dbl_mb(qc0(1)+it-1).gt.0.0d0) then
c                     dbl_mb(qc0(1)+it-1) = scal*dbl_mb(qc0(1)+it-1)
c                  end if
c               end do
c            end if
c         end if
!$OMP END MASTER

         end if !** if Qset **



c      do it=1,5000
c         call nwpw_cosmo_calc_residual(dbl_mb(qcoord_ptr),
c     >                                 dbl_mb(qs(1)),
c     >                                 dbl_mb(BQc(1)),
c     >                                 dbl_mb(qc0(1)),
c     >                                 dbl_mb(qc1(1)))
c         call Parallel_shared_vector_copy(.true.,nefc,
c     >                                    dbl_mb(qc1(1)),
c     >                                    dbl_mb(qc0(1)))
c      end do

c      call nwpw_cosmo_calc_Aq(dbl_mb(qcoord_ptr),
c     >                        dbl_mb(qs(1)),
c     >                        dbl_mb(qc0(1)),
c     >                        dbl_mb(qc1(1)))

c      esum = 0.0d0
c      esum1 = 0.0d0
c      do it=1,nefc
c         write(*,*) "BQ,Aq,q,BQ+Aq=",dbl_mb(BQc(1)+it-1),
c     >                             dbl_mb(qc1(1)+it-1),
c     >                             dbl_mb(qc0(1)+it-1),
c     >         dbl_mb(qc1(1)+it-1)+dbl_mb(BQc(1)+it-1)
c         esum  = esum + dbl_mb(qc0(1)+it-1)*dbl_mb(qc1(1)+it-1)
c         esum1 = esum1 + dbl_mb(qc0(1)+it-1)*dbl_mb(BQc(1)+it-1)
c      end do
c      write(*,*) "-0.5*qAq=",-0.5d0*esum
c      write(*,*) "qBQ + 0.5*qAq=",esum1+0.5d0*esum
c      write(*,*) "qBQ + qAq=",esum1+esum

         if (cosmo1_on) call ion_cosmo_reset_katm(dbl_mb(qc0(1)))


      !** if Qt**
      else

         if (firstrestart.and.cosmo1_on) then
            call ion_cosmo_reset_katm(dbl_mb(qc0(1)))
         end if 

      end if !** else Qt**

      firstrestart = .false.

      return
      end

      subroutine nwpw_cosmo_lagrange(charge,invAsum,N,Ainv,q)
      implicit none
      real*8  charge,invAsum
      integer N
      real*8 Ainv(N,N),q(N)

*     **** local variables ****
      integer i,j
      real*8 lmd,tsum

      lmd = -charge
      do i =1,N
        lmd = lmd + q(i)
      end do
      lmd = lmd/invAsum

      do i =1,N
         tsum = 0.0d0
         do j=1,N
            tsum = tsum + Ainv(i,j)
         end do
         q(i) = q(i) - lmd*tsum
      end do

      return
      end 

      subroutine nwpw_cosmo_generate_Ainv(rion)
      implicit none
      real*8 rion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_cosmo.fh"

*     **** local variables ****
      integer i,j,k
      integer tid,nthreads
      real*8  x,y,z,r

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()


      call Parallel_shared_vector_zero(.true.,nefc*nefc,dbl_mb(AA(1)))
      do k=tid+1,nindx,nthreads
         i = int_mb(i_indx(1)+k-1)
         j = int_mb(j_indx(1)+k-1)
         if (i.eq.j) then
            dbl_mb(AA(1)+(i-1)+(j-1)*nefc)
     >         = 3.8d0/dsqrt(dbl_mb(qs(1)+i-1))
         else
            x = rion(1,i)-rion(1,j)
            y = rion(2,i)-rion(2,j)
            z = rion(3,i)-rion(3,j)
            r = dsqrt(x*x + y*y + z*z)
            dbl_mb(AA(1)+(i-1)+(j-1)*nefc) = 1.0d0/r
            dbl_mb(AA(1)+(j-1)+(i-1)*nefc) = 1.0d0/r
         end if
      end do
      call Parallel_Vector_SumAll(nefc*nefc,dbl_mb(AA(1)))

      call nwpw_matrix_invert(nefc,dbl_mb(AA(1)),nefc)

!$OMP MASTER
      do i=1,nefc*nefc
         invAsum = invAsum + dbl_mb(AA(1)+i-1)
      end do
!$OMP END MASTER

      return
      end



