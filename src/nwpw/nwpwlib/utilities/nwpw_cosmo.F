
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

*     **** local variables ****
      logical  omaster,oprint,value
      integer  taskid,np,i,j,pcount
      integer  MASTER
      parameter(MASTER=0)

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,nefc,rtdb
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

*     ***** external functions *****
      logical  control_print
      external control_print
      integer  Parallel_threadid
      external Parallel_threadid

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

      if (.not.btdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) 
     >   cosmo_on = .false.

      if(cosmo_on) then
         call cosmo_initialize(rtdb,geom,-1,oprint)
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

         screen = (dielec - 1.0d0)/(dielec+0.5d0)

         if(.not.btdb_get(rtdb,'cosmo:nefc',mt_int,1,nefc)) nefc = 0

         value = BA_alloc_get(mt_dbl,nefc,'qc0',qc0(2),qc0(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'qc1',qc1(2),qc1(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'BQc',BQc(2),BQc(1))
         value = value.and.BA_alloc_get(mt_dbl,nefc,'qs',  qs(2), qs(1))
         value = value.and.
     >           BA_alloc_get(mt_dbl,nefc*nefc,'AA',AA(2),AA(1))
         if (.not. value) 
     >     call errquit('nwpw_cosmo_init:out of heap memory',0,MA_ERR)

         if(.not.btdb_get(rtdb,'cosmo:efcz',mt_dbl,nefc,dbl_mb(qc0(1))))
     >     call errquit('nwpw_cosmo_qs: failed',0,rtdb_err)

         if(.not.btdb_get(rtdb,'cosmo:efcs',mt_dbl,nefc,dbl_mb(qs(1))))
     >     call errquit('nwpw_cosmo_qsurfaces: failed',0,rtdb_err)


         Qset = .false.
         Qion = 0.0d0
         Qelc = 0.0d0

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
     >          BA_alloc_get(mt_int,nindx,'j_indx',j_indx(2),j_indx(1))
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
#include "errquit.fh"

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

      logical value

c     == Clean up cosmo stuff, if necessary ==
      if (cosmo_on) then
         call cosmo_tidy(rtdb)
         cosmo_on = .false.
         cosmo_phase = 1
         nefc = 0
         gcds = 0.0d0
         value = BA_free_heap(i_indx(2))
         value = value.and.BA_free_heap(j_indx(2))
         value = value.and.BA_free_heap(qc0(2))
         value = value.and.BA_free_heap(qc1(2))
         value = value.and.BA_free_heap(BQc(2))
         value = value.and.BA_free_heap(qs(2))
         value = value.and.BA_free_heap(AA(2))
         if (.not. value) 
     >      call errquit('nwpw_cosmo_end:heap memory',0,MA_ERR)
      endif
      return
      end

      real*8 function nwpw_cosmo_screen_qc0(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

      nwpw_cosmo_screen_qc0 = screen*dbl_mb(qc0(1)+i-1)
      return
      end


      real*8 function nwpw_cosmo_qc0(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

      nwpw_cosmo_qc0 = dbl_mb(qc0(1)+i-1)
      return
      end

      real*8 function nwpw_cosmo_qc1(i)
      implicit none
      integer i

#include "bafdecls.fh"
#include "errquit.fh"

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

      nwpw_cosmo_qc1 = dbl_mb(qc1(1)+i-1)
      return
      end


      integer function nwpw_cosmo_BQ_ptr()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

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

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

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

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_phase = 0
      if (cosmo_on) then
         nwpw_cosmo_phase = cosmo_phase
      endif
      return
      end


*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_smd               *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_smd()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_smd = .false.
      if (cosmo_on) then
         nwpw_cosmo_smd = do_cosmo_smd
      end if
      return
      end

*     *******************************************
*     *                                         *
*     *            nwpw_cosmo_on                *
*     *                                         *
*     *******************************************
      logical function nwpw_cosmo_on()
      implicit none

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      nwpw_cosmo_on = cosmo_on
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

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      if (cosmo_on) then
         if(.not.btdb_get(rtdb,'cosmo:efcc',mt_dbl,3*nefc,qcoord))
     >     call errquit('nwpw_cosmo_qcoord: failed',0,rtdb_err)
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

*     **** nwpw_cosmo common block ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

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

*     *** local variables ****
      integer i,j,k
      real*8  e,x,y,z,r

      real*8 sa,sb
      common /pack_conjg_tcc_indot_sasb/ sa,sb

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd

      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

*     **** external functions ****
      real*8   ion_rion
      external ion_rion

      
!$OMP MASTER
            sb = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sb)
      do k=1,nindx
         i = int_mb(i_indx(1)+k-1)
         j = int_mb(j_indx(1)+k-1)
         if (i.eq.j) then
            sb = sb + 0.5d0*3.8d0
     >                /dsqrt(dbl_mb(qs(1)+i-1))
     >                *dbl_mb(qc0(1)+i-1)**2
         else
            x = ion_rion(1,i)-ion_rion(1,j)
            y = ion_rion(2,i)-ion_rion(2,j)
            z = ion_rion(3,i)-ion_rion(3,j)
            r = dsqrt(x*x + y*y + z*z)
            sb = sb + dbl_mb(qc0(1)+i-1)*dbl_mb(qc0(1)+j-1)/r
         end if
      end do
!$OMP END DO
      call Parallel_SumAll(sb)

      nwpw_cosmo_Eqq = sb*screen
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

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd
      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

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

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd
      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

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

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd
      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

*     **** local variables ****
      integer it,qcoord_ptr,ia,ne(2),ispin
      real*8  esum,esum1,esum2,qsum1,qsum2,scal

*     **** external functions ****
      integer  ion_rion_ptr,ion_nion,ion_nkatm0,ion_natm
      external ion_rion_ptr,ion_nion,ion_nkatm0,ion_natm
      real*8   psp_zv
      external psp_zv
      


      if (.not.Qset) then
         Qset = .true.
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
     >              (-1.0d0),
     >              dbl_mb(AA(1)),nefc,
     >              dbl_mb(BQc(1)),nefc,
     >              (0.0d0),
     >              dbl_mb(qc0(1)),nefc)

c*     !**** rescale the total surface charge ****
c!$OMP MASTER
c      qsum1 = 0.0d0
c      qsum2 = 0.0d0
c      do it=1,nefc
c         if (dbl_mb(qc0(1)+it-1).gt.0.0d0) then
c            qsum1 = qsum1 + dbl_mb(qc0(1)+it-1)
c         else
c            qsum2 = qsum2 + dbl_mb(qc0(1)+it-1)
c         end if
c      end do
c      if (dabs(Qion+Qelc).gt.1.0d-6) then
c         scal = -(Qelc+Qion)/(qsum1+qsum2)
c         call dscal(nefc,scal,dbl_mb(qc0(1)),1)
c      else
c         if (dabs(qsum1).gt.dabs(qsum2)) then
c            scal = -qsum1/qsum2
c            do it=1,nefc
c               if (dbl_mb(qc0(1)+it-1).lt.0.0d0) then
c                  dbl_mb(qc0(1)+it-1) = scal*dbl_mb(qc0(1)+it-1)
c               end if
c            end do
c         else
c            scal = -qsum2/qsum1
c            do it=1,nefc
c               if (dbl_mb(qc0(1)+it-1).gt.0.0d0) then
c                  dbl_mb(qc0(1)+it-1) = scal*dbl_mb(qc0(1)+it-1)
c               end if
c            end do
c         end if
c      end if
c      qsum1 = 0.0d0
c      do it=1,nefc
c         qsum1 = qsum1 + dbl_mb(qc0(1)+it-1)
c      end do
c!$OMP END MASTER



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

      call ion_cosmo_reset_katm(dbl_mb(qc0(1)))

      return
      end

      subroutine nwpw_cosmo_generate_Ainv(rion)
      implicit none
      real*8 rion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_cosmo common blocks ****
      logical cosmo_on,do_cosmo_smd
      integer cosmo_phase,rtdb,nefc
      real*8  screen,dielec,gcds
      common /nwpw_cosmo_blk/ screen,dielec,gcds,cosmo_phase,rtdb,nefc,
     >                        cosmo_on,do_cosmo_smd
      integer i_indx(2),j_indx(2),nindx
      common /nwpw_cosmo_blk2/ i_indx,j_indx,nindx

      integer qc0(2),qc1(2),BQc(2),qs(2),AA(2)
      real*8  Qion,Qelc
      logical Qset
      common /nwpw_cosmo_blk3/ Qion,Qelc,qc0,qc1,BQc,qs,AA,Qset

*     **** local variables ****
      integer i,j,k
      integer tid,nthreads
      real*8  x,y,z,r

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid      = Parallel_threadid()
      nthreads = Parallel_nthreads()

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

c      write(*,*) "A matrix="
c      do i=1,nefc
c        write(*,*) (dbl_mb(AA(1)+(i-1)+(j-1)*nefc),j=1,nefc)
c        write(*,*)
c      end do

      call nwpw_matrix_invert(nefc,dbl_mb(AA(1)),nefc)

c      write(*,*) "inverse A matrix="
c      do i=1,nefc
c        write(*,*) (dbl_mb(AA(1)+(i-1)+(j-1)*nefc),j=1,nefc)
c        write(*,*)
c      end do

      return
      end

