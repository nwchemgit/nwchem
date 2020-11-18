*
* $Id$
*

*     ***********************************
*     *					*
*     *		geodesic_init		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic_init()
      implicit none

#include "errquit.fh"
#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      integer npack1
      

*     **** external functions **** 
      integer  psi_ne,psi_neq
      external psi_ne,psi_neq

      logical  Dneall_m_allocate
      external Dneall_m_allocate

      call Pack_npack(1,npack1)
c      nemax = psi_ne(1)+psi_ne(2)
c      nelc1 = psi_ne(1)

      value = BA_alloc_get(mt_dcpl,npack1*(psi_neq(1)+psi_neq(2)),
     >                     'U',U(2),U(1))

      value = value.and.Dneall_m_allocate(0,Vt)
c      value = value.and.
c     >        BA_alloc_get(mt_dbl,2*nelc1*nelc1,
c     >                     'Vt',Vt(2),Vt(1))

      value = value.and.
     >        BA_alloc_get(mt_dbl,(psi_ne(1)+psi_ne(2)),
     >                     'S',S(2),S(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *		geodesic_finalize	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      logical  Dneall_m_free
      external Dneall_m_free

      value =           BA_free_heap(S(2))
      value = value.and.Dneall_m_free(Vt)
c      value = value.and.BA_free_heap(Vt(2))
      value = value.and.BA_free_heap(U(2))
      if (.not. value) call errquit('error freeing of heap memory',0,
     &       MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		geodesic_start		*
*     *					*
*     ***********************************
*
*     This routine initializes the geodesic module
* for a linesearch.  Basically this routine just
* calculates the SVD decomposition of the search direction, 
* A=HY-Y(Y^tHY) or A=(determined from CG). The only requirement
* of the search direction is that it is tangent to the direction
* spanned by Y.   It returns the maximum value in the diagonal 
* Sigma matrix, and it also returns the linegradient determined
* by the direction A.
* 
*     Entry - A: gradient
*             SA: S*gradient if paw
*     Exit  - max_sigma:
*             dE:
*             SA: S*U if paw
*     Uses - geodesic common block
*

      subroutine geodesic_start(A,max_sigma,dE)
      implicit none
      complex*16 A(*)
      real*8     max_sigma,dE

#include "bafdecls.fh"
#include "errquit.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

      integer spsi1(2),spsi2(2)
      common / psi_paw_block / spsi1,spsi2

*     **** local variables ****
      integer i,npack1,V(2),ispin,neq(2)
      real*8 de_private

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack,psp_pawexist
      integer  psi_ispin,psi_ne,psi_neq
      real*8   electron_eorbit_noocc
      external Dneall_m_push_get,Dneall_m_pop_stack,psp_pawexist
      external psi_ispin,psi_ne,psi_neq
      external electron_eorbit_noocc

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)

*     **** allocate tmp space ****
      if (.not.Dneall_m_push_get(0,V))
     >   call errquit('geodesic_start:out of stack memory',0,MA_ERR)
      
*     **** HomeGrown SVD ****
      if (psp_pawexist()) then
         ispin = psi_ispin()
         neq(1)= psi_neq(1)
         neq(2)= psi_neq(2)
         call psp_overlap_S(ispin,neq,A,dcpl_mb(spsi1(1)))
         call Dneall_f_SVD_ASA1(0,A,
     >                          dcpl_mb(spsi1(1)),
     >                          dcpl_mb(U(1)),npack1,
     >                          dbl_mb(S(1)),dbl_mb(V(1)))
         call psp_overlap_S(ispin,neq,dcpl_mb(U(1)),dcpl_mb(spsi1(1)))
         call Dneall_f_SVD_ASA2(0,dcpl_mb(U(1)),
     >                          dcpl_mb(spsi1(1)),npack1)

      else
         call Dneall_f_SVD(0,A,dcpl_mb(U(1)),npack1,
     >                     dbl_mb(S(1)),dbl_mb(V(1)))
      end if

      max_sigma = 0.0d0
      do i=1,(psi_ne(1)+psi_ne(2))
        if (dabs(dbl_mb(S(1)+i-1)).gt.max_sigma) 
     >      max_sigma = dabs(dbl_mb(S(1)+i-1))
      end do

*     **** calculate Vt ****
      call Dneall_mm_transpose(0,dbl_mb(V(1)),dbl_mb(Vt(1)))

*     **** calculate 2*<A|H|psi> ****
      de_private = 2.0d0*electron_eorbit_noocc(A)

      dE = de_private

*     **** deallocate tmp space ****   
      if (.not.Dneall_m_pop_stack(V))  
     >   call errquit('geodesic_start:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end

      subroutine pspw_calc_Vt(n,A,B)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 B(n,n)
      integer i,j

      do j=1,n
      do i=1,n
        A(i,j) = B(j,i)
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *		geodesic_get		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

      integer spsi1(2),spsi2(2)
      common / psi_paw_block / spsi1,spsi2

*     **** local variables ****
      logical    value
      integer    npack1,nemax,ispin,ne(2),neq(2),shift,ms
      real*8     zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)
c      real*8     sum1,sum2,sum3
      real*8     sum1

      real*8 sum2(2)
      common /geodescic_sum2/ sum2

      integer    taskid, MASTER
      parameter  (MASTER=0)

*     **** external functions ****
      integer  psi_ispin,psi_ne,psi_neq
      external psi_ispin,psi_ne,psi_neq
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      logical  psp_pawexist
      external psp_pawexist

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0
      call Pack_npack(1,npack1)
      ispin  = psi_ispin()
      ne(1)  = psi_ne(1)
      ne(2)  = psi_ne(2)
      neq(1) = psi_neq(1)
      neq(2) = psi_neq(2)
      nemax  = ne(1) + ne(2)

*     **** allocate tmp space ****
      value =           Dneall_m_push_get(0,tmp1)
      value = value.and.Dneall_m_push_get(0,tmp2)
      value = value.and.Dneall_m_push_get(0,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) call errquit('geodesic_get:out of stack',0,MA_ERR)


      call Dneall_mm_SCtimesVtrans(0,t,dbl_mb(S(1)),
     >                                dbl_mb(Vt(1)),
     >                                dbl_mb(tmp1(1)),
     >                                dbl_mb(tmp3(1)),
     >                                dbl_mb(tmpC(1)),
     >                                dbl_mb(tmpS(1)))


      call Dneall_mmm_Multiply2(0,dbl_mb(Vt(1)),
     >                           dbl_mb(tmp1(1)),
     >                           dbl_mb(tmp2(1)))

      call Dneall_fmf_Multiply(0,Yold,npack1,
     >                          dbl_mb(tmp2(1)),1.0d0,
     >                          Ynew,0.0d0)

      call Dneall_fmf_Multiply(0,dcpl_mb(U(1)),npack1,
     >                          dbl_mb(tmp3(1)),1.0d0,
     >                          Ynew,1.0d0)



!$OMP BARRIER
*     **** Orthonormality Check ****
      if (psp_pawexist()) then
         call psp_overlap_S(ispin,neq,Ynew,dcpl_mb(spsi1(1)))
         do ms=1,ispin
            shift = 1 + (ms-1)*neq(1)*npack1
            call Grsm_gg_itrace(npack1,neq(ms),
     >                       Ynew(shift),
     >                       dcpl_mb(spsi1(1)+shift-1),sum2(ms))
         end do
!$OMP BARRIER
         call Parallel_Vector_SumAll(ispin,sum2)
         do ms=1,ispin
            sum1 = dble(ne(ms))
            if (dabs(sum2(ms)-sum1).gt.1.0d-10) then
                shift = 1 + (ms-1)*neq(1)*npack1
               call Dneall_f_Sortho(ms,Ynew,dcpl_mb(spsi1(1)),npack1)
            end if
         end do
      else
         do ms=1,ispin
            shift = 1 + (ms-1)*neq(1)*npack1
            call Grsm_gg_itrace(npack1,neq(ms),
     >                          Ynew(shift),Ynew(shift),sum2(ms))
         end do
!$OMP BARRIER
         call Parallel_Vector_SumAll(ispin,sum2)
         do ms=1,ispin
            sum1 = dble(ne(ms))
            !write(*,*) "sum1,sum2=",sum1,sum2(ms),dabs(sum2(ms)-sum1)
            if (dabs(sum2(ms)-sum1).gt.1.0d-10) then
                shift = 1 + (ms-1)*neq(1)*npack1
            !write(*,*) "INTO GRAMSCHMIDT"
               call Dneall_f_GramSchmidt(ms,Ynew,npack1)
            !write(*,*) "OUT GRAMSCHMIDT"
            end if
         end do
      end if
!$OMP BARRIER


*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Dneall_m_pop_stack(tmp3)
      value = value.and.Dneall_m_pop_stack(tmp2)
      value = value.and.Dneall_m_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('geodesic:get:error popping stack memory',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end




*     ***********************************
*     *					*
*     *		SCtimesVtrans		*
*     *					*
*     ***********************************

      subroutine SCtimesVtrans(t,n,S,Vt,A,B,scal1,scal2)
      implicit none
      real*8 t
      integer n
      real*8  S(n),Vt(n,n)
      real*8  A(n,n),B(n,n)
      real*8 scal1(n),scal2(n)

      integer j,k

      do j=1,n
        scal1(j) = dcos(S(j)*t)
        scal2(j) = dsin(S(j)*t)
      end do

      do k=1,n
      do j=1,n
          A(j,k) = scal1(j)*Vt(j,k)
          B(j,k) = scal2(j)*Vt(j,k)
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *		geodesic_transport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic_transport(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      logical    value
      integer    npack1,nemax
      real*8     zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)

*     **** allocate tmp space ****
      value =           Dneall_m_push_get(0,tmp1)
      value = value.and.Dneall_m_push_get(0,tmp2)
      value = value.and.Dneall_m_push_get(0,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) 
     >   call errquit('geodesic_transport: out of stack',0,MA_ERR)


      call Dneall_mm_SCtimesVtrans2(0,t,dbl_mb(S(1)),
     >                                dbl_mb(Vt(1)),
     >                                dbl_mb(tmp1(1)),
     >                                dbl_mb(tmp3(1)),
     >                                dbl_mb(tmpC(1)),
     >                                dbl_mb(tmpS(1)))

      call Dneall_mmm_Multiply2(0,dbl_mb(Vt(1)),
     >                           dbl_mb(tmp1(1)),
     >                           dbl_mb(tmp2(1)))

      call Dneall_fmf_Multiply(0,Yold,npack1,
     >                          dbl_mb(tmp2(1)),-1.0d0,
     >                          Ynew,0.0d0)

      call Dneall_fmf_Multiply(0,dcpl_mb(U(1)),npack1,
     >                          dbl_mb(tmp3(1)),1.0d0,
     >                          Ynew,1.0d0)

*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Dneall_m_pop_stack(tmp3)
      value = value.and.Dneall_m_pop_stack(tmp2)
      value = value.and.Dneall_m_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('geodesic_transport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end

*     ***********************************
*     *                                 *
*     *         SCtimesVtrans2           *
*     *                                 *
*     ***********************************

      subroutine SCtimesVtrans2(t,n,S,Vt,A,B,scal1,scal2)
      implicit none
      real*8 t
      integer n
      real*8  S(n),Vt(n,n)
      real*8  A(n,n),B(n,n)
      real*8 scal1(n),scal2(n)

      integer j,k

      do j=1,n
        scal1(j) = S(j)*dsin(S(j)*t)
        scal2(j) = S(j)*dcos(S(j)*t)
      end do

      do k=1,n
      do j=1,n
          A(j,k) = scal1(j)*Vt(j,k)
          B(j,k) = scal2(j)*Vt(j,k)
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *		geodesic_Gtransport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic_Gtransport(t,Yold,tG)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 tG(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S


*     **** local variables ****
      logical    value
      integer    npack1,nemax
      real*8     zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)

*     **** allocate tmp space ****
      value =           Dneall_m_push_get(0,tmp1)
      value = value.and.Dneall_m_push_get(0,tmp2)
      value = value.and.Dneall_m_push_get(0,tmp3)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not. value) 
     >   call errquit('geodesic_Gtransport:out of stack',0,MA_ERR)

         
      call Dneall_ffm_Multiply(0,dcpl_mb(U(1)),tG,npack1,
     >                           dbl_mb(tmp2(1)))

      call Dneall_mm_SCtimesVtrans3(0,t,dbl_mb(S(1)),
     >                                dbl_mb(tmp2(1)),
     >                                dbl_mb(tmp1(1)),
     >                                dbl_mb(tmp3(1)),
     >                                dbl_mb(tmpC(1)),
     >                                dbl_mb(tmpS(1)))

      call Dneall_mmm_Multiply2(0,dbl_mb(Vt(1)),
     >                           dbl_mb(tmp1(1)),
     >                           dbl_mb(tmp2(1)))

      call Dneall_fmf_Multiply(0,Yold,npack1,
     >                          dbl_mb(tmp2(1)),(-1.0d0),
     >                          tG,1.0d0)

      call Dneall_fmf_Multiply(0,dcpl_mb(U(1)),npack1,
     >                          dbl_mb(tmp3(1)),(-1.0d0),
     >                          tG,1.0d0)

*     **** deallocate tmp space ****     
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Dneall_m_pop_stack(tmp3)
      value = value.and.Dneall_m_pop_stack(tmp2)
      value = value.and.Dneall_m_pop_stack(tmp1)
      if (.not. value) 
     > call errquit('geodesic_gtransport:error popping stack',0,MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *                                 *
*     *         SCtimesVtrans3	        *
*     *                                 *
*     ***********************************

      subroutine SCtimesVtrans3(t,n,S,Vt,A,B,scal1,scal2)
      implicit none
      real*8 t
      integer n
      real*8  S(n),Vt(n,n)
      real*8  A(n,n),B(n,n)
      real*8  scal1(n),scal2(n)

      integer j,k

      do j=1,n
        scal1(j) = dsin(S(j)*t)
        scal2(j) = 1.0d0-dcos(S(j)*t)
      end do

      do k=1,n
      do j=1,n
          A(j,k) = scal1(j)*Vt(j,k)
          B(j,k) = scal2(j)*Vt(j,k)
      end do
      end do

      return
      end



