*
* $Id$
*

*     ***********************************
*     *                                 *
*     *         paw_geodesic_init       *
*     *                                 *
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine paw_geodesic_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax,nelc1
      

*     **** external functions **** 
      integer  paw_psi_ne
      external paw_psi_ne


      call Pack_npack(1,npack1)
      nemax = paw_psi_ne(1)+paw_psi_ne(2)
      nelc1 = paw_psi_ne(1)

      value = BA_alloc_get(mt_dcpl,npack1*nemax,
     >                     'U',U(2),U(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,2*nelc1*nelc1,
     >                     'Vt',Vt(2),Vt(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,nemax,
     >                     'S',S(2),S(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      return
      end

*     ***********************************
*     *                                 *
*     *         paw_geodesic_finalize   *
*     *                                 *
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine paw_geodesic_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S

*     **** local variables ****
      logical value

      value = BA_free_heap(S(2))
      value = value.and.
     >        BA_free_heap(Vt(2))
      value = value.and.
     >        BA_free_heap(U(2))
      if (.not. value) call errquit('error freeing of heap memory',0,
     &       MA_ERR)

      return
      end



*     ***********************************
*     *                                 *
*     *         paw_geodesic_start      *
*     *                                 *
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
*     Exit  - max_sigma:
*             dE:
*     Uses - geodesic common block
*

      subroutine paw_geodesic_start(A,max_sigma,dE)
      implicit none
#include "errquit.fh"
      complex*16 A(*)
      real*8     max_sigma,dE

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ms,i,j,n
c      integer ierr
      integer shift1,shift2
      integer tmp1(2),tmp2(2),V(2)
      integer indx1,indx2,indx3
      real*8 sum1

*     **** external functions **** 
      integer  paw_psi_ispin,paw_psi_ne
      real*8   paw_electron_eorbit,paw_psi_CheckOrtho
      external paw_psi_ispin,paw_psi_ne
      external paw_electron_eorbit,paw_psi_CheckOrtho

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)
      nemax = paw_psi_ne(1) + paw_psi_ne(2)


*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(4*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'V',V(2),V(1))

      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

*     ***********************
*     **** HomeGrown SVD ****
*     ***********************
      
      max_sigma = 0.0d0
      call dcopy(nemax,0.0d0,0,dbl_mb(S(1)),1)
      do ms=1,paw_psi_ispin()
         shift1 = (ms-1)*paw_psi_ne(1)
         shift2 = (ms-1)*paw_psi_ne(1)*paw_psi_ne(1)
         n = paw_psi_ne(ms)
         if (n.le.0) go to 30

         call paw_Pack_SVD(npack1,n,
     >                  A(1+shift1*npack1),
     >                  dcpl_mb(U(1)+shift1*npack1),
     >                  dbl_mb(S(1)+shift1),
     >                  dbl_mb(V(1)+shift2))


         !!call paw_psi_CheckOrtho2(npack1,n,dcpl_mb(U(1)+shift1*npack1))
         !sum1=paw_psi_CheckOrtho(npack1,n,dcpl_mb(U(1)+shift1*npack1))
         !write(*,*) "<U|O|U>:",sum1



         do i=1,n
            if (dabs(dbl_mb(S(1)+i-1+shift1)).gt.max_sigma) 
     >        max_sigma = dabs(dbl_mb(S(1)+i-1+shift1))

         end do

         call paw_calc_Vt(n,dbl_mb(Vt(1)+shift2),dbl_mb(V(1)+shift2))

  30     continue         
      end do

      !write(*,*) "<psi|O|U>:"
      !call paw_psi_1check_Tangent(dcpl_mb(U(1)))
      !write(*,*)


*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*paw_electron_eorbit(A)

*     **** deallocate tmp space ****     
      value = value.and.
     >        BA_pop_stack(V(2))
      value = value.and.
     >        BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end

      subroutine paw_calc_Vt(n,A,B)
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
*     *                                 *
*     *         paw_geodesic_get        *
*     *                                 *
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine paw_geodesic_get(t,Yold,Ynew)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S

*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2,shift
      real*8     scal1,scal2
      real*8     zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)
      real*8     sum1,sum2

      integer    taskid, MASTER
      parameter  (MASTER=0)

*     **** external functions ****
      integer  paw_psi_ispin,paw_psi_ne
      real*8   paw_psi_CheckOrtho
      external paw_psi_ispin,paw_psi_ne
      external paw_psi_CheckOrtho

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0
      call Pack_npack(1,npack1)
      nemax = paw_psi_ne(1) + paw_psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpC',
     >                    tmpC(2),tmpC(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpS',
     >                    tmpS(2),tmpS(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)


      do ms = 1,paw_psi_ispin()
        n      = paw_psi_ne(ms)
        if (n.eq.0) go to 100  !*** ferromagnetic case ***
        shift1 = (ms-1)*paw_psi_ne(1)
        shift2 = (ms-1)*paw_psi_ne(1)*paw_psi_ne(1)

        call SCtimesVtrans(t,n,dbl_mb(S(1)+shift1),
     >                    dbl_mb(Vt(1)+shift2),
     >                    dbl_mb(tmp1(1)+shift2),
     >                    dbl_mb(tmp3(1)+shift2),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))

        call dgemm('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2),n,
     >              dbl_mb(tmp1(1)+shift2),n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2),n)



        call dgemm('N','N',2*npack1,n,n,
     >              one,
     >              Yold(1+shift1*npack1),2*npack1,
     >              dbl_mb(tmp2(1)+shift2),n,
     >              zero,
     >              Ynew(1+shift1*npack1),2*npack1)

        call dgemm('N','N',2*npack1,n,n,
     >              one,
     >              dcpl_mb(U(1)+shift1*npack1),2*npack1,
     >              dbl_mb(tmp3(1)+shift2),n,
     >              one,
     >              Ynew(1+shift1*npack1),2*npack1)



*     **** Orthonormality Check ****     
      !call paw_psi_CheckOrtho2(npack1,n,Ynew(1+shift1*npack1))
      sum1=paw_psi_CheckOrtho(npack1,n,Ynew(1+shift1*npack1))
      if (dabs(sum1).gt.1.0d-10) then
c         call Parallel_taskid(taskid)
c         if (taskid.eq.MASTER) then
c           write(*,*) "Gram-Schmidt performed in paw_geodesic_get"
c           write(*,*) "sum1:",sum1
c         end if
         call paw_psi_MakeOrtho(npack1,n,Ynew(1+shift1*npack1))
      end if

 100    continue

      end do



*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *                                 *
*     *         paw_geodesic_transport  *
*     *                                 *
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine paw_geodesic_transport(t,Yold,Ynew)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S

*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2,shift
      real*8     zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)

*     **** external functions ****
      integer  paw_psi_ispin,paw_psi_ne
      external paw_psi_ispin,paw_psi_ne

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = paw_psi_ne(1) + paw_psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpC',
     >                    tmpC(2),tmpC(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpS',
     >                    tmpS(2),tmpS(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)


      do ms = 1,paw_psi_ispin()
        n      = paw_psi_ne(ms)
        if (n.eq.0) go to 200
        shift1 = (ms-1)*paw_psi_ne(1)
        shift2 = (ms-1)*paw_psi_ne(1)*paw_psi_ne(1)

        call SCtimesVtrans2(t,n,dbl_mb(S(1)+shift1),
     >                    dbl_mb(Vt(1)+shift2),
     >                    dbl_mb(tmp1(1)+shift2),
     >                    dbl_mb(tmp3(1)+shift2),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))



        call dgemm('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2), n,
     >              dbl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2), n)




        call dgemm('N','N',2*npack1,n,n,
     >              (-1.0d0),
     >              Yold(1+shift1*npack1),2*npack1,
     >              dbl_mb(tmp2(1)+shift2),n,
     >              zero,
     >              Ynew(1+shift1*npack1),2*npack1)

        call dgemm('N','N',2*npack1,n,n,
     >              one,
     >              dcpl_mb(U(1)+shift1*npack1),2*npack1,
     >              dbl_mb(tmp3(1)+shift2),n,
     >              one,
     >              Ynew(1+shift1*npack1),2*npack1)

 200    continue

      end do

*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end






*     ***********************************
*     *                                 *
*     *         paw_geodesic_Gtransport *
*     *                                 *
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine paw_geodesic_Gtransport(t,Yold,tG)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 tG(*)

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / paw_geodesic_block / U,Vt,S


*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2
      real*8     zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)
      integer    tmpC(2),tmpS(2)

*     **** external functions ****
      integer  paw_psi_ispin,paw_psi_ne
      external paw_psi_ispin,paw_psi_ne

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = paw_psi_ne(1) + paw_psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpC',
     >                    tmpC(2),tmpC(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(nemax),'tmpS',
     >                    tmpS(2),tmpS(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      do ms = 1,paw_psi_ispin()
        n      = paw_psi_ne(ms)
        if (n.eq.0) go to 300   !*** ferromagnetic case ***
        shift1 = (ms-1)*paw_psi_ne(1)
        shift2 = (ms-1)*paw_psi_ne(1)*paw_psi_ne(1)

        call Pack_ccm_dot(1,n,
     >             dcpl_mb(U(1)+shift1*npack1),
     >             tG(1+shift1*npack1),
     >             dbl_mb(tmp2(1)+shift2))

        call SCtimesVtrans3(t,n,dbl_mb(S(1)+shift1),
     >                    dbl_mb(tmp2(1)+shift2),
     >                    dbl_mb(tmp1(1)+shift2),
     >                    dbl_mb(tmp3(1)+shift2),
     >                    dbl_mb(tmpC(1)),
     >                    dbl_mb(tmpS(1)))


        call dgemm('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2), n,
     >              dbl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2), n)



        call dgemm('N','N',2*npack1,n,n,
     >              (-1.0d0),
     >              Yold(1+shift1*npack1),2*npack1,
     >              dbl_mb(tmp2(1)+shift2),n,
     >              one,
     >              tG(1+shift1*npack1),2*npack1)

        call dgemm('N','N',2*npack1,n,n,
     >              (-1.0d0),
     >              dcpl_mb(U(1)+shift1*npack1),2*npack1,
     >              dbl_mb(tmp3(1)+shift2),n,
     >              one,
     >              tG(1+shift1*npack1),2*npack1)

 300    continue

      end do !*** ms ***

*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *                                 *
*     *            paw_Pack_SVD         *
*     *                                 *
*     ***********************************

      subroutine paw_Pack_SVD(npack,ne,A,U,S,V)
      implicit none
#include "errquit.fh"
      integer npack,ne
      complex*16 A(npack,ne)
      complex*16 U(npack,ne)
      real*8     S(ne)
      real*8     V(ne,ne)

#include "bafdecls.fh"

*     **** local variables ****
      logical value
      integer i,j,ierr
      integer tmp1(2),tmp2(2)
c     real*8 sum1,sumv

      real*8 ddot
      external ddot

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(4*ne*ne),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(ne),'tmp2',tmp2(2),tmp2(1))

      if (.not. value)
     >   call errquit('out of stack memory',0, MA_ERR)


*     ***** generate V and Sigma^2 ****
      !call Grsm_ggm_sym_dot(npack,ne,A,A,V)
      !call Pack_ccm_sym_dot(1,ne,A,A,V)
      call paw_overlap_matrix_gen(ne,ne,A,A,V)


      !call eigen(ne,ne,V,S,dbl_mb(tmp1(1)),ierr)
      call dsyev('V','U',ne,
     >           V,ne,
     >           S,
     >           dbl_mb(tmp1(1)),(4*ne*ne),
     >           ierr)


*     **** generate U*Sigma ****
      call dgemm('N','N',2*npack,ne,ne,
     >              (1.0d0),
     >              A,2*npack,
     >              V,ne,
     >              (0.0d0),
     >              U,2*npack)
      !call Grsm_gmg_Mul(npack,ne,A,V,U)


*     **** normalize U*sigma, and sqrt(S^2) ****
      do i=1,ne
         !call Pack_cc_idot(1,U(1,i),U(1,i),dbl_mb(tmp2(1)+i-1))
        call paw_overlap_matrix_gen(1,1,U(1,i),U(1,i),
     >                              dbl_mb(tmp2(1)+i-1))
      end do
      !call D3dB_Vector_SumAll(ne,dbl_mb(tmp2(1)))

c     do i=1,ne
c        dbl_mb(tmp2(1)+i-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+i-1))
c        call Pack_c_SMul(1,dbl_mb(tmp2(1)+i-1),
c    >                    U(1,i),U(1,i))
c
c     end do
      do i=1,ne
         dbl_mb(tmp2(1)+i-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+i-1))
      end do

      do i=1,ne
         call dscal(2*npack,dbl_mb(tmp2(1)+i-1),U(1,i),1)
      end do


      do i=1,ne
        if (S(i).lt.0.0d0) S(i) = dabs(S(i))
        S(i) = dsqrt(S(i))
      end do


c      sumv = 0.0d0
c      do i=1,ne
c          sumv = sumv + ddot(ne,V(1,i),1,V(1,i),1)
c      end do
c
c       call Grsm_gg_trace(npack,ne,U,U,sum1)
c       write(*,*) "Tr(U), Tr(V):",sum1,sumv



      value = value.and.BA_pop_stack(tmp2(2))
      value = value.and.BA_pop_stack(tmp1(2))
      if (.not. value)
     > call errquit('error popping stack memory',0, MA_ERR)

      return
      end
