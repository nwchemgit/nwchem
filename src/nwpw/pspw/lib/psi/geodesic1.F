*
* $Id$
*

*     ***********************************
*     *					*
*     *		geodesic1_init		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic1_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic1_block / U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax,nelc1
      

*     **** external functions **** 
      integer  psi_ne
      external psi_ne


      call Pack_npack(1,npack1)
      nemax = psi_ne(1)+psi_ne(2)
      nelc1 = psi_ne(1)

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
*     *					*
*     *		geodesic1_finalize	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic1_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic1_block / U,Vt,S

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
*     *					*
*     *		geodesic1_start		*
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
*     Exit  - max_sigma:
*             dE:
*     Uses - geodesic common block
*

      subroutine geodesic1_start(A,max_sigma,dE)
      implicit none
#include "errquit.fh"
      complex*16 A(*)
      real*8     max_sigma,dE

#include "bafdecls.fh"

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic1_block / U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ms,i,j,n
c      integer ierr
      integer shift1,shift2
      integer tmp1(2),tmp2(2),V(2)

*     **** external functions **** 
      integer  psi_ispin,psi_ne
      real*8   electron_eorbit
      external psi_ispin,psi_ne
      external electron_eorbit

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)


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
      do ms=1,psi_ispin()
         shift1 = (ms-1)*psi_ne(1)
         shift2 = (ms-1)*psi_ne(1)*psi_ne(1)
         n = psi_ne(ms)

         call Pack_SVD(npack1,n,
     >                  A(1+shift1*npack1),
     >                  dcpl_mb(U(1)+shift1*npack1),
     >                  dbl_mb(S(1)+shift1),
     >                  dbl_mb(V(1)+shift2))

         do i=1,n
            if (dabs(dbl_mb(S(1)+i-1+shift1)).gt.max_sigma) 
     >        max_sigma = dabs(dbl_mb(S(1)+i-1+shift1))

         end do

*        *** calculate Vt ***
         do i=1,n
         do j=1,n
             dbl_mb(Vt(1)+j-1+(i-1)*n + shift2) 
     >     =  dbl_mb(V(1)+i-1+(j-1)*n + shift2) 
         end do
         end do
         
      end do

*     **** calculate 2*<A|H|psi> ****
      dE = 2.0d0*electron_eorbit(A)

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


*     ***********************************
*     *                                 *
*     *            Pack_SVD             *
*     *                                 *
*     ***********************************

      subroutine Pack_SVD(npack,ne,A,U,S,V)
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
      integer i,ierr
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
      call Grsm_ggm_dot(npack,ne,A,A,V)
      call eigen(ne,ne,V,S,dbl_mb(tmp1(1)),ierr)

         
*     **** generate U*Sigma ****
      call Grsm_gmg_Mul(npack,ne,A,V,U)


*     **** normalize U*sigma, and sqrt(S^2) ****
      do i=1,ne
         call Pack_cc_idot(1,U(1,i),U(1,i),dbl_mb(tmp2(1)+i-1))
      end do
      call D3dB_Vector_SumAll(ne,dbl_mb(tmp2(1)))
      
      do i=1,ne
         dbl_mb(tmp2(1)+i-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+i-1))

         call Pack_c_SMul(1,dbl_mb(tmp2(1)+i-1),
     >                    U(1,i),U(1,i))
 
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



      value = value.and.
     >        BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)
       
      return
      end

*     ***********************************
*     *					*
*     *		geodesic1_get		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic1_get(t,Yold,Ynew)
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
      common / geodesic1_block / U,Vt,S

*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2
      real*8     scal1,scal2
      real*8     zero,one
      integer    tmp1(2),tmp2(2),tmp3(2)
      real*8     sum1,sum2

      integer    taskid, MASTER
      parameter  (MASTER=0)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0
      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)


      do ms = 1,psi_ispin()
        n      = psi_ne(ms)
        if (n.eq.0) go to 100  !*** ferromagnetic case ***
        shift1 = (ms-1)*psi_ne(1)
        shift2 = (ms-1)*psi_ne(1)*psi_ne(1)

        do j=1,n
          scal1 = dcos(dbl_mb(S(1)+j-1+shift1)*t)
          scal2 = dsin(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,n
           dbl_mb(tmp1(1)+j-1+(k-1)*n+shift2) 
     >          = scal1 * dbl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2) 
     >          = scal2 * dbl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do

        end do

        call DGEMM('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2),n,
     >              dbl_mb(tmp1(1)+shift2),n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2),n)


        do k=1,n
          call dcopy(2*npack1,zero,0,Ynew(1+(k-1+shift1)*npack1),1)
          do j=1,n
            call Pack_cc_daxpy(1,
     >                         dbl_mb(tmp2(1)+j-1+(k-1)*n+shift2),
     >                         Yold(1+(j-1+shift1)*npack1),
     >                         Ynew(1+(k-1+shift1)*npack1))
            call Pack_cc_daxpy(1,
     >                         dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2),
     >                         dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >                         Ynew(1+(k-1+shift1)*npack1))

          end do
        end do

 100    continue

      end do


*     **** Orthonormality Check ****     
      call Grsm_gg_trace(npack1,nemax,Ynew,Ynew,sum2)
  
      sum1 = dble(nemax)
      if (dabs(sum2-sum1).gt.1.0d-10) then
         call Parallel_taskid(taskid)
c         if (taskid.eq.MASTER) 
c    >      write(*,*) 
c    >      "Warning: Gram-Schmidt being performed in geodesic_get"
         call Grsm_g_MakeOrtho(npack1,psi_ne(1),Ynew)
         if (psi_ispin().gt.1) then
           call Grsm_g_MakeOrtho(npack1,psi_ne(2),
     >                           Ynew(1 + psi_ne(1)*npack1))
         end if
      end if



*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmp3(2))
      value = value.and.
     >        BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end

*     ***********************************
*     *					*
*     *		geodesic1_transport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic1_transport(t,Yold,Ynew)
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
      common / geodesic1_block / U,Vt,S

*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2
      real*8     zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)


      do ms = 1,psi_ispin()
        n      = psi_ne(ms)
        if (n.eq.0) go to 200
        shift1 = (ms-1)*psi_ne(1)
        shift2 = (ms-1)*psi_ne(1)*psi_ne(1)

        do j=1,n
          scal1 = dbl_mb(S(1)+j-1+shift1)
     >            *dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = dbl_mb(S(1)+j-1+shift1)
     >            *dcos(dbl_mb(S(1)+j-1+shift1)*t)

          do k=1,n
           dbl_mb(tmp1(1)+j-1+(k-1)*n+shift2) 
     >       = scal1 * dbl_mb(Vt(1)+j-1+(k-1)*n+shift2)
           dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2) 
     >       = scal3 * dbl_mb(Vt(1)+j-1+(k-1)*n+shift2)
          end do
        end do

        call DGEMM('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2), n,
     >              dbl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2), n)


        do k=1,n
          call dcopy(2*npack1,zero,0,Ynew(1+(k-1+shift1)*npack1),1)
          do j=1,n
            call Pack_cc_daxpy(1,
     >              (-1.0d0*dbl_mb(tmp2(1)+j-1+(k-1)*n+shift2)),
     >              Yold(1+(j-1+shift1)*npack1),
     >              Ynew(1+(k-1+shift1)*npack1))
            call Pack_cc_daxpy(1,
     >              dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2),
     >              dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >              Ynew(1+(k-1+shift1)*npack1))

          end do
        end do

 200    continue

      end do

*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmp3(2))
      value = value.and.
     >        BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end


*     ***********************************
*     *					*
*     *		geodesic1_Gtransport	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic1_Gtransport(t,Yold,tG)
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
      common / geodesic1_block / U,Vt,S


*     **** local variables ****
      logical    value
      integer    npack1,nemax
      integer    ms,j,k,n,shift1,shift2
      real*8     zero,one
      real*8     scal1,scal3
      integer    tmp1(2),tmp2(2),tmp3(2)

*     **** external functions ****
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne

      call nwpw_timing_start(10)
      zero = 0.0d0
      one  = 1.0d0

      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)

*     **** allocate tmp space ****
      value = BA_push_get(mt_dbl,(2*nemax*nemax),'tmp1',
     >                    tmp1(2),tmp1(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_dbl,(2*nemax*nemax),'tmp3',
     >                    tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('out of stack memory',0, MA_ERR)

      do ms = 1,psi_ispin()
        n      = psi_ne(ms)
        if (n.eq.0) go to 300   !*** ferromagnetic case ***
        shift1 = (ms-1)*psi_ne(1)
        shift2 = (ms-1)*psi_ne(1)*psi_ne(1)

        call Grsm_ggm_dot(npack1,n,
     >             dcpl_mb(U(1)+shift1*npack1),
     >             tG(1+shift1*npack1),
     >             dbl_mb(tmp1(1)))
        call Grsm_ggm_dot(npack1,n,
     >            dcpl_mb(U(1)+shift1*npack1),
     >            tG(1+shift1*npack1),
     >            dbl_mb(tmp3(1)))

        do j=1,n
          scal1 =       dsin(dbl_mb(S(1)+j-1+shift1)*t)
          scal3 = 1.0d0-dcos(dbl_mb(S(1)+j-1+shift1)*t)
          do k=1,n
           dbl_mb(tmp1(1)+j-1+(k-1)*n+shift2) 
     >        = scal1 * dbl_mb(tmp1(1)+j-1+(k-1)*n + shift2)
           dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2) 
     >        = scal3 * dbl_mb(tmp3(1)+j-1+(k-1)*n + shift2)
          end do
        end do
        
        call DGEMM('T','N',n,n,n,one,
     >              dbl_mb(  Vt(1)+shift2), n,
     >              dbl_mb(tmp1(1)+shift2), n,
     >              zero,
     >              dbl_mb(tmp2(1)+shift2), n)


        do k=1,n
          do j=1,n
            call Pack_cc_daxpy(1,
     >             (-1.0d0*dbl_mb(tmp2(1)+j-1+(k-1)*n+shift2)),
     >             Yold(1+(j-1+shift1)*npack1),
     >               tG(1+(k-1+shift1)*npack1))
            call Pack_cc_daxpy(1,
     >             (-1.0d0*dbl_mb(tmp3(1)+j-1+(k-1)*n+shift2)),
     >             dcpl_mb(U(1)+(j-1+shift1)*npack1),
     >             tG(1+(k-1+shift1)*npack1))
          end do
        end do

 300    continue

      end do !*** ms ***

*     **** deallocate tmp space ****     
      value = BA_pop_stack(tmp3(2))
      value = value.and.
     >        BA_pop_stack(tmp2(2))
      value = value.and.
     >        BA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack memory',0, MA_ERR)

      call nwpw_timing_end(10)

      return
      end
