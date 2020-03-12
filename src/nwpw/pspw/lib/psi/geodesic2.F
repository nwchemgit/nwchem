*
* $Id$
*

*     ***********************************
*     *					*
*     *		geodesic2_init		*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic2_init()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1

*     **** external functions **** 
      integer  psi_ne,psi_neq
      external psi_ne,psi_neq
      logical  Dneall_m_allocate,Dneall_4m_allocate
      external Dneall_m_allocate,Dneall_4m_allocate

      call Pack_npack(1,npack1)
c      nemax = psi_neq(1)+psi_neq(2)
c      ne(1) = psi_neq(1)
c      ne(2) = psi_neq(2)

      value = BA_alloc_get(mt_dcpl,npack1*(psi_neq(1)+psi_neq(2)),
     >                     'Hold',Hold(2),Hold(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*(psi_neq(1)+psi_neq(2)),
     >                     'Q',Q(2),Q(1))
      value = value.and.Dneall_m_allocate(0,R)
      value = value.and.Dneall_m_allocate(0,A)
      value = value.and.Dneall_4m_allocate(0,V)
      value = value.and.Dneall_4m_allocate(0,W)

      value = value.and.
     >        BA_alloc_get(mt_dbl,2*(psi_ne(1)+psi_ne(2)),
     >                     'S',S(2),S(1))
      if (.not. value) call errquit('geodesic2_init:allocating heap',0,
     &       MA_ERR)

      return
      end

*     ***********************************
*     *					*
*     *		geodesic2_finalize	*
*     *					*
*     ***********************************
*
*     Uses - geodesic common block
*

      subroutine geodesic2_finalize()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value

      logical  Dneall_m_free
      external Dneall_m_free

      value =           BA_free_heap(S(2))   
      value = value.and.Dneall_m_free(W)
      value = value.and.Dneall_m_free(V)
      value = value.and.Dneall_m_free(A)
      value = value.and.Dneall_m_free(R)
      value = value.and.BA_free_heap(Q(2))
      value = value.and.BA_free_heap(Hold(2))
      if (.not. value) 
     >  call errquit('geodesic2_finalize:freeing heap memory',0,
     &       MA_ERR)

      return
      end



*     ***********************************
*     *					*
*     *		geodesic2_start		*
*     *					*
*     ***********************************
*
*     This routine determines the pxp matrices R and YA, and
* the orthogonal nxp matrix Q.   Q and R are determined from 
* the QR decomposition of the projected direction (I-YY^t)H, and
* YH is defined as the Lagrange Multiplier pxp matrix Y^tH.  
*
*     Uses - geodesic2 common block
*
      subroutine geodesic2_start(Y,H,max_sigma,dE)
      implicit none
      complex*16 Y(*)
      complex*16 H(*)
      real*8     max_sigma,dE

#include "bafdecls.fh"
#include "errquit.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ms,n,ispin,ne(2)
      integer shift,shift2,i,j
      integer T(2)

*     **** external functions **** 
      logical  Dneall_4m_allocate,Dneall_m_free
      integer  psi_ispin,psi_neq
      real*8   electron_eorbit_noocc
      external Dneall_4m_allocate,Dneall_m_free
      external psi_ispin,psi_neq
      external electron_eorbit_noocc

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)
      nemax = psi_neq(1) + psi_neq(2)
      ispin = psi_ispin()
      ne(1) = psi_neq(1)
      ne(2) = psi_neq(2)

*     **** allocate tmp space ****
      value = Dneall_4m_allocate(0,T)
      if (.not. value) 
     >   call errquit('geodesic2_start: pushing stack',0, MA_ERR)

*     **** Hold <-- H ****
      !call dcopy(2*npack1*nemax,H,1,dcpl_mb(Hold(1)),1)
      call Parallel_shared_vector_copy(.true.,2*npack1*nemax,
     >                                 H,dcpl_mb(Hold(1)))

*     **** calculate A=<Y|H> ****
      call Dneall_ffm_Multiply(0,Y,H,npack1,dbl_mb(A(1)))

*     **** calculate Q=(I-YYt)H - should not be necessary but just in case ****
      call Dneall_fmf_Multiply(0,Y,npack1,
     >                         dbl_mb(A(1)),1.0d0,
     >                         dcpl_mb(Q(1)),0.0d0)
      call DAXPY_OMP(2*npack1*nemax,(-1.0d0),H,1,dcpl_mb(Q(1)),1)
      call DSCAL_OMP(2*npack1*nemax,(-1.0d0),dcpl_mb(Q(1)),1)

*     **** calculate QR using Modified Gram-Schmidt ****
      call Dneall_fm_QR(0,dcpl_mb(Q(1)),npack1,dbl_mb(R(1)))

*     **** generate T from A and R ****
*       -     -
*  T = |A, -R^t|
*      |R,  0  |
*       -     -
      call Dneall_AR_to_4m(0,dbl_mb(A(1)),dbl_mb(R(1)),dbl_mb(T(1)))

*     **** Factor T--> V,W,and S ****
      call Dneall_4m_FactorSkew(0, dbl_mb(T(1)),
     >                             dbl_mb(V(1)),
     >                             dbl_mb(W(1)),
     >                             dbl_mb(S(1)))

*     **** calculate dE ****
      dE = 2.0d0*electron_eorbit_noocc(H)


*     **** deallocate tmp space ****     
      value = Dneall_m_free(T)
      if (.not. value) 
     > call errquit('geodesic2_start:popping stack',1, MA_ERR)

      call nwpw_timing_end(10)
      return
      end

*     ***********************************
*     *					*
*     *	      geodesic2_generate_T	*
*     *					*
*     ***********************************
*
*     This routine determines T.  T is defined
* to be a 2nx2n skew symmetric matrix.
*
*       -     -
*  T = |A, -R^t|
*      |R,  0  |
*       -     -
*
*     Entry - n: dimension of matrices A and R
*             A: an nxn skew symmetric matrix
*             R: an nxn matrix
*     Exit - T: a 2nx2n skew symmetric matrix

      subroutine geodesic2_generate_T(n,A,R,T)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 R(n,n)
      real*8 T(2*n,2*n)

*     **** local variables ****
      integer i,j

      call dcopy(4*n*n,0.0d0,0,T,1)
      
*     **** copy A to upper-left of T ****
      do j=1,n
      do i=1,n
         T(i,j) = A(i,j)
      end do
      end do

*     **** copy R to lower-left of T ****
      do j=1,n
      do i=1,n
         T(i+n,j) = R(i,j)
      end do
      end do

*     **** copy -R^t to upper-right of T ****
      do j=1,n
      do i=1,n
         T(i,j+n) = -R(j,i)
      end do
      end do

      return
      end 


*     ***********************************
*     *					*
*     *		geodesic2_get		*
*     *					*
*     ***********************************
*
*     This routine calculates 
*
*   Ynew = Yold*M(t) + Q*N(t)
*
*   where
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*
      subroutine geodesic2_get(t,Yold,Ynew)
      implicit none
      real*8     t
      complex*16 Yold(*)
      complex*16 Ynew(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "geodesic2.fh"

*     **** local variables ****
      logical value
      integer npack1,ispin,ne(2)
      integer ms,n,j,k,shift1,shift2
      integer MM(2),NN(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack


      call Pack_npack(1,npack1)

*     **** allocate tmp space ****
      value =           Dneall_m_push_get(0,MM)
      value = value.and.Dneall_m_push_get(0,NN)
      if (.not. value) 
     >   call errquit('geodesic2_get: pushing stack',0, MA_ERR)

      call geodesic2_get_MandN(t,dbl_mb(MM(1)),dbl_mb(NN(1)))

      call Dneall_fmf_Multiply(0,
     >                         Yold,npack1,
     >                         dbl_mb(MM(1)),
     >                         1.0d0,
     >                         Ynew,
     >                         0.0d0)
      call Dneall_fmf_Multiply(0,
     >                         dcpl_mb(Q(1)),npack1,
     >                         dbl_mb(NN(1)),
     >                         1.0d0,
     >                         Ynew,
     >                         1.0d0)

*     **** deallocate tmp space ****     
      value =           Dneall_m_pop_stack(NN)
      value = value.and.Dneall_m_pop_stack(MM)
      if (.not. value) 
     > call errquit('geodesic2_get:popping stack',1, MA_ERR)

      return
      end


*     ***********************************
*     *					*
*     *		geodesic2_transport	*
*     *					*
*     ***********************************
*
*     This routine calculates 
*
*   Hnew = Hold*M(t)    + Yold*R^t*N(t)
*
*   where
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*
      subroutine geodesic2_transport(t,Yold,Hnew)
      implicit none
#include "errquit.fh"
      real*8     t
      complex*16 Yold(*)
      complex*16 Hnew(*)

#include "bafdecls.fh"
#include "geodesic2.fh"


*     **** local variables ****
      logical value
      integer npack1
      integer MM(2),NN(2),TT(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack


      call Pack_npack(1,npack1)

*     **** allocate tmp space ****
      value =           Dneall_m_push_get(0,MM)
      value = value.and.Dneall_m_push_get(0,NN)
      value = value.and.Dneall_m_push_get(0,TT)
      if (.not. value) 
     >   call errquit('geodesic2_transport: pushing stack',0, MA_ERR)

      call geodesic2_get_MandN(t,dbl_mb(MM(1)),dbl_mb(NN(1)))

*     **** TT(t) = -R^t*NN(t) ****
      call Dneall_mmm_Multiply2ab(0,
     >                          dbl_mb(R(1)),
     >                          dbl_mb(NN(1)),
     >                          (-1.0d0),
     >                          dbl_mb(TT(1)),
     >                          0.0d0)


*     *** Hnew <-- Hold*M(t) + Yold*TT(t) ***
      call Dneall_fmf_Multiply(0,
     >                         dcpl_mb(Hold(1)),npack1,
     >                         dbl_mb(MM(1)),
     >                         1.0d0,
     >                         Hnew,
     >                         0.0d0)
      call Dneall_fmf_Multiply(0,
     >                         Yold,npack1,
     >                         dbl_mb(TT(1)),
     >                         1.0d0,
     >                         Hnew,
     >                         1.0d0)

*     **** deallocate tmp space ****    
      value =           Dneall_m_pop_stack(TT)
      value = value.and.Dneall_m_pop_stack(NN)
      value = value.and.Dneall_m_pop_stack(MM)
      if (.not. value) 
     > call errquit('geodesic2_transport:popping stack',1, MA_ERR)


      return
      end

      subroutine geodesic2_checkMN(name,n,MM,NN)
      implicit none
      character*(*) name
      integer n
      real*8 MM(n,n),NN(n,n)

*     **** local variables ****
      integer i,j,k
      real*8 tmp(n,n),sum

      do j=1,n
      do i=1,n
        sum = 0.0d0
        do k=1,n
          sum = sum + MM(k,i)*NN(k,j)
        end do
        tmp(i,j) = sum
      end do
      end do

      write(*,*)
      write(*,*) "checkMN:",name
      do i=1,n
        write(*,'(22F8.3)') (tmp(i,j),j=1,n)
      end do
      write(*,*)
      return
      end

*     ***********************************
*     *					*
*     *	      geodesic2_get_MandN	*
*     *					*
*     ***********************************
*
*   This routine returns 
*        -    -               - -
*       | M(t) | = Exp(t*T)* | I |
*       | N(t) |             | 0 |
*        -    -               - -
*   where 
* 
*      T =  U*Sigma*U^H, with U=(V+iW)
*
*      is a skew matrix that is decomposed into V,W,and Sigma
*      
      subroutine geodesic2_get_MandN(t,M,N)
      implicit none
      real*8 t
      real*8 M(*),N(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "geodesic2.fh"


*     **** local variables ****
      logical value
      integer ispin,ne(2)
      integer ms,i,j,nn,shift,shift2,shift3
      integer RR(2)

*     **** external functions ****
      logical  Dneall_4m_push_get,Dneall_m_pop_stack
      external Dneall_4m_push_get,Dneall_m_pop_stack

*     **** allocate tmp space ****
      value = Dneall_4m_push_get(0,RR)
      if (.not. value) 
     >   call errquit('geodesic2_get_MandN: pushing stack',0, MA_ERR)

      call Dneall_4m_RotateSkew(0,t,
     >                          dbl_mb(V(1)),
     >                          dbl_mb(W(1)),
     >                          dbl_mb(S(1)),
     >                          dbl_mb(RR(1)) )

      call Dneall_4m_to_MN(0,dbl_mb(RR(1)),M,N)

*     **** deallocate tmp space ****     
      value = Dneall_m_pop_stack(RR)
      if (.not. value) 
     > call errquit('geodesic2_get_MandN:popping stack',1, MA_ERR)

      return
      end
