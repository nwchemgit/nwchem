*
* $Id: geodesic2.F,v 1.1 2002-06-12 00:59:53 bylaska Exp $
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

#include "mafdecls.fh"

*     **** geodesic2 common block ***
      integer Q(2)  ! npack1xne(1) matrix        and npack*ne(2) matrix
      integer U(2)  ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer Vt(2) ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer S(2)  ! diagonal matrices of SVD decomposition 
      common / geodesic2_block / Q,U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax,nelc1
      

*     **** external functions **** 
      integer  psi_ne
      external psi_ne

      call Pack_npack(1,npack1)
      nemax = psi_ne(1)+psi_ne(2)
      nelc1 = psi_ne(1)

      value = MA_alloc_get(mt_dcpl,npack1*nemax,
     >                     'Q',Q(2),Q(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,8*nelc1*nelc1,
     >                     'U',U(2),U(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,8*nelc1*nelc1,
     >                     'Vt',Vt(2),Vt(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*nemax,
     >                     'S',S(2),S(1))
      if (.not. value) call errquit('geodesic2_init:allocating heap',0)

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

#include "mafdecls.fh"

*     **** geodesic2 common block ***
      integer Q(2)  ! npack1xne(1) matrix        and npack*ne(2) matrix
      integer U(2)  ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer Vt(2) ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer S(2)  ! diagonal matrices of SVD decomposition 
      common / geodesic2_block / Q,U,Vt,S

*     **** local variables ****
      logical value

      value =           MA_free_heap(S(2))   
      value = value.and.MA_free_heap(Vt(2))
      value = value.and.MA_free_heap(U(2))
      value = value.and.MA_free_heap(Q(2))
      if (.not. value) 
     >  call errquit('geodesic2_finalize:freeing heap memory',0)

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

#include "mafdecls.fh"

*     **** geodesic2 common block ***
      integer Q(2)  ! npack1xne(1) matrix        and npack*ne(2) matrix
      integer U(2)  ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer Vt(2) ! (2*ne(1))*(2*ne(1)) matrix and (2*ne(2))*(2*ne(2)) matrix
      integer S(2)  ! diagonal matrices of SVD decomposition 
      common / geodesic2_block / Q,U,Vt,S

*     **** local variables ****
      logical value
      integer npack1,nemax
      integer ms,i,j,n,ispin,ne(2)
c      integer ierr
      integer shift1,shift2
      integer A(2)

*     **** external functions **** 
      integer  psi_ispin,psi_ne
      real*8   electron_eorbit
      external psi_ispin,psi_ne
      external electron_eorbit

      call nwpw_timing_start(10)
      call Pack_npack(1,npack1)
      nemax = psi_ne(1) + psi_ne(2)
      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)

*     **** allocate tmp space ****
      value = MA_push_get(mt_dbl,(ne(1)*ne(1)+ne(2)*ne(2)),
     >                    'A',A(2),A(1))
      if (.not. value) 
     >   call errquit('geodesic2_start: pushing stack',0)


*     **** calculate A=<Y|H> ****
      do ms=1,ispin
         call Grsm_ggm_dot(npack1,ne(ms),
     >                     Y(1+(ms-1)*ne(1)*npack1),
     >                     H(1+(ms-1)*ne(1)*npack1),
     >                     dbl_mb(A(1)+(ms-1)*ne(1)*ne(1)))
      end do
      stop


*     **** calculate QR using Modified Gram-Schmidt ****

*     **** generate T ****
*     **** SVD T--> U*S*Vt ****


*     **** deallocate tmp space ****     
      value = value.and.
     >        MA_pop_stack(A(2))
      if (.not. value) 
     > call errquit('geodesic2_start:popping stack',1)

      
      call nwpw_timing_end(10)

      return
      end

