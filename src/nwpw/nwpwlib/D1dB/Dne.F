* $Id: Dne.F,v 1.2 2006-09-20 19:18:48 bylaska Exp $
*


*     ***********************************
*     *					*
*     *	        Dne_init		*	
*     *					*
*     ***********************************

      subroutine Dne_init(ispin_in,ne_in,map_in)
      implicit none
      integer ispin_in
      integer ne_in(2)
      integer map_in

#include "Dne.fh"

      ispin   = ispin_in
      ne(1)   = ne_in(1)
      ne(2)   = ne_in(2)

      neq(1) = 0
      neq(2) = 0
      call D1dB_init(1,ne(1),map_in)
      call D1dB_nq(1,neq(1))
      if (ispin.eq.2) then
        call D1dB_init(2,ne(2),map_in)
        call D1dB_nq(2,neq(2))
      end if

      return
      end


*     ***********************************
*     *					*
*     *	          Dne_end   		*	
*     *					*
*     ***********************************

      subroutine Dne_end()
      implicit none

#include "Dne.fh"

      call D1dB_end(1)
      if (ispin.eq.2) call D1dB_end(2)

      return
      end
         

*     ***********************************
*     *					*
*     *	          Dneall_ntoqp 		*	
*     *					*
*     ***********************************

      subroutine Dneall_ntoqp(n,q,p)      
      implicit none
      integer n,q,p

#include "Dne.fh"

      if (n.le.ne(1)) then
        call D1dB_ktoqp(1,n,q,p)
      else
        call D1dB_ktoqp(2,n,q,p)
        q = q + neq(1)
      end if

      return
      end


*     ***********************************
*     *					*
*     *	        Dneall_neq		*	
*     *					*
*     ***********************************

      subroutine Dneall_neq(nqtmp)
      implicit none
      integer nqtmp(2)

#include "Dne.fh"

      nqtmp(1) = neq(1)
      nqtmp(2) = neq(2)
      return 
      end


*     ***********************************
*     *					*
*     *	        Dneall_ne 		*	
*     *					*
*     ***********************************

      subroutine Dneall_ne(netmp)
      implicit none
      integer netmp(2)

#include "Dne.fh"

      netmp(1) = ne(1)
      netmp(2) = ne(2)
      return 
      end

c
c
c
c     ****************************************
c     *                                      *
c     *        Dneall_m_allocate             *
c     *                                      *
c     ****************************************
      logical function Dneall_m_allocate(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "mafdecls.fh"

      integer size
     
      if (mb.eq.0) then 
         size = ne(1)*ne(1) + ne(2)*ne(2)
      else
         size = ne(mb)*ne(mb)
      end if
     
       Dneall_m_allocate = MA_alloc_get(mt_dbl,size,'hml',hml(2),hml(1))
      return
      end
c
c     ****************************************
c     *                                      *
c     *        Dneall_m_free                 *
c     *                                      *
c     ****************************************
      logical function Dneall_m_free(hml)
      implicit none
      integer hml(2)

#include "mafdecls.fh"

      Dneall_m_free = MA_free_heap(hml(2))
      return
      end

c
c     ****************************************
c     *                                      *
c     *        Dneall_m_push_get             *
c     *                                      *
c     ****************************************
      logical function Dneall_m_push_get(mb,hml)
      implicit none
      integer mb
      integer hml(2)

#include "Dne.fh"
#include "mafdecls.fh"

      integer size
      if (mb.eq.0) then
         size = ne(1)*ne(1) + ne(2)*ne(2)
      else
         size = ne(mb)*ne(mb)
      end if

      Dneall_m_push_get = MA_push_get(mt_dbl,size,'hml',hml(2),hml(1))
      return
      end

c     ****************************************
c     *                                      *
c     *        Dneall_m_pop_stack            *
c     *                                      *
c     ****************************************
      logical function Dneall_m_pop_stack(hml)
      implicit none
      integer hml(2)

#include "mafdecls.fh"

      Dneall_m_pop_stack = MA_pop_stack(hml(2))
      return
      end



c     ****************************************
c     *                                      *
c     *        DMatrix_fmf_Multiply          *
c     *                                      *
c     ****************************************

      subroutine DMatrix_fmf_Multiply(mb,Ain,npack1,hml,alpha,Aout,beta)
      implicit none
      integer    mb
      complex*16 Ain(*)
      integer    npack1
      real*8     hml(*)
      real*8     alpha
      complex*16 Aout(*)
      real*8     beta
           
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2

c      call Pack_npack(1,npack1)

      if (mb.eq.0) then
         ms1 = 1    
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
         n     = ne(ms)
         if (n.le.0) go to 30
         shift  = 1 + (ms-1)*ne(1)*npack1
         shift2 = 1 + (ms-1)*ishift2
         call DGEMM('N','N',2*npack1,n,n,
     >             (alpha),
     >             Ain(shift),  2*npack1,
     >             hml(shift2),    n,
     >             (beta),
     >             Aout(shift),2*npack1)
   30    continue
      end do

      return
      end


c     ****************************************
c     *                                      *
c     *        DMatrix_ffm_sym_Multiply      *
c     *                                      *
c     ****************************************

      subroutine DMatrix_ffm_sym_Multiply(mb,A1,A2,npack1,hml)
      implicit none
      integer    mb
      complex*16 A1(*),A2(*)
      integer    npack1
      real*8     hml(*)
           
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,shift2,ishift2

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
         shift  = 1 + (ms-1)*ne(1)*npack1
         shift2 = 1 + (ms-1)*ishift2
         n     = ne(ms)
         if (n.le.0) go to 30

         call Pack_ccm_sym_dot(1,n,
     >                     A1(shift),
     >                     A2(shift),
     >                     hml(shift2))
  30     continue
      end do

      return
      end


c     ****************************************
c     *                                      *
c     *        DMatrix_ffm_Multiply          *
c     *                                      *
c     ****************************************

      subroutine DMatrix_ffm_Multiply(mb,A1,A2,npack1,hml)
      implicit none
      integer    mb
      complex*16 A1(*),A2(*)
      integer    npack1
      real*8     hml(*)
           
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift2,ishift2,shift

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
         shift  = 1+(ms-1)*ne(1)*npack1
         shift2 = 1+(ms-1)*ishift2
         n     = ne(ms)
         if (n.le.0) go to 30

         call Pack_ccm_dot(1,n,
     >                     A1(shift),
     >                     A2(shift),
     >                     hml(shift2))
  30     continue
      end do

      return
      end



c     ****************************************
c     *                                      *
c     *        DMatrix_m_diagonalize         *
c     *                                      *
c     ****************************************

      subroutine DMatrix_m_diagonalize(mb,hml,eig,assending)
      implicit none
      integer    mb
      real*8     hml(*),eig(*)
      logical    assending

#include "mafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     ***** local variables ****
      logical value
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,ierr,i,j,indx
      integer tmp1(2)

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift1 = ne(1)
         ishift2 = ne(1)*ne(1)
         call dcopy((ne(1)+ne(2)),0.0d0,0,eig,1)
      else
         ms1 = mb
         ms2 = mb
         ishift1 = 0
         ishift2 = 0
         call dcopy(ne(mb),0.0d0,0,eig,1)
      end if

      value = MA_push_get(mt_dbl,(2*ne(1)*ne(1)),'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('DMatrix_m_diagonalize:out of stack',0,MA_ERR)

*     ***** diagonalize the matrix *****
      do ms=ms1,ms2
         shift1 = 1+(ms-1)*ishift1
         shift2 = 1+(ms-1)*ishift2
         if (ne(ms).le.0) go to 30

         call DSYEV('V','U',ne(ms),
     >              hml(shift2),ne(ms), 
     >              eig(shift1),
     >              dbl_mb(tmp1(1)),2*ne(1)*ne(1),
     >              ierr)

         if (.not.assending)
     >    call eigsrt(eig(shift1),
     >                hml(shift2),
     >                ne(ms),ne(ms))

  30    continue
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('error popping stack in DMatrix_m_diagonalize',0,
     &       MA_ERR)

      return
      end




c     ****************************************
c     *                                      *
c     *        DMatrix_f_SVD                 *
c     *                                      *
c     ****************************************

      subroutine DMatrix_f_SVD(mb,A,U,npack1,S,V)
      implicit none
      integer    mb
      complex*16 A(*),U(*)
      integer    npack1      
      real*8     S(*),V(*)

#include "mafdecls.fh"           
#include "errquit.fh"
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift,tmp2(2)
      integer nj   

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         nj  = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         nj  = ne(mb)
      end if

*     **** generate V and Sigma^2 ****
      call Dmatrix_ffm_sym_Multiply(mb,A,A,npack1,V)
      call Dmatrix_m_diagonalize(mb,V,S,.false.)

*     **** generate U*Sigma ****
      call Dmatrix_fmf_Multiply(mb,A,npack1,V,1.0d0,U,0.0d0)

      
*     **** normalize U*sigma ****
      if (.not. MA_push_get(mt_dbl,neq(1),'tmp2',tmp2(2),tmp2(1)))
     >   call errquit('DMatrix_f_SVD:out of stack memory',0,MA_ERR)

      do ms=ms1,ms2
         if (neq(ms).le.0) go to 30
         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call Pack_cc_idot(1,U(shift),U(shift),dbl_mb(tmp2(1)+n-1))
            shift = shift + npack1
         end do
         call D3dB_Vector_SumAll(neq(ms),dbl_mb(tmp2(1)))

         do n=1,neq(ms)
            dbl_mb(tmp2(1)+n-1) = 1.0d0/dsqrt(dbl_mb(tmp2(1)+n-1))
         end do

         shift = 1+(ms-1)*neq(1)*npack1
         do n=1,neq(ms)
            call dscal(2*npack1,dbl_mb(tmp2(1)+n-1),U(shift),1)
            shift = shift + npack1
         end do

 30     continue
      end do

      if (.not.MA_pop_stack(tmp2(2))) 
     >  call errquit('DMatrix_f_SVD: popping stack memory',0,MA_ERR)


*     **** calculated sqrt(S^2) ****
      do n=1,nj
         if (S(n).lt.0.0d0) S(n) = dabs(S(n))
         S(n) = dsqrt(S(n))
      end do
      return
      end

c     ****************************************
c     *                                      *
c     *        DMatrix_mm_transpose          *
c     *                                      *
c     ****************************************

      subroutine DMatrix_mm_transpose(mb,Min,Mout)
      implicit none
      integer mb
      real*8  Min(*),Mout(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift2,ishift2

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
        if (ne(ms).le.0) goto 30
        shift2 = 1 + (ms-1)*ishift2
        call DMatrix_mm_transpose_sub(ne(ms),Min(shift2),Mout(shift2))
 30     continue
      end do

      return
      end

      subroutine DMatrix_mm_transpose_sub(n,A,B)
      implicit none
      integer n
      real*8 A(n,n)
      real*8 B(n,n)
      integer i,j
      do j=1,n
         do i=1,n
            B(i,j) = A(j,i)
         end do
      end do
      return
      end


c     ****************************************
c     *                                      *
c     *        DMatrix_mm_SCtimesVtrans       *
c     *                                      *
c     ****************************************

      subroutine DMatrix_mm_SCtimesVtrans(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift1 = ne(1)
         ishift2 = ne(1)*ne(1)
         nj = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         ishift1 = 0
         ishift2 = 0
         nj = ne(mb)
      end if

      do j=1,nj
         SA(j) = dcos(S(j)*t)
         SB(j) = dsin(S(j)*t)
      end do

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift1 = 1 + (ms-1)*ishift1
        shift2 = 1 + (ms-1)*ishift2

        indx2 = shift2
        do k=1,ne(ms)
           indx1 = shift1
           do j=1,ne(ms)
              A(indx2) = SA(indx1)*Vt(indx2)
              B(indx2) = SB(indx1)*Vt(indx2)
              indx2 = indx2 + 1
              indx1 = indx1 + 1
           end do
        end do

 30     continue
      end do


      return
      end


c     ****************************************
c     *                                      *
c     *        DMatrix_mm_SCtimesVtrans2     *
c     *                                      *
c     ****************************************

      subroutine DMatrix_mm_SCtimesVtrans2(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj
 
      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift1 = ne(1)
         ishift2 = ne(1)*ne(1)
         nj = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         ishift1 = 0
         ishift2 = 0
         nj = ne(mb)
      end if

      do j=1,nj
         SA(j) = S(j)*dsin(S(j)*t)
         SB(j) = S(j)*dcos(S(j)*t)
      end do

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30

        shift1 = 1 + (ms-1)*ishift1
        shift2 = 1 + (ms-1)*ishift2

        indx2 = shift2
        do k=1,ne(ms)
           indx1 = shift1
           do j=1,ne(ms)
              A(indx2) = SA(indx1)*Vt(indx2)
              B(indx2) = SB(indx1)*Vt(indx2)
              indx2 = indx2 + 1
              indx1 = indx1 + 1
           end do
        end do

 30     continue
      end do

      return
      end



c     ****************************************
c     *                                      *
c     *        DMatrix_mm_SCtimesVtrans3     *
c     *                                      *
c     ****************************************

      subroutine DMatrix_mm_SCtimesVtrans3(mb,t,S,Vt,A,B,SA,SB)
      implicit none
      integer mb
      real*8 t
      real*8 S(*)
      real*8 Vt(*)
      real*8 A(*),B(*)
      real*8 SA(*),SB(*)

#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,shift1,shift2,ishift1,ishift2,k,j,indx1,indx2
      integer nj

      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift1 = ne(1)
         ishift2 = ne(1)*ne(1)
         nj = ne(1)+ne(2)
      else
         ms1 = mb
         ms2 = mb
         ishift1 = 0
         ishift2 = 0
         nj = ne(mb)
      end if

      do j=1,nj
         SA(j) = dsin(S(j)*t)
         SB(j) = 1.0d0-dcos(S(j)*t)
      end do

      do ms=ms1,ms2
        if (ne(ms).le.0) go to 30
        shift1 = 1 + (ms-1)*ishift1
        shift2 = 1 + (ms-1)*ishift2

        indx2 = shift2
        do k=1,ne(ms)
           indx1 = shift1
           do j=1,ne(ms)
              A(indx2) = SA(indx1)*Vt(indx2)
              B(indx2) = SB(indx1)*Vt(indx2)
              indx2 = indx2 + 1
              indx1 = indx1 + 1
           end do
        end do

 30     continue
      end do

      return
      end


c     ****************************************
c     *                                      *
c     *        DMatrix_mmm_Multiply2          *
c     *                                      *
c     ****************************************

      subroutine DMatrix_mmm_Multiply2(mb,A,B,C)
      implicit none
      integer mb
      real*8 A(*),B(*),C(*)
           
#include "Dne.fh"

*     **** local variables ****
      integer ms,ms1,ms2,n,shift2,ishift2
  
      if (mb.eq.0) then
         ms1 = 1
         ms2 = ispin
         ishift2 = ne(1)*ne(1)
      else
         ms1 = mb
         ms2 = mb
         ishift2 = 0
      end if

      do ms=ms1,ms2
         n     = ne(ms)
         if (n.le.0) go to 30
         shift2 = 1 + (ms-1)*ishift2
         call DGEMM('T','N',n,n,n,1.0d0,
     >             A(shift2), n,
     >             B(shift2), n,
     >             0.0d0,
     >             C(shift2), n)
   30    continue
      end do

      return
      end
