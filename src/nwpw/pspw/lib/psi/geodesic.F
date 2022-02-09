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

#include "bafdecls.fh"
#include "errquit.fh"

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
*     *		geodesic_get_U_ptr      *
*     *					*
*     ***********************************
      integer function geodesic_get_U_ptr()
      implicit none

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

      geodesic_get_U_ptr = U(1)
      return
      end


*     ***********************************
*     *					*
*     *	       geodesic_ABC_init        *
*     *					*
*     ***********************************
      subroutine geodesic_ABC_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** common block ****
      integer ABC(2)
      real*8  Ermdr
      common / geodesic_ABC_block / ABC,Ermdr

*     **** external functions ****
      logical  Dneall_m_allocate_block
      external Dneall_m_allocate_block

      if (.not.Dneall_m_allocate_block(0,4,ABC))
     >   call errquit('geodesic_ABC_init:alloc ABC',0,MA_ERR)
      return
      end

*     ***********************************
*     *					*
*     *	       geodesic_ABC_finalize    *
*     *					*
*     ***********************************
      subroutine geodesic_ABC_finalize()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / geodesic_ABC_block / ABC,Ermdr

*     **** external functions ****
      logical  Dneall_m_free
      external Dneall_m_free

      if (.not.Dneall_m_free(ABC))
     >   call errquit('geodesic_ABC_end:deallocate ABC',0,MA_ERR)
      return
      end



*     ***********************************
*     *					*
*     *	       geodesic_ABC_start       *
*     *					*
*     ***********************************
      subroutine geodesic_ABC_start(ispin,ne,npack1,Ermdr0,Y,Z,HY,HZ)
      implicit none
      integer    ispin,ne(2),npack1
      real*8     Ermdr0
      complex*16 Y(*),Z(*),HY(*),HZ(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / geodesic_ABC_block / ABC,Ermdr

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      integer ms,ishift,A,Bdag,C,tmp1(2)

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      integer  psi_ne
      external psi_ne
      real*8   Dneall_m_trace
      external Dneall_m_trace

      Ermdr = Ermdr0

      call Dneall_m_size(0,ishift)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      if (.not.Dneall_m_push_get(0,tmp1))
     >   call errquit('geodesic_ABC_start:stack',0,MA_ERR)

*     **** Copy Hpsi_k to HY then ****
*     ****   generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
      call electron_get_Hpsi_k(HY)
      call electron_run0(Z)
      call electron_get_Hpsi_k(HZ)

*     **** generate A=<Y|HY>, Bdag=<U|HY>, C=<U|HU> matrices ****
      call Dneall_ffm_sym_Multiply(0,Y,HY,npack1,dbl_mb(A))
      call Dneall_ffm_Multiply(0,Z,HY,npack1,dbl_mb(Bdag))
      call Dneall_ffm_sym_Multiply(0,Z,HZ,npack1,dbl_mb(C))

*     **** compute V'*A*V ****
      call Dneall_mmm_Multiply3(0,dbl_mb(A),
     >                            dbl_mb(Vt(1)),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(Vt(1)),
     >                           dbl_mb(tmp1(1)),
     >                           1.0d0,
     >                           dbl_mb(A),0.0d0)

*     **** compute B'*V ****
      call Dneall_mmm_Multiply3(0,dbl_mb(Bdag),
     >                           dbl_mb(Vt(1)),
     >                           dbl_mb(tmp1(1)))
      call Parallel_shared_vector_copy(.true.,ishift,
     >                                 dbl_mb(tmp1(1)),dbl_mb(Bdag))

      if (.not.Dneall_m_pop_stack(tmp1))
     >   call errquit('geodesic_ABC_start:pop',1,MA_ERR)

      return
      end

      subroutine junk_mat_print(n,A)
      implicit none
      integer n
      real*8 A(n,n)
      integer i,j

      write(*,*) "mat_print:"
      do i=1,n
         write(*,*) (A(i,j),j=1,n)
      end do
      
      return
      end

*     ***********************************
*     *					*
*     *	       geodesic_ABC_energy      *
*     *					*
*     ***********************************
      real*8 function geodesic_ABC_energy(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / geodesic_ABC_block / ABC,Ermdr

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      integer ishift,A,Bdag,C,nemax
      integer tmp1(2),tmp2(2),tmp3(2),tmpC(2),tmpS(2),tmp4(2)
      real*8  energy

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne
      real*8   Dneall_m_trace
      external Dneall_m_trace


      call Dneall_m_size(0,ishift)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      nemax = psi_ne(1)+psi_ne(2)

*     **** push stack space ****
      value =           Dneall_m_push_get(0,tmp1)
      value = value.and.Dneall_m_push_get(0,tmp2)
      value = value.and.Dneall_m_push_get(0,tmp3)
      value = value.and.Dneall_m_push_get(0,tmp4)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) call errquit('geodesic_ABC_energy:stack',0,MA_ERR)

      call Dneall_mm_SCtimesVtrans(0,t,dbl_mb(S(1)),
     >                                dbl_mb(Vt(1)),
     >                                dbl_mb(tmp1(1)),
     >                                dbl_mb(tmp3(1)),
     >                                dbl_mb(tmpC(1)),
     >                                dbl_mb(tmpS(1)))

*     **** V*cos*B*sin*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(Bdag),
     >                            dbl_mb(tmp3(1)),
     >                            dbl_mb(tmp2(1)))
      call Dneall_mmm_Multiply2(0,dbl_mb(tmp1(1)),
     >                            dbl_mb(tmp2(1)),
     >                            dbl_mb(tmp4(1)))

*     **** V*cos*A*cos*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(tmp1(1)),
     >                            dbl_mb(A),
     >                            dbl_mb(tmp2(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp2(1)),
     >                           dbl_mb(tmp1(1)),1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)

*     **** V*sin*B'*cos*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(tmp3(1)),
     >                            dbl_mb(Bdag),
     >                            dbl_mb(tmp2(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp2(1)),
     >                           dbl_mb(tmp1(1)),1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)

*     **** V*sin*C*sin*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(tmp3(1)),
     >                            dbl_mb(C),
     >                            dbl_mb(tmp2(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp2(1)),
     >                           dbl_mb(tmp3(1)),1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)

      !tmp4 = tmp1'*A*tmp1
      !tmp4 += tmp1'*B*tmp3
      !tmp4 += tmp3'*B'*tmp1
      !tmp4 += tmp3'*C*tmp3

      energy = Dneall_m_trace(0,dbl_mb(tmp4(1)))
      if (psi_ispin().eq.1) energy = energy+energy

*     **** pop stack space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.Dneall_m_pop_stack(tmp4)
      value = value.and.Dneall_m_pop_stack(tmp3)
      value = value.and.Dneall_m_pop_stack(tmp2)
      value = value.and.Dneall_m_pop_stack(tmp1)
      if (.not.value) call errquit('geodesic_ABC_energy:pop',1,MA_ERR)

      geodesic_ABC_energy = energy + Ermdr
      return
      end

*     ***********************************
*     *					*
*     *	      geodesic_ABC_denergy      *
*     *					*
*     ***********************************
      real*8 function geodesic_ABC_denergy(t)
      implicit none
      real*8 t

#include "bafdecls.fh"
#include "errquit.fh"

*     **** ABC common block ****
      integer ABC(2)
      real*8  Ermdr
      common / geodesic_ABC_block / ABC,Ermdr

*     **** geodesic common block ***
      integer U(2)
      integer Vt(2)
      integer S(2)
      common / geodesic_block / U,Vt,S

*     **** local variables ****
      logical value
      integer nemax,ishift,A,Bdag,C
      integer CV(2),SV(2),SSV(2),CSV(2),tmpC(2),tmpS(2)
      integer tmp1(2),tmp4(2)
      real*8  denergy,denergyB,denergyA,denergyC,denergyBdag

*     **** external functions ****
      logical  Dneall_m_push_get,Dneall_m_pop_stack
      external Dneall_m_push_get,Dneall_m_pop_stack
      integer  psi_ispin,psi_ne
      external psi_ispin,psi_ne
      real*8   Dneall_m_trace
      external Dneall_m_trace


      call Dneall_m_size(0,ishift)
      A    = ABC(1)
      Bdag = ABC(1) +   ishift
      C    = ABC(1) + 2*ishift

      nemax = psi_ne(1)+psi_ne(2)

*     **** push stack space ****
      value =           Dneall_m_push_get(0,CV)
      value = value.and.Dneall_m_push_get(0,SV)
      value = value.and.Dneall_m_push_get(0,SSV)
      value = value.and.Dneall_m_push_get(0,CSV)
      value = value.and.Dneall_m_push_get(0,tmp1)
      value = value.and.Dneall_m_push_get(0,tmp4)
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpC',tmpC(2),tmpC(1))
      value = value.and.BA_push_get(mt_dbl,nemax,'tmpS',tmpS(2),tmpS(1))
      if (.not.value) call errquit('geodesic_ABC_denergy:push',0,MA_ERR)

      call Dneall_mm_SCtimesVtrans(0,t,dbl_mb(S(1)),
     >                               dbl_mb(Vt(1)),
     >                               dbl_mb(CV(1)),
     >                               dbl_mb(SV(1)),
     >                               dbl_mb(tmpC(1)),
     >                               dbl_mb(tmpS(1)))
      call Dneall_mm_SCtimesVtrans2(0,t,dbl_mb(S(1)),
     >                                dbl_mb(Vt(1)),
     >                                dbl_mb(SSV(1)),
     >                                dbl_mb(CSV(1)),
     >                                dbl_mb(tmpC(1)),
     >                                dbl_mb(tmpS(1)))


*     **** -V*S*sin*B*sin*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(Bdag),
     >                            dbl_mb(SV(1)),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply2(0,dbl_mb(SSV(1)),
     >                            dbl_mb(tmp1(1)),
     >                            dbl_mb(tmp4(1)))

      denergy = -Dneall_m_trace(0,dbl_mb(tmp4(1)))

*     ****  V*cos*B*cos*S*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(Bdag),
     >                            dbl_mb(CSV(1)),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply2(0,dbl_mb(CV(1)),
     >                            dbl_mb(tmp1(1)),
     >                            dbl_mb(tmp4(1)))

      denergy = denergy + Dneall_m_trace(0,dbl_mb(tmp4(1)))


*     **** -V*S*sin*A*cos*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(SSV(1)),
     >                            dbl_mb(A),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(CV(1)), -1.0d0,
     >                           dbl_mb(tmp4(1)),0.0d0)
*     **** -V*cos*A*sin*S*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(CV(1)),
     >                            dbl_mb(A),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(SSV(1)),-1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)

*     **** V*S*cos* B' *cos*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(CSV(1)),
     >                            dbl_mb(Bdag),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(CV(1)),  1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)
*     **** -V*sin*B'*sin*S*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(SV(1)),
     >                            dbl_mb(Bdag),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(SSV(1)), -1.0d0,
     >                           dbl_mb(tmp4(1)), 1.0d0)

*     ****  V*S*cos*C*sin*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(CSV(1)),
     >                            dbl_mb(C),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(SV(1)),  1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)
*     ****  V*sin*C*cos*S*V' ****
      call Dneall_mmm_Multiply2(0,dbl_mb(SV(1)),
     >                            dbl_mb(C),
     >                            dbl_mb(tmp1(1)))
      call Dneall_mmm_Multiply(0,dbl_mb(tmp1(1)),
     >                           dbl_mb(CSV(1)), 1.0d0,
     >                           dbl_mb(tmp4(1)),1.0d0)

      denergy = denergy + Dneall_m_trace(0,dbl_mb(tmp4(1)))
      if (psi_ispin().eq.1) denergy = denergy+denergy

*     **** pop stack space ****
      value =           BA_pop_stack(tmpS(2))
      value = value.and.BA_pop_stack(tmpC(2))
      value = value.and.BA_pop_stack(tmp4(2))
      value = value.and.BA_pop_stack(tmp1(2))
      value = value.and.Dneall_m_pop_stack(CSV)
      value = value.and.Dneall_m_pop_stack(SSV)
      value = value.and.Dneall_m_pop_stack(SV)
      value = value.and.Dneall_m_pop_stack(CV)
      if (.not.value) call errquit('geodesic_ABC_denergy:pop',1,MA_ERR)

      geodesic_ABC_denergy = denergy
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



