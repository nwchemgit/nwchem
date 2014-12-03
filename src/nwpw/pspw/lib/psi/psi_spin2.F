c $Id$


*     ******************************
*     *                            *
*     *    Calculate_psi_spin2     *
*     *                            *
*     ******************************

      subroutine Calculate_psi_spin2(ispin,ne,nfft3d,psi,
     >                               fractional,occ,Sab)
      implicit none
      integer ispin,ne(2)
      integer nfft3d
      complex*16 psi(nfft3d,*)
      logical fractional
      real*8 occ(*)
      real*8 Sab

#include "util.fh"
#include "bafdecls.fh"
#include "errquit.fh"
#include "stdio.fh"
ccc#include "frac_occ.fh"

*     **** local variables ****
      logical oprint
      integer i,j,ms,sm(2)
      real*8  Sexact,SS,en(2)
       
      integer taskid,MASTER
      parameter (MASTER=0)

*     **** external functions ****
      logical  control_print,Dneall_sm_push_get,Dneall_sm_pop_stack
      external control_print,Dneall_sm_push_get,Dneall_sm_pop_stack
      real*8   Dneall_sm_value
      external Dneall_sm_value


      call Parallel_taskid(taskid)
      oprint = (taskid.eq.MASTER).and.control_print(print_medium)

*     ****************************
*     ***** calculate <S^2> ******
*     ****************************

      if (ispin.eq.2) then

         en(1)     = ne(1)
         en(ispin) = ne(ispin)
         if (fractional) then
            do ms=1,ispin
            en(ms) =0.0
            do i=1,ne(ms)
              en(ms) = en(ms) + occ(i+(ms-1)*ne(1))
            end do
            end do
         end if

         Sexact = ((en(1)-en(2))/2.0d0)
     >          * (((en(1)-en(2))/2.0d0)+1.0d0)
         Sab = 0.0d0
         if (ne(2).gt.0) then
            if (.not.Dneall_sm_push_get(sm)) then
              call errquit('Calculate_psi_spin2:out of stack',0,MA_ERR)
            end if
            call Dneall_ffsm_Multiply(psi,psi,nfft3d,dbl_mb(sm(1)))
            do i=1,ne(1)
            do j=1,ne(2)
*             *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
               SS = Dneall_sm_value(i,j,dbl_mb(sm(1)))
               Sab = Sab + SS*SS
            end do
            end do

            if (.not.Dneall_sm_pop_stack(sm)) then
              call errquit('Calculate_psi_spin2:popping stack',0,MA_ERR)
            end if
         end if
         Sab = Sexact + en(2) - Sab

c         do i=1,ne(1)
c         do j=1,ne(2)
c*           *** SS =  GCDOTC(nida,nidb,c2(1,i),c2(1,j+ne(1))) ****
c            call Pack_cc_dot(1,psi(1,i),psi(1,j+ne(1)),SS)
c            Sab = Sab + SS*SS
c         end do
c         end do
c         Sab = Sexact + en(2) - Sab
      else
         Sexact = 0.0d0
         Sab    = 0.0d0
      end if

      if (oprint) then
         write(LuOut,'(//A)') '=== Spin Contamination ==='
         write(LuOut,*)
         write(LuOut,*) "<Sexact^2> = ", Sexact
         write(LuOut,*) "<S^2>      = ", Sab
         write(LuOut,*)
      end if

      return
      end
