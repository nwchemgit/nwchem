c
c $Id$
c
 
*     ***********************************
*     *					*
*     *        Orb_Analysis		*
*     *					*
*     ***********************************

      subroutine Orb_Analysis(iunit,flag,
     >                        ispin,ne,npack1,nemax,c2,
     >                        lmax,b0,a,sum)
      
      implicit none
      integer    iunit,flag
      integer    ispin,ne(2)
      integer    npack1,nemax
      complex*16 c2(npack1,*)

      integer lmax(*)
      real*8  B0(0:5,*)
      real*8  a(36,nemax,*)
      real*8  sum(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** parallel variables ****
      integer  taskid
      integer  MASTER
      parameter(MASTER=0)

*     **** local variables ****
      logical value
      integer spin
      integer i,k,l,n,ll,l1,l2
      integer llmax,nfft3d
      real*8  omega,w,pi,twopi
      real*8  f0,f1,f21,f22,f23,f24,f25
      integer n1(2),n2(2)



*     **** external functions ****
      logical  aorbs_init,aorbs_readall,control_mulliken_kawai
      external aorbs_init,aorbs_readall,control_mulliken_kawai
      integer  ion_nion_qm,ion_katm_qm
      external ion_nion_qm,ion_katm_qm
      real*8   ion_rion
      external ion_rion

      call Parallel_taskid(taskid)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)


*     **** calculated POPULATION ANALYSIS ****

      if (control_mulliken_kawai()) then
        call Orb_pop_kawai(flag,ispin,ne,npack1,nemax,c2,lmax,b0,a,sum)
      else
        call Orb_pop_aorb(flag,ispin,ne,npack1,nemax,c2,lmax,a,sum)
      end if

*     **** write out POPULATION ANALYSIS ****
      if (taskid.eq.MASTER) then
         write(iunit,*) n2(ispin),ion_nion_qm()
         write(iunit,911) (ion_rion(1,i),ion_rion(2,i),ion_rion(3,i),
     >                     (lmax(ion_katm_qm(i))+1)**2,
     >                     i=1,ion_nion_qm())
 911     format(f10.4, f10.4, f10.4, i3)

         do spin=1,ispin
         do n=n1(spin),n2(spin)

            do i=1,ion_nion_qm()
            do l=0,lmax(ion_katm_qm(i))
               l1 =l**2+1
               l2 = (l+1)**2
               write(iunit,1531) (A(ll,n,i),ll=l1,l2)
            end do
            end do
         end do
         end do
         write(iunit,*)
 1531    format(8f11.5)
         call util_flush(iunit)
      endif

      return 
      end
 
