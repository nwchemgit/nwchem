c      program test
c      implicit none
c
c      integer nion,n,na,nqm,nb,lr
c      real*8  a_in,b_in,aqm,da,a_out,b_out
c
c      write(*,*) "Enter a,b,aqm,nion,lr:"
c      read(*,*) a_in,b_in,aqm,nion,lr
c 
c      call FMM_Generate_Box(a_in,b_in,aqm,nion,lr,
c     >                      n,na,nqm,nb,da,a_out,b_out)
c
c      write(*,*)
c      write(*,*) "da,n,n*da      = ",n,da,da*n
c      write(*,*) "a_in,a_out,na  = ",na,a_in,a_out
c      write(*,*) "b_in,b_out,nb  = ",nb,b_in,b_out
c      write(*,*) "aqm,nqm        = ",nqm,aqm
c      write(*,*) 
c      stop
c      end

*     **********************************************
*     *                                            *
*     *             FMM_start                      *
*     *                                            *
*     **********************************************
*
      subroutine FMM_start()
      implicit none

#include "util.fh"
#include "stdio.fh"
#include "fmm.fh"

*     **** local variables ****
      integer MASTER,taskid,np
      parameter (MASTER=0)

      logical oprint
      integer l,m,i,j,k,lm

*     **** external functions ****
      logical  control_fmm,control_print
      integer  control_fmm_lmax,control_fmm_lr,ion_nion,control_version
      external control_fmm,control_print
      external control_fmm_lmax,control_fmm_lr,ion_nion,control_version

      fmm      = control_fmm().and.(control_version().eq.4)
      if (fmm) then
         lmax = control_fmm_lmax()
         lr   = control_fmm_lr()
         levels = nint(dlog(dble(ion_nion()))/dlog(8.0d0)) - lr
         nboxes = 8**levels
         tboxes = 0
         do l=0,levels
            tboxes = tboxes + 8**l
         end  do

         call Parallel_taskid(taskid)
         call Parallel_np(np)
         oprint = (taskid.eq.MASTER).and.control_print(print_medium)

c         value = MA_alloc_get(mt_dbl,26*(lmax+1)**2,'first_tlm_fmm',
c     >                       first_tlm_fmm(2),first_tlm_fmm(1))
c         value = value.and.
c     >           MA_alloc_get(mt_dbl,98*(lmax+1)**2,'second_tlm_fmm',
c     >                        second_tlm_fmm(2),second_tlm_fmm(1))
c         value = value.and.
c     >           MA_alloc_get(mt_dbl,7216,'interaction_tlm_fmm',
c     >                        interaction_tlm_fmm(2),interaction_tlm_fmm(1))
c
c*        **** first neigbors of the current box      ****
c*        **** first_count = 3*3*3 - 1*1*1 = 27-1 = 26 ****
c         first_count = 0
c         do k=-1,1
c            do j=-1,1
c               do i=-1,1
c                  if ((i.ne.0).or.(j.ne.0).or.(k.eq.0)) then
c                     do l=0,fmm_lmax
c                     do m=-l,l
c                        dbl_mb(first_tlm_fmm(1)+first_count) 
c     >                     = Tesseral3_lm(l,m,dble(i),dble(j),dble(k))
c                        first_count = first_count + 1
c                     end do
c                     end do
c                  end if
c               end do
c            end do
c         end do
c
c*        **** second nearest neigbors of the current box ****
c*        **** second_count = 5*5*5 - 3*3*3 = 125-27 = 98 ****
c         second_count = 0
c         do k=-2,2
c            do j=-2,2
c               do i=-2,2
c                  if ((i*i+j*j+k*k).gt.2) then
c                     do l=0,fmm_lmax
c                     do m=-l,l
c                        dbl_mb(second_tlm_fmm(1)+second_count) 
c     >                     = Tesseral3_lm(l,m,dble(i),dble(j),dble(k))
c                        second_count = second_count + 1
c                     end do
c                     end do
c                  end if
c               end do
c            end do
c         end do

c*        **** interaction neigbors - set of boxes which are children of            ****
c*        **** the nearest and second nearest neighbors of the current boxes        ****
c*        **** parent which are not nearest or second nearest neighbors of          ****
c*        **** the current box                                                      ****
c*        **** interaction_count=8*(8*5*5*5-98)=8*(8*125-98)=8*(1000-98)=8*902=7216 ****

         if (oprint) then
            write(luout,100) lmax,lr,
     >                       levels,
     >                       nboxes,tboxes
         end if
      end if
      return
  100 format(/" FMM - Fast Multipole Algorithm",
     >    /"     - fmm_lmax = ",i3," fmm_lr = ",i3,
     >    /"     - number of refinement levels = ",i2,
     >    /"     - nboxes (finest level) = ",i8,      
     >    /"     - nboxes (total)        = ",i8/)      
      end


*     **********************************************
*     *                                            *
*     *                FMM_end                     *
*     *                                            *
*     **********************************************

      subroutine FMM_end()
      implicit none

#include "fmm.fh"

      if (fmm) then
      end if
      return
      end


      logical function FMM_fmm()
      implicit none

#include "fmm.fh"

      FMM_fmm = fmm
      return
      end

      integer function FMM_lmax()
      implicit none

#include "fmm.fh"
     
      FMM_lmax = lmax
      return
      end


*     **********************************************
*     *                                            *
*     *             FMM_Generate_Box               *
*     *                                            *
*     **********************************************
*
*   This routine generates a box gridding for FMM calculations.
*
*   Entry - a_in - guess for start of box 
*           b_in - guess for end of box 
*           aqm  - simple cubic side length
*           nion - number of ions
*           lr   - level reduction factor
* 
*   Exit - levels - number of tree levels
*          n      - number of boxes from a to b
*          na     - number of boxes from a to (-aqm/2)
*          nqm    - number of boxes from (-aqm/2) to (aqm/2)
*          nb     - number of boxes from (aqm/2) to b
*          da     - length of small boxes
*          a_out  - box start
*          b_out  - box end
*
      subroutine FMM_Generate_Box(a_in,b_in,aqm,nion,lr,
     >                            levels,n,na,nqm,nb,da,a_out,b_out)
      implicit none
      real*8 a_in,b_in,aqm
      integer nion,lr
      integer levels,n,na,nqm,nb
      real*8 da,a_out,b_out

      levels = idint(dlog(1.0d0*nion)/dlog(8.0d0))-lr
      n = 2**levels

      nqm = idint(n*aqm/(b_in-a_in))
      da = aqm/dble(nqm)

      na = idint(n*(-0.5d0*aqm-a_in)/(b_in-a_in))
      nb = n - nqm - na
      a_out = -0.5d0*aqm - na*da
      b_out =  0.5d0*aqm + nb*da
      return
      end



*     **********************************************
*     *                                            *
*     *              FMM_rion_Llm                  *
*     *                                            *
*     **********************************************
*
      subroutine FMM_rion_Llm(lmax,nion,katm,rion,qatm,rsphere2,Llm)
      implicit none
      integer lmax
      integer nion,katm(*)
      real*8 rion(3,*),qatm(*),rsphere2,Llm(*)

*     **** local variables ****
      integer taskid,np
      integer ii,l,m,idx
      real*8 fourpi,r2

*     **** external functions ****
      real*8   Tesseral3_lm_over_rl
      external Tesseral3_lm_over_rl

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      fourpi = 16.0d0*datan(1.0d0)
      call dcopy((lmax+1)**2,0.0d0,0,Llm,1)
      do ii=1,nion
         if (mod(ii,np).eq.taskid) then
            r2 = rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2
            if (r2.gt.rsphere2) then
               idx = 1
               do l=0,lmax
               do m=-l,l
                  Llm(idx) = Llm(idx) 
     >                     + (fourpi/dble(2*l+1))
     >                      *qatm(katm(ii))
     >                      *Tesseral3_lm_over_rl(l,m,rion(1,ii),
     >                                                rion(2,ii),
     >                                                rion(3,ii))
                  idx = idx + 1
               end do
               end do
            end if
         end if
      end do
      call Parallel_Vector_SumAll((lmax+1)**2,Llm)
      return
      end


c*     **********************************************
c*     *                                            *
c*     *              FMM_rgrid_Mlm                 *
c*     *                                            *
c*     **********************************************
c
c      subroutine FMM_rgrid_Mlm(lmax,ngrid,rgrid,Mlm)
c      implicit none
c      integer lmax
c      integer ngrid
c      real*8  rgrid(3,*)
c      real*8  Mlm(*)
c
c*     **** local variables ****
c      integer i,l,m,idx
c      real*8 x,y,z
c
c*     **** external functions ****
c      real*8   Tesseral3_lm_rl
c      external Tesseral3_lm_rl
c
c      do i=1,ngrid
c         x = rgrid(1,i)
c         y = rgrid(2,i)
c         z = rgrid(3,i)
c         vl(i) = 0.0d0
c         idx = 1
c         do l=0,lmax
c         do m=-l,l
c            vl(i) = vl(i) + Llm(idx)*Tesseral3_lm_over_rl(l,m,x,y,z)
c            idx = idx + 1
c         end do
c         end do
c      end do
c      return
c      end 

