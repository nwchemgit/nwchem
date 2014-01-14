
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
      integer MASTER,taskid
      parameter (MASTER=0)

      logical oprint
      integer l

*     **** external functions ****
      logical  control_fmm,control_print
      integer  control_fmm_lmax,control_fmm_lr,ion_nion,control_version
      external control_fmm,control_print
      external control_fmm_lmax,control_fmm_lr,ion_nion,control_version

      fmm      = control_fmm().and.(control_version().eq.4)
      if (fmm) then
         fmm_lmax = control_fmm_lmax()
         fmm_lr   = control_fmm_lr()
         fmm_levels = idint(dlog(dble(ion_nion()))/dlog(8.0d0)) - fmm_lr
         fmm_nboxes = 8**fmm_levels
         fmm_tboxes = 0
         do l=0,fmm_levels
            fmm_tboxes = fmm_tboxes + 8**l
         end  do

         oprint = (taskid.eq.MASTER).and.control_print(print_medium)

         if (oprint) then
            write(luout,100) fmm_lmax,fmm_lr,
     >                       fmm_levels,
     >                       fmm_nboxes,fmm_tboxes
         end if
      end if
      return
  100 format(/"FMM - Fast Multipole Algorithm",
     >    /"   - fmm_lmax   =",i3," fmm_lr = ",i3,
     >    /"   - number of refinement levels = ",i3,
     >    /"   - nboxes (finest level) = ",i8,      
     >    /"   - nboxes (total)        = ",i8)      
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



c*     **********************************************
c*     *                                            *
c*     *              FMM_rion_Llm                  *
c*     *                                            *
c*     **********************************************
c
c      subroutine FMM_rion_Llm(lmax,nion,katm,rion,qatm,rsphere,Llm)
c      implicit none
c      integer lmax
c      integer nion,katm(*)
c      real*8 rion(3,*),qatm(*),rsphere,Llm(*)
c
c*     **** local variables ****
c      integer taskid,np
c      integer ii,l,m,idx
c      real*8 fourpi,r
c
c*     **** external functions ****
c      real*8   Tesseral3_lm_over_rl
c      external Tesseral3_lm_over_rl
c
c      call Parallel_np(np)
c      call Parallel_taskid(taskid)
c
c      fourpi = 16.0d0*datan(1.0d0)
c      call dcopy((lmax+1)**2,0.0d0,0,Llm,1)
c      do ii=1,nion
c         if (mod(ii,np).eq.taskid) then
c            r = dsqrt(rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2)
c            if (r.ge.rsphere) then
c               idx = 1
c               do l=0,lmax
c               do m=-l,l
c                  Llm(idx) = Llm(idx) 
c     >                     + (fourpi/dble(2*l+1))
c     >                      *qatm(katm(ii))
c     >                      *Tesseral3_lm_over_rl(l,m,rion(1,ii),
c     >                                                rion(2,ii),
c     >                                                rion(3,ii))
c                  idx = idx +1
c               end do
c               end do
c            end if
c         end if
c      end do
c      call Parallel_Vector_SumAll((lmax+1)**2,Llm)
c
c      return
c      end
c
c*     **********************************************
c*     *                                            *
c*     *              FMM_Llm_vlocal                *
c*     *                                            *
c*     **********************************************
c
c      subroutine FMM_Llm_vlocal(lmax,Llm,rgrid,vl)
c      implicit none
c      integer lmax
c      real*8  Llm(*)
c      integer ngrid
c      real*8  rgrid(3,*),vl(*)
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
c            vl(i) = vl(i) + Llm(idx)*Tesseral3_lm_rl(l,m,x,y,z)
c            idx = idx + 1
c         end do
c         end do
c      end do
c      return
c      end 

