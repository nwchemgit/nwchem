
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

      subroutine FMM_rion_Llm(lmax,nion,katm,rion,qatm,rsphere,Llm)
      implicit none
      integer lmax
      integer nion,katm(*)
      real*8 rion(3,*),qatm(*),rsphere,Llm(*)

*     **** local variables ****
      integer taskid,np
      integer ii,l,m,idx
      real*8 fourpi,r

*     **** external functions ****
      real*8   Tesseral3_lm_over_rl
      external Tesseral3_lm_over_rl

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      fourpi = 16.0d0*datan(1.0d0)
      call dcopy((lmax+1)**2,0.0d0,0,Llm,1)
      do ii=1,nion
         if (mod(ii,np).eq.taskid) then
            r = dsqrt(rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2)
            if (r.ge.rsphere) then
               idx = 1
               do l=0,lmax
               do m=-l,l
                  Llm(idx) = Llm(idx) 
     >                     + (fourpi/dble(2*l+1))
     >                      *qatm(katm(ii))
     >                      *Tesseral3_lm_over_rl(l,m,rion(1,ii),
     >                                                rion(2,ii),
     >                                                rion(3,ii))
                  idx = idx +1
               end do
               end do
            end if
         end if
      end do
      call Parallel_Vector_SumAll((lmax+1)**2,Llm)

      return
      end

*     **********************************************
*     *                                            *
*     *              FMM_Llm_vlocal                *
*     *                                            *
*     **********************************************

      subroutine FMM_Llm_vlocal(lmax,Llm,rgrid,vl)
      implicit none
      integer lmax
      real*8  Llm(*)
      integer ngrid
      real*8  rgrid(3,*),vl(*)

*     **** local variables ****
      integer i,l,m,idx
      real*8 x,y,z

*     **** external functions ****
      real*8   Tesseral3_lm_rl
      external Tesseral3_lm_rl

      do i=1,ngrid
         x = rgrid(1,i)
         y = rgrid(2,i)
         z = rgrid(3,i)
         vl(i) = 0.0d0
         idx = 1
         do l=0,lmax
         do m=-l,l
            vl(i) = vl(i) + Llm(idx)*Tesseral3_lm_rl(l,m,x,y,z)
            idx = idx + 1
         end do
         end do
      end do
      return
      end 



