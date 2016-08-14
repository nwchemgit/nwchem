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

c         value = BA_alloc_get(mt_dbl,26*(lmax+1)**2,'first_tlm_fmm',
c     >                       first_tlm_fmm(2),first_tlm_fmm(1))
c         value = value.and.
c     >           BA_alloc_get(mt_dbl,98*(lmax+1)**2,'second_tlm_fmm',
c     >                        second_tlm_fmm(2),second_tlm_fmm(1))
c         value = value.and.
c     >           BA_alloc_get(mt_dbl,7216,'interaction_tlm_fmm',
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
*     *             FMM_Setup_Tree                 *
*     *                                            *
*     **********************************************

*   This routine generates a box gridding and tree for FMM calculations.
*
*   Entry - nion - number of ions
*           rion - atom positions
*           levels - number of tree levels
* 
*   Exit - a_out  - box length
*          ion_boxindx - box indexes for the ions

      subroutine FMM_Setup_Tree(nion,rion,levels,a_out,ion_boxindx)
      implicit none
      integer nion
      real*8  rion(3,*)
      integer levels

      real*8 a_out
      integer ion_boxindx(3,*)

*     **** local variables ****
      integer ii,n,n2
      real*8  aover2

      aover2 = 0.0d0
      do ii=1,nion
         if (dabs(rion(1,ii)).gt.aover2) aover2=dabs(rion(1,ii))
         if (dabs(rion(2,ii)).gt.aover2) aover2=dabs(rion(2,ii))
         if (dabs(rion(3,ii)).gt.aover2) aover2=dabs(rion(3,ii))
      end do
      aover2 = aover2 + 0.5d0
      a_out = 2*aover2

      n = 2**levels
      n2 = n*n
      do ii=1,nion
         ion_boxindx(1,ii) =  idint(n*(rion(1,ii)+aover2)/a_out)
         ion_boxindx(2,ii) =  idint(n*(rion(2,ii)+aover2)/a_out)
         ion_boxindx(3,ii) =  idint(n*(rion(3,ii)+aover2)/a_out)
      end do

      return
      end

*     **********************************************
*     *                                            *
*     *             FMM_Setup_Mlm                  *
*     *                                            *
*     **********************************************

      subroutine FMM_Setup_Mlm(nion,katm,rion,qatm,
     >                         levels,a,ion_boxindx,
     >                         lmax,Mlm)
      implicit none
      integer nion,katm(*)
      real*8  rion(3,*)
      real*8  qatm(*)
      integer levels
      real*8  a
      integer ion_boxindx(3,*)
      integer lmax
      complex*16 Mlm((lmax+1)**2,*)

*     **** local variables ****
      integer taskid,np
      integer ii,i,j,k,n,n2,ibox,idx,l,m
      real*8 aover2,da,xc,yc,zc,q

*     **** external functions ****
      complex*16 rSphericalHarmonic3
      external   rSphericalHarmonic3

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      n = 2**levels
      n2 = n*n
      da = a/dble(a)
      aover2 = a/2.0d0
      call dcopy(2*((lmax+1)**2)*(n*n*n),0.0d0,0,Mlm,1)
      do ii=1,nion
         if (mod(ii,np).eq.taskid) then
            q = qatm(katm(ii))
            xc = -aover2 + (ion_boxindx(1,ii)+0.5d0)*da
            yc = -aover2 + (ion_boxindx(2,ii)+0.5d0)*da
            zc = -aover2 + (ion_boxindx(3,ii)+0.5d0)*da
            ibox = 1 + ion_boxindx(1,ii) 
     >           + ion_boxindx(2,ii)*n 
     >           + ion_boxindx(3,ii)*n2
            idx = 1
            do l=0,lmax
            do m=-l,l
               Mlm(idx,ibox) = Mlm(idx,ibox) 
     >                       + q*rSphericalHarmonic3(l,-m,
     >                                               xc-rion(1,ii),
     >                                               yc-rion(2,ii),
     >                                               zc-rion(3,ii))
               idx = idx + 1
            end do
            end do
         end if
      end do
      call Parallel_Vector_SumAll(((lmax+1)**2)*(n*n*n),Mlm)

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
      subroutine FMM_rion_Llm(lmax,nion,nion_qm,katm,rion,qatm,
     >                        lmbda0,rsphere2,Llm)
      implicit none
      integer lmax
      integer nion,nion_qm,katm(*)
      real*8 rion(3,*),qatm(*),lmbda0,rsphere2,Llm(*)

*     **** local variables ****
      integer taskid,np
      integer ii,l,m,idx
      real*8 fourpi,r2,lmbda

*     **** external functions ****
      real*8   Tesseral3_lm_over_rl
      external Tesseral3_lm_over_rl

      call Parallel_np(np)
      call Parallel_taskid(taskid)

      fourpi = 16.0d0*datan(1.0d0)
      !call dcopy((lmax+1)**2,0.0d0,0,Llm,1)
      call Parallel_shared_vector_zero(.true.,(lmax+1)**2,Llm)
      do ii=1,nion
         if (mod(ii,np).eq.taskid) then
            r2 = rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2
            if (r2.gt.rsphere2) then
               if (ii.gt.nion_qm) then
                  lmbda = lmbda0
               else
                  lmbda = 1.0d0
               end if
               !idx = 1
!$OMP DO
               do l=0,lmax
               do m=-l,l
                  idx = l*l+l+m+1
                  Llm(idx) = Llm(idx) 
     >                     + (fourpi/dble(2*l+1))
     >                      *qatm(katm(ii))*lmbda
     >                      *Tesseral3_lm_over_rl(l,m,rion(1,ii),
     >                                                rion(2,ii),
     >                                                rion(3,ii))
                  !idx = idx + 1
               end do
               end do
!$OMP END DO
            end if
         end if
      end do
      call Parallel_Vector_SumAll((lmax+1)**2,Llm)
      return
      end

*     **********************************************
*     *                                            *
*     *              FMM_rion_Llm_mm               *
*     *                                            *
*     **********************************************
*
      subroutine FMM_rion_Llm_mm(lmax,nion,nion_qm,katm,rion,qatm,
     >                           rsphere2,Llm)
      implicit none
      integer lmax
      integer nion,nion_qm,katm(*)
      real*8 rion(3,*),qatm(*),lmbda0,rsphere2,Llm(*)

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
      !call dcopy((lmax+1)**2,0.0d0,0,Llm,1)
      call Parallel_shared_vector_zero(.true.,(lmax+1)**2,Llm)
      do ii=nion_qm+1,nion
         if (mod(ii,np).eq.taskid) then
            r2 = rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2
            if (r2.gt.rsphere2) then

               !idx = 1
!$OMP DO
               do l=0,lmax
               do m=-l,l
                  idx = l*l+l+m+1
                  Llm(idx) = Llm(idx)
     >                     + (fourpi/dble(2*l+1))
     >                      *qatm(katm(ii))
     >                      *Tesseral3_lm_over_rl(l,m,rion(1,ii),
     >                                                rion(2,ii),
     >                                                rion(3,ii))
                  !idx = idx + 1
               end do
               end do
!$OMP END DO
            end if
         end if
      end do
      call Parallel_Vector_SumAll((lmax+1)**2,Llm)
      return
      end




*     **********************************************
*     *                                            *
*     *              FMM_fion_Mlm                  *
*     *                                            *
*     **********************************************
*
      subroutine FMM_fion_Mlm(lmax,nion,nion_qm,katm,rion,qatm,
     >                        lmbda0,rsphere2,
     >                        Mlm,fion)
      implicit none
      integer lmax
      integer nion,nion_qm,katm(*)
      real*8 rion(3,*),qatm(*),lmbda0,rsphere2,Mlm(*),fion(3,nion)

*     **** local variables ****
      integer taskid,np
      integer ii,l,m,idx
      real*8 fourpi,r2,coef,dtx,dty,dtz,lmbda


      call Parallel_np(np)
      call Parallel_taskid(taskid)

      fourpi = 16.0d0*datan(1.0d0)
      do ii=1,nion
         if (mod(ii,np).eq.taskid) then
            r2 = rion(1,ii)**2+rion(2,ii)**2+rion(3,ii)**2
            if (r2.gt.rsphere2) then
               if (ii.gt.nion_qm) then
                  lmbda = lmbda0
               else
                  lmbda = 1.0d0
               end if
c               idx = 1
!$OMP DO REDUCTION(+:fion)
               do l=0,lmax
                  coef = (fourpi/dble(2*l+1))*qatm(katm(ii))
                  do m=-l,l
                     idx = l*l+l+m+1
                     call dTesseral3_lm_over_rl(l,m,rion(1,ii),
     >                                              rion(2,ii),
     >                                              rion(3,ii),
     >                                              dtx,dty,dtz)
                     fion(1,ii) = fion(1,ii) + dtx*coef*Mlm(idx)*lmbda
                     fion(2,ii) = fion(2,ii) + dty*coef*Mlm(idx)*lmbda
                     fion(3,ii) = fion(3,ii) + dtz*coef*Mlm(idx)*lmbda
c                     idx = idx + 1
                  end do
               end do
!$OMP END DO
            end if
         end if
      end do
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

*     ***************************************
*     *                                     *
*     *           FMM_A_over_i              *
*     *                                     *
*     ***************************************
      complex*16 function FMM_A_over_i(l,m)
      implicit none
      integer l,m

*     **** local variables ****
      integer    i,mod_m
      real*8     A
      complex*16 img
       
 
      img  = dcmplx(0.0d0,1.0d0)
      mod_m = abs(m)
      A = 1
      do i=1,(l-mod_m)
        A = A*dble(i)
      end do
      A = A*A
      do i=(l-mod_m+1),(l+mod_m)
        A = A*dble(i)
      end do
      A = ((-1)**l)/dsqrt(A)

      FMM_A_over_i = A/img**mod_m
      return
      end

*     **********************************************
*     *                                            *
*     *              FMM_Translate_Mlm             *
*     *                                            *
*     **********************************************

      subroutine FMM_Translate_Mlm(lmax,Olm,Nlm,x,y,z)
      implicit none
      integer lmax
      complex*16 Olm(*)
      complex*16 Nlm(*)
      real*8 x,y,z

#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer l,m,idx,n
      integer rYlm(2),A(2),B(2),AB(2),tmp22(2),tmp42(2)

*     **** external functions ****
      complex*16 rSphericalHarmonic3
      external   rSphericalHarmonic3

      n = (2*lmax+2)*(4*lmax+2)

*     **** allocate stack ****
      value = BA_push_get(mt_dcpl,(lmax+1)**2,'rYlm',rYlm(2),rYlm(1))
      value=value.and.BA_push_get(mt_dcpl,n,'A',  A(2), A(1))
      value=value.and.BA_push_get(mt_dcpl,n,'B',  B(2), B(1))
      value=value.and.BA_push_get(mt_dcpl,n,'AB',AB(2),AB(1))
      value=value.and.
     >      BA_push_get(mt_dcpl,(2*lmax+2),'tmp22',tmp22(2),tmp22(1))
      value=value.and.
     >      BA_push_get(mt_dcpl,(4*lmax+2),'tmp42',tmp42(2),tmp42(1))
      if (.not.value)
     >   call errquit('FMM_Translate_Mlm:out of stack',0,MA_ERR)

      call dcffti(2*lmax+2,dcpl_mb(tmp22(1)))
      call dcffti(4*lmax+2,dcpl_mb(tmp42(1)))

      idx = 0
      do l=0,lmax
         do m=-l,l
            dcpl_mb(rYlm(1)+idx) = rSphericalHarmonic3(l,-m,x,y,z)
         end do
      end do

      call FMM_Translate_Mlm_sub1(lmax,Olm,dcpl_mb(rYlm(1)),
     >                           dcpl_mb(A(1)),dcpl_mb(B(1)))

      call FMM_Translate_sub2(lmax,dcpl_mb(A(1)),dcpl_mb(B(1)),
     >                        dcpl_mb(AB(1)),
     >                        dcpl_mb(tmp22(1)),dcpl_mb(tmp42(1)))

      call FMM_Translate_Mlm_sub3(lmax,dcpl_mb(AB(1)),Nlm)
            

*     **** deallocate stack ****
      value = BA_pop_stack(tmp42(2))
      value = value.and.BA_pop_stack(tmp22(2))
      value = value.and.BA_pop_stack(AB(2))
      value = value.and.BA_pop_stack(B(2))
      value = value.and.BA_pop_stack(A(2))
      value = value.and.BA_pop_stack(rYlm(2))
      if (.not.value)
     >   call errquit('FMM_Translate_Mlm:popping stack',0,MA_ERR)
      return
      end

*     ***************************************
*     *                                     *
*     *        FMM_Translate_Mlm_sub1       *
*     *                                     *
*     ***************************************
      subroutine FMM_Translate_Mlm_sub1(lmax,Olm,rYlm,A,B)
      implicit none
      integer lmax
      complex*16 Olm(*),rYlm(*)
      complex*16 A(2*lmax+2,4*lmax+2)
      complex*16 B(2*lmax+2,4*lmax+2)


*     **** local variables ****
      integer idx,i,j,l,m
      complex*16 ctmp

*     **** external functions ****
      complex*16 FMM_A_over_i
      external   FMM_A_over_i

      call dcopy(2*(2*lmax+2)*(4*lmax+2),0.0d0,0,A,1)
      call dcopy(2*(2*lmax+2)*(4*lmax+2),0.0d0,0,B,1)
      idx = 1
      j   = 1
      do l=0,lmax
         i=1
         do m=-l,l
            ctmp = FMM_A_over_i(l,m)
            A(i,j) = Olm(idx) *ctmp
            B(i,j) = rYlm(idx)*ctmp
             
            i = i + 1
            idx = idx + 1
         end do
         j = j + 1
      end do
      return
      end 

*     ***************************************
*     *                                     *
*     *        FMM_Translate_sub2           *
*     *                                     *
*     ***************************************
      subroutine FMM_Translate_sub2(lmax,A,B,AB,tmp22,tmp42)
      implicit none
      integer lmax
      complex*16  A(2*lmax+2,4*lmax+2)
      complex*16  B(2*lmax+2,4*lmax+2)
      complex*16 AB(2*lmax+2,4*lmax+2)
      complex*16 tmp22(*),tmp42(*)

*     **** local variables ****
      integer i,j,n

      n = (2*lmax+2)*(4*lmax+2)
      do j=1,(4*lmax+2)
         call dcfftf(2*lmax+2,A(1,j),tmp22)
         call dcfftf(2*lmax+2,B(1,j),tmp22)
      end do
      do i=1,(2*lmax+2)
         call dcfftf(4*lmax+2,A(i,1),tmp42)
         call dcfftf(4*lmax+2,B(i,1),tmp42)
      end do
      do j=1,(4*lmax+2)
         do i=1,(2*lmax+2)
            AB(i,j) = A(i,j)*B(i,j)
         end do
      end do

      do i=1,(2*lmax+2)
         call dcfftb(4*lmax+2,AB(i,1),tmp42)
      end do
      do j=1,(4*lmax+2)
         call dcfftb(2*lmax+2,AB(1,j),tmp22)
      end do
      call dscal(2*n,1.0d0/dble(n),AB,1)

      return
      end

*     ***************************************
*     *                                     *
*     *        FMM_Translate_Mlm_sub3       *
*     *                                     *
*     ***************************************
      subroutine FMM_Translate_Mlm_sub3(lmax,AB,Mlm)
      implicit none
      integer lmax
      complex*16 AB(2*lmax+2,4*lmax+2)
      complex*16 Mlm(*)

*     **** local variables ****
      integer idx,i,j,l,m

*     **** external functions ****
      complex*16 FMM_A_over_i
      external   FMM_A_over_i

* Jeff added these because they are implicit with Intel Fortran 
* and either this code is never used, is fine to give undefined behavior
* or assumes initialization to zero.
      idx = 0
      j   = 0
      i   = 0

      do l=0,lmax
         do m=-l,l
            Mlm(idx) = AB(i,j)/FMM_A_over_i(l,m)
            i = i + 1
            idx = idx + 1
         end do
         j = j + 1
      end do
      return
      end


*     **********************************************
*     *                                            *
*     *              FMM_Translate_Llm             *
*     *                                            *
*     **********************************************

      subroutine FMM_Translate_Llm(lmax,Olm,Nlm,x,y,z)
      implicit none
      integer lmax
      complex*16 Olm(*)
      complex*16 Nlm(*)
      real*8 x,y,z

#include "bafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      logical value
      integer l,m,idx,n
      integer rYlm(2),A(2),B(2),AB(2),tmp22(2),tmp42(2)

*     **** external functions ****
      complex*16 rSphericalHarmonic3
      external   rSphericalHarmonic3

      n = (2*lmax+2)*(4*lmax+2)

*     **** allocate stack ****
      value = BA_push_get(mt_dcpl,(lmax+1)**2,'rYlm',rYlm(2),rYlm(1))
      value=value.and.BA_push_get(mt_dcpl,n,'A',  A(2), A(1))
      value=value.and.BA_push_get(mt_dcpl,n,'B',  B(2), B(1))
      value=value.and.BA_push_get(mt_dcpl,n,'AB',AB(2),AB(1))
      value=value.and.
     >      BA_push_get(mt_dcpl,2*lmax+2,'tmp22',tmp22(2),tmp22(1))
      value=value.and.
     >      BA_push_get(mt_dcpl,4*lmax+2,'tmp42',tmp42(2),tmp42(1))
      if (.not.value)
     >   call errquit('FMM_Translate_Llm:out of stack',0,MA_ERR)

      call dcffti(2*lmax+2,dcpl_mb(tmp22(1)))
      call dcffti(4*lmax+2,dcpl_mb(tmp42(1)))

      idx = 0
      do l=0,lmax
         do m=-l,l
            dcpl_mb(rYlm(1)+idx) = rSphericalHarmonic3(l,m,x,y,z)
         end do
      end do

      call FMM_Translate_Llm_sub1(lmax,Olm,dcpl_mb(rYlm(1)),
     >                           dcpl_mb(A(1)),dcpl_mb(B(1)))

      call FMM_Translate_sub2(lmax,dcpl_mb(A(1)),dcpl_mb(B(1)),
     >                        dcpl_mb(AB(1)),
     >                        dcpl_mb(tmp22(1)),dcpl_mb(tmp42(1)))

      call FMM_Translate_Llm_sub3(lmax,dcpl_mb(AB(1)),Nlm)


*     **** deallocate stack ****
      value = BA_pop_stack(tmp42(2))
      value = value.and.BA_pop_stack(tmp22(2))
      value = value.and.BA_pop_stack(AB(2))
      value = value.and.BA_pop_stack(B(2))
      value = value.and.BA_pop_stack(A(2))
      value = value.and.BA_pop_stack(rYlm(2))
      if (.not.value)
     >   call errquit('FMM_Translate_Llm:popping stack',0,MA_ERR)
      return
      end


*     ***************************************
*     *                                     *
*     *        FMM_Translate_Llm_sub1       *
*     *                                     *
*     ***************************************
      subroutine FMM_Translate_Llm_sub1(lmax,Olm,rYlm,A,B)
      implicit none
      integer lmax
      complex*16 Olm(*),rYlm(*)
      complex*16 A(2*lmax+2,4*lmax+2)
      complex*16 B(2*lmax+2,4*lmax+2)


*     **** local variables ****
      integer idx,i,j,l,m,sgn
      complex*16 ctmp

*     **** external functions ****
      complex*16 FMM_A_over_i
      external   FMM_A_over_i

      call dcopy(2*(2*lmax+2)*(4*lmax+2),0.0d0,0,A,1)
      call dcopy(2*(2*lmax+2)*(4*lmax+2),0.0d0,0,B,1)
      idx = 1
      j   = 1
      sgn = 1
      do l=0,lmax
         i=1
         do m=-l,l
            ctmp = FMM_A_over_i(l,m)
            A(i,j) = sgn*Olm(idx)/ctmp
            B(i,j) = rYlm(idx)*ctmp

            i = i + 1
            idx = idx + 1
         end do
         j = j + 1
         sgn = -sgn
      end do
      return
      end



*     ***************************************
*     *                                     *
*     *        FMM_Translate_Llm_sub3       *
*     *                                     *
*     ***************************************
      subroutine FMM_Translate_Llm_sub3(lmax,AB,Llm)
      implicit none
      integer lmax
      complex*16 AB(2*lmax+2,4*lmax+2)
      complex*16 Llm(*)

*     **** local variables ****
      integer idx,i,j,l,m,sgn

*     **** external functions ****
      complex*16 FMM_A_over_i
      external   FMM_A_over_i

* Jeff: idx was not initialized before and defaults to zero with Intel,
* which appears to be the primary compiler for testing.
      idx = 0

      sgn = 1
      j=1
      do l=0,lmax
         i=1
         do m=-l,l
            Llm(idx) = sgn*AB(i,j)/FMM_A_over_i(l,m)
            i = i + 1
            idx = idx + 1
         end do
         j = j + 1
         sgn = -sgn
      end do
      return
      end


*     ***************************************
*     *                                     *
*     *        FMM_reduce_levels            *
*     *                                     *
*     ***************************************
*
*    This function translates index values between different levels of an octree grid

      integer function FMM_reduce_levels(nl,l,lr)
      implicit none
      integer nl,l,lr

      FMM_reduce_levels = ishft(ishft(ibits(nl,2*l,l),-lr),2*(l-lr))
     >                  + ishft(ishft(ibits(nl,l,l),-lr),(l-lr)) 
     >                  + ishft(ibits(nl,0,l),-lr)
      return
      end


