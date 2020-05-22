********************************* real Gaussian Integrals **************************************

*     *************************************************
*     *                                               *
*     *         nwpw_gintegrals_set_gcount            *
*     *                                               *
*     *************************************************

      subroutine nwpw_gintegrals_set_gcount()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer taskid,np,pcount,gcount,nshl3d,tcount
      integer l1,m1,l2,m2,iii,jjj,iia,jja
      integer tid,nthr

*     **** external functions ****
      integer  control_version,ewald_nshl3d
      external control_version,ewald_nshl3d
      integer  Parallel_maxthreads
      external Parallel_maxthreads

      call Parallel_taskid(taskid)
      call Parallel_np(np)
      nthr = Parallel_maxthreads()
      do l1 = 1,nthr
         int_mb(tgauss(1)+l1-1) = 0
      end do

      periodic = (control_version().eq.3)

      if (periodic) then
         nshl3d = ewald_nshl3d()
      else
         nshl3d = 1
      end if

      pcount = 0
      gcount = 0
      tcount = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)

*        **** calculate on-site integrals ****
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
               tid = mod(gcount,nthr)
               int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + 1
               gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     tid = mod(gcount,nthr)
                     int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + 1
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end if
         end do
         end do

*        **** calculate IJ integrals ****
         do jjj=iii+1,nion_paw
            jja = int_mb(katm_paw(1)+jjj-1)

            do l1=0,int_mb(mult_l(1)+iia-1)
            do m1=-l1,l1
               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                     tid = mod(gcount,nthr)
                     int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + 1
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end do
            end do
         end do
      end do
      ngauss_max = gcount
      ngauss     = gcount


      tcount = 0
      do l1 = 1,nthr
         int_mb(tgauss_shift(1)+l1-1) = tcount
         tcount = tcount + int_mb(tgauss(1)+l1-1)
      end do

      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_gintegrals_init                *
*     *                                               *
*     *************************************************

      subroutine nwpw_gintegrals_init()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value
      integer nthr

*     **** external functions ****
      integer  Parallel_maxthreads
      external Parallel_maxthreads

      nthr = Parallel_maxthreads()
      value =           BA_alloc_get(mt_int,nthr,"tgauss",
     >                               tgauss(2),tgauss(1))
      value = value.and.BA_alloc_get(mt_int,nthr,"tgauss_shift",
     >                               tgauss_shift(2),tgauss_shift(1))

      call nwpw_gintegrals_set_gcount()

      value = value.and.BA_alloc_get(mt_int,ngauss_max,"lm1_gauss",
     >                              lm1_gauss(2),lm1_gauss(1))
      value = value.and.BA_alloc_get(mt_int,ngauss_max,"lm2_gauss",
     >                              lm2_gauss(2),lm2_gauss(1))
      value = value.and.BA_alloc_get(mt_int,ngauss_max,"iii1_gauss",
     >                              iii1_gauss(2),iii1_gauss(1))
      value = value.and.BA_alloc_get(mt_int,ngauss_max,"iii2_gauss",
     >                              iii2_gauss(2),iii2_gauss(1))
      value = value.and.BA_alloc_get(mt_dbl,ngauss_max,"e_gauss",
     >                              e_gauss(2),e_gauss(1))
      value = value.and.BA_alloc_get(mt_dbl,3*ngauss_max,"f_gauss",
     >                              f_gauss(2),f_gauss(1))
      if (.not.value)
     > call errquit("nwpw_gintegrals_init:cannot allocate memory",
     >             0,MA_ERR)

      return
      end


*     *************************************************
*     *                                               *
*     *           nwpw_gintegrals_end                 *
*     *                                               *
*     *************************************************
      subroutine nwpw_gintegrals_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      logical value

      value =           BA_free_heap(tgauss(2))
      value = value.and.BA_free_heap(tgauss_shift(2))
      value = value.and.BA_free_heap(lm1_gauss(2))
      value = value.and.BA_free_heap(lm2_gauss(2))
      value = value.and.BA_free_heap(iii1_gauss(2))
      value = value.and.BA_free_heap(iii2_gauss(2))
      value = value.and.BA_free_heap(e_gauss(2))
      value = value.and.BA_free_heap(f_gauss(2))
      if (.not.value)
     > call errquit("nwpw_gintegrals_end:cannot allocate memory",
     >             0,MA_ERR)

      return
      end



*     *************************************************
*     *                                               *
*     *             nwpw_gintegrals_set               *
*     *                                               *
*     *************************************************
c
c  The logic of this routine needs to be completely reworked for threading.
c  It's well designed for MPI parallelism, so one option is to expand all the
c  data structures over tasks and threads instead of just tasks.
c  Another option is to define thread shifts for indx.... However the threshold
c  check with tole would have to be eliminated.
c
      subroutine nwpw_gintegrals_set(move)
      implicit none
      logical move

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     ***** local variables ****
      real*8 tole
      parameter (tole=1.0d-25)

      logical value
      integer taskid,np,pcount,gcount
      integer ii,jj,ia,ja,indx,shft,inds
      integer iii,jjj,iia,jja
      integer lm1,l1,m1,lm2,l2,m2
      integer l,nshl3d,rcell,rcell_hndl
      real*8 R1(3),R12(3),s1,s2,Rab(3),Rba(3),R,ss
      real*8 W1,W2,W3,W4,dW1(3),dW2(3),dW3(3),dW4(3)
      real*8 e1,de1(3)

      integer tid,nthr
      integer lm1_tauss(2),lm2_tauss(2),iii1_tauss(2),iii2_tauss(2)
      integer e_tauss(2),f_tauss(2)

*     **** external functions ****
      real*8   ion_rion,nwpw_WGaussian,nwpw_UGaussian
      external ion_rion,nwpw_WGaussian,nwpw_UGaussian
      integer  nwpw_doublefactorial
      external nwpw_doublefactorial
      integer  ewald_nshl3d,ewald_rcell_ptr
      external ewald_nshl3d,ewald_rcell_ptr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(34)
      call Parallel_taskid(taskid)
      call Parallel_np(np)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      shft = int_mb(tgauss_shift(1)+tid)

c     **** allocate temporary memory ****
      value =           BA_push_get(mt_int,ngauss_max,"lm1_tauss",
     >                              lm1_tauss(2),lm1_tauss(1))
      value = value.and.BA_push_get(mt_int,ngauss_max,"lm2_tauss",
     >                              lm2_tauss(2),lm2_tauss(1))
      value = value.and.BA_push_get(mt_int,ngauss_max,"iii1_tauss",
     >                              iii1_tauss(2),iii1_tauss(1))
      value = value.and.BA_push_get(mt_int,ngauss_max,"iii2_tauss",
     >                              iii2_tauss(2),iii2_tauss(1))
      value = value.and.BA_push_get(mt_dbl,ngauss_max,"e_tauss",
     >                              e_tauss(2),e_tauss(1))
      value = value.and.BA_push_get(mt_dbl,3*ngauss_max,"f_tauss",
     >                              f_tauss(2),f_tauss(1))

      if (periodic) then
         nshl3d = ewald_nshl3d()
         rcell  = ewald_rcell_ptr()
      else
         if (.not. BA_push_get(mt_dbl,3,"rcellflm",rcell_hndl,rcell))
     >   call errquit("nwpw_compcharge_set_gintegrals:stack",1,MA_ERR)

         nshl3d = 1
!$OMP MASTER
         dbl_mb(rcell)   = 0.0d0
         dbl_mb(rcell+1) = 0.0d0
         dbl_mb(rcell+2) = 0.0d0
!$OMP END MASTER
!$OMP BARRIER
      end if


      !**** these should not need to be called!!!! ****
!$OMP MASTER
      call dcopy(ngauss_max,0.0d0,0,dbl_mb(e_tauss(1)),1)
      call dcopy(3*ngauss_max,0.0d0,0,dbl_mb(f_tauss(1)),1)
      call dcopy(ngauss_max,0.0d0,0,dbl_mb(e_gauss(1)),1)
      call dcopy(3*ngauss_max,0.0d0,0,dbl_mb(f_gauss(1)),1)
!$OMP END MASTER
!$OMP BARRIER

      pcount = 0
      gcount = 0
      indx   = 0
      do iii=1,nion_paw
         iia = int_mb(katm_paw(1)+iii-1)
         s1  = dbl_mb(sigma_paw(1)+iia-1)

*        **** calculate on-site integrals ****
         lm1 = 0
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=-l1,l1
            if (mod(pcount,np).eq.taskid) then
            if (mod(gcount,nthr).eq.tid) then
                W1=nwpw_UGaussian(l1,m1,s1,l1,m1,s1)
                W2=nwpw_UGaussian(l1,m1,s1,l1,m1,sigma_smooth)
                W4=nwpw_UGaussian(l1,m1,sigma_smooth,l1,m1,sigma_smooth)
                e1 = 0.5d0*W1 + 0.5d0*W4 - W2
c               !if (dabs(e1).gt.tole) then
                inds = indx + shft
                dbl_mb(e_tauss(1)+inds) = e1
                int_mb(lm1_tauss(1)+inds) = (iii-1)*2*lm_size_max+lm1
                int_mb(lm2_tauss(1)+inds) = (iii-1)*2*lm_size_max+lm1
                int_mb(iii1_tauss(1)+inds) = iii
                int_mb(iii2_tauss(1)+inds) = iii

                indx = indx + 1
c               !end if
            end if
            gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               lm2 = 0
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                  if (mod(gcount,nthr).eq.tid) then
                     e1 = 0.0d0
                     do l=2,nshl3d
                        Rab(1) = dbl_mb(rcell+l-1)
                        Rab(2) = dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                           W1 = nwpw_WGaussian(l1,m1,s1,l2,m2,s1,Rab)
                           W2 = nwpw_WGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
c                           W3 = nwpw_WGaussian(l1,m1,sigma_smooth,
c     >                                         l2,m2,s1,Rab)
                           W4 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           !e1=e1+0.5d0*W1+0.5d0*W4-0.5d0*W2-0.5d0*W3
                           !e1=e1+(W1+W4-2.0d0*W2)
                           !e1=e1+(W1+W4-W2-W3)
                           e1=e1+0.5d0*W1+0.5d0*W4-W2
                        end if
                     end do

                     !if (dabs(e1).gt.tole) then
                     inds = indx + shft

                     dbl_mb(e_tauss(1)+inds) = e1
                     int_mb(lm1_tauss(1)+inds) 
     >                     = (iii-1)*2*lm_size_max+lm1
                     int_mb(lm2_tauss(1)+inds)
     >                     = (iii-1)*2*lm_size_max+lm2
                     int_mb(iii1_tauss(1)+inds) = iii
                     int_mb(iii2_tauss(1)+inds) = iii

                     indx = indx + 1
                     !end if
                  end if
                  gcount = gcount + 1
                  end if
                  pcount = pcount + 1
                  lm2 = lm2 + 1
               end do !*** end m2 ***
               end do !*** end l2 ***
            end if

            lm1 = lm1 + 1
         end do !*** end m1***
         end do !*** end l1***


*        **** calculate IJ integrals ****
         ii  = int_mb(ion_pawtoion(1)+iii-1)
         R1(1) = ion_rion(1,ii)
         R1(2) = ion_rion(2,ii)
         R1(3) = ion_rion(3,ii)
         do jjj=iii+1,nion_paw
            jja = int_mb(katm_paw(1)+jjj-1)
            s2  = dbl_mb(sigma_paw(1)+jja-1)

            jj  = int_mb(ion_pawtoion(1)+jjj-1)
            R12(1) = R1(1) - ion_rion(1,jj)
            R12(2) = R1(2) - ion_rion(2,jj)
            R12(3) = R1(3) - ion_rion(3,jj)

            lm1 = 0
            do l1=0,int_mb(mult_l(1)+iia-1)
            do m1=-l1,l1

               lm2 = 0
               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=-l2,l2
                  if (mod(pcount,np).eq.taskid) then
                  if (mod(gcount,nthr).eq.tid) then
                     e1 = 0.0d0
                     de1(1) = 0.0d0
                     de1(2) = 0.0d0
                     de1(3) = 0.0d0
                     do l=1,nshl3d
                        Rab(1) = R12(1) + dbl_mb(rcell+l-1)
                        Rab(2) = R12(2) + dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = R12(3) + dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                        if (move) then
                           call nwpw_dWGaussian(l1,m1,s1,l2,m2,s2,Rab,
     >                                          W1,dW1)
                           call nwpw_dWGaussian(l1,m1,s1,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W2,dW2)
                           call nwpw_dWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,s2,Rab,
     >                                          W3,dW3)
                           call nwpw_dWGaussian(l1,m1,sigma_smooth,
     >                                          l2,m2,sigma_smooth,Rab,
     >                                          W4,dW4)
                           e1 = e1 + (W1 + W4 - W2 - W3)
                           de1(1)=de1(1)+(dW1(1)+dW4(1)-dW2(1)-dW3(1))
                           de1(2)=de1(2)+(dW1(2)+dW4(2)-dW2(2)-dW3(2))
                           de1(3)=de1(3)+(dW1(3)+dW4(3)-dW2(3)-dW3(3))

                        else
                           W1 = nwpw_WGaussian(l1,m1,s1,l2,m2,s2,Rab)
                           W2 = nwpw_WGaussian(l1,m1,s1,
     >                                         l2,m2,sigma_smooth,Rab)
                           W3 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,s2,Rab)
                           W4 = nwpw_WGaussian(l1,m1,sigma_smooth,
     >                                         l2,m2,sigma_smooth,Rab)
                           e1 = e1 + (W1 + W4 - W2 - W3)
                        end if
                        end if
                     end do
                     !if (dabs(e1).gt.tole) then
                     inds = indx + shft

                     dbl_mb(e_tauss(1)+inds) = e1
                     if (move) then
                        dbl_mb(f_tauss(1)+3*inds)   = de1(1)
                        dbl_mb(f_tauss(1)+3*inds+1) = de1(2)
                        dbl_mb(f_tauss(1)+3*inds+2) = de1(3)
                     end if
                     int_mb(lm1_tauss(1)+inds)
     >                  = (iii-1)*2*lm_size_max+lm1
                     int_mb(lm2_tauss(1)+inds)
     >                  = (jjj-1)*2*lm_size_max+lm2
                     int_mb(iii1_tauss(1)+inds) = iii
                     int_mb(iii2_tauss(1)+inds) = jjj

                     indx = indx + 1
                     !end if
                  end if
                  gcount = gcount + 1
                  end if
                  pcount = pcount + 1
                  lm2 = lm2 + 1
               end do
               end do

               lm1 = lm1 + 1
            end do
            end do

         end do
      end do

!$OMP BARRIER
!$OMP MASTER
      call nwpw_gintegral_stripper(ngauss_max,
     >                             int_mb(iii1_tauss(1)),
     >                             int_mb(iii2_tauss(1)),
     >                             int_mb(lm1_tauss(1)),
     >                             int_mb(lm2_tauss(1)),
     >                             dbl_mb(e_tauss(1)),
     >                             dbl_mb(f_tauss(1)),
     >                             ngauss,
     >                             int_mb(iii1_gauss(1)),
     >                             int_mb(iii2_gauss(1)),
     >                             int_mb(lm1_gauss(1)),
     >                             int_mb(lm2_gauss(1)),
     >                             dbl_mb(e_gauss(1)),
     >                             dbl_mb(f_gauss(1)))
!$OMP END MASTER
!$OMP BARRIER
c Need to have barrier before deallocation below.  This barrier could be removed
c if the tauss variables were on the heap and not deallocated rather than stack.


c     **** deallocate stack memory ****
      if (.not.periodic) then
        if (.not.BA_pop_stack(rcell_hndl))
     >   call errquit("nwpw_compcharge_set_gintegrals:stack",2,MA_ERR)
      end if
      value =           BA_pop_stack(f_tauss(2))
      value = value.and.BA_pop_stack(e_tauss(2))
      value = value.and.BA_pop_stack(iii2_tauss(2))
      value = value.and.BA_pop_stack(iii1_tauss(2))
      value = value.and.BA_pop_stack(lm2_tauss(2))
      value = value.and.BA_pop_stack(lm1_tauss(2))
        if (.not.value)
     >   call errquit("nwpw_compcharge_set_gintegrals:stack",3,MA_ERR)

      call nwpw_timing_end(34)
      return
      end

*     *************************************************
*     *                                               *
*     *             nwpw_gintegrals_stripper          *
*     *                                               *
*     *************************************************
c
c  This routine is used to remove unecessary integrals
c
      subroutine nwpw_gintegral_stripper(ng_in,
     >                                   iii1_in,iii2_in,
     >                                   lm1_in,lm2_in,e_in,f_in,
     >                                   ng_out,
     >                                   iii1_out,iii2_out,
     >                                   lm1_out,lm2_out,e_out,f_out)
      implicit none
      integer ng_in,iii1_in(*),iii2_in(*),lm1_in(*),lm2_in(*)
      real*8  e_in(*),f_in(*)
      integer ng_out,iii1_out(*),iii2_out(*),lm1_out(*),lm2_out(*)
      real*8  e_out(*),f_out(*)

c     **** local variables ****
      real*8 tole
      parameter (tole=1.0d-25)
      integer i

      ng_out = 0
      do i=1,ng_in
         if (dabs(e_in(i)).gt.tole) then
            ng_out = ng_out + 1
            iii1_out(ng_out) = iii1_in(i)
            iii2_out(ng_out) = iii2_in(i)
            lm1_out(ng_out)  = lm1_in(i)
            lm2_out(ng_out)  = lm2_in(i)
            e_out(ng_out)    = e_in(i)
            f_out(3*(ng_out-1)+1) = f_in(3*(i-1)+1)
            f_out(3*(ng_out-1)+2) = f_in(3*(i-1)+2)
            f_out(3*(ng_out-1)+3) = f_in(3*(i-1)+3)
         end if
      end do
      return
      end


********************************* real Gaussian Integrals **************************************

