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
      integer taskid,np,pcount,gcount,nshl3d,tcount,gcount0
      integer l1,m1,l2,m2,iii,jjj,iia,jja
      integer tid,nthr,nn

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
         do m1=0,l1
            if (m1.eq.0) nn = 1
            if (m1.gt.0) nn = 2
            if (mod(pcount,np).eq.taskid) then
               tid = mod(gcount,nthr)
               int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + nn
               tcount = tcount + nn
               gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=0,l2
                  if ((m1.eq.0).and.(m2.eq.0)) nn = 1
                  if ((m1.eq.0).and.(m2.gt.0)) nn = 2
                  if ((m1.gt.0).and.(m2.eq.0)) nn = 2
                  if ((m1.gt.0).and.(m2.gt.0)) nn = 4
                  if (mod(pcount,np).eq.taskid) then
                     tid = mod(gcount,nthr)
                     int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + nn
                     tcount = tcount + nn
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
            do m1=0,l1
               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=0,l2
                  if ((m1.eq.0).and.(m2.eq.0)) nn = 1
                  if ((m1.eq.0).and.(m2.gt.0)) nn = 2
                  if ((m1.gt.0).and.(m2.eq.0)) nn = 2
                  if ((m1.gt.0).and.(m2.gt.0)) nn = 4
                  if (mod(pcount,np).eq.taskid) then
                     tid = mod(gcount,nthr)
                     int_mb(tgauss(1)+tid) = int_mb(tgauss(1)+tid) + nn
                     tcount = tcount + nn
                     gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end do
            end do
         end do
      end do
      ngauss_max = tcount
      ngauss     = tcount

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
      integer l1,m1,l2,m2
      integer l,nshl3d,rcell,rcell_hndl
      integer lm1(4),lm2(4),n,i,nn
      real*8  R1(3),R12(3),s1,s2,Rab(3),Rba(3),R,ss
      real*8  W1,W2,W3,W4,dW1(3),dW2(3),dW3(3),dW4(3)
      real*8  W(4),dW(3,4)
      real*8  e1(4),de1(3,4)

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
         do l1=0,int_mb(mult_l(1)+iia-1)
         do m1=0,l1
            if (m1.eq.0) nn = 1
            if (m1.gt.0) nn = 2
            if (mod(pcount,np).eq.taskid) then
            if (mod(gcount,nthr).eq.tid) then
                W1=nwpw_UGaussian(l1,m1,s1,l1,m1,s1)
                W2=nwpw_UGaussian(l1,m1,s1,l1,m1,sigma_smooth)
                W4=nwpw_UGaussian(l1,m1,sigma_smooth,l1,m1,sigma_smooth)
                e1(1) = 0.5d0*W1 + 0.5d0*W4 - W2
                lm1(1) = l1*(l1+1) + m1
                if (nn.gt.1) then
                   W1=nwpw_UGaussian(l1,-m1,s1,l1,-m1,s1)
                   W2=nwpw_UGaussian(l1,-m1,s1,l1,-m1,sigma_smooth)
                   W4=nwpw_UGaussian(l1,-m1,sigma_smooth,
     >                               l1,-m1,sigma_smooth)
                   e1(2) = 0.5d0*W1 + 0.5d0*W4 - W2
                   lm1(2) = l1*(l1+1) - m1 
                end if

c               !if (dabs(e1).gt.tole) then
                do i=1,nn
                   inds = indx + shft
                   dbl_mb(e_tauss(1)+inds) = e1(i)
                   int_mb(lm1_tauss(1)+inds) 
     >             = (iii-1)*2*lm_size_max+lm1(i)
                   int_mb(lm2_tauss(1)+inds) 
     >             = (iii-1)*2*lm_size_max+lm1(i)
                   int_mb(iii1_tauss(1)+inds) = iii
                   int_mb(iii2_tauss(1)+inds) = iii
                   indx = indx + 1
                end do
c               !end if
            end if
            gcount = gcount + 1
            end if
            pcount = pcount + 1

            if (nshl3d.gt.1) then
               do l2=0,int_mb(mult_l(1)+iia-1)
               do m2=0,l2
                  if ((m1.eq.0).and.(m2.eq.0)) nn = 1
                  if ((m1.eq.0).and.(m2.gt.0)) nn = 2
                  if ((m1.gt.0).and.(m2.eq.0)) nn = 2
                  if ((m1.gt.0).and.(m2.gt.0)) nn = 4
                  if (mod(pcount,np).eq.taskid) then
                  if (mod(gcount,nthr).eq.tid) then
                     do i=1,4
                        e1(i) = 0.0d0
                     end do
                     do l=2,nshl3d
                        Rab(1) = dbl_mb(rcell+l-1)
                        Rab(2) = dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                           call nwpw_WGaussian_block(l1,m1,s1,l2,m2,s1,
     >                                               sigma_smooth,Rab,
     >                                               n,lm1,lm2,W)
                           do i=1,n
                              e1(i) = e1(i) + 0.5d0*W(i)
                           end do
                        end if
                     end do

                     !if (dabs(e1).gt.tole) then
                     do i=1,nn
                        inds = indx + shft
                        dbl_mb(e_tauss(1)+inds) = e1(i)
                        int_mb(lm1_tauss(1)+inds) 
     >                        = (iii-1)*2*lm_size_max+lm1(i)
                        int_mb(lm2_tauss(1)+inds)
     >                        = (iii-1)*2*lm_size_max+lm2(i)
                        int_mb(iii1_tauss(1)+inds) = iii
                        int_mb(iii2_tauss(1)+inds) = iii

                        indx = indx + 1
                     end do
                     !end if
                  end if
                  gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do
            end if
         end do
         end do

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

            do l1=0,int_mb(mult_l(1)+iia-1)
            do m1=0,l1

               do l2=0,int_mb(mult_l(1)+jja-1)
               do m2=0,l2
                  if ((m1.eq.0).and.(m2.eq.0)) nn = 1
                  if ((m1.eq.0).and.(m2.gt.0)) nn = 2
                  if ((m1.gt.0).and.(m2.eq.0)) nn = 2
                  if ((m1.gt.0).and.(m2.gt.0)) nn = 4
                  if (mod(pcount,np).eq.taskid) then
                  if (mod(gcount,nthr).eq.tid) then
                     do i=1,4
                        e1(i)    = 0.0d0
                        de1(1,i) = 0.0d0
                        de1(2,i) = 0.0d0
                        de1(3,i) = 0.0d0
                     end do
                     do l=1,nshl3d
                        Rab(1) = R12(1) + dbl_mb(rcell+l-1)
                        Rab(2) = R12(2) + dbl_mb(rcell+l-1+nshl3d)
                        Rab(3) = R12(3) + dbl_mb(rcell+l-1+2*nshl3d)
                        R = dsqrt(Rab(1)**2 + Rab(2)**2 + Rab(3)**2)
                        if (R.lt.(4*sigma_smooth)) then
                        if (move) then
                           call nwpw_dWGaussian_block(l1,m1,s1,l2,m2,s2,
     >                                               sigma_smooth,Rab,
     >                                               n,lm1,lm2,W,dW)
                           do i=1,n
                              e1(i)  = e1(i)  + W(i)
                              de1(1,i) = de1(1,i) + dW(1,i)
                              de1(2,i) = de1(2,i) + dW(2,i)
                              de1(3,i) = de1(3,i) + dW(3,i)
                           end do
                        else
                           call nwpw_WGaussian_block(l1,m1,s1,l2,m2,s2,
     >                                               sigma_smooth,Rab,
     >                                               n,lm1,lm2,W)
                           do i=1,n
                              e1(i) = e1(i) + W(i)
                           end do
                        end if
                        end if
                     end do
                     !if (dabs(e1).gt.tole) then
                     do i=1,nn
                        inds = indx + shft
                        dbl_mb(e_tauss(1)+inds) = e1(i)
                        if (move) then
                           dbl_mb(f_tauss(1)+3*inds)   = de1(1,i)
                           dbl_mb(f_tauss(1)+3*inds+1) = de1(2,i)
                           dbl_mb(f_tauss(1)+3*inds+2) = de1(3,i)
                        end if
                        int_mb(lm1_tauss(1)+inds)
     >                     = (iii-1)*2*lm_size_max+lm1(i)
                        int_mb(lm2_tauss(1)+inds)
     >                     = (jjj-1)*2*lm_size_max+lm2(i)
                        int_mb(iii1_tauss(1)+inds) = iii
                        int_mb(iii2_tauss(1)+inds) = jjj

                        indx = indx + 1
                     end do
                     !end if
                  end if
                  gcount = gcount + 1
                  end if
                  pcount = pcount + 1
               end do
               end do

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


*     *******************************************************
*     *                                                     *
*     *             nwpw_WGaussian_block                    *
*     *                                                     *
*     *******************************************************
 
      subroutine nwpw_WGaussian_block(l1,mod_m1,sigma1,
     >                                l2,mod_m2,sigma2,
     >                                sigma_smooth,Rab,
     >                                n,lm1,lm2,W)
      implicit none
      integer l1,mod_m1
      real*8  sigma1
      integer l2,mod_m2
      real*8  sigma2
      real*8  sigma_smooth
      real*8  Rab(3)
      integer n
      integer lm1(*),lm2(*)
      real*8  W(*)

*     **** local variables ****
      complex*16 CW4,CW4p,CW4m,CW4pp,CW4pm,CW4mp,CW4mm


*     **** external functions ****
      complex*16 nwpw_CWGaussian
      external   nwpw_CWGaussian
      
      if ((mod_m1.eq.0).and.(mod_m2.eq.0)) then
         n  = 1
         lm1(1) = l1*(l1+1) 
         lm2(1) = l2*(l2+1) 
 
         CW4 = nwpw_CWGaussian(l1,mod_m1,sigma1,l2,mod_m2,sigma2,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma_smooth,l2,mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1,mod_m1,sigma_smooth,
     >                     l2,mod_m2,sigma_smooth,Rab)

         W(1) = dble(CW4)

      else if (mod_m1.eq.0) then
         n  = 2
         lm1(1) = l1*(l1+1)
         lm2(1) = l2*(l2+1) - mod_m2

         lm1(2) = l1*(l1+1)
         lm2(2) = l2*(l2+1) + mod_m2

         CW4p = nwpw_CWGaussian(l1,mod_m1,sigma1,l2,mod_m2,sigma2,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma_smooth,l2,mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1,mod_m1,sigma_smooth,
     >                     l2,mod_m2,sigma_smooth,Rab)

         CW4m = nwpw_CWGaussian(l1,mod_m1,sigma1,l2,-mod_m2,sigma2,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab)
     >   - nwpw_CWGaussian(l1,mod_m1,sigma_smooth,l2,-mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1,mod_m1,sigma_smooth,
     >                     l2,-mod_m2,sigma_smooth,Rab)

         !** m2<0 **
         CW4 = (CW4m - (-1.0d0)**mod_m2 * CW4p)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
         W(1) = dble(CW4)

         !** m2>0 **
         CW4 = (CW4m + (-1.0d0)**mod_m2 * CW4p)/dsqrt(2.0d0)
         W(2) = dble(CW4)

      else if (mod_m2.eq.0) then
         n  = 2
         lm1(1) = l1*(l1+1) - mod_m1
         lm2(1) = l2*(l2+1) 

         lm1(2) = l1*(l1+1) + mod_m1
         lm2(2) = l2*(l2+1) 

         CW4p = nwpw_CWGaussian(l1, mod_m1,sigma1,l2,mod_m2,sigma2,Rab)
     >   - nwpw_CWGaussian(l1, mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab)
     >   - nwpw_CWGaussian(l1, mod_m1,sigma_smooth,l2,mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1, mod_m1,sigma_smooth,
     >                     l2,mod_m2,sigma_smooth,Rab)

         CW4m = nwpw_CWGaussian(l1,-mod_m1,sigma1,l2,mod_m2,sigma2,Rab)
     >   - nwpw_CWGaussian(l1,-mod_m1,sigma1,l2,mod_m2,sigma_smooth,Rab)
     >   - nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,l2,mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,
     >                     l2,mod_m2,sigma_smooth,Rab)

         !** m1<0 **
         CW4 = (CW4m - (-1.0d0)**mod_m1 * CW4p)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
         W(1) = dble(CW4)

         !** m1>0 **
         CW4 = (CW4m + (-1.0d0)**mod_m1 * CW4p)/dsqrt(2.0d0)
         W(2) = dble(CW4)

      else
         n = 4
         lm1(1) = l1*(l1+1) - mod_m1
         lm2(1) = l2*(l2+1) - mod_m2

         lm1(2) = l1*(l1+1) - mod_m1
         lm2(2) = l2*(l2+1) + mod_m2

         lm1(3) = l1*(l1+1) + mod_m1
         lm2(3) = l2*(l2+1) - mod_m2

         lm1(4) = l1*(l1+1) + mod_m1
         lm2(4) = l2*(l2+1) + mod_m2

         CW4pp =nwpw_CWGaussian(l1, mod_m1,sigma1,l2, mod_m2,sigma2,Rab)
     >  - nwpw_CWGaussian(l1, mod_m1,sigma1,l2, mod_m2,sigma_smooth,Rab)
     >  - nwpw_CWGaussian(l1, mod_m1,sigma_smooth,l2, mod_m2,sigma2,Rab)
     >  + nwpw_CWGaussian(l1, mod_m1,sigma_smooth,
     >                    l2, mod_m2,sigma_smooth,Rab)

         CW4pm =nwpw_CWGaussian(l1, mod_m1,sigma1,l2,-mod_m2,sigma2,Rab)
     >  - nwpw_CWGaussian(l1, mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab)
     >  - nwpw_CWGaussian(l1, mod_m1,sigma_smooth,l2,-mod_m2,sigma2,Rab)
     >  + nwpw_CWGaussian(l1, mod_m1,sigma_smooth,
     >                    l2,-mod_m2,sigma_smooth,Rab)

         CW4mp =nwpw_CWGaussian(l1,-mod_m1,sigma1,l2, mod_m2,sigma2,Rab)
     >   -nwpw_CWGaussian(l1,-mod_m1,sigma1,l2, mod_m2,sigma_smooth,Rab)
     >   -nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,l2, mod_m2,sigma2,Rab)
     >   + nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,
     >                     l2, mod_m2,sigma_smooth,Rab)

         CW4mm =nwpw_CWGaussian(l1,-mod_m1,sigma1,l2,-mod_m2,sigma2,Rab)
     >  - nwpw_CWGaussian(l1,-mod_m1,sigma1,l2,-mod_m2,sigma_smooth,Rab)
     >  - nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,l2,-mod_m2,sigma2,Rab)
     >  + nwpw_CWGaussian(l1,-mod_m1,sigma_smooth,
     >                    l2,-mod_m2,sigma_smooth,Rab)

         !** m1<0 and m2<0 **
         CW4 = -(CW4mm
     >          + (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >          - (-1.0d0)**mod_m1          * CW4pm
     >          - (-1.0d0)**mod_m2          * CW4mp)/2.0d0
         W(1) = dble(CW4)

         !** m1<0 and m2>0 **
         CW4 = (CW4mm
     >         - (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >         - (-1.0d0)**mod_m1          * CW4pm
     >         + (-1.0d0)**mod_m2          * CW4mp)
     >         *dcmplx(0.0d0,1.0d0/2.0d0)
         W(2) = dble(CW4)

         !** m1>0 and m2<0 **
         CW4 = (CW4mm
     >         - (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >         + (-1.0d0)**mod_m1      * CW4pm
     >         - (-1.0d0)**mod_m2      * CW4mp)
     >        *dcmplx(0.0d0,1.0d0/2.0d0)
         W(3) = dble(CW4)

         !** m1>0 and m2>0 **
         CW4 =  (CW4mm
     >          + (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >          + (-1.0d0)**mod_m1      * CW4pm
     >          + (-1.0d0)**mod_m2      * CW4mp)/2.0d0
         W(4) = dble(CW4)

      end if

      return
      end


*     *******************************************************
*     *                                                     *
*     *             nwpw_dWGaussian_block                   *
*     *                                                     *
*     *******************************************************

      subroutine nwpw_dWGaussian_block(l1,mod_m1,sigma1,
     >                                 l2,mod_m2,sigma2,
     >                                 sigma_smooth,Rab,
     >                                 n,lm1,lm2,W,dW)
      implicit none
      integer l1,mod_m1
      real*8  sigma1
      integer l2,mod_m2
      real*8  sigma2
      real*8  sigma_smooth
      real*8  Rab(3)
      integer n
      integer lm1(*),lm2(*)
      real*8  W(*),dW(3,*)

*     **** local variables ****
      integer i
      complex*16 etmp,detmp(3)
      complex*16 CW4,CW4p,CW4m,CW4pp,CW4pm,CW4mp,CW4mm
      complex*16 dCW4(3),dCW4p(3),dCW4m(3)
      complex*16 dCW4pp(3),dCW4pm(3)
      complex*16 dCW4mp(3),dCW4mm(3)

      if ((mod_m1.eq.0).and.(mod_m2.eq.0)) then
         n  = 1
         lm1(1) = l1*(l1+1) 
         lm2(1) = l2*(l2+1) 
 
         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4     = etmp
         dCW4(1) = detmp(1)
         dCW4(2) = detmp(2)
         dCW4(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4     = CW4 - etmp
         dCW4(1) = dCW4(1) - detmp(1)
         dCW4(2) = dCW4(2) - detmp(2)
         dCW4(3) = dCW4(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4     = CW4 - etmp
         dCW4(1) = dCW4(1) - detmp(1)
         dCW4(2) = dCW4(2) - detmp(2)
         dCW4(3) = dCW4(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4     = CW4 + etmp
         dCW4(1) = dCW4(1) + detmp(1)
         dCW4(2) = dCW4(2) + detmp(2)
         dCW4(3) = dCW4(3) + detmp(3)


         W(1) = dble(CW4)
         dW(1,1) = dble(dCW4(1))
         dW(2,1) = dble(dCW4(2))
         dW(3,1) = dble(dCW4(3))

      else if (mod_m1.eq.0) then
         n  = 2
         lm1(1) = l1*(l1+1)
         lm2(1) = l2*(l2+1) - mod_m2

         lm1(2) = l1*(l1+1)
         lm2(2) = l2*(l2+1) + mod_m2

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4p     = etmp
         dCW4p(1) = detmp(1)
         dCW4p(2) = detmp(2)
         dCW4p(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p - etmp
         dCW4p(1) = dCW4p(1) - detmp(1)
         dCW4p(2) = dCW4p(2) - detmp(2)
         dCW4p(3) = dCW4p(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p - etmp
         dCW4p(1) = dCW4p(1) - detmp(1)
         dCW4p(2) = dCW4p(2) - detmp(2)
         dCW4p(3) = dCW4p(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p + etmp
         dCW4p(1) = dCW4p(1) + detmp(1)
         dCW4p(2) = dCW4p(2) + detmp(2)
         dCW4p(3) = dCW4p(3) + detmp(3)




         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                        l2,-mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4m     = etmp
         dCW4m(1) = detmp(1)
         dCW4m(2) = detmp(2)
         dCW4m(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                        l2,-mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m - etmp
         dCW4m(1) = dCW4m(1) - detmp(1)
         dCW4m(2) = dCW4m(2) - detmp(2)
         dCW4m(3) = dCW4m(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                        l2,-mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m - etmp
         dCW4m(1) = dCW4m(1) - detmp(1)
         dCW4m(2) = dCW4m(2) - detmp(2)
         dCW4m(3) = dCW4m(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                        l2,-mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m + etmp
         dCW4m(1) = dCW4m(1) + detmp(1)
         dCW4m(2) = dCW4m(2) + detmp(2)
         dCW4m(3) = dCW4m(3) + detmp(3)

         !** m2<0 **
         CW4 = (CW4m - (-1.0d0)**mod_m2 * CW4p)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
         W(1) = dble(CW4)
         do i=1,3
            CW4 = (dCW4m(i) - (-1.0d0)**mod_m2 * dCW4p(i))
     >            *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
            dW(i,1) = dble(CW4)
         end do

         !** m2>0 **
         CW4 = (CW4m + (-1.0d0)**mod_m2 * CW4p)/dsqrt(2.0d0)
         W(2) = dble(CW4)
         do i=1,3
            CW4 = (dCW4m(i) + (-1.0d0)**mod_m2 * dCW4p(i))/dsqrt(2.0d0)
            dW(i,2) = dble(CW4)
         end do


      else if (mod_m2.eq.0) then
         n  = 2
         lm1(1) = l1*(l1+1) - mod_m1
         lm2(1) = l2*(l2+1) 

         lm1(2) = l1*(l1+1) + mod_m1
         lm2(2) = l2*(l2+1) 


         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                        l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4p     = etmp
         dCW4p(1) = detmp(1)
         dCW4p(2) = detmp(2)
         dCW4p(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                        l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p - etmp
         dCW4p(1) = dCW4p(1) - detmp(1)
         dCW4p(2) = dCW4p(2) - detmp(2)
         dCW4p(3) = dCW4p(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                        l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p - etmp
         dCW4p(1) = dCW4p(1) - detmp(1)
         dCW4p(2) = dCW4p(2) - detmp(2)
         dCW4p(3) = dCW4p(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                        l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4p     = CW4p + etmp
         dCW4p(1) = dCW4p(1) + detmp(1)
         dCW4p(2) = dCW4p(2) + detmp(2)
         dCW4p(3) = dCW4p(3) + detmp(3)

        call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                        l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4m     = etmp
         dCW4m(1) = detmp(1)
         dCW4m(2) = detmp(2)
         dCW4m(3) = detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                        l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m - etmp
         dCW4m(1) = dCW4m(1) - detmp(1)
         dCW4m(2) = dCW4m(2) - detmp(2)
         dCW4m(3) = dCW4m(3) - detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                        l2,mod_m2,sigma2,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m - etmp
         dCW4m(1) = dCW4m(1) - detmp(1)
         dCW4m(2) = dCW4m(2) - detmp(2)
         dCW4m(3) = dCW4m(3) - detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                        l2,mod_m2,sigma_smooth,
     >                        Rab,etmp,detmp)
         CW4m     = CW4m + etmp
         dCW4m(1) = dCW4m(1) + detmp(1)
         dCW4m(2) = dCW4m(2) + detmp(2)
         dCW4m(3) = dCW4m(3) + detmp(3)


         !** m1<0 **
         CW4 = (CW4m - (-1.0d0)**mod_m1 * CW4p)
     >         *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
         W(1) = dble(CW4)
         do i=1,3
            CW4 = (dCW4m(i) - (-1.0d0)**mod_m1 * dCW4p(i))
     >            *dcmplx(0.0d0,1.0d0/dsqrt(2.0d0))
            dW(i,1) = dble(CW4)
         end do

         !** m1>0 **
         CW4 = (CW4m + (-1.0d0)**mod_m1 * CW4p)/dsqrt(2.0d0)
         W(2) = dble(CW4)
         do i=1,3
            CW4 = (dCW4m(i) + (-1.0d0)**mod_m1 * dCW4p(i))/dsqrt(2.0d0)
            dW(i,2) = dble(CW4)
         end do

      else
         n = 4
         lm1(1) = l1*(l1+1) - mod_m1
         lm2(1) = l2*(l2+1) - mod_m2

         lm1(2) = l1*(l1+1) - mod_m1
         lm2(2) = l2*(l2+1) + mod_m2

         lm1(3) = l1*(l1+1) + mod_m1
         lm2(3) = l2*(l2+1) - mod_m2

         lm1(4) = l1*(l1+1) + mod_m1
         lm2(4) = l2*(l2+1) + mod_m2

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4pp = etmp
         dCW4pp(1) = detmp(1)
         dCW4pp(2) = detmp(2)
         dCW4pp(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4pp = CW4pp - etmp
         dCW4pp(1) = dCW4pp(1) - detmp(1)
         dCW4pp(2) = dCW4pp(2) - detmp(2)
         dCW4pp(3) = dCW4pp(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4pp = CW4pp - etmp
         dCW4pp(1) = dCW4pp(1) - detmp(1)
         dCW4pp(2) = dCW4pp(2) - detmp(2)
         dCW4pp(3) = dCW4pp(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4pp = CW4pp + etmp
         dCW4pp(1) = dCW4pp(1) + detmp(1)
         dCW4pp(2) = dCW4pp(2) + detmp(2)
         dCW4pp(3) = dCW4pp(3) + detmp(3)


         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,-mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4pm = etmp
         dCW4pm(1) = detmp(1)
         dCW4pm(2) = detmp(2)
         dCW4pm(3) = detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma1,
     >                         l2,-mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4pm = CW4pm - etmp
         dCW4pm(1) = dCW4pm(1) - detmp(1)
         dCW4pm(2) = dCW4pm(2) - detmp(2)
         dCW4pm(3) = dCW4pm(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,-mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4pm = CW4pm - etmp
         dCW4pm(1) = dCW4pm(1) - detmp(1)
         dCW4pm(2) = dCW4pm(2) - detmp(2)
         dCW4pm(3) = dCW4pm(3) - detmp(3)

         call nwpw_dCWGaussian(l1,mod_m1,sigma_smooth,
     >                         l2,-mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4pm = CW4pm + etmp
         dCW4pm(1) = dCW4pm(1) + detmp(1)
         dCW4pm(2) = dCW4pm(2) + detmp(2)
         dCW4pm(3) = dCW4pm(3) + detmp(3)



         call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4mp = etmp
         dCW4mp(1) = detmp(1)
         dCW4mp(2) = detmp(2)
         dCW4mp(3) = detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4mp = CW4mp - etmp
         dCW4mp(1) = dCW4mp(1) - detmp(1)
         dCW4mp(2) = dCW4mp(2) - detmp(2)
         dCW4mp(3) = dCW4mp(3) - detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4mp = CW4mp - etmp
         dCW4mp(1) = dCW4mp(1) - detmp(1)
         dCW4mp(2) = dCW4mp(2) - detmp(2)
         dCW4mp(3) = dCW4mp(3) - detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                         l2,mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4mp = CW4mp + etmp
         dCW4mp(1) = dCW4mp(1) + detmp(1)
         dCW4mp(2) = dCW4mp(2) + detmp(2)
         dCW4mp(3) = dCW4mp(3) + detmp(3)


         call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                         l2,-mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4mm = etmp
         dCW4mm(1) = detmp(1)
         dCW4mm(2) = detmp(2)
         dCW4mm(3) = detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma1,
     >                         l2,-mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4mm = CW4mm - etmp
         dCW4mm(1) = dCW4mm(1) - detmp(1)
         dCW4mm(2) = dCW4mm(2) - detmp(2)
         dCW4mm(3) = dCW4mm(3) - detmp(3)

         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                         l2,-mod_m2,sigma2,
     >                         Rab,etmp,detmp)
         CW4mm = CW4mm - etmp
         dCW4mm(1) = dCW4mm(1) - detmp(1)
         dCW4mm(2) = dCW4mm(2) - detmp(2)
         dCW4mm(3) = dCW4mm(3) - detmp(3)
     
         call nwpw_dCWGaussian(l1,-mod_m1,sigma_smooth,
     >                         l2,-mod_m2,sigma_smooth,
     >                         Rab,etmp,detmp)
         CW4mm = CW4mm + etmp
         dCW4mm(1) = dCW4mm(1) + detmp(1)
         dCW4mm(2) = dCW4mm(2) + detmp(2)
         dCW4mm(3) = dCW4mm(3) + detmp(3)

         !** m1<0 and m2<0 **
         CW4 = -(CW4mm
     >          + (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >          - (-1.0d0)**mod_m1          * CW4pm
     >          - (-1.0d0)**mod_m2          * CW4mp)/2.0d0
         W(1) = dble(CW4)
         do i=1,3
            CW4 = -(dCW4mm(i)
     >          + (-1.0d0)**(mod_m1+mod_m2) * dCW4pp(i)
     >          - (-1.0d0)**mod_m1          * dCW4pm(i)
     >          - (-1.0d0)**mod_m2          * dCW4mp(i))/2.0d0
            dW(i,1) = dble(CW4)
         end do

         !** m1<0 and m2>0 **
         CW4 = (CW4mm
     >         - (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >         - (-1.0d0)**mod_m1          * CW4pm
     >         + (-1.0d0)**mod_m2     * CW4mp)*dcmplx(0.0d0,1.0d0/2.0d0)
         W(2) = dble(CW4)
         do i=1,3
            CW4 = (dCW4mm(i)
     >         - (-1.0d0)**(mod_m1+mod_m2) * dCW4pp(i)
     >         - (-1.0d0)**mod_m1          * dCW4pm(i)
     >         + (-1.0d0)**mod_m2 * dCW4mp(i))*dcmplx(0.0d0,1.0d0/2.0d0)
            dW(i,2) = dble(CW4)
         end do

         !** m1>0 and m2<0 **
         CW4 = (CW4mm
     >         - (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >         + (-1.0d0)**mod_m1     * CW4pm
     >         - (-1.0d0)**mod_m2     * CW4mp)*dcmplx(0.0d0,1.0d0/2.0d0)
         W(3) = dble(CW4)
         do i=1,3
            CW4 = (dCW4mm(i)
     >         - (-1.0d0)**(mod_m1+mod_m2) * dCW4pp(i)
     >         + (-1.0d0)**mod_m1 * dCW4pm(i)
     >         - (-1.0d0)**mod_m2 * dCW4mp(i))*dcmplx(0.0d0,1.0d0/2.0d0)
            dW(i,3) = dble(CW4)
         end do

         !** m1>0 and m2>0 **
         CW4 =  (CW4mm
     >          + (-1.0d0)**(mod_m1+mod_m2) * CW4pp
     >          + (-1.0d0)**mod_m1      * CW4pm
     >          + (-1.0d0)**mod_m2      * CW4mp)/2.0d0
         W(4) = dble(CW4)
         do i=1,3
            CW4 =  (dCW4mm(i)
     >          + (-1.0d0)**(mod_m1+mod_m2) * dCW4pp(i)
     >          + (-1.0d0)**mod_m1      * dCW4pm(i)
     >          + (-1.0d0)**mod_m2      * dCW4mp(i))/2.0d0
            dW(i,4) = dble(CW4)
         end do

      end if

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
      integer i
      real*8 tole
      parameter (tole=1.0d-25)

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
