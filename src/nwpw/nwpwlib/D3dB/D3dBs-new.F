ccccccccccc  #define TCGMSG
#define NBLOCKS 3



*     ***********************************
*     *					*
*     *	        D3dBs_(c,r,t)_Zero	*	
*     *					*
*     ***********************************
      subroutine D3dBs_c_Zero(nb,A)     
      implicit none 
      integer nb
      complex A(*)

#include "D3dB.fh"

      call Parallel_shared_vector_szero(.true.,n2ft3d_map(nb),A)
      return
      end

      subroutine D3dBs_c_nZero(nb,n,A)
      implicit none
      integer nb,n
      complex A(*)

#include "D3dB.fh"

      call Parallel_shared_vector_szero(.true.,n*n2ft3d_map(nb),A)
      return
      end


      subroutine D3dBs_r_Zero(nb,A)     
      implicit none 
      integer nb
      real    A(*)

#include "D3dB.fh"

      call Parallel_shared_vector_szero(.true.,n2ft3d_map(nb),A)
      return
      end


      subroutine D3dBs_r_nZero(nb,n,A)
      implicit none
      integer nb,n
      real    A(*)

#include "D3dB.fh"

      call Parallel_shared_vector_szero(.true.,n*n2ft3d_map(nb),A)
      return
      end


*     ***********************************
*     *					*
*     *	        D3dBs_(c,r,t)_Copy	*	
*     *					*
*     ***********************************
      subroutine D3dBs_c_Copy(nb,A,B)     
      implicit none 
      integer nb
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,2*nfft3d_map(nb),A,B)
      return
      end

      subroutine D3dBs_c_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,n*2*nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dBs_r_Copy(nb,A,B)     
      implicit none 
      integer nb
      real A(*)
      real B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,n2ft3d_map(nb),A,B)
      return
      end


      subroutine D3dBs_r_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      real A(*)
      real B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,n*n2ft3d_map(nb),A,B)
      return
      end


      subroutine D3dBs_t_Copy(nb,A,B)     
      implicit none 
      integer nb
      real A(*)
      real B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dBs_t_nCopy(nb,n,A,B)
      implicit none
      integer nb,n
      real A(*)
      real B(*)

#include "D3dB.fh"

      call Parallel_shared_vector_scopy(.true.,n*nfft3d_map(nb),A,B)
      return
      end


      subroutine D3dBs_tc_Copy(nb,A,B)
      implicit none
      integer nb
      real    A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = cmplx(A(i),0.0)
      end do
!$OMP END DO

      return
      end

      subroutine D3dBs_ct_Copy(nb,A,B)
      implicit none
      integer nb
      complex A(*)
      real    B(*)

#include "D3dB.fh"

      integer i
!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = real(A(i))
      end do
!$OMP END DO

      return
      end



*     ***********************************
*     *					*
*     *	        D3dBs_fft_init		*
*     *					*
*     ***********************************

      subroutine D3dBs_fft_init(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz

#ifdef FFTW3
#include <fftw3.f>
       integer Atest(2),nxh,nxhz,nxhy
#endif

      logical value
      integer i,mthr, Parallel_maxthreads
      external Parallel_maxthreads

      mthr = Parallel_maxthreads()
      

      !call D3dB_nfft3d(nb,nfft3d)

      value = BA_alloc_get(mt_scpl,mthr*(2*nx(nb)+15),
     >        'sfttmpx',tmpx(2,nb),tmpx(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_scpl,mthr*(2*ny(nb)+15),
     >        'sfttmpy',tmpy(2,nb),tmpy(1,nb))
      value = value.and.
     >        BA_alloc_get(mt_scpl,mthr*(2*nz(nb)+15),
     >        'sfttmpz',tmpz(2,nb),tmpz(1,nb))
      if (.not. value) 
     >   call errquit('D3dBs_fft_init:out of heap memory',0, MA_ERR)

      call Parallel_shared_vector_szero(.true.,
     >                       (2*nx(nb)+15),scpl_mb(tmpx(1,nb)))

      do i=1,mthr
        !write(*,*) "DEBUG init fft arrays of thread ", i-1
        call srffti(nx(nb),scpl_mb(tmpx(1,nb)+(i-1)*(2*nx(nb)+15)))
        call scffti(ny(nb),scpl_mb(tmpy(1,nb)+(i-1)*(2*ny(nb)+15)))
        call scffti(nz(nb),scpl_mb(tmpz(1,nb)+(i-1)*(2*nz(nb)+15)))
      end do


#ifdef FFTW3


       iforward  = FFTW_FORWARD
       ibackward = FFTW_BACKWARD
c       iestimate = FFTW_PATIENT
c       iestimate = FFTW_MEASURE
c       iestimate = FFTW_ESTIMATE
c       iestimate = FFTW_EXHAUSTIVE
       call icopy(nplans*NBLOCKS,0,0,splans,1)
      if (mapping.eq.1) then
         nxh = (nx(nb)/2+1)
         nxhz = nxh*nz(nb)
         nxhy = nxh*ny(nb)
         if (.not.BA_alloc_get(mt_scpl,nx(nb)*ny(nb)*nq(nb),
     >                       'Atest',Atest(2),Atest(1)))
     >     call errquit('D3dBs_fft_init:out of heap memory',0,MA_ERR)

         call sfftw_plan_many_dft(splans(1,nb),1,nz(nb),nxh,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(2,nb),1,ny(nb),nxh,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1, 
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft_c2r(splans(3,nb),1,nx(nb),
     >        ny(nb)*nq(nb),
     >        scpl_mb(Atest(1)),nxh       *ny(nb)*nq(nb),1,nxh,
     >        scpl_mb(Atest(1)),(nx(nb)+2)*ny(nb)*nq(nb),1,nx(nb)+2, 
     >        FFTW_EXHAUSTIVE)

         call sfftw_plan_many_dft_r2c(splans(4,nb),1,nx(nb),
     >        ny(nb)*nq(nb),
     >        scpl_mb(Atest(1)),(nx(nb)+2)*ny(nb)*nq(nb),1,nx(nb)+2,
     >        scpl_mb(Atest(1)),nxh       *ny(nb)*nq(nb),1,nxh,
     >        FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(5,nb),1,ny(nb),nxh,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1, 
     >                        iforward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(6,nb),1,nz(nb),nxh,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)

         call sfftw_plan_many_dft(splans(7,nb),1,nz(nb),1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(8,nb),1,ny(nb),1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        ibackward,FFTW_EXHAUSTIVE)

          call sfftw_plan_many_dft(splans(9,nb),1,ny(nb),1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        scpl_mb(Atest(1)),nxhy,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)
          call sfftw_plan_many_dft(splans(10,nb),1,nz(nb),1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        scpl_mb(Atest(1)),nxhz,nxh,1,
     >                        iforward,FFTW_EXHAUSTIVE)


         if (.not.BA_free_heap(Atest(2)))
     >   call errquit('D3dBs_fft_init:freeing heap',0,MA_ERR)

      else

         nxh  = (nx(nb)/2+1)
         if (.not.BA_alloc_get(mt_scpl,nfft3d(nb),'Atest',
     >                         Atest(2),Atest(1)))
     >     call errquit('D3dBs_fft_init:out of heap memory',0,MA_ERR)

         call sfftw_plan_many_dft(splans(11,nb),1,nz(nb),nq3(nb),
     >                        scpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        scpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(12,nb),1,ny(nb),nq2(nb),
     >                        scpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        scpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft_c2r(splans(13,nb),1,nx(nb),
     >        nq1(nb),
     >        scpl_mb(Atest(1)),nfft3d(nb),1,nxh,
     >        scpl_mb(Atest(1)),n2ft3d(nb),1,(nx(nb)+2), 
     >        FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft_r2c(splans(14,nb),1,nx(nb),
     >        nq1(nb),
     >        scpl_mb(Atest(1)),n2ft3d(nb),1,nx(nb)+2,
     >        scpl_mb(Atest(1)),nfft3d(nb),1,nxh,
     >        FFTW_EXHAUSTIVE)

         call sfftw_plan_many_dft(splans(15,nb),1,ny(nb),nq2(nb),
     >                        scpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        scpl_mb(Atest(1)),ny(nb)*nq2(nb),1,ny(nb),
     >                        iforward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(16,nb),1,nz(nb),nq3(nb),
     >                        scpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        scpl_mb(Atest(1)),nz(nb)*nq3(nb),1,nz(nb),
     >                        iforward,FFTW_EXHAUSTIVE)

         call sfftw_plan_many_dft(splans(17,nb),1,nz(nb),1,
     >                        scpl_mb(Atest(1)),nz(nb),1,1,
     >                        scpl_mb(Atest(1)),nz(nb),1,1,
     >                        ibackward,FFTW_EXHAUSTIVE)
         call sfftw_plan_many_dft(splans(18,nb),1,ny(nb),1,
     >                        scpl_mb(Atest(1)),ny(nb),1,1,
     >                        scpl_mb(Atest(1)),ny(nb),1,1,
     >                        ibackward,FFTW_EXHAUSTIVE)

          call sfftw_plan_many_dft(splans(19,nb),1,ny(nb),1,
     >                        scpl_mb(Atest(1)),ny(nb),1,1,
     >                        scpl_mb(Atest(1)),ny(nb),1,1,
     >                        iforward,FFTW_EXHAUSTIVE)
          call sfftw_plan_many_dft(splans(20,nb),1,nz(nb),1,
     >                        scpl_mb(Atest(1)),nz(nb),1,1,
     >                        scpl_mb(Atest(1)),nz(nb),1,1,
     >                        iforward,FFTW_EXHAUSTIVE)

         if (.not.BA_free_heap(Atest(2)))
     >   call errquit('D3dBs_fft_init:freeing heap',0,MA_ERR)

      end if
#endif

      return
      end

*     ***********************************
*     *                                 *
*     *         D3dBs_fft_end           *
*     *                                 *
*     ***********************************

      subroutine D3dBs_fft_end(nb)
      implicit none
      integer nb

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz
#ifdef USE_OPENMP
      logical is_par
      common    / Debug_openmp / is_par
#endif

      logical value
      integer i

#ifdef FFTW3
      do i=1,nplans
         if (splans(i,nb).ne.0) call sfftw_destroy_plan(splans(i,nb))
      end do
#endif

      value =           BA_free_heap(tmpx(2,nb))
      value = value.and.BA_free_heap(tmpy(2,nb))
      value = value.and.BA_free_heap(tmpz(2,nb))
      if (.not.value)
     >   call errquit(
     >   'D3dBs_fft_end:error deallocatingof heap memory',0,MA_ERR)

      return
      end





*     ***********************************
*     *					*
*     *	        D3dBs_cr_fft3b		*
*     *					*
*     ***********************************

      subroutine D3dBs_cr_fft3b(nb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex       *
*      inverse fft                                  *
*           A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]   * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed and the imaginary   *
*              part of A is set to zero             *
*       uses - D3dB_c_transpose_jk, dcopy           *
*                                                   *
*****************************************************

      implicit none
      integer nb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz

*     *** local variables ***
      integer i,j,k,q,indx
      integer nxh,nxhy,nxhz,indx0,indx1

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value

      integer  tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_scpl,(nfft3d(nb)),'ffttmp2',
     >                    tmp2(2),tmp2(1))
      value = value.and.
     >      BA_push_get(mt_real,(n2ft3d(nb)),'ffttmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

       nxh = (nx(nb)/2+1)
       nxhz = nxh*nz(nb)
       nxhy = nxh*ny(nb)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,kz,ky) <- A(kx,ky,kz)      ***
*     ********************************************
c     call D3dBs_c_transpose_jk(nb,A,dcpl_mb(tmp2(1)),dbl_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************


#ifdef FFTW3
      do q=1,nq(nb)
        indx = 1+(q-1)*nxhz
        call sfftw_execute_dft(splans(1,nb),A(indx),A(indx))
      end do
#else
      indx0 = 0
      do q=1,nq(nb)
      do i=1,nxh

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            scpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call scfftb(nz(nb),scpl_mb(tmp2(1)),scpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(indx1) = dcpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do

      end do
      indx0 = indx0 + nxhz
      end do
#endif

*     ****************************************************
*     ***         Do a transpose of A                  ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ****************************************************
      call D3dBs_c_transpose_jk(nb,A,scpl_mb(tmp2(1)),real_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************


#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhy
         call sfftw_execute_dft(splans(2,nb),A(indx),A(indx))
      end do
#else
      indx0 = 0
      do q=1,nq(nb)
      do i=1,(nx(nb)/2+1)
      
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            scpl_mb(tmp2(1)+j-1) = A(indx)
            indx = indx + nxh
         end do
         call scfftb(ny(nb),scpl_mb(tmp2(1)),scpl_mb(tmpy(1,nb)))
         do j=1,ny(nb)
            A(indx1) = scpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do

      end do
      indx0 = indx0 + nxhy
      end do
#endif

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *********************************************************************

#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(3,nb),A,A)

#else

      call cshift1_sfftb(nx(nb),ny(nb),nq(nb),1,A)
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         !indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call srfftb(nx(nb),A(indx),scpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call zeroend_sfftb(nx(nb),ny(nb),nq(nb),1,A)
#endif


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


*     *****************************************************
*     ***     do fft along kz dimension                 ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *****************************************************

#ifdef FFTW3
      call sfftw_execute_dft(splans(11,nb),A,A)

#else

      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nq3(nb),nthr
        call scfftb(nz(nb),A(1+(i-1)*nz(nb)),scpl_mb(tmpz(1,nb)+offset))
      end do
!$OMP BARRIER

#endif

      call D3dBs_c_transpose_ijk(nb,3,A,scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************************
#ifdef FFTW3
      call sfftw_execute_dft(splans(12,nb),A,A)

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nq2(nb),nthr
        call scfftb(ny(nb),A(1+(i-1)*ny(nb)),scpl_mb(tmpy(1,nb)+offset))
      end do
!$OMP BARRIER
#endif

      call D3dBs_c_transpose_ijk(nb,4,A,scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *********************************************************************

#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(13,nb),A,A)

#else
      offset=tid*(2*nx(nb)+15)
      call cshift1_sfftb(nx(nb),nq1(nb),1,1,A)
      do i=tid+1,nq1(nb),nthr
        call srfftb(nx(nb),A(1+(i-1)*nxh),scpl_mb(tmpx(1,nb)+offset))
      end do
!$OMP BARRIER

      call zeroend_sfftb(nx(nb),nq1(nb),1,1,A)

#endif

      end if

    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = value.and.BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))
      if (.not. value) call errquit('popping stack memory',0,MA_ERR)

      call nwpw_timing_end(1)

      return
      end

      subroutine D3dBs_fftbx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real    tmpx(*)
      complex A(nxh,n)
      integer i

      do i=1,n
         call srfftb(nx,A(1,i),tmpx)
      end do
      return
      end

      subroutine D3dBs_fftby_sub2(n,ny,tmpy,A)
      implicit none
      integer n,ny
      real    tmpy(*)
      complex A(ny,n)
      integer i

      do i=1,n
         call scfftb(ny,A(1,i),tmpy)
      end do
      return
      end

      subroutine D3dBs_fftbz_sub2(n,nz,tmpz,A)
      implicit none
      integer n,nz
      real    tmpz(*)
      complex A(nz,n)
      integer i

      do i=1,n
         call scfftb(nz,A(1,i),tmpz)
      end do
      return
      end

      subroutine cshift1_sfftb(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real A(*)

      integer i,j,indx

!$OMP DO
      do j=1,(ny*nq*ne)
        indx = 1+(j-1)*(nx+2)
         do i=2,nx
            A(indx+i-1) = A(indx+i)
         end do
      end do
!$OMP END DO
         
      return
      end

      subroutine cshift1_sfftb1(nx,A)
      implicit none
      integer nx
      real A(*)
      integer i
      do i=2,nx
         A(i) = A(i+1)
      end do
      return
      end

      subroutine zeroend_sfftb1(nx,A)
      implicit none
      integer nx
      real A(*)
      integer i
      A(nx+1) = 0.0
      A(nx+2) = 0.0
      return
      end

      subroutine zeroend_sfftb(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real A(*)

      integer i,indx

!$OMP DO
      do i=1,(ny*nq*ne)
         indx = nx+1+(i-1)*(nx+2)
         A(indx)   = 0.0
         A(indx+1) = 0.0
      end do
!$OMP END DO

      return
      end




*     ***********************************
*     *					*
*     *	        D3dBs_rc_fft3f		*
*     *					*
*     ***********************************

      subroutine D3dBs_rc_fft3f(nb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional single precision complex *
*      to complex fft                               *
*   A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]    * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*       uses - transpose1 subroutine                *
*                                                   *
*****************************************************

      implicit none
      integer nb
      complex A(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,indx1
      integer nxh,nxhy,nxhz

      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value

      integer  tid,nthr,offset,nnn
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      call nwpw_timing_start(1)


*     ***** allocate temporary space ****
      !call D3dB_nfft3d(nb,nfft3d)
      nnn = nfft3d(nb)
      if ((nthr*ny(nb)).gt.nnn) nnn = nthr*ny(nb)
      if ((nthr*nz(nb)).gt.nnn) nnn = nthr*nz(nb)
      value = BA_push_get(mt_scpl,nnn,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        BA_push_get(mt_real,(n2ft3d(nb)),'tmp3',tmp3(2),tmp3(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      nxh = (nx(nb)/2+1)
      nxhz = nxh*nz(nb)
      nxhy = nxh*ny(nb)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ****************************************************************
*     ***     do fft along nx(nb) dimension                        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ****************************************************************
#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(4,nb),A,A)

#else
      offset=tid*(2*nx(nb)+15)
      do j=tid+1,ny(nb)*nq(nb),nthr
         call srfftf(nx(nb),A((j-1)*nxh+1),scpl_mb(tmpx(1,nb)+offset))
      end do
      call cshift_sfftf(nx(nb),ny(nb),nq(nb),1,A)
#endif


*     ********************************************************
*     ***     do fft along ny(nb) dimension                ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************************
#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhy
         call sfftw_execute_dft(splans(5,nb),A(indx),A(indx))
      end do

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nxh,nthr
      indx = i
      indx1= i
      do q=1,nq(nb)
         do j=1,ny(nb)
            !indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            scpl_mb(tmp2(1)+j-1+tid*ny(nb)) = A(indx)
            indx = indx + nxh
         end do

         call scfftf(ny(nb),scpl_mb(tmp2(1)+tid*ny(nb)),
     >                     scpl_mb(tmpy(1,nb)+offset))

         do j=1,ny(nb)
            !indx = i + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx1) = scpl_mb(tmp2(1)+j-1+tid*ny(nb))
            indx1 = indx1 + nxh
         end do
      end do
      end do

#endif

*     ****************************************************
*     ***         Do a transpose of A                  ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ****************************************************
      call D3dBs_c_transpose_jk(nb,A,scpl_mb(tmp2(1)),real_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
#ifdef FFTW3
      do q=1,nq(nb)
         indx = 1+(q-1)*nxhz
         call sfftw_execute_dft(splans(6,nb),A(indx),A(indx))
      end do

#else
      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nxh,nthr
      indx  = i
      indx1 = i
      do q=1,nq(nb)

         do k=1,nz(nb)
            scpl_mb(tmp2(1)+k-1+tid*nz(nb)) = A(indx)
            indx = indx + nxh
         end do
         call scfftf(nz(nb),scpl_mb(tmp2(1)+tid*nz(nb)),
     >                     scpl_mb(tmpz(1,nb)+offset))
         do k=1,nz(nb)
            !indx = i + (k-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
            A(indx1) = scpl_mb(tmp2(1)+k-1+tid*nz(nb))
            indx1 = indx1 + nxh
         end do

      end do
      end do

#endif

*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call D3dBs_c_transpose_jk(nb,A,scpl_mb(tmp2(1)),real_mb(tmp3(1)))



      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ****************************************************************
*     ***     do fft along nx(nb) dimension                        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ****************************************************************
#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(14,nb),A,A)
#else
      offset=tid*(2*nx(nb)+15)
      do i=tid+1,nq1(nb),nthr
        call srfftf(nx(nb),A(1+(i-1)*nxh),scpl_mb(tmpx(1,nb)+offset))
      end do
      call cshift_sfftf(nx(nb),nq1(nb),1,1,A)
#endif
      
      call D3dBs_c_transpose_ijk(nb,1,A,scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     ********************************************
*     ***     do fft along ny(nb) dimension        ***
*     ***   A(ky,nz(nb),kx) <- fft1d[A(ny(nb),nz(nb),kx)]  ***
*     ********************************************
#ifdef FFTW3
      call sfftw_execute_dft(splans(15,nb),A,A)

#else
      offset=tid*(2*ny(nb)+15)
      do i=tid+1,nq2(nb),nthr
        call scfftf(ny(nb),A(1+(i-1)*ny(nb)),scpl_mb(tmpy(1,nb)+offset))
      end do
!$OMP BARRIER
#endif
      call D3dBs_c_transpose_ijk(nb,2,A,scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))


*     ************************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(nb),kx,ky)]  ***
*     ************************************************
#ifdef FFTW3
      call sfftw_execute_dft(splans(16,nb),A,A)
#else
      offset=tid*(2*nz(nb)+15)
      do i=tid+1,nq3(nb),nthr
        call scfftf(nz(nb),A(1+(i-1)*nz(nb)),scpl_mb(tmpz(1,nb)+offset))
      end do
!$OMP BARRIER

#endif

      end if

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)

      return
      end

      subroutine D3dBs_fftfx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real    tmpx(2*nx+15)
      complex A(nxh,n)
      integer i

      do i=1,n
         call srfftf(nx,A(1,i),tmpx)
      end do
      return
      end

      subroutine  D3dBs_fftfy_sub2(n,ny,tmpy,A)
      implicit none
      integer n,ny
      real    tmpy(4*ny+15)
      complex A(ny,n)
      integer i,indx

      do i=1,n
         call scfftf(ny,A(1,i),tmpy)
      end do
      return
      end

      subroutine  D3dBs_fftfz_sub2(n,nz,tmpz,A)
      implicit none
      integer n,nz
      real    tmpz(4*nz+15)
      complex A(nz,n)
      integer i

      do i=1,n
         call scfftf(nz,A(1,i),tmpz)
      end do
      return
      end


      subroutine cshift_sfftf(nx,ny,nq,ne,A)
      implicit none
      integer nx,ny,nq,ne
      real A(*)

      integer i,j,indx

!$OMP BARRIER
!$OMP DO 
      do j=1,(ny*nq*ne)
        indx = 1+(j-1)*(nx+2)

         do i=nx,2,-1
            A(indx+i) = A(indx+i-1)
         end do
         A(indx+1)    = 0.0
         A(indx+nx+1) = 0.0
!         indx = indx + (nx+2)
      end do
!$OMP END DO

      return
      end 


      subroutine cshift_sfftf_ab(nx,ny,nq,ne,A,B)
      implicit none
      integer nx,ny,nq,ne
      real A(*)
      real B(*)

      integer i,j,indx

      indx = 1
      do j=1,(ny*nq*ne)
CDIR$ NOVECTOR
         do i=nx,2,-1
            B(indx+i) = A(indx+i-1)
         end do
         B(indx+1)    = 0.0
         B(indx+nx+1) = 0.0
         indx = indx + (nx+2)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	       D3dBs_(c,r,t)_SMul 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = scale * A
* where scale is a real*8 number.

      subroutine D3dBs_c_SMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      real       scale
      complex A(*)
      complex C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end


      subroutine D3dBs_c_SMul1(nb,scale,A)
      implicit none 
      integer    nb
      real    scale
      complex A(*)
      
#include "D3dB.fh"

      integer i
      
      do i=1,nfft3d_map(nb)
         A(i) = scale*A(i)
      end do 
      return 
      end



      subroutine D3dBs_r_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real scale
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = scale*A(i)
      end do
!$OMP END DO
      return
      end 

      subroutine D3dBs_r_SMul1(nb,scale,A)
      implicit none
      integer nb
      real scale
      real A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = scale*A(i)
      end do
!$OMP END DO
      return
      end


      subroutine D3dBs_t_SMul(nb,scale,A,C)     
      implicit none 
      integer nb
      real scale
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end 


      subroutine D3dBs_t_SMul1(nb,scale,A)
      implicit none
      integer nb
      real scale
      real A(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         A(i) = scale*A(i)
      end do
      return
      end



      subroutine D3dBs_c_ZMul(nb,scale,A,C)     
      implicit none 
      integer    nb
      complex scale
      complex A(*)
      complex C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = scale*A(i)
      end do
      return
      end

      subroutine D3dBs_r_Power1(nb,y,A)
      implicit none
      integer nb
      real  y
      real A(*)

#include "D3dB.fh"

      integer ii

      do ii=1,n2ft3d_map(nb)
         A(ii) = A(ii)**y
      end do
      
      return
      end


      subroutine D3dBs_rr_Power(nb,y,A,B)
      implicit none
      integer nb
      real  y
      real A(*)
      real B(*)

#include "D3dB.fh"

      integer i

      do i=1,n2ft3d_map(nb)
         B(i) = A(i)**y
      end do
      return
      end


*     ***********************************
*     *					*
*     *	       D3dBs_ct_Sqr	 	*	
*     *					*
*     ***********************************

*  This routine performs the operation	C = A * A

      subroutine D3dBs_ct_Sqr(nb,A,C)     
      implicit none 
      integer    nb
      complex A(*)
      real    C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = real(A(i))**2 + aimag(A(i))**2
      end do
      return
      end

*     ***********************************
*     *					*
*     *	       D3dBs_rr_Sqr	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_Sqr(nb,A,C)     
      implicit none 
      integer nb
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i)**2
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *                                 *
*     *        D3dBs_rr_SqrAdd          *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rr_SqrAdd(nb,A,C)
      implicit none
      integer nb
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = C(i) + A(i)**2
      end do
!$OMP END DO
      return
      end




*     ***********************************
*     *                                 *
*     *        D3dBs_rr_Sqr1            *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rr_Sqr1(nb,A)
      implicit none
      integer nb
      real A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = A(i)**2
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *					*
*     *	       D3dBs_rr_Sqrt	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_Sqrt(nb,A,C)     
      implicit none 
      integer nb
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

      do i=1,n2ft3d_map(nb)
         C(i) = sqrt(A(i))
      end do
      return
      end



*     ***********************************
*     *                                 *
*     *        D3dB_rr_Sqrt1            *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rr_Sqrt1(nb,A)
      implicit none
      integer nb
      real A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         A(i) = sqrt(A(i))
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *					*
*     *	       D3dBs_tt_Sqr	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_tt_Sqr(nb,A,C)     
      implicit none 
      integer nb
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = A(i)**2
      end do
      return
      end

             
*     ***********************************
*     *					*
*     *	         D3dBs_cc_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_dot(nb,A,B,sumall)     
      implicit none
      integer nb 
      complex A(*)
      complex B(*)
      real    sumall


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      real    sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0d0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ***** kx!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nx(nb)/2+1)
            index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >            + (j-1)*(nx(nb)/2+1) + i
            sum = sum + real(A(index)) * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
         end do
         end do
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nx(nb)/2+1)*ny(nb) + (j-1)*(nx(nb)/2+1) + 1
            sum = sum + real(A(index))  * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
         end do
      end do

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + real(A(index))  * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - real(A(index))  * real(B(index))
     >             - aimag(A(index)) * aimag(B(index))
         end if
      end do
      end do

      end if
      

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dBs_SumAll(sum)
      end if

      call nwpw_timing_end(2)

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_cc_idot 	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_idot(nb,A,B,sumall)     
      implicit none
      integer nb 
      complex A(*)
      complex B(*)
      real    sumall


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      real    sum


      call nwpw_timing_start(2)

c      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ***** kx!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nx(nb)/2+1)
            index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >            + (j-1)*(nx(nb)/2+1) + i
            sum = sum + real(A(index))  * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
         end do
         end do
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nx(nb)/2+1)*ny(nb) + (j-1)*(nx(nb)/2+1) + 1
            sum = sum + real(A(index))  * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
         end do
      end do
      

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + real(A(index))  * real(B(index))
     >                + aimag(A(index)) * aimag(B(index))
      end do
      sum = sum*2.0d0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - real(A(index))  * real(B(index))
     >             - aimag(A(index)) * aimag(B(index))
         end if
      end do
      end do

      end if


*     **** do not add up sums from other nodes ****
       
      call nwpw_timing_end(2)

      sumall = sum
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_tt_dot  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_tt_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real A(*)
      real B(*)
      real sumall

#include "D3dB.fh"

      integer i,j,k,q,index,np,nxh,taskid,p
      real    sum

      nxh=nx(nb)/2
      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ***** k!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
            sum = sum + A(index)*B(index)
         end do
         end do
      end do
      sum = sum*2.0

*     **** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
            sum = sum + A(index)*B(index)
         end do
      end do

      
      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)*B(index)
      end do
      sum = sum*2.0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)*B(index)
         end if
      end do
      end do

      end if
      

*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dBs_SumAll(sum)
      end if

      sumall = sum
      return
      end


*     ***********************************
*     *					*
*     *	         D3dBs_tt_idot 	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_tt_idot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real A(*)
      real B(*)
      real sumall


#include "D3dB.fh"

      integer i,j,k,q,index,np,nxh,taskid,p
      real    sum


      call nwpw_timing_start(2)

      nxh=nx(nb)/2
c      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ***** k!=0 plane, so double count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
            sum = sum + A(index)*B(index)
         end do
         end do
      end do
      sum = sum*2.0

*     **** kx==0 plane, so single count *****
      do q=1,nq(nb)
         do j=1,ny(nb)
            i=1
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
            sum = sum + A(index)*B(index)
         end do
      end do
      


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)*B(index)
      end do
      sum = sum*2.00

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(1,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)*B(index)
         end if
      end do
      end do

      end if


*     **** !!!! do not add up sums from other nodes ****

      call nwpw_timing_end(2)

      sumall = sum
      return
      end



*     ***********************************
*     *					*
*     *	         D3dBs_rr_dot  	 	*	
*     *					*
*     ***********************************
*     shared memory output
*     - sumall

      subroutine D3dBs_rr_dot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real  A(*)
      real  B(*)
      real  sumall

#include "D3dB.fh"

      integer i,np
      real    sum
      common /D3dBs_RR_TSUM/ sum

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
!$OMP MASTER
      sum = 0.0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO REDUCTION(+:sum)
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)*B(i)
      end do
!$OMP END DO


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dBs_SumAll(sum)
      end if

!$OMP MASTER
      sumall = sum
!$OMP END MASTER
!$OMP BARRIER
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_rr_idot   	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_idot(nb,A,B,sumall)     
      implicit none 
      integer nb
      real  A(*)
      real  B(*)
      real  sumall

#include "D3dB.fh"

      integer i,np
      real    sum
      common /D3dBs_RR_TSUM/ sum


*     **** sum up dot product on this node ****
!$OMP MASTER
      sum = 0.0
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO reduction(+:sum)
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)*B(i)
      end do
!$OMP END DO

*     **** add up sums from other nodes ****
*     if (np.gt.1) then
*        call D3dB_SumAll(sum)
*     end if

!$OMP MASTER
      sumall = sum
!$OMP END MASTER
!$OMP BARRIER
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_rrm_sym_dot      *       
*     *                                 *
*     ***********************************

      subroutine D3dBs_rrm_sym_dot(nb,n,A,B,matrix)
      implicit none
      integer nb,n
      real  A(*)
      real  B(*)
      real  matrix(n,n)

#include "D3dB.fh"

*     **** local variables ****
      integer j,k
      integer np

      call nwpw_timing_start(2)
      call Parallel2d_np_i(np)

      do k=1,n
        call SGEMM_OMP('T','N',k,1,n2ft3d(nb),
     >             1.00,
     >             A,n2ft3d(nb),
     >             B(1+(k-1)*n2ft3d(nb)),n2ft3d(nb),
     >             0.00,
     >             matrix(1,k),k)
      end do

!$OMP DO
      do k=1,n
      do j=k+1,n
        matrix(j,k) = matrix(k,j)
      end do
      end do
!$OMP END DO

      if (np.gt.1) call D3dBs_Vector_SumAll(n*n,matrix)
      call nwpw_timing_end(2)
      return
      end





*     ***********************************
*     *					*
*     *	         D3dBs_cc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_Mul(nb,A,B,C)     
      implicit none 
      integer    nb
      complex A(*)
      complex B(*)
      complex C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = conjg(A(i)) * B(i)
      end do

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_cc_Mul2          *
*     *                                 *
*     ***********************************

      subroutine D3dBs_cc_Mul2(nb,A,B)
      implicit none
      integer    nb
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         B(i) = A(i) * B(i)
      end do

      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_lc_Mask 	 	*
*     *					*
*     ***********************************

      subroutine D3dBs_lc_Mask(nb,masker,A)
      implicit none 
      integer nb
      logical masker(*)
      complex A(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         if (masker(i)) A(i) = cmplx(0.0,0.0)
      end do
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_lr_Mask  	*	
*     *					*
*     ***********************************

      subroutine D3dBs_lr_Mask(nb,masker,A)
      implicit none 
      integer nb
      logical masker(*)
      real    A(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         if (masker(i)) A(i) = 0.0
      end do
      return
      end


*     ***********************************
*     *					*
*     *	         D3dBs_tc_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_tc_Mul(nb,A,B,C)     
      implicit none 
      integer nb
      real    A(*)
      complex B(*)
      complex C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO
      return
      end



*     ***********************************
*     *                                 *
*     *          D3dBs_tc_Mul2          *
*     *                                 *
*     ***********************************

      subroutine D3dBs_tc_Mul2(nb,A,B)
      implicit none
      integer nb
      real    A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) * A(i)
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_rr_Mul  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_Mul(nb,A,B,C)     
      implicit none 
      integer nb
      real A(*)
      real B(*)
      real C(*)
      integer i,n

#include "D3dB.fh"

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) * B(i)
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Mul2          *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rr_Mul2(nb,A,B)
      implicit none
      integer nb
      real A(*)
      real B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) * A(i)
      end do
!$OMP END DO

      return
      end



*     ***********************************
*     *					*
*     *	         D3dBs_cc_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_Sum(nb,A,B,C)     
      implicit none 
      integer nb
      complex A(*)
      complex B(*)
      complex C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         C(i) = A(i) + B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          D3dBs_cc_Sum2          *
*     *                                 *
*     ***********************************
      subroutine D3dBs_cc_Sum2(nb,A,B)
      implicit none
      integer    nb
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + A(i)
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *					*
*     *	         D3dBs_rr_Sum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_Sum(nb,A,B,C)     
      implicit none 
      integer nb
      real  A(*)
      real  B(*)
      real  C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = B(i)+A(i)
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Sum2          *
*     *                                 *
*     ***********************************
      subroutine D3dBs_rr_Sum2(nb,A,B)
      implicit none
      integer nb
      real  A(*)
      real  B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) + A(i)
      end do
!$OMP END DO
      return
      end


*     ***********************************
*     *					*
*     *	         D3dBs_cc_Sub  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_Sub(nb,A,B,C)     
      implicit none 
      integer nb
      complex A(*)
      complex B(*)
      complex C(*)

#include "D3dB.fh"

      integer i

      do i=1,nfft3d_map(nb)
         C(i) = A(i) - B(i)
      end do

      return
      end


*     ***********************************
*     *					*
*     *	         D3dBs_rr_Sub   	*	
*     *					*
*     ***********************************

      subroutine D3dBs_rr_Sub(nb,A,B,C)     
      implicit none 

#include "D3dB.fh"

      integer nb
      real  A(n2ft3d_map(nb))
      real  B(n2ft3d_map(nb))
      real  C(n2ft3d_map(nb))

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) - B(i)
      end do
!$OMP END DO
      return
      end



*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Sub2          *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rr_Sub2(nb,A,B)
      implicit none

#include "D3dB.fh"

      integer nb
      real    A(n2ft3d_map(nb))
      real    B(n2ft3d_map(nb))

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) - A(i)
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *         D3dBs_rr_Multiply2      *       
*     *                                 *
*     ***********************************
      subroutine D3dBs_rr_Multiply2(nb,A,B)
      implicit none 
      integer nb
      real  A(*)
      real  B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i)*A(i)
      end do
!$OMP END DO
      return
      end



*     ***********************************
*     *                                 *
*     *       D3dBs_rrr_MultiplyAdd     *       
*     *                                 *
*     ***********************************
      subroutine D3dBs_rrr_MultiplyAdd(nb,A,B,C)
      implicit none
      integer nb
      real  A(*)
      real  B(*)
      real  C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = C(i) + B(i)*A(i)
      end do
!$OMP END DO
      return
      end



*     ***********************************
*     *					*
*     *	        D3dBs_cc_zaxpy 	 	*	
*     *					*
*     ***********************************
      subroutine D3dBs_cc_zaxpy(nb,alpha,A,B)     
      implicit none 
      integer nb
      complex alpha
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + alpha*A(i)
      end do
!$OMP END DO
      return
      end



*     ***********************************
*     *					*
*     *	         D3dBs_cc_daxpy	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_daxpy(nb,alpha,A,B)     
      implicit none 
      integer    nb
      real    alpha
      complex A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = B(i) + alpha*A(i)
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_rr_daxpy  	*	
*     *					*
*     ***********************************
      subroutine D3dBs_rr_daxpy(nb,alpha,A,B)     
      implicit none 
      integer nb
      real  alpha
      real  A(*)
      real  B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         B(i) = B(i) + alpha* A(i)
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Divide        *
*     *                                 *
*     ***********************************
      subroutine D3dBs_rr_Divide(nb,A,B,C)
      implicit none
      integer nb
      real A(*)
      real B(*)
      real C(*)

#include "D3dB.fh"

      real eta
      parameter (eta=1.0e-9)

      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
CDIR$ NOVECTOR
!$OMP DO
      do index = 1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         if (abs(B(index)) .le. eta) then
           C(index) = 0.00
         else
           C(index) = A(index) / B(index)
         end if
      end do
!$OMP END DO
      !end do
      !end do
      !end do
      return
      end



*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Divide2       *
*     *                                 *
*     ***********************************
      subroutine D3dBs_rr_Divide2(nb,A,B)
      implicit none
      integer nb
      real A(*)
      real B(*)

#include "D3dB.fh"

      real eta
      parameter (eta=1.0e-9)

      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
CDIR$ NOVECTOR
!$OMP DO
      do index = 1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         if (abs(A(index)) .le. eta) then
           B(index) = 0.0
         else
           B(index) = B(index) / A(index)
         end if
      end do
!$OMP END DO
      !end do
      !end do
      !end do

      return
      end



*     ***********************************
*     *                                 *
*     *          D3dBs_r_ABS            *
*     *                                 *
*     ***********************************
      subroutine D3dBs_r_ABS(nb,A,C)
      implicit none
      integer nb
      real A(*)
      real C(*)

#include "D3dB.fh"

      integer index

      !do q=1,nq(nb)
      !do j=1,ny(nb)
      !do i=1,nx(nb)
!$OMP DO
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         C(index) = abs(A(index))
      end do
!$OMP END DO
      !end do
      !end do
      !end do
      return
      end

      subroutine D3dBs_r_abs1(nb,A)
      implicit none
      integer nb
      real A(*)

#include "D3dB.fh"

      integer index

!$OMP DO
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         A(index) = abs(A(index))
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_r_ABSMAX         *
*     *                                 *
*     ***********************************
      subroutine D3dBs_r_ABSMAX(nb,A,amax)
      implicit none
      integer nb
      real A(*)
      real amax

#include "D3dB.fh"

      integer index
      real aa

      amax = 0.0
      do index=1,n2ft3d_map(nb)
         !index = i + (j-1)*(nx(nb)+2) + (q-1)*(nx(nb)+2)*ny(nb)
         aa = abs(A(index))
         if (aa.gt.amax) amax = aa
      end do
      call D3dBs_MaxAll(amax)

      return
      end


*     ***********************************
*     *                                 *
*     *         D3dBs_r_ZeroNegative    *
*     *                                 *
*     ***********************************

      subroutine D3dBs_r_ZeroNegative(nb,A)
      implicit none
      integer nb
      real A(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         if (A(i).lt.0.0d0) A(i) = 0.0
      end do
!$OMP END DO
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_rr_Minus         *
*     *                                 *
*     ***********************************
      subroutine D3dBs_rr_Minus(nb,A,B,C)
      implicit none
      integer nb
      real  A(*)
      real  B(*)
      real  C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,n2ft3d_map(nb)
         C(i) = A(i) - B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *					*
*     *          D3dBs_r_Zero_Ends0	*
*     *					*
*     ***********************************
      subroutine D3dBs_r_Zero_Ends0(nb,A)
      integer nb
      real A(*) 

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2) 
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))

            A(index)   = 0.0
            A(index+1) = 0.0
         end do
      end do


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0

        end do
        end do
      end if

      if (n2ft3d_map(nb).lt.n2ft3d(nb)) then
         call scopy((n2ft3d(nb)-n2ft3d_map(nb)),
     >              0.0,0,A(n2ft3d_map(nb)+1),1)
      end if

      return
      end





*     ***********************************
*     *					*
*     *          D3dBs_r_Zero_Ends	*
*     *					*
*     ***********************************

      subroutine D3dBs_r_Zero_Ends(nb,A)
      integer nb
      real A(*) 

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
!$OMP DO
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2) 
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))

            A(index)   = 0.0
            A(index+1) = 0.0
         end do
      end do
!$OMP END DO


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
!$OMP DO
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 0.0

        end do
        end do
!$OMP end DO
      end if

!$OMP MASTER
      if (n2ft3d_map(nb).lt.n2ft3d(nb)) then
         call scopy((n2ft3d(nb)-n2ft3d_map(nb)),
     >              0.00,0,A(n2ft3d_map(nb)+1),1)
      end if
!$OMP END MASTER
!$OMP BARRIER

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dBs_r_notZero_Ends    *
*     *                                 *
*     ***********************************
      subroutine D3dBs_r_notZero_Ends(nb,A)
      integer nb
      real A(*) 

#include "D3dB.fh"

      integer j,k,q,index,taskid,p

      !**** slab mapping ****
      if (mapping.eq.1) then
      do q=1,nq(nb)
         do j=1,ny(nb)
            index = (nx(nb)+1) + (j-1)*(nx(nb)+2)
     >                         + (q-1)*(nx(nb)+2)*(ny(nb))
            A(index)   = 1.0
            A(index+1) = 1.0
         end do
      end do


      !**** hilbert mapping ****
      else
        call Parallel2d_taskid_i(taskid)
        do k=1,nz(nb)
        do j=1,ny(nb)

         call D3dB_ijktoindex2p(nb,(nx(nb)+1),j,k,index,p)
         if (p.eq.taskid) A(index) = 1.0

         call D3dB_ijktoindex2p(nb,(nx(nb)+2),j,k,index,p)
         if (p.eq.taskid) A(index) = 1.0

        end do
        end do
      end if

      return
      end




*     ***********************************
*     *					*
*     *	         D3dBs_r_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_r_dsum(nb,A,sumall)     
      implicit none
      integer nb 
      real  A(*)
      real  sumall

#include "D3dB.fh"

      integer i,np
      real    sum

      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0
      do i=1,n2ft3d_map(nb)
         sum = sum + A(i)
      end do

*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call D3dB_SumAll(sum)
      end if

      sumall = sum

      return
      end

*     ***********************************
*     *					*
*     *	         D3dBs_t_dsum  	 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_t_dsum(nb,A,sumall)     
      implicit none 
      integer nb
      real  A(*)
      real  sumall

#include "D3dB.fh"

      integer i,j,k,q,np,nxh,index,taskid,p
      real sum

      nxh = nx(nb)/2
      call Parallel2d_np_i(np)

*     **** sum up dot product on this node ****
      sum = 0.0

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     ***** k!=0 plane so double count *****
      do q=1,nq(nb)
      do j=1,ny(nb)
         do i=2,(nxh+1)
            index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + i
            sum = sum + A(index)
         end do
      end do
      end do
      sum = sum*2.0

*     ***** k==0 plane, so single count *****
      do q=1,nq(nb)
      do j=1,ny(nb)
         index = (q-1)*(nxh+1)*ny(nb) + (j-1)*(nxh+1) + 1
         sum = sum + A(index)
      end do
      end do


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     ***** kx!=0 plane, so double count *****
      do index=1,nfft3d_map(nb)
            sum = sum + A(index)
      end do
      sum = sum*2.0

*     ***** kx==0 plane, so single count *****
      do k=1,nz(nb)
      do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         sum = sum - A(index)
         end if
      end do
      end do

      end if


*     **** add up sums from other nodes ****
      if (np.gt.1) then
        call D3dBs_SumAll(sum)
      end if

      sumall = sum

      return
      end


*     ***********************************
*     *					*
*     *	     D3dBs_cc_Vector_dot 	*	
*     *					*
*     ***********************************
      subroutine D3dBs_cc_Vector_dot(nb,nnfft3d,nn,ne,A,B,sumall) 
      implicit none 
      integer    nb
      integer    nnfft3d,nn,ne
      complex A(*)
      complex B(*)
      real    sumall(nn,nn)


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      integer index1,index2
      integer n,m,shift1,shift2
      real    sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = (n-1)*nnfft3d
        shift2 = (m-1)*nnfft3d
        sum    = 0.0

*       ***** kx!=0 plane, so double count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
           do i=2,(nx(nb)/2+1)
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + i
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + real(A(index1))  * real(B(index2))
     >                  + aimag(A(index1)) * aimag(B(index2))
           end do
           end do
        end do
        sum = sum*2.0

*       ***** kx==0 plane, so single count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
              i=1
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + 1
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + real(A(index1))  * real(B(index2))
     >                  + aimag(A(index1)) * aimag(B(index2))
           end do
        end do

        sumall(n,m) = sum
        sumall(m,n) = sum
      end do
      end do


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     **** sum up dot product on this node ****
      do n=1,ne
      do m=n,ne

        shift1 = (n-1)*nnfft3d
        shift2 = (m-1)*nnfft3d
        sum    = 0.0

*       ***** kx!=0 plane, so double count *****
        do index=1,nfft3d_map(nb)
            index1 = index+shift1
            index2 = index+shift2
            sum = sum + real(A(index1))  * real(B(index2))
     >                + aimag(A(index1)) * aimag(B(index2))
        end do
        sum = sum*2.0

*       ***** kx==0 plane, so single count *****
        do k=1,nz(nb)
        do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         index1 = index+shift1
         index2 = index+shift2
         sum = sum - real(A(index1))  * real(B(index2))
     >             - aimag(A(index1)) * aimag(B(index2))
         end if
        end do
        end do

        sumall(n,m) = sum
        sumall(m,n) = sum
      end do
      end do

      end if


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dBs_Vector_SumAll(nn*ne,sumall)
      end if

      call nwpw_timing_end(2)

      return
      end



*     ***********************************
*     *					*
*     *	     D3dBs_cc_Vector_ndot 	*	
*     *					*
*     ***********************************

      subroutine D3dBs_cc_Vector_ndot(nb,nnfft3d,ne,A,B,sumall)
      implicit none 
      integer    nb
      integer    nnfft3d,ne
      complex A(*)
      complex B(*)
      real    sumall(ne)


#include "D3dB.fh"

      integer i,j,k,q,index,np,taskid,p
      integer index1,index2
      integer n,shift1,shift2
      real    sum


      call nwpw_timing_start(2)

      call Parallel2d_np_i(np)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = (n-1)*nnfft3d
        shift2 = 0
        sum    = 0.0
*       ***** kx!=0 plane, so double count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
           do i=2,(nx(nb)/2+1)
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + i
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + real(A(index1))  * real(B(index2))
     >                  + aimag(A(index1)) * aimag(B(index2))
           end do
           end do
        end do
        sum = sum*2.0

*       ***** kx==0 plane, so single count *****
        do q=1,nq(nb)
           do j=1,ny(nb)
              i=1
              index = (q-1)*(nx(nb)/2+1)*ny(nb) 
     >              + (j-1)*(nx(nb)/2+1) + 1
              index1 = index+shift1
              index2 = index+shift2
              sum = sum + real(A(index1))  * real(B(index2))
     >                  + aimag(A(index1)) * aimag(B(index2))
           end do
        end do

        sumall(n) = sum
      end do


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call Parallel2d_taskid_i(taskid)
*     **** sum up dot product on this node ****
      do n=1,ne

        shift1 = (n-1)*nnfft3d
        shift2 = 0
        sum    = 0.0

*       ***** kx!=0 plane, so double count *****
        do index=1,nfft3d_map(nb)
            index1 = index+shift1
            index2 = index+shift2
            sum = sum + real(A(index1))  * real(B(index2))
     >                + aimag(A(index1)) * aimag(B(index2))
        end do
        sum = sum*2.0

*       ***** kx==0 plane, so single count *****
        do k=1,nz(nb)
        do j=1,ny(nb)
         i=1
         call D3dB_ijktoindexp(nb,i,j,k,index,p)
         if (p.eq.taskid) then
         index1 = index+shift1
         index2 = index+shift2
         sum = sum - real(A(index1))  * real(B(index2))
     >             - aimag(A(index1)) * aimag(B(index2))
         end if
        end do
        end do

        sumall(n) = sum
      end do

      end if


*     **** add up sums from other nodes ****
      if (np.gt.1) then
         call D3dBs_Vector_SumAll(ne,sumall)
      end if

      call nwpw_timing_end(2)

      return
      end






*     ***********************************
*     *                                 *
*     *          D3dBs_ic_Mul            *
*     *                                 *
*     ***********************************
cpgi$r opt=1
      subroutine D3dBs_ic_Mul(nb,A,B,C)
      implicit none
      integer    nb
      real    A(*)
      complex B(*)
      complex C(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
            C(i) = cmplx(0.0,A(i)) * B(i)
      end do
!$OMP END DO

      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_ic_Mul2          *
*     *                                 *
*     ***********************************
cpgi$r opt=1
      subroutine D3dBs_ic_Mul2(nb,A,B)
      implicit none
      integer    nb
      real    A(*)
      complex B(*)

#include "D3dB.fh"

      integer i

!$OMP DO
      do i=1,nfft3d_map(nb)
         B(i) = cmplx(0.0,A(i)) * B(i)
      end do
!$OMP END DO

      return
      end


*     ***********************************
*     *                                 *
*     *          D3dBs_r_Expand         *
*     *                                 *
*     ***********************************

      subroutine D3dBs_r_Expand(nb1,A,nb2,B)
      implicit none
      integer nb1
      real    A(*)
      integer nb2
      real    B(*)

#include "D3dB.fh"

      integer i,j,q,index1,index2

      if (mapping.eq.1) then
      call Parallel_shared_vector_szero(.true.,
     >                                 nq(nb2)*ny(nb2)*(nx(nb2)+2),B)
!$OMP DO
      do j=1,ny(nb1)
      do q=1,nq(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         B(index2) = A(index1)
      end do
      end do
      end do
!$OMP END DO

      else
      call Parallel_shared_vector_szero(.true.,n2ft3d(nb2),B)
!$OMP DO
      do q=1,nq1(nb1)
      do i=1,nx(nb1)
         index1 = i + (q-1)*(nx(nb1)+2)
         index2 = i + (q-1)*(nx(nb2)+2)
         B(index2) = A(index1)
      end do
      end do
!$OMP END DO
      end if
      return
      end

*     ***********************************
*     *                                 *
*     *          D3dBs_r_Contract       *
*     *                                 *
*     ***********************************
      subroutine D3dBs_r_Contract(nb2,B,nb1,A)
      implicit none
      integer nb2
      real    B(*)
      integer nb1
      real    A(*)

#include "D3dB.fh"

      integer i,j,q,index1,index2

      if (mapping.eq.1) then
      call Parallel_shared_vector_szero(.true.,
     >                                 nq(nb1)*ny(nb1)*(nx(nb1)+2),A)
!$OMP DO
      do j=1,ny(nb1)
      do q=1,nq(nb1)
      do i=1,nx(nb1)
         index1 = i + (j-1)*(nx(nb1)+2) + (q-1)*(nx(nb1)+2)*ny(nb1)
         index2 = i + (j-1)*(nx(nb2)+2) + (q-1)*(nx(nb2)+2)*ny(nb2)
         A(index1) = B(index2)
      end do
      end do
      end do
!$OMP END DO

      else
      call Parallel_shared_vector_szero(.true.,n2ft3d(nb1),A)
!$OMP DO
      do q=1,nq1(nb1)
      do i=1,nx(nb1)
         index1 = i + (q-1)*(nx(nb1)+2)
         index2 = i + (q-1)*(nx(nb2)+2)
         A(index1) = B(index2)
      end do
      end do
!$OMP END DO
      end if
      return
      end


      subroutine D3dBs_pfft_index1_copy(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex a(*),b(*)
      integer i

#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        b(i) = a(index(i))
      end do
!$OMP END DO
      return
      end

      subroutine D3dBs_pfft_index2_copy(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex  a(*),b(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        b(index(i)) = a(i)
      end do
!$OMP END DO NOWAIT
      return
      end


      subroutine D3dBs_pfft_index2_zero(n,index,a)
      implicit none
      integer n
      integer index(*)
      complex a(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        a(index(i)) = 0.0
      end do
!$OMP END DO
      return
      end



      subroutine D3dBs_pfft_index2_copy_conjg(n,index,a,b)
      implicit none
      integer n
      integer index(*)
      complex  a(*),b(*)
      integer i
#ifndef CRAY
!DIR$ ivdep
#endif
!$OMP DO
      do i=1,n
        b(index(i)) = conjg(a(i))
      end do
!$OMP END DO NOWAIT
      return
      end

