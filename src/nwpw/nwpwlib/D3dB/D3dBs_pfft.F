
#define NBLOCKS 3



*     ***********************************
*     *					*
*     *	        D3dBs_cr_pfft3b		*
*     *					*
*     ***********************************

      subroutine D3dBs_cr_pfft3b(nb,nbb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional single precision complex *
*      to complex inverse fft                       *
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
*       uses - D3dBs_c_ptranspose_jk, scopy         *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,ierr
      integer nxh,nxhy,nxhz,indx0,indx1,indx2

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(1)
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

*     ***** allocate temporary space ****
      !call D3dBs_nfft3d(nb,nfft3d)
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
c     call D3dBs_c_stranspose_jk(nb,A,scpl_mb(tmp2(1)),real_mb(tmp3(1)))

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *************************************************

      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call sfftw_execute_dft(splans(7,nb),A(indx),A(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            scpl_mb(tmp2(1)+k-1) = A(indx)
            indx = indx + nxh
         end do
         call scfftb(nz(nb),scpl_mb(tmp2(1)),scpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(indx1) = scpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
#endif

       end if
      end do
      indx0 = indx0 + nxhz
      end do

*     ****************************************************
*     ***         Do a ptranspose of A                 ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ****************************************************
      call D3dBs_c_ptranspose1_jk(nbb,A,scpl_mb(tmp2(1)),
     >                                  real_mb(tmp3(1)))

*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************************
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call sfftw_execute_dft(splans(8,nb),A(indx),A(indx))
#else
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
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *************************************************

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


*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(17,nb),A(indx),A(indx))
         end if
         indx = indx + nz(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(nb)+15)
      do q=tid+1,nq3(nb),nthr
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call scfftb(nz(nb),A(1+(q-1)*nz(nb)),
     >                 scpl_mb(tmpz(1,nb)+offset))
         end if
         !indx = indx + nz(nb)
      end do
!$OMP BARRIER

#endif

      call D3dBs_c_ptranspose_ijk(nbb,3,A,
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(18,nb),A(indx),A(indx))
         end if
         indx = indx + ny(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(nb)+15)
      do q=tid+1,nq2(nb),nthr
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call scfftb(ny(nb),A(1+(q-1)*ny(nb)),
     >                  scpl_mb(tmpy(1,nb)+offset))
         end if
         !indx = indx + ny(nb)
      end do
!$OMP BARRIER

#endif

      call D3dBs_c_ptranspose_ijk(nbb,4,A,
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *********************************************************************

#ifdef FFTW3
         call sfftw_execute_dft_c2r(splans(13,nb),A,A)
#else

      !indx = 1
      offset = tid*(2*nx(nb)+15)
      call cshift1_sfftb(nx(nb),nq1(nb),1,1,A)
      do q=tid+1,nq1(nb),nthr
         call srfftb(nx(nb),A(1+(q-1)*nxh),
     >              scpl_mb(tmpx(1,nb)+offset))
         !indx = indx + nxh
      end do
!$OMP BARRIER
      call zeroend_sfftb(nx(nb),nq1(nb),1,1,A)

#endif

      end if
    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end




*     ***********************************
*     *					*
*     *	        D3dBs_rc_pfft3f		*
*     *					*
*     ***********************************

      subroutine D3dBs_rc_pfft3f(nb,nbb,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional a single precision       *
*     complex to complex fft                        *
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
*       uses - ptranspose1 subroutine               *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb
      complex A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,k,q,indx,indx1,ierr,indx2,indx0
      integer nxh,nxhy,nxhz

      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(1)

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

*     ***** allocate temporary space ****
      !call D3dBs_nfft3d(nb,nfft3d)
      value = BA_push_get(mt_scpl,(nfft3d(nb)),'tmp2',tmp2(2),tmp2(1))
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
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         call srfftf(nx(nb),A(indx),scpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call cshift_sfftf(nx(nb),ny(nb),nq(nb),1,A)
#endif


*     ********************************************************
*     ***     do fft along ny(nb) dimension                ***
*     ***   A(kx,ky,nz(nb)) <- fft1d[A(kx,ny(nb),nz(nb))]  ***
*     ********************************************************

       indx0 = 0
       indx2 = 0
       do q=1,nq(nb)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
          indx  = i + indx0
          call sfftw_execute_dft(splans(9,nb),A(indx),A(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,ny(nb)
             scpl_mb(tmp2(1)+k-1) = A(indx)
             indx = indx + nxh
          end do
          call scfftf(ny(nb),scpl_mb(tmp2(1)),scpl_mb(tmpy(1,nb)))
          do k=1,ny(nb)
             A(indx1) = scpl_mb(tmp2(1)+k-1)
             indx1 = indx1 + nxh
          end do
#endif
        end if
       end do
       indx0 = indx0 + nxhy
       end do



*     ****************************************************
*     ***         Do a ptranspose of A                 ***
*     ***      A(ky,nz(nb),ky) <- A(kx,ky,nz(nb))      ***
*     ****************************************************
      call D3dBs_c_ptranspose2_jk(nbb,A,scpl_mb(tmp2(1)),
     >                                  real_mb(tmp3(1)))


*     ********************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(nb),ky)]  ***
*     ********************************************
      
       indx0 = 0
       indx2 = 0
       do q=1,nq(nb)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then
 
#ifdef FFTW3
          indx  = i + indx0
          call sfftw_execute_dft(splans(10,nb),A(indx),A(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,nz(nb)
             scpl_mb(tmp2(1)+k-1) = A(indx)
             indx = indx + nxh
          end do
          call scfftf(nz(nb),scpl_mb(tmp2(1)),scpl_mb(tmpz(1,nb)))
          do k=1,nz(nb)
             A(indx1) = scpl_mb(tmp2(1)+k-1)
             indx1 = indx1 + nxh
          end do
#endif
        end if
       end do
       indx0 = indx0 + nxhz
       end do


*     ********************************************
*     ***         Do a transpose of A          ***
*     ***      A(kx,ky,kz) <- A(kx,kz,ky)      ***
*     ********************************************
c     call D3dBs_c_transpose_jk(nb,A,scpl_mb(tmp2(1)),real_mb(tmp3(1)))



      !*************************
      !**** hilbert mapping ****
      !*************************
      else

*     ********************************************
*     ***     do fft along nx(nb) dimension        ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d[A(nx(nb),ny(nb),nz(nb))]  ***
*     ********************************************

#ifdef FFTW3
          call sfftw_execute_dft_r2c(splans(14,nb),A,A)

#else
      !indx = 1
      offset = tid*(2*nx(nb)+15)
      do q=tid+1,nq1(nb),nthr
         call srfftf(nx(nb),A(1+(q-1)*nxh),scpl_mb(tmpx(1,nb)+offset))
         !indx = indx + nxh
      end do
      call cshift_sfftf(nx(nb),nq1(nb),1,1,A)

#endif

      call D3dBs_c_ptranspose_ijk(nbb,1,A,
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     ********************************************************
*     ***     do fft along ny(nb) dimension                ***
*     ***   A(ky,nz(nb),kx) <- fft1d[A(ny(nb),nz(nb),kx)]  ***
*     ********************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(19,nb),A(indx),A(indx))
         end if
         indx = indx + ny(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(nb)+15)
      do q=tid+1,nq2(nb),nthr
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call scfftf(ny(nb),A(1+(q-1)*ny(nb)),
     >                 scpl_mb(tmpy(1,nb)+offset))
         end if
         !indx = indx + ny(nb)
      end do
!$OMP BARRIER

#endif

      call D3dBs_c_ptranspose_ijk(nbb,2,A,
     >                            scpl_mb(tmp2(1)),
     >                            real_mb(tmp3(1)))

*     ************************************************
*     ***     do fft along nz(nb) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(nb),kx,ky)]  ***
*     ************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(20,nb),A(indx),A(indx))
         end if
         indx = indx + nz(nb)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(nb)+15)
      do q=tid+1,nq3(nb),nthr
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call scfftf(nz(nb),A(1+(q-1)*nz(nb)),
     >                 scpl_mb(tmpz(1,nb)+offset))
         end if
         !indx = indx + nz(nb)
      end do
!$OMP BARRIER

#endif

      end if

*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))

      call nwpw_timing_end(1)
      return
      end






*     ***********************************
*     *					*
*     *	        D3dBs_cr_mpfft3b	*
*     *					*
*     ***********************************

      subroutine D3dBs_cr_mpfft3b(nb,nbb,m,A)

*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional single precision complex *
*      to complex inverse fft                       *
*    A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]  * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed and the imaginary   *
*              part of A is set to zero             *
*       uses - D3dBs_c_ptranspose_jk, dcopy          *
*                                                   *
*****************************************************

      implicit none
      integer nb,nbb,m
      complex  A(*)

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer i,j,k,q,indx,ierr,ii,shift,shift0
      integer nxh,nxhy,nxhz,indx0,indx1,indx2

      
      !integer tmp1(2),tmp2(2),tmp3(2)
      integer tmp2(2),tmp3(2)
      logical value


      call nwpw_timing_start(1)

*     ***** allocate temporary space ****
      !call D3dBs_nfft3d(nb,nfft3d)
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

      shift = 1
      do ii=1,m
      shift0=shift-1
*     *****************************************************
*     ***     do fft along kz dimension                 ***
*     ***   A(kx,nz(nb),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *****************************************************

      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx = i+indx0+shift0
         call sfftw_execute_dft(splans(7,nb),A(indx),A(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do k=1,nz(nb)
            scpl_mb(tmp2(1)+k-1) = A(shift0+indx)
            indx = indx + nxh
         end do
         call scfftb(nz(nb),scpl_mb(tmp2(1)),scpl_mb(tmpz(1,nb)))
         do k=1,nz(nb)
            A(shift0+indx1) = scpl_mb(tmp2(1)+k-1)
            indx1 = indx1 + nxh
         end do
#endif

       end if
      end do
      indx0 = indx0 + nxhz
      end do


*     ****************************************************
*     ***         Do a ptranspose of A                 ***
*     ***      A(kx,ky,nz(nb)) <- A(kx,nz(nb),ky)      ***
*     ****************************************************
      call D3dBs_c_ptranspose1_jk(nbb,A(shift),
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))
      shift = shift + nfft3d(nb)
      end do

      shift = 1
      do ii=1,m
      shift0 = shift-1
*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(kx,ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ky,nz(nb))]  ***
*     *************************************************************
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(nb)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx = i+indx0+shift0
         call sfftw_execute_dft(splans(8,nb),A(indx),A(indx))
#else
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(nb)
            scpl_mb(tmp2(1)+j-1) = A(shift0+indx)
            indx = indx + nxh
         end do
         call scfftb(ny(nb),scpl_mb(tmp2(1)),scpl_mb(tmpy(1,nb)))
         do j=1,ny(nb)
            A(shift0+indx1) = scpl_mb(tmp2(1)+j-1)
            indx1 = indx1 + nxh
         end do
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *********************************************************************

#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(3,nb),A(shift0),A(shift0))

#else
      call cshift1_sfftb(nx(nb),ny(nb),nq(nb),1,A(shift))
      indx = 1
      do q=1,nq(nb)
      do j=1,ny(nb)
         !indx = 1 + (j-1)*(nx(nb)/2+1) + (q-1)*(nx(nb)/2+1)*ny(nb)
         call srfftb(nx(nb),A(shift0+indx),scpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      end do
      call zeroend_sfftb(nx(nb),ny(nb),nq(nb),1,A(shift))
#endif


      shift = shift + nfft3d(nb)
      end do

      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      shift = 1
      do ii=1,m
      shift0=shift-1

*     *****************************************************
*     ***     do fft along kz dimension                 ***
*     ***   A(nz(nb),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *****************************************************

      indx = 1
      do q=1,nq3(nb)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then

#ifdef FFTW3
         call sfftw_execute_dft(splans(17,nb),
     >                         A(indx+shift0),A(indx+shift0))
#else
         call scfftb(nz(nb),A(shift0+indx),scpl_mb(tmpz(1,nb)))
#endif

         end if
         indx = indx + nz(nb)
      end do

      call D3dBs_c_ptranspose_ijk(nbb,3,A(shift),
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(ny(nb),nz(nb),kx) <- fft1d^(-1)[A(ky,nz(nb),kx)]  ***
*     *************************************************************

      indx = 1
      do q=1,nq2(nb)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then

#ifdef FFTW3
         call sfftw_execute_dft(splans(18,nb),
     >                          A(indx+shift0),A(indx+shift0))
#else
         call scfftb(ny(nb),A(shift0+indx),scpl_mb(tmpy(1,nb)))
#endif

         end if
         indx = indx + ny(nb)
      end do

      call D3dBs_c_ptranspose_ijk(nbb,4,A(shift),
     >                           scpl_mb(tmp2(1)),
     >                           real_mb(tmp3(1)))

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(nb),ny(nb),nz(nb)) <- fft1d^(-1)[A(kx,ny(nb),nz(nb))]  ***
*     *********************************************************************

#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(13,nb),A(shift0),A(shift0))

#else
      call cshift1_sfftb(nx(nb),nq1(nb),1,1,A(shift))
      indx = 1
      do q=1,nq1(nb)
         !indx = 1 + (q-1)*(nx(nb)/2+1)
         call srfftb(nx(nb),A(shift0+indx),scpl_mb(tmpx(1,nb)))
         indx = indx + nxh
      end do
      call zeroend_sfftb(nx(nb),nq1(nb),1,1,A(shift))
#endif

      shift = shift + nfft3d(nb)
      end do


      end if
    
*     **** deallocate temporary space  ****
      value = BA_pop_stack(tmp3(2))
      value = BA_pop_stack(tmp2(2))
      !value = BA_pop_stack(tmp1(2))

      call nwpw_timing_end(1)
      return
      end




*     ***********************************
*     *                                 *
*     *         D3dBs_pfftfx            *
*     *                                 *
*     ***********************************
*
*     do fft along nx(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dBs_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex A(*)
      complex tmp1(*)
      complex tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz

      !*** local variables ***
      integer j,q,indx,nxh,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      nxh = (nx(id)/2+1)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)

#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(4,id),A,A)

#else

#ifdef NUSE_OPENMP
      call D3dBs_pfftfx_sub(ny(id)*nq(id),nx(id),nxh,
     >                    scpl_mb(tmpx(1,id)),A)
      call cshift_sfftf(nx(id),ny(id),nq(id),1,A)

#else
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call srfftf(nx(id),A(indx),scpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call cshift_fftf(nx(id),ny(id),nq(id),1,A)
#endif
#endif
      !call dcopy((nx(id)+2)*ny(id)*nq(id),A,1,tmp1,1)
      call Parallel_shared_vector_scopy(.true.,
     >                                  (nx(id)+2)*ny(id)*nq(id),A,tmp1)
      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(31)


#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(14,id),A,A)

#else

      offset = tid*(2*nx(id)+15)
      do q=tid+1,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call srfftf(nx(id),A(indx),scpl_mb(tmpx(1,id)+offset))
      end do
      call cshift_sfftf(nx(id),nq1(id),1,1,A)

#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_start(nbb,1,A,tmp1,tmp2,
     >                                 request,reqcnt,40)
      call nwpw_timing_end(32)

      end if

      return
      end

      subroutine  D3dBs_pfftfx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real tmpx(2*nx+15)
      real A(2*nxh,n)
      integer i,indx

      do i=1,n
         call srfftf(nx,A(1,i),tmpx)
      end do
      return
      end





*     ***********************************
*     *                                 *
*     *         D3dBs_pfftfx0           *
*     *                                 *
*     ***********************************
*
*     do fft along nx(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dBs_pfftfx0(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex tmp1(*)
      complex tmp2(*)
      integer    request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz

      !*** local variables ***
      integer j,q,indx,nxh,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      nxh = (nx(id)/2+1)

      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)

#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(4,id),tmp1,tmp1)

#else

#ifdef NUSE_OPENMP
      call D3dBs_pfftfx_sub(ny(id)*nq(id),nx(id),nxh,
     >                      scpl_mb(tmpx(1,id)),tmp1)
      call cshift_sfftf(nx(id),ny(id),nq(id),1,tmp1)

#else
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call srfftf(nx(id),tmp1(indx),scpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call cshift_sfftf(nx(id),ny(id),nq(id),1,tmp1)
#endif
#endif
      call nwpw_timing_end(31)

      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(31)

#ifdef FFTW3
      call sfftw_execute_dft_r2c(splans(14,id),tmp1,tmp1)

#else

      offset = tid*(2*nx(id)+15)
      do q=tid+1,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call srfftf(nx(id),tmp1(indx),scpl_mb(tmpx(1,id)+offset))
      end do
      call cshift_sfftf(nx(id),nq1(id),1,1,tmp1)

#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_start(nbb,1,tmp1,tmp2,tmp1,
     >                                 request,reqcnt,40)
      call nwpw_timing_end(32)

      end if

      return
      end





*     ***********************************
*     *                                 *
*     *         D3dBs_pfftfy            *
*     *                                 *
*     ***********************************
*
*     do fft along ny(id) dimension       
*
*        A(kx,ny(id),nz(id)) <- fft1d[A(nx(id),ny(id),nz(id))] 
*     

      subroutine D3dBs_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,j,k,indx,indx0,indx1,indx2,q,nxh,nxhy,id
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      nxhy = nxh*ny(id)


      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then
      call nwpw_timing_start(31)

*     ********************************************************
*     ***     do fft along ny(id) dimension                ***
*     ***   A(kx,ky,nz(id)) <- fft1d[A(kx,ny(id),nz(id))]  ***
*     ********************************************************

#ifdef FFTW3
       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then
          indx  = i + indx0
          call sfftw_execute_dft(splans(9,1),tmp1(indx),tmp1(indx))
        end if
       end do
       indx0 = indx0 + nxhy
       end do

#else
       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then
          indx  = i + indx0
          indx1 = indx
          do k=1,ny(id)
             tmp2(k) = tmp1(indx)
             indx = indx + nxh
          end do
          call scfftf(ny(id),tmp2,scpl_mb(tmpy(1,id)))
          do k=1,ny(id)
             tmp1(indx1) = tmp2(k)
             indx1 = indx1 + nxh
          end do
        end if
       end do
       indx0 = indx0 + nxhy
       end do
#endif


*     ****************************************************
*     ***         Do a ptranspose of A                 ***
*     ***      A(ky,nz(id),ky) <- A(kx,ky,nz(id))      ***
*     ****************************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose2_jk_start(nbb,tmp1,tmp2,tmp1,
     >                                  request,reqcnt,41)
      call nwpw_timing_end(32)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_end(nbb,1,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ********************************************************
*     ***     do fft along ny(id) dimension                ***
*     ***   A(ky,nz(id),kx) <- fft1d[A(ny(id),nz(id),kx)]  ***
*     ********************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
         call sfftw_execute_dft(splans(19,1),tmp1(indx),tmp1(indx))
         end if
         indx = indx + ny(id)
      end do
#else

      !indx = 1
      offset = tid*(2*ny(id)+15)
      do q=tid+1,nq2(id),nthr
         indx = 1+(q-1)*ny(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call scfftf(ny(id),tmp1(indx),scpl_mb(tmpy(1,id)+offset))
         end if
         !indx = indx + ny(id)
      end do
!$OMP BARRIER

#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_start(nbb,2,tmp1,tmp2,tmp1,
     >                                  request,reqcnt,42)
      call nwpw_timing_end(32)


      end if

      return
      end

      subroutine  D3dBs_pfftfy_sub2(n,ny,row2,tmpy,A)
      implicit none
      integer n,ny
      logical    row2(n)
      real     tmpy(4*ny+15)
      real     A(2*ny,n)
      integer i
      do i=1,n
         if (.not.row2(i)) then
            call scfftf(ny,A(1,i),tmpy)
         end if
      end do
      return
      end




*     ***********************************
*     *                                 *
*     *         D3dBs_pfftfz            *
*     *                                 *
*     ***********************************

      subroutine D3dBs_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


      !*** local variables ***
      integer i,k,q,nxh,nxhz,id
      integer indx,indx0,indx1,indx2
      integer tid,nthr,offset
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      nxhz = nxh*nz(id)



      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose2_jk_end(nbb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ************************************************
*     ***     do fft along nz(id) dimension        ***
*     ***   A(kx,kz,ky) <- fft1d[A(kx,nz(id),ky)]  ***
*     ************************************************


       indx0 = 0
       indx2 = 0
       do q=1,nq(id)
       do i=1,nxh
        indx2 = indx2 + 1
        if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
          indx  = i + indx0
          call sfftw_execute_dft(splans(10,1),tmp2(indx),tmp2(indx))

#else
          indx  = i + indx0
          indx1 = indx
          do k=1,nz(id)
             tmp1(k) = tmp2(indx)
             indx = indx + nxh
          end do
          call scfftf(nz(id),tmp1,scpl_mb(tmpz(1,id)))
          do k=1,nz(id)
             tmp2(indx1) = tmp1(k)
             indx1 = indx1 + nxh
          end do
#endif

        end if
       end do
       indx0 = indx0 + nxhz
       end do
      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_end(nbb,2,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     ************************************************
*     ***     do fft along nz(id) dimension        ***
*     ***   A(kz,kx,ky) <- fft1d[A(nz(id),kx,ky)]  ***
*     ************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(20,1),tmp2(indx),tmp2(indx))
         end if
         indx = indx + nz(id)
      end do
#else

      !indx = 1
      offset = tid*(2*nz(id)+15)
      do q=tid+1,nq3(id),nthr
         indx = 1+(q-1)*nz(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call scfftf(nz(id),tmp2(indx),scpl_mb(tmpz(1,id)+offset))
         end if
         !indx = indx + nz(id)
      end do
!$OMP BARRIER

#endif
      call nwpw_timing_end(31)

      end if

      call nwpw_timing_start(32)
      call Packs_c_pack_start(nbb,tmp2,tmp1,request,reqcnt,43)
      call nwpw_timing_end(32)

      return
      end

      subroutine  D3dBs_pfftfz_sub2(n,nz,row3,tmpz,A)
      implicit none
      integer n,nz
      logical row3(n)
      real tmpz(4*nz+15)
      real A(2*nz,n)
      integer i
      do i=1,n
         if (.not.row3(i)) then
            call scfftf(nz,A(1,i),tmpz)
         end if
      end do
      return
      end




*     ***********************************
*     *                                 *
*     *        D3dBs_pfftf_final        *
*     *                                 *
*     ***********************************
*

      subroutine D3dBs_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      call nwpw_timing_start(32)
      call Packs_c_pack_end(nbb,tmp2,request,reqcnt)
      call nwpw_timing_end(32)

      return
      end




*     ***********************************
*     *                                 *
*     *         D3dBs_pfftf_step        *
*     *                                 *
*     ***********************************
      subroutine D3dBs_pfftf_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex A(*)
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt

     
      if (step.eq.1) then
         call D3dBs_pfftfx(nbb,A,tmp1,tmp2,request,reqcnt)
      else if (step.eq.2) then
         call D3dBs_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.3) then
         call D3dBs_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.4) then
         call D3dBs_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      end if
!$OMP BARRIER

      return
      end


*     ***********************************
*     *                                 *
*     *     D3dBs_pfftf_step_center     *
*     *                                 *
*     ***********************************
      subroutine D3dBs_pfftf_step_center(center,step,nbb,A,tmp1,tmp2,
     >                                  request,reqcnt)
      implicit none
      integer center(3)
      integer step
      integer nbb
      complex A(*)
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt

      if (step.eq.1) then
         call D3dBs_r_1to3_Contract_start3(center,A,tmp2,request,reqcnt)
      else if (step.eq.2) then
         call D3dBs_r_1to3_Contract_mid3(tmp2,tmp1,request,reqcnt)
      else if (step.eq.3) then
         call D3dBs_r_1to3_Contract_final3(tmp1,tmp2,request,reqcnt)

      else if (step.eq.4) then
         call D3dBs_pfftfx0(nbb,tmp2,tmp1,request,reqcnt)
      else if (step.eq.5) then
         call D3dBs_pfftfy(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dBs_pfftfz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.7) then
         call D3dBs_pfftf_final(nbb,tmp1,tmp2,request,reqcnt)
      end if

!$OMP BARRIER

      return
      end






*     ***********************************
*     *					*
*     *	        D3dBs_pfftbz		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of      
*      a three dimensional single precision complex to complex     
*      inverse fft                                
*           A(nx,ny(id),nz(id)) <- FFT3^(-1)[A(kx,ky,kz)]  
*                                                  
*      Entry -                                     
*              A: a column distribuded 3d block    
*              tmp: tempory work space must be at  
*                    least the size of (complex)   
*                    (nfft*nfft + 1) + 10*nfft     
*                                               
*       Exit - A is transformed and the imaginary
*              part of A is set to zero          
*       uses - D3dBs_c_ptranspose_jk, dcopy        
*                                                  


      subroutine D3dBs_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex  tmp1(*)
      complex  tmp2(*)
      integer  request(*),reqcnt
      

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer i,k,q,indx,ierr,id
      integer nxh,nxhz,indx0,indx1,indx2,tid,nthr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

       nxh = (nx(id)/2+1)
       nxhz = nxh*nz(id)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *****************************************************
*     ***     do fft along kz dimension                 ***
*     ***   A(kx,nz(id),ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
*     *****************************************************

      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
       indx2 = indx2 + 1
       if (.not.log_mb(zero_row3(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call sfftw_execute_dft(splans(7,1),tmp1(indx),tmp1(indx))
#else

         indx  = i + indx0
         indx1 = indx
         do k=1,nz(id)
            tmp2(k) = tmp1(indx)
            indx = indx + nxh
         end do
         call scfftb(nz(id),tmp2,scpl_mb(tmpz(1,id)))
         do k=1,nz(id)
            tmp1(indx1) = tmp2(k)
            indx1 = indx1 + nxh
         end do
#endif
       end if
      end do
      indx0 = indx0 + nxhz
      end do

*     ****************************************************
*     ***         Do a ptranspose of A                 ***
*     ***      A(kx,ky,nz(id)) <- A(kx,nz(id),ky)      ***
*     ****************************************************
      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose1_jk_start(nbb,tmp1,tmp2,tmp1,
     >                                  request,reqcnt,44)
      call nwpw_timing_end(32)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else


      call nwpw_timing_start(31)
*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(nz(id),kx,ky) <- fft1d^(-1)[A(kz,kx,ky)]  ***
*     *************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq3(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(17,1),tmp1(indx),tmp1(indx))
         end if
         indx = indx + nz(id)
      end do
#else


      !indx = 1
      do q=1+tid,nq3(id),nthr
         indx = 1 + (q-1)*nz(id)
         if (.not.log_mb(zero_row3(1,nbb)+q-1)) then
            call scfftb(nz(id),tmp1(indx),
     >                 scpl_mb(tmpz(1,id)+tid*(2*nz(id)+15)))
         end if
         !indx = indx + nz(id)
      end do
!$OMP BARRIER

#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_start(nbb,3,tmp1,tmp2,tmp1,
     >                                  request,reqcnt,45)
      call nwpw_timing_end(32)

      end if
      return
      end

      subroutine  D3dBs_pfftbz_sub2(n,nz,row3,tmpz,A)
      implicit none
      integer n,nz
      logical row3(n)
      real    tmpz(4*nz+15)
      real    A(2*nz,n)
      integer i
      do i=1,n
         if (.not.row3(i)) then
            call scfftb(nz,A(1,i),tmpz)
         end if
      end do
      return
      end


*     ***********************************
*     *					*
*     *	        D3dBs_pfftby		*
*     *					*
*     ***********************************

      subroutine D3dBs_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz


*     *** local variables ***
      integer i,j,q,indx,ierr,id
      integer nxh,nxhy,indx0,indx1,indx2,tid,nthr
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

       if (nbb.lt.2) then
          id = 1
       else
          id = 3
       end if
      
       tid  = Parallel_threadid()
       nthr = Parallel_nthreads()

       nxh = (nx(id)/2+1)
       nxhy = nxh*ny(id)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then


*     ********************************************
*     ***         Do a ptranspose of A          ***
*     ***      A(kx,ky,nz(id)) <- A(kx,nz(id),ky)      ***
*     ********************************************
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose1_jk_end(nbb,tmp2,tmp1,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(kx,ny(id),nz(id)) <- fft1d^(-1)[A(kx,ky,nz(id))]  ***
*     *************************************************************
    
      indx0 = 0
      indx2 = 0
      do q=1,nq(id)
      do i=1,nxh
        indx2 = indx2 + 1

        if (.not.log_mb(zero_row2(1,nbb)+indx2-1)) then

#ifdef FFTW3
         indx  = i + indx0
         call sfftw_execute_dft(splans(8,id),tmp2(indx),tmp2(indx))

#else
         indx  = i + indx0
         indx1 = indx
         do j=1,ny(id)
            tmp1(j) = tmp2(indx)
            indx = indx + nxh
         end do
         call scfftb(ny(id),tmp1,scpl_mb(tmpy(1,id)))
         do j=1,ny(id)
            tmp2(indx1) = tmp1(j)
            indx1 = indx1 + nxh
         end do
#endif

        end if
      end do
      indx0 = indx0 + nxhy
      end do

      call nwpw_timing_end(31)


      !*************************
      !**** hilbert mapping ****
      !*************************
      else

      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_end(nbb,3,tmp2,tmp1,request,reqcnt)

      call nwpw_timing_end(32)
      call nwpw_timing_start(31)
*     *************************************************************
*     ***     do fft along ky dimension                         ***
*     ***   A(ny(id),nz(id),kx) <- fft1d^(-1)[A(ky,nz(id),kx)]  ***
*     *************************************************************

#ifdef FFTW3
      indx = 1
      do q=1,nq2(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call sfftw_execute_dft(splans(18,id),tmp2(indx),tmp2(indx))
         end if
         indx = indx + ny(id)
      end do
#else
      do q=1+tid,nq2(id),nthr
         indx = 1 + (q-1)*ny(id)
         if (.not.log_mb(zero_row2(1,nbb)+q-1)) then
            call scfftb(ny(id),tmp2(indx),
     >                 scpl_mb(tmpy(1,id)+tid*(2*ny(id)+15)))
         end if
      end do
!$OMP BARRIER

#endif

      call nwpw_timing_end(31)
      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_start(nbb,4,tmp2,tmp1,tmp2,
     >                                  request,reqcnt,46)
      call nwpw_timing_end(32)


      end if

      return
      end

      subroutine  D3dBs_pfftby_sub2(n,ny,row2,tmpy,A)
      implicit none
      integer n,ny
      logical row2(n)
      real    tmpy(4*ny+15)
      real    A(2*ny,n)
      integer i
      do i=1,n
         if (.not.row2(i)) then
            call scfftb(ny,A(1,i),tmpy)
         end if
      end do
      return
      end





*     ***********************************
*     *					*
*     *	        D3dBs_pfftbx		*
*     *					*
*     ***********************************

      subroutine D3dBs_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      implicit none
      integer nbb
      complex  tmp1(*)
      complex  tmp2(*)
      integer  request(*),reqcnt

#include "bafdecls.fh"
#include "errquit.fh"
#include "D3dB.fh"
#include "D3dB_pfft.fh"

      integer tmpx(2,NBLOCKS),tmpy(2,NBLOCKS),tmpz(2,NBLOCKS)
      common    / D3dBs_fft / tmpx,tmpy,tmpz



*     *** local variables ***
      integer j,q,indx,ierr,id
      integer nxh,tid,nthr

*     **** external functions ****
      integer  Parallel_threadid,Parallel_nthreads
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      tid  = Parallel_threadid()
      nthr = Parallel_nthreads()

      nxh = (nx(id)/2+1)
      
      !**********************
      !**** slab mapping ****
      !**********************
      if (mapping.eq.1) then

      call nwpw_timing_start(31)
*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(id),ny(id),nz(id)) <- fft1d^(-1)[A(kx,ny(id),nz(id))]  ***
*     *********************************************************************


#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(3,id),tmp2,tmp2)

#else
#ifdef NUSE_OPENMP
      call cshift1_sfftb(nx(id),ny(id),nq(id),1,tmp2)
      call D3dBs_pfftbx_sub(nq(id)*ny(id),nx(id),nxh,
     >                      scpl_mb(tmpx(1,id)),tmp2)
      call zeroend_sfftb(nx(id),ny(id),nq(id),1,tmp2)
#else
      call cshift1_sfftb(nx(id),ny(id),nq(id),1,tmp2)
      indx = 1
      do q=1,nq(id)
      do j=1,ny(id)
         call srfftb(nx(id),tmp2(indx),scpl_mb(tmpx(1,id)))
         indx = indx + nxh
      end do
      end do
      call zeroend_sfftb(nx(id),ny(id),nq(id),1,tmp2)
#endif
#endif
      call Parallel_shared_vector_scopy(.true.,
     >                                  (nx(id)+2)*ny(id)*nq(id),
     >                                  tmp2,tmp1)

      call nwpw_timing_end(31)
*     *************************************************
      !*************************
      !**** hilbert mapping ****
      !*************************
      else



      call nwpw_timing_start(32)
      call D3dBs_c_ptranspose_ijk_end(nbb,4,tmp1,tmp2,request,reqcnt)
      call nwpw_timing_end(32)
      call nwpw_timing_start(31)

*     *********************************************************************
*     ***     do fft along kx dimension                                 ***
*     ***   A(nx(id),ny(id),nz(id)) <- fft1d^(-1)[A(kx,ny(id),nz(id))]  ***
*     *********************************************************************

#ifdef FFTW3
      call sfftw_execute_dft_c2r(splans(13,id),tmp1,tmp1)

#else
      indx = 1
      do q=1+tid,nq1(id),nthr
         indx = 1+(q-1)*nxh
         call cshift1_sfftb1(nx(id),tmp1(indx))
         call srfftb(nx(id),tmp1(indx),
     >               scpl_mb(tmpx(1,id)+tid*(2*nx(id)+15))) 
         call zeroend_sfftb1(nx(id),tmp1(indx))
         !indx = indx + nxh
      end do
#endif

      call nwpw_timing_end(31)
      end if
    
      return
      end


      subroutine  D3dBs_pfftbx_sub(n,nx,nxh,tmpx,A)
      implicit none
      integer n,nx,nxh
      real    tmpx(2*nx+15)
      real    A(2*nxh,n)
      integer i
      do i=1,n
         call srfftb(nx,A(1,i),tmpx)
      end do

      return
      end



*     ***********************************
*     *                                 *
*     *         D3dBs_pfftb_step        *
*     *                                 *
*     ***********************************
*

      subroutine D3dBs_pfftb_step(step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer step
      integer nbb
      complex A(*)
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"

      integer tmp1z(2),tmp2z(2)
      common / pffts_queuez_common / tmp1z,tmp2z


      if (step.eq.1) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_start(nbb,A,tmp1,request,reqcnt,47)
         call nwpw_timing_end(32)
      else if (step.eq.2) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_mid(nbb,tmp1,tmp2,
     >                          scpl_mb(tmp1z(1)),scpl_mb(tmp2z(1)),
     >                          request,reqcnt,48)
         call nwpw_timing_end(32)
      else if (step.eq.3) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_end(nbb,tmp1,
     >                          scpl_mb(tmp1z(1)),scpl_mb(tmp2z(1)),
     >                          request,reqcnt)
         call nwpw_timing_end(32)
      else if (step.eq.4) then
         call D3dBs_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.5) then
         call D3dBs_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dBs_pfftbx(nbb,tmp1,tmp2,request,reqcnt)
      end if

      return
      end



*     ***********************************
*     *                                 *
*     *      D3dBs_pfftb_step_center    *
*     *                                 *
*     ***********************************

      subroutine D3dBs_pfftb_step_center(center,
     >                            step,nbb,A,tmp1,tmp2,request,reqcnt)
      implicit none
      integer center(3)
      integer step
      integer nbb
      complex A(*)
      complex tmp1(*)
      complex tmp2(*)
      integer request(*),reqcnt


#include "bafdecls.fh"

      integer tmp1z(2),tmp2z(2)
      common / pffts_queuez_common / tmp1z,tmp2z

      if (step.eq.1) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_start(nbb,A,tmp1,request,reqcnt,47)
         call nwpw_timing_end(32)
      else if (step.eq.2) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_mid(nbb,tmp1,tmp2,
     >                          scpl_mb(tmp1z(1)),scpl_mb(tmp2z(1)),
     >                          request,reqcnt,48)
         call nwpw_timing_end(32)
      else if (step.eq.3) then
         call nwpw_timing_start(32)
         call Packs_c_unpack_end(nbb,tmp1,
     >                          scpl_mb(tmp1z(1)),scpl_mb(tmp2z(1)),
     >                          request,reqcnt)
         call nwpw_timing_end(32)
      else if (step.eq.4) then
         call D3dBs_pfftbz(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.5) then
         call D3dBs_pfftby(nbb,tmp1,tmp2,request,reqcnt)
      else if (step.eq.6) then
         call D3dBs_pfftbx(nbb,tmp1,tmp2,request,reqcnt)

      else if (step.eq.7) then
         call D3dBs_r_3to1_Expand_start3(center,tmp1,tmp2,
     >                                   request,reqcnt)
      else if (step.eq.8) then
         call D3dBs_r_3to1_Expand_mid3(tmp2,tmp1,request,reqcnt)
      else if (step.eq.9) then
         call D3dBs_r_3to1_Expand_final3(tmp1,tmp2,request,reqcnt)
      end if

      return
      end


*     ***********************************
*     *                                 *
*     *     D3dBs_pfft3_queue_init      *
*     *                                 *
*     ***********************************

      subroutine D3dBs_pfft3_queue_init(qmax_in)
      implicit none
      integer qmax_in


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),acenter(2)
      integer aqmax,aqsize,alast_index
      common / pffts_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                               areqcnt,acenter,
     >                               aqmax,aqsize,alast_index

      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bcenter(2)
      integer bqmax,bqsize,blast_index
      common / pffts_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                               breqcnt,bcenter,
     >                               bqmax,bqsize,blast_index
      
      integer tmp1z(2),tmp2z(2)
      common / pffts_queuez_common / tmp1z,tmp2z

      logical value
      integer np,size,zplane_sizetmp

      call Parallel2d_np_i(np)

c     **** allocate aqueue ****
      aqmax       = qmax_in
      aqsize      = 0
      alast_index = aqmax
      size  = nfft3d(1)*aqmax*2
      value = BA_alloc_get(mt_scpl,size,'atmp',atmp(2),atmp(1))
      size  = np*aqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,2*size,'arequest',
     >                   arequest(2),arequest(1))
      size  = aqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'aqindx',aqindx(2),aqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'aqstatus',aqstatus(2),aqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'areqcnt',areqcnt(2),areqcnt(1))
      size  = 3*aqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'acenter',acenter(2),acenter(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      !call dcopy(4*nfft3d(1)*aqmax,0.0d0,0,dcpl_mb(atmp(1)),1)
      call Parallel_shared_vector_szero(.true.,4*nfft3d(1)*aqmax,
     >                                  scpl_mb(atmp(1)))


c     **** allocate bqueue ****
      bqmax       = qmax_in
      bqsize      = 0
      blast_index = bqmax
      size  = nfft3d(1)*bqmax*2
      value = BA_alloc_get(mt_scpl,size,'btmp',btmp(2),btmp(1))
      size  = np*bqmax*2
      value = value.and.
     >      BA_alloc_get(mt_int,2*size,'brequest',
     >                   brequest(2),brequest(1))
      size  = bqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bqindx',bqindx(2),bqindx(1))
      value = value.and.
     >      BA_alloc_get(mt_int,size,'bqstatus',bqstatus(2),bqstatus(1))
      value = value.and.
     >        BA_alloc_get(mt_int,size,'breqcnt',breqcnt(2),breqcnt(1))
      size  = 3*bqmax
      value = value.and.
     >        BA_alloc_get(mt_int,size,'bcenter',bcenter(2),bcenter(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      !call dcopy(4*nfft3d(1)*bqmax,0.0d0,0,dcpl_mb(btmp(1)),1)
      call Parallel_shared_vector_szero(.true.,4*nfft3d(1)*bqmax,
     >                                  scpl_mb(btmp(1)))



      call D3dB_zplane_size(1,zplane_sizetmp)
      value = value.and.
     >        BA_alloc_get(mt_scpl,(zplane_sizetmp),
     >                    'tmp1z',tmp1z(2),tmp1z(1))
      value = value.and.
     >        BA_alloc_get(mt_scpl,(zplane_sizetmp),
     >                    'tmp2z',tmp2z(2),tmp2z(1))
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)
      

      return
      end




*     ***********************************
*     *                                 *
*     *     D3dBs_pfft3_queue_end       *
*     *                                 *
*     ***********************************

      subroutine D3dBs_pfft3_queue_end()
      implicit none


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queues ****
      integer aqindx(2),aqstatus(2)
      integer atmp(2),arequest(2),areqcnt(2),acenter(2)
      integer aqmax,aqsize,alast_index
      common / pffts_aqueue_common / aqindx,aqstatus,atmp,arequest,
     >                               areqcnt,acenter,
     >                               aqmax,aqsize,alast_index

      integer bqindx(2),bqstatus(2)
      integer btmp(2),brequest(2),breqcnt(2),bcenter(2)
      integer bqmax,bqsize,blast_index
      common / pffts_bqueue_common / bqindx,bqstatus,btmp,brequest,
     >                               breqcnt,bcenter,
     >                               bqmax,bqsize,blast_index

      integer tmp1z(2),tmp2z(2)
      common / pffts_queuez_common / tmp1z,tmp2z

      logical value

      value =           BA_free_heap(atmp(2))
      value = value.and.BA_free_heap(arequest(2))
      value = value.and.BA_free_heap(aqindx(2))
      value = value.and.BA_free_heap(aqstatus(2))
      value = value.and.BA_free_heap(areqcnt(2))
      value = value.and.BA_free_heap(acenter(2))

      value = value.and.BA_free_heap(btmp(2))
      value = value.and.BA_free_heap(brequest(2))
      value = value.and.BA_free_heap(bqindx(2))
      value = value.and.BA_free_heap(bqstatus(2))
      value = value.and.BA_free_heap(breqcnt(2))
      value = value.and.BA_free_heap(bcenter(2))

      value = value.and.BA_free_heap(tmp1z(2))
      value = value.and.BA_free_heap(tmp2z(2))
      if (.not. value) 
     > call errquit('D3dBs_pfft3_queue_end:error freeing heap',0,MA_ERR)

      return
      end







*     ***********************************
*     *                                 *
*     *    D3dBs_cr_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function D3dBs_cr_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      D3dBs_cr_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end



*     ***********************************
*     *                                 *
*     *    D3dBs_rc_pfft3_queue_filled  *
*     *                                 *
*     ***********************************

      logical function D3dBs_rc_pfft3_queue_filled()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      D3dBs_rc_pfft3_queue_filled = (qsize.ge.qmax)
      return
      end




*     ***********************************
*     *                                 *
*     *     D3dBs_rc_pfft3f_queuein     *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rc_pfft3f_queuein(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,status,np

      call nwpw_timing_start(30)
      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call Parallel2d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
     
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(id)*(2*(indx-1))
         shift2=nfft3d(id)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call D3dBs_pfftf_step(status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(id)*(2*(last_index-1))
      shift2=nfft3d(id)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call D3dBs_pfftf_step(status,nbb,A,
     >                     scpl_mb(tmp(1)+shift1),
     >                     scpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

      call nwpw_timing_end(30)
      return
      end




*     ***********************************
*     *                                 *
*     *     D3dBs_rc_pfft3f_queueout    *
*     *                                 *
*     ***********************************

      subroutine D3dBs_rc_pfft3f_queueout(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,indx1,status,np

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)


      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.4)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
           
            status = int_mb(qstatus(1)+indx-1)+1
            shift1=nfft3d(id)*(2*(indx-1))
            shift2=nfft3d(id)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call D3dBs_pfftf_step(status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do


      shift2=nfft3d(id)*(2*(indx1-1)+1)
      call Pack_c_sCopy(nbb,scpl_mb(tmp(1)+shift2),A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      call nwpw_timing_end(30)
      return
      end



*     ***********************************
*     *                                 *
*     *     D3dBs_cr_pfft3b_queuein      *
*     *                                 *
*     ***********************************

      subroutine D3dBs_cr_pfft3b_queuein(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3,id
      integer q,indx,status,np
      integer  Parallel_threadid,Parallel_nthreads,tid,nthr
      external Parallel_threadid,Parallel_nthreads

      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if

      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)


      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
        
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(id)*(2*(indx-1))
         shift2=nfft3d(id)*(2*(indx-1)+1)
         shift3=2*np*(indx-1)
         call D3dBs_pfftb_step(status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      
        call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do


      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(id)*(2*(last_index-1))
      shift2=nfft3d(id)*(2*(last_index-1)+1)
      shift3=2*np*(last_index-1)

      call D3dBs_pfftb_step(status,nbb,A,
     >                     scpl_mb(tmp(1)+shift1),
     >                     scpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

!$OMP BARRIER
      call nwpw_timing_end(30)

      return
      end




*     ***********************************
*     *                                 *
*     *     D3dBs_cr_pfft3b_queueout    *
*     *                                 *
*     ***********************************

      subroutine D3dBs_cr_pfft3b_queueout(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter, qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np,id


      if (nbb.lt.2) then
         id = 1
      else
         id = 3
      end if
c      tid = Parallel_threadid()
      call Parallel2d_np_i(np)


      call nwpw_timing_start(30)

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.6)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1) + 1
            shift1=nfft3d(id)*(2*(indx-1))
            shift2=nfft3d(id)*(2*(indx-1)+1)
            shift3=2*np*(indx-1)
            call D3dBs_pfftb_step(status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))

            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do
!$OMP BARRIER


      shift1=nfft3d(id)*(2*(indx1-1))
      call Parallel_shared_vector_scopy(.true.,2*nfft3d(id),
     >                                 scpl_mb(tmp(1)+shift1),
     >                                 A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      call nwpw_timing_end(30)
      return
      end





cccccccccccccccccccccc  queue center routines   ccccccccccccccccccccccccccccccccccccc



*     ******************************************
*     *                                        *
*     *     D3dBs_rc_pfft3f_queuein_center      *
*     *                                        *
*     ******************************************

      subroutine D3dBs_rc_pfft3f_queuein_center(nbb,center,A)
      implicit none
      integer nbb
      integer center(3)
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np

      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)

      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)
         !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=4*np*(indx-1)
         call D3dBs_pfftf_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
         call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do

      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_vector_icopy(.false.,3,center,
     >                           int_mb(qcenter(1)+3*(last_index-1)))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=4*np*(last_index-1)

      call D3dBs_pfftf_step_center(int_mb(qcenter(1)+3*(last_index-1)),
     >                     status,nbb,A,
     >                     scpl_mb(tmp(1)+shift1),
     >                     scpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

      call nwpw_timing_end(30)
      return
      end




*     ******************************************
*     *                                        *
*     *     D3dBs_rc_pfft3f_queueout_center    *
*     *                                        *
*     ******************************************

      subroutine D3dBs_rc_pfft3f_queueout_center(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_bqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np

      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.7)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1)+1
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=4*np*(indx-1)
            call D3dBs_pfftf_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do

      shift2=nfft3d(1)*(2*(indx1-1)+1)
      call Pack_c_sCopy(nbb,scpl_mb(tmp(1)+shift2),A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      call nwpw_timing_end(30)
      return
      end





*     ******************************************
*     *                                        *
*     *     D3dBs_cr_pfft3b_queuein_center     *
*     *                                        *
*     ******************************************

      subroutine D3dBs_cr_pfft3b_queuein_center(nbb,center,A)
      implicit none
      integer nbb
      integer center(3)
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                               qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,status,np
      integer  Parallel_threadid,Parallel_nthreads,tid,nthr
      external Parallel_threadid,Parallel_nthreads


      call nwpw_timing_start(30)
      call Parallel2d_np_i(np)
c      tid  = Parallel_threadid()
c      nthr = Parallel_nthreads()


      do q=1,qsize
         indx   = int_mb(qindx(1)+q-1)

         status = int_mb(qstatus(1)+indx-1)+1
         shift1=nfft3d(1)*(2*(indx-1))
         shift2=nfft3d(1)*(2*(indx-1)+1)
         shift3=4*np*(indx-1)
         call D3dBs_pfftb_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))
      
        call Parallel_shared_iinc(.true.,1,int_mb(qstatus(1)+indx-1))
      end do


      call Parallel_shared_irot(.true.,qmax,last_index)
      call Parallel_shared_iinc(.false.,1,qsize)
      call Parallel_shared_icopy(.false.,last_index,
     >                           int_mb(qindx(1)+qsize-1))
      call Parallel_shared_vector_icopy(.false.,3,center,
     >                           int_mb(qcenter(1)+3*(last_index-1)))
      call Parallel_shared_icopy(.true.,1,
     >                           int_mb(qstatus(1)+last_index-1))

      status = 1
      shift1=nfft3d(1)*(2*(last_index-1))
      shift2=nfft3d(1)*(2*(last_index-1)+1)
      shift3=4*np*(last_index-1)

      call D3dBs_pfftb_step_center(int_mb(qcenter(1)+3*(last_index-1)),
     >                     status,nbb,A,
     >                     scpl_mb(tmp(1)+shift1),
     >                     scpl_mb(tmp(1)+shift2),
     >                     int_mb(request(1)+shift3),
     >                     int_mb(reqcnt(1)+last_index-1))

!$OMP BARRIER
      call nwpw_timing_end(30)

      return
      end




*     ******************************************
*     *                                        *
*     *     D3dBs_cr_pfft3b_queueout_center    *
*     *                                        *
*     ******************************************

      subroutine D3dBs_cr_pfft3b_queueout_center(nbb,A)
      implicit none
      integer nbb
      complex A(*)


#include "bafdecls.fh"
#include "errquit.fh"

#include "D3dB.fh"
#include "D3dB_pfft.fh"


      !**** pfft_queue ****
      integer qindx(2),qstatus(2)
      integer tmp(2),request(2),reqcnt(2),qcenter(2)
      integer qmax,qsize,last_index
      common / pffts_aqueue_common / qindx,qstatus,tmp,request,reqcnt,
     >                              qcenter,qmax,qsize,last_index

      !*** local variables ***
      integer shift1,shift2,shift3
      integer q,indx,indx1,status,np
c      integer  Parallel_threadid,tid
c      external Parallel_threadid


c      tid = Parallel_threadid()
      call Parallel2d_np_i(np)

      call nwpw_timing_start(30)

      indx1   = int_mb(qindx(1))
      do while (int_mb(qstatus(1)+indx1-1).lt.9)
         do q=1,qsize
            indx   = int_mb(qindx(1)+q-1)
            !int_mb(qstatus(1)+indx-1) = int_mb(qstatus(1)+indx-1) + 1
            status = int_mb(qstatus(1)+indx-1) + 1
            shift1=nfft3d(1)*(2*(indx-1))
            shift2=nfft3d(1)*(2*(indx-1)+1)
            shift3=4*np*(indx-1)
            call D3dBs_pfftb_step_center(int_mb(qcenter(1)+3*(indx-1)),
     >                        status,nbb,A,
     >                        scpl_mb(tmp(1)+shift1),
     >                        scpl_mb(tmp(1)+shift2),
     >                        int_mb(request(1)+shift3),
     >                        int_mb(reqcnt(1)+indx-1))

            call Parallel_shared_iinc(.true.,1,
     >                                int_mb(qstatus(1)+indx-1))
         end do
      end do
!$OMP BARRIER


      shift2=nfft3d(1)*(2*(indx1-1)+1)
      call Parallel_shared_vector_scopy(.true.,2*nfft3d(1),
     >                                 scpl_mb(tmp(1)+shift2),
     >                                 A)

      call Parallel_shared_iinc(.true.,-1,qsize)
      call Parallel_shared_lshift(.true.,qsize,int_mb(qindx(1)))

      call nwpw_timing_end(30)
      return
      end


