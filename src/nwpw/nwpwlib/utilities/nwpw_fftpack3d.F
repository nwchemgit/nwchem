*     ***********************************
*     *					*
*     *	        nwpw_cr_fft3b		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of     
*      a three dimensional complex to real
*      inverse fft                                
*           A(nx,ny,nz) <- FFT3^(-1)[A(kx,ky,kz)] 
*                                                 
*      Entry - 					  
*              tmp: tempory work space must be at 
*                    least the size of (complex)  
*                    (nfft*nfft + 1) + 10*nfft    
*                                                 
*       Exit - A is transformed and the imaginary 
*              part of A is set to zero           
*                                                 

      subroutine nwpw_cr_fft3b(A,
     >                    inc2c,inc3c,
     >                    inc2r,inc3r,
     >                    nx,ny,nz,tmp,ntmp)

      implicit none
      
      real*8  A(*)
      integer inc2c,inc3c
      integer inc2r,inc3r
      integer nx,ny,nz
      real*8  tmp(*)
      integer  ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tz,ty,tx

      tz = 2*nz+1
      ty = 2*ny+1
      tx = nx+3

*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,ky,nz) <- fft1d^(-1)[A(kx,ky,kz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call dcffti(nz,tmp(tz))
      do j=1,ny
      do i=1,(nx/2+1)
         indx = (2*i-1) + (j-1)*(inc2r) 
         call zcopy(nz,A(indx),inc3c,tmp,1)
         call dcfftb(nz,tmp,tmp(tz))
         call zcopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do


*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call dcffti(ny,tmp(ty))
      do k=1,nz
      do i=1,(nx/2+1)
         indx = (2*i-1) + (k-1)*(inc3r) 
         call zcopy(ny,A(indx),inc2c,tmp,1)
         call dcfftb(ny,tmp,tmp(ty))
         call zcopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)]  ***
*     ***   1-d complex to real transforms          ***
*     *************************************************
      call drffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2r) + (k-1)*(inc3r)
         call dcopy((nx+2),A(indx),1,tmp,1)
         do i=2,nx
            tmp(i) = tmp(i+1)
         end do
         call drfftb(nx,tmp,tmp(tx))
         tmp(nx+1) = 0.0d0
         tmp(nx+2) = 0.0d0
         call dcopy((nx+2),tmp,1,A(indx),1)
      end do
      end do
      return
      end




*     ***********************************
*     *					*
*     *	       nwpw_rc_fft3f		*
*     *					*
*     ***********************************

      subroutine nwpw_rc_fft3f(A,
     >                    inc2c,inc3c,
     >                    inc2r,inc3r,
     >                    nx,ny,nz,tmp,ntmp)


*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx,ny,nz)]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*                                                   *
*****************************************************
      implicit none
      real*8   A(*)
      integer inc2c,inc3c
      integer inc2r,inc3r
      integer nx,ny,nz
      real*8  tmp(*)
      integer ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tx,ty,tz

      tx = nx   + 3
      ty = 2*ny + 1
      tz = 2*nz + 1


*     ********************************************
*     ***     do fft along nx dimension        ***
*     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
*     ********************************************
      call drffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2r) + (k-1)*(inc3r)
         call dcopy((nx+2),A(indx),1,tmp,1)
         call drfftf(nx,tmp,tmp(tx))
         do i=nx,2,-1
            tmp(i+1) = tmp(i)
         end do
         tmp(2)    = 0.0d0
         tmp(nx+2) = 0.0d0
         call dcopy((nx+2),tmp,1,A(indx),1)
      end do
      end do

*     ********************************************
*     ***     do fft along ny dimension        ***
*     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
*     ********************************************
      call dcffti(ny,tmp(ty))
      do k=1,nz
      do i=1,(nx/2+1)
         indx = (2*i-1) + (k-1)*(inc3r) 
         call zcopy(ny,A(indx),inc2c,tmp,1)
         call dcfftf(ny,tmp,tmp(ty))
         call zcopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do


*     ********************************************
*     ***     do fft along nz dimension        ***
*     ***   A(kx,ky,kz) <- fft1d[A(kx,ky,nz)]  ***
*     ********************************************
      call dcffti(nz,tmp(tz))
      do j=1,ny
      do i=1,(nx/2+1)
         indx = (2*i-1) + (j-1)*(inc2r) 
         call zcopy(nz,A(indx),inc3c,tmp,1)
         call dcfftf(nz,tmp,tmp(tz))
         call zcopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do

      return
      end


*     ***********************************
*     *					*
*     *	        nwpw_cc_fft3b		*
*     *					*
*     ***********************************

*                                                  
*      This routine performs the operation of     
*      a three dimensional complex to complex     
*      inverse fft                                
*           A(nx,ny,nz) <- FFT3^(-1)[A(kx,ky,kz)] 
*                                                 
*      Entry - 					  
*              tmp: tempory work space must be at 
*                    least the size of (complex)  
*                    (nfft*nfft + 1) + 10*nfft    
*                                                 
*       Exit - A is transformed and the imaginary 
*              part of A is set to zero           
*                                                 

      subroutine nwpw_cc_fft3b(A,inc2c,inc3c,
     >                    nx,ny,nz,tmp,ntmp)

      implicit none
      
      complex*16 A(*)
      integer    inc2c,inc3c
      integer    nx,ny,nz
      complex*16 tmp(*)
      integer    ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tz,ty,tx

      tz = nz+1
      ty = ny+1
      tx = nx+1


*     *************************************************
*     ***     do fft along kz dimension             ***
*     ***   A(kx,ky,nz) <- fft1d^(-1)[A(kx,ky,kz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call dcffti(nz,tmp(tz))
      do j=1,ny
      do i=1,nx
         indx = i + (j-1)*(inc2c) 
         call zcopy(nz,A(indx),inc3c,tmp,1)
         call dcfftb(nz,tmp,tmp(tz))
         call zcopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do


*     *************************************************
*     ***     do fft along ky dimension             ***
*     ***   A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call dcffti(ny,tmp(ty))
      do k=1,nz
      do i=1,nx
         indx = i + (k-1)*(inc3c) 
         call zcopy(ny,A(indx),inc2c,tmp,1)
         call dcfftb(ny,tmp,tmp(ty))
         call zcopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do

*     *************************************************
*     ***     do fft along kx dimension             ***
*     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)]  ***
*     ***   1-d complex to complex transforms       ***
*     *************************************************
      call dcffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2c) + (k-1)*(inc3c)
         call zcopy(nx,A(indx),1,tmp,1)
         call dcfftb(nx,tmp,tmp(tx))
         call zcopy(nx,tmp,1,A(indx),1)
      end do
      end do
      return
      end




*     ***********************************
*     *					*
*     *	       nwpw_cc_fft3f		*
*     *					*
*     ***********************************

      subroutine nwpw_cc_fft3f(A,inc2c,inc3c,
     >                    nx,ny,nz,tmp,ntmp)


*****************************************************
*                                                   *
*      This routine performs the operation of       *
*      a three dimensional complex to complex fft   *
*           A(kx,ky,kz) <- FFT3[A(nx,ny,nz)]        * 
*                                                   *
*      Entry - 					    *
*              A: a column distribuded 3d block     *
*              tmp: tempory work space must be at   *
*                    least the size of (complex)    *
*                    (nfft*nfft + 1) + 10*nfft      * 
*                                                   *
*       Exit - A is transformed                     *
*                                                   *
*                                                   *
*****************************************************
      implicit none
      complex*16 A(*)
      integer    inc2c,inc3c
      integer    nx,ny,nz
      complex*16 tmp(*)
      integer    ntmp


*     *** local variables ***
      integer i,j,k,indx
      integer tx,ty,tz

      tx = nx + 1
      ty = ny + 1
      tz = nz + 1


*     ********************************************
*     ***     do fft along nx dimension        ***
*     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
*     ********************************************
      call dcffti(nx,tmp(tx))
      do k=1,nz
      do j=1,ny
         indx = 1 + (j-1)*(inc2c) + (k-1)*(inc3c)
         call zcopy(nx,A(indx),1,tmp,1)
         call dcfftf(nx,tmp,tmp(tx))
         call zcopy(nx,tmp,1,A(indx),1)
      end do
      end do

*     ********************************************
*     ***     do fft along ny dimension        ***
*     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
*     ********************************************
      call dcffti(ny,tmp(ty))
      do k=1,nz
      do i=1,nx
         indx = i + (k-1)*(inc3c) 
         call zcopy(ny,A(indx),inc2c,tmp,1)
         call dcfftf(ny,tmp,tmp(ty))
         call zcopy(ny,tmp,1,A(indx),inc2c)
      end do
      end do


*     ********************************************
*     ***     do fft along nz dimension        ***
*     ***   A(kx,ky,kz) <- fft1d[A(kx,ky,nz)]  ***
*     ********************************************
      call dcffti(nz,tmp(tz))
      do j=1,ny
      do i=1,nx
         indx = i + (j-1)*(inc2c) 
         call zcopy(nz,A(indx),inc3c,tmp,1)
         call dcfftf(nz,tmp,tmp(tz))
         call zcopy(nz,tmp,1,A(indx),inc3c)
      end do
      end do

      return
      end


