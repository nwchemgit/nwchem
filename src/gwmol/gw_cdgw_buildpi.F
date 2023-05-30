#include "gw_simd.fh"
      subroutine gw_cdgw_buildpi(which,ovup,ovdw,pi,factor,ldpi,nri,
     $                           npoles,main,ipol,maxpoles)
C
C     Build the polarizability matrix \Pi (Equation (39)) and transform
C     it to the dielectric matrix (inline Equation just below).
C
C     Two cases should be distinguished: 
C
C     When the argument of Equation (39) is purely imaginary, the factor 
C     in square brackets is real and positive. In this case, we can 
C     compute the square root of the factor and use DSFRK (RFP version of DSYRK), 
C     in order to obtain \Pi.
C
C     When the argument is a general complex number, the factor is no
C     longer guaranteed to be real nor positive. In thi case, we could
C     either work in complex arithmetics and take the square root of the
C     factor, or we can work in real arithmetics and assemble \Pi via 
C     DSYR2K. Here, the latter approach is implemented. 
C
C     Note that we actually compute -\Pi.

      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
      logical,intent(in) :: main
      integer,intent(in) :: ldpi, nri, maxpoles, ipol
      character(len=1),intent(in) :: which
      integer,intent(in),dimension(ipol) :: npoles
      double precision,intent(inout),dimension(nri,*) :: ovup,ovdw

      double precision,intent(in),dimension(maxpoles,ipol) :: factor
      double precision,intent(out),dimension(ldpi,*) :: pi

      double precision scal, t0, scal2
      integer ipole,iri,ma_dbl,l_scr,k_scr,mynpoles,llpole,mymax
      integer ulpole,nbatch,ibatch,idx,idy
      logical ok,dorfp
      integer nthreads
#ifdef USE_OPENMP
      integer, external :: omp_get_max_threads
      nthreads = omp_get_max_threads()
#else
      nthreads = 1
#endif

      !
      ! Scaling factor for closed-shells
      !
      scal = 2d0
      if (ipol.eq.1) scal = 2d0*scal

      !
      ! Case 1. Purely imaginary omegas
      !
      if (which.eq.'w') then

        !
        ! Alpha spin channel
        !
        if (npoles(1).gt.0) then

          !
          ! Scale ERIs with factor in square brackets of Equation (39)
          !
!$omp     parallel do MYSIMD schedule(static) collapse(2)
          do ipole=1,npoles(1)
          do iri=1,nri
            ovup(iri,ipole) = ovup(iri,ipole)*factor(ipole,1)
          enddo
          enddo
!$omp     end parallel do MYSIMD

          !
          ! DSFRK is the RFP version of DSYRK
          !
          call dsfrk('n','l','n',nri,npoles(1),scal,ovup,nri,0d0,pi)

          !
          ! Back-transform ERIs, since we are going to need them again
          !
!$omp     parallel do MYSIMD schedule(static) collapse(2)
          do ipole=1,npoles(1)
          do iri=1,nri
            ovup(iri,ipole) = ovup(iri,ipole)/factor(ipole,1)
          enddo
          enddo
!$omp     end parallel do MYSIMD

        endif

        !
        ! Beta spin channel
        !
        if (ipol.gt.1.and.npoles(2).gt.0) then

          !
          ! Scale ERIs with factor in square brackets of Equation (39)
          !
!$omp     parallel do MYSIMD schedule(static) collapse(2)
          do ipole=1,npoles(2)
          do iri=1,nri
            ovdw(iri,ipole) = ovdw(iri,ipole)*factor(ipole,2)
          enddo
          enddo
!$omp     end parallel do MYSIMD

          !
          ! DSFRK is the RFP version of DSYRK
          !
          call dsfrk('n','l','n',nri,npoles(2),scal,ovdw,nri,1d0,pi)

          !
          ! Back-transform ERIs, since we are going to need them again
          !
!$omp     parallel do MYSIMD schedule(static) collapse(2)
          do ipole=1,npoles(2)
          do iri=1,nri
            ovdw(iri,ipole) = ovdw(iri,ipole)/factor(ipole,2)
          enddo
          enddo
!$omp     end parallel do MYSIMD

        endif

        !
        ! Transform \Pi into the dielectric matrix \epsilon. Since
        ! contributions from all MPI ranks are going to be summed, only
        ! one rank needs to add 1.
        !
        ! Note: Here we employ the rigt RFP format since the matrix
        ! was obtained with DSFRK and LDPI will have the right value.
        !
        if (main) then
          if ( mod(nri,2) .eq. 0 ) then
            do iri=1,nri/2
              pi(iri,iri) = pi(iri,iri) + 1d0
              pi(iri+1,iri) = pi(iri+1,iri) + 1d0
            enddo
          else
            pi(1,1) = pi(1,1) + 1d0
            do iri=2,nri/2+1
              pi(iri-1,iri) = pi(iri-1,iri) + 1d0
              pi(iri,iri) = pi(iri,iri) + 1d0
            enddo
          endif
        endif


      !
      ! Case 2. Indefinite factor
      !
      else if (which.eq.'r') then

        !
        ! DSYR2K essentially adds the missing 2d0
        !
        scal = scal/2d0

        !
        ! See how many MO pairs we can do at a time
        !
        ma_dbl = nint(0.9*ma_inquire_avail(mt_dbl))/nri

        if (ma_dbl.eq.0) 
     $    call errquit('gw_buildpi: not enough memory',0,MA_ERR)

        !
        ! If an MPI rank has no work, just zero-out the field and skip
        !
        if (npoles(1).eq.0) then
          call dfill(nri**2,0d0,pi,1)
          goto 200
        endif

        !
        ! Get the number of batches needed to accomodate all info and
        ! get upper and lower limits
        !
        nbatch = npoles(1)/ma_dbl
        if (mod(npoles(1),ma_dbl).ne.0) nbatch = nbatch + 1

        mymax = npoles(1)/nbatch
        if (mod(npoles(1),nbatch).ne.0) mymax = mymax + 1

        mymax = min(npoles(1),mymax)
        ok = ma_push_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

        !
        ! If dorfp=.true., and in house version of DSFRK for 2-rank
        ! updates (DSFR2K_DMR) will be used. The numbers below were obtained
        ! empirically and might not lead to an optimal algorithm in all
        ! computational architechtures.
        !
        dorfp = (nri.le.2048) .and. (mymax.gt.512)

        !
        ! Alpha spin channel
        !
        do ibatch=1,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(1))
          mynpoles = ulpole - llpole + 1

          call dfill(nri*mynpoles,0d0,dbl_mb(k_scr),1)

          !
          ! The combination of OpenMP and OpenBLAS with pthreads might
          ! lead to slow downs when BLAS routines are called inside
          ! OMP regions. For now, we will set BLAS threads to 1.
          !
#ifdef USE_OPENMP
          call util_blas_set_num_threads(1)
#endif
!$omp     parallel do schedule(static) private(idx,idy)
          do ipole=1,mynpoles
            idx = k_scr + (ipole-1)*nri
            idy = llpole + ipole - 1
            call daxpy(nri,factor(idy,1),ovup(1,idy),1,dbl_mb(idx),1)
          enddo
!$omp     end parallel do
#ifdef USE_OPENMP
          call util_blas_set_num_threads(nthreads)
#endif

          !
          ! Use either DSYR2K, DSFR2K, or DGEMM. One of these should be
          ! faster depending on the matrix shape. The current values
          ! have not been optimized across different platforms.
          !
          !
          scal2 = 1d0
          if (ibatch.eq.1) scal2 = 0d0

          if (nri.gt.2048) then
            call dsyr2k('l','n',nri,mynpoles,scal,ovup(1,llpole),nri,
     $                   dbl_mb(k_scr),nri,scal2,pi,nri)
          elseif (dorfp) then
            call dsfr2k_dmr(ovup(1,llpole),dbl_mb(k_scr),pi,nri,
     $                      mynpoles,scal,scal2)
          else
            call dgemm('n','t',nri,nri,mynpoles,2d0*scal,ovup(1,llpole),
     $                  nri,dbl_mb(k_scr),nri,scal2,pi,nri)
          endif

        enddo
        ok = ma_chop_stack(l_scr)

        !
        ! Beta spin channel
        !
  200   if (ipol.gt.1 .and. npoles(2).gt.0) then

          !
          ! Get the number of batches needed to accomodate all info and
          ! get upper and lower limits
          !
          nbatch = npoles(2)/ma_dbl
          if (mod(npoles(2),ma_dbl).ne.0) nbatch = nbatch + 1

          mymax = npoles(2)/nbatch
          if (mod(npoles(2),nbatch).ne.0) mymax = mymax + 1

          mymax = min(npoles(2),mymax)
          ok = ma_push_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

          do ibatch=1,nbatch
            llpole = (ibatch-1)*mymax + 1
            ulpole = min(ibatch*mymax,npoles(2))
            mynpoles = ulpole - llpole + 1

            call dfill(nri*mynpoles,0d0,dbl_mb(k_scr),1)

            !
            ! The combination of OpenMP and OpenBLAS with pthreads might
            ! lead to slow downs when BLAS routines are called inside
            ! OMP regions. For now, we will set BLAS threads to 1.
            !
#ifdef USE_OPENMP
            call util_blas_set_num_threads(1)
#endif
!$omp       parallel do schedule(static)
            do ipole=1,mynpoles
              idx = k_scr + (ipole-1)*nri
              idy = llpole + ipole - 1
              call daxpy(nri,factor(idy,2),ovdw(1,idy),1,dbl_mb(idx),1)
            enddo
!$omp       end parallel do
#ifdef USE_OPENMP
            call util_blas_set_num_threads(nthreads)
#endif

            if (nri.gt.2048) then
              call dsyr2k('l','n',nri,mynpoles,scal,ovdw(1,llpole),nri,
     $                     dbl_mb(k_scr),nri,1d0,pi,nri)
            elseif (dorfp) then
              call dsfr2k_dmr(ovdw(1,llpole),dbl_mb(k_scr),pi,nri,
     $                        mynpoles,scal,1d0)
            else
              call dgemm('n','t',nri,nri,mynpoles,2d0*scal,
     $                    ovdw(1,llpole),nri,dbl_mb(k_scr),nri,1d0,
     $                    pi,nri)
            endif
          enddo
          ok = ma_chop_stack(l_scr)
        endif

        !
        ! Put the matrix in RFP format
        !
        if (.not.dorfp) call std2rfp(pi,nri)

        !
        ! Transform \Pi into the dielectric matrix \epsilon. Since
        ! contributions from all MPI ranks are going to be summed, only
        ! one rank needs to add 1.
        !
        ! Note: Here we employ the shifted RFP format (see std2rfp).
        !
        if (main) then
          if ( mod(nri,2) .eq. 0 ) then
            do iri=1,nri/2
              pi(iri,iri) = pi(iri,iri) + 1d0
              pi(iri,iri+1) = pi(iri,iri+1) + 1d0
            enddo
          else
            pi(1,1) = pi(1,1) + 1d0
            do iri=2,nri/2+1
              pi(iri-1,iri) = pi(iri-1,iri) + 1d0
              pi(iri,iri) = pi(iri,iri) + 1d0
            enddo
          endif
        endif

      endif


      end subroutine

! *********************************************************************

      subroutine gw_cdgw_buildpi_cmplx(which,ovup,ovdw,pi,factor,ldpi,
     $                                 nri,npoles,main,ipol,maxpoles)
C
C     Complex arithmetic version of gw_cdgw_buildpi. Here, we do not
C     need to distinguish between the two cases since they can be easily
C     handled in complex arithmetics. However, we do need to copy the
C     data from doubles to complex arrays.
C
      implicit none
#include "mafdecls.fh"      
      logical,intent(in) :: main
      integer,intent(in) :: ldpi, nri, maxpoles, ipol
      character(len=1),intent(in) :: which
      integer,intent(in),dimension(ipol) :: npoles
      double precision,intent(in),dimension(nri,*) :: ovup,ovdw

      double complex,intent(in),dimension(maxpoles,ipol) :: factor
      double complex,intent(out),dimension(ldpi,*) :: pi

      double complex :: scal, zero, one, fac
      integer ipole,iri,ma_dbl,l_scr,k_scr,mynpoles,llpole,mymax
      integer ulpole,nbatch,ibatch,isp
      logical ok

      zero = dcmplx(0d0,0d0)
      scal = dcmplx(2d0,0d0)
      one  = dcmplx(1d0,0d0)

      if (ipol.eq.1) scal = 2d0*scal

      ma_dbl = ma_inquire_avail(mt_dcpl)/nri
      fac = -one

      do isp=1,ipol

        if (npoles(isp).eq.0) cycle

        !
        ! See how many batches are needed depending on the available
        ! memory
        !
        nbatch = npoles(isp)/ma_dbl + 1
        mymax = min(npoles(isp)/nbatch + 1,npoles(isp))
        ok = ma_alloc_get(mt_dcpl,nri*mymax,'scratch',l_scr,k_scr)
        fac = fac + one

        do ibatch=1,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(isp))
          mynpoles = ulpole - llpole + 1

          !
          ! Copy data from doubles to double complex arrays
          !
          if (isp.eq.1) then
            call zlacp2('A',nri,mynpoles,ovup(1,llpole),nri,
     $                   dcpl_mb(k_scr),nri)
          else
            call zlacp2('A',nri,mynpoles,ovdw(1,llpole),nri,
     $                   dcpl_mb(k_scr),nri)
          endif

          !
          ! Scale ERIs with square roots of the square bracket in Equation
          ! (39).
          !
          do ipole=1,mynpoles
            call zscal(nri,factor(ipole+llpole-1,isp),
     $                 dcpl_mb(k_scr+(ipole-1)*nri),1)
          enddo

          call zsyrk('l','n',nri,mynpoles,scal,dcpl_mb(k_scr),nri,1d0,
     $                pi,nri)

        enddo
        ok = ma_free_heap(l_scr)
      enddo

      !
      ! Transform into dielectric matrix
      !
      if (main) then
        do iri=1,nri
          pi(iri,iri) = pi(iri,iri) + one
        enddo
      endif


      end subroutine

! *********************************************************************

      subroutine dsfr2k_dmr(a,b,c,m,n,alpha,beta)
      !
      ! In-house RFP version of DSYR2K
      !
      ! Note: when the size of the matrix is even, in reality we are
      ! using a shifted RFP. For example, for N=4 the full format and
      ! RFP formats are given by:
      !
      ! [  a11  ---  ---  ---  ]        [  a33  a34  ]
      ! [  a12  a22  ---  ---  ]        [  a11  a44  ]
      ! [  a13  a23  a33  ---  ]        [  a12  a22  ]
      ! [  a14  a24  a34  a44  ]        [  a13  a23  ]
      !                                 [  a14  a24  ]
      !
      ! The format we are obtaining here is
      !
      ! [  a11  a33  a34 ]
      ! [  a12  a22  a44 ]
      ! [  a13  a23  --- ]
      ! [  a14  a24  --- ]

      implicit none
      integer m,n
      double precision,intent(in) :: a(m,n),b(m,n),alpha,beta
      double precision,intent(out) :: c(m,m)

      if (mod(m,2).eq.0) then
        call dsyr2k('l','n',m/2,n,alpha,a,m,b,m,beta,c,m)
        call dsyr2k('u','n',m/2,n,alpha,a(m/2+1,1),m,b(m/2+1,1),m,beta,
     $              c(1,2),m)
        call dgemm('n','t',m/2,m/2,n,2d0*alpha,a(m/2+1,1),m,b,m,beta,
     $              c(m/2+1,1),m)
      else
        call dsyr2k('l','n',m/2+1,n,alpha,a,m,b,m,beta,c,m)
        call dsyr2k('u','n',m/2,n,alpha,a(m/2+2,1),m,b(m/2+2,1),m,beta,
     $              c(1,2),m)
        call dgemm('n','t',m/2,m/2+1,n,2d0*alpha,a(m/2+2,1),m,b,m,beta,
     $             c(m/2+2,1),m)
      endif

      end subroutine

! *********************************************************************

      subroutine std2rfp(matrix,n)

      ! Transform a matrix from full format to RFP format
      !
      ! Note: when the size of the matrix is even, in reality we are
      ! using a shifted RFP. For example, for N=4 the full format and
      ! RFP formats are given by:
      !
      ! [  a11  ---  ---  ---  ]        [  a33  a34  ]
      ! [  a12  a22  ---  ---  ]        [  a11  a44  ]
      ! [  a13  a23  a33  ---  ]        [  a12  a22  ]
      ! [  a14  a24  a34  a44  ]        [  a13  a23  ]
      !                                 [  a14  a24  ]
      !
      ! The format we are using here is
      !
      ! [  a11  a33  a34 ]
      ! [  a12  a22  a44 ]
      ! [  a13  a23  --- ]
      ! [  a14  a24  --- ]

      implicit none
      integer,intent(in) :: n
      double precision   :: matrix(n,n)
      integer i

      if (mod(n,2).eq.0) then
        do i=n/2+1,n
          matrix(i-n/2,i-n/2+1:n/2+1) = matrix(i:,i)
        enddo
      else
        do i=n/2+2,n
          matrix(i-n/2-1,i-n/2:n/2+1) = matrix(i:,i)
        enddo
      endif

      end subroutine
