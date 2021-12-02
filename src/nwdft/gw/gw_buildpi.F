      subroutine gw_buildpi(which,ovup,ovdw,pi,factor,ldpi,nri,npoles,
     $                      main,ipol,maxpoles)
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

      double precision scal, t0
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
      scal = 2d0
      if (ipol.eq.1) scal = 2d0*scal

      if (which.eq.'w') then

        if (npoles(1).gt.0) then

          t0 = util_wallsec()
!$omp     parallel do simd schedule(static) collapse(2)
          do ipole=1,npoles(1)
          do iri=1,nri
c            call dscal(nri,factor(ipole,1),ovup(1,ipole),1)
            ovup(iri,ipole) = ovup(iri,ipole)*factor(ipole,1)
          enddo
          enddo
!$omp     end parallel do simd

          call dsfrk('n','l','n',nri,npoles(1),scal,ovup,nri,0d0,pi)

!$omp     parallel do simd schedule(static) collapse(2)
          do ipole=1,npoles(1)
          do iri=1,nri
c            call dscal(nri,1d0/factor(ipole,1),ovup(1,ipole),1)
            ovup(iri,ipole) = ovup(iri,ipole)/factor(ipole,1)
          enddo
          enddo
!$omp     end parallel do simd

        endif

        if (ipol.gt.1.and.npoles(2).gt.0) then

!$omp     parallel do schedule(static)
          do ipole=1,npoles(2)
            call dscal(nri,factor(ipole,2),ovdw(1,ipole),1)
          enddo
!$omp     end parallel do 

          call dsfrk('n','l','n',nri,npoles(2),scal,ovdw,nri,1d0,pi)

!$omp     parallel do schedule(static)
          do ipole=1,npoles(2)
            call dscal(nri,1d0/factor(ipole,2),ovdw(1,ipole),1)
          enddo
!$omp     end parallel do

        endif

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

      else if (which.eq.'r') then

        scal = scal/2d0
        ma_dbl = nint(0.9*ma_inquire_avail(mt_dbl))/nri

        if (ma_dbl.eq.0) 
     $    call errquit('gw_buildpi: not enough memory',0,MA_ERR)

        if (npoles(1).eq.0) then
          call dfill(nri**2,0d0,pi,1)
          goto 200
        endif

        nbatch = npoles(1)/ma_dbl
        if (mod(npoles(1),ma_dbl).ne.0) nbatch = nbatch + 1

        mymax = npoles(1)/nbatch
        if (mod(npoles(1),nbatch).ne.0) mymax = mymax + 1

        mymax = min(npoles(1),mymax)
        ok = ma_push_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

        dorfp = (nri.le.2048) .and. (mymax.gt.512)

!$omp   parallel do schedule(static) private(idx)
        do ipole=1,mymax
          idx = k_scr + (ipole-1)*nri
          call dfill(nri,0d0,dbl_mb(idx),1)
          call daxpy(nri,factor(ipole,1),ovup(1,ipole),1,dbl_mb(idx),1)
        enddo
!$omp   end parallel do

        if (nri.gt.2048) then
          call dsyr2k('l','n',nri,mymax,scal,ovup,nri,dbl_mb(k_scr),
     $                 nri,0d0,pi,nri)     
        elseif (dorfp) then
          call dsfr2k_(ovup,dbl_mb(k_scr),pi,nri,mymax,scal,0d0)
        else
          call dgemm('n','t',nri,nri,mymax,2d0*scal,ovup,nri,
     $                dbl_mb(k_scr),nri,0d0,pi,nri)     
        endif


        do ibatch=2,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(1))
          mynpoles = ulpole - llpole + 1

!$omp     parallel do schedule(static) private(idx,idy)
          do ipole=1,mynpoles
            idx = k_scr + (ipole-1)*nri
            idy = llpole + ipole - 1
            call dfill(nri,0d0,dbl_mb(idx),1)
            call daxpy(nri,factor(idy,1),ovup(1,idy),1,dbl_mb(idx),1)
          enddo
!$omp     end parallel do

          if (nri.gt.2048) then
            call dsyr2k('l','n',nri,mynpoles,scal,ovup(1,llpole),nri,
     $                   dbl_mb(k_scr),nri,1d0,pi,nri)     
          elseif (dorfp) then
            call dsfr2k_(ovup(1,llpole),dbl_mb(k_scr),pi,nri,mynpoles,
     $                  scal,1d0)
          else
            call dgemm('n','t',nri,nri,mynpoles,2d0*scal,ovup(1,llpole),
     $                  nri,dbl_mb(k_scr),nri,1d0,pi,nri)
          endif

        enddo

        ok = ma_chop_stack(l_scr)

  200   if (ipol.gt.1 .and. npoles(2).gt.0) then

          nbatch = npoles(2)/ma_dbl + 1
          mymax = min(npoles(2)/nbatch + 1,npoles(2))
          ok = ma_push_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

!$omp     parallel do schedule(static)
          do ipole=1,mymax
            call dcopy(nri,ovdw(1,ipole),1,
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
            call dscal(nri,factor(ipole,2),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
!$omp     end parallel do

          if (nri.gt.2048) then
            call dsyr2k('l','n',nri,mymax,scal,ovdw,nri,dbl_mb(k_scr),
     $                   nri,1d0,pi,nri)
          elseif (dorfp) then
            call dsfr2k_(ovdw,dbl_mb(k_scr),pi,nri,mymax,scal,1d0)
          else
            call dgemm('n','t',nri,nri,mymax,2d0*scal,ovdw,nri,
     $                  dbl_mb(k_scr),nri,1d0,pi,nri)
          endif

          do ibatch=2,nbatch
            llpole = (ibatch-1)*mymax + 1
            ulpole = min(ibatch*mymax,npoles(2))
            mynpoles = ulpole - llpole + 1

!$omp       parallel do schedule(static)
            do ipole=1,mynpoles
              call dcopy(nri,ovdw(1,llpole+ipole-1),1,
     $                   dbl_mb(k_scr+(ipole-1)*nri),1)
              call dscal(nri,factor(llpole+ipole-1,2),
     $                   dbl_mb(k_scr+(ipole-1)*nri),1)
            enddo
!$omp       end parallel do

            if (nri.gt.2048) then
              call dsyr2k('l','n',nri,mynpoles,scal,ovdw(1,llpole),nri,
     $                     dbl_mb(k_scr),nri,1d0,pi,nri)     
            elseif (dorfp) then
              call dsfr2k_(ovdw(1,llpole),dbl_mb(k_scr),pi,nri,mynpoles,
     $                    scal,1d0)
            else
              call dgemm('n','t',nri,nri,mynpoles,2d0*scal,
     $                    ovdw(1,llpole),nri,dbl_mb(k_scr),nri,1d0,
     $                    pi,nri)
            endif

          enddo

          ok = ma_chop_stack(l_scr)

        endif

        if (.not.dorfp) call std2rfp(pi,nri)

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



      subroutine gw_buildpi_cmplx(which,ovup,ovdw,pi,factor,ldpi,nri,
     $                      npoles,main,ipol,maxpoles)
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

        nbatch = npoles(isp)/ma_dbl + 1
        mymax = min(npoles(isp)/nbatch + 1,npoles(isp))
        ok = ma_alloc_get(mt_dcpl,nri*mymax,'scratch',l_scr,k_scr)
        fac = fac + one

        if (isp.eq.1) then
          call zlacp2('A',nri,mymax,ovup,nri,dcpl_mb(k_scr),nri)
        else
          call zlacp2('A',nri,mymax,ovdw,nri,dcpl_mb(k_scr),nri)
        endif

        do ipole=1,mymax
          call zscal(nri,factor(ipole,isp),dcpl_mb(k_scr+(ipole-1)*nri),
     $               1)
        enddo
        call zsyrk('l','n',nri,mymax,scal,dcpl_mb(k_scr),nri,fac,pi,nri)

        do ibatch=2,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(isp))
          mynpoles = ulpole - llpole + 1

          if (isp.eq.1) then
            call zlacp2('A',nri,mynpoles,ovup(1,llpole),nri,
     $                  dcpl_mb(k_scr),nri)
          else
            call zlacp2('A',nri,mynpoles,ovdw(1,llpole),nri,
     $                   dcpl_mb(k_scr),nri)
          endif

          do ipole=1,mynpoles
              call zscal(nri,factor(llpole+ipole-1,isp),
     $               dcpl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call zsyrk('l','n',nri,mynpoles,scal,dcpl_mb(k_scr),nri,one,
     $               pi,nri)
        enddo
        ok = ma_free_heap(l_scr)
      enddo

      if (main) then
        do iri=1,nri
          pi(iri,iri) = pi(iri,iri) + one
        enddo
      endif


      end subroutine


      subroutine gw_scal(n,from,to,scal)
      integer n
      double precision from(n),to(n),scal

      to(1:n) = scal*from(1:n)


      end subroutine

      subroutine dsfr2k_(a,b,c,m,n,alpha,beta)
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

      subroutine std2rfp(matrix,n)
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
