      subroutine gw_buildpi(which,ovup,ovdw,pi,factor,ldpi,nri,npoles,
     $                      main,ipol,maxpoles)
      implicit none
#include "mafdecls.fh"      
#include "errquit.fh"
      logical,intent(in) :: main
      integer,intent(in) :: ldpi, nri, maxpoles, ipol
      character(len=1),intent(in) :: which
      integer,intent(in),dimension(ipol) :: npoles
      double precision,intent(inout),dimension(nri,*) :: ovup,ovdw

      double precision,intent(in),dimension(maxpoles,ipol) :: factor
      double precision,intent(out),dimension(ldpi,*) :: pi

      double precision scal
      integer ipole,iri,ma_dbl,l_scr,k_scr,mynpoles,llpole,mymax
      integer ulpole,nbatch,ibatch
      logical ok

      scal = 2d0
      if (ipol.eq.1) scal = 2d0*scal

      if (which.eq.'w') then

        if (npoles(1).gt.0) then

!$omp     parallel do
          do ipole=1,npoles(1)
            ovup(1:nri,ipole) = ovup(1:nri,ipole)*factor(ipole,1)
          enddo
!$omp     end parallel do

          call dsfrk('n','l','n',nri,npoles(1),scal,ovup,nri,0d0,pi)

!$omp     parallel do
          do ipole=1,npoles(1)
            ovup(1:nri,ipole) = ovup(1:nri,ipole)/factor(ipole,1)
          enddo
!$omp     end parallel do

        endif

        if (ipol.gt.1.and.npoles(2).gt.0) then

!$omp     parallel do
          do ipole=1,npoles(2)
            ovdw(1:nri,ipole) = ovdw(1:nri,ipole)*factor(ipole,2)
          enddo
!$omp     end parallel do

          call dsfrk('n','l','n',nri,npoles(2),scal,ovdw,nri,1d0,pi)

!$omp     parallel do
          do ipole=1,npoles(2)
            ovdw(1:nri,ipole) = ovdw(1:nri,ipole)/factor(ipole,2)
          enddo
!$omp     end parallel do

        endif

        if (main) then
          if ( mod(nri,2) .eq. 0 ) then
!$omp       parallel do
            do iri=1,nri/2
              pi(iri,iri) = pi(iri,iri) + 1d0
              pi(iri+1,iri) = pi(iri+1,iri) + 1d0
            enddo
!$omp       end parallel do
          else
            pi(1,1) = pi(1,1) + 1d0
!$omp       parallel do
            do iri=2,nri/2+1
              pi(iri-1,iri) = pi(iri-1,iri) + 1d0
              pi(iri,iri) = pi(iri,iri) + 1d0
            enddo
!$omp       end parallel do
          endif
        endif

      else if (which.eq.'r') then

        scal = scal/2d0

        if (npoles(1).eq.0) goto 200

        ma_dbl = ma_inquire_avail(mt_dbl)/nri

        if (ma_dbl.eq.0) 
     $    call errquit('gw_buildpi: not enough memory',0,MA_ERR)

        nbatch = npoles(1)/ma_dbl
        if (mod(npoles(1),ma_dbl).ne.0) nbatch = nbatch + 1

        mymax = npoles(1)/nbatch
        if (mod(npoles(1),nbatch).ne.0) mymax = mymax + 1

        mymax = min(npoles(1),mymax)
        ok = ma_alloc_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

        call dcopy(nri*mymax,ovup,1,dbl_mb(k_scr),1)
        do ipole=1,mymax
          call dscal(nri,factor(ipole,1),dbl_mb(k_scr+(ipole-1)*nri),1)
        enddo
        call dgemm('n','t',nri,nri,mymax,2d0*scal,dbl_mb(k_scr),nri,
     $              ovup,nri,0d0,pi,nri)


        do ibatch=2,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(1))
          mynpoles = ulpole - llpole + 1
          call dcopy(nri*mynpoles,ovup(1,llpole),1,dbl_mb(k_scr),1)
          do ipole=1,mynpoles
            call dscal(nri,factor(llpole+ipole-1,1),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call dgemm('n','t',nri,nri,mynpoles,2d0*scal,dbl_mb(k_scr),
     $                nri,ovup(1,llpole),nri,1d0,pi,nri)
          
        enddo

        ok = ma_free_heap(l_scr)

  200   if (ipol.gt.1 .and. npoles(2).gt.0) then
          nbatch = npoles(2)/ma_dbl + 1
          mymax = min(npoles(2)/nbatch + 1,npoles(2))
          ok = ma_alloc_get(mt_dbl,nri*mymax,'scratch',l_scr,k_scr)

          call dcopy(nri*mymax,ovdw,1,dbl_mb(k_scr),1)
          do ipole=1,mymax
            call dscal(nri,factor(ipole,2),
     $                 dbl_mb(k_scr+(ipole-1)*nri),1)
          enddo
          call dsyr2k('l','n',nri,mymax,scal,ovdw,nri,dbl_mb(k_scr),nri,
     $                 1d0,pi,nri)

          do ibatch=2,nbatch
            llpole = (ibatch-1)*mymax + 1
            ulpole = min(ibatch*mymax,npoles(2))
            mynpoles = ulpole - llpole + 1
            call dcopy(nri*mynpoles,ovdw(1,llpole),1,dbl_mb(k_scr),1)
            do ipole=1,mynpoles
              call dscal(nri,factor(llpole+ipole-1,2),
     $                   dbl_mb(k_scr+(ipole-1)*nri),1)
            enddo
            call dsyr2k('l','n',nri,mynpoles,scal,ovdw(1,llpole),nri,
     $                   dbl_mb(k_scr),nri,1d0,pi,nri)     
          enddo

          ok = ma_free_heap(l_scr)

        endif

        if (main) then
!$omp     parallel do
          do iri=1,nri
            pi(iri,iri) = pi(iri,iri) + 1d0
          enddo
!$omp     end parallel do
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
      double precision,intent(inout),dimension(nri,*) :: ovup,ovdw

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
          call zscal(nri,sqrt(factor(ipole,isp)),
     $               dcpl_mb(k_scr+(ipole-1)*nri),1)
        enddo
        call zsyrk('l','n',nri,mymax,scal,dcpl_mb(k_scr),nri,fac,pi,nri)

        do ibatch=2,nbatch
          llpole = (ibatch-1)*mymax + 1
          ulpole = min(ibatch*mymax,npoles(isp))
          mynpoles = ulpole - llpole + 1

          if (isp.eq.1) then
            call zlacp2('A',nri,mynpoles,ovup,nri,dcpl_mb(k_scr),nri)
          else
            call zlacp2('A',nri,mynpoles,ovdw,nri,dcpl_mb(k_scr),nri)
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
!$omp   parallel do
        do iri=1,nri
          pi(iri,iri) = pi(iri,iri) + one
        enddo
!$omp   end parallel do
      endif


      end subroutine
