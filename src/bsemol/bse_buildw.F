      subroutine bse_buildW(pars,wia,w,nmo,nri,maxpoles,ipol)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bse.fh"

      type(bse_params_t) :: pars
      integer nmo, nri, maxpoles, ipol

      character(*), parameter :: pname = 'bse_buildw: '
      
      double precision wia(maxpoles,ipol)
      double precision w(nri,nri)

      integer :: ovlo,ovhi,kov(2),isp
      integer :: ldpi, sizepi
      integer :: alloc, info, ld, ipole
      double precision,allocatable :: factor(:,:)

      allocate(factor(maxpoles,ipol),stat=alloc)
      if (alloc.ne.0)
     &  call errquit(pname//'allocation failed',0,MA_ERR) 
C
C     Decide the leading dimension of the \Pi array according to
C     the rectangular full-packed (RFP) format
C     
      if (mod(nri,2).eq.0) then
        ldpi = nri + 1
      else
        ldpi = nri
      endif
      sizepi = (nri*(nri+1))/2

      ! Get pointer to local integrals
      do isp=1,pars%ipol
        kov(isp) = 1
        if (pars%mynpoles(isp).lt.1) cycle

        ovlo = pars%ovlo(isp)
        ovhi = pars%ovhi(isp)
        call ga_access(pars%g_eriov(1),1,nri,ovlo,ovhi,kov(isp),ld)

!$omp   parallel do
        do ipole=1,pars%mynpoles(isp)
          factor(ipole,isp) = dsqrt(1.0d0/wia(ipole,isp))
        enddo
!$omp   end parallel do
      enddo

      ! Build the Polarizability matrix and transform it to the
      ! dielectric matrix \epsilon = 1 - \pi
      call gw_cdgw_buildpi('w',dbl_mb(kov(1)),dbl_mb(kov(2)),w,factor,
     $                      ldpi,nri,pars%mynpoles,pars%me.eq.0,
     $                      ipol,maxpoles) 
      call ga_dgop(1038,w,sizepi,'+')

      !Factorize dieletric matrix for further use
      call dpftrf('n','l',nri,w,info)
      
      !Deallocate
      deallocate(factor,stat=alloc)
      if (alloc.ne.0)
     &  call errquit(pname//'deallocation failed',0,MA_ERR) 

      !Release GAs
      do isp=1,ipol
        if (pars%mynpoles(isp).lt.1) cycle
        ovlo = pars%ovlo(isp)
        ovhi = pars%ovhi(isp)
        call ga_release_update(pars%g_eriov(isp),1,nri,ovlo,ovhi)
      enddo

      end
