*
* $Id: c_semicore.F,v 1.1 2001-12-13 03:53:31 bylaska Exp $
*

#define TCGMSG

      subroutine c_semicore_init()      
      implicit none

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / c_ccore2 / ncore_sum

*     **** semicore3 common block ****
      integer dncore(2)
      common / c_ccore3 / dncore

*     *** local variables ****
      integer npack0,npsp,nfft3d
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      npsp = ion_nkatm()

      value = MA_alloc_get(mt_dbl,(npsp*npack0*5),
     >                    'ncore',ncore(2),ncore(1))
      value = value.and.
     >         MA_alloc_get(mt_dbl,(npsp),
     >                    'rcore',rcore(2),rcore(1))
      value = value.and.
     >        MA_alloc_get(mt_log,(npsp+1),
     >                    'semicore',semicore(2),semicore(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),
     >                    'ncore_sum',ncore_sum(2),ncore_sum(1))

      value = value.and.
     >        MA_alloc_get(mt_dbl,(nfft3d),
     >                    'dncore',dncore(2),dncore(1))

      call dcopy(npsp*npack0*5, 0.0d0,0,dbl_mb(ncore(1)), 1)
      call dcopy(npsp,        0.0d0,0,dbl_mb(rcore(1)),1)
      call dcopy(npsp,        0.0d0,0,dbl_mb(ncore_sum(1)),1)
      call dcopy(nfft3d,      0.0d0,0,dbl_mb(dncore(1)),1)

      return
      end


      subroutine c_semicore_end()      
      implicit none

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / c_ccore2 / ncore_sum

*     **** semicore3 common block ****
      integer dncore(2)
      common / c_ccore3 / dncore

      logical value

      value = MA_free_heap(ncore(2))
      value = MA_free_heap(rcore(2))
      value = MA_free_heap(semicore(2))
      value = MA_free_heap(ncore_sum(2))
      value = MA_free_heap(dncore(2))

      return
      end


*     ***************************
*     *                         *
*     *      cpsp_semicore       *
*     *                         *
*     ***************************

      logical function cpsp_semicore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore common block ****
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore


      cpsp_semicore = log_mb(semicore(1)+ia)
      return
      end

*     ***************************
*     *                         *
*     *        cpsp_rcore        *
*     *                         *
*     ***************************

      real*8 function cpsp_rcore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore


      cpsp_rcore = dbl_mb(rcore(1)+ia-1)
      return
      end

*     ***************************
*     *                         *
*     *        cpsp_ncore        *
*     *                         *
*     ***************************

      real*8 function cpsp_ncore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / c_ccore2 / ncore_sum


      cpsp_ncore = dbl_mb(ncore_sum(1)+ia-1)
      return
      end


*     **********************************
*     *                         	   *
*     *    c_semicore_density_update     *
*     *                                *
*     **********************************

      subroutine c_semicore_density_update()
      implicit none


#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore

*     **** semicore3 common block ****
      integer dncore(2)
      common / c_ccore3 / dncore

*     *** local variables ***
      logical value
      real*8 scal2
      integer ii,ia
      integer npack0,nfft3d
      integer exi(2),tmp(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega
      external ion_nion,ion_katm
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()
      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)


      call dcopy(nfft3d, 0.0d0,0, dbl_mb(dncore(1)),1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

        if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call cstrfac(ii,dcpl_mb(exi(1)))
           call Cram_c_pack(0,dcpl_mb(exi(1)))

*          **** put sqrt(core-density) at atom position ****
           call Cram_rc_MUL(0,dbl_mb(ncore(1)+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp(1)))
c          call Cram_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Cram_c_unpack(0,dcpl_mb(tmp(1)))
           call C3dB_cr_fft3b(1,dcpl_mb(tmp(1)))

*          ****  square it and add to density ****
           call C3dB_rr_Sqr(1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
           call C3dB_rr_Sum(1,dbl_mb(dncore(1)),
     >                        dcpl_mb(tmp(1)),
     >                        dbl_mb(dncore(1)))
        end if

      end do
      call C3dB_r_SMul(1,scal2*scal2,dbl_mb(dncore(1)),
     >                               dbl_mb(dncore(1)))

      value = MA_pop_stack(tmp(2))
      value = MA_pop_stack(exi(2))

      return
      end


*     ***************************
*     *                         *
*     *    c_semicore_check       *
*     *                         *
*     ***************************

      subroutine c_semicore_check()
      implicit none

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logical semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / c_ccore2 / ncore_sum

*     *** local variables ***
      logical value
      real*8 scal2,dV,sum
      integer ia
      integer npack0,nfft3d
      integer tmp(2)
      integer nx,ny,nz,k

*     **** external functions ****
      integer  ion_nkatm
      real*8   lattice_omega
      external ion_nkatm
      external lattice_omega

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      dV = lattice_omega()/dble(nx*ny*nz)
      
      scal2 = 1.0d0/lattice_omega()
      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)


      do ia=1,ion_nkatm()

        dbl_mb(ncore_sum(1)+ia-1) = 0.0d0
        if (log_mb(semicore(1)+ia)) then

*          **** put sqrt(core-density) at atom position ****
           do k=1,npack0
              dcpl_mb(tmp(1)+k-1) 
     >          = dcmplx(dbl_mb(ncore(1)+(ia-1)*npack0*5+k-1),0.0d0)
           end do
           call Cram_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Cram_c_unpack(0,dcpl_mb(tmp(1)))
           call C3dB_cr_fft3b(1,dcpl_mb(tmp(1)))

*          ****  square it  ****
           call C3dB_rr_Sqr(1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))

*          ****  integrate it ****
           call C3dB_r_dsum(1,dcpl_mb(tmp(1)),sum)
           dbl_mb(ncore_sum(1)+ia-1) = sum*dV
        end if

      end do

      value = MA_pop_stack(tmp(2))
      return
      end



*     **********************************
*     *                         	   *
*     *    c_semicore_density            *
*     *                                *
*     **********************************

      subroutine c_semicore_density(dns)
      implicit none
      real*8 dns(*)


#include "mafdecls.fh"

*     **** semicore3 common block ****
      integer dncore(2)
      common / c_ccore3 / dncore

*     **** local variables ****
      integer nfft3d

      call C3dB_nfft3d(1,nfft3d)
      call dcopy(nfft3d,dbl_mb(dncore(1)),1,dns,1)
      return
      end


