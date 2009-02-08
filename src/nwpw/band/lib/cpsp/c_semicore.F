*
* $Id: c_semicore.F,v 1.8 2009-02-08 03:26:25 bylaska Exp $
*

#define TCGMSG

      subroutine c_semicore_init()      
      implicit none

#include "mafdecls.fh"
#include "c_semicore_common.fh"

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
#include "c_semicore_common.fh"


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
#include "c_semicore_common.fh"


      cpsp_semicore = log_mb(semicore(1)+ia)
      return
      end

*     ***************************
*     *                         *
*     *        cpsp_rcore       *
*     *                         *
*     ***************************

      real*8 function cpsp_rcore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"
#include "c_semicore_common.fh"

      cpsp_rcore = dbl_mb(rcore(1)+ia-1)
      return
      end

*     ***************************
*     *                         *
*     *        cpsp_ncore       *
*     *                         *
*     ***************************

      real*8 function cpsp_ncore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"
#include "c_semicore_common.fh"

      cpsp_ncore = dbl_mb(ncore_sum(1)+ia-1)
      return
      end


*     **********************************
*     *                                *
*     *    c_semicore_density_update   *
*     *                                *
*     **********************************

      subroutine c_semicore_density_update()
      implicit none
#include "errquit.fh"


#include "mafdecls.fh"
#include "c_semicore_common.fh"

*     *** local variables ***
      logical value
      real*8 scal2,scal
      integer ii,ia
      integer npack0,nfft3d
      integer exi(2),tmp(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega
      external ion_nion,ion_katm
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()
      scal = scal2*scal2
      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      call dcopy(nfft3d, 0.0d0,0, dbl_mb(dncore(1)),1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

        if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call cstrfac(ii,dcpl_mb(exi(1)))
           call Cram_c_pack(0,dcpl_mb(exi(1)))

*          **** put sqrt(core-density) at atom position ****
           call Cram_rc_Mul(0,dbl_mb(ncore(1)+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp(1)))
c          call Cram_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Cram_c_unpack(0,dcpl_mb(tmp(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(tmp(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(tmp(1)))

*          ****  square it and add to density ****
c          call C3dB_cr_Sqr(1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
c          call C3dB_rr_Sum(1,dbl_mb(dncore(1)),
c    >                        dcpl_mb(tmp(1)),
c    >                        dbl_mb(dncore(1)))
           call C3dB_cr_aSqrpy(1,scal,dcpl_mb(tmp(1)),
     >                                 dbl_mb(dncore(1)))
        end if

      end do
c     call C3dB_r_SMul(1,scal2*scal2,dbl_mb(dncore(1)),
c    >                               dbl_mb(dncore(1)))

      value = MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) call errquit('c_semicore:pop stack',1,MA_ERR)

      return
      end


*     ***************************
*     *                         *
*     *    c_semicore_check     *
*     *                         *
*     ***************************

      subroutine c_semicore_check()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"
#include "c_semicore_common.fh"

*     *** local variables ***
      logical value
      real*8 scal2,dV,sum
      integer ia
      integer npack0,nfft3d
      integer tmp(2),tmp2(2)
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
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)


      do ia=1,ion_nkatm()

        dbl_mb(ncore_sum(1)+ia-1) = 0.0d0
        if (log_mb(semicore(1)+ia)) then

*          **** put sqrt(core-density) at atom position ****
           do k=1,npack0
              dcpl_mb(tmp(1)+k-1) 
     >          = dcmplx(dbl_mb(ncore(1)+(ia-1)*npack0*5+k-1),0.0d0)
           end do
c           call Cram_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
           call Cram_c_SMul1(0,scal2,dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Cram_c_unpack(0,dcpl_mb(tmp(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(tmp(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(tmp(1)))

*          ****  square it  ****
           call C3dB_cr_Sqr(1,dcpl_mb(tmp(1)),dbl_mb(tmp2(1)))

*          ****  integrate it ****
           call C3dB_r_dsum(1,dbl_mb(tmp2(1)),sum)
           dbl_mb(ncore_sum(1)+ia-1) = sum*dV
        end if

      end do

      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp(2))
      if (.not.value) call errquit('c_semicore_check:pop',1,MA_ERR)
      return
      end



*     **********************************
*     *                                *
*     *    c_semicore_density          *
*     *                                *
*     **********************************

      subroutine c_semicore_density(dns)
      implicit none
      real*8 dns(*)


#include "mafdecls.fh"
#include "c_semicore_common.fh"

*     **** local variables ****
      integer nfft3d

      call C3dB_nfft3d(1,nfft3d)
      call dcopy(nfft3d,dbl_mb(dncore(1)),1,dns,1)
      return
      end




*     **********************************
*     *                                *
*     *        c_semicore_euv          *
*     *                                *
*     **********************************

      subroutine c_semicore_euv(stress)
      implicit none
      real*8 stress(3,3)


#include "mafdecls.fh"
#include "errquit.fh"
#include "c_semicore_common.fh"


*     **** local variables ****
      logical value
      integer u,v,s,nfft3d
      real*8  VV,WW(3,3),hm(3,3),tstress(3,3)
      real*8  wus,pi,scal
      integer Aus(2)

*     **** external functions ****
      real*8   c_electron_pxc_rho,lattice_unitg
      external c_electron_pxc_rho,lattice_unitg


      call C3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dbl,nfft3d,'Aus',Aus(2),Aus(1))
      if (.not. value) 
     >  call errquit('c_semicore_euv:out of stack memory',0, MA_ERR)


*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     **** compute V ****
      VV = -c_electron_pxc_rho(dbl_mb(dncore(1)))

*     **** compute W ****
      do u=1,3
      do s=u,3
        call c_semicore_gen_Aus(s,u,dbl_mb(Aus(1)))

        wus = -c_electron_pxc_rho(dbl_mb(Aus(1)))
        WW(u,s) = wus
        WW(s,u) = wus 
      end do
      end do

*     **** stress(u,v) = V*hm(u,v) + Sum(s){W(u,s)*hm(s,v) }  ****   
      do v=1,3
      do u=1,3
        stress(u,v) = VV*hm(u,v)
        do s=1,3
           stress(u,v) = stress(u,v) + WW(u,s)*hm(s,v)
        end do
      end do
      end do

      value = MA_pop_stack(Aus(2))
      if (.not.value) call errquit('c_semicore_euv: pop stack',1,MA_ERR)
      return
      end


*     **********************************
*     *                                *
*     *      c_semicore_gen_Aus        *
*     *                                *
*     **********************************

      subroutine c_semicore_gen_Aus(u,s,Aus)
      implicit none
      integer u,s
      real*8 Aus(*)

#include "mafdecls.fh" 
#include "errquit.fh"
#include "c_semicore_common.fh"


*     **** coulomb common block ****
      integer vc_indx,vc_hndl
      common / c_vc_block / vc_indx,vc_hndl


*     *** local variables ***
      logical value
      real*8 scal2,fourpi,ss,sum
      integer ii,ia
      integer npack0,nfft3d
      integer exi(2),tmp1(2),tmp2(2),tmp3(2),Gu,Gs,rmp1(2),rmp2(2)

*     **** external functions ****
      integer  ion_nion,ion_katm,c_G_indx
      real*8   lattice_omega
      external ion_nion,ion_katm,c_G_indx
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()
      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)

  
*     **** allocate memory off of MA stack ****
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'rmp1',rmp1(2),rmp1(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nfft3d,'rmp2',rmp2(2),rmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,npack0,'tmp3',tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('semicore_gen_Aus:out of stack memory',0, MA_ERR)

*     **** define Gu and Gs ****
      Gu = tmp1(1)
      Gs = tmp2(1)
      call dcopy(nfft3d,dbl_mb(c_G_indx(u)),1,dcpl_mb(Gu),1)
      call dcopy(nfft3d,dbl_mb(c_G_indx(s)),1,dcpl_mb(Gs),1)
      call Cram_r_pack(0,dcpl_mb(Gu))
      call Cram_r_pack(0,dcpl_mb(Gs))


*     **** tmp3(G) = Gu*Gs/G **** 
      fourpi = 16.0d0*datan(1.0d0)
      ss     = 1.0d0/fourpi
      call Cram_r_SMul(0,ss,dbl_mb(vc_indx),dbl_mb(tmp3(1)))
c      call Cram_rr_Sqrt(0,dbl_mb(tmp3(1)),dbl_mb(tmp3(1)))
c      call Cram_rr_Mul(0,dbl_mb(tmp3(1)),
c     >                   dcpl_mb(Gu),
c     >                   dbl_mb(tmp3(1)))
c      call Cram_rr_Mul(0,dbl_mb(tmp3(1)),
c     >                   dcpl_mb(Gs),
c     >                   dbl_mb(tmp3(1)))
      call Cram_rr_Sqrt1(0,dbl_mb(tmp3(1)))
      call Cram_rr_Mul2(0,dcpl_mb(Gu),dbl_mb(tmp3(1)))
      call Cram_rr_Mul2(0,dcpl_mb(Gs),dbl_mb(tmp3(1)))


*     **** calculate Aus ****  
      call dcopy(nfft3d,0.0d0,0,Aus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

        if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call cstrfac_pack(0,ii,dcpl_mb(exi(1)))

*          **** put phi and dphi/dG  at atom position ****
           call Cram_rc_Mul(0,dbl_mb(ncore(1)+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp1(1)))

           call Cram_rc_Mul(0,dbl_mb(ncore(1)+npack0+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp2(1)))
c           call Cram_rc_Mul(0,dbl_mb(tmp3(1)),
c     >                        dcpl_mb(tmp2(1)),
c     >                        dcpl_mb(tmp2(1)))
           call Cram_rc_Mul2(0,dbl_mb(tmp3(1)),dcpl_mb(tmp2(1)))


*          **** Put put tmp1 and tmp2 into real space ****
           call Cram_c_unpack(0,dcpl_mb(tmp1(1)))
           call Cram_c_unpack(0,dcpl_mb(tmp2(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
           !call C3dB_cr_fft3b(1,dcpl_mb(tmp2(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(tmp1(1)))
           call C3dB_cr_pfft3b(1,0,dcpl_mb(tmp2(1)))

           call C3dB_cr_real(1,dcpl_mb(tmp1(1)),dbl_mb(rmp1(1)))
           call C3dB_cr_real(1,dcpl_mb(tmp2(1)),dbl_mb(rmp2(1)))


*          ****  add tmp1*tmp2=phi*dphi*dphi/dG density ****
c           call C3dB_rr_Mul(1,dbl_mb(rmp1(1)),
c     >                        dbl_mb(rmp2(1)),
c     >                        dbl_mb(rmp1(1)))
c           call C3dB_rr_Sum(1,Aus,
c     >                        dbl_mb(rmp1(1)),
c     >                        Aus)
           call C3dB_rr_Mul2(1,dbl_mb(rmp2(1)),dbl_mb(rmp1(1)))
           call C3dB_rr_Sum2(1,dbl_mb(rmp1(1)),Aus)
        end if

      end do
c      call C3dB_r_SMul(1,scal2*scal2,Aus,Aus)
      call C3dB_r_SMul1(1,scal2*scal2,Aus)


*     **** pop stack memory ****
      value =           MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(rmp2(2))
      value = value.and.MA_pop_stack(rmp1(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) 
     > call errquit('c_semicore_gen_Aus:popping stack memory',2,MA_ERR)

      return
      end
