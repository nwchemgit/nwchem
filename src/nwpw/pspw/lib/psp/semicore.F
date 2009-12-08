*
* $Id: semicore.F,v 1.10 2009-02-07 03:50:58 bylaska Exp $
*

#define TCGMSG

*     **********************************
*     *                                *
*     *        semicore_init           *
*     *                                *
*     **********************************

      subroutine semicore_init()      
      implicit none

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / ccore2 / ncore_sum

*     **** semicore3 common block ****
      integer dncore(2)
      common / ccore3 / dncore

*     *** local variables ****
      integer npack0,npsp,n2ft3d
      logical value

*     **** external functions *****
      integer  ion_nkatm
      external ion_nkatm
 

      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)
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
     >        MA_alloc_get(mt_dbl,(n2ft3d),
     >                    'dncore',dncore(2),dncore(1))

      call dcopy(npsp*npack0*5, 0.0d0,0,dbl_mb(ncore(1)), 1)
      call dcopy(npsp,        0.0d0,0,dbl_mb(rcore(1)),1)
      call dcopy(npsp,        0.0d0,0,dbl_mb(ncore_sum(1)),1)
      call dcopy(n2ft3d,      0.0d0,0,dbl_mb(dncore(1)),1)

      return
      end

*     **********************************
*     *                                *
*     *        semicore_end            *
*     *                                *
*     **********************************

      subroutine semicore_end()      
      implicit none

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / ccore2 / ncore_sum

*     **** semicore3 common block ****
      integer dncore(2)
      common / ccore3 / dncore

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
*     *      psp_semicore       *
*     *                         *
*     ***************************

      logical function psp_semicore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore


      psp_semicore = log_mb(semicore(1)+ia)
      return
      end

*     ***************************
*     *                         *
*     *        psp_rcore        *
*     *                         *
*     ***************************

      real*8 function psp_rcore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore


      psp_rcore = dbl_mb(rcore(1)+ia-1)
      return
      end

*     ***************************
*     *                         *
*     *        psp_ncore        *
*     *                         *
*     ***************************

      real*8 function psp_ncore(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / ccore2 / ncore_sum


      psp_ncore = dbl_mb(ncore_sum(1)+ia-1)
      return
      end


*     **********************************
*     *                                *
*     *    semicore_density_update     *
*     *                                *
*     **********************************

      subroutine semicore_density_update()
      implicit none
#include "errquit.fh"


#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** semicore3 common block ****
      integer dncore(2)
      common / ccore3 / dncore

*     *** local variables ***
      logical value
      real*8 scal2
      integer ii,ia
      integer npack0,nfft3d,n2ft3d
      integer exi(2),tmp(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      real*8   lattice_omega
      external ion_nion,ion_katm
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) 
     >  call errquit('semicore_density_update: push stack',0, MA_ERR)

      
      call dcopy(n2ft3d, 0.0d0,0, dbl_mb(dncore(1)),1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

        if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
c           call strfac(ii,dcpl_mb(exi(1)))
c           call Pack_c_pack(0,dcpl_mb(exi(1)))
           call strfac_pack(0,ii,dcpl_mb(exi(1)))


*          **** put sqrt(core-density) at atom position ****
           call Pack_tc_Mul(0,dbl_mb(ncore(1)+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp(1)))
c          call Pack_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Pack_c_unpack(0,dcpl_mb(tmp(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp(1)))

*          ****  square it and add to density ****
c           call D3dB_rr_Sqr(1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
c           call D3dB_rr_Sum(1,dbl_mb(dncore(1)),
c     >                        dcpl_mb(tmp(1)),
c     >                        dbl_mb(dncore(1)))
           call D3dB_rr_Sqr1(1,dcpl_mb(tmp(1)))
           call D3dB_rr_Sum2(1,dcpl_mb(tmp(1)),dbl_mb(dncore(1)))
        end if

      end do
      call D3dB_r_Zero_Ends(1,dbl_mb(dncore(1)))
c      call D3dB_r_SMul(1,scal2*scal2,dbl_mb(dncore(1)),
c     >                               dbl_mb(dncore(1)))
      call D3dB_r_SMul1(1,scal2*scal2,dbl_mb(dncore(1)))

      value =           MA_pop_stack(tmp(2))
      value = value.and.MA_pop_stack(exi(2))
      if (.not.value) 
     >  call errquit('semicore_density_update: pop stack',1, MA_ERR)
      return
      end


*     ***************************
*     *                         *
*     *    semicore_check       *
*     *                         *
*     ***************************

      subroutine semicore_check()
      implicit none
#include "errquit.fh"

#include "mafdecls.fh"

*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logical semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** semicore2 common block ****
      integer ncore_sum(2)
      common / ccore2 / ncore_sum

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


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dV = lattice_omega()/dble(nx*ny*nz)
      
      scal2 = 1.0d0/lattice_omega()
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))
      if (.not. value) 
     >  call errquit('semicore_check:out of stack memory',0, MA_ERR)

      do ia=1,ion_nkatm()

        dbl_mb(ncore_sum(1)+ia-1) = 0.0d0
        if (log_mb(semicore(1)+ia)) then

*          **** put sqrt(core-density) at atom position ****
           do k=1,npack0
              dcpl_mb(tmp(1)+k-1) 
     >          = dcmplx(dbl_mb(ncore(1)+(ia-1)*npack0*5+k-1),0.0d0)
           end do
c           call Pack_c_SMul(0,scal2,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
           call Pack_c_SMul1(0,scal2,dcpl_mb(tmp(1)))

*          **** Put put tmp into real space ****
           call Pack_c_unpack(0,dcpl_mb(tmp(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp(1)))

*          ****  square it  ****
c           call D3dB_rr_Sqr(1,dcpl_mb(tmp(1)),dcpl_mb(tmp(1)))
           call D3dB_rr_Sqr1(1,dcpl_mb(tmp(1)))

*          ****  integrate it ****
           call D3dB_r_dsum(1,dcpl_mb(tmp(1)),sum)
           dbl_mb(ncore_sum(1)+ia-1) = sum*dV

        end if

      end do

      value = MA_pop_stack(tmp(2))
      if (.not.value) call errquit('semicore_check: pop stack',1,
     &       MA_ERR)
      return
      end



*     **********************************
*     *                                *
*     *    semicore_density            *
*     *                                *
*     **********************************

      subroutine semicore_density(dns)
      implicit none
      real*8 dns(*)


#include "mafdecls.fh"

*     **** semicore3 common block ****
      integer dncore(2)
      common / ccore3 / dncore

*     **** local variables ****
      integer n2ft3d

      call D3dB_n2ft3d(1,n2ft3d)
      call dcopy(n2ft3d,dbl_mb(dncore(1)),1,dns,1)
      return
      end



*     **********************************
*     *                                *
*     *        semicore_euv            *
*     *                                *
*     **********************************

      subroutine semicore_euv(stress)
      implicit none
#include "errquit.fh"
      real*8 stress(3,3)


#include "mafdecls.fh"

*     **** semicore3 common block ****
      integer dncore(2)
      common / ccore3 / dncore

*     **** local variables ****
      logical value
      integer u,v,s,n2ft3d
      real*8  VV,WW(3,3),hm(3,3),tstress(3,3)
      real*8  wus,pi,scal
      integer Aus(2)

*     **** external functions ****
      real*8   electron_pxc_rho,lattice_unitg
      external electron_pxc_rho,lattice_unitg


      call D3dB_n2ft3d(1,n2ft3d)
      value = MA_push_get(mt_dbl,n2ft3d,'Aus',Aus(2),Aus(1))
      if (.not. value) 
     >  call errquit('semicore_euv:out of stack memory',0, MA_ERR)



*     *** define hm ****
      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)
      do v=1,3
      do u=1,3
         hm(u,v) = scal*lattice_unitg(u,v)
      end do
      end do

*     **** compute V ****
      VV = -electron_pxc_rho(dbl_mb(dncore(1)))


*     **** compute W ****
      do u=1,3
      do s=u,3
        call semicore_gen_Aus(s,u,dbl_mb(Aus(1)))

        wus = -electron_pxc_rho(dbl_mb(Aus(1)))
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
      if (.not.value) call errquit('semicore_euv: pop stack',1, MA_ERR)
      return
      end


*     **********************************
*     *                                *
*     *        semicore_gen_Aus        *
*     *                                *
*     **********************************

      subroutine semicore_gen_Aus(u,s,Aus)
      implicit none
#include "errquit.fh"
      integer u,s
      real*8 Aus(*)

#include "mafdecls.fh" 

*     **** semicore common block ****
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / ccore / ncore,rcore,semicore

*     **** coulomb common block ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     *** local variables ***
      logical value
      real*8 scal2,fourpi,ss,sum
      integer ii,ia
      integer npack0,nfft3d,n2ft3d
      integer exi(2),tmp1(2),tmp2(2),tmp3(2),Gu,Gs

*     **** external functions ****
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega
      external ion_nion,ion_katm,G_indx
      external lattice_omega

      scal2 = 1.0d0/lattice_omega()
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

  
*     **** allocate memory off of MA stack ****
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp2',tmp2(2),tmp2(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'exi',exi(2),exi(1))
      value = value.and.
     >        MA_push_get(mt_dbl,npack0,'tmp3',tmp3(2),tmp3(1))
      if (.not. value) 
     >   call errquit('semicore_gen_Aus:out of stack memory',0, MA_ERR)

*     **** define Gu and Gs ****
      Gu = tmp1(1)
      Gs = tmp2(1)
      call dcopy(nfft3d,dbl_mb(G_indx(u)),1,dcpl_mb(Gu),1)
      call dcopy(nfft3d,dbl_mb(G_indx(s)),1,dcpl_mb(Gs),1)
      call Pack_t_pack(0,dcpl_mb(Gu))
      call Pack_t_pack(0,dcpl_mb(Gs))


*     **** tmp3(G) = Gu*Gs/G **** 
      fourpi = 16.0d0*datan(1.0d0)
      ss     = 1.0d0/fourpi
      call Pack_t_SMul(0,ss,dbl_mb(vc_indx),dbl_mb(tmp3(1)))
c      call Pack_tt_Sqrt(0,dbl_mb(tmp3(1)),dbl_mb(tmp3(1)))
c      call Pack_tt_Mul(0,dbl_mb(tmp3(1)),
c     >                   dcpl_mb(Gu),
c     >                   dbl_mb(tmp3(1)))
c      call Pack_tt_Mul(0,dbl_mb(tmp3(1)),
c     >                   dcpl_mb(Gs),
c     >                   dbl_mb(tmp3(1)))
      call Pack_tt_Sqrt1(0,dbl_mb(tmp3(1)))
      call Pack_tt_Mul2(0,dcpl_mb(Gu),dbl_mb(tmp3(1)))
      call Pack_tt_Mul2(0,dcpl_mb(Gs),dbl_mb(tmp3(1)))
      


*     **** calculate Aus ****  
      call dcopy(n2ft3d,0.0d0,0,Aus,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

        if (log_mb(semicore(1)+ia)) then

*          **** structure factor and local pseudopotential ****
           call strfac_pack(0,ii,dcpl_mb(exi(1)))

*          **** put phi and dphi/dG  at atom position ****
           call Pack_tc_MUL(0,dbl_mb(ncore(1)+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp1(1)))

           call Pack_tc_Mul(0,dbl_mb(ncore(1)+npack0+(ia-1)*npack0*5),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(tmp2(1)))
c           call Pack_tc_Mul(0,dbl_mb(tmp3(1)),
c     >                        dcpl_mb(tmp2(1)),
c     >                        dcpl_mb(tmp2(1)))
           call Pack_tc_Mul2(0,dbl_mb(tmp3(1)),
     >                        dcpl_mb(tmp2(1)))

*          **** Put put tmp1 and tmp2 into real space ****
           call Pack_c_unpack(0,dcpl_mb(tmp1(1)))
           call Pack_c_unpack(0,dcpl_mb(tmp2(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp1(1)))
           call D3dB_cr_fft3b(1,dcpl_mb(tmp2(1)))


*          ****  add tmp1*tmp2=phi*dphi*dphi/dG density ****
c           call D3dB_rr_Mul(1,dcpl_mb(tmp1(1)),
c     >                        dcpl_mb(tmp2(1)),
c     >                        dcpl_mb(tmp1(1)))
c           call D3dB_rr_Sum(1,Aus,
c     >                        dcpl_mb(tmp1(1)),
c     >                        Aus)
           call D3dB_rr_Mul2(1,dcpl_mb(tmp2(1)),dcpl_mb(tmp1(1)))
           call D3dB_rr_Sum2(1,dcpl_mb(tmp1(1)),Aus)
        end if

      end do
      call D3dB_r_Zero_Ends(1,Aus)
c      call D3dB_r_SMul(1,scal2*scal2,Aus,Aus)
      call D3dB_r_SMul1(1,scal2*scal2,Aus)


*     **** pop stack memory ****
      value =           MA_pop_stack(tmp3(2))
      value = value.and.MA_pop_stack(exi(2))
      value = value.and.MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      if (.not. value) 
     >   call errquit('semicore_gen_Aus:popping stack memory',2, MA_ERR)

      return
      end
