

*     ***********************************
*     *					*
*     *		electron_init		*
*     *					*
*     ***********************************
      subroutine electron_init()
      implicit none

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter


*     **** local variables ****
      logical value
      integer n2ft3d

*     **** external functions ****
      integer  psi_ispin,psi_ne,control_version
      external psi_ispin,psi_ne,control_version

      counter = 0

      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)

*     **** get nfft3d, and n2ft3d ****
      call Pack_npack(1,npack1)
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

*     **** allocate memory ****
      value = MA_alloc_get(mt_dcpl,npack1*(ne(1)+ne(2)),
     >                     'Hpsi_k',Hpsi_k(2),Hpsi_k(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d*(ne(1)+ne(2)),
     >                     'psi_r',psi_r(2),psi_r(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'vl2',vl(2),vl(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,n2ft3d,
     >                     'vl_lr',vl_lr(2),vl_lr(1))

      if (control_version().eq.3) then
        value = value.and.
     >        MA_alloc_get(mt_dcpl,npack0,
     >                     'vc',vc(2),vc(1))
      end if
      if (control_version().eq.4) then
        value = value.and.
     >        MA_alloc_get(mt_dcpl,nfft3d,
     >                     'vc',vc(2),vc(1))
      end if

      value = value.and.
     >        MA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xcp',xcp(2),xcp(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,2*n2ft3d,
     >                     'xce',xce(2),xce(1))
      if (.not. value) 
     >  call errquit('electron_init: out of heap memory',0)

      return
      end


*     ***********************************
*     *					*
*     *		electron_finalize       *
*     *					*
*     ***********************************
      subroutine electron_finalize()
      implicit none

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical value

*     **** free heap  memory ****
      value = MA_free_heap(Hpsi_k(2))
      value = value.and.
     >        MA_free_heap(psi_r(2))
      value = value.and.
     >        MA_free_heap(vl(2))
      value = value.and.
     >        MA_free_heap(vl_lr(2))
      value = value.and.
     >        MA_free_heap(vc(2))
      value = value.and.
     >        MA_free_heap(xcp(2))
      value = value.and.
     >        MA_free_heap(xce(2))
      if (.not. value) 
     >  call errquit('electron_init: error freeing heap memory',0)

      return
      end

*     ***********************************
*     *					*
*     *		electron_count		*
*     *					*
*     ***********************************
      integer function electron_count()
      implicit none

*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter

      electron_count = counter
      return
      end


*     ***********************************
*     *									*
*     *		electron_run				*
*     *									*
*     ***********************************
      subroutine electron_run(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

*     **** electron_counter common block ****
      integer counter
      common / electron_counter / counter

      counter = counter+1

      call electron_gen_psi_r(psi_k)
      call electron_gen_densities(psi_k,dn,dng,dnall)
      call electron_gen_scf_potentials(dn,dng,dnall)
      call electron_gen_Hpsi_k(psi_k)

      return
      end


*     ***********************************
*     *									*
*     *		electron_sd_update  		*
*     *									*
*     ***********************************
      subroutine electron_sd_update(psi1,psi2,dte)
      implicit none
      complex*16 psi1(*),psi2(*)
      real*8     dte

#include "mafdecls.fh"

      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

      call electron_sd_subupdate(npack1,(ne(1)+ne(2)),
     >                           psi1,psi2,dcpl_mb(Hpsi_k(1)),
     >                           dte)

      return
      end


      subroutine electron_sd_subupdate(npack1,nn,
     >                                 psi1,psi2,Hpsi,dte)
      implicit none
      integer    npack1,nn
      complex*16 psi1(npack1,nn)
      complex*16 psi2(npack1,nn)
      complex*16 Hpsi(npack1,nn)
      real*8     dte

      integer n

*     ************************************
*     **** do a steepest descent step ****
*     ************************************
      do n=1,nn
        call Pack_c_SMul(1,(-dte),Hpsi(1,n),psi2(1,n))
        call Pack_cc_Sum(1,psi2(1,n),psi1(1,n),psi2(1,n))
      end do

      return
      end

*     ***********************************
*     *					*
*     *		electron_energy		*
*     *					*
*     ***********************************
      real*8 function electron_energy(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0


*     **** local variables ****
      integer n2ft3d
      integer ii,ms,n1(2),n2(2),nx,ny,nz
      real*8  sum,eorbit,ehartr,exc,pxc,exc2,pxc2,dv

*     **** external functions *****
      integer  control_version
      real*8   lattice_omega,coulomb_e,electron_ehartree2
      external control_version
      external lattice_omega,coulomb_e,electron_ehartree2
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
     

      dv = lattice_omega()/dble(nx*ny*nz)

      n2ft3d = 2*nfft3d
      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)


*     *** get orbital energies ****
      eorbit = 0.0d0
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
c           call D3dB_cc_idot(1,
c     >                       psi_k(1+(ii-1)*nfft3d),
c     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*nfft3d),
c     >                       sum)
           call Pack_cc_idot(1,
     >                       psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       sum)

           eorbit = eorbit + sum

         end do
      end do
      call D3dB_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit+eorbit

           
*     **** get coulomb energy ****
      if (control_version().eq.3) ehartr = coulomb_e(dng)
      if (control_version().eq.4) ehartr = electron_ehartree2(dn)
        

*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dnall(1),
     >                dbl_mb(xce(1)),
     >                exc)
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     .                 pxc)
      if (ispin.eq.1) then
         exc= exc + exc 
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dnall(1+n2ft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
         call D3dB_rr_dot(1,dn(1+n2ft3d),
     >                    dbl_mb(xcp(1)+n2ft3d),
     >                    pxc2)
         exc= exc + exc2
         pxc= pxc + pxc2
      end if
      exc = exc*dv
      pxc = pxc*dv
          

*     **** total energy ****
      electron_energy = eorbit + exc - ehartr - pxc

      return
      end




*     ***********************************
*     *					*
*     *		electron_eorbit		*
*     *					*
*     ***********************************
      real*8 function electron_eorbit(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0


*     **** local variables ****
      integer ii,ms,n1(2),n2(2)
      real*8  sum,eorbit

      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)


*     *** get orbital energies ****
      eorbit = 0.0d0
      do ms=1,ispin
         do ii=n1(ms),n2(ms)
c           call D3dB_cc_idot(1,psi_k(1+(ii-1)*nfft3d),
c     >                      dcpl_mb(Hpsi_k(1)+(ii-1)*nfft3d),
c     >                       sum)
           call Pack_cc_idot(1,psi_k(1+(ii-1)*npack1),
     >                       dcpl_mb(Hpsi_k(1)+(ii-1)*npack1),
     >                       sum)

           eorbit = eorbit + sum

         end do
      end do
      call D3dB_SumAll(eorbit)
      if (ispin.eq.1) eorbit = eorbit+eorbit

           
      electron_eorbit = eorbit 

      return
      end
 

*     ***********************************
*     *					*
*     *		electron_ehartree	*
*     *					*
*     ***********************************
      real*8 function electron_ehartree(dng)
      implicit none
      complex*16 dng(*)


*     **** external functions ****
      real*8   coulomb_e
      external coulomb_e

      electron_ehartree = coulomb_e(dng)

      return
      end

*     ***********************************
*     *									*
*     *		electron_ehartree2			*
*     *									*
*     ***********************************
      real*8 function electron_ehartree2(dn)
      implicit none
      real*8     dn(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      real*8 ehartr, ehart1,ehart2,dv
      integer nx,ny,nz

*     ***** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      dv = lattice_omega()/dble(nx*ny*nx)

      call D3dB_rr_dot(1,dn(1),
     >                   dcpl_mb(vc(1)),
     >                   ehart1)
      call D3dB_rr_dot(1,dn(1+(ispin-1)*2*nfft3d),
     >                   dcpl_mb(vc(1)),
     >                   ehart2)
      ehartr = 0.5d0*(ehart1+ehart2)*dv

      electron_ehartree2 = ehartr

      return
      end

*     ***********************************
*     *					*
*     *		electron_exc		*
*     *					*
*     ***********************************
      real*8 function electron_exc(dnall)
      implicit none
      real*8 dnall(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      integer nx,ny,nz
      real*8  exc,exc2,dv

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dnall,
     >                 dbl_mb(xce(1)),
     >                 exc)
      if (ispin.eq.1) then
         exc= exc + exc 
      else
         call D3dB_rr_dot(1,dnall(1+2*nfft3d),
     >                    dbl_mb(xce(1)),
     >                    exc2)
         exc= exc + exc2
      end if
      exc = exc*dv
        
  
      electron_exc =  exc 

      return
      end


*     ***********************************
*     *					*
*     *		electron_pxc		*
*     *					*
*     ***********************************
      real*8 function electron_pxc(dn)
      implicit none
      real*8 dn(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      integer nx,ny,nz
      real*8  pxc,pxc2,dv

*     **** external functions *****
      real*8   lattice_omega
      external lattice_omega

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      
      dv = lattice_omega()/dble(nx*ny*nz)


*     **** get exchange-correlation energy ****
      call D3dB_rr_dot(1,dn(1),
     >                 dbl_mb(xcp(1)),
     >                 pxc)
      if (ispin.eq.1) then
         pxc= pxc + pxc 
      else
         call D3dB_rr_dot(1,dn(1+2*nfft3d),
     >                    dbl_mb(xcp(1)+2*nfft3d),
     >                    pxc2)
         pxc= pxc + pxc2
      end if
      pxc = pxc*dv
        
  
      electron_pxc =  pxc

      return
      end


 

*     ***********************************
*     *					*
*     *		electron_get_Hpsi_k	*
*     *					*
*     ***********************************
      subroutine electron_get_Hpsi_k(Hpsi_k_new)
      implicit none
      complex*16 Hpsi_k_new(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

      call dcopy(2*npack1*(ne(1)+ne(2)),
     >           dcpl_mb(Hpsi_k(1)),1,
     >           Hpsi_k_new,1)
      return
      end



*     ***************************
*     *				*
*     *	   electron_ispin	*
*     *				*
*     ***************************
      integer function electron_ispin()
      implicit none

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

      electron_ispin = ispin
      return
      end


*     ***************************
*     *				*
*     *	     electron_ne	*
*     *				*
*     ***************************
      integer function electron_ne(ms)
      implicit none
      integer ms

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

      electron_ne = ne(ms)
      return
      end


*     ***********************************
*     *					*
*     *	    electron_get_Tgradient 	*
*     *					*
*     ***********************************

      subroutine electron_get_Tgradient(psi_k,hml,THpsi_k)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)
      complex*16 THpsi_k(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     ***** local variables ****
      integer ms,n,n1(2),shift
     

      n1(1) = 1
      n1(2) = ne(1)+1


*     **** no preconditioning ***
      do ms=1,ispin
            shift = (ms-1)*ne(1)*ne(1)
            n     = ne(ms)
            call Grsm_gmg_Mul(npack1,n,
     >                        psi_k(1+(ms-1)*ne(1)*npack1),
     >                        hml(shift+1),
     >                        THpsi_k(1+(ms-1)*ne(1)*npack1))
      end do
      call Grsm_ggg_Sub(npack1,(ne(1)+ne(2)),
     >                  dcpl_mb(Hpsi_k(1)),
     >                  THpsi_k,
     >                  THpsi_k)
      call Grsm_gg_dScale(npack1,(ne(1)+ne(2)),
     >                    (-1.0d0),
     >                    THpsi_k,
     >                    THpsi_k)
      
      return
      end

*     ***************************
*     *				*
*     *	    electron_gen_hml 	*
*     *				*
*     ***************************

      subroutine electron_gen_hml(psi_k,hml)
      implicit none
      complex*16 psi_k(*)
      real*8     hml(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      integer ms,n,n1(2),shift
     

      n1(1) = 1
      n1(2) = ne(1) + 1


      do ms=1,ispin
         shift = (ms-1)*ne(1)*ne(1)
         n     = ne(ms)
         call Grsm_ggm_dot(npack1,n,
     >                     psi_k(1+(ms-1)*ne(1)*npack1),
     >                     dcpl_mb(Hpsi_k(1)+(ms-1)*ne(1)*npack1),
     >                     hml(shift+1))
      end do


      return
      end



**************************************************************************
**************************************************************************
*******    routines below this line are for internal use only    *********
**************************************************************************
**************************************************************************

*     ***********************************
*     *					*
*     *		electron_gen_Hpsi_k	*
*     *					*
*     ***********************************

      subroutine electron_gen_Hpsi_k(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical move
      real*8  fion(3,1)

*     **** external functions ****
      integer  control_version
      external control_version

      move = .false.
*     ******************
*     **** get Hpsi ****
*     ******************
      if (control_version().eq.3)
     >  call psi_H(ispin,ne,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dcpl_mb(vl(1)),
     >             dcpl_mb(vc(1)),
     >             dbl_mb(xcp(1)),
     >             dcpl_mb(Hpsi_k(1)),
     >             move,
     >             fion)

      if (control_version().eq.4)
     >  call psi_Hv4(ispin,ne,psi_k,
     >             dbl_mb(psi_r(1)),
     >             dcpl_mb(vl(1)),
     >             dbl_mb(vl_lr(1)),
     >             dcpl_mb(vc(1)),
     >             dbl_mb(xcp(1)),
     >             dcpl_mb(Hpsi_k(1)),
     >             move,
     >             fion)

      call Grsm_gg_dscale(npack1,(ne(1)+ne(2)),(-1.0d0),
     >                    dcpl_mb(Hpsi_k(1)),
     >                    dcpl_mb(Hpsi_k(1)))

      return
      end

*     ***************************
*     *				*
*     *	    electron_gen_psi_r	*
*     *				*
*     ***************************

      subroutine electron_gen_psi_r(psi_k)
      implicit none
      complex*16 psi_k(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      integer n,nemax,n2ft3d

      n2ft3d = 2*nfft3d
      nemax = ne(1) + ne(2)

c     call Grsm_gg_Copy(npack1,nemax,psi_k,dbl_mb(psi_r(1)))
      do n=1,nemax
         call Pack_c_Copy(1,psi_k(1+(n-1)*npack1),
     >                    dbl_mb(psi_r(1)+(n-1)*n2ft3d))
      end do

      call Grsm_gh_fftb(nfft3d,nemax,dbl_mb(psi_r(1)))
      call Grsm_h_Zero_Ends(nfft3d,nemax,dbl_mb(psi_r(1)))

      return
      end


*     ***************************
*     *				*
*     *	 electron_gen_densities	*
*     *				*
*     ***************************

      subroutine electron_gen_densities(psi_k,dn,dng,dnall)
      implicit none
      complex*16 psi_k(*)
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0
     
*     **** local variables ****
      integer i,ms,nx,ny,nz,n2ft3d
      integer n,n1(2),n2(2)
      real*8  scal1,scal2
      integer tmp1(2)
      logical value

*     ***** external functions *****
      logical  psp_semicore
      real*8   lattice_omega
      external psp_semicore
      external lattice_omega

      n1(1) = 1
      n1(2) = ne(1) + 1
      n2(1) = ne(1)
      n2(2) = ne(1) + ne(2)

      n2ft3d = 2*nfft3d
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()


*     *********************
*     **** generate dn ****
*     *********************
      call dcopy(2*n2ft3d,0.0d0,0,dn,1)
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            do i=1,n2ft3d
               dn(i+(ms-1)*n2ft3d) 
     >            = dn(i+(ms-1)*n2ft3d) 
     >            + scal2*(dbl_mb(psi_r(1)+i-1+(n-1)*n2ft3d)**2)
            end do
         end do
         call D3dB_r_Zero_Ends(1,dn(1+(ms-1)*n2ft3d))
      end do
        
*     **********************
*     **** generate dng ****
*     **********************
      value = MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_densities: out of stack memory',0)

      call D3dB_rr_Sum(1,dn,dn(1+(ispin-1)*n2ft3d),dbl_mb(tmp1(1)))
      call D3dB_r_SMul(1,scal1,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
      call D3dB_rc_fft3f(1,dbl_mb(tmp1(1)))
      call Pack_c_pack(0,dbl_mb(tmp1(1)))
      call Pack_c_Copy(0,dbl_mb(tmp1(1)),dng)


*       ********************************************************
*       **** generate dnall - used for semicore corrections ****
*       ********************************************************
        if (psp_semicore(0)) then
           call semicore_density(dbl_mb(tmp1(1)))
           call D3dB_r_SMul(1,0.5d0,dbl_mb(tmp1(1)),dbl_mb(tmp1(1)))
        else
           call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
        end if
        do ms=1,ispin
          call D3dB_rr_Sum(1,dn(1+(ms-1)*n2ft3d),
     >                     dbl_mb(tmp1(1)),
     >                     dnall(1+(ms-1)*n2ft3d))
        end do



      value = MA_pop_stack(tmp1(2)) 
      if (.not. value) call errquit(
     >     'electron_gen_densities: error popping stack memory',0)


      return
      end



*     ***********************************
*     *			 		*
*     *   electron_gen_scf_potentials	*
*     *					*
*     ***********************************

      subroutine electron_gen_scf_potentials(dn,dng,dnall)

      implicit none
      real*8     dn(*)
      complex*16 dng(*)
      real*8     dnall(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     ***** local variables ****
      logical value
      integer n2ft3d
      integer tmp1(2)

*     **** external functions ****
      integer  control_gga,control_version
      external control_gga,control_version
     
      n2ft3d = 2*nfft3d

      if (control_version().eq.3) then
         call coulomb_v(dng,dcpl_mb(vc(1)))
      end if

      if (control_version().eq.4)  then
         value = MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',0)
      
         call D3dB_rr_Sum(1,dn(1),dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(tmp1(1)))
         call coulomb2_v(dbl_mb(tmp1(1)),dcpl_mb(vc(1)))

         value = MA_pop_stack(tmp1(2)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',0)
      end if

      if (control_gga().eq.0) then
         value = MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
         if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: out of stack memory',0)
        call vxc(n2ft3d,ispin,dnall,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)),
     >                    dbl_mb(tmp1(1)))
        value = MA_pop_stack(tmp1(2)) 
        if (.not. value) call errquit(
     >     'electron_gen_scf_potentials: error popping stack memory',0)

      else
        call v_bwexc(n2ft3d,ispin,dnall,
     >                    dbl_mb(xcp(1)),
     >                    dbl_mb(xce(1)))
      end if

      return
      end




*     ***********************************
*     *			 		*
*     *   electron_gen_vl_potential 	*
*     *					*
*     ***********************************

      subroutine electron_gen_vl_potential()
      implicit none

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical move,value
      integer n2ft3d
      integer tmp1(2)
      integer tmp2(2)
      integer r_grid(2)

*     **** external functions *****
      logical  pspw_qmmm_found
      integer  control_version
      external pspw_qmmm_found
      external control_version

      value = MA_push_get(mt_dcpl,(nfft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl,(3),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call 
     >   errquit('electron_gen_vl_potential: out of stack memory',0)

      move = .false.
      call v_local(dcpl_mb(vl(1)),
     >               move,
     >               dcpl_mb(tmp1(1)),
     >               dbl_mb(tmp2(1)))

      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >  'electron_gen_vl_potential: error popping stack memory',0)


*      *** long-range psp for charge systems ***
      if (control_version().eq.4) then
         value = MA_push_get(mt_dbl,(6*nfft3d),'r_grid',
     >                       r_grid(2),r_grid(1)) 

         call lattice_r_grid(dbl_mb(r_grid(1)))
         call v_lr_local(dbl_mb(r_grid(1)),
     >                   dbl_mb(vl_lr(1)))

*        **** generate QM/MM potential ****
         if (pspw_qmmm_found()) then
           n2ft3d = 2*nfft3d
           call pspw_qmmm_Generate_V(n2ft3d,dbl_mb(r_grid(1)),
     >                                     dbl_mb(vl_lr(1)))
         end if

         value = MA_pop_stack(r_grid(2))
         if (.not. value) call errquit(
     >   'electron_gen_vl_potential: error popping stack memory',0)

      end if
 


      return
      end



*     ***********************************
*     *			 		*
*     *   electron_psi_vl_ave	 	*
*     *					*
*     ***********************************

      real*8 function electron_psi_vl_ave(psi1,dn)
      implicit none
      complex*16 psi1(*)
      real*8     dn(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2)
      integer nx,ny,nz,n2ft3d,np
      real*8 elocal,sum,scal1,scal2,dv
      integer tmp1(2),tmp2(2)

*     **** external functions ***
      integer  control_version
      real*8   lattice_omega
      external control_version
      external lattice_omega


      call Parallel_np(np)

      n2ft3d = 2*nfft3d
      value = MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      scal1 = 1.0d0/dble(nx*ny*nz)
      scal2 = 1.0d0/lattice_omega()
      dv    = scal1/scal2


*     **** average Kohn-Sham v_local energy ****
       call Pack_c_Copy(0,dcpl_mb(vl(1)),dbl_mb(tmp1(1)))
      call Pack_c_unpack(0,dbl_mb(tmp1(1)))
      call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
      elocal = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call D3dB_rr_Mul(1,
     >                       dbl_mb(tmp1(1)),
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)))

c           call D3dB_rc_fft3f(1,dbl_mb(tmp2(1)))
c           call Pack_c_pack(1,dbl_mb(tmp2(1)))
c           call Pack_cc_dot(1,psi1(1+(n-1)*npack1),
c    >                         dbl_mb(tmp2(1)),
c    >                         sum)
            
            call D3dB_rr_idot(1,
     >                       dbl_mb(psi_r(1)+(n-1)*n2ft3d),
     >                       dbl_mb(tmp2(1)),
     >                       sum)

            elocal = elocal + sum*scal1*scal2
         end do
      end do
      if (np.gt.1) call D3dB_SumAll(elocal)
      if (ispin.eq.1) elocal = 2.0d0*elocal

*     *** add in long range part ****
      if (control_version().eq.4) then
       call D3dB_rr_dot(1,dn(1),dbl_mb(vl_lr(1)),sum)
       elocal = elocal + sum*dv
       call D3dB_rr_dot(1,dn(1+(ispin-1)*n2ft3d),
     >                    dbl_mb(vl_lr(1)),sum)
       elocal = elocal + sum*dv
      end if


 
      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0)

      electron_psi_vl_ave = elocal
      return
      end



*     ***********************************
*     *			 						*
*     *   electron_psi_vnl_ave	 		*
*     *									*
*     ***********************************

      real*8 function electron_psi_vnl_ave(psi1)
      implicit none
      complex*16 psi1(*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

*     **** local variables ****
      logical value
      integer n,ms,n1(2),n2(2),np
      integer nee(2)
      integer n2ft3d
      real*8 enlocal,sum
      integer tmp1(2),tmp2(2)


      call Parallel_np(np)

      n2ft3d = 2*nfft3d
      value = MA_push_get(mt_dbl,(n2ft3d),'tmp1',tmp1(2),tmp1(1)) 
      value = value.and.
     >        MA_push_get(mt_dbl,(n2ft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit(
     >            'electron_psi_vl_ave: out of stack memory',0)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)



*     **** average Kohn-Sham v_nonlocal energy ****
      nee(1) = 1
      nee(2) = 0
      enlocal = 0.0d0
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            call dcopy(n2ft3d,0.0d0,0,dbl_mb(tmp1(1)),1)
            call v_nonlocal(ispin,nee,psi1(1+(n-1)*npack1),
     >                      dbl_mb(tmp1(1)),
     >                      .false.,dbl_mb(tmp2(1)))
            call Pack_cc_idot(1,psi1(1+(n-1)*npack1),
     >                         dbl_mb(tmp1(1)),
     >                         sum)
            enlocal = enlocal - sum
         end do
      end do
      if (np.gt.1) call D3dB_SumAll(enlocal)
      if (ispin.eq.1) enlocal = 2.0d0*enlocal

 
      value = MA_pop_stack(tmp2(2))
      value = value.and.
     >        MA_pop_stack(tmp1(2))
      if (.not. value) call errquit(
     >           'electron_psi_vl_ave: error popping stack memory',0)

      electron_psi_vnl_ave = enlocal
      return
      end




*     ***********************************
*     *			 		*
*     *   electron_semicoreforce 	*
*     *					*
*     ***********************************

      subroutine electron_semicoreforce(fion)
      implicit none
      real*8 fion(3,*)

#include "mafdecls.fh"

*     ***** electron common block ****
      integer ispin,ne(2),nfft3d,npack1,npack0
      integer Hpsi_k(2)
      integer psi_r(2)
      integer vl(2)
      integer vl_lr(2)
      integer vc(2)
      integer xcp(2)
      integer xce(2)
      common / electron_block / Hpsi_k,psi_r,
     >                          vl,vl_lr,vc,xcp,xce,
     >                          ne,ispin,nfft3d,npack1,npack0

     
      call semicore_xc_F(ispin,dbl_mb(xcp(1)),fion)

      return
      end




