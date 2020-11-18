*
* $Id$
*

*     *************************************
*     *                                   *
*     *             psi_H                 *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,neq           - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vl                  - local psp in k-space
*           v_field,field_exist - external fields (e.g. QM/MM)
*           vc                  - coulomb potential in k-space
*           xcp                 - xc potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*
      subroutine psi_H(ispin,neq,
     >                 psi_k,psi_r,
     >                 vl,
     >                 v_field,field_exist,
     >                 vc,xcp,Hpsi_k,
     >                 move,fion,
     >                 fractional,occ)
      implicit none
      integer    ispin,neq(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vl(*)
      real*8     v_field(*)
      logical    field_exist
      complex*16 vc(*)
      real*8     xcp(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)
      logical fractional
      real*8 occ(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dosic,dohfx,dodftu,meta,tamd,meta_gga
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1,indx1,indx2
      integer vall(2),vpsi(2),Hpsi_r(2)
      logical value,done


*     **** external functions ****
      logical  pspw_SIC,pspw_SIC_relaxed,D3dB_rc_pfft3_queue_filled
      logical  pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm,meta_found
      logical  psp_pawexist,nwpw_meta_gga_on,tamd_found
      real*8   lattice_omega
      external pspw_SIC,pspw_SIC_relaxed,D3dB_rc_pfft3_queue_filled
      external pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm,meta_found
      external psp_pawexist,nwpw_meta_gga_on,tamd_found
      external lattice_omega


      dosic = (pspw_SIC().and.pspw_SIC_relaxed()) 
      dohfx = (pspw_HFX().and.pspw_HFX_relaxed()) 
      dodftu = psp_U_psputerm()
      meta = meta_found()
      tamd = tamd_found()
      meta_gga = nwpw_meta_gga_on()
      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
     
      if ((dosic).or.(dohfx)) then
         value = value.and.
     >           BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,
     >                        'Hpsi_r',Hpsi_r(2),Hpsi_r(1))
         !call dcopy((ne(1)+ne(2))*n2ft3d,0.0d0,0,dbl_mb(Hpsi_r(1)),1)
         call Parallel_shared_vector_zero(.false.,
     >                                    (neq(1)+neq(2))*n2ft3d,
     >                                    dbl_mb(Hpsi_r(1)))
      end if
      if (.not. value) call errquit('psi_H:out of stack memory',0,
     &       MA_ERR)
      !call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vall(1)),1)
      !call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vpsi(1)),1)

      !call D3dB_c_Zero(1,dcpl_mb(vall(1)))
      !call D3dB_c_Zero(1,dcpl_mb(vpsi(1)))
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1) + 1
      n2(2) = neq(1) + neq(2)

*     *********************************************************
*     **** add kinetic energy - Hpsi_k = Ke*psi_k ****
*     *********************************************************
      call ke(ispin,neq,psi_k,Hpsi_k)

*     *********************************************************
*     **** add non-local psp - Hpsi_k = Hpsi_k + Vnl*psi_k ****
*     *********************************************************
      if (psp_pawexist()) call psp_set_vc(vc)  !*** this should be located where vc is generated ***
      call v_nonlocal(ispin,neq,psi_k,Hpsi_k,move,fion,fractional,occ)


*     *********************************************************
*     **** add DFT+U correction -                             *
*     *********************************************************
      if (dodftu) then
         call psp_U_v_nonlocal(ispin,neq,psi_k,Hpsi_k,move,fion,
     >                         fractional,occ,.true.)
      end if

*     *********************************************************
*     **** add Metadynamics potential                         *
*     *********************************************************
      if (meta) then
         call meta_force(ispin,neq,psi_k,Hpsi_k,move,fion)
      end if

*     *********************************************************
*     **** add TAMD potential                                 *
*     *********************************************************
      if (tamd) then
         call tamd_force(ispin,neq,psi_k,Hpsi_k,move,fion)
      end if

*     ****************************************************************************
*     **** add extra Meta GGA terms Hpsi_k -= 0.5 * grad dot (dfdtau*gradpsi) ****
*     ****************************************************************************
      if (meta_gga) then
         call nwpw_meta_gga_Hpsik(ispin,neq,psi_k,Hpsi_k)
      end if


*     *******************************************************
*     **** SIC corrections- Hpsi_r = Vsic*psi_r ****
*     *******************************************************
      if (dosic) then
         call pspw_potential_SIC(ispin,psi_r,dbl_mb(Hpsi_r(1)))
      end if

*     *******************************************************
*     **** HF exchange corrections- Hpsi_r = Vhfx*psi_r ****
*     *******************************************************
      if (dohfx) then
         call pspw_potential_HFX(ispin,psi_r,dbl_mb(Hpsi_r(1)))
      end if


*     **************************************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  (+ vcmp if paw) ****
*     **************************************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      call Pack_cc_Sum2(0,vc,dcpl_mb(vall(1)))
      if (psp_pawexist()) call psp_add_vcmp(dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     ******************************
*     **** add v_field to vall  ****
*     ******************************
c      if (field_exist)
c     >  call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_field,dcpl_mb(vall(1)))
      if (field_exist)
     >  call D3dB_rr_Sum2(1,v_field,dcpl_mb(vall(1)))


*     **************************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r +Vsic*psi_r] **** 
*     ************************************************************** 


      indx1 = 1
      indx2 = 1
      done  = .false.
      do while (.not.done)
         if (indx1.le.n2(ispin)) then
            ms = 1
            if (indx1.gt.n2(1)) ms = 2

*           **** vpsi_r = (vall_r+xcp)*psi_r + Hpsi_r****
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))

c            call D3dB_rr_Mul(1,dcpl_mb(vpsi(1)),
c     >                       psi_r(1+(indx1-1)*n2ft3d),
c     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul2(1,psi_r(1+(indx1-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))

c            if ((dosic).or.(dohfx))
c     >        call D3dB_rr_Sum(1,dcpl_mb(vpsi(1)),
c     >                         dbl_mb(Hpsi_r(1)+(indx1-1)*n2ft3d),
c     >                         dcpl_mb(vpsi(1)))
            if ((dosic).or.(dohfx))
     >        call D3dB_rr_Sum2(1,dbl_mb(Hpsi_r(1)+(indx1-1)*n2ft3d),
     >                            dcpl_mb(vpsi(1)))

*           **** put vpsi_r in k-space ****
            call D3dB_rc_pfft3f_queuein(1,dcpl_mb(vpsi(1)))
            indx1 = indx1 + 1

         end if

         if ((D3dB_rc_pfft3_queue_filled()).or.
     >       (indx1.gt.n2(ispin))) then
            call D3dB_rc_pfft3f_queueout(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(indx2-1)*npack1))
            indx2 = indx2 + 1
         end if

         done = ((indx1.gt.n2(ispin)).and.(indx2.gt.n2(ispin)))
      end do
 

      call nwpw_timing_start(12)

      value = .true.
      if ((dosic).or.(dohfx)) 
     >  value = value.and.BA_pop_stack(Hpsi_r(2))
      value = value.and.BA_pop_stack(vpsi(2))
      value = value.and.BA_pop_stack(vall(2))
      if (.not. value) call errquit('psi_H:popping stack memory',1,0)

      call nwpw_timing_end(12)
      return
      end




*     *************************************
*     *                                   *
*     *             psi_H_vall            *
*     *                                   *
*     *************************************
*
*  This subroutine is used by block cg and bfgs minimizers.
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,neq           - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vall                - local psp in k-space
*   Exit - Hpsi_k - gradient in k-space
*
      subroutine psi_H_vall(ispin,neq,
     >                 psi_k,psi_r,
     >                 vall,
     >                 Hpsi_k)
      implicit none
      integer    ispin,neq(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      real*8     vall(*)
      complex*16 Hpsi_k(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dosic,dohfx,dodftu,meta,tamd,meta_gga
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      !real*8  scal1,scal2
      real*8  scal1,fion(3,1)

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1,indx1,indx2
      integer vpsi(2),Hpsi_r(2)
      logical value,done

*     **** external functions ****
      logical  pspw_SIC,pspw_SIC_relaxed,D3dB_rc_pfft3_queue_filled
      logical  pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm,meta_found
      logical  psp_pawexist,nwpw_meta_gga_on,tamd_found
      real*8   lattice_omega
      external pspw_SIC,pspw_SIC_relaxed,D3dB_rc_pfft3_queue_filled
      external pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm,meta_found
      external psp_pawexist,nwpw_meta_gga_on,tamd_found
      external lattice_omega

      dosic = (pspw_SIC().and.pspw_SIC_relaxed()) 
      dohfx = (pspw_HFX().and.pspw_HFX_relaxed()) 
      dodftu = psp_U_psputerm()
      meta = meta_found()
      tamd = tamd_found()
      meta_gga = nwpw_meta_gga_on()
      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d

      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if ((dosic).or.(dohfx)) then
         value = value.and.
     >           BA_push_get(mt_dbl,(neq(1)+neq(2))*n2ft3d,
     >                        'Hpsi_r',Hpsi_r(2),Hpsi_r(1))
         call Parallel_shared_vector_zero(.false.,
     >                                    (neq(1)+neq(2))*n2ft3d,
     >                                    dbl_mb(Hpsi_r(1)))
      end if
      if (.not.value) 
     >   call errquit('psi_H_vall:out of stack memory',0,MA_ERR)
      call nwpw_timing_end(12)
      

      !scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = neq(1)
      n1(2) = neq(1) + 1
      n2(2) = neq(1) + neq(2)

*     *********************************************************
*     **** add kinetic energy - Hpsi_k = Ke*psi_k ****
*     *********************************************************
      call ke(ispin,neq,psi_k,Hpsi_k)


*     *********************************************************
*     **** add non-local psp - Hpsi_k = Hpsi_k + Vnl*psi_k ****
*     *********************************************************
      !if (psp_pawexist()) call psp_set_vc(vc) !*** this should be where vc is generated ***
      call v_nonlocal(ispin,neq,psi_k,Hpsi_k,.false.,fion,.false.,fion)


*     *********************************************************
*     **** add DFT+U correction -                             *
*     *********************************************************
      if (dodftu) then
         call psp_U_v_nonlocal(ispin,neq,psi_k,Hpsi_k,.false.,fion,
     >                         .false.,fion,.false.)
      end if

*     *********************************************************
*     **** add Metadynamics potential                         *
*     *********************************************************
      if (meta) then
         call meta_force(ispin,neq,psi_k,Hpsi_k,.false.,fion)
      end if

*     *********************************************************
*     **** add TAMD potential                                 *
*     *********************************************************
      if (tamd) then
         call tamd_force(ispin,neq,psi_k,Hpsi_k,.false.,fion)
      end if

*     ****************************************************************************
*     **** add extra Meta GGA terms Hpsi_k -= 0.5 * grad dot (dfdtau*gradpsi) ****
*     ****************************************************************************
      if (meta_gga) then
         call nwpw_meta_gga_Hpsik(ispin,neq,psi_k,Hpsi_k)
      end if


*     *******************************************************
*     **** SIC corrections- Hpsi_r = Vsic*psi_r ****
*     *******************************************************
      if (dosic) then
         call pspw_potential_SIC(ispin,psi_r,dbl_mb(Hpsi_r(1)))
      end if

*     *******************************************************
*     **** HF exchange corrections- Hpsi_r = Vhfx*psi_r ****
*     *******************************************************
      if (dohfx) then
         call pspw_potential_HFX(ispin,psi_r,dbl_mb(Hpsi_r(1)))
      end if


*     **************************************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  (+ vcmp if paw) ****
*     **************************************************************************
      !call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      !call Pack_cc_Sum2(0,vc,dcpl_mb(vall(1)))
      !if (psp_pawexist()) call psp_add_vcmp(dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      !call Pack_c_unpack(0,dcpl_mb(vall(1)))
      !call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     ******************************
*     **** add v_field to vall  ****
*     ******************************
c      if (field_exist)
c     >  call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_field,dcpl_mb(vall(1)))
c      if (field_exist)
c     >  call D3dB_rr_Sum2(1,v_field,dcpl_mb(vall(1)))


*     *******************************************************
*     **** Hpsi_k = Hpsi_k + FFT[(vall)*psi_r + Hpsi_r]  **** 
*     ******************************************************* 
      indx1 = 1
      indx2 = 1
      done  = .false.
      do while (.not.done)
         if (indx1.le.n2(ispin)) then
            ms = 1
            if (indx1.gt.n2(1)) ms = 2

*           **** vpsi_r = (vall_r)*psi_r + Hpsi_r****
            call D3dB_rr_Mul(1,vall(1+(ms-1)*n2ft3d),
     >                       psi_r(1+(indx1-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))

            if ((dosic).or.(dohfx))
     >        call D3dB_rr_Sum2(1,dbl_mb(Hpsi_r(1)+(indx1-1)*n2ft3d),
     >                            dcpl_mb(vpsi(1)))

*           **** put vpsi_r in k-space ****
            call D3dB_rc_pfft3f_queuein(1,dcpl_mb(vpsi(1)))
            indx1 = indx1 + 1

         end if

         if ((D3dB_rc_pfft3_queue_filled()).or.
     >       (indx1.gt.n2(ispin))) then
            call D3dB_rc_pfft3f_queueout(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(indx2-1)*npack1))
            indx2 = indx2 + 1
         end if

         done = ((indx1.gt.n2(ispin)).and.(indx2.gt.n2(ispin)))
      end do
 

      call nwpw_timing_start(12)
      value = .true.
      if ((dosic).or.(dohfx)) value=value.and.BA_pop_stack(Hpsi_r(2))
      value = value.and.BA_pop_stack(vpsi(2))
      if (.not.value) 
     >  call errquit('psi_H_vall:popping stack memory',1,MA_ERR)
      call nwpw_timing_end(12)

      return
      end







*     *************************************
*     *                                   *
*     *             psi_Horb              *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number 
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*   Exit - Hpsi_k - gradient in k-space
*
      subroutine psi_Horb(hfxon,ispin,ms,
     >                 psi_r,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      logical    hfxon
      integer    ispin,ms
      real*8     psi_r(*)
      real*8     vall(*)      
      complex*16 orb_k(*)
      real*8     orb_r(*)
      complex*16 Horb_k(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dohfx,meta_gga
      integer nne(2)
      integer nx,ny,nz
      real*8  scal1,fion(3,1)

*     **** tempory space ****
      integer nfft3d,n2ft3d,npack1
      integer vpsi(2),Horb_r(2)
      logical value


*     **** external functions ****
      logical  pspw_SIC,pspw_SIC_relaxed,nwpw_meta_gga_on
      logical  pspw_HFX,pspw_HFX_relaxed
      real*8   lattice_omega
      external pspw_SIC,pspw_SIC_relaxed,nwpw_meta_gga_on
      external pspw_HFX,pspw_HFX_relaxed
      external lattice_omega


      dohfx = pspw_HFX().and.pspw_HFX_relaxed().and.hfxon
      meta_gga = nwpw_meta_gga_on()
      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (dohfx) then
         value = value.and.
     >           BA_push_get(mt_dbl,n2ft3d,'Horb_r',Horb_r(2),Horb_r(1))
         call dcopy(n2ft3d,0.0d0,0,dbl_mb(Horb_r(1)),1)
      end if
      if (.not. value) 
     >  call errquit('psi_Horb:out of stack memory',0, MA_ERR)
      !call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vpsi(1)),1)
      call D3dB_c_Zero(1,dcpl_mb(vpsi(1)))
      call nwpw_timing_end(12)
      

      if (ms.eq.1) then
        nne(1) = 1
        nne(2) = 0
      else
        nne(1) = 0
        nne(2) = 1
      end if
*
      scal1 = 1.0d0/dble(nx*ny*nz)


*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call ke(ispin,nne,orb_k,Horb_k)

*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call v_nonlocal(ispin,nne,orb_k,Horb_k,
     >                .false.,fion,.false.,fion)


*     **********************************************************************
*     **** add extra Meta GGA terms Hpsi_k += grad dot (dfdtau*gradpsi) ****
*     **********************************************************************
      if (meta_gga) then
         call nwpw_meta_gga_Hpsik(ispin,nne,orb_k,Horb_k)
      end if

*     ***************************
*     **** SIC corrections?? ****
*     ***************************
c     if (pspw_SIC().and.pspw_SIC_relaxed()) then
c        call pspw_potential_SIC(ispin,psi_r,Hpsi_k)
c     end if

*     **********************
*     **** HFX energy   ****
*     **********************
      if (dohfx) then
         call pspw_potential_HFX_orb(ms,
     >                               psi_r,
     >                               orb_r,
     >                               dbl_mb(Horb_r(1)))
      end if


*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] **** 
*     **************************************************        

*     **** vpsi_r = (vall_r)*orb_r ****
      call D3dB_rr_Mul(1,vall(1+(ms-1)*n2ft3d),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))

*     **** vpsi_r = vpsi_r + Hpsi_r ****
c      if (dohfx)
c     >   call D3dB_rr_Sum(1,dcpl_mb(vpsi(1)),
c     >                      dbl_mb(Horb_r(1)),
c     >                      dcpl_mb(vpsi(1)))
      if (dohfx)
     >   call D3dB_rr_Sum2(1,dbl_mb(Horb_r(1)),
     >                      dcpl_mb(vpsi(1)))

*     **** put vpsi_r in k-space ****
      call D3dB_rc_pfft3f(1,1,dcpl_mb(vpsi(1)))
      call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*     **** add to Horb_k ****
      call Pack_cc_daxpy(1,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)

      call nwpw_timing_start(12)
      value = .true.
      if (dohfx) 
     >  value = value.and.BA_pop_stack(Horb_r(2))
      value = value.and.BA_pop_stack(vpsi(2))
      if (.not.value) 
     >  call errquit('psi_Horb: popping stack',1, MA_ERR)
      call nwpw_timing_end(12)

      return
      end


*     *************************************
*     *                                   *
*     *         psi_Horb_replicated       *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number 
*           psi_k,psi_r         - distributed orbitals in k-space and r-space
*           orb_k,orb_h         - replicated orbital in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*           notvirtual
*   Exit - Horb_k - replicated gradient in k-space
*
      subroutine psi_Horb_replicated(hfxon,ispin,ms,
     >                 psi_r,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      logical    hfxon
      integer    ispin,ms
      real*8     psi_r(*)
      real*8     vall(*)      
      complex*16 orb_k(*)
      real*8     orb_r(*)
      complex*16 Horb_k(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dohfx,dodftu
      integer nne(2)
      integer nx,ny,nz
      real*8  scal1,fion(3,1)

*     **** tempory space ****
      integer nfft3d,n2ft3d,npack1
      integer vpsi(2),Horb_r(2)
      logical value,meta_gga


*     **** external functions ****
      logical  pspw_SIC,pspw_SIC_relaxed,nwpw_meta_gga_on
      logical  pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm
      real*8   lattice_omega
      external pspw_SIC,pspw_SIC_relaxed,nwpw_meta_gga_on
      external pspw_HFX,pspw_HFX_relaxed,psp_U_psputerm
      external lattice_omega


      dohfx = pspw_HFX().and.pspw_HFX_relaxed().and.hfxon
      dodftu = psp_U_psputerm().and.hfxon
      meta_gga = nwpw_meta_gga_on()
      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (dohfx) then
         value = value.and.
     >           BA_push_get(mt_dbl,n2ft3d,'Horb_r',Horb_r(2),Horb_r(1))
         call dcopy(n2ft3d,0.0d0,0,dbl_mb(Horb_r(1)),1)
      end if
      if (.not. value) 
     >  call errquit('psi_Horb:out of stack memory',0, MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vpsi(1)),1)
      call nwpw_timing_end(12)
      
      if (ms.eq.1) then
        nne(1) = 1
        nne(2) = 0
      else
        nne(1) = 0
        nne(2) = 1
      end if
*
      scal1 = 1.0d0/dble(nx*ny*nz)


*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call ke(ispin,nne,orb_k,Horb_k)

*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call v_nonlocal(ispin,nne,orb_k,Horb_k,
     >                .false.,fion,.false.,fion)

*     *********************************************************
*     **** add DFT+U term                                  ****
*     *********************************************************
      if (dodftu) then
         call psp_U_v_nonlocal(ispin,nne,orb_k,Horb_k,
     >                .false.,fion,.false.,fion,.false.)
      end if

*     **********************************************************************
*     **** add extra Meta GGA terms Hpsi_k += grad dot (dfdtau*gradpsi) ****
*     **********************************************************************
      if (meta_gga) then
         call nwpw_meta_gga_Hpsik(ispin,nne,orb_k,Horb_k)
      end if


*     ***************************
*     **** SIC corrections?? ****
*     ***************************
c     if (pspw_SIC().and.pspw_SIC_relaxed()) then
c        call pspw_potential_SIC(ispin,psi_r,Hpsi_k)
c     end if

*     **********************
*     **** HFX energy   ****
*     **********************
      if (dohfx) then
         call pspw_potential_HFX_orb_replicated(ms,
     >                               psi_r,
     >                               orb_r,
     >                               dbl_mb(Horb_r(1)))
      end if


*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] **** 
*     **************************************************        

*     **** vpsi_r = (vall_r)*orb_r ****
      call D3dB_rr_Mul(1,vall(1+(ms-1)*n2ft3d),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))

*     **** vpsi_r = vpsi_r + Hpsi_r ****
      if (dohfx)
     >   call D3dB_rr_Sum2(1,dbl_mb(Horb_r(1)),
     >                      dcpl_mb(vpsi(1)))

*     **** put vpsi_r in k-space ****
      call D3dB_rc_pfft3f(1,1,dcpl_mb(vpsi(1)))
      call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*     **** add to Horb_k ****
      call Pack_cc_daxpy(1,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)

      call nwpw_timing_start(12)
      value = .true.
      if (dohfx) 
     >  value = value.and.BA_pop_stack(Horb_r(2))
      value = value.and.BA_pop_stack(vpsi(2))
      if (.not.value) 
     >  call errquit('psi_Horb: popping stack',1, MA_ERR)
      call nwpw_timing_end(12)

      return
      end



*     *************************************
*     *                                   *
*     *             psi_H1                *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k + FFT[(vall)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,ne            - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vl                  - local psp in k-space
*           v_field,field_exist - external fields (e.g. QM/MM)
*   Exit - Hpsi_k - gradient in k-space
*
      subroutine psi_H1(ispin,ne,
     >                 psi_k,psi_r,
     >                 vl,
     >                 v_field,field_exist,
     >                 Hpsi_k)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vl(*)
      real*8     v_field(*)
      logical    field_exist
      complex*16 Hpsi_k(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1,indx1,indx2
      integer vall(2),vpsi(2)
      logical value,done,meta_gga
      real*8 fion(3),occ(3)


*     **** external functions ****
      logical  D3dB_rc_pfft3_queue_filled,nwpw_meta_gga_on
      real*8   lattice_omega
      external D3dB_rc_pfft3_queue_filled,nwpw_meta_gga_on
      external lattice_omega


      meta_gga = nwpw_meta_gga_on()
      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('psi_H:out of stack memory',0,
     &       MA_ERR)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vall(1)),1)
      call dcopy(2*nfft3d,0.0d0,0,dcpl_mb(vpsi(1)),1)
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

*     *********************************************************
*     **** add kinetic energy - Hpsi_k = Ke*psi_k ****
*     *********************************************************
      call ke(ispin,ne,psi_k,Hpsi_k)

*     *********************************************************
*     **** add non-local psp - Hpsi_k = Hpsi_k + Vnl*psi_k ****
*     *********************************************************
      call v_nonlocal(ispin,ne,psi_k,Hpsi_k,.false.,fion,.false.,occ)



*     **********************************************************************
*     **** add extra Meta GGA terms Hpsi_k += grad dot (dfdtau*gradpsi) ****
*     **********************************************************************
      if (meta_gga) then
         call nwpw_meta_gga_Hpsik(ispin,ne,psi_k,Hpsi_k)
      end if


*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl       ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     ******************************
*     **** add v_field to vall  ****
*     ******************************
      if (field_exist)
     >  call D3dB_rr_Sum2(1,v_field,dcpl_mb(vall(1)))


*     **********************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r]*psi_r  **** 
*     ********************************************** 
      indx1 = 1
      indx2 = 1
      done  = .false.
      do while (.not.done)
         if (indx1.le.n2(ispin)) then
            ms = 1
            if (indx1.gt.n2(1)) ms = 2

*           **** vpsi_r = (vall_r+xcp)*psi_r ****
            call D3dB_r_Copy(1,dcpl_mb(vall(1)),dcpl_mb(vpsi(1)))

            call D3dB_rr_Mul2(1,psi_r(1+(indx1-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))


*           **** put vpsi_r in k-space ****
            call D3dB_rc_pfft3f_queuein(1,dcpl_mb(vpsi(1)))
            indx1 = indx1 + 1

         end if

         if ((D3dB_rc_pfft3_queue_filled()).or.
     >       (indx1.gt.n2(ispin))) then
            call D3dB_rc_pfft3f_queueout(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(indx2-1)*npack1))
            indx2 = indx2 + 1
         end if

         done = ((indx1.gt.n2(ispin)).and.(indx2.gt.n2(ispin)))
      end do
 

      call nwpw_timing_start(12)

      value = .true.
      value = value.and.BA_pop_stack(vpsi(2))
      value = value.and.BA_pop_stack(vall(2))
      if (.not. value) call errquit('psi_H:popping stack memory',1,0)

      call nwpw_timing_end(12)
      return
      end


