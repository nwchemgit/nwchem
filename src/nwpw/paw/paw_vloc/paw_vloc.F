*
* $Id: paw_vloc.F,v 1.4 2004-08-13 01:16:43 bylaska Exp $
*


!**************************************************
!
!       Name: paw_vloc_init
!
!       Purpose: 
!
!       Created: Wed Feb 12 15:18:24 PST 2003
!**************************************************
      subroutine paw_vloc_init()
      implicit none

#include "mafdecls.fh"
#include "paw_ma.fh"
#include "bessel_transform.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_vloc_data.fh"

      !**** local variables ****
      real*8 small
      parameter (small=1.0d-9)

      logical ok
      integer i,ia,npack0,nfft3d,nr,nkatm
      integer ps_ptr,i_rgrid
      integer Gx(2),Gy(2),Gz(2),f(2)
      real*8 log_amesh,gg,fourpi

      !*** external functions ****
      integer  G_indx
      external G_indx

      fourpi = 16.0d0*datan(1.0d0)

      !*** allocate storage for local pseudopotential ***
      nkatm = ion_nkatm()
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      ok = my_alloc(mt_dbl,(nkatm*npack0),'vl',vl)
      if (.not.ok) 
     > call errquit('paw_vloc_init:out of heap memory',0,1)

      !*** allocate temporary space ***
      ok = MA_push_get(mt_dbl,nfft3d,'Gx',Gx(2),Gx(1))
      ok = ok.and.
     >        MA_push_get(mt_dbl,nfft3d,'Gy',Gy(2),Gy(1))
      ok = ok.and.
     >        MA_push_get(mt_dbl,nfft3d,'Gz',Gz(2),Gz(1))
      if (.not.ok)
     > call errquit('paw_vloc_init: out of stack memory',0,0)

      call D3dB_t_Copy(1,dbl_mb(G_indx(1)),dbl_mb(Gx(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(2)),dbl_mb(Gy(1)))
      call D3dB_t_Copy(1,dbl_mb(G_indx(3)),dbl_mb(Gz(1)))
      call Pack_t_pack(0,dbl_mb(Gx(1)))
      call Pack_t_pack(0,dbl_mb(Gy(1)))
      call Pack_t_pack(0,dbl_mb(Gz(1)))  


      !*** format the pseudopotentials onto the reciprocal space ***
      do ia=1,nkatm
        ps_ptr          = paw_basis_i_v_ps(ia)
        i_rgrid         = paw_basis_i_rgrid(ia)
        nr              = int_mb(paw_basis_i_ngrid(ia))
        log_amesh       = dbl_mb(paw_basis_i_log_amesh(ia))

        ok = MA_push_get(mt_dbl,nr,'f',f(2),f(1))
        if (.not.ok)
     >  call errquit('paw_vloc_init: out of stack memory',0,1)
        do i=1,nr
          dbl_mb(f(1)+i-1) = dbl_mb(i_rgrid+i-1)
     >                      *dbl_mb(ps_ptr+i-1)
        end do

        do i=1,npack0
          gg  = dsqrt( dbl_mb(Gx(1)+i-1)*dbl_mb(Gx(1)+i-1)
     >               + dbl_mb(Gy(1)+i-1)*dbl_mb(Gy(1)+i-1)
     >               + dbl_mb(Gz(1)+i-1)*dbl_mb(Gz(1)+i-1) )


          if (gg.gt.small) then
          dbl_mb(vl(1)+(ia-1)*npack0 + i-1)
     >       = fourpi
     >        *spher_bessel_transform(0,nr,log_amesh,
     >                                dbl_mb(i_rgrid),
     >                                dbl_mb(f(1)),
     >                                gg)
          else
          dbl_mb(vl(1)+(ia-1)*npack0 + i-1)
     >       = fourpi
     >        *spher_bessel0_transform(0,nr,log_amesh,
     >                                dbl_mb(i_rgrid),
     >                                dbl_mb(f(1)))

          end if
        end do !*i*
       
        ok = MA_pop_stack(f(2))
        if (.not.ok)
     >  call errquit('paw_vloc_init: error popping stack',0,2)

      end do !*ia*

      !*** deallocate temporary space ***
      ok =        MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      if (.not.ok)
     >  call errquit('paw_vloc_init: error popping stack',0,3)

      return
      end


!**************************************************
!
!       Name: paw_vloc_end
!
!       Purpose: 
!
!       Created: Wed Feb 12 15:18:24 PST 2003
!**************************************************
      subroutine paw_vloc_end()
      implicit none

#include "mafdecls.fh"
#include "paw_ma.fh"
#include "paw_vloc_data.fh"

      !**** local variables ****
      logical ok

      !*** deallocate storage for local pseudopotential ***
      ok = my_dealloc(vl)
      if (.not.ok) 
     > call errquit('paw_vloc_end:out of heap memory',0,1)

      return
      end

!**************************************************
!
!       Name: paw_vloc
!
!       Purpose: 
!
!       Created: Wed Feb 12 15:18:24 PST 2003
!**************************************************
      subroutine paw_vloc(vl_out,move,dng,fion)
      implicit none
      complex*16 vl_out(*)
      logical    move
      complex*16 dng(*)
      real*8     fion(3,*)

#include "mafdecls.fh"
#include "paw_ma.fh"
#include "paw_vloc_data.fh"


*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical ok

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      ok = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. ok) call errquit('paw_vloc:out of stack memory',0,0)
      if (move) then
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        if (.not. ok) call errquit('paw_vloc:out of stack memory',0,1)
        G(1)  = G_indx(1)
        G(2)  = G_indx(2)
        G(3)  = G_indx(3)

        !**** define Gx,Gy and Gz in packed space ****
        call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Pack_t_pack(0,dbl_mb(Gx(1)))
        call Pack_t_pack(0,dbl_mb(Gy(1)))
        call Pack_t_pack(0,dbl_mb(Gz(1)))
      end if

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac_pack(0,ii,dcpl_mb(exi(1)))

*       **** add to local psp ****
        call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))
        call Pack_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)


        if (move) then

          do i=1,npack0
             dbl_mb(xtmp(1)+i-1)
     >          = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >           - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
          end do
         call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
         call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
         call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))

       end if


      end do
      
      ok = .true.
      if (move) then
      ok = ok.and.MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      ok = ok.and.MA_pop_stack(xtmp(2))
      end if
      ok = ok.and.MA_pop_stack(vtmp(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not. ok) 
     > call errquit('paw_vloc:error popping stack ',0,2)

      call nwpw_timing_end(5)


      return
      end

!**************************************************
!
!       Name: paw_vloc_force
!
!       Purpose: 
!
!       Created: Wed Feb 12 15:18:24 PST 2003
!**************************************************
      subroutine paw_vloc_force(dng,fion)
      implicit none
      complex*16 dng(*)
      real*8     fion(3,*)

#include "mafdecls.fh"
#include "paw_ma.fh"
#include "paw_vloc_data.fh"


*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical ok

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call nwpw_timing_start(5)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      ok = MA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      ok = ok.and.
     >     MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. ok) call errquit('paw_vloc:out of stack memory',0,0)

        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        ok = ok.and.
     >       MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
        if (.not. ok) call errquit('paw_vloc:out of stack memory',0,1)
        G(1)  = G_indx(1)
        G(2)  = G_indx(2)
        G(3)  = G_indx(3)

        !**** define Gx,Gy and Gz in packed space ****
        call D3dB_t_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
        call D3dB_t_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
        call D3dB_t_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
        call Pack_t_pack(0,dbl_mb(Gx(1)))
        call Pack_t_pack(0,dbl_mb(Gy(1)))
        call Pack_t_pack(0,dbl_mb(Gz(1)))

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call strfac_pack(0,ii,dcpl_mb(exi(1)))

*       **** add to local psp ****
        call Pack_tc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))


*       **** force evaluation ****
        do i=1,npack0
           dbl_mb(xtmp(1)+i-1)
     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
        end do
        call Pack_tt_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
        call Pack_tt_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
        call Pack_tt_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))


      end do
      
      ok = .true.
      ok = ok.and.MA_pop_stack(Gz(2))
      ok = ok.and.MA_pop_stack(Gy(2))
      ok = ok.and.MA_pop_stack(Gx(2))
      ok = ok.and.MA_pop_stack(xtmp(2))
      ok = ok.and.MA_pop_stack(vtmp(2))
      ok = ok.and.MA_pop_stack(exi(2))
      if (.not. ok) 
     > call errquit('paw_vloc:error popping stack ',0,2)

      call nwpw_timing_end(5)


      return
      end


