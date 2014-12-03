*
* $Id$
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

#include "bafdecls.fh"
#include "paw_ma.fh"
#include "paw_basis.fh"
#include "paw_geom.fh"
#include "paw_vloc_data.fh"

      logical ok
      integer npack0,nfft3d,nkatm


      !*** allocate storage for local pseudopotential ***
      nkatm = ion_nkatm()
      call Pack_npack(0,npack0)
      call D3dB_nfft3d(1,nfft3d)
      ok = my_alloc(mt_dbl,(nkatm*npack0),'vl',vl)
      if (.not.ok) 
     > call errquit('paw_vloc_init:out of heap memory',0,1)


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

#include "bafdecls.fh"
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

#include "bafdecls.fh"
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

      ok = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. ok) call errquit('paw_vloc:out of stack memory',0,0)
      if (move) then
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
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
c        call Pack_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)
        call Pack_cc_Sum2(0,dcpl_mb(vtmp(1)),vl_out)


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
      ok = ok.and.BA_pop_stack(Gz(2))
      ok = ok.and.BA_pop_stack(Gy(2))
      ok = ok.and.BA_pop_stack(Gx(2))
      ok = ok.and.BA_pop_stack(xtmp(2))
      end if
      ok = ok.and.BA_pop_stack(vtmp(2))
      ok = ok.and.BA_pop_stack(exi(2))
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

#include "bafdecls.fh"
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

      ok = BA_push_get(mt_dcpl,npack0,'exi', exi(2), exi(1))
      ok = ok.and.
     >     BA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      if (.not. ok) call errquit('paw_vloc:out of stack memory',0,0)

        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
        ok = ok.and.
     >       BA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
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
      ok = ok.and.BA_pop_stack(Gz(2))
      ok = ok.and.BA_pop_stack(Gy(2))
      ok = ok.and.BA_pop_stack(Gx(2))
      ok = ok.and.BA_pop_stack(xtmp(2))
      ok = ok.and.BA_pop_stack(vtmp(2))
      ok = ok.and.BA_pop_stack(exi(2))
      if (.not. ok) 
     > call errquit('paw_vloc:error popping stack ',0,2)

      call nwpw_timing_end(5)


      return
      end


      subroutine paw_vloc_ptr(ptr)
      implicit none
      integer ptr
#include "paw_vloc_data.fh"

      ptr = vl(1)
      return
      end
