*
* $Id: cpsp.F,v 1.3 2001-12-20 00:30:50 bylaska Exp $
*

#define TCGMSG

      subroutine cpsp_init()      
      implicit none

#include "mafdecls.fh"



*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      integer nfft3d,npack0,nbrillioun
      logical value

*     **** external functions *****
      integer  ion_nkatm, brillioun_nbrillioun
      external ion_nkatm, brillioun_nbrillioun
 


      npsp = ion_nkatm()
      nbrillioun = brillioun_nbrillioun()

      call Cram_npack(0,npack0)
      call C3dB_nfft3d(1,nfft3d)

      value = MA_alloc_get(mt_dbl,(npsp*npack0),
     >                    'vl',vl(2),vl(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(lmmax_max*npsp*nfft3d*nbrillioun),
     >                    'vnl',vnl(2),vnl(1))
      value = value.and.
     >         MA_alloc_get(mt_dbl,(lmmax_max*npsp),
     >                    'vnlnrm',vnlnrm(2),vnlnrm(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'zv',zv(2),zv(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp),'amass',amass(2),amass(1))
      value = value.and.
     >        MA_alloc_get(mt_dbl,(npsp*(lmax_max+1)),'rc',rc(2),rc(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmmax',lmmax(2),lmmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'lmax',lmax(2),lmax(1))
      value = value.and.
     >        MA_alloc_get(mt_int,(npsp),'locp',locp(2),locp(1))

      if (.not. value) call errquit('out of heap memory',0)

      call dcopy(npsp*npack0,          0.0d0,0,dbl_mb(vl(1)), 1)
      call dcopy(lmmax_max*npsp*nfft3d*nbrillioun,
     >           0.0d0,0,dbl_mb(vnl(1)),1)
      call dcopy(lmmax_max*npsp,       0.0d0,0,dbl_mb(vnlnrm(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(zv(1)),1)
      call dcopy(npsp,                 0.0d0,0,dbl_mb(amass(1)),1)
      call dcopy(npsp*(lmax_max+1),    0.0d0,0,dbl_mb(rc(1)),1)

*     **** allocate semicore data ****
      call c_semicore_init()

      return
      end




      subroutine cpsp_end()      
      implicit none

#include "mafdecls.fh"


*     **** common blocks ****
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      logical value

*     **** external functions ****

*     **** deallocate semicore data ****
      call c_semicore_end()

      value = MA_free_heap(vl(2))
      value = MA_free_heap(vnl(2))
      value = MA_free_heap(vnlnrm(2))
      value = MA_free_heap(zv(2))
      value = MA_free_heap(amass(2))
      value = MA_free_heap(rc(2))
      value = MA_free_heap(lmmax(2))
      value = MA_free_heap(lmax(2))
      value = MA_free_heap(locp(2))

      return
      end



*     ***********************************
*     *					*
*     *	 	   cpsp_zv		*
*     *					*
*     ***********************************

      real*8 function cpsp_zv(ia)
      implicit none
      integer ia

#include "mafdecls.fh"


*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      cpsp_zv = dbl_mb(zv(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_amass		*
*     *					*
*     ***********************************

      real*8 function cpsp_amass(ia)
      implicit none
      integer ia

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      cpsp_amass = dbl_mb(amass(1)+ia-1)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_rc		*
*     *					*
*     ***********************************

      real*8 function cpsp_rc(i,ia)
      implicit none
      integer i,ia

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


c     cpsp_rc = rc(i,ia)
      cpsp_rc = dbl_mb(rc(1) + i + (lmax_max+1)*(ia-1))
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_atom		*
*     *					*
*     ***********************************

      character*2 function cpsp_atom(ia)
      implicit none
      integer  ia

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      cpsp_atom = atom(ia)
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_lmmax		*
*     *					*
*     ***********************************

      integer function cpsp_lmmax(ia)
      implicit none
      integer  ia

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx),lmax(2),locp(2)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      cpsp_lmmax = int_mb(lmmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_lmax		*
*     *					*
*     ***********************************

      integer function cpsp_lmax(ia)
      implicit none
      integer  ia

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      cpsp_lmax = int_mb(lmax(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_locp		*
*     *					*
*     ***********************************

      integer function cpsp_locp(ia)
      implicit none
      integer  ia

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

      cpsp_locp = int_mb(locp(1)+ia-1)
      return
      end

*     ***********************************
*     *					*
*     *	 	   cpsp_npsp		*
*     *					*
*     ***********************************

      integer function cpsp_npsp()
      implicit none


*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


      cpsp_npsp = npsp
      return
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_v_local  		*
*     *					*
*     ***********************************

      subroutine cpsp_v_local(vl_out,move,dng,fion)
      implicit none
      complex*16 vl_out(*)
      logical    move
      complex*16 dng(*)
      real*8     fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - cstrfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call current_second(tim1)  
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Cram_r_pack(0,dbl_mb(Gx(1)))
      call Cram_r_pack(0,dbl_mb(Gy(1)))
      call Cram_r_pack(0,dbl_mb(Gz(1)))

      call dcopy((2*npack0),0.0d0,0,vl_out,1)
      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call cstrfac(ii,dcpl_mb(exi(1)))
        call Cram_c_pack(0,dcpl_mb(exi(1)))
        
*       **** add to local psp ****
        call Cram_rc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))
        call Cram_cc_Sum(0,vl_out,dcpl_mb(vtmp(1)),vl_out)


        if (move) then
       
          do i=1,npack0
             dbl_mb(xtmp(1)+i-1) 
     >          = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >           - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
          end do
         call Cram_rr_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
         call Cram_rr_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
         call Cram_rr_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
   
       end if
        

      end do
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call current_second(tim2)
      times(5) = times(5) + (tim2-tim1)
      return 
      end


*     ***********************************
*     *					*
*     *	 	   cpsp_f_vlocal  		*
*     *					*
*     ***********************************

      subroutine cpsp_f_vlocal(dng,fion)
      implicit none
      complex*16 dng(*)
      real*8     fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - cstrfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


*     *** local variables ***
      integer nfft3d,npack0
      integer i,ii,ia
      integer exi(2),vtmp(2),xtmp(2),G(3)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      integer  G_indx,ion_nion,ion_katm
      external G_indx,ion_nion,ion_katm

      call current_second(tim1)  
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Cram_r_pack(0,dbl_mb(Gx(1)))
      call Cram_r_pack(0,dbl_mb(Gy(1)))
      call Cram_r_pack(0,dbl_mb(Gz(1)))

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call cstrfac(ii,dcpl_mb(exi(1)))
        call Cram_c_pack(0,dcpl_mb(exi(1)))
        
*       **** add to local psp ****
        call Cram_rc_Mul(0,dbl_mb(vl(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(vtmp(1)))

         do i=1,npack0
           dbl_mb(xtmp(1)+i-1) 
     >        = dimag(dng(i))* dble(dcpl_mb(vtmp(1)+i-1))
     >         - dble(dng(i))*dimag(dcpl_mb(vtmp(1)+i-1))
        end do

        call Cram_rr_dot(0,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),fion(1,ii))
        call Cram_rr_dot(0,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),fion(2,ii))
        call Cram_rr_dot(0,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),fion(3,ii))
 
      end do
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call current_second(tim2)
      times(5) = times(5) + (tim2-tim1)
      return 
      end




*     ***********************************
*     *					*
*     *	 	   cpsp_v_nonlocal  		*
*     *					*
*     ***********************************

      subroutine cpsp_v_nonlocal(ispin,ne,
     >                           psi1,psi2,move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)
      logical move
      real*8 fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - cstrfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


*     *** local variables ***
      integer nfft3d,G(3),npack1
      integer i,ii,ia,l,n,nn,nb,neall,nbrill
      real*8  omega,weight
      complex*16 ctmp,cxr
      integer exi(2),tmp4(2),vtmp(2),xtmp(2),zsw1(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,G_indx,brillioun_nbrillioun
      real*8   lattice_omega,brillioun_weight
      external is_sORd
      external ion_nion,ion_katm,G_indx,brillioun_nbrillioun
      external lattice_omega,brillioun_weight

      call current_second(tim1)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'tmp4', tmp4(2), tmp4(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nn,'zsw1',zsw1(2),zsw1(1))
      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)

        call cstrfac(ii,dcpl_mb(tmp4(1)))

        do nb=1,nbrill
           call Cram_npack(nb,npack1)
           call C3dB_c_Copy(1,dcpl_mb(tmp4(1)),dcpl_mb(exi(1)))
           call Cram_c_pack(nb,dcpl_mb(exi(1)))
           call cstrfac_k(ii,nb,cxr)
           call Cram_c_ZMul(nb,cxr,dcpl_mb(exi(1)),dcpl_mb(exi(1)))


*       **** structure factor and local pseudopotential ****

        do l=1,int_mb(lmmax(1)+ia-1)

*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
              call Cram_rc_Mul(nb,dbl_mb(vnl(1)+(l-1)*nfft3d
     >                           +(ia-1)*nfft3d*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              call Cram_irc_Mul(nb,dbl_mb(vnl(1)+(l-1)*nfft3d
     >                           +(ia-1)*nfft3d*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))
           end if


*         **** do kleinman-bylander multiplication ****
          call Cram_cc_nzdot(nb,nn,
     >                      psi1(1+(nb-1)*neall*nfft3d),
     >                      dcpl_mb(vtmp(1)),
     >                      dcpl_mb(zsw1(1)))

          do n=1,nn
          
             dcpl_mb(zsw1(1)+n-1) = dcpl_mb(zsw1(1)+n-1)
     >                      /(omega*dbl_mb(vnlnrm(1)
     >                         +(l-1)+(ia-1)*lmmax_max))
             call Cram_cc_zaxpy(nb,(-dcpl_mb(zsw1(1)+n-1)),
     >                          dcpl_mb(vtmp(1)),
     >                psi2(1+(n-1)*nfft3d +(nb-1)*neall*nfft3d)) 

          end do

          if (move) then
             weight = brillioun_weight(nb)
             do n=1,nn
                if (ispin.eq.1) 
     >            dcpl_mb(zsw1(1)+n-1)=2.0d0*dcpl_mb(zsw1(1)+n-1)

                do i=1,npack1
                   ctmp = psi1(i+(n-1)*nfft3d)
     >                  *dconjg(dcpl_mb(vtmp(1)+i-1))
     >                  *dconjg(dcpl_mb(zsw1(1)+n-1))
                   dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
                end do

*               **** define Gx,Gy and Gz in packed space ****
                call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
                call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
                call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
                call Cram_r_pack(nb,dbl_mb(Gx(1)))
                call Cram_r_pack(nb,dbl_mb(Gy(1)))
                call Cram_r_pack(nb,dbl_mb(Gz(1)))
                call Cram_rr_idot(nb,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
                call Cram_rr_idot(nb,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+1+3*(n-1)))
                call Cram_rr_idot(nb,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+2+3*(n-1)))
             end do

             call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
             do n=1,nn
                fion(1,ii) = fion(1,ii) 
     >                     + 2.0d0*weight
     >                            *dbl_mb(sum(1)+3*(n-1))
                fion(2,ii) = fion(2,ii) 
     >                     + 2.0d0*weight
     >                            *dbl_mb(sum(1)+1+3*(n-1))
                fion(3,ii) = fion(3,ii) 
     >                     + 2.0d0*weight
     >                            *dbl_mb(sum(1)+2+3*(n-1))
             end do
          end if
       
        
        end do

      end do
      end do

      value = MA_pop_stack(sum(2))
      value = MA_pop_stack(zsw1(2))
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(tmp4(2))
      value = MA_pop_stack(exi(2))

      call current_second(tim2)
      times(6) = times(6) + (tim2-tim1)

      return 
      end



*     ***********************************
*     *					*
*     *	 	   cpsp_f_vnonlocal 		*
*     *					*
*     ***********************************

      subroutine cpsp_f_vnonlocal(ispin,ne,psi1,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      real*8 fion(3,*)

#include "mafdecls.fh"

*     **** common blocks ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - cstrfac_time
*     9 - maskC_time
      real*8 tim1,tim2
      real*8 times(50)
      common / timing / times


*     *** local variables ***
      integer nfft3d,G(3),npack1
      integer i,ii,ia,l,n,nn
      real*8  omega
      complex*16 ctmp
      integer exi(2),vtmp(2),xtmp(2),sw1(2),sum(2)
      integer Gx(2),Gy(2),Gz(2)
      logical value

*     **** external functions ****
      logical  is_sORd
      integer  ion_nion,ion_katm,G_indx
      real*8   lattice_omega
      external is_sORd
      external ion_nion,ion_katm,G_indx
      external lattice_omega

      call current_second(tim1)
      
*     **** allocate local memory ****
      nn = ne(1)+ne(2)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(1,npack1)
      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vtmp',vtmp(2),vtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'xtmp',xtmp(2),xtmp(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gx',Gx(2),Gx(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gy',Gy(2),Gy(1))
      value = value.and.
     >        MA_push_get(mt_dbl, nfft3d,'Gz',Gz(2),Gz(1))
      value = value.and.
     >        MA_push_get(mt_dbl,nn,'sw1',sw1(2),sw1(1))
      value = value.and.
     >      MA_push_get(mt_dbl,3*nn,'sum',sum(2),sum(1))
      if (.not. value) call errquit('out of stack memory',0)
      G(1)  = G_indx(1)
      G(2)  = G_indx(2)
      G(3)  = G_indx(3)

*     **** define Gx,Gy and Gz in packed space ****
      call C3dB_r_Copy(1,dbl_mb(G(1)),dbl_mb(Gx(1)))
      call C3dB_r_Copy(1,dbl_mb(G(2)),dbl_mb(Gy(1)))
      call C3dB_r_Copy(1,dbl_mb(G(3)),dbl_mb(Gz(1)))
      call Cram_r_pack(1,dbl_mb(Gx(1)))
      call Cram_r_pack(1,dbl_mb(Gy(1)))
      call Cram_r_pack(1,dbl_mb(Gz(1)))

      omega = lattice_omega()

      do ii=1,ion_nion()
        ia=ion_katm(ii)

*       **** structure factor and local pseudopotential ****
        call cstrfac(ii,dcpl_mb(exi(1)))
        call Cram_c_pack(1,dcpl_mb(exi(1)))

        do l=1,int_mb(lmmax(1)+ia-1)

*          **** phase factor does not matter therefore ****
*          **** (-i)^l is the same as (i)^l in the     ****
*          **** Rayleigh scattering formula            ****

*          *** current function is s or d ****
           if (is_sORd(l,int_mb(lmax(1)+ia-1),
     >                   int_mb(locp(1)+ia-1))
     >        ) then
              call Cram_rc_Mul(1,dbl_mb(vnl(1)+(l -1)*npack1
     >                           +(ia-1)*npack1*lmmax_max),
     >                        dcpl_mb(exi(1)),
     >                        dcpl_mb(vtmp(1)))

*          *** current function is p or f ****
           else
              do i=1,npack1
                 dcpl_mb(vtmp(1)+i-1) 
     >                     = dbl_mb(vnl(1) + (i-1)
     >                                     + (l-1)*npack1
     >                          + (ia-1)*npack1*lmmax_max)
     >                * dcmplx(-dimag(dcpl_mb(exi(1)+i-1)),
     >                          dble(dcpl_mb(exi(1)+i-1)))
              end do
           end if

*         **** do kleinman-bylander multiplication ****
          call Cram_cc_nzdot(1,nn,
     >                      psi1,
     >                      dcpl_mb(vtmp(1)),
     >                      dbl_mb(sw1(1)))

          do n=1,nn
             dbl_mb(sw1(1)+n-1) = dbl_mb(sw1(1)+n-1)
     >                      /(omega*dbl_mb(vnlnrm(1)
     >                         +(l-1)+(ia-1)*lmmax_max))
 
          end do

          
          do n=1,nn
             if (ispin.eq.1) 
     >         dbl_mb(sw1(1)+n-1)=2.0d0*dbl_mb(sw1(1)+n-1)

               do i=1,npack1
                  ctmp = psi1(i+(n-1)*npack1)
     >                 *dconjg(dcpl_mb(vtmp(1)+i-1))
                  dbl_mb(xtmp(1)+i-1) = dimag(ctmp)
               end do

               call Cram_rr_idot(1,dbl_mb(Gx(1)),dbl_mb(xtmp(1)),
     >                            dbl_mb(sum(1)+3*(n-1)))
               call Cram_rr_idot(1,dbl_mb(Gy(1)),dbl_mb(xtmp(1)),
     >                           dbl_mb(sum(1)+1+3*(n-1)))
               call Cram_rr_idot(1,dbl_mb(Gz(1)),dbl_mb(xtmp(1)),
     >                           dbl_mb(sum(1)+2+3*(n-1)))

          end do

          call C3dB_Vector_Sumall(3*(nn),dbl_mb(sum(1)))
  
          do n=1,nn
             fion(1,ii) = fion(1,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+3*(n-1))
             fion(2,ii) = fion(2,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+1+3*(n-1))
             fion(3,ii) = fion(3,ii) 
     >                  + 2.0d0*dbl_mb(sw1(1)+n-1)
     >                         *dbl_mb(sum(1)+2+3*(n-1))
          end do
       
        
        end do

      end do

      value = MA_pop_stack(sum(2))
      value = MA_pop_stack(sw1(2))
      value = MA_pop_stack(Gz(2))
      value = MA_pop_stack(Gy(2))
      value = MA_pop_stack(Gx(2))
      value = MA_pop_stack(xtmp(2))
      value = MA_pop_stack(vtmp(2))
      value = MA_pop_stack(exi(2))

      call current_second(tim2)
      times(6) = times(6) + (tim2-tim1)

      return 
      end

      subroutine cpsp_read(fname,
     >                       version,
     >                       nfft,unita,
     >                       atom,amass,zv,lmmax,lmax,locp,
     >                       rc,vnlnrm,
     >                       nfft3d,npack0,lmmax_max,
     >                       vl,vnl,
     >                       semicore,rcore,ncore,
     >                       tmp,tmp2,
     >                       ierr)
      implicit none 
      character*20 fname
      integer version
      integer nfft(3)
      real*8  unita(3,3)
      character*2 atom
      real*8 amass,zv
      integer lmmax
      integer lmax
      integer locp
      real*8 rc(*),vnlnrm(*)
      integer nfft3d,npack0,lmmax_max
      real*8 vl(*)
      real*8 vnl(nfft3d,lmmax_max,*)
      logical semicore
      real*8  rcore
      real*8  ncore(*)
      complex*16 tmp(*)
      real*8     tmp2(*)
      integer ierr

#ifdef TCGMSG
#include "tcgmsg.fh"
#include "msgtypesf.h"
#endif

*    *** local variables ***
      integer MASTER,taskid
      parameter(MASTER=0)
      integer n,l,nb
      integer msglen
      integer iatom(2)
      character*255 full_filename
      real*8 kv(3)
      integer nbrillioun

      integer  brillioun_nbrillioun
      real*8   brillioun_k
      external brillioun_nbrillioun
      external brillioun_k


      call Parallel_taskid(taskid)
  
*     **** open fname binary file ****
      if (taskid.eq.MASTER) then
         call util_file_name_noprefix(fname,.false.,
     >                             .false.,
     >                       full_filename)
         l = index(full_filename,' ') - 1
         call openfile(5,full_filename,l,'r',l)
         call iread(5,version,1)
         call iread(5,nfft,3)
         call dread(5,unita,9)
         call cread(5,atom,2)
         call dread(5,amass,1)
         call dread(5,zv,1)
         call iread(5,lmax,1)
         call iread(5,locp,1)
         lmmax=(lmax+1)**2 - (2*locp+1)
         amass = amass*1822.89d0
         call dread(5,rc,lmax+1)
         call dread(5,vnlnrm,lmmax)
         call dread(5,rcore,1)
         call iread(5,nbrillioun,1)
         ierr = 0
         if (nbrillioun.eq.brillioun_nbrillioun()) then
            do nb=1,nbrillioun
               call dread(5,kv,3)
               if ((brillioun_k(1,nb).ne.kv(1)).or.
     >             (brillioun_k(2,nb).ne.kv(2)).or.
     >             (brillioun_k(3,nb).ne.kv(3))) ierr = 1
            end do
         else
            ierr = 1
         end if
      end if

      msglen = 1
      call BRDCST(9+MSGINT,ierr,mitob(msglen),MASTER)
      if (ierr.ne.0) then
         return
      end if

*     **** send header data to all processors ****
#ifdef TCGMSG
      msglen = 1
      call BRDCST(9+MSGINT,version,mitob(msglen),MASTER)
      msglen = 3
      call BRDCST(9+MSGINT,nfft,mitob(msglen),MASTER)
      msglen = 9
      call BRDCST(9+MSGDBL,unita,mdtob(msglen),MASTER)

      iatom(1) = ichar(atom(1:1))
      iatom(2) = ichar(atom(2:2))
      msglen = 2
      call BRDCST(9+MSGCHR,iatom,mitob(msglen),MASTER)
      atom(1:1) = char(iatom(1))
      atom(2:2) = char(iatom(2))

      msglen = 1
      call BRDCST(9+MSGDBL,amass,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGDBL,zv,mdtob(msglen),MASTER)
      msglen = 1
      call BRDCST(9+MSGINT,lmax,mitob(msglen),MASTER)
      call BRDCST(9+MSGINT,locp,mitob(msglen),MASTER)
      lmmax=(lmax+1)**2 - (2*locp+1)

      msglen=lmax+1
      call BRDCST(9+MSGDBL,rc,mdtob(msglen),MASTER)
      msglen=lmmax
      call BRDCST(9+MSGDBL,vnlnrm,mdtob(msglen),MASTER)
      msglen=1
      call BRDCST(9+MSGDBL,rcore,mdtob(msglen),MASTER)

#endif

*     **** determine semicore value ****
      if (rcore.gt.0.0d0) then
         semicore = .true.
      else
         semicore = .false.
      end if


 
*     *** read in vl 3d block ***
      call C3dB_r_read(1,5,tmp2,tmp)
      call Cram_r_pack(0,tmp2)
      call Cram_r_Copy(0,tmp2,vl)

*     **** read in vnl 3d blocks ****
      do nb=1,brillioun_nbrillioun()
      do n=1,lmmax
         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(nb,tmp2)
         call Cram_r_Copy(1,tmp2,vnl(1,n,nb))
      end do
      end do

*     **** read in semicore density block ****
      if (semicore) then
         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+2*npack0))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+3*npack0))

         call C3dB_r_read(1,5,tmp2,tmp)
         call Cram_r_pack(0,tmp2)
         call Cram_r_Copy(0,tmp2,ncore(1+4*npack0))
      end if

*     *** close fname binary file ***
      if (taskid.eq.MASTER) then
c       close(11)
         call closefile(5)
      end if

      ierr = 0
      return
      end

*     ***********************************
*     *									*
*     *	 	  cpsp_readall  				*
*     *									*
*     ***********************************

      subroutine cpsp_readall()
      implicit none
      
#include "mafdecls.fh"


*     **** semicore common block ****
c     real*8  ncore(nfft3d,nkatmx),rcore(nkatmx)
c     logocal semicore(0:nkatmx)
      integer ncore(2),rcore(2)
      integer semicore(2)
      common / c_ccore / ncore,rcore,semicore

*     **** psp common block ****
c     real*8     vl(nfft3d,nkatmx)
c     real*8     vnl(nfft3d,16,nkatmx)
c     real*8     vnlnrm(16,nkatmx)
c     real*8     zv(nkatmx),amass(nkatmx),rc(0:4,nkatmx)
c     integer    lmmax(nkatmx)
      integer    vl(2)
      integer    vnl(2)
      integer    vnlnrm(2)
      integer    zv(2),amass(2),rc(2)
      integer    lmmax(2),lmax(2),locp(2)
      integer    npsp
   
      integer nkatmx,lmmax_max,lmax_max
      parameter (nkatmx=50,lmax_max=4,lmmax_max=16)
      character*2 atom(nkatmx)
      common / psp_block / vl,vnl,vnlnrm,zv,amass,rc,
     >                     lmmax,lmax,locp,npsp,atom


*     **** local variables ****
      integer ngp(3),version,nfft3d,npack0
      integer ia,l
      real*8 unita(3,3)
      character*12 boundry
      integer tmp(2),tmp2(2),ierr
      logical value,found,correct_box
      character*4  element
      character*20 fname

*     **** parallel i/o variable ****
      integer MASTER,taskid
      parameter(MASTER=0)

*     **** external functions ****
      logical      pspw_filefind
      integer      control_ngrid
      real*8       control_unita
      character*12 control_boundry
      character*2  ion_atom
      external     pspw_filefind
      external     control_ngrid
      external     control_unita
      external     control_boundry
      external     ion_atom
      

      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)
      call Parallel_taskid(taskid)

*     *** set semicore(0) *****
      log_mb(semicore(1)) = .false.

      value = MA_push_get(mt_dbl,(2*nfft3d),'tmp',tmp(2),tmp(1))
      if (.not. value) call errquit('out of stack memory',0)

      value = MA_push_get(mt_dbl,(nfft3d),'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0)

*     **** read pseudopotentials ****
      do ia=1,npsp

*      **** define formatted psp name ****
       element = '    '
       element = ion_atom(ia)
       l = index(element,' ') - 1
       fname = element(1:l)//'.cpp'
        

       found = .false.
       do while (.not.found)

         if (pspw_filefind(fname)) then
            call cpsp_read(fname,
     >                  version,
     >                  ngp,unita,
     >                  atom(ia),
     >                  dbl_mb(amass(1)+ia-1),
     >                  dbl_mb(zv(1)+ia-1),
     >                  int_mb(lmmax(1)+ia-1),
     >                  int_mb(lmax(1)+ia-1),
     >                  int_mb(locp(1)+ia-1),
     >                  dbl_mb(rc(1) + (ia-1)*(lmax_max+1)),
     >                  dbl_mb(vnlnrm(1) + (ia-1)*lmmax_max),
     >                  nfft3d,npack0,lmmax_max,
     >                  dbl_mb(vl(1) + (ia-1)*npack0),
     >                  dbl_mb(vnl(1)+ (ia-1)*nfft3d*lmmax_max),
     >                  log_mb(semicore(1)+ia),
     >                  dbl_mb(rcore(1)+ia-1),
     >                  dbl_mb(ncore(1)+ (ia-1)*npack0*5),
     >                  dbl_mb(tmp(1)),dbl_mb(tmp2(1)),
     >                  ierr)


*          **** set semicore(0) ****
           if (log_mb(semicore(1)+ia)) log_mb(semicore(1)) = .true.
           if (ierr.gt.1) go to 9000

*          **************************************************************
*          ***** logic for finding out if psp is correctly formatted ****
*          **************************************************************
           correct_box = .true.
           boundry = control_boundry()
           l =index(boundry,' ') - 1
           if ( (ngp(1).ne.control_ngrid(1)) .or.
     >       (ngp(2).ne.control_ngrid(2)) .or.
     >       (ngp(3).ne.control_ngrid(3)) .or. 
     >       (unita(1,1).ne.control_unita(1,1)) .or.
     >       (unita(2,1).ne.control_unita(2,1)) .or.
     >       (unita(3,1).ne.control_unita(3,1)) .or.
     >       (unita(1,2).ne.control_unita(1,2)) .or.
     >       (unita(2,2).ne.control_unita(2,2)) .or.
     >       (unita(3,2).ne.control_unita(3,2)) .or.
     >       (unita(1,3).ne.control_unita(1,3)) .or.
     >       (unita(2,3).ne.control_unita(2,3)) .or.
     >       (unita(3,3).ne.control_unita(3,3)) .or.
     >       ((boundry(1:l).eq.'periodic').and.(version.ne.3)).or.
     >       ((boundry(1:l).eq.'aperiodic').and.(version.ne.4))) then
              correct_box = .false.
              if (taskid.eq.MASTER) then
              write(6,*) "pseudopotential is not correctly formatted:",
     >                    fname
              end if
           end if
           if (correct_box) found = .true.
           if (ierr.eq.1)   then
              found = .false.
              if (taskid.eq.MASTER) then
              write(6,*) "pseudopotential is not correctly formatted-",
     >                   "bad brillioun zone:",fname
              end if
           end if

         end if

*        **** generate formatted pseudopotential atom.cpp *****
         if (.not.found) then
             call cpsp_formatter_auto(ion_atom(ia))
         end if

       end do !***do while ****


      end do
 9000 value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp(2))

      return
      end

