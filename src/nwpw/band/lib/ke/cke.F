*
* $Id: cke.F,v 1.2 2002-01-11 18:26:00 bylaska Exp $
*
      subroutine cke_init()
      implicit none 

#include "mafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer nfft3d,G(3)
      integer i,nb,nbrillioun
      real*8  gg,gx,gy,gz,kx,ky,kz
      logical value
 
      integer tmp1(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx, brillioun_nbrillioun
      real*8   brillioun_k
      external G_indx, brillioun_nbrillioun
      external brillioun_k

      nbrillioun = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      G(1)= G_indx(1)
      G(2)= G_indx(2)
      G(3)= G_indx(3)

      value = MA_alloc_get(mt_dbl,nfft3d*nbrillioun,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) call errquit('out of heap memory',0)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)


      do nb=1,nbrillioun
         kx = brillioun_k(1,nb)
         ky = brillioun_k(2,nb)
         kz = brillioun_k(3,nb)
         do i = 1,nfft3d
            gx  =  dbl_mb(G(1)+i-1)+kx
            gy  =  dbl_mb(G(2)+i-1)+ky
            gz  =  dbl_mb(G(3)+i-1)+kz
            gg  = (gx*gx + gy*gy + gz*gz)
            dbl_mb(tmp1(1)+i-1) = -0.5d0*gg
         end do

         call Cram_r_pack(nb,dbl_mb(tmp1(1)))
         call Cram_r_Copy(nb,dbl_mb(tmp1(1)),
     >                       dbl_mb(tg_indx+(nb-1)*nfft3d))
      end do
      

      value = MA_pop_stack(tmp1(2))
      return
      end

      subroutine cke_end()
      implicit none 
#include "mafdecls.fh"
#include "cke_common.fh"

      logical value

      value = MA_free_heap(tg_hndl)
      return
      end


      subroutine cke(ispin,ne,psi1,psi2)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi1(*)
      complex*16 psi2(*)

#include "mafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer nfft3d,nbrillioun
      integer n,nb,indx

*     **** external functions ****
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      
      call C3dB_nfft3d(1,nfft3d)
      nbrillioun = brillioun_nbrillioun()

      do nb=1,nbrillioun
      do n=1,(ne(1)+ne(2))
         indx = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2))

         call Cram_rc_Mul(nb,
     >                    dbl_mb(tg_indx+(nb-1)*nfft3d),
     >                    psi1(indx),
     >                    psi2(indx))
      end do
      end do
      

      return 
      end

      subroutine cke_ave(ispin,ne,psi1,ave)
      implicit none
      integer ispin,ne(2)
      complex*16 psi1(*)
      real*8     ave

#include "mafdecls.fh"
#include "cke_common.fh"


*     **** local variables ****
      integer nfft3d,np,indx,nb
      integer ms,n,n1(2),n2(2)
      real*8  sum

      integer tmp1(2)
      logical value

      integer  brillioun_nbrillioun
      real*8   brillioun_weight
      external brillioun_nbrillioun
      external brillioun_weight
     
      call Parallel_np(np)
      call C3dB_nfft3d(1,nfft3d)
      value = MA_push_get(mt_dcpl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0)


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      ave = 0.0d0
      do nb=1,brillioun_nbrillioun()
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            indx =  1+ (n-1)*nfft3d + (nb-1)*nfft3d*(ne(1)+ne(2))
            call Cram_rc_Mul(nb,
     >                       dbl_mb(tg_indx+(nb-1)*nfft3d),
     >                       psi1(indx),
     >                       dcpl_mb(tmp1(1)))
            call Cram_cc_idot(nb,psi1(indx),
     >                       dcpl_mb(tmp1(1)),
     >                       sum)
            ave = ave + sum*brillioun_weight(nb)
         end do
      end do
      end do
      if (np.gt.1) call C3dB_SumAll(ave)
      if (ispin.eq.1) ave = 2.0d0*ave
      ave = -ave

      value = MA_pop_stack(tmp1(2))
      return 
      end 







