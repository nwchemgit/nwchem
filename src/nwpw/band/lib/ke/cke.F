*
* $Id: cke.F,v 1.11 2007-08-23 23:49:50 bylaska Exp $
*
      subroutine cke_init()
      implicit none 
#include "errquit.fh"

#include "mafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer nfft3d,npack1,G(3)
      integer i,nb,nbrillioun
      real*8  gg,gx,gy,gz,kx,ky,kz
      logical value
 
      integer tmp1(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  c_G_indx, brillioun_nbrillioun
      real*8   brillioun_k
      external c_G_indx, brillioun_nbrillioun
      external brillioun_k

      nbrillioun = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      G(1)= c_G_indx(1)
      G(2)= c_G_indx(2)
      G(3)= c_G_indx(3)

      value = MA_alloc_get(mt_dbl,npack1*nbrillioun,
     >                     'tg',tg_hndl,tg_indx)
      if (.not. value) 
     > call errquit('cke_init:out of heap memory',0,MA_ERR)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('cke_init:out of stack memory',0,MA_ERR)


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
     >                       dbl_mb(tg_indx+(nb-1)*npack1))
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
      integer npack1,nbrillioun
      integer n,nb,indx

*     **** external functions ****
      integer  brillioun_nbrillioun
      external brillioun_nbrillioun
      
      call Cram_max_npack(npack1)
      nbrillioun = brillioun_nbrillioun()

      do nb=1,nbrillioun
      do n=1,(ne(1)+ne(2))
         indx = 1 + (n-1)*npack1 + (nb-1)*npack1*(ne(1)+ne(2))

         call Cram_rc_Mul(nb,
     >                    dbl_mb(tg_indx+(nb-1)*npack1),
     >                    psi1(indx),
     >                    psi2(indx))
      end do
      end do
      

      return 
      end



      subroutine cke_orb(nb,orb1,orb2)
      implicit none
      integer    nb
      complex*16 orb1(*)
      complex*16 orb2(*)

#include "mafdecls.fh"
#include "cke_common.fh"

*     **** local variables ****
      integer npack1
      integer ne1,shifts
      integer cpsi_ne
      logical cpsi_spin_orbit
      external cpsi_spin_orbit
      external cpsi_ne

      call Cram_max_npack(npack1)
      call Cram_rc_Mul(nb,
     >                 dbl_mb(tg_indx+(nb-1)*npack1),
     >                 orb1,
     >                 orb2)
      if (cpsi_spin_orbit()) then
        ne1=cpsi_ne(1)
        shifts=ne1*npack1
        call Cram_rc_Mul(nb,
     >                 dbl_mb(tg_indx+(nb-1)*npack1),
     >                 orb1(shifts+1),
     >                 orb2(shifts+1))
      end if      
      return 
      end


      subroutine cke_ave(ispin,ne,psi1,ave)
      implicit none
#include "errquit.fh"
      integer ispin,ne(2)
      complex*16 psi1(*)
      real*8     ave

#include "mafdecls.fh"
#include "cke_common.fh"


*     **** local variables ****
      integer npack1,np,indx,nb
      integer ms,n,n1(2),n2(2)
      real*8  sum

      integer tmp1(2)
      logical value

      integer  brillioun_nbrillioun
      real*8   brillioun_weight
      external brillioun_nbrillioun
      external brillioun_weight
     
      call Parallel_np(np)
      call Cram_max_npack(npack1)
      value = MA_push_get(mt_dcpl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('cke_ave:out of stack memory',0,MA_ERR)


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      ave = 0.0d0
      do nb=1,brillioun_nbrillioun()
      do ms=1,ispin
         do n=n1(ms),n2(ms)
            indx =  1+ (n-1)*npack1 + (nb-1)*npack1*(ne(1)+ne(2))
            call Cram_rc_Mul(nb,
     >                       dbl_mb(tg_indx+(nb-1)*npack1),
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



      subroutine cke_euv(ispin,ne,psi,euv)
      implicit none 
      integer ispin,ne(2)
      complex*16 psi(*)
      real*8 euv(3,3)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer npack1,nfft3d,nbrillioun,G(2,3)
      integer nb,i,j,ms,n,n1(2),n2(2),indx
      integer u,v,s
      logical value

      real*8 pi,scal,sum,kx,ky,kz,weight
      real*8 hm(3,3),Aus(3,3)
      integer tmp1(2),tmp2(2)

*     **** external functions ****
      integer  c_G_indx, brillioun_nbrillioun
      real*8   brillioun_k,brillioun_weight
      external c_G_indx, brillioun_nbrillioun
      external brillioun_k,brillioun_weight

      real*8   lattice_unitg,lattice_omega,lattice_unita
      external lattice_unitg,lattice_omega,lattice_unita


      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)

      pi   = 4.0d0*datan(1.0d0)
      scal = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      



      nbrillioun = brillioun_nbrillioun()
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) 
     > call errquit('cke_euv:out of stack  memory',0,MA_ERR)

      value = MA_push_get(mt_dbl,npack1,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)

      value = MA_push_get(mt_dcpl,npack1,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) 
     > call errquit('cke_euv:out of stack memory',0,MA_ERR)

    

*     **** calculate Aus = Sum(k)Sum(n)Sum(G) psi(G+k,n)**2 (G(u)+k(u))(G(s)+k(s)) ****
      call dcopy(9,0.0d0,0,Aus,1)    
      do nb=1,nbrillioun
        weight = brillioun_weight(nb)
        kx = brillioun_k(1,nb)
        ky = brillioun_k(2,nb)
        kz = brillioun_k(3,nb)
        call dcopy(nfft3d,dbl_mb(c_G_indx(1)),1,dbl_mb(G(1,1)),1)
        call dcopy(nfft3d,dbl_mb(c_G_indx(2)),1,dbl_mb(G(1,2)),1)
        call dcopy(nfft3d,dbl_mb(c_G_indx(3)),1,dbl_mb(G(1,3)),1)
        do i = 1,nfft3d
          dbl_mb(G(1,1)+i-1) =  dbl_mb(G(1,1)+i-1)+kx
          dbl_mb(G(1,2)+i-1) =  dbl_mb(G(1,2)+i-1)+ky
          dbl_mb(G(1,3)+i-1) =  dbl_mb(G(1,3)+i-1)+kz
        end do
        call Cram_r_pack(nb,dbl_mb(G(1,1)))
        call Cram_r_pack(nb,dbl_mb(G(1,2)))
        call Cram_r_pack(nb,dbl_mb(G(1,3)))

        do u=1,3
        do s=u,3
          call Cram_rr_Mul(nb,dbl_mb(G(1,u)),
     >                        dbl_mb(G(1,s)),
     >                        dbl_mb(tmp1(1)))
          do ms=1,ispin
          do n=n1(ms),n2(ms)    
            indx = 1 + (n-1)*npack1 + (nb-1)*npack1*(ne(1)+ne(2))
            call Cram_rc_Mul(nb,dbl_mb(tmp1(1)),
     >                       psi(indx),
     >                       dcpl_mb(tmp2(1)))
            call Cram_cc_dot(nb,psi(indx),
     >                       dcpl_mb(tmp2(1)),
     >                       sum)
           
             Aus(u,s) = Aus(u,s) + sum*weight
          end do
          end do
        end do
        end do
      end do

      do u=1,3
      do s=u+1,3
         Aus(s,u) = Aus(u,s)
      end do
      end do
      if (ispin.eq.1) call dscal(9,2.0d0,Aus,1)

*     *** calculate euv = -Sum(s) hm(s,v)*Aus(u,s)
      call dcopy(9,0.0d0,0,euv,1)
      do v=1,3
      do u=1,3
         do s=1,3
            euv(u,v) = euv(u,v) - Aus(u,s)*hm(s,v) 
         end do
      end do
      end do

      value = MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      value = value.and.MA_pop_stack(G(2,3)) 
      value = value.and.MA_pop_stack(G(2,2)) 
      value = value.and.MA_pop_stack(G(2,1))
      if (.not. value) call errquit('error poping stack memory',0,
     &       MA_ERR)
      return
      end




