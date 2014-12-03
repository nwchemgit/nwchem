

*     *************************************************
*     *                                               *
*     *          nwpw_compcharge_atomic_coulomb       *
*     *                                               *
*     *************************************************

*     This routine computes the atomic coulomb interactions

      subroutine nwpw_compcharge_atomic_coulomb(ii,ia,ispin,ne,nprj,
     >                                          sw1,sw2,eh_atom)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  sw2(ne(1)+ne(2),nprj)
      real*8  eh_atom

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

      return
      end


*     *************************************************
*     *                                               *
*     *              nwpw_compcharge_Watom            *
*     *                                               *
*     *************************************************
*
*
      real*8 function nwpw_compcharge_Watom(ii,ia,ispin,ne,nprj,sw1)
      implicit none
      integer ii,ia
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  dmatrix(nprj,nprj)

#include "bafdecls.fh"
#include "errquit.fh"
#include "nwpw_compcharge.fh"

*     **** local variables ****
      integer iii,iia,n,k,ms,lm,iprj,jprj,shift,indx,indx1,n1(2),n2(2)
      real*8  scal,coeff,w

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega


      iii = int_mb(ion_iontopaw(1) +ii-1)
      iia = int_mb(katm_iontopaw(1)+ia-1)

      indx = (iii-1)*2*lm_size_max
      call dcopy(2*lm_size_max,0.0d0,0,dbl_mb(Qlm(1)+indx),1)
      n1(1) = 1
      n1(2) = ne(1)+1
      n2(1) = ne(1)
      n2(2) = ne(1)+ne(2)
      scal = 1.0d0/lattice_omega()

      shift = int_mb(shift_Tndiff(1)+iia-1)
      do k=1,int_mb(nindx_Tndiff(1)+iia-1)
         lm    = int_mb(lm_Tndiff(1)+shift+k-1)
         iprj  = int_mb(iprj_Tndiff(1)+shift+k-1)
         jprj  = int_mb(jprj_Tndiff(1)+shift+k-1)
         coeff = scal*dbl_mb(coeff_Tndiff(1)+shift+k-1)
         do ms=1,ispin
            w = 0.0d0
            do n=n1(ms),n2(ms)
               w = w + sw1(n,iprj)*sw1(n,jprj)
            end do
c            indx1 = (iii-1)*2*lm_size_max+(ms-1)*lm_size_max + lm
            indx1 = indx + (ms-1)*lm_size_max + lm
            dbl_mb(Qlm(1)+indx1) = dbl_mb(Qlm(1)+indx1) + coeff*w
         end do
      end do
      call D1dB_Vector_SumAll(2*lm_size_max,dbl_mb(Qlm(1)+indx))

      return
      end 

*     ***************************************************
*     *                                                 *
*     *         nwpw_compcharge_density_matrix          *
*     *                                                 *
*     ***************************************************
*
*     This routine computes a spin dependent density matrix using sw1
*
      subroutine nwpw_compcharge_density_matrix(ispin,ne,nprj,sw1,wmatrix)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  wmatrix(nprj,nprj,2)

*     **** local variables ****
      integer i,j,ms,n,n1(2),n2(2)

      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1)+1
      n2(2) = ne(1)+ne(2)

      call dcopy(2*nprj*nprj,0.0d0,0,wmatrix,1)
      do j=1,nprj
         do ms=1,ispin
         do n=n1(ms),n2(ms)
            wmatrix(j,j,ms) = wmatrix(j,j,ms) + sw1(n,j)*sw1(n,j)
         end do
         end do
      end do

      do j=1,nprj
         do i=j+1,nprj
            do ms=1,ispin
            do n=n1(ms),n2(ms)
               tmp = sw1(n,i)*sw1(n,j)
               wmatrix(i,j,ms) = wmatrix(i,j,ms) + tmp
               wmatrix(j,i,ms) = wmatrix(j,i,ms) + tmp
            end do
            end do
         end do
      end do
      return
      end


*    **********************************************************
*    *                                                        *
*    *              nwpw_compcharge_Watom1_sub                *
*    *                                                        *
*    **********************************************************
      real*8 function nwpw_compcharge_Watom1_sub(ispin,nprj,wmatrix,
     >                                            nindx1,
     >                                            iprj_indx1,jprj_indx1,
     >                                            kprj_indx1,lprj_indx1,
     >                                            vh_indx1)
      implicit none
      integer ispin,ne(2),nprj
      real*8  wmatrix(nprj,nprj,2)
      integer nindx1
      integer iprj_indx1(*),jprj_indx1(*)
      integer kprj_indx1(*),lprj_indx1(*)
      real*8  vh_indx1(*)

*     **** local variables ****
      integer k,ms,lm,iprj,jprj,kprj,lprj
      real*8 vh,W

      W = 0.0d0
      do k=1,nindx1
         iprj = iprj_indx1(k)
         jprj = jprj_indx1(k)
         kprj = kprj_indx1(k)
         lprj = lprj_indx1(k)
         vh   = vh_indx1(k)
         do ms=1,ispin
            W = W + 0.5d0*vh*
     >                   *wdensity(iprj,jprj,ms)
     >                   *wdensity(kprj,lprj,ms)
         end do
      end do
      if (ispin.eq.1) W = W + W

      nwpw_compcharge_Watom1_sub = W
      return 
      end

*    **********************************************************
*    *                                                        *
*    *              nwpw_compcharge_Watom2_sub                *
*    *                                                        *
*    **********************************************************
      real*8 function nwpw_compcharge_Watom2_sub(ispin,nprj,wmatrix,
     >                                           lmsize_max,Qlm,
     >                                           nindx2,lm_indx2(*),
     >                                           iprj_indx2,jprj_indx2,
     >                                           vcmp_indx2)
      implicit none
      integer ispin,ne(2),nprj
      real*8  wmatrix(nprj,nprj,2)
      integer nindx2,lm_indx2(*)
      integer iprj_indx2(*),jprj_indx2(*)
      real*8  vcmp_indx2(*)

*     **** local variables ****
      integer k,lm,iprj,jprj,ms
      real*8 W,vg

      W = 0.0d0
      do k=1,nindx2
         lm   = lm_indx2(k)
         iprj = iprj_indx2(k)
         jprj = jprj_indx2(k)
         vcmp = vcmp_indx2(k)
         do ms=1,ispin
            W = W + (Qlm(lm,1)+Qlm(lm,ispin))*vcmp*wdensity(iprj,jprj,ms)
         end do
      end do
      if (ispin.eq.1) W = W + W

      nwpw_compcharge_Watom2_sub = W
      return
      end

*    **********************************************************
*    *                                                        *
*    *              nwpw_compcharge_Watom3_sub                *
*    *                                                        *
*    **********************************************************
      real*8 function nwpw_compcharge_Watom3_sub(ispin,lmsize_max,Qlm,
     >                                           nindx3,lm_indx3,
     >                                           vg_indx3)
      implicit none
      integer ispin,ne(2),nprj
      real*8  wmatrix(nprj,nprj,2)
      integer nindx3,lm_indx3(*)
      real*8  vg_indx3(*)

*     **** local variables ****
      integer k,lm
      real*8 W,vg

      W = 0.0d0
      do k=1,nindx3
         lm = lm_indx3(k)
         vg = vg_indx3(k)
         W  = W + 0.5d0*vg*(Qlm(lm,1)+Qlm(lm,ispin))
      end do
      if (ispin.eq.1) W = W + W
     
      nwpw_compcharge_Watom3_sub = W
      return
      end

      subroutine nwpw_compcharge_Vatom1_sub(ispin,ne,nprj,sw1,wmatrix,
     >                                      nindx1,
     >                                      iprj_indx1,jprj_indx1,
     >                                      kprj_indx1,lprj_indx1,
     >                                      vh_indx1)
      implicit none
      integer ispin,ne(2),nprj
      real*8  sw1(ne(1)+ne(2),nprj)
      real*8  wmatrix(nprj,nprj,2)
      integer nindx1
      integer iprj_indx1(*),jprj_indx1(*)
      integer kprj_indx1(*),lprj_indx1(*)
      real*8  vh_indx1(*)

      do k=1,nindx1
         iprj = iprj_indx1(k)
         jprj = jprj_indx1(k)
         kprj = kprj_indx1(k)
         lprj = lprj_indx1(k)
         vh   = vh_indx1(k)
         wab = 0.0d0
         wcd = 0.0d0
         do n=1,ne(1)+ne(2)
            wcd = wcd + sw1(n,kprj)*sw1(n,lprj)
         end do
         if (ispin.eq.1) wcd = wcd+wcd
         do n=1,ne(1)+ne(2)
            sw2(n,iprj) = vh*wcd*sw1(n,jprj)
         end do
      end do

      return
      end


c $Id$
