c
c $Id: truncating_function.F,v 1.1 2002-03-03 17:42:41 bylaska Exp $
c

***********************************************************
*                                                         *
*            Truncating Function module                   *
*                                                         *
*          Interfaced to nwchem-PSPW code                 *
*                                                         *
*   -- developed by Eric J. Bylaska on February 27, 2002  *
*                                                         *
***********************************************************
*
*
*

*     **********************************
*     *                                *
*     *   Trunctating_Function_init    *
*     *                                *
*     **********************************

      subroutine Truncating_Function_init(rtdb)
      implicit none
      integer rtdb

#include "mafdecls.fh"
#include "rtdb.fh"
#include "truncating_function.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)

      integer npack0,nfft3d,G(3)
      integer i,ii,ia
      real*8  gg,rr
      logical value
      integer tmp1(2)

*     **** external functions ****
      integer  G_indx,ion_nkatm,ion_atn
      real*8   Truncating_Function_transform
      external G_indx,ion_nkatm,ion_atn
      external Truncating_Function_transform

*
*     ----- van der waals radii in angstroms -Needs to be updated   -----
*
      real*8 vander(106)
      data (vander(i),i=1,106)
     1 /1.20d0,1.20d0,
     2 1.37d0,1.45d0,1.45d0,1.50d0,1.50d0,1.40d0,1.35d0,1.30d0,
     3 1.57d0,1.36d0,1.24d0,1.17d0,1.80d0,1.75d0,1.70d0,1.80d0,
     4  18*1.80d0,
     5  18*1.80d0,
     6  32*1.80d0,
     7  20*1.80d0/



      call Parallel_taskid(taskid)

      if (.not.rtdb_get(rtdb,'nwpw:truncating_rcut_scaling',
     >                  mt_dbl,1,rcut))
     >   rcut = 1.2d0

      if (.not.rtdb_get(rtdb,'nwpw:truncating_power',mt_dbl,1,pp))
     >   pp = 6.0d0

      if (taskid.eq.MASTER) then
          write(*,1160)
          write(*,1161) rcut
          write(*,1162) pp
      end if
 1160 FORMAT(/' Truncating Function definition:')
 1161 FORMAT(5X, '- Truncating rcut scaling  = (', E11.3,' )  ')
 1162 FORMAT(5X, '- Truncating power         = (', E11.3,' )  ')


*     **** define Truncating Functions ****      
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

      nt = ion_nkatm()

      value = MA_alloc_get(mt_dbl,(nt*npack0),'Trnc',Trnc(2),Trnc(1))
      if (.not. value) 
     >  call errquit('Truncating_Function_init: out of heap memory',0)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) 
     >   call errquit('Truncating_Function_init: pushing stack',1)

*     ***** Generate Fourier transform of Truncating Functions *****
      do ia=1,nt
        ii = ion_atn(ia)
        rr = (vander(ii)/0.529177d0)*rcut

        do i = 1,nfft3d
          gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >          + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >          + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )
          dbl_mb(tmp1(1)+i-1) = 
     >       Truncating_Function_transform(gg,rr,pp)
        end do
        call Pack_t_pack(0,dbl_mb(tmp1(1)))
        call Pack_t_Copy(0,dbl_mb(tmp1(1)),
     >                     dbl_mb(Trnc(1)+(ia-1)*npack0))
      end do

      value = MA_pop_stack(tmp1(2))
      if (.not. value) 
     >   call errquit('Truncating_Function_init: popping stack',3)

      return
      end

*     **********************************
*     *                                *
*     *   Truncating_Function_end      *
*     *                                *
*     **********************************

      subroutine Truncating_Function_end()
      implicit none

#include "mafdecls.fh"
#include "truncating_function.fh"

      if (.not. MA_free_heap(Trnc(2))) 
     >  call errquit('Truncating_Function_end: freeing heap',0)
      return
      end

*     ****************************
*     *                          *
*     *   Truncating_Function    *
*     *                          *
*     ****************************
      subroutine Truncating_Function(T)
      implicit none
      complex*16 T(*)

#include "mafdecls.fh"
#include "truncating_function.fh"

*     ***** local variables ****
      logical value
      integer ii,ia
      integer nion,nfft3d,npack0
      integer exi(2)

*     **** external functions ****
      integer  ion_nion,ion_katm
      external ion_nion,ion_katm

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = MA_push_get(mt_dcpl,nfft3d,'exi', exi(2), exi(1))
      if (.not. value) 
     >  call errquit('Truncating_Function: pushing stack',0)

      call Pack_c_Zero(0,T)

      nion = ion_nion()
      do ii=1,nion
        ia = ion_katm(ii)

*       **** structure factor ****
        call strfac(ii,dcpl_mb(exi(1)))
        call Pack_c_pack(0,dcpl_mb(exi(1)))

*       **** add to Total Truncating Function ****
        call Pack_tc_Mul(0,dbl_mb(Trnc(1)+npack0*(ia-1)),
     >                   dcpl_mb(exi(1)),
     >                   dcpl_mb(exi(1)))
        call Pack_cc_Sum(0,T,dcpl_mb(exi(1)),T)

      end do

      value = MA_pop_stack(exi(2))
      if (.not. value) 
     >  call errquit('Truncating_Function: popping stack',1)
      return
      end

*     ********************************
*     *   	                         *
*     *   grad_Truncating_Function   *
*     *                              *
*     ********************************
      subroutine grad_Truncating_Function(dng,V_mm,fion)
      implicit none
      complex*16 dng(*),V_mm(*)
      real*8 fion(3,*)

      return
      end 


*     *************************************
*     *                                   *
*     *   Truncating_Function_transform   *
*     *                                   *
*     *************************************

*     This function returns the fourier transform of
*
*            T(r) = (1.0d0-(1.0d0-dexp(-(r/rcut)**(pp+2)))**pp)
*
*      May want to use:
*            T(r) = exp(-(r/rcut)**pp)
*
*     Entry - gg: g squared
*             rcut:
*             pp:
*
*     Exit - returns
*                      /infty
*                     |
*      T(g) = (4*pi)* | r**2 * T(r)* j0(gr) dr
*                     |
*                     / 0

      real*8 function Truncating_Function_transform(gg,rcut,pp)
      implicit none
      real*8 gg,rcut,pp

*     **** local variables ****
      integer nrho
      parameter (nrho=5000)

      integer i
      real*8 q,fourpi,r,sum,drho,pp2,transform

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

*     **** gg > 0 ****
      if (gg.gt.1.0d-9) then
         q  = dsqrt(gg)
         sum = 0.0d0
         do i=2,nrho-1
            r = (i-1)*drho
            sum = sum + r*sin(q*r)
     >           *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         end do
         r = drho*(nrho-1)
         sum = sum + r*sin(q*r)
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         transform = (fourpi/q)*sum*drho

*     **** gg == 0 *****
      else
         sum = 0.0d0
         do i=2,nrho-1
            r = (i-1)*drho
            sum = sum + r*r*
     >           *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         end do
         r = drho*(nrho-1)
         sum = sum + r*r*
     >        *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
         transform = (fourpi)*sum*drho
      end if

      Truncating_Function_transform = transform
      return
      end




