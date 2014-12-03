*     **********************************
*     *	                               *
*     *          BLCJ_SR_init          *
*     *                                *
*     **********************************

*   This reads in Gaussian Basis and initializes
*  variables for the SR contribution 
*
      subroutine BLCJ_SR_init(rtdb)
      implicit none
#include "errquit.fh"
      integer rtdb

#include "bafdecls.fh"
#include "btdb.fh"

*     **** BLCJ_SR common block ****
      integer Nbasis
      integer iatom(2)
      integer c(2),g(2),coeff(2)
      common / BLCJ_SR_blck / c,coeff,g,iatom,Nbasis

*     **** local variables ***
      logical value
      integer taskid
      integer MASTER
      parameter(MASTER=0)

      integer i
      real*8 sqrtpi

      call Parallel_taskid(taskid)

      value = btdb_get(rtdb,'water_pseudopotential:num_basis',
     >                 mt_int,1,Nbasis)

      if (Nbasis.gt.0) then

        value = BA_alloc_get(mt_int,Nbasis,'iatom',
     >                       iatom(2),iatom(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'c',c(2),c(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'coeff',coeff(2),coeff(1))
        value = value.and.
     >          BA_alloc_get(mt_dbl,Nbasis,'g',g(2),g(1))
        if (.not. value) 
     >   call errquit('init_BLCJ_SR:out of heap memory',0, MA_ERR)

         value = btdb_get(rtdb,'water_pseudopotential:sr_ion_number',
     >                   mt_int,Nbasis,int_mb(iatom(1)))
         value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:sr_c',
     >                   mt_dbl,Nbasis,dbl_mb(c(1)))
         value = value.and.
     >          btdb_get(rtdb,'water_pseudopotential:sr_g',
     >                   mt_dbl,Nbasis,dbl_mb(g(1)))
         if (.not.value)
     >    call errquit('init_BLCJ_SR:rtdb error',0, RTDB_ERR)



         if (taskid.eq.MASTER) then
           write(*,*)
           write(*,*) "Gaussian Basis for SR Water Pseudopotential"
           write(*,*) "-------------------------------------------"
           write(*,*)
           write(*,'  (A15,2x,A15,2x,A4)') "coeff","exponent","atom"
           do i=1,Nbasis
             if (int_mb(iatom(1)+i-1).eq.1) 
     >         write(*,'  (E15.6,2x,E15.6,2x,A4)') dbl_mb(c(1)+i-1),
     >                                          dbl_mb(g(1)+i-1),"O"
             if (int_mb(iatom(1)+i-1).eq.2) 
     >         write(*,'  (E15.6,2x,E15.6,2x,A4)') dbl_mb(c(1)+i-1),
     >                                            dbl_mb(g(1)+i-1),"H1"
             if (int_mb(iatom(1)+i-1).eq.3) 
     >         write(*,'  (E15.6,2x,E15.6,2x,A4)') dbl_mb(c(1)+i-1), 
     >                                            dbl_mb(g(1)+i-1),"H2"
           end do
           write(*,*)
         end if

         sqrtpi = dsqrt(4.0d0*datan(1.0d0))
         do i=1,Nbasis
            dbl_mb(coeff(1)+i-1) = dbl_mb(c(1)+i-1)
     >                           * (1.0d0/(sqrtpi*dbl_mb(g(1)+i-1)))**3
            dbl_mb(g(1)+i-1) = 1.0d0/dbl_mb(g(1)+i-1)**2
         end do

      end if

      return
      end

*     **********************************
*     *	                               *
*     *          BLCJ_SR_end	       *
*     *                                *
*     **********************************

      subroutine BLCJ_SR_end()

#include "bafdecls.fh"

*     **** BLCJ_SR common block ****
      integer Nbasis
      integer iatom(2)
      integer c(2),g(2),coeff(2)
      common / BLCJ_SR_blck / c,coeff,g,iatom,Nbasis

*     **** local variables ***
      logical value

      if (Nbasis.gt.0) then
         value = BA_free_heap(iatom(2))
         value = value.and.BA_free_heap(c(2))
         value = value.and.BA_free_heap(coeff(2))
         value = value.and.BA_free_heap(g(2))
         if (.not.value) call errquit('error BA_free_heap',0, MA_ERR)
      end if

      return
      end



*     **********************************
*     *	                               *
*     *          BLCJ_SR               *
*     *                                *
*     **********************************

*   This routine calculate the SR contribution 
*  to the BLCJ water psp.
*

      subroutine BLCJ_SR(RO,R1,R2,
     >                   n2ft3d,rgrid,
     >                   Vwpsp)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      integer n2ft3d
      real*8 rgrid(3,n2ft3d)
      real*8 Vwpsp(n2ft3d)

#include "bafdecls.fh"

*     **** BLCJ_SR common block ****
      integer Nbasis
      integer iatom(2)
      integer c(2),g(2),coeff(2)
      common / BLCJ_SR_blck / c,coeff,g,iatom,Nbasis

*     **** local variables ****
      integer i,k
      real*8 d2
      real*8 Rgw(3,3)

      Rgw(1,1) = RO(1)
      Rgw(2,1) = RO(2)
      Rgw(3,1) = RO(3)
      Rgw(1,2) = R1(1)
      Rgw(2,2) = R1(2)
      Rgw(3,2) = R1(3)
      Rgw(1,3) = R2(1)
      Rgw(2,3) = R2(2)
      Rgw(3,3) = R2(3)

      do k=1,n2ft3d
         do i=1,Nbasis
            d2 = (rgrid(1,k)-Rgw(1,int_mb(iatom(1)+i-1)))**2 
     >         + (rgrid(2,k)-Rgw(2,int_mb(iatom(1)+i-1)))**2
     >         + (rgrid(3,k)-Rgw(3,int_mb(iatom(1)+i-1)))**2 

            Vwpsp(k) = Vwpsp(k) 
     >               + dbl_mb(coeff(1)+i-1)*dexp(-d2*dbl_mb(g(1)+i-1))
         end do
      end do

      return 
      end

*     **********************************
*     *	                               *
*     *          BLCJ_SR_Fwater        *
*     *                                *
*     **********************************

*   This routine calculates the water force for
*  the SR contribution to the BLCJ water psp.
*

      subroutine BLCJ_SR_Fwater(RO,R1,R2,
     >                   n2ft3d,rgrid,rho,dv,
     >                   fo,f1,f2)
      implicit none
      real*8 RO(3),R1(3),R2(3)
      integer n2ft3d
      real*8 rgrid(3,*),rho(*)
      real*8 dv
      real*8 fo(3),f1(3),f2(3)

#include "bafdecls.fh"

*     **** BLCJ_SR common block ****
      integer Nbasis
      integer iatom(2)
      integer c(2),g(2),coeff(2)
      common / BLCJ_SR_blck / c,coeff,g,iatom,Nbasis

*     **** local variables ****
      integer i,k
      real*8 x,y,z,d2,desr
      real*8 Rgw(3,3), ff(3,3)

      Rgw(1,1) = RO(1)
      Rgw(2,1) = RO(2)
      Rgw(3,1) = RO(3)
      Rgw(1,2) = R1(1)
      Rgw(2,2) = R1(2)
      Rgw(3,2) = R1(3)
      Rgw(1,3) = R2(1)
      Rgw(2,3) = R2(2)
      Rgw(3,3) = R2(3)

      call dcopy(9,0.0d0,0,ff,1)
      do k=1,n2ft3d
         do i=1,Nbasis
            x  = rgrid(1,k) - Rgw(1,int_mb(iatom(1)+i-1))
            y  = rgrid(2,k) - Rgw(2,int_mb(iatom(1)+i-1))
            z  = rgrid(3,k) - Rgw(3,int_mb(iatom(1)+i-1))
            d2 = (x*x+y*y+z*z) 

            desr = dv*rho(k)
     >           * dbl_mb(coeff(1)+i-1)
     >           * dexp(-d2*dbl_mb(g(1)+i-1))
     >           *  (-2.0d0*dbl_mb(g(1)+i-1))
            ff(1,int_mb(iatom(1)+i-1)) = 
     >      ff(1,int_mb(iatom(1)+i-1)) + (x)*desr

            ff(2,int_mb(iatom(1)+i-1)) = 
     >      ff(2,int_mb(iatom(1)+i-1)) + (y)*desr

            ff(3,int_mb(iatom(1)+i-1)) = 
     .      ff(3,int_mb(iatom(1)+i-1)) + (z)*desr
         end do
      end do
      fo(1) = ff(1,1)
      fo(2) = ff(2,1)
      fo(3) = ff(3,1)
      f1(1) = ff(1,2)
      f1(2) = ff(2,2)
      f1(3) = ff(3,2)
      f2(1) = ff(1,3)
      f2(2) = ff(2,3)
      f2(3) = ff(3,3)
      call D3dB_Vector_SumAll(3,fo)
      call D3dB_Vector_SumAll(3,f1)
      call D3dB_Vector_SumAll(3,f2)
      return 
      end

c $Id$
