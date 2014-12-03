c
c $Id$
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

#include "bafdecls.fh"
#include "btdb.fh"
#include "truncating_function.fh"
#include "errquit.fh"

*     **** local variables ****
      integer MASTER,taskid
      parameter (MASTER=0)
c      integer nrho
c      integer rho(2),f(2)

      logical value,found
      integer npack0,nfft3d,G(3),nfft(3),n2ft3d
      integer i,ii,ia,l,j,k,zero,pzero,qzero,nxh,nyh
      real*8  gg,rr,rcut_in,pp_in
      real*8  unita(3,3),r,f,x,y,z
      integer tmp1(2),tmp2(2),rgrid(2)
      character*4  element
      character*50 fname
      character*255 full_filename
      integer nx,ny,nz
      real*8  scal1




*     **** external functions ****
      logical      nwpw_filefind
      character*2 ion_atom
      integer  G_indx,ion_nkatm,ion_atn,control_ngrid
      real*8   Truncating_Function_transform,control_unita
      real*8   Truncating_Function_transform0
      external nwpw_filefind
      external ion_atom
      external G_indx,ion_nkatm,ion_atn,control_ngrid
      external Truncating_Function_transform,control_unita
      external Truncating_Function_transform0


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

      if (.not.btdb_get(rtdb,'nwpw:truncating_rcut_scaling',
     >                  mt_dbl,1,rcut))
     >   rcut = 1.2d0

      if (.not.btdb_get(rtdb,'nwpw:truncating_power',mt_dbl,1,pp))
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

      nt = ion_nkatm()

      value = BA_alloc_get(mt_dbl,(nt),'Trnc',Trnc(2),Trnc(1))
      if (.not. value) 
     >  call errquit('Truncating_Function_init: out of heap memory',0,
     &       MA_ERR)


*     ***** Generate Fourier transform of Truncating Functions *****
      do ia=1,nt
        ii = ion_atn(ia)
        rr = (vander(ii)/0.529177d0)*rcut
        if (taskid.eq.MASTER) then
          write(*,1163) ion_atom(ia),rr,(vander(ii)/0.529177d0)
        end if
 1163   FORMAT(5X,A2,' radii = ',E11.3,' (', E11.3,' van der Waal)')

        dbl_mb(Trnc(1)+ia-1) = rr

      end do

      return
      end

*     **********************************
*     *                                *
*     *   Truncating_Function_end      *
*     *                                *
*     **********************************

      subroutine Truncating_Function_end()
      implicit none
#include "errquit.fh"

#include "bafdecls.fh"
#include "truncating_function.fh"

      if (.not. BA_free_heap(Trnc(2))) 
     >  call errquit('Truncating_Function_end: freeing heap',0, MA_ERR)
      return
      end

*     ****************************
*     *                          *
*     *   Truncating_Function    *
*     *                          *
*     ****************************
      subroutine Truncating_Function(rgrid,T)
      implicit none
#include "errquit.fh"
      real*8 rgrid(3,*)
      real*8 T(*)

#include "bafdecls.fh"
#include "truncating_function.fh"

*     ***** local variables ****
      logical value
      integer ii
      integer nion,n2ft3d
      integer tmp(2)
     

*     **** external functions ****
      integer  ion_nion
      external ion_nion


      call D3dB_n2ft3d(1,n2ft3d)
      value = BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1))
      if (.not. value)
     >  call errquit(' Switching_Charge_Generate_V:pushing stack',0,
     &       MA_ERR)

      call Truncating_Function_atom(n2ft3d,1,rgrid,T)
      nion = ion_nion()
      do ii=2,nion
        call Truncating_Function_atom(n2ft3d,ii,rgrid,dbl_mb(tmp(1)))
        call Truncating_Function_Union(n2ft3d,T,dbl_mb(tmp(1)),T)
      end do

      value = BA_pop_stack(tmp(2))
      if (.not. value) 
     >  call errquit('Truncating_Function: popping stack',1, MA_ERR)
      return
      end

*     **********************************
*     *                                *
*     *   Truncating_Function_index    *
*     *                                *
*     **********************************

      subroutine Truncating_Function_index(nindx,indx,rgrid,T)
      implicit none
      integer nindx,indx(*)
      real*8 rgrid(3,*)
      real*8 T(*)

#include "bafdecls.fh"
#include "truncating_function.fh"
#include "errquit.fh"

*     ***** local variables ****
      integer ii
      integer n2ft3d
      integer tmp(2)

      call D3dB_n2ft3d(1,n2ft3d)
c      if (.not.BA_push_get(mt_dbl,n2ft3d,'tmp',tmp(2),tmp(1)))
c     >  call errquit('Truncating_Function_index:pushing stack',0,
c     >       MA_ERR)

      call dcopy(n2ft3d,0.0d0,0,T,1)
      do ii=1,nindx
        call Truncating_Function_atom(n2ft3d,indx(ii),rgrid,T)
      end do

c      if (.not.BA_pop_stack(tmp(2)))
c     >  call errquit('Truncating_Function_index:pop stack',1,MA_ERR)

      return
      end



      subroutine Truncating_Function_Union(n2ft3d,T1,T2,T12)
      implicit none
      integer n2ft3d
      real*8 T1(*),T2(*),T12(*)

      integer i

      do i=1,n2ft3d
        T12(i) = 0.5d0*(T1(i) + T2(i) + dabs(T1(i)-T2(i)))
      end do

      return
      end

      subroutine Truncating_Function_atom(n2ft3d,ii,rgrid,T)
      implicit none
      integer n2ft3d,ii
      real*8 rgrid(3,*)
      real*8 T(*)

#include "bafdecls.fh"
#include "truncating_function.fh"

*     **** local variables ****
      integer ia,i,j,n1,n2,n3
      real*8 xii,yii,zii,r,rr,pp2
      real*8 x,y,z,x0,y0,z0,f,a(3,3)

*     **** external functions ****
      integer  ion_katm
      real*8   ion_rion,lattice_unita
      external ion_katm
      external ion_rion,lattice_unita

      do j=1,3
         do i=1,3
            a(i,j) = lattice_unita(i,j)
         end do
      end do

      ia = ion_katm(ii)
      pp2 = pp+2.0d0
      rr = dbl_mb(Trnc(1)+ia-1)


      xii = ion_rion(1,ii)
      yii = ion_rion(2,ii)
      zii = ion_rion(3,ii)
      do i=1,n2ft3d 

         x0 = (xii-rgrid(1,i))
         y0 = (yii-rgrid(2,i))
         z0 = (zii-rgrid(3,i))

*        **** take into account images ****
         do n3=-1,1
         do n2=-1,1
         do n1=-1,1
            x = x0 + n1*a(1,1) + n2*a(1,2) + n3*a(1,3)
            y = y0 + n1*a(2,1) + n2*a(2,2) + n3*a(2,3)
            z = z0 + n1*a(3,1) + n2*a(3,2) + n3*a(3,3)
            r = dsqrt(x*x + y*y + z*z)
            if (r.lt.rr) T(i) = 1.0d0
         end do
         end do
         end do

      end do

      return
      end



*     ********************************
*     *   	                         *
*     *   grad_Truncating_Function   *
*     *                              *
*     ********************************
      subroutine grad_Truncating_Function(Amm,fion)
      implicit none
      complex*16 Amm(*)
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
      parameter (nrho=5001)

      integer i,no,ne
      real*8 q,fourpi,r,sum,drho,pp2,f

      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho-1)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

      q  = dsqrt(gg)
c      sum = 0.0d0
c      do i=2,nrho-1
c        r = (i-1)*drho
c        sum = sum 
c     >    + r*sin(q*r)*(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c      end do
c      r = drho*(nrho-1)
c      sum = sum + 0.5d0*r*sin(q*r)
c     >     *(1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)


*     **** simpson integration ****
c      no = nrho/2
c      ne = no+1
c
c      r = drho*(nrho-1)
c      f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c      sum = -0.5d0*r*r*f*sin(q*r)
c      do i=1,no
c         r = 2*(i-1)*drho
c         f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c         sum = sum + 2.0d0*r*r*f*sin(q*r)
c      end do
c      do i=1,ne
c         r = (2*(i-1)+1)*drho
c         f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c         sum = sum + 4.0d0*r*r*f*sin(q*r)
c      end do
c      sum = sum/3.0d0
c      Truncating_Function_transform = (fourpi/q)*sum*drho

      Truncating_Function_transform = 
     >  (fourpi/q**2)*(sin(q*rcut)/q - rcut*cos(q*rcut))
      return
      end




*     *************************************
*     *                                   *
*     *   Truncating_Function_transform0  *
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
*                        /infty
*                       |
*      T(g=0) = (4*pi)* | r**2 * T(r)  dr
*                       |
*                       / 0

      real*8 function Truncating_Function_transform0(rcut,pp)
      implicit none
      real*8 rcut,pp

*     **** local variables ****
      integer nrho,no,ne
      parameter (nrho=5001)

      integer i
      real*8 fourpi,r,sum,drho,pp2,f


      pp2    = pp + 2.0d0
      drho   = 2.0d0*rcut/dble(nrho-1)
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))

c      sum = 0.0d0
c      do i=2,nrho-1
c            r = (i-1)*drho
c            f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c            sum = sum + r*r*f
c      end do

*     **** simpson integration ****
c     no = nrho/2
c     ne = no+1
c
c     r = drho*(nrho-1)
c     f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c     sum = -0.5d0*r*r*f
c     do i=1,no
c        r = 2*(i-1)*drho
c        f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c        sum = sum + 2.0d0*r*r*f
c     end do
c     do i=1,ne
c        r = (2*(i-1)+1)*drho
c        f = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)
c        sum = sum + 4.0d0*r*r*f
c     end do
c     sum = sum/3.0d0
c      Truncating_Function_transform0 = (fourpi)*sum*drho
       Truncating_Function_transform0 = (fourpi/3.0d0)*(rcut**3)

      return
      end




