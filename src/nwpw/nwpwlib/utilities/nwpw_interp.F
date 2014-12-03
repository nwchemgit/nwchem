*
*     $Id$
*

*     *********************************
*     *                               *
*     *      nwpw_interp_init         *
*     *                               *
*     *********************************
      subroutine nwpw_interp_init(ndim0,norder0)
      implicit none
      integer ndim0,norder0

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variable ****
      logical value
      integer i,k
      real*8  tmp

      ndim_max   = ndim0
      norder_max = norder0
      value = BA_alloc_get(mt_dbl,norder_max,'root',root(2),root(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,norder_max,'coef',coeff(2),coeff(1))
      value = value.and.
     >        BA_alloc_get(mt_dbl,norder_max*ndim_max*2,
     >                     'pbasis',pbasis(2),pbasis(1))
      if (.not.value) 
     >   call errquit('nwpw_interp_init:out of heap',0,MA_ERR)

      do i=0,norder_max-1
         dbl_mb(root(1)+i) = -1.0d0 + i*2.0d0/dble(norder_max-1)
      end do

      do i=0,norder_max-1
         tmp = 1.0d0
         do k=i+1,(norder_max+i-1)
            tmp = tmp
     >           *(dbl_mb(root(1)+i)-dbl_mb(root(1)+mod(k,norder_max)))
         end do
         dbl_mb(coeff(1)+i) = 1.0d0/tmp
      end do

      return
      end

*     *********************************
*     *                               *
*     *      nwpw_interp_end          *
*     *                               *
*     *********************************
      subroutine nwpw_interp_end()
      implicit none

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variables ****
      logical value

      value =           BA_free_heap(pbasis(2))
      value = value.and.BA_free_heap(coeff(2))
      value = value.and.BA_free_heap(root(2))
      if (.not.value)
     >   call errquit('nwpw_interp_end:freeing heap',0,MA_ERR)
      return
      end


*     *********************************
*     *                               *
*     *      nwpw_interp_basis        *
*     *                               *
*     *********************************
      real*8 function nwpw_interp_basis(i,x)
      implicit none
      integer i
      real*8 x

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variables ****
      integer k
      real*8  f

      f = 1.0d0
      do k=i+1,(i+norder_max-1)
         f = f*(x-dbl_mb(root(1)+mod(k,norder_max)))
      end do
      f = f*dbl_mb(coeff(1)+i)

      nwpw_interp_basis = f
      return
      end


*     *********************************
*     *                               *
*     *      nwpw_interp_dbasis       *
*     *                               *
*     *********************************
      real*8 function nwpw_interp_dbasis(i,x)
      implicit none
      integer i
      real*8 x

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variables ****
      integer k,kk
      real*8  tmp,df

      df = 0.0d0
      do k=i+1,(i+norder_max-1)
         tmp = 1.0d0
         do kk=i+1,(i+norder_max-1)
            if (kk.ne.k) tmp=tmp*(x-dbl_mb(root(1)+mod(kk,norder_max)))
         end do
         df = df + tmp
      end do
      df = df*dbl_mb(coeff(1)+i)

      nwpw_interp_dbasis = df
      return
      end


*     *********************************
*     *                               *
*     *          nwpw_interp          *
*     *                               *
*     *********************************
      real*8 function nwpw_interp(ndim,N,a,b,p,mesh,x)
      implicit none
      integer ndim,N(*),p(*)
      real*8  a(*),b(*),mesh(*),x(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variables ****
      logical failed
      integer d,n1,nn1,i,ii,im,ip,ishift,jm(10),j(10)
      real*8  xm,xp,xtilde,f,fi

*     **** external functions ****
      real*8   nwpw_interp_basis
      external nwpw_interp_basis

      failed = .false.
      do d=1,ndim
         failed = failed.or.(x(d).lt.a(d))
         failed = failed.or.(x(d).gt.b(d))
      end do
      if (failed) then
         nwpw_interp = 0.0d0
         return
      end if


      ishift = norder_max/2 - 1 + mod(norder_max,2)

      do d=1,ndim
         im = (N(d)-1)*(x(d)-a(d))/(b(d)-a(d)) - ishift
         ip = im+norder_max-1
         if (p(d).eq.0) then
            if (im<0) then
               im=0 
               ip=im+norder_max-1
            end if
            if (ip>(N(d)-1)) then
               ip=(N(d)-1)
               im=ip-norder_max+1
            end if
            xm = a(d) + im*(b(d)-a(d))/dble(N(d)-1)
            xp = a(d) + ip*(b(d)-a(d))/dble(N(d)-1)
            xtilde = (2.0d0*x(d) - xp - xm)/(xp-xm)
         else
            xm = a(d) + im*(b(d)-a(d))/dble(N(d))
            xp = a(d) + ip*(b(d)-a(d))/dble(N(d))
            xtilde = (2.0d0*x(d) - xp - xm)/(xp-xm)
            do while (xtilde.ge.1.0d0)
               xtilde = xtilde-2.0d0
            end do
            do while (xtilde.lt.-1.0d0)
               xtilde = xtilde+2.0d0
            end do
         end if
         jm(d) = im
         do i=0,norder_max-1
            dbl_mb(pbasis(1)+i+(d-1)*norder_max) 
     >         = nwpw_interp_basis(i,xtilde)
         end do
      end do

      do d=1,ndim
         j(d) = 0
      end do
      n1=1
      do d=1,ndim
         n1 = n1*norder_max
      end do

      f = 0.0d0
      do i=1,n1
         nn1 = N(1)
         ii = mod(jm(1)+j(1)+N(1),N(1))
         do d=2,ndim
            ii = ii + mod(jm(d)+j(d)+N(d),N(d))*nn1
            nn1 = nn1*N(d)
         end do
         fi=1.0d0
         do d=1,ndim
            fi = fi*dbl_mb(pbasis(1)+j(d)+(d-1)*norder_max)
         end do
         f = f + mesh(ii+1)*fi

         j(1) = j(1) + 1
         do d=1,ndim-1
            if (j(d).ge.norder_max) then
               j(d) = 0
               j(d+1) = j(d+1)+1
            end if
         end do
      end do

      nwpw_interp = f
      return
      end


*     *********************************
*     *                               *
*     *          nwpw_dinterp         *
*     *                               *
*     *********************************
      subroutine nwpw_dinterp(ndim,N,a,b,p,mesh,x,df)
      implicit none
      integer ndim,N(*),p(*)
      real*8  a(*),b(*),mesh(*),x(*),df(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** nwpw_interp_common block ****
      integer ndim_max,norder_max
      integer root(2),coeff(2),pbasis(2)
      common /nwpw_interp_common/ root,coeff,pbasis,ndim_max,norder_max

*     **** local variables ****
      logical failed
      integer d,dd,n1,nn1,i,ii,im,ip,ishift,k,jm(10),j(10)
      real*8  xm,xp,xtilde,dx,fi

*     **** external functions ****
      real*8   nwpw_interp_basis,nwpw_interp_dbasis
      external nwpw_interp_basis,nwpw_interp_dbasis

      failed = .false.
      do d=1,ndim
         failed = failed.or.(x(d).lt.a(d))
         failed = failed.or.(x(d).gt.b(d))
      end do
      if (failed) then
         return
      end if

      ishift = norder_max/2 - 1 + mod(norder_max,2)

      do d=1,ndim
         im = (N(d)-1)*(x(d)-a(d))/(b(d)-a(d)) - ishift
         ip = im+norder_max-1
         if (p(d).eq.0) then
            if (im<0) then
               im=0
               ip=im+norder_max-1
            end if
            if (ip>(N(d)-1)) then
               ip=(N(d)-1)
               im=ip-norder_max+1
            end if
            xm = a(d) + im*(b(d)-a(d))/dble(N(d)-1)
            xp = a(d) + ip*(b(d)-a(d))/dble(N(d)-1)
            xtilde = (2.0d0*x(d) - xp - xm)/(xp-xm)
         else
            xm = a(d) + im*(b(d)-a(d))/dble(N(d))
            xp = a(d) + ip*(b(d)-a(d))/dble(N(d))
            xtilde = (2.0d0*x(d) - xp - xm)/(xp-xm)
            do while (xtilde.ge.1.0d0)
               xtilde = xtilde-2.0d0
            end do
            do while (xtilde.lt.-1.0d0)
               xtilde = xtilde+2.0d0
            end do
         end if
         dx = 2.0d0/(xp-xm)
         jm(d) = im
         do i=0,norder_max-1
            dbl_mb(pbasis(1)+2*i+  (d-1)*2*norder_max)
     >         = nwpw_interp_basis(i,xtilde)
            dbl_mb(pbasis(1)+2*i+1+(d-1)*2*norder_max)
     >         = dx*nwpw_interp_dbasis(i,xtilde)
         end do
      end do

      do d=1,ndim
         j(d) = 0
      end do
      n1=1
      do d=1,ndim
         n1 = n1*norder_max
      end do
      do d=1,ndim
         df(d) = 0.0d0
      end do

      do i=1,n1
         nn1 = N(1)
         ii = mod(jm(1)+j(1)+N(1),N(1))
         do d=2,ndim
            ii = ii + mod(jm(d)+j(d)+N(d),N(d))*nn1
            nn1 = nn1*N(d)
         end do
         do d=0,ndim-1
            fi = dbl_mb(pbasis(1)+2*j(d+1)+1+d*2*norder_max)
            do dd=d+1,(d+ndim-1)
               k = mod(dd,ndim)
               fi = fi*dbl_mb(pbasis(1)+2*j(k+1)+k*2*norder_max)
            end do
            df(d+1) = df(d+1) + mesh(ii+1)*fi
         end do


         j(1) = j(1) + 1
         do d=1,ndim-1
            if (j(d).ge.norder_max) then
               j(d) = 0
               j(d+1) = j(d+1)+1
            end if
         end do
      end do


      return
      end
