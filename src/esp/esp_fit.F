      subroutine esp_fit(g,v,x,q,a,b,c,akeep,qkeep,n,qf)
c
c $Id: esp_fit.F,v 1.4 1997-10-07 23:47:48 d3j191 Exp $
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
c
      real*8 g(3,mxgrid),v(mxgrid)
      real*8 x(3,natoms),q(natoms),qf(ndim,4)
      real*8 a(ndim,ndim),b(ndim),c(ndim)
      real*8 akeep(ndim,ndim),qkeep(natoms)
      integer n(natoms)
c
      integer i,j,k,jcon,lcon
      integer l_it0,l_it1,l_it2,l_dt,ii0ptr,ii1ptr,ii2ptr,idtptr
      real*8 sum,dif,difm
c
      nf=1
c
c     set up matrix of linear coefficients
c
      do 1 i=1,natoms
      do 2 j=i,natoms
      sum=zero
      do 3 k=1,ngrid
      sum=sum+one/dsqrt( ( (x(1,i)-g(1,k))*(x(1,i)-g(1,k))+
     +                       (x(2,i)-g(2,k))*(x(2,i)-g(2,k))+
     +                       (x(3,i)-g(3,k))*(x(3,i)-g(3,k)) )*
     +                     ( (x(1,j)-g(1,k))*(x(1,j)-g(1,k))+
     +                       (x(2,j)-g(2,k))*(x(2,j)-g(2,k))+
     +                       (x(3,j)-g(3,k))*(x(3,j)-g(3,k)) ) )
    3 continue
      a(i,j)=sum
      a(j,i)=sum
    2 continue
      a(i,natoms+1)=one
      do 4 j=2,nconst
      a(i,natoms+j)=zero
    4 continue
    1 continue
      do 5 j=1,natoms
      a(natoms+1,j)=one
      do 6 i=2,nconst
      a(natoms+i,j)=zero
    6 continue
    5 continue
      do 7 i=1,nconst
      do 8 j=1,nconst
      a(natoms+i,natoms+j)=zero
    8 continue
    7 continue
c
c     construct column vector c
c
      do 9 i=1,natoms
      sum=zero
      do 10 k=1,ngrid
      sum=sum+v(k)/dsqrt((x(1,i)-g(1,k))*(x(1,i)-g(1,k))+
     +                   (x(2,i)-g(2,k))*(x(2,i)-g(2,k))+
     +                   (x(3,i)-g(3,k))*(x(3,i)-g(3,k)))
   10 continue
      c(i)=sum
    9 continue
      c(natoms+1)=charge
c
      jcon=0
      do 11 i=2,nconst
      c(natoms+i)=rcon(i-1)
      jcon=jcon+1
      lcon=icon(jcon)
      do 12 j=1,lcon
      jcon=jcon+1
      k=icon(jcon)
      if(k.gt.0) then
      a(natoms+i,k)=one
      a(k,natoms+i)=one
      else
      a(natoms+i,-k)=-one
      a(-k,natoms+i)=-one
      endif
   12 continue
   11 continue
      do 13 i=1,ndim
      do 14 j=1,ndim
      akeep(i,j)=a(i,j)
   14 continue
   13 continue
c
c     allocate work arrays for the matrix inversion
c
      if(.not.ma_push_get(mt_int,2*ndim,'it0',l_it0,ii0ptr))
     + call errquit('esp: ma_push_get it0 failed',me)
      if(.not.ma_push_get(mt_int,ndim,'it1',l_it1,ii1ptr))
     + call errquit('esp: ma_push_get it1 failed',me)
      if(.not.ma_push_get(mt_int,ndim,'it2',l_it2,ii2ptr))
     + call errquit('esp: ma_push_get it2 failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'dt',l_dt,idtptr))
     + call errquit('esp: ma_push_get dt failed',me)
c
c     unconstrained esp fitting
c
      labelf(nf)='      ESP   '
      labelc(nf)='            '
c
c     invert the matrix a
c
      call esp_matinv(a,natoms+1,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)
c
c     matrix multiplication b = A_inv * c
c
      do 15 i=1,natoms+1
      sum=zero
      do 16 j=1,natoms+1
      sum=sum+a(i,j)*c(j)
   16 continue
      qf(i,nf)=sum
   15 continue
      do 17 i=1,ndim
      do 18 j=1,ndim
      a(i,j)=akeep(i,j)
   18 continue
   17 continue
c
c     constrained esp fitting
c
      if(nconst.gt.1) then
      nf=nf+1
      labelf(nf)='      ESP   '
      labelc(nf)='     constr '
c
c     invert the matrix a
c
      call esp_matinv(a,ndim,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)
c
c     matrix multiplication b = A_inv * c
c
      do 19 i=1,ndim
      sum=zero
      do 20 j=1,ndim
      sum=sum+a(i,j)*c(j)
      a(i,j)=akeep(i,j)
   20 continue
      qf(i,nf)=sum
   19 continue
      endif
c
c     restraint fitting
c
      if(irest.gt.0) then
      nf=nf+1
      labelf(nf)='     RESP   '
      labelc(nf)='            '
c
   21 continue
c
      if(irest.eq.1) then
      do 22 i=1,natoms
      if(ihfree.eq.0.or.n(i).ne.1) a(i,i)=akeep(i,i)+resa
   22 continue
      endif
c
      if(irest.eq.2) then
      do 23 i=1,natoms
      q(i)=zero
   23 continue
      endif
c
      niter=0
   24 continue
      niter=niter+1
c
      if(irest.eq.2) then
      do 25 i=1,ndim
      do 26 j=1,ndim
      a(i,j)=akeep(i,j)
   26 continue
   25 continue
      do 27 i=1,natoms
      if(ihfree.eq.0.or.n(i).ne.1)
     + a(i,i)=a(i,i)+resa/dsqrt(q(i)*q(i)+resb*resb)
      qkeep(i)=q(i)
   27 continue
      endif
      if(nf.ne.4) then
c
c     invert the matrix a
c
      call esp_matinv(a,natoms+1,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)
c
c     matrix multiplication b = A_inv * c
c
      do 28 i=1,natoms+1
      sum=zero
      do 29 j=1,natoms+1
      sum=sum+a(i,j)*c(j)
   29 continue
      q(i)=sum
      qf(i,nf)=sum
   28 continue
      else
c
c     invert the matrix a
c
      call esp_matinv(a,ndim,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)
c
c     matrix multiplication b = A_inv * c
c
      do 30 i=1,ndim
      sum=zero
      do 31 j=1,ndim
      sum=sum+a(i,j)*c(j)
   31 continue
      q(i)=sum
      qf(i,nf)=sum
   30 continue
      endif
c
      if(irest.eq.2) then
      difm=zero
      do 32 i=1,natoms
      dif=(q(i)-qkeep(i))*(q(i)-qkeep(i))
      if(difm.lt.dif) difm=dif
   32 continue
      dif=dsqrt(dif)
      if(dif.gt.toler.and.niter.lt.maxit) goto 24
      endif
c
      if(nf.eq.3) then
      nf=nf+1
      labelf(nf)='      RESP  '
      labelc(nf)='     constr '
      do 33 i=1,ndim
      do 34 j=1,ndim
      a(i,j)=akeep(i,j)
   34 continue
   33 continue
      goto 21
      endif
c
      endif
c
c     deallocate work arrays
c
      if(.not.ma_pop_stack(l_dt))
     + call errquit('esp: ma_pop_stack dt failed',me)
      if(.not.ma_pop_stack(l_it2))
     + call errquit('esp: ma_pop_stack it2 failed',me)
      if(.not.ma_pop_stack(l_it1))
     + call errquit('esp: ma_pop_stack it1 failed',me)
      if(.not.ma_pop_stack(l_it0))
     + call errquit('esp: ma_pop_stack it0 failed',me)
c
      return
      end
