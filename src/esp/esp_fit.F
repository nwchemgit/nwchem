      subroutine esp_fit(g,v,x,q,a,b,bkeep,akeep,qkeep,n,qf,l,tag)
c
c $Id$
c
      implicit none
c
#include "esp_params.fh"
#include "esp_common.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "global.fh"
c
      logical esp_matinv
      external esp_matinv
c
      real*8 g(3,mxgrid),v(mxgrid)
      real*8 x(3,natoms),q(natoms),qf(ndim,6)
      real*8 a(ndim,ndim),b(ndim),bkeep(ndim)
      real*8 akeep(ndim,ndim),qkeep(natoms)
      integer n(natoms),l(natoms)
      character*16 tag(natoms)
      logical luncon
c
      integer i,j,k,jcon,lcon,ndima,istage
      integer l_it0,l_it1,l_it2,l_dt,ii0ptr,ii1ptr,ii2ptr,idtptr
      real*8 sum,dif,difm,d,q0
c
      nf=1
      luncon=.true.
      istage=0
      q0=0.0d0
c
c     set up matrix of linear coefficients
c     ------------------------------------
c
      do 1 i=1,ndim
        do 2 j=1,ndim
          a(i,j)=zero
    2   continue
    1 continue
c
      do 3 i=me+1,natoms,np
        do 4 j=i,natoms
          sum=zero
          do 5 k=1,ngrid
            sum=sum+one/dsqrt( ( (x(1,i)-g(1,k))*(x(1,i)-g(1,k))+
     +                           (x(2,i)-g(2,k))*(x(2,i)-g(2,k))+
     +                           (x(3,i)-g(3,k))*(x(3,i)-g(3,k)) )*
     +                         ( (x(1,j)-g(1,k))*(x(1,j)-g(1,k))+
     +                           (x(2,j)-g(2,k))*(x(2,j)-g(2,k))+
     +                           (x(3,j)-g(3,k))*(x(3,j)-g(3,k)) ) )
    5     continue
          a(i,j)=sum
          a(j,i)=sum
    4   continue
        a(i,natoms+1)=one
        do 6 j=2,nconst
          a(i,natoms+j)=zero
    6   continue
    3 continue
c
      if(np.gt.1) call ga_dgop(mre_007,a,ndim*ndim,'+')
c
      do 7 j=1,natoms
        a(natoms+1,j)=one
    7 continue
c
c     construct column vector b
c     -------------------------
c
      do 8 i=1,natoms
        sum=zero
        do 9 k=1,ngrid
          sum=sum+v(k)/dsqrt((x(1,i)-g(1,k))*(x(1,i)-g(1,k))+
     +                       (x(2,i)-g(2,k))*(x(2,i)-g(2,k))+
     +                       (x(3,i)-g(3,k))*(x(3,i)-g(3,k)))
    9   continue
        b(i)=sum
    8 continue
      b(natoms+1)=charge
c
      ndima=natoms+1
      jcon=0
      do 10 i=2,nconst
        b(natoms+i)=rcon(i-1)
        jcon=jcon+1
        lcon=icon(jcon)
        do 11 j=1,lcon
          jcon=jcon+1
          k=icon(jcon)
          if(k.gt.0) then
            a(natoms+i,k)=one
            a(k,natoms+i)=one
          else
            a(natoms+i,-k)=-one
            a(-k,natoms+i)=-one
          endif
   11   continue
        ndima=ndima+1
   10 continue
c
c     save a copy in akeep and bkeep
c     ------------------------------
c
      do 12 i=1,ndim
        do 13 j=1,ndim
          akeep(i,j)=a(i,j)
   13   continue
        bkeep(i)=b(i)
   12 continue
c
c     allocate work arrays for the matrix inversion
c     ---------------------------------------------
c
      if(.not.ma_push_get(mt_int,2*ndim,'it0',l_it0,ii0ptr))
     + call md_abort('esp: ma_push_get it0 failed',me)
      if(.not.ma_push_get(mt_int,ndim,'it1',l_it1,ii1ptr))
     + call md_abort('esp: ma_push_get it1 failed',me)
      if(.not.ma_push_get(mt_int,ndim,'it2',l_it2,ii2ptr))
     + call md_abort('esp: ma_push_get it2 failed',me)
      if(.not.ma_push_get(mt_dbl,ndim,'dt',l_dt,idtptr))
     + call md_abort('esp: ma_push_get dt failed',me)
c
c     zeroth stage
c     ------------
c
c================================
c     unconstrained esp fitting =
c================================
c
      labelf(nf)='      ESP   '
      labelc(nf)='            '
c
c     invert the matrix a
c     -------------------
c
      if(.not.esp_matinv(a,natoms+1,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim))
     + call md_abort('Singular matrix in esp_fit',1)
c
c     matrix multiplication b = A_inv * c
c     -----------------------------------
c
      do 14 i=1,natoms+1
      sum=zero
      do 15 j=1,natoms+1
      sum=sum+a(i,j)*b(j)
   15 continue
      qf(i,nf)=sum
   14 continue
c
c==============================
c     constrained esp fitting =
c==============================
c
c     set matrix a back
c     -----------------
c
      do 16 i=1,ndim
      do 17 j=1,ndim
      a(i,j)=akeep(i,j)
   17 continue
   16 continue
c
      if(ndima.gt.natoms+1) then
      nf=nf+1
      labelf(nf)='      ESP   '
      labelc(nf)='     constr '
c
c     invert the matrix a
c     -------------------
c
      if(.not.esp_matinv(a,ndima,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim))
     + call md_abort('Singular matrix in esp_fit',2)
c
c     matrix multiplication b = A_inv * c
c     -----------------------------------
c
      do 18 i=1,ndim
      sum=zero
      do 19 j=1,ndim
      sum=sum+a(i,j)*b(j)
      a(i,j)=akeep(i,j)
   19 continue
      qf(i,nf)=sum
   18 continue
      endif
c
c========================
c     restraint fitting =
c========================
c
c     first stage
c     -----------
      istage=1
c
      if(irest.gt.0) then
      nf=nf+1
      luncon=.true.
      labelf(nf)='      RESP  '
      labelc(nf)='            '
c
c     initialize charges with zero
c     ----------------------------
c
   20 continue
c
      do 21 i=1,natoms
      qkeep(i)=zero
   21 continue
c
      niter=0
c
   22 continue
c
      niter=niter+1
c
c     copy stored matrix a
c     --------------------
c
      do 23 i=1,ndim
      do 24 j=1,ndim
      a(i,j)=akeep(i,j)
   24 continue
   23 continue
c
c     harmonic restraints
c     -------------------
c
      if(irest.eq.1) then
      do 26 i=1,natoms
      if(ihfree.eq.0.or.n(i).ne.1) a(i,i)=akeep(i,i)+resa
   26 continue
      if(q0.gt.0.0d0) then
      do 126 i=1,natoms
      do 27 j=i+1,natoms
      if(.not.(tag(i)(1:2).eq.'C '.and.tag(j)(1:2).eq.'H ')) goto 126
      a(j,j)=akeep(j,j)+resa*(q(j)-q0)/q(j)
      b(j)=bkeep(j)+q0*resa*(q(j)-q0)/q(j)
c      write(*,'(a,i5,2f12.6)') 'q',j,q(j),q0
   27 continue
  126 continue
      endif
      endif
c
c     hyperbolic restraints
c     ---------------------
c
      if(irest.eq.2) then
      do 28 i=1,natoms
      if(ihfree.eq.0.or.n(i).ne.1)
     + a(i,i)=akeep(i,i)+resa/dsqrt(q(i)*q(i)+resb*resb)
   28 continue
      if(q0.gt.0.0d0) then
      do 128 i=1,natoms
      do 29 j=i+1,natoms
      if(.not.(tag(i)(1:2).eq.'C '.and.tag(j)(1:2).eq.'H ')) goto 128
      a(j,j)=akeep(j,j)+resa/dsqrt((q(j)-q0)**2+resb*resb)
      b(j)=bkeep(j)+q0*resa/dsqrt((q(j)-q0)**2+resb*resb)
c      write(*,'(a,i5,4f12.6)') 'q',j,q(j),q0
   29 continue
  128 continue
      endif
      endif
c
c     invert matrix a and perform matrix multiplication b = A_inv * c
c     ---------------------------------------------------------------
c
      if(luncon) then
c
      if(.not.esp_matinv(a,natoms+1,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)) then
      if(me.eq.0) write(lfnout,1999)
 1999 format(/,' Possible fit problem: singular constraint matrix',/)
      endif
c
      do 30 i=1,natoms+1
      sum=zero
      do 31 j=1,natoms+1
      sum=sum+a(i,j)*b(j)
   31 continue
      q(i)=sum
      qf(i,nf)=sum
   30 continue
c
      else
c
      if(.not.esp_matinv(a,ndima,int_mb(ii0ptr),int_mb(ii1ptr),
     + int_mb(ii2ptr),dbl_mb(idtptr),ndim)) then
      if(me.eq.0) write(lfnout,1999)
      endif
c
      do 32 i=1,ndim
      sum=zero
      do 33 j=1,ndim
      sum=sum+a(i,j)*b(j)
   33 continue
      q(i)=sum
      qf(i,nf)=sum
   32 continue
      endif
c
c     check convergence
c     -----------------
c
      if(irest.eq.2.or.q0.gt.0.0d0) then
      difm=zero
      do 34 i=1,natoms
      dif=(q(i)-qkeep(i))*(q(i)-qkeep(i))
      if(difm.lt.dif) difm=dif
   34 continue
      difm=dsqrt(difm)
c      write(*,'(3i4,5f12.6)') nf,niter,maxit,difm,toler,q(1),resa,resb
      do 25 i=1,natoms
      qkeep(i)=q(i)
   25 continue
      if(difm.gt.toler.and.niter.lt.maxit) goto 22
      endif
c
      if(istage.eq.1) then
      istage=2
      if(ndima.gt.natoms+1.and.luncon) then
      nf=nf+1
      luncon=.false.
      labelf(nf)='      RESP  '
      labelc(nf)='     constr '
      goto 20
      endif
      endif
c
c     second and third stage RESP
c
      if(istage.eq.2.or.(istage.eq.3.and.qh.gt.0.0d0)) then
      nf=nf+1
      luncon=.false.
      if(istage.eq.2) then
      labelf(nf)='      RESP2 '
      else
      labelf(nf)='      RESP3 '
      q0=qh
      endif
      labelc(nf)='     constr '
      istage=istage+1
      resa=resc
      resb=resd
c
      do 37 i=1,natoms
      l(i)=1
   37 continue
      do 38 i=1,natoms-1
      do 39 j=i+1,natoms
      if((tag(i)(1:2).eq.'C '.and.tag(j)(1:2).eq.'H ').or.
     + (tag(i)(1:2).eq.'H '.and.tag(j)(1:2).eq.'C ')) then
      d=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2+(x(3,i)-x(3,j))**2)
      if(cau2nm*d.lt.0.12) then
      l(i)=0
      l(j)=0
      endif
      endif
   39 continue
   38 continue
      do 40 i=natoms+2,ndim
      do 41 j=1,ndim
      akeep(i,j)=zero
      akeep(j,i)=zero
   41 continue
   40 continue
      ndima=natoms+1
      jcon=0
      do 42 i=2,nconst
      jcon=jcon+1
      lcon=icon(jcon)
      if(lcon.lt.0) lcon=-lcon
      do 143 j=1,lcon
      if(l(iabs(icon(jcon+j))).eq.0) goto 144
 143  continue
      jcon=jcon+lcon
      goto 42
 144  continue
      ndima=ndima+1
      do 43 j=1,lcon
      k=icon(jcon+j)
      if(k.gt.0) then
      akeep(ndima,k)=one
      akeep(k,ndima)=one
      else
      akeep(ndima,-k)=-one
      akeep(-k,ndima)=-one
      endif
   43 continue
      if(me.eq.0)
     + write(lfnout,1043) rcon(i-1),
     + (icon(jcon+j),tag(iabs(icon(jcon+j)))(1:2),j=1,lcon)
 1043 format(' Stage 2 constraint charge ',f12.6,' on atoms',
     + 8(i5,a2),/,(47x,8(i5,a2)))
      jcon=jcon+lcon
      b(ndima)=rcon(i-1)
   42 continue
      do 45 i=1,natoms
      if(tag(i)(1:2).eq.'C ') then
      k=0
      do 46 j=1,natoms
      if(tag(j)(1:2).eq.'H ') then
      d=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2+(x(3,i)-x(3,j))**2)
      if(cau2nm*d.lt.0.12) then
      if(k.gt.0) then
      ndima=ndima+1
      akeep(ndima,k)=one
      akeep(k,ndima)=one
      akeep(ndima,j)=-one
      akeep(j,ndima)=-one
      b(ndima)=zero
      endif
      k=j
      endif
      endif
   46 continue
      endif
   45 continue
      do 44 i=1,natoms
      if(l(i).eq.1) then
      ndima=ndima+1
      akeep(ndima,i)=one
      akeep(i,ndima)=one
      b(ndima)=q(i)
      endif
   44 continue
      goto 20
      endif
c
      endif
 9999 continue
c
c     deallocate work arrays
c
      if(.not.ma_pop_stack(l_dt))
     + call md_abort('esp: ma_pop_stack dt failed',me)
      if(.not.ma_pop_stack(l_it2))
     + call md_abort('esp: ma_pop_stack it2 failed',me)
      if(.not.ma_pop_stack(l_it1))
     + call md_abort('esp: ma_pop_stack it1 failed',me)
      if(.not.ma_pop_stack(l_it0))
     + call md_abort('esp: ma_pop_stack it0 failed',me)
c
      return
      end
