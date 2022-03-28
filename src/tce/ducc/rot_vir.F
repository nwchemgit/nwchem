      subroutine rot_vir(rtdb,q_orb,q_ele_a,q_ele_b,nacto,nrotv,
     &                   nrot,xrot)
      IMPLICIT NONE
#include "mafdecls.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
      integer rtdb             ! Run-time database
      integer q_orb,q_ele_a,q_ele_b
      integer nacto,nrotv
      integer nrot
      integer info
      integer i,j,k,l
      integer frozc
      double precision xrot(nrot,nrot) !XR mat. (parameters from input)
      double precision xrot2(nrot,nrot) !XR square 
      double precision full_r(nrot,nrot)
      double precision vr(nrot,nrot),vl(nrot,nrot)
      double precision w(nrot,nrot),wt(nrot,nrot)
      double precision eigenre(nrot),eigenim(nrot)
      double precision work(4*nrot)
      double precision blockm(nbf,nrot),blockm_r(nbf,nrot)
      double precision gmov(nbf,nmo(1))
      double precision del_k,xxx
      logical nodezero
      nodezero=(ga_nodeid().eq.0)
c
c
c   valid for RHF reference
c
c
      frozc=nocc(1)-q_ele_a
c
      do j=1,nmo(1)
      do i=1,nbf
       gmov(i,j) = 0.0d0
      enddo
      enddo
c zeroing blocks
      do i=1,nbf
      do j=1,nrot
       blockm(i,j)=gmov(i,j+frozc)
       blockm_r(i,j)=0.0d0
      enddo
      enddo
c
c ---- Xrot ---
      if(nodezero) then 
       write(6,*)'--- Xrot matrix ---'
       do i=1,nrot
       do j=1,nrot
        write(6,22) i,j,xrot(i,j)
       enddo
       enddo
       write(6,*)'------------------'
       call util_flush(6)
      endif
 22   format(2i5,2x,f12.6)
c --------------
c
c
       call ga_get(g_movecs(1),1,nbf,1,nmo(1),
     1         gmov,nbf)
c
c carving a block
      do i=1,nbf
      do j=1,nrot
       blockm(i,j)=gmov(i,j+frozc+q_ele_a-nacto)
      enddo
      enddo
c
c
      do i=1,nrot
       eigenre(i) = 0.0d0
       eigenim(i) = 0.0d0
      do j=1,nrot
       xrot2(i,j)  = 0.0d0
       full_r(i,j) = 0.0d0
       vr(i,j)     = 0.0d0
       vl(i,j)     = 0.0d0
       w(i,j)      = 0.0d0
       wt(i,j)     = 0.0d0
      enddo
      enddo
c forming xrot2
      do i=1,nrot
      do j=1,nrot
      do k=1,nrot
       xrot2(i,j)=xrot2(i,j)+xrot(i,k)*xrot(k,j)
      enddo
      enddo 
      enddo 
c
      call util_dgeev('v','v',nrot,xrot2,nrot,
     1  eigenre,eigenim,vl,nrot,vr,nrot,
     2  work,4*nrot,info)
c
      if (info.ne.0) call errquit('rot_vir: dgeev failed',info,
     1  unknown_err)
c fixing small values
      do i=1,nrot
       if(dabs(eigenre(i)).lt.1.0d-10) then
        xxx=dabs(eigenre(i))
        eigenre(i)=-xxx 
       endif
      enddo
c
      do i=1,nrot
      do j=1,nrot
       w(i,j)=vr(i,j)
       wt(i,j)=vr(j,i)
      enddo
      enddo
c R-matrix R=Wcos(tau)W^T + W (tau)^-1 sin(tau) W^T X^R -----
      do i=1,nrot
      do j=1,nrot
      do k=1,nrot
       full_r(i,j)=full_r(i,j)+w(i,k)*dcos(dsqrt(-eigenre(k)))*wt(k,j)
      enddo
      enddo
      enddo
c
      do i=1,nrot
      do j=1,nrot
      do k=1,nrot
       if(dabs(eigenre(k)).lt.1.0d-10) then
        del_k=1.0d0
       else
        del_k=(dsin(dsqrt(-eigenre(k)))/(dsqrt(-eigenre(k))))
       endif
      do l=1,nrot
       full_r(i,j)=full_r(i,j)+w(i,k)
     &       *del_k
     &       *wt(k,l)*xrot(l,j)
      enddo
      enddo
      enddo
      enddo
c rotating blockm --> blockm_r
      do i=1,nbf
      do j=1,nrot
      do k=1,nrot
       blockm_r(i,j)=blockm_r(i,j)+blockm(i,k)*full_r(k,j)
      enddo
      enddo
      enddo
c insert the block
      do i=1,nbf
      do j=1,nrot
       gmov(i,j+frozc+q_ele_a-nacto)=blockm_r(i,j)
      enddo
      enddo
c
      call util_flush(6)
c ---------------------------------
c puting block into g_movecs
c
      if(nodezero) then 
       call ga_put(g_movecs(1),1,nbf,1,nmo(1),
     1         gmov,nbf)
      endif
c
      return 
      end
c
c
c
c
c
c
      subroutine xrot_ini(nrot,xrot)
#include "tce_prop.fh"
#include "errquit.fh"
      integer nrot
      integer i,j,k
      double precision xrot(nrot,nrot)
c
      do i=1,nrot
      do j=1,nrot
       xrot(i,j)=0.0d0
      enddo
      enddo 
c
c check if nrot*(nrot-1)/2 < maxnumrot
      if((nrot*(nrot-1)/2).gt.maxnumrot) 
     &  call errquit('xrto_ini: maxnumrot limit met',info,
     & unknown_err)
c
       k=1
       do i=1,nrot
       do j=i+1,nrot
        xrot(i,j) =  arot(k)
        xrot(j,i) = -arot(k)
        k=k+1
       enddo
       enddo
c
      return
      end 


