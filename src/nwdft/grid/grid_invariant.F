      Subroutine grid_invariant
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
c
      logical oprint
      integer lcoord, icoord, lcharge, icharge, ltags, itags
      double precision matalr(3,3),c(3,3)
      integer numgr3,n
      logical dogr3rot
      parameter (numgr3=6)
      character*10 grrname,gr3name(numgr3)
      data gr3name/
     1'C3        ',
     3'D3        ',
     5'C3v       ',
     7'C3h       ',
     9'D3h       ',
     1'D3d       '/
cedo      parameter (numgr3=13)
cedo      character*10 grrname,gr3name(numgr3)
cedo      data gr3name/
cedo     1'C3        ','C6        ',
cedo     3'D3        ','D6        ',
cedo     5'C3v       ','C6v       ',
cedo     7'C3h       ','C6h       ',
cedo     9'D3h       ','D6h       ',
cedo     1'D3d       ','D6d       ',' S6       '/
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
c
c     Define the rotation matrix that will be used to ensure that the quadrature
c     grid is rotationally invariant.
c
      if (iqrot.eq.0)then
c
c        Turn off the rotational invariance by defining a unit matrix.
c
         qrot(1,1) = 1.D0
         qrot(2,1) = 0.D0
         qrot(3,1) = 0.D0
         qrot(1,2) = 0.D0
         qrot(2,2) = 1.D0
         qrot(3,2) = 0.D0
         qrot(1,3) = 0.D0
         qrot(2,3) = 0.D0
         qrot(3,3) = 1.D0
c
      else
c
c        allocate space for atomic coordinates and charges
c
         if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &      icoord))call errquit(
     .        'grid_invariant: failed to alloc coordinates',0, MA_ERR)
         if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &      icharge))call errquit(
     .        'grid_invariant: failed to alloc charges',0, MA_ERR)
         if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &      ltags, itags))call errquit(
     /        'grid_invariant: failed to alloc center tags',0, MA_ERR)
c     
         if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                           Dbl_MB(icoord), Dbl_MB(icharge)))
     &        call errquit('grid_invariant: geom_cart_get failed',74,
     &       GEOM_ERR)
c
         call grid_invar_mat(dbl_mb(icoord), dbl_mb(icharge))
         if (.not. MA_Pop_Stack(ltags))
     &      call errquit('grid_invariant: pop stack failed.',0, MA_ERR)
         if (.not. MA_Pop_Stack(lcharge))
     &      call errquit('grid_invariant: pop stack failed.',0, MA_ERR)
         if (.not. MA_Pop_Stack(lcoord))
     &      call errquit('grid_invariant: pop stack failed.',0, MA_ERR)
      endif
c
c     rotate by Ahlrichs prescription for groups with 3-fold axis
c
      call sym_group_name(geom, grrname)
!      write(LuOut,'("GRGR",A,"!")') grrname
      dogr3rot=.false.
      do n=1,numgr3
         if(grrname.eq.gr3name(n)) dogr3rot=.true.
      enddo
      if(dogr3rot) then
      if(ga_nodeid().eq.0) write(LuOut,*) ' Rotation of axis '
         matalr(1,1)= sqrt(6d0)/3d0
         matalr(2,1)=0d0
         matalr(3,1)= sqrt(3d0)/3d0
         matalr(1,2)=-sqrt(6d0)/6d0
         matalr(2,2)= sqrt(2d0)/2d0
         matalr(3,2)= sqrt(3d0)/3d0
         matalr(1,3)=-sqrt(6d0)/6d0
         matalr(2,3)=-sqrt(2d0)/2d0
         matalr(3,3)= sqrt(3d0)/3d0
         call dcopy(9,qrot,1,c,1)
         call dgemm('N','N',3,3,3,1d0,matalr,3,c,3,
     ,        0d0,qrot,3)
      endif
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint)then
         write(LuOut,*)' qrot: '
         call output(qrot, 1, 3, 1, 3, 3, 3, 1)
      endif
      return
      end
      Subroutine grid_invar_mat(coord, charge)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c
#include "stdio.fh"
#include "cdft.fh"
#include "global.fh"
#include "util.fh"

c      
      double precision eps
      parameter (eps = 1.d-9)
      double precision coord(3,ncenters), charge(ncenters)
c
c     allocate space for M and T (rotational invariance; 
c     see B.G. Johnson et al, CPL, 220, 377 (1994), and auxiliary 
c     species R_a - T and (R_a -T)(R_a - T)^T
c
      double precision tvec(3), ramt(3), ramtm(3,3), eigs(3)
      double precision ztot, lramt,coss,dett
c
      integer ictr, icart, jcart
      double precision deter3,qin(3,3),ddiag,dasum,
     ,     valratio
      external deter3,dasum
      logical oprint
c
      oprint=util_print('quadrature', print_high)
      call dcopy(9,0d0,0,qrot, 1)
      call dcopy(3,0d0,0,tvec, 1)
c
c     atomic case 
c
      if(ncenters.eq.1) goto 111
c
c     calculate the T vector
c
      ztot = 0.d0
      do ictr = 1, ncenters
         ztot = ztot+charge(ictr)
         do icart = 1, 3
            tvec(icart) = tvec(icart) + 
     &                    charge(ictr)*coord(icart,ictr)
         enddo
      enddo
c
c     normalize T
c
      if (abs(ztot).gt.EPS) then
         lramt=0d0
         do icart = 1, 3
            tvec(icart) = tvec(icart)/ztot
            lramt=lramt+abs(tvec(icart))
         enddo
      else
         write(LuOut,*)
         write(LuOut,*) 'grid_invariant:  ',
     &                  'The total nuclear charge is zero.'
         call errquit('Exiting in grid_invariant.',4, UNKNOWN_ERR)
      endif
c
c     calculate the M matrx
c
      do ictr = 1, ncenters
c
c        Form R_a - T
c
         call dcopy(3,0d0,0,ramt, 1)
         lramt = 0.d0
         do icart = 1, 3
            ramt(icart) = coord(icart,ictr) - tvec(icart)
            lramt = lramt + ramt(icart)*ramt(icart)
         enddo
c
c        form (R_a -T)(R_a - T)^T
c
         do icart = 1, 3
            do jcart = 1, 3
               ramtm(icart,jcart) = ramt(icart)*ramt(jcart)
            enddo
         enddo
c
c        form M in qrot
c
         do icart = 1, 3
            qrot(icart,icart) = qrot(icart,icart) + 
     &                          charge(ictr)*lramt
            do jcart = 1, 3
               qrot(icart,jcart) = qrot(icart,jcart) - 
     &                             ramtm(icart,jcart)*charge(ictr)
            enddo
         enddo
      enddo
      valratio=(abs(qrot(1,2))+abs(qrot(1,3))+abs(qrot(2,3)))+
     +     (abs(qrot(1,1))+abs(qrot(2,2))+abs(qrot(3,3)))
      if(valratio.gt.1d-20) then
      valratio=(abs(qrot(1,2))+abs(qrot(1,3))+abs(qrot(2,3)))/
     /     (abs(qrot(1,1))+abs(qrot(2,2))+abs(qrot(3,3)))
      endif
      if(ga_nodeid().eq.0.and.oprint) then
         write(LuOut,*)' ratio offdiag/diag ',valratio
      endif
      if(valratio.lt.1d-3) then
         if(ga_nodeid().eq.0.and.oprint) then
            write(luout,*)
         endif
         call dfill(9,0d0,qrot,1)
         goto 111
      endif
c
c     diagonalize M
c
      call dcopy(9,qrot,1,qin,1)
      call dcopy(3,0d0,0,eigs, 1)
      call util_jacobinosort(3, qrot,3, eigs)
      if(ga_nodeid().eq.0.and.oprint) then
         write(LuOut,*)' nuclear charge qrot: '
         call output(qin, 1, 3, 1, 3, 3, 3, 1)
         write(LuOut,*)'eig: ',eigs
      endif
      
      ddiag=abs(eigs(1)-qin(1,1))+abs(eigs(2)-qin(2,2))+
     +     abs(eigs(3)-qin(3,3))
      if(abs(ddiag).lt.eps) then
         if(ga_nodeid().eq.0.and.oprint) then
            write(LuOut,*) ' qrot input diag '
         endif
         call dfill(9,0d0,qrot,1)
         goto 111
      endif
      dett=deter3(qrot)
      coss=(qrot(1,1)+qrot(2,2)+qrot(3,3)-1d0)*.5d0
      if ((dett-coss*coss).lt.0) then
         if(ga_nodeid().eq.0.and.oprint) then
            write(LuOut,*)'wqrot: '
            call output(qrot, 1, 3, 1, 3, 3, 3, 1)
            write(LuOut,*)'eig: ',eigs
         
            write(LuOut,*) ' WARNING: nonsense rotational matrix'
            write(LuOut,*) ' cosine      ',coss
            write(LuOut,*) ' sine square ',dett-coss*coss
            write(LuOut,*) ' Identity matrix is used '
            write(LuOut,*) '     '
         endif
         call dfill(9,0d0,qrot,1)
         goto 111
      endif
      return
 111  qrot(1,1)=1d0
      qrot(2,2)=1d0
      qrot(3,3)=1d0
      return
      end
