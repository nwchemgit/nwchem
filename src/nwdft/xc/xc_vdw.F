c     $Id: xc_vdw.F,v 1.2 2005-03-28 23:45:46 edo Exp $
      subroutine xc_vdw(rtdb,geom,exc,force,what,s6)
      implicit none
      character *(*) what
      integer geom,rtdb
      double precision exc,force(*),s6
#include "geom.fh"
#include "mafdecls.fh"
#include "errquit.fh"
#include "util.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer n
      integer itags,ltags,i_xyz,l_xyz,icharge,lcharge,
     I     l_fvdw,i_fvdw
      external xc_vdw_e
      double precision xc_vdw_e,evdw,scalea
      integer iz,lz,i
      logical xc_vdw_init
      external xc_vdw_init
      logical oprint,oprinth
c
      oprint = util_print('vdw', print_medium)
      oprinth = util_print('vdw high', print_high)

c
      if (.not. geom_ncent(geom, n)) 
     &   call errquit('xcvdw: geom_ncent failed',geom, GEOM_ERR)
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_xyz,i_xyz))
     &   call errquit('xcvdw: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,n,'atns',lz,iz))
     &   call errquit('xcvdw: cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,n,'charge',lcharge,icharge))
     &   call errquit('xcvdw: cannot allocate charge',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,n*16,'tags',ltags,itags))
     &   call errquit('xcvdw: cannot allocate tags',0, MA_ERR)
      if (.not. geom_cart_get(geom, n, Byte_MB(itags),
     &                        Dbl_MB(i_xyz), Dbl_MB(icharge)))
     &   call errquit('xcvdw: geom_cart_get failed',74, GEOM_ERR)
      if (.not.ma_pop_stack(ltags))
     &   call errquit('xcvdw: cannot pop stack',2, MA_ERR)
c
c     get atns
c
      do i=1,n
          int_mb(iz+i-1)=dbl_mb(icharge+i-1)
      enddo
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('xcvdw: cannot pop stack',3, MA_ERR)
c
c     conversion factor angs 2 au
c
       if(.not.geom_get_ang2au(geom, scalea)) call
     S     errquit('xcvdw: gang2au failed',0,0)
      if(.not.xc_vdw_init(scalea))
     &   call errquit('xcvdw: vwdinit failed',0, 0)
      if(what.eq.'energy') then
         evdw=xc_vdw_e(s6,n,dbl_mb(i_xyz),int_mb(iz))
         if(oprint.and.ga_nodeid().eq.0) 
     W        write(luout,*) ' vdW contrib ',evdw
         Exc=Exc+evdw
      elseif(what.eq.'forces') then
      if (.not.MA_push_get(MT_Dbl,n*3,'xyz',l_fvdw,i_fvdw))
     &   call errquit('xcvdw: cannot allocate forcev',0, MA_ERR)
         call xc_vdw_der(s6,n,dbl_mb(i_xyz),int_mb(iz),
     D     dbl_mb(i_fvdw),scalea)
         if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' vdW contrib '
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               dbl_mb(i_fvdw+3*(i-1)),
     Y               dbl_mb(i_fvdw+3*(i-1)+1),
     Z        dbl_mb(i_fvdw+3*(i-1)+2)
             enddo
             write(luout,*) ' before vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo
          endif
          call daxpy(3*n,1d0,dbl_mb(i_fvdw),1,force,1)
          if(oprinth.and.ga_nodeid().eq.0) then
             write(luout,*) ' after vdw contr @@@@@'
             do i=1,n
                write(luout,'(I2,3F10.7," F = ",3(1PE13.5))')
     Z               int_mb(iz+i-1),
     X               dbl_mb(i_xyz+3*(i-1)),
     Y               dbl_mb(i_xyz+3*(i-1)+1),
     Z               dbl_mb(i_xyz+3*(i-1)+2),
     X               force(1+3*(i-1)),
     Y               force(1+3*(i-1)+1),
     Z               force(1+3*(i-1)+2)
             enddo
          endif
      endif
      if (.not.ma_chop_stack(l_xyz)) 
     C   call errquit('xcvdw: cannot pop stack',4, MA_ERR)
      return
      end
      subroutine xc_vdw_der(s6,n,x,z,force)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2694 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
      double precision s6
      integer n
      double precision x(3,*),force(3,*)
      integer z(*)
c
      integer i,j,A
      double precision c6ij_sk,rAj
      external c6ij_sk
      double precision drajdxa,r0aj
      double precision ff1,rr
#include "xc_vdw.fh"
c
c     deriv
c
      do A=1,n
         force(1,A)=0d0
         force(2,A)=0d0
         force(3,A)=0d0
         do j=1,n
            if(A.ne.j) then
               rAj=sqrt(
     +              (x(1,A)-x(1,j))**2 +
     +              (x(2,A)-x(2,j))**2 +
     +              (x(3,A)-x(3,j))**2)
               r0aj=r0(z(A))+r0(z(j))
               ff= fdmp(rAj,r0aj)
               ff1= f1dmp(rAj,r0aj,ff)
               rr=c6ij_sk(A,j,z)/(rAj**6)*
     *             (-6d0*ff/rAj+ff1)
               do i=1,3
                  drAjdxa=1d0/rAj*(x(i,A)-x(i,j))
                  force(i,A)=force(i,A)-rr*drAjdxa
               enddo
            endif
         enddo
      enddo
c
c     fudge factor if needed
c      
      if(abs(s6-1d0).gt.1d-9) 
     F     call dscal(3*n,s6,force,1)
      return
      end
      double precision function c6ij_sk(i,j,z)
      implicit none
      integer i,j
      integer z(*)
#include "xc_vdw.fh"
c     
      integer zi,zj
      double precision ni,nj
      external xc_neff
      double precision xc_neff,ci,cj
c
      zi=z(i)
      zj=z(j)
      ni=xc_neff(zi)
      nj=xc_neff(zj)
      if(c6(zi).eq.-1d0.or.c6(zj).eq.-1d0) then
         write(6,*) ' i j',i,j,' c6 ',c6(zi),c6(zj)
         stop 'cra'
      endif
      ci=c6(zi)
      cj=c6(zj)
      c6ij_sk=2d0*(ci*ci*cj*cj*ni*nj)**(1d0/3d0)/
     /     ((nj*nj*ci)**(1d0/3d0)+(ni*ni*cj)**(1d0/3d0))
      c6ij_sk=c6ij_sk
      return
      end
      double precision function xc_neff(z)
      implicit none
      integer z
c
      integer nv
c
      if(z.eq.1) then
         xc_neff=0.8d0
         return
      endif
      if(z.le.2) nv=z
      if(z.le.10) nv=z-2
      if(z.gt.10) stop 'neff stop '
      xc_neff=1.17d0+0.33d0*nv
      return
      end
      double precision function xc_vdw_e(s6,n,x,z)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2694 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
#include "xc_vdw.fh"
      double precision s6
      integer n
      double precision x(3,*)
      integer z(*)
c
      integer i,j
      double precision c6ij_sk,rij
      external c6ij_sk
c
      xc_vdw_e=0d0
      do i=1,n-1
         do j=i+1,n
            rij=sqrt(
     +           (x(1,i)-x(1,j))**2 +
     +           (x(2,i)-x(2,j))**2 +
     +           (x(3,i)-x(3,j))**2)
            xc_vdw_e=xc_vdw_e-c6ij_sk(i,j,z)*
     *           fdmp(rij,r0(z(i))+r0(z(j)))*
     *           (rij)**(-6)
         enddo
      enddo
c
c     fudge factor if needed
c      
      xc_vdw_e=xc_vdw_e*s6
!      write(0,*) ' xc_vdw in kcal/mol ',xc_vdw_e
!      xc_vdw_e=xc_vdw_e/au2kcalmol
      return
      end
      logical function xc_vdw_init(angs2au)
      implicit none
      double precision angs2au
#include "xc_vdw.fh"
c
      integer i
c
      do i=1,120
         c6(i)=-1
         r0(i)=-1
      enddo
cH
      c6(1)=0.163d0
      r0(1)=111d0
cC
      c6(6)=1.881d0
      r0(6)=161d0
cN
      c6(7)=1.11d0
      r0(7)=155d0
cO
      c6(8)=0.684d0
      r0(8)=149d0
cF
      c6(9)=0.57d0
      r0(9)=143d0
cNe
      c6(10)=0.45d0
      r0(10)=138d0
c
c     convert to au
c
      do i=1,120
c     j nm6/mol to  kcal nm6/mol  to au nm6 to au
         c6(i)=c6(i)*joule2kcal/au2kcalmol/((angs2nm/angs2au)**6)
c picom to angs to au
         r0(i)=r0(i)/angs2pm*angs2au
      enddo
      xc_vdw_init=.true.
      return
      end
