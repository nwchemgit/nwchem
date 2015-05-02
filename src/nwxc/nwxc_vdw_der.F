c
c     empirical dispersion: derivatives
c
      subroutine nwxc_vdw_der(s6,s8,sr6,sr8,a1,a2,n,x,z,force)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "nwxc_vdw.fh"
#include "global.fh"
c
      double precision s6,s8,sr6,sr8,a1,a2
      integer n
      double precision x(3,n),force(3,n)
      integer z(n)
c
      integer i,j,k,A,l_cnij,k_cnij,l_cnijk,k_cnijk
      double precision nwxc_c6ij_sk
      external nwxc_c6ij_sk
      double precision drajdxa
      double precision ff1,rr,ff
      double precision nwxc_fdmp,nwxc_f1dmp,cnA,cnj
      external nwxc_c6cn,nwxc_crd_nr,nwxc_fdmp,nwxc_f1dmp
      double precision nwxc_c6cn,nwxc_crd_nr
      double precision fac6,fac8,fdmp6,fdmp6a,fdmp8,fdmp8a,Qfac
      double precision rAj,rAk,rjk,r0aj,r0ak,r0jk,c6Aj,grad_c6(3)
      double precision dxAj,dyAj,dzAj,dxAk,dyAk,dzAk,dxjk,dyjk,dzjk
      double precision tmp6,tmp6a,tmp8,tmp8a
      double precision nwxc_fdmpbj, nwxc_fdmpbj_d1
      external         nwxc_fdmpbj, nwxc_fdmpbj_d1
c
c     Derivatives of Grimme dispersion term
c
c  DFT-D1 / DFT-D2
c
      if (ivdw.le.2) then
         do A=1,n
            force(1,A)=0d0
            force(2,A)=0d0
            force(3,A)=0d0
            if (Z(A).ne.0) then
              do j=1,n
                 if(A.ne.j.and.Z(j).ne.0) then
                    rAj=sqrt(
     +                 (x(1,A)-x(1,j))**2 +
     +                 (x(2,A)-x(2,j))**2 +
     +                 (x(3,A)-x(3,j))**2)
                    r0aj=r0(z(A))+r0(z(j))
                    ff= nwxc_fdmp(rAj,r0aj)
                    ff1= nwxc_f1dmp(rAj,r0aj,ff)
                    rr=nwxc_c6ij_sk(A,j,z)/(rAj**6)*
     *               ((-6d0*ff/rAj)+ff1)
                    do i=1,3
                       drAjdxa=(x(i,A)-x(i,j))/rAj
                       force(i,A)=force(i,A)-rr*drAjdxa
                    enddo
                 endif
              enddo
            endif
         enddo
         if(abs(s6-1d0).gt.1d-9) 
     F        call dscal(3*n,s6,force,1)
c
c DFT-D3
c
      else if (ivdw.eq.3) then
c
c        Precompute coordinate derivatives C6 dependency
c
         if (.not.ma_push_get(mt_dbl,3*n,'xcvdw cnij',l_cnij,k_cnij))
     &      call errquit('xcvdw cnij: cannot allocate cnij',0, MA_ERR)
         if (.not.ma_push_get(mt_dbl,3*n*n,'vdw cnijk',l_cnijk,k_cnijk))
     &      call errquit('vdw cnijk: cannot allocate cnijk',0, MA_ERR)
c
         call nwxc_crd_nr_der(n,x,z,dbl_mb(k_cnij),dbl_mb(k_cnijk))
c
         do A=1,n
           force(1,A)=0.0d0
           force(2,A)=0.0d0
           force(3,A)=0.0d0
           if (Z(A).ne.0) then
             do j=1,n
               if(A.ne.j.and.Z(j).ne.0) then      
                  dxAj=x(1,A)-x(1,j)
                  dyAj=x(2,A)-x(2,j)
                  dzAj=x(3,A)-x(3,j)
                  rAj=dxAj**2+dyAj**2+dzAj**2
c
c                 Two center derivatives. Grimme uses screening to reduce 
c                 computational work
c
c                 Screening r^2 distance vs threshold of 20000.0
c
                  if (rAj.gt.20000.d0) goto 901
c
c                 Factors
c
                  r0aj=r0AB(z(A),z(j))
                  Qfac=Qatom(z(A))*Qatom(z(j))
                  fac6=(dsqrt(rAj)/(sr6*r0aj))**(-alpha)
                  fac8=(dsqrt(rAj)/(sr8*r0aj))**(-(alpha+2.0d0))
                  fdmp6=1.0d0/(1.0d0+6.0d0*fac6)
                  fdmp8=1.0d0/(1.0d0+6.0d0*fac8)
c
c                 Coordination dependent C6_AB value
c
                  cnA=nwxc_crd_nr(A,n,x,z)
                  cnj=nwxc_crd_nr(j,n,x,z)
                  c6Aj=nwxc_c6cn(z(A),z(j),cnA,cnj)
c
c                 Get gradient for coordination number dependent C6
c
                  call nwxc_c6_grad(grad_c6,A,j,A,x,z,n,
     &                              dbl_mb(k_cnij),dbl_mb(k_cnijk))
c
                  tmp6=6.0d0*fdmp6*s6*c6Aj/(rAj**4.0d0)
                  tmp8=6.0d0*fdmp8*s8*c6Aj*Qfac/(rAj**5.0d0)
c
c                 dx contribution to A
c
                  tmp6a=tmp6*dxAj
                  tmp8a=tmp8*dxAj
                  force(1,A)=force(1,A)
     $              +(1.0d0-fdmp6*fac6*alpha)*tmp6a
     $              -fdmp6*s6*grad_c6(1)/(rAj**3.0d0)
     $              +(4.0d0-3.0d0*fdmp8*fac8*(alpha+2.0d0))*tmp8a
     $              -3.0d0*fdmp8*s8*grad_c6(1)*Qfac/(rAj**4.0d0)
c
c                 dy contribution to A
c
                  tmp6a=tmp6*dyAj
                  tmp8a=tmp8*dyAj
                  force(2,A)=force(2,A)
     $              +(1.0d0-fdmp6*fac6*alpha)*tmp6a
     $              -fdmp6*s6*grad_c6(2)/(rAj**3.0d0)
     $              +(4.0d0-3.0d0*fdmp8*fac8*(alpha+2.0d0))*tmp8a
     $              -3.0d0*fdmp8*s8*grad_c6(2)*Qfac/(rAj**4.0d0)
c
c                 dz contribution to A
c
                  tmp6a=tmp6*dzAj
                  tmp8a=tmp8*dzAj
                  force(3,A)=force(3,A)
     $              +(1.0d0-fdmp6*fac6*alpha)*tmp6a
     $              -fdmp6*s6*grad_c6(3)/(rAj**3.0d0)
     $              +(4.0d0-3.0d0*fdmp8*fac8*(alpha+2.0d0))*tmp8a
     $              -3.0d0*fdmp8*s8*grad_c6(3)*Qfac/(rAj**4.0d0)
 901              continue
               endif
             enddo
c
c            Three center derivatives. Grimme uses aggressive screening
c            to get this N^3 contribution back to N^2
c
             do j=2,n
               if(A.ne.j.and.Z(j).ne.0) then      
                  rAj=sqrt(
     +                 (x(1,A)-x(1,j))**2 +
     +                 (x(2,A)-x(2,j))**2 +
     +                 (x(3,A)-x(3,j))**2)
                  r0aj=r0AB(z(A),z(j))
c
c                 Screening per Grimme
c
                  if (rAj.gt.1600d0*r0aj/r0AB(1,1)) goto 910
c
c                 Third center involved
c
                  do k=1,j-1
                     if(A.ne.k.and.Z(k).ne.0) then      
                       dxAk=x(1,A)-x(1,k)
                       dyAk=x(2,A)-x(2,k)
                       dzAk=x(3,A)-x(3,k)
                       rAk=dxAk**2+dyAk**2+dzAk**2
                       r0ak=r0AB(z(A),z(k))
                       dxjk=x(1,j)-x(1,k)
                       dyjk=x(2,j)-x(2,k)
                       dzjk=x(3,j)-x(3,k)
                       rjk=dxjk**2+dyjk**2+dzjk**2
                       r0jk=r0AB(z(j),z(k))
c
c                      Screening r^2 distance vs threshold of 1600.0*(radii Ak)
c
                       if ((rAk.gt.1600.0d0*r0ak/r0AB(1,1)).or.
     $                     (rjk.gt.1600.0d0*r0jk/r0AB(1,1))) goto 911
c
c                      Get gradient for coordination number dependent C6 for three centers
c
                       call nwxc_c6_grad(grad_c6,j,k,A,x,z,n,
     &                                   dbl_mb(k_cnij),dbl_mb(k_cnijk))
                       fac6=(sr6*r0jk/dsqrt(rjk))**(alpha)
                       fac8=(sr8*r0jk/dsqrt(rjk))**(alpha+2.0d0)
                       fdmp6=1.0d0/(1.0d0+6.0d0*fac6)
                       fdmp8=1.0d0/(1.0d0+6.0d0*fac8)
c
c                      dx, dy, and dz contribution to A
c
                       Qfac=Qatom(z(j))*Qatom(z(k))
                       force(1,A)=force(1,A)
     $                      -fdmp6*s6*grad_c6(1)/(rjk**3.0d0)
     $                      -3.0d0*fdmp8*s8*grad_c6(1)*Qfac/(rjk**4.0d0)
                       force(2,A)=force(2,A)
     $                      -fdmp6*s6*grad_c6(2)/(rjk**3.0d0)
     $                      -3.0d0*fdmp8*s8*grad_c6(2)*Qfac/(rjk**4.0d0)
                       force(3,A)=force(3,A)
     $                      -fdmp6*s6*grad_c6(3)/(rjk**3.0d0)
     $                      -3.0d0*fdmp8*s8*grad_c6(3)*Qfac/(rjk**4.0d0)
                     endif
 911              continue
                  enddo
 910           continue
               endif
             enddo
           endif
         enddo
         if (.not.ma_pop_stack(l_cnijk))
     $      call errquit('xcvdw cnijk: cannot pop cnijk',4, MA_ERR)
         if (.not.ma_pop_stack(l_cnij))
     $      call errquit('xcvdw cnij: cannot pop cnij',4, MA_ERR)
c
c DFT-D3BJ
c
      else if (ivdw.eq.4) then
c
c        Precompute coordinate derivatives C6 dependency
c
         if (.not.ma_push_get(mt_dbl,3*n,'xcvdw cnij',l_cnij,k_cnij))
     &      call errquit('xcvdw cnij: cannot allocate cnij',0, MA_ERR)
         if (.not.ma_push_get(mt_dbl,3*n*n,'vdw cnijk',l_cnijk,k_cnijk))
     &      call errquit('vdw cnijk: cannot allocate cnijk',0, MA_ERR)
c
         call crd_nr_der(n,x,z,dbl_mb(k_cnij),dbl_mb(k_cnijk))
c
         do A=1,n
           force(1,A)=0.0d0
           force(2,A)=0.0d0
           force(3,A)=0.0d0
           if (Z(A).ne.0) then
             do j=1,n
               if(A.ne.j.and.Z(j).ne.0) then      
                  dxAj=x(1,A)-x(1,j)
                  dyAj=x(2,A)-x(2,j)
                  dzAj=x(3,A)-x(3,j)
                  rAj=dxAj**2+dyAj**2+dzAj**2
c
c                 Two center derivatives. Grimme uses screening to reduce 
c                 computational work
c
c                 Screening r^2 distance vs threshold of 20000.0
c
                  if (rAj.gt.20000.d0) goto 941
c
c                 Factors
c
                  rAj = dsqrt(rAj)
                  Qfac=Qatom(z(A))*Qatom(z(j))
c
c                 Coordination dependent C6_AB value
c
                  cnA=nwxc_crd_nr(A,n,x,z)
                  cnj=nwxc_crd_nr(j,n,x,z)
                  c6Aj=nwxc_c6cn(z(A),z(j),cnA,cnj)
                  c8=3.0d0*c6Aj*Qfac
                  r0Aj=dsqrt(3.0d0*Qfac)
c
c                 Get gradient for coordination number dependent C6
c
                  call nwxc_c6_grad(grad_c6,A,j,A,x,z,n,
     &                              dbl_mb(k_cnij),dbl_mb(k_cnijk))
c
c                 dx contribution to A
c
                  force(1,A)=force(1,A)
     $              -s6*c6Aj*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,6)
     $               *dxAj/rAj
     $              -s8*c8*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,8)
     $               *dxAj/rAj
     $              -s6*grad_c6(1)*nwxc_fdmpbj(rAj,r0Aj,a1,a2,6)
     $              -s8*3.0d0*grad_c6(1)*Qfac*nwxc_fdmpbj(rAj,r0Aj,
     $                a1,a2,8)
c
c                 dy contribution to A
c
                  force(2,A)=force(2,A)
     $              -s6*c6Aj*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,6)
     $               *dyAj/rAj
     $              -s8*c8*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,8)
     $               *dyAj/rAj
     $              -s6*grad_c6(2)*nwxc_fdmpbj(rAj,r0Aj,a1,a2,6)
     $              -s8*3.0d0*grad_c6(2)*Qfac*nwxc_fdmpbj(rAj,r0Aj,
     $                a1,a2,8)
c
c                 dz contribution to A
c
                  force(3,A)=force(3,A)
     $              -s6*c6Aj*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,6)
     $               *dzAj/rAj
     $              -s8*c8*nwxc_fdmpbj_d1(rAj,r0Aj,a1,a2,8)
     $               *dzAj/rAj
     $              -s6*grad_c6(3)*nwxc_fdmpbj(rAj,r0Aj,a1,a2,6)
     $              -s8*3.0d0*grad_c6(3)*Qfac*nwxc_fdmpbj(rAj,r0Aj,
     $                a1,a2,8)
 941              continue
               endif
             enddo
c
c            Three center derivatives. Grimme uses aggressive screening
c            to get this N^3 contribution back to N^2
c
             do j=1,n
               if(A.ne.j.and.z(j).ne.0) then      
                  rAj=sqrt(
     +                 (x(1,A)-x(1,j))**2 +
     +                 (x(2,A)-x(2,j))**2 +
     +                 (x(3,A)-x(3,j))**2)
                  r0aj=r0AB(z(A),z(j))
c
c                 Screening per Grimme
c
                  if (rAj.gt.1600d0*r0aj/r0AB(1,1)) goto 950
c
c                 Third center involved
c
                  do k=1,n
                     if(A.ne.k.and.k.ne.j.and.z(k).ne.0) then      
                       dxAk=x(1,A)-x(1,k)
                       dyAk=x(2,A)-x(2,k)
                       dzAk=x(3,A)-x(3,k)
                       rAk=dxAk**2+dyAk**2+dzAk**2
                       r0ak=r0AB(z(A),z(k))
                       dxjk=x(1,j)-x(1,k)
                       dyjk=x(2,j)-x(2,k)
                       dzjk=x(3,j)-x(3,k)
                       rjk=dxjk**2+dyjk**2+dzjk**2
                       r0jk=r0AB(z(j),z(k))
c
c                      Screening r^2 distance vs threshold of 1600.0*(radii Ak)
c
                       if ((rAk.gt.1600.0d0*r0ak/r0AB(1,1)).or.
     $                     (rjk.gt.1600.0d0*r0jk/r0AB(1,1))) goto 951
c
c                      Get gradient for coordination number dependent C6 for three centers
c
                       call nwxc_c6_grad(grad_c6,j,k,A,x,z,n,
     &                                   dbl_mb(k_cnij),dbl_mb(k_cnijk))
                       rjk=dsqrt(rjk)
c
c                      dx, dy, and dz contribution to A
c
                       Qfac=Qatom(z(j))*Qatom(z(k))
                       fdmp6=nwxc_fdmpbj(rjk,dsqrt(3.0d0*Qfac),a1,a2,6)
                       fdmp8=nwxc_fdmpbj(rjk,dsqrt(3.0d0*Qfac),a1,a2,8)
                       force(1,A)=force(1,A)
     $                           -fdmp6*s6*grad_c6(1)
     $                           -3.0d0*fdmp8*s8*grad_c6(1)*Qfac
                       force(2,A)=force(2,A)
     $                           -fdmp6*s6*grad_c6(2)
     $                           -3.0d0*fdmp8*s8*grad_c6(2)*Qfac
                       force(3,A)=force(3,A)
     $                           -fdmp6*s6*grad_c6(3)
     $                           -3.0d0*fdmp8*s8*grad_c6(3)*Qfac
                     endif
 951              continue
                  enddo
 950           continue
               endif
             enddo
           endif
         enddo
         if (.not.ma_pop_stack(l_cnijk))
     $      call errquit('xcvdw cnijk: cannot pop cnijk',4, MA_ERR)
         if (.not.ma_pop_stack(l_cnij))
     $      call errquit('xcvdw cnij: cannot pop cnij',4, MA_ERR)
      endif
c
#ifdef DEBUG
      write(6,*) ' gradient vdw called'
#endif
      return
      end
c $Id$
