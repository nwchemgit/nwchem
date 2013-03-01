c
c     empirical dispersion: energy
c
      double precision function nwxc_vdw_e(s6,s8,sr6,sr8,n,x,z)
c
c     S. Grimme J Comp Chem 25, 1463 (2004)
c     U. Zimmerli, M Parrinello and P. Koumoutsakos, JCP. 120, 2693 (2004)
c     Q. Wu and W. Yang, JCP. 116, 515 (2002)
c
      implicit none
c
#include "nwxc_vdw.fh"
c
      double precision s6,s8,sr6,sr8
      integer n
      double precision x(3,n)
      integer z(n)
c
      integer i,j
      double precision fdmp, fdmp3, cni, cnj,c6d3
      double precision c6ij_sk,rij,c6cn,crd_nr,e6,e8
      external c6ij_sk,c6cn,nxtask,crd_nr,fdmp3
c
      nwxc_vdw_e=0d0
      e6=0.0d0
      e8=0.0d0
c
c DFT-D1 / DFT-D2
c
      if (ivdw.le.2) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
               nwxc_vdw_e=nwxc_vdw_e-c6ij_sk(i,j,z)*
     *            fdmp(rij,r0(z(i))+r0(z(j)))*
     *            (rij)**(-6.0d0)
            enddo
          endif
        enddo
        nwxc_vdw_e=nwxc_vdw_e*s6
c
c DFT-D3
c
c As off August, 2011 Grimme states: "Adding three-body corrections is
c currently not recommended, as very little is known about the three-
c body behaviour of common DFs in overlapping density regions."
c http://toc.uni-muenster.de/DFTD3/data/man.pdf, section 1.3.
c Hence the three-body terms have not been implemented.
c
c The reference to three-center derivatives in the gradient code
c refers to contributions that come from differentiating the
c coordination dependent dispersion coefficients.
c
      else if (ivdw.eq.3) then
        do i=1,n-1
          if (Z(i).ne.0) then
            do j=i+1,n
               rij=dsqrt(
     +            (x(1,i)-x(1,j))**2 +
     +            (x(2,i)-x(2,j))**2 +
     +            (x(3,i)-x(3,j))**2)
               cni=crd_nr(i,n,x,z)
               cnj=crd_nr(j,n,x,z)
               c6d3=c6cn(z(i),z(j),cni,cnj)
               c8=3.0d0*c6d3*Qatom(z(i))*Qatom(z(j))
               e6=e6-c6d3*fdmp3(rij,r0AB(z(i),z(j))*sr6,alpha)*
     *            (rij)**(-6.0d0)
               e8=e8-c8*fdmp3(rij,r0AB(z(i),z(j))*sr8,alpha+2.0d0)*
     *            (rij)**(-8.0d0)
            enddo
          endif
        enddo
        nwxc_vdw_e=e6*s6+e8*s8
      endif
c
      return
      end
