c
c     calculate the oscillator strength including higher order contributions
c
      subroutine tddft_oscstr(tmom, eexcit, mag_transdip, 
     &                        oscstr, oscstrtot)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
c
      double precision tmom(20), eexcit, oscstr(3),oscstrtot
      double precision mag_transdip(3)
c
c     Local variables:
c
      double precision edip, equad, mdip, edip_eocta, edip_mquad
      double precision trace, qmat(3,3), rmat(3,3), tmat(3,3)
      integer i,j, k, l, m, n
c
c     Define translation table for quadrupole incices in
c     packed storage in tmom(i) below
c     XX=5,XY=YX=6, XZ=ZX=7,YY=8,YZ=ZY=9,ZZ=10
      integer qindex(3,3)
      data qindex(1,1),qindex(1,2),qindex(1,3)/5,6,7/
      data qindex(2,1),qindex(2,2),qindex(2,3)/6,8,9/
      data qindex(3,1),qindex(3,2),qindex(3,3)/7,9,10/
c
      double precision cl, rau2cgs
      parameter (cl = 137.0359895d0, rau2cgs = 235.726327d0)
      double precision zero, one, two, three, half, third
      parameter (zero = 0d0, one=1d0, two=2d0, three=3d0,
     &   half=one/two, third=one/three)
c
      logical nodezero
c
      nodezero=(ga_nodeid().eq.0)
c
c     Electric dipole contribution
      edip = tmom(2)**2+tmom(3)**2+tmom(4)**2
      edip = (2.d0/3.d0)*eexcit*edip
c
c     Electric quadrupole contribution
      do i = 1,3
        do j = 1,3
          qmat(i,j) = tmom(qindex(i,j))
          if (j.ne.i) qmat(j,i) = tmom(qindex(j,i))
        end do
      end do
c
c     traceless quadrupole tensor
      trace = zero
      do i = 1,3
        trace = trace + tmom(qindex(i,i))
      end do
      do i = 1,3
        qmat(i,i) = qmat(i,i) - third*trace
      end do
c
      equad = 0.d0
      do i = 1,3
       do j = 1,3
         equad = equad + qmat(i,j)*qmat(i,j)
       end do
      end do
      equad = (1.d0/(20.d0*cl*cl))*eexcit*eexcit*eexcit*equad
c
c     debug
c     trace = zero
c     do i = 1,3
c        trace = trace + qmat(i,i)
c     end do
c     if (nodezero) then
c       write (luout,*) 'traceless quadrupole moment tensor'
c       call output(qmat,1,3,1,3,3,3,1)
c       write (luout,'(1X,a,f15.5)') 'Trace = ', trace
c     end if
c
c     Magnetic dipole contribution (length)
      mdip = mag_transdip(1)*mag_transdip(1)
     &     + mag_transdip(2)*mag_transdip(2)
     &     + mag_transdip(3)*mag_transdip(3)
      mdip = (2.d0/(3.d0*cl*cl))*eexcit*mdip
c
c     el_dipole*el_octapole
c     tmom(2)  = x
c     tmom(3)  = y
c     tmom(4)  = z
c     tmom(11) = xxx
c     tmom(12) = xxy = yxx = xyz
c     tmom(13) = xxz = zxx = xzx
c     tmom(14) = xyy = yxy = yyx
c     tmom(15) = xyz = zxy = yzx = zyx = yxz = xzy
c     tmom(16) = xzz = zxz = zzx
c     tmom(17) = yyy 
c     tmom(18) = yyz = zyy = yzy
c     tmom(19) = yzz = zyz = zzy
c     tmom(20) = zzz
c      edip_eocta = tmom(2)*tmom(11)+tmom(3)*tmom(12)+tmom(4)*tmom(13)
c     &           + tmom(2)*tmom(14)+tmom(3)*tmom(17)+tmom(4)*tmom(18)
c     &           + tmom(2)*tmom(16)+tmom(3)*tmom(19)+tmom(4)*tmom(20)
c      edip_eocta = -(2.d0/(45.d0*cl*cl))*eexcit*eexcit*eexcit*edip_eocta
c
c     el_dipole*mag_quadrupole
c     edip_mquad = 0.d0  ! not implemented yet
c
c     assemble oscillator strength
      oscstr(1) = edip
      oscstr(2) = equad
      oscstr(3) = mdip
      oscstrtot = oscstr(1) + oscstr(2) + oscstr(3)
c
c     debug
c     if (nodezero) then
c        write (luout,*) 'oscstrtot :', oscstrtot
c        write (luout,*) 'edip      :', edip
c        write (luout,*) 'equad     :', equad
c        write (luout,*) 'mdip      :', mdip
c        write (luout,*) 'edip_eocta:', edip_eocta
c     end if
c
      end
