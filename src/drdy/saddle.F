C $Id: saddle.F,v 1.2 2001-01-08 15:52:01 windus Exp $
      Subroutine drdy_saddle
c
C Calculate potential and hessian for saddle point geometry
C  and write out to file 30      
C
      Implicit None ! Double Precision (A-H,O-Z)
#include "drdyP.fh"
#include "stdio.fh"
c
      Integer fu6
      Integer i, j, nbar
C
C
      fu6 = luout
      do i = 1, natom
         iatom(i) = i
      enddo
      do i = 1, n3
         x(i) = xr(i,5)
      enddo
c
      write (fu6,600)
      write (fu6,610) (j+1,(x(3*j+i),i=1,3),j=0,natom-1)
      s = 0.0d00
C  Get potential and first and second derivatives
      nbar = (n3*(n3+1))/2
      Call drdy_pot2(x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
      write (fu6,620) v,v*627.5095d00
      write (fu6,610) (j+1,(dx(3*j+i),i=1,3),j=0,natom-1)
c
      vspsv = v
      do i = 1,n3
         xspsv(i) = x(i)
         dxspsv(i) = dx(i)
         do j = 1,n3
            fspsv(j,i) = f(j,i)
         enddo
      enddo
      call drdy_rphwrt(5)
      return
c
600   Format(//1x,10(1h*),' Saddle point or starting point',//,
     * ' Geometry in mass-scaled cartesians (bohrs):',
     *   //,17x,1hx,15x,1hy,15x,1hz)
610   Format(1x,i3,4x,1p,3e16.6)
620   Format(/,' Potential energy=',1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians',
     *   ' (au):')
      end                                           
