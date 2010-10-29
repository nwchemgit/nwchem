C $Id$
      Subroutine drdy_saddle
c
C Calculate potential and hessian for saddle or starting point geometry
C
      Implicit None 
#include "drdyP.fh"
#include "stdio.fh"
c
      Integer i,ierr,ii,ip,j,nbar
      Double Precision swf
C
C  Geometry for saddle or starting point read in subroutine input
C
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
      vspsv = v
      swf = s
      if(lgs(38).gt.0) call drdy_wfiles(swf,0)
      if (lgs(39).eq.0) then
         write (fu6,620) v,v*ckcal
      else
         call drdy_potsp(x,vspc,scr1,amass,natom,n3)
         vspspsv = vspc
         write (fu6,621) v,v*ckcal,vspc,vspc*ckcal
         if(lgs(38).eq.4) call drdy_wfiles(swf,1)
      endif
      write (fu6,610) (j+1,(dx(3*j+i),i=1,3),j=0,natom-1)
c
      do i = 1,n3
         xspsv(i) = x(i)
         dxspsv(i) = dx(i)
         do j = 1,n3
            fspsv(j,i) = f(j,i)
         enddo
      enddo
C
C If no saddle point, then project out gradient      
      ierr = 0
      if (lgs(1).eq.0) then
         call drdy_projct(ierr)
      else
         ii = 0
         do i = 1,n3
            do ip = 1,i
               ii = ii + 1
               hess(ii) = f(ip,i)
            enddo
         enddo
      endif
C  Diagonalize hessian and write out frequencies and eigenvectors
      if (ierr.eq.0) call drdy_fdiag(n3,nf(5)+1,ierr,1)
      if (ierr.eq.0) then
C Store frequencies in freqr array
         do i = 1,nf(5)+1
            freqr(i,5) = scr2(i)
         enddo
      endif
C      
      return
600   Format(//1x,10(1h*),' Saddle point or starting point',//,
     * ' Geometry in mass-scaled cartesians (bohrs):',
     *   //,17x,1hx,15x,1hy,15x,1hz)
610   Format(1x,i3,4x,1p,3e16.6)
620   Format(/,' Potential energy=',T20,1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians',
     *   ' (au):')
621   Format(/,' Potential energy=',T20,1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',/,5x,'Single point',T20,1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians (au):')
      end                                           
