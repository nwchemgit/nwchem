C $Id: react.F,v 1.2 2001-01-08 15:52:00 windus Exp $
      Subroutine drdy_react(jop)
C
C Calculate potential and hessian for reactant and product geometries
C  and write out to file 30      
C
      Implicit Double Precision (A-H,O-Z)
C
#include "drdyP.fh"
C
C  Geometries for reactants or products read in subroutine input
C
      Do iloop = 1,2
         iop = 2*iloop - 1
         if (iop.gt.jop) then
            do i = 1,n3
               x(i) = xr(i,iop)
            end do
            if (iop.eq.1) then
               write (6,600)
               write (6,610)
            else
               write (6,620)
               write (6,610)
            endif
            write (6,640) (j+1,(x(3*j+i),i=1,3),j=0,natom-1)
C  Get potential and first and second derivatives
            nbar = (n3*(n3+1))/2
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
            write (6,650) v,v*627.5095
            write (6,640) (j+1,(dx(3*j+i),i=1,3),j=0,natom-1)
C  Write out information to file 30
            Call drdy_react2(iop)
         endif
         if (iop.eq.1 .and. lgs(6).le.2 .or.
     *    iop.eq.3 .and. lgs(6).eq.1.or.lgs(6).eq.3) then
            iop = iop+1
            if (iop.gt.jop) then
               v = 0.0
               call drdy_react2(iop)
            endif
         end if
      enddo
      RETURN
C
600   Format(//1x,10(1h*),' Reactants',/)
610   Format(1x,'Geometry in mass-scaled cartesians (bohrs):',
     *   //,17x,1hx,15x,1hy,15x,1hz)
620   Format(//1x,10(1h*),' Products',/)
640   Format(1x,i3,4x,1p,3e16.6)
650   Format(/,' Potential energy=',1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians',
     *   ' (au):')
      END
