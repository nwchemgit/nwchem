C $Id$
      Subroutine drdy_react(jop)
C
C Calculate potential and hessian for reactant and product geometries
C  and write out to file 30      
C
      Implicit none
      Integer jop
      Integer i,iloop,iop,j,nbar
      Double Precision swf
C
#include "drdyP.fh"
C
C  Geometries for reactants or products read in subroutine input
C
      swf = sreact
      Do iloop = 1,2
         iop = 2*iloop - 1
         if (iop.gt.jop) then
            do i = 1,n3
               x(i) = xr(i,iop)
            end do
            if (iop.eq.1) then
               write (fu6,600)
               write (fu6,610)
            else
               write (fu6,620)
               write (fu6,610)
            endif
            write (fu6,640) (j+1,(x(3*j+i),i=1,3),j=0,natom-1)
C  Get potential and first and second derivatives
            nbar = (n3*(n3+1))/2
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
            if (iloop.eq.1) then
               vr = v
            else
               vp = v
            endif
            if(lgs(38).gt.0) call drdy_wfiles(swf,0)
            if (lgs(39).eq.0) then
               write (fu6,650) v,v*ckcal
            else
               call drdy_potsp(x,vspc,scr1,amass,natom,n3)
               if (iloop.eq.1) then
                  vrsp = vspc
               else
                  vpsp = vspc
               endif
               write (fu6,651) v,v*ckcal,vspc,vspc*ckcal
               if(lgs(38).eq.4) call drdy_wfiles(swf,1)
            endif
            write (fu6,640) (j+1,(dx(3*j+i),i=1,3),j=0,natom-1)
C  put correct parts of hessian matrix into fsv and call rphwrt
            Call drdy_react2(iop)
         endif
         if ((iop.eq.1 .and. lgs(6).le.2) .or.
     *    (iop.eq.3 .and. (lgs(6).eq.1.or.lgs(6).eq.3))) then
            iop = iop+1
            if (iop.gt.jop) then
               v = 0.0
C  put correct parts of hessian matrix into fsv and call rphwrt
               call drdy_react2(iop)
            endif
         end if
         swf = sprod
      enddo
      RETURN
C
600   Format(//1x,10(1h*),' Reactants',/)
610   Format(1x,'Geometry in mass-scaled cartesians (bohrs):',
     *   //,17x,1hx,15x,1hy,15x,1hz)
620   Format(//1x,10(1h*),' Products',/)
640   Format(1x,i3,4x,1p,3e16.6)
650   Format(/,' Potential energy=',T20,1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians',
     *   ' (au):')
651   Format(/,' Potential energy=',T20,1pe18.10,' Hartree=',e13.5,
     *   ' kcal/mol',/,5x,'Single point',T20,1pe18.10,' Hartree=',
     *   e13.5,
     *   ' kcal/mol',//,' Derivatives in mass-scaled cartesians (au):')
      END
