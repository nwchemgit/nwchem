      Subroutine drdy_path2(iopt,delss,sold,smax,xsign)
C
C  Calculate steepest descent path from saddle point
C     iopt  -  option for following MEP
C              =1, use Euler integrator
C              =2, use Page-McIver LQA algorithm
C              =3, use Page-McIver CLQA algorithm
C              =4, use Page-McIver CUBE algorithm
C
      Implicit Double Precision (A-H,O-Z)
c
#include "drdyP.fh"
C
      data eps/1.d-6/
C
      nbar = (n3*(n3+1))/2
      do while (s*xsign.lt.smax*xsign-dels*eps)
         if (iopt.eq.1) then
            call drdy_euler (n3tm,n3,x,dx,dels)
         elseif (iopt.eq.2) then
            call drdy_pmlqa (n3tm,n3,nbar,x,dx,f,hess,amass,dels,
     *       vec0,vec1,u0,scr1,vec2,ierr)
         elseif (iopt.eq.3) then
            delsf = xsign*(s - sold)
            call drdy_pmclqa (n3tm,n3,nbar,x,dx,f,fold,hess,amass,
     &            dels,delsf,vec0,vec1,u0,scr1,vec2,ierr)
         elseif (iopt.ge.4) then
            delsf = xsign*(s - sold)
            call drdy_pmcube (n3tm,n3,x,dx,f,fold,amass,dels,delsf,
     *       vec0,vec1,vec2)
         endif
         call drdy_center (natom,n3,amass,x,scr1)
         if (iopt.ge.3) then
C Store previous hessian matrix in fold
            do i = 1,n3
               do j = 1,n3
                  fold(j,i) = f(j,i)
               enddo
            enddo
            sold = s
         endif
         s = s + dels*xsign
         delss = delss + dels
         if (iopt.eq.1 .and. delss.lt.delsv*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,
     &           natom,n3,n3tm,nbar)
         endif
         write (6,610) s,v-vzero,(x(i),i=1,n3)
         write (6,611) (dx(i),i=1,n3)
c
c         call calcs(xspsv,x,sval)
c         write (6,*) 'sval = ',sval        
c
         if (delss.ge.delsv*(1.0d00-eps)) then
            call drdy_rphwrt(7)
            delss = 0.0d00
         endif
      enddo
      return
610   Format(1x,0pf10.5,1pe13.5,(t26,1p6e13.5))
611   Format(t26,1p6e13.5)
      end
