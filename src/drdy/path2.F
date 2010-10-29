C $Id$
      Subroutine drdy_path2(iopt,delss,delsh,sold,smax,xsign)
C
C  Calculate steepest descent path from saddle point
C     iopt  -  option for following MEP
C              =1, use Euler integrator
C              =2, use Page-McIver LQA algorithm
C              =3, use Page-McIver CLQA algorithm
C              =4, use Page-McIver CUBE algorithm
C
      Implicit none
#include "errquit.fh"
      Integer iopt
      Double Precision delss,delsh,sold,smax,xsign
      Integer i,ierr,j,nbar,nosc
      Double Precision delsf,dotp,dxmag,shess,swf
C
      Double Precision eps
c
#include "drdyP.fh"
C
      data eps/1.d-6/
C
      delsf = xsign*(s-sold)
      shess = s
      nosc = 0
      nbar = (n3*(n3+1))/2
      do while (s*xsign.lt.smax*xsign-dels*eps)
         if (iopt.le.1) then
            call drdy_euler (n3tm,n3,x,dx,dels)
         elseif (iopt.eq.2) then
            ierr = 0
            call drdy_pmlqa (n3tm,n3,nbar,x,dx,f,hess,amass,dels,
     *       vec0,vec1,u0,scr1,vec2,ierr)
            if (ierr.ne.0) then
               write(fu6,*) ' error in pmlqa'
               call errquit('error in pmlqa',555, UNKNOWN_ERR)
            endif
         elseif (iopt.eq.3) then
            ierr = 0
            call drdy_pmclqa (n3tm,n3,nbar,x,dx,f,fold,hess,amass,
     &            dels,delsf,vec0,vec1,u0,scr1,vec2,ierr)
            if (ierr.ne.0) then
               write(fu6,*) ' error in pmclqa'
               call errquit('error in pmlqa',555, UNKNOWN_ERR)
            endif
         elseif (iopt.ge.4) then
            call drdy_pmcube (n3tm,n3,x,dx,f,fold,amass,dels,delsf,
     *       vec0,vec1,vec2)
         endif
         call drdy_center (natom,n3,amass,x,scr1)
         s = s + dels*xsign
         delss = delss + dels
         delsh = delsh + dels
         if (delss.lt.delsv*(1.0d00-eps) .and.
     *    delsh.lt.delhss*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            if (iopt.ge.3) then
C Store previous hessian matrix in fold
               do i = 1,n3
                  do j = 1,n3
                     fold(j,i) = f(j,i)
                  enddo
               enddo
               sold = shess
               delsf = xsign*(s-sold)
               shess = s
            endif
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,
     *                 nbar)
            delsh = 0.0d00
         endif
         write(fu6,610) s,v-vzero,(x(i),i=1,n3)
         dxmag = 0.0d00
         do i = 1,n3
            dxmag = dxmag + dx(i)*dx(i)
         enddo
         dxmag = sqrt(dxmag)
         write(fu6,611) dxmag,(dx(i)/dxmag,i=1,n3)
C Check for oscillating gradient
         dotp = 0.0d00
         do i = 1,n3
            dotp = dotp + dxold(i)*dx(i)/dxmag
            dxold(i) = dx(i)/dxmag
         enddo
         if (dotp.lt. 0.0d00) then
            nosc = nosc + 1
            if (nosc .ge. 5) then
               write (fu6,6000)
               call errquit('The gradient is unstable',555, INPUT_ERR)
            endif 
         else
            nosc = 0
         endif
         if (delss.ge.delsv*(1.0d00-eps)) then
            swf = s
            if(lgs(38).gt.0) call drdy_wfiles(swf,0)
            if (lgs(39).ne.0) then
               call drdy_potsp(x,vspc,scr1,amass,natom,n3)
               write (fu6,612) vspc-vzerosp
               if(lgs(38).eq.4) call drdy_wfiles(swf,1)
            endif
            call drdy_rphwrt(7)
            if(lgs(2).ne.0) then
               ierr = 0
               call drdy_projct(ierr)
               if (ierr.eq.0) call drdy_fdiag(n3,nf(5),ierr,lgs(2))
            endif
            delss = 0.0d00
         endif
      enddo
      return
610   Format (t4,0pf12.5,1pe13.5,T33,'x=',(T35,0p9f9.4))
611   Format (15x,1pe13.5,T32,'dx=',(T35,0p9f9.4))
612   Format(11x,'Vsp=',1pe13.5)
6000  Format (//,' ********** The direction of the gradient vector',
     * ' has changed by more than 90 degrees for the last 5 steps.',/,
     * 12x,'Either the step size should be decreased, or another',
     * ' integrator should be tried.')
      end
