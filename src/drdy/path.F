C $Id$
      Subroutine drdy_path(ns)
C
C  Calculate steepest descent path from saddle point
C     lgs(31)  -  option for following MEP
C              =1, use Euler integrator
C              =2, use Page-McIver LQA algorithm
C              =3, use Page-McIver CLQA algorithm
C              =4, use Page-McIver CUBE algorithm
C
      Implicit none
      Integer ns
      Integer i,ierr,ii,iopt,j,nbar
      Double Precision delss,delsh,dxmag,smax,sold,swf,xsign

      Double Precision eps
c
#include "drdyP.fh"
C
      Data eps/1.d-6/
C
      sold = 0.0d0
      iopt = lgs(31)
      nbar = (n3*(n3+1))/2
      write (fu6,600) dels,delsv,delhss,dir,slm,slp
      if (iopt.le.1) then
         write (fu6,601)
      else if (iopt.eq.2) then
         write (fu6,602)
      else if (iopt.eq.3) then
         write (fu6,603)
      else
         write (fu6,604)
      end if
C
      if (lgs(1).ne.0) then
C Set up for calculation on reactant side (if lgs(1)=0 then no saddle
C    point and reaction path is in product direction only)
         xsign = -1.0d00
         smax = slm
         if (ns.gt.0 .and. sgrid(1).lt.0.0d00) then
C Restart from grid point with most negative s value
            ii = 0
            dxmag = 0.0d00
            do i = 1,n3
               x(i) = xgrid(i,1)
               dx(i) = dxgrid(i,1)
               dxmag = dxmag + dx(i)*dx(i)
               do j = 1,i
                  ii = ii + 1
                  f(j,i) = hgrid(ii,1)
                  f(i,j) = hgrid(ii,1)
               enddo
            enddo
            dxmag = sqrt(dxmag)
            if (iopt.ge.3) then
C For PG CLQA and CUBE algorithms we need the Hessian at a previous point
               if (ns.ge.2) then
                  sold = sgrid(2)
                  ii = 0
                  do i = 1,n3
                     do j = 1,i
                        ii = ii + 1
                        fold(j,i) = hgrid(ii,2)
                        fold(i,j) = hgrid(ii,2)
                     enddo
                  enddo
               else
C Use saddle point Hessian as previous one
                  sold = 0.0
                  do i = 1,n3
                     do j = 1,n3
                        fold(j,i) = fspsv(j,i)
                     enddo
                  enddo
               endif
            endif
            write(fu6,605)
            write(fu6,606)
            write(fu6,607) sgrid(1),vgrid(1),(x(i),i=1,n3)
            write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
            s = sgrid(1)
            delss = 0.0d00
            delsh = 0.0d00
            do i = 1,n3
               dxold(i) = dx(i)/dxmag
            enddo
C
         else
C Take step off saddle point in reactant direction
            write(fu6,609)
            write(fu6,606)
            write(fu6,607) 0.,vspsv-vzero,(xspsv(i),i=1,n3)
            write(fu6,608) 0.,(vec0sv(i),i=1,n3)
            if (lgs(39).ne.0) write (fu6,6001) vspspsv-vzerosp
            do i = 1,n3
               x(i) = xspsv(i)-dels*(vec0sv(i)-0.5d00*dels*vec1sv(i))
               do j = 1,n3
                  f(j,i) = fspsv(j,i)
               enddo
            enddo
            if (iopt.ge.3) then
C For CLQA and CUBE, store hessian matrix in fold
               sold = 0.0d00
               do i = 1,n3
                  do j = 1,n3
                     fold(j,i) = f(j,i)
                  enddo
               enddo
            endif
            s = -dels
            delss = dels
            delsh = dels
C Calculate g and optionally f at new geometry
            if (iopt.lt.3 .and. delss.lt.delsv*(1.0d00-eps) .and.
     *       delsh.lt.delhss*(1.0d00-eps)) then
               Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
            else
               Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,
     *            n3tm,nbar)
               delsh = 0.0d00
            endif
            write(fu6,607) s,v-vzero,(x(i),i=1,n3)
            dxmag = 0.0d00
            do i = 1,n3
               dxmag = dxmag + dx(i)*dx(i)
            enddo
            dxmag = sqrt(dxmag)
            write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
            do i = 1,n3
               dxold(i) = dx(i)/dxmag
            enddo
            if (delss.ge.delsv*(1.0d00-eps)) then
               swf = s
               if(lgs(38).gt.0) call drdy_wfiles(swf,0)
               if (lgs(39).ne.0) then
                  call drdy_potsp(x,vspc,scr1,amass,natom,n3)
                  write (fu6,6001) vspc-vzerosp
                  if(lgs(38).eq.4) call drdy_wfiles(swf,1)
               endif
               call drdy_rphwrt (7)
               if(lgs(2).ne.0) then
                  ierr = 0
                  call drdy_projct(ierr)
                  if (ierr.eq.0) call drdy_fdiag(n3,nf(5),ierr,lgs(2))
               endif
               delss = 0.0d00
            endif
         endif
C
         Call drdy_path2(iopt,delss,delsh,sold,smax,xsign)
C
      endif
C
C Calculation on product side
      xsign = 1.0d00
      smax = slp
      if (((lgs(1).ne.0.and.ns.gt.0) .or. ns.gt.1) .and. 
     *   sgrid(ns).gt.0.0d00) then
C Restart from grid point with most positive s value
         ii = 0
         dxmag = 0.0
         do i = 1,n3
            x(i) = xgrid(i,ns)
            dx(i) = dxgrid(i,ns)
            dxmag = dxmag + dx(i)*dx(i)
            do j = 1,i
               ii = ii + 1
               f(j,i) = hgrid(ii,ns)
               f(i,j) = hgrid(ii,ns)
            enddo
         enddo
         dxmag = sqrt(dxmag)
         if (iopt.ge.3) then
C For PM CLQA and CUBE algorithms we need the Hessian at a previous point
            if (ns.ge.2) then
               sold = sgrid(ns-1)
               ii = 0
               do i = 1,n3
                  do j = 1,i
                     ii = ii + 1
                     fold(j,i) = hgrid(ii,ns-1)
                     fold(i,j) = hgrid(ii,ns-1)
                  enddo
               enddo
            else
C Use saddle point Hessian as previous one (note that for lgs(1)=0,
C  ns>1 so this section will not be called
               sold = 0.0d00
               do i = 1,n3
                  do j = 1,n3
                     fold(j,i) = fspsv(j,i)
                  enddo
               enddo
            endif
         endif
         write(fu6,610)
         write(fu6,606)
         write(fu6,607) sgrid(ns),vgrid(ns),(x(i),i=1,n3)
         write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
         do i = 1,n3
            dxold(i) = dx(i)/dxmag
         enddo
         s = sgrid(ns)
         delss = 0.0d00
         delsh = 0.0d00
C 
      elseif (lgs(1).ne.0) then
C Calculation for step off saddle point
         write(fu6,611)
         write(fu6,606)
         write(fu6,607) 0.,vspsv-vzero,(xspsv(i),i=1,n3)
         write(fu6,608) 0.,(vec0sv(i),i=1,n3)
C Take step off saddle point in product direction
         do i = 1,n3
            x(i) = xspsv(i) + dels*(vec0sv(i) +
     &           0.5d0*dels*vec1sv(i))
            do j = 1,n3
               f(j,i) = fspsv(j,i)
            enddo
         enddo
         if (iopt.ge.3) then
C For CLQA and CUBE, store hessian matrix in fold
            sold = 0.0d00
            do i = 1,n3
               do j = 1,n3
                  fold(j,i) = f(j,i)
               enddo
            enddo
         endif
         s = dels
         delss = dels
         delsh = dels
C Calculate g and optionally f at new geometry
         if (iopt.lt.3 .and. delss.lt.delsv*(1.0d00-eps) .and.
     *    delsh.lt.delhss*(1.0-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,
     *         n3tm,nbar)
            delsh = 0.0d00
         endif
         write(fu6,607) s,v-vzero,(x(i),i=1,n3)
         dxmag = 0.0
         do i = 1,n3
            dxmag = dxmag + dx(i)*dx(i)
         enddo
         dxmag = sqrt(dxmag)
         write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
         do i = 1,n3
            dxold(i) = dx(i)/dxmag
         enddo
         if (delss.ge.delsv*(1.0d00-eps)) then
            swf = s
            if(lgs(38).gt.0) call drdy_wfiles(swf,0)
            if (lgs(39).ne.0) then
               call drdy_potsp(x,vspc,scr1,amass,natom,n3)
               write (fu6,6001) vspc-vzerosp
               if(lgs(38).eq.4) call drdy_wfiles(swf,1)
            endif
            call drdy_rphwrt (7)
            if(lgs(2).ne.0) then
               ierr = 0
               call drdy_projct(ierr)
               if (ierr.eq.0) call drdy_fdiag(n3,nf(5),ierr,lgs(2))
            endif
            delss = 0.0d00
         endif
      else
C No saddle point, take step from starting geometry
         dxmag = 0.0d00
         do i = 1,n3
            x(i) = xspsv(i)
            dx(i) = dxspsv(i)
            dxmag = dxmag + dx(i)*dx(i)
            do j = 1,n3
               f(j,i) = fspsv(j,i)
C For CLQA and CUBE, store hessian matrix in fold
               fold(j,i) = fspsv(j,i)
            enddo
         enddo
         dxmag = sqrt(dxmag)
         write(fu6,612)
         write(fu6,606)
         write(fu6,607) 0.,vspsv-vzero,(x(i),i=1,n3)
         write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
         if (iopt.le.1) then
            call drdy_euler (n3tm,n3,x,dx,dels)
         else
            ierr = 0
            call drdy_pmlqa (n3tm,n3,nbar,x,dx,f,hess,amass,
     *       dels,vec0,vec1,u0,scr1,vec2,ierr)
            sold = 0.0
         endif
         s = dels
         delss = dels
         delsh = dels
C Calculate g and f at new geometry
         if (iopt.lt.3 .and. delss.lt.delsv*(1.0d00-eps) .and.
     *    delsh.lt.delhss*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,
     *        n3tm,nbar)
            delsh = 0.0d00
         endif
         write(fu6,607) s,v-vzero,(x(i),i=1,n3)
         dxmag = 0.0d00
         do i = 1,n3
            dxmag = dxmag + dx(i)*dx(i)
         enddo
         dxmag = sqrt(dxmag)
         write(fu6,608) dxmag,(dx(i)/dxmag,i=1,n3)
         do i = 1,n3
            dxold(i) = dx(i)/dxmag
         enddo
         if (delss.ge.delsv*(1.0d00-eps)) then
            swf = s
            if(lgs(38).gt.0) call drdy_wfiles(swf,0)
            if (lgs(39).ne.0) then
               call drdy_potsp(x,vspc,scr1,amass,natom,n3)
               write (fu6,6001) vspc-vzerosp
               if(lgs(38).eq.4) call drdy_wfiles(swf,1)
            endif
            call drdy_rphwrt (7)
            if(lgs(2).ne.0) then
               ierr = 0
               call drdy_projct(ierr)
               if (ierr.eq.0) call drdy_fdiag(n3,nf(5),ierr,lgs(2))
            endif
            delss = 0.0
         endif
      endif
C
      Call drdy_path2(iopt,delss,delsh,sold,smax,xsign)
C
      Return
600   Format(//,' ********** Calculate points along reaction path',//,
     *   5x,' dels=',T20,1pe13.5,/,5x,' delsv=',T20,1pe13.5,/,
     *   5x,' delhss=',T20,1pe13.5,/,5x,' dir=',T20,1pe13.5,/,
     *   5x,' slm=', T20,1pe13.5,/,5x,' slp=',T20,1pe13.5)
601   Format(/,' Path followed using Euler integrator')
602   Format(/,' Path followed using Page-McIver LQA algorithm')
603   Format(/,' Path followed using Page-McIver CLQA algorithm')
604   Format(/,' Path followed using Page-McIver CUBE algorithm')
605   Format(/,' Take step from restart geometry in reactant',
     *   ' direction')
606   Format(/,'  All quantities in atomic units, x is vector of ',
     *   'mass-weighted coordinates,',/,'  dx is the normalized mass-',
     *   'weighted gradient vector, and |dV/x| is the magnitude of',
     *   ' the gradient',/,9x,'s',7x,'V,|dV/dx|')
607   Format (t4,0pf12.5,1pe13.5,T33,'x=',(T35,0p9f9.4))
608   Format (15x,1pe13.5,T32,'dx=',(T35,0p9f9.4))
6001  Format(11x,'Vsp=',1pe13.5)
609   Format(/,' Take step from saddle point in reactant direction')
610   Format(/,' Take step from restart geometry in product',
     *   ' direction')
611   Format(/,' Take step from saddle point in product direction')
612   Format(/,' Take step from initial geometry along gradient')
      End
