      Subroutine drdy_path(ns)
C
C  Calculate steepest descent path from saddle point
C     lgs(31)  -  option for following MEP
C              =1, use Euler integrator
C              =2, use Page-McIver LQA algorithm
C              =3, use Page-McIver CLQA algorithm
C              =4, use Page-McIver CUBE algorithm
C
      Implicit Double Precision (A-H,O-Z)
c
#include "drdyP.fh"
C
      Data eps/1.d-6/
C
      iopt = lgs(31)
      nbar = (n3*(n3+1))/2
      write (6,600) dels,delsv,dir,slm,slp
      if (iopt.le.1) then
         write (6,601)
      else if (iopt.eq.2) then
         write (6,602)
      else if (iopt.eq.3) then
         write (6,603)
      else
         write (6,604)
      end if
C
C Set up for calculation on reactant side          
      xsign = -1.0d00
      smax = slm
      if (ns.gt.0 .and. sgrid(1).lt.0.0d00) then
C Restart from grid point with most negative s value
         write (6,608)
         write (6,610) sgrid(1),vgrid(1),(xgrid(i,1),i=1,n3)
         write (6,611) (dxgrid(i,1),i=1,n3)
         ii = 0
         do i = 1,n3
            x(i) = xgrid(i,1)
            dx(i) = dxgrid(i,1)
            do j = 1,i
               ii = ii + 1
               f(j,i) = hgrid(ii,1)
            enddo
         enddo
         do i = 2,n3
            do j = 1,i-1
               f(i,j) = f(j,i)
            enddo
         enddo
         if (iopt.gt.3) then
            if (ns.ge.2) then
               sold = sgrid(2)
               ii = 0
               do i = 1,n3
                  do j = 1,i
                     ii = ii + 1
                     fold(j,i) = hgrid(ii,2)
                  enddo
               enddo
               do i = 2,n3
                  do j = 1,i-1
                     fold(i,j) = fold(j,i)
                  enddo
               enddo
            else
               sold = sgrid(1) + dels
               do i = 1,n3
                  do j = 1,n3
                     fold(i,j) = f(i,j)
                  enddo
               enddo
            endif
         endif
         s = sgrid(1)
         delss = 0.0d00
C
         Call drdy_path2(iopt,delss,sold,smax,xsign)
C
      elseif(lgs(1).ne.0) then
C Calculation for step off saddle point
         write (6,605)
         write (6,610) 0,vspsv-vzero,(xspsv(i),i=1,n3)
         write (6,611) (vec0sv(i),i=1,n3)
C Take step off saddle point in reactant direction
         do i = 1,n3
            x(i) = xspsv(i) - dels*(vec0sv(i) -
     &           0.5d00*dels*vec1sv(i))
            do j = 1,n3
               f(j,i) = fspsv(j,i)
            enddo
         enddo
C Store hessian matrix in fold
         if (iopt.ge.3) then
            sold = 0.0d00
            do i = 1,n3
               do j = 1,n3
                  fold(j,i) = f(j,i)
               enddo
            enddo
         endif
         s = -dels
         delss = dels
C Calculate g and optionally f at new geometry
         if (iopt.eq.1 .and. delss.lt.delsv*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
         endif
         write (6,610) s,v-vzero,(x(i),i=1,n3)
         write (6,611) (dx(i),i=1,n3)
         if (delss.ge.delsv*(1.0-eps)) then
            call drdy_rphwrt (7)
            delss = 0.0d00
         endif
C
         Call drdy_path2(iopt,delss,sold,smax,xsign)
C
      endif
C
C Calculation on product side
      xsign = 1.0d00
      smax = slp
      if (ns.gt.0 .and. sgrid(ns).gt.0.0d00) then
C Restart from grid point with most positive s value
         write (6,609)
         write (6,610) sgrid(ns),vgrid(ns),(xgrid(i,ns),i=1,n3)
         write (6,611) (dxgrid(i,ns),i=1,n3)
         ii = 0
         do i = 1,n3
            x(i) = xgrid(i,ns)
            dx(i) = dxgrid(i,ns)
            do j = 1,i
               ii = ii + 1
               f(j,i) = hgrid(ii,ns)
            enddo
         enddo
         do i = 2,n3
            do j = 1,i-1
               f(i,j) = f(j,i)
            enddo
         enddo
         if (iopt.gt.3) then
            if (ns.ge.2) then
               sold = sgrid(ns-1)
               ii = 0
               do i = 1,n3
                  do j = 1,i
                     ii = ii + 1
                     fold(j,i) = hgrid(ii,ns-1)
                  enddo
               enddo
               do i = 2,n3
                  do j = 1,i-1
                     fold(i,j) = fold(j,i)
                  enddo
               enddo
            else
               sold = sgrid(ns) - dels
               do i = 1,n3
                  do j = 1,n3
                     fold(i,j) = f(i,j)
                  enddo
               enddo
            endif
         endif
         s = sgrid(ns)
         delss = 0.0d00
         Call drdy_path2(iopt,delss,sold,smax,xsign)
C 
      elseif(lgs(1).ne.0) then
C Calculation for step off saddle point
         write (6,607)
         write (6,610) 0,vspsv-vzero,(xspsv(i),i=1,n3)
         write (6,611) (vec0sv(i),i=1,n3)
C Take step off saddle point in reactant direction
         do i = 1,n3
            x(i) = xspsv(i) + dels*(vec0sv(i) +
     &           0.5d0*dels*vec1sv(i))
            do j = 1,n3
               f(j,i) = fspsv(j,i)
            enddo
         enddo
C Store hessian matrix in fold
         if (iopt.ge.3) then
            sold = 0.0d00
            do i = 1,n3
               do j = 1,n3
                  fold(j,i) = f(j,i)
               enddo
            enddo
         endif
         s = dels
         delss = dels
C Calculate g and optionally f at new geometry
         if (iopt.eq.1 .and. delss.lt.delsv*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
         endif
         write (6,610) s,v-vzero,(x(i),i=1,n3)
         write (6,611) (dx(i),i=1,n3)
         if (delss.ge.delsv*(1.0-eps)) then
            call drdy_rphwrt (7)
            delss = 0.0d00
         endif
         Call drdy_path2(iopt,delss,sold,smax,xsign)
      else
C No saddle point, take step from geometry in product valley along MEP
         write (6,606)
         write (6,610) 0,vspsv-vzero,(xspsv(i),i=1,n3)
         write (6,611) (dxspsv(i),i=1,n3)
         if (iopt.eq.1) then
            call drdy_euler (n3tm,n3,xspsv,dxspsv,dels)
         elseif (iopt.ge.2) then
            ndbar = (n3*(n3+1))/2
            call drdy_pmlqa (n3tm,n3,ndbar,xspsv,
     &            dxspsv,fspsv,hess,amass,
     &            dels,vec0,vec1,u0,scr1,vec2,ierr)
         endif
         s = dels
         delss = dels
C Store hessian matrix in fold
         if (iopt.ge.3) then
            sold = 0.0d00
            do i = 1,n3
               do j = 1,n3
                  fold(j,i) = fspsv(j,i)
               enddo
            enddo
         endif
C Calculate g and f at new geometry
         if (iopt.eq.1 .and. delss.lt.delsv*(1.0d00-eps)) then
            Call drdy_pot1 (x,v,dx,scr1,amass,natom,n3)
         else
            Call drdy_pot2 (x,v,dx,f,hess,scr1,amass,natom,n3,n3tm,nbar)
         endif
         write (6,610) s,v-vzero,(x(i),i=1,n3)
         write (6,611) (dx(i),i=1,n3)
         if (delss.ge.delsv*(1.0-eps)) then
            call drdy_rphwrt (7)
            delss = 0.0d00
         endif
         Call drdy_path2(iopt,delss,sold,smax,xsign)
      endif
C
c     Call path2(iopt,delss,sold,smax,xsign)
C
      Return
600   Format(//,' Calculate points along reaction path',//,
     *   5x,' dels=',T20,1pe13.5,/,5x,' delsv=',T20,1pe13.5,/,
     *   5x,' dir=',T20,1pe13.5,/,
     *   5x,' slm=', T20,1pe13.5,/,5x,' slp=',T20,1pe13.5)
601   Format(/,' Path followed using Euler integrator')
602   Format(/,' Path followed using Page-McIver LQA algorithm')
603   Format(/,' Path followed using Page-McIver CLQA algorithm')
604   Format(/,' Path followed using Page-McIver CUBE algorithm')
605   Format(/,' Take step from saddle point in reactant direction',
     *   /,7x,'s',9x,'v',10x,'x/dx=')
606   Format(/,' Take step from initial geometry along gradient',
     *   /,7x,'s',9x,'v',10x,'x/dx=')
607   Format(/,' Take step from saddle point in product direction',
     *   /,7x,'s',9x,'v',10x,'x/dx=')
608   Format(/,' Take step from restart geometry in reactant',
     *   ' direction',/,7x,'s',9x,'v',10x,'x/dx=')
609   Format(/,' Take step from restart geometry in product',
     *   ' direction',/,7x,'s',9x,'v',10x,'x/dx=')
610   Format(1x,0pf10.5,1pe13.5,(t26,1p6e13.5))
611   Format(t26,1p6e13.5)
      End
