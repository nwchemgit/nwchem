      Subroutine dftg_beckew(grid_written,
     ,     d_qwght, qwght, qxyz, xyz, Rij, rq, p,
     &                       zeta, dzeta, d_p, ictr, bcenters, 
     .     nq_orig,nq,
     &                       iandex)
c
C     $Id$
c
      implicit none
#include "errquit.fh"
c
#include "cdft.fh"
c      
      logical grid_written ! [in]
      integer bcenters, nq,nq_orig
      double precision qxyz(3,nq), xyz(3,*)
      double precision Rij(*)
      integer ictr(*)
      integer iandex(bcenters)
      double precision p(nq_orig,*), rq(nq_orig,*), 
     ,     qwght(nq_orig), zeta(nq_orig),
     &                 dzeta(3,nq_orig,bcenters)
      double precision d_p(3,nq_orig,*)
      double precision d_qwght(3,nq_orig,*)
C     
      integer i, j, n, ij
      integer iind, jind, istep, iti
      double precision mu, nu
      double precision x,  y, z, x0i, y0i, z0i
      double precision sk
      double precision aaij, rrij
C     
      double precision damuab(3), dbmuba(3)
      double precision dBPA(3), dAPA(3)
      double precision xi, yi, zi, dmu1, dmu2, dmu3
      integer ictr_reduc
      double precision tmu, tmu1
      integer A, B
      double precision fac, tfac, TOLL,eps
      parameter (tfac=-27.d0/16.d0,eps=1.d-16)
      parameter(TOLL=0.925d0)
      double precision ppp
      ppp(x)=0.5D0*x*(3.D0-x*x)
c     
c     Evaluate derivative of Becke space partitioning weight. 
c     
      call dfill(nq*3*bcenters, 0.d0, d_p, 1)
      call dfill(nq*3*bcenters, 0.d0, dzeta, 1)
      call dfill(nq*3*bcenters, 0.d0, d_qwght, 1)
C     if (TOLL.lt.1000.) return
C     
C     see where ictr is in the reduced list
C     
      ictr_reduc = 1        ! take care of compiler warnings
      do i = 1, bcenters
         if (iandex(i).eq.ictr(1))then
            ictr_reduc=i
            goto 123
         endif
      enddo
      call errquit('dftg_beckew: unknown center', 1234, UNKNOWN_ERR)
  123 continue
      do 20 i = 1, bcenters
c
         x0i=xyz(1, i)
         y0i=xyz(2, i)
         z0i=xyz(3, i)
c
         do 10 n = 1, nq
c
            x = qxyz(1,n) - x0i
            y = qxyz(2,n) - y0i
            z = qxyz(3,n) - z0i
c          
            rq(n,i) = sqrt(x*x + y*y + z*z) 
c
            p(n,i) = 1.D0
c
   10    continue
   20 continue
c
      do 50 i = 1, bcenters
         iind = iandex(i)
         iti = iatype(iind)
         istep = (iind*(iind-1))/2
         do 40 j = 1, i-1
            jind = iandex(j)
            ij = istep + jind
            rrij = Rij(ij)
c          
            aaij=aratio(iti,iatype(jind))
c          
            do 30 n = 1, nq
c 
               mu = (rq(n,i) - rq(n,j))*rrij
c
               if (mu.gt.TOLL)then
                  p(n,i)=0.d0
               else
c
c                 i>j
c
c                 atomic size adjustment 
c             
                 nu = mu + aaij*(1.D0 - mu*mu)
                  sk = 0.5D0 - 0.5D0*ppp(ppp(ppp(nu)))
c
                  p(n,i) = p(n,i)*sk
                  p(n,j)=p(n,j)*(1.d0 - sk)
c
               endif
   30       continue
   40    continue
   50 continue
c     
c     compute sum of Becke weights for normalization
c     
      call dfill(nq, 0.d0, zeta, 1)
      do j = 1, bcenters
         do n = 1, nq
            zeta(n) = zeta(n) + p(n,j)
         enddo
      enddo
c
      do A = 1, bcenters
         iind=iandex(A)
         iti = iatype(iind)
c       
c        derivation variable B
c       
         do B = 1, bcenters
            if (A.ne.B)then
               jind = iandex(B)
c
               if (A.ge.B)then
                  ij = (iind*(iind-1))/2 + jind
               else
                  ij = (jind*(jind-1))/2 + iind
               endif
c
               aaij = aratio(iti,iatype(jind))
               rrij = Rij(ij)
c
               dmu1 = rrij*(xyz(1,A)-xyz(1,B))
               dmu2 = rrij*(xyz(2,A)-xyz(2,B))
               dmu3 = rrij*(xyz(3,A)-xyz(3,B))
c            
               do n = 1, nq
c
                  mu = (rq(n,A) - rq(n,B))*rrij
c
c                    atomic size adjustment 
c               
                    nu = mu + aaij*(1.D0 - mu*mu)
c
c                    t(mu) = -27/16*(1-mu*mu)
c
                     tmu = tfac*(1.D0-nu*nu)
c               
c                    compute p1(nu)
c               
                     nu = ppp(nu)
c
c                    t(mu) = t(mu)*(1-p1(mu)*p1(mu))
c
                     tmu = tmu*(1.D0-nu*nu)
c               
c                    compute p2(nu)
c               
                     nu = ppp(nu)
c
c                    t(mu) = t(mu)*(1-p2(mu)*p2(mu))
c
                     tmu = tmu*(1.D0-nu*nu)
c               
c                    compute p3(nu)
c               
                     nu = ppp(nu)
c                
                     sk = 0.5D0 - 0.5D0*nu
c                     tmu=(1.d0-ppp(ppp(nu))**2)*(1.-ppp(nu)**2)*
c     *                    (1.-nu*nu)*tfac  
c                     sk = 0.5D0 - 0.5D0*ppp(ppp(ppp(nu)))
                     if(abs(sk).gt.eps) then
                     tmu = tmu/sk
c               
c                    compute D(B)mu(AB)
c               
                     xi = qxyz(1,n) - xyz(1,B)
                     yi = qxyz(2,n) - xyz(2,B)
                     zi = qxyz(3,n) - xyz(3,B)
c               
c                    atomic size adjustment derivative
c               
                    fac=(1.-aaij*2.d0*mu)
                     dbmuba(1) = -(xi/rq(n,B) + mu*dmu1)*rrij*fac
                     dbmuba(2) = -(yi/rq(n,B) + mu*dmu2)*rrij*fac
                     dbmuba(3) = -(zi/rq(n,B) + mu*dmu3)*rrij*fac
c
                     tmu1=tmu*p(n,A)
c               
c                    term \Delta_B PA
c               
                     dBPA(1)= -tmu1*dbmuba(1)
                     dBPA(2)= -tmu1*dbmuba(2)
                     dBPA(3)= -tmu1*dbmuba(3)
c
                     dzeta(1,n,B) = dzeta(1,n,B)+ dBPA(1)
                     dzeta(2,n,B) = dzeta(2,n,B)+ dBPA(2)
                     dzeta(3,n,B) = dzeta(3,n,B)+ dBPA(3)
c
                     if (A.ne.ictr_reduc)then
c               
c                    term \Delta_A PA (partial)
c               
c               
c                    compute D(A)mu(AB)
c               
                     xi = qxyz(1,n) - xyz(1,A)
                     yi = qxyz(2,n) - xyz(2,A)
                     zi = qxyz(3,n) - xyz(3,A)
                     fac=(1.-aaij*2.d0*mu)
                     damuab(1) = -(xi/rq(n,A)+mu*dmu1)*rrij*fac
                     damuab(2) = -(yi/rq(n,A)+mu*dmu2)*rrij*fac
                     damuab(3) = -(zi/rq(n,A)+mu*dmu3)*rrij*fac
                     dAPA(1) = tmu1*damuab(1)
                     dAPA(2) = tmu1*damuab(2)
                     dAPA(3) = tmu1*damuab(3)
                     dzeta(1,n,A) = dzeta(1,n,A)+ dAPA(1)
                     dzeta(2,n,A) = dzeta(2,n,A)+ dAPA(2)
                     dzeta(3,n,A) = dzeta(3,n,A)+ dAPA(3)
                   else
                     d_p(1,n,B) =  dBPA(1)
                     d_p(2,n,B) =  dBPA(2)
                     d_p(3,n,B) =  dBPA(3)
                   endif
                 endif
               enddo
             endif
           enddo
         enddo
      if(.not.grid_written) then
         do n=1,nq
            if(p(n,ictr_reduc).gt.eps)then
               qwght(n) = (p(n,ictr_reduc)*qwght(n))/zeta(n)
            else
               qwght(n)=0d0
            endif
         enddo
      endif

c     
c     compute \Delta_i W_ictr
c     
c     \Delta_B PA -\delta_B Z*PA/Z
c     
      do B = 1, bcenters
         if (B.ne.ictr_reduc)then
            do n = 1, nq
               if(p(n,ictr_reduc).gt.1d-16) then
                  d_qwght(1,n,B)= (d_p(1,n,B)/p(n,ictr_reduc) -
     &                 dzeta(1,n,B)/zeta(n))*qwght(n)
                  d_qwght(2,n,B)= (d_p(2,n,B)/p(n,ictr_reduc) -
     &                 dzeta(2,n,B)/zeta(n))*qwght(n)
                  d_qwght(3,n,B)= (d_p(3,n,B)/p(n,ictr_reduc) -
     &                 dzeta(3,n,B)/zeta(n))*qwght(n)
               endif
            enddo
         endif
      enddo
      n=1
 111  if(n.le.nq) then
      if(abs(qwght(n)).lt.1d-16) then
         do j=n,nq
            qwght(j)=qwght(j+1)
            qxyz(1,j)=qxyz(1,j+1)
            qxyz(2,j)=qxyz(2,j+1)
            qxyz(3,j)=qxyz(3,j+1)
         enddo
         do i=1,bcenters
            do j=n,nq
               d_qwght(1,j,i)=d_qwght(1,j+1,i)
               d_qwght(2,j,i)=d_qwght(2,j+1,i)
               d_qwght(3,j,i)=d_qwght(3,j+1,i)
            enddo
         enddo
         nq=nq-1
      else
         n=n+1
      endif
      goto 111
      endif
      return
      end
