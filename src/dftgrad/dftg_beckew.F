      Subroutine dftg_beckew(d_qwght,qwght,qxyz,xyz,Rij,aratio,rq,p,
     &     zeta,dzeta,d_p,ictr,ncenters,iatype,nq,
     &     iandex)

C     $Id: dftg_beckew.F,v 1.1 1995-12-22 22:41:15 og845 Exp $
      implicit none
#include "bas_staticP.fh"
      
      integer ncenters,nq

      double precision qxyz(3,nq),xyz(3,*)

      double precision aratio(ntags_bsmx,*)
      double precision Rij(*)
      integer iatype(*)
      integer ictr(*)
      integer iandex(ncenters)
      double precision p(nq,*),rq(nq,*),qwght(nq),zeta(nq),
     &     dzeta(3,nq,ncenters)
      double precision d_p(3,nq,*)
      double precision d_qwght(3,nq,*)
C     
      integer i,j,n,ij
      integer iind,jind,istep,iti
      double precision mu,nu
      double precision x,y,z,x0i,y0i,z0i
      double precision sk
      double precision aaij,rrij
C     
      double precision damuab(3),dbmuab(3)
      double precision dBPA(3),dAPA(3)
      double precision xi,yi,zi,dmu1,dmu2,dmu3
      integer ictr_reduc
      double precision tmu,tmu1
      integer A,B
      double precision fac,tfac,TOLL
      parameter (tfac=-27./16.)
      parameter(TOLL=0.99d0)
c     
c     Evaluate derivative of Becke space partitioning weight. 
c     
c*******************************************************************************
      
      call dfill(nq*3*ncenters,0.d0,d_p,1)
      call dfill(nq*3*ncenters,0.d0,dzeta,1)
      call dfill(nq*3*ncenters,0.d0,d_qwght,1)
C     if(TOLL.lt.1000.) return
C     
C     see where is ictr in the reduced list
C     
      do i=1,ncenters
        if(iandex(i).eq.ictr(1)) then
          ictr_reduc=i
          goto 123
        endif
      enddo
      call errquit(' bummer ! :',1234)
  123 continue
      do 20 i = 1,ncenters

        x0i=xyz(1,i)
        y0i=xyz(2,i)
        z0i=xyz(3,i)

        do 10 n = 1,nq

          x = qxyz(1,n) - x0i
          y = qxyz(2,n) - y0i
          z = qxyz(3,n) - z0i
          

          rq(n,i) = sqrt(x*x + y*y + z*z) 

          p(n,i) = 1.D0

   10   continue
   20 continue

      do 50 i=1,ncenters
        iind=iandex(i)
        iti = iatype(iind)
        istep=(iind*(iind-1))/2
        do 40 j = 1,i-1

          jind=iandex(j)

          ij=istep+jind
          rrij=Rij(ij)
          
          aaij=aratio(iti,iatype(jind))
          
          do 30 n = 1,nq

            mu = (rq(n,i) - rq(n,j))*rrij

            if(mu.gt.TOLL) then

              p(n,i)=0.d0


            else

c             i>j

c             
c             atomic size adjustment 
c             
C             nu = mu + aaij*(1.D0 - mu*mu)
              nu=mu
C             
C             compute p1(nu)
C             
              nu = 0.5D0*nu*(3.D0 - nu*nu)
C             
C             compute p2(nu)
C             
              nu = 0.5D0*nu*(3.D0 - nu*nu)
C             
C             compute p3(nu)
C             
              nu = 0.5D0*nu*(3.D0 - nu*nu)
              
              sk = 0.5D0 - 0.5D0*nu

              p(n,i)=p(n,i)*sk
              p(n,j)=p(n,j)*(1.d0 - sk)

            endif
   30     continue

   40   continue

   50 continue
C     
C     compute sum of Becke weights for normalization
C     
      call dfill(nq,0.d0,zeta,1)
      do j=1,ncenters
        do  n = 1,nq
          zeta(n) = zeta(n) + p(n,j)
        enddo
      enddo
      do  A= 1,ncenters
        iind=iandex(A)
        iti = iatype(iind)
C       
C       derivation variable B
C       
        do  B=1,ncenters
          if(A.ne.B) then
C        write(0,*) ' A ',A,' B ',B,' istep ',istep
            
            jind=iandex(B)
            if(A.ge.B) then
              ij=(iind*(iind-1))/2+jind
            else
              ij=(jind*(jind-1))/2+iind
            endif

            aaij=aratio(iti,iatype(jind))
            rrij=Rij(ij)
            
            dmu1 = rrij*(xyz(1,A)-xyz(1,B))
            dmu2 = rrij*(xyz(2,A)-xyz(2,B))
            dmu3 = rrij*(xyz(3,A)-xyz(3,B))
            
            do  n = 1,nq

              mu = (rq(n,A) - rq(n,B))*rrij
              if(mu.le.TOLL) then

c               
c               atomic size adjustment 
c               
C               nu = mu + aaij*(1.D0 - mu*mu)
                nu =mu
C
C               t(mu)=-27/16*(1-mu*mu)
C
                tmu=tfac*(1.D0-nu*nu)
C               
C               compute p1(nu)
C               
                nu = 0.5D0*nu*(3.D0 - nu*nu)
C
C               t(mu)=t(mu)*(1-p1(mu)*p1(mu))
C
                tmu=tmu*(1.D0-nu*nu)
C               
C               compute p2(nu)
C               
                nu = 0.5D0*nu*(3.D0 - nu*nu)
C
C               t(mu)=t(mu)*(1-p2(mu)*p2(mu))
C
                tmu=tmu*(1.D0-nu*nu)
C               
C               compute p3(nu)
C               
                nu = 0.5D0*nu*(3.D0 - nu*nu)
                
                sk = 0.5D0 - 0.5D0*nu
                tmu=tmu/sk

C               
C               atomic size adjustment derivative
C               
C               fac=(1.-aaij*2.*mu)
                fac=1.d0
C               
C               compute D(A)mu(AB)
C               

                xi = qxyz(1,n) - xyz(1,A)
                yi = qxyz(2,n) - xyz(2,A)
                zi = qxyz(3,n) - xyz(3,A)
 
                damuab(1) = (xi/rq(n,A) - mu*dmu1)*rrij*fac
                damuab(2) = (yi/rq(n,A) - mu*dmu2)*rrij*fac
                damuab(3) = (zi/rq(n,A) - mu*dmu3)*rrij*fac
C               
C               compute D(B)mu(AB)
C               
                xi = qxyz(1,n) - xyz(1,B)
                yi = qxyz(2,n) - xyz(2,B)
                zi = qxyz(3,n) - xyz(3,B)


                dbmuab(1) = -(xi/rq(n,B) - mu*dmu1)*rrij*fac
                dbmuab(2) = -(yi/rq(n,B) - mu*dmu2)*rrij*fac
                dbmuab(3) = -(zi/rq(n,B) - mu*dmu3)*rrij*fac

                tmu1=tmu*p(n,A)
C               
C               term \Delta_B PA
C               
                dBPA(1)= tmu1*dbmuab(1)
                dBPA(2)= tmu1*dbmuab(2)
                dBPA(3)= tmu1*dbmuab(3)
C               
C               term \Delta_A PA (partial)
C               
                dAPA(1) = tmu1*damuab(1)
                dAPA(2) = tmu1*damuab(2)
                dAPA(3) = tmu1*damuab(3)

                dzeta(1,n,B) = dzeta(1,n,B)+ dBPA(1)
                dzeta(2,n,B) = dzeta(2,n,B)+ dBPA(2)
                dzeta(3,n,B) = dzeta(3,n,B)+ dBPA(3)

                if(A.eq.ictr_reduc) then
                  d_p(1,n,B) = d_p(1,n,B)+ dBPA(1)
                  d_p(2,n,B) = d_p(2,n,B)+ dBPA(2)
                  d_p(3,n,B) = d_p(3,n,B)+ dBPA(3)
                else
                  dzeta(1,n,A) = dzeta(1,n,A)+ dAPA(1)
                  dzeta(2,n,A) = dzeta(2,n,A)+ dAPA(2)
                  dzeta(3,n,A) = dzeta(3,n,A)+ dAPA(3)
                endif

              endif
            enddo
          endif
        enddo
      enddo
      
C     
C     compute \Delta_i W_ictr
C     

C     
C       \Delta_B PA -\delta_B Z*PA/Z
C     
      do B=1,ncenters
        if(B.ne.ictr_reduc) then
          do n=1,nq
            d_qwght(1,n,B)= (d_p(1,n,B) -
     -           dzeta(1,n,B)*p(n,ictr_reduc)/zeta(n))*qwght(n)/zeta(n)
            d_qwght(2,n,B)= (d_p(2,n,B) -
     -           dzeta(2,n,B)*p(n,ictr_reduc)/zeta(n))*qwght(n)/zeta(n)
            d_qwght(3,n,B)= (d_p(3,n,B) -
     -           dzeta(3,n,B)*p(n,ictr_reduc)/zeta(n))*qwght(n)/zeta(n)
          enddo
        endif
      enddo
      do  n = 1,nq
        qwght(n) = (p(n,ictr_reduc)*qwght(n))/zeta(n)
      enddo
                
C        write(0,1235) qxyz(1,n),qxyz(2,n),qxyz(3,n),
C     &               dAPA(1),dAPA(2),dAPA(3),
C     &               dBPA(1),dBPA(2),dBPA(3)
C 1235   format(1x,3f12.6,'A',3f10.6,'B',3f10.6,i4)
C        write(0,1234) B,qxyz(1,n),qxyz(2,n),qxyz(3,n),
C     & d_qwght(1,n,B),d_qwght(2,n,B),d_qwght(3,n,B)
C 1234   format(i4,3f12.6,' | ',3f12.6)

      return
      end
