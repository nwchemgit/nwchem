      Subroutine dftg_gridssw(grid_written,
     ,     d_qwght,qwght, qxyz, xyz, Rij, rq, p, 
     ,     zeta, dzeta, d_p,
     .     ictr, nctrs_pruned, nq_orig,nq, 
     ,     lscreen)
c
C$Id$
c
      implicit none
#include "errquit.fh"
      logical grid_written
      integer nctrs_pruned ! [in] natoms after signf
      integer nq           ! [in] no. grid pts
      integer nq_orig      ! [in/out] no. grid pts after compression
      logical lscreen ! [in] screen weights
c
      double precision qxyz(3,nq)! grid points [input]
      double precision xyz(3,*)  ! atom coordinates [input] 
c
      double precision Rij(*) !interatomic distances [input]
      integer ictr ! [in] ctr of quadr
      double precision p(nctrs_pruned)
      double precision rq(nq_orig,*)  ! sum of p(n) [ output]
      double precision qwght(nq_orig) ! weights [output]
      double precision d_qwght(3,nq_orig,*) ! weight deriv [output]
c
      double precision d_p(3,nq_orig,*), zeta(nq_orig),
     &                 dzeta(3,nq_orig,*)
      integer i, j,  ij
      integer iind,jind,   AA,NN
      double precision mu, radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision  toll,rag 
      parameter(toll=1d-13)

      logical inotA,jnota
c
      integer n
      double precision ass,distnc,distnn,wsum,dist2
      double precision damuab(3), dbmuba(3)
      double precision dBPA(3)
      double precision xi, yi, zi, dmu1, dmu2, dmu3
      double precision tmu, tmu1,mu1,dmu1dmu,dskdmu1,asqrtpi
      double precision xa,ya,za,alpha_erf,murij
      integer A,B,n1,n2
c 20 digits ln(10)=2.3025
c 0.5d0 factor because of mu^2??
      double precision undovl
      parameter(undovl=20d0*2.3025d0*0.5d0)
      external fastexp
      double precision fastexp,asse
#include "xc_erftab.fh"
c     
c     RE Stratmann, GE Scuseria, MJ Frisch, Chem Phys Lett 257, 213 (1996)
c     Evaluate Stratman space partitioning weight. Then, incorporate it
c     with weights from the single-center quadratures to form the total 
c     multi-center quadrature weight.
c     The following 2 lines are to satisfy compiler warnings.
c     
      NN = 1
      asqrtpi=1d0/sqrt(4*datan(1d0))
      call dfill(nq*3*nctrs_pruned, 0.d0, dzeta, 1)
      call dfill(nq, 0.d0, zeta, 1)
      call dfill(nq*3*nctrs_pruned, 0.d0, d_p, 1)
      call dfill(nq*3*nctrs_pruned, 0.d0, d_qwght, 1)
        ass=ass_erf1
        asse=ass-eps
        alpha_erf=alpha_erf1
      do  i = 1, nctrs_pruned

        x0i = xyz(1,i)
        y0i = xyz(2,i)
        z0i = xyz(3,i)
       
        do  n = 1,nq
          x = qxyz(1,n) - x0i
          y = qxyz(2,n) - y0i
          z = qxyz(3,n) - z0i
         
          rq(n,i) = sqrt(x*x + y*y + z*z) 
        enddo
      enddo
      call a_dist(xyz, Rij, nctrs_pruned,.false.)
      AA=ictr
c
c     find nearest neighb
c
      distnn=1.d+10
      x0i=xyz(1,AA)
      y0i=xyz(2,AA)
      z0i=xyz(3,AA)
      do i=1,nctrs_pruned
        if(i.ne.AA) then
          distnc=(xyz(1,i)-x0i)*(xyz(1,i)-x0i)+
     +         (xyz(2,i)-y0i)*(xyz(2,i)-y0i)+
     +         (xyz(3,i)-z0i)*(xyz(3,i)-z0i)
          if(distnc.lt.distnn) then
            distnn=distnc
            NN=i
          endif
        endif
      enddo
          
      
      radwgh=(1.d0-ass)*sqrt(distnn)*0.5d0
      do n = 1,nq
c
c       check if grid point is within sphere where w=1
c
         if(rq(n,AA).ge.radwgh+eps) then
            n1=n
            goto 31
         endif
      enddo
c     all inside
      return
 31   continue

      dist2=asse*sqrt(distnn)*4d0
      n2=nq
      do n=nq,n1,-1
         if ((rq(n,AA)-rq(n,nn)).gt.dist2) then
            qwght(n)=0d0
         else
            n2=n
            goto 32
         endif
      enddo
32     continue
#ifdef USE_OPENMP
!$omp parallel do private(n,p)
#endif
      do n=n1,n2
c
c       compute mu_AN
c
        mu=(rq(n,AA)-rq(n,nn))/sqrt(distnn)
        if (mu.gt.asse) then
          p(AA)=0d0
          zeta(n)=1d0
          goto 1100
        endif

        call dfill(nctrs_pruned,1.d0,p,1)
        do  i = 2, nctrs_pruned
          inota=i.ne.AA
          rag=rq(n,i)
          ij = (i*(i-1))/2
          do  j = 1, i-1

            jnota=j.ne.AA
c           
            ij=ij+1
              mu = (rag - rq(n,j))*Rij(ij)
              if (mu.ge.asse) then
                p(i)=0.d0

              elseif (mu.le.-asse) then
                p(j)=0.d0

              else
                 if(inota.and.jnota) then
c                 
c                 use interpolation for erfs
c
                    sk=erf1c(mu)
                    if(mu.lt.0d0) sk=1d0-sk
                 else
                    sk=erf1(mu)
                 endif
                 p(i) = p(i)*sk
                 p(j) = p(j)*(1d0 - sk)
              endif
            enddo ! end loop over j
          enddo   ! end loop over i
c       
c       compute sum of partitioning weights for normalization
c       
c       
        wsum=0.d0
        do i = 1, nctrs_pruned
          wsum=wsum+p(i)
        enddo
        if(abs(wsum).lt.toll) goto 300
        zeta(n) = 1d0/wsum
 1100   continue
      do A = 1, nctrs_pruned
         if(abs(p(A)).gt.toll) then
         iind=A
          inota=A.ne.AA
          xA = (qxyz(1,n) - xyz(1,A))/rq(n,A)
          yA = (qxyz(2,n) - xyz(2,A))/rq(n,A)
          zA = (qxyz(3,n) - xyz(3,A))/rq(n,A)
c       
c        derivation variable B
c       
         do B = 1, nctrs_pruned
          jnota=B.ne.AA
            if (A.ne.B)then
               jind = B
c
               if (A.ge.B)then
                  ij = (A*(A-1))/2 + B
               else
                  ij = (B*(B-1))/2 + A
               endif
c
c            
                  mu = (rq(n,A) - rq(n,B))*Rij(ij)
                  if(abs(mu).lt.asse) then
                     mu1=mu/(1d0-mu*mu)*alpha_erf
                     sk=erf1c(mu)
                     if(mu.lt.0d0) sk=1d0-sk
                     if(mu1.gt.undovl) then
                        dskdmu1=0d0
                        tmu=0d0
                     else
                        dmu1 = Rij(ij)*(xyz(1,A)-xyz(1,B))
                        dmu2 = Rij(ij)*(xyz(2,A)-xyz(2,B))
                        dmu3 = Rij(ij)*(xyz(3,A)-xyz(3,B))
                        dmu1dmu=((mu*mu+1d0)/(1d0-mu*mu)**2)*alpha_erf
                        dskdmu1=-exp(-mu1*mu1)*asqrtpi
                        tmu = dskdmu1*dmu1dmu
                     endif
                  else
                     tmu=0d0
                     sk=0d0
                  endif
c               
                
                  if(abs(sk).gt.toll.and.abs(tmu).gt.0d0) then
c               
c                    compute D(B)mu(AB)
c               
                     xi = qxyz(1,n) - xyz(1,B)
                     yi = qxyz(2,n) - xyz(2,B)
                     zi = qxyz(3,n) - xyz(3,B)
c               
c                    atomic size adjustment derivative
c               
                     dbmuba(1) = -(xi/rq(n,B) + mu*dmu1)*Rij(ij)
                     dbmuba(2) = -(yi/rq(n,B) + mu*dmu2)*Rij(ij)
                     dbmuba(3) = -(zi/rq(n,B) + mu*dmu3)*Rij(ij)
c
                     tmu1=tmu*p(A)/sk
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

                     if (inota)then
c               
c                    term \Delta_A PA (partial)
c               
c               
c                    compute D(A)mu(AB)
c               
                        damuab(1) = -(xA+mu*dmu1)*Rij(ij)
                        damuab(2) = -(yA+mu*dmu2)*Rij(ij)
                        damuab(3) = -(zA+mu*dmu3)*Rij(ij)
                        dzeta(1,n,A) = dzeta(1,n,A)+tmu1*damuab(1)
                        dzeta(2,n,A) = dzeta(2,n,A)+tmu1*damuab(2)
                        dzeta(3,n,A) = dzeta(3,n,A)+tmu1*damuab(3)
                     else
                        d_p(1,n,B) =  dBPA(1)
                        d_p(2,n,B) =  dBPA(2)
                        d_p(3,n,B) =  dBPA(3)
                     endif
                 endif
             endif
           enddo ! B loop
          endif
         enddo   ! A loop
 300     continue
      if(.not.grid_written) then
            if(abs(p(AA)).gt.toll)then
               qwght(n) = (p(AA)*qwght(n))*zeta(n)
            else
               qwght(n)=0d0
            endif
      endif
c     
c     compute \Delta_i W_ictr
c     
c     \Delta_B PA -\delta_B Z*PA/Z
c     
      do B = 1, nctrs_pruned
         if (B.ne.AA)then
               if(abs(p(AA)).gt.toll) then
               d_qwght(1,n,b)= (d_p(1,n,B)/p(AA) -
     &           dzeta(1,n,B)*zeta(n))*qwght(n)
               d_qwght(2,n,b)= (d_p(2,n,B)/p(AA) -
     &           dzeta(2,n,B)*zeta(n))*qwght(n)
               d_qwght(3,n,b)= (d_p(3,n,B)/p(AA) -
     &           dzeta(3,n,B)*zeta(n))*qwght(n)
              endif
         endif
      enddo
      enddo ! n loop
#ifdef USE_OPENMP
!$omp end parallel do
#endif
      return
      end

