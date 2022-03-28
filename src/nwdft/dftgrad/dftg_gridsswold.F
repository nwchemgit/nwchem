      Subroutine dftg_gridsswold(grid_written,
     ,     d_qwght,qwght, qxyz, xyz, Rij, rq, p, 
     ,     dzeta, d_p,
     .     ictr, nctrs_pruned, nq_orig,nq, 
     ,     whichfunc,lscreen)
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
      character*8 whichfunc ! [input]
      double precision p(*)
      double precision rq(nq_orig,*)  ! sum of p(n) [ output]
      double precision qwght(nq_orig) ! weights [output]
      double precision d_qwght(3,nq_orig,*) ! weight deriv [output]
c
      double precision d_p(3,*), 
     &                 dzeta(3,*)
      integer i, j,  ij
      integer iind,jind,   AA,NN
      double precision mu, nu,radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision  toll,rag 
      parameter(toll=1d-13)

      logical inotA,jnota
c
      integer iw,n
      double precision ass,distnc,distnn,wsum
      double precision damuab(3), dbmuba(3)
      double precision dBPA(3)
      double precision xi, yi, zi, dmu1, dmu2, dmu3
      double precision tmu, tmu1,mu1,mu14,dmu1dmu,dskdmu1,asqrtpi
      double precision xa,ya,za,alpha_erf
      double precision  zetan
      integer A,B,n1
#include "xc_erftab.fh"
c     
c     RE Stratmann, GE Scuseria, MJ Frisch, Chem Phys Lett 257, 213 (1996)
c     Evaluate Stratman space partitioning weight. Then, incorporate it
c     with weights from the single-center quadratures to form the total 
c     multi-center quadrature weight.
c     The following 2 lines are to satisfy compiler warnings.
c     
      NN = 1
      iw = 0
      asqrtpi=1d0/sqrt(4*datan(1d0))
      call dfill(nq*3*nctrs_pruned, 0.d0, d_qwght, 1)
      if(whichfunc.eq.'erf1') then
        ass=ass_erf1
        alpha_erf=alpha_erf1
        iw=1
      elseif(whichfunc.eq.'ssf ') then
        ass=ass_ssf
        iw=0
      elseif(whichfunc.eq.'erf2') then
        ass=ass_erf2
        iw=2
        alpha_erf=alpha_erf2
      else
        ass=0d0
        alpha_erf=0d0
        call errquit(' dftg_gridssw: wrong family',0, UNKNOWN_ERR)
      endif
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
      do n=n1,nq
         call dfill(3*nctrs_pruned, 0.d0, dzeta, 1)
         call dfill(3*nctrs_pruned, 0.d0, d_p, 1)
         zetan=0d0
c
c       compute mu_AN
c
        mu=(rq(n,AA)-rq(n,nn))/sqrt(distnn)
        if (mu.gt.ass-eps) then
          p(AA)=0d0
          zetan=1d0
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
              if (mu.ge.(ass-eps)) then
                p(i)=0.d0

              elseif (mu.le.-(ass-eps)) then
                p(j)=0.d0

              else
                if(iw.eq.0) then
                  nu=mu/ass
                  mu=nu*nu*nu*nu
                  sk=0.5d0-(7.d0*(5.d0+3.d0*mu)-
     -                 5.d0*nu*nu*(7.d0+mu))*nu/32.d0
                  
                elseif(inota.and.jnota) then
c                 
c                 use interpolation for erfs
c
                     if(iw.eq.1) then
                        sk=erf1c(mu)
                     else
                        sk=erf2c(mu)
                     endif
                  if(mu.lt.0d0) sk=1d0-sk
               else
                  if(iw.eq.1) then
                     sk=erf1(mu)
                  else
                     sk=erf2(mu)
                  endif
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
        zetan = 1d0/wsum
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
               dmu1 = Rij(ij)*(xyz(1,A)-xyz(1,B))
               dmu2 = Rij(ij)*(xyz(2,A)-xyz(2,B))
               dmu3 = Rij(ij)*(xyz(3,A)-xyz(3,B))
c            
                  mu = (rq(n,A) - rq(n,B))*Rij(ij)
                  if(abs(mu).lt.ass-eps) then
                     if(iw.eq.0) then
                        mu1=mu/ass
                        mu14=mu1*mu1*mu1*mu1
                        sk=0.5d0-(7.d0*(5.d0+3.d0*mu14)-
     -                       5.d0*mu1*mu1*(7.d0+mu14))*mu1/32.d0
                        dmu1dmu=1d0/ass
                        dskdmu1=-((1.d0+3.d0*mu14)-
     -                       mu1*mu1*(3.d0+mu14))*35d0/32.d0
                     else
                        if(iw.eq.1) then
                           mu1=mu/(1d0-mu*mu)
                           dmu1dmu=(mu*mu+1d0)/(1d0-mu*mu)**2
                        elseif(iw.eq.2) then
                           mu1=mu/(1d0-mu*mu)**2
                           dmu1dmu=(3d0*mu*mu+1d0)/(1d0-mu*mu)**3
                        endif
                        mu1=mu1*alpha_erf
                        dmu1dmu=dmu1dmu*alpha_erf
                        if(inota.and.jnota) then
                           if(iw.eq.1) then
                              sk=erf1c(mu)
                           else
                              sk=erf2c(mu)
                           endif
                           if(mu.lt.0d0) sk=1d0-sk
                        else
#if defined (WIN32) || defined(LINUX) || defined(CRAY)
! DERFC is a GNU extension
                           sk=0.5d0*derfc(mu1)
#else
! ERFC is Fortran 2008
                           sk=0.5d0*erfc(mu1)
#endif
                        endif
                        dskdmu1=-exp(-mu1*mu1)*asqrtpi
                     endif
                     tmu = dskdmu1*dmu1dmu
                  else
                     tmu=0d0
                     sk=0d0
                  endif
c               
                
                     if(abs(sk).gt.toll) then
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
                     dzeta(1,B) = dzeta(1,B)+ dBPA(1)
                     dzeta(2,B) = dzeta(2,B)+ dBPA(2)
                     dzeta(3,B) = dzeta(3,B)+ dBPA(3)
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
                     dzeta(1,A) = dzeta(1,A)+tmu1*damuab(1)
                     dzeta(2,A) = dzeta(2,A)+tmu1*damuab(2)
                     dzeta(3,A) = dzeta(3,A)+tmu1*damuab(3)
                   else
                     d_p(1,B) =  dBPA(1)
                     d_p(2,B) =  dBPA(2)
                     d_p(3,B) =  dBPA(3)
                   endif
                 endif
             endif
           enddo ! B loop
          endif
         enddo   ! A loop
 300     continue
      if(.not.grid_written) then
!         do n=1,nq
            if(abs(p(AA)).gt.toll)then
               qwght(n) = (p(AA)*qwght(n))*zetan
            else
               qwght(n)=0d0
            endif
!         enddo
      endif
c     
c     compute \Delta_i W_ictr
c     
c     \Delta_B PA -\delta_B Z*PA/Z
c     
      do B = 1, nctrs_pruned
         if (B.ne.AA)then
!            do n = 1, nq
               if(abs(p(AA)).gt.toll) then
               d_qwght(1,n,b)= (d_p(1,B)/p(AA) -
     &           dzeta(1,B)*zetan)*qwght(n)
               d_qwght(2,n,b)= (d_p(2,B)/p(AA) -
     &           dzeta(2,B)*zetan)*qwght(n)
               d_qwght(3,n,b)= (d_p(3,B)/p(AA) -
     &           dzeta(3,B)*zetan)*qwght(n)
              endif
!            enddo
         endif
      enddo
      enddo ! n loop
      n=0
      do j=1,nq
         if(abs(qwght(j)).lt.toll.and.lscreen) then
            n=n+1
         else
            if(n.ne.0) then
               qwght(j-n)=qwght(j)
               qxyz(1,j-n)=qxyz(1,j)
               qxyz(2,j-n)=qxyz(2,j)
               qxyz(3,j-n)=qxyz(3,j)
               do i=1,nctrs_pruned
                  d_qwght(1,j-n,i)=d_qwght(1,j,i)
                  d_qwght(2,j-n,i)=d_qwght(2,j,i)
                  d_qwght(3,j-n,i)=d_qwght(3,j,i)
               enddo
            endif
         endif
      enddo
      return
      end
