      Subroutine grid_ssw(xyzw, xyz, Rij, rq, p, ictr, 
     &                     nctrs_pruned, nq, 
     &                     iscr, whichfunc)
c
C$Id: grid_ssw.F,v 1.2 1999-07-10 01:06:09 edo Exp $
c
      implicit none
#include "xc_erftab.fh"

c
c      
      integer nctrs_pruned, nq![input]
c
      double precision xyzw(4,nq)! grid and weight points [input]
      double precision xyz(3,*)  ! atom coordinates [input] 
c
      double precision Rij(*) !interatomic distances [input]
      integer ictr
      integer iscr(*)
      character*4 whichfunc ! [input]
      double precision p(*)
      double precision rq(nq,*)  ! sum of p(n) [ output]

c
      integer i, j, n,  ij
      integer  AA,NN,ipoint
      double precision mu, nu,radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision  toll_scr,eps,rag 
      logical inotA,jnota
c
      integer iw
      double precision ass,distnc,distnn,wsum,
     ,ass_erf1,ass_erf2,ass_ssf,ass_erf3

c
c     SSF original 0.64d0 value for polynomial function
c
#if defined(SGI) || defined(DECOSF)
      double precision derfc
#else
      double precision erfc
#endif
      parameter(ass_ssf=0.64d0,ass_erf1=.875d0,ass_erf2=0.77d0,
     ,ass_erf3=0.70d0,
     ,     eps=0.02d0,
     ,     toll_scr=1d-14)
c     
c     RE Stratmann, GE Scuseria, MJ Frisch, Chem Phys Lett 257, 213 (1996)
c     Evaluate Stratman space partitioning weight. Then, incorporate it
c     with weights from the single-center quadratures to form the total 
c     multi-center quadrature weight.
c     
      if(whichfunc.eq.'ssf ') then
        ass=ass_ssf
        iw=0
      elseif(whichfunc.eq.'erf1') then
        ass=ass_erf1
        iw=1
      elseif(whichfunc.eq.'erf2') then
        ass=ass_erf2
        iw=2
      elseif(whichfunc.eq.'erf3') then
        ass=ass_erf3
        iw=3
      endif
c      toll=max((tol_rho*1d-5),1d-33)
      do  i = 1, nctrs_pruned

        x0i = xyz(1,i)
        y0i = xyz(2,i)
        z0i = xyz(3,i)
       
        do  n = 1,nq
          x = xyzw(1,n) - x0i
          y = xyzw(2,n) - y0i
          z = xyzw(3,n) - z0i
         
          rq(n,i) = sqrt(x*x + y*y + z*z) 
        enddo
      enddo
      call a_dist(xyz, Rij, nctrs_pruned)
      AA=ictr
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
      do n=1,nq
c
c       check if grid point is within sphere where w=1
c
        if(rq(n,AA).lt.radwgh+eps) goto 30
c
c       compute mu_AN
c
        mu=(rq(n,AA)-rq(n,nn))/sqrt(distnn)
        if (mu.gt.ass-eps) then
          xyzw(4,n)=0
          goto 30
        endif

        call dfill(nctrs_pruned,1.d0,p,1)
        do  i = 2, nctrs_pruned
          inota=i.ne.AA
          ij = (i*(i-1))/2
          rag=rq(n,i)

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
                  ipoint=int(scale*abs(mu))
                  nu = abs(mu) - table(0,ipoint)
                  sk= (table(1,ipoint) + nu*(table(2,ipoint) + 
     $                 nu*table(3,ipoint)))
                  if(mu.lt.0d0) sk=1d0-sk
#if defined(SGI) || defined(DECOSF)
                elseif(iw.eq.1) then
                  mu=mu/(1d0-mu*mu)
                  mu=mu/(0.34d0)
                  sk=0.5d0*derfc(mu)
                elseif(iw.eq.2) then
                  mu=mu/(1d0-mu*mu)**2
                  mu=mu/(0.34d0)
                  sk=0.5d0*derfc(mu)
                elseif(iw.eq.3) then
                  mu=mu/(1d0-mu*mu)**3
                  mu=mu/(0.34d0)
                  sk=0.5d0*derfc(mu)
                endif
#else
                elseif(iw.eq.1) then
                  mu=mu/(1d0-mu*mu)
                  mu=mu/(0.34d0)
                  sk=0.5d0*erfc(mu)
                elseif(iw.eq.2) then
                  mu=mu/(1d0-mu*mu)**2
                  mu=mu/(0.34d0)
                  sk=0.5d0*erfc(mu)
                elseif(iw.eq.3) then
                  mu=mu/(1d0-mu*mu)**3
                  mu=mu/(0.34d0)
                  sk=0.5d0*erfc(mu)
                endif
#endif
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
        xyzw(4,n) = p(AA)*xyzw(4,n)/wsum
   30   continue
      enddo
      call grid_compress(nq,xyzw,iscr,toll_scr)
        
      return
      end
