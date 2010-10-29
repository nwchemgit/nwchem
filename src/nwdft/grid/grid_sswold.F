      Subroutine grid_sswold(xyzw, xyz, Rij, rq, p, ictr, 
     &                     nctrs_pruned, nq, 
     &                     iscr, indx, whichfunc)
c
C$Id$
c
      implicit none
#include "errquit.fh"
#include "cdft.fh"

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
      integer indx(*)
      character*8 whichfunc ! [input]
      double precision p(*)
      double precision rq(nq,*)  ! sum of p(n) [ output]

c
      integer i, j,   ij
      integer  AA,NN
      double precision mu, nu,radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision  toll_scr,rag 
      logical inotA,jnota
c
      integer iw,n
      double precision ass,distnc,distnn,wsum
      parameter (toll_scr=1d-13)
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
      if(whichfunc.eq.'ssf ') then
        ass=ass_ssf
        iw=0
      elseif(whichfunc.eq.'erf1') then
        ass=ass_erf1
        iw=1
      elseif(whichfunc.eq.'erf2') then
        ass=ass_erf2
        iw=2
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
      call a_dist(xyz, Rij, nctrs_pruned,.false.)
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
        if(wsum.ne.0) then
           xyzw(4,n) = p(AA)*xyzw(4,n)/wsum
        else
           call errquit(' ssw : wsum .eq.0 ',0, UNKNOWN_ERR)
!           xyzw(4,n) =0d0
        endif
   30   continue
      enddo
      call grid_compress(nq,xyzw,iscr,toll_scr)
        
      return
      end
