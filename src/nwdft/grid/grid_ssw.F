cedo#define USE_VEC 1
      Subroutine grid_ssw(xyzw, xyz, Rij, rq, p, ictr, 
     &                     nctrs_pruned, nq, 
     &                     iscr, indx, whichfunc)
c
C$Id$
c
      implicit none
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
      double precision p(nctrs_pruned)
      double precision rq(nq,*)  ! sum of p(n) [ output]

c
      integer i
      integer  AA,NN
      double precision radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision  toll_scr
c
      integer n,n1,n2
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
      ass=ass_erf1
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
        if(rq(n,AA).ge.radwgh+eps) then
          n1=n
          goto 31
        endif
      enddo
c     all inside
      return
31    continue  
c
c       compute mu_AN
c
      distnn=(ass-eps)*sqrt(distnn)
      n2=nq
      do n=nq,n1,-1
         if ((rq(n,AA)-rq(n,nn)).gt.distnn) then
            xyzw(4,n)=0
         else
            n2=n
            goto 32
         endif
      enddo
32     continue
c
#ifdef USE_OPENMP
!$omp parallel do 
!$omp&  private(n,p)  shared(xyzw)
c1!$omp& default(shared) private(n,p) 
cc!$omp& reduction(+:wsum)
#endif
      do n=n1,n2
         call grid_ssw0(nq,n,nctrs_pruned,
     ,     p,rq,rij,ass)
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
            call errquit(' ssw : wsum .eq.0 ',0,0)
                                !           xyzw(4,n) =0d0
         endif
      enddo
#ifdef USE_OPENMP
!$omp end parallel do
#endif
      call grid_compress(nq,xyzw,iscr,toll_scr)
        
      return
      end
#if 1
      subroutine grid_ssw0(nq,n,nctrs_pruned,
     ,     p,rq,rij,ass)
      implicit none
      integer nq ! complete no. grid pts
      integer n ! grid point
      integer nctrs_pruned ! no. centers
      double precision p(*)
      double precision rq(nq,*)
      double precision rij(*)
      double precision ass ! erfc spread
      integer i,j,ij
      double precision rag,mu,sk
      double precision leps,toll_scr,x
      parameter(leps=1d-6,toll_scr=1d-10)
#include "xc_erftab.fh"
c
      call dcopy(nctrs_pruned,1.d0,0,p,1)
      do  i = 2, nctrs_pruned
         ij = (i*(i-1))/2
         rag=rq(n,i)
         do  j = 1, i-1
            ij=ij+1
               if(
     $           abs(p(i)).ge.toll_scr.or.
     $           abs(p(j)).ge.toll_scr) then
            mu = (rag - rq(n,j))*rij(ij)
            if (mu.gt.ass) then
               p(i)=0.d0
            elseif (mu.lt.-ass) then
               p(j)=0.d0
            else
c     
c     use interpolation for erfs
c
               if(abs(mu).lt.leps) then
                  sk=-1.88063178008d0*mu+0.5d0
               else
                  sk=erf1c(mu)
               endif
               if(mu.gt.0d0) sk = 1d0 - sk
               p(j) = p(j)*sk
               p(i) = p(i)*(1d0 - sk)
            endif
        endif
      enddo                     ! end loop over j
      enddo                     ! end loop over i

      return
      end
#else      
      subroutine grid_ssw0(nq,n,nctrs_pruned,
     ,     p,rq,rij,ass)
      implicit none
      integer nq ! complete no. grid pts
      integer n ! grid point
      integer nctrs_pruned ! no. centers
      double precision p(*)
      double precision rq(nq,*)
      double precision rij(*)
      double precision ass ! erfc spread
      integer i,j,ij
      double precision rag,mu,sk
      double precision leps,toll_scr,x
      parameter(leps=1d-6,toll_scr=1d-10)
#include "xc_erftab.fh"
c
#ifdef USE_VEC
      call dfill(nctrs_pruned,1.d0,p,1)
      do  i = 1, nctrs_pruned
#else
      call dcopy(nctrs_pruned,1.d0,0,p,1)
      do  i = 2, nctrs_pruned
#endif
         ij = (i*(i-1))/2
         rag=rq(n,i)
#ifdef USE_VEC
CDIR$IVDEP
         do  j = 1, nctrs_pruned
            if(i.ne.j) then
               if(j.lt.i) then
                  ij=(i*(i-1))/2+j
               else
                  ij=(j*(j-1))/2+i
               endif
#else         
         do  j = 1, i-1
            ij=ij+1
#endif
               if(
c#ifndef USE_VEC
#if 1
     $           abs(p(i)).ge.toll_scr.or.
#endif
     $           abs(p(j)).ge.toll_scr) then
            mu = (rag - rq(n,j))*rij(ij)
#ifdef USE_VEC
                  if (mu.lt.-ass) then
                     p(j)=0.d0
                  elseif(mu.lt.ass) then
#else
            if (mu.gt.ass) then
               p(i)=0.d0
            elseif (mu.lt.-ass) then
               p(j)=0.d0
            else
#endif
c     
c     use interpolation for erfs
c
               if(abs(mu).lt.leps) then
                  sk=-1.88063178008d0*mu+0.5d0
               else
                  sk=erf1c(mu)
               endif
               if(mu.gt.0d0) sk = 1d0 - sk
               p(j) = p(j)*sk
#ifndef USE_VEC
               p(i) = p(i)*(1d0 - sk)
#endif
            endif
        endif
#ifdef USE_VEC
        endif
#endif
      enddo                     ! end loop over j
      enddo                     ! end loop over i

      return
      end
#endif
