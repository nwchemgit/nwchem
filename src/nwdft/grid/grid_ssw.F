      Subroutine grid_ssw(xyzw, xyz, Rij, rq, p, ictr, 
     &                     nctrs_pruned, nq, 
     &                     iscr, whichfunc,indx)
c
C$Id: grid_ssw.F,v 1.19 2003-10-17 22:55:37 carlfahl Exp $
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
      character*4 whichfunc ! [input]
      double precision p(*)
      double precision rq(nq,*)  ! sum of p(n) [ output]

c
      integer i, j,   ij
      integer  AA,NN
      double precision mu, radwgh
      double precision x, y, z, x0i, y0i, z0i
      double precision sk
      double precision  toll_scr,rag 
      logical inotA
c
      integer iw,n,n1,n2
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
        do n=n1,n2
           call dcopy(nctrs_pruned,1.d0,0,p,1)
           do  i = 2, nctrs_pruned
              inota=i.ne.AA
              ij = (i*(i-1))/2
              rag=rq(n,i)
              
              do  j = 1, i-1
                 ij=ij+1
                 mu = (rag - rq(n,j))*Rij(ij)
                 if (mu.gt.(ass-eps)) then
                    p(i)=0.d0
                 elseif (mu.lt.-(ass-eps)) then
                    p(j)=0.d0
                 else
#if !(defined(IFCLINUX) && defined(LINUX64))
                    if(inota.and. j.ne.AA) then
c     
c     use interpolation for erfs
c
                       sk=erf1c(mu)
                       if(mu.lt.0d0) then
                          p(j) = p(j)*sk
                          p(i) = p(i)*(1d0 - sk)
                       else
                          p(i) = p(i)*sk
                          p(j) = p(j)*(1d0 - sk)
                       endif
                    else
#endif
                       sk=erf1(mu)
                       p(i) = p(i)*sk
                       p(j) = p(j)*(1d0 - sk)
#if !(defined(IFCLINUX) && defined(LINUX64))
                    endif
#endif
                 endif
              enddo             ! end loop over j
           enddo                ! end loop over i

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
        enddo
      call grid_compress(nq,xyzw,iscr,toll_scr)
        
      return
      end
