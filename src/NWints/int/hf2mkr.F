      Subroutine hf2mkr(Axyz,Bxyz,Cxyz,Dxyz,alpha_p,alpha_q,R0,IJK,
     &    Pxyz,Qxyz,PQ,ff,R,ffscr,rscr,NPP,NPQ,Lr,Lr3)
c
c $Id$

      implicit none
#include "case.fh"
c
      double precision PI,P1
      Parameter (PI=3.1415926535898D0,P1=4.D0/PI)

c--> ints
      integer npp,npq
      integer Lr, Lr3
      integer IJK(0:Lr,0:Lr,0:Lr)

c--> Cartesian Coordinates
      double precision Axyz(3),Bxyz(3),Cxyz(3),Dxyz(3)

c--> Exponents

      double precision alpha_p(2,NPP),alpha_q(2,NPQ)

c--> Auxiliary Function Integrals & Index

      double precision R0((NPP*NPQ),Lr3)

c--> Scratch Space

      double precision Pxyz(3,NPP),Qxyz(3,NPQ),PQ((NPP*NPQ),3),
     &  ff(2,(NPP*NPQ)), R((NPP*NPQ),0:Lr,Lr3)
c
      double precision ffscr(2,(NPP*NPQ)), Rscr((NPP*NPQ),0:Lr,Lr3),
     & rhoscaled,radd

      integer mp, mq, mr, j, n
      double precision a,b,c,d,f1,f2,p,q,rho
      double precision PQx, PQy, PQz 
c
c Evaluate the auxiliary function integrals. These integrals are scaled by a 
c factor appropriate for ERIs, defined as
c
c                          1/2
c  ES = ((4/PI)*(pq/(p+q)))           where p = a + b and q = c + d.
c
c******************************************************************************
c
c Define the center "P".

      do 100 mp = 1,NPP

       a = alpha_p(1,mp)
       b = alpha_p(2,mp)

       f1 = a/(a+b)
       f2 = b/(a+b)

       Pxyz(1,mp) = f1*Axyz(1) + f2*Bxyz(1)
       Pxyz(2,mp) = f1*Axyz(2) + f2*Bxyz(2)
       Pxyz(3,mp) = f1*Axyz(3) + f2*Bxyz(3)

  100 continue

c Define the center "Q".

      do 110 mq = 1,NPQ

       c = alpha_q(1,mq)
       d = alpha_q(2,mq)

       f1 = c/(c+d)
       f2 = d/(c+d)

       Qxyz(1,mq) = f1*Cxyz(1) + f2*Dxyz(1)
       Qxyz(2,mq) = f1*Cxyz(2) + f2*Dxyz(2)
       Qxyz(3,mq) = f1*Cxyz(3) + f2*Dxyz(3)

  110 continue

c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.

      mr = 0
      do 125 mp = 1,NPP
        do 120 mq = 1,NPQ
          mr = mr + 1

          p = alpha_p(1,mp) + alpha_p(2,mp)
          q = alpha_q(1,mq) + alpha_q(2,mq)

          rho = p*q/(p+q)

          PQx = Pxyz(1,mp) - Qxyz(1,mq)
          PQy = Pxyz(2,mp) - Qxyz(2,mq)
          PQz = Pxyz(3,mp) - Qxyz(3,mq)

          PQ(mr,1) = PQx
          PQ(mr,2) = PQy
          PQ(mr,3) = PQz
c
          R(mr,0,1) = rho*(PQx**2 + PQy**2 + PQz**2)
          ff(1,mr) = 2.d0*dsqrt(rho/pi)
          ff(2,mr) = -2.D0*rho
c
          if (doscreen) then
           rhoscaled = rho
           call case_md(rhoscaled)
           Rscr(mr,0,1) = rhoscaled*(PQx**2 + PQy**2 + PQz**2)
           ffscr(1,mr) = 2.d0*dsqrt(rhoscaled/pi)
           ffscr(2,mr) = -2.D0*rhoscaled
          end if
c
  120   continue
  125 continue

c Evaluate the incomplete gamma function.

      call igamma(R,(NPP*NPQ),Lr)
      if (doscreen) call igamma(Rscr,(NPP*NPQ),Lr)

c Define the initial auxiliary functions (i.e., R000j, j=1,Lr).

      do 135 j = 0,Lr
        do 130 mr = 1,(NPP*NPQ)
             R(mr,j,1) = ff(1,mr)*R(mr,j,1)
             ff(1,mr) = ff(1,mr)*ff(2,mr)
  130   continue
             if (doscreen) then
                do  mr = 1,(NPP*NPQ)
                   Rscr(mr,j,1) = ffscr(1,mr)*Rscr(mr,j,1)
                   ffscr(1,mr) = ffscr(1,mr)*ffscr(2,mr)
                enddo
             endif
  135 continue

c Recursively build the remaining auxiliary functions (i.e., RIJKj, j=0).
        
      call hfmkr(R,IJK,PQ,(NPP*NPQ),Lr,Lr3)
      if (doscreen) call hfmkr(Rscr,IJK,PQ,(NPP*NPQ),Lr,Lr3)
       
c Transfer to R0 array.
        
c       write(6,*) "cam_alpha:",cam_alpha
c       write(6,*) "cam_beta:",cam_beta
c       write(6,*) "doscreen:",doscreen
c       write(6,*) "cam_srhf:",cam_srhf
c
      if (doscreen) then
         do  n = 1,Lr3
            if (cam_srhf) then  ! for pure short-range HF (1-erf(r)/r)
               do  mr = 1,(NPP*NPQ)
                  R0(mr,n) = R(mr,0,n)-Rscr(mr,0,n)
               enddo
            else
               do  mr = 1,(NPP*NPQ)
                  R0(mr,n) = cam_alpha*R(mr,0,n)+cam_beta*Rscr(mr,0,n)
               enddo
            endif 
         enddo
      else
         do  n = 1,Lr3
            do  mr = 1,(NPP*NPQ)
               R0(mr,n) = R(mr,0,n)
            enddo
         enddo
      endif
        
      end
