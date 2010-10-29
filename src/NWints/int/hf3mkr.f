      Subroutine hf3mkr(Axyz,Bxyz,Cxyz,alpha,Gxyz,
     &    RS,GC,ff,R,R0,IJK,Nabc,Lg,Lg3)
c
c $Id$
c
      Implicit none 
c::passed
      integer Nabc, Lg, Lg3
c--> Cartesian Coordinates for centers a, b, c
      
      Double Precision Axyz(3),Bxyz(3),Cxyz(3) ! [input]
      
c--> Exponents (1:3,*) and ES prefactors (4:4,*)
      
      Double Precision alpha(4,Nabc) ! [input]
      
c--> Auxiliary Function Integrals & Index
      
      Double Precision R0(Nabc,Lg3) ! [output]
      Integer IJK(0:Lg,0:Lg,0:Lg) ! [output]
      
c--> Scratch Space
      
      Double Precision Gxyz(3,Nabc), GC(Nabc,3)
      Double Precision RS(Nabc), ff(2,Nabc), R(Nabc,0:Lg,Lg3)
c::local
      double precision PI, PI4
      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)
c
      double precision a, b, c, abci
*rak: double precision ab, abi
      double precision Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz 
*rak:      Double Precision Px, Py, Pz
      double precision GCx, GCy, GCz 
      double precision alpha_t
*acc_debug:      double precision accy
*acc_debug:      integer accy_cnt
*acc_debug:      logical reached
      integer mp, j, m, n

c
c Define the auxiliary function integrals necessary to compute the three 
c center nuclear attraction integrals (NAIs). These integrals are scaled 
c by an appropriate factor, RS, defined as
c
c         / a + b + c \ 1/2
c   RS = | ----------- |
c         \    PI/4   /
c
c
c******************************************************************************
*      call dfill((Nabc*(Lg+1)*Lg3), 0.0d00 , R, 1)
*      call dfill((Nabc*Lg3),0.0d00, r0, 1)
      
c Define the center "P" plus C to get "G" center.
      
      Ax = Axyz(1)
      Ay = Axyz(2)
      Az = Axyz(3)
      
      Bx = Bxyz(1)
      By = Bxyz(2)
      Bz = Bxyz(3)
      
      Cx = Cxyz(1)
      Cy = Cxyz(2)
      Cz = Cxyz(3)
      
      do 00100 mp = 1,Nabc
        
        a = alpha(1,mp)
        b = alpha(2,mp)
        c = alpha(3,mp)

*rak:        ab = a + b
*rak:        abi = 1/ab
*rak:
*rak:        px = abi*(a*Ax + b*Bx)
*rak:        py = abi*(a*Ay + b*By)
*rak:        pz = abi*(a*Az + b*Bz)
*rak:   abci = 1/(ab+c)        

        abci = 1/(a+b+c)
        
*rak:        Gxyz(1,mp) = abci*(ab*px + c*Cx)
*rak:        Gxyz(2,mp) = abci*(ab*py + c*Cy)
*rak:        Gxyz(3,mp) = abci*(ab*pz + c*Cz)
        
        Gxyz(1,mp) = abci*(a*Ax + b*Bx + c*Cx)
        Gxyz(2,mp) = abci*(a*Ay + b*By + c*Cy)
        Gxyz(3,mp) = abci*(a*Az + b*Bz + c*Cz)

c Define the scaling factor.
        
        RS(mp) = sqrt((a+b+c)*PI4)
        
00100 continue
      
c Define factors necessary to compute incomplete gamma function and the
c auxiliary functions.
      
       do 00200 m = 1,Nabc

         alpha_t = alpha(1,m) + alpha(2,m) + alpha(3,m)
         
         ff(1,m) = RS(m)
         ff(2,m) = -2.D0*alpha_t
         
         GCx = Gxyz(1,m) - Cx
         GCy = Gxyz(2,m) - Cy
         GCz = Gxyz(3,m) - Cz
         
         R(m,0,1) = alpha_t*(GCx*GCx + GCy*GCy + GCz*GCz)
         
         GC(m,1) = GCx
         GC(m,2) = GCy
         GC(m,3) = GCz
         
00200  continue
       
c Evaluate the incomplete gamma function.

       call igamma(R,Nabc,Lg)

*acc_debug:       accy = 1.0d-30
*acc_debug:       accy_cnt = 0
*acc_debug:       reached = .false.
*acc_debug:00001  continue
*acc_debug:       call igamma_acc(R,Nabc,Lg,accy,reached)
*acc_debug:       if (.not.reached) then
*acc_debug:         accy_cnt = accy_cnt + 1
*acc_debug:         accy = accy/5.0d00
*acc_debug:         write(6,*)' accy reset to ',accy,' count is ',accy_cnt
*acc_debug:         goto 00001
*acc_debug:       endif

c Define the initial auxiliary functions (i.e., R000j, j=1,Lg).
       
       do 00300 j = 0,Lg
         do 00400 m = 1,Nabc
           R(m,j,1) = ff(1,m)*R(m,j,1)
           ff(1,m) = ff(1,m)*ff(2,m)
00400    continue
00300  continue
       
c Recursively build the remaining auxiliary functions (i.e., RIJKj, j=0).
       
       call hfmkr(R,IJK,GC,Nabc,Lg,Lg3)
       
c Transfer to R0 array.
       
       do 00500 n = 1,Lg3
         do 00600 m = 1,Nabc
           R0(m,n) = R(m,0,n)
00600    continue
00500  continue
       
       end
