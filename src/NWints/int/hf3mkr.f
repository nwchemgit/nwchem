      Subroutine hf3mkr(Axyz,Bxyz,Cxyz,alpha,Gxyz,
     &    RS,PC,ff,R,R0,IJK,Nabc,Lg,Lg3)
c
c $Id: hf3mkr.f,v 1.1 1995-10-30 20:56:35 d3e129 Exp $
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
      
      Double Precision Gxyz(3,Nabc), PC(Nabc,3)
      Double Precision RS(Nabc), ff(2,Nabc), R(Nabc,0:Lg,Lg3)
c::local
      double precision PI, PI4
      Parameter (PI=3.1415926535898D0,PI4=4.D0/PI)
c
      double precision a, b, c, abci
      double precision Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz 
      double precision PCx, PCy, PCz 
      double precision alpha_t
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
        
        abci = 1/(a+b+c)
        
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
         
         PCx = Gxyz(1,m) - Cx
         PCy = Gxyz(2,m) - Cy
         PCz = Gxyz(3,m) - Cz
         
         R(m,0,1) = alpha_t*(PCx*PCx + PCy*PCy + PCz*PCz)
         
         PC(m,1) = PCx
         PC(m,2) = PCy
         PC(m,3) = PCz
         
00200  continue
       
c Evaluate the incomplete gamma function.
       
       call igamma(R,Nabc,Lg)
       
c Define the initial auxiliary functions (i.e., R000j, j=1,Lr).
       
       do 00300 j = 0,Lg
         do 00400 m = 1,Nabc
           R(m,j,1) = ff(1,m)*R(m,j,1)
           ff(1,m) = ff(1,m)*ff(2,m)
00400    continue
00300  continue
       
c Recursively build the remaining auxiliary functions (i.e., RIJKj, j=0).
       
       call hfmkr(R,IJK,PC,Nabc,Lg,Lg3)
       
c Transfer to R0 array.
       
       do 00500 n = 1,Lg3
         do 00600 m = 1,Nabc
           R0(m,n) = R(m,0,n)
00600    continue
00500  continue
       
       end
