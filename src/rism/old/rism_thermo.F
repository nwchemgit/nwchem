
       subroutine thermo(nu,nv,nvv,ngr,icl,ims,den,t,rgrid,qqu,qvv,
     * hnew,cnew,solperm,tau) 
       implicit none
       real*8 muh, mugf, muuh(1:nu,1:nvv),muugf(1:nu,1:nvv)        
       integer  nu, nv,nvv, ngr, i,j1,j,icl        
       real*8 rgrid(1:ngr)  
       real*8 pi,dr,kb,t,kcal,theta,ck,echar,solperm,nav,tau
       real*8 den(1:nv)
       integer ims(1:nv)
       real*8 cnew(1:nu,1:nvv,1:ngr), hnew(1:nu,1:nvv,1:ngr)
       real*8 ul(1:nu,1:nvv,1:ngr),hq(1:nu,1:nvv,1:ngr)
       real*8 ulr(1:nu,1:nvv,1:ngr),qvv(1:nvv),qqu(1:nu)
       real*8 erf
#ifndef GFORTRAN
       external erf
#endif
c
       pi=2*asin(1.0)
       kb=8.13441e-3
       kcal=0.00198717 
       mugf=0 
       muh=0
       do j=1,nu
        do j1=1,nvv
         muugf(j,j1)=0
         muuh(j,j1)=0
        enddo
       enddo
       dr=rgrid(2)-rgrid(1)
       nav=6.02214179e+23
       echar=4.803e-10
       ck=nav*echar**2/solperm*0.01
c  ulr electrostatic part of the interaction potential in real space
       do j=1,nu
        do j1=1,nvv     
         do i=1,ngr 
          ulr(j,j1,i)=ck*qqu(j)*qvv(j1)/rgrid(i)*erf(tau*rgrid(i))/4.1868
         enddo    
        enddo
       enddo 
c calculations of partial contibutions 
       do j=1,nu
        do j1=1,nvv  
         do i=1,ngr
          if(icl.eq.1) then 
            hq(j,j1,i)=hnew(j,j1,i)
          endif
          if(icl.eq.2) then 
           hq(j,j1,i)= -theta(-hnew(j,j1,i))
          endif           
          muugf(j,j1)=muugf(j,j1)-(2*cnew(j,j1,i)+hnew(j,j1,i)*
     *    (cnew(j,j1,i)-ulr(j,j1,i)/kcal/t))*rgrid(i)**2*dr*2*pi
          muuh(j,j1)=muuh(j,j1)-(2*cnew(j,j1,i)+hnew(j,j1,i)*
     *    (cnew(j,j1,i)-ulr(j,j1,i)/kcal/t))*rgrid(i)**2*dr*2*pi
     *    +hnew(j,j1,i)*hq(j,j1,i)*rgrid(i)**2*dr*2*pi   
         enddo    
        enddo
       enddo 
       do j=1,nu
        do j1=1,nv
         mugf=mugf+den(j1)*muugf(j,ims(j1))*kcal*t 
         muh=muh+den(j1)*muuh(j,ims(j1))*kcal*t       
        enddo
       enddo             
c           print*, mugf,muh
c  mugf chem.potential in gaussian approximation
c muh   chem.potential in hnc approximation         
       open(3, file='tout.data', status='unknown')         
        write(3,108) muh,mugf
       close(3)
  108  format(2(2x,e10.4))
       return 
       end subroutine
c
      real*8 function theta(x)
      real*8 x
       if(x.le.0.0)then
        theta=x
        else
        theta=0
       endif
       return
       end      
c
c $Id: rism_thermo.F 21176 2011-10-10 06:35:49Z d3y133 $
