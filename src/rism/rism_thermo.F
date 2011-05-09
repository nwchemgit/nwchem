       subroutine thermo(nu,nv,nvv,ngr,icl,ims,den,t,rgrid,
     * hnew,cnew,ul) 
       implicit none
       real*8 muh, mugf, muuh(1:nu,1:nvv),muugf(1:nu,1:nvv)        
       integer  nu, nv,nvv, ngr, i,j1,j,icl        
       real*8 rgrid(1:ngr)  
       real*8 pi,dr,kb,t,kcal,theta
       real*8 den(1:nv)
       integer ims(1:nv)
       real*8 cnew(1:nu,1:nvv,1:ngr), hnew(1:nu,1:nvv,1:ngr)
       real*8 ul(1:nu,1:nvv,1:ngr),hq(1:nu,1:nvv,1:ngr)
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
     *    (cnew(j,j1,i)-1.0/kb/t*ul(j,j1,i)))*rgrid(i)**2*dr*2*pi
          muuh(j,j1)=muuh(j,j1)-(2*cnew(j,j1,i)+hnew(j,j1,i)*
     *    (cnew(j,j1,i)-1.0/kb/t*ul(j,j1,i)))*rgrid(i)**2*dr*2*pi
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
