c
       subroutine chicreat(rdffile,nd,nv,ngr,ims,nvv,kgrid,tv,
     * isv,den,mv,xv,yv,zv,chi) 
       implicit none
       character*(*) rdffile
       real*8 wv(1:nvv,1:nvv,1:ngr), chi(1:nvv,1:nvv,1:ngr),hd(1:ngr)
       integer  nv, ngr, isv(1:nv), mv(1:nv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv),sigv(1:nv)    
       integer  i,j,j1,j2,icr,nvv,ims(1:nv),nd       
       real*8 r1, r2, dr, dk, rgrid(1:ngr),kgrid(1:ngr),pi
       real*8 hsol(1:nd,1:ngr), rr(1:ngr)
       real*8 h1(1:nvv,1:nvv,1:ngr),h2(1:ngr),h3(1:nvv,1:nvv,1:ngr)
       character (4) tv(1:nv)
       pi=2*asin(1.0)
       
c 
c reading of solvent rdfs 
        open(3, file=rdffile, status='old')      
        do i=1,ngr
         read(3,104) rr(i), (hsol(j,i), j=1,nd)
        enddo
       close(3)
104    format(11(1x,f12.7)) 
       dr=rr(2)-rr(1) 
!      do i=1,ngr,40
!       print*, (hsol(j,i), j=1,9)
!      enddo
c 
c      sin fft transform of rdfs with proper arrangement in array
       do j=1,nvv
        do j1=j,nvv
         do i=2,ngr 
          h1(j,j1,i)=(hsol(nvv*(j-1)-j*(j-1)/2+j1,i)-1)*rr(i) 
          h2(i)=h1(j,j1,i)
         enddo
         call sinft(h2,ngr)
c        normalization of sin-fft with excluding the zeropoint (x=0)
         do i=1,ngr-1
          h3(j,j1,i)=h2(i+1)/kgrid(i)
         enddo
         h3(j,j1,ngr)=h3(j,j1,ngr-1)
        enddo
       enddo
c      symmetric indeces
       do j=1,nvv
        do j1=j,nvv
         do i=1,ngr 
          h3(j1,j,i)=h3(j,j1,i)
         enddo           
        enddo
       enddo
!      do i=1,ngr,40       
!       print*, (h3(j,3,i), j=1,4) 
!      enddo 
c
c wv(i,j,k) intramolecular solvent function
c
       call wcreat(nv,ngr,nvv,ims,kgrid,isv,xv,yv,zv,wv) 
!      do i=1,ngr,40
!       print*, (wv(1,j1,i), j1=1,nvv)
!      enddo
c
c chi(i,j,k) solvent susceptibility
c      
       dr=rr(2)-rr(1)   
       do i=1,ngr
        do j1=1,nv
         do j2=1,nv
          chi(ims(j1),ims(j2),i)=wv(ims(j1),ims(j2),i)
     *    +4*pi*dr*mv(j1)*mv(j2)*den(j1)*h3(ims(j1),ims(j2),i)
         enddo
        enddo      
       enddo
!      do i=1,ngr,40
!       print*, (chi(2,j,i), j=1,4) 
!      enddo
       return
       end subroutine 
c
c
c caculation of reduced solvent matrix
c
       subroutine wcreat(nv,ngr,nvv,ims,kgrid,isv,xv,yv,zv,wv) 
       implicit none
       integer  nv,nvv, ngr
       integer  isv(nv),ims(nv)
       double precision kgrid(1:ngr)
       double precision xv(1:nv),yv(1:nv),zv(1:nv)
       double precision wv(nvv,nvv,ngr)
c 
       integer  i, j1,j2       
       double precision pi
       double precision dist,co
       double precision sinr
       external sinr
c
       pi=2*asin(1.0)  
c
c      initialization to zero
c      ----------------------
       call dfill(0.0d0,wv,1,nvv*nvv*ngr)
c 
c      dist(i,j) distance between i and j sites
c      co(i,j,k) sinr between i and j sites
       do i=1,ngr
        do j1=1,nv       
         do j2=1,nv
          if(isv(j1).eq.isv(j2)) then        
           dist=sqrt(((xv(j1)-xv(j2))**2+
     *     (yv(j1)-yv(j2))**2+(zv(j1)-zv(j2))**2))
           co=sinr(kgrid(i)*dist)
           wv(ims(j1),ims(j2),i)= wv(ims(j1),ims(j2),i)+co
          endif
         enddo  
        enddo 
       enddo
c
c     testing
c     -------
      do i=1,ngr,40
       write(43,*) ((wv(j1,j2,i), j1=1,nvv),j2=1,nvv)
      enddo
      return
      end subroutine 

