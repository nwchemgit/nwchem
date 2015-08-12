       subroutine chicreat(rdffile,nd,nv,ngr,ims,nvv,kgrid,tv,
     * isv,den,mv,xv,yv,zv,chi,qvv,tau,okspace) 
       implicit none
       character*(*) rdffile
       real*8 wv(1:nvv,1:nvv,1:ngr), chi(1:nvv,1:nvv,1:ngr),hd(1:ngr)
       integer  nv, ngr, isv(1:nv), mv(1:nv),mvt(1:nvv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv),sigv(1:nv)    
       integer  i,j,j1,j2,icr,nvv,ims(1:nv),nd,ngt,i1
       logical okspace
       real*8 r1, r2, dr, dk, rgrid(1:ngr),kgrid(1:ngr),pi
       real*8 hsol(1:nd,1:ngr), rr(1:ngr),huk(1:nd,1:ngr)
       real*8 h1(1:nvv,1:nvv,1:ngr),h2(1:ngr),h3(1:nvv,1:nvv,1:ngr)
       character (4) tv(1:nv)
       real*8 qvv(1:nvv),qo(1:nvv),qm(1:nvv),tau
       real*8 dentv(1:nvv),del(1:nvv),sqm 
       integer iq
       integer fn
c
       external util_get_io_unit
       logical  util_get_io_unit
c
       pi=2*asin(1.0)  
       if(.not.util_get_io_unit(fn))
     >   call errquit("cannot get file number",0,0)
       open(unit=fn,status="old",form="formatted",file=rdffile)
c      reading of solvent rdfs 

        do i=1,ngr
         read(fn,*) rr(i), (hsol(j,i), j=1,nd)        
        enddo
c      reading and then symmetric indeces
       do j=1,nvv
        do j1=j,nvv
         do i=1,ngr 
          h3(j,j1,i)=hsol(nvv*(j-1)-j*(j-1)/2+j1,i)   
          h3(j1,j,i)=h3(j,j1,i)
         enddo           
        enddo
       enddo
c      wv(i,j,k) intramolecular solvent function
       call wcreat(nv,ngr,nvv,ims,kgrid,isv,xv,yv,zv,wv) 
c      chi(i,j,k) solvent susceptibility 
       do i=1,ngr
        do j1=1,nv
         do j2=1,nv
          chi(ims(j1),ims(j2),i)=wv(ims(j1),ims(j2),i)
     *    +mv(j1)*mv(j2)*h3(ims(j1),ims(j2),i)
         enddo
        enddo      
       enddo
       close(fn) 
       return
       end subroutine 
c
C>
C>      computes solvent susceptibility
C>      -------------------------------
C>    @param[in] rdffile name of the file containing rdf
C>    @param[in] nd total number of different rdfs (e.g. 2 for water),
C>     or another words number of unique solvent site-site pairs
C>    @param[in] nv total number of atom in solvent molecule (e.g. 3 for
C>     water)
C>    @param[in] ngr total number of grid points in rdf
C>    @param[in] ims total number of grid points in rdf
      
       subroutine chicreat_old(rdffile,nd,nv,ngr,ims,nvv,kgrid,tv,
     * isv,den,mv,xv,yv,zv,chi) 
       implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
       character*(*) rdffile
       real*8 wv(1:nvv,1:nvv,1:ngr), chi(1:nvv,1:nvv,1:ngr),hd(1:ngr)
       integer  nv, ngr, isv(1:nv), mv(1:nv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv),sigv(1:nv)    
       integer  i,j,j1,j2,icr,nvv,ims(1:nv),nd       
       real*8 r1, r2, dr, dk, rgrid(1:ngr),kgrid(1:ngr),pi
       real*8 hsol(1:nd,1:ngr), rr(1:ngr)
       real*8 h1(1:nvv,1:nvv,1:ngr),h2(1:ngr),h3(1:nvv,1:nvv,1:ngr)
       character (4) tv(1:nv)
c
       character*50 pname,message
       integer fn
       pname = "rism_solvent_chi"
       pi=2*asin(1.0)
c       
c     reading of solvent rdfs 
c     ------------------------
      message = "reading solvent rdf"
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=rdffile)
      do i=1,ngr
        read(fn,*,ERR=30,END=30) rr(i), (hsol(j,i), j=1,nd)
      enddo
      close(fn)
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
30    call errquit(pname//message,0,0)
       end subroutine 
c
C>     Caculation of reduced intra-molecular correlation matrix in k-space
C>     \f$\int \psi dx\f$       
C>     @param[in]  nv    number of solvent atoms      
C>     @param[in]  nvv   number of unique solvent atoms      
C>     @param[in]  ngr   number of grid points
C>     @param[in]  ims   atom type array
C>     @param[in]  kgrid kspace array
C>     @param[in]  isv residue/segment index of solvent atoms      
C>     @param[in]  xv  x-coord
C>     @param[in]  yv  y-coord
C>     @param[in]  zv  z-coord
C>     @param[out] wv  intra-molecular correlation matrix in k-space
C>     \psi
C>    
       subroutine wcreat(nv,ngr,nvv,ims,kgrid,isv,xv,yv,zv,wv) 
       implicit none
C>     @param[in] input panrameter       
       integer  nv,nvv, ngr
C>     input parameter       
       integer  isv(nv),ims(nv)
       double precision kgrid(1:ngr)  !< our fine level matrix
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
c      do i=1,ngr,40
c       write(43,*) ((wv(j1,j2,i), j1=1,nvv),j2=1,nvv)
c      enddo
      return
      end subroutine 

C>     Function rism_util_power_2
C>     returns true if number can be represented as power of 2
C>      
       function rism_util_power_2(n)
       implicit none
       logical rism_util_power_2
       integer n
       integer i,nsb,imsb

       nsb = 0
       do i=0,bit_size(n)-1
        if(btest(n,i)) then
          nsb = nsb+1
          imsb = i
        end if
       end do

       rism_util_power_2 = nsb.eq.1
       return
       end

      subroutine rism_message(buffer)
      implicit none
#include "stdio.fh"
#include "global.fh"
#include "inp.fh"
      character*(*) buffer
c
      character*12 cnode
      write(cnode,'("node ",I3,": ")') ga_nodeid()
      write(luout,*) cnode,buffer(1:inp_strlen(buffer))
      call util_flush(luout) 
      return
      end
c $Id$
