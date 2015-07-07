c program qqq 
c with diis implementation
       implicit none
       integer  nu, nv, ngr ,icr,icl,i,dd 
       real*8 r1, r2, tau, solperm,t, tol,lambd,dr  
c parameters(nv number of solute sites  nu number of solute sites ngr number of grid points)
c solperm solvent relative permittivity
c tau(aa-1) separation between short and long range functions in 1/angstroms
c icr type of combination rule: 1-arithmetic; 2-geometric 
c t(k)  temperature in kelvins
c tol residue tolerance
c lamd mixing parameter
c icl closure type: 1-hnc; 2-kh
c tau debye parameter for spliting potnetial into short and long part
c solperm dielectric permitivity of the solvent
c dd number of diis basis vectors
c      ngr=4096
c      nu=12 
c      nv=4 
c      tau=1
c      solperm=1
c      icr=1
c      t=298
c      tol=1e-5
c      lambd=0.5 
c      icl=2 
       open(3, file='uvpar.data', status='old')
       read(3,102) nu,nv,icr,icl,ngr,tau,solperm,t,tol,lambd,dd 
       print*, nu,nv,icr,icl,ngr,tau,solperm,t,tol,lambd,dd 
       close(3)
c parameters(r1 first grid point r2 second grid point)
c      reading of r1 and r2  
       open(3, file='full.data', status='old')      
        read(3,104) r1,r2          
       close(3)
104    format(1x,f12.7) 
       dr=r2-r1
!      print*, r1,r2,dr
       call risminput(nv,nu,ngr,dr,r1,tau,
     * solperm,icr,t,tol,lambd,icl,dd) 
 102   format(4(1x,i4),1x,i8,5(2x,e9.4),1x,i4) 
       stop
       end 
cccc
       subroutine risminput(nv,nu,ngr,dr,r1,tau,
     * solperm,icr,t,tol,lambd,icl,dd) 
       implicit none
       integer  nu, isu(1:nu), mu(1:nu)
       real*8 wu(1:nu,1:nu,1:ngr),sgvv(1:nv)
       real*8 xu(1:nu),yu(1:nu),zu(1:nu),sigu(1:nu)
       real*8 epsiu(1:nu),qqu(1:nu) 
       integer  nv, ngr, isv(1:nv), mv(1:nv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv),sigv(1:nv)
       real*8 epsiv(1:nv),qqv(1:nv)    
       real*8 sigvv(1:nv), epsvv(1:nv),qvv(1:nv)
       integer  i, icr, icl,nvv,ims(1:nv),nd,k,dd       
       real*8 r1,dr, dk, rgrid(1:ngr),kgrid(1:ngr),pi 
       character (4) tv(1:nv), tu(1:nu)
       real*8 tau,solperm,t, tol,lambd,kb
       integer j1,j
       real*8 c3(1:nu,1:nv,1:ngr),hu(1:nu,1:nv,1:ngr) 
       real*8 gam(1:nu,1:nv,1:ngr)
c
c      solvent input
c tv type of cite isv type of species mv multiplicity  den density (1/a3) xvyvzv[xyz](a) 
c sigv  sigma(a) epsiv epsilon(kj/mol) qqv charge(e)
c      number of lines in solvt..data should be equal nv
       open(3, file='solvent3.data', status='old')      
        do i=1,nv
         read(3,103)tv(i),isv(i),mv(i),den(i),  
     *   xv(i),yv(i),zv(i),sigv(i), epsiv(i),qqv(i) 
!        print *, tv(i),isv(i),mv(i),den(i),  
!    *   xv(i),yv(i),zv(i),sigv(i), epsiv(i),qqv(i)  
        enddo
       close(3)
 103   format(a4,2x,i4,2x,i4,7(2x,e10.4))
c
c      sorting
       call sort(nv,tv,isv,mv,ims,nvv) 
!      print*, (ims(i), i=1,nv)
c      calculations interaction parameters for reduced matrix
       call vpot(nv,ims,nvv,sigv,epsiv,qqv,sgvv,epsvv,qvv)
!      do i=1,nvv 
!       print*, sgvv(i),epsvv(i),qvv(i),i,ims(i)
!      enddo
c
c      solute input
c tu type of cite isu type of species mu multiplicity  xuyuzu[xyz](a) 
c      sigu sigma(a) epsiu epsilon(kj/mol) qqu charge(e)
       open(3, file='solute2.data', status='old')      
        do i=1,nu
         read(3,102)tu(i),isu(i),mu(i),  
     *   xu(i),yu(i),zu(i),sigu(i),epsiu(i),qqu(i) 
         print *,tu(i),isu(i),mu(i),  
     *   xu(i),yu(i),zu(i),sigu(i),epsiu(i),qqu(i)
        enddo
       close(3)
 102   format(a4,2x,i4,2x,i4,6(2x,e10.4))                
c
c      creat rgrid rgrd(i) creat kgrid kgrid(i)
       pi=2*asin(1.0)
c      the zero point is excluded
       r1=dr   
       dk=pi/(ngr+1)/dr
       do i=1,ngr
        rgrid(i)=r1+dr*(i-1) 
        kgrid(i)=i*dk 
       enddo 
c
c      wu(i,j,k) intramolecular solute function
       call wucreat(nu,ngr,kgrid,tu,isu,xu,yu,zu,wu) 
!      do k=1,ngr,40
!       print*, (wu(1,i,k), i=1,9)
!      enddo
c       
       call rism(nu,nv,nvv,ngr,icl,kgrid,rgrid,tv,isv,den,mv,xv,yv,zv,
     * icr,ims,tau,solperm,sgvv,epsvv,qvv,sigu,epsiu,qqu,t,lambd,
     * tol,wu,dd)          
       return
       end subroutine 
cccc
c
c subroutine for rearrangement of the solvent data which excludes the same sites
c defined by paramters tv and mv
c
       subroutine sort(nv,tv,isv,mv,ims,nvv)
       implicit none
       integer  nv, isv(1:nv), mv(1:nv),ims(1:nv),mvt(1:nv)
       integer  i,i1,j1,j2,icc,imm(1:nv),max,nvv      
       character (4) tv(1:nv)
c
c      initialization of connectivity vector      
       do i=1,nv
        ims(i)=i 
        imm(i)=1
        mvt(i)=mv(i)  
       enddo
c      sorting 
       do j1=1,nv-1
        do j2=j1+1,nv
         if((tv(j1).eq.tv(j2)).and.(isv(j1).eq.isv(j2))
     *    .and.(mv(j2).ne.1)) then 
c         replace current value of ims by the first found
          ims(j2)=ims(j1)
c         change indicator of replacement
          imm(j2)=2
c         change multiplicity
          mv(j2)=1 
c         shift all others except replaced by 1
          if(j2.ne.nv) then
           do i1=j2+1,nv
            if(imm(i1).eq.1) then  
             ims(i1)=ims(i1)-1
             else 
             ims(i1)=ims(i1)
            endif 
           enddo
          endif  
          else 
          ims(j2)=ims(j2)
         endif
!       print*, j2,ims(j2),mv(j2),imm(j2)
        enddo
       enddo
c      evaluation of maximal value of different sites
       max=ims(1) 
       do i=2,nv
        if(ims(i).ge.max) then
         max=ims(i)
         else
         max=max 
        endif
       enddo
       nvv=max
c      replace changed mv by the initial one 
       do i=1,nv
        mv(i)=mvt(i)
       enddo
       return
       end subroutine 
c
c      subroutine for calculation potential parameters for reduced matrix
c
       subroutine vpot(nv,ims,nvv,sigv,epsiv,qqv,sgvv,epsvv,qvv) 
       implicit none
       integer  nv,nvv 
       real*8 sgvv(1:nvv),epsvv(1:nvv),qvv(1:nvv)        
       real*8 sigv(1:nv),epsiv(1:nv),qqv(1:nv)    
       integer  i,ims(1:nv),j1,j2,icr       
       do i=1,nv  
        sgvv(ims(i))=sigv(i)
        epsvv(ims(i))=epsiv(i)
        qvv(ims(i))=qqv(i)
!       print*, sgvv(ims(i)),sigv(i),i,ims(i)
       enddo 
       return
       end subroutine
c
c      prepration of wu(i,j,k) intramoplecular solute matrix
c
       subroutine wucreat(nu,ngr,kgrid,tu,isu,xu,yu,zu,wu)
       implicit none
       integer  nu, ngr, isu(1:nu), mu(1:nu)
       real*8 wu(1:nu,1:nu,1:ngr), dist(1:nu,1:nu)
       real*8 xu(1:nu),yu(1:nu),zu(1:nu)   
       integer  i, j1,j2,icr       
       real*8 kgrid(1:ngr),sik  
       character (4) tu(1:nu) 
c      dis(i,j) intramolecular distances      
       do i=1,ngr
        do j1=1,nu
         do j2=1,nu
          dist(j1,j2)=(((xu(j1)-xu(j2))**2+(yu(j1)-yu(j2))**2+
     *    (zu(j1)-zu(j2))**2))**(1.0/2)
          wu(j1,j2,i)=sik(kgrid(i)*dist(j1,j2))
         enddo
        enddo
       enddo 
       return
       end subroutine 
c
       real*8 function sik(x)
       real*8 x
       if(x.lt.1d-22)then
        sik=1
       else
       sik=sin(x)/x
       endif
       return
       end
c
c creat susceptibility of solvent
c
       subroutine chicreat(nd,nv,ngr,ims,nvv,kgrid,tv,
     * isv,den,mv,xv,yv,zv,chi) 
       implicit none
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
        open(3, file='full.data', status='old')      
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
       call wcreat(nv,ngr,nvv,ims,kgrid,tv,isv,xv,yv,zv,wv) 
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
c caculation of reduced solvent matrix
c
       subroutine wcreat(nv,ngr,nvv,ims,kgrid,tv,isv,xv,yv,zv,wv) 
       implicit none
       real*8 wv(1:nvv,1:nvv,1:ngr),dist(1:nv,1:nv)
       real*8 ws(1:nv,1:nvv,1:ngr)
       integer  nv, ngr, isv(1:nv),mv(1:nv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv),sigv(1:nv)    
       integer  i, j1,j2,nvv,ims(1:nv),t       
       real*8 kgrid(1:ngr),pi,sik,co(1:nv,1:nv,1:ngr) 
       character (4) tv(1:nv)
       pi=2*asin(1.0)  
c
c      initialization of wv(i,j,k) & calculating of the 11 element
c
       do i=1, ngr         
        do j1=1,nv
         do j2=1,nv          
           wv(ims(j1),ims(j2),i)=0
           ws(j1,ims(j2),i)=0  
         enddo
        enddo
       enddo            
c      dist(i,j) distance between i and j sites
c      co(i,j,k) sik between i and j sites
c      reducing number of coloumns
       do i=1,ngr
        do j1=1,nv       
         do j2=1,nv
c         calculations all other columns
c         if they belong to the same molecule
          if(isv(j1).eq.isv(j2)) then        
           dist(j1,j2)=(((xv(j1)-xv(j2))**2+
     *     (yv(j1)-yv(j2))**2+(zv(j1)-zv(j2))**2))**(1.0/2)           
           co(j1,j2,i)=sik(kgrid(i)*dist(j1,j2))
           ws(j1,ims(j2),i)= ws(j1,ims(j2),i)+co(j1,j2,i)
          endif
         enddo  
        enddo 
       enddo
c      reducing the number of rows 
       do i=1,ngr        
        do j2=1,nvv       
         do j1=1,nv   
          wv(ims(j1),j2,i)=wv(ims(j1),j2,i)+ws(j1,j2,i)        
         enddo  
        enddo
       enddo 
!      do i=1,ngr,40
!       print*, (wv(1,j1,i), j1=1,nvv)
!      enddo
       return
       end subroutine 
c
c  creat solute-solvent potentials
c
       subroutine potcreat(nvv,nu,ngr,icr,kgrid,rgrid,tau,solperm,
     * sgvv,epsvv,qvv,sigu,epsiu,qqu,plj,ul) 
       implicit none
       integer  nu,nv, ngr,nvv 
       real*8 sigu(1:nu),epsiu(1:nu),qqu(1:nu)        
       real*8 sgvv(1:nvv),epsvv(1:nvv),qvv(1:nvv)    
       integer  i, j1,j2,icr       
       real*8 r1, r2, dr, dk, rgrid(1:ngr),kgrid(1:ngr),pi 
       real*8 sigvv(1:nu,1:nvv), epsivv(1:nu,1:nvv)     
       real*8 ul(1:nu,1:nvv,1:ngr),plj(1:nu,1:nvv,1:ngr)
       real*8 nav,tau, echar,solperm,uthre,ck 
c      constants for potential
       pi=2*asin(1.0)
       nav=6.02214179e+23
       echar=4.803e-10
       ck=nav*echar**2/solperm*0.01
       uthre=1000   
c
c calculation lj parameters
       call combrule(nu,nvv,icr,sgvv,sigu,epsiu,epsvv,sigvv,epsivv)
c
c plj(j1,j2,r) lj potential+short part from coulomb
c ul(j1,j2,k) is the fourier trasnform of long range interactions
       do i=1,ngr 
        do j1=1,nu
         do j2=1,nvv
          ul(j1,j2,i)=4*pi*ck*qqu(j1)*qvv(j2)*exp(-kgrid(i)**2/4/tau**2)
          ul(j1,j2,i)=ul(j1,j2,i)/kgrid(i)**2
          plj(j1,j2,i)=4*epsivv(j1,j2)*
     *    ((sigvv(j1,j2)/rgrid(i))**(12)-(sigvv(j1,j2)/rgrid(i))**(6))
          plj(j1,j2,i)=plj(j1,j2,i)+ck
     *    *qqu(j1)*qvv(j2)/rgrid(i)*(1-erf(tau*rgrid(i)))
          if (plj(j1,j2,i).ge. uthre) then 
           plj(j1,j2,i)=uthre
          endif 
         enddo
        enddo
       enddo
!      do i=1,40
!       print*, (plj(1,j1,i), j1=1,4)
!      enddo 
       return
       end subroutine 
c
c   calculations of parameters for solute-solvent interactions
c
       subroutine combrule(nu,nvv,icr,sgvv,sigu,
     * epsiu,epsvv,sigvv,epsivv) 
       implicit none
       integer  nu, nv, icr,nvv
       real*8 sigu(1:nu), epsiu(1:nu), sgvv(1:nvv), epsvv(1:nvv)          
       real*8 sigvv(1:nu,1:nvv), epsivv(1:nu,1:nvv)
       integer  i, j1,j2   
       do j1=1,nu
        do j2=1,nvv
         epsivv(j1,j2)=(epsiu(j1)*epsvv(j2))**(1.0/2)
         if (icr.eq.1) then
          sigvv(j1,j2)=(sigu(j1)+sgvv(j2))/2
          else
          sigvv(j1,j2)=(sigu(j1)*sgvv(j2))**(1.0/2)
         endif
!        print *, sigvv(j1,j2), epsivv(j1,j2), j1,j2
        enddo
       enddo
       return
       end subroutine  
c 
c site-site ornstein-zernike in k-space
c
       subroutine ssoz(nvv,nu,ngr,nv,ims,mv,wu,c3,chi,hu)
       implicit none
       integer  nu, nvv, nv, ngr, i,j1,j,k1,ims(1:nv),mv(1:nv)        
       real*8 rgrid(1:ngr),kgrid(1:ngr)  
       real*8 pi,ct
       real*8 c3(1:nu,1:nvv,1:ngr),hu(1:nu,1:nvv,1:ngr)
       real*8 ctt(1:nu,1:nvv,1:ngr),h1(1:nu,1:nvv,1:ngr)
       real*8 wu(1:nu,1:nu,1:ngr), chi(1:nvv,1:nvv,1:ngr)  
c      calculations of right product of matrix via summation over solvent sites         
       do i=1,ngr
        do j1=1,nu 
         do j=1,nvv
          ct=0
          do k1=1,nvv
           ct=c3(j1,k1,i)*chi(k1,j,i)+ct
          enddo 
          ctt(j1,j,i)=ct
         enddo
        enddo
       enddo
!      do i=1,ngr,40       
!       print*, (ctt(1,j,i), j=1,4) 
!      enddo 
c      calculations of left product of matrix via summation over solute sites        
       do i=1,ngr
        do j1=1,nu 
         do j=1,nvv
          ct=0
          do k1=1,nu
           ct=wu(j1,k1,i)*ctt(k1,j,i)+ct
          enddo 
          hu(j1,j,i)=ct
          h1(j1,j,i)=hu(j1,j,i)
         enddo
        enddo
       enddo
       do i=1,ngr
        do j1=1,nu 
         do j=1,nv
          hu(j1,ims(j),i)=h1(j1,ims(j),i)/mv(j)
         enddo
        enddo
       enddo
!      do i=1,ngr,40       
!       print*, (hu(1,j,i),h1(1,j,i), j=1,4)  
!      enddo 
       return
       end subroutine
c
c  subroutine for closure
c                       
       subroutine clos(nvv,nu,ngr,t,plj,cr,gam,icl)
       implicit none
       integer  nu, nvv, ngr, i,j1,j,icl
       real*8 plj(1:nu,1:nvv,1:ngr)          
       real*8 t, kb,pi,lexp
       real*8 cr(1:nu,1:nvv,1:ngr),gam(1:nu,1:nvv,1:ngr)           
       kb=8.13441e-3 
       do i=1,ngr
        do j=1,nu
         do j1=1,nvv
          if(icl.eq.1) then
           cr(j,j1,i)=exp(-1.0/kb/t*plj(j,j1,i)+gam(j,j1,i))
     *     -gam(j,j1,i)-1
          endif  
          if(icl.eq.2) then
           cr(j,j1,i)=lexp(-1.0/kb/t*plj(j,j1,i)+gam(j,j1,i))
     *     -gam(j,j1,i)-1
          endif
         enddo
        enddo
       enddo
       return
       end subroutine 
c
       real*8 function lexp(x)
       real*8 x
       if(x.le.0.0)then
        lexp=exp(x)
        else
        lexp=1+x
       endif
       return
       end      
c
ccc
c    
       subroutine rism(nu,nv,nvv,ngr,icl,kgrid,rgrid,tv,isv,
     * den,mv,xv,yv,zv,
     * icr,ims,tau,solperm,sgvv,epsvv,qvv,sigu,epsiu,qqu,t,
     * lambd,tol,wu,dd)       
       implicit none
c      parameters
       integer  nu, nv,nd,nvv, ngr, i,j1,j,k1,icl,k        
       real*8 rgrid(1:ngr),kgrid(1:ngr)  
       real*8 pi,dk,t,tau,solperm,tol,lambd,del0,kb,dr
c      solute
       real*8 wu(1:nu,1:nu,1:ngr), chi(1:nvv,1:nvv,1:ngr)   
       real*8 sigu(1:nu),epsiu(1:nu),qqu(1:nu)
c      solvent 
       integer  icr,isv(1:nv),mv(1:nv),ims(1:nv)
       real*8 den(1:nv),xv(1:nv),yv(1:nv),zv(1:nv)  
       real*8 wv(1:nvv,1:nvv,1:ngr)
       real*8 sgvv(1:nvv),epsvv(1:nvv),qvv(1:nvv)
       character (4) tv(1:nv)
c      potential
       real*8 plj(1:nu,1:nvv,1:ngr), ul(1:nu,1:nvv,1:ngr)
c      functions
       real*8 cr(1:nu,1:nvv,1:ngr),c2(1:ngr),tt(1:nu,1:nvv)  
       real*8 ht(1:ngr),c3(1:nu,1:nvv,1:ngr),cf3(1:nu,1:nvv,1:ngr)
       real*8 cold(1:nu,1:nvv,1:ngr),cnew(1:nu,1:nvv,1:ngr)
       real*8 gfold(1:nu,1:nvv,1:ngr),hu(1:nu,1:nvv,1:ngr) 
       real*8 gold(1:nu,1:nvv,1:ngr), gnew(1:nu,1:nvv,1:ngr)   
       real*8 hold(1:nu,1:nvv,1:ngr), hnew(1:nu,1:nvv,1:ngr)             
       real*8 del(1:nu,1:nvv),mmaxi,mmax(1:nu),mmmax 
c      diis functions and counts
       real*8 ggo(1:dd,1:nu,1:nvv,1:ngr),dgg(1:dd,1:nu,1:nvv,1:ngr)
       integer kd,dd 
c
       pi=2*asin(1.0)       
c      nd total number of different rdfs
       nd=(nvv+1)*nvv/2     
c      chi(i,j,k) solvent susceptibility
c
       call chicreat(nd,nv,ngr,ims,nvv,kgrid,tv,isv,den,mv,xv,yv,zv,chi) 
!      do i=1,ngr,40 
!       print*, (chi(2,j,i), j=1,4) 
!      enddo
c
c      plj(j1,j2,r) lj potential+short part from coulomb
c      ul(j1,j2,k) is the fourier trasnform of long range interactions
c
       call potcreat(nvv,nu,ngr,icr,kgrid,rgrid,tau,solperm,
     * sgvv,epsvv,qvv,sigu,epsiu,qqu,plj,ul) 
!      do k=1,40
!       print*, (plj(i,1,k), i=1,9)
!      enddo            
cc     initial value of c_short(r)       
       kb=8.13441e-3
       do i=1,ngr
        do j=1,nu
         do j1=1,nvv
          gold(j,j1,i)=0
          cold(j,j1,i)=0
         enddo
        enddo
       enddo  
c      starting counts for diis (kd) and for iterations k1
       k1=1
       kd=1
  1    continue  
c      calculating c_new(r) and h_new(r) from closure
       call clos(nvv,nu,ngr,t,plj,cnew,gold,icl)
       do j=1,nu
        do j1=1,nvv
         do i=1,ngr
          hold(j,j1,i)=gold(j,j1,i)+cnew(j,j1,i) 
         enddo
        enddo
       enddo 
!      do i=1,ngr,40       
!       print*, rgrid(i),(cnew(1,j,i), j=1,4)  
!      enddo          
c      fast fourier transform
       dr=rgrid(2)-rgrid(1)      
       do j=1,nu
        do j1=1,nvv
         c2(1)=0
         do i=1,ngr-1 
          c2(i+1)=cnew(j,j1,i)*rgrid(i) 
         enddo
         call sinft(c2,ngr)
c        normalization of sin-fft with excluding the zeropoint (x=0)
         do i=1,ngr-1
          cf3(j,j1,i)=4*pi*c2(i+1)/kgrid(i)*dr 
         enddo
         cf3(j,j1,ngr)=cf3(j,j1,ngr-1) 
        enddo
       enddo       
!      do i=1,40       
!       print*, kgrid(i), (cf3(1,j1,i), j1=1,4) 
!      enddo        
       dk=kgrid(2)-kgrid(1)
       pi=2*asin(1.0)       
c       adding the long-range potential to c_short 
c
       do i=1,ngr
        do j=1,nu
         do j1=1,nvv
          cf3(j,j1,i)=cf3(j,j1,i)-1.0/kb/t*ul(j,j1,i)
         enddo
        enddo
       enddo 
c      calculations of fourier transform of h(i,j,k) by
c      site-site ornstein-zerinike equations
c       
       call ssoz(nvv,nu,ngr,nv,ims,mv,wu,cf3,chi,hu)
!      do i=1,ngr,40       
!       print*, (hu(1,j,i), j=1,4)  
!      enddo 
c      inverse fourier transform of h(i,j,k) & calculations of gamma(i,i,r) 
       do j=1,nu
        do j1=1,nvv
         ht(1)=0
         do i=1,ngr-1 
          ht(i+1)=hu(j,j1,i)*kgrid(i)         
         enddo
         call sinft(ht,ngr)
         do i=1,ngr-1
          hnew(j,j1,i)=ht(i+1)*dk/2/rgrid(i)/pi**2
          gnew(j,j1,i)=hnew(j,j1,i)-cnew(j,j1,i) 
         enddo   
         gnew(j,j1,ngr)=gnew(j,j1,ngr-1)
         hnew(j,j1,ngr)=hnew(j,j1,ngr-1) 
        enddo
       enddo 
!      do i=1,40       
!       print*, (gnew(1,j,i), j=1,4)  
!      enddo                      
c      intialization del
       do j=1,nu
        do j1=1,nvv
         del(j,j1)=0
        enddo
       mmax(j)=0 
       enddo
       mmaxi=0 
c      evaluation of accuracy
       do j=1,nu
        do j1=1,nvv
         do i=1,ngr
          del(j,j1)=del(j,j1)+(gnew(j,j1,i)-gold(j,j1,i))**2                                
         enddo
         mmax(j)=del(j,j1)+mmax(j)                 
        enddo
        mmaxi=mmaxi+mmax(j)  
       enddo  
!      do j=1,nu
!      print*, (del(j,j1), j1=1,4),mmmax  
!      enddo
c      normalization
       del0=(mmaxi/ngr/nu/nvv)**(1.0/2)
c      diis procedure
       call  diis(nu,nvv,ngr,lambd,kd,dd,gold,gnew,ggo,dgg)
c      old functions as new functions for new cycle
       do j=1,nu
        do j1=1,nvv
         do i=1,ngr
          gold(j,j1,i)=gnew(j,j1,i)
          cold(j,j1,i)=cnew(j,j1,i) 
         enddo    
        enddo
       enddo 
  107  format(9(1x,e10.4))
       k1=k1+1 
       print*, k1,kd,del0       
       if(k1.ge.700) goto 2         
        if(del0.ge.tol) goto 1
 2     continue  
       open(3, file='g1.data', status='old')  
        do k=1,ngr        
         write(3,107)(gnew(1,i,k)+cnew(1,i,k), i=1,4)  
        enddo
       close(3)
       return 
       end subroutine
c 
c      subroutine for evaluation diis residues dgg and matrix coefficients ss
c
       subroutine diis(nu,nvv,ngr,lam,k1,dd,hold,hnew,ggo,dgg)
       implicit none
       real*8 ggo(1:dd,1:nu,1:nvv,1:ngr),dgg(1:dd,1:nu,1:nvv,1:ngr)
       real*8 ss(1:dd+1,1:dd+1),hold(1:nu,1:nvv,1:ngr)
       real*8 s0(1:dd+1),s1,lam,hnew(1:nu,1:nvv,1:ngr)
       integer k1,dd,jd,jd1
       integer j,j1,i,nu,nvv,ngr
       integer ipiv(1:dd+1),info
c      initialization of overlap matrix for diis
       do jd=1,dd
        do jd1=1,dd
         ss(jd,jd1)=0
        enddo
        s0(jd)=0
        ss(jd,dd+1)=-1
        ss(dd+1,jd)=-1
       enddo
       ss(dd+1,dd+1)=0 
       s0(dd+1)=-1
c      calculation of basis
       do j=1,nu
        do j1=1,nvv
         do i=1,ngr     
          ggo(k1,j,j1,i)=hnew(j,j1,i) 
          dgg(k1,j,j1,i)=hnew(j,j1,i)-hold(j,j1,i)
c         calculation of elements of overlaping  matrix
          if(k1.eq.dd) then
           do jd=1,dd
            do jd1=jd,dd        
             ss(jd,jd1)=(dgg(jd,j,j1,i)*dgg(jd1,j,j1,i))/ngr+ss(jd,jd1)
c            symmerization
             ss(jd1,jd)=ss(jd,jd1)
            enddo
           enddo
          endif
         enddo
        enddo 
       enddo
c      solution of diis equation and change of baisis
       if(k1.eq.dd) then
        do jd=1,dd 
c        print*, (ss(jd,jd1), jd1=1,dd)
        enddo     
        call dgesv(dd+1,1,ss,dd+1,ipiv,s0,dd+1,info)            
        do i=1,ngr
         do j=1,nu
          do j1=1,nvv
           s1=0
c          calculation sum for new solution
           do jd=1,dd          
            s1=s1+s0(jd)*ggo(jd,j,j1,i)+lam*s0(jd)*dgg(jd,j,j1,i)  
           enddo
c          change of  diis basis
           do jd=1,dd-1           
            ggo(jd,j,j1,i)=ggo(jd+1,j,j1,i)
            dgg(jd,j,j1,i)=dgg(jd+1,j,j1,i)  
           enddo
           hnew(j,j1,i)=s1
          enddo
         enddo
        enddo 
       endif
c      matrix singular or has illegal value see dgeesv
       if((k1.eq.dd).and.(info.ne.0)) then 
        k1=1
        print*, info
        go to 1
       endif        
c      change count for diis                  
       if(k1.lt.dd) then
        k1=k1+1
        else 
        k1=dd
       endif  
 1     continue    
       return
       end subroutine
c
c uses realft
c calculates the sine transform of a set of n real-valued data points stored in array y(1:n).
c the number n must be a power of 2. on exit y is replaced by its transform. this program,
c without changes, also calculates the inverse sine transform, but in this case the output array
c should be multiplied by 2/n.
c
      subroutine sinft(y,n)
      integer n
      real*8 y(n)
      integer j
      real*8 sum,y1,y2
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n)
      wr=1.0d0
      wi=0.0d0
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      y(1)=0.0
      do 11 j=1,n/2
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
        y1=wi*(y(j+1)+y(n-j+1))
        y2=0.5*(y(j+1)-y(n-j+1))
        y(j+1)=y1+y2
        y(n-j+1)=y1-y2
11    continue
      call realft(y,n,+1)
      sum=0.0
      y(1)=0.5*y(1)
      y(2)=0.0
      do 12 j=1,n-1,2
        sum=sum+y(j)
        y(j)=y(j+1)
        y(j+1)=sum
12    continue
      return
      end subroutine
c
cc uses four1
c calculates the fourier transform of a set of n real-valued data points. replaces this data
c (which is stored in array data(1:n)) by the positive frequency half of its complex fourier
c transform. the real-valued rst and last components of the complex transform are returned
c as elements data(1) and data(2), respectively. n must be a power of 2. this routine
c also calculates the inverse transform of a complex data array if it is the transform of real
c data. (result in this case must be multiplied by 2/n.)
c
      subroutine realft(data,n,isign)
      integer isign,n
      real*8 data(n)
      integer i,i1,i2,i3,i4,n2p3
      real*8 c1,c2,h1i,h1r,h2i,h2r,wis,wrs
      double precision theta,wi,wpi,wpr,wr,wtemp
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if (isign.eq.1) then
        c2=-0.5
        call four1(data,n/2,+1)
      else
        c2=0.5
        theta=-theta
      endif
      wpr=-2.0d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.0d0+wpr
      wi=wpi
      n2p3=n+3
      do 11 i=2,n/4
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(data(i1)+data(i3))
        h1i=c1*(data(i2)-data(i4))
        h2r=-c2*(data(i2)+data(i4))
        h2i=c2*(data(i1)-data(i3))
        data(i1)=h1r+wrs*h2r-wis*h2i
        data(i2)=h1i+wrs*h2i+wis*h2r
        data(i3)=h1r-wrs*h2r+wis*h2i
        data(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=data(1)
        data(1)=h1r+data(2)
        data(2)=h1r-data(2)
      else
        h1r=data(1)
        data(1)=c1*(h1r+data(2))
        data(2)=c1*(h1r-data(2))
        call four1(data,n/2,-1)
      endif
      return
      end subroutine
cc replaces data(1:2*nn) by its discrete fourier transform, if isign is input as 1; or replaces
c data(1:2*nn) by nn times its inverse discrete fourier transform, if isign is input as -1.
c data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn
c must be an integer power of 2 (this is not checked for!).
      subroutine four1(data,nn,isign)
      integer isign,nn
      real*8 data(2*nn)
      integer i,istep,j,m,mmax,n
      real*8 tempi,tempr
      double precision theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*sin(0.5d0*theta)**2
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
            tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      goto 2
      endif
      return
      end subroutine
ccc lapack subroutines
ccc
      subroutine dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )
*
*  -- lapack driver routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*
*     .. scalar arguments ..
      integer            info, lda, ldb, n, nrhs
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  dgesv computes the solution to a real system of linear equations
*     a * x = b,
*  where a is an n-by-n matrix and x and b are n-by-nrhs matrices.
*
*  the lu decomposition with partial pivoting and row interchanges is
*  used to factor a as
*     a = p * l * u,
*  where p is a permutation matrix, l is unit lower triangular, and u is
*  upper triangular.  the factored form of a is then used to solve the
*  system of equations a * x = b.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of linear equations, i.e., the order of the
*          matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the n-by-n coefficient matrix a.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (output) integer array, dimension (n)
*          the pivot indices that define the permutation matrix p;
*          row i of the matrix was interchanged with row ipiv(i).
*
*  b       (input/output) double precision array, dimension (ldb,nrhs)
*          on entry, the n-by-nrhs matrix of right hand side matrix b.
*          on exit, if info = 0, the n-by-nrhs solution matrix x.
*
*  ldb     (input) integer
*          the leading dimension of the array b.  ldb >= max(1,n).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero.  the factorization
*                has been completed, but the factor u is exactly
*                singular, so the solution could not be computed.
*
*  =====================================================================
*
*     .. external subroutines ..
      external           dgetrf, dgetrs, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if( n.lt.0 ) then
         info = -1
      else if( nrhs.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, n ) ) then
         info = -4
      else if( ldb.lt.max( 1, n ) ) then
         info = -7
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgesv ', -info )
         return
      end if
*
*     compute the lu factorization of a.
*
      call dgetrf( n, n, a, lda, ipiv, info )
      if( info.eq.0 ) then
*
*        solve the system a*x = b, overwriting b with x.
*
         call dgetrs( 'no transpose', n, nrhs, a, lda, ipiv, b, ldb,
     $                info )
      end if
      return
*
*     end of dgesv
*
      end
ccc
      subroutine dgetrf( m, n, a, lda, ipiv, info )
*
*  -- lapack routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*
*     .. scalar arguments ..
      integer            info, lda, m, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dgetrf computes an lu factorization of a general m-by-n matrix a
*  using partial pivoting with row interchanges.
*
*  the factorization has the form
*     a = p * l * u
*  where p is a permutation matrix, l is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and u is upper
*  triangular (upper trapezoidal if m < n).
*
*  this is the right-looking level 3 blas version of the algorithm.
*
*  arguments
*  =========
*
*  m       (input) integer
*          the number of rows of the matrix a.  m >= 0.
*
*  n       (input) integer
*          the number of columns of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the m-by-n matrix to be factored.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,m).
*
*  ipiv    (output) integer array, dimension (min(m,n))
*          the pivot indices; for 1 <= i <= min(m,n), row i of the
*          matrix was interchanged with row ipiv(i).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*          > 0:  if info = i, u(i,i) is exactly zero. the factorization
*                has been completed, but the factor u is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. local scalars ..
      integer            i, iinfo, j, jb, nb
*     ..
*     .. external subroutines ..
      external           dgemm, dgetf2, dlaswp, dtrsm, xerbla
*     ..
*     .. external functions ..
      integer            ilaenv
      external           ilaenv
*     ..
*     .. intrinsic functions ..
      intrinsic          max, min
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetrf', -info )
         return
      end if
*
*     quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
*     determine the block size for this environment.
*
      nb = ilaenv( 1, 'dgetrf', ' ', m, n, -1, -1 )
      if( nb.le.1 .or. nb.ge.min( m, n ) ) then
*
*        use unblocked code.
*
         call dgetf2( m, n, a, lda, ipiv, info )
      else
*
*        use blocked code.
*
         do 20 j = 1, min( m, n ), nb
            jb = min( min( m, n )-j+1, nb )
*
*           factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            call dgetf2( m-j+1, jb, a( j, j ), lda, ipiv( j ), iinfo )
*
*           adjust info and the pivot indices.
*
            if( info.eq.0 .and. iinfo.gt.0 )
     $         info = iinfo + j - 1
            do 10 i = j, min( m, j+jb-1 )
               ipiv( i ) = j - 1 + ipiv( i )
   10       continue
*
*           apply interchanges to columns 1:j-1.
*
            call dlaswp( j-1, a, lda, j, j+jb-1, ipiv, 1 )
*
            if( j+jb.le.n ) then
*
*              apply interchanges to columns j+jb:n.
*
               call dlaswp( n-j-jb+1, a( 1, j+jb ), lda, j, j+jb-1,
     $                      ipiv, 1 )
*
*              compute block row of u.
*
               call dtrsm( 'left', 'lower', 'no transpose', 'unit', jb,
     $                     n-j-jb+1, one, a( j, j ), lda, a( j, j+jb ),
     $                     lda )
               if( j+jb.le.m ) then
*
*                 update trailing submatrix.
*
                  call dgemm( 'no transpose', 'no transpose', m-j-jb+1,
     $                        n-j-jb+1, jb, -one, a( j+jb, j ), lda,
     $                        a( j, j+jb ), lda, one, a( j+jb, j+jb ),
     $                        lda )
               end if
            end if
   20    continue
      end if
      return
*
*     end of dgetrf
*
      end
ccc
      subroutine dgetf2( m, n, a, lda, ipiv, info )
*
*  -- lapack routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*
*     .. scalar arguments ..
      integer            info, lda, m, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dgetf2 computes an lu factorization of a general m-by-n matrix a
*  using partial pivoting with row interchanges.
*
*  the factorization has the form
*     a = p * l * u
*  where p is a permutation matrix, l is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and u is upper
*  triangular (upper trapezoidal if m < n).
*
*  this is the right-looking level 2 blas version of the algorithm.
*
*  arguments
*  =========
*
*  m       (input) integer
*          the number of rows of the matrix a.  m >= 0.
*
*  n       (input) integer
*          the number of columns of the matrix a.  n >= 0.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the m by n matrix to be factored.
*          on exit, the factors l and u from the factorization
*          a = p*l*u; the unit diagonal elements of l are not stored.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,m).
*
*  ipiv    (output) integer array, dimension (min(m,n))
*          the pivot indices; for 1 <= i <= min(m,n), row i of the
*          matrix was interchanged with row ipiv(i).
*
*  info    (output) integer
*          = 0: successful exit
*          < 0: if info = -k, the k-th argument had an illegal value
*          > 0: if info = k, u(k,k) is exactly zero. the factorization
*               has been completed, but the factor u is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      double precision   sfmin 
      integer            i, j, jp
*     ..
*     .. external functions ..
      double precision   dlamch      
      integer            idamax
      external           dlamch, idamax
*     ..
*     .. external subroutines ..
      external           dger, dscal, dswap, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max, min
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      if( m.lt.0 ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( lda.lt.max( 1, m ) ) then
         info = -4
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetf2', -info )
         return
      end if
*
*     quick return if possible
*
      if( m.eq.0 .or. n.eq.0 )
     $   return
*
*     compute machine safe minimum 
* 
      sfmin = dlamch('s')  
*
      do 10 j = 1, min( m, n )
*
*        find pivot and test for singularity.
*
         jp = j - 1 + idamax( m-j+1, a( j, j ), 1 )
         ipiv( j ) = jp
         if( a( jp, j ).ne.zero ) then
*
*           apply the interchange to columns 1:n.
*
            if( jp.ne.j )
     $         call dswap( n, a( j, 1 ), lda, a( jp, 1 ), lda )
*
*           compute elements j+1:m of j-th column.
*
            if( j.lt.m ) then 
               if( abs(a( j, j )) .ge. sfmin ) then 
                  call dscal( m-j, one / a( j, j ), a( j+1, j ), 1 ) 
               else 
                 do 20 i = 1, m-j 
                    a( j+i, j ) = a( j+i, j ) / a( j, j ) 
   20            continue 
               end if 
            end if 
*
         else if( info.eq.0 ) then
*
            info = j
         end if
*
         if( j.lt.min( m, n ) ) then
*
*           update trailing submatrix.
*
            call dger( m-j, n-j, -one, a( j+1, j ), 1, a( j, j+1 ), lda,
     $                 a( j+1, j+1 ), lda )
         end if
   10 continue
      return
*
*     end of dgetf2
*
      end
ccc
      subroutine dger(m,n,alpha,x,incx,y,incy,a,lda)
*     .. scalar arguments ..
      double precision alpha
      integer incx,incy,lda,m,n
*     ..
*     .. array arguments ..
      double precision a(lda,*),x(*),y(*)
*     ..
*
*  purpose
*  =======
*
*  dger   performs the rank 1 operation
*
*     a := alpha*x*y' + a,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and a is an m by n matrix.
*
*  arguments
*  ==========
*
*  m      - integer.
*           on entry, m specifies the number of rows of the matrix a.
*           m must be at least zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry, n specifies the number of columns of the matrix a.
*           n must be at least zero.
*           unchanged on exit.
*
*  alpha  - double precision.
*           on entry, alpha specifies the scalar alpha.
*           unchanged on exit.
*
*  x      - double precision array of dimension at least
*           ( 1 + ( m - 1 )*abs( incx ) ).
*           before entry, the incremented array x must contain the m
*           element vector x.
*           unchanged on exit.
*
*  incx   - integer.
*           on entry, incx specifies the increment for the elements of
*           x. incx must not be zero.
*           unchanged on exit.
*
*  y      - double precision array of dimension at least
*           ( 1 + ( n - 1 )*abs( incy ) ).
*           before entry, the incremented array y must contain the n
*           element vector y.
*           unchanged on exit.
*
*  incy   - integer.
*           on entry, incy specifies the increment for the elements of
*           y. incy must not be zero.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, n ).
*           before entry, the leading m by n part of the array a must
*           contain the matrix of coefficients. on exit, a is
*           overwritten by the updated matrix.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. lda must be at least
*           max( 1, m ).
*           unchanged on exit.
*
*  further details
*  ===============
*
*  level 2 blas routine.
*
*  -- written on 22-october-1986.
*     jack dongarra, argonne national lab.
*     jeremy du croz, nag central office.
*     sven hammarling, nag central office.
*     richard hanson, sandia national labs.
*
*  =====================================================================
*
*     .. parameters ..
      double precision zero
      parameter (zero=0.0d+0)
*     ..
*     .. local scalars ..
      double precision temp
      integer i,info,ix,j,jy,kx
*     ..
*     .. external subroutines ..
      external xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic max
*     ..
*
*     test the input parameters.
*
      info = 0
      if (m.lt.0) then
          info = 1
      else if (n.lt.0) then
          info = 2
      else if (incx.eq.0) then
          info = 5
      else if (incy.eq.0) then
          info = 7
      else if (lda.lt.max(1,m)) then
          info = 9
      end if
      if (info.ne.0) then
          call xerbla('dger  ',info)
          return
      end if
*
*     quick return if possible.
*
      if ((m.eq.0) .or. (n.eq.0) .or. (alpha.eq.zero)) return
*
*     start the operations. in this version the elements of a are
*     accessed sequentially with one pass through a.
*
      if (incy.gt.0) then
          jy = 1
      else
          jy = 1 - (n-1)*incy
      end if
      if (incx.eq.1) then
          do 20 j = 1,n
              if (y(jy).ne.zero) then
                  temp = alpha*y(jy)
                  do 10 i = 1,m
                      a(i,j) = a(i,j) + x(i)*temp
   10             continue
              end if
              jy = jy + incy
   20     continue
      else
          if (incx.gt.0) then
              kx = 1
          else
              kx = 1 - (m-1)*incx
          end if
          do 40 j = 1,n
              if (y(jy).ne.zero) then
                  temp = alpha*y(jy)
                  ix = kx
                  do 30 i = 1,m
                      a(i,j) = a(i,j) + x(ix)*temp
                      ix = ix + incx
   30             continue
              end if
              jy = jy + incy
   40     continue
      end if
*
      return
*
*     end of dger  .
*
      end
ccc
      subroutine dscal(n,da,dx,incx)
*     .. scalar arguments ..
      double precision da
      integer incx,n
*     ..
*     .. array arguments ..
      double precision dx(*)
*     ..
*
*  purpose
*  =======
*
*     dscal scales a vector by a constant.
*     uses unrolled loops for increment equal to one.
*
*  further details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      integer i,m,mp1,nincx
*     ..
*     .. intrinsic functions ..
      intrinsic mod
*     ..
      if (n.le.0 .or. incx.le.0) return
      if (incx.eq.1) go to 20
*
*        code for increment not equal to 1
*
      nincx = n*incx
      do 10 i = 1,nincx,incx
          dx(i) = da*dx(i)
   10 continue
      return
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
   20 m = mod(n,5)
      if (m.eq.0) go to 40
      do 30 i = 1,m
          dx(i) = da*dx(i)
   30 continue
      if (n.lt.5) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
          dx(i) = da*dx(i)
          dx(i+1) = da*dx(i+1)
          dx(i+2) = da*dx(i+2)
          dx(i+3) = da*dx(i+3)
          dx(i+4) = da*dx(i+4)
   50 continue
      return
      end
ccc
      subroutine dswap(n,dx,incx,dy,incy)
*     .. scalar arguments ..
      integer incx,incy,n
*     ..
*     .. array arguments ..
      double precision dx(*),dy(*)
*     ..
*
*  purpose
*  =======
*
*     interchanges two vectors.
*     uses unrolled loops for increments equal one.
*
*  further details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      double precision dtemp
      integer i,ix,iy,m,mp1
*     ..
*     .. intrinsic functions ..
      intrinsic mod
*     ..
      if (n.le.0) return
      if (incx.eq.1 .and. incy.eq.1) go to 20
*
*       code for unequal increments or equal increments not equal
*         to 1
*
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
          dtemp = dx(ix)
          dx(ix) = dy(iy)
          dy(iy) = dtemp
          ix = ix + incx
          iy = iy + incy
   10 continue
      return
*
*       code for both increments equal to 1
*
*
*       clean-up loop
*
   20 m = mod(n,3)
      if (m.eq.0) go to 40
      do 30 i = 1,m
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
   30 continue
      if (n.lt.3) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
          dtemp = dx(i)
          dx(i) = dy(i)
          dy(i) = dtemp
          dtemp = dx(i+1)
          dx(i+1) = dy(i+1)
          dy(i+1) = dtemp
          dtemp = dx(i+2)
          dx(i+2) = dy(i+2)
          dy(i+2) = dtemp
   50 continue
      return
      end
ccc
      subroutine xerbla( srname, info )
*
*  -- lapack auxiliary routine (preliminary version) --
*     univ. of tennessee, univ. of california berkeley and nag ltd..
*     november 2006
*
*     .. scalar arguments ..
      character*(*)      srname
      integer            info
*     ..
*
*  purpose
*  =======
*
*  xerbla  is an error handler for the lapack routines.
*  it is called by an lapack routine if an input parameter has an
*  invalid value.  a message is printed and execution stops.
*
*  installers may consider modifying the stop statement in order to
*  call system-specific exception-handling facilities.
*
*  arguments
*  =========
*
*  srname  (input) character*(*)
*          the name of the routine which called xerbla.
*
*  info    (input) integer
*          the position of the invalid parameter in the parameter list
*          of the calling routine.
*
* =====================================================================
*
*     .. intrinsic functions ..
      intrinsic          len_trim
*     ..
*     .. executable statements ..
*
      write( *, fmt = 9999 )srname( 1:len_trim( srname ) ), info
*
      stop
*
 9999 format( ' ** on entry to ', a, ' parameter number ', i2, ' had ',
     $      'an illegal value' )
*
*     end of xerbla
*
      end
ccc
      subroutine dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
*
*  -- lapack routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*
*     .. scalar arguments ..
      character          trans
      integer            info, lda, ldb, n, nrhs
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * ), b( ldb, * )
*     ..
*
*  purpose
*  =======
*
*  dgetrs solves a system of linear equations
*     a * x = b  or  a' * x = b
*  with a general n-by-n matrix a using the lu factorization computed
*  by dgetrf.
*
*  arguments
*  =========
*
*  trans   (input) character*1
*          specifies the form of the system of equations:
*          = 'n':  a * x = b  (no transpose)
*          = 't':  a'* x = b  (transpose)
*          = 'c':  a'* x = b  (conjugate transpose = transpose)
*
*  n       (input) integer
*          the order of the matrix a.  n >= 0.
*
*  nrhs    (input) integer
*          the number of right hand sides, i.e., the number of columns
*          of the matrix b.  nrhs >= 0.
*
*  a       (input) double precision array, dimension (lda,n)
*          the factors l and u from the factorization a = p*l*u
*          as computed by dgetrf.
*
*  lda     (input) integer
*          the leading dimension of the array a.  lda >= max(1,n).
*
*  ipiv    (input) integer array, dimension (n)
*          the pivot indices from dgetrf; for 1<=i<=n, row i of the
*          matrix was interchanged with row ipiv(i).
*
*  b       (input/output) double precision array, dimension (ldb,nrhs)
*          on entry, the right hand side matrix b.
*          on exit, the solution matrix x.
*
*  ldb     (input) integer
*          the leading dimension of the array b.  ldb >= max(1,n).
*
*  info    (output) integer
*          = 0:  successful exit
*          < 0:  if info = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. parameters ..
      double precision   one
      parameter          ( one = 1.0d+0 )
*     ..
*     .. local scalars ..
      logical            notran
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. external subroutines ..
      external           dlaswp, dtrsm, xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic          max
*     ..
*     .. executable statements ..
*
*     test the input parameters.
*
      info = 0
      notran = lsame( trans, 'n' )
      if( .not.notran .and. .not.lsame( trans, 't' ) .and. .not.
     $    lsame( trans, 'c' ) ) then
         info = -1
      else if( n.lt.0 ) then
         info = -2
      else if( nrhs.lt.0 ) then
         info = -3
      else if( lda.lt.max( 1, n ) ) then
         info = -5
      else if( ldb.lt.max( 1, n ) ) then
         info = -8
      end if
      if( info.ne.0 ) then
         call xerbla( 'dgetrs', -info )
         return
      end if
*
*     quick return if possible
*
      if( n.eq.0 .or. nrhs.eq.0 )
     $   return
*
      if( notran ) then
*
*        solve a * x = b.
*
*        apply row interchanges to the right hand sides.
*
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, 1 )
*
*        solve l*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'lower', 'no transpose', 'unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        solve u*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'upper', 'no transpose', 'non-unit', n,
     $               nrhs, one, a, lda, b, ldb )
      else
*
*        solve a' * x = b.
*
*        solve u'*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'upper', 'transpose', 'non-unit', n, nrhs,
     $               one, a, lda, b, ldb )
*
*        solve l'*x = b, overwriting b with x.
*
         call dtrsm( 'left', 'lower', 'transpose', 'unit', n, nrhs, one,
     $               a, lda, b, ldb )
*
*        apply row interchanges to the solution vectors.
*
         call dlaswp( nrhs, b, ldb, 1, n, ipiv, -1 )
      end if
*
      return
*
*     end of dgetrs
*
      end
ccc
      subroutine dgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
*     .. scalar arguments ..
      double precision alpha,beta
      integer k,lda,ldb,ldc,m,n
      character transa,transb
*     ..
*     .. array arguments ..
      double precision a(lda,*),b(ldb,*),c(ldc,*)
*     ..
*
*  purpose
*  =======
*
*  dgemm  performs one of the matrix-matrix operations
*
*     c := alpha*op( a )*op( b ) + beta*c,
*
*  where  op( x ) is one of
*
*     op( x ) = x   or   op( x ) = x',
*
*  alpha and beta are scalars, and a, b and c are matrices, with op( a )
*  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
*
*  arguments
*  ==========
*
*  transa - character*1.
*           on entry, transa specifies the form of op( a ) to be used in
*           the matrix multiplication as follows:
*
*              transa = 'n' or 'n',  op( a ) = a.
*
*              transa = 't' or 't',  op( a ) = a'.
*
*              transa = 'c' or 'c',  op( a ) = a'.
*
*           unchanged on exit.
*
*  transb - character*1.
*           on entry, transb specifies the form of op( b ) to be used in
*           the matrix multiplication as follows:
*
*              transb = 'n' or 'n',  op( b ) = b.
*
*              transb = 't' or 't',  op( b ) = b'.
*
*              transb = 'c' or 'c',  op( b ) = b'.
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry,  m  specifies  the number  of rows  of the  matrix
*           op( a )  and of the  matrix  c.  m  must  be at least  zero.
*           unchanged on exit.
*
*  n      - integer.
*           on entry,  n  specifies the number  of columns of the matrix
*           op( b ) and the number of columns of the matrix c. n must be
*           at least zero.
*           unchanged on exit.
*
*  k      - integer.
*           on entry,  k  specifies  the number of columns of the matrix
*           op( a ) and the number of rows of the matrix op( b ). k must
*           be at least  zero.
*           unchanged on exit.
*
*  alpha  - double precision.
*           on entry, alpha specifies the scalar alpha.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, ka ), where ka is
*           k  when  transa = 'n' or 'n',  and is  m  otherwise.
*           before entry with  transa = 'n' or 'n',  the leading  m by k
*           part of the array  a  must contain the matrix  a,  otherwise
*           the leading  k by m  part of the array  a  must contain  the
*           matrix a.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program. when  transa = 'n' or 'n' then
*           lda must be at least  max( 1, m ), otherwise  lda must be at
*           least  max( 1, k ).
*           unchanged on exit.
*
*  b      - double precision array of dimension ( ldb, kb ), where kb is
*           n  when  transb = 'n' or 'n',  and is  k  otherwise.
*           before entry with  transb = 'n' or 'n',  the leading  k by n
*           part of the array  b  must contain the matrix  b,  otherwise
*           the leading  n by k  part of the array  b  must contain  the
*           matrix b.
*           unchanged on exit.
*
*  ldb    - integer.
*           on entry, ldb specifies the first dimension of b as declared
*           in the calling (sub) program. when  transb = 'n' or 'n' then
*           ldb must be at least  max( 1, k ), otherwise  ldb must be at
*           least  max( 1, n ).
*           unchanged on exit.
*
*  beta   - double precision.
*           on entry,  beta  specifies the scalar  beta.  when  beta  is
*           supplied as zero then c need not be set on input.
*           unchanged on exit.
*
*  c      - double precision array of dimension ( ldc, n ).
*           before entry, the leading  m by n  part of the array  c must
*           contain the matrix  c,  except when  beta  is zero, in which
*           case c need not be set on entry.
*           on exit, the array  c  is overwritten by the  m by n  matrix
*           ( alpha*op( a )*op( b ) + beta*c ).
*
*  ldc    - integer.
*           on entry, ldc specifies the first dimension of c as declared
*           in  the  calling  (sub)  program.   ldc  must  be  at  least
*           max( 1, m ).
*           unchanged on exit.
*
*  further details
*  ===============
*
*  level 3 blas routine.
*
*  -- written on 8-february-1989.
*     jack dongarra, argonne national laboratory.
*     iain duff, aere harwell.
*     jeremy du croz, numerical algorithms group ltd.
*     sven hammarling, numerical algorithms group ltd.
*
*  =====================================================================
*
*     .. external functions ..
      logical lsame
      external lsame
*     ..
*     .. external subroutines ..
      external xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic max
*     ..
*     .. local scalars ..
      double precision temp
      integer i,info,j,l,ncola,nrowa,nrowb
      logical nota,notb
*     ..
*     .. parameters ..
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
*     ..
*
*     set  nota  and  notb  as  true if  a  and  b  respectively are not
*     transposed and set  nrowa, ncola and  nrowb  as the number of rows
*     and  columns of  a  and the  number of  rows  of  b  respectively.
*
      nota = lsame(transa,'n')
      notb = lsame(transb,'n')
      if (nota) then
          nrowa = m
          ncola = k
      else
          nrowa = k
          ncola = m
      end if
      if (notb) then
          nrowb = k
      else
          nrowb = n
      end if
*
*     test the input parameters.
*
      info = 0
      if ((.not.nota) .and. (.not.lsame(transa,'c')) .and.
     +    (.not.lsame(transa,'t'))) then
          info = 1
      else if ((.not.notb) .and. (.not.lsame(transb,'c')) .and.
     +         (.not.lsame(transb,'t'))) then
          info = 2
      else if (m.lt.0) then
          info = 3
      else if (n.lt.0) then
          info = 4
      else if (k.lt.0) then
          info = 5
      else if (lda.lt.max(1,nrowa)) then
          info = 8
      else if (ldb.lt.max(1,nrowb)) then
          info = 10
      else if (ldc.lt.max(1,m)) then
          info = 13
      end if
      if (info.ne.0) then
          call xerbla('dgemm ',info)
          return
      end if
*
*     quick return if possible.
*
      if ((m.eq.0) .or. (n.eq.0) .or.
     +    (((alpha.eq.zero).or. (k.eq.0)).and. (beta.eq.one))) return
*
*     and if  alpha.eq.zero.
*
      if (alpha.eq.zero) then
          if (beta.eq.zero) then
              do 20 j = 1,n
                  do 10 i = 1,m
                      c(i,j) = zero
   10             continue
   20         continue
          else
              do 40 j = 1,n
                  do 30 i = 1,m
                      c(i,j) = beta*c(i,j)
   30             continue
   40         continue
          end if
          return
      end if
*
*     start the operations.
*
      if (notb) then
          if (nota) then
*
*           form  c := alpha*a*b + beta*c.
*
              do 90 j = 1,n
                  if (beta.eq.zero) then
                      do 50 i = 1,m
                          c(i,j) = zero
   50                 continue
                  else if (beta.ne.one) then
                      do 60 i = 1,m
                          c(i,j) = beta*c(i,j)
   60                 continue
                  end if
                  do 80 l = 1,k
                      if (b(l,j).ne.zero) then
                          temp = alpha*b(l,j)
                          do 70 i = 1,m
                              c(i,j) = c(i,j) + temp*a(i,l)
   70                     continue
                      end if
   80             continue
   90         continue
          else
*
*           form  c := alpha*a'*b + beta*c
*
              do 120 j = 1,n
                  do 110 i = 1,m
                      temp = zero
                      do 100 l = 1,k
                          temp = temp + a(l,i)*b(l,j)
  100                 continue
                      if (beta.eq.zero) then
                          c(i,j) = alpha*temp
                      else
                          c(i,j) = alpha*temp + beta*c(i,j)
                      end if
  110             continue
  120         continue
          end if
      else
          if (nota) then
*
*           form  c := alpha*a*b' + beta*c
*
              do 170 j = 1,n
                  if (beta.eq.zero) then
                      do 130 i = 1,m
                          c(i,j) = zero
  130                 continue
                  else if (beta.ne.one) then
                      do 140 i = 1,m
                          c(i,j) = beta*c(i,j)
  140                 continue
                  end if
                  do 160 l = 1,k
                      if (b(j,l).ne.zero) then
                          temp = alpha*b(j,l)
                          do 150 i = 1,m
                              c(i,j) = c(i,j) + temp*a(i,l)
  150                     continue
                      end if
  160             continue
  170         continue
          else
*
*           form  c := alpha*a'*b' + beta*c
*
              do 200 j = 1,n
                  do 190 i = 1,m
                      temp = zero
                      do 180 l = 1,k
                          temp = temp + a(l,i)*b(j,l)
  180                 continue
                      if (beta.eq.zero) then
                          c(i,j) = alpha*temp
                      else
                          c(i,j) = alpha*temp + beta*c(i,j)
                      end if
  190             continue
  200         continue
          end if
      end if
*
      return
*
*     end of dgemm .
*
      end
ccc
      subroutine dtrsm(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
*     .. scalar arguments ..
      double precision alpha
      integer lda,ldb,m,n
      character diag,side,transa,uplo
*     ..
*     .. array arguments ..
      double precision a(lda,*),b(ldb,*)
*     ..
*
*  purpose
*  =======
*
*  dtrsm  solves one of the matrix equations
*
*     op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,
*
*  where alpha is a scalar, x and b are m by n matrices, a is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
*
*     op( a ) = a   or   op( a ) = a'.
*
*  the matrix x is overwritten on b.
*
*  arguments
*  ==========
*
*  side   - character*1.
*           on entry, side specifies whether op( a ) appears on the left
*           or right of x as follows:
*
*              side = 'l' or 'l'   op( a )*x = alpha*b.
*
*              side = 'r' or 'r'   x*op( a ) = alpha*b.
*
*           unchanged on exit.
*
*  uplo   - character*1.
*           on entry, uplo specifies whether the matrix a is an upper or
*           lower triangular matrix as follows:
*
*              uplo = 'u' or 'u'   a is an upper triangular matrix.
*
*              uplo = 'l' or 'l'   a is a lower triangular matrix.
*
*           unchanged on exit.
*
*  transa - character*1.
*           on entry, transa specifies the form of op( a ) to be used in
*           the matrix multiplication as follows:
*
*              transa = 'n' or 'n'   op( a ) = a.
*
*              transa = 't' or 't'   op( a ) = a'.
*
*              transa = 'c' or 'c'   op( a ) = a'.
*
*           unchanged on exit.
*
*  diag   - character*1.
*           on entry, diag specifies whether or not a is unit triangular
*           as follows:
*
*              diag = 'u' or 'u'   a is assumed to be unit triangular.
*
*              diag = 'n' or 'n'   a is not assumed to be unit
*                                  triangular.
*
*           unchanged on exit.
*
*  m      - integer.
*           on entry, m specifies the number of rows of b. m must be at
*           least zero.
*           unchanged on exit.gchuev-874 
*
*  n      - integer.
*           on entry, n specifies the number of columns of b.  n must be
*           at least zero.
*           unchanged on exit.
*
*  alpha  - double precision.
*           on entry,  alpha specifies the scalar  alpha. when  alpha is
*           zero then  a is not referenced and  b need not be set before
*           entry.
*           unchanged on exit.
*
*  a      - double precision array of dimension ( lda, k ), where k is m
*           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
*           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
*           upper triangular part of the array  a must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           a is not referenced.
*           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
*           lower triangular part of the array  a must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           a is not referenced.
*           note that when  diag = 'u' or 'u',  the diagonal elements of
*           a  are not referenced either,  but are assumed to be  unity.
*           unchanged on exit.
*
*  lda    - integer.
*           on entry, lda specifies the first dimension of a as declared
*           in the calling (sub) program.  when  side = 'l' or 'l'  then
*           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
*           then lda must be at least max( 1, n ).
*           unchanged on exit.
*
*  b      - double precision array of dimension ( ldb, n ).
*           before entry,  the leading  m by n part of the array  b must
*           contain  the  right-hand  side  matrix  b,  and  on exit  is
*           overwritten by the solution matrix  x.
*
*  ldb    - integer.
*           on entry, ldb specifies the first dimension of b as declared
*           in  the  calling  (sub)  program.   ldb  must  be  at  least
*           max( 1, m ).
*           unchanged on exit.
*
*  further details
*  ===============
*
*  level 3 blas routine.
*
*
*  -- written on 8-february-1989.
*     jack dongarra, argonne national laboratory.
*     iain duff, aere harwell.
*     jeremy du croz, numerical algorithms group ltd.
*     sven hammarling, numerical algorithms group ltd.
*
*  =====================================================================
*
*     .. external functions ..
      logical lsame
      external lsame
*     ..
*     .. external subroutines ..
      external xerbla
*     ..
*     .. intrinsic functions ..
      intrinsic max
*     ..
*     .. local scalars ..
      double precision temp
      integer i,info,j,k,nrowa
      logical lside,nounit,upper
*     ..
*     .. parameters ..
      double precision one,zero
      parameter (one=1.0d+0,zero=0.0d+0)
*     ..
*
*     test the input parameters.
*
      lside = lsame(side,'l')
      if (lside) then
          nrowa = m
      else
          nrowa = n
      end if
      nounit = lsame(diag,'n')
      upper = lsame(uplo,'u')
*
      info = 0
      if ((.not.lside) .and. (.not.lsame(side,'r'))) then
          info = 1
      else if ((.not.upper) .and. (.not.lsame(uplo,'l'))) then
          info = 2
      else if ((.not.lsame(transa,'n')) .and.
     +         (.not.lsame(transa,'t')) .and.
     +         (.not.lsame(transa,'c'))) then
          info = 3
      else if ((.not.lsame(diag,'u')) .and. (.not.lsame(diag,'n'))) then
          info = 4
      else if (m.lt.0) then
          info = 5
      else if (n.lt.0) then
          info = 6
      else if (lda.lt.max(1,nrowa)) then
          info = 9
      else if (ldb.lt.max(1,m)) then
          info = 11
      end if
      if (info.ne.0) then
          call xerbla('dtrsm ',info)
          return
      end if
*
*     quick return if possible.
*
      if (m.eq.0 .or. n.eq.0) return
*
*     and when  alpha.eq.zero.
*
      if (alpha.eq.zero) then
          do 20 j = 1,n
              do 10 i = 1,m
                  b(i,j) = zero
   10         continue
   20     continue
          return
      end if
*
*     start the operations.
*
      if (lside) then
          if (lsame(transa,'n')) then
*
*           form  b := alpha*inv( a )*b.
*
              if (upper) then
                  do 60 j = 1,n
                      if (alpha.ne.one) then
                          do 30 i = 1,m
                              b(i,j) = alpha*b(i,j)
   30                     continue
                      end if
                      do 50 k = m,1,-1
                          if (b(k,j).ne.zero) then
                              if (nounit) b(k,j) = b(k,j)/a(k,k)
                              do 40 i = 1,k - 1
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   40                         continue
                          end if
   50                 continue
   60             continue
              else
                  do 100 j = 1,n
                      if (alpha.ne.one) then
                          do 70 i = 1,m
                              b(i,j) = alpha*b(i,j)
   70                     continue
                      end if
                      do 90 k = 1,m
                          if (b(k,j).ne.zero) then
                              if (nounit) b(k,j) = b(k,j)/a(k,k)
                              do 80 i = k + 1,m
                                  b(i,j) = b(i,j) - b(k,j)*a(i,k)
   80                         continue
                          end if
   90                 continue
  100             continue
              end if
          else
*
*           form  b := alpha*inv( a' )*b.
*
              if (upper) then
                  do 130 j = 1,n
                      do 120 i = 1,m
                          temp = alpha*b(i,j)
                          do 110 k = 1,i - 1
                              temp = temp - a(k,i)*b(k,j)
  110                     continue
                          if (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  120                 continue
  130             continue
              else
                  do 160 j = 1,n
                      do 150 i = m,1,-1
                          temp = alpha*b(i,j)
                          do 140 k = i + 1,m
                              temp = temp - a(k,i)*b(k,j)
  140                     continue
                          if (nounit) temp = temp/a(i,i)
                          b(i,j) = temp
  150                 continue
  160             continue
              end if
          end if
      else
          if (lsame(transa,'n')) then
*
*           form  b := alpha*b*inv( a ).
*
              if (upper) then
                  do 210 j = 1,n
                      if (alpha.ne.one) then
                          do 170 i = 1,m
                              b(i,j) = alpha*b(i,j)
  170                     continue
                      end if
                      do 190 k = 1,j - 1
                          if (a(k,j).ne.zero) then
                              do 180 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  180                         continue
                          end if
  190                 continue
                      if (nounit) then
                          temp = one/a(j,j)
                          do 200 i = 1,m
                              b(i,j) = temp*b(i,j)
  200                     continue
                      end if
  210             continue
              else
                  do 260 j = n,1,-1
                      if (alpha.ne.one) then
                          do 220 i = 1,m
                              b(i,j) = alpha*b(i,j)
  220                     continue
                      end if
                      do 240 k = j + 1,n
                          if (a(k,j).ne.zero) then
                              do 230 i = 1,m
                                  b(i,j) = b(i,j) - a(k,j)*b(i,k)
  230                         continue
                          end if
  240                 continue
                      if (nounit) then
                          temp = one/a(j,j)
                          do 250 i = 1,m
                              b(i,j) = temp*b(i,j)
  250                     continue
                      end if
  260             continue
              end if
          else
*
*           form  b := alpha*b*inv( a' ).
*
              if (upper) then
                  do 310 k = n,1,-1
                      if (nounit) then
                          temp = one/a(k,k)
                          do 270 i = 1,m
                              b(i,k) = temp*b(i,k)
  270                     continue
                      end if
                      do 290 j = 1,k - 1
                          if (a(j,k).ne.zero) then
                              temp = a(j,k)
                              do 280 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  280                         continue
                          end if
  290                 continue
                      if (alpha.ne.one) then
                          do 300 i = 1,m
                              b(i,k) = alpha*b(i,k)
  300                     continue
                      end if
  310             continue
              else
                  do 360 k = 1,n
                      if (nounit) then
                          temp = one/a(k,k)
                          do 320 i = 1,m
                              b(i,k) = temp*b(i,k)
  320                     continue
                      end if
                      do 340 j = k + 1,n
                          if (a(j,k).ne.zero) then
                              temp = a(j,k)
                              do 330 i = 1,m
                                  b(i,j) = b(i,j) - temp*b(i,k)
  330                         continue
                          end if
  340                 continue
                      if (alpha.ne.one) then
                          do 350 i = 1,m
                              b(i,k) = alpha*b(i,k)
  350                     continue
                      end if
  360             continue
              end if
          end if
      end if
*
      return
*
*     end of dtrsm .
*
      end
ccc
      subroutine dlaswp( n, a, lda, k1, k2, ipiv, incx )
*
*  -- lapack auxiliary routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*
*     .. scalar arguments ..
      integer            incx, k1, k2, lda, n
*     ..
*     .. array arguments ..
      integer            ipiv( * )
      double precision   a( lda, * )
*     ..
*
*  purpose
*  =======
*
*  dlaswp performs a series of row interchanges on the matrix a.
*  one row interchange is initiated for each of rows k1 through k2 of a.
*
*  arguments
*  =========
*
*  n       (input) integer
*          the number of columns of the matrix a.
*
*  a       (input/output) double precision array, dimension (lda,n)
*          on entry, the matrix of column dimension n to which the row
*          interchanges will be applied.
*          on exit, the permuted matrix.
*
*  lda     (input) integer
*          the leading dimension of the array a.
*
*  k1      (input) integer
*          the first element of ipiv for which a row interchange will
*          be done.
*
*  k2      (input) integer
*          the last element of ipiv for which a row interchange will
*          be done.
*
*  ipiv    (input) integer array, dimension (k2*abs(incx))
*          the vector of pivot indices.  only the elements in positions
*          k1 through k2 of ipiv are accessed.
*          ipiv(k) = l implies rows k and l are to be interchanged.
*
*  incx    (input) integer
*          the increment between successive values of ipiv.  if ipiv
*          is negative, the pivots are applied in reverse order.
*
*  further details
*  ===============
*
*  modified by
*   r. c. whaley, computer science dept., univ. of tenn., knoxville, usa
*
* =====================================================================
*
*     .. local scalars ..
      integer            i, i1, i2, inc, ip, ix, ix0, j, k, n32
      double precision   temp
*     ..
*     .. executable statements ..
*
*     interchange row i with row ipiv(i) for each of rows k1 through k2.
*
      if( incx.gt.0 ) then
         ix0 = k1
         i1 = k1
         i2 = k2
         inc = 1
      else if( incx.lt.0 ) then
         ix0 = 1 + ( 1-k2 )*incx
         i1 = k2
         i2 = k1
         inc = -1
      else
         return
      end if
*
      n32 = ( n / 32 )*32
      if( n32.ne.0 ) then
         do 30 j = 1, n32, 32
            ix = ix0
            do 20 i = i1, i2, inc
               ip = ipiv( ix )
               if( ip.ne.i ) then
                  do 10 k = j, j + 31
                     temp = a( i, k )
                     a( i, k ) = a( ip, k )
                     a( ip, k ) = temp
   10             continue
               end if
               ix = ix + incx
   20       continue
   30    continue
      end if
      if( n32.ne.n ) then
         n32 = n32 + 1
         ix = ix0
         do 50 i = i1, i2, inc
            ip = ipiv( ix )
            if( ip.ne.i ) then
               do 40 k = n32, n
                  temp = a( i, k )
                  a( i, k ) = a( ip, k )
                  a( ip, k ) = temp
   40          continue
            end if
            ix = ix + incx
   50    continue
      end if
*
      return
*
*     end of dlaswp
*
      end
ccc
      integer function ilaenv( ispec, name, opts, n1, n2, n3, n4 )
*
*  -- lapack auxiliary routine (version 3.2.1)                        --
*
*  -- april 2009                                                      --
*
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*
*     .. scalar arguments ..
      character*( * )    name, opts
      integer            ispec, n1, n2, n3, n4
*     ..
*
*  purpose
*  =======
*
*  ilaenv is called from the lapack routines to choose problem-dependent
*  parameters for the local environment.  see ispec for a description of
*  the parameters.
*
*  ilaenv returns an integer
*  if ilaenv >= 0: ilaenv returns the value of the parameter specified by ispec
*  if ilaenv < 0:  if ilaenv = -k, the k-th argument had an illegal value.
*
*  this version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  this routine will not function correctly if it is converted to all
*  lower case.  converting it to all upper case is allowed.
*
*  arguments
*  =========
*
*  ispec   (input) integer
*          specifies the parameter to be returned as the value of
*          ilaenv.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for n less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines (deprecated)
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ilaenv(2,...) and m by ilaenv(5,...)
*          = 6: the crossover point for the svd (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a qr factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift qr method
*               for nonsymmetric eigenvalue problems (deprecated)
*          = 9: maximum size of the subproblems at the bottom of the
*               computation tree in the divide-and-conquer algorithm
*               (used by xgelsd and xgesdd)
*          =10: ieee nan arithmetic can be trusted not to trap
*          =11: infinity arithmetic can be trusted not to trap
*          12 <= ispec <= 16:
*               xhseqr or one of its subroutines,
*               see iparmq for detailed explanation
*
*  name    (input) character*(*)
*          the name of the calling subroutine, in either upper case or
*          lower case.
*
*  opts    (input) character*(*)
*          the character options to the subroutine name, concatenated
*          into a single character string.  for example, uplo = 'u',
*          trans = 't', and diag = 'n' for a triangular routine would
*          be specified as opts = 'utn'.
*
*  n1      (input) integer
*  n2      (input) integer
*  n3      (input) integer
*  n4      (input) integer
*          problem dimensions for the subroutine name; these may not all
*          be required.
*
*  further details
*  ===============
*
*  the following conventions have been used when calling ilaenv from the
*  lapack routines:
*  1)  opts is a concatenation of all of the character options to
*      subroutine name, in the same order that they appear in the
*      argument list for name, even if they are not used in determining
*      the value of the parameter specified by ispec.
*  2)  the problem dimensions n1, n2, n3, n4 are specified in the order
*      that they appear in the argument list for name.  n1 is used
*      first, n2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  the parameter value returned by ilaenv is checked for validity in
*      the calling subroutine.  for example, ilaenv is used to retrieve
*      the optimal blocksize for strtri as follows:
*
*      nb = ilaenv( 1, 'strtri', uplo // diag, n, -1, -1, -1 )
*      if( nb.le.1 ) nb = max( 1, n )
*
*  =====================================================================
*
*     .. local scalars ..
      integer            i, ic, iz, nb, nbmin, nx
      logical            cname, sname
      character          c1*1, c2*2, c4*2, c3*3, subnam*6
*     ..
*     .. intrinsic functions ..
      intrinsic          char, ichar, int, min, real
*     ..
*     .. external functions ..
      integer            ieeeck, iparmq
      external           ieeeck, iparmq
*     ..
*     .. executable statements ..
*
      go to ( 10, 10, 10, 80, 90, 100, 110, 120,
     $        130, 140, 150, 160, 160, 160, 160, 160 )ispec
*
*     invalid value for ispec
*
      ilaenv = -1
      return
*
   10 continue
*
*     convert name to upper case if the first character is lower case.
*
      ilaenv = 1
      subnam = name
      ic = ichar( subnam( 1: 1 ) )
      iz = ichar( 'z' )
      if( iz.eq.90 .or. iz.eq.122 ) then
*
*        ascii character set
*
         if( ic.ge.97 .and. ic.le.122 ) then
            subnam( 1: 1 ) = char( ic-32 )
            do 20 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ic.ge.97 .and. ic.le.122 )
     $            subnam( i: i ) = char( ic-32 )
   20       continue
         end if
*
      else if( iz.eq.233 .or. iz.eq.169 ) then
*
*        ebcdic character set
*
         if( ( ic.ge.129 .and. ic.le.137 ) .or.
     $       ( ic.ge.145 .and. ic.le.153 ) .or.
     $       ( ic.ge.162 .and. ic.le.169 ) ) then
            subnam( 1: 1 ) = char( ic+64 )
            do 30 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ( ic.ge.129 .and. ic.le.137 ) .or.
     $             ( ic.ge.145 .and. ic.le.153 ) .or.
     $             ( ic.ge.162 .and. ic.le.169 ) )subnam( i:
     $             i ) = char( ic+64 )
   30       continue
         end if
*
      else if( iz.eq.218 .or. iz.eq.250 ) then
*
*        prime machines:  ascii+128
*
         if( ic.ge.225 .and. ic.le.250 ) then
            subnam( 1: 1 ) = char( ic-32 )
            do 40 i = 2, 6
               ic = ichar( subnam( i: i ) )
               if( ic.ge.225 .and. ic.le.250 )
     $            subnam( i: i ) = char( ic-32 )
   40       continue
         end if
      end if
*
      c1 = subnam( 1: 1 )
      sname = c1.eq.'s' .or. c1.eq.'d'
      cname = c1.eq.'c' .or. c1.eq.'z'
      if( .not.( cname .or. sname ) )
     $   return
      c2 = subnam( 2: 3 )
      c3 = subnam( 4: 6 )
      c4 = c3( 2: 3 )
*
      go to ( 50, 60, 70 )ispec
*
   50 continue
*
*     ispec = 1:  block size
*
*     in these examples, separate code is provided for setting nb for
*     real and complex.  we assume that nb will take the same value in
*     single or double precision.
*
      nb = 1
*
      if( c2.eq.'ge' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or.
     $            c3.eq.'qlf' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nb = 32
            else
               nb = 32
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'po' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nb = 32
         else if( sname .and. c3.eq.'gst' ) then
            nb = 64
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trf' ) then
            nb = 64
         else if( c3.eq.'trd' ) then
            nb = 32
         else if( c3.eq.'gst' ) then
            nb = 64
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nb = 32
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nb = 32
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nb = 32
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nb = 32
            end if
         end if
      else if( c2.eq.'gb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n4.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'pb' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            else
               if( n2.le.64 ) then
                  nb = 1
               else
                  nb = 32
               end if
            end if
         end if
      else if( c2.eq.'tr' ) then
         if( c3.eq.'tri' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( c2.eq.'la' ) then
         if( c3.eq.'uum' ) then
            if( sname ) then
               nb = 64
            else
               nb = 64
            end if
         end if
      else if( sname .and. c2.eq.'st' ) then
         if( c3.eq.'ebz' ) then
            nb = 1
         end if
      end if
      ilaenv = nb
      return
*
   60 continue
*
*     ispec = 2:  minimum block size
*
      nbmin = 2
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or. c3.eq.
     $       'qlf' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         else if( c3.eq.'tri' ) then
            if( sname ) then
               nbmin = 2
            else
               nbmin = 2
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( c3.eq.'trf' ) then
            if( sname ) then
               nbmin = 8
            else
               nbmin = 8
            end if
         else if( sname .and. c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nbmin = 2
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nbmin = 2
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nbmin = 2
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nbmin = 2
            end if
         else if( c3( 1: 1 ).eq.'m' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nbmin = 2
            end if
         end if
      end if
      ilaenv = nbmin
      return
*
   70 continue
*
*     ispec = 3:  crossover point
*
      nx = 0
      if( c2.eq.'ge' ) then
         if( c3.eq.'qrf' .or. c3.eq.'rqf' .or. c3.eq.'lqf' .or. c3.eq.
     $       'qlf' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'hrd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         else if( c3.eq.'brd' ) then
            if( sname ) then
               nx = 128
            else
               nx = 128
            end if
         end if
      else if( c2.eq.'sy' ) then
         if( sname .and. c3.eq.'trd' ) then
            nx = 32
         end if
      else if( cname .and. c2.eq.'he' ) then
         if( c3.eq.'trd' ) then
            nx = 32
         end if
      else if( sname .and. c2.eq.'or' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nx = 128
            end if
         end if
      else if( cname .and. c2.eq.'un' ) then
         if( c3( 1: 1 ).eq.'g' ) then
            if( c4.eq.'qr' .or. c4.eq.'rq' .or. c4.eq.'lq' .or. c4.eq.
     $          'ql' .or. c4.eq.'hr' .or. c4.eq.'tr' .or. c4.eq.'br' )
     $           then
               nx = 128
            end if
         end if
      end if
      ilaenv = nx
      return
*
   80 continue
*
*     ispec = 4:  number of shifts (used by xhseqr)
*
      ilaenv = 6
      return
*
   90 continue
*
*     ispec = 5:  minimum column dimension (not used)
*
      ilaenv = 2
      return
*
  100 continue
*
*     ispec = 6:  crossover point for svd (used by xgelss and xgesvd)
*
      ilaenv = int( real( min( n1, n2 ) )*1.6e0 )
      return
*
  110 continue
*
*     ispec = 7:  number of processors (not used)
*
      ilaenv = 1
      return
*
  120 continue
*
*     ispec = 8:  crossover point for multishift (used by xhseqr)
*
      ilaenv = 50
      return
*
  130 continue
*
*     ispec = 9:  maximum size of the subproblems at the bottom of the
*                 computation tree in the divide-and-conquer algorithm
*                 (used by xgelsd and xgesdd)
*
      ilaenv = 25
      return
*
  140 continue
*
*     ispec = 10: ieee nan arithmetic can be trusted not to trap
*
*     ilaenv = 0
      ilaenv = 1
      if( ilaenv.eq.1 ) then
         ilaenv = ieeeck( 1, 0.0, 1.0 )
      end if
      return
*
  150 continue
*
*     ispec = 11: infinity arithmetic can be trusted not to trap
*
*     ilaenv = 0
      ilaenv = 1
      if( ilaenv.eq.1 ) then
         ilaenv = ieeeck( 0, 0.0, 1.0 )
      end if
      return
*
  160 continue
*
*     12 <= ispec <= 16: xhseqr or one of its subroutines. 
*
      ilaenv = iparmq( ispec, name, opts, n1, n2, n3, n4 )
      return
*
*     end of ilaenv
*
      end
ccc
      integer          function ieeeck( ispec, zero, one )
*
*  -- lapack auxiliary routine (version 3.2.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     june 2010
*
*     .. scalar arguments ..
      integer            ispec
      real               one, zero
*     ..
*
*  purpose
*  =======
*
*  ieeeck is called from the ilaenv to verify that infinity and
*  possibly nan arithmetic is safe (i.e. will not trap).
*
*  arguments
*  =========
*
*  ispec   (input) integer
*          specifies whether to test just for inifinity arithmetic
*          or whether to test for infinity and nan arithmetic.
*          = 0: verify infinity arithmetic only.
*          = 1: verify infinity and nan arithmetic.
*
*  zero    (input) real
*          must contain the value 0.0
*          this is passed to prevent the compiler from optimizing
*          away this code.
*
*  one     (input) real
*          must contain the value 1.0
*          this is passed to prevent the compiler from optimizing
*          away this code.
*
*  return value:  integer
*          = 0:  arithmetic failed to produce the correct answers
*          = 1:  arithmetic produced the correct answers
*
*     .. local scalars ..
      real               nan1, nan2, nan3, nan4, nan5, nan6, neginf,
     $                   negzro, newzro, posinf
*     ..
*     .. executable statements ..
      ieeeck = 1
*
      posinf = one / zero
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if
*
      neginf = -one / zero
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if
*
      negzro = one / ( neginf+one )
      if( negzro.ne.zero ) then
         ieeeck = 0
         return
      end if
*
      neginf = one / negzro
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if
*
      newzro = negzro + zero
      if( newzro.ne.zero ) then
         ieeeck = 0
         return
      end if
*
      posinf = one / newzro
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if
*
      neginf = neginf*posinf
      if( neginf.ge.zero ) then
         ieeeck = 0
         return
      end if
*
      posinf = posinf*posinf
      if( posinf.le.one ) then
         ieeeck = 0
         return
      end if
*
*
*
*
*     return if we were only asked to check infinity arithmetic
*
      if( ispec.eq.0 )
     $   return
*
      nan1 = posinf + neginf
*
      nan2 = posinf / neginf
*
      nan3 = posinf / posinf
*
      nan4 = posinf*zero
*
      nan5 = neginf*negzro
*
      nan6 = nan5*zero
*
      if( nan1.eq.nan1 ) then
         ieeeck = 0
         return
      end if
*
      if( nan2.eq.nan2 ) then
         ieeeck = 0
         return
      end if
*
      if( nan3.eq.nan3 ) then
         ieeeck = 0
         return
      end if
*
      if( nan4.eq.nan4 ) then
         ieeeck = 0
         return
      end if
*
      if( nan5.eq.nan5 ) then
         ieeeck = 0
         return
      end if
*
      if( nan6.eq.nan6 ) then
         ieeeck = 0
         return
      end if
*
      return
      end
ccc
      integer function idamax(n,dx,incx)
*     .. scalar arguments ..
      integer incx,n
*     ..
*     .. array arguments ..
      double precision dx(*)
*     ..
*
*  purpose
*  =======
*
*     idamax finds the index of element having max. absolute value.
*
*  further details
*  ===============
*
*     jack dongarra, linpack, 3/11/78.
*     modified 3/93 to return if incx .le. 0.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*  =====================================================================
*
*     .. local scalars ..
      double precision dmax
      integer i,ix
*     ..
*     .. intrinsic functions ..
      intrinsic dabs
*     ..
      idamax = 0
      if (n.lt.1 .or. incx.le.0) return
      idamax = 1
      if (n.eq.1) return
      if (incx.eq.1) go to 20
*
*        code for increment not equal to 1
*
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
          if (dabs(dx(ix)).le.dmax) go to 5
          idamax = i
          dmax = dabs(dx(ix))
    5     ix = ix + incx
   10 continue
      return
*
*        code for increment equal to 1
*
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
          if (dabs(dx(i)).le.dmax) go to 30
          idamax = i
          dmax = dabs(dx(i))
   30 continue
      return
      end
ccc
      logical function lsame(ca,cb)
*
*  -- lapack auxiliary routine (version 3.1) --
*     univ. of tennessee, univ. of california berkeley and nag ltd..
*     november 2006
*
*     .. scalar arguments ..
      character ca,cb
*     ..
*
*  purpose
*  =======
*
*  lsame returns .true. if ca is the same letter as cb regardless of
*  case.
*
*  arguments
*  =========
*
*  ca      (input) character*1
*
*  cb      (input) character*1
*          ca and cb specify the single characters to be compared.
*
* =====================================================================
*
*     .. intrinsic functions ..
      intrinsic ichar
*     ..
*     .. local scalars ..
      integer inta,intb,zcode
*     ..
*
*     test if the characters are equal
*
      lsame = ca .eq. cb
      if (lsame) return
*
*     now test for equivalence if both characters are alphabetic.
*
      zcode = ichar('z')
*
*     use 'z' rather than 'a' so that ascii can be detected on prime
*     machines, on which ichar returns a value with bit 8 set.
*     ichar('a') on prime machines returns 193 which is the same as
*     ichar('a') on an ebcdic machine.
*
      inta = ichar(ca)
      intb = ichar(cb)
*
      if (zcode.eq.90 .or. zcode.eq.122) then
*
*        ascii is assumed - zcode is the ascii code of either lower or
*        upper case 'z'.
*
          if (inta.ge.97 .and. inta.le.122) inta = inta - 32
          if (intb.ge.97 .and. intb.le.122) intb = intb - 32
*
      else if (zcode.eq.233 .or. zcode.eq.169) then
*
*        ebcdic is assumed - zcode is the ebcdic code of either lower or
*        upper case 'z'.
*
          if (inta.ge.129 .and. inta.le.137 .or.
     +        inta.ge.145 .and. inta.le.153 .or.
     +        inta.ge.162 .and. inta.le.169) inta = inta + 64
          if (intb.ge.129 .and. intb.le.137 .or.
     +        intb.ge.145 .and. intb.le.153 .or.
     +        intb.ge.162 .and. intb.le.169) intb = intb + 64
*
      else if (zcode.eq.218 .or. zcode.eq.250) then
*
*        ascii is assumed, on prime machines - zcode is the ascii code
*        plus 128 of either lower or upper case 'z'.
*
          if (inta.ge.225 .and. inta.le.250) inta = inta - 32
          if (intb.ge.225 .and. intb.le.250) intb = intb - 32
      end if
      lsame = inta .eq. intb
*
*     return
*
*     end of lsame
*
      end
ccc
      integer function iparmq( ispec, name, opts, n, ilo, ihi, lwork )
*
*  -- lapack auxiliary routine (version 3.2) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     november 2006
*     
*     .. scalar arguments ..
      integer            ihi, ilo, ispec, lwork, n
      character          name*( * ), opts*( * )
*
*  purpose
*  =======
*
*       this program sets problem and machine dependent parameters
*       useful for xhseqr and its subroutines. it is called whenever 
*       ilaenv is called with 12 <= ispec <= 16
*
*  arguments
*  =========
*
*       ispec  (input) integer scalar
*              ispec specifies which tunable parameter iparmq should
*              return.
*
*              ispec=12: (inmin)  matrices of order nmin or less
*                        are sent directly to xlahqr, the implicit
*                        double shift qr algorithm.  nmin must be
*                        at least 11.
*
*              ispec=13: (inwin)  size of the deflation window.
*                        this is best set greater than or equal to
*                        the number of simultaneous shifts ns.
*                        larger matrices benefit from larger deflation
*                        windows.
*
*              ispec=14: (inibl) determines when to stop nibbling and
*                        invest in an (expensive) multi-shift qr sweep.
*                        if the aggressive early deflation subroutine
*                        finds ld converged eigenvalues from an order
*                        nw deflation window and ld.gt.(nw*nibble)/100,
*                        then the next qr sweep is skipped and early
*                        deflation is applied immediately to the
*                        remaining active diagonal block.  setting
*                        iparmq(ispec=14) = 0 causes ttqre to skip a
*                        multi-shift qr sweep whenever early deflation
*                        finds a converged eigenvalue.  setting
*                        iparmq(ispec=14) greater than or equal to 100
*                        prevents ttqre from skipping a multi-shift
*                        qr sweep.
*
*              ispec=15: (nshfts) the number of simultaneous shifts in
*                        a multi-shift qr iteration.
*
*              ispec=16: (iacc22) iparmq is set to 0, 1 or 2 with the
*                        following meanings.
*                        0:  during the multi-shift qr sweep,
*                            xlaqr5 does not accumulate reflections and
*                            does not use matrix-matrix multiply to
*                            update the far-from-diagonal matrix
*                            entries.
*                        1:  during the multi-shift qr sweep,
*                            xlaqr5 and/or xlaqraccumulates reflections and uses
*                            matrix-matrix multiply to update the
*                            far-from-diagonal matrix entries.
*                        2:  during the multi-shift qr sweep.
*                            xlaqr5 accumulates reflections and takes
*                            advantage of 2-by-2 block structure during
*                            matrix-matrix multiplies.
*                        (if xtrmm is slower than xgemm, then
*                        iparmq(ispec=16)=1 may be more efficient than
*                        iparmq(ispec=16)=2 despite the greater level of
*                        arithmetic work implied by the latter choice.)
*
*       name    (input) character string
*               name of the calling subroutine
*
*       opts    (input) character string
*               this is a concatenation of the string arguments to
*               ttqre.
*
*       n       (input) integer scalar
*               n is the order of the hessenberg matrix h.
*
*       ilo     (input) integer
*       ihi     (input) integer
*               it is assumed that h is already upper triangular
*               in rows and columns 1:ilo-1 and ihi+1:n.
*
*       lwork   (input) integer scalar
*               the amount of workspace available.
*
*  further details
*  ===============
*
*       little is known about how best to choose these parameters.
*       it is possible to use different values of the parameters
*       for each of chseqr, dhseqr, shseqr and zhseqr.
*
*       it is probably best to choose different parameters for
*       different matrices and different parameters at different
*       times during the iteration, but this has not been
*       implemented --- yet.
*
*
*       the best choices of most of the parameters depend
*       in an ill-understood way on the relative execution
*       rate of xlaqr3 and xlaqr5 and on the nature of each
*       particular eigenvalue problem.  experiment may be the
*       only practical way to determine which choices are most
*       effective.
*
*       following is a list of default values supplied by iparmq.
*       these defaults may be adjusted in order to attain better
*       performance in any particular computational environment.
*
*       iparmq(ispec=12) the xlahqr vs xlaqr0 crossover point.
*                        default: 75. (must be at least 11.)
*
*       iparmq(ispec=13) recommended deflation window size.
*                        this depends on ilo, ihi and ns, the
*                        number of simultaneous shifts returned
*                        by iparmq(ispec=15).  the default for
*                        (ihi-ilo+1).le.500 is ns.  the default
*                        for (ihi-ilo+1).gt.500 is 3*ns/2.
*
*       iparmq(ispec=14) nibble crossover point.  default: 14.
*
*       iparmq(ispec=15) number of simultaneous shifts, ns.
*                        a multi-shift qr iteration.
*
*                        if ihi-ilo+1 is ...
*
*                        greater than      ...but less    ... the
*                        or equal to ...      than        default is
*
*                                0               30       ns =   2+
*                               30               60       ns =   4+
*                               60              150       ns =  10
*                              150              590       ns =  **
*                              590             3000       ns =  64
*                             3000             6000       ns = 128
*                             6000             infinity   ns = 256
*
*                    (+)  by default matrices of this order are
*                         passed to the implicit double shift routine
*                         xlahqr.  see iparmq(ispec=12) above.   these
*                         values of ns are used only in case of a rare
*                         xlahqr failure.
*
*                    (**) the asterisks (**) indicate an ad-hoc
*                         function increasing from 10 to 64.
*
*       iparmq(ispec=16) select structured matrix multiply.
*                        (see ispec=16 above for details.)
*                        default: 3.
*
*     ================================================================
*     .. parameters ..
      integer            inmin, inwin, inibl, ishfts, iacc22
      parameter          ( inmin = 12, inwin = 13, inibl = 14,
     $                   ishfts = 15, iacc22 = 16 )
      integer            nmin, k22min, kacmin, nibble, knwswp
      parameter          ( nmin = 75, k22min = 14, kacmin = 14,
     $                   nibble = 14, knwswp = 500 )
      real               two
      parameter          ( two = 2.0 )
*     ..
*     .. local scalars ..
      integer            nh, ns
*     ..
*     .. intrinsic functions ..
      intrinsic          log, max, mod, nint, real
*     ..
*     .. executable statements ..
      if( ( ispec.eq.ishfts ) .or. ( ispec.eq.inwin ) .or.
     $    ( ispec.eq.iacc22 ) ) then
*
*        ==== set the number simultaneous shifts ====
*
         nh = ihi - ilo + 1
         ns = 2
         if( nh.ge.30 )
     $      ns = 4
         if( nh.ge.60 )
     $      ns = 10
         if( nh.ge.150 )
     $      ns = max( 10, nh / nint( log( real( nh ) ) / log( two ) ) )
         if( nh.ge.590 )
     $      ns = 64
         if( nh.ge.3000 )
     $      ns = 128
         if( nh.ge.6000 )
     $      ns = 256
         ns = max( 2, ns-mod( ns, 2 ) )
      end if
*
      if( ispec.eq.inmin ) then
*
*
*        ===== matrices of order smaller than nmin get sent
*        .     to xlahqr, the classic double shift algorithm.
*        .     this must be at least 11. ====
*
         iparmq = nmin
*
      else if( ispec.eq.inibl ) then
*
*        ==== inibl: skip a multi-shift qr iteration and
*        .    whenever aggressive early deflation finds
*        .    at least (nibble*(window size)/100) deflations. ====
*
         iparmq = nibble
*
      else if( ispec.eq.ishfts ) then
*
*        ==== nshfts: the number of simultaneous shifts =====
*
         iparmq = ns
*
      else if( ispec.eq.inwin ) then
*
*        ==== nw: deflation window size.  ====
*
         if( nh.le.knwswp ) then
            iparmq = ns
         else
            iparmq = 3*ns / 2
         end if
*
      else if( ispec.eq.iacc22 ) then
*
*        ==== iacc22: whether to accumulate reflections
*        .     before updating the far-from-diagonal elements
*        .     and whether to use 2-by-2 block structure while
*        .     doing it.  a small amount of work could be saved
*        .     by making this choice dependent also upon the
*        .     nh=ihi-ilo+1.
*
         iparmq = 0
         if( ns.ge.kacmin )
     $      iparmq = 1
         if( ns.ge.k22min )
     $      iparmq = 2
*
      else
*        ===== invalid value of ispec =====
         iparmq = -1
*
      end if
*
*     ==== end of iparmq ====
*
      end
ccc
      double precision function dlamch( cmach )
*
*  -- lapack auxiliary routine (version 3.3.0) --
*  -- lapack is a software package provided by univ. of tennessee,    --
*  -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
*     based on lapack dlamch but with fortran 95 query functions
*     see: http://www.cs.utk.edu/~luszczek/lapack/lamch.html
*     and  http://www.netlib.org/lapack-dev/lapack-coding/program-style.html#id2537289
*     july 2010
*
*     .. scalar arguments ..
      character          cmach
*     ..
*
*  purpose
*  =======
*
*  dlamch determines double precision machine parameters.
*
*  arguments
*  =========
*
*  cmach   (input) character*1
*          specifies the value to be returned by dlamch:
*          = 'e' or 'e',   dlamch := eps
*          = 's' or 's ,   dlamch := sfmin
*          = 'b' or 'b',   dlamch := base
*          = 'p' or 'p',   dlamch := eps*base
*          = 'n' or 'n',   dlamch := t
*          = 'r' or 'r',   dlamch := rnd
*          = 'm' or 'm',   dlamch := emin
*          = 'u' or 'u',   dlamch := rmin
*          = 'l' or 'l',   dlamch := emax
*          = 'o' or 'o',   dlamch := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*
* =====================================================================
*
*     .. parameters ..
      double precision   one, zero
      parameter          ( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. local scalars ..
      double precision   rnd, eps, sfmin, small, rmach
*     ..
*     .. external functions ..
      logical            lsame
      external           lsame
*     ..
*     .. intrinsic functions ..
      intrinsic          digits, epsilon, huge, maxexponent,
     $                   minexponent, radix, tiny
*     ..
*     .. executable statements ..
*
*
*     assume rounding, not chopping. always.
*
      rnd = one
*
      if( one.eq.rnd ) then
         eps = epsilon(zero) * 0.5
      else
         eps = epsilon(zero)
      end if
*
      if( lsame( cmach, 'e' ) ) then
         rmach = eps
      else if( lsame( cmach, 's' ) ) then
         sfmin = tiny(zero)
         small = one / huge(zero)
         if( small.ge.sfmin ) then
*
*           use small plus a bit, to avoid the possibility of rounding
*           causing overflow when computing  1/sfmin.
*
            sfmin = small*( one+eps )
         end if
         rmach = sfmin
      else if( lsame( cmach, 'b' ) ) then
         rmach = radix(zero)
      else if( lsame( cmach, 'p' ) ) then
         rmach = eps * radix(zero)
      else if( lsame( cmach, 'n' ) ) then
         rmach = digits(zero)
      else if( lsame( cmach, 'r' ) ) then
         rmach = rnd
      else if( lsame( cmach, 'm' ) ) then
         rmach = minexponent(zero)
      else if( lsame( cmach, 'u' ) ) then
         rmach = tiny(zero)
      else if( lsame( cmach, 'l' ) ) then
         rmach = maxexponent(zero)
      else if( lsame( cmach, 'o' ) ) then
         rmach = huge(zero)
      else
         rmach = zero
      end if
*
      dlamch = rmach
      return
*
*     end of dlamch
*
      end
************************************************************************
*
      double precision function dlamc3( a, b )
*
*  -- lapack auxiliary routine (version 3.3.0) --
*     univ. of tennessee, univ. of california berkeley and nag ltd..
*     november 2010
*
*     .. scalar arguments ..
      double precision   a, b
*     ..
*
*  purpose
*  =======
*
*  dlamc3  is intended to force  a  and  b  to be stored prior to doing
*  the addition of  a  and  b ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  arguments
*  =========
*
*  a       (input) double precision
*  b       (input) double precision
*          the values a and b.
*
* =====================================================================
*
*     .. executable statements ..
*
      dlamc3 = a + b
*
      return
*
*     end of dlamc3
*
      end
*
************************************************************************
c $Id: rism1d.F 21176 2011-10-10 06:35:49Z d3y133 $
