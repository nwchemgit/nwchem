C> \ingroup task
C> @{
       function task_rism(rtdb)
       implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
       integer rtdb
       logical  task_rism
       integer n
       logical rism_util_power_2
       external rism_util_power_2
c
       call rism_print_header()
c       
c      load data and allocate memory
c      -----------------------------
       call rism_prepare(rtdb)
#ifdef RISM_DEBUG       
       call rism_message("finished with rism_prepare")
#endif       
       call rism_print_params(rtdb)
       call rism_print_solute_configuration(rtdb)
c       stop
c       
c      main rism routine call
c      ----------------------
       call rism_message("calling rism wrapper")
       call rism_wrapper(rtdb)
       call rism_cleanup() 
       task_rism = .true.
       call rism_message("completed task rism")
       return
      end 
C> @}
c       
       subroutine rism_wrapper(rtdb)
       implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
      integer nu,nv,nvv,ngr
      integer icl,icr
      integer i_rgrid,i_kgrid
      integer i_tv,i_den,i_isv,i_mv
      integer i_xv,i_yv,i_zv
      integer i_ims
      logical result
      character*32 sname,dname,pname
      double precision t,tau,lambd,tol,solperm
      integer dd
      integer i_sgvv,i_epsvv,i_qvv
      integer i_sigu,i_epsiu,i_qqu,i_wu
      character*72 rdffile,tag
      logical okspace
c       
      pname = "rism_wrapper"
c      
c
c     solute data
c     -------------- 
      sname = "solute"
      dname = "natoms"
      call db_data_get_int(sname,dname,1,nu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get"//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "sigma"
      call db_data_get_index(sname,dname,i_sigu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "epsilon"
      call db_data_get_index(sname,dname,i_epsiu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "charge"
      call db_data_get_index(sname,dname,i_qqu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "struct_factor"
      call db_data_get_index(sname,dname,i_wu,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)


#ifdef RISM_DEBUG      
      call rism_message("finished with solute data")
#endif      
c
c     solvent data
c     -------------- 
      sname = "solvent"
      dname = "natoms"
      call db_data_get_int(sname,dname,1,nv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get"//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "natoms_reduced"
      call db_data_get_int(sname,dname,1,nvv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get"//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "atom_name"
      call db_data_get_index(sname,dname,i_tv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "residue_index"
      call db_data_get_index(sname,dname,i_isv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "density"
      call db_data_get_index(sname,dname,i_den,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "multiplicity"
      call db_data_get_index(sname,dname,i_mv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "xcoord"
      call db_data_get_index(sname,dname,i_xv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "ycoord"
      call db_data_get_index(sname,dname,i_yv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "zcoord"
      call db_data_get_index(sname,dname,i_zv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "map_reduced"
      call db_data_get_index(sname,dname,i_ims,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)


      dname = "sigma_reduced"
      call db_data_get_index(sname,dname,i_sgvv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "epsilon_reduced"
      call db_data_get_index(sname,dname,i_epsvv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "charge_reduced"
      call db_data_get_index(sname,dname,i_qvv,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      call rism_message("rism_wrapper 1")
c
c     grid data
c     -------------- 
      sname = "grid"
      dname = "npoints"
      call db_data_get_int(sname,dname,1,ngr,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "rgrid"
      call db_data_get_index(sname,dname,i_rgrid,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "okspace"
      call db_data_get_log(sname,dname,1,okspace,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      dname = "kgrid"
      call db_data_get_index(sname,dname,i_kgrid,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
      call rism_message("rism_wrapper 2")
c
c     parameters
c     -------------- 
      sname = "parameters"
      dname = "closure"
      call db_data_get_int(sname,dname,1,icl,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "vdw_rule"
      call db_data_get_int(sname,dname,1,icr,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "solvent_permittivity"
      call db_data_get_dbl(sname,dname,1,solperm,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "tau"
      call db_data_get_dbl(sname,dname,1,tau,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)


      dname = "tolerance"
      call db_data_get_dbl(sname,dname,1,tol,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "mixing"
      call db_data_get_dbl(sname,dname,1,lambd,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "temperature"
      call db_data_get_dbl(sname,dname,1,t,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)

      dname = "diis"
      call db_data_get_int(sname,dname,1,dd,result)
      if(.not.result)
     >  call errquit(pname//"cannot get "//sname(1:inp_strlen(sname))
     >               //dname(1:inp_strlen(dname)),0,0)
c      
c     get filename for solvent g(r)
c     -----------------------------
      tag = "rism:solvent:rdf"
      if(.not.rtdb_cget(rtdb,tag,1,rdffile)) 
     >  call errquit("cannot get "//tag,0,0)
      call rism_message("rism_wrapper 3")

       call rism_message("getting ready for main rism")

       call rism(rtdb,rdffile,nu,nv,nvv,ngr,icl,
     +           dbl_mb(i_kgrid),dbl_mb(i_rgrid),
     +           byte_mb(i_tv),int_mb(i_isv),dbl_mb(i_den),
     +           int_mb(i_mv),
     +           dbl_mb(i_xv),dbl_mb(i_yv),dbl_mb(i_zv),
     +           icr,int_mb(i_ims),tau,solperm,
     +           dbl_mb(i_sgvv),dbl_mb(i_epsvv),
     +           dbl_mb(i_qvv),dbl_mb(i_sigu),
     +           dbl_mb(i_epsiu),dbl_mb(i_qqu),
     +           t,lambd,
     +           tol,dbl_mb(i_wu),dd,okspace)          

c       call rism(nu,nv,nvv,ngr,icl,kgrid,rgrid,tv,isv,den,mv,xv,yv,zv,
c     * icr,ims,tau,solperm,sgvv,epsvv,qvv,sigu,epsiu,qqu,t,lambd,
c     * tol,wu,dd)          

       return
      end 

c creat susceptibility of solvent
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
       external util_erf
       double precision util_erf
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
     *    *qqu(j1)*qvv(j2)/rgrid(i)*(1-util_erf(tau*rgrid(i)))
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
       subroutine rism(rtdb,rdffile,nu,nv,nvv,ngr,
     * icl,kgrid,rgrid,tv,isv,
     * den,mv,xv,yv,zv,
     * icr,ims,tau,solperm,sgvv,epsvv,qvv,sigu,epsiu,qqu,t,
     * lambd,tol,wu,dd,okspace)       
       implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
c      parameters
       integer rtdb
       integer kd,dd 
       character*(*) rdffile
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
       logical okspace
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
       double precision muh,mugf
       write(76,*) nu,nv,nvv,ngr,icl,icr,tau,solperm
       write(76,*) T,lambd,tol

c
       pi=2*asin(1.0)       
c      nd total number of different rdfs
       nd=(nvv+1)*nvv/2     
c      chi(i,j,k) solvent susceptibility
c
c      WORK HERE
       call chicreat(rdffile,nd,nv,ngr,ims,nvv,
     +               kgrid,tv,isv,den,mv,xv,yv,zv,chi,
     +               qvv,tau,okspace) 
c       
c       call chicreat(rdffile,nd,nv,ngr,ims,nvv,
c     +               kgrid,tv,isv,den,mv,xv,yv,zv,chi)

c
c      plj(j1,j2,r) lj potential+short part from coulomb
c      ul(j1,j2,k) is the fourier trasnform of long range interactions
       call potcreat(nvv,nu,ngr,icr,kgrid,rgrid,tau,solperm,
     * sgvv,epsvv,qvv,sigu,epsiu,qqu,plj,ul) 
c
cc     initial value of c_short(r)       
       kb=8.3144e-3
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
       write(*,*) "starting iterations"
1      continue 
c      calculating c_new(r) and h_new(r) from closure
       call clos(nvv,nu,ngr,t,plj,cnew,gold,icl)
       do j=1,nu
        do j1=1,nvv
         do i=1,ngr
          hold(j,j1,i)=gold(j,j1,i)+cnew(j,j1,i) 
         enddo
        enddo
       enddo 
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
       dk=kgrid(2)-kgrid(1)
       pi=2*asin(1.0)       
c       adding the long-range potential to c_short 
       do i=1,ngr
        do j=1,nu
         do j1=1,nvv
          cf3(j,j1,i)=cf3(j,j1,i)-1.0/kb/t*ul(j,j1,i)
         enddo
        enddo
       enddo 
c      calculations of fourier transform of h(i,j,k) by
c      site-site ornstein-zerinike equations       
       call ssoz(nvv,nu,ngr,nv,ims,mv,wu,cf3,chi,hu)
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
       k1=k1+1 
       write(*,*) k1,kd,del0       
       if(k1.ge.500) goto 2         
        if(del0.ge.tol) goto 1
 2     continue  
       call thermo(nu,nv,nvv,ngr,icl,ims,den,t,rgrid,qqu,qvv,
     * hnew,cnew,solperm,tau,muh,mugf)
       if (.not. rtdb_put(rtdb, 'rism:energy', mt_dbl, 1, mugf))
     $           call errquit('failed rism',0,0)
       if (.not. rtdb_put(rtdb, 'rism:muh', mt_dbl, 1, muh))
     $           call errquit('failed rism',0,0)
       if (.not. rtdb_put(rtdb, 'rism:mugf', mt_dbl, 1, mugf))
     $           call errquit('failed rism',0,0)


       return 
       end subroutine
cc
c 
c      subroutine for evaluation diis residues dgg and matrix coefficients ss
c
       subroutine diis(nu,nvv,ngr,lam,k1,dd,hold,hnew,ggo,dgg)
       implicit none
       integer k1,dd,jd,jd1
       real*8 ggo(1:dd,1:nu,1:nvv,1:ngr),dgg(1:dd,1:nu,1:nvv,1:ngr)
       real*8 ss(1:dd+1,1:dd+1),hold(1:nu,1:nvv,1:ngr)
       real*8 s0(1:dd+1),s1,lam,hnew(1:nu,1:nvv,1:ngr)
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
c
       subroutine thermo(nu,nv,nvv,ngr,icl,ims,den,t,rgrid,qqu,qvv,
     * hnew,cnew,solperm,tau,muh,mugf) 
       implicit none
#include "stdio.fh"
       real*8 muh, mugf, muuh(1:nu,1:nvv),muugf(1:nu,1:nvv)        
       integer  nu, nv,nvv, ngr, i,j1,j,icl,nd        
       real*8 rgrid(1:ngr),dentv(1:nvv),cu(1:nu*nvv,1:ngr) 
       real*8 pi,dr,kb,t,kcal,theta,ck,echar,solperm,nav,tau
       real*8 den(1:nv),gu(1:nu*nvv,1:ngr),nc(1:nu*nvv,1:ngr)
       integer ims(1:nv)
       real*8 cnew(1:nu,1:nvv,1:ngr), hnew(1:nu,1:nvv,1:ngr)
       real*8 ul(1:nu,1:nvv,1:ngr),hq(1:nu,1:nvv,1:ngr)
       real*8 ulr(1:nu,1:nvv,1:ngr),qvv(1:nvv),qqu(1:nu)
       external util_erf
       double precision util_erf
       do j=1,nv
        dentv(ims(j))=den(ims(j))
       enddo      
       pi=2*asin(1.0)
       kb=8.31441e-3
       kcal=0.0019859
       nd=(nvv+1)*nvv/2 
       dr=rgrid(2)-rgrid(1) 
       do j=1,nu
        do j1=1,nvv
         gu(nvv*(j-1)+j1,1)=0
         nc(nvv*(j-1)+j1,1)=0
         cu(nvv*(j-1)+j1,1)=cnew(j,j1,1) 
         do i=2,ngr 
          cu(nvv*(j-1)+j1,i)=cnew(j,j1,i) 
          gu(nvv*(j-1)+j1,i)=hnew(j,j1,i-1)+1 
          nc(nvv*(j-1)+j1,i)=nc(nvv*(j-1)+j1,i-1)
     *    +4*pi*dr*dentv(j1)*gu(nvv*(j-1)+j1,i)*rgrid(i)**2.0
         enddo
        enddo
       enddo 
       open(3, file='rdf_out.data', status='unknown')  
        do i=1,ngr        
         write(3,*) rgrid(i),(gu(j,i), j=1,nu*nvv)          
        enddo
       close(3)
       open(3, file='nc_out.data', status='unknown')  
        do i=1,ngr        
         write(3,*) rgrid(i),(nc(j,i), j=1,nu*nvv)          
        enddo
       close(3)
       open(3, file='c_out.data', status='unknown')  
        do i=1,ngr        
         write(3,*) rgrid(i),(cu(j,i), j=1,nu*nvv)          
        enddo
       close(3)
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
c      ulr electrostatic part of the interaction potential in real space
       do j=1,nu
        do j1=1,nvv      
         do i=2,ngr 
          ulr(j,j1,i)=ck*qqu(j)*qvv(j1)/rgrid(i)*util_erf(tau*rgrid(i))
         enddo 
         ulr(j,j1,1)=ulr(j,j1,2)
        enddo
       enddo 
c      calculations of partial contibutions 
       do j=1,nu
        do j1=1,nvv  
         do i=2,ngr
          if(icl.eq.1) then 
            hq(j,j1,i)=hnew(j,j1,i)
          endif
          if(icl.eq.2) then 
           hq(j,j1,i)= -theta(-hnew(j,j1,i))
          endif          
          muugf(j,j1)=muugf(j,j1)-(2*cnew(j,j1,i)+hnew(j,j1,i)*
     *    (cnew(j,j1,i)-ulr(j,j1,i)/kb/t))*rgrid(i)**2*dr*2*pi
          muuh(j,j1)=muuh(j,j1)-(2*cnew(j,j1,i)+hnew(j,j1,i)*
     *    (cnew(j,j1,i)-ulr(j,j1,i)/kb/t))*rgrid(i)**2*dr*2*pi
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
      
c      mugf chem.potential in gaussian approximation
c      muh   chem.potential in hnc approximation         
       write(luout,108) "(hnc approximation)",muh
       write(luout,108) "(gaussian approximation)",mugf
  108  format("Chemical potential",A,T43,F10.4)
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

c $Id$
