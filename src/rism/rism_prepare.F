      subroutine rism_prepare(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer rtdb
      logical master

      master = ga_nodeid().eq.0
c
c     order may be important here
c     ---------------------------
      if(master) then
       call rism_prepare_params(rtdb)
       call rism_prepare_grid(rtdb)
       call rism_prepare_solvent(rtdb)
       call rism_prepare_solute(rtdb)
c       call db_data_print(6)
      end if
c
c     send data to other processors
c     -----------------------------
      call ga_sync()
      call smd_data_brdcst_create(master)
      return
 
      end

      subroutine rism_prepare_params(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
#include "rism_params.fh"
      integer rtdb
c
      character*(72) filename
      integer i
      integer fn
      character*72 tag
      character*72 buffer
      character*1 b1
      integer master
      integer nu
      character*32 sname
      double precision dr,r1,pi,dk
      integer i_kgrid,i_rgrid,ngrid
      integer icl
      integer icr
      double precision f
      integer ipar
c 
      sname = "parameters"
      call db_namespace_create(sname)
c      
c     closure type
c     -----------------------------
      tag = "rism:closure"
      if(.not.rtdb_cget(rtdb,tag,1,buffer)) 
     >  tag = "hnc"
      if (inp_compare(.false.,"hnc", buffer)) then
        icl=t_hnc
      else if (inp_compare(.false.,"kh", buffer)) then
        icl=t_kh
      else
        call errquit("unknown closure "//buffer,0,0)
      end if

      tag = "closure"
      call db_data_put_int(sname,tag,1,icl)
c      
c     vdw combination rule
c     -----------------------------
      tag = "rism:vdw:rule"
      if(.not.rtdb_cget(rtdb,tag,1,buffer)) 
     >  tag="geometric"
      if (inp_compare(.false.,"arithmetic", buffer)) then
        icr=t_arithmetic
      else if (inp_compare(.false.,"geometric", buffer)) then
        icr=t_geometric
      else
        call errquit("unknown rule "//buffer,0,0)
      end if

      tag = "vdw_rule"
      call db_data_put_int(sname,tag,1,icr)
c      
c     solute permittivity
c     -----------------------------
      tag = "rism:permittivity"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) 
     >  f = 1.0d0
      tag = "solute_permittivity"
      call db_data_put_dbl(sname,tag,1,f)
c      
c     coulomb tail separation
c     -----------------------------
      tag = "rism:tau"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) 
     >  f=1.0d0
      tag = "tau"
      call db_data_put_dbl(sname,tag,1,f)
c      
c     tolerance
c     -----------------------------
      tag = "rism:tol"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) 
     >  f=1.0d-5
      tag = "tolerance"
      call db_data_put_dbl(sname,tag,1,f)
c      
c     mixing
c     -----------------------------
      tag = "rism:lambda"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) 
     >  f=0.0d0
      tag = "mixing"
      call db_data_put_dbl(sname,tag,1,f)
c      
c     temperature
c     -----------------------------
      tag = "rism:temp"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) 
     >  f=298.15d0
      tag = "temperature"
      call db_data_put_dbl(sname,tag,1,f)
c      
c     diis
c     -----------------------------
      tag = "rism:diis"
      if(.not.rtdb_get(rtdb,tag,mt_int,1,ipar)) 
     >  ipar = 9
      tag = "diis"
      write(*,*) "ipar=",ipar
      call db_data_put_int(sname,tag,1,ipar)
c
      call db_data_namespace_print(sname,55)
      return
      end

      subroutine rism_prepare_grid(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer rtdb
c
      character*(72) filename
      integer i
      integer fn
      character*72 tag
      character*72 buffer
      character*1 b1
      integer master
      integer nu
      character*32 sname
      double precision dr,r1,pi,dk
      integer i_kgrid,i_rgrid,ngrid
c 
      sname = "grid"
      call db_namespace_create(sname)
c      
      tag="rism:ngrid"
      if (.not.rtdb_get(rtdb,tag,mt_int,1,ngrid)) 
     >  call errquit("cannot get "//tag,0,0)
      tag = "npoints"
      call db_data_put_int(sname,tag,1,ngrid)
c      
c     get filename for solvent g(r)
c     -----------------------------
      tag = "rism:solvent:rdf"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
      call rism_gr_grid(filename,dr)
      tag = "dr"
      call db_data_put_dbl(sname,tag,1,dr)
c      
      pi=2*asin(1.0)
      dk=pi/(ngrid+1)/dr
      tag = "dk"
      call db_data_put_dbl(sname,tag,1,dk)
c      
      call db_data_create_get(sname,"rgrid",ngrid,mt_dbl,i_rgrid)
      call db_data_create_get(sname,"kgrid",ngrid,mt_dbl,i_kgrid)
c     the zero point is excluded
      r1=dr   
      do i=1,ngrid
       dbl_mb(i_rgrid+i-1)=r1+dr*(i-1)
       dbl_mb(i_kgrid+i-1)=i*dk
c       rgrid(i)=r1+dr*(i-1) 
c       kgrid(i)=i*dk 
      enddo 
      return
      end

      subroutine rism_prepare_solvent(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer rtdb
c
      character*(72) filename
      integer i
      integer fn
      character*72 tag
      character*72 buffer
      character*1 b1
      integer nv,nvv
      character*32 sname
      integer  i_epsiv
      integer  i_qqv
      integer  i_den
      integer  i_xv
      integer  i_yv
      integer  i_zv
      integer  i_sigv
      integer  i_isv
      integer  i_mv
      integer  i_tv
      integer  i_ims
      integer  i_sgvv
      integer  i_epsvv
      integer  i_qvv
c 
c     create container for solvent data
c     --------------------------------- 
      sname = "solvent"
      call db_namespace_create(sname)
c      
c     get solvent filename
c     --------------------
      tag = "rism:solvent:model"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
c      
c     get number of solvent atoms
c     ---------------------------
      call rism_nsolvent(nv,filename)
      tag = "natoms"
      call db_data_put_int(sname,tag,1,nv)
c      
      call db_data_create_get(sname,"atom_name",nv,MT_STR,i_tv)
      call db_data_create_get(sname,"residue_index",nv,mt_int,i_isv)
      call db_data_create_get(sname,"multiplicity",nv,mt_int,i_mv)
      call db_data_create_get(sname,"density",nv,mt_dbl,i_den)
      call db_data_create_get(sname,"xcoord",nv,mt_dbl,i_xv)
      call db_data_create_get(sname,"ycoord",nv,mt_dbl,i_yv)
      call db_data_create_get(sname,"zcoord",nv,mt_dbl,i_zv)
      call db_data_create_get(sname,"sigv",nv,mt_dbl,i_sigv)
      call db_data_create_get(sname,"epsiv",nv,mt_dbl,i_epsiv)
      call db_data_create_get(sname,"qqv",nv,mt_dbl,i_qqv)
c
      call rism_read_solvent_file(filename,nv,
     +                            byte_mb(i_tv),int_mb(i_isv),
     +                            int_mb(i_mv),dbl_mb(i_den),
     +                            dbl_mb(i_xv),dbl_mb(i_yv),
     +                            dbl_mb(i_zv),dbl_mb(i_sigv),
     +                            dbl_mb(i_epsiv),dbl_mb(i_qqv))

c     array that maps to reduced solvent parameters
      call db_data_create_get(sname,"map_reduced",nv,mt_int,i_ims)
c     calculate size and mapping for reduced solvent parameters
      call rism_sort(nv,byte_mb(i_tv),int_mb(i_isv),int_mb(i_mv),
     +               int_mb(i_ims),nvv)
      tag = "natoms_reduced"
      call db_data_put_int(sname,tag,1,nvv)
c     fill in reduced solvent parameters
      call db_data_create_get(sname,"sigma_reduced",nvv,mt_dbl,i_sgvv)
      call db_data_create_get(sname,"epsilon_reduced",
     +                        nvv,mt_dbl,i_epsvv)
      call db_data_create_get(sname,"charge_reduced",nvv,mt_dbl,i_qvv)
      call rism_vpot(nv,int_mb(i_ims),nvv,
     +          dbl_mb(i_sigv),dbl_mb(i_epsiv),
     +          dbl_mb(i_qqv),dbl_mb(i_sgvv),
     +          dbl_mb(i_epsvv),dbl_mb(i_qvv))
       do i=1,nvv
         write(33,*) dbl_mb(i_sgvv+i-1),dbl_mb(i_epsvv+i-1),
     +               dbl_mb(i_qvv+i-1),i,int_mb(i_ims+i-1)
       end do
      end

      subroutine rism_prepare_solute(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer rtdb
c
      character*(72) filename
      integer i
      integer fn
      character*72 tag
      character*72 buffer
      character*1 b1
      integer master
      integer nu
      character*32 sname
      logical result
c 
      integer  i_epsiu
      integer  i_qqu
      integer  i_xu
      integer  i_yu
      integer  i_zu
      integer  i_sigu
      integer  i_isu
      integer  i_mu
      integer  i_tu
      integer  i_wu
      integer  i_kgrid,ngrid
c
      sname = "solute"
      call db_namespace_create(sname)
c      
c     get solute filename
c     --------------------
      tag = "rism:solute"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
c      
c     get number of solute atoms
c     ---------------------------
      call rism_natoms(nu,filename)

      tag = "natoms"
      call db_data_put_int(sname,tag,1,nu)

      call db_data_create_get(sname,"tu",nu,MT_STR,i_tu)
      call db_data_create_get(sname,"isu",nu,mt_int,i_isu)
      call db_data_create_get(sname,"mu",nu,mt_int,i_mu)
      call db_data_create_get(sname,"xu",nu,mt_dbl,i_xu)
      call db_data_create_get(sname,"yu",nu,mt_dbl,i_yu)
      call db_data_create_get(sname,"zu",nu,mt_dbl,i_zu)
      call db_data_create_get(sname,"sigma",nu,mt_dbl,i_sigu)
      call db_data_create_get(sname,"epsilon",nu,mt_dbl,i_epsiu)
      call db_data_create_get(sname,"charge",nu,mt_dbl,i_qqu)
 
c
      call rism_read_solute_file(filename,nu,
     +                            byte_mb(i_tu),int_mb(i_isu),
     +                            int_mb(i_mu),
     +                            dbl_mb(i_xu),dbl_mb(i_yu),
     +                            dbl_mb(i_zu),dbl_mb(i_sigu),
     +                            dbl_mb(i_epsiu),dbl_mb(i_qqu))
c
c     compute solute structure factor
c     -------------------------------
      call db_data_get_ind_dim("grid","kgrid",
     +                           i_kgrid,ngrid,result)
      if(.not.result)
     >  call errquit("cannot get kgrid",0,0)

      call db_data_create_get(sname,"struct_factor",
     +                        nu*nu*ngrid,mt_dbl,i_wu)
      call rism_compute_solute_sfac(nu,ngrid,
     +                            dbl_mb(i_kgrid),
     +                            dbl_mb(i_xu),dbl_mb(i_yu),
     +                            dbl_mb(i_zu),
     +                            dbl_mb(i_wu))


      return
      end

      subroutine rism_vpot(nv,ims,nvv,sigv,epsiv,qqv,sgvv,epsvv,qvv) 
      implicit none
      integer  nv,nvv 
      real*8 sgvv(nvv),epsvv(nvv),qvv(nvv)        
      real*8 sigv(nv),epsiv(nv),qqv(nv)    
      integer  i,ims(nv),j1,j2,icr       
      do i=1,nv  
        sgvv(ims(i))=sigv(i)
        epsvv(ims(i))=epsiv(i)
        qvv(ims(i))=qqv(i)
!       print*, sgvv(ims(i)),sigv(i),i,ims(i)
      enddo 
      return
      end subroutine
c
      subroutine rism_compute_solute_sfac(nu,ngr,
     +                            kgrid,
     +                            xu,yu,
     +                            zu,
     +                            wu)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      integer nu,ngr
      double precision kgrid(ngr)
      double precision xu(nu),yu(nu),zu(nu),sigu(nu)
      double precision wu(nu,nu,ngr)
c
      integer i,j1,j2
      double precision dist
c      
      double precision sinr
      external sinr
c
      do i=1,ngr
       do j1=1,nu
        do j2=1,nu
         dist=((xu(j1)-xu(j2))**2+(yu(j1)-yu(j2))**2+
     *   (zu(j1)-zu(j2))**2)
         dist =SQRT(dist)
         wu(j1,j2,i)=sinr(kgrid(i)*dist)
        enddo
       enddo
      enddo 
      return
      end
c
      subroutine rism_sort(nv,tv,isv,mv,ims,nvv)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      character(db_string_size) tv(nv)
      integer  nv, isv(nv), mv(nv),ims(nv)
      integer  i,i1,j1,j2,icc,maxn,nvv      
c     automatic array change it ASAP
      integer mvt(nv),imm(1:nv)
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
       maxn=ims(1) 
       do i=2,nv
        if(ims(i).ge.maxn) then
         maxn=ims(i)
        endif
       enddo
       nvv=maxn
c      replace changed mv by the initial one 
       do i=1,nv
        mv(i)=mvt(i)
       enddo
       return
       end subroutine 
c
      subroutine rism_read_solvent_file(filename,nv,
     +                            tv,isv,
     +                            mv,den,
     +                            xv,yv,
     +                            zv,sigv,
     +                            epsiv,qqv)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      character*(*) filename
      character(db_string_size) tv(nv)
      integer  nv, isv(nv), mv(nv)
      real*8 den(nv),xv(nv),yv(nv),zv(nv),sigv(nv)
      real*8 epsiv(nv),qqv(nv) 
c
      integer i
      integer fn
      character*180 buffer
      character*1 b1
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      i=0
10    continue
      read(fn,'(A180)',ERR=40,END=40) buffer
      if(buffer(1:1).eq."#") goto 10
      i = i+1
      read(buffer,*,ERR=40,END=40) tv(i),isv(i),mv(i),den(i),  
     +   xv(i),yv(i),zv(i),sigv(i), epsiv(i),qqv(i) 
      if(i.lt.nv) goto 10 
30    close(fn)
      return
40    continue
      call errquit("error reading solvent file",0,0)
      end

      subroutine rism_read_solute_file(filename,nu,
     +                            tu,isu,
     +                            mu,
     +                            xu,yu,
     +                            zu,sigu,
     +                            epsiu,qqu)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      character*(*) filename
      character(db_string_size) tu(nu)
      integer  nu, isu(nu), mu(nu)
      real*8 xu(nu),yu(nu),zu(nu),sigu(nu)
      real*8 epsiu(nu),qqu(nu) 
c
      integer i
      integer fn
      character*180 buffer
      character*1 b1
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      i=0
10    continue
      read(fn,'(A180)',ERR=40,END=40) buffer
      if(buffer(1:1).eq."#") goto 10
      i = i+1
      read(buffer,*,ERR=40,END=40) tu(i),isu(i),mu(i),
     +   xu(i),yu(i),zu(i),sigu(i), epsiu(i),qqu(i) 
      if(i.lt.nu) goto 10 
30    close(fn)
      return
40    continue
      call errquit("error reading solute file",0,0)
      end

      function sinr(x)
      implicit none
      double precision sinr
      double precision x
      if(x.lt.1d-22)then
       sinr=1
      else
      sinr=sin(x)/x
      endif
      return
      end

      subroutine rism_nsolvent(nv,filename)
c find number of atoms in the solvent molecule
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer nv
      character*(*) filename
c
      integer i
      integer fn
      character*72 buffer
      character*1 b1
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      nv=0
10    continue        
      read(fn,*,ERR=30,END=30) b1
      if(b1.eq.' ') goto 10
      if(b1.eq."#") goto 10
      nv = nv+1
      goto 10
30    close(fn)
      return
      end

      subroutine rism_natoms(nt,filename)
c find number of atoms from the file
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      integer nt
      character*(*) filename
c
      integer i
      integer fn
      character*72 buffer
      character*1 b1
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      nt=0
10    continue        
      read(fn,*,ERR=30,END=30) b1
      if(b1.eq.' ') goto 10
      if(b1.eq."#") goto 10
      nt = nt+1
      goto 10
30    close(fn)
      return
      end

      subroutine rism_gr_grid(filename,dr)
c find number of atoms from the file
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      character*(*) filename
      double precision dr
c
      double precision r1,r2
      integer fn
      character*72 buffer
      character*1 b1
      character*30 pname
c      
      pname = "rism_gr_grid"
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      read(fn,*,ERR=30,END=30) r1
      read(fn,*,ERR=30,END=30) r2
      dr = r2-r1
      close(fn)
      return
30    call errquit(pname//"errors",0,0)
      end
