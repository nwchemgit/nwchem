#define RISM_PREP_DEBUG1
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
      logical omaster
      integer master
      character*80 pname

      master = 0
      pname = "rism_prepare"

#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif      
      omaster = ga_nodeid().eq.master
c
c     order may be important here
c     ---------------------------
      if(omaster) then
       call rism_prepare_params(rtdb)
       call rism_prepare_grid(rtdb)
       call rism_prepare_solvent(rtdb)
       call rism_prepare_solute(rtdb)
c       call db_data_print(6)
      end if
#ifdef RISM_PREP_DEBUG      
      call rism_message("1 "//pname)
#endif
c
c     send data to other processors
c     -----------------------------
#ifdef RISM_PREP_DEBUG      
      call rism_message("2 "//pname)
#endif      
      call ga_sync()
#ifdef RISM_PREP_DEBUG      
      call rism_message("3 "//pname)
#endif      
      call db_data_brdcst_create(master)
#ifdef RISM_PREP_DEBUG      
      call rism_message("4 "//pname)
#endif      
#ifdef RISM_PREP_DEBUG      
      call rism_message("out "//pname)
#endif      
      return
 
      end

      subroutine rism_cleanup()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "inp.fh"
#include "util.fh"
#include "rtdb.fh"
#include "global.fh"
#include "db_const_data.fh"
      integer master
      character*80 pname

      pname = "rism_cleanup"
      call db_namespace_destroy("solute") 
      call db_namespace_destroy("solvent") 
      call db_namespace_destroy("grid") 
      call db_namespace_destroy("parameters") 
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
      character*80 pname
c
      pname = "rism_prepare_params"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif
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
     >  tag="arithmetic"
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
c     solvent permittivity
c     -----------------------------
      tag = "rism:solvent:permittivity"
      if(.not.rtdb_get(rtdb,tag,mt_dbl,1,f)) then
        f = 1.0d0
#ifdef RISM_PREP_DEBUG      
        call rism_message(pname//
     +       " setting default for solvent permittivity")
#endif      
      end if
      tag = "solvent_permittivity"
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
      call db_data_put_int(sname,tag,1,ipar)
c
      call db_data_namespace_print(sname,55)
#ifdef RISM_PREP_DEBUG      
      call rism_message("out "//pname)
#endif      
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
      double precision d
      character*80 pname
c
      logical okspace
      character*256 compiled_name
      logical rism_util_power_2
      external rism_util_power_2
c 
      pname = "rism_prepare_grid"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif      
      call util_nwchem_srcdir(compiled_name)
      compiled_name=compiled_name(1:inp_strlen(compiled_name))
     $     //"/rism/data/"

      sname = "grid"
      call db_namespace_create(sname)
c      
c      tag="rism:ngrid"
c      if (.not.rtdb_get(rtdb,tag,mt_int,1,ngrid)) 
c     >  call errquit("cannot get "//tag,0,0)
c      tag = "npoints"
c      call db_data_put_int(sname,tag,1,ngrid)
c      
c     get filename for solvent g(r)
c     -----------------------------
      tag = "rism:solvent:rdf"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) then
C     >  call errquit("cannot get "//tag,0,0)
        filename=compiled_name(1:inp_strlen(compiled_name))//"water.rdf"
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,filename)) 
     >  call errquit("cannot put "//tag,0,0)
      end if

      write(*,*) "filename=",filename
      call rism_grid_parse(filename,ngrid,d)
c     
c     check if grid size is power of 2
c     --------------------------------
      if(rism_util_power_2(ngrid)) then
#ifdef RISM_DEBUG              
        call rism_message("grid is power of 2")
#endif        
      else
        call errquit(pname//"grid size has to be power of 2",0,0)
      end if
      tag = "npoints"
      call db_data_put_int(sname,tag,1,ngrid)
c
c     check if rdf is given in kspace
c     -------------------------------
      tag="rism:solvent:kspace"
      if (.not.rtdb_get(rtdb,tag,mt_log,1,okspace)) then
        okspace = .false.
      end if
      tag = "okspace"
      call db_data_put_log(sname,tag,1,okspace)
      call db_data_create_get(sname,"rgrid",ngrid,mt_dbl,i_rgrid)
      call db_data_create_get(sname,"kgrid",ngrid,mt_dbl,i_kgrid)
      pi=2*asin(1.0)
      if(okspace) then
         dk = d
         dr = pi/ngrid/dk
      else
         dr = d
         dk = pi/ngrid/dr
      end if
      tag = "dr"
      call db_data_put_dbl(sname,tag,1,dr)
      tag = "dk"
      call db_data_put_dbl(sname,tag,1,dk)
c      
c     the zero point is excluded
      r1=dr   
      do i=1,ngrid
       dbl_mb(i_rgrid+i-1)=r1+dr*(i-1)
       dbl_mb(i_kgrid+i-1)=i*dk
       write(44,*) i,dbl_mb(i_rgrid+i-1),dbl_mb(i_kgrid+i-1)
      enddo 
#ifdef RISM_PREP_DEBUG      
      call rism_message("out "//pname)
#endif      
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
      character*80 pname 
      character*256 compiled_name
c
      pname = "rism_prepare_solvent"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif
c 
c     create container for solvent data
c     --------------------------------- 
      sname = "solvent"
      call db_namespace_create(sname)
c      
      call util_nwchem_srcdir(compiled_name)
      compiled_name=compiled_name(1:inp_strlen(compiled_name))
     $     //"/rism/data/"
c      
c     get solvent filename
c     --------------------
      tag = "rism:solvent:model"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) then
c     >  call errquit("cannot get "//tag,0,0)
        filename=compiled_name(1:inp_strlen(compiled_name))
     $//"water.conf"
        if (.not.rtdb_cput(rtdb,tag(1:inp_strlen(tag)),1,filename)) 
     >  call errquit("cannot put "//tag,0,0)
      end if
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
#ifdef RISM_PREP_DEBUG      
      call rism_message("out "//pname)
#endif      
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
      integer  i_tt
      integer  i_wu
      integer  i_kgrid,ngrid
      character*80 pname 
c
      pname = "rism_prepare_solute"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif
c
      sname = "solute"
      call db_namespace_create(sname)
c      
c     get solute filename
c     --------------------
      tag = "rism:solute:configuration"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
c      
c     get number of solute atoms
c     ---------------------------
c      call rism_natoms(nu,filename)
      call rism_natoms_xyz(nu,filename)

      tag = "natoms"
      call db_data_put_int(sname,tag,1,nu)

      call db_data_create_get(sname,"tu",nu,MT_STR,i_tu)
      call db_data_create_get(sname,"type",nu,MT_STR,i_tt)
c      isu and mu arrays are not needed
c      keeping them only for compatibility with original version
      call db_data_create_get(sname,"isu",nu,mt_int,i_isu)
      call db_data_create_get(sname,"mu",nu,mt_int,i_mu)
c      
      call db_data_create_get(sname,"xu",nu,mt_dbl,i_xu)
      call db_data_create_get(sname,"yu",nu,mt_dbl,i_yu)
      call db_data_create_get(sname,"zu",nu,mt_dbl,i_zu)
      call db_data_create_get(sname,"sigma",nu,mt_dbl,i_sigu)
      call db_data_create_get(sname,"epsilon",nu,mt_dbl,i_epsiu)
      call db_data_create_get(sname,"charge",nu,mt_dbl,i_qqu)
 
c     legacy format
c      call rism_read_solute_file(filename,nu,
c     +                            byte_mb(i_tu),int_mb(i_isu),
c     +                            int_mb(i_mu),
c     +                            dbl_mb(i_xu),dbl_mb(i_yu),
c     +                            dbl_mb(i_zu),dbl_mb(i_sigu),
c     +                            dbl_mb(i_epsiu),dbl_mb(i_qqu))
      call rism_read_solute_file_xyz(filename,nu,
     +                            byte_mb(i_tu),
     +                            dbl_mb(i_xu),dbl_mb(i_yu),
     +                            dbl_mb(i_zu),
     +                            dbl_mb(i_qqu),
     +                            byte_mb(i_tt))
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

c
c     fill in vdw parameters
c     ---------------------
      tag = "rism:vdw:parameters"
      if(.not.rtdb_cget(rtdb,tag,1,filename)) 
     >  call errquit("cannot get "//tag,0,0)

      call rism_read_vdw_params(filename,nu,
     +                            byte_mb(i_tt),
     +                            dbl_mb(i_sigu),
     +                            dbl_mb(i_epsiu))
c
      call db_data_destroy(sname,"isu",result)
      if(.not.result)
     >  call errquit("cannot destroy isu",0,0)
      call db_data_destroy(sname,"mu",result)
      if(.not.result)
     >  call errquit("cannot destroy mu",0,0)
c      
#ifdef RISM_PREP_DEBUG      
      call rism_message("out "//pname)
#endif

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
C>    Calculates number of unique atoms (nvv) 
C>    and enumerates them (ism)
C>    @param[in] nv number of solvent atoms      
C>    @param[in] tv names of solvent atoms      
C>    @param[in] isv residue/segment index of solvent atoms      
C>    @param[in] multiplicites/populations of solvent atoms
C>    @param[out] enumerated array of unique solvent atoms
C>    in the order of occurence
C>    @param[out] total number of unique solvent atoms
      subroutine rism_sort(nv,tv,isv,mv,ims,nvv)
      implicit none
#include "errquit.fh"
#include "db_const_data.fh"
      character(db_string_size) tv(nv)
      integer  nv, isv(nv), mv(nv)
      integer  ims(nv),nvv      
c     local variables
      integer  j1,j2
      integer ims_cur,mvtemp
      character*82 message
c      
      do j1=1,nv
        ims(j1)=0 
      enddo
c
      nvv = 0
      ims_cur=0
      do j1=1,nv
        if(ims(j1).eq.0) then
          mvtemp = 1
          ims_cur = ims_cur + 1
          ims(j1) = ims_cur
          nvv = nvv + 1
          do j2=j1+1,nv
            if((tv(j1).eq.tv(j2)).and.(isv(j1).eq.isv(j2))
     *         .and.(mv(j2).eq.mv(j1))) then 
              ims(j2)=ims_cur
              mvtemp = mvtemp + 1
            end if
          end do
          if(mvtemp.ne.mv(j1)) then
            write(message,'("inconsistent multiplicity ",A4,2(I4,2X))') 
     >           tv(j1),isv(j1), mv(j1)
            call errquit(message,0,0)
          end if
        end if
      end do

c      write(*,*) "my ims",nvv
c      do j1=1,nv
c        write(*,*) j1,ims(j1)
c      enddo
c
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
     +                            tu,
     +                            xu,yu,
     +                            zu,
     +                            qqu)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      character*(*) filename
      character(db_string_size) tu(nu)
      integer  nu
      real*8 xu(nu),yu(nu),zu(nu)
      real*8 qqu(nu) 
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
      read(buffer,*,ERR=40,END=40) tu(i),
     +   xu(i),yu(i),zu(i),qqu(i) 
      if(i.lt.nu) goto 10 
30    close(fn)
      return
40    continue
      call errquit("error reading solute file",0,0)
      end

      subroutine rism_read_solute_file_xyz(filename,nu,
     +                            tu,
     +                            xu,yu,
     +                            zu,
     +                            qqu,
     +                            tt)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      character*(*) filename
      character(db_string_size) tu(nu)
      integer  nu
      real*8 xu(nu),yu(nu),zu(nu)
      real*8 qqu(nu) 
      character(db_string_size) tt(nu)
c
      integer i
      integer nu1
      integer fn
      character*180 buffer
      character*1 b1
      character*72 message
      character*30 pname
c
      pname = "rism_read_solute_file_xyz"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif      
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      message = "reading number of atoms"
      read(fn,'(A180)',ERR=40,END=40) buffer
      read(buffer,*,ERR=40,END=40) nu1
      if(nu.ne.nu1) then
        message = "incompatible number of atoms:"
        goto 40
      end if
c      
      message = "reading title line"
      read(fn,'(A180)',ERR=40,END=40) buffer
c      
      message = "reading coordinates"
      do i=1,nu
        read(fn,'(A180)',ERR=40,END=40) buffer
        read(buffer,*,ERR=40,END=40) tu(i),
     +     xu(i),yu(i),zu(i),qqu(i),tt(i)
#ifdef RISM_DEBUG      
        write(*,*) tu(i),
     +     xu(i),yu(i),zu(i),qqu(i),tt(i)
#endif      
      end do
      close(fn)
      return
40    continue
      call errquit(message(1:inp_strlen(message))//pname,0,0)
      end

      subroutine rism_read_solute_file_legacy(filename,nu,
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

      subroutine rism_read_vdw_params(filename,nu,
     +                            tu,
     +                            sigu,
     +                            epsiu)

      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "db_const_data.fh"
      character*(*) filename
      character(db_string_size) tu(nu)
      integer  nu 
      real*8 epsiu(nu),sigu(nu) 
c
      integer i,ia
      integer fn
      character*180 buffer
      character*5 tag
      double precision sigma, eps
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
c     loop over atom names
      do ia=1,nu
10      continue
        read(fn,'(A180)',ERR=40,END=40) buffer
        if(buffer(1:1).eq."#") goto 10
        read(buffer,*,ERR=40,END=40) tag,sigma,eps
        if (inp_compare(.false.,tag, tu(ia))) then
           sigu(ia) = sigma
           epsiu(ia) =eps
        else
           goto 10
        end if
        rewind(fn)
      end do
30    close(fn)
      return
40    continue
      call errquit("error reading vdw parameters",0,0)
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

      subroutine rism_natoms_xyz(nt,filename)
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
      character*30 pname
c
      pname = "rism_natoms_xyz"
#ifdef RISM_PREP_DEBUG      
      call rism_message("in "//pname)
#endif      
c      
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
      read(fn,*,ERR=30,END=30) nt
      close(fn)
      return
30    continue
      call errquit("cannot get number of solute atoms",0,0)
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

      subroutine rism_grid_parse(filename,n,d)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
      character*(*) filename
      integer n
      double precision d
c
      double precision g1,g2
      integer fn
      character*72 buffer
      character*1 b1
      character*30 pname
c      
      pname = "rism_grid_parse"
      if(.not.util_get_io_unit(fn))
     >  call errquit("cannot get file number",0,0)
      open(unit=fn,status="old",form="formatted",file=filename)
c  
c     first determine spacing
c     -----------------------
      read(fn,*,ERR=30,END=31) g1
      read(fn,*,ERR=30,END=31) g2
      d = g2-g1
      n=2
c  
c     now determine number of points
c     -----------------------------
28    continue
      read(fn,*,ERR=30,END=29) g1
      n=n+1
      goto 28
29    continue
      close(fn)
      return
30    call errquit(pname//"error reading file "//filename,0,0)
31    call errquit(pname//"premature end of file "//filename,0,0)
      end
c $Id$
