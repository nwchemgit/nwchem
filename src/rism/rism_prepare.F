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
      integer master
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
      sname = "rism"
      master = 0
c      
c     get solvent filename
c     --------------------
      tag = "rism:solvent"
      if(.not.rtdb_cget(rtdb,"rism:solvent",1,filename)) 
     >  call errquit("cannot get "//tag,0,0)
c      
c     get number of solvent atoms
c     ---------------------------
      if(ga_nodeid().eq.master) then
        call rism_nsolvent(nv,filename)
      end if
      call ga_brdcst(msg_smd,nv,
     >     ma_sizeof(mt_int,1,mt_byte),master)

      tag = "nsolvent"
      call db_data_put_int(sname,tag,1,nv)
c      
      call db_data_create_get(sname,"tv",nv,MT_STR,i_tv)
      call db_data_create_get(sname,"isv",nv,mt_int,i_isv)
      call db_data_create_get(sname,"mv",nv,mt_int,i_mv)
      call db_data_create_get(sname,"den",nv,mt_dbl,i_den)
      call db_data_create_get(sname,"xv",nv,mt_dbl,i_xv)
      call db_data_create_get(sname,"yv",nv,mt_dbl,i_yv)
      call db_data_create_get(sname,"zv",nv,mt_dbl,i_zv)
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
      call db_data_create_get(sname,"ims",nv,mt_int,i_ims)
c     calculate size and mapping for reduced solvent parameters
      call rism_sort(nv,byte_mb(i_tv),int_mb(i_isv),int_mb(i_mv),
     +               int_mb(i_ims),nvv)
      tag = "nsolvent_red"
      call db_data_put_int(sname,tag,1,nvv)
c     fill in reduced solvent parameters
      call db_data_create_get(sname,"sgvv",nvv,mt_dbl,i_sgvv)
      call db_data_create_get(sname,"epsvv",nvv,mt_dbl,i_epsvv)
      call db_data_create_get(sname,"qvv",nvv,mt_dbl,i_qvv)
      call rism_vpot(nv,int_mb(i_ims),nvv,
     +          dbl_mb(i_sigv),dbl_mb(i_epsiv),
     +          dbl_mb(i_qqv),dbl_mb(i_sgvv),
     +          dbl_mb(i_epsvv),dbl_mb(i_qvv))
c      call vpot(nv,ims,nvv,sigv,epsiv,qqv,sgvv,epsvv,qvv)
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
