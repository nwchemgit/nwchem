      subroutine smd_rdinp()
c
      implicit none
#include "errquit.fh"
#include "stdio.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "srtdb.fh"
c
      integer srtdb,rtdb
      character*32 tag
      character*32 pname
      character*255 token
      character*80 mtoken(10)
      integer  itoken(10)
      double precision  ftoken(10)
      integer ip,np
c
      pname = "smd_rdinp: "
c
c      write(luout,*) "in ",pname
c
      call smd_srtdb_get_handle(srtdb)
      call smd_rtdb_get_handle(rtdb)
c
c     coordinates
c     ----------
      tag="smd:coordfile"
      if (.not.rtdb_cget(rtdb,tag,1,mtoken(1))) 
     >     call errquit(pname//'failed to store'//tag,0,
     >      RTDB_ERR)
      
      if (.not.srtdb_cput(srtdb,tag,1,mtoken(1)))
     >     call errquit(pname//'failed to store'//tag,0,
     >      RTDB_ERR)
      

c
c     parameters
c     ----------
      tag="smd:paramfile"
      if (.not.rtdb_cget(rtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_cput(srtdb,tag,1,mtoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     kvec
c     -----
      tag="smd:kvec"
      if (.not.rtdb_get(rtdb,tag,mt_int,3,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_int,3,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     lat_a
c     -----
      tag="smd:lat_a"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     lat_b
c     -----
      tag="smd:lat_b"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     lat_c
c     -----
      tag="smd:lat_c"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,3,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     ndata
c     -----
      np = 1
      tag="smd:ndata"
      if (.not.rtdb_get(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     nequil
c     -----
      np = 1
      tag="smd:nequil"
      if (.not.rtdb_get(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     nprint
c     -----
      np = 1
      tag="smd:nprint"
      if (.not.rtdb_get(rtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_int,np,itoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     step
c     -----
      np = 1
      tag="smd:step"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     rcut
c     ----
      np = 1
      tag="smd:rcut"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     temp
c     -----
      np = 1
      tag="smd:temp_target"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     ewald
c     -----
      np = 1
      tag="smd:ewald"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
c
c     verlet
c     ------
      np = 1
      tag="smd:rcut_verlet"
      if (.not.rtdb_get(rtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
      
      if (.not.srtdb_put(srtdb,tag,mt_dbl,np,ftoken(1))) 
     >      call errquit(pname//'failed to store'//tag,0,
     >       RTDB_ERR)
          
c      if(.not.srtdb_print(srtdb,.true.))      
c     >      call errquit(pname//'failed to print srtdb',0,
c     >       RTDB_ERR)
c
      end


c $Id$
