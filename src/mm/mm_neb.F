      subroutine mm_neb_init(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
#include "msgids.fh"
#include "tcgmsg.fh"
      
      integer rtdb

      logical orestrt
      logical ofile
      integer nbeads
      integer i, j
      character*255 filename
      character*255 crdparmfilename
c     character*255 qmcoordfilename
c     character*255 mmcoordfilename
c     character*255 bndparmfilename
c     character*255 parmfilename
c     character*255 bqfilename

      character*255 crdparmfile(2)
c     character*255 qmcoordfile(2)
c     character*255 mmcoordfile(2)
c     character*255 bndparmfile(2)
c     character*255 parmfile(2)
c     character*255 bqfile(2)
      character*30 pname
      logical custom_path

      pname = "mm_neb_init"
c     write(*,*) pname


      if(.not.rtdb_get(rtdb,'mm:neb:restart',mt_log,1,orestrt))
     >       orestrt = .false.
       
      if (.not.orestrt) then
        if(.not.rtdb_cget(rtdb,"mm:neb_path_limits",
     &                   2,crdparmfile)) goto 911
c       if(.not.rtdb_cget(rtdb,"mm:neb_path_limits_qm",
c    &                   2,qmcoordfile)) goto 911
c       if(.not.rtdb_cget(rtdb,"mm:neb_path_limits_mm",
c    &                   2,mmcoordfile)) goto 911
c       if(.not.rtdb_cget(rtdb,"mm:bndparmfile",
c    &                   2,bndparmfile)) goto 911
c       if(.not.rtdb_cget(rtdb,"mm:parmfile",
c    &                   2,parmfile)) goto 911
c       if(.not.rtdb_cget(rtdb,"mm:bqfile",
c    &                   2,bqfile)) goto 911
      
        call mm_neb_create_geom(rtdb,"neb_start",
     &                crdparmfile(1))
c    &                qmcoordfile(1),mmcoordfile(1),
c    &                bndparmfile(1),parmfile(1),
c    &                bqfile(1))
        call mm_neb_create_geom(rtdb,"neb_end",
     &                crdparmfile(2))
c    &                qmcoordfile(2),mmcoordfile(2),
c    &                bndparmfile(2),parmfile(2),
c    &                bqfile(2))
      else
c       write(*,*) 'This is a restart run'
c       custom_path = .true.
c       write(*,*) 'custom_path = ', custom_path
        if(.not.rtdb_put(rtdb,'neb:custom_path',mt_log,1,custom_path)) 
     >        call errquit(pname,0,RTDB_ERR)
        filename = ' '
        if(.not.rtdb_get(rtdb,'neb:nbeads',mt_int,1,nbeads))
     >    call errquit("error getting nbeads",0,-1)
        do i=1,nbeads
          write(filename,10) 'nwchem_',i
c         write(*,*) 'filename = ', filename(1:inp_strlen(filename))
          inquire(file=filename,exist=ofile)
          if(ofile) then
            if(ga_nodeid().eq.0) then
              write(*,*) 'mm:neb: found ', 
     >                  filename(1:inp_strlen(filename))
            end if
          else 
            call errquit("error: not found! "//
     >                    filename(1:inp_strlen(filename)),0,-1)
          end if

c         write(filename,11) 'nwchem_',i
c         write(*,*) 'filename = ', filename(1:inp_strlen(filename))
c         inquire(file=filename,exist=ofile)
c         if(ofile) then
c           write(*,*) 'mm:neb: found ', 
c    >                  filename(1:inp_strlen(filename))
c         else 
c           call errquit("error: not found! "//
c    >                    filename(1:inp_strlen(filename)),0,-1)
c         end if

c         write(filename,12) 'nwchem_',i
c         write(*,*) 'filename = ', filename(1:inp_strlen(filename))
c         inquire(file=filename,exist=ofile)
c         if(ofile) then
c           write(*,*) 'mm:neb: found ', 
c    >                  filename(1:inp_strlen(filename))
c         else 
c           call errquit("error: not found! "//
c    >                    filename(1:inp_strlen(filename)),0,-1)
c         end if

c         write(filename,13) 'nwchem_',i
c         write(*,*) 'filename = ', filename(1:inp_strlen(filename))
c         inquire(file=filename,exist=ofile)
c         if(ofile) then
c           write(*,*) 'mm:neb: found ', 
c    >                  filename(1:inp_strlen(filename))
c         else 
c           call errquit("error: not found! "//
c    >                    filename(1:inp_strlen(filename)),0,-1)
c         end if

c         write(filename,14) 'nwchem_',i
c         write(*,*) 'filename = ', filename(1:inp_strlen(filename))
c         inquire(file=filename,exist=ofile)
c         if(ofile) then
c           write(*,*) 'mm:neb: found ', 
c    >                  filename(1:inp_strlen(filename))
c         else 
c           call errquit("error: not found! "//
c    >                    filename(1:inp_strlen(filename)),0,-1)
c         end if

        end do
        crdparmfilename = ' '
c       qmcoordfilename = ' '
c       mmcoordfilename = ' '
c       bndparmfilename = ' '
c       parmfilename    = ' '
c       bqfilename      = ' '
        write(crdparmfilename,10) 'nwchem_', 1
c       write(qmcoordfilename,10) 'nwchem_', 1
c       write(mmcoordfilename,11) 'nwchem_', 1
c       write(parmfilename,   12) 'nwchem_', 1
c       write(bndparmfilename,13) 'nwchem_', 1
c       write(bqfilename,     14) 'nwchem_', 1

        call mm_neb_data_load(rtdb,crdparmfilename)
c    &                            qmcoordfilename,mmcoordfilename,
c    &                            bndparmfilename,parmfilename,
c    &                            bqfilename)

      end if
   
10    format(a,i3.3,'.crdparms')
c10    format(a,i3.3,'.qmcoords')
c11    format(a,i3.3,'.mmcoords')
c12    format(a,i3.3,'.parms')
c13    format(a,i3.3,'.bndparms')
c14    format(a,i3.3,'.bq')
   
      return
911   call errquit("error"//pname,0,-1)  

      end

      subroutine mm_neb_create_geom(rtdb,geomname,crdparmfile)
c    &                              qmcoordfile,mmcoordfile,
c    &                              bndparmfile,parmfile,
c    &                              bqfile)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer rtdb
      character*(*) geomname
      character*(*) crdparmfile
c     character*(*) qmcoordfile
c     character*(*) mmcoordfile
c     character*(*) bndparmfile
c     character*(*) parmfile
c     character*(*) bqfile
 
      character*30 pname

      pname = "mm_neb_create_geom"
c     write(*,*) pname

      if(.not.rtdb_cput(rtdb,"mm:crdparms:load",1,
     &                  crdparmfile(1:inp_strlen(crdparmfile))))
     & call errquit("error "//"mm:crdparmfile",0,-1)

c     if(.not.rtdb_cput(rtdb,"mm:qmcoord",1,
c    &                  qmcoordfile(1:inp_strlen(qmcoordfile))))
c    & call errquit("error "//"mm:qmcoordfile",0,-1)

c     if(.not.rtdb_cput(rtdb,"mm:mmcoord",1,
c    &                  mmcoordfile(1:inp_strlen(mmcoordfile))))
c    & call errquit("error "//"mm:mmcoordfile",0,-1)

c     if(.not.rtdb_cput(rtdb,"mm:bndparm",1,
c    &                  bndparmfile(1:inp_strlen(bndparmfile))))
c    & call errquit("error "//"mm:bndparmfile",0,-1)

c     if(.not.rtdb_cput(rtdb,"mm:parm",1,
c    &                  parmfile(1:inp_strlen(parmfile))))
c    & call errquit("error "//"mm:parmfile",0,-1)

c     if(.not.rtdb_cput(rtdb,"mm:bq",1,
c    &                  bqfile(1:inp_strlen(bqfile))))
c    & call errquit("error "//"mm:bqfile",0,-1)

      call mm_coords_init(rtdb)
      call mm_bndparm_init(rtdb)
      call mm_bond_coords_init()  
      call mm_vdw_init(rtdb)
      call mm_vdw_coords_init()
      call mm_links_init(rtdb)      
      call mm_geom_init(rtdb,geomname)      
      call mm_bq_init(rtdb)
c     call mm_vdw_init(rtdb)
c     call mm_vdw_coords_init()

      end 

      subroutine mm_neb_update_ibead(rtdb)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"
      integer rtdb

      logical orestrt
      integer ibead      
      integer nbeads
      character*255 crdparm(2)
      character*255 crdparmfile
c     character*255 qmcoordfile
c     character*255 mmcoordfile
c     character*255 parmfile
c     character*255 bndparmfile
c     character*255 bqfile
      character*30 pname

      pname = "mm_neb_update_ibead"
c     write(*,*) pname

      if(.not.rtdb_get(rtdb,'mm:neb:restart',mt_log,1,orestrt))
     >       orestrt = .false.

      if(.not.rtdb_get(rtdb,'neb:nbeads',mt_int,1,nbeads))
     >  call errquit('neb:nbeads get',0,-1)
c     write(*,*) 'nbeads = ', nbeads

      if(.not.rtdb_get(rtdb,'neb:ibead',mt_int,1,ibead))
     >  call errquit('neb:ibead get',0,-1)
c     write(*,*) 'ibead = ', ibead

      crdparmfile = ' '
c     qmcoordfile = ' '
c     mmcoordfile = ' '
c     parmfile    = ' '
c     bndparmfile = ' '
c     bqfile      = ' '

      if(.not.orestrt) then 
        if(.not.rtdb_cget(rtdb,"mm:neb_path_limits",
     &                   2,crdparm))
     &    call errquit("error getting crdparm "//pname,0 -1)
        if(ibead.ne.nbeads) then
          crdparmfile = crdparm(1)
        else
          crdparmfile = crdparm(2)
        end if

      else       
        write(crdparmfile,10) 'nwchem_',ibead
      end if

10    format (a,i3.3,'.crdparms') 
c     write(*,*) 'crdparmfile= ', crdparmfile

c     write(qmcoordfile,10) 'nwchem_',ibead
c10    format (a,i3.3,'.qmcoords') 
c      write(*,*) 'qmcoordfile= ', qmcoordfile

c      write(mmcoordfile,11) 'nwchem_',ibead
c11    format (a,i3.3,'.mmcoords') 
c      write(*,*) 'mmcoordfile= ', mmcoordfile

c      write(parmfile,12) 'nwchem_',ibead
c12    format (a,i3.3,'.parms') 
c      write(*,*) 'parmfile= ', parmfile

c      write(bndparmfile,13) 'nwchem_',ibead
c13    format (a,i3.3,'.bndparms') 
c      write(*,*) 'bndparmfile= ', bndparmfile

c      write(bqfile,14) 'nwchem_',ibead
c14    format (a,i3.3,'.bq') 
c      write(*,*) 'bqfile= ', bqfile
  
      call mm_neb_data_load(rtdb,crdparmfile)
c     call mm_neb_data_load(rtdb,
c    >                      qmcoordfile,mmcoordfile,
c    >                      bndparmfile,parmfile,
c    >                      bqfile)

c     update active atom coordinates
      call mm_geom_push_active(rtdb)

      end

      subroutine mm_neb_data_load(rtdb,crdparmfile)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer rtdb
      character*(*) crdparmfile
      logical scnb_adj
 
      character*30 pname

      pname = "mm_neb_data_load"
c     write(*,*) pname      

      if(.not.rtdb_cput(rtdb,"mm:crdparms:load",1,
     &                  crdparmfile(1:inp_strlen(crdparmfile))))
     & call errquit("error "//"mm:crdparmfile",0,-1)

c     write(*,*) 'crdparmfile= ', crdparmfile
      call mm_coords_init(rtdb)
      call mm_bndparm_init(rtdb)
      call mm_bond_coords_init()  
      call mm_vdw_init(rtdb)
      call mm_vdw_coords_init()
      call mm_links_init(rtdb)

c     adjust scnb if true
      scnb_adj = .true.
      if(.not.rtdb_get(rtdb,'mm:scnb_adjust',mt_log,1,scnb_adj)) 
     >  scnb_adj = .false.
      if(scnb_adj) call mm_vdw_adj_scnb()

      call mm_geom_init(rtdb,"none")      
      call mm_bq_init(rtdb)
c     call mm_vdw_init(rtdb)
c     call mm_vdw_coords_init()

      end 

      subroutine mm_neb_data_load2(rtdb,
     &                            qmcoordfile,mmcoordfile,
     &                            bndparmfile,parmfile,
     &                            bqfile)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "global.fh"
#include "inp.fh"

      integer rtdb
      character*(*) qmcoordfile
      character*(*) mmcoordfile
      character*(*) bndparmfile
      character*(*) parmfile
      character*(*) bqfile
 
      character*30 pname

      pname = "mm_neb_data_load2"
c     write(*,*) pname      

      if(.not.rtdb_cput(rtdb,"mm:qmcoord",1,
     &                  qmcoordfile(1:inp_strlen(qmcoordfile))))
     & call errquit("error "//"mm:qmcoordfile",0,-1)

      if(.not.rtdb_cput(rtdb,"mm:mmcoord",1,
     &                  mmcoordfile(1:inp_strlen(mmcoordfile))))
     & call errquit("error "//"mm:mmcoordfile",0,-1)

      if(.not.rtdb_cput(rtdb,"mm:bndparm",1,
     &                  bndparmfile(1:inp_strlen(bndparmfile))))
     & call errquit("error "//"mm:bndparmfile",0,-1)

      if(.not.rtdb_cput(rtdb,"mm:parm",1,
     &                  parmfile(1:inp_strlen(parmfile))))
     & call errquit("error "//"mm:parmfile",0,-1)

      if(.not.rtdb_cput(rtdb,"mm:bq",1,
     &                  bqfile(1:inp_strlen(bqfile))))
     & call errquit("error "//"mm:bqfile",0,-1)

      call mm_coords_init(rtdb)
      call mm_bndparm_init(rtdb)
      call mm_bond_coords_init()  
      call mm_links_init(rtdb)
      call mm_geom_init(rtdb,"none")      
      call mm_bq_init(rtdb)
      call mm_vdw_init(rtdb)
      call mm_vdw_coords_init()

      end

      subroutine mm_neb_create_xyz_file(luout,ibead)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "util.fh"
#include "global.fh"

      integer luout
      integer ibead

      logical value
      integer taglen,i,ii,geom,geomlen,nbeads,nion
      real*8  energy,q,rxyz(3)
      character*255 rtdb_name,geom_name
      character*2   symbol
      character*16  t,name
      character*50  tag

      integer bead_rtdb
      common /bead_list/ bead_rtdb
      
*     **** external functions ****
      integer     inp_strlen,size_bead_list
      character*7 bead_index_name
      external    inp_strlen,size_bead_list
      external    bead_index_name

      character*30 pname
      double precision scale_length
 
      pname = "mm_neb_create_xyz_file"
c     write(*,*) pname
    
      tag       = 'bead_list'
      taglen    = inp_strlen(tag)
      nbeads    = size_bead_list(tag)
      value = geom_create(geom,'bead_tmp') 

      rtdb_name = tag(1:taglen)//bead_index_name(ibead)//':geom_name'
c     write(*,*) 'rtdb_name = ', rtdb_name
      if (.not.rtdb_cget(bead_rtdb,rtdb_name,1,geom_name))
     >   geom_name = 'bead'//bead_index_name(ibead)//':geom'
      geomlen = inp_strlen(geom_name)      
      value = value.and.geom_rtdb_load(bead_rtdb,
     >                          geom,geom_name(1:geomlen))
      value = value.and.geom_ncent(geom,nion)

      call util_convert_units("au","ang",scale_length)

      if (ga_nodeid().eq.0) then
      write(luout,*) nion
      write(luout,*) 
      end if
      do ii=1,nion       
         value = value.and.geom_cent_get(geom,ii,t,rxyz,q)
         value = value.and. geom_tag_to_element(t,symbol,name,q)
c        rxyz(1)= rxyz(1)*scale_length
c        rxyz(2)= rxyz(2)*scale_length
c        rxyz(3)= rxyz(3)*scale_length
         rxyz(1)= rxyz(1)*0.529177d0
         rxyz(2)= rxyz(2)*0.529177d0
         rxyz(3)= rxyz(3)*0.529177d0
         if (ga_nodeid().eq.0) then
         write(luout,'(A2,6x,3F12.6)') symbol,rxyz
         end if
      end do
      value = value.and.geom_destroy(geom)

c     if (ga_nodeid().eq.0) then
c        write(luout,*) 
c     end if

      if (.not.value) 
     >  call errquit('create_xyz_file_bead_list failed',0,0)

      return
      end 
