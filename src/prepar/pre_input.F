      subroutine pre_input(irtdb)
c
c $Id: pre_input.F,v 1.8 1998-04-29 18:35:17 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item
c
      character*3 slvnam
      character*80 source
      character*10 slvmdl
      character*80 ffield,sysnam
      character*255 dir_s,dir_x,dir_u,dir_t,commnd
      integer newtop,newseq,newrst
c
      integer len,number,mgrid,mnoe,inoe(2),maxscf
      real*8 touch,xpndw,rgrid,edge,edges(3),dnoe(5),qscale
      character*4  cnoe(2)
c
c      if(ga_nodeid().eq.0) then
c      if(util_print('prepare',print_medium)) then
c      write(6,*)
c      write(6,*)
c      call util_print_centered(6,'PREPARE Input Module',40,.true.)
c      write(6,*)
c      write(6,*)
c      endif
c      endif
c
c     set the defaults
c     ----------------
c
      ffield=' '
      dir_s=' '
      dir_x=' '
      dir_u=' '
      dir_t=' '
      sysnam=' '
      slvnam='HOH'
      slvmdl='spce      '
      newtop=0
      newseq=0
      newrst=0
      source=' '
      mgrid=24
      rgrid=0.2d0
      mnoe=0
      maxscf=0
      qscale=-1.2d0
c
      touch=0.25
      xpndw=0.10
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,commnd))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      commnd=commnd(1:index(commnd,' ')-1)//'.cmd '
      open(unit=68,file=commnd(1:index(commnd,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=68)
      number=0
c
c     read the prepare input block
c     ----------------------------
c
    1 continue
c
      if(.not.inp_read()) call errquit('pre_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
c     print & noprint
c     ---------------
c
      if(inp_compare(.false.,'print',item).or.
     + inp_compare(.false.,'noprint',item)) then
      call util_print_input(irtdb,'prepare')
      goto 1
      endif
c
c
c     force field
c     -----------
c
      if(inp_compare(.false.,'amber',item)) then
      ffield='amber '
      goto 2
      endif
      if(inp_compare(.false.,'charmm',item)) then
      ffield='charmm '
      goto 2
      endif
      if(inp_compare(.false.,'gromos',item)) then
      ffield='gromos '
      goto 2
      endif
      if(inp_compare(.false.,'oplsa',item)) then
      ffield='oplsa '
      goto 2
      endif
      if(inp_compare(.false.,'oplsg',item)) then
      ffield='oplsg '
      goto 2
      endif
c
c     directories
c     -----------
c
      if(inp_compare(.false.,'standard',item)) then
      if(.not.inp_a(dir_s)) call errquit('pre_input: error dir_s',0)
      goto 2
      endif
      if(inp_compare(.false.,'extensions',item)) then
      if(.not.inp_a(dir_x)) call errquit('pre_input: error dir_x',0)
      goto 2
      endif
      if(inp_compare(.false.,'user',item)) then
      if(.not.inp_a(dir_u)) call errquit('pre_input: error dir_u',0)
      goto 2
      endif
      if(inp_compare(.false.,'temporary',item)) then
      if(.not.inp_a(dir_t)) call errquit('pre_input: error dir_t',0)
      goto 2
      endif
c
c     sysnam
c     ------
c
      if(inp_compare(.false.,'system',item)) then
      if(.not.inp_a(sysnam)) call errquit('pre_input: error sysnam',0)
      goto 2
      endif
c
c     source
c     ------
c
      if(inp_compare(.false.,'source',item)) then
      if(.not.inp_a(source)) call errquit('pre_input: error source',0)
      goto 2
      endif
c
c     max atoms in scf and charge scaling
c     -----------------------------------
c
      if(inp_compare(.false.,'maxscf',item)) then
      if(.not.inp_i(maxscf)) call errquit('pre_input: error maxscf',0)
      goto 2
      endif
      if(inp_compare(.false.,'qscale',item)) then
      if(.not.inp_f(qscale)) call errquit('pre_input: error qscale',0)
      goto 2
      endif
c
c     slvnam & slvmdl
c     ---------------
c
      if(inp_compare(.false.,'solvent',item)) then
   20 continue
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      if(inp_compare(.false.,'name',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      slvnam=item(1:3)
      goto 20
      else
      if(inp_compare(.false.,'model',item)) then
      if(.not.inp_a(slvmdl)) call errquit('pre_input: solvent',9999)
      goto 20
      endif
      endif
      endif
c
c     topology
c     --------
c
      if(inp_compare(.false.,'topology',item)) then
      newtop=1
   21 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'sequence',item)) then
      newseq=1
      goto 2
      endif
      goto 3
      endif
c
c     newrst
c     ------
c
      if(inp_compare(.false.,'newrst',item)) then
      newrst=1
      goto 2
      endif
c
c     counterion grid parameters
c     --------------------------
c
      if(inp_compare(.false.,'grid',item)) then
      if(.not.inp_i(mgrid)) call errquit('pre_input: grid',9999)
      if(.not.inp_f(rgrid)) call errquit('pre_input: grid',9999)
      goto 2
      endif
c
c     noe's
c     -----
c
      if(inp_compare(.false.,'noe',item)) then
      if(.not.inp_i(inoe(1))) call errquit('pre_input: noe',9999)
      if(.not.inp_a(cnoe(1))) call errquit('pre_input: noe',9999)
      if(.not.inp_i(inoe(2))) call errquit('pre_input: noe',9999)
      if(.not.inp_a(cnoe(2))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(1))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(2))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(3))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(4))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(5))) call errquit('pre_input: noe',9999)
      mnoe=mnoe+1
      write(68,6820) inoe,cnoe,dnoe
 6820 format('disres',2i5,2a4,3f12.6,2e12.5)
      number=number+1
      goto 2
      endif
c
c     touch
c     -----
c
      if(inp_compare(.false.,'touch',item)) then
      if(.not.inp_f(touch)) call errquit('pre_input: touch',9999)
      write(68,6801) touch
 6801 format('touch ',f12.6)
      number=number+1
      goto 2
      endif
c
c     periodic
c     --------
c
      if(inp_compare(.false.,'periodic',item)) then
      write(68,6802)
 6802 format('period ')
      number=number+1
      goto 2
      endif
c
c     vacuum
c     ------
c
      if(inp_compare(.false.,'vacuo',item)) then
      write(68,6803)
 6803 format('vacuum ')
      number=number+1
      goto 2
      endif
c
c     center
c     ------
c
      if(inp_compare(.false.,'center',item)) then
      write(68,6818)
 6818 format('center ')
      number=number+1
      goto 2
      endif
c
c     orient
c     ------
c
      if(inp_compare(.false.,'orient',item)) then
      write(68,6819)
 6819 format('orient ')
      number=number+1
      goto 2
      endif
c
c     solvate
c     ------
c
      if(inp_compare(.false.,'solvate',item)) then
      write(68,6804) touch
 6804 format('touch ',f12.6)
      if(.not.inp_a(item)) then
      write(68,6805) xpndw
 6805 format('rdcoor',/,'period ',/,'solvnt',/,'solvat cube ',/,
     + 'expand with ',f12.6,/,'wr_rst')
      number=number+1
      goto 1
      else
      if(inp_compare(.false.,'cube',item)) then
      if(.not.inp_f(edge)) then
      write(68,6811)
 6811 format('period',/,'solvat cube')
      else
      write(68,6812) edge
 6812 format('period',/,'solvat cube',f12.6)
      endif
      elseif(inp_compare(.false.,'box',item)) then
      if(.not.inp_f(edges(1))) then
      write(68,6813)
 6813 format('period',/,'solvat box')
      else
      if(.not.inp_f(edges(2))) edges(2)=edges(1)
      if(.not.inp_f(edges(2))) edges(3)=edges(2)
      write(68,6814) edges
 6814 format('period',/,'solvat box',3f12.6)
      endif
      elseif(inp_compare(.false.,'sphere',item)) then
      if(.not.inp_f(edge)) then
      write(68,6815)
 6815 format('vacuum',/,'solvat sphere')
      else
      write(68,6816) edge
 6816 format('vacuum',/,'solvat sphere',f12.6)
      endif
      endif
      endif
      number=number+1
      goto 2
      endif
c
c     write
c     -----
c
      if(inp_compare(.false.,'write',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: write',9999)
      if(inp_compare(.false.,'solute',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: write',9999)
      if(.not.inp_compare(.false.,'pdb',item))
     + call errquit('pre_input: write solute',9999)
      if(.not.inp_a(item)) then
      write(68,6806) ' '
      else
      write(68,6806) item(1:index(item,' ')-1)
 6806 format('wr_pdb solute ',a)
      endif
      elseif(inp_compare(.false.,'pdb',item)) then
      if(.not.inp_a(item)) then
      write(68,6807) ' '
      else
      write(68,6807) item(1:index(item,' ')-1)
 6807 format('wr_pdb ',a)
      endif
      elseif(inp_compare(.false.,'rst',item)) then
      if(.not.inp_a(item)) then
      write(68,6808) ' '
      else
      write(68,6808) item(1:index(item,' ')-1)
 6808 format('wr_rst ',a)
      endif
      endif
      number=number+1
      goto 2
      endif
c
c     read
c     ----
c
      if(inp_compare(.false.,'read',item)) then
      if(.not.inp_a(item)) then
      write(68,6817)
 6817 format('rdcoor')
      else
      if(inp_compare(.false.,'pdb',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: read',9999)
      write(68,6809) item(1:index(item,' ')-1)
 6809 format('rd_pdb ',a)
      elseif(inp_compare(.false.,'rst',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: read',9999)
      write(68,6810) item(1:index(item,' ')-1)
 6810 format('rd_rst ',a)
      endif
      endif
      number=number+1
      goto 2
      endif
c
c     end
c     ---
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,1000) item
 1000 format('pre_input: unknown field ',a)
      call errquit('pre_input: error reading input file',0)
      endif
c
      len=index(dir_s,' ')-1
      if(len.gt.0.and.dir_s(len:len).ne.'/') dir_s(len+1:len+2)='/ '
      len=index(dir_x,' ')-1
      if(len.gt.0.and.dir_x(len:len).ne.'/') dir_x(len+1:len+2)='/ '
      len=index(dir_u,' ')-1
      if(len.gt.0.and.dir_s(len:len).ne.'/') dir_u(len+1:len+2)='/ '
      len=index(dir_t,' ')-1
      if(len.gt.0.and.dir_t(len:len).ne.'/') dir_t(len+1:len+2)='/ '
c
c     put on rtdb
c     -----------
c
      if(ffield(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:ffield',1,ffield))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(dir_s(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:dir_s',1,dir_s))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(dir_x(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:dir_x',1,dir_x))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(dir_u(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:dir_u',1,dir_u))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(dir_t(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:dir_t',1,dir_t))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(sysnam(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:sysnam',1,sysnam))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(source(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:source',1,source))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(slvnam(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:slvnam',1,slvnam))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(slvmdl(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'prep:slvmdl',1,slvmdl))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(newtop.ne.0) then
      if(.not.rtdb_put(irtdb,'prep:newtop',mt_int,1,newtop))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(newseq.ne.0) then
      if(.not.rtdb_put(irtdb,'prep:newseq',mt_int,1,newseq))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(newrst.ne.0) then
      if(.not.rtdb_put(irtdb,'prep:newrst',mt_int,1,newrst))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(mgrid.gt.0) then
      if(.not.rtdb_put(irtdb,'prep:mgrid',mt_int,1,mgrid))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:rgrid',mt_dbl,1,rgrid))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(maxscf.gt.0) then
      if(.not.rtdb_put(irtdb,'prep:maxscf',mt_int,1,maxscf))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(qscale.gt.0.0d0) then
      if(.not.rtdb_put(irtdb,'prep:qscale',mt_dbl,1,qscale))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(number.gt.0) then
      close(unit=68)
      else
      close(unit=68,status='delete')
      endif
c
      return
c
 9999 continue
      call errquit('pre_input failed',9999)
      return
      end
