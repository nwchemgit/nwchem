      subroutine ana_input(irtdb)
c
c $Id: ana_input.F,v 1.7 1999-10-25 17:44:26 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item,select,filtag,fpref,project
      character*255 filref,filtrj,filcop,filsup
      integer lcop,lsup,lrms
      character*3 fmtcop,fmtsup
c
      integer irmsd
      integer isel,jsel,numsel,iatag,jatag,numtag
      character*255 atom
      real*8 timoff
      integer ifrfr,ifrto,mcopf,msupf,ifrst,ilast
c
      real*8 rcut,spac,rtag,rsel
      integer lplt,itag
      character*255 filchg,filplt,filxyz
      character*24 tag
c
      irmsd=0
      lcop=0
      lsup=0
      ifrfr=1
      ifrto=0
      mcopf=0
      msupf=0
      ifrst=0
      ilast=0
      timoff=0.0d0
c
      lplt=0
      rcut=0.3d0
      spac=0.05d0
      rtag=0.0d0
      rsel=0.0d0
c
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,'Analysis Input Module',40,.true.)
      endif
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,fpref))
     + call errquit('rtdb_cget failed on file_prefix',9999)
c
      select=fpref(1:index(fpref,' ')-1)//'.sel '
      filtag=fpref(1:index(fpref,' ')-1)//'.tag '
c
c     default reference taken from ${file_prefix}.rst
c
      filref=fpref(1:index(fpref,' ')-1)//'.rst '
c
c     default trajectory taken from ${file_prefix}.trj
c
      filtrj=fpref(1:index(fpref,' ')-1)//'.trj '
c
c     default esp charges take from ${file_prefix}.q
c
      filchg=' '
      filplt=' '
      filxyz=' '
c
      filcop=' '
      fmtcop='trj'
c
      filsup=' '
      fmtsup='trj'
c
      open(unit=69,file=select(1:index(select,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=69)
      numsel=0
c
      open(unit=68,file=filtag(1:index(filtag,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=68)
      numtag=0
c
    1 continue
c
      if(.not.inp_read()) call errquit('ana_input: inp_read failed',0)
c
    2 continue
c
      if(.not.inp_a(item)) goto 1
c
    3 continue
c
      if(inp_compare(.false.,'system',item)) then
      if(.not.inp_a(project)) call errquit('ana_input: error system',0)
      if(filchg(1:1).eq.' ') then
      filchg=project(1:index(project,' ')-1)//'.q '
      endif
      if(filplt(1:1).eq.' ') then
      filplt=project(1:index(project,' ')-1)//'.plt '
      endif
      if(filxyz(1:1).eq.' ') then
      filxyz=project(1:index(project,' ')-1)//'.xyz '
      endif
      goto 2
      endif
c
c     solute root mean square devation
c
      if(inp_compare(.false.,'rmsdev',item)) then
      irmsd=1
      goto 2
      endif
c
      if(inp_compare(.false.,'reference',item)) then
      if(.not.inp_a(filref))
     + call errquit('ana_input: reference failed',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'file',item)) then
      if(.not.inp_a(filtrj)) call errquit('ana_input: file failed',0)
      if(.not.inp_i(ifrst)) ifrst=0
      if(.not.inp_i(ilast)) ilast=0
      goto 2
      endif
c
      if(inp_compare(.false.,'super',item)) then
      if(.not.inp_i(lsup)) lsup=1
      if(.not.inp_i(msupf)) msupf=0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trj',item)) then
      fmtsup='trj'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'arc',item)) then
      fmtsup='arc'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'bam',item)) then
      fmtsup='bam'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'amb',item)) then
      fmtsup='amb'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'mvm',item)) then
      fmtsup='mvm'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'xyz',item)) then
      fmtsup='xyz'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'eci',item)) then
      fmtsup='eci'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'frm',item)) then
      fmtsup='frm'
      if(.not.inp_a(item)) goto 1
      endif
      filsup=item
      goto 2
      endif
c
      if(inp_compare(.false.,'copy',item)) then
      if(.not.inp_i(lcop)) lcop=1
      if(.not.inp_i(mcopf)) mcopf=0
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'trj',item)) then
      fmtcop='trj'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'arc',item)) then
      fmtcop='arc'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'bam',item)) then
      fmtcop='bam'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'amb',item)) then
      fmtcop='amb'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'mvm',item)) then
      fmtcop='mvm'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'xyz',item)) then
      fmtcop='xyz'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'eci',item)) then
      fmtcop='eci'
      if(.not.inp_a(item)) goto 1
      endif
      if(inp_compare(.false.,'frm',item)) then
      fmtcop='frm'
      if(.not.inp_a(item)) goto 1
      endif
      filcop=item
      goto 2
      endif
c
      if(inp_compare(.false.,'select',item)) then
      rsel=0.0d0
      isel=0
      jsel=0
      atom=' '
      if(inp_i(isel)) then
      if(.not.inp_i(jsel)) jsel=isel
      else
      if(inp_f(rsel)) then
      if(inp_i(isel)) then
      if(.not.inp_i(jsel)) jsel=isel
      endif
      endif
      endif
      if(.not.inp_a(atom)) goto 4
    4 continue
      numsel=numsel+1
      write(69,6900) isel,jsel,rsel,atom(1:index(atom,' ')-1)
 6900 format(2i7,f12.6,a)
      if(inp_a(atom)) goto 4
      goto 2
      endif
c
      if(inp_compare(.false.,'label',item)) then
      rtag=0.0d0
      if(.not.inp_i(itag)) call errquit('ana_input: label failed',0)
      if(.not.inp_a(tag)) call errquit('ana_input: label failed',0)
      iatag=0
      jatag=0
      atom=' '
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      else
      if(inp_f(rtag)) then
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      endif
      endif
      endif
      if(.not.inp_a(atom)) goto 44
   44 continue
      numtag=numtag+1
      write(68,6800) itag,tag,iatag,jatag,rtag,atom(1:index(atom,' ')-1)
 6800 format(i3,a24,2i7,f12.6,a)
      if(inp_a(atom)) goto 44
      goto 2
      endif
c
      if(inp_compare(.false.,'time',item)) then
      if(.not.inp_f(timoff)) call errquit('ana_input: time failed',0)
      goto 2
      endif
c
      if(inp_compare(.false.,'frames',item)) then
      if(.not.inp_i(ifrfr)) call errquit('ana_input: maxfrm failed',0)
      if(.not.inp_i(ifrto)) then
      ifrto=ifrfr
      ifrfr=1
      endif
      goto 2
      endif
c
      if(inp_compare(.false.,'charges',item)) then
      if(.not.inp_a(filchg)) call errquit('ana_input: filchg failed',0)
      goto 2
      endif
      if(inp_compare(.false.,'plt',item)) then
      if(.not.inp_a(filplt)) call errquit('ana_input: filplt failed',0)
      goto 2
      endif
      if(inp_compare(.false.,'xyz',item)) then
      if(.not.inp_a(filxyz)) call errquit('ana_input: filxyz failed',0)
      goto 2
      endif
      if(inp_compare(.false.,'espgrid',item)) then
      if(.not.inp_i(lplt)) lplt=1
    5 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'spacing',item)) then
      if(.not.inp_f(spac)) call errquit('ana_input: spac failed',0)
      goto 5
      endif
      if(inp_compare(.false.,'range',item)) then
      if(.not.inp_f(rcut)) call errquit('ana_input: rcut failed',0)
      goto 5
      endif
      goto 3
      endif
c
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,1999) item
 1999 format('ana_input: unknown field ',a)
      call errquit('ana_input: error reading input file',0)
      endif
c
      fpref=filtrj
      if(index(filtrj,'.coo').gt.0) then
      fpref=filtrj(1:index(filtrj,'.coo')-1)
      elseif(index(filtrj,'.sco').gt.0) then
      fpref=filtrj(1:index(filtrj,'.sco')-1)
      elseif(index(filtrj,'.trj').gt.0) then
      fpref=filtrj(1:index(filtrj,'.trj')-1)
      else
      filtrj=fpref(1:index(fpref,' ')-1)//'.trj'
      endif
c
      if(filcop(1:1).eq.' ') then
      filcop=fpref(1:index(fpref,' ')-1)//'-copy.trj'
      endif
c
      if(filsup(1:1).eq.' ') then
      filsup=fpref(1:index(fpref,' ')-1)//'-super.trj'
      endif
c
      if(filchg(1:1).eq.' ') then
      filchg=fpref(1:index(fpref,' ')-1)//'.q '
      endif
c
c     default plt data to  ${file_prefix}.plt
c
      if(filplt(1:1).eq.' ') then
      filplt=fpref(1:index(fpref,' ')-1)//'.plt '
      endif
c
c     default xyz data to  ${file_prefix}.plt
c
      if(filxyz(1:1).eq.' ') then
      filxyz=fpref(1:index(fpref,' ')-1)//'.xyz '
      endif
c
c     put input in rtdb
c
      if(project(1:1).ne.' ') then
      if(.not.rtdb_cput(irtdb,'ana:project',1,project))
     + call errquit('ana_input: rtdb_put failed',0)
      endif
c
c     root mean square deviation
c
      if(.not.rtdb_put(irtdb,'ana:irmsd',mt_int,1,irmsd))
     + call errquit('ana_input: rtdb_put irmsdfailed',0)
c
c     reference input file
c
      if(.not.rtdb_cput(irtdb,'ana:filref',1,filref))
     + call errquit('ana_input: rtdb_put failed',0)
c
c     trajectory input file
c
      if(.not.rtdb_cput(irtdb,'ana:filtrj',1,filtrj))
     + call errquit('ana_input: rtdb_put failed',0)
c
c     copy trajectory
c
      if(.not.rtdb_cput(irtdb,'ana:filcop',1,filcop))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'ana:fmtcop',1,fmtcop))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:lcop',mt_int,1,lcop))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:mcopf',mt_int,1,mcopf))
     + call errquit('ana_input: rtdb_put failed',0)
c
c     superimpose trajectory
c
      if(.not.rtdb_cput(irtdb,'ana:filsup',1,filsup))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'ana:fmtsup',1,fmtsup))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:lsup',mt_int,1,lsup))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:msupf',mt_int,1,msupf))
     + call errquit('ana_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'ana:timoff',mt_dbl,1,timoff))
     + call errquit('ana_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'ana:ifrfr',mt_int,1,ifrfr))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:ifrto',mt_int,1,ifrto))
     + call errquit('ana_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'ana:ifrst',mt_int,1,ifrst))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:ilast',mt_int,1,ilast))
     + call errquit('ana_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'ana:lplt',mt_int,1,lplt))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:spac',mt_dbl,1,spac))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'ana:rcut',mt_dbl,1,rcut))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'ana:filchg',1,filchg))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'ana:filplt',1,filplt))
     + call errquit('ana_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'ana:filxyz',1,filxyz))
     + call errquit('ana_input: rtdb_put failed',0)
c
      if(numsel.gt.0) then
      close(unit=69)
      else
      close(unit=69,status='delete')
      endif
c
      if(numtag.gt.0) then
      close(unit=68)
      else
      close(unit=68,status='delete')
      endif
c
      return
c
 9999 continue
      call errquit('ana_input: error',0)
c
      return
      end
