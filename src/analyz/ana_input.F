      subroutine ana_input(irtdb)
c
c $Id: ana_input.F,v 1.14 2000-03-18 04:29:48 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
c
      character*255 item,prefix,atom,tag,atomi,atomj,atomk,atoml
      character*255 filcmd,filchg,filplt,filxyz,filref,filtrj,filcop
      integer lfncmd,numcmd,len,ivec,nspac,ibond,iangl,itors,iplan
      integer ifrfr,ifrto,ifrsk,ifrst,ilast,isel,jsel,itag,iatag,jatag
      real*8 timoff,rsel,rtag,rcut
      logical lref,lfil
c
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,'Analysis Input Module',40,.true.)
      endif
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,prefix))
     + call errquit('rtdb_cget failed on file_prefix',9999)
c
      len=0
      lfncmd=69
      numcmd=0
      lref=.false.
      lfil=.false.
c
      timoff=0.0d0
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
c     system command defines initial filenames and opens command file
c     ---------------------------------------------------------------
c
      if(inp_compare(.false.,'system',item)) then
      if(len.gt.0) call errquit('ana_input: duplicate system',0)
      if(.not.inp_a(prefix)) call errquit('ana_input: error system',0)
      len=index(prefix,' ')-1
      filcmd=prefix(1:len)//'.cmd '
      filchg=prefix(1:len)//'.q '
      filplt=prefix(1:len)//'.plt '
      filxyz=prefix(1:len)//'.xyz '
      filref=prefix(1:len)//'.rst '
      filtrj=prefix(1:len)//'.trj '
c
      open(unit=lfncmd,file=filcmd(1:index(filcmd,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=lfncmd)
      numcmd=0
c
      goto 2
      endif
c
c     define reference file
c     ---------------------
c
      if(inp_compare(.false.,'reference',item)) then
      if(.not.inp_a(filref)) call errquit('Error in reference file',0)
      write(lfncmd,1000) filref(1:index(filref,' ')-1)
 1000 format('refer  ',a)
      numcmd=numcmd+1
      lref=.true.
      goto 2
      endif
      if(.not.lref) then
      write(lfncmd,1000) filref(1:index(filref,' ')-1)
      numcmd=numcmd+1
      lref=.true.
      endif
c
c     define trajectory file(s)
c     -------------------------
c
      if(inp_compare(.false.,'file',item)) then
      if(.not.inp_a(filtrj)) call errquit('Error in trajectory file',0)
      if(index(filtrj,'.trj').le.0)
     + filtrj=filtrj(1:index(filtrj,' ')-1)//'.trj '
      if(.not.inp_i(ifrst)) ifrst=0
      if(.not.inp_i(ilast)) ilast=0
      write(lfncmd,1001) ifrst,ilast,filtrj(1:index(filtrj,' ')-1)
 1001 format('file   ',2i5,a)
      numcmd=numcmd+1
      lfil=.true.
      goto 2
      endif
c
c     define simulation time offset
c     -----------------------------
c
      if(inp_compare(.false.,'time',item)) then
      if(.not.inp_f(timoff)) call errquit('Error in time offset',0)
      write(lfncmd,1002) timoff
 1002 format('time   ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     define frequency processing frames
c     ----------------------------------
c
      if(inp_compare(.false.,'frames',item)) then
      if(.not.inp_i(ifrfr)) call errquit('Error in frames',0)
      if(.not.inp_i(ifrto)) then
      ifrto=ifrfr
      ifrfr=1
      ifrsk=1
      else
      if(.not.inp_i(ifrsk)) ifrsk=1
      endif
      write(lfncmd,1003) ifrfr,ifrto,ifrsk
 1003 format('frame  ',3i10)
      numcmd=numcmd+1
      goto 2
      endif
c
c     copy frames to file
c     -------------------
c
      if(inp_compare(.false.,'copy',item)) then
      if(.not.inp_a(filcop)) call errquit('Error in copy',0)
      if(index(filcop,'.').le.0)
     + filcop=filcop(1:index(filcop,' ')-1)//'.trj '
      if(filcop.eq.filtrj) call errquit('File name conflict',0)
      write(lfncmd,1004) filcop(1:index(filcop,' ')-1)
 1004 format('copy   ',a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     superimpose frames and write to file
c     ------------------------------------
c
      if(inp_compare(.false.,'super',item)) then
      if(.not.inp_a(filcop)) call errquit('Error in super',0)
      if(index(filcop,'.').le.0)
     + filcop=filcop(1:index(filcop,' ')-1)//'.trj '
      if(filcop.eq.filtrj) call errquit('File name conflict',0)
      write(lfncmd,1005) filcop(1:index(filcop,' ')-1)
 1005 format('super  ',a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     scan frames and write to file
c     -----------------------------
c
      if(inp_compare(.false.,'scan',item)) then
      if(.not.inp_a(filcop)) call errquit('Error in scan',0)
      write(lfncmd,1016) filcop(1:index(filcop,' ')-1)
 1016 format('scan   ',a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     select atoms
c     ------------
c
      if(inp_compare(.false.,'select',item)) then
      rsel=0.0d0
      isel=0
      jsel=0
      if(inp_i(isel)) then
      if(.not.inp_i(jsel)) jsel=isel
      else
      if(inp_f(rsel)) then
      if(inp_i(isel)) then
      if(.not.inp_i(jsel)) jsel=isel
      endif
      endif
      endif
      if(.not.inp_a(atom)) atom=' '
    4 continue
      write(lfncmd,1006) isel,jsel,rsel,atom(1:index(atom,' ')-1)
 1006 format('select ',2i7,f12.6,a)
      numcmd=numcmd+1
      if(atom(1:1).eq.' ') goto 2
      if(inp_a(atom)) goto 4
      goto 2
      endif
c
c     define bonds
c     ------------
c
      if(inp_compare(.false.,'bond',item)) then
      if(.not.inp_i(ibond)) call errquit('Error in bond',0)
      if(.not.inp_a(atomi)) call errquit('Error in bond',0)
      if(.not.inp_a(atomj)) call errquit('Error in bond',0)
      write(lfncmd,1012) ibond
 1012 format('bond   ',i5)
      write(lfncmd,1013) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1013) atomj(1:index(atomj,' ')-1)
 1013 format(a)
      goto 2
      endif

c
c     define angles
c     -------------
c
      if(inp_compare(.false.,'angle',item)) then
      if(.not.inp_i(iangl)) call errquit('Error in angle',0)
      if(.not.inp_a(atomi)) call errquit('Error in angle',0)
      if(.not.inp_a(atomj)) call errquit('Error in angle',0)
      if(.not.inp_a(atomk)) call errquit('Error in angle',0)
      write(lfncmd,1014) iangl
 1014 format('angle  ',i5)
      write(lfncmd,1013) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1013) atomj(1:index(atomj,' ')-1)
      write(lfncmd,1013) atomk(1:index(atomk,' ')-1)
      goto 2
      endif
c
c     define torsions
c     ---------------
c
      if(inp_compare(.false.,'torsion',item)) then
      if(.not.inp_i(itors)) call errquit('Error in torsion',0)
      if(.not.inp_a(atomi)) call errquit('Error in torsion',0)
      if(.not.inp_a(atomj)) call errquit('Error in torsion',0)
      if(.not.inp_a(atomk)) call errquit('Error in torsion',0)
      if(.not.inp_a(atoml)) call errquit('Error in torsion',0)
      write(lfncmd,1015) itors
 1015 format('torsio ',i5)
      write(lfncmd,1013) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1013) atomj(1:index(atomj,' ')-1)
      write(lfncmd,1013) atomk(1:index(atomk,' ')-1)
      write(lfncmd,1013) atoml(1:index(atoml,' ')-1)
      goto 2
      endif
c
c     define impropers
c     ----------------
c
      if(inp_compare(.false.,'improper',item)) then
      if(.not.inp_i(iplan)) call errquit('Error in improper',0)
      if(.not.inp_a(atomi)) call errquit('Error in improper',0)
      if(.not.inp_a(atomj)) call errquit('Error in improper',0)
      if(.not.inp_a(atomk)) call errquit('Error in improper',0)
      if(.not.inp_a(atoml)) call errquit('Error in improper',0)
      write(lfncmd,1017) iplan
 1017 format('improp ',i5)
      write(lfncmd,1013) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1013) atomj(1:index(atomj,' ')-1)
      write(lfncmd,1013) atomk(1:index(atomk,' ')-1)
      write(lfncmd,1013) atoml(1:index(atoml,' ')-1)
      goto 2
      endif
c
c     label atoms
c     -----------
c
      if(inp_compare(.false.,'label',item)) then
      rtag=0.0d0
      if(.not.inp_i(itag)) call errquit('Error in label',0)
      if(.not.inp_a(tag)) call errquit('Error in label',0)
      iatag=0
      jatag=0
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      else
      if(inp_f(rtag)) then
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      endif
      endif
      endif
      if(.not.inp_a(atom)) atom=' '
    5 continue
      write(lfncmd,1007) itag,tag,iatag,jatag,rtag,
     + atom(1:index(atom,' ')-1)
 1007 format('label  ',i3,a24,2i7,f12.6,a)
      numcmd=numcmd+1
      if(atom(1:1).eq.' ') goto 2
      if(inp_a(atom)) goto 5
      goto 2
      endif
c
c     root mean square deviation analysis
c     -----------------------------------
c
      if(inp_compare(.false.,'rmsdev',item)) then
      write(lfncmd,1008)
 1008 format('rmsdev ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     essential dynamics analysis
c     ---------------------------
c
      if(inp_compare(.false.,'essential',item)) then
      if(.not.inp_a(item)) item=' ' 
      write(lfncmd,1009) item
 1009 format('xaver  ',/,'covar  ',a)
      numcmd=numcmd+2
      goto 2
      endif
c
c     essential dynamics analysis
c     ---------------------------
c
      if(inp_compare(.false.,'project',item)) then
      if(.not.inp_i(ivec)) call errquit('No vector specified',0)
      if(.not.inp_a(item)) call errquit('No projection file',0)
      write(lfncmd,1010) ivec,item
 1010 format('projec ',i7,a)
      numcmd=numcmd+1
      goto 1
      endif
c
c     electrostatic potential
c     -----------------------     
c
      if(inp_compare(.false.,'esp',item)) then
      if(.not.inp_a(filchg)) call errquit('ana_input: filchg failed',0)
      if(index(filchg,'.xyzq').le.0)
     + filchg=filchg(1:index(filchg,' ')-1)//'.xyzq '
      if(.not.inp_i(nspac)) nspac=10
      if(.not.inp_f(rcut)) rcut=1.0d0
      if(.not.inp_a(filplt)) then
      filplt=filchg(1:index(filchg,'.xyzq'))
      endif
      if(index(filplt,'.plt').le.0)
     + filplt=filchg(1:index(filchg,'.xyzq')-1)//'.plt '
      write(lfncmd,1011) nspac,rcut,filchg(1:index(filchg,' ')-1),
     + filplt(1:index(filplt,' ')-1)
 1011 format('espplt ',i4,f12.6,a,' ',a)
      numcmd=numcmd+1
      goto 1
      endif
c
c      if(inp_compare(.false.,'charges',item)) then
c      if(.not.inp_a(filchg)) call errquit('ana_input: filchg failed',0)
c      goto 2
c      endif
c      if(inp_compare(.false.,'plt',item)) then
c      if(.not.inp_a(filplt)) call errquit('ana_input: filplt failed',0)
c      goto 2
c      endif
c      if(inp_compare(.false.,'xyz',item)) then
c      if(.not.inp_a(filxyz)) call errquit('ana_input: filxyz failed',0)
c      goto 2
c      endif
c      if(inp_compare(.false.,'espgrid',item)) then
c      if(.not.inp_i(lplt)) lplt=1
c    5 continue
c      if(.not.inp_a(item)) goto 1
c      if(inp_compare(.false.,'spacing',item)) then
c      if(.not.inp_f(spac)) call errquit('ana_input: spac failed',0)
c      goto 5
c      endif
c      if(inp_compare(.false.,'range',item)) then
c      if(.not.inp_f(rcut)) call errquit('ana_input: rcut failed',0)
c      goto 5
c      endif
c      goto 3
c      endif
c
c 
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,1999) item
 1999 format('Unknown field ',a)
      call errquit('Error reading input file',0)
      endif
c
c     close command file
c
      write(lfncmd,1099)
 1099 format('end    ')
      close(unit=lfncmd)
      if(numcmd.le.0) call errquit('Error: no analysis commands',0)
c
c     put command file name on rtdb
c     -----------------------------
c
      if(.not.rtdb_cput(irtdb,'ana:filcmd',1,filcmd))
     + call errquit('Error rtdb_put filcmd failed',0)
c
      return
c
 9999 call errquit('Error opening command file',0)
c
      return
      end
