      subroutine ana_input(irtdb)
c
c $Id$
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
      character*255 item,prefix,atom,tag,atomi,atomj,atomk,atoml,filtop
      character*255 filcmd,filchg,filplt,filxyz,filref,filtrj,filcop
      character*255 filpin,filhol,filpou,filavg
      integer lfncmd,numcmd,len,ivec,nspac,ibond,iangl,itors,iplan,icd
      integer ifrfr,ifrto,ifrsk,ifrst,ilast,isel,jsel,ksel,lsel
      integer iopt,ivctr,iorder
      integer itag,iatag,jatag,ipbc,nrdf,mcopf
      real*8 timoff,rsel,rtag,rcut,scale,cpk,stick,rval,rmin,rmax,rang
      real*8 arota,rrdf
      logical lref,lfil,lsol,lsuper,lselo
      integer iesppb,lent,idpdb,indx,icent,jcent,nwhb,nclean
      character*8 option
c
      if(ga_nodeid().eq.0) then
      call util_print_centered(6,'Analysis Input Module',40,.true.)
      endif
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,prefix))
     + call md_abort('rtdb_cget failed on file_prefix',9999)
c
      len=0
      lfncmd=69
      numcmd=0
      lref=.false.
      lfil=.false.
      nrdf=0
      rrdf=1.0d0
c
      timoff=0.0d0
      scale=1.0d0
      cpk=1.0d0
      stick=1.0d0
      icent=0
c
    1 continue
c
      if(.not.inp_read()) call md_abort('ana_input: inp_read failed',0)
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
      if(len.gt.0) call md_abort('ana_input: duplicate system',0)
      if(.not.inp_a(prefix)) call md_abort('ana_input: error system',0)
      len=index(prefix,' ')-1
      filcmd=prefix(1:len)//'.cmd '
      filchg=prefix(1:len)//'.q '
      filplt=prefix(1:len)//'.plt '
      filxyz=prefix(1:len)//'.xyz '
      filref=prefix(1:len)//'.rst '
      filtrj=prefix(1:len)//'.trj '
      lent=index(prefix,'_')-1
      if(lent.le.0) lent=len
      filtop=prefix(1:lent)//'.top '
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
      if(.not.inp_a(filref)) call md_abort('Error in reference file',0)
      lent=index(filref,'_')-1
      if(lent.le.0) lent=index(filref,'.')-1
      if(lent.le.0) call md_abort('Error in topology file',0)
      filtop=filref(1:lent)//'.top '
      write(lfncmd,1000) filref(1:index(filref,' ')-1)
      write(lfncmd,2000) filtop(1:index(filtop,' ')-1)
 1000 format('refer  ',a)
 2000 format('topol  ',a)
      numcmd=numcmd+1
      lref=.true.
      goto 2
      endif
      if(.not.lref) then
      write(lfncmd,1000) filref(1:index(filref,' ')-1)
      numcmd=numcmd+1
      write(lfncmd,2000) filtop(1:index(filtop,' ')-1)
      numcmd=numcmd+1
      lref=.true.
      endif
c
c     define trajectory file(s)
c     -------------------------
c
      if(inp_compare(.false.,'file',item)) then
      if(.not.inp_a(filtrj)) call md_abort('Error in trajectory file',0)
      if(index(filtrj,'.trj').le.0.and.index(filtrj,'.sco').le.0.and.
     + index(filtrj,'.coo').le.0)
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
      if(.not.inp_f(timoff)) call md_abort('Error in time offset',0)
      write(lfncmd,1002) timoff
 1002 format('time   ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     define centering fraction
c     -------------------------
c
      if(inp_compare(.false.,'center',item)) then
      if(.not.inp_i(icent)) call md_abort('Error in center',0)
      if(.not.inp_i(jcent)) jcent=icent
      write(lfncmd,1306) icent,jcent
 1306 format('center ',2i5)
      numcmd=numcmd+1
      goto 2
      endif
c
c     define frequency processing frames
c     ----------------------------------
c
      if(inp_compare(.false.,'frames',item)) then
      if(.not.inp_i(ifrfr)) call md_abort('Error in frames',0)
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
c     define scale (for povray output only)
c     -------------------------------------
c
      if(inp_compare(.false.,'scale',item)) then
      if(.not.inp_f(scale)) call md_abort('Error in scale',0)
      write(lfncmd,1303) scale
 1303 format('scale ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     define scale (for povray output only)
c     -------------------------------------
c
      if(inp_compare(.false.,'cpk',item)) then
      if(.not.inp_f(cpk)) cpk=1.0d0
      write(lfncmd,1304) cpk
 1304 format('cpk   ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
      if(inp_compare(.false.,'stick',item)) then
      if(.not.inp_f(stick)) stick=1.0d0
      write(lfncmd,1305) stick
 1305 format('stick ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     write pdb file of selected frame
c     --------------------------------
c
      if(inp_compare(.false.,'write',item)) then
      if(.not.inp_i(idpdb)) idpdb=0
      lsol=.false.
      lsuper=.false.
 3301 continue
      if(.not.inp_a(filcop)) call md_abort('Error in copy',0)
      if(filcop(1:6).eq.'solute') then
      lsol=.true.
      goto 3301
      endif
      if(filcop(1:5).eq.'super') then
      lsuper=.true.
      goto 3301
      endif
      option='        '
      if(lsuper) option='sup     '
      if(lsol) option=option(1:3)//' sol '
      if(index(filcop,'.').le.0)
     + filcop=filcop(1:index(filcop,' ')-1)//'.pdb '
      if(filcop.eq.filtrj) call md_abort('File name conflict',0)
      write(lfncmd,3004) idpdb,option,filcop(1:index(filcop,' ')-1)
 3004 format('wrtpdb',i8,a8,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     copy frames to file
c     -------------------
c
      if(inp_compare(.false.,'copy',item)) then
      if(.not.inp_i(mcopf)) mcopf=0
      if(.not.inp_a(filcop)) call md_abort('Error in copy',0)
      lsol=.false.
      if(filcop(1:6).eq.'solute') then
      if(.not.inp_a(filcop)) call md_abort('Error in super',0)
      lsol=.true.
      endif
      rang=0.0d0
      if(filcop(1:6).eq.'rotate') then
      if(.not.inp_f(rang)) call md_abort('Error in rotation time',0)
      if(.not.inp_a(filcop)) call md_abort('Error in super',0)
      endif
      if(index(filcop,'.').le.0)
     + filcop=filcop(1:index(filcop,' ')-1)//'.trj '
      if(filcop.eq.filtrj) call md_abort('File name conflict',0)
      if(lsol) then
      write(lfncmd,1004) mcopf,rang,filcop(1:index(filcop,' ')-1)
 1004 format('copys  ',i5,f12.6,a)
      else
      write(lfncmd,1024) mcopf,rang,filcop(1:index(filcop,' ')-1)
 1024 format('copy   ',i5,f12.6,a)
      endif
      numcmd=numcmd+1
      goto 2
      endif
c
c     clean up frames to file
c     -----------------------
c
      if(inp_compare(.false.,'clean',item)) then
      if(.not.inp_i(nclean)) nclean=0
      write(lfncmd,1034) nclean
 1034 format('clean  ',i6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     superimpose frames and write to file
c     ------------------------------------
c
      if(inp_compare(.false.,'super',item)) then
      if(.not.inp_a(filcop)) call md_abort('Error in super',0)
      lsol=.false.
      if(filcop(1:6).eq.'solute') then
      if(.not.inp_a(filcop)) call md_abort('Error in super',0)
      lsol=.true.
      endif
      rang=0.0d0
      if(filcop(1:6).eq.'rotate') then
      if(.not.inp_f(rang)) call md_abort('Error in rotation time',0)
      if(.not.inp_a(filcop)) call md_abort('Error in super',0)
      endif
      if(index(filcop,'.').le.0)
     + filcop=filcop(1:index(filcop,' ')-1)//'.trj '
      if(filcop.eq.filtrj) call md_abort('File name conflict',0)
      if(lsol) then
      write(lfncmd,1005) rang,filcop(1:index(filcop,' ')-1)
 1005 format('supers ',f12.6,a)
      else
      write(lfncmd,1025) rang,filcop(1:index(filcop,' ')-1)
 1025 format('super  ',f12.6,a)
      endif
      numcmd=numcmd+1
      goto 2
      endif
c
c     scan frames and write to file
c     -----------------------------
c
      if(inp_compare(.false.,'scan',item)) then
      if(.not.inp_a(item)) call md_abort('Error in scan',0)
      if(inp_compare(.false.,'super',item)) then
      if(.not.inp_a(filcop)) call md_abort('Error in scan',0)
      write(lfncmd,1016) filcop(1:index(filcop,' ')-1)
 1016 format('scans  ',a)
      else
      filcop=item
      write(lfncmd,1023) filcop(1:index(filcop,' ')-1)
 1023 format('scan   ',a)
      endif
      numcmd=numcmd+1
      goto 2
      endif
c
c     select atoms
c     ------------
c
      if(inp_compare(.false.,'select',item)) then
      lsuper=.false.
      lselo=.false.
      rsel=0.0d0
 1255 continue
      if(.not.inp_a(atom)) call md_abort('Error in select',0)
      if(atom(1:5).eq.'super') then
      lsuper=.true.
      lselo=.false.
      goto 1255
      endif
      if(atom(1:6).eq.'output') then
      lsuper=.false.
      lselo=.true.
      goto 1255
      endif
      if(atom(1:5).eq.'range') then
      if(.not.inp_f(rsel)) call md_abort('Error in select',1)
      goto 1255
      endif
      if(atom(1:7).eq.'solvent') then
      if(lsuper) call md_abort('Can not select super solvent',0)
      if(.not.inp_f(rsel)) call md_abort('No range specified',0)
      write(lfncmd,1206) rsel
 1206 format('wselct ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
      if(atom(1:4).eq.'save') then
      if(.not.inp_a(atom)) call md_abort('No file specified',0)
      if(lsuper) then
      write(lfncmd,1258) atom(1:index(atom,' ')-1)
 1258 format('savsls ',a)
      else
      write(lfncmd,1208) atom(1:index(atom,' ')-1)
 1208 format('savsel ',a)
      endif
      numcmd=numcmd+1
      goto 2
      endif
      if(atom(1:4).eq.'read') then
      if(.not.inp_a(atom)) call md_abort('No file specified',0)
      if(lsuper) then
      write(lfncmd,1259) atom(1:index(atom,' ')-1)
 1259 format('rdsels ',a)
      else
      write(lfncmd,1209) atom(1:index(atom,' ')-1)
 1209 format('rdsel  ',a)
      endif
      numcmd=numcmd+1
      goto 2
      endif
    4 continue
      if(lsuper) then
      write(lfncmd,1056) rsel,atom(1:index(atom,' ')-1)
 1056 format('selecs ',f12.6,a)
      elseif(lselo) then
      write(lfncmd,1057) rsel,atom(1:index(atom,' ')-1)
 1057 format('seleco ',f12.6,a)
      else
      write(lfncmd,1006) rsel,atom(1:index(atom,' ')-1)
 1006 format('select ',f12.6,a)
      endif
      numcmd=numcmd+1
      if(inp_a(atom)) goto 4
      goto 2
      endif
c
c     define groups of solute atoms
c     -----------------------------
c
      if(inp_compare(.false.,'define',item)) then
      if(.not.inp_i(isel)) call md_abort('ana_input: define',0)
      if(.not.inp_f(rsel)) rsel=0.0d0
      if(.not.inp_a(atom)) call md_abort('ana_input: define',0)
      if(atom(1:7).eq.'solvent') then
      if(.not.inp_f(rsel)) rsel=0.0d0
      if(.not.inp_a(atom)) atom=' '
      write(lfncmd,1207) isel,rsel,atom(1:index(atom,' ')-1)
 1207 format('wdefin ',i7,f12.6,a)
      numcmd=numcmd+1
      goto 2
      endif
    6 continue
      write(lfncmd,1020) isel,rsel,atom(1:index(atom,' ')-1)
 1020 format('define ',i7,f12.6,a)
      numcmd=numcmd+1
      if(inp_a(atom)) goto 6
      goto 2
      endif
c
c     do things with groups of atoms
c     ------------------------------
c
      if(inp_compare(.false.,'groups',item)) then
      if(.not.inp_i(isel)) call md_abort('No igroup specified',0)
      if(.not.inp_i(jsel)) jsel=-1
      if(.not.inp_i(ksel)) ksel=-1
      if(.not.inp_i(lsel)) lsel=-1
      if(.not.inp_a(item)) call md_abort('No group function',0)
      ipbc=0
      if(inp_compare(.false.,'periodic',item)) then
      if(.not.inp_i(ipbc)) ipbc=3
      if(.not.inp_a(item)) call md_abort('No group function',0)
      endif
      iopt=0
      if(inp_compare(.false.,'distance',item)) iopt=1
      if(inp_compare(.false.,'distances',item)) iopt=2
      if(inp_compare(.false.,'zcoordinate',item)) iopt=3
      if(inp_compare(.false.,'angle',item)) iopt=4
      if(inp_compare(.false.,'vectors',item)) iopt=5
      if(.not.inp_f(rsel)) rsel=0.0d0
      if(.not.inp_f(rval)) rval=rsel
      if(.not.inp_a(item)) item=' '
      if(iopt.eq.1.and.jsel.le.0)
     + call md_abort('Error in groups definition',4)
      if(iopt.eq.2.and.jsel.le.0)
     + call md_abort('Error in groups definition',5)
      if(iopt.eq.4.and.(jsel.le.0.or.ksel.le.0))
     + call md_abort('Error in groups definition',6)
c      if(iopt.eq.5.and.(jsel.le.0.or.ksel.le.0.or.lsel.le.0))
c     + call md_abort('Error in groups definition',7)
      write(lfncmd,1022) isel,jsel,ksel,lsel,iopt,ipbc,rsel,rval,
     + item(1:index(item,' ')-1)
 1022 format('groups ',4i7,i5,i3,2f12.6,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     do things with groups of atoms on the trajectory
c     ------------------------------------------------
c
      if(inp_compare(.false.,'histogram',item)) then
      if(.not.inp_i(isel)) call md_abort('No group specified',0)
      if(.not.inp_i(jsel)) jsel=0
      if(.not.inp_a(item)) call md_abort('No group function',0)
      ipbc=0
      if(inp_compare(.false.,'periodic',item)) then
      if(.not.inp_i(ipbc)) ipbc=3
      if(.not.inp_a(item)) call md_abort('No group function',0)
      endif
      ksel=0
      if(inp_compare(.false.,'zcoordinate',item)) ksel=1
      if(.not.inp_f(rsel)) rsel=0.0d0
      if(.not.inp_f(rval)) rval=rsel
      if(.not.inp_a(item)) item=' '
      write(lfncmd,1021) isel,jsel,ksel,ipbc,rsel,rval,
     + item(1:index(item,' ')-1)
 1021 format('histo  ',2i7,i5,i3,2f12.6,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     electrostatic potentials at specific positions
c     ----------------------------------------------
c
      if(inp_compare(.false.,'esp_points',item)) then
      if(.not.inp_a(filpin)) filpin='espp.inp '
      if(.not.inp_a(filhol)) filpin='espp.hol '
      if(.not.inp_a(filpou)) filpou='espp.dat '
      if(.not.inp_a(filavg)) filpou='espp.avg '
      write(lfncmd, 1055) filpin(1:index(filpin,' ')),
     + filhol(1:index(filhol,' ')),filpou(1:index(filpou,' ')),
     + filavg(1:index(filavg,' '))
 1055 format('esppnt ',a,a,a,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     radial distribution function
c     ----------------------------
c
      if(inp_compare(.false.,'rdf',item)) then
      if(.not.inp_i(nrdf))
     + call md_abort('No rdf number of bins specified',0)
      if(.not.inp_f(rrdf)) call md_abort('No rdf range specified',0)
      write(lfncmd,1032) nrdf,rrdf
 1032 format('rdf    ',i7,f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     group analysis
c     --------------
c
      if(inp_compare(.false.,'group',item)) then
      if(.not.inp_i(isel)) call md_abort('No group specified',0)
      if(.not.inp_a(item)) call md_abort('No group function',0)
      ipbc=0
      if(inp_compare(.false.,'periodic',item)) then
      if(.not.inp_i(ipbc)) ipbc=3
      if(.not.inp_a(item)) call md_abort('No group function',0)
      endif
      ksel=0
      if(inp_compare(.false.,'local',item)) ksel=1
      if(.not.inp_f(rsel)) rsel=0.0d0
      if(.not.inp_f(rval)) rval=rsel
      if(.not.inp_a(item)) item=' '
      if(ksel.eq.0) call md_abort('Unknown group function specified',0)
      write(lfncmd,1030) isel,ksel,ipbc,rsel,rval,
     + item(1:index(item,' ')-1)
 1030 format('group  ',i7,i5,i3,2f12.6,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     calculate order parameter
c     -------------------------
c
      if(inp_compare(.false.,'order',item)) then
      if(.not.inp_i(iorder)) call md_abort('No iorder specified',0)
      if(.not.inp_i(isel)) call md_abort('No igroup specified',0)
      if(.not.inp_i(jsel)) call md_abort('No jgroup specified',0)
      if(.not.inp_a(atomi)) call md_abort('No iatom specified',0)
      if(.not.inp_a(atomj)) call md_abort('No jatom specified',0)
      if(isel.eq.jsel) call md_abort('Order groups identical',isel)
      write(lfncmd,1028) iorder,isel,jsel 
 1028 format('order  ',i5,2i7)
      write(lfncmd,1029) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1029) atomj(1:index(atomj,' ')-1)
 1029 format(a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     define vectors
c     --------------
c
      if(inp_compare(.false.,'vector',item)) then
      if(.not.inp_i(ivctr)) call md_abort('Error in vector',0)
      if(.not.inp_a(atomi)) call md_abort('Error in vector',0)
      if(.not.inp_a(atomj)) call md_abort('Error in vector',0)
      write(lfncmd,1026) ivctr
 1026 format('vector ',i5)
      write(lfncmd,1027) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1027) atomj(1:index(atomj,' ')-1)
 1027 format(a)
      goto 2
      endif
c
c     define bonds
c     ------------
c
      if(inp_compare(.false.,'bond',item).or.
     + inp_compare(.false.,'distance',item)) then
      if(.not.inp_i(ibond)) call md_abort('Error in bond',0)
      if(.not.inp_a(atomi)) call md_abort('Error in bond',0)
      if(.not.inp_a(atomj)) call md_abort('Error in bond',0)
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
      if(.not.inp_i(iangl)) call md_abort('Error in angle',0)
      if(.not.inp_a(atomi)) call md_abort('Error in angle',0)
      if(.not.inp_a(atomj)) call md_abort('Error in angle',0)
      if(.not.inp_a(atomk)) call md_abort('Error in angle',0)
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
      if(.not.inp_i(itors)) call md_abort('Error in torsion',0)
      if(.not.inp_a(atomi)) call md_abort('Error in torsion',0)
      if(.not.inp_a(atomj)) call md_abort('Error in torsion',0)
      if(.not.inp_a(atomk)) call md_abort('Error in torsion',0)
      if(.not.inp_a(atoml)) call md_abort('Error in torsion',0)
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
      if(.not.inp_i(iplan)) call md_abort('Error in improper',0)
      if(.not.inp_a(atomi)) call md_abort('Error in improper',0)
      if(.not.inp_a(atomj)) call md_abort('Error in improper',0)
      if(.not.inp_a(atomk)) call md_abort('Error in improper',0)
      if(.not.inp_a(atoml)) call md_abort('Error in improper',0)
      write(lfncmd,1017) iplan
 1017 format('improp ',i5)
      write(lfncmd,1013) atomi(1:index(atomi,' ')-1)
      write(lfncmd,1013) atomj(1:index(atomj,' ')-1)
      write(lfncmd,1013) atomk(1:index(atomk,' ')-1)
      write(lfncmd,1013) atoml(1:index(atoml,' ')-1)
      goto 2
      endif
c
c     charge distribution
c     -------------------
c
      if(inp_compare(.false.,'charge_distribution',item)) then
      if(.not.inp_i(icd)) call md_abort('Error in improper',0)
      write(lfncmd,1018) icd
 1018 format('chgdst ',i5)
      goto 2
      endif
c
c     label atoms
c     -----------
c
      if(inp_compare(.false.,'label',item).or.
     + inp_compare(.false.,'render',item).or.
     + inp_compare(.false.,'color',item)) then
      rtag=0.0d0
      rval=1.0d0
      if(inp_compare(.false.,'render',item)) then
      itag=1
      elseif(inp_compare(.false.,'color',item)) then
      itag=2
      else
      if(.not.inp_i(itag)) call md_abort('Error in label',0)
      endif
      if(.not.inp_a(tag)) call md_abort('Error in label',0)
      iatag=-1
      jatag=-1
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      if(.not.inp_f(rtag)) rtag=0.0d0
      else
      if(inp_f(rval)) then
      if(inp_i(iatag)) then
      if(.not.inp_i(jatag)) jatag=iatag
      if(.not.inp_f(rtag)) rtag=0.0d0
      endif
      endif
      endif
      if(.not.inp_a(atom)) atom=' '
    5 continue
      write(lfncmd,1007) itag,tag,rval,iatag,jatag,rtag,
     + atom(1:index(atom,' ')-1)
 1007 format('label  ',i3,a24,f12.6,2i7,f12.6,a)
      numcmd=numcmd+1
      if(atom(1:1).eq.' ') goto 2
      if(inp_a(atom)) goto 5
      goto 2
      endif
c
c     define rotation
c     ---------------
c
      if(inp_compare(.false.,'rotate',item)) then
      if(.not.inp_a(item)) call md_abort('Error in rotate',0)
      if(inp_compare(.false.,'off',item)) then
      write(lfncmd,1588)
 1588 format('rotat0 ')
      numcmd=numcmd+1
      goto 2
      endif
      if(inp_compare(.false.,'x',item)) then
      if(.not.inp_f(arota)) call md_abort('Error in rotate',0)
      write(lfncmd,1587) arota
 1587 format('rotatx ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
      if(inp_compare(.false.,'y',item)) then
      if(.not.inp_f(arota)) call md_abort('Error in rotate',0)
      write(lfncmd,1586) arota
 1586 format('rotaty ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
      if(inp_compare(.false.,'z',item)) then
      if(.not.inp_f(arota)) call md_abort('Error in rotate',0)
      write(lfncmd,1585) arota
 1585 format('rotatz ',f12.6)
      numcmd=numcmd+1
      goto 2
      endif
      endif
c
c     ramachandran analysis
c     ---------------------
c
      if(inp_compare(.false.,'dielectric',item)) then
      write(lfncmd,1089)
 1089 format('diel   ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     ramachandran analysis
c     ---------------------
c
      if(inp_compare(.false.,'ramachandran',item)) then
      write(lfncmd,1088)
 1088 format('rama   ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     root mean square deviation analysis
c     -----------------------------------
c
      if(inp_compare(.false.,'rmsd',item)) then
      write(lfncmd,1008)
 1008 format('rmsdev ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     set hydrogen bond criteria
c     --------------------------
c
      if(inp_compare(.false.,'hbond',item)) then
      if(inp_a(item)) then
 6644 continue
      if(inp_compare(.false.,'distance',item)) then
      if(.not.inp_f(rmin)) call md_abort('Error in hbond',0)
      if(.not.inp_f(rmax)) then
      rmax=rmin
      rmin=0.0d0
      endif
      write(lfncmd,1066) rmin,rmax
 1066 format('hbdist ',2f12.6)
      if(inp_a(item)) goto 6644
      elseif(inp_compare(.false.,'donorangle',item).or.
     + inp_compare(.false.,'angle',item)) then
      if(.not.inp_f(rmin)) call md_abort('Error in hbond',0)
      if(.not.inp_f(rmax)) rmax=3.1416d0
      write(lfncmd,1067) rmin,rmax
 1067 format('hbdang ',2f12.6)
      if(inp_a(item)) goto 6644
      elseif(inp_compare(.false.,'acceptorangle',item)) then
      if(.not.inp_f(rmin)) call md_abort('Error in hbond',0)
      if(.not.inp_f(rmax)) then
      rmax=3.1416d0
      endif
      write(lfncmd,1068) rmin,rmax
 1068 format('hbaang ',2f12.6)
      if(inp_a(item)) goto 6644
      elseif(inp_compare(.false.,'solvent',item)) then
      if(.not.inp_i(nwhb)) nwhb=1
      write(lfncmd,1168) nwhb
 1168 format('numwhb ',i5)
      if(inp_a(item)) goto 6644
      endif
      endif
      write(lfncmd,1069)
 1069 format('hbond  ')
      goto 2
      endif
c
c     average structure
c     -----------------
c
      if(inp_compare(.false.,'average',item)) then
      if(.not.inp_a(item)) call md_abort('No average file',0)
      if(item(1:5).eq.'super') then
      if(.not.inp_a(item)) call md_abort('No average file',0)
      write(lfncmd,1070) item
 1070 format('saver ',a)
      numcmd=numcmd+1
      else
      write(lfncmd,1071) item
 1071 format('aver  ',a)
      numcmd=numcmd+1
      endif
      goto 2
      endif
c
c     essential dynamics analysis
c     ---------------------------
c
      if(inp_compare(.false.,'essential',item)) then
      if(.not.inp_a(item)) item=' ' 
      write(lfncmd,1009) item(1:index(item,' ')-1)
 1009 format('saver  ',/,'covar  ',a)
      numcmd=numcmd+2
      goto 2
      endif
c
c     essential dynamics analysis
c     ---------------------------
c
      if(inp_compare(.false.,'project',item)) then
      if(.not.inp_i(ivec)) call md_abort('No vector specified',0)
      if(.not.inp_a(item)) call md_abort('No projection file',0)
      write(lfncmd,1010) ivec,item(1:index(item,' ')-1)
 1010 format('projec ',i7,a)
      numcmd=numcmd+1
      goto 1
      endif
c
c     electrostatic potential
c     -----------------------     
c
      if(inp_compare(.false.,'esp',item)) then
      if(.not.inp_i(nspac)) nspac=10
      if(.not.inp_f(rcut)) rcut=1.0d0
      iesppb=0
      filplt=' '
      filchg=' '
      if(inp_a(item)) then
      if(inp_compare(.false.,'periodic',item)) then
      if(.not.inp_i(iesppb)) iesppb=3
      if(inp_a(filchg)) then
      if(.not.inp_a(filplt)) filplt=' '
      endif
      else
      filchg=item
      if(.not.inp_a(filplt)) filplt=' '
      endif
      endif
      if(filchg(1:1).ne.' ') then
      if(index(filchg,'.xyzq').le.0)
     + filchg=filchg(1:index(filchg,' ')-1)//'.xyzq '
      if(filplt(1:1).eq.' ')
     + filplt=filchg(1:index(filchg,'.xyzq')-1)//'.plt '
      endif
      if(filplt(1:1).ne.' ') then
      if(index(filplt,'.plt').le.0)
     + filplt=filchg(1:index(filchg,'.xyzq')-1)//'.plt '
      endif
      write(lfncmd,1011) nspac,iesppb,rcut,
     + filchg(1:index(filchg,' ')-1),
     + filplt(1:index(filplt,' ')-1)
 1011 format('espplt ',i4,i4,f12.6,a,' ',a)
      numcmd=numcmd+1
      goto 1
      endif
c
      if(inp_compare(.false.,'index',item)) then
      if(.not.inp_i(indx)) indx=1
      write(lfncmd,1019) indx
 1019 format('index ',i10)
      numcmd=numcmd+1
      goto 2
      endif
c
c     make report
c     -----------
c
      if(inp_compare(.false.,'report',item)) then
      if(.not.inp_a(item)) call md_abort('No report file',0)
      if(.not.inp_a(tag)) tag=' '
      option='        '
      if(inp_compare(.false.,'local',tag)) option='local   '
      write(lfncmd,1031) tag,item(1:index(item,' ')-1)
 1031 format('report ',a8,a)
      numcmd=numcmd+1
      goto 2
      endif
c 
      if(.not.inp_compare(.false.,'end',item)) then
      write(6,1999) item
 1999 format('Unknown field ',a)
      call md_abort('Error reading input file',0)
      endif
c
c     close command file
c
      write(lfncmd,1099)
 1099 format('end    ')
      close(unit=lfncmd)
      if(numcmd.le.0) call md_abort('Error: no analysis commands',0)
c
c     put command file name on rtdb
c     -----------------------------
c
      if(.not.rtdb_cput(irtdb,'ana:filcmd',1,filcmd))
     + call md_abort('Error rtdb_put filcmd failed',0)
c
      return
c
 9999 call md_abort('Error opening command file',0)
c
      return
      end
