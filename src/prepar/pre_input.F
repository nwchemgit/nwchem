      subroutine pre_input(irtdb)
c
c $Id: pre_input.F,v 1.57 2000-06-09 23:27:00 d3j191 Exp $
c
      implicit none
c
#include "rtdb.fh"
#include "inp.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer irtdb
      character*255 item,atom,type,atomi,atomj,atomk,atoml
c
      integer mcount
      parameter(mcount=10)
      integer mfract
      parameter(mfract=250)
c
      character*3 slvnam
      character*80 source
      character*10 slvmdl
      character*80 ffield,sysnam
      character*255 dir_s,dir_x,dir_u,dir_t,commnd,modify
      integer newtop,newseq,newrst,icount(10),multip,ncount
      integer nfract,ifract(mfract)
c
      integer len,numcmd,nummod,mgrid,mnoe,maxscf,num,mset(2),nset
      integer isgm,icyren,model,nxlnk,iconst,npmfr
      real*8 touch,xpnd(3),rgrid,edge,edges(3),dnoe(5),qscale,charge
      real*8 polar,value,forcon,forc2,box(3)
      character*11 cfix
      real*8 rfix
      character*1 altloc,chain
      integer ires,jres,kres,i,ignore,igroup,jgroup,kgroup,lgroup
      integer nxrep,nyrep,nzrep,itran,jtran,j,nkeepw
      real*8 rrep,rdist,xtran(3)
      real*8 value1,value2
      character*4 scount(mcount)
      real*8 xrmin,xrmax,yrmin,yrmax,zrmin,zrmax
c
      logical ltouch,lrdcoo,lslvnt,needwt,lwrcoo
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
      altloc='A'
      chain='A'
      newtop=0
      newseq=0
      newrst=0
      source=' '
      mgrid=24
      rgrid=0.2d0
      mnoe=0
      maxscf=0
      qscale=-1.0d0
      ncount=0
      nfract=0
      icyren=0
      iconst=0
      model=0
      nxlnk=0
      ignore=0
      nxrep=0
      nyrep=0
      nzrep=0
      rrep=0.25d0
      rdist=0.0d0
c
      touch=0.25d0
      xpnd(1)=0.00d0
      xpnd(2)=0.00d0
      xpnd(3)=0.00d0
      box(1)=0.0d0
      box(2)=0.0d0
      box(3)=0.0d0
c
      ltouch=.false.
      lrdcoo=.false.
      lwrcoo=.false.
      lslvnt=.false.
      needwt=.false.
c
      if(.not.rtdb_cget(irtdb,'file_prefix',1,commnd))
     + call errquit('rtdb_cget failed on file_prefix',9999)
      modify=commnd(1:index(commnd,' ')-1)//'.mod '
      open(unit=67,file=modify(1:index(modify,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=67)
      nummod=0
      commnd=commnd(1:index(commnd,' ')-1)//'.cmd '
      open(unit=68,file=commnd(1:index(commnd,' ')-1),
     + form='formatted',status='unknown',err=9999)
      rewind(unit=68)
      write(68,7801)
 7801 format('rdcoor ',/,'center',/,'boxsiz')
      lrdcoo=.true.
      numcmd=0
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
c     end
c     ---
c
      if(inp_compare(.false.,'end',item)) then
      goto 9
      endif
c
      lwrcoo=.false.
c
c     print & noprint
c     ---------------
c
      if(inp_compare(.false.,'print',item).or.
     + inp_compare(.false.,'noprint',item)) then
      call util_print_input(irtdb,'prepar')
      goto 1
      endif
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
c     ignore not found bonded parameters for quantum atoms
c
      if(inp_compare(.false.,'ignore',item)) then
      ignore=1
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
      num=0
   20 continue
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      if(inp_compare(.false.,'name',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: solvent',9999)
      slvnam=item(1:3)
      if(num.eq.2) goto 2
      num=1
      goto 20
      else
      if(inp_compare(.false.,'model',item)) then
      if(.not.inp_a(slvmdl)) call errquit('pre_input: solvent',9999)
      if(num.eq.1) goto 2
      num=2
      goto 20
      endif
      endif
      endif
c
c     alternate location
c     ------------------
c
      if(inp_compare(.false.,'altloc',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: altloc',9999)
      altloc=item(1:1)
      goto 2
      endif
c
c     chain
c     -----
c
      if(inp_compare(.false.,'chain',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: chain',9999)
      chain=item(1:1)
      goto 2
      endif
c
c     topology
c     --------
c
      if(inp_compare(.false.,'new_top',item)) then
      newtop=1
   21 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'new_seq',item)) then
      newseq=1
      goto 2
      endif
      goto 3
      endif
c
c     model
c     -----
c
      if(inp_compare(.false.,'model',item)) then
      if(.not.inp_i(model)) call errquit('pre_input: model',9999)
      goto 2
      endif
c
c     newrst
c     ------
c
      if(inp_compare(.false.,'new_rst',item)) then
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
c     rename commands
c     ---------------
c
      if(inp_compare(.false.,'sscyx',item)) then
      icyren=1
      goto 2
      endif
c
c     completion
c     ----------
c
      if(inp_compare(.false.,'construct',item)) then
      iconst=1
      goto 2
      endif
      if(inp_compare(.false.,'hbuild',item)) then
      iconst=1
      goto 2
      endif
c
c     divide molecules into fractions
c     -------------------------------
c
      if(inp_compare(.false.,'fraction',item)) then
    4 continue
      if(.not.inp_i(i)) goto 2
      nfract=nfract+1
      if(nfract.gt.mfract) call errquit('Too many fractions',nfract)
      ifract(nfract)=i
      goto 4
      endif
c
c     add counter ions to sequence file
c     ---------------------------------
c
      if(inp_compare(.false.,'counter',item)) then
      ncount=ncount+1
      if(ncount.gt.mcount) call errquit('Too many counters',ncount)
      if(.not.inp_i(icount(ncount)))
     + call errquit('pre_input: counter',9999)
      if(.not.inp_a(item)) call errquit('pre_input: counter',9999)
      scount(ncount)=item(1:4)
      goto 2
      endif
c
c     fix commands
c     ------------
c
      if(inp_compare(.false.,'fix',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: fix',9999)
      if(inp_compare(.false.,'atoms',item)) then
      cfix(8:11)='atm '
      else if(inp_compare(.false.,'segments',item)) then
      cfix(8:11)='sgm '
      else
      call errquit('pre_input: fix',9999)
      endif
      if(.not.inp_a(item)) call errquit('pre_input: fix',9999)
      if(inp_compare(.false.,'beyond',item)) then
      cfix(1:7)='fixout '
      else if(inp_compare(.false.,'within',item)) then
      cfix(1:7)='fixwin '
      else
      call errquit('pre_input: fix',9999)
      endif
      if(.not.inp_f(rfix)) call errquit('pre_input: fix',9999)
      if(.not.inp_a(item)) call errquit('pre_input: fix',9999)
      write(68,6841) cfix,rfix,item
 6841 format(a11,f12.6,1x,a)
      numcmd=numcmd+1
      goto 2
      endif
c
c     topology modification commands
c     ------------------------------
c
      if(inp_compare(.false.,'modify',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: modify',9999)
c
      if(inp_compare(.false.,'atom',item)) then
      if(.not.inp_a(atom)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
 3312 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3312
      endif
      if(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3312
      endif
      if(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3312
      endif
      if(inp_compare(.false.,'charge',item)) then
      if(.not.inp_f(charge)) call errquit('pre_input: charge',9999)
      write(67,6701) (atom(1:index(atom,' ')-1),mset(i),charge,i=1,nset)
 6701 format('atom ',a,t19,i1,' chg ',f12.6)
      nummod=nummod+1
      goto 3312
      endif
      if(inp_compare(.false.,'polar',item)) then
      if(.not.inp_f(polar)) call errquit('pre_input: polar',9999)
      write(67,6702) (atom(1:index(atom,' ')-1),mset(i),polar,i=1,nset)
 6702 format('atom ',a,t19,i1,' pol ',f12.9)
      nummod=nummod+1
      goto 3312
      endif
      if(inp_compare(.false.,'type',item)) then
      if(.not.inp_a(type)) call errquit('pre_input: type',9999)
      write(67,6703) (atom(1:index(atom,' ')-1),mset(i),
     + type(1:index(type,' ')-1),i=1,nset)
 6703 format('atom ',a,t19,i1,' typ ',a)
      nummod=nummod+1
      goto 3312
      endif
      if(inp_compare(.false.,'dummy',item)) then
      write(67,6704) (atom(1:index(atom,' ')-1),mset(i),i=1,nset)
 6704 format('atom ',a,t19,i1,' dum ')
      nummod=nummod+1
      goto 3312
      endif
      if(inp_compare(.false.,'self',item)) then
      write(67,6705) (atom(1:index(atom,' ')-1),mset(i),i=1,nset)
 6705 format('atom ',a,t19,i1,' ego ')
      nummod=nummod+1
      goto 3312
      endif
      if(inp_compare(.false.,'quantum',item)) then
      write(67,6706) (atom(1:index(atom,' ')-1),mset(i),i=1,nset)
 6706 format('atom ',a,t19,i1,' qtm ')
      nummod=nummod+1
      goto 3312
      else
      call errquit('pre_input: atom',9999)
      endif
c
      else if(inp_compare(.false.,'segment',item)) then
      if(.not.inp_i(isgm)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
 3313 continue
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3313
      endif
      if(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3313
      endif
      if(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3313
      endif
      if(inp_compare(.false.,'dummy',item)) then
      write(67,6707) (isgm,mset(i),i=1,nset)
 6707 format('segment ',i5,t19,i1,' dum ')
      nummod=nummod+1
      goto 3313
      endif
      if(inp_compare(.false.,'self',item)) then
      write(67,6708) (isgm,mset(i),i=1,nset)
 6708 format('segment ',i5,t19,i1,' ego ')
      nummod=nummod+1
      goto 3313
      endif
      if(inp_compare(.false.,'quantum',item)) then
      write(67,6709) (isgm,mset(i),i=1,nset)
 6709 format('segment ',i5,t19,i1,' qtm ')
      nummod=nummod+1
      goto 3313
      else
      call errquit('pre_input: segment',9999)
      endif
c
      else if(inp_compare(.false.,'bond',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
 3314 continue
      if(inp_f(value)) then
      if(.not.inp_f(forcon)) call errquit('pre_input: modify',9999)
      write(67,6710) (mset(i),value,forcon,atomi(1:index(atomi,' ')-1),
     + atomj(1:index(atomj,' ')-1),i=1,nset)
 6710 format('bond ',i1,1x,f12.6,e12.5,1x,a,1x,a)
      nummod=nummod+1
      endif
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3314
      endif
      if(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3314
      endif
      if(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3314
      else
      call errquit('pre_input: bond',9999)
      endif
c
      else if(inp_compare(.false.,'angle',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomk)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
 3315 continue
      if(inp_f(value)) then
      if(.not.inp_f(forcon)) call errquit('pre_input: modify',9999)
      write(67,6711) (mset(i),value,forcon,atomi(1:index(atomi,' ')-1),
     + atomj(1:index(atomj,' ')-1),atomk(1:index(atomk,' ')-1),i=1,nset)
 6711 format('angle ',i1,1x,f12.6,e12.5,1x,a,1x,a,1x,a)
      nummod=nummod+1
      endif
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3315
      elseif(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3315
      elseif(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3315
      else
      call errquit('pre_input: angle',9999)
      endif
c
      else if(inp_compare(.false.,'torsion',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomk)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atoml)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
      multip=0
 3316 continue
      if(inp_f(value)) then
      if(.not.inp_f(forcon)) call errquit('pre_input: modify',9999)
      write(67,6712) (mset(i),multip,value,forcon,
     + atomi(1:index(atomi,' ')-1),atomj(1:index(atomj,' ')-1),
     + atomk(1:index(atomk,' ')-1),atoml(1:index(atoml,' ')-1),i=1,nset)
 6712 format('torsion ',i1,1x,i1,1x,f12.6,e12.5,1x,a,1x,a,1x,a,1x,a)
      nummod=nummod+1
      endif
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3316
      endif
      if(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3316
      endif
      if(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3316
      endif
      if(inp_compare(.false.,'multiplicity',item)) then
      if(.not.inp_i(multip)) call errquit('pre_input: multip',9999)
      if(multip.lt.0.or.multip.gt.3)
     + call errquit('pre_input: multip',9999)
      goto 3316
      else
      call errquit('pre_input: torsion',9999)
      endif
c
      else if(inp_compare(.false.,'plane',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atomk)) call errquit('pre_input: modify',9999)
      if(.not.inp_a(atoml)) call errquit('pre_input: modify',9999)
      mset(1)=0
      nset=1
 3317 continue
      if(inp_f(value)) then
      if(.not.inp_f(forcon)) call errquit('pre_input: modify',9999)
      write(67,6713) (mset(i),value,forcon,
     + atomi(1:index(atomi,' ')-1),atomj(1:index(atomj,' ')-1),
     + atomk(1:index(atomk,' ')-1),atoml(1:index(atoml,' ')-1),i=1,nset)
 6713 format('plane ',i1,1x,f12.6,e12.5,1x,a,1x,a,1x,a,1x,a)
      nummod=nummod+1
      endif
      if(.not.inp_a(item)) goto 1
      if(inp_compare(.false.,'set',item)) then
      if(.not.inp_i(mset(1))) call errquit('pre_input: mset',9999)
      if(mset(1).lt.0.or.mset(1).gt.3)
     + call errquit('pre_input: mset',9999)
      goto 3317
      endif
      if(inp_compare(.false.,'initial',item)) then
      mset(1)=1
      mset(2)=2
      nset=2
      goto 3317
      endif
      if(inp_compare(.false.,'final',item)) then
      mset(1)=3
      nset=1
      goto 3317
      else
      call errquit('pre_input: plane',9999)
      endif
c
      elseif(inp_compare(.false.,'sequence',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: sequence',9999)
 3318 continue
      write(67,6716) item
 6716 format('sequence ',a)
      nummod=nummod+1
      if(inp_a(item)) goto 3318
      goto 1
c
      else
      call errquit('pre_input: modify',9999)
      endif
      goto 2
      endif
c
c     topology link commands
c     ----------------------
c
      if(inp_compare(.false.,'link',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: link atom 1',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: link atom 2',9999)
      i=index(atomi,':')
      if(i.le.1) call errquit('pre_input: link atom 1',9999)
      read(atomi(1:i-1),*) ires
      atomi=atomi(i+1:80)
      i=index(atomj,':')
      if(i.le.1) call errquit('pre_input: link atom 2',9999)
      read(atomj(1:i-1),*) jres
      atomj=atomj(i+1:80)
      if(atomi(1:1).eq.'_') atomi(1:1)=' '
      if(atomj(1:1).eq.'_') atomj(1:1)=' '
      write(67,6715) ires,atomi(1:4),jres,atomj(1:4)
 6715 format('link ',i5,a4,i5,a4)
      nummod=nummod+1
      nxlnk=nxlnk+1
      goto 1
      endif
c
c     update
c
      if(inp_compare(.false.,'update',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: update',9999)
c
      if(inp_compare(.false.,'lists',item)) then
      write(67,6714)
 6714 format('update')
      nummod=nummod+1
      else
      call errquit('pre_input: update',9999)
      endif
      goto 2
      endif
c
c     noe's
c     -----
c
      if(inp_compare(.false.,'noe',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: noe',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(1))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(2))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(3))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(4))) call errquit('pre_input: noe',9999)
      if(.not.inp_f(dnoe(5))) call errquit('pre_input: noe',9999)
      mnoe=mnoe+1
      write(68,6820) atomi(1:index(atomi,' ')-1),
     + atomj(1:index(atomj,' ')-1),dnoe
 6820 format('disres',1x,a,1x,a,/,3f12.6,2e12.5)
      numcmd=numcmd+1
      goto 2
      endif
c
c     select
c     ------
c
      if(inp_compare(.false.,'select',item)) then
      if(.not.inp_i(igroup)) call errquit('pre_input: select',9999)
      if(.not.inp_a(atomi)) call errquit('pre_input: noe',9999)
 6830 continue
      write(68,6831) igroup,atomi(1:index(atomi,' ')-1)
 6831 format('select',i5,1x,a)
      numcmd=numcmd+1
      if(.not.inp_a(atomi)) goto 1
      goto 6830
      endif
c
c     pmf
c     ---
c
      if(inp_compare(.false.,'pmf',item)) then
      if(.not.inp_a(type)) call errquit('pre_input: pmf',9999)
      if(.not.inp_i(igroup)) call errquit('pre_input: pmf',9999)
      if(type(1:8).eq.'distance') then
      if(.not.inp_i(jgroup)) call errquit('pre_input: pmf',9999)
      endif
      if(type(1:5).eq.'angle') then
      if(.not.inp_i(jgroup)) call errquit('pre_input: pmf',9999)
      if(.not.inp_i(kgroup)) call errquit('pre_input: pmf',9999)
      endif
      if(type(1:7).eq.'torsion') then
      if(.not.inp_i(kgroup)) call errquit('pre_input: pmf',9999)
      if(.not.inp_i(lgroup)) call errquit('pre_input: pmf',9999)
      endif
      if(type(1:5).ne.'align'.and.type(1:6).ne.'planar') then
      if(.not.inp_f(value1)) call errquit('pre_input: pmf',9999)
      if(.not.inp_f(value2)) call errquit('pre_input: pmf',9999)
      else
      value1=0.0d0
      value2=0.0d0
      endif
      if(.not.inp_f(forcon)) call errquit('pre_input: pmf',9999)
      if(.not.inp_f(forc2)) call errquit('pre_input: pmf',9999)
      if(type(1:5).eq.'align') then
      write(68,6844) igroup,forcon,forc2
 6844 format('pmf line',i5,15x,24x,2e12.5)
      elseif(type(1:6).eq.'planar') then
      write(68,6845) igroup,forcon,forc2
 6845 format('pmf plan',i5,15x,24x,2e12.5)
      elseif(type(1:8).eq.'distance') then
      write(68,6840) igroup,jgroup,value1,value2,forcon,forc2
 6840 format('pmf dist',2i5,10x,2f12.6,2e12.5)
      elseif(type(1:5).eq.'angle') then
      write(68,6842) igroup,jgroup,kgroup,value1,value2,forcon,forc2
 6842 format('pmf angl',3i5,5x,2f12.6,2e12.5)
      elseif(type(1:7).eq.'torsion') then
      write(68,6843) igroup,jgroup,kgroup,lgroup,value1,value2,
     + forcon,forc2
 6843 format('pmf tors',4i5,2f12.6,2e12.5)
      endif
      goto 2
      endif
c
c     repeat selections and pmfs
c     --------------------------
c
      if(inp_compare(.false.,'pmf_repeat',item)) then
      if(.not.inp_i(npmfr)) call errquit('pre_input: pmf_repeat',9999)
      write(68,6846) npmfr
 6846 format('rep_pmf',i5)
      numcmd=numcmd+1
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
      ltouch=.true.
      numcmd=numcmd+1
      goto 2
      endif
c
c     periodic
c     --------
c
      if(inp_compare(.false.,'periodic',item)) then
      write(68,6802)
 6802 format('period ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     vacuum
c     ------
c
      if(inp_compare(.false.,'vacuo',item)) then
      write(68,6803)
 6803 format('vacuum ')
      numcmd=numcmd+1
      goto 2
      endif
c
c     center
c     ------
c
      if(inp_compare(.false.,'center',item)) then
      write(68,6818)
 6818 format('center ',/,'boxsiz')
      numcmd=numcmd+1
      goto 2
      endif
c
c     orient
c     ------
c
      if(inp_compare(.false.,'orient',item)) then
      write(68,6819)
 6819 format('orient ',/,'boxsiz')
      numcmd=numcmd+1
      goto 2
      endif
c
c     translate
c     ---------
c
      if(inp_compare(.false.,'translate',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: translate',9999)
      if(inp_compare(.false.,'molecule',item)) then
      if(.not.inp_i(itran)) call errquit('pre_input: translate',9999)
      if(.not.inp_i(jtran)) jtran=itran
      if(.not.inp_f(xtran(1))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(2))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(3))) call errquit('pre_input: translate',9999)
      write(68,6852) itran,jtran,xtran
 6852 format('tramol ',2i5,3f12.6)
      numcmd=numcmd+1
      elseif(inp_compare(.false.,'segment',item)) then
      if(.not.inp_i(itran)) call errquit('pre_input: translate',9999)
      if(.not.inp_i(jtran)) jtran=itran
      if(.not.inp_f(xtran(1))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(2))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(3))) call errquit('pre_input: translate',9999)
      write(68,6853) itran,jtran,xtran
 6853 format('trasgm ',2i5,3f12.6)
      numcmd=numcmd+1
      elseif(inp_compare(.false.,'atom',item)) then
      if(.not.inp_i(itran)) call errquit('pre_input: translate',9999)
      if(.not.inp_i(jtran)) jtran=itran
      if(.not.inp_f(xtran(1))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(2))) call errquit('pre_input: translate',9999)
      if(.not.inp_f(xtran(3))) call errquit('pre_input: translate',9999)
      write(68,6854) itran,jtran,xtran
 6854 format('traatm ',2i5,3f12.6)
      numcmd=numcmd+1
      else
      call errquit('pre_input: translate',9999)
      endif
      goto 2
      endif
c
c     remove solvent
c     --------------
c
      if(inp_compare(.false.,'remove',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: translate',9999)
      if(inp_compare(.false.,'solvent',item)) then
      xrmin=0.0d0
      xrmax=0.0d0
      yrmin=0.0d0
      yrmax=0.0d0
      zrmin=0.0d0
      zrmax=0.0d0
      if(.not.inp_a(item)) call errquit('pre_input: translate',9999)
      j=-1
      if(inp_compare(.false.,'inside',item)) j=1
      if(inp_compare(.false.,'outside',item)) j=0
      if(j.eq.-1) call errquit('pre_input: translate',9999)
      do 6855 i=1,3
      if(inp_a(item)) then
      if(inp_compare(.false.,'x',item)) then
      if(.not.inp_f(xrmin)) call errquit('pre_input remove',0)
      if(.not.inp_f(xrmax)) call errquit('pre_input remove',0)
      endif
      if(inp_compare(.false.,'y',item)) then
      if(.not.inp_f(yrmin)) call errquit('pre_input remove',0)
      if(.not.inp_f(yrmax)) call errquit('pre_input remove',0)
      endif
      if(inp_compare(.false.,'z',item)) then
      if(.not.inp_f(zrmin)) call errquit('pre_input remove',0)
      if(.not.inp_f(zrmax)) call errquit('pre_input remove',0)
      endif
      else
      goto 6857
      endif
 6855 continue
 6857 continue
      write(68,6856) xrmin,xrmax,yrmin,yrmax,zrmin,zrmax,j
 6856 format('removw ',6f12.6,i1)
      numcmd=numcmd+1
      goto 2
      endif
      endif
c
c     expand with
c     -----------
c
      if(inp_compare(.false.,'expand',item)) then
      if(.not.inp_f(xpnd(1))) call errquit('pre_input: xpndw',9999)
      if(.not.inp_f(xpnd(2))) xpnd(2)=xpnd(1)
      if(.not.inp_f(xpnd(3))) xpnd(3)=xpnd(2)
      write(68,6849) xpnd
 6849 format('expand with ',3f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     expand-set with
c     ---------------
c
      if(inp_compare(.false.,'envelope',item)) then
      if(.not.inp_f(xpnd(1))) call errquit('pre_input: xpndw',9999)
      if(.not.inp_f(xpnd(2))) xpnd(2)=xpnd(1)
      if(.not.inp_f(xpnd(3))) xpnd(3)=xpnd(2)
      goto 2
      endif
c
c     solvate
c     -------
c
      if(inp_compare(.false.,'solvate',item)) then
      if(.not.ltouch) write(68,6804) touch
 6804 format('touch ',f12.6)
      ltouch=.true.
      if(.not.lrdcoo) write(68,7804)
 7804 format('rdcoor ')
      lrdcoo=.true.
      if(.not.lslvnt) write(68,7805)
 7805 format('solvnt ')
      lslvnt=.true.
      if(.not.inp_a(item)) then
      write(68,6805) xpnd
 6805 format('period ',/,'solvat cube ',/,
     + 'expand with ',3f12.6,/,'wr_rst')
      lwrcoo=.true.
      numcmd=numcmd+1
      goto 1
      else
      if(inp_compare(.false.,'cube',item)) then
      if(.not.inp_f(edge)) then
      write(68,6811)
 6811 format('period',/,'solvat cube')
      write(68,7807) xpnd
 7807 format('expand with ',3f12.6)
      else
      write(68,6812) edge
 6812 format('period',/,'solvat cube',f12.6)
      write(68,7807) xpnd
      endif
      elseif(inp_compare(.false.,'box',item)) then
      if(.not.inp_f(edges(1))) then
      write(68,6813)
 6813 format('period',/,'solvat box')
      write(68,7807) xpnd
      else
      if(.not.inp_f(edges(2))) edges(2)=edges(1)
      if(.not.inp_f(edges(3))) edges(3)=edges(2)
      write(68,6814) edges
 6814 format('period',/,'solvat box',3f12.6)
      write(68,7807) xpnd
      endif
      elseif(inp_compare(.false.,'troct',item)) then
      if(.not.inp_f(edge)) then
      write(68,6821)
 6821 format('period',/,'solvat troct')
      write(68,7807) xpnd
      else
      write(68,6822) edge
 6822 format('period',/,'solvat troct',f12.6)
      write(68,7807) xpnd
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
      needwt=.true.
      numcmd=numcmd+1
      goto 2
      endif
c
c     box
c     ---
c
      if(inp_compare(.false.,'box',item)) then
      if(.not.inp_f(box(1))) call errquit('pre_input: box',9999)
      if(.not.inp_f(box(2))) call errquit('pre_input: box',9999)
      if(.not.inp_f(box(3))) call errquit('pre_input: box',9999)
      write(68,8686) box
 8686 format('setbox',3f12.6)
      numcmd=numcmd+1
      goto 2
      endif
c
c     write
c     -----
c
      if(inp_compare(.false.,'write',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: write',9999)
      if(inp_compare(.false.,'solute',item)) then
      if(.not.inp_i(nkeepw)) nkeepw=-1
      if(.not.inp_a(item)) call errquit('pre_input: write',9999)
      if(inp_compare(.false.,'pdb',item)) then
      if(.not.inp_a(item)) then
      write(68,6806) nkeepw,' '
      else
      write(68,6806) nkeepw,item(1:index(item,' ')-1)
 6806 format('wr_pdb solute ',i10,a)
      endif
      elseif(inp_compare(.false.,'xyz',item)) then
      if(.not.inp_a(item)) then
      write(68,6824) nkeepw,' '
      else
      write(68,6824) nkeepw,item(1:index(item,' ')-1)
 6824 format('wr_xyz solute ',i10,a)
      endif
      elseif(index(item,'.pdb').gt.0) then
      write(68,6806) nkeepw,item(1:index(item,' ')-1)
      elseif(index(item,'.xyz').gt.0) then
      write(68,6824) nkeepw,item(1:index(item,' ')-1)
      else
      call errquit('pre_input: write solute',9999)
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
      lwrcoo=.true.
      endif
      elseif(inp_compare(.false.,'rst_old',item)) then
      if(.not.inp_a(item)) then
      write(68,6808) ' '
      else
      write(68,6848) item(1:index(item,' ')-1)
 6848 format('wr_old ',a)
      lwrcoo=.true.
      endif
      elseif(inp_compare(.false.,'xyz',item)) then
      if(.not.inp_a(item)) then
      write(68,6823) ' '
      else
      write(68,6823) item(1:index(item,' ')-1)
 6823 format('wr_xyz ',a)
      lwrcoo=.true.
      endif
      elseif(index(item,'.pdb').gt.0) then
      write(68,6807) item(1:index(item,' ')-1)
      elseif(index(item,'.rst').gt.0) then
      write(68,6808) item(1:index(item,' ')-1)
      elseif(index(item,'.rst_old').gt.0) then
      write(68,6848) item(1:index(item,'.rst_old')+3)
      elseif(index(item,'.xyz').gt.0) then
      write(68,6823) item(1:index(item,' ')-1)
      else
      call errquit('Problem with write file',0)
      endif
      needwt=.false.
      numcmd=numcmd+1
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
      if(numcmd.eq.0) rewind(68)
      write(68,6809) item(1:index(item,' ')-1)
 6809 format('rd_pdb ',a)
      lrdcoo=.true.
      elseif(inp_compare(.false.,'rst',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: read',9999)
      if(numcmd.eq.0) rewind(68)
      write(68,6810) item(1:index(item,' ')-1)
 6810 format('rd_rst ',a)
      lrdcoo=.true.
      elseif(inp_compare(.false.,'rst_old',item)) then
      if(.not.inp_a(item)) call errquit('pre_input: read',9999)
      if(numcmd.eq.0) rewind(68)
      write(68,6825) item(1:index(item,' ')-1)
 6825 format('rd_old ',a)
      lrdcoo=.true.
      elseif(index(item,'.pdb').gt.0) then
      if(numcmd.eq.0) rewind(68)
      write(68,6809) item(1:index(item,' ')-1)
      lrdcoo=.true.
      elseif(index(item,'.rst').gt.0) then
      if(numcmd.eq.0) rewind(68)
      write(68,6810) item(1:index(item,' ')-1)
      lrdcoo=.true.
      elseif(index(item,'.qrs').gt.0) then
      if(numcmd.eq.0) rewind(68)
      write(68,6810) item(1:index(item,' ')-1)
      lrdcoo=.true.
      elseif(index(item,'.rst_old').gt.0) then
      if(numcmd.eq.0) rewind(68)
      write(68,6825) item(1:index(item,'.rst_old')+3)
      lrdcoo=.true.
      endif
      endif
      numcmd=numcmd+1
      goto 2
      endif
c
c     align
c     -----
c
      if(inp_compare(.false.,'align',item)) then
      if(.not.inp_a(atomi)) call errquit('pre_input: align',9999)
      if(.not.inp_a(atomj)) call errquit('pre_input: align',9999)
      if(.not.inp_a(atomk)) call errquit('pre_input: align',9999)
      i=index(atomi,':')
      if(i.le.1) call errquit('pre_input: align atom 1',9999)
      read(atomi(1:i-1),*) ires
      atomi=atomi(i+1:80)
      i=index(atomj,':')
      if(i.le.1) call errquit('pre_input: align atom 2',9999)
      read(atomj(1:i-1),*) jres
      atomj=atomj(i+1:80)
      i=index(atomk,':')
      if(i.le.1) call errquit('pre_input: align atom 3',9999)
      read(atomk(1:i-1),*) kres
      atomk=atomk(i+1:80)
      if(atomi(1:1).eq.'_') atomi(1:1)=' '
      if(atomj(1:1).eq.'_') atomj(1:1)=' '
      if(atomk(1:1).eq.'_') atomk(1:1)=' '
      write(68,6850) ires,atomi(1:4),jres,atomj(1:4),kres,atomk(1:4)
 6850 format('align  ',i5,a4,i5,a4,i5,a4)
      numcmd=numcmd+1
      goto 1
      endif
c
c     repeat
c     ------
c
      if(inp_compare(.false.,'repeat',item)) then
      if(.not.inp_i(nxrep)) call errquit('pre_input: repeat',9999)
      if(.not.inp_i(nyrep)) call errquit('pre_input: repeat',9999)
      if(.not.inp_i(nzrep)) call errquit('pre_input: repeat',9999)
      if(.not.inp_f(rrep)) then
      endif
      if(.not.inp_f(rdist)) then
      endif
      if(nxrep.lt.0) call errquit('pre_input: repeat',9999)
      if(nyrep.lt.0) call errquit('pre_input: repeat',9999)
      if(nzrep.lt.-2) call errquit('pre_input: repeat',9999)
      write(68,6851) nxrep,nyrep,nzrep,rrep,rdist
 6851 format('repeat ',3i5,2f12.6)
      numcmd=numcmd+1
      goto 1
      endif
c
      write(6,1000) item
 1000 format('pre_input: unknown field ',a)
      call errquit('pre_input: error reading input file',0)
c
c     end
c     ---
c
    9 continue
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
      if(altloc(1:1).ne.'A') then
      if(.not.rtdb_cput(irtdb,'prep:altloc',1,altloc))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(chain(1:1).ne.'A') then
      if(.not.rtdb_cput(irtdb,'prep:chain',1,chain))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(model.ne.0) then
      if(.not.rtdb_put(irtdb,'prep:model',mt_int,1,model))
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
      if(.not.rtdb_put(irtdb,'prep:newtop',mt_int,1,newtop))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:newseq',mt_int,1,newseq))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:newrst',mt_int,1,newrst))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(.not.rtdb_put(irtdb,'prep:nfract',mt_int,1,nfract))
     + call errquit('pre_input: rtdb_put failed',0)
      if(nfract.gt.0) then
      if(.not.rtdb_put(irtdb,'prep:ifract',mt_int,nfract,ifract))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'prep:iconst',mt_int,1,iconst))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:icyren',mt_int,1,icyren))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:ncount',mt_int,1,ncount))
     + call errquit('pre_input: rtdb_put failed',0)
      if(ncount.gt.0) then
      if(.not.rtdb_put(irtdb,'prep:icount',mt_int,ncount,icount))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_cput(irtdb,'prep:scount',ncount,scount))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
      if(mgrid.gt.0) then
      if(.not.rtdb_put(irtdb,'prep:mgrid',mt_int,1,mgrid))
     + call errquit('pre_input: rtdb_put failed',0)
      if(.not.rtdb_put(irtdb,'prep:rgrid',mt_dbl,1,rgrid))
     + call errquit('pre_input: rtdb_put failed',0)
      endif
c
      if(.not.rtdb_put(irtdb,'prep:ignore',mt_int,1,ignore))
     + call errquit('pre_input: rtdb_put failed',0)
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
      if(.not.rtdb_put(irtdb,'prep:nxlnk',mt_int,1,nxlnk))
     + call errquit('pre_input: rtdb_put failed',0)
c
      if(nummod.gt.0) then
      close(unit=67)
      else
      close(unit=67,status='delete')
      endif
c
      if(numcmd.gt.0) then
      if(.not.lwrcoo) then
      write(68,7809)
 7809 format('wr_rst')
      endif
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
