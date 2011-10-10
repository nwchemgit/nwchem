      subroutine argos_ti()
c
      implicit none
c
#include "argos.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "msgids.fh"
c
      logical frequency,argos_prop_mcti_step,argos_prop_rdmri
      logical argos_space_rdmri
      integer argos_prop_dfr
      external frequency,argos_prop_mcti_step,argos_prop_rdmri
      external argos_space_rdmri,argos_prop_dfr
c
      logical done
      integer npp,i,irunp,jequi,jdacq,ndec,ndum
      real*8 rdum
      character*256 filrun
c
      if(nserie.eq.0) then
      if(mropt.ge.2) then
      if(me.eq.0) then
      open(unit=lfnmri,file=filmri(1:index(filmri,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmri) npp
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      endif
      krun=0
      endif
      if(me.eq.0) then
      open(unit=lfnmro,file=filmro(1:index(filmro,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      write(lfnmro) np,mrun,mequi,mdacq
      endif
      else
c
      if(me.eq.0) then
      open(unit=lfnmro,file=filmro(1:index(filmro,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmro) npp,mrun,mequi,mdacq
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      endif
      krun=0
      do 11 irun=1,mrun
      if(me.eq.0) then
      read(lfnmro,end=12,err=12) irunp,kequi,kdacq
      endif
      if(.not.argos_prop_rdmri(lfnmro,ndec,mropt)) goto 12
      if(.not.argos_space_rdmri(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs))) goto 12
      krun=irun
   11 continue
   12 continue
c
      if(me.eq.0) then
      open(unit=lfnmri,file=filmri(1:index(filmri,' ')-1),
     + form='unformatted',status='unknown',err=9999)
      read(lfnmri) npp
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      rewind(lfnmro)
      read(lfnmro) npp,mrun,mequi,mdacq
      if(npp.ne.np) call md_abort('Number of nodes changed',npp)
      rewind(lfngib)
      endif
c
      do 13 irun=1,mrun
      if(me.eq.0) then
      read(lfnmri) irunp,kequi,kdacq
      endif
      if(.not.argos_prop_rdmri(lfnmri,ndec,mropt))
     + call md_abort('Error in mri file',0)
      if(.not.argos_space_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs)))
     + call md_abort('Error in mri file',0)
      read(lfnmro) irunp,kequi,kdacq
      if(.not.argos_prop_rdmri(lfnmro,ndec,mropt))
     + call md_abort('Error in mro file',0)
      if(.not.argos_space_rdmri(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs)))
     + call md_abort('Error in mro file',0)
      if(me.eq.0) then
      read(lfngib,2000) ndum
 2000 format(7x,i7)
      read(lfngib,2001) (rdum,i=1,24)
 2001 format(4e20.12)
      read(lfngib,2001) (rdum,i=1,ndum)
      read(lfngib,2001) (rdum,i=1,ndum)
      read(lfngib,2002) ndum
 2002 format(i10)
      read(lfngib,2003) rdum
      read(lfngib,2003) rdum
 2003 format(e20.12)
      if(ndec.gt.0) read(lfngib,2001) (rdum,i=1,nsa)
      endif
   13 continue
c
      endif
c
      call ga_brdcst(mrg_d30,krun,ma_sizeof(mt_int,1,mt_byte),0)
c
      do 1 irun=krun+1,mrun
c
      if(npg.gt.1) then
      if(irun.ne.meg+1) goto 1
      endif
c
      if(me.eq.0) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filtrj(1:index(filtrj,'.trj')-1),'-',
     + irun,'.trj '
      open(unit=lfntrj,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      call argos_cafe_trjhdr(lfntrj)
      endif
      if(nfprop.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filprp(1:index(filprp,'.prp')-1),'-',
     + irun,'.prp '
      open(unit=lfnprp,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      endif
      if(nfrdf.gt.0) then
      write(filrun,'(a,a,i5.5,a)') filrdf(1:index(filprp,'.prp')-1),'-',
     + irun,'.rdf '
      open(unit=lfnrdf,file=filrun(1:index(filrun,' ')-1),
     + form='formatted',status='unknown')
      endif
      endif
c
      if(irun.eq.1.and.iand(ivopt,2).eq.2) nfgaus=ivreas
      if(irun.eq.maxlam.and.iand(ivopt,4).eq.4) nfgaus=ivreas
c
      lfirst=.true.
c
c     initialize parameters
c
      call argos_cafe_lambda(lamtyp,irun,maxlam,elam,lfnout,lfnpmf,
     + rlambd,dlambd,filnam)
c
c     property initialization
c
      call argos_prop_init()
c
      if(mropt.ge.2) then
      if(me.eq.0) then
      read(lfnmri,end=399,err=399) irunp,kequi,kdacq
      if(irunp.ne.irun) call md_abort('Number of run changed',irunp)
      if(.not.argos_prop_rdmri(lfnmri,ndec,mropt))
     + call md_abort('Error in mri',0)
      if(kequi+kdacq.lt.mequi) then
      kequi=kequi+kdacq
      kdacq=0
      elseif(kequi.lt.mequi) then
      if(kdacq+kequi-mequi.gt.0) then
      kdacq=argos_prop_dfr(kdacq+kequi-mequi)
      else
      kdacq=0
      endif
      kequi=mequi
      endif
      goto 398
  399 continue
      kequi=0
      kdacq=0
  398 continue
      endif
      call ga_brdcst(mrg_d31,kequi,ma_sizeof(mt_int,1,mt_byte),0)
      call ga_brdcst(mrg_d32,kdacq,ma_sizeof(mt_int,1,mt_byte),0)
c
c     kequi will be 0 if no records could be read from lfnmri, and
c     the coordinates and velocities will be used from the previous
c     lambda run.
c     if this is the first run there must be something wrong
c
      if(kequi.eq.0.and.kdacq.eq.0) then
      if(irun.eq.1) call md_abort('No records found on mri',me)
      else
      if(.not.argos_space_rdmri(lfnmri,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs)))
     + call md_abort('Error reading mri in argos_space_rdmri',0)
      endif
      endif
c
      lxw=.false.
      lvw=.false.
      lxs=.false.
      lvs=.false.
      lesp=.false.
c
      if(mropt.eq.2) then
      kequi=0
      kdacq=0
      call argos_prop_init
      endif
c
c     equilibration
c
      lequi=.true.
      jequi=kequi
      lprpmf=.false.
      iprpmf=0
      do 2 iequi=kequi+1,mequi
c
      mdstep=mdstep+1
      lpmfc=npmf.gt.1.and.iequi.gt.npmf
      call timer_start(201)
      call argos_newton()
      call timer_stop(201)
      jequi=iequi
      stime=stime+tstep
    2 continue
      lpmfc=.true.
c
c     data gathering
c
      mdstep=kdacq
      if(kdacq.eq.0) stime=zero
c
      ndec=0
      if(npgdec.gt.1) call argos_cafe_dera_init()
c
      lequi=.false.
      jdacq=kdacq
      lprpmf=.true.
      iprpmf=-1
      do 3 idacq=kdacq+1,mdacq
c
      mdstep=mdstep+1
c
      lxw=frequency(mdstep,nfcoor)
      lvw=frequency(mdstep,nfvelo)
      lfw=frequency(mdstep,nfforc)
      lxs=frequency(mdstep,nfscoo)
      lvs=frequency(mdstep,nfsvel)
      lfs=frequency(mdstep,nfsfor)
c
      call argos_timer_init()
c
      call timer_start(201)
c
      call argos_newton()
c
      if(npgdec.gt.1) ndec=ndec+1
c
      call argos_cafe_mcti_kin(int_mb(i_is+(lsatt-1)*msa),
     + int_mb(i_is+(lsgan-1)*msa),dbl_mb(i_vs),nsaloc)
c
      done=argos_prop_mcti_step(idacq,ldacq)
c
      call argos_prop_proper(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme,
     + epot,epotw,epotsw,epots,volume,dwr,dbl_mb(i_dsr),ekin,etot,
     + npolit,dbl_mb(i_gsm),dbl_mb(i_esa),box,dbl_mb(i_xsmp))
      call argos_prop_step(mdstep,stime,eww,dbl_mb(i_esw),
     + dbl_mb(i_ess),dbl_mb(i_fss),dbl_mb(i_esk),epme)
c
      if(lfw.or.lfs) then
      call argos_space_gaputf(me,dbl_mb(i_fw),nwmloc,dbl_mb(i_fs),
     + nsaloc)
      endif
c
      write(projct,4000) nserie,irun,mrun,mequi,idacq,tmpext,
     + filnam(1:32)
 4000 format(i2,' ti ',i5,'/',i5,' :',i7,'+ ',i7,' @ ',f7.2,' K ',a)
      if(frequency(mdstep,nfrest))
     + call argos_wrtrst(lfnrst,rfile,.true.)
      if(frequency(mdstep,nftime)) call argos_wrtime
c
      call timer_stop(201)
c
      jdacq=idacq
      if(idacq.gt.ldacq.and.done) goto 4
c
      stime=stime+tstep
c
    3 continue
    4 continue
c
      if(me.eq.0) then
      write(lfnmro) irun,jequi,jdacq
      endif
      call argos_prop_wrtmro(lfnmro,ndec)
      if(me.eq.0) then
      call argos_space_wrtmro(lfnmro,stime,pres,temp,tempw,temps,
     + int_mb(i_iw),dbl_mb(i_xw),dbl_mb(i_vw),dbl_mb(i_xwcr),
     + int_mb(i_is),dbl_mb(i_xs),dbl_mb(i_vs),projct)
      endif
c
      call argos_prop_mcti_run(rlambd,dlambd,ndec)
c
      if(ndec.gt.2) call argos_cafe_print_dera(lfnout,ndec)
c
      if(me.eq.0) then
      if(nfcoor.gt.0.or.nfscoo.gt.0.or.nfvelo.gt.0.or.nfsvel.gt.0) then
      close(unit=lfntrj,status='keep')
      endif
      if(nfprop.gt.0) then
      close(unit=lfnprp,status='keep')
      endif
      endif
c
    1 continue
c
      call argos_prop_mcti(npgdec,filnam)
c
      if(me.eq.0) then
      if(nserie.eq.0) then
      if(mropt.ge.2) then
      close(unit=lfnmri,status='keep')
      endif
      close(unit=lfnmro,status='keep')
      else
      close(unit=lfnmro,status='keep')
      endif
      endif
c
      return
c
 9999 continue
      call md_abort('Failed to open file mro',0)
      return
      end
c $Id$
