      subroutine argos_cafe_start(irtdbi,lfnout,lfntop,filtop,
     + ndistr,
     + npmf,nipmf,nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai,
     + mdalg,npbt,nbxt,rcs,rcl,bx,
     + jpme,jorder,jgx,jgy,jgz,nodpme,spmet,step,tols,mshw,mshs,nosh,
     + iipolt,iitmps,iiprss,iipopt,rtmpx,rprsx,rtmpw,rtmps,rpres,
     + sclq,fpmf,iislow,tempi,tempwi,tempsi,compr,ntyp,idset,ipset1,
     + ipset2,issscl,delta,nfanal,lpbc,npgi,fldi,fvec,ffrq,npenrg,ictrl,
     + nbiasi,mropti,incl,ltwn,nseqi,i_lseqi,nfhopi,rhopi,thopi,ndumsi,
     + ipbtpi,lfnhopi,iradgi,nbgeti,npreci)
c $Id$
      implicit none
c
#include "argos_cafe_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer lfnout,lfntop,irtdbi,lfnhopi
      character*(*) filtop
      integer ndistr,nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai
      integer jpme,jorder,jgx,jgy,jgz,mshw,mshs,nosh,npmf,nipmf,npgi
      real*8 spmet,step,tols,bx(3),fldi,fvec(3),ffrq,sclq,rhopi,thopi
      integer mdalg,npbt,nbxt,iipolt,nfanal,ictrl,incl
      integer nfhopi,ipbtpi,iradgi,nbgeti,npreci
      real*8 rcs,rcl,compr,delta,fpmf
      integer iitmps,iiprss,iislow,iipopt,nseqi,i_lseqi
      integer nodpme,ntyp,idset,ipset1,ipset2,ndumsi
      integer issscl,npenrg,nbiasi,mropti
      real*8 rtmpx,rprsx,rtmpw,rtmps,rpres,tempi,tempwi,tempsi
      logical lpbc,ltwn
c
      integer i,itemp(2)
      character*3 string
c
      call ga_sync()
c
      if(.not.ma_verify_allocator_stuff()) then
      call md_abort('ERROR IN MEMORY USE AT START ARGOS_CAFE_START',0)
      endif
c
      me=ga_nodeid()
      np=ga_nnodes()
c
c     dimensions initially set to zero indicating non-allocated
c
      irtdb=irtdbi
      lpress=lpbc
c
      itscal=iitmps
      ipscal=iiprss
      ipopt=iipopt
      tmpext=rtmpx
      prsext=rprsx
      tmwrlx=rtmpw
      tmsrlx=rtmps
      prsrlx=rpres
      scaleq=sclq
      facpmf=fpmf
      nbget=nbgeti
      nprec=npreci
c
      ntype=ntyp
      mdalgo=mdalg
      includ=incl
      ipbtyp=ipbtpi
c
      lfnhop=lfnhopi
      nhops=0
c
      lfree=ntype.eq.3
c
      nbs=0
      mscr=0
      lscr=.false.
      llst=.false.
      lpair=.true.
      llist=.false.
      lpmf=npmf.gt.0
      ndxp=0
      nlda=0
      maxl=0
      lanal=nfanal.gt.0
      npener=npenrg
      icntrl=ictrl
      iradgy=iradgi
c
      npgdec=npgi
c
      mwm=mwmi
      mwa=mwai
      msf=msfi
      msm=msmi
      msa=msai
      nwm=nwmi
      nwa=nwai
      nsf=nsfi
      nsm=nsmi
      nsa=nsai
      mscr=max(mwm+1,msa+1)
      nwmtot=nwm
      nsatot=nsa
c
      rshrt=rcs
      rlong=rcl
      rshrt2=rshrt*rshrt
      rlong2=rlong*rlong
      ltwin=ltwn
      lssscl=issscl.ne.0
c
      nbxtyp=nbxt
      npbtyp=npbt
c
      do 1 i=1,3
      box(i)=bx(i)
      boxh(i)=half*bx(i)
    1 continue
c
      ipolt=iipolt
      islow=iislow
c
      lstype=1
c
      ngc=1
      ngl=1
      nfrdf=99999
      ifstep=1
      ngrww=0
      ngrsw=0
      ngrss=0
      ireact=0
      iset=idset
      issscl=0
      nrwrec=0
      isolvo=0
c
      ndums=ndumsi
c
      lpww=1
      lpsw=1
      lpss=1
c
      npww=1
      npsw=1
      npss=1
      if(ltwin) then
      npww=2
      npsw=2
      npss=2
      endif
c
      mgc=1
      mgl=1
      mgr=1
c
      rffww=0.0d0
      rffsw=0.0d0
      rffss=0.0d0
c
      tstep=step
      tstepi=one/tstep
      tolsha=tols
      mshitw=mshw
      mshits=mshs
      noshak=nosh
c
      temp=tempi
      tempw=tempwi
      temps=tempsi
c
      q14fac=0.833333d0
      facpsc=compr*tstep/prsrlx
c
      field=fldi
      fvect(1)=fvec(1)
      fvect(2)=fvec(2)
      fvect(3)=fvec(3)
      ffreq=ffrq
c
      pi=four*atan(one)
      twopi=two*pi
c
      wbox=zero
c
      shift0(1)=zero
      shift0(2)=zero
      shift0(3)=delta
      shift0(4)=delta
      shift0(5)=zero
      shift0(6)=delta
      shift1(1)=delta
      shift1(2)=delta
      shift1(3)=zero
      shift1(4)=-delta
      shift1(5)=delta
      shift1(6)=zero
c
      numpmf=0
c
      rhop=rhopi
      rhop2=rhop*rhop
      thop=thopi
      nfhop=nfhopi
      nseq=nseqi
      mseq=nseqi
      i_lseq=i_lseqi
c
      lqhop=nfhop.ne.0
c
      if(lqhop.and.lfree)
     + call md_abort('Proton hopping thermodynamics not allowed',0)
c
      ithint=ntype.eq.3
      ipert2=ithint.or.(iset.eq.1.and.(ipset1.eq.2.or.ipset2.eq.2))
      ipert3=ithint.or.(iset.eq.1.and.(ipset1.eq.2.or.ipset2.eq.2))
      do 2 i=1,24
      ith(i)=.false.
      ip2(i)=.false.
      ip3(i)=.false.
    2 continue
c
      if(.not.ma_push_get(mt_byte,16*nsatot,'snam',l_snam,i_snam))
     + call md_abort('Failed to allocate snam',me)
c
      if(lqhop) then
      if(.not.ma_push_get(mt_int,mseq,'mprot',l_mprot,i_mprot))
     + call md_abort('Failed to allocate mprot',me)
      endif
c
      if(npgdec.gt.1) then
      if(.not.ma_push_get(mt_dbl,6*nsatot,'sti',l_sti,i_sti))
     + call md_abort('Failed to allocate sti',me)
      else
      if(.not.ma_push_get(mt_dbl,1,'sti',l_sti,i_sti))
     + call md_abort('Failed to allocate sti',me)
      endif
c
      call argos_cafe_rdtop(lfntop,filtop,byte_mb(i_snam))
      call argos_cafe_topol_init(lfnout)
c
c     distance restraints
c     -------------------
c
      if(ndistr.gt.0) then
      if(me.eq.0) then
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=9999)
      rewind(lfntop)
    5 continue
      read(lfntop,100,end=9999,err=9999) string
  100 format(a3)
      if(string.ne.'noe') goto 5
      read(lfntop,1000) ndrs
 1000 format(i5)
      endif
      call ga_brdcst(mcf_55,ndrs,ma_sizeof(mt_int,1,mt_byte),0)
      if(ndrs.gt.0) then
      if(.not.ma_push_get(mt_int,2*ndrs,'idrs',l_idrs,i_idrs))
     + call md_abort('Failed to allocate idrs',0)
      if(.not.ma_push_get(mt_dbl,6*ndrs,'rdrs',l_rdrs,i_rdrs))
     + call md_abort('Failed to allocate rdrs',0)
      if(.not.ma_push_get(mt_dbl,6*ndrs,'xdrs',l_xdrs,i_xdrs))
     + call md_abort('Failed to allocate xdrs',0)
      endif
      call argos_cafe_rddrs(lfntop,int_mb(i_idrs),dbl_mb(i_rdrs))
      if(me.eq.0) then
      close(unit=lfntop)
      endif
      endif
c
c     proton hopping: donor-acceptor pair list allocation
c     ---------------------------------------------------
c
      if(nhop.gt.0) then
      if(.not.ma_push_get(mt_int,16*nhop*3,'lda',l_lda,i_lda))
     + call md_abort('Failed to allocate lda',0)
      if(.not.ma_push_get(mt_dbl,11*nhop*3,'rda',l_rda,i_rda))
     + call md_abort('Failed to allocate rda',0)
      if(.not.ma_push_get(mt_dbl,4*nhop*3,'uda',l_uda,i_uda))
     + call md_abort('Failed to allocate uda',0)
      if(.not.ma_push_get(mt_dbl,nhop*3,'pda',l_pda,i_pda))
     + call md_abort('Failed to allocate pda',0)
      if(.not.ma_push_get(mt_int,nhop*30,'lsthop',l_lsthop,i_lsthop))
     + call md_abort('Failed to allocate lsthop',0)
      if(.not.ma_push_get(mt_dbl,nhop*15,'timhop',l_timhop,i_timhop))
     + call md_abort('Failed to allocate timhop',0)
      endif
c
c     potential of mean force
c     -----------------------
c
      if(lpmf) then
      if(me.eq.0) then
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='old',err=9998)
      rewind(lfntop)
    6 continue
      read(lfntop,100,end=9998,err=9998) string
      if(string.ne.'pmf') goto 6
      read(lfntop,2000) itemp
 2000 format(2i5)
      endif
      call ga_brdcst(mcf_59,itemp,ma_sizeof(mt_int,2,mt_byte),0)
      numpmf=itemp(1)
      npmfa=itemp(2)
      if(numpmf.gt.0) then
      if(.not.ma_push_get(mt_int,8*numpmf,'ipmf',l_ipmf,i_ipmf))
     + call md_abort('Failed to allocate ipmf',0)
      if(.not.ma_push_get(mt_int,4*numpmf*npmfa,'jpmf',l_jpmf,i_jpmf))
     + call md_abort('Failed to allocate jpmf',npmfa)
      if(.not.ma_push_get(mt_dbl,18*numpmf,'rpmf',l_rpmf,i_rpmf))
     + call md_abort('Failed to allocate rpmf',0)
      if(.not.ma_push_get(mt_dbl,16*numpmf,'xpmf',l_xpmf,i_xpmf))
     + call md_abort('Failed to allocate xpmf',0)
      if(.not.ma_push_get(mt_dbl,12*numpmf,'ypmf',l_ypmf,i_ypmf))
     + call md_abort('Failed to allocate ypmf',0)
      if(.not.ma_push_get(mt_dbl,4*numpmf,'wpmf',l_wpmf,i_wpmf))
     + call md_abort('Failed to allocate wpmf',0)
      if(.not.ma_push_get(mt_dbl,numpmf,'upmf',l_upmf,i_upmf))
     + call md_abort('Failed to allocate upmf',0)
      endif
      call argos_cafe_rdpmf(lfnout,lfntop,int_mb(i_ipmf),int_mb(i_jpmf),
     + dbl_mb(i_rpmf))
      call ga_brdcst(mcf_75,nbias,ma_sizeof(mt_int,1,mt_byte),0)
      if(me.eq.0) then
      close(unit=lfntop)
      endif
      endif
c
      if(ithint) then
      do 3 i=1,24
      ip2(i)=ith(i)
      ip3(i)=ith(i)
    3 continue
      endif
c
c     particle-mesh Ewald initialization
c     ----------------------------------
c
      ipme=jpme
      morder=jorder
      ngx=jgx
      ngy=jgy
      ngz=jgz
      ngmax=max(ngx,ngy,ngz)
      ngrx=ngx+morder
      ngry=ngy+morder
      ngrz=ngz
      pmetol=spmet
      if(ipme.gt.0) then
      if(morder.gt.25) call md_abort('morder too large',0)
      call argos_cafe_alpha
      call argos_pme_start(alpha,morder,1,nodpme,
     + ngx,ngy,ngz,mwm,mwa,msa,icntrl,nbget)
      endif
c
      call argos_cafe_pardif(dbl_mb(i_mas),dbl_mb(i_vdw),dbl_mb(i_chg),
     + int_mb(i_iwa),int_mb(i_iwq),
     + mbt(1),numb(1),mbp(1),dbl_mb(i_bnd(1)),
     + mht(1),numh(1),mhp(1),dbl_mb(i_ang(1)),
     + mdt(1),numd(1),mdp(1),dbl_mb(i_dih(1)),
     + mit(1),numi(1),mip(1),dbl_mb(i_imp(1)),
     + mbt(2),mbp(2),dbl_mb(i_bnd(2)),mht(2),mhp(2),dbl_mb(i_ang(2)),
     + mdt(2),mdp(2),dbl_mb(i_dih(2)),mit(2),mip(2),dbl_mb(i_imp(2)))
c
      if(me.eq.0.and.nsf.ne.nsfi) then
      write(*,'(a,a)') ' Number of fractions differs on topology ',
     + ' and restart files'
      endif
c
c     initialize QHOP parameters
c
      if(nfhop.gt.0) call qhop_init(lfntop,filtop,lfnhop,me)
c
      nsfi=nsf
      nbiasi=nbias
      mropt=mropti
      nipmf=npmfi
c
      return
 
 9998 continue
      call md_abort('Potentials of mean force input not found',0)
 9999 continue
      call md_abort('Distance restraints file not found',0)
c
      return
      end
