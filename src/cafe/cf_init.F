      subroutine cf_start(lfntop,filtop,
     + nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai,
     + mdalg,npbt,nbxt,rcs,rcl,rqm,bx,
     + jpme,jorder,jgx,jgy,jgz,nodpme,spmet,step,tols,mshw,mshs,
     + iiqmmm,iipolt,iitmps,iiprss,rtmpx,rprsx,rtmpw,rtmps,rpres,
     + iislow,tempi,tempwi,tempsi,compr,ntyp,idset,ipset1,ipset2,
     + issscl,delta)
c
      implicit none
c
#include "cf_common.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      integer lfntop
      character*(*) filtop
      integer nwmi,mwmi,nwai,mwai,nsfi,msfi,nsmi,msmi,nsai,msai
      integer jpme,jorder,jgx,jgy,jgz,mshw,mshs
      real*8 spmet,step,tols,bx(3)
      integer mdalg,npbt,nbxt,iiqmmm,iipolt
      real*8 rcs,rcl,rqm,compr,delta
      integer iitmps,iiprss,iislow,nodpme,ntyp,idset,ipset1,ipset2
      integer issscl
      real*8 rtmpx,rprsx,rtmpw,rtmps,rpres,tempi,tempwi,tempsi
c
      integer i
c
      me=ga_nodeid()
      np=ga_nnodes()
c
c     dimensions initially set to zero indicating non-allocated
c
      itscal=iitmps
      ipscal=iiprss
      tmpext=rtmpx
      prsext=rprsx
      tmwrlx=rtmpw
      tmsrlx=rtmps
      prsrlx=rpres
c
      ntype=ntyp
      mdalgo=mdalg
c
      nbs=0
      mscr=0
      lscr=.false.
      llst=.false.
      lpair=.true.
      llist=.false.
      ndxp=0
      maxl=0
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
      rqmmm=rqm
      rshrt2=rshrt*rshrt
      rlong2=rlong*rlong
      rqmmm2=rqmmm*rqmmm
      ltwin=rlong.gt.rshrt
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
      iqmmm=iiqmmm
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
      npener=0
      issscl=0
      nrwrec=0
      isolvo=0
c
      lpww=1
      lpsw=1
      lpss=1
c
      npww=1
      npsw=1
      npss=1
      if(rlong.gt.rshrt) then
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
c
      temp=tempi
      tempw=tempwi
      temps=tempsi
c
      q14fac=0.833333d0
      facpsc=compr*tstep/prsrlx
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
      ithint=ntype.eq.3
      ipert2=ithint.or.(iset.eq.1.and.(ipset1.eq.2.or.ipset2.eq.2))
      ipert3=ithint.or.(iset.eq.1.and.(ipset1.eq.2.or.ipset2.eq.2))
      do 2 i=1,24
      ith(i)=.false.
      ip2(i)=.false.
      ip3(i)=.false.
    2 continue
c
      call cf_rdtop(lfntop,filtop)
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
      pmetol=spmet
      if(ipme.gt.0) then
      if(morder.gt.25) call errquit('morder too large',0)
      call cf_alpha
      call pme_start(alpha,morder,1,nodpme,ngx,ngy,ngz,mwm,mwa,msa)
      endif
c
      call cf_pardif(dbl_mb(i_mas),dbl_mb(i_vdw),dbl_mb(i_chg),
     + int_mb(i_iwa),int_mb(i_iwq),
     + mbt(1),numb(1),mbp(1),dbl_mb(i_bnd(1)),
     + mht(1),numh(1),mhp(1),dbl_mb(i_ang(1)),
     + mdt(1),numd(1),mdp(1),dbl_mb(i_dih(1)),
     + mit(1),numi(1),mip(1),dbl_mb(i_imp(1)),
     + mbt(2),mbp(2),dbl_mb(i_bnd(2)),mht(2),mhp(2),dbl_mb(i_ang(2)),
     + mdt(2),mdp(2),dbl_mb(i_dih(2)),mit(2),mip(2),dbl_mb(i_imp(2)))
c
      return
      end
      subroutine cf_inita(nat,nap,nqt,nqp)
c
c     cf_init initializes atom type and charge type information
c
c     in:  nat, integer number of atom types
c          nap, integer number of parameters per atom type
c          nqt, integer number of charge types
c          nqp, integer number of parameters per charge type
c
c     this routine allocates the memory necessary to store all
c     atom type based force field parameters
c
      implicit none
c
      integer nat,nap,nqt,nqp
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      mat=max(1,nat)
      map=max(1,nap)
      mqt=max(1,nqt)
      mqp=max(1,nqp)
c
c     allocate memory
c
c     for the 4 van der Waals parameters (c6,c12,c6t,c12t) 1 : set 1
c                                                          2 : set 2
c                                                          3 : set 3
c                                                          4 : derivative
c                                                          5 : set 1
c                                                          6 : set 2
c
      if(.not.ma_push_get(mt_dbl,mat*mat*map*6,'vdw',l_vdw,i_vdw))
     + call errquit('Allocation failed for vdw',0)
c
c     for atomic masses 1 : mass set 1
c                       2 : mass set 2
c                       3 : mass set 3
c                       4 : mass derivative
c                       5 : mass set 1
c                       6 : mass set 2
c                           
      if(.not.ma_push_get(mt_dbl,mat*6,'mas',l_mas,i_mas))
     + call errquit('Allocation failed for mas',0)
c
c     for atom type names 1 : name set 1
c                         2 : name set 2
c                         3 : name set 3
c
      if(.not.ma_push_get(mt_byte,6*mat*3,'nam',l_nam,i_nam))
     + call errquit('Allocation failed for nam',0)
c
c     for atomic numbers  1 : atomic number set 1
c                         2 : atomic number set 2
c                         3 : atomic number set 3
c
      if(.not.ma_push_get(mt_int,mat*3,'num',l_num,i_num))
     + call errquit('Allocation failed for num',0)
c
      if(.not.ma_push_get(mt_int,mat*mat,'ias',l_ias,i_ias))
     + call errquit('Allocation failed for ias',0)
c
c     for charge types
c
      if(.not.ma_push_get(mt_dbl,mqt*mqp*6,'chg',l_chg,i_chg))
     + call errquit('Allocation failed for chg',0)
c
      return
      end
      subroutine cf_para(ia,nam,wgt,num)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia
      character*6 nam(3)
      real*8 wgt(3)
      integer num(3)
c
      if(ia.lt.1.or.ia.gt.mat) call errquit('Error in para',0)
c
      call cf_copya(ia,nam,wgt,num,byte_mb(i_nam),
     + dbl_mb(i_mas),int_mb(i_num))
c
      if(abs(wgt(1)-wgt(2)).gt.tiny) ip2(13)=.true.
      if(abs(wgt(1)-wgt(3)).gt.tiny) ip3(13)=.true.
      if(abs(wgt(2)-wgt(3)).gt.tiny) ith(13)=.true.
c
      return
      end
c
      subroutine cf_copya(ia,nami,wgti,numa,nam,wgt,num)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia
      character*6 nami(3),nam(mat,3)
      real*8 wgti(3),wgt(mat,6)
      integer numa(3),num(mat,3)
c
      integer i
c
      do 1 i=1,3
      nam(ia,i)=nami(i)
      wgt(ia,i)=wgti(i)
      num(ia,i)=numa(i)
    1 continue
      wgt(ia,4)=wgti(3)-wgti(2)
      wgt(ia,5)=wgti(2)
      wgt(ia,6)=wgti(3)
c
      return
      end
      subroutine cf_parv(ia,ja,vdwai,vdwri)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia,ja
      real*8 vdwai(6),vdwri(6)
c
c     van der Waals parameters are entered here as:
c
c     vdwai(1:6) : c6a (set 1,2,3),  c6b (set 1,2,3)
c     vdwri(1:6) : c12a (set 1,2,3), c12b (set 1,2,3)
c
      if(ia.lt.1.or.ia.gt.mat) call errquit('Error in para',0)
      if(ja.lt.1.or.ja.gt.mat) call errquit('Error in para',0)
c
      call cf_copyv(ia,ja,vdwai,vdwri,dbl_mb(i_vdw))
c
      if(abs(vdwai(1)-vdwai(2)).gt.tiny) ip2(14)=.true.
      if(abs(vdwai(1)-vdwai(3)).gt.tiny) ip3(14)=.true.
      if(abs(vdwai(2)-vdwai(3)).gt.tiny) ith(14)=.true.
      if(abs(vdwri(1)-vdwri(2)).gt.tiny) ip2(14)=.true.
      if(abs(vdwri(1)-vdwri(3)).gt.tiny) ip3(14)=.true.
      if(abs(vdwri(2)-vdwri(3)).gt.tiny) ith(14)=.true.
c
      return
      end
      subroutine cf_parq(ia,chgi)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia
      real*8 chgi(6)
c
      if(ia.lt.1.or.ia.gt.mqt) call errquit('Error in parq',0)
c
      call cf_copyq(ia,chgi,dbl_mb(i_chg))
c
      if(abs(chgi(1)-chgi(3)).gt.tiny) ip2(16)=.true.
      if(abs(chgi(1)-chgi(5)).gt.tiny) ip3(16)=.true.
      if(abs(chgi(3)-chgi(5)).gt.tiny) ith(16)=.true.
      if(abs(chgi(2)-chgi(4)).gt.tiny) ip2(17)=.true.
      if(abs(chgi(2)-chgi(6)).gt.tiny) ip3(17)=.true.
      if(abs(chgi(4)-chgi(6)).gt.tiny) ith(17)=.true.
c
      return
      end
      subroutine cf_copyq(ia,chgi,chg)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia
      real*8 chgi(6),chg(mqt,mqp,6)
c
      integer i
c
      do 1 i=1,3
      chg(ia,1,i)=qfac*chgi(2*i-1)
      chg(ia,2,i)=chgi(2*i)
    1 continue
c
      do 2 i=1,2
      chg(ia,i,4)=chg(ia,i,3)-chg(ia,i,2)
      chg(ia,i,5)=chg(ia,i,2)
      chg(ia,i,6)=chg(ia,i,3)
    2 continue
c
      return
      end
      subroutine cf_copyv(ia,ja,vdwai,vdwri,vdw)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer ia,ja
      real*8 vdwai(6),vdwri(6),vdw(mat,mat,map,6)
c
      integer i,j
c
      do 1 i=1,3
      vdw(ia,ja,1,i)=vdwai(i)
      vdw(ia,ja,2,i)=vdwai(i+3)
      vdw(ia,ja,3,i)=vdwri(i)
      vdw(ia,ja,4,i)=vdwri(i+3)
    1 continue
      do 2 i=1,4
      vdw(ia,ja,i,4)=vdw(ia,ja,i,3)-vdw(ia,ja,i,2)
      vdw(ia,ja,i,5)=vdw(ia,ja,i,2)
      vdw(ia,ja,i,6)=vdw(ia,ja,i,3)
    2 continue
      if(ia.ne.ja) then
      do 3 i=1,4
      do 4 j=1,6
      vdw(ja,ia,i,j)=vdw(ia,ja,i,j)
    4 continue
    3 continue
      endif
c
      return
      end
      subroutine cf_initb(id,ntot,nbt,nbp,nht,nhp,ndt,ndp,nit,nip,
     + ntt,ntp,nxt,nxp,na)
c
c     cf_init initializes the classical forces API routines
c
c     in: nbt, integer number of bond types
c          nbp, integer number of parameters per bond type
c          nht, integer number of angle types
c          nhp, integer number of parameters per angle type
c          ndt, integer number of dihedral types
c          ndp, integer number of parameters per dihedral type
c          nit, integer number of improper dihedral types
c          nip, integer number of parameters per improper dihedral type
c
c     this routine allocates the memory necessary to store all
c     force field parameters
c
      implicit none
c
      integer id,ntot,nbt,nbp,nht,nhp,ndt,ndp,nit,nip,ntt,ntp,nxt,nxp
      integer na
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      nbs=nbs+1
      if(nbs.gt.mbs) call errquit('Error 1 in initb',0)
      if(id.ne.nbs) call errquit('Error 2 in initb',0)
c
      numb(id)=nbt
      numh(id)=nht
      numd(id)=ndt
      numi(id)=nit
      numt(id)=ntt
      numx(id)=nxt
c
      mbt(id)=max(1,nbt)
      mbp(id)=max(1,nbp)
      mht(id)=max(1,nht)
      mhp(id)=max(1,nhp)
      mdt(id)=max(1,ndt)
      mdp(id)=max(1,ndp)
      mit(id)=max(1,nit)
      mip(id)=max(1,nip)
      mtt(id)=max(ntot,ntt)
      mtp(id)=max(2,ntp)
      mxt(id)=max(ntot,nxt)
      mxp(id)=max(2,nxp)
c
c     allocate memory
c
c     for bond types
c
      if(.not.ma_push_get(mt_dbl,mbt(id)*mbp(id)*6,'bnd',
     + l_bnd(id),i_bnd(id)))
     + call errquit('Allocation failed for bnd',0)
      if(.not.ma_push_get(mt_int,mbt(id)*mbp(id)*3,'ibnd',
     + l_ibnd(id),i_ibnd(id)))
     + call errquit('Allocation failed for ibnd',0)
      if(.not.ma_push_get(mt_dbl,mbt(id)*2,'rbnd',
     + l_rbnd(id),i_rbnd(id)))
     + call errquit('Allocation failed for rbnd',0)
c
c     for angle types
c
      if(.not.ma_push_get(mt_dbl,mht(id)*mhp(id)*6,'ang',
     + l_ang(id),i_ang(id)))
     + call errquit('Allocation failed for ang',0)
      if(.not.ma_push_get(mt_int,mht(id)*mhp(id)*4,'iang',
     + l_iang(id),i_iang(id)))
     + call errquit('Allocation failed for iang',0)
      if(.not.ma_push_get(mt_dbl,mht(id)*2,'rang',
     + l_rang(id),i_rang(id)))
     + call errquit('Allocation failed for rang',0)
c
c     for dihedral types
c
      if(.not.ma_push_get(mt_dbl,mdt(id)*mdp(id)*6,'dih',
     + l_dih(id),i_dih(id)))
     + call errquit('Allocation failed for dih',0)
      if(.not.ma_push_get(mt_int,mdt(id)*mdp(id)*5,'idih',
     + l_idih(id),i_idih(id)))
     + call errquit('Allocation failed for idih',0)
      if(.not.ma_push_get(mt_dbl,mdt(id)*2,'rdih',
     + l_rdih(id),i_rdih(id)))
     + call errquit('Allocation failed for rdih',0)
c
c     for improper dihedral types
c
      if(.not.ma_push_get(mt_dbl,mit(id)*mip(id)*6,'imp',
     + l_imp(id),i_imp(id)))
     + call errquit('Allocation failed for imp',0)
      if(.not.ma_push_get(mt_int,mit(id)*mip(id)*5,'iimp',
     + l_iimp(id),i_iimp(id)))
     + call errquit('Allocation failed for iimp',0)
      if(.not.ma_push_get(mt_dbl,mit(id)*2,'rimp',
     + l_rimp(id),i_rimp(id)))
     + call errquit('Allocation failed for rimp',0)
c
c     for third neighbor lists
c
      if(.not.ma_push_get(mt_int,(mtt(id)+1)*mtp(id),'trd',
     + l_itrd(id),i_itrd(id)))
     + call errquit('Allocation failed for trd',0)
c
c     for excluded or non-bonded lists
c
      if(.not.ma_push_get(mt_int,(mxt(id)+1)*mxp(id),'xcl',
     + l_ixcl(id),i_ixcl(id)))
     + call errquit('Allocation failed for xcl',0)
c
      if(id.eq.1.and.mwa.gt.0) then
      if(.not.ma_push_get(mt_int,mwa,'iwa',l_iwa,i_iwa))
     + call errquit('Allocation failed for iw',0)
      if(.not.ma_push_get(mt_int,mwa,'iwq',l_iwq,i_iwq))
     + call errquit('Allocation failed for iwq',0)
      endif
c
      if(id.eq.2.and.msm.gt.0) then
      if(.not.ma_push_get(mt_dbl,msa,'wsm',l_wsm,i_wsm))
     + call errquit('Allocation failed forwsm',0)
      endif
c
      return
      end
      subroutine cf_parwiq(id,ia,iq)
c
c     cf_parwiq
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,ia,iq
c
      if(id.lt.1.or.id.gt.mwa) call errquit('Error in parwiq',0)
      call cf_copwiq(id,ia,iq,int_mb(i_iwa),int_mb(i_iwq),
     + dbl_mb(i_mas),dbl_mb(i_chg),dbl_mb(i_vdw))
c
      return
      end
      subroutine cf_copwiq(id,ia,iq,iwa,iwq,wgt,chg,vdw)
c
      implicit none
c
#include "cf_common.fh"
c
      integer i,j,id,ia,iq,iwa(mwa),iwq(mwa)
      real*8 wgt(mat,6),chg(mqt,mqp,6),vdw(mat,mat,map,6)
c
      iwa(id)=ia
      iwq(id)=iq
c
      if(abs(wgt(ia,1)-wgt(ia,2)).gt.tiny) ip2(1)=.true.
      if(abs(wgt(ia,1)-wgt(ia,3)).gt.tiny) ip3(1)=.true.
      if(abs(wgt(ia,2)-wgt(ia,3)).gt.tiny) ith(1)=.true.
      do 1 j=1,map
      do 2 i=1,mat
      if(abs(vdw(ia,i,j,1)-vdw(ia,i,j,2)).gt.tiny) ip2(2)=.true.
      if(abs(vdw(ia,i,j,1)-vdw(ia,i,j,3)).gt.tiny) ip3(2)=.true.
      if(abs(vdw(ia,i,j,2)-vdw(ia,i,j,3)).gt.tiny) ith(2)=.true.
    2 continue
    1 continue
      if(abs(chg(iq,1,1)-chg(iq,1,2)).gt.tiny) ip2(4)=.true.
      if(abs(chg(iq,1,1)-chg(iq,1,3)).gt.tiny) ip3(4)=.true.
      if(abs(chg(iq,1,2)-chg(iq,1,3)).gt.tiny) ith(4)=.true.
      if(abs(chg(iq,2,1)-chg(iq,2,2)).gt.tiny) ip2(5)=.true.
      if(abs(chg(iq,2,1)-chg(iq,2,3)).gt.tiny) ip3(5)=.true.
      if(abs(chg(iq,2,2)-chg(iq,2,3)).gt.tiny) ith(5)=.true.
c
      return
      end
      subroutine cf_parbnd(id,idp,ip,p)
c
c     cf_parbnd
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(4)
      real*8 p(*)
c
      if(id.le.0.or.id.gt.nbs) call errquit('Error in parbnd',0)
      if(idp.lt.1.or.idp.gt.mbt(id)) call errquit('Error in parbnd',0)
      if(id.eq.2.and.ip(3).lt.0) nmult(1)=nmult(1)+1
c
      call cf_pcopy(idp,ip,p,4,mbp(id),int_mb(i_ibnd(id)),
     + dbl_mb(i_bnd(id)),mbt(id),
     + ip2(6+(id-1)*12),ip3(6+(id-1)*12),ith(6+(id-1)*12))
c
      return
      end
      subroutine cf_parang(id,idp,ip,p)
c
c     cf_parang
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(5)
      real*8 p(*)
c
      if(id.le.0.or.id.gt.nbs) call errquit('Error 1 in parang',0)
      if(idp.lt.1.or.idp.gt.mht(id)) call errquit('Error 2 in parang',0)
      if(id.eq.2.and.ip(4).lt.0) nmult(2)=nmult(2)+1
c
      call cf_pcopy(idp,ip,p,5,mhp(id),int_mb(i_iang(id)),
     + dbl_mb(i_ang(id)),mht(id),
     + ip2(6+(id-1)*12),ip3(6+(id-1)*12),ith(6+(id-1)*12))
c
      return
      end
      subroutine cf_pardih(id,idp,ip,p)
c
c     cf_pardih
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(6)
      real*8 p(*)
c
      if(id.le.0.or.id.gt.nbs) call errquit('Error 1 in pardih',0)
      if(idp.lt.1.or.idp.gt.mdt(id)) call errquit('Error 2 in pardih',0)
      if(id.eq.2.and.ip(5).lt.0) nmult(3)=nmult(3)+1
c
      call cf_pcopy(idp,ip,p,6,mdp(id),int_mb(i_idih(id)),
     + dbl_mb(i_dih(id)),mdt(id),
     + ip2(6+(id-1)*12),ip3(6+(id-1)*12),ith(6+(id-1)*12))
c
      return
      end
      subroutine cf_parimp(id,idp,ip,p)
c
c     cf_parimp
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,idp,ip(6)
      real*8 p(*)
c
      if(id.le.0.or.id.gt.nbs) call errquit('Error 1 in parimp',0)
      if(idp.lt.1.or.idp.gt.mit(id)) call errquit('Error 2 in parimp',0)
      if(id.eq.2.and.ip(5).lt.0) nmult(4)=nmult(4)+1
c
      call cf_pcopy(idp,ip,p,6,mip(id),int_mb(i_iimp(id)),
     + dbl_mb(i_imp(id)),mit(id),
     + ip2(6+(id-1)*12),ip3(6+(id-1)*12),ith(6+(id-1)*12))
c
      return
      end
      subroutine cf_pcopy(id,ipari,pari,idimi,idimp,ipar,par,idimt,
     + ip2log,ip3log,ithlog)
c
      integer id,idimi,idimp,idimt
      integer ipari(idimi),ipar(idimt,idimi)
      real*8 pari(idimp,3),par(idimt,idimp,6)
      logical ip2log,ip3log,ithlog
c
      do 1 i=1,idimi
      ipar(id,i)=ipari(i)
    1 continue
      do 2 i=1,idimp
      par(id,i,1)=pari(i,1)
      par(id,i,2)=pari(i,2)
      par(id,i,3)=pari(i,3)
      par(id,i,4)=pari(i,3)-pari(i,2)
      par(id,i,5)=pari(i,2)
      par(id,i,6)=pari(i,3)
      if(abs(pari(i,1)-pari(i,2)).gt.tiny) ip2log=.true.
      if(abs(pari(i,1)-pari(i,3)).gt.tiny) ip3log=.true.
      if(abs(pari(i,2)-pari(i,3)).gt.tiny) ithlog=.true.
    2 continue
c
      return
      end
      subroutine cf_ndxtrd(id,idx,jdx,n)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,n
      integer idx(n),jdx(n)
c
      if(n.gt.mtt(id)) call errquit('Error in ndxtrd',0)
c
      if(id.eq.1) then
      call cf_icopy(idx,jdx,n,int_mb(i_itrd(id)),mtt(id))
      else
      call cf_jcopy(idx,jdx,n,int_mb(i_itrd(id)),mtt(id))
      call cf_index(int_mb(i_itrd(id)),mtt(id))
      endif
c
      return
      end
      subroutine cf_ndxxcl(id,idx,jdx,n)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer id,n
      integer idx(n),jdx(n)
c
      if(n.gt.mxt(id)) call errquit('Error in ndxxcl',0)
c
      if(id.eq.1) then
      call cf_icopy(idx,jdx,n,int_mb(i_ixcl(id)),mxt(id))
      else
      call cf_jcopy(idx,jdx,n,int_mb(i_ixcl(id)),mxt(id))
      call cf_index(int_mb(i_ixcl(id)),mxt(id))
      endif
c
      return
      end
      subroutine cf_icopy(idx,jdx,n,ip,nip)
c
      implicit none
c
      integer n,nip
      integer idx(n),jdx(n),ip(0:nip,2)
c
      integer i
c
      do 1 i=1,n
      ip(i,1)=idx(i)
      ip(i,2)=jdx(i)
    1 continue
      if(n.lt.nip) then
      do 2 i=n+1,nip
      ip(i,1)=0
      ip(i,2)=0
    2 continue
      endif
c
      return
      end
      subroutine cf_jcopy(idx,jdx,n,ip,nip)
c
      implicit none
c
      integer n,nip
      integer idx(n),jdx(n),ip(0:nip,2)
c
      integer i
c
      do 1 i=1,n
      ip(i,1)=-idx(i)
      ip(i,2)=jdx(i)
    1 continue
      if(n.lt.nip) then
      do 2 i=n+1,nip
      ip(i,1)=0
    2 continue
      endif
c
      return
      end
      subroutine cf_index(ip,np)
c
      implicit none
c
      integer np
      integer ip(0:np,2)
c
      integer i,j
c
      ip(0,1)=0
      do 1 i=1,np
      j=iabs(ip(i,1))
      if(j.gt.0) ip(j,1)=i
    1 continue
      do 2 i=1,np
      if(ip(i,1).lt.0) ip(i,1)=ip(i-1,1)
    2 continue
c
      return
      end
      subroutine cf_salloc()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lenscr
c
      if(lscr) call errquit('Error 1 in cf_salloc',0)
      if(mscr.le.0) call errquit('Error 2 in cf_salloc',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i1',l_s1i1,i_s1i1))
     + call errquit('Failed to allocate scratch array s1i1',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i2',l_s1i2,i_s1i2))
     + call errquit('Failed to allocate scratch array s1i2',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i3',l_s1i3,i_s1i3))
     + call errquit('Failed to allocate scratch array s1i3',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i4',l_s1i4,i_s1i4))
     + call errquit('Failed to allocate scratch array s1i4',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i5',l_s1i5,i_s1i5))
     + call errquit('Failed to allocate scratch array s1i5',0)
c
      if(.not.ma_alloc_get(mt_int,mscr,'s1i6',l_s1i6,i_s1i6))
     + call errquit('Failed to allocate scratch array s1i6',0)
c
      lenscr=max(nsatot,3*mscr)
      if(.not.ma_alloc_get(mt_int,lenscr,'s2i1',l_s2i1,i_s2i1))
     + call errquit('Failed to allocate scratch array s2i1',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r1',l_s1r1,i_s1r1))
     + call errquit('Failed to allocate scratch array s1r1',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r2',l_s1r2,i_s1r2))
     + call errquit('Failed to allocate scratch array s1r2',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r3',l_s1r3,i_s1r3))
     + call errquit('Failed to allocate scratch array s1r3',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r4',l_s1r4,i_s1r4))
     + call errquit('Failed to allocate scratch array s1r4',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r5',l_s1r5,i_s1r5))
     + call errquit('Failed to allocate scratch array s1r5',0)
c
      if(.not.ma_alloc_get(mt_dbl,mscr,'s1r6',l_s1r6,i_s1r6))
     + call errquit('Failed to allocate scratch array s1r6',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mscr,'s3r1',l_s3r1,i_s3r1))
     + call errquit('Failed to allocate scratch array s3r1',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mscr,'s3r2',l_s3r2,i_s3r2))
     + call errquit('Failed to allocate scratch array s3r2',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr1',l_smr1,i_smr1))
     + call errquit('Failed to allocate scratch array smr1',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr2',l_smr2,i_smr2))
     + call errquit('Failed to allocate scratch array smr2',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr3',l_smr3,i_smr3))
     + call errquit('Failed to allocate scratch array smr3',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr4',l_smr4,i_smr4))
     + call errquit('Failed to allocate scratch array smr4',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr5',l_smr5,i_smr5))
     + call errquit('Failed to allocate scratch array smr5',0)
c
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mscr,'smr6',l_smr6,i_smr6))
     + call errquit('Failed to allocate scratch array smr6',0)
c
      if(ipme.gt.0) then
      if(.not.ma_alloc_get(mt_dbl,3*msa,'pmes',l_pmes,i_pmes))
     + call errquit('Failed to allocate scratch array pmes',0)
      if(.not.ma_alloc_get(mt_dbl,3*mwa*mwm,'pmew',l_pmew,i_pmew))
     + call errquit('Failed to allocate scratch array pmew',0)
      lenscr=3*morder*(mwa*mwm+msa)
      if(.not.ma_alloc_get(mt_dbl,lenscr,'theta',l_theta,i_theta))
     + call errquit('Failed to allocate scratch array theta',0)
      if(.not.ma_alloc_get(mt_dbl,lenscr,'dtheta',l_dtheta,i_dtheta))
     + call errquit('Failed to allocate scratch array dtheta',0)
      endif
c
      lscr=.true.
c
      return
      end
      subroutine cf_sfree()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
c      if(.not.lscr) call errquit('Error 1 in cf_sfree',0)
      if(.not.lscr) return
c
      if(ipme.gt.0) then
      if(.not.ma_free_heap(l_dtheta))
     + call errquit('Failed to free scratch array dtheta',0)
      if(.not.ma_free_heap(l_theta))
     + call errquit('Failed to free scratch array theta',0)
      if(.not.ma_free_heap(l_pmew))
     + call errquit('Failed to free scratch array pmew',0)
      if(.not.ma_free_heap(l_pmes))
     + call errquit('Failed to free scratch array pmes',0)
      endif
c
      if(.not.ma_free_heap(l_smr6))
     + call errquit('Failed to free scratch array smr6',0)
c
      if(.not.ma_free_heap(l_smr5))
     + call errquit('Failed to free scratch array smr5',0)
c
      if(.not.ma_free_heap(l_smr4))
     + call errquit('Failed to free scratch array smr4',0)
c
      if(.not.ma_free_heap(l_smr3))
     + call errquit('Failed to free scratch array smr3',0)
c
      if(.not.ma_free_heap(l_smr2))
     + call errquit('Failed to free scratch array smr2',0)
c
      if(.not.ma_free_heap(l_smr1))
     + call errquit('Failed to free scratch array smr1',0)
c
      if(.not.ma_free_heap(l_s3r2))
     + call errquit('Failed to free scratch array s3r2',0)
c
      if(.not.ma_free_heap(l_s3r1))
     + call errquit('Failed to free scratch array s3r1',0)
c
      if(.not.ma_free_heap(l_s1r6))
     + call errquit('Failed to free scratch array s1r6',0)
c
      if(.not.ma_free_heap(l_s1r5))
     + call errquit('Failed to free scratch array s1r5',0)
c
      if(.not.ma_free_heap(l_s1r4))
     + call errquit('Failed to free scratch array s1r4',0)
c
      if(.not.ma_free_heap(l_s1r3))
     + call errquit('Failed to free scratch array s1r3',0)
c
      if(.not.ma_free_heap(l_s1r2))
     + call errquit('Failed to free scratch array s1r2',0)
c
      if(.not.ma_free_heap(l_s1r1))
     + call errquit('Failed to free scratch array s1r1',0)
c
      if(.not.ma_free_heap(l_s2i1))
     + call errquit('Failed to free scratch array s2i1',0)
c
      if(.not.ma_free_heap(l_s1i6))
     + call errquit('Failed to free scratch array s1i6',0)
c
      if(.not.ma_free_heap(l_s1i5))
     + call errquit('Failed to free scratch array s1i5',0)
c
      if(.not.ma_free_heap(l_s1i4))
     + call errquit('Failed to free scratch array s1i4',0)
c
      if(.not.ma_free_heap(l_s1i3))
     + call errquit('Failed to free scratch array s1i3',0)
c
      if(.not.ma_free_heap(l_s1i2))
     + call errquit('Failed to free scratch array s1i2',0)
c
      if(.not.ma_free_heap(l_s1i1))
     + call errquit('Failed to free scratch array s1i1',0)
c
      lscr=.false.
c
      return
      end
      subroutine cf_lalloc()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      if(llst) call errquit('Error 1 in cf_lalloc',0)
c
      maxl=ma_inquire_avail(mt_int)/2
      if(.not.ma_alloc_get(mt_int,maxl,'list',l_list,i_list))
     + call errquit('Failed to allocate list',0)
c
      llst=.true.
      llist=.false.
c
      return
      end
      subroutine cf_lfree()
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
c      if(.not.llst) call errquit('Error 1 in cf_lfree',0)
c
      if(llst) then
      if(.not.ma_free_heap(l_list))
     + call errquit('Failed to free list',0)
      endif
c
      llst=.false.
      llist=.false.
c
      return
      end
      subroutine cf_free()
c
      implicit none
c
      call cf_lfree()
      call cf_sfree()
c
      return
      end
      subroutine cf_init(lp,llng,bx,vl,vli,zwi,zs,eww,esw,ess)
c     + vw,vwt,vs,vst,zwi,zs,eww,esw,ess)
c
c     in r*8 : bx(3) : box dimensions
c     in log : lp    : flag to force pairlist recalculation
c
      implicit none
c
#include "cf_common.fh"
c
c      real*8 vw(mwm,3,mwa),vwt(mwm,3,mwa),vs(msa,3),vst(msa,3)
      real*8 bx(3),vl(3,3),vli(3,3),zwi(3,2),zs(msf,3,2)
      real*8 eww(mpe,2),esw(msf,mpe,2),ess(msf,msf,mpe,2)
      logical lp,llng
c
      integer i,j,k,l,m
c
      lpair=lp
      llong=llng
c
      ntest=0
      m=1
      if(llong) m=2
c
      do 1 i=1,3
      box(i)=bx(i)
      boxh(i)=half*bx(i)
      do 2 j=1,3
      vlat(i,j)=vl(i,j)
      vlati(i,j)=vli(1,3)
    2 continue
    1 continue
      volume=box(1)*box(2)*box(3)
c
      do 8 l=1,m
      do 3 i=1,mpe
      eww(i,l)=zero
      do 4 j=1,msf
      esw(j,i,l)=zero
      do 5 k=1,msf
      ess(k,j,i,l)=zero
    5 continue
    4 continue
    3 continue
      do 6 i=1,3
      zw(i,l)=zero
      zwi(i,l)=zero
      do 7 j=1,msf
      zs(j,i,l)=zero
    7 continue
    6 continue
      if(ntype.eq.3) then
      do 9 i=1,24
      deriv(i,l)=zero
    9 continue
      endif
      ep2(l)=zero
      ep3(l)=zero
    8 continue
c
      lpww=1
      lpsw=1
      lpss=1
      if(llong) then
      lpww=2
      lpsw=2
      lpss=2
      endif
c
      if(lpair) then
      llww=0
      lsww=0
      llsw=0
      lssw=0
      llss=0
      lsss=0
      endif
c
      if(.not.llist) lpair=.true.
c
c     allocate memory for the scratch arrays
c
      if(.not.lscr) call cf_salloc()
c
c     allocate memory for the pairlists
c
      if(.not.llst) call cf_lalloc()
c
      ndxp=0
c
      if(ipme.gt.0) call pme_init()
c
      return
      end
      subroutine cf_lstmul(ixmul,imul,
     + msb,ibnd,msh,iang,msd,idih,mso,iimp)
c
      implicit none
c
#include "cf_common.fh"
c
      integer ixmul(mmuli),imul(mmult,4)
      integer msb,msh,msd,mso
      integer ibnd(msb,3),iang(msh,4),idih(msd,5),iimp(mso,5)
c
      integer i,j,k,m
c
      k=0
      m=0
c
      if(nmult(1).gt.0) then
      do 1 i=1,numb(2)
      if(ibnd(i,3).lt.0) then
      m=m+1
      ixmul(m)=i
      do 2 j=1,k
      if(ibnd(i,1).eq.imul(k,1)) goto 3
    2 continue
      k=k+1
      imul(k,1)=ibnd(i,1)
    3 continue
      do 4 j=1,k
      if(ibnd(i,2).eq.imul(k,1)) goto 1
    4 continue
      k=k+1
      imul(k,1)=ibnd(i,2)
      endif
    1 continue
      endif
c
      if(nmult(2).gt.0) then
      do 5 i=1,numh(2)
      if(iang(i,4).lt.0) then
      m=m+1
      ixmul(m)=i
      do 6 j=1,k
      if(iang(i,1).eq.imul(k,1)) goto 7
    6 continue
      k=k+1
      imul(k,1)=iang(i,1)
    7 continue
      do 8 j=1,k
      if(iang(i,2).eq.imul(k,1)) goto 9
    8 continue
      k=k+1
      imul(k,1)=iang(i,2)
    9 continue
      do 10 j=1,k
      if(iang(i,3).eq.imul(k,1)) goto 5
   10 continue
      k=k+1
      imul(k,1)=iang(i,3)
      endif
    5 continue
      endif
c
      if(nmult(3).gt.0) then
      do 11 i=1,numd(2)
      if(idih(i,5).lt.0) then
      m=m+1
      ixmul(m)=i
      do 12 j=1,k
      if(idih(i,1).eq.imul(k,1)) goto 13
   12 continue
      k=k+1
      imul(k,1)=idih(i,1)
   13 continue
      do 14 j=1,k
      if(idih(i,2).eq.imul(k,1)) goto 15
   14 continue
      k=k+1
      imul(k,1)=idih(i,2)
   15 continue
      do 16 j=1,k
      if(idih(i,3).eq.imul(k,1)) goto 17
   16 continue
      k=k+1
      imul(k,1)=idih(i,3)
   17 continue
      do 18 j=1,k
      if(idih(i,4).eq.imul(k,1)) goto 11
   18 continue
      k=k+1
      imul(k,1)=idih(i,4)
      endif
   11 continue
      endif
c
      if(nmult(4).gt.0) then
      do 19 i=1,numi(2)
      if(iimp(i,5).lt.0) then
      m=m+1
      ixmul(m)=i
      do 20 j=1,k
      if(iimp(i,1).eq.imul(k,1)) goto 21
   20 continue
      k=k+1
      imul(k,1)=iimp(i,1)
   21 continue
      do 22 j=1,k
      if(iimp(i,2).eq.imul(k,1)) goto 23
   22 continue
      k=k+1
      imul(k,1)=iimp(i,2)
   23 continue
      do 24 j=1,k
      if(iimp(i,3).eq.imul(k,1)) goto 25
   24 continue
      k=k+1
      imul(k,1)=iimp(i,3)
   25 continue
      do 26 j=1,k
      if(iimp(i,4).eq.imul(k,1)) goto 19
   26 continue
      k=k+1
      imul(k,1)=iimp(i,4)
      endif
   19 continue
      endif
c
      nmul=k
c
      return
      end
      subroutine cf_pardif(wgt,vdw,chg,iwat,iwqt,mwb,nwb,nbp,bndw,
     + mwh,nwh,nhp,angw,mwd,nwd,ndp,dihw,mwo,nwo,nop,oopw,
     + msb,npb,bnds,msh,nph,angs,msd,npd,dihs,mso,npo,oops)
c
      implicit none
c
#include "cf_common.fh"
c
      real*8 wgt(mat,6),vdw(mat,mat,map,6),chg(mqt,mqp,6)
      integer iwat(mwa),iwqt(mwa)
      integer mwb,mwh,mwd,mwo,nbp,nhp,ndp,nop
      integer nwb,nwh,nwd,nwo
      integer msb,msh,msd,mso,npb,nph,npd,npo
      real*8 bndw(mwb,nbp,6),angw(mwh,nhp,6)
      real*8 dihw(mwd,ndp,6),oopw(mwo,nop,6)
      real*8 bnds(msb,npb,6),angs(msh,nph,6)
      real*8 dihs(msd,npd,6),oops(mso,npo,6)
c
      return
      end
      subroutine cf_print_top(lfnout)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lfnout
c
      if(me.eq.0) then
      write(lfnout,1000)
 1000 format(/,' TOPOLOGY')
c
      call cf_print_atypes(lfnout,byte_mb(i_nam),int_mb(i_num),
     + dbl_mb(i_mas),dbl_mb(i_vdw))
      call cf_print_ctypes(lfnout,dbl_mb(i_chg))
c      call cf_print_solvent(lfnout)
c      call cf_print_solute(lfnout)
      endif
c
      return
      end
      subroutine cf_print_atypes(lfnout,nam,num,wgt,vdw)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfnout
      character*6 nam(mat,3)
      integer num(mat,3)
      real*8 wgt(mat,6),vdw(mat,mat,map,6)
c
      integer i,j
c
      write(lfnout,1000)
 1000 format(/,' ATOM TYPES',//,t8,'Type',t27,'Atomic number',
     + t44,'Mass',/)
c
      do 1 i=1,mat
      write(lfnout,1001) i,(nam(i,j)(1:5),j=1,3),(num(i,j),j=1,3),
     + (wgt(i,j),j=1,3)
 1001 format(i5,1x,3(1x,a5),3i5,3f12.6)
    1 continue
c
      return
      end
      subroutine cf_print_ctypes(lfnout,chg)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfnout
      real*8 chg(mqt,mqp,6)
c
      integer i,j
c
      write(lfnout,1000)
 1000 format(/,' CHARGE TYPES',//,t10,'Charge',t49,'Polarization',/)
c
      do 1 i=1,mqt
      write(lfnout,1001) i,(chg(i,1,j)/qfac,j=1,3),(chg(i,2,j),j=1,3)
 1001 format(i5,3f12.6,3x,3f12.6)
    1 continue
c
      return
      end
      subroutine cf_print_solvent(lfnout)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfnout
c
      write(lfnout,1000)
 1000 format(/,' SOLVENT',//)
      return
      end
      subroutine cf_print_solute(lfnout)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lfnout
c
      write(lfnout,1000)
 1000 format(/,' SOLUTE',//)
      return
      end
      subroutine cf_lambda(lamtyp,lambda,maxlam,explam,lfnout,
     + rlambd,dlambd)
c
      implicit none
c
#include "cf_common.fh"
#include "mafdecls.fh"
c
      integer lamtyp,lambda,maxlam,lfnout
      real*8 explam,rlambd,dlambd
c
      call cf_lam(lamtyp,lambda,maxlam,explam,lfnout,rlambd,dlambd,
     + dbl_mb(i_mas),dbl_mb(i_vdw),dbl_mb(i_chg),
     + mbt(1),numb(1),mbp(1),dbl_mb(i_bnd(1)),
     + mht(1),numh(1),mhp(1),dbl_mb(i_ang(1)),
     + mdt(1),numd(1),mdp(1),dbl_mb(i_dih(1)),
     + mit(1),numi(1),mip(1),dbl_mb(i_imp(1)),
     + mbt(2),mbp(2),dbl_mb(i_bnd(2)),mht(2),mhp(2),dbl_mb(i_ang(2)),
     + mdt(2),mdp(2),dbl_mb(i_dih(2)),mit(2),mip(2),dbl_mb(i_imp(2)))
c
      return
      end
      subroutine cf_lam(lamtyp,lambda,maxlam,explam,lfnout,
     + rlambd,dlambd,
     + wgt,vdw,chg,mwb,nwb,nbp,bndw,
     + mwh,nwh,nhp,angw,mwd,nwd,ndp,dihw,mwo,nwo,nop,oopw,
     + msb,npb,bnds,msh,nph,angs,msd,npd,dihs,mso,npo,oops)
c
      implicit none
c
#include "cf_common.fh"
c
      integer lamtyp,lambda,maxlam,lfnout
      real*8 explam,rlambd,dlambd
      real*8 wgt(mat,6),vdw(mat,mat,map,6),chg(mqt,mqp,6)
      integer mwb,mwh,mwd,mwo,nbp,nhp,ndp,nop
      integer nwb,nwh,nwd,nwo
      integer msb,msh,msd,mso,npb,nph,npd,npo
      real*8 bndw(mwb,nbp,6),angw(mwh,nhp,6)
      real*8 dihw(mwd,ndp,6),oopw(mwo,nop,6)
      real*8 bnds(msb,npb,6),angs(msh,nph,6)
      real*8 dihs(msd,npd,6),oops(mso,npo,6)
c
      real*8 tlam(3),rlam(3),slam(3)
      integer i,j,k,l
      character*10 pdate,ptime
c
      tlam(1)=dble(lambda-1)/dble(maxlam-1)
      tlam(2)=(dble(lambda-1)-half)/dble(maxlam-1)
      tlam(3)=(dble(lambda-1)+half)/dble(maxlam-1)
      if(tlam(2).lt.zero) tlam(2)=zero
      if(tlam(3).gt.one) tlam(3)=one
c
      if(iabs(lamtyp).eq.1) then
      rlam(1)=tlam(1)
      rlam(2)=tlam(2)
      rlam(3)=tlam(3)
      elseif(iabs(lamtyp).eq.2) then
      rlam(1)=tlam(1)**explam
      rlam(2)=tlam(2)**explam
      rlam(3)=tlam(3)**explam
      elseif(iabs(lamtyp).eq.3) then
      rlam(1)=one-(one-tlam(1))**explam
      rlam(2)=one-(one-tlam(2))**explam
      rlam(3)=one-(one-tlam(3))**explam
      endif
      if(lamtyp.lt.0) then
      rlam(1)=one-rlam(1)
      rlam(2)=one-rlam(2)
      rlam(3)=one-rlam(3)
      endif
      slam(1)=one-rlam(1)
      slam(2)=one-rlam(2)
      slam(3)=one-rlam(3)
c
      do 1 k=1,3
      do 2 i=1,mat
      wgt(i,k)=rlam(k)*wgt(i,6)+slam(k)*wgt(i,5)
      do 3 l=1,map
      do 4 j=1,mat
      vdw(i,j,l,k)=rlam(k)*vdw(i,j,l,6)+slam(k)*vdw(i,j,l,5)
    4 continue
    3 continue
    2 continue
      do 5 j=1,mqp
      do 6 i=1,mqt
      chg(i,j,k)=rlam(k)*chg(i,j,6)+slam(k)*chg(i,j,5)
    6 continue
    5 continue
      do 7 j=1,nbp
      do 8 i=1,mwb
      bndw(i,j,k)=rlam(k)*bndw(i,j,6)+slam(k)*bndw(i,j,5)
    8 continue
    7 continue
      do 9 j=1,nhp
      do 10 i=1,mwh
      angw(i,j,k)=rlam(k)*angw(i,j,6)+slam(k)*angw(i,j,5)
   10 continue
    9 continue
      do 11 j=1,ndp
      do 12 i=1,mwd
      dihw(i,j,k)=rlam(k)*dihw(i,j,6)+slam(k)*dihw(i,j,5)
   12 continue
   11 continue
      do 13 j=1,nop
      do 14 i=1,mwo
      oopw(i,j,k)=rlam(k)*oopw(i,j,6)+slam(k)*oopw(i,j,5)
   14 continue
   13 continue
      do 15 j=1,nbp
      do 16 i=1,msb
      bnds(i,j,k)=rlam(k)*bnds(i,j,6)+slam(k)*bnds(i,j,5)
   16 continue
   15 continue
      do 17 j=1,nhp
      do 18 i=1,msh
      angs(i,j,k)=rlam(k)*angs(i,j,6)+slam(k)*angs(i,j,5)
   18 continue
   17 continue
      do 19 j=1,ndp
      do 20 i=1,msd
      dihs(i,j,k)=rlam(k)*dihs(i,j,6)+slam(k)*dihs(i,j,5)
   20 continue
   19 continue
      do 21 j=1,nop
      do 22 i=1,mso
      oops(i,j,k)=rlam(k)*oops(i,j,6)+slam(k)*oops(i,j,5)
   22 continue
   21 continue
      shift0(k)=rlam(k)*shift0(6)+slam(k)*shift0(5)
      shift1(k)=rlam(k)*shift1(6)+slam(k)*shift1(5)
    1 continue
c
      if(me.eq.0) then
      call swatch(pdate,ptime)
      write(lfnout,1000) rlam(1),ith,pdate,ptime
 1000 format(/' THERMODYNAMIC INTEGRATION, LAMBDA=',f10.5,
     + 5x,24l2,t110,2a10)
      endif
c
      rlambd=rlam(1)
      dlambd=rlam(3)-rlam(2)
c
      return
      end
