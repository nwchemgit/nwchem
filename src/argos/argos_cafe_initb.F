      subroutine argos_cafe_initb(id,ntot,nbt,nbp,nht,nhp,ndt,
     + ndp,nit,nip,ntt,ntp,nxt,nxp,na)
c
c     argos_cafe_init initializes the classical forces API routines
c
c     in:  nbt, integer number of bond types
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
#include "argos_cafe_common.fh"
#include "mafdecls.fh"
c
      nbs=nbs+1
      if(nbs.gt.mbs) call md_abort('Error 1 in initb',0)
      if(id.ne.nbs) call md_abort('Error 2 in initb',0)
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
      if(.not.ma_push_get(mt_dbl,mbt(id)*mbp(id)*mset,'bnd',
     + l_bnd(id),i_bnd(id)))
     + call md_abort('Allocation failed for bnd',0)
      if(.not.ma_push_get(mt_int,mbt(id)*mbp(id)*3,'ibnd',
     + l_ibnd(id),i_ibnd(id)))
     + call md_abort('Allocation failed for ibnd',0)
      if(.not.ma_push_get(mt_dbl,mbt(id)*2,'rbnd',
     + l_rbnd(id),i_rbnd(id)))
     + call md_abort('Allocation failed for rbnd',0)
c
c     for angle types
c
      if(.not.ma_push_get(mt_dbl,mht(id)*mhp(id)*mset,'ang',
     + l_ang(id),i_ang(id)))
     + call md_abort('Allocation failed for ang',0)
      if(.not.ma_push_get(mt_int,mht(id)*mhp(id)*4,'iang',
     + l_iang(id),i_iang(id)))
     + call md_abort('Allocation failed for iang',0)
      if(.not.ma_push_get(mt_dbl,mht(id)*2,'rang',
     + l_rang(id),i_rang(id)))
     + call md_abort('Allocation failed for rang',0)
c
c     for Urey-Bradley terms
c
      if(.not.ma_push_get(mt_dbl,mht(id)*2,'rub',
     + l_rub(id),i_rub(id)))
     + call md_abort('Allocation failed for rub',0)
c
c     for dihedral types
c
      if(.not.ma_push_get(mt_dbl,mdt(id)*mdp(id)*mset,'dih',
     + l_dih(id),i_dih(id)))
     + call md_abort('Allocation failed for dih',0)
      if(.not.ma_push_get(mt_int,mdt(id)*mdp(id)*5,'idih',
     + l_idih(id),i_idih(id)))
     + call md_abort('Allocation failed for idih',0)
      if(.not.ma_push_get(mt_dbl,mdt(id)*2,'rdih',
     + l_rdih(id),i_rdih(id)))
     + call md_abort('Allocation failed for rdih',0)
c
c     for improper dihedral types
c
      if(.not.ma_push_get(mt_dbl,mit(id)*mip(id)*mset,'imp',
     + l_imp(id),i_imp(id)))
     + call md_abort('Allocation failed for imp',0)
      if(.not.ma_push_get(mt_int,mit(id)*mip(id)*5,'iimp',
     + l_iimp(id),i_iimp(id)))
     + call md_abort('Allocation failed for iimp',0)
      if(.not.ma_push_get(mt_dbl,mit(id)*2,'rimp',
     + l_rimp(id),i_rimp(id)))
     + call md_abort('Allocation failed for rimp',0)
c
c     for third neighbor lists
c
      if(.not.ma_push_get(mt_int,(mtt(id)+1)*mtp(id),'trd',
     + l_itrd(id),i_itrd(id)))
     + call md_abort('Allocation failed for trd',0)
c
c     for excluded or non-bonded lists
c
      if(.not.ma_push_get(mt_int,(mxt(id)+1)*mxp(id),'xcl',
     + l_ixcl(id),i_ixcl(id)))
     + call md_abort('Allocation failed for xcl',0)
c
      if(id.eq.1.and.mwa.gt.0) then
      if(.not.ma_push_get(mt_int,mwa,'iwa',l_iwa,i_iwa))
     + call md_abort('Allocation failed for iw',0)
      if(.not.ma_push_get(mt_int,mwa,'iwq',l_iwq,i_iwq))
     + call md_abort('Allocation failed for iwq',0)
      endif
c
      if(id.eq.2.and.msm.gt.0) then
      if(.not.ma_push_get(mt_dbl,msm,'wsm',l_wsm,i_wsm))
     + call md_abort('Allocation failed for wsm',0)
c
      endif
c
      return
      end
c $Id$
