      logical function pre_mkfrg(irtdb,source,iunit,lfnpar,filpar,
     + dir_s,dir_x,dir_u,dir_t,
     + lfnout,lprint,
     + lseq,cseq,mseq,nseq,lsgm,csgm,msgm,nsgm,latm,catm,xatm,qatm,matm,
     + natm,lbnd,mbnd,lang,mang,ldih,mdih,limp,mimp)
c
c     function to prepare missing fragment definitions
c
c     in  : iunit     = dbase logical file number
c           dir_t     = dbase directory name
c           lfnout    = output file logical file number
c           lprint    = integer for print to output
c           mseq      = dimension of the sequence list
c
c     out : lseq(1,*) = segment numbers
c           lseq(2,*) = number of atoms
c           lseq(3,*) = index to unique segment
c           lsgm(1,i) = number of segments of type i
c           lsgm(2,i) = source: 0=not found; 1=s; 2=x; 3=u; 4=t;
c           cseq      = segment names
c
c
      implicit none
c
#include "mafdecls.fh"
c
      logical pre_hnames
      external pre_hnames
c
      integer irtdb
      character*80 source
      integer mring
      parameter (mring=100)
      integer lring(6,mring)
      logical aring(mring)
c
      integer pre_atnum
      logical pre_short,pre_bnd,pre_ang,pre_dih,pre_impctr
      external pre_atnum
      external pre_short,pre_bnd,pre_ang,pre_dih,pre_impctr
      logical pre_atype,pre_charge
      external pre_atype,pre_charge
c
      integer iunit,lfnpar,lfnout
      integer mseq,msgm,matm
      integer nseq,nsgm,natm
      integer mbnd,mang,mdih,mimp
      integer lseq(4,mseq),lsgm(3,msgm),latm(5,matm)
      character*5 cseq(2,mseq),csgm(msgm)
      character*6 catm(3,matm)
      real*8 xatm(3,matm),qatm(matm)
      integer lbnd(2,mbnd),lang(3,mang),ldih(4,mdih),limp(4,mimp)
      character*255 filpar,filnam
      character*255 dir_s,dir_x,dir_u,dir_t
      integer lprint
c
      integer len_t,length,len_f
c
      integer mtyp,i_ltyp,l_ltyp
      integer nring3,nring4,nring5,nring6
c
      integer i,j,natoms,ilo,ihi,iseq,jlo,jhi,ilist
      integer ibnd,nbnd,nang,ndih,nconn
      logical found
c
      if(lprint.ge.5) then
      write(lfnout,2000)
 2000 format(/,'pre_mkfrg ')
      endif
c
      len_t=index(dir_t,' ')-1
      ilist=0
c
      do 1 i=1,nsgm
      if(lsgm(2,i).eq.0) then
c
c     determine index first segment of this type in sequence list
c     -----------------------------------------------------------
c
      iseq=0
      do 2 j=1,mseq
      if(lseq(2,j).eq.i) then
      iseq=j
      goto 3
      endif
    2 continue
      pre_mkfrg=.false.
      return
    3 continue
c
c     find indices ilo : first atom of this segment
c     ------------ ihi : last atom for this segment
c                  jlo : first atom
c                  jhi : last atom
c
      ilo=lseq(3,iseq)
      ihi=lseq(3,iseq+1)-1
c
      if(ilist.eq.0) then
c
      jlo=1
      jhi=natm
c      if(csgm(i)(5:5).eq.'N'.or.csgm(i)(5:5).eq.'M') jlo=ilo
c      if(csgm(i)(5:5).eq.'C'.or.csgm(i)(5:5).eq.'M') jhi=ihi
c
c     determine the number of 'real' atoms
c     ------------------------------------
c
      natoms=0
      do 5 j=jlo,jhi
      latm(2,j)=pre_atnum(catm(2,j)(1:2))
      if(latm(2,j).gt.0) natoms=natoms+1
    5 continue
c
c     complete the list of bonds for the fragment
c     -------------------------------------------
c
      nbnd=0
      if(.not.pre_bnd(xatm,latm,matm,jlo,jhi,lbnd,mbnd,nbnd))
     + call errquit('pre_bnd failed',9999)
c
c     redetermine hydrogen names
c     --------------------------
      if(source(1:3).ne.'pdb') then
      if(.not.pre_hnames(latm,catm,matm,natm,lbnd,mbnd,nbnd))
     + call errquit('pre_hnames failed',9999)
      endif
c
c     complete the list of angles for the fragment
c     --------------------------------------------
c
      nang=0
      if(.not.pre_ang(lbnd,mbnd,nbnd,lang,mang,nang))
     + call errquit('pre_ang failed',9999)
c
c     complete the list of torsions for the fragment
c     --------------------------------------------
c
      ndih=0
      if(.not.pre_dih(lang,mang,nang,ldih,mdih,ndih))
     + call errquit('pre_dih failed',9999)
c
c     count number of bonds for each atom
c     -----------------------------------
c
      do 13 j=1,natm
      latm(3,j)=0
   13 continue
      do 14 ibnd=1,nbnd
      latm(3,lbnd(1,ibnd))=latm(3,lbnd(1,ibnd))+1
      latm(3,lbnd(2,ibnd))=latm(3,lbnd(2,ibnd))+1
   14 continue
c
      ilist=1
      endif
c
c     determine center types
c     ----------------------
c
      if(.not.pre_impctr(xatm,latm,matm,jlo,ilo,ihi,jhi,
     + lbnd,mbnd,nbnd,lang,mang,nang,ldih,mdih,ndih,
     + lring,aring,mring,nring3,nring4,nring5,nring6))
     + call errquit('pre_impctr failed',9999)
c
c     memory allocation for work array ltyp
c     -------------------------------------
c
      mtyp=jhi-jlo+1
      if(.not.ma_push_get(mt_int,15*mtyp,"ltyp",l_ltyp,i_ltyp))
     + call errquit('Memory allocation failed for ltyp',9999)
c
c     determine atom types
c     --------------------
c
      if(.not.pre_atype(lfnout,lprint,lfnpar,filpar,
     + dir_s,dir_x,dir_u,dir_t,latm,catm,matm,lbnd,mbnd,nbnd,
     + jlo,ilo,ihi,jhi,int_mb(i_ltyp),mtyp,lring,aring,mring,
     + nring3,nring4,nring5,nring6))
     + call errquit('pre_atype failed',9999)
c
c     memory deallocation
c     -------------------
c
      if(.not.ma_pop_stack(l_ltyp))
     + call errquit('Memory deallocation failed for ltyp',9999)
c
c
c     guestimate partial charges
c     --------------------------
c
      if(.not.pre_charge(irtdb,source,jlo,ilo,ihi,jhi,
     + latm,catm,xatm,qatm,matm,natm,lbnd,mbnd,nbnd))
     + call errquit('pre_charge failed',9999)
c
c     write the fragment file
c     -----------------------
c
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
      filnam=dir_t(1:len_t)//csgm(i)(1:length)//'.frg '
      len_f=index(filnam,' ')-1
c
      open(unit=iunit,file=filnam(1:len_f),form='formatted',
     + status='new',err=9999)
c
      write(iunit,2001)
 2001 format('# This is an automatically generated fragment file',/,
     + '# Atom types and connectivity were derived from coordinates',/,
     + '# Atomic partial charges are crude guestimations',/,
     + '#')
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
      write(iunit,2002) csgm(i)(1:length),ihi-ilo+1
 2002 format('$',a,/,i5)
      write(lfnout,3002) csgm(i)(1:length)
 3002 format(/,'FRAGMENT ',a,/,
     + '  num name type   link cntr  grp pgrp    charge     polarizab',
     + /)
      found=.true.
      do 6 j=ilo,ihi
      if(latm(2,j).gt.0) then
c
      write(lfnout,2003) j+1-ilo,catm(2,j),catm(3,j),
     + latm(5,j),latm(4,j),1,1,qatm(j),0.0
      write(iunit,2003) j+1-ilo,catm(2,j),catm(3,j),
     + latm(5,j),latm(4,j),1,1,qatm(j),0.0
 2003 format(i5,a4,2x,a6,4i5,2f12.6)
      if(catm(3,j)(1:1).eq.' ') found=.false.
      endif
    6 continue
c
      nconn=0 
      do 15 ibnd=1,nbnd
      if(lbnd(1,ibnd).ge.ilo.and.lbnd(1,ibnd).le.ihi.and.
     + lbnd(2,ibnd).ge.ilo.and.lbnd(2,ibnd).le.ihi) then
      write(iunit,2004) lbnd(1,ibnd)+1-ilo,lbnd(2,ibnd)+1-ilo
 2004 format(2i5)
      nconn=nconn+1
      lbnd(1,nconn)=lbnd(1,ibnd)+1-ilo
      lbnd(2,nconn)=lbnd(2,ibnd)+1-ilo
      endif
   15 continue
c
      if(nconn.gt.0) then
      write(lfnout,3004) (lbnd(1,ibnd),lbnd(2,ibnd),ibnd=1,nconn)
 3004 format(/,'CONNECTIVITY',/,(12(i3,'-',i3,2x)))
      endif
c
      close(unit=iunit)
c
      write(lfnout,3005) filnam(1:len_f)
 3005 format(/,'CREATED FRAGMENT FILE ',a)
c
      if(.not.found) then
      write(lfnout,3006) csgm(i)(1:length)
 3006 format(/,'ATOM TYPES COULD NOT BE DETERMINED FOR FRAGMENT ',a)
      goto 9999
      endif
c
      lsgm(2,i)=-4
c
      endif
    1 continue
c
      pre_mkfrg=.true.
      return
c
 9999 continue
      pre_mkfrg=.false.
      return
      end
