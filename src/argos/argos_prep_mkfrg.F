      logical function argos_prep_mkfrg(irtdb,source,iunit,lfnpar,
     + lfnout,iconst,
     + lseq,cseq,mseq,nseq,lsgm,csgm,msgm,nsgm,latm,catm,xatm,qatm,matm,
     + natm,lato,cato,xato,qato,mato,nato,lbnd,mbnd,lbndt,mbndt,
     + lang,mang,ldih,mdih,limp,mimp,
     + maxscf,qscale,llnk,clnk,mlnk,nlnk)
c
c $Id$
c
c     function to prepare missing fragment definitions
c
c     in  : iunit     = dbase logical file number
c           dir_t     = dbase directory name
c           lfnout    = output file logical file number
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
#include "util.fh"
#include "argos_prep_common.fh"
c
      logical argos_prep_hnames,argos_prep_hybrid
      external argos_prep_hnames,argos_prep_hybrid
c
      integer irtdb
      character*80 source
      integer mring,maxscf,iconst
      parameter (mring=2000)
      integer lring(6,mring)
      logical aring(mring)
c
      integer argos_prep_atnum
      logical argos_prep_short,argos_prep_bnd,argos_prep_ang,
     + argos_prep_dih,argos_prep_impctr,argos_prep_bonds
      external argos_prep_atnum
      external argos_prep_short,argos_prep_bnd,argos_prep_ang,
     +argos_prep_dih,argos_prep_impctr,argos_prep_bonds
      logical argos_prep_atype,argos_prep_charge
      external argos_prep_atype,argos_prep_charge
c
      integer iunit,lfnpar,lfnout
      integer mseq,msgm,matm,mato
      integer nseq,nsgm,natm,nato
      integer mbnd,mang,mdih,mimp
      integer mbndt,nbndt
      integer lseq(6,mseq),lsgm(3,msgm),latm(5,matm),lato(5,mato)
      character*5 cseq(2,mseq),csgm(msgm)
      character*6 catm(3,matm),cato(3,mato)
      real*8 xatm(3,matm),qatm(matm),xato(3,mato),qato(mato),qscale,qsum
      integer lbnd(2,mbnd),lang(3,mang),ldih(4,mdih),limp(4,mimp)
      integer lbndt(2,mbndt)
      character*255 filnam
      integer mlnk,nlnk
      integer llnk(4,mlnk)
      character*4 clnk(2,mlnk)
c
      integer lentmp,length,len_f
c
      integer mtyp,i_ltyp,l_ltyp
      integer nring3,nring4,nring5,nring6
c
      integer i,j,natoms,ilo,ihi,iseq,jlo,jhi,ilist
      integer ibnd,nbnd,nang,ndih,nconn,linkid
      logical found
      character*10 date,time
c
      integer i_wrk,l_wrk
c
      argos_prep_mkfrg=.true.
c
      if(util_print('where',print_debug)) then
      write(lfnout,2000)
 2000 format(/,'argos_prep_mkfrg ')
      endif
c
      lentmp=index(dirpar(mdirpar),' ')-1
      ilist=0
c
      nbndt=0
      if(.not.argos_prep_bonds(xatm,latm,matm,1,natm,lbndt,mbndt,nbndt,
     + llnk,clnk,mlnk,nlnk))
     + call md_abort('Error in argos_prep_bonds',0)
c
      if(iconst.gt.0) then
      if(.not.argos_prep_hybrid(xatm,latm,matm,natm,lbndt,mbndt,nbndt))
     + call md_abort('Error in argos_prep_hybrid',0)
      endif
c
      do 1 i=1,nsgm
      if(lsgm(2,i).eq.0) then
      ilist=0
c
c      if(util_print('where',print_debug)) then
      write(lfnout,2005) csgm(i)
 2005 format(/,' Creating fragment for residue ',a5,/)
c      endif
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
      argos_prep_mkfrg=.false.
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
      linkid=3
      do 55 j=ilo,ihi
      if(latm(5,j).ge.3) then
      latm(5,j)=linkid
      linkid=linkid+1
      endif
   55 continue
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
      if(latm(2,j).le.0) latm(2,j)=argos_prep_atnum(catm(2,j)(1:2))
      if(latm(2,j).gt.0) natoms=natoms+1
    5 continue
c
c     complete the list of bonds for the fragment
c     -------------------------------------------
c
      nbnd=0
      if(.not.ma_push_get(mt_int,12*mato,'wrk',l_wrk,i_wrk))
     + call md_abort('Unable to allocate wrk array',0)
      if(.not.argos_prep_bnd(xatm,latm,catm,matm,natm,ilo,ihi,
     + xato,lato,cato,mato,nato,lbnd,mbnd,nbnd,maxscf,
     + llnk,clnk,mlnk,nlnk,iconst,mang,lang,mdih,ldih,
     + int_mb(i_wrk)))
     + call md_abort('argos_prep_bnd failed',9999)
      if(.not.ma_pop_stack(l_wrk))
     + call md_abort('Unable to deallocate wrk array',0)
c
c     redetermine hydrogen names
c     --------------------------
      if(source(1:3).ne.'pdb') then
      if(.not.argos_prep_hnames(lato,cato,mato,nato,lbnd,mbnd,nbnd))
     + call md_abort('argos_prep_hnames failed',9999)
      endif
c
c     complete the list of angles for the fragment
c     --------------------------------------------
c
      nang=0
      if(.not.argos_prep_ang(lbnd,mbnd,nbnd,lang,mang,nang))
     + call md_abort('argos_prep_ang failed',9999)
c
c     complete the list of torsions for the fragment
c     --------------------------------------------
c
      ndih=0
      if(.not.argos_prep_dih(lang,mang,nang,ldih,mdih,ndih))
     + call md_abort('argos_prep_dih failed in argos_prep_mkfrg',9999)
c
c     count number of bonds for each atom
c     -----------------------------------
c
      do 13 j=1,nato
      lato(1,j)=lato(3,j)
      lato(3,j)=0
      if(lato(5,j).gt.0) lato(3,j)=1
   13 continue
      do 14 ibnd=1,nbnd
      lato(3,lbnd(1,ibnd))=lato(3,lbnd(1,ibnd))+1
      lato(3,lbnd(2,ibnd))=lato(3,lbnd(2,ibnd))+1
   14 continue
c
      ilist=1
      endif
c
c     determine center types
c     ----------------------
c
      if(.not.argos_prep_impctr(lfnout,lato,mato,1,1,nato,nato,
     + lbnd,mbnd,nbnd,lang,mang,nang,ldih,mdih,ndih,
     + lring,aring,mring,nring3,nring4,nring5,nring6))
     + call md_abort('argos_prep_impctr failed',9999)
c
c     memory allocation for work array ltyp
c     -------------------------------------
c
      mtyp=nato+50
      if(.not.ma_push_get(mt_int,15*mtyp,'ltyp',l_ltyp,i_ltyp))
     + call md_abort('Memory allocation failed for ltyp',9999)
c
c     determine atom types
c     --------------------
c
      if(.not.argos_prep_atype(lfnout,lfnpar,
     + lato,cato,mato,lbnd,mbnd,nbnd,
     + 1,1,nato,nato,int_mb(i_ltyp),mtyp,lring,aring,mring,
     + nring3,nring4,nring5,nring6,
     + latm,matm,natm,lbndt,mbndt,nbndt)) then
      call md_abort('argos_prep_atype failed',9999)
      endif
c
c     memory deallocation
c     -------------------
c
      if(.not.ma_pop_stack(l_ltyp))
     + call md_abort('Memory deallocation failed for ltyp',9999)
c
c     check if all atom types could be determined
c     -------------------------------------------
c
      found=.true.
      do 16 j=1,nato
      if(lato(2,j).gt.0) then
      if(cato(3,j)(1:1).eq.' ') found=.false.
      endif
   16 continue
c
c     guestimate partial charges
c     --------------------------
c
      if(found) then
      if(.not.argos_prep_charge(irtdb,lfnout,lfnpar,
     + source,jlo,ilo,ihi,jhi,
     + lato,cato,xato,qato,mato,nato,lbnd,mbnd,nbnd,maxscf,qscale))
     + call md_abort('argos_prep_charge failed',9999)
      else
      if(.not.argos_prep_charge(irtdb,lfnout,lfnpar,
     + source,jlo,ilo,ihi,jhi,
     + lato,cato,xato,qato,mato,nato,lbnd,mbnd,nbnd,0,qscale))
     + call md_abort('argos_prep_charge failed',9999)
      endif
c
c     write the fragment file
c     -----------------------
c
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
      filnam=dirpar(mdirpar)(1:lentmp)//csgm(i)(1:length)//'.frg '
      len_f=index(filnam,' ')-1
      if(.not.found) filnam=filnam(1:len_f)//'_TMP'
      len_f=index(filnam,' ')-1
c
      if(util_print('where',print_debug)) then
      write(lfnout,'(a,a)') 'Writing new fragment file ',
     + filnam(1:len_f)
      endif
      open(unit=iunit,file=filnam(1:len_f),form='formatted',
     + status='unknown',err=9999)
c
      call swatch(date,time)
c
      write(iunit,2001) date,time
 2001 format('# This is an automatically generated fragment file',/,
     + '# Atom types and connectivity were derived from coordinates',/,
     + '# Atomic partial charges are crude estimates',/,
     + '# ',2a10,/,'#')
      length=index(csgm(i),' ')-1
      if(length.le.0) length=5
c
c     write: csgm(i) : fragment name preceded by $
c            nato    : number of atoms
c            1       : number of parameter sets / protonation states
c            1       : default protonation state
c            0       : number of Z-matrix definitions
c            csgm(i) : file name
c
      write(iunit,2002) csgm(i)(1:length),nato,1,1,0,csgm(i)(1:length)
 2002 format('$',a,/,4i5,/,a)
      if(util_print('sequence',print_medium)) then
      write(lfnout,3002) csgm(i)(1:length)
 3002 format(/,' Fragment ',t40,a,//,
     + '  num name  type   link cntr  grp pgrp    charge     polarizab',
     + /)
      endif
      found=.true.
      qsum=0.0d0
      do 6 j=1,nato
      if(lato(2,j).gt.0) then
c
      if(util_print('sequence',print_medium)) then
      write(lfnout,3003) j,cato(2,j),cato(3,j),
     + lato(5,j),lato(4,j),0,1,1,qato(j),0.0
 3003 format(i5,1x,a4,2x,a6,5i5,2f12.6)
      qsum=qsum+qato(j)
      endif
      write(iunit,2003) j,cato(2,j),cato(3,j),
     + lato(5,j),lato(4,j),0,1,1,qato(j),0.0
 2003 format(i5,a4,2x,a6,5i5,2f12.6)
      if(cato(3,j)(1:1).eq.' ') found=.false.
      endif
    6 continue
      if(util_print('sequence',print_medium)) then
      write(lfnout,3009) qsum
 3009 format(43x,'------------',/,25x,'total charge ',f12.6)
      endif
c
      nconn=0
      do 15 ibnd=1,nbnd
      if(lbnd(1,ibnd).ge.1.and.lbnd(1,ibnd).le.nato.and.
     + lbnd(2,ibnd).ge.1.and.lbnd(2,ibnd).le.nato) then
      write(iunit,2004) lbnd(1,ibnd),lbnd(2,ibnd)
 2004 format(2i5)
      nconn=nconn+1
      endif
   15 continue
c
      if(nconn.gt.0.and.util_print('sequence',print_medium)) then
      write(lfnout,3004)
      do 18 ibnd=1,nbnd
      if(lbnd(1,ibnd).ge.1.and.lbnd(1,ibnd).le.nato.and.
     + lbnd(2,ibnd).ge.1.and.lbnd(2,ibnd).le.nato)
     + write(lfnout,3005) lbnd(1,ibnd),lbnd(2,ibnd)
   18 continue
 3004 format(/,' Connectivity',/)
 3005 format(5x,i3,'-',i3)
      endif
c
      close(unit=iunit)
      if(util_print('sequence',print_medium)) then
      write(lfnout,3006)
 3006 format(' ')
      endif
c
      if(util_print('sequence',print_medium)) then
      write(lfnout,3007) filnam(1:len_f)
 3007 format(' Created fragment',t40,a,/)
      endif
c
      if(.not.found.and.util_print('sequence',print_none)) then
      write(lfnout,3008) csgm(i)(1:length)
 3008 format(' Unresolved atom types in fragment ',a,/)
      argos_prep_mkfrg=.false.
c      call md_abort('Unresolved atom types',0)
      endif
c
      lsgm(2,i)=-6
c
      endif
    1 continue
c
c      argos_prep_mkfrg=.true.
      return
c
 9999 continue
      argos_prep_mkfrg=.false.
      return
      end
