      subroutine pre_wrtpmf(lfntop,filtop,numpmf,maxpmf,ipmf,rpmf,iopt,
     + maxgrp,maxatm,igroup,lgroup,mgroup,npmfr,ipmfr,
     + nsarep,msa,nsa,isegm,ismol,csa,iequi,msm,nsm,itemp,jtemp)
c
      implicit none
c
      integer lfntop
      character*(*) filtop
c
      integer pre_atom_index,pre_base
      external pre_atom_index,pre_base
c
      integer numpmf,maxpmf,maxgrp,maxatm,npmfr,ipmfr,nsarep,msa,nsa
      integer msm,nsm
      integer ipmf(maxpmf,5),iopt(maxpmf,3)
      real*8 rpmf(maxpmf,4)
      integer igroup(maxgrp,maxatm),lgroup(maxgrp),mgroup(maxgrp)
      integer isegm(msa),ismol(msa),iequi(msm,4),itemp(msa),jtemp(msa)
      character*16 csa(msa)
      character*80 target
c
      integer i,j,k,l,m,n,io,ibase,jbase,idatom,nume,nump,im
      integer natoms,natomsj
      character*3 string
      real*8 pi
c
      pi=4.0d0*atan(1.0d0)
c
      if(npmfr.le.0) npmfr=1
      if(ipmfr.le.0) ipmfr=0
c
c     ipmf(i,1) : type of function in pmf i
c
c                 1 : distance
c                 2 : angle
c                 3 : torsion
c                 4 : improper
c                 5 : alignment
c                 6 : planar
c                 7 : base pair
c                 8 : z-axis alignment each atom
c                 9 : z-distance
c                10 : z-axis alignment cog
c
c     ipmf(i,n) : group(s) involved in pmf i
c
      m=0
      nume=0
      nump=0
      do 5 i=1,numpmf
c
c     alignments in more than one equivalent molecule
c
      if(iopt(i,1).eq.1) then
      if(ipmf(i,1).eq.5) then
      nump=nump+lgroup(ipmf(i,2))
      nume=1
      if(mgroup(ipmf(i,2)).gt.0) then
      l=iequi(mgroup(ipmf(i,2)),3)
      if(l.ne.0) call md_abort('pmf not on first equivalent molecule',0)
      do 51 j=1,nsm
      if(iequi(j,3).eq.mgroup(ipmf(i,2))) then
      nump=nump+lgroup(ipmf(i,2))
      nume=nume+1
      endif
   51 continue
      m=max(m,nume*lgroup(ipmf(i,2)))
      endif
      endif
c
c     alignments in combination of more than one equivalent molecule
c
      elseif(iopt(i,1).eq.2) then
      if(ipmf(i,1).eq.6.or.ipmf(i,1).eq.8.or.
     + ipmf(i,1).eq.9.or.ipmf(i,1).eq.10) then
      nump=nump+1
      nume=1
      if(mgroup(ipmf(i,2)).gt.0) then
      l=iequi(mgroup(ipmf(i,2)),3)
      if(l.ne.0) call md_abort('pmf not on first equivalent molecule',0)
      do 52 j=1,nsm
      if(iequi(j,3).eq.mgroup(ipmf(i,2))) then
      nume=nume+1
      endif
   52 continue
      m=max(m,nume*lgroup(ipmf(i,2)))
      endif
      endif
c
      else
      if(ipmf(i,1).eq.5) then
      nump=nump+lgroup(ipmf(i,2))
      elseif(ipmf(i,1).eq.7) then
      nump=nump+3
      else
      nump=nump+1
      endif
      endif
      n=-1
      if(ipmf(i,1).eq.1) n=1
      if(ipmf(i,1).eq.2) n=2
      if(ipmf(i,1).eq.3) n=3
      if(ipmf(i,1).eq.4) n=3
      if(ipmf(i,1).eq.5) n=0
      if(ipmf(i,1).eq.6) n=0
      if(ipmf(i,1).eq.7) m=max(m,1)
      if(ipmf(i,1).eq.8) n=0
      if(ipmf(i,1).eq.9) n=1
      if(ipmf(i,1).eq.10) n=0
      if(n.ge.0) then
      do 6 j=2,n+2
      m=max(m,lgroup(ipmf(i,j)))
    6 continue
      endif
    5 continue
c
c     write pmf
c
      if(numpmf.gt.0) then
      open(unit=lfntop,file=filtop(1:index(filtop,' ')-1),
     + form='formatted',status='unknown',err=9999)
    9 continue
      read(lfntop,100,end=10,err=9999) string
  100 format(a3)
      if(string.ne.'pmf') goto 9
      goto 11
   10 continue
c     position right before EOF
c     to avoid problems on chinook 
c     (MV 11/12/2008 2:50 PM )
      backspace(lfntop)
      write(lfntop,101)
  101 format('pmf')
   11 continue
      write(lfntop,1000) nump,m
 1000 format(2i5)
c      do 12 l=1,npmfr
      l=1
      io=0
      do 1 i=1,numpmf
      if(i.le.ipmfr.or.l.eq.1) then
      if(ipmf(i,1).eq.0) call md_abort('Unknown pmf type',ipmf(i,1))
c
c     distance pmf
c
      if(ipmf(i,1).eq.1.or.ipmf(i,1).eq.9) then
      if(iopt(i,1).ne.2) then
      write(lfntop,1001) ipmf(i,1),iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,3),0,0,i
 1001 format(3i5,2f12.6,2e12.5,5i5)
      do 2 k=2,3
      if(lgroup(ipmf(i,k)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,k))
      write(lfntop,1002) (igroup(ipmf(i,k),j)+io,j=1,lgroup(ipmf(i,k)))
 1002 format(10i7)
    2 continue
      else
      natoms=0
      do 111 k=1,lgroup(ipmf(i,2))
      natoms=natoms+1
      itemp(natoms)=igroup(ipmf(i,2),k)
      if(iequi(ismol(itemp(natoms)),3).ne.0)
     + call md_abort('pmf combine not on first equivalent molecule',0)
  111 continue
      do 112 im=1,nsm
      do 113 k=1,lgroup(ipmf(i,2))
      if(iequi(im,3).eq.ismol(igroup(ipmf(i,2),k))) then
      natoms=natoms+1
      itemp(natoms)=igroup(ipmf(i,2),k)+
     + iequi(im,4)-iequi(ismol(igroup(ipmf(i,2),k)),4)
      endif
  113 continue
  112 continue
      if(natoms.le.0)
     + call md_abort('No atoms selected in group',ipmf(i,2))
      natomsj=0
      do 211 k=1,lgroup(ipmf(i,3))
      natomsj=natomsj+1
      jtemp(natomsj)=igroup(ipmf(i,3),k)
      if(iequi(ismol(jtemp(natomsj)),3).ne.0)
     + call md_abort('pmf combine not on first equivalent molecule',0)
  211 continue
      do 212 im=1,nsm
      do 213 k=1,lgroup(ipmf(i,3))
      if(iequi(im,3).eq.ismol(igroup(ipmf(i,3),k))) then
      natomsj=natomsj+1
      jtemp(natomsj)=igroup(ipmf(i,3),k)+
     + iequi(im,4)-iequi(ismol(igroup(ipmf(i,3),k)),4)
      endif
  213 continue
  212 continue
      if(natoms.le.0)
     + call md_abort('No atoms selected in group',ipmf(i,3))
      write(lfntop,1001) ipmf(i,1),iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),natoms,natomsj,0,0,i
      write(lfntop,1002) (itemp(im),im=1,natoms)
      write(lfntop,1002) (jtemp(im),im=1,natomsj)
      endif
      endif
c
c     angle pmf
c
      if(ipmf(i,1).eq.2) then
      write(lfntop,1001) 2,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,4),0,i
      do 3 k=2,4
      if(lgroup(ipmf(i,k)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,k))
      write(lfntop,1002) (igroup(ipmf(i,k),j)+io,j=1,lgroup(ipmf(i,k)))
    3 continue
      endif
c
c     torsion pmf
c
      if(ipmf(i,1).eq.3) then
      write(lfntop,1001) 3,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,5),i
      do 4 k=2,5
      if(lgroup(ipmf(i,k)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,k))
      write(lfntop,1002) (igroup(ipmf(i,k),j)+io,j=1,lgroup(ipmf(i,k)))
    4 continue
      endif
c
c     torsion pmf
c
      if(ipmf(i,1).eq.4) then
      write(lfntop,1001) 4,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),(lgroup(ipmf(i,j)),j=2,5),i
      do 13 k=2,5
      if(lgroup(ipmf(i,k)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,k))
      write(lfntop,1002) (igroup(ipmf(i,k),j)+io,j=1,lgroup(ipmf(i,k)))
   13 continue
      endif
c
c     alignment pmf
c
      if(ipmf(i,1).eq.5) then
      if(lgroup(ipmf(i,2)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,2))
      do 71 im=1,nsm
      if(im.eq.1) then
      io=0
      elseif(iequi(im,3).eq.mgroup(ipmf(i,2))) then
      io=iequi(im,4)-iequi(mgroup(ipmf(i,2)),4)
      else
      goto 71
      endif
      do 7 k=1,lgroup(ipmf(i,2))
      write(lfntop,1001) 5,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),1,lgroup(ipmf(i,2)),0,0,i
      write(lfntop,1002) igroup(ipmf(i,2),k)+io
      write(lfntop,1002) (igroup(ipmf(i,2),j)+io,j=1,lgroup(ipmf(i,2)))
    7 continue
      if(iopt(i,1).ne.1) goto 711
   71 continue
  711 continue
      endif
c
c      6: planar pmf
c      8: z-axis pmf on each atom
c     10: z-axis pmf on COG
c
      if(ipmf(i,1).eq.6.or.ipmf(i,1).eq.8.or.ipmf(i,1).eq.10) then
      if(iopt(i,1).ne.2) then
      write(lfntop,1001) ipmf(i,1),iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),lgroup(ipmf(i,2)),0,0,0,i
      if(lgroup(ipmf(i,2)).le.0)
     + call md_abort('No atoms selected in group',ipmf(i,2))
      write(lfntop,1002) (igroup(ipmf(i,2),j),j=1,lgroup(ipmf(i,2)))
      else
      natoms=0
      do 81 k=1,lgroup(ipmf(i,2))
      natoms=natoms+1
      itemp(natoms)=igroup(ipmf(i,2),k)
      if(iequi(ismol(itemp(natoms)),3).ne.0)
     + call md_abort('pmf combine not on first equivalent molecules',0)
   81 continue
      do 82 im=1,nsm
      do 83 k=1,lgroup(ipmf(i,2))
      if(iequi(im,3).eq.ismol(igroup(ipmf(i,2),k))) then
      natoms=natoms+1
      itemp(natoms)=igroup(ipmf(i,2),k)+
     + iequi(im,4)-iequi(ismol(igroup(ipmf(i,2),k)),4)
      endif
   83 continue
   82 continue
      if(natoms.le.0)
     + call md_abort('No atoms selected in group',ipmf(i,2))
      write(lfntop,1001) ipmf(i,1),iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),natoms,0,0,0,i
      write(lfntop,1002) (itemp(im),im=1,natoms)
      endif
      endif
c
c     base pair pmf
c
      if(ipmf(i,1).eq.7) then
      ibase=pre_base(ipmf(i,2),isegm,csa,msa,nsa)
      jbase=pre_base(ipmf(i,3),isegm,csa,msa,nsa)
      if(ibase.eq.0.or.jbase.eq.0)
     + call md_abort('Error in pmf basepair: not a base',0)
      if((ibase.eq.1.and.jbase.eq.5).or.
     + (ibase.eq.2.and.jbase.eq.3)) then
c
      write(lfntop,1001) 1,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),1,1,0,0,i
      write(target,'(i5,a)') ipmf(i,2),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      write(lfntop,1001) 4,iopt(i,2),iopt(i,3),
     + pi,pi,(rpmf(i,j),j=3,4),1,1,1,1,i
      write(target,'(i5,a)') ipmf(i,2),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_C5 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      write(lfntop,1001) 4,iopt(i,2),iopt(i,3),
     + pi,pi,(rpmf(i,j),j=3,4),1,1,1,1,i
      write(target,'(i5,a)') ipmf(i,3),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_C5 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      endif
      if((ibase.eq.5.and.jbase.eq.1).or.
     + (ibase.eq.3.and.jbase.eq.2)) then
      write(lfntop,1001) 1,iopt(i,2),iopt(i,3),
     + (rpmf(i,j),j=1,4),1,1,0,0,i
      write(target,'(i5,a)') ipmf(i,2),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      write(lfntop,1001) 4,iopt(i,2),iopt(i,3),
     + pi,pi,(rpmf(i,j),j=3,4),1,1,1,1,i
      write(target,'(i5,a)') ipmf(i,2),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_C5 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      write(lfntop,1001) 4,iopt(i,2),iopt(i,3),
     + pi,pi,(rpmf(i,j),j=3,4),1,1,1,1,i
      write(target,'(i5,a)') ipmf(i,3),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_C5 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,3),':_N1 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
      write(target,'(i5,a)') ipmf(i,2),':_N3 '
      idatom=pre_atom_index(target,isegm,csa,msa,nsa)
      write(lfntop,1002) idatom
c
      endif
c
      endif
c
      endif
    1 continue
c   12 continue
c
      close(unit=lfntop)
      endif
c
      return
c
 9999 continue
      call md_abort('Could not read filtop to write pmf',0)
      return
      end
c $Id$
