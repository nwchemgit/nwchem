c
c $Id$
c

      SUBROUTINE md_inpt(inputfile,
     $                   iseed,tstep,nstep,nequl,nprnt,ntype,ncons,
     $    consatm,nbond,bondatm,nshel,shelatm,lveloc)

      implicit none
      character*(*) inputfile

      include 'p_input.inc'
      include 'p_array.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'
      include 'cm_temp.inc'
      include 'cm_ewld.inc'
      include 'cm_cuto.inc'
      include 'cm_pote.inc'
      include 'cm_cons.inc'
      include 'cm_bond.inc'
      include 'cm_shel.inc'

      integer i,j,k
      integer ierr,ifield,istart,inum_char,infile
      integer iseed,ntype,nequl,nstep,nprnt
      integer ncons,consatm,nbond,bondatm,nshel,shelatm
      integer npote,potindex,pkey,patm1,patm2,imin,imax,idif

      real*8 tstep,para1,para2,para3

      character word*4,irecord*1,filename*20

      logical lend,lcoord,lveloc,lvectr

      dimension consatm(mxcons,2),bondatm(mxbond,2),shelatm(mxshel,2)
      dimension istart(irec_len),inum_char(irec_len),irecord(irec_len)

      open(unit=input,file=inputfile)
      infile=input

      nstep=0
      nequl=0
      nprnt=1
      ntype=0
      npote=0
      ncons=0
      nbond=0
      nshel=0
      natms=0
      iseed=620419483
      kmaxx=0
      kmaxy=0
      kmaxz=0
      tstep=0.0
      temp=0.0
      rcut=0.0
      vcut=0.0
      alpha=0.0
      lcoord=.false.
      lveloc=.false.
      lvectr=.false.

!##########################################################################
 100  call read_lin(infile,ifield,istart,inum_char,irecord,lend)
      if(lend)goto 200
      if(ifield.gt.0)then
       call ex_4char(istart(1),irecord,word)
!##########################################################################
       if(word.eq.'INCL')then
        j=0
        do k=istart(2),istart(2)+inum_char(2)
         j=j+1
         filename(j:j)=irecord(k)
        enddo
        write(output,"(/,1x,'Openning include file :',a20)")
     $   filename(1:j)
        open(unit=iinclude,file=filename(1:j))
        infile=iinclude
        goto 100
!##########################################################################
       elseif(word.eq.'VECT')then
        lvectr=.true.
        read(infile,*)latt(1,1),latt(2,1),latt(3,1)
        read(infile,*)latt(1,2),latt(2,2),latt(3,2)
        read(infile,*)latt(1,3),latt(2,3),latt(3,3)
        write(output,"(/,1x,'Simulation cell vectors')")
        write(output,'(3g20.10)')latt(1,1),latt(2,1),latt(3,1)
        write(output,'(3g20.10)')latt(1,2),latt(2,2),latt(3,2)
        write(output,'(3g20.10)')latt(1,3),latt(2,3),latt(3,3)
!##########################################################################
       elseif(word.eq.'COOR')then
        call ex_inter(inum_char(2),istart(2),irecord,natms,ierr)
        lcoord=.true.
        if(ierr.eq.1)then
         write(output,*)'Invalid number of particles'
        endif
        if(natms.gt.mxatms)then
         write(output,*)'Increase mxatms to ',natms
         stop
        endif
        do j=1,natms
         read(infile,*,end=300)atmtype(j),ccc(j,1),ccc(j,2),ccc(j,3)
        enddo
        write(output,"(/,1x,i9,' atom coordinates found')")natms
!##########################################################################
       elseif(word.eq.'VELO')then
        call ex_inter(inum_char(2),istart(2),irecord,natms,ierr)
        lveloc=.true.
        if(ierr.eq.1)then
         write(output,*)'Invalid number of particles'
         stop
        endif
        if(natms.gt.mxatms)then
         write(output,*)'Increase mxatms to ',natms
         stop
        endif
        do j=1,natms
         read(infile,*,end=300)vvv(j,1),vvv(j,2),vvv(j,3)
        enddo
        write(output,"(/,1x,i9,' atom velocities found')")natms
!##########################################################################
       elseif(word.eq.'ATOM')then
        call ex_inter(inum_char(2),istart(2),irecord,ntype,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of particle types'
         stop
        endif
        if(ntype.gt.mxtype)then
         write(output,*)'Increase mxtype to ',ntype
         stop
        endif
        write(output,"(/,1x,'Atom types: ',i6)")ntype
        potindex=0
        do j=1,ntype
         read(infile,*,end=300)typname(j),typmass(j),typchge(j),
     $    typmol(j)
         write(output,'(a8,1x,2(g20.10,1x),i6)')typname(j),typmass(j),
     $      typchge(j),typmol(j)
         do k=j,ntype
          potindex=potindex+1
          potkey(potindex)=0
         enddo
        enddo
!##########################################################################
       elseif(word.eq.'POTE')then
        call ex_inter(inum_char(2),istart(2),irecord,npote,ierr)
        if(npote.gt.mxpote)then
         write(output,*)'Increase mxpote to ',npote
         stop
        endif
        if(ntype.eq.0)then
         write(output,*)'Atom type info needs to be provided first'
         stop
        endif
        write(output,"(/,1x,'Potential parameters:')")
        do j=1,npote
         read(infile,*,end=300)pkey,patm1,patm2,para1,para2,para3
         imin=min(patm1,patm2)
         imax=max(patm1,patm2)
         idif=imax-imin
         potindex=0
         do k=1,imin-1
          potindex=potindex+(ntype-k+1)
         enddo
         potindex=potindex+idif+1
         potkey(potindex)=pkey 
         potpar(potindex,1)=para1
         potpar(potindex,2)=para2
         potpar(potindex,3)=para3
        enddo
        potindex=0
        do j=1,ntype
         do k=j,ntype
         potindex=potindex+1
         if(potkey(potindex).gt.0)then
          write(output,'(3(i3,1x),3(1x,g20.10))')j,k,potkey(potindex),
     $     potpar(potindex,1),potpar(potindex,2),potpar(potindex,3)
          endif
         enddo
        enddo
!##########################################################################
       elseif(word.eq.'CONS')then
        call ex_inter(inum_char(2),istart(2),irecord,ncons,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of constraints'
         stop
        endif
        if(ncons.gt.mxcons)then
         write(output,*)'Increase mxcons to ',ncons
         stop
        endif
        write(output,"(/,1x,'Distance constraints: ',i6)")ncons
        do j=1,ncons
         read(infile,*,end=300)consatm(j,1),consatm(j,2),consdist(j)
         write(output,'(2(i6,1x),g20.10)')consatm(j,1),consatm(j,2),
     $      consdist(j)
        enddo
!##########################################################################
       elseif(word.eq.'BOND')then
        call ex_inter(inum_char(2),istart(2),irecord,nbond,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of bonds'
         stop
        endif
        if(nbond.gt.mxbond)then
         write(output,*)'Increase mxbond to ',nbond
         stop
        endif
        write(output,"(/,1x,'Bonded interactions: ',i6)")nbond
        do j=1,nbond
         read(infile,*,end=300)bondatm(j,1),bondatm(j,2),bondfrce(j),
     $      bonddist(j)
         write(output,'(2(i6,1x),2(g20.10,1x))')bondatm(j,1),
     $      bondatm(j,2),bondfrce(j),bonddist(j)
        enddo
!##########################################################################
       elseif(word.eq.'SHEL')then
        call ex_inter(inum_char(2),istart(2),irecord,nshel,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of shells'
         stop
        endif
        if(nshel.gt.mxshel)then
         write(output,*)'Increase mxshel to ',nshel
         stop
        endif
        write(output,"(/,1x,'Number of core-shell units: ',i6)")nshel
        do j=1,nshel
         read(infile,*,end=300)shelatm(j,1),shelatm(j,2),shelfrce(j)
         write(output,'(2(i6,1x),g20.10)')shelatm(j,1),shelatm(j,2),
     $      shelfrce(j)
        enddo
!##########################################################################
       elseif(word.eq.'SEED')then
        call ex_inter(inum_char(2),istart(2),irecord,iseed,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid seed for random number generator'
         stop
        endif
        write(output,"(/,1x,'Seed for random number generator :',i18)")
     $   iseed
!##########################################################################
       elseif(word.eq.'TEMP')then
        call ex_1real(inum_char(2),istart(2),irecord,temp,ierr)
        if(ierr.eq.2)then
         write(output,*)'Invalid temperature'
         stop
        endif
        write(output,"(/,1x,'Temperature :',3x,f12.4)")temp
!##########################################################################
       elseif(word.eq.'CUTO')then
        call ex_1real(inum_char(2),istart(2),irecord,rcut,ierr)
        if(ierr.eq.2)then
         write(output,*)'Invalid cutoff distance'
         stop
        endif
        write(output,"(/,1x,'Cutoff distance :',3x,f12.4)")rcut
        rcutsq=rcut**2
!##########################################################################
       elseif(word.eq.'VERL')then
        call ex_1real(inum_char(2),istart(2),irecord,vcut,ierr)
        if(ierr.eq.2)then
         write(output,*)'Invalid verlet shell cutoff'
         stop
        endif
        write(output,"(/,1x,'Verlet shell cutoff :',3x,f12.4)")vcut
        vcutsq=vcut**2
!##########################################################################
       elseif(word.eq.'EWAL')then
        call ex_1real(inum_char(2),istart(2),irecord,alpha,ierr)
        if(ierr.eq.2)then
         write(output,*)'Invalid ewald parameter'
         stop
        endif
        write(output,"(/,1x,'Ewald parameter :',3x,f12.4)")alpha
!##########################################################################
       elseif(word.eq.'KVEC')then
        call ex_inter(inum_char(2),istart(2),irecord,kmaxx,ierr)
        call ex_inter(inum_char(3),istart(3),irecord,kmaxy,ierr)
        call ex_inter(inum_char(4),istart(4),irecord,kmaxz,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of k-vectors'
         stop
        endif
        if(kmaxx.gt.mxkvect)then
         write(output,*)'Increase mxkvect to ',kmaxx
         stop
        endif
        if(kmaxy.gt.mxkvect)then
         write(output,*)'Increase mxkvect to ',kmaxy
         stop
        endif
        if(kmaxz.gt.mxkvect)then
         write(output,*)'Increase mxkvect to ',kmaxz
         stop
        endif
        write(output,"(/,1x,'k-vectors :',3i4)")kmaxx,kmaxy,kmaxz
!##########################################################################
       elseif(word.eq.'STEP')then
        call ex_inter(inum_char(2),istart(2),irecord,nstep,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of MD steps'
         stop
        endif
        write(output,"(/,1x,'Number of MD steps :',i12)")nstep
!##########################################################################
       elseif(word.eq.'EQUI')then
        call ex_inter(inum_char(2),istart(2),irecord,nequl,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid number of equilibration steps'
         stop
        endif
        write(output,"(/,1x,'Number of equilibration steps :',i12)")
     $ nequl
!##########################################################################
       elseif(word.eq.'PRIN')then
        call ex_inter(inum_char(2),istart(2),irecord,nprnt,ierr)
        if(ierr.eq.1)then
         write(output,*)'Invalid printing frequency'
         stop
        endif
        write(output,"(/,1x,'Print every ',i9,' steps')")nprnt
!##########################################################################
       elseif(word.eq.'TIME')then
        call ex_1real(inum_char(2),istart(2),irecord,tstep,ierr)
        if(ierr.eq.2)then
         write(output,*)'Invalid integration timestep'
         stop
        endif
        write(output,"(/,1x,'Integration timestep :',3x,f12.4)")tstep
!##########################################################################
       endif
      endif

      goto 100

 200  continue

      if(infile.eq.input)then
       if(nstep.eq.0)then
        write(output,*)'Number of MD steps has not been defined'
        stop
       elseif(ntype.eq.0)then
        write(output,*)'Atom types have not been defined'
        stop
       elseif(natms.eq.0)then
        write(output,*)'Number of atoms has not been defined'
        stop
       elseif(kmaxx.eq.0)then
        write(output,*)'Number of k-vectors in x has not been defined'
        stop
       elseif(kmaxy.eq.0)then
        write(output,*)'Number of k-vectors in y has not been defined'
        stop
       elseif(kmaxz.eq.0)then
        write(output,*)'Number of k-vectors in z has not been defined'
        stop
       elseif(tstep.eq.0.0)then
        write(output,*)'Integration timestep has not been defined'
        stop
       elseif(temp.eq.0.0)then
        write(output,*)'Temperature has not been defined'
        stop
       elseif(rcut.eq.0.0)then
        write(output,*)'Short-range cutoff has not been defined'
        stop
       elseif(vcut.eq.0.0)then
        write(output,*)'Verlet list cutoff has not been defined'
        stop
       elseif(alpha.eq.0.0)then
        write(output,*)'Ewald parameter has not been defined'
        stop
       elseif(.not.lcoord)then
        write(output,*)'Could not find atomic coordinates'
        stop
       elseif(.not.lvectr)then
        write(output,*)'Could not find simulation cell vectors'
        stop
       endif
       write(output,*)
       return
      elseif(infile.eq.iinclude)then
       infile=input
       goto 100
      endif

 300  continue
      write(output,"(/,1x,'EOF reached in keyword ',a4)")word

      END
