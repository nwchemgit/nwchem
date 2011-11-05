***********************************************************************
*    This program converts NWChem trajectory file
*    into xyz format.
*        
*    Date:   07/11/2011
*    Author: Marat Valiev (PNNL)
***********************************************************************
       program main
       implicit none
       logical ofile
       character*(180) buffer
       integer itype
       integer i,n
       character*180 file_in
       character*180 file_out
       integer l
       integer istatus
       character*80 message
       integer nvm,nva,nu,nprec,nv
       logical lv(4),lu(4)
       integer fn_in,fn_out
       character*5 tvr
       character*5 , dimension(:), allocatable :: tva,tua,tur
       double precision, dimension(:,:), allocatable :: cv,cu
       integer, dimension(:), allocatable :: iur
       integer ifr,nfr
       logical ok
       logical overb
       logical ohelp
c
c       n = command_argument_count()
c       write(*,*) "number of command line args",n
c      intial values for command line args
       file_in  = " "
       file_out = " "
       overb = .false.
       ohelp = .false.
c      --------------------------------------------      
c      beging parsing command line arguments if any
c      --------------------------------------------      
       i = 0
16     continue
       i = i+1
       call get_command_argument(i,buffer,l,istatus)
c       write(*,*) "argument ",i,buffer
       if(istatus.ne.0) goto 18
       if(buffer.eq."-frames") then
          i = i+1
          call get_command_argument(i,buffer,l,istatus)
          if(istatus.ne.0) goto 18
          read(buffer,*) nfr
          go to 16
       else if(buffer.eq."-v") then
          overb=.true.
          go to 16
       else if(buffer.eq."-help") then
          ohelp=.true.
          write(*,1000)
          go to 14
       else 
          if(file_in .eq. ' ') then
            file_in=buffer
            go to 16
          else if(file_out.eq." ") then
            file_out=buffer
            go to 16
          end if
       end if 
c      ---------------------------      
c      end of command line parsing
c      ---------------------------      
18     continue
c      -----------------------
c      finialize io
c      ----------------------
       if(file_in.eq." ") then
         file_in="nwchem.trj"
       end if
       if(file_out.eq." ") then
         file_out="nwchem.xyz"
       end if
       inquire(file=file_in,exist=ofile)
       fn_in = 10
       if(ofile) then
          open(fn_in,file=file_in,
     $            form='formatted',status='old',err=911)
       else
           message = "no file found: "//file_in
           goto 911
       end if 

       fn_out = 11
       if(ofile) then
               message = "opening output file "//file_out
          open(fn_out,file=file_out,
     $            form='formatted',status='unknown',err=911)
       else
           message = "no file found: "//file_in
           goto 911
       end if 
c      testing
1000   format(
     > "NAME:",/
     > " traj2xyz converts of NWChem trajectory file",//,
     >  "SYNOPSIS",/,
     > " traj2xyz [-help] [-v ] [-frames number of frames]",
     > " input output ",//,
     >  "DESCRIPTION",/,
     > " The traj2xyz utility converts NWChem trajectory file",
     > " into xyz format ",/,
     > " The option are as follows:",/,
     > " -help  prints out this message",
     > " -v  generate more verbose output",
     > " -frames  number of frames to process"
     > )
       call trj_read_header(fn_in,nvm,nva,nu,nprec)
       nv = nva*nvm
       if(overb) then
          write(*,*) "number of solvent molecules",nvm
          write(*,*) "number of solvent atoms per molecule",nva
          write(*,*) "number of solvent atoms ",nv
          write(*,*) "number of solute molecules",nu
          write(*,*) "trajectory precision",nprec
       end if
       allocate(tva(nva))
       allocate(tua(nu))
       allocate(tur(nu))
       allocate(iur(nu))
       allocate(cv(3,nvm*nva))
       allocate(cu(3,nu))
       call trj_read_solvent_specs(fn_in,nva,tvr,tva)
       call trj_read_solute_specs(fn_in,nu,tur,tua,iur)
       ifr = 0
12     ifr = ifr+1
       if(ifr.gt.nfr) goto 14
       call trj_read_coords(fn_in,nv,nu,cv,cu,ok)
       if(ok) then
         call trj_write_coords_xyz(fn_out,ifr,nva,nvm,nu,tva,tua,cv,cu)
         goto 12
       end if
c       call trj_read_dimensions(fn_in,lv,lu,nvm,nva,nu)
14     continue
       return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "Emergency STOP"
       write(*,*) message
       stop
       end 

      subroutine trj_read_solvent_specs(fn,nva,tvr,tva)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: nva
      character*5, intent(out) :: tvr
      character*5, intent(out) :: tva(nva)
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     skip one line      
      read(fn,1000,err=911,end=911) card
c      
      sep = " "
      do i=1,nva
        read(fn,1000,err=911,end=911) card
        i0 = 1
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tvr
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tva(i)
      end do
c      do i=1,nva
c        write(*,*) tvr,tva(i)
c      end do
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_solute_specs(fn,n,tur,tua,iur)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: n
      character*5, intent(out) :: tur(n)
      character*5, intent(out) :: tua(n)
      integer, intent(out) :: iur(n)
c      
      integer nu,nva
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
      sep = " "
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     find how many lines to skip 
      read(fn,1000,err=911,end=911) card
      i0 = 1
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nva
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nu
      if(nu.gt.n) then
        message = "insufficient solute size "
        goto 911
      end if
c     skip solvent block      
      do i=1,nva
        read(fn,1000,err=911,end=911) card
      end do
      do i=1,nu
        read(fn,1000,err=911,end=911) card
        i0 = 1
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tur(i)
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) tua(i)
        call get_next_token(i0,card,sep,token,ostatus)
        read(token,*) iur(i)
      end do
c      do i=1,nu
c        write(*,*) tur(i),tua(i),iur(i)
c      end do
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_header(fn,nvm,nva,nu,nprec)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(out)   :: nvm
      integer, intent(out)   :: nva
      integer, intent(out)   :: nu
      integer, intent(out)   :: nprec
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      rewind(fn)
c      
c     look for the "header"
c     ---------------------
  100 continue
      read(fn,1000,end=911) card
c      write(*,*) card
 1000 format(a)
      if(card(1:6).ne.'header') goto 100
c     parse the line      
      read(fn,1000,err=911,end=911) card
      i0 = 1
      sep = " "
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nva
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nu
      call get_next_token(i0,card,sep,token,ostatus)
      call get_next_token(i0,card,sep,token,ostatus)
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nvm
      call get_next_token(i0,card,sep,token,ostatus)
      read(token,*) nprec
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_write_coords_xyz(fn,ifr,nva,nvm,nu,tva,tua,cv,cu)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: ifr
      integer, intent(in)   :: nva
      integer, intent(in)   :: nvm
      integer, intent(in)   :: nu
      character*5, intent(in) :: tva(nva)
      character*5, intent(in) :: tua(nu)
      double precision, intent(in)   :: cv(3,nva*nvm)
      double precision, intent(in)   :: cu(3,nu)
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c      
      integer nv
      integer j,i1
      nv = nva*nvm
      write(fn,*) nu+nva*nvm
      write(fn,'("frame ",I6)') ifr
      i1 = 0
      do i=1,nvm
      do j=1,nva
        i1 = i1+1
        write(fn,'(A5,3f8.3)') tva(j),(cv(k,i1)*10,k=1,3)
      end do
      enddo
      do i=1,nu
        write(fn,'(A5,3f8.3)') tua(i),(cu(k,i)*10,k=1,3)
      enddo
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_coords(fn,nv,nu,cv,cu,ok)
      implicit none
      integer, intent(in)    :: fn  
      integer, intent(in)   :: nv
      integer, intent(in)   :: nu
      logical, intent(out)  :: ok
      double precision, intent(out)   :: cv(3,nv)
      double precision, intent(out)   :: cu(3,nu)
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c
      ok = .false.
      sep = " "
c      
c     look for the start of the fame record
c     -------------------------------------
  100 continue
      read(fn,1000,end=11,err=911) card
c      write(*,*) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 100
c
c     skip five lines here
c     --------------------
      do i=1,5
        read(fn,1000,err=911,end=911) card
      end do
      do i=1,nv
        read(fn,1000,err=911,end=911) card
        i0 = 1
        do k=1,3
          call get_next_token(i0,card,sep,token,ostatus)
          read(token,*,ERR=911,END=911) cv(k,i)
        end do
      end do
      do i=1,nu
        read(fn,1000,err=911,end=911) card
        i0 = 1
        do k=1,3
          call get_next_token(i0,card,sep,token,ostatus)
          read(token,*,ERR=911,END=911) cu(k,i)
        end do
      end do
      ok = .true.
11    continue
      return
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine trj_read_dimensions(fn,lv,lu,nvm,nva,nu)
      implicit none
      integer, intent(in)    :: fn  
      logical, intent(out)   :: lv(4)
      logical, intent(out)   :: lu(4)
      integer, intent(out)   :: nvm
      integer, intent(out)   :: nva
      integer, intent(out)   :: nu
c      
      character*80 card
      character*80 message
      character*30 token
      character*1 sep
      character*1 a1
      integer i0
      integer i,k,lb,lt
      logical ostatus
c      
c     look for the start of the fame record
c     -------------------------------------
  100 continue
      read(fn,1000,end=911) card
c      write(*,*) card
 1000 format(a)
      if(card(1:5).ne.'frame') goto 100
c
c     skip four lines here
c     --------------------
      do i=1,4
        read(fn,1000,err=911,end=911) card
c        write(*,*) card
      end do
      read(fn,1000,err=911,end=911) card
      i0 = 1
      sep = " "
      call get_next_token(i0,card,sep,token,ostatus)
c      write(*,*) "token=",token
      read(token,1003) (lv(i),i=1,4),(lu(i),i=1,4)
1003  format(8l1,3i10)
      if(lv(4).or.lu(4)) then
         message = "found pressure field record?"
         goto 911
      end if
      call get_next_token(i0,card,sep,token,ostatus)
c      write(*,*) "token=",token
      read(token,*) nvm
      call get_next_token(i0,card,sep,token,ostatus)
c      write(*,*) "token=",token
      read(token,*) nva
      call get_next_token(i0,card,sep,token,ostatus)
c      write(*,*) "token=",token
      read(token,*) nu
      return
cc
c      if(card(8:8).eq.'T'.or.card(8:8).eq.'F') then
c      read(card,1003) lxw,lvw,lfw,lpw,lxs,lvs,lfs,lps,nwm,nwa,nsa
c 1003 format(8l1,3i10)
c      elseif(card(6:6).eq.'T'.or.card(6:6).eq.'F') then
c      read(card,1023) lxw,lvw,lfw,lxs,lvs,lfs,nwm,nwa,nsa
c 1023 format(6l1,3i10)
c      lpw=.false.
c      lps=.false.
c      else
c      read(card,1033) lxw,lvw,lxs,lvs,nwm,nwa,nsa
c 1033 format(4l1,3i10)
911    continue       
c      if you reach this you are in trouble
       write(*,*) "found error"
       write(*,*) message
       stop

      end subroutine

      subroutine get_next_token(i0,buffer,sep,token,ostatus)
      implicit none
      integer ,         intent(inout) :: i0
      character(LEN=*), intent(in)    :: buffer
      character(LEN=*), intent(out)   :: token
      character*1,      intent(in)    :: sep
      logical,          intent(out)   :: ostatus
c      
      character*1 a1
      integer i,k,lb,lt
c      
      token = " "
      ostatus = .false.
      lb=len(buffer)
      lt=len(token)
      k=0
      do i=i0,lb
       a1=buffer(i:i)
       if(a1.eq." ") then
          if(k.gt.0.and.sep.eq." ") goto 30
       else if(a1.eq.sep.and.sep.ne." ") then
          goto 30
       else
         k = k+1
         if(k.gt.lt) goto 30
         token(k:k) = a1
       end if
      end do
30    continue
      i0=i+1
      ostatus = k.gt.0
      end subroutine
c $Id$
