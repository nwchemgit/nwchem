      program read_genr
      implicit none
      integer max_rank
      parameter (max_rank=1000)
      double precision grad0(max_rank), hess(max_rank*max_rank)
      integer iatom, rank, ixyz, iflag0
      integer nat
      open(unit=69,file='gen_at_restart',
     &    form='unformatted',
     &    access='sequential',
     &    status='old',
     &    err=09990)
      read(69) iatom, ixyz, rank, iflag0
c
c     compute nat <=> number of atoms
c
      if (mod(rank,3).ne.0) then
        write(6,*)' rank is not a multiple of 3?'
        stop ' rank error '
      endif
      nat = rank/3
      if (rank.gt.max_rank) then
        write(6,*)' rank     = ',rank
        write(6,*)' max_rank = ',max_rank
        write(6,*)' increase max_rank in source and recompile/link'
        stop ' max_rank error '
      endif
      write(6,*)' restart file completed atom  :', iatom
      write(6,*)' restart file completed xyz   :', ixyz
      write(6,*)' restart file rank of hessian :', rank
      if (iflag0.eq.1) then
        write(6,*)' restart file does have grad0 :',' true'
      else
        write(6,*)' restart file does have grad0 :',' false'
      endif 
      write(6,*)' number of atoms              :',nat
      if (iflag0.eq.1) call read_prt_grad0(69,grad0,rank,nat)
      call read_prt_hess(69,hess,rank,nat)
      stop ' done '
09990 write(6,*)' error opening restart file '
      stop '09990'
      end
      subroutine read_prt_hess(iu,hess,rank,nat)
      implicit none
c
      integer iu
      integer rank, nat
      double precision hess(rank,rank)
      integer icnt, jcnt, iatom, ixyz, jatom, jxyz
      character*1 labxyz(3)
      data labxyz /'X','Y','Z'/
c
      read(iu) hess
      write(6,*)' nonzero hessian values '
      jcnt = 0
      do 00100 jatom=1,nat
        do 00200 jxyz = 1,3
          jcnt = jcnt + 1
          icnt = 0
          do 00300 iatom=1,nat
            do 00400 ixyz = 1,3
              icnt = icnt + 1
              if (abs(hess(icnt,jcnt)).gt.1.0d-07) then
                write(6,10000)iatom,labxyz(ixyz),
     &                        jatom,labxyz(jxyz),
     &                        hess(icnt,jcnt)
              endif
00400       continue
00300     continue
00200   continue
00100 continue
      write(6,*)' '
      write(6,*)' hessian matrix '
      write(6,*)' '
      call output(hess,1,rank,1,rank,rank,rank,1)
10000 format('atom',i4,1x,a1,1x,'atom',i4,1x,a1,1x,f10.6)
      end
      subroutine read_prt_grad0(iu,grad0,rank,nat)
      implicit none
c
      integer iu
      integer rank, nat
      double precision grad0(rank)
      integer icnt, iatom, ixyz
      character*1 labxyz(3)
      data labxyz /'X','Y','Z'/
c
      read(iu) grad0        
      write(6,*)' nonzero grad0 values '
      icnt = 0
      do 00100 iatom = 1,nat
        do 00200 ixyz = 1,3
          icnt = icnt + 1
          if (abs(grad0(icnt)).gt.1.0d-07) then
            write(6,10000)iatom,labxyz(ixyz),grad0(icnt)
          endif
00200   continue
00100 continue
10000 format(' atom ',i4,1x,a1,1x,f10.6)
      end
      subroutine output (z,rowlow,rowhi,collow,colhi,rowdim,coldim,
     $     nctl)
c.......................................................................
c output prints a real*8 matrix in formatted form with numbered rows
c and columns.  the input is as follows;
c        matrix(*,*).........matrix to be output
c        rowlow..............row number at which output is to begin
c        rowhi...............row number at which output is to end
c        collow..............column number at which output is to begin
c        colhi...............column number at which output is to end
c        rowdim..............row dimension of matrix(*,*)
c        coldim..............column dimension of matrix(*,*)
c        nctl................carriage control flag; 1 for single space
c                                                   2 for double space
c                                                   3 for triple space
c the parameters that follow matrix are all of type integer*4.  the
c program is set up to handle 5 columns/page with a 1p5d24.15 format for
c the columns.  if a different number of columns is required, change
c formats 1000 and 2000, and initialize kcol with the new number of
c columns.
c author;  nelson h.f. beebe, quantum theory project, university of
c          florida, gainesville
c.......................................................................
C$Id$
      implicit none
      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
      integer nctl, i, j, last, k
      double precision z(rowdim,coldim), zero
      character*8 asa(3), column, ctl, blank
      data column/'column'   /,asa/' ','00000000'  ,
     1     '--------'  /,blank/' '/
      data kcol/8/
      data zero/0.d00/
      do 11 i=rowlow,rowhi
         do 10 j=collow,colhi
            if (z(i,j).ne.zero) go to 15
 10      continue
 11   continue
      write (6,3000)
 3000 format (/' zero matrix'/)
      go to 3
 15   continue
      ctl = blank
      if ((nctl.le.3).and.(nctl.gt.0)) ctl = asa(nctl)
      if (rowhi.lt.rowlow) go to 3
      if (colhi.lt.collow) go to 3
      last = min(colhi,collow+kcol-1)
      do 2 begin = collow,colhi,kcol
*         write (6,1000) (column,i,i = begin,last)
         write (6,1000) (i,i = begin,last)
         do 1 k = rowlow,rowhi
            do 4 i=begin,last
               if (z(k,i).ne.zero) go to 5
 4          continue
            go to 1
 5          write (6,2000) ctl,k,(z(k,i), i = begin,last)
 1       continue
         last = min(last+kcol,colhi)
 2    continue
 3    return
* kcol = 4
* 1000 format (/1h ,16x,3(a6,i3,2x),(a6,i3))
* 2000 format (a1,3hrow,i4,2x,4f17.11)
* kcol = 8
*
* if U like having rows and columns labelled with row and col
* use these
*
* 1000 format (/1h ,11x,7(a3,i3,3x),(a3,i3))
* 2000 format (a1,'row',i4,1x,8f9.4)
c
 1000 format (/1h ,8x,7('   ',i3,3x),('   ',i3))
 2000 format (a1,i4,1x,8f9.4)
      end
