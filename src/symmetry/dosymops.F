ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c   subroutine dosymops applies symops to the assymetric unit atoms
c      input to it, outputs a list of unique atoms defining the 
c      molecule, surface or crystal. it also sorts them by name
c
c                                       3/1/94 A.C. Hess
c
c
c   some variables:
c       coords = input coordinates of atoms in assymetric unit
c        nata = number of atoms in assymtric unit
c       coord = full coordinate list of all atoms
c         nat = total number of atoms
c      ipcord = pointer before the first empty spot in coord
c        xnew = array holds atoms generated for a given assym atoms 
c               for the list of entire list of symops, 4th colum is
c               used in an in place comparison of atoms to determine
c               uniqueness
c iasym_order = keeps track of where the original assymetric unit
c               atoms are in the final (full) coordinate list.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dosymops(symops,nops,coords,tags,nata,itype,
     &     tags_new,coords_new,nareal,nat,cdist,cang,iasym_order,
     $     charge, charge_new, mass, mass_new, 
     $     invnucexp, invnucexp_new, threquiv)
C$Id$
*rak typed everyting
      implicit none
#include "util.fh"
c
c has to be a .F file since print_debug is in the include definition and
c with the implicit real*8 it was of type real not integer
#include "errquit.fh"
c
      double precision tol, threquiv
      integer maxops, iodat
      parameter (maxops=192,iodat=11)
c
      integer nata, nat
      integer nops, itype, nareal
      character*16 tags(nata),tags_new(nat)
      double precision symops(maxops*3,4),coords_new(3,nat),
     &    xnew(maxops,4),
     &    coords(3,nata),cdist(3),cang(3)
      integer iasym_order(nata)
      logical oprint
      double precision charge(*), charge_new(*), mass(*),
     $     mass_new(*), invnucexp(*), invnucexp_new(*)
c
      integer i, j, k, l, ii, ll, nata_old,
     &    ipcord, ipt, nauniq, iopnum, im, iop
      double precision sum, dx, dy, dz, x, y, z 
c
c     set tolerance
c
      tol = threquiv
      nata_old = nata
c
c     Modified below (RJH) so that can send in a list with symmetry
c     equivalent atoms and not generate redundant ones and also
c     identify the unique ones correctly
c
c--> application of symops to atoms in assym. unit
c-->
c--> outer loops go: 1)over all atoms, 2) over the number of symops
c--> 3) mult. current symop times current atom coordinates. Done for
c--> all symops for a given atom, 4) the list is then sorted, and the
c--> unique atoms added to the master list.
c
c--> if its a fractional system move all atoms into cell before application
c--> of symops. still no 1-d systems
c
c bulk
      oprint = util_print('dosymops',print_debug)
c
      if(itype.eq.3) then
        do 900 i=1,nata
          do 910 j=1,3
 2000       if(coords(j,i).lt.(0.0d+00)) then
              coords(j,i)=coords(j,i)+1.0d+00
              goto 2000
            elseif (coords(j,i).ge.1.0d+00) then
              coords(j,i)=coords(j,i)-1.0d+00
              goto 2000
            endif
            if(dabs(coords(j,i)) .lt. 1.0d-10) then
              coords(j,i) = 0.0d+00
            endif
            if(dabs(dabs(coords(j,i))-1.0d+00) .lt. 1.0d-10) then
              coords(j,i) = 0.0d+00
            endif
  910     continue
  900   continue
c surfaces
      elseif(itype.eq.2) then
         do 902 i=1,nata
            do 912 j=1,2
 2002          if(coords(j,i).lt.0.0d+00) then
                  coords(j,i)=coords(j,i)+1.0d+00
                  goto 2002
               elseif (coords(j,i).ge.1.0d+00) then
                  coords(j,i)=coords(j,i)-1.0d+00
                  goto 2002
               endif
               if(dabs(coords(j,i)) .lt. 1.0d-10) then
                 coords(j,i) = 0.0d+00
               endif
               if(dabs(dabs(coords(j,i))-1.0d+00) .lt. 1.0d-10) then
                 coords(j,i) = 0.0d+00
               endif
 912        continue
 902     continue
c         
      elseif(itype.eq.1) then
         do 904 i=1,nata
            j = 3
 2004          if(coords(j,i).lt.0.0d+00) then
                  coords(j,i)=coords(j,i)+1.0d+00
                  goto 2004
               elseif (coords(j,i).ge.1.0d+00) then
                  coords(j,i)=coords(j,i)-1.0d+00
                  goto 2004
               endif
               if(dabs(coords(j,i)) .lt. 1.0d-10) then
                 coords(j,i) = 0.0d+00
               endif
               if(dabs(dabs(coords(j,i))-1.0d+00) .lt. 1.0d-10) then
                 coords(j,i) = 0.0d+00
               endif
 904     continue
c         
      endif
c
c symops*coordinates
c
      ipcord=0
      nareal=0
      sum=0.0d+00
      ipt=nops+1
      nauniq = 0
c
      do 300 i=1,nata
c
c     Check to see if the atom is already present
c
         do l=1,ipcord
            dx=dabs(coords(1,i)-coords_new(1,l))
            dy=dabs(coords(2,i)-coords_new(2,l))
            dz=dabs(coords(3,i)-coords_new(3,l))
            if(dx+dy+dz.le.3*tol) then
c              If the tags don't match, we are probably in trouble
               if(tags(i) .eq. tags_new(l)) then
                 if (oprint)write(6,*) ' ATOM ',i,' is already present '
                 goto 300
               endif
            endif
         enddo
         if (oprint) write(6,544) 
     $        ' doing atom ',i,coords(1,i),coords(2,i),coords(3,i)
 543     format(a,2i4,3f12.6)
 544     format(a,i4,3f12.6)
         iopnum=0
c
c--> reset flag for atom uniqueness
c
         do 570 im=1,nops
            xnew(im,4)=0.0
570      continue
         do 305 iop=1,nops
            do 310 j=1,3
               iopnum=iopnum+1
               do 320 k=1,3
                  sum=sum+symops(iopnum,k)*coords(k,i)
320            continue
               xnew(iop,j)=sum+symops(iopnum,4)
               sum=0.0d+00
310         continue
            if (oprint) write(6,544) 
     $           ' op ', iop,xnew(iop,1),xnew(iop,2),xnew(iop,3)
c
c--> shift it into the home cell. This operation depends on sysytem type
c--> in place for 3-d,2-d and molecules now.
c
            if(itype.eq.3) then ! bulk
 3000         do 380 ii=1,3
                if(dabs(xnew(iop,ii)).lt.1.0d-10) then
                  xnew(iop,ii)=0.0d0
                endif
                if(dabs((dabs(xnew(iop,ii))-1.0d0)).lt.1.0d-10) then
                   xnew(iop,ii)=0.0d0
                endif
                if(xnew(iop,ii).lt.(0.0d0)) then
                  xnew(iop,ii)=xnew(iop,ii)+1.0d0
                  goto 3000
                elseif (xnew(iop,ii).ge.1.0d0) then
                  xnew(iop,ii)=xnew(iop,ii)-1.0d0
                  goto 3000
                endif
  380         continue
            elseif(itype.eq.2) then ! surfaces
 3020         do 382 ii=1,2
                if(dabs(xnew(iop,ii)).lt.1.0d-10) then
                  xnew(iop,ii)=0.0d0
                endif
                if(dabs((dabs(xnew(iop,ii))-1.0d0)).lt.1.0d-10) then
                   xnew(iop,ii)=0.0d0
                endif
                if(xnew(iop,ii).lt.0.0d0) then
                  xnew(iop,ii)=xnew(iop,ii)+1.0d0
                  goto 3020
                elseif (xnew(iop,ii).ge.1.0d0) then
                  xnew(iop,ii)=xnew(iop,ii)-1.0d0
                  goto 3020
                endif
 382         continue
            endif
  305     continue
c
c--> place this assym. atom into final coord list
c
         ipcord=ipcord + 1
         nauniq = nauniq + 1
         nareal = nareal + 1
         do 330 ll=1,3
            coords_new(ll,ipcord)=coords(ll,i)
330      continue
         tags_new(ipcord)=tags(i)
         charge_new(ipcord)=charge(i)
         invnucexp_new(ipcord)=invnucexp(i)
         mass_new(ipcord)=mass(i)
         iasym_order(nauniq)=ipcord
*
*rak: changed K -> i  was that right?
*
         if (oprint) write(6,543) ' keepingX ', ipcord,nareal,
     $        xnew(i,1),xnew(i,2),xnew(i,3)
c
c--> remove copies of the assymetric atom from the xnew list
c
         x=coords(1,i)
         y=coords(2,i)
         z=coords(3,i)
         if (oprint) write(6,544) ' comparing ', i, x,y,x
         do 390 l=1,nops
               dx=dabs(x-xnew(l,1))
               dy=dabs(y-xnew(l,2))
               dz=dabs(z-xnew(l,3))
               if (oprint) write(6,544)
     $              '    with ',l,xnew(l,1),xnew(l,2),xnew(l,3)
               if (oprint) write(6,544)
     $              '    diff ',l,dx,dy,dz
               if(dx.le.tol.and.dy.le.tol.and.dz.le.tol) then
                  xnew(l,4)=1.0d0
               endif
390      continue
c
c--> sift the xnew list for unique atoms
c
         do 400 k=1,nops
            if(xnew(k,4).ne.1.0d0) then
               x=xnew(k,1)
               y=xnew(k,2)
               z=xnew(k,3)
               do 410 l=k+1,nops
                  dx=dabs(x-xnew(l,1))
                  dy=dabs(y-xnew(l,2))
                  dz=dabs(z-xnew(l,3))
                  if(dx.le.tol.and.dy.le.tol.and.dz.le.tol) then
                     if (oprint) write(6,*) ' KEEP IT',l
                     xnew(l,4)=1.0d0
                  endif
410            continue
            endif
400      continue
c
c--> put the new unique atoms into the final list, keeping track of 
c--> what kind they are.
c
         do 420 k=1,nops
            if(xnew(k,4).ne.1.0d0) then
               ipcord=ipcord+1
               nareal=nareal+1
               do 430 ll=1,3
                  coords_new(ll,ipcord)=xnew(k,ll)
430            continue
               if (oprint) write(6,543) ' keeping ', ipcord,nareal,
     $              xnew(k,1),xnew(k,2),xnew(k,3)
               tags_new(ipcord)=tags(i)
               charge_new(ipcord)=charge(i)
               invnucexp_new(ipcord)=invnucexp(i)
               mass_new(ipcord)=mass(i)
            endif
420      continue
300   continue
c
      nata = nauniq
      if (oprint) then 
         write(*,25) nareal
         do 510 i=1,nareal
            write(*,2122) i,tags_new(i),(coords_new(j,i),j=1,3),
     $           charge_new(i)
 510     continue
 2122    format(i5,2x,a16,2x,3f12.6,2x,f10.2)
      endif
*
*     RJH.  If the user input all of the atoms don't want
*     to reorder them since it may break zmatrix/autoz
*
      if (itype.eq.0 .and. nareal.eq.nata_old) then
c
c     Map unique atoms back into old ordering
c
         do ii = 1, nauniq
            i = iasym_order(ii)
            do j = 1, nareal
               dx = abs(coords_new(1,i) - coords(1,j))
               dy = abs(coords_new(2,i) - coords(2,j))
               dz = abs(coords_new(3,i) - coords(3,j))
               if(dx.le.tol.and.dy.le.tol.and.dz.le.tol) then
                  iasym_order(ii) = j
                  goto 9876
               endif
            enddo
            call errquit('dosymops: failed matching uniq',i, INPUT_ERR)
 9876       continue
         enddo
         do i = 1, nareal
            tags_new(i) = tags(i)
         enddo
         call dcopy(nareal, charge, 1, charge_new, 1)
         call dcopy(nareal, invnucexp, 1, invnucexp_new, 1)
         call dcopy(nareal, mass, 1, mass_new, 1)
         call dcopy(nareal*3, coords, 1, coords_new, 1)
      endif
c     
c
c--> write stuff out
c
c--> coordinate output stuff
c
c     if(itype.ne.0) then
c        write(*,13)
c        if(itype.eq.3) then
c           write(*,28) 
c           write(*,29)
c        elseif(itype.eq.2) then
c           write(*,32)
c           write(*,31)
c        elseif(itype.eq.1) then
c           write(*,33)
c           write(*,34)
c        endif
c        write(*,14)
c        call wrcell(cdist,cang,itype)
c     
c--   > list of atoms in assymetric unit
c     
c        write(*,10) nata
c        do 500 i=1,nata
c           write(*,20) i,tags(i),(coords(j,i),j=1,3)
c500     continue
c
c--> list of total number of atoms unit cell of the solid
c
c        write(*,25) nareal
c        do 510 i=1,nareal
c           write(*,20) i,tags_new(i),(coord(j,i),j=1,3)
c510     continue
c     endif
c
10    format(/,22x,i5,6x,'SYMMETRY DISTINCT ATOMS',/)
13    format(///,'--------------------- PERIODIC SYSTEM GEOMETRY INFORMA
     $TION',' -------------------',/)
 14   format(/,'------------------------------------------------------',
     &'-----------------------',/)
 20   format(10x,i5,2x,a16,3f14.8)
 25   format(//,22x,i5,6x,'ATOMS IN UNIT CELL',/)
 28   format(28X,'SYSTEM TYPE = BULK CRYSTAL',/)
 29   format(25x,' FRACTIONAL ATOMIC COORDINATES ')
 30   format(/,15x,'-------------- ATOMIC COORDINATES --------------',/)
 31   format(6x,'FRACTIONAL & CARTESIAN COORDINATES: Z-COMPONENT UNITS =
     $ANGSTROMS')
 32   format(29x,'SYSTEM TYPE = SURFACE',/)
 33   format(29x,'SYSTEM TYPE = POLYMER',/)
 34   format(5x,'FRACTIONAL & CARTESIAN COORDINATES: XY-COMPONENTS IN UN
     $ITS = ANGSTROMS')
c
      if (oprint) call util_flush(6)
      end
