      subroutine sym_map(itype, tags, symops, nops, nareal, map, xnew,
     &     acord, oprint, ilbl, threquiv)
C$Id$
      implicit none
#include "errquit.fh"
c
      integer itype
      integer nareal
      character*16 tags(nareal)
      integer nops
      integer maxops
      double precision tol
      double precision threquiv
      parameter(maxops=192)
      double precision symops(maxops*3,4),xnew(nops,3),acord(3,nareal)
c     
      integer map(nops,nareal)
      integer ilbl(nareal)
      logical oprint
      integer i, j, k, ii, iop, iopnum, num
      double precision sum, dx, dy, dz
c     
c     creates an equiv. table for mapping of atoms under sym ops
c     stores the matrix representations of the group operators, number
c     of operators and equiv.
c     
c     table format:
c     atm #1    atm #2  ..........   atm #N
c     op#1
c     
c     op#2
c     .
c     .
c     .
c     
c     op#N   
c     
c     therefore, searching down a col. will give you all the atoms equivalent
c     to atom N under all operations of the group.
c     
c     ach 10/14/94
c     
c     Initialize map to invalid center number so that can check
c     all entries are set at the end
c     
      do i = 1, nareal
         do iop = 1, nops
            map(iop,i) = 0
         enddo
      enddo
c
c     set tolerance
c
      tol=threquiv
c     
c     symops*coordinates
c     
      sum=0.0d+00
      do 300 i=1,nareal
         iopnum=0
         do 305 iop=1,nops
            do 310 j=1,3
               iopnum=iopnum+1
               do 320 k=1,3
                  sum=sum+symops(iopnum,k)*acord(k,i)
 320           continue
               xnew(iop,j)=sum+symops(iopnum,4)
               sum=0.0d+00
 310        continue
c     
c--   > shift it into the home cell. This operation depends on sysytem type
c--   > in place for 3-d and molecules now.
c     
            if(itype.eq.3) then
 3000          do 380 ii=1,3
                  if(dabs(xnew(iop,ii)).lt.1.0d-10) then
                     xnew(iop,ii)=0.0d0
                  endif
                  if(dabs((dabs(xnew(iop,ii))-1.0d0)).lt.1.0d-10) then
                     xnew(iop,ii)=0.0d0
                  endif
                  if (xnew(iop,ii).lt.(0.0d0)) then
                      xnew(iop,ii)=xnew(iop,ii)+1.0d0
                      goto 3000
                  elseif (xnew(iop,ii).ge.1.0d0) then
                     xnew(iop,ii)=xnew(iop,ii)-1.0d0
                     goto 3000
                  endif
 380           continue
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
c--   > load tags equiv. to the tag of atom i for search
         num=0
         do 400 j=1,nareal
            if(tags(i).eq.tags(j)) then
               num=num+1
               ilbl(num)=j
            endif
 400     continue
c--   > search those tags to identify atoms in xnew
         do 410 j=1,nops
            do 420 k=1,num
               dx=dabs(acord(1,ilbl(k))-xnew(j,1))
               if(dx.lt.tol) then
                  dy=dabs(acord(2,ilbl(k))-xnew(j,2))
                  if(dy.lt.tol) then
                     dz=dabs(acord(3,ilbl(k))-xnew(j,3))
                     if(dz.lt.tol) then
                        map(j,i)=ilbl(k)
                        goto 421
                     endif
                  endif
               endif
 420        continue
            write(6,*) ' atom = ', i
            write(6,*) ' op   = ', j
            write(6,*) ' xnew   ', (xnew(j,k),k=1,3)
            call util_flush(6)
            call errquit('sym_map: no match for 1000*op+atom ',1000*j+i,
     &       INPUT_ERR)
 421        continue
 410     continue
 300  continue
c
c     Check the map has all entries in it.
c
      do i = 1, nareal
         do j = 1, nops
            if (map(j,i).le.0 .or. map(j,i).gt.nareal) call errquit
     $           ('sym_map: invalid map for 1000*op+atom ',1000*j+i,
     &       INPUT_ERR)
         enddo
      enddo
c     
c--   > test prints (prints map)
c     
c     write(*,*) 'map printed out nops rows nareal columns'
c     do 500 i=1,nops
c     write(*,'(8(i2,3x))') (map(i,j),j=1,nareal)
c     500 continue
c     
      if (oprint) write(*,10)
 10   format(//,18x,'<<< SYMMETRY OPERATION/ATOM MAPPING BUILT >>>')
      return
      end




