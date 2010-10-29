      subroutine opprint(symops,rotoop,maxops,nops,itype)
C$Id$
      implicit none
      integer maxops, nops, itype
      character*2 rotoop(maxops)
      double precision symops(maxops*3,4)
      integer i, j, k, indx
c***********************************************************************
c     
c     print the matrix representations of the operators
c     
c***********************************************************************
        write(*,9)
        write(*,12)
        write(*,13) nops
        do 100 i=1,nops
          indx=(i-1)*3
c         
c         Add next block to avoid printing out bad rotoop data for pt. grps
c         
          if(itype.eq.0) then
            rotoop(i)='pt'
            write(*,8) rotoop(i)
          elseif(itype.eq.3) then
            write(*,8) rotoop(i)
          endif
          do 150 j=1,3
            write(*,10) (symops((indx+j),k), k=1,4)
  150     continue
  100   continue
    8 format(/,25x,a2,' fold Rotoinversion Operator')
    9 format(//,20x,'---------- FULL LISTING OF GROUP ----------')
   10 format(19x,4(f10.6))
   12 format(/,18x,'MATRIX REPRESENTATIONS OF THE GROUP OPERATORS',/)
   13 format(5x,'EXCLUDING THE IDENTITY OPERATOR THERE ARE',I3,' OPERATO
     &RS IN THIS GROUP')
      return
      end
