      subroutine selci_output(z,rowlow,rowhi,collow,colhi,rowdim,coldim,
     $     nctl)
*
* $Id: selci_output.f,v 1.1 2003-04-07 23:48:56 windus Exp $
*
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
      implicit double precision (a-h,o-z)
      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
      dimension z(rowdim,coldim)
      dimension asa(3)
      data column/8hcolumn   /,asa/8h          ,8h00000000  ,
     1     8h--------  /,blank/8h          /
      data kcol/4/
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
      last = min0(colhi,collow+kcol-1)
      do 2 begin = collow,colhi,kcol
         write (6,1000) (column,i,i = begin,last)
         do 1 k = rowlow,rowhi
            do 4 i=begin,last
               if (z(k,i).ne.zero) go to 5
 4          continue
            go to 1
 5          write (6,2000) ctl,k,(z(k,i), i = begin,last)
 1       continue
         last = min0(last+kcol,colhi)
 2    continue
 3    return
 1000 format (/1h ,16x,3(a6,i4,7x),(a6,i4))
 2000 format (a1,3hrow,i4,2x,4f17.11)
      end
