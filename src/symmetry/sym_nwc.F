c***********************************************************************
c  Subroutine sym_nwc: a 'friend' of geom
c
c Driver routine:Attributes
c
c   calls gensym chain:  (creates matrix representaions of the 
c                         symmetry operators)
c   calls dosymops:      (applies operatorst to symmetry distinct atoms
c                         prints and stores expanded lists)
c   This routine has direct access to geom private include files and
c   thus common blocks.
c
c   calls sym_map: creates atom pair mapping relations based upon the 
c                  symmetry operations and deposits symmetry operations
c                  and mapping tables to rtdb for future retrieval.
c                         
c
c                                          A.C.H. 5/12/94
c                 Interfaces to NWCHEM     A.C.H. 10/10/94 + RJH
c  
c***********************************************************************
      subroutine sym_nwc(geom,rtdb,nata,oprint)
      implicit real*8 (a-h,o-z)
      parameter(maxops=192)
      Parameter (EPS=1.D-14)

#include "mafdecls.h"
#include "geom.fh"
#include "geomP.fh"
#include "rtdb.fh"

      logical LResult, oprint
      integer geom,rtdb
      dimension cdist(3),cang(3)
      external pertbl
c
c--> parameters needed by sym chain
c
      itype  = isystype(geom)
      numgrp = group_number(geom)
      numset = setting_number(geom)
c
c--> lattice vectors & angles (cell constants)
c
      do 1000 i=1,3
        cdist(i)=lattice_vectors(i,geom)
        cang(i) =lattice_angles(i,geom)
 1000 continue
c
c--> check for system type (dimensionality) - leave dummy if structure.
c
      if (itype.eq.3) then
         call gensym(itype,numgrp,numset,sym_ops(1,1,geom),
     $       nops,oprint, group_name(geom))
      elseif(itype.eq.0) then
         call gensym(itype,numgrp,numset,sym_ops(1,1,geom),
     $       nops,oprint, group_name(geom))
      endif
c
c     scratch space needed by dosymops: tags, coordinates
c     space on stack: pointer list, points to the assymetric (distinct) atoms
c     expanded coord list
c
      LResult = MA_Push_Get(MT_Byte,nata*(nops+1)*16,'scratch labels',
     &     latags, iatags)

      LResult = MA_Push_Get(MT_Dbl,nata*(nops+1)*3,'scratch atoms',
     &     lcoord, icoord) .and. lresult

      Lresult = MA_Push_Get(MT_Int,nata,'distinct atoms',ldatom,idatom)
     $     .and. lresult
      if (.not. lresult) call errquit('sym_nwc: dosymops ma failed', 0)
c
      call dosymops(sym_ops(1,1,geom),nops,coords(1,1,geom),tags,nata,
     &     itype,Byte_MB(iatags),Dbl_MB(icoord),nareal,nata*(nops+1),
     $     cdist,cang,Int_MB(idatom), oprint)

c
      if(nareal.gt.max_cent) then
        write(iw,*)' Too many atoms, increase max_cent in ',
     &       'geomP.fh error in sym_nwc : ', nareal
        call errquit('sym: too many atoms', 1)
      endif
c
c     Allocate memory for the center map on the heap since it persists
c
      if (nops .gt. 0) then
      if (.not. MA_Alloc_Get(MT_Int,nops*nareal,'atom sym map',
     &     lmscr,imscr)) call errquit
     $     ('sym_nwc: ma_alloc_get of sym map failed', nops*nareal)
      else
       lmscr = -1
       imscr = 99999999
      endif
c
c--> derive mapping table for atom pairs permuted by sym operations
c--> NOTE: at this point coords will contain cartesian coordinates
c          and dbl_mb(icoord) has the fractional coordinates
c          (unless system is molecular in which case is also cartesian)
c
      if (nops .gt. 0) then
      if (.not. MA_Push_Get(MT_Dbl,nops*3,'xnew scratch',
     &     lxnew,ixnew)) call errquit
     $     ('sym_nwc: ma for xnew scratch', nops*3)
      endif
      if (.not. MA_Push_Get(MT_INT,nareal,'ilbl',
     &     llbl,ilbl)) call errquit
     $     ('sym_nwc: ma for ilbl', nareal)
c
      call sym_map(itype, Byte_MB(iatags), sym_ops(1,1,geom),
     $     nops, nareal, Int_Mb(imscr),Dbl_Mb(ixnew),Dbl_Mb(icoord),
     $     oprint,INT_MB(ilbl))
c
      LResult = MA_Pop_Stack(llbl)
      if (nops .gt. 0) then
      LResult = MA_Pop_Stack(lxnew)
      endif
c
c     Insert new information into the geometry data structures
c     Geom will hold the cartesian coords ... fractional ones
c     remain in icoord.
c     Note that the memory for the map (lmscr) is now freed in geom_destroy
c
      call sym_put_in_geom(geom, nata, nareal,
     $     Byte_MB(iatags), int_mb(idatom), dbl_mb(icoord),
     $     lmscr, imscr, nops)

c
c---> print gmat,amat,vol evaluated by geom_nd routines
c
      if(isystype(geom).ne.0 .and. oprint) then
c
c--> Print A matrix.
c
         write(*,9)
         do 513 i=1,3
            write(*,15) (amatrix(i,j,geom), j=1,3)
 513     continue
c     
c--> Print G matrix.
c
        write(*,10)
        do 130 i=1,3
          write(*,15) (metric_matrix(i,j,geom), j=1,3)
  130   continue
c     
c--> Print volume of unit cell.
c     
        write(*,17) volume_direct(geom)
        write(*,14)
c
c--> Print cartesian coordinates of unit cell.
c
        write(*,30)
        l=0
        do 510 i=1,nareal
          write(*,20) i,tags(i,geom),(coords(j,i,geom), j=1,3)
          l=l+3
  510   continue
c
      endif
c
c--> final clean up of memory: gensym and dosymops activities
c
      LResult = MA_Pop_Stack(ldatom)
      LResult = MA_Pop_Stack(lcoord)
      LResult = MA_Pop_Stack(latags)
c
c--> formats
c
    9 format(//15x,
     & 'Crystallographic to Cartesian Transformation Matrix',/)
   10 format(//,21x,'THE METRIC MATRIX FOR THE DIRECT LATTICE ',/)
   15 format(20x,3(f12.6,2x))      
   17 format(/,19x,'UNIT CELL VOLUME = ',f12.6,3x,'ANGSTOMS**3')
   14 format(/,80('-'),/)
c   20 format(15x,i5,5x,a16,3(f10.6,2x))
   20 format(10x,i5,2x,a16,3(f10.6,2x))
   30 format(/,15x,'--------- CARTESIAN COORDINATES (A.U.) ---------',
     &/)
      return
      end
