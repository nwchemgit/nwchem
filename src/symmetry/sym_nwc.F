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
c   call sym_init_inv_op: to make list of inverse operations
c
c                                          A.C.H. 5/12/94
c                 Interfaces to NWCHEM     A.C.H. 10/10/94 + RJH
c  
c***********************************************************************
      subroutine sym_nwc(geom,rtdb,nata,oprint,scale,threquiv,nops)
C$Id$
      implicit real*8 (a-h,o-z)
#include "errquit.fh"
      parameter(maxops=192)
      Parameter (EPS=1.D-14)

#include "mafdecls.fh"
#include "geom.fh"
#include "nwc_const.fh"
#include "geomP.fh"
#include "rtdb.fh"

      logical LResult, oprint
      integer geom,rtdb
      dimension cdist(3),cang(3)
      double precision mass_new(max_cent), charge_new(max_cent)
      double precision invnucexp_new(max_cent)
C      external pertbl
c
c--> parameters needed by sym chain
c
      itype  = isystype(geom)
      numgrp = group_number(geom)
      numset = setting_number(geom)
c
c   lattice vectors & angles (cell constants)
c   values should have been read in and initialized properly for
c   each system dimension
c
      do 1000 i=1,3
        cdist(i)=lattice_vectors(i,geom)
        cang(i) =lattice_angles(i,geom)
 1000 continue
c
c    generate symmetry operators for whatever group (crystals, surfaces
c    polymers or molecules) have been requested.
c
      if ((itype.eq.3).and.(numset.gt.2)) then
          !write(*,*) "generating unsual groups"
          call gensym_extra(itype,group_number(geom),numset,
     $                      sym_ops(1,1,geom),
     $                      nops,oprint, group_name(geom))
      else
         call gensym(itype,numgrp,numset,sym_ops(1,1,geom),
     $               nops,oprint, group_name(geom),geom,rtdb)
      end if

c
c     scratch space needed by dosymops: tags, coordinates
c     space on stack: pointer list, points to the assymetric (distinct) atoms
c     apply symmops to coordinate list
c
      LResult = MA_Push_Get(MT_Byte,nata*(nops+1)*16,'scratch labels',
     &     latags, iatags)
      LResult = MA_Push_Get(MT_Dbl,nata*(nops+1)*3,'scratch atoms',
     &     lcoord, icoord) .and. lresult
      Lresult = MA_Push_Get(MT_Int,nata,'distinct atoms',ldatom,idatom)
     $     .and. lresult
      if (.not. lresult) call errquit('sym_nwc: dosymops ma failed', 0,
     &       MA_ERR)
c
      call dosymops(sym_ops(1,1,geom),nops,coords(1,1,geom),
     $     tags(1,geom),nata,
     &     itype,Byte_MB(iatags),Dbl_MB(icoord),nareal,nata*(nops+1),
     $     cdist,cang,Int_MB(idatom), 
     $     charge(1,geom), charge_new, geom_mass(1,geom),
     $     mass_new, geom_invnucexp(1,geom),invnucexp_new,threquiv)
c
      if(nareal.gt.max_cent) then
        write(6,*)' Too many atoms, increase nw_max_atom in ',
     &       'util/nwc_const.fh error in sym_nwc : ', nareal
        call errquit('sym: too many atoms', 1, INPUT_ERR)
      endif
c
c     Allocate memory for the center map on the heap since it persists
c
      if (nops .gt. 0) then
      if (.not. MA_Alloc_Get(MT_Int,nops*nareal,'atom sym map',
     &     lmscr,imscr)) call errquit
     $     ('sym_nwc: ma_alloc_get of sym map failed', nops*nareal,
     &       MA_ERR)
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
     &       lxnew,ixnew)) call errquit
     $       ('sym_nwc: ma for xnew scratch', nops*3, MA_ERR)
      else
        ixnew = 0
      endif
      if (.not. MA_Push_Get(MT_INT,nareal,'ilbl',
     &     llbl,ilbl)) call errquit
     $     ('sym_nwc: ma for ilbl', nareal, MA_ERR)
c
      call sym_map(itype, Byte_MB(iatags), sym_ops(1,1,geom),
     $     nops, nareal, Int_Mb(imscr),Dbl_Mb(ixnew),Dbl_Mb(icoord),
     $     oprint,INT_MB(ilbl), threquiv)
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
     $     lmscr, imscr, nops, charge_new, mass_new, invnucexp_new)
c
c ... sloppy way for the moment to fix units for 2-d case
c
      if(itype.eq.2) then
         do i=1,nareal
            coords(3,i,geom)=coords(3,i,geom)*scale
         enddo
      endif
c
c---> print gmat,amat,vol evaluated by geom_nd routines
c
c
c---> print gmat,amat,vol evaluated by geom_nd routines
c
c
c--> Print Amatrix, G matrix, volume for 3-d and 2-d cells.
c
c      if(itype.eq.3) then
c         write(*,9)
c         do 513 i=1,3
c            write(*,15) (amatrix(i,j,geom), j=1,3)
c 513     continue
c
c         write(*,21)
c         do 517 i=1,3
c            write(*,15) (amatrix_inv(i,j,geom), j=1,3)
c 517     continue
c
c         write(*,10)
c         do 130 i=1,3
c            write(*,15) (metric_matrix(i,j,geom), j=1,3)
c 130     continue
c
c         write(*,17) volume_direct(geom)
c
c      elseif(itype.eq.2) then
c         write(*,9)
c         do 630 i=1,3
c            write(*,15) (amatrix(i,j,geom), j=1,3)
c 630     continue
c         write(*,21)
c         do 631 i=1,3
c            write(*,15) (amatrix_inv(i,j,geom), j=1,3)
c 631     continue
c     
c         write(*,10)
c         do 640 i=1,3
c            write(*,15) (metric_matrix(i,j,geom), j=1,3)
c 640     continue
c     
c         write(*,18) volume_direct(geom)
c      endif
c
c--> Print cartesian coordinates of unit cell.
c
c      if(itype.eq.3) then
c        write(*,30)
c        l=0
c        do 510 i=1,nareal
c           write(*,20) i,tags(i,geom),(coords(j,i,geom), j=1,3)
c           l=l+3
c 510    continue
c         write(*,14)
c   convert z-coords of slab to bohr from angstroms         
c      elseif(itype.eq.2) then
c         do i=1,nareal
c            coords(3,i,geom)=coords(3,i,geom)*angstrom_to_au
c         enddo
c         write(*,30)
c         l=0
c         do 512 i=1,nareal
c            write(*,20) i,tags(i,geom),(coords(j,i,geom), j=1,3)
c            l=l+3
c 512     continue
c         write(*,14)
c      endif
c
c--> final clean up of memory: gensym and dosymops activities
c
      LResult = MA_Pop_Stack(ldatom)
      LResult = MA_Pop_Stack(lcoord)
      LResult = MA_Pop_Stack(latags)
c
      call sym_init_inv_op(geom)
c
c--> formats
c
 9    format(//20x,
     &     'CRYSTALLOGRAPHIC TO CARTESIAN MATRIX (A.U.)',/)
 10   format(//,23x,'METRIC MATRIX:DIRECT LATTICE (A.U.)',/)
 14   format(/,80('-'),//)
 15   format(20x,3(f12.6,2x))      
 17   format(/,19x,'UNIT CELL VOLUME = ',f12.6,3x,'BOHR**3')
 18   format(/,19x,'UNIT CELL AREA = ',f12.6,3x,'BOHR**3')
 20   format(10x,i5,2x,a16,3(f10.6,2x))
 21   format(//20x,
     & 'CARTESIAN TO CRYSTALLOGRAPHIC MATRIX (A.U.)',/)
 30   format(/,15x,'--------- CARTESIAN COORDINATES (A.U.) ---------',
     &/)
      return
      end



