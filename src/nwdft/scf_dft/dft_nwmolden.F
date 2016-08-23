c
c     $Id: dft_scf.F 28430 2016-07-27 18:47:23Z edo $
c
c     interface to nwmolden
c     need to print un-normalized basis exp to avoid troubles
c     Before using JANPA,  run the .molden file through 
c     molden2molden.jar -NormalizeBF to get the correct Molden normalization
c     https://sourceforge.net/p/janpa/wiki/nwchem2molden/
      subroutine dft_nwmolden(rtdb,geom,
     I     ipol,eval,occ,ir,
     N     nbf, nmo,
     I     g_movecs,scftype)
      implicit none
#include "errquit.fh"
#include "bas.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
      integer rtdb
      integer geom 
      integer ipol
      double precision eval(*),occ(*)
      integer ir(*)
      integer nbf, nmo
      integer g_movecs(*)
      character*(*) scftype
c
      character*20 pname
      integer ispin
      integer bas
      logical nwmolden_movecs,nwmolden_basis,nwmolden_geom
      external nwmolden_movecs,nwmolden_basis,nwmolden_geom
      logical nwmolden_openfile,nwmolden_closefile
      external nwmolden_openfile,nwmolden_closefile
      logical do_io
c
c
      do_io=ga_nodeid().eq.0
      pname='dft_nwmolden: '
      if(do_io) then
         if(.not.nwmolden_openfile())
     E        call errquit(pname//'nwmolden_openfile failed ',0,0)
      if(.not.nwmolden_geom(geom))
     E     call errquit(pname//'nwmolden_geom failed ',0,0)
      endif
c     create internal basis handle
      if (.not. bas_create(bas, 'ao basis'))
     &   call errquit(pname//'bas_create?', 2, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, bas, 'ao basis'))
     &   call errquit(pname//'no ao basis ', 3, BASIS_ERR)
      if(do_io) then
         if(.not.nwmolden_basis(bas))
     E        call errquit(pname//'nwmolden_basis failed ',0,0)
c
         do ispin=1,ipol
            if(.not.nwmolden_movecs(bas,nbf,nmo,
     1           eval((ispin-1)*nbf+1),occ((ispin-1)*nbf+1),
     S           ir((ispin-1)*nbf+1),
     G           g_movecs(ispin),scftype))
     E           call errquit(pname//'nwmolden_movecs failed ',0,0)
         enddo
      endif
c
         if (.not. bas_destroy(bas))
     &      call errquit(pname//'not able to destroy CD_bas:',33,
     &       BASIS_ERR)
         if(do_io) then
         if(.not.nwmolden_closefile())
     E        call errquit(pname//'nwmolden_closefile failed ',0,0)
         write(6,*) ' '
         write(6,*) ' Molden file created '
         write(6,*) ' If you plan to use JANPA '
         write(6,*) '   run the .molden file through '
         write(6,*) '   molden2molden.jar -NormalizeBF to get the'
         write(6,*) '   correct Molden normalization '
         write(6,*) '  '
         endif
         call ga_sync()
         return
         end
