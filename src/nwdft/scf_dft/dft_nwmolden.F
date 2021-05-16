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
     I     g_movecs,scftype,rw)
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
      character*(*) rw ! [in] read - write
c
      character*20 pname
      integer ispin
      integer bas,io_res
      character*255 filename
      integer inntsize
      logical nwmolden_movecs,nwmolden_basis,nwmolden_geom
      external nwmolden_movecs,nwmolden_basis,nwmolden_geom
      logical nwmolden_openfile,nwmolden_closefile
      external nwmolden_openfile,nwmolden_closefile
      logical do_io,do_norm_janpa,do_norm_nwchem
      logical int_norm_janpa,int_normalize
      external int_norm_janpa,int_normalize
c
c
      inntsize = MA_sizeof(MT_INT,1,MT_BYTE)
      do_io=ga_nodeid().eq.0
      pname='dft_nwmolden: '
      if(rw.eq.'write') then
         call util_file_name('molden', .false., .false.,filename)
      else
         if (.not. rtdb_cget(rtdb, 'vectors:molden', 1, filename))
     $        call errquit('dftnwmolden: no file',0, 0)
      endif
      io_res=1
      if (do_io) then
         if(.not.nwmolden_openfile(rw,filename))
     E        io_res=0
      endif
      call ga_brdcst(2020, io_res, inntsize, 0)
      if(io_res.eq.0) 
     E        call errquit(pname//'nwmolden_openfile failed ',0,0)
      if(.not.nwmolden_geom(rw,geom))
     E     call errquit(pname//'nwmolden_geom failed ',0,0)
c     create internal basis handle
      if (.not. bas_create(bas, 'ao basis'))
     &   call errquit(pname//'bas_create?', 2, BASIS_ERR)
      if (.not. bas_rtdb_load(rtdb, geom, bas, 'ao basis'))
     &   call errquit(pname//'no ao basis ', 3, BASIS_ERR)
      if (.not. rtdb_get(rtdb,'molden:do_norm_janpa', mt_log, 1, 
     A     do_norm_janpa))  do_norm_janpa=.false.
      if (.not. rtdb_get(rtdb,'molden:do_norm_nwchem', mt_log, 1, 
     A     do_norm_nwchem))  do_norm_nwchem=.false.
      if(ga_nodeid().eq.0) then
         write(6,*) ' do_norm_janpa ',do_norm_janpa
         write(6,*) ' do_norm_nwchem ',do_norm_nwchem
      endif
      if(do_norm_janpa) then
         if(.not.int_norm_janpa(rtdb,bas))
     &        call errquit(pname//'norm_janpa failed ', 5, BASIS_ERR)
      endif
      if(do_norm_nwchem) then
         if(.not.int_normalize(rtdb,bas))
     &        call errquit(pname//'int_normalize failed ', 5, BASIS_ERR)
      endif
      if(do_io) then
         if(.not.nwmolden_basis(rw,bas))
     E        call errquit(pname//'nwmolden_basis failed ',0,0)
c

         if(.not.nwmolden_movecs(rw,bas,nbf,nmo,
     1        eval,occ,ir,
     G        g_movecs,scftype))
     E        call errquit(pname//'nwmolden_movecs failed ',0,0)
         
      endif
c     broadcast evals when reading
      if(rw.eq.'read') then
         ispin=ma_sizeof(MT_DBL,nbf,MT_BYTE)
         if(scftype.eq.'UHF') ispin=ispin*2
         call ga_brdcst(1984, eval, ispin, 0)
      endif
      
c
         if (.not. bas_destroy(bas))
     &      call errquit(pname//'not able to destroy bas:',33,
     &       BASIS_ERR)
         if(do_io) then
         if(.not.nwmolden_closefile())
     E        call errquit(pname//'nwmolden_closefile failed ',0,0)
      if(.not.do_norm_janpa.and.rw.eq.'write') then
         write(6,*) ' '
         write(6,*) ' Molden file created '
         write(6,*) ' If you plan to use JANPA '
         write(6,*) '   run the .molden file through '
         write(6,*) '   molden2molden.jar -NormalizeBF to get the'
         write(6,*) '   correct Molden normalization '
         write(6,*) '  '
      endif
         endif
         call ga_sync()
         return
         end
      logical function int_norm_janpa(rtdb,basisin)
c $Id: int_norm.F 24111 2013-04-21 06:55:44Z d3y133 $
      implicit none
#include "apiP.fh"
#include "errquit.fh"
#include "bas.fh"
#include "nwc_const.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "bas_exndcf_dec.fh"
#include "stdio.fh"
      integer basisin  !< [Input] The basis set handle
      integer rtdb     !< [Input] The runtime database handle
      integer basis, nucont, iucont, type, nprim, ngen, ept, cpt, igen
      integer norm_typ
      logical norm_ok
      logical force_norm_rel
      character*30 pname
c
#include "bas_exndcf_sfn.fh"
      pname = "int_norm_janpa"
c
c     Normalize the contraction coefficients in this basis.
c     Eventually we need to store both the normalized and unnormalized
c     coefficients so that the user can print out something recognizable
c
      int_norm_janpa = bas_check_handle(basisin, 'int_norm_janpa')
      if (.not. int_norm_janpa) return
c
      basis = basisin + BASIS_HANDLE_OFFSET
c
c
      nucont = infbs_head(HEAD_NCONT,basis)
      do iucont = 1, nucont
         type  = infbs_cont(CONT_TYPE, iucont,basis)
         nprim = infbs_cont(CONT_NPRIM,iucont,basis) 
         ngen  = infbs_cont(CONT_NGEN, iucont,basis)
         ept   = infbs_cont(CONT_IEXP, iucont,basis)
         cpt   = infbs_cont(CONT_ICFP, iucont,basis)
         if (type .ge. 0) then
            do igen = 0,ngen-1
               call nmcoef_janpa2(dbl_mb(mb_exndcf(ept, basis)),
     &              dbl_mb(mb_exndcf(cpt+igen*nprim,basis)),
     &              type, nprim)
            enddo
c....... sp
         else if (type .eq. -1) then
            call errquit(pname//' not ready for itype ',type,0)
         else if (type .eq. -2) then
c........ spd
            call errquit(pname//' not ready for itype ',type,0)
         else
            call errquit(pname//'invalid contraction type', 0,
     &       BASIS_ERR)
         endif
      enddo
c
      end
c     normalization compatible with JANPA softwares requirements
C     T.Y.Nikolaienko, L.A.Bulavin, D.M.Hovorun; Comput.Theor.Chem. (2014)
C     V. 1050, P. 15-22, DOI: 10.1016/j.comptc.2014.10.002
      subroutine nmcoef_janpa2(z,cc,l,np)

      Implicit none
      Integer l,np ! [in]
      Double precision z(*)  ! [in]
      Double precision cc(*) ! [in/out] normalized in exit
c local
      Double precision f     
#include "basdeclsP.fh"
      double precision pi,sum
      integer l_lim
      integer i,j
      double precision ovlp_prim,z1
c     Overlap integral for a gaussian function
c     see M.Ho and J.M.Hernandez-Perez, The Mathematica Journal, 14(3), 2012. doi:10.3888/tmj.14-3.
c     see Libint programmers manual at
c     https://github.com/evaleev/libint/blob/master/doc/progman/progman.tex
      ovlp_prim(z1,pi,l,f)=(pi/z1)**1.5d0*f/(2.d0*z1)**l
c******************************************************************************
      pi=acos(-1d0)

c compute double factorial of (2*L-1)
      f = 1.0d0
      if(l.gt.1) then
         do i = (2*l-1),2,-2
            f = f*dble(i)
         end do
      endif
#ifdef SLOW
c
c     Normalize each primitive 
c
      do i = 1,np
        cc(i)= cc(i)/sqrt(ovlp_prim(2d0*z(i),pi,l,f))
      enddo

c Normalize the contractions.

      sum = 0.D0
      do i=1,np
         sum = sum + cc(i)*cc(i)*ovlp_prim(2d0*z(i),pi,l,f)
         do j= i+1,np
            sum = sum + 2d0*cc(i)*cc(j)*
     *           ovlp_prim(z(i)+z(j),pi,l,f)
         enddo
      enddo
#else
c
c Normalize the contractions.
c
      sum = 0.D0
      do i = 1,np
         sum = sum + cc(i)*cc(i)
c     Normalize each primitive 
         cc(i)= cc(i)/sqrt(ovlp_prim(2d0*z(i),pi,l,f))
      enddo
      do i=1,np
         do j= i+1,np
            sum = sum + 2d0*cc(i)*cc(j)*
     *           ovlp_prim(z(i)+z(j),pi,l,f)
         enddo
      enddo
#endif
c..   
      sum = 1d0/sqrt(sum)
      
      do  i = 1,np
         cc(i) = cc(i)*sum*
c     take out primitive normalization factor previously used
     *        sqrt(ovlp_prim(2d0*z(i),pi,l,f))
      enddo
      return
      end

