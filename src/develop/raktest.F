      subroutine raktest(rtdb)
      implicit none
c $Id: raktest.F,v 1.14 1994-06-23 00:18:15 vg038 Exp $
      integer rtdb
c
      write(6,*)' run rjh test '
c      call rjhtest(rtdb)
      write(6,*)' rjh test done '
      write(6,*) ' rak test'
      call dointsum(rtdb,'geometry', 'mo basis')
      call derivtest(rtdb,'geometry', 'mo basis')
c
      end
      subroutine dointsum(rtdb,geom_name, basis_name)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer basis, geom
      character*(*) geom_name, basis_name
c
      logical int_normalize
      double precision ddot
      external int_normalize
      external ddot
c
      call hf_print_set(0)
c
      write(6,*)' rtdb      :',rtdb
      write(6,*)' geom_name : <',geom_name,'>'
      write(6,*)' basis_name: <',basis_name,'>'
c
      if (.not. geom_create(geom, geom_name))
     $     call errquit('dointsum: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('dointsum: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('dointsum: geom_print failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('dointsum: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $      call errquit('dointsum: basis load failed', 0)
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
      if (.not. gbs_map_print(basis))
     $     call errquit('dointsum: gbs_map_print failed', 0)
c
c normalize basis set
c
      if (.not.int_normalize(basis))
     &       call errquit('dointsum: basis norm. failed', 0)
      write(6,*)' after normalization'
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
c
c      write(6,*)' with integrals  printed '
c      call int_chk_sum(basis,.true.)
      write(6,*)' without integrals  printed '
      call int_chk_sum(basis,.false.)
c
      if (.not.(
     &       (bas_destroy(basis)).and.(geom_destroy(geom))
     &       ))
     &       call errquit('rak:error destroying gbs pair',0)
      call MA_summarize_allocated_blocks()
c
      end
      subroutine derivtest(rtdb,geom_name, basis_name)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
c
      integer rtdb
      integer basis, geom
      integer nscr, hscr, iscr
      integer nbf, nshells, nat, maxg, size2e,size2ed
      integer h2e, i2e, h2efd, i2efd, h2ep, i2ep, h2em, i2em 
      integer h2ediff, i2ediff
      integer hbuf2e, ibuf2e
      integer buf1esz, hbuf1e, ibuf1e
      integer sa_size, hsa, isa, hsafd, isafd
      integer hta, ita, htafd, itafd
      integer threen
      integer hxyz, ixyz
      character*(*) geom_name, basis_name
      logical int_normalize
      double precision ddot
      external int_normalize
      external ddot
c
      call hf_print_set(0)
c
      write(6,*)' rtdb      :',rtdb
      write(6,*)' geom_name : <',geom_name,'>'
      write(6,*)' basis_name: <',basis_name,'>'
c
      if (.not. geom_create(geom, geom_name))
     $     call errquit('dointsum: geom_create failed?', 0)
      if (.not. geom_rtdb_load(rtdb, geom, geom_name))
     $     call errquit('dointsum: geom_load failed', 0)
      if (.not. geom_print(geom))
     $     call errquit('dointsum: geom_print failed', 0)
      if (.not. bas_create(basis, basis_name))
     $     call errquit('dointsum: basis create failed', 0)
      if (.not. bas_rtdb_load(rtdb, geom, basis, basis_name))
     $      call errquit('dointsum: basis load failed', 0)
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
      if (.not. gbs_map_print(basis))
     $     call errquit('dointsum: gbs_map_print failed', 0)
c
c normalize basis set
c
      if (.not.int_normalize(basis))
     &       call errquit('dointsum: basis norm. failed', 0)
      write(6,*)' after normalization'
      call intd_init(1,basis)
      if (.not. bas_print(basis))
     $      call errquit('dointsum: basis print failed', 0)
c
c      write(6,*)' with integrals  printed '
c      call int_chk_sum(basis,.true.)
      write(6,*)' without integrals  printed '
      call int_chk_sum(basis,.false.)

c
      write(6,*)' derivative test '
c
*      call hf_print_set(1)
*c
*      nscr = 1 000 000
*      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
*     &       call errquit('raktest: allocation scr failed',nscr)
*      call dfill(nscr,0.0,dbl_mb(iscr),1)
*      call onedint_fd(basis,geom,dbl_mb(iscr),nscr)
*      if (.not.MA_Pop_Stack(hscr))
*     &       write(6,*)'raktest: pop on hscr failed'
c
      
      if(.not.bas_numbf(basis,nbf))
     &       stop 'raktest:bas_numbf failed'
      if(.not.bas_numcont(basis,nshells))
     &       stop 'raktest:bas_numcont failed'
      if(.not.geom_ncent(geom,nat))
     &       stop 'raktest:geom_ncent failed'
c
      sa_size = (nbf*nbf)*(3*nat)
c..... space for sa
      if (.not.MA_Push_Get(MT_dbl,sa_size,'sa buffer',hsa,isa))
     &       call errquit('raktest: allocation for sa failed')
c..... space for safd
      if (.not.MA_Push_Get(MT_dbl,sa_size,'safd buffer',hsafd,isafd))
     &       call errquit('raktest: allocation for safd failed')
      threen = 3*nat
      if (.not.MA_Push_Get(MT_dbl,2*threen,'xyzxyz buffer',hxyz,ixyz))
     &       call errquit('raktest: allocation for safd failed')
      buf1esz = nbf*nbf*3
      if (.not.MA_Push_Get(MT_dbl,buf1esz,'safd buffer',hbuf1e,ibuf1e))
     &       call errquit('raktest: allocation for buf1e failed')
c..... scratch space
      nscr = 100 000
      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
     &       call errquit('raktest: allocation scr failed',nscr)
c..... space for Ta&Tafd so to sum it to Va&Fafd
      if (.not.MA_Push_Get(MT_Dbl,sa_size,'Ta',hta,ita))
     &       call errquit('raktest: allocation scr failed',sa_size)
      if (.not.MA_Push_Get(MT_Dbl,sa_size,'Ta fd',htafd,itafd))
     &       call errquit('raktest: allocation scr failed',sa_size)
      call dfill(sa_size,0.0,dbl_mb(ita),1)
      call dfill(sa_size,0.0,dbl_mb(itafd),1)
      call dfill(sa_size,0.0,dbl_mb(isa),1)
      call dfill(sa_size,0.0,dbl_mb(isafd),1)
      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call buildSa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,.true.)
#if defined(SUN)
      call flush(6)
#endif
      call buildSanew(dbl_mb(isafd),nat,nshells,nbf,basis,
     &       dbl_mb(iscr),nscr,.true.)
      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
      write(6,*)' Sanew difference norm =',
     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
#if defined(SUN)
      call flush(6)
#endif
      call buildSafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
     &       dbl_mb(ibuf1e+2*nbf*nbf),
     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
      write(6,*)' Sa difference norm =',
     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
#if defined(SUN)
      call flush(6)
#endif
c
c.... do T
      call dfill(sa_size,0.0,dbl_mb(isa),1)
      call dfill(sa_size,0.0,dbl_mb(isafd),1)
      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call buildTa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,.true.)
      call dcopy(sa_size,dbl_mb(isa),1,dbl_mb(ita),1) 
#if defined(SUN)
      call flush(6)
#endif
      call buildTafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
     &       dbl_mb(ibuf1e+2*nbf*nbf),
     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
      call dcopy(sa_size,dbl_mb(isafd),1,dbl_mb(itafd),1) 
      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
      write(6,*)' Ta difference norm =',
     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
#if defined(SUN)
      call flush(6)
#endif
c.. do V
      call dfill(sa_size,0.0,dbl_mb(isa),1)
      call dfill(sa_size,0.0,dbl_mb(isafd),1)
      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
      call dfill(buf1esz,0.0,dbl_mb(ibuf1e),1)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call buildVa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,.true.)
      call daxpy(sa_size,1.0d00,dbl_mb(isa),1,dbl_mb(ita),1)
#if defined(SUN)
      call flush(6)
#endif
      call buildVafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
     &       dbl_mb(ibuf1e+2*nbf*nbf),
     &       dbl_mb(ixyz),dbl_mb(ixyz+threen),.true.)
      call daxpy(sa_size,1.0d00,dbl_mb(isafd),1,dbl_mb(itafd),1)
      call daxpy(sa_size,-1.0d00,dbl_mb(isa),1,dbl_mb(isafd),1)
      call printdmat('Va difference matrix',dbl_mb(isafd),nat,nbf)
      write(6,*)' Va difference norm =',
     &       ddot(sa_size,dbl_mb(isafd),1,dbl_mb(isafd),1)
#if defined(SUN)
      call flush(6)
#endif
c... print H derivs
      call printdmat(' Ha ',dbl_mb(ita),nat,nbf)
      call printdmat(' Ha <<<FD>>> ',dbl_mb(itafd),nat,nbf)
      call daxpy(sa_size,-1.0d00,dbl_mb(ita),1,dbl_mb(itafd),1)
      write(6,*)' Ha difference norm =',
     &       ddot(sa_size,dbl_mb(itafd),1,dbl_mb(itafd),1)
c
      call buildHa(dbl_mb(itafd),nat,nshells,nbf,basis,
     &       dbl_mb(iscr),nscr,.true.)
      call printdmat(' Ha New',dbl_mb(itafd),nat,nbf)
      call daxpy(sa_size,-1.0d00,dbl_mb(ita),1,dbl_mb(itafd),1)
      write(6,*)' Ha new difference norm =',
     &       ddot(sa_size,dbl_mb(itafd),1,dbl_mb(itafd),1)
#if defined(SUN)
      call flush(6)
#endif
c
      if (.not.MA_Pop_Stack(htafd))
     &       call errquit('raktest: ma_pop_stack htafd failed',911)
      if (.not.MA_Pop_Stack(hta))
     &       call errquit('raktest: ma_pop_stack hta failed',911)
      if (.not.MA_Pop_Stack(hscr))
     &       call errquit('raktest: ma_pop_stack hscr failed',911)
      if (.not.MA_Pop_Stack(hbuf1e))
     &       call errquit('raktest: ma_pop_stack hbu1e failed',911)
      if (.not.MA_Pop_Stack(hxyz))
     &       call errquit('raktest: ma_pop_stack hxyz failed',911)
      if (.not.MA_Pop_Stack(hsafd))
     &       call errquit('raktest: ma_pop_stack hsafd failed',911)
      if (.not.MA_Pop_Stack(hsa))
     &       call errquit('raktest: ma_pop_stack hsa failed',911)
c
c
c
c... do two electron integral stuff
c
      call int_mem_2e4c(maxg,nscr)
      write(6,*)' nscr    => ',nscr
      write(6,*)' maxg    => ',maxg
      nscr = nscr*12    ! 8x scratch for derivatives
      maxg = 12*maxg    ! quadruple buffer size for derivatives
      size2e = nbf*(nbf+1)/2
      size2e = size2e*(size2e+1)/2
      size2ed= size2e*3*nat
      write(6,*)' nscr    => ',nscr
      write(6,*)' maxg    => ',maxg
      write(6,*)' nbf     => ',nbf
      write(6,*)' nat     => ',nat
      write(6,*)' size2e  => ',size2e
      write(6,*)' size2ed => ',size2ed
#if defined(SUN)
      call flush(6)
#endif
      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e d',h2ediff,i2ediff))
     &       call errquit('raktest: allocation h2ediff failed',size2ed)
      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e d',h2e,i2e))
     &       call errquit('raktest: allocation h2e failed',size2ed)
      if (.not.MA_Push_Get(MT_Dbl,size2ed,'2e fd',h2efd,i2efd))
     &       call errquit('raktest: allocation h2efd failed',size2ed)
      if (.not.MA_Push_Get(MT_Dbl,size2e,'2e +buf',h2ep,i2ep))
     &       call errquit('raktest: allocation h2ep failed',size2e)
      if (.not.MA_Push_Get(MT_Dbl,size2e,'2e -buf',h2em,i2em))
     &       call errquit('raktest: allocation h2em failed',size2e)
      if (.not.MA_Push_Get(MT_Dbl,nscr,'scr 2ed',hscr,iscr))
     &       call errquit('raktest: allocation hscr failed',nscr)
      if (.not.MA_Push_Get(MT_Dbl,maxg,'buf2e 2ed',hbuf2e,ibuf2e))
     &       call errquit('raktest: allocation hbuf2e failed',maxg)
      threen = 3*nat
      if (.not.MA_Push_Get(MT_dbl,2*threen,'xyzxyz buffer',hxyz,ixyz))
     &       call errquit('raktest: allocation for safd failed')
      call dfill(size2ed,0.0,dbl_mb(i2e),1)
      call dfill(size2ed,0.0,dbl_mb(i2efd),1)
      call dfill(size2e,0.0,dbl_mb(i2ep),1)
      call dfill(size2e,0.0,dbl_mb(i2em),1)
      call dfill(2*threen,0.0,dbl_mb(ixyz),1)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
      
      call build_d2e(geom,basis,nbf,nshells,nat,size2e,dbl_mb(i2e),
     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
*      call print2de(nat,size2e,dbl_mb(i2e),1.0d-05,nbf,'anal:d2e')
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
      call build_d2efd(geom,basis,nbf,nshells,nat,size2e,
     &       dbl_mb(i2efd),
     &       dbl_mb(i2ep),dbl_mb(i2em),dbl_mb(ixyz),
     &       dbl_mb(ixyz+threen),
     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
*      call print2de(nat,size2e,dbl_mb(i2efd),1.0d-05,nbf,'fd:d2e')
      call dcopy(size2ed,dbl_mb(i2efd),1,dbl_mb(i2ediff),1)
      call daxpy(size2ed,-1.0d00,dbl_mb(i2e),1,dbl_mb(i2ediff),1)
      write(6,*)' 2ed difference norm = ',
     &       ddot(size2ed,dbl_mb(i2ediff),1,dbl_mb(i2ediff),1)
*      call print2de_diff(nat,size2e,dbl_mb(i2ediff),
*     &       dbl_mb(i2e),dbl_mb(i2efd),1.0d-05,nbf)
      call dfill(size2ed,0.0,dbl_mb(i2efd),1)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call dfill(maxg,0.0,dbl_mb(ibuf2e),1)
      call build_d2e_new(basis,nbf,nshells,nat,size2e,dbl_mb(i2efd),
     &       maxg,dbl_mb(ibuf2e),nscr,dbl_mb(iscr))
      call dcopy(size2ed,dbl_mb(i2efd),1,dbl_mb(i2ediff),1)
      call daxpy(size2ed,-1.0d00,dbl_mb(i2e),1,dbl_mb(i2ediff),1)
      write(6,*)' new 2ed difference norm = ',
     &       ddot(size2ed,dbl_mb(i2ediff),1,dbl_mb(i2ediff),1)

      call print2de_diff(nat,size2e,dbl_mb(i2ediff),
     &       dbl_mb(i2e),dbl_mb(i2efd),1.0d-05,nbf)

      if (.not.MA_Pop_Stack(hxyz))
     &       call errquit('raktest: ma_pop_stack hxyz failed',911)
      if (.not.MA_Pop_Stack(hbuf2e))
     &       call errquit('raktest: ma_pop_stack hbuf2e failed',911)
      if (.not.MA_Pop_Stack(hscr))
     &       call errquit('raktest: ma_pop_stack hscr failed',911)
      if (.not.MA_Pop_Stack(h2em))
     &       call errquit('raktest: ma_pop_stack h2em failed',911)
      if (.not.MA_Pop_Stack(h2ep))
     &       call errquit('raktest: ma_pop_stack h2ep failed',911)
      if (.not.MA_Pop_Stack(h2efd))
     &       call errquit('raktest: ma_pop_stack h2efd failed',911)
      if (.not.MA_Pop_Stack(h2e))
     &       call errquit('raktest: ma_pop_stack h2e failed',911)
      if (.not.MA_Pop_Stack(h2ediff))
     &       call errquit('raktest: ma_pop_stack h2ediff failed',911)
      call int_terminate()
      if (.not.(
     &       (bas_destroy(basis)).and.(geom_destroy(geom))
     &       ))
     &       call errquit('rak:error destroying gbs pair',0)
      call MA_summarize_allocated_blocks()
      end
c.......................................................................
      Subroutine onedint_fd(basisin,geom,W0,maxW0)

      implicit none 

      Logical canAB,spherical
      integer basisin, basis, geom

c      double precision EPS
c      Parameter (EPS=1.D-15)

#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c----#include <dimensions.fh>
c----#include <atoms.fh>
c----#include <AO_basis_set.fh>
c----#include <timer.fh>

c--> Scratch Space

      integer ltime
      parameter (ltime= 69)
      integer maxW0, itime, i, nshells, i_top
      integer ucont, nshbfa, nshbfb, nint, ncenters
      integer iscr, nscr, i2, i1
      integer i_s0, i_s2, i_s1, i_sa
      integer i_t0, i_t2, i_t1, i_ta
      integer i_v0, i_v2, i_v1, i_va
      integer isha, ictra, igeoa, ipa, npa, ica, nca, la
      integer ishb, ictrb, igeob, ipb, npb, icb, ncb, lb
      double precision time(10)
      Double Precision W0(maxW0),Txyz(3,2)
      double precision delta
      double precision adifv, adifs, adift
c
c Compute analytic 1-electron integrals derivatives and compare to finite
c difference evaluation.
c
c******************************************************************************

      basis = basisin + BASIS_HANDLE_OFFSET
c Set timer.

      itime = 1

      do 10 i = 1,10
       time(i) = 0.D0
   10 continue

      spherical = .FALSE.

c Compute the 1-electron integrals (S, T, V).

c      open(1,file='1e')

      if (.not.bas_numcont(basisin,nshells))
     &       stop ' bas_numcont:error'
      
      ncenters = ncenter(geom)
      do 190 isha = 1,nshells

       do 180 ishb = 1,isha
  
c         write(6,*)' coords at <isha:ishb> <',isha,':',ishb,'>'
c         if (.not. geom_print(geom)) stop 'onedint_fd:geom_print'
c Set the canonicalization switches.

        if( isha.eq.ishb )then
         canAB = .false.
        else
         canAB = .FALSE.
        end if

c Define the basis sets on each of the centers.

        ictra = ibs_cn2ce(isha,basis)
        igeoa = ibs_geom(basis)
        ucont = ibs_cn2ucn(isha,basis)
        ipa   = infbs_cont(CONT_IEXP,ucont,basis)
        npa   = infbs_cont(CONT_NPRIM,ucont,basis)
        ica   = infbs_cont(CONT_ICFP,ucont,basis)
        nca   = infbs_cont(CONT_NGEN,ucont,basis)
        la    = infbs_cont(CONT_TYPE,ucont,basis)

        ictrb = ibs_cn2ce(ishb,basis)
        igeob = igeoa                       ! always same basis
        ucont = ibs_cn2ucn(ishb,basis)
        ipb   = infbs_cont(CONT_IEXP,ucont,basis)
        npb   = infbs_cont(CONT_NPRIM,ucont,basis)
        icb   = infbs_cont(CONT_ICFP,ucont,basis) 
        ncb   = infbs_cont(CONT_NGEN,ucont,basis) 
        lb    = infbs_cont(CONT_TYPE,ucont,basis) 

c Define the number of integrals in this block of shells.
 
        if( spherical )then
         nshbfa = 2*la+1
         nshbfb = 2*lb+1
        else
         nshbfa = ((la+1)*(la+2))/2
         nshbfb = ((lb+1)*(lb+2))/2
       end if
 
        if( canAB )then
         nint = (nshbfa*(nshbfa+1))/2
        else
         nint = nshbfa*nshbfb
        end if
 
c Assign pointers to scratch space.

        i_S0  = 1
        i_S1  = i_S0 + nint*6
        i_S2  = i_S1 + nint
        i_top = i_S2 + nint - 1

        i_T0  = i_top + 1
        i_T1  = i_T0 + nint*6
        i_T2  = i_T1 + nint
        i_top = i_T2 + nint - 1

        i_V0   = i_top  + 1
        i_V1   = i_V0   + nint*(3*ncenters)
        i_V2   = i_V1   + nint
        i_top  = i_V2   + nint

        if( i_top.gt.maxW0 )then

         write(*,*) 'ONEDINT:  Insufficient scratch space.'
         write(*,*) '          needed    ',i_top
         write(*,*) '          allocated ',maxW0

         write(*,*) 'S0:  ',i_S0
         write(*,*) 'S1:  ',i_S1
         write(*,*) 'S2:  ',i_S2

         write(*,*) 'T0:  ',i_T0
         write(*,*) 'T1:  ',i_T1
         write(*,*) 'T2:  ',i_T2

         write(*,*) 'V0:  ',i_V0
         write(*,*) 'V1:  ',i_V1
         write(*,*) 'V2:  ',i_V2

         stop

        end if

        iscr = i_top + 1
        nscr = maxW0 - i_top

c Analytic derivatives.

        i2  = 1
        i1  = 1

        i_Sa = i_S0 + ( ((i2-1)*3 + i1) - 1 )*nint
        i_Ta = i_T0 + ( ((i2-1)*3 + i1) - 1 )*nint
        i_Va = i_V0 + ( ((i2-1)*3 + i1) - 1 )*nint

                     itime = 1

c        call hf1d(xyz(1,ictra),zprims(ipa),ccoef(ica),npa,nca,la,
c     &            xyz(1,ictrb),zprims(ipb),ccoef(icb),npb,ncb,lb,
c     &            xyz,zan,ncenters,ictra,ictrb,
c     &            W0(i_Sa),W0(i_Ta),W0(i_Va),nint,
c     &            .TRUE.,.TRUE.,.TRUE.,canAB,W0(iscr),nscr)
c
c        write(6,*)' hfd1: <1> nscr:',nscr
        call hf1d(
     &   coords(1,ictra,geom),exndcf(ipa,basis),
     &                      exndcf(ica,basis),npa,nca,la,ictra,
     &   coords(1,ictrb,geom),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,ictrb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &          W0(i_Sa),W0(i_Ta),W0(i_Va),nint,
     &          .TRUE.,.TRUE.,.TRUE.,canAB,.false.,W0(iscr),nscr)
        
                     itime = 0

c Numeric derivatives of 2-ctr OIs and KEIs.

c Set finite difference perturbation.

        DELTA = 1.D-06

c Set coordinates for centers.

        Txyz(1,1) = coords(1,ictra,geom) ! xyz(1,ictra)
        Txyz(2,1) = coords(2,ictra,geom) ! xyz(2,ictra)
        Txyz(3,1) = coords(3,ictra,geom) ! xyz(3,ictra)

        Txyz(1,2) = coords(1,ictrb,geom) ! xyz(1,ictrb)
        Txyz(2,2) = coords(2,ictrb,geom) ! xyz(2,ictrb)
        Txyz(3,2) = coords(3,ictrb,geom) ! xyz(3,ictrb)

c Loop over the centers.

        do 120 i2 = 1,2

c Cartesian components of the 1st derivative.

         do 110 i1 = 1,3

c Plus delta.

          Txyz(i1,i2) = Txyz(i1,i2) + DELTA

c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
c       &             xyz,zan,ncenters,
c       &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
c       &             .TRUE.,.TRUE.,.FALSE.,canAB,W0(iscr),nscr)
c          write(6,*)' hf1: <1> nscr:',nscr
          call hf1(
     &   txyz(1,1),exndcf(ipa,basis),
     &                      exndcf(ica,basis),npa,nca,la,
     &   txyz(1,2),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
     &             .TRUE.,.TRUE.,.FALSE.,canAB,
     &           .false.,W0(iscr),nscr)

          Txyz(i1,i2) = Txyz(i1,i2) - DELTA

c Minus delta.

          Txyz(i1,i2) = Txyz(i1,i2) - DELTA

c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
c       &             xyz,zan,ncenters,
c       &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
c       &             .TRUE.,.TRUE.,.FALSE.,canAB,W0(iscr),nscr)
c          write(6,*)' hf1: <2> nscr:',nscr
          call hf1(
     &   txyz(1,1),exndcf(ipa,basis),
     &                      exndcf(ica,basis),npa,nca,la,
     &   txyz(1,2),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &           W0(i_S2),W0(i_T2),W0(i_V2),nint,
     &           .TRUE.,.TRUE.,.FALSE.,canAB,
     &           .false.,W0(iscr),nscr)
          Txyz(i1,i2) = Txyz(i1,i2) + DELTA

c Compute the 1st derivative of the integrals via a central difference formula.

          call cdif(W0(i_S1),W0(i_S2),DELTA,nint)
          call cdif(W0(i_T1),W0(i_T2),DELTA,nint)

c Check accuracy of derivatives.

          i_Sa = i_S0 + ( ((i2-1)*3 + i1) - 1 )*nint
          i_Ta = i_T0 + ( ((i2-1)*3 + i1) - 1 )*nint

          adifS = 0.D0
          adifT = 0.D0
          do 100 i = 1,nint
           adifS = max( adifS, abs( W0(i_Sa+(i-1)) - W0(i_S1+(i-1)) ) )
           adifT = max( adifT, abs( W0(i_Ta+(i-1)) - W0(i_T1+(i-1)) ) )
  100     continue

          if( adifS.gt.1.D-08)then
           write(21,'(a,e16.8)') 'max. abs. dif. ',adifS
           write(21,'(3i5)') isha,ictra,la
           write(21,'(3i5)') ishb,ictrb,lb
           write(21,'(3i5)') i1,i2
           write(21,'(i10,3e22.14)')
     &                 (i, W0(i_Sa+(i-1)),  W0(i_S1+(i-1)),
     &                    (W0(i_Sa+(i-1)) - W0(i_S1+(i-1))), i=1,nint)
           write(6,'(a,e16.8)') 'S:max. abs. dif. ',adifS
           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
           write(6,'(a,3i5)') 'i1,i2',i1,i2
           write(6,'(i10,3e22.14)')
     &                 (i, W0(i_Sa+(i-1)),  W0(i_S1+(i-1)),
     &                    (W0(i_Sa+(i-1)) - W0(i_S1+(i-1))), i=1,nint)
c          stop
          end if

          if( adifT.gt.1.D-08)then
           write(31,'(a,e16.8)') 'max. abs. dif. ',adifT
           write(31,'(3i5)') isha,ictra,la
           write(31,'(3i5)') ishb,ictrb,lb
           write(31,'(3i5)') i1,i2
           write(31,'(i10,3e22.14)')
     &                 (i, W0(i_Ta+(i-1)),  W0(i_T1+(i-1)),
     &                    (W0(i_Ta+(i-1)) - W0(i_T1+(i-1))), i=1,nint)
           write(6,'(a,e16.8)') 'T max. abs. dif. ',adifT
           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
           write(6,'(a,3i5)') 'i1,i2',i1,i2
           write(6,'(i10,3e22.14)')
     &                 (i, W0(i_Ta+(i-1)),  W0(i_T1+(i-1)),
     &                    (W0(i_Ta+(i-1)) - W0(i_T1+(i-1))), i=1,nint)
c          stop
          end if

  110    continue

  120   continue


c Numeric derivatives of NAIs.

c Loop over the centers.

        do 150 i2 = 1,ncenters

c Cartesian components of the 1st derivative.

         do 140 i1 = 1,3

c         if( (i2.eq.ictra) .or. (i2.eq.ictrb) ) go to 150

c Plus delta.

          if( i2.eq.ictra )then
           Txyz(i1,1) = Txyz(i1,1) + DELTA
          else if( i2.eq.ictrb )then
           Txyz(i1,2) = Txyz(i1,2) + DELTA
          end if

c          xyz(i1,i2) = xyz(i1,i2) + DELTA
          coords(i1,i2,geom) = coords(i1,i2,geom) + DELTA

c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
c       &             xyz,zan,ncenters,
c       &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
c       &             .FALSE.,.FALSE.,.TRUE.,canAB,W0(iscr),nscr)
c          write(6,*)' hf1: <3> nscr:',nscr
          call hf1(
     &   txyz(1,1),exndcf(ipa,basis),
     &                      exndcf(ica,basis),npa,nca,la,
     &   txyz(1,2),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &             W0(i_S1),W0(i_T1),W0(i_V1),nint,
     &             .FALSE.,.FALSE.,.TRUE.,canAB,
     &           .false.,W0(iscr),nscr)
          if( i2.eq.ictra )then
           Txyz(i1,1) = Txyz(i1,1) - DELTA
          else if( i2.eq.ictrb )then
           Txyz(i1,2) = Txyz(i1,2) - DELTA
          end if

c          xyz(i1,i2) = xyz(i1,i2) - DELTA
          coords(i1,i2,geom) = coords(i1,i2,geom) - DELTA

c Minus delta.

          if( i2.eq.ictra )then
           Txyz(i1,1) = Txyz(i1,1) - DELTA
          else if( i2.eq.ictrb )then
           Txyz(i1,2) = Txyz(i1,2) - DELTA
          end if

c          xyz(i1,i2) = xyz(i1,i2) - DELTA
          coords(i1,i2,geom) = coords(i1,i2,geom) - DELTA

c            call hf1(Txyz(1,1),zprims(ipa),ccoef(ica),npa,nca,la,
c       &             Txyz(1,2),zprims(ipb),ccoef(icb),npb,ncb,lb,
c       &             xyz,zan,ncenters,
c       &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
c       &             .FALSE.,.FALSE.,.TRUE.,canAB,W0(iscr),nscr)
c          write(6,*)' hf1: <4> nscr:',nscr
          call hf1(
     &   txyz(1,1),exndcf(ipa,basis),
     &                      exndcf(ica,basis),npa,nca,la,
     &   txyz(1,2),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &             W0(i_S2),W0(i_T2),W0(i_V2),nint,
     &             .FALSE.,.FALSE.,.TRUE.,canAB,
     &           .false.,W0(iscr),nscr)

          if( i2.eq.ictra )then
           Txyz(i1,1) = Txyz(i1,1) + DELTA
          else if( i2.eq.ictrb )then
           Txyz(i1,2) = Txyz(i1,2) + DELTA
          end if

c          xyz(i1,i2) = xyz(i1,i2) + DELTA
          coords(i1,i2,geom) = coords(i1,i2,geom) + DELTA

c Compute the 1st derivative of the integrals via a central difference formula.

          call cdif(W0(i_V1),W0(i_V2),DELTA,nint)

c Check accuracy of derivatives.

          i_Va = i_V0 + ( ((i2-1)*3 + i1) - 1 )*nint

          adifV = 0.D0
          do 130 i = 1,nint
           adifV = max( adifV, abs( W0(i_Va+(i-1)) - W0(i_V1+(i-1)) ) )
  130     continue

          if( adifV.gt.1.D-08)then
           write(41,'(a,e16.8)') 'max. abs. dif. ',adifV
           write(41,'(3i5)') isha,ictra,la
           write(41,'(3i5)') ishb,ictrb,lb
           write(41,'(3i5)') i1,i2
           write(41,'(i10,3e22.14)')
     &                 (i, W0(i_Va+(i-1)),  W0(i_V1+(i-1)),
     &                    (W0(i_Va+(i-1)) - W0(i_V1+(i-1))), i=1,nint)
           write(6,'(a,e16.8)') 'V max. abs. dif. ',adifV
           write(6,'(a,3i5)') 'isha,ictra,la',isha,ictra,la
           write(6,'(a,3i5)') 'ishb,ictrb,lb',ishb,ictrb,lb
           write(6,'(a,3i5)') 'i1,i2',i1,i2
           write(6,'(i10,3e22.14)')
     &                 (i, W0(i_Va+(i-1)),  W0(i_V1+(i-1)),
     &                    (W0(i_Va+(i-1)) - W0(i_V1+(i-1))), i=1,nint)
c          stop
          end if

  140    continue

  150   continue

c  c Determine the indices for this block of integrals.
c  
c          ioffseta = kloc(isha) - 1
c          ioffsetb = kloc(ishb) - 1
c  
c          nn = -1
c  
c          do 170 ma = 1,nshbfa
c  
c           if( canAB )then
c            mb_limit = ma
c           else
c            mb_limit = nshbfb
c           end if
c  
c           do 160 mb = 1,mb_limit
c  
c            nn = nn + 1
c  
c            i = ioffseta + ma
c            j = ioffsetb + mb
c  
c            ij = (i*(i-1))/2 + j
c  
c    160    continue
c  
c    170   continue

  180  continue

  190 continue

      close(1)

      itime = 1
c      if( itime.eq.1 )then
c       write(ltime,9999) 'hfset    ',time(1)
c       write(ltime,9999) 'hfmke    ',time(2)
c       write(ltime,9999) 'O2Is     ',time(3)
c       write(ltime,9999) 'KEIs     ',time(4)
c       write(ltime,9999) 'NAIs     ',time(5)
c      end if

 9999 format('  time spent in ',a10,20x,f20.6)

      end
c.......................................................................
      Subroutine cdif(f1,f2,delta,nfcns)

      Implicit real*8 (a-h,o-z)
      Implicit integer (i-n)

      Dimension f1(nfcns),f2(nfcns)
c
c Compute 1st derivatives via a central difference formula.
c
c******************************************************************************

      do 10 n = 1,nfcns
       f1(n) = (f1(n) - f2(n))/(2.D0*delta)
   10 continue

      end
c.......................................................................
      subroutine buildSa(Sa,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision sa(nbf,nbf,3*nat),scr(nscr)
      logical printit
c
      integer ish,jsh, basis, id, ids, ide, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF, compute
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildSa: geom does not match'
c
      FT = .true.
      FF = .false.
      call dfill(nbf*nbf*3*nat,0.0,sa,1)
      call dfill(nscr,0.0,scr,1)
c
      do 00100 ish = 1,ncont
        do 00200 jsh = 1,ish
          iatom = ibs_cn2ce(ish,basis)
          jatom = ibs_cn2ce(jsh,basis)
          compute = iatom.ne.jatom
          if (compute) then
            ucont = ibs_cn2ucn(ish,basis)
            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
            igen  = infbs_cont(CONT_NGEN,ucont,basis)
            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
            icf   = infbs_cont(CONT_ICFP,ucont,basis)
            itype = infbs_cont(CONT_TYPE,ucont,basis) 
            
            ucont = ibs_cn2ucn(jsh,basis)
            jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
            jgen  = infbs_cont(CONT_NGEN,ucont,basis)
            jexp  = infbs_cont(CONT_IEXP,ucont,basis)
            jcf   = infbs_cont(CONT_ICFP,ucont,basis)
            jtype = infbs_cont(CONT_TYPE,ucont,basis) 
                
            nint = (itype+1)*(itype+2)/2
            nint = nint * (jtype+1)*(jtype+2)/2
            used = nint*6
            iscr = used + 1
            call hf1d(
     &             coords(1,iatom,geom),
     &             exndcf(iexp,basis),
     &             exndcf(icf,basis),
     &             inp,igen,itype,iatom,
c
     &             coords(1,jatom,geom),
     &             exndcf(jexp,basis),
     &             exndcf(jcf,basis),
     &             jnp,jgen,jtype,jatom,
c
     &             coords(1,1,geom),charge(1,geom),nat,
     &             scr,scr,scr,nint,
     &             FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
            icount = 0
            ids = (iatom-1)*3+1
            ide = ids + 2
            do 00300 id = ids,ide
              do 00400 ibf = ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00500 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  sa(ibf,jbf,id) = scr(icount)
                  sa(jbf,ibf,id) = scr(icount)
00500           continue
00400         continue
00300       continue
            ids = (jatom-1)*3+1
            ide = ids + 2
            do 00600 id = ids,ide
              do 00700 ibf = ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
     &               ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  sa(ibf,jbf,id) = scr(icount)
                  sa(jbf,ibf,id) = scr(icount)
00800           continue
00700         continue
00600       continue
          endif
00200   continue
00100 continue

      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Sa <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Sa ',1.0d-05)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine buildSafd(Sa,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision sa(nbf,nbf,3*nat),scr(nscr)
      double precision bplus(nbf,nbf),bminus(nbf,nbf)
      double precision b(nbf,nbf)
      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
      double precision delta
      logical printit
c
      integer ish,jsh, basis, id, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildSafd: geom does not match'
c
      FT = .true.
      FF = .false.
      delta = 1.0d-06
      call dfill(nbf*nbf*3*nat,0.0,sa,1)
      call dfill(nscr,0.0,scr,1)
c
      id = 0
      do 00050 datom = 1,nat
        do 00060 dxyz = 1,3
          id = id + 1
          write(6,*)' Sa:fd <atom:',datom,'> <xyz:',dxyz,'>'
#if defined(SUN)
      call flush(6)
#endif
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*          call print_coords('   coords',coords(1,1,geom),nat)
*          call print_coords(' + coords',xyzxyzp,nat)
*          call print_coords(' - coords',xyzxyzm,nat)
          do 00100 ish = 1,ncont
            do 00200 jsh = 1,ncont
              ucont = ibs_cn2ucn(ish,basis)
              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
              igen  = infbs_cont(CONT_NGEN,ucont,basis)
              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
              icf   = infbs_cont(CONT_ICFP,ucont,basis)
              itype = infbs_cont(CONT_TYPE,ucont,basis) 
              iatom = ibs_cn2ce(ish,basis)
              
              ucont = ibs_cn2ucn(jsh,basis)
              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
              jatom = ibs_cn2ce(jsh,basis)
              
              nint = (itype+1)*(itype+2)/2
              nint = nint * (jtype+1)*(jtype+2)/2
              used = nint
              iscr = used + 1
c
c... do center
              call hf1(
     &               coords(1,iatom,geom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               coords(1,jatom,geom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               coords(1,1,geom),charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  b(jbf,ibf) = scr(icount)
00401           continue
00301         continue
c
c... do plus
              call dfill(nscr,0.0,scr,1)
              call hf1(
     &               xyzxyzp(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzp(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzp,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bplus(jbf,ibf) = scr(icount)
00400           continue
00300         continue
c
c........ do minus
              call dfill(nscr,0.0,scr,1)
c
              call hf1(
     &               xyzxyzm(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzm(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzm,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bminus(jbf,ibf) = scr(icount)
00600           continue
00500         continue
00200       continue
00100     continue
*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
          call cdif(bplus,bminus,delta,nbf*nbf)
          call dcopy(nbf*nbf,bplus,1,sa(1,1,id),1)
00060   continue
00050 continue
      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Sa <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Sa FD ',1.0d-05)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine print_coords(msg,xyzxyz,nat)
      implicit none
      character*(*) msg
      integer nat, i, j
      double precision xyzxyz(3,nat)
c
      write(6,*)' '
      write(6,*)' '
      write(6,10000)
      write(6,'(1x,a)')msg
      write(6,10000)
      do 00100 i=1,nat
        write(6,10001)i,(xyzxyz(j,i),j=1,3)
00100 continue
      write(6,10000)
      return
10000 format(100('-'))
10001 format(1x,'<atom:',i3,'>  <x:',1pd20.10,'>  <y:',1pd20.10,
     &       '>  <z:',1pd20.10,'>')
      end
c.......................................................................
      subroutine myoutput (z,rowlow,rowhi,collow,colhi,rowdim,coldim,
     $     nctl,msg,zerotol)
c.......................................................................
c myoutput prints a real*8 matrix in formatted form with numbered rows
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
c        msg.................informative message
c        zerotol.............information
c the parameters that follow matrix are all of type integer.  the
c program is set up to handle 5 columns/page with a 1p5d24.15 format for
c the columns.  if a different number of columns is required, change
c formats 1000 and 2000, and initialize kcol with the new number of
c columns.
c author;  nelson h.f. beebe, quantum theory project, university of
c          florida, gainesville
c modified: RJ Harrison - 8 columns 
c modified: Rick A. Kendall, HPCC group, Pacific Northwest Laboratory
c.......................................................................
      implicit none
      integer rowlow,rowhi,collow,colhi,rowdim,coldim,begin,kcol
      integer nctl, i, j, last, k
      double precision z(rowdim,coldim)
      double precision zerotol
      character*(*) msg
      character*8 asa(3), column, ctl, blank
      data column/'column'   /,asa/' ','00000000'  ,
     1     '--------'  /,blank/' '/
      data kcol/8/
c
      write(6,4000)msg
 4000 format(1x,80('-'),/,2x,a,/,1x,80('-'),/)
c
      do 11 i=rowlow,rowhi
         do 10 j=collow,colhi
            if (abs(z(i,j)).ge.zerotol) go to 15
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
         write (6,1000) (column,i,i = begin,last)
         do 1 k = rowlow,rowhi
            do 4 i=begin,last
               if (abs(z(k,i)).ge.zerotol) go to 5
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
 1000 format (/1h ,11x,7(a3,i3,3x),(a3,i3))
 2000 format (a1,'row',i4,1x,8f9.4)
      end
c.......................................................................
      subroutine buildTa(Ta,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision Ta(nbf,nbf,3*nat),scr(nscr)
      logical printit
c
      integer ish,jsh, basis, id, ids, ide, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF, compute
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildTa: geom does not match'
c
      FT = .true.
      FF = .false.
      call dfill(nbf*nbf*3*nat,0.0,Ta,1)
      call dfill(nscr,0.0,scr,1)
c
      do 00100 ish = 1,ncont
        do 00200 jsh = 1,ish
          iatom = ibs_cn2ce(ish,basis)
          jatom = ibs_cn2ce(jsh,basis)
          compute = iatom.ne.jatom
          if (compute) then
            ucont = ibs_cn2ucn(ish,basis)
            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
            igen  = infbs_cont(CONT_NGEN,ucont,basis)
            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
            icf   = infbs_cont(CONT_ICFP,ucont,basis)
            itype = infbs_cont(CONT_TYPE,ucont,basis) 
            
            ucont = ibs_cn2ucn(jsh,basis)
            jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
            jgen  = infbs_cont(CONT_NGEN,ucont,basis)
            jexp  = infbs_cont(CONT_IEXP,ucont,basis)
            jcf   = infbs_cont(CONT_ICFP,ucont,basis)
            jtype = infbs_cont(CONT_TYPE,ucont,basis) 
                
            nint = (itype+1)*(itype+2)/2
            nint = nint * (jtype+1)*(jtype+2)/2
            used = nint*6
            iscr = used + 1
            call hf1d(
     &             coords(1,iatom,geom),
     &             exndcf(iexp,basis),
     &             exndcf(icf,basis),
     &             inp,igen,itype,iatom,
c
     &             coords(1,jatom,geom),
     &             exndcf(jexp,basis),
     &             exndcf(jcf,basis),
     &             jnp,jgen,jtype,jatom,
c
     &             coords(1,1,geom),charge(1,geom),nat,
     &             scr,scr,scr,nint,
     &             FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
            icount = 0
            ids = (iatom-1)*3+1
            ide = ids + 2
            do 00300 id = ids,ide
              do 00400 ibf = ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00500 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  Ta(ibf,jbf,id) = scr(icount)
                  Ta(jbf,ibf,id) = scr(icount)
00500           continue
00400         continue
00300       continue
            ids = (jatom-1)*3+1
            ide = ids + 2
            do 00600 id = ids,ide
              do 00700 ibf = ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  Ta(ibf,jbf,id) = scr(icount)
                  Ta(jbf,ibf,id) = scr(icount)
00800           continue
00700         continue
00600       continue
          endif
00200   continue
00100 continue

      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Ta <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(Ta(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Ta ',1.0d-05)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine buildTafd(Ta,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision Ta(nbf,nbf,3*nat),scr(nscr)
      double precision bplus(nbf,nbf),bminus(nbf,nbf)
      double precision b(nbf,nbf)
      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
      double precision delta
      logical printit
c
      integer ish,jsh, basis, id, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildTafd: geom does not match'
c
      FT = .true.
      FF = .false.
      delta = 1.0d-06
      call dfill(nbf*nbf*3*nat,0.0,Ta,1)
      call dfill(nscr,0.0,scr,1)
c
      id = 0
      do 00050 datom = 1,nat
        do 00060 dxyz = 1,3
          id = id + 1
          write(6,*)' Ta:fd <atom:',datom,'> <xyz:',dxyz,'>'
#if defined(SUN)
      call flush(6)
#endif
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*          call print_coords('   coords',coords(1,1,geom),nat)
*          call print_coords(' + coords',xyzxyzp,nat)
*          call print_coords(' - coords',xyzxyzm,nat)
          do 00100 ish = 1,ncont
            do 00200 jsh = 1,ncont
              ucont = ibs_cn2ucn(ish,basis)
              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
              igen  = infbs_cont(CONT_NGEN,ucont,basis)
              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
              icf   = infbs_cont(CONT_ICFP,ucont,basis)
              itype = infbs_cont(CONT_TYPE,ucont,basis) 
              iatom = ibs_cn2ce(ish,basis)
              
              ucont = ibs_cn2ucn(jsh,basis)
              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
              jatom = ibs_cn2ce(jsh,basis)
              
              nint = (itype+1)*(itype+2)/2
              nint = nint * (jtype+1)*(jtype+2)/2
              used = nint
              iscr = used + 1
c
c... do center
              call hf1(
     &               coords(1,iatom,geom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               coords(1,jatom,geom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               coords(1,1,geom),charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  b(jbf,ibf) = scr(icount)
00401           continue
00301         continue
c
c... do plus
              call dfill(nscr,0.0,scr,1)
              call hf1(
     &               xyzxyzp(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzp(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzp,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bplus(jbf,ibf) = scr(icount)
00400           continue
00300         continue
c
c........ do minus
              call dfill(nscr,0.0,scr,1)
c
              call hf1(
     &               xyzxyzm(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzm(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzm,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FT,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bminus(jbf,ibf) = scr(icount)
00600           continue
00500         continue
00200       continue
00100     continue
*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
          call cdif(bplus,bminus,delta,nbf*nbf)
          call dcopy(nbf*nbf,bplus,1,Ta(1,1,id),1)
00060   continue
00050 continue
      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Ta <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(Ta(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Ta FD ',1.0d-05)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine buildVa(Va,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision Va(nbf,nbf,3*nat),scr(nscr)
      logical printit
c
      integer ish,jsh, basis, id, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildVa: geom does not match'
c
      FT = .true.
      FF = .false.
      call dfill(nbf*nbf*3*nat,0.0,Va,1)
      call dfill(nscr,0.0,scr,1)
c
      do 00100 ish = 1,ncont
        do 00200 jsh = 1,ish
          iatom = ibs_cn2ce(ish,basis)
          jatom = ibs_cn2ce(jsh,basis)

          ucont = ibs_cn2ucn(ish,basis)
          inp   = infbs_cont(CONT_NPRIM,ucont,basis)
          igen  = infbs_cont(CONT_NGEN,ucont,basis)
          iexp  = infbs_cont(CONT_IEXP,ucont,basis)
          icf   = infbs_cont(CONT_ICFP,ucont,basis)
          itype = infbs_cont(CONT_TYPE,ucont,basis) 
          
          ucont = ibs_cn2ucn(jsh,basis)
          jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
          jgen  = infbs_cont(CONT_NGEN,ucont,basis)
          jexp  = infbs_cont(CONT_IEXP,ucont,basis)
          jcf   = infbs_cont(CONT_ICFP,ucont,basis)
          jtype = infbs_cont(CONT_TYPE,ucont,basis) 
          
          nint = (itype+1)*(itype+2)/2
          nint = nint * (jtype+1)*(jtype+2)/2
          used = nint*3*nat
          iscr = used + 1
          call hf1d(
     &           coords(1,iatom,geom),
     &           exndcf(iexp,basis),
     &           exndcf(icf,basis),
     &           inp,igen,itype,iatom,
c
     &           coords(1,jatom,geom),
     &           exndcf(jexp,basis),
     &           exndcf(jcf,basis),
     &           jnp,jgen,jtype,jatom,
c
     &           coords(1,1,geom),charge(1,geom),nat,
     &           scr,scr,scr,nint,
     &           FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
          icount = 0
          id = 0
          do 00300 datom = 1,nat
            do 00400 dxyz = 1,3
              id = id + 1
              do 00500 ibf = ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  Va(ibf,jbf,id) = scr(icount)
                  Va(jbf,ibf,id) = scr(icount)
00600           continue
00500         continue
00400       continue
00300     continue
00200   continue
00100 continue

      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Va <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(Va(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Va ',1.0d-05)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine buildVafd(Va,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,b,bplus,bminus,xyzxyzp,xyzxyzm,printit)
      implicit none
c
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c
      integer nat,ncont,nbf,basisin,nscr,geom
      double precision Va(nbf,nbf,3*nat),scr(nscr)
      double precision bplus(nbf,nbf),bminus(nbf,nbf)
      double precision b(nbf,nbf)
      double precision xyzxyzp(3,nat),xyzxyzm(3,nat)
      double precision delta
      logical printit
c
      integer ish,jsh, basis, id, icount
      integer testgeom, ucont, nint, iscr, used, jbf, ibf
      integer datom,dxyz
      integer iatom,inp,igen,iexp,icf,itype
      integer jatom,jnp,jgen,jexp,jcf,jtype
      logical FT, FF
c
      basis = basisin + BASIS_HANDLE_OFFSET
      testgeom = ibs_geom(basis)
      if (testgeom.ne.geom) stop 'buildVafd: geom does not match'
c
      FT = .true.
      FF = .false.
      delta = 1.0d-06
      call dfill(nbf*nbf*3*nat,0.0,Va,1)
      call dfill(nscr,0.0,scr,1)
c
      id = 0
      do 00050 datom = 1,nat
        do 00060 dxyz = 1,3
          id = id + 1
          write(6,*)' Va:fd <atom:',datom,'> <xyz:',dxyz,'>'
#if defined(SUN)
      call flush(6)
#endif
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzp,1)
          xyzxyzp(dxyz,datom) = xyzxyzp(dxyz,datom) + delta
          call dcopy(3*nat,coords(1,1,geom),1,xyzxyzm,1)
          xyzxyzm(dxyz,datom) = xyzxyzm(dxyz,datom) - delta
*          call print_coords('   coords',coords(1,1,geom),nat)
*          call print_coords(' + coords',xyzxyzp,nat)
*          call print_coords(' - coords',xyzxyzm,nat)
          do 00100 ish = 1,ncont
            do 00200 jsh = 1,ncont
              ucont = ibs_cn2ucn(ish,basis)
              inp   = infbs_cont(CONT_NPRIM,ucont,basis)
              igen  = infbs_cont(CONT_NGEN,ucont,basis)
              iexp  = infbs_cont(CONT_IEXP,ucont,basis)
              icf   = infbs_cont(CONT_ICFP,ucont,basis)
              itype = infbs_cont(CONT_TYPE,ucont,basis) 
              iatom = ibs_cn2ce(ish,basis)
              
              ucont = ibs_cn2ucn(jsh,basis)
              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
              jatom = ibs_cn2ce(jsh,basis)
              
              nint = (itype+1)*(itype+2)/2
              nint = nint * (jtype+1)*(jtype+2)/2
              used = nint
              iscr = used + 1
c
c... do center
              call hf1(
     &               coords(1,iatom,geom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               coords(1,jatom,geom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               coords(1,1,geom),charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00301 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00401 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  b(jbf,ibf) = scr(icount)
00401           continue
00301         continue
c
c... do plus
              call dfill(nscr,0.0,scr,1)
              call hf1(
     &               xyzxyzp(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzp(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzp,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00300 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bplus(jbf,ibf) = scr(icount)
00400           continue
00300         continue
c
c........ do minus
              call dfill(nscr,0.0d00,scr,1)
c
              call hf1(
     &               xyzxyzm(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyzm(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyzm,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FF,FF,FT,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00500 ibf=ibs_cn2bfr(1,ish,basis),
     &               ibs_cn2bfr(2,ish,basis)
                do 00600 jbf = ibs_cn2bfr(1,jsh,basis),
     &                 ibs_cn2bfr(2,jsh,basis)
                  icount = icount + 1
                  bminus(jbf,ibf) = scr(icount)
00600           continue
00500         continue
00200       continue
00100     continue
*          write(6,*)' b  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(b,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
*          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
*          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
          call cdif(bplus,bminus,delta,nbf*nbf)
          call dcopy(nbf*nbf,bplus,1,Va(1,1,id),1)
00060   continue
00050 continue
      if (.not.printit) return
      id = 0
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,i2,a,i2,a)')
     &           ' Va <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(Va(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           ' Va FD ',1.0d-05)
01000   continue
00900 continue
      end
      subroutine printdmat(msg,Xa,nat,nbf)
      implicit none
c
      character*(*) msg
      integer nat,nbf
      double precision Xa(nbf,nbf,3*nat)
c
      integer datom, dxyz, id
c
      id = 0
      do 00100 datom = 1,nat
        do 00200 dxyz = 1,3
          id = id + 1
          write(6,'(//,a,a,i2,a,i2,a)')
     &           msg,' <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(Xa(1,1,id),1,nbf,1,nbf,nbf,nbf,1,
     &           msg,1.0d-05)
00200   continue
00100 continue
      return
      end
      subroutine build_d2e(geom,basisin,nbf,ncont,nat,intsize,d2e,
     &       maxg,gbuf,nscr,scr)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::passed
      integer geom,basisin,nbf,nat,intsize,maxg,nscr,ncont
      double precision d2e(intsize,3,nat)
      double precision gbuf(maxg)
      double precision scr(nscr)
c::local
      integer ucont,basis,nint,iduse,ipfrom,ipto,icount,id
      integer idatom(4)
      integer index,ixyz
c
C     SET-UP FUNCTIONS
C
      integer ibf, ish, iatom, inp, igen, iexp, icf, itype
      integer jbf, jsh, jatom, jnp, jgen, jexp, jcf, jtype
      integer kbf, ksh, katom, knp, kgen, kexp, kcf, ktype
      integer lbf, lsh, latom, lnp, lgen, lexp, lcf, ltype
c
c...  statement function
C
C     SET-UP FUNCTIONS
      integer min,max
      intrinsic min,max
      integer i,j,k,l,isym2,isym4
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
     &               (max(isym2(i,j),isym2(k,l))-1)/2+
     &               min(isym2(i,j),isym2(k,l))
c
      basis = basisin + BASIS_HANDLE_OFFSET
c      
      do 00100 ish = 1,ncont
c.
        iatom = ibs_cn2ce(ish,basis)
        ucont = ibs_cn2ucn(ish,basis)
        inp   = infbs_cont(CONT_NPRIM,ucont,basis)
        igen  = infbs_cont(CONT_NGEN,ucont,basis)
        iexp  = infbs_cont(CONT_IEXP,ucont,basis)
        icf   = infbs_cont(CONT_ICFP,ucont,basis)
        itype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
        do 00200 jsh= 1,ish
c.
          jatom = ibs_cn2ce(jsh,basis)
          ucont = ibs_cn2ucn(jsh,basis)
          jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
          jgen  = infbs_cont(CONT_NGEN,ucont,basis)
          jexp  = infbs_cont(CONT_IEXP,ucont,basis)
          jcf   = infbs_cont(CONT_ICFP,ucont,basis)
          jtype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
          do 00300 ksh = 1,jsh
c.
            katom = ibs_cn2ce(ksh,basis)
            ucont = ibs_cn2ucn(ksh,basis)
            knp   = infbs_cont(CONT_NPRIM,ucont,basis)
            kgen  = infbs_cont(CONT_NGEN,ucont,basis)
            kexp  = infbs_cont(CONT_IEXP,ucont,basis)
            kcf   = infbs_cont(CONT_ICFP,ucont,basis)
            ktype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
            do 00400 lsh = 1,ksh
c.
              latom = ibs_cn2ce(lsh,basis)
              ucont = ibs_cn2ucn(lsh,basis)
              lnp   = infbs_cont(CONT_NPRIM,ucont,basis)
              lgen  = infbs_cont(CONT_NGEN,ucont,basis)
              lexp  = infbs_cont(CONT_IEXP,ucont,basis)
              lcf   = infbs_cont(CONT_ICFP,ucont,basis)
              ltype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
              nint = (itype+1)*(itype+2)/2
              nint = nint*(jtype+1)*(jtype+2)/2
              nint = nint*(ktype+1)*(ktype+2)/2
              nint = nint*(ltype+1)*(ltype+2)/2
c
              if (
     &               iatom.eq.jatom.and.
     &               jatom.eq.katom.and.
     &               katom.eq.latom
     &               ) then
                continue
              else
                call hf2d(
     &                 coords(1,iatom,geom),exndcf(iexp,basis),
     &                 exndcf(icf,basis),inp,igen,itype,iatom,
c
     &                 coords(1,jatom,geom),exndcf(jexp,basis),
     &                 exndcf(jcf,basis),jnp,jgen,jtype,jatom,
c
     &                 coords(1,katom,geom),exndcf(kexp,basis),
     &                 exndcf(kcf,basis),knp,kgen,ktype,katom,
c
     &                 coords(1,latom,geom),exndcf(lexp,basis),
     &                 exndcf(lcf,basis),lnp,lgen,ltype,latom,
c
     &                 gbuf,nint,.false.,.false.,.false.,.false.,
     &                 scr,nscr)
c 
                call logic_atom(idatom,iatom,jatom,katom,latom)
c do sums
                do 00500 iduse=2,4
                  if(idatom(iduse).gt.0) then 
                    continue
                  else
                    ipto   = (abs(idatom(iduse)) -1)*nint*3 + 1
                    ipfrom = (iduse-1)*nint*3 + 1
                    call daxpy(nint*3,1.0d00,
     &                     gbuf(ipfrom),1,gbuf(ipto),1)
                  endif
00500           continue
c
c ... fill integral buffers              
                icount = 0
                do 00600 iduse = 1,4
                  if (idatom(iduse).lt.1) then
                    icount = icount + nint*3
                  else
                    id = idatom(iduse)
                    do 00700 ixyz = 1,3
                      do 00800 ibf = ibs_cn2bfr(1,ish,basis),
     &                       ibs_cn2bfr(2,ish,basis)
                        do 00900 jbf = ibs_cn2bfr(1,jsh,basis),
     &                         ibs_cn2bfr(2,jsh,basis)
                          do 01000 kbf = ibs_cn2bfr(1,ksh,basis),
     &                           ibs_cn2bfr(2,ksh,basis)
                            do 01100 lbf = ibs_cn2bfr(1,lsh,basis),
     &                             ibs_cn2bfr(2,lsh,basis)
                              index = isym4(ibf,jbf,kbf,lbf)
                              icount = icount + 1
                              d2e(index,ixyz,id) = gbuf(icount)
01100                       continue
01000                     continue
00900                   continue
00800                 continue
00700               continue
                  endif
00600           continue
c
              endif
c
00400       continue
00300     continue
00200   continue
00100 continue
      end
      subroutine logic_atom(idat,iat,jat,kat,lat)
      implicit none
      integer iat,jat,kat,lat
      integer idat(4)
c
      idat(1) = iat
      idat(2) = jat
      idat(3) = kat
      idat(4) = lat
      if (iat.eq.jat) idat(2) = -1
      if (iat.eq.kat) idat(3) = -1
      if (iat.eq.lat) idat(4) = -1
      if (jat.eq.kat) then
        if(idat(2).gt.0) then
          idat(3) = -2
        else
          idat(3) = idat(2)
        endif
      endif
      if (jat.eq.lat) then
        if(idat(2).gt.0) then
          idat(4) = -2
        else
          idat(4) = idat(2)
        endif
      endif
      if (kat.eq.lat) then
        if(idat(3).gt.0) then
          idat(4) = -3
        else
          idat(4) = idat(3)
        endif
      endif
      end
      subroutine build_d2efd(geom,basisin,nbf,ncont,nat,intsize,d2e,
     &       bp,bm,xyzp,xyzm,
     &       maxg,gbuf,nscr,scr)
      implicit none
#include "bas.fh"
#include "basP.fh"
#include "basdeclsP.fh"
#include "geom.fh"
#include "geomP.fh"
#include "geobasmapP.fh"
c::passed
      integer geom,basisin,nbf,nat,intsize,maxg,nscr,ncont
      double precision d2e(intsize,3,nat)
      double precision gbuf(maxg)
      double precision scr(nscr)
      double precision bp(intsize),bm(intsize)
      double precision xyzp(3,nat),xyzm(3,nat)
c:: local
      integer basis
      double precision delta
      integer dxyz,datom,nint,ucont,index,icount
      integer ibf, ish, iatom, inp, igen, iexp, icf, itype
      integer jbf, jsh, jatom, jnp, jgen, jexp, jcf, jtype
      integer kbf, ksh, katom, knp, kgen, kexp, kcf, ktype
      integer lbf, lsh, latom, lnp, lgen, lexp, lcf, ltype
c
c
c...  statement function
C
C     SET-UP FUNCTIONS
      integer min,max
      intrinsic min,max
      integer i,j,k,l,isym2,isym4
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
     &               (max(isym2(i,j),isym2(k,l))-1)/2+
     &               min(isym2(i,j),isym2(k,l))
c
      delta = 1.0d-06
      basis = basisin + BASIS_HANDLE_OFFSET
c      
      do 00100 datom = 1,nat
        do 00200 dxyz = 1,3
          call dcopy(3*nat,coords(1,1,geom),1,xyzp,1)
          call dcopy(3*nat,coords(1,1,geom),1,xyzm,1)
          xyzp(dxyz,datom) =  xyzp(dxyz,datom) + delta
          xyzm(dxyz,datom) =  xyzm(dxyz,datom) - delta
          call print_coords('2e +',xyzp,nat)
          call print_coords('2e -',xyzm,nat)
c
          do 00300 ish = 1,ncont
c.
            iatom = ibs_cn2ce(ish,basis)
            ucont = ibs_cn2ucn(ish,basis)
            inp   = infbs_cont(CONT_NPRIM,ucont,basis)
            igen  = infbs_cont(CONT_NGEN,ucont,basis)
            iexp  = infbs_cont(CONT_IEXP,ucont,basis)
            icf   = infbs_cont(CONT_ICFP,ucont,basis)
            itype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
            do 00400 jsh = 1,ish
c.
              jatom = ibs_cn2ce(jsh,basis)
              ucont = ibs_cn2ucn(jsh,basis)
              jnp   = infbs_cont(CONT_NPRIM,ucont,basis)
              jgen  = infbs_cont(CONT_NGEN,ucont,basis)
              jexp  = infbs_cont(CONT_IEXP,ucont,basis)
              jcf   = infbs_cont(CONT_ICFP,ucont,basis)
              jtype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
              do 00500 ksh = 1,jsh
c.
                katom = ibs_cn2ce(ksh,basis)
                ucont = ibs_cn2ucn(ksh,basis)
                knp   = infbs_cont(CONT_NPRIM,ucont,basis)
                kgen  = infbs_cont(CONT_NGEN,ucont,basis)
                kexp  = infbs_cont(CONT_IEXP,ucont,basis)
                kcf   = infbs_cont(CONT_ICFP,ucont,basis)
                ktype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
                do 00600 lsh = 1,ksh
c.
                  latom = ibs_cn2ce(lsh,basis)
                  ucont = ibs_cn2ucn(lsh,basis)
                  lnp   = infbs_cont(CONT_NPRIM,ucont,basis)
                  lgen  = infbs_cont(CONT_NGEN,ucont,basis)
                  lexp  = infbs_cont(CONT_IEXP,ucont,basis)
                  lcf   = infbs_cont(CONT_ICFP,ucont,basis)
                  ltype = infbs_cont(CONT_TYPE,ucont,basis) 
c.
                  nint = (itype+1)*(itype+2)/2
                  nint = nint*(jtype+1)*(jtype+2)/2
                  nint = nint*(ktype+1)*(ktype+2)/2
                  nint = nint*(ltype+1)*(ltype+2)/2
c
c... do plus
                  call dfill(nscr,0.0d00,scr,1)
                  call dfill(maxg,0.0d00,gbuf,1)
                  call hf2(
     &                   xyzp(1,iatom),exndcf(iexp,basis),
     &                   exndcf(icf,basis),inp,igen,itype,
c
     &                   xyzp(1,jatom),exndcf(jexp,basis),
     &                   exndcf(jcf,basis),jnp,jgen,jtype,
c
     &                   xyzp(1,katom),exndcf(kexp,basis),
     &                   exndcf(kcf,basis),knp,kgen,ktype,
c
     &                   xyzp(1,latom),exndcf(lexp,basis),
     &                   exndcf(lcf,basis),lnp,lgen,ltype,
c
     &                   gbuf,nint,.false.,.false.,.false.,.false.,
     &                   scr,nscr)
c
                  icount = 0
                  do 00700 ibf = ibs_cn2bfr(1,ish,basis),
     &                   ibs_cn2bfr(2,ish,basis)
                    do 00800 jbf = ibs_cn2bfr(1,jsh,basis),
     &                     ibs_cn2bfr(2,jsh,basis)
                      do 00900 kbf = ibs_cn2bfr(1,ksh,basis),
     &                       ibs_cn2bfr(2,ksh,basis)
                        do 01000 lbf = ibs_cn2bfr(1,lsh,basis),
     &                         ibs_cn2bfr(2,lsh,basis)
                          index = isym4(ibf,jbf,kbf,lbf)
                          icount = icount + 1
                          bp(index) = gbuf(icount)
01000                   continue
00900                 continue
00800               continue
00700             continue
c
c... do minus
                  call dfill(nscr,0.0d00,scr,1)
                  call dfill(maxg,0.0d00,gbuf,1)
                  call hf2(
     &                   xyzm(1,iatom),exndcf(iexp,basis),
     &                   exndcf(icf,basis),inp,igen,itype,
c
     &                   xyzm(1,jatom),exndcf(jexp,basis),
     &                   exndcf(jcf,basis),jnp,jgen,jtype,
c
     &                   xyzm(1,katom),exndcf(kexp,basis),
     &                   exndcf(kcf,basis),knp,kgen,ktype,
c
     &                   xyzm(1,latom),exndcf(lexp,basis),
     &                   exndcf(lcf,basis),lnp,lgen,ltype,
c
     &                   gbuf,nint,.false.,.false.,.false.,.false.,
     &                   scr,nscr)
c
                  icount = 0
                  do 00701 ibf = ibs_cn2bfr(1,ish,basis),
     &                   ibs_cn2bfr(2,ish,basis)
                    do 00801 jbf = ibs_cn2bfr(1,jsh,basis),
     &                     ibs_cn2bfr(2,jsh,basis)
                      do 00901 kbf = ibs_cn2bfr(1,ksh,basis),
     &                       ibs_cn2bfr(2,ksh,basis)
                        do 01001 lbf = ibs_cn2bfr(1,lsh,basis),
     &                         ibs_cn2bfr(2,lsh,basis)
                          index = isym4(ibf,jbf,kbf,lbf)
                          icount = icount + 1
                          bm(index) = gbuf(icount)
01001                   continue
00901                 continue
00801               continue
00701             continue
c
00600           continue
00500         continue
00400       continue
00300     continue
          call cdif(bp,bm,delta,intsize)
          call dcopy(intsize,bp,1,d2e(1,dxyz,datom),1)
00200   continue
00100 continue
c
      end
      subroutine print2de_diff(nat,size2e,diff,d2e,d2efd,thresh,nbf)
      implicit none 
      integer nat,size2e,index
      double precision diff (size2e,3,nat)
      double precision d2e  (size2e,3,nat)
      double precision d2efd(size2e,3,nat)
      double precision thresh, ra2fd,  rfd2a
      integer ibf,jbf,kbf,lbf
      integer nbf
c
      integer dxyz,datom

      do 00100 datom = 1,nat
        do 00200 dxyz = 1,3
          do 00300 index = 1,size2e
            if (abs(diff(index,dxyz,datom)).gt.thresh) then
              call getfromi(index,nbf,ibf,jbf,kbf,lbf)
              if (abs(d2efd(index,dxyz,datom)).gt.1.0d-05) then
                ra2fd = d2e(index,dxyz,datom)/d2efd(index,dxyz,datom)
              else
                ra2fd = -565.6589d00
              endif
              if (abs(d2e(index,dxyz,datom)).gt.1.0d-05) then
                rfd2a = d2efd(index,dxyz,datom)/d2e(index,dxyz,datom)
              else
                rfd2a = -6.0201023d00
              endif
              write(6,10000)
     &               ibf,jbf,kbf,lbf,
     &               index,dxyz,datom,diff(index,dxyz,datom),
     &               d2e(index,dxyz,datom),d2efd(index,dxyz,datom),
     &               ra2fd,rfd2a
#if defined(SUN)
              call flush(6)
#endif              
            endif
00300     continue
00200   continue
00100 continue
10000 format(' <',i2,i2,'|',i2,i2,'||',i5,'|xyz:',i2,'|a:',i2,'>  ',
     &       '<diff:',1pd14.6,'> ','<anal:',1pd14.6,'> ',
     &       '<fd:',1pd14.6,'> <ra2fd:',1pd14.6,'> <rfd2a:',
     &       1pd14.6,'>')
      end
      subroutine getfromi(index,nbf,i,j,k,l)
      implicit none
      integer index,nbf,i,j,k,l,lhi
c
      integer icount
c
      icount = 1
      do 00100 i=1,nbf
        do 00200 j=1,i
          do 00300 k=1,i
            lhi = k
            if (k.eq.i) lhi = j
            do 00400 l=1,lhi
              if (icount.eq.index) return
              icount = icount + 1
00400       continue
00300     continue
00200   continue
00100 continue
      write(6,*)' index/nbf =',index,nbf
      stop ' getfromi fucked '
      end
      subroutine print2de(nat,size2e,d2e,thresh,nbf,msg)
      implicit none 
      integer nat,size2e,index
      double precision d2e  (size2e,3,nat)
      double precision thresh
      integer ibf,jbf,kbf,lbf
      character*(*)msg
      integer nbf
c
      integer dxyz,datom

      do 00100 datom = 1,nat
        do 00200 dxyz = 1,3
          do 00300 index = 1,size2e
            if (abs(d2e(index,dxyz,datom)).gt.thresh) then
              call getfromi(index,nbf,ibf,jbf,kbf,lbf)
              write(6,10000) ibf,jbf,kbf,lbf,
     &               index,dxyz,datom,msg,d2e(index,dxyz,datom)
#if defined(SUN)
              call flush(6)
#endif              
            endif
00300     continue
00200   continue
00100 continue
10000 format(' <',i2,i2,'|',i2,i2,'||',i5,'|xyz:',i2,'|a:',i2,'>  ',
     &       '<',a,':',1pd14.6,'> ')
      end
c.......................................................................
      subroutine buildSanew(Sa,nat,ncont,nbf,basis,
     &       scr,nscr,printit)
      implicit none
c
#include "bas.fh"
c
      integer nat,ncont,nbf,basis,nscr,idatom(2)
      double precision sa(nbf,nbf,3,nat),scr(nscr)
      logical printit
c
      integer used, iscr, newnscr, nbfmax
      integer ish, ibf, ibflo, ibfhi
      integer jsh, jbf, jbflo, jbfhi
      integer id, dxyz, datom, icount
c
      call dfill(nbf*nbf*3*nat,0.0,sa,1)
      call dfill(nscr,0.0,scr,1)
c
      if(.not.bas_nbf_cn_max(basis,nbfmax))
     &       call errquit('raktest:',922)
c
      call hf_print_set(1)
      used = nbfmax*nbfmax*3*2
      iscr = used + 1
      newnscr = nscr - used

      do 00100 ish = 1,ncont
        do 00200 jsh = 1,ish

          call intd_1eov(basis,ish,basis,jsh,newnscr,scr(iscr),
     &           used,scr,idatom)

          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &           call errquit('error',911)
          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &           call errquit('error',911)
          
          icount = 0
          do 00300 id = 1,2
            if (idatom(id).gt.0) then
              do 00400 dxyz = 1,3
                do 00500 ibf = ibflo,ibfhi
                  do 00600 jbf = jbflo,jbfhi
                    icount = icount + 1
                    sa(ibf,jbf,dxyz,idatom(id)) = scr(icount)
                    sa(jbf,ibf,dxyz,idatom(id)) = scr(icount)
00600             continue
00500           continue
00400         continue
            else
              icount = icount + 3*(ibfhi-ibflo+1)*(jbfhi-jbflo+1)
            endif
00300     continue
00200   continue
00100 continue

      if (.not.printit) return
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          write(6,'(//,a,i2,a,i2,a)')
     &           ' SaNEW <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(sa(1,1,dxyz,datom),1,nbf,1,nbf,nbf,nbf,1,
     &           ' SaNEW ',1.0d-05)
01000   continue
00900 continue
      call hf_print_set(0)
      end
c.......................................................................
      subroutine buildHa(Ha,nat,ncont,nbf,basis,
     &       scr,nscr,printit)
      implicit none
c
#include "bas.fh"
c
      integer nat,ncont,nbf,basis,nscr
      double precision ha(nbf,nbf,3,nat),scr(nscr)
      logical printit
c
      integer used, iscr, newnscr, nbfmax
      integer ish, ibf, ibflo, ibfhi
      integer jsh, jbf, jbflo, jbfhi
      integer dxyz, datom, icount
c
      call dfill(nbf*nbf*3*nat,0.0,ha,1)
      call dfill(nscr,0.0,scr,1)
c
      if(.not.bas_nbf_cn_max(basis,nbfmax))
     &       call errquit('raktest:',922)
c
      used = nbfmax*nbfmax*3*nat
      iscr = used + 1
      newnscr = nscr - used
c
      do 00100 ish = 1,ncont
        do 00200 jsh = 1,ish
c
          call intd_1eh1(basis,ish,basis,jsh,newnscr,scr(iscr),
     &           used,scr)
c
          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &           call errquit('error',911)
          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &           call errquit('error',911)
c          
          icount = 0
          do 00300 datom = 1,nat
              do 00400 dxyz = 1,3
                do 00500 ibf = ibflo,ibfhi
                  do 00600 jbf = jbflo,jbfhi
                    icount = icount + 1
                    ha(ibf,jbf,dxyz,datom) = scr(icount)
                    ha(jbf,ibf,dxyz,datom) = scr(icount)
00600             continue
00500           continue
00400         continue
00300     continue
00200   continue
00100 continue

      if (.not.printit) return
      do 00900 datom = 1,nat
        do 01000 dxyz = 1,3
          write(6,'(//,a,i2,a,i2,a)')
     &           ' HaNEW <atom:',datom,'>  <xyz:',dxyz,'>'
          call myoutput(ha(1,1,dxyz,datom),1,nbf,1,nbf,nbf,nbf,1,
     &           ' HaNEW ',1.0d-05)
01000   continue
00900 continue
      end
c......................................................................
      subroutine build_d2e_new(basis,nbf,ncont,nat,intsize,d2e,
     &       maxg,gbuf,nscr,scr)
      implicit none
#include "bas.fh"
c::passed
      integer basis,nbf,nat,intsize,maxg,nscr,ncont
      double precision d2e(intsize,3,nat)
      double precision gbuf(maxg)
      double precision scr(nscr)
c::local
      integer idatom(4)
      integer dxyz,id,index,icount
      integer ish,ibf,ibflo,ibfhi
      integer jsh,jbf,jbflo,jbfhi
      integer ksh,kbf,kbflo,kbfhi
      integer lsh,lbf,lbflo,lbfhi
c
C     SET-UP FUNCTIONS
c
c...  statement function
C
C     SET-UP FUNCTIONS
      integer min,max
      intrinsic min,max
      integer i,j,k,l,isym2,isym4
      isym2(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
      isym4(i,j,k,l)=max(isym2(i,j),isym2(k,l))*
     &               (max(isym2(i,j),isym2(k,l))-1)/2+
     &               min(isym2(i,j),isym2(k,l))
c
c      
      do 00100 ish = 1,ncont
        do 00200 jsh= 1,ish
          do 00300 ksh = 1,jsh
            do 00400 lsh = 1,ksh
c
              call intd_2e4c(basis,ish,jsh,basis,ksh,lsh,
     &               nscr,scr,maxg,gbuf,idatom)
c
          if(.not.bas_cn2bfr(basis,ish,ibflo,ibfhi))
     &           call errquit('error',911)
          if(.not.bas_cn2bfr(basis,jsh,jbflo,jbfhi))
     &           call errquit('error',911)
          if(.not.bas_cn2bfr(basis,ksh,kbflo,kbfhi))
     &           call errquit('error',911)
          if(.not.bas_cn2bfr(basis,lsh,lbflo,lbfhi))
     &           call errquit('error',911)
c
c ... fill integral buffers              
                icount = 0
                do 00600 id = 1,4
                  if (idatom(id).gt.0) then
                    do 00700 dxyz = 1,3
                      do 00800 ibf = ibflo,ibfhi
                        do 00900 jbf = jbflo,jbfhi
                          do 01000 kbf = kbflo,kbfhi
                            do 01100 lbf = lbflo,lbfhi
                              index = isym4(ibf,jbf,kbf,lbf)
                              icount = icount + 1
                              d2e(index,dxyz,idatom(id)) = gbuf(icount)
01100                       continue
01000                     continue
00900                   continue
00800                 continue
00700               continue
                  else
                    icount = icount +
     &                     3*(ibfhi-ibflo+1)*(jbfhi-jbflo+1)*
     &                     (kbfhi-kbflo+1)*(lbfhi-lbflo+1)
                  endif
00600           continue
c
00400       continue
00300     continue
00200   continue
00100 continue
      end
