      subroutine raktest(rtdb)
      implicit none
      integer rtdb
c
      write(6,*)' run rjh test '
c      call rjhtest(rtdb)
      write(6,*)' rjh test done '
      write(6,*) ' rak test'
      call dointsum(rtdb,'geometry', 'mo basis')
c
      end
      subroutine dointsum(rtdb,geom_name, basis_name)
      implicit none
#include "geom.fh"
#include "bas.fh"
#include "mafdecls.fh"
      logical int_normalize
      external int_normalize
c
      integer rtdb
      integer basis, geom
      integer nscr, hscr, iscr
      integer nbf, nshells, nat
      integer buf1esz, hbuf1e, ibuf1e
      integer sa_size, hsa, isa, hsafd, isafd
      integer threen
      integer hxyz, ixyz
      integer ta_size
      integer va_size
      character*(*) geom_name, basis_name
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
c      write(6,*)' without integrals  printed '
      call int_chk_sum(basis,.false.)

c
      write(6,*)' derivative test '
c
      call hf_print_set(1)
c
      nscr = 1 000 000
      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
     &       call errquit('raktest: allocation scr failed',nscr)
      call dfill(nscr,0.0,dbl_mb(iscr),1)
      call onedint_fd(basis,geom,dbl_mb(iscr),nscr)
      if (.not.MA_Pop_Stack(hscr))
     &       write(6,*)'raktest: pop on hscr failed'
c
      
      if(.not.bas_numbf(basis,nbf))
     &       stop 'raktest:bas_numbf failed'
      if(.not.bas_numcont(basis,nshells))
     &       stop 'raktest:bas_numcont failed'
      if(.not.geom_ncent(geom,nat))
     &       stop 'raktest:geom_ncent failed'
c
      sa_size = (nbf*nbf)*(3*nat)
      ta_size = sa_size
      va_size = va_size
c..... space for sa
      if (.not.MA_Push_Get(MT_dbl,sa_size,'sa buffer',hsa,isa))
     &       call errquit('raktest: allocation for sa failed')
c..... space for safd
      if (.not.MA_Push_Get(MT_dbl,sa_size,'safd buffer',hsafd,isafd))
     &       call errquit('raktest: allocation for safd failed')
      threen = 3*nat
      if (.not.MA_Push_Get(MT_dbl,threen,'xyzxyz buffer',hxyz,ixyz))
     &       call errquit('raktest: allocation for safd failed')
      buf1esz = nbf*nbf*2
      if (.not.MA_Push_Get(MT_dbl,buf1esz,'safd buffer',hbuf1e,ibuf1e))
     &       call errquit('raktest: allocation for buf1e failed')
c..... scratch space
      nscr = 100 000
      if (.not.MA_Push_Get(MT_Dbl,nscr,'int buffer',hscr,iscr))
     &       call errquit('raktest: allocation scr failed',nscr)
      call buildSa(dbl_mb(isa),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,.true.)
      call buildSafd(dbl_mb(isafd),nat,nshells,nbf,basis,geom,
     &       dbl_mb(iscr),nscr,dbl_mb(ibuf1e),dbl_mb(ibuf1e+nbf*nbf),
     &       dbl_mb(ixyz),.true.)
      if (.not.MA_Pop_Stack(hscr))
     &       write(6,*)'raktest: pop on hscr failed'
      if (.not.MA_Pop_Stack(hbuf1e))
     &       write(6,*)'raktest: pop on hbuf1e failed'
      if (.not.MA_Pop_Stack(hxyz))
     &       write(6,*)'raktest: pop on hxyz failed'
      if (.not.MA_Pop_Stack(hsafd))
     &       write(6,*)'raktest: pop on hsafd failed'
      if (.not.MA_Pop_Stack(hsa))
     &       write(6,*)'raktest: pop on hsa failed'
c
      if (.not.(
     &       (bas_destroy(basis)).and.(geom_destroy(geom))
     &       ))
     &       call errquit('rak:error destroying gbs pair',0)
      end
