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

      open(1,file='1e')

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
     &                      exndcf(ica,basis),npa,nca,la,
     &   coords(1,ictrb,geom),exndcf(ipb,basis),
     &                      exndcf(icb,basis),npb,ncb,lb,
     &   coords(1,1,geom),charge(1,geom),ncenters,
     &                      ictra,ictrb,
     &          W0(i_Sa),W0(i_Ta),W0(i_Va),nint,
     &          .TRUE.,.TRUE.,.TRUE.,canAB,W0(iscr),nscr)
        
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
      if( itime.eq.1 )then
       write(ltime,9999) 'hfset    ',time(1)
       write(ltime,9999) 'hfmke    ',time(2)
       write(ltime,9999) 'O2Is     ',time(3)
       write(ltime,9999) 'KEIs     ',time(4)
       write(ltime,9999) 'NAIs     ',time(5)
      end if

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
     &             inp,igen,itype,
c
     &             coords(1,jatom,geom),
     &             exndcf(jexp,basis),
     &             exndcf(jcf,basis),
     &             jnp,jgen,jtype,
     &             iatom,jatom,
c
     &             coords(1,1,geom),charge(1,geom),nat,
     &             scr,scr,scr,nint,
     &             FT,FF,FF,FF,scr(iscr),(nscr-used))
            icount = 0
            ids = (iatom-1)*3+1
            ide = ids + 2
            do 00300 id = ids,ide
              do 00400 jbf = ibs_cn2bfr(1,jsh,basis),
     &               ibs_cn2bfr(2,jsh,basis)
                do 00500 ibf = ibs_cn2bfr(1,ish,basis),
     &                 ibs_cn2bfr(2,ish,basis)
                  icount = icount + 1
                  sa(ibf,jbf,id) = scr(icount)
00500           continue
00400         continue
00300       continue
            ids = (jatom-1)*3+1
            ide = ids + 2
            do 00600 id = ids,ide
              do 00700 jbf = ibs_cn2bfr(1,jsh,basis),
     &               ibs_cn2bfr(2,jsh,basis)
                do 00800 ibf = ibs_cn2bfr(1,ish,basis),
     &                 ibs_cn2bfr(2,ish,basis)
                  icount = icount + 1
                  sa(ibf,jbf,id) = scr(icount)
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
     &           ' <atom:',datom,'>  <xyz:',dxyz,'>'
          call output(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1)
01000   continue
00900 continue
      end
c.......................................................................
      subroutine buildSafd(Sa,nat,ncont,nbf,basisin,geom,
     &       scr,nscr,bplus,bminus,xyzxyz,printit)
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
      double precision xyzxyz(3,nat)
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
      delta = 1.0d-05
      call dfill(nbf*nbf*3*nat,0.0,sa,1)
      call dfill(nscr,0.0,scr,1)
c
      call dcopy(3*nat,coords(1,1,geom),1,xyzxyz,1)
c
      id = 0
      do 00050 datom = 1,nat
        do 00060 dxyz = 1,3
          id = id + 1
          do 00100 ish = 1,ncont
            do 00200 jsh = 1,ish
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
              used = nint
              iscr = used + 1
              call dcopy(3*nat,coords(1,1,geom),1,xyzxyz,1)
              xyzxyz(dxyz,datom) = xyzxyz(dxyz,datom) + delta
              call print_coords(' + coords',xyzxyz,nat)
              call hf1(
     &               xyzxyz(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyz(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyz,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00300 jbf = ibs_cn2bfr(1,jsh,basis),
     &               ibs_cn2bfr(2,jsh,basis)
                do 00400 ibf=ibs_cn2bfr(1,ish,basis),
     &                 ibs_cn2bfr(2,ish,basis)
                  icount = icount + 1
                  bplus(jbf,ibf) = scr(icount)
00400           continue
00300         continue
              call dfill(nscr,0.0,scr,1)
              call dcopy(3*nat,coords(1,1,geom),1,xyzxyz,1)
              xyzxyz(dxyz,datom) = xyzxyz(dxyz,datom) - delta
              call print_coords(' - coords',xyzxyz,nat)
              call hf1(
     &               xyzxyz(1,iatom),
     &               exndcf(iexp,basis),
     &               exndcf(icf,basis),
     &               inp,igen,itype,
c
     &               xyzxyz(1,jatom),
     &               exndcf(jexp,basis),
     &               exndcf(jcf,basis),
     &               jnp,jgen,jtype,
c
     &               xyzxyz,charge(1,geom),nat,
     &               scr,scr,scr,nint,
     &               FT,FF,FF,FF,FF,scr(iscr),(nscr-used))
              icount = 0
              do 00500 jbf = ibs_cn2bfr(1,jsh,basis),
     &               ibs_cn2bfr(2,jsh,basis)
                do 00600 ibf=ibs_cn2bfr(1,ish,basis),
     &                 ibs_cn2bfr(2,ish,basis)
                  icount = icount + 1
                  bminus(jbf,ibf) = scr(icount)
00600           continue
00500         continue
00200       continue
00100     continue
          write(6,*)' bplus  <nat:xyz>  <',datom,':',dxyz,'>'
          call output(bplus,1,nbf,1,nbf,nbf,nbf,1)
          write(6,*)' bminus <nat:xyz>  <',datom,':',dxyz,'>'
          call output(bminus,1,nbf,1,nbf,nbf,nbf,1)
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
     &           ' <<<FD>>><atom:',datom,'>  <xyz:',dxyz,'>'
          call output(sa(1,1,id),1,nbf,1,nbf,nbf,nbf,1)
01000   continue
00900 continue
      end
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
10000 format(80('-'))
10001 format(1x,'<atom:',i3,'>  <x:',1pd10.5,'>  <y:',1pd10.5,
     &       '>  <z:',1pd10.5,'>')
      end
c.......................................................................
