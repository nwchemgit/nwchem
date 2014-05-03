c
C> \brief Setup the COSMO cavity surface
C>
      subroutine hnd_cosset(rtdb,nat,c,radius)
      implicit none
c
c              ----- starting from -icosahedron- -----
c
c     pass, napex, nface, error =   0      12      20      20
c     pass, napex, nface, error =   1      42      80     100    0.4982
c     pass  napex, nface, error =   2     162     320     420    0.1848
c     pass  napex, nface, error =   3     642    1280    1700    0.0523
c     pass  napex, nface, error =   4    2562    5120    6820    0.0135
c     pass  napex, nface, error =   5   10242   20480   27300    0.0034
c
c              ----- starting from -octahedron-  -----
c
c     pass, napex, nface, error =   0       6       8       8
c     pass, napex, nface, error =   1      18      32      40    0.8075
c     pass  napex, nface, error =   2      66     128     168    0.4557
c     pass  napex, nface, error =   3     258     512     680    0.1619
c     pass  napex, nface, error =   4    1026    2048    2728    0.0451
c     pass  napex, nface, error =   5    4098    8192   10920    0.0116
c     pass  napex, nface, error =   6   16386   32768   43688    0.0029
c
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "stdio.fh"
c
      integer rtdb                 !< [Input] The RTDB handle
      integer nat                  !< [Input] The number of atoms
      double precision c(3,nat)    !< [Input] The atomic coordinates
      double precision radius(nat) !< [Input] The atomic radii
c
      double precision dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision adiag,dsurf,dvol,srfmol,volmol,ptspatm
      integer          lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
c
      integer   mxface, mxapex
      parameter (mxface=43688)
      parameter (mxapex=16386)
      logical     dbug, stat
      integer l_i10,  i10
      integer l_i20,  i20
      integer l_i30,  i30
      integer l_i40,  i40
      integer l_i50,  i50
      integer l_i60,  i60
      integer l_i70,  i70
      integer l_i80,  i80
      integer l_i90,  i90
      integer l_i100, i100
      integer l_i110, i110
      integer l_i120, i120
      integer l_i130, i130
      integer need
c
      dbug=.false.
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9999)
      endif
c
      if(ificos.eq.0.and.maxbem.gt.6) then
         write(luout,*) '-maxbem- too large for parameters in -cosset-'
         call errquit('hnd_cosset, -maxbem- too large',911,0)
      elseif(ificos.ne.0.and.maxbem.gt.7) then
         write(luout,*) '-maxbem- too large for parameters in -cosset-'
         call errquit('hnd_cosset, -maxbem- too large',911,0)
      endif
c
c     ----- partition memory -----
c
      need = 6*nat + 7*mxface + 7*mxface*nat + 3*mxapex
c
c     ----- allocate memory block -----
c
c     if(.not.ma_push_get(mt_dbl,need,'mem init:cosmo:hnd_cosset:1',
c    &    i_init,init))
c    & call errquit('hnd_cosset, malloc of init  failed',911,MA_ERR)
c
      stat = .true.
      stat = stat.and.ma_push_get(mt_dbl,3*nat,"xyzatm",l_i10,i10)
      stat = stat.and.ma_push_get(mt_dbl,  nat,"ratm",l_i20,i20)
      stat = stat.and.ma_push_get(mt_int,  nat,"nspa",l_i30,i30)
      stat = stat.and.ma_push_get(mt_int,  nat,"nppa",l_i40,i40)
      stat = stat.and.ma_push_get(mt_int,3*mxface,"ijkfac",l_i50,i50)
      stat = stat.and.ma_push_get(mt_dbl,3*mxface,"xyzseg",l_i60,i60)
      stat = stat.and.ma_push_get(mt_int,  mxface,"ijkseg",l_i70,i70)
      stat = stat.and.ma_push_get(mt_log,  mxface*nat,"insseg",
     &                                                l_i80,i80)
      stat = stat.and.ma_push_get(mt_dbl,3*mxface*nat,"xyzspa",
     &                                                l_i90,i90)
      stat = stat.and.ma_push_get(mt_int,  mxface*nat,"ijkspa",
     &                                                l_i100,i100)
      stat = stat.and.ma_push_get(mt_int,  mxface*nat,"numpps",
     &                                                l_i110,i110)
      stat = stat.and.ma_push_get(mt_dbl,3*mxapex    ,"apex",
     &                                                l_i120,i120)
      stat = stat.and.ma_push_get(mt_dbl,  mxface*nat,"xyzff",
     &                                                l_i130,i130)
c     i10 =init                    ! xyzatm(3,nat)
c     i20 =i10 +3*nat              !   ratm(  nat)
c     i30 =i20 +  nat              !   nspa(  nat)
c     i40 =i30 +  nat              !   nppa(  nat)
c     i50 =i40 +  nat              ! ijkfac(3,mxface)
c     i60 =i50 +3*mxface             ! xyzseg(3,mxface)
c     i70 =i60 +3*mxface             ! ijkseg(  mxface)
c     i80 =i70 +  mxface             ! insseg(  mxface,nat)
c     i90 =i80 +  mxface*nat         ! xyzspa(3,mxface,nat)
c     i100=i90 +3*mxface*nat         ! ijkspa(  mxface,nat)
c     i110=i100+  mxface*nat         ! numpps(  mxface,nat)
c     i120=i110+  mxface*nat         ! apex(3,mxapex)
c
c     ----- get -cosmo- surface -----
c
      call hnd_cossrf(nat,c,radius,nat,mxface,mxapex,
     1                dbl_mb(i10),dbl_mb(i20),int_mb(i30),int_mb(i40),
     2                int_mb(i50),dbl_mb(i60),int_mb(i70),
     3                log_mb(i80),dbl_mb(i90),int_mb(i100),int_mb(i110),
     4                dbl_mb(i120),dbl_mb(i130),rtdb)
c
c     ----- release memory block -----
c
      if(.not.ma_chop_stack(l_i10)) call
     &  errquit('hnd_cosset, ma_pop_stack of init failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- surface',
     2       /,10x,15(1h-))
      end
c
C> \brief Generate the COSMO cavity surface
C>
      subroutine hnd_cossrf(nat,c,radius,mxatm,mxfac,mxapx,
     1                  xyzatm,ratm,nspa,nppa,
     2                  ijkfac,xyzseg,ijkseg,insseg,
     3                  xyzspa,ijkspa,numpps,apex,xyzff,rtdb)
      implicit none
c
#include "nwc_const.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
      integer rtdb, nat
      integer mxatm
      integer mxfac
      integer mxapx
      double precision      c(3,nat  )
      double precision radius(    nat)
      double precision xyzatm(3,mxatm)
      double precision   ratm(  mxatm)
      integer            nspa(  mxatm)
      integer            nppa(  mxatm)
      integer          ijkfac(3,mxfac)
      double precision xyzseg(3,mxfac)
      integer          ijkseg(  mxfac)
      logical          insseg(  mxfac,mxatm)
      double precision xyzspa(3,mxfac,mxatm)
      integer          ijkspa(  mxfac,mxatm)
      integer          numpps(  mxfac,mxatm)
      double precision   apex(3,mxapx)
      double precision  xyzff(  mxfac,mxatm)
c
      double precision bohr
      parameter (bohr=0.529177249d+00)
      logical    some
      logical    out
      logical    more
      logical    dbug
      double precision dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision adiag,dsurf,dvol,srfmol,volmol,ptspatm
      integer          lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
c
      integer i, iat, lfac, lseg, ndiv, nfac, nseg
      integer mfac
c
      dbug=.false.
      more=.false.
      more=more.or.dbug
      out =.false. 
      out =out.or.more
      some=.false.
      some=some.or.out
c
c     ----- approximate sphere with segments and points -----
c
      do iat = 1, mxatm
        nspa(iat) = 0
        nppa(iat) = 0
      enddo
      nseg = 0
      nfac = 0
      ndiv = 0
      call hnd_cossph(nseg,nfac,ndiv,
     1                ijkfac,xyzseg,ijkseg,mxfac,apex,mxapx,
     2                dsurf,dvol,adiag)
      ptspatm = dble(nseg)
c
c     ----- debug printing -----
c
      if(out.and.ga_nodeid().eq.0) then
         write(luout,9999) nseg,nfac,ndiv,dsurf,dvol
         write(luout,9995) adiag
         if(more) then
            write(luout,9998)
            do lseg=1,nseg
               write(luout,9997) lseg,xyzseg(1,lseg),xyzseg(2,lseg),
     1                             xyzseg(3,lseg),ijkseg(  lseg)
            enddo
         endif
         if(dbug) then
            write(luout,9996)
            do lfac=1,nfac
               mfac=lfac+nseg
               write(luout,9997) mfac,xyzseg(1,mfac),xyzseg(2,mfac),
     1                             xyzseg(3,mfac),ijkseg(  mfac)
            enddo
         endif
      endif
c
c     ----- set molecule -----
c
      do iat=1,nat
         do i=1,3
            xyzatm(i,iat)=c(i,iat)
         enddo
      enddo
      do iat=1,nat
         if(radius(iat).eq.0.0d0) then
            ratm(iat)=0.0d0
         else
            ratm(iat)=radius(iat)/bohr
         endif
      enddo
c
c     ----- create -solvent accessible surface- of the molecule -----
c
      call hnd_cossas(nat,xyzatm,ratm,mxatm,
     1                nspa,nppa,xyzspa,ijkspa,
     2                nseg,nfac,xyzseg,ijkseg,insseg,
     3                numpps,xyzff,mxfac,rtdb)
c
      return
 9999 format(' nseg,nfac,ndiv=nfac/nseg,dsurf,dvol = ',3i7,2f10.6)
 9998 format('  seg  ','      x     ','      y     ','      z     ',
     1       ' seg ',/,1x,47(1h-))
 9997 format(i7,3f12.8,i5,f12.8)
 9996 format('  fac  ','      x     ','      y     ','      z     ',
     1       ' seg ',/,1x,47(1h-))
 9995 format(' adiag           = ',f12.6)
      end
C>
C> \brief Construct the Solvent Accessible Surface (SAS) from the
C> triangulated spheres
C>
C> ## The legacy of Klamt and Sch&uuml;&uuml;rmann ##
C>
C> This subroutine was originally written to implement the algorithm
C> to construct the Solvent Accessible Surface as proposed by 
C> Klamt and Sch&uuml;&uuml;rmann [1]. This algorithm worked as follows:
C>
C> If two spheres partially overlap then parts of the surface need
C> to be eliminated. Segments that are contained entirely within the
C> sphere of another atom will be eliminated completely. Segments
C> that straddle the boundary between two spheres will have their
C> surface reduced proportional to the fraction that resides within the
C> sphere of the other atom. This fraction is established by counting
C> the number of faces that fall within the sphere of the other atom. 
C> In addition the location of the charge representing a segment should
C> be calculated as the center of the remaining points representing the
C> faces (see [1] page 802, 2nd column, 2nd paragraph), but currently
C> that is not done.
C>
C> To understand the approach suggested above it is important to know
C> the concepts "segments" and "faces". 
C> - "Segments" are parts of the surface of the sphere that will be
C>   represented by a single COSMO charge.
C> - "Faces" are further refinements of segments. I.e. segments have
C>   been partitioned into a number of faces. The faces are used to
C>   eliminate parts of a segment that are within the sphere of another
C>   atom. By counting the remaining faces the surface area of a segment
C>   can be adjusted.
C> The trick with segments and faces is needed to create a smoother
C> boundary between neighboring spheres without introducing large 
C> numbers of COSMO charges. The smooth boundary is needed to keep
C> discretization errors small, whereas "small" numbers of COSMO charges
C> are needed to keep the cost of calculating the COSMO charges low.
C>
C> The segments have been generated in `hnd_cossph` by partitioning the
C> triangles of the original polyhedron `minbem` times. The faces have
C> generated by partitioning the segments an additional `maxbem-minbem`
C> times.
C>
C> ## The new smooth approach of York and Karplus ##
C>
C> The approach by Klamt and Sch&uuml;&uuml;rmann [1] led to problems
C> because the corresponding potential energy surface was not 
C> continuous. York and Karplus [2] proposed a method that provides
C> a smooth potential energy surface and this subroutine was changed
C> to implement this new approach. This meant that some things stayed
C> the same as before (for example minbem still works the same way
C> to generate the surface charge positions), other things changed
C> significantly (maxbem and rsolve are now entirely obsolete and not
C> used anymore, also the elimination of point charges is no longer
C> based on reducing the segment surface area until it vanishes, instead
C> the surface charge of a segment is quenched with a switching function
C> to eliminate the contribution of a surface point).
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="http://dx.doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="http://dx.doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
      subroutine hnd_cossas(nat,xyzatm,ratm,mxatom,
     1                      nspa,nppa,xyzspa,ijkspa,
     2                      nseg,nfac,xyzseg,ijkseg,insseg,
     3                      numpps,xyzff,mxface,rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bq.fh"
#include "prop.fh"
c
      integer rtdb    !< [Input] The RTDB handle
      integer nat     !< [Input] The actual number of atoms
      integer mxface  !< [Input] The maximum number of faces
      integer mxatom  !< [Input] The maximum number of atoms
      integer nseg    !< [Input] The actual number of segments
      integer nfac    !< [Input] The actual number of faces
c
      double precision bohr
      parameter  (bohr=0.529177249d+00)
      logical     dbug
      double precision dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision adiag,dsurf,dvol,srfmol,volmol,ptspatm
      integer lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
      double precision xyzatm(3,mxatom) !< [Input] The atom positions
      double precision   ratm(  mxatom) !< [Input] The atom radii
      integer            nspa(  mxatom) !< [Output] The number of
                                        !< segments remaining for
                                        !< each atom
      integer            nppa(  mxatom) !< [Output] The number of faces
                                        !< remaining for each atom
      double precision xyzseg(3,mxface) !< [Input] The coordinates of 
                                        !< the segment and face points
                                        !< on the unit sphere
      integer          ijkseg(  mxface) !< [Input] List for every
                                        !< face what the corresponding
                                        !< segment is, if ijkseg(ii) is
                                        !< 0 then face ii should be 
                                        !< ignored (has been eliminated)
      logical          insseg(  mxface,mxatom) !< [Output] If .false.
                                               !< keep the segment or
                                               !< face, discard it
                                               !< otherwise
      double precision xyzspa(3,mxface,mxatom)
      integer          ijkspa(  mxface,mxatom)
      integer          numpps(  mxface,mxatom)
      double precision  xyzff(  mxface,mxatom)
      double precision zero, one
      data zero    /0.0d+00/
      data one     /1.0d+00/
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer l_efclb, k_efclb, k_efciat, l_efciat
      double precision ratm_real,dij,dum,cavdsp,pi,zetai,zetaii
      integer m,mfac,mseg
      integer nefc,iat,jat,npp,i,iseg,ifac,ief,ipp
c
      double precision cosff
      external         cosff
c
      double precision dist, xi, yi, zi, xj, yj, zj, rin, rout, alphai
      parameter (alphai = 0.5d0)
      dist(xi,yi,zi,xj,yj,zj)=sqrt((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
      rin(iat)=ratm(iat)*(1.0d0-alphai*gammas*sqrt(0.25d0**minbem))
      rout(iat)=ratm(iat)*(1.0d0+(1.0d0-alphai)*gammas*
     &                     sqrt(0.25d0**minbem))
c
      dbug=.false.
      pi = acos(-1.0d0)
c
      if(ga_nodeid().eq.0) then
         write(luout,9999)
      endif
c
c     ----- print atomic centers -----
c
      if(ga_nodeid().eq.0) then
         write(luout,9998)
         do iat=1,nat
            write(luout,9997) iat,xyzatm(1,iat),xyzatm(2,iat),
     1                                       xyzatm(3,iat),
     2                    (ratm(iat)*bohr)
         enddo
      endif
c
c     ----- clear arrays ..... -----
c
      do iat=1,nat
         do i=1,mxface
            ijkspa(i,iat)=0
            numpps(i,iat)=0
            xyzff(i,iat)=zero
         enddo
      enddo
c
c     ----- sift through atomic centers and decide if a face -----
c           belongs to the -sas- or is inside the molecule.
c
      do iat=1,nat
c
         if(ratm(iat).ne.zero) then
            do iseg=1,nseg
               ijkspa(iseg,iat)=ijkseg(iseg)
               xyzff(iseg,iat)=one
               do m=1,3
                  xyzspa(m,iseg,iat)=xyzseg(m,iseg)*ratm(iat)
     1                              +xyzatm(m,iat)
               enddo
            enddo
            do ifac=1,nfac
               ijkspa(ifac+nseg,iat)=ijkseg(ifac+nseg)
               do m=1,3
                  xyzspa(m,ifac+nseg,iat)=xyzseg(m,ifac+nseg)*ratm(iat)
     1                                   +xyzatm(m,iat)
               enddo
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9996) iat
               write(luout,9995) (ijkspa(ifac+nseg,iat),ifac=1,nfac)
            endif
            do jat=1,nat
               dij=dist(xyzatm(1,iat),
     1                  xyzatm(2,iat),
     2                  xyzatm(3,iat),
     3                  xyzatm(1,jat),
     4                  xyzatm(2,jat),
     5                  xyzatm(3,jat))
               if((jat.ne.iat).and.(ratm(jat).ne.zero)
     1                        .and.(dij.lt.(ratm(iat)+rout(jat)))) then
                  do iseg=1,nseg
                     dum=dist(xyzspa(1,iseg,iat),
     1                        xyzspa(2,iseg,iat),
     2                        xyzspa(3,iseg,iat),
     3                        xyzatm(1,jat),
     4                        xyzatm(2,jat),
     5                        xyzatm(3,jat))
                     xyzff(iseg,iat) = xyzff(iseg,iat) *
     1                 cosff((dum-rin(jat))/(rout(jat)-rin(jat)))
                  enddo
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9996) iat
               write(luout,9995) (ijkspa(ifac+nseg,iat),ifac=1,nfac)
            endif
c
c     ----- check segments belonging to -sas- -----
c
            do iseg=1,nseg
               insseg(iseg,iat)=.not.(xyzff(iseg,iat).ge.swtol)
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9994) iat
               write(luout,9993) (insseg(iseg,iat),iseg=1,nseg)
            endif
            mseg=0
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) mseg=mseg+1
            enddo
            mfac=0
            nspa(iat)=mseg
            nppa(iat)=mfac
c
         endif
c
      enddo
c
      if(ga_nodeid().eq.0) then
         write(luout,9985) nseg,nfac
         write(luout,9992)
         do iat=1,nat
            npp=0
            do iseg=1,nseg
               npp=npp+numpps(iseg,iat)
            enddo
            write(luout,9991) iat,nspa(iat),nppa(iat),npp
         enddo
      endif
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9987)
         do iat=1,nat
            do iseg=1,nseg
               write(luout,9986) iat,iseg,numpps(iseg,iat)
            enddo
         enddo
      endif
c
c    Count the number of surface points, i.e. number of point charges
c    and generate memory to store them
c
      nefc = 0
      do iat=1,nat
         if(ratm(iat).ne.zero) then
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) nefc = nefc+1
            enddo
         endif
      enddo
c
c     Allocate memory for point charges
c
      if(.not.ma_push_get(mt_dbl,nefc*3,'cosmo efcc',l_efcc,k_efcc))
     & call errquit('cosmo_cossas malloc k_efcc failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcs',l_efcs,k_efcs))
     & call errquit('cosmo_cossas malloc k_efcs failed',911,MA_ERR)
      if(.not.ma_push_get(mt_int,nefc,'cosmo efciat',l_efciat,k_efciat))
     & call errquit('cosmo_cossas malloc k_efciat failed',911,MA_ERR)
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo efcz',l_efcz,k_efcz))
     & call errquit('cosmo_cossas malloc k_efcz failed',911,MA_ERR)
      if(.not.ma_push_get(mt_byte,nefc*8,'cosmo tags',l_efclb,k_efclb))
     & call errquit('cosmo_cossas malloc k_tag failed',911,MA_ERR)
c
c     ----- save coordinates of surface points -----
c           save segment surfaces
c           save segment to atom mapping
c
      srfmol=zero
      volmol=zero
      ief   =0
      do iat=1,nat
         if(ratm(iat).ne.zero) then
            ratm_real=ratm(iat)
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) then
                  ief=ief+1
                  do i=1,3
                     dbl_mb(k_efcc+3*(ief-1)+i-1)=xyzatm(i,iat)
     1                          +xyzseg(i,iseg)*ratm_real
                  enddo
c
c                 --- eval eq.(67) from [2] ---
c
                  ipp=numpps(iseg,iat)
                  dum=4.00d0**(maxbem-minbem)
                  dum=16.0d0 ! MAXBEM is obsolete in York and Karplus approach
                  zetai=zeta*sqrt(nseg*dum)/(ratm_real*sqrt(2.0d0*pi))
                  zetaii=zetai/sqrt(2.0d0)
                  dbl_mb(k_efcs+ief-1) = (1.0d0/xyzff(iseg,iat))
     1                                 * 2.0d0*zetaii/sqrt(pi)
                  srfmol   = srfmol + xyzff(iseg,iat)*dsurf*ratm_real**2
                  volmol   = volmol + xyzff(iseg,iat)*dvol *ratm_real**3
                  int_mb(k_efciat+ief-1)=iat
               endif
            enddo
         endif
      enddo
      srfmol=srfmol*(bohr**2)
      volmol=volmol*(bohr**3)
c
      if(ga_nodeid().eq.0) then
         write(luout,9990) nefc
         write(luout,9984) srfmol 
         write(luout,9983) volmol 
      endif 
c
c     ----- Cavity/Dispersion free energy ---
c           Sitkoff, Sharp, and Honig,
c           J.Phys.Chem. 98, 1978 (1994)
c
      cavdsp=0.860+0.005*srfmol
      if(ga_nodeid().eq.0) then
         write(luout,9981) cavdsp
      endif
c
c     ----- print segment surfaces -----
c
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9989)
         do ief=1,nefc
            write(luout,9988) ief,dbl_mb(k_efcs+ief-1),
     &                        int_mb(k_efciat+ief-1)
         enddo
      endif
c
      do ief=1,nefc
         dbl_mb(k_efcz+ief-1)=zero
      enddo
      do ief=1,nefc
         byte_mb(k_efclb+(ief-1)*8)='        '
      enddo
c
c     ----- write out to -rtdb- -----
c
      if(.not.rtdb_put(rtdb,'cosmo:nefc',mt_int,1     ,nefc))
     $   call errquit('hnd_cossas: rtdb put failed for nefc  ',911,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcc',mt_dbl,3*nefc,dbl_mb(k_efcc)))
     $   call errquit('hnd_cossas: rtdb put failed for efcc  ',912,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,dbl_mb(k_efcz)))
     $   call errquit('hnd_cossas: rtdb put failed for efcz  ',913,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcs',mt_dbl,  nefc,dbl_mb(k_efcs)))
     $   call errquit('hnd_cossas: rtdb put failed for efcs  ',914,
     &       rtdb_err)
c
c     ----- reset cosmo:rawt to avoid trouble in cosmo charge 
c           calculation -----
c
      if(.not.rtdb_put(rtdb,'cosmo:rawt',mt_dbl,  nefc,dbl_mb(k_efcz)))
     $   call errquit('hnd_cossas: rtdb put failed for rawt  ',915,
     &       rtdb_err)
c
c     We will need the next bit of information to calculate the analytic
c     COSMO gradients. This table describes the relationship between
c     the COSMO charges and the associated atoms. So we better save this
c     info.
c
      if(.not.rtdb_put(rtdb,'cosmo:efciat',mt_int,nefc,
     $                 int_mb(k_efciat)))
     $   call errquit('hnd_cossas: rtdb put failed for iatefc',916,
     &       rtdb_err)
c     if(.not.rtdb_cput(rtdb,'char variable',nefc,byte_mb(k_efclb)))
c    $   call errquit('hnd_cossas: rtdb put failed for efclab',917,
c    &       rtdb_err)
c
      if(.not.ma_pop_stack(l_efclb)) call
     &   errquit('cosmo_cossas chop stack k_efclb failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcz)) call
     &   errquit('cosmo_cossas chop stack k_efcz failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efciat)) call
     &   errquit('cosmo_cossas chop stack k_efciat failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcs)) call
     &   errquit('cosmo_cossas chop stack k_efcs failed',911,MA_ERR)
      if(.not.ma_pop_stack(l_efcc)) call
     &   errquit('cosmo_cossas chop stack k_efcc failed',911,MA_ERR)
c
      return
 9999 format(/,1x,'solvent accessible surface',/,1x,26(1h-))
 9998 format(/,1x,'---------- ATOMIC COORDINATES (A.U.) ----------',
     1            '-- VDWR(ANG.) --')
 9997 format(  1x,i5,3f14.8,f10.3)
 9996 format(/,1x,'---------- SEGMENTS FOR -IAT- = ',i5)
 9995 format(16i4)
 9994 format(/,1x,'-INSSEG- FACES FOR IAT = ',i5)
 9993 format(16l4)
 9992 format(  1x,'atom',' ( ','  nspa',',','  nppa',' )',/,1x,22(1h-))
 9991 format(  1x,   i4 ,' ( ',     i6 ,',',     i6 ,' )',i8)
 9990 format(  1x,'number of -cosmo- surface points = ',i8)
 9989 format(  1x,'SEGMENT SURFACES =',/,1x,18(1h-))
 9988 format(i8,f10.5,i5)
 9987 format(  1x,'NUMBER OF FACES / SEGMENT =',/,1x,27(1h-))
 9986 format(3i5)
 9985 format(' number of segments per atom = ',i10,/,
     1       ' number of   points per atom = ',i10)
 9984 format(' molecular surface = ',f10.3,' angstrom**2')
 9983 format(' molecular volume  = ',f10.3,' angstrom**3')
 9981 format(' G(cav/disp)       = ',f10.3,' kcal/mol')
      end
c
C> \brief Triangulate a sphere using the Boundary Element Method (BEM)
C>
C> This routine approximates a sphere starting from either an
C> octahedron or an icosahedron and partitioning the triangles that
C> make up these polyhedra. Each triangle is partitioned into four
C> triangles at each level in the recursion. The procedure is starting
C> from an equal sided triangle, select the midpoints of all three
C> sides, and draw a triangle through the three midpoints. This way four
C> triangles of equal size are obtained. However, the midpoints of the
C> original sides do not lie on the surface of the sphere and therefore
C> they need to be projected outwards. This outwards projection changes
C> the size of the central triangle more than that of the other three.
C> So in the final triangulation the triangles are not all of the same
C> size, but this is ignored in the COSMO formalism.
C>
C> Ultimately we are interested only in the triangles at 2 levels of
C> granularity:
C>
C> - minbem: these triangles are referred to as "segments" and
C>   represent the sphere and their centers become the positions for
C>   the COSMO charges.
C>
C> - maxbem: these triangles are referred to as "faces" and they are 
C>   used to adjust the surface of the segments in regions where two
C>   atomic spheres meet and the segments straddle the boundary between
C>   both spheres.
C>
C> All other triangles are reduced to mere artefacts of the triangle
C> generation algorithm. The array `ijkseg` administrates what the 
C> status of a triangle is. It lists for each face which segment it is
C> part of.
C>
C> In addition this routine computes `adiag` of [1] Eq.(B1). The 
C> expression in this routine can be obtained from Eq.(B1) as
C> \f{eqnarray*}{
C>   \frac{1}{2R}
C>   &=&\frac{M}{2}\sum_{\nu=2}^M\frac{M^{-2}}{||t_1-t_\nu||}
C>    + MM^{-2}\frac{a_{\mathrm{diag}}}{2} \\\\
C>   \frac{a_{\mathrm{diag}}}{M}
C>   &=&\frac{1}{R}-\sum_{\nu=2}^M\frac{M^{-1}}{||t_1-t_\nu||} \\\\
C>   a_{\mathrm{diag}}
C>   &=&\frac{M}{R}-\sum_{\nu=2}^M\frac{1}{||t_1-t_\nu||}
C> \f}
C> The expression implemented in this routine can be mapped onto this
C> by
C> \f{eqnarray*}{
C>   a_{\mathrm{diag}}' &=& \left(\frac{4\pi}{M}\right)^{1/2}
C>         \left(M-\sum_{\nu=2}\frac{1}{||t_1'-t_\nu'||}\right)
C> \f}
C> where \f$a_{\mathrm{diag}}'\f$ is `adiag` as calculated in this
C> routine, the \f$t'\f$ are points in the unit sphere (as opposed to
C> \f$t\f$ which are points on the sphere with radius \f$R\f$).
C> In the routines `hnd_coschg`, `hnd_cosaxd` and `hnd_cosxad` this is
C> multiplied with 
C> \f$|S_\mu|^{-1/2} = \left(\frac{4\pi R^2}{M}\right)^{-1/2}\f$ to give
C> the proper \f$a_{\mathrm{diag}}\f$
C> \f{eqnarray*}{
C>   a_{\mathrm{diag}} &=& \left(\frac{4\pi}{M}\right)^{1/2}
C>         \left(M-\sum_{\nu=2}\frac{1}{||t_1'-t_\nu'||}\right)
C>         \left(\frac{M}{4\pi R^2}\right)^{1/2} \\\\
C>   &=& \frac{1}{R}
C>       \left(M-\sum_{\nu=2}\frac{1}{||t_1'-t_\nu'||}\right) \\\\
C>   &=& \left(\frac{M}{R}-\sum_{\nu=2}\frac{1}{||t_1-t_\nu||}\right)
C> \f}
C> In ref.[1] Eq.(B2) is wrong because of a spurious scale factor
C> \f$4\pi R^2/M\f$.
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="http://dx.doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
      subroutine hnd_cossph(nseg,nfac,ndiv,
     1                  ijkfac,xyzseg,ijkseg,mxface,apex,mxapex,
     2                  dsurf,dvol,adiag)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
c              ----- starting from -icosahedron- -----
c
c     pass, napex, nface, error =   0      12      20      20
c     pass, napex, nface, error =   1      42      80     100    0.4982
c     pass  napex, nface, error =   2     162     320     420    0.1848
c     pass  napex, nface, error =   3     642    1280    1700    0.0523
c     pass  napex, nface, error =   4    2562    5120    6820    0.0135
c     pass  napex, nface, error =   5   10242   20480   27300    0.0034
c
c              ----- starting from -octahedron-  -----
c
c     pass, napex, nface, error =   0       6       8       8
c     pass, napex, nface, error =   1      18      32      40    0.8075
c     pass  napex, nface, error =   2      66     128     168    0.4557
c     pass  napex, nface, error =   3     258     512     680    0.1619
c     pass  napex, nface, error =   4    1026    2048    2728    0.0451
c     pass  napex, nface, error =   5    4098    8192   10920    0.0116
c     pass  napex, nface, error =   6   16386   32768   43688    0.0029
c
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      dimension xyzseg(3,*)
c
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
c
      parameter (mxpass=    7)
      dimension minfac(mxpass)
      dimension maxfac(mxpass)
      dimension minico(mxpass)
      dimension maxico(mxpass)
      dimension minoct(mxpass)
      dimension maxoct(mxpass)
      dimension ijknew(3)
      dimension ijkold(3)
      equivalence (ijkold(1),iold),(ijkold(2),jold),(ijkold(3),kold)
      equivalence (ijknew(1),inew),(ijknew(2),jnew),(ijknew(3),knew)
      logical icos
      logical octa
      logical some,out,dbug
      data minico /    1,   21,  101,  421, 1701, 6821,    0/
      data maxico /   20,  100,  420, 1700, 6820,27300,    0/
      data minoct /    1,    9,   41,  169,  681, 2729,10921/
      data maxoct /    8,   40,  168,  680, 2728,10920,43688/
      data zero  /0.0d+00/
      data one   /1.0d+00/
      data two   /2.0d+00/
      data three /3.0d+00/
      data four  /4.0d+00/
c
      dist(xi,yi,zi,xj,yj,zj)=sqrt((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
c
      dbug=.false.
      out =.false.
      out =out.or.dbug
      some=.false.
      some=some.or.out
c
      pi=four*atan(one)
      rad=one
      cir= two*pi*rad
      srf=four*pi*rad**2
      vol=four*pi*rad**3/three
c
      npass=maxbem
c
c     ----- define  hedron  -----
c           define -minfac- 
c           define -maxfac- 
c
      icos=ificos.ne.0
      octa=.not.icos
      if(icos) then
         call hnd_sphico(apex,napex,ijkfac,ijkseg,nface)
         do ipass=1,mxpass
            minfac(ipass)=minico(ipass)
            maxfac(ipass)=maxico(ipass)
         enddo
      endif
      if(octa) then
         call hnd_sphoct(apex,napex,ijkfac,ijkseg,nface)
         do ipass=1,mxpass
            minfac(ipass)=minoct(ipass)
            maxfac(ipass)=maxoct(ipass)
         enddo
      endif
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         if(icos) then
            write(luout,9994)
         endif
         if(octa) then
            write(luout,9982)
         endif
         if(out) then
            write(luout,*) '-minbem- = ',minbem
            write(luout,*) '-maxbem- = ',maxbem
            write(luout,*) '-minfac- = ',minfac
            write(luout,*) '-maxfac- = ',maxfac
            write(luout,*) '-npass - = ',npass
            write(luout,9999)
            do iapex=1,napex
               write(luout,9998) iapex,apex(1,iapex),
     1                              apex(2,iapex),
     2                              apex(3,iapex)
            enddo
         endif
      endif
c
c     ----- loop over divisions to create sphere -----
c
      mxfac=0
      ipass=1
  100 ipass=ipass+1
         mnfac=mxfac+1
         mxfac=nface
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9996) ipass,napex,nface,mnfac,mxfac
         endif
c
         dmin =one
         mapex=napex
         mface=nface
         do lface=mnfac,mxfac
            iold=ijkfac(1,lface)
            jold=ijkfac(2,lface)
            kold=ijkfac(3,lface)
            call hnd_sphapx(apex,mapex,ijkfac,ijkseg,mface,lface,
     1                      ijkold,ijknew,dijk)
            dmin=min(dmin,dijk)
         enddo
         napex=mapex
         nface=mface
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9995) napex,nface
         endif
c
c     ----- print out new apeces -----
c
         if(dbug.and.ga_nodeid().eq.0) then
            do iapex=1,napex
               write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                              apex(3,iapex)
            enddo
         endif
c
c     ----- print approximate volume -----
c
         radapp=    dmin
         radrat=    dmin
         raderr=one-radrat
         srfapp=srf*dmin**2
         srfrat=    dmin**2
         srferr=one-srfrat
         volapp=vol*dmin**3
         volrat=    dmin**3
         volerr=one-volrat
         if(out.and.ga_nodeid().eq.0) then
            write(luout,9997) vol,volapp,volrat,volerr
            write(luout,9992) srf,srfapp,srfrat,srferr
            write(luout,9991) rad,radapp,radrat,raderr
         endif
c
c     ----- assign early segment to each face -----
c
         if(ipass.gt.(minbem+1)) then
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9981) ipass
               write(luout,9980) (minfac(i),i=1,ipass)
               write(luout,9979) (maxfac(i),i=1,ipass)
            endif
            maxseg=maxfac(minbem)
            lfacmn=minfac(ipass)
            lfacmx=maxfac(ipass)
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9990) ipass
               write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
            endif
            do lface=lfacmn,lfacmx
               ijkseg(lface)=ijkseg(ijkseg(lface))
               if(ijkseg(lface).gt.maxseg.and.ga_nodeid().eq.0) then
                  write(luout,9987) lface,ijkseg(lface)
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9989) ipass
               write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
            endif
         endif
c
      if(ipass.lt.npass) go to 100
c
c     ----- end of loop over tessalating passes -----
c
      if(dbug.and.ga_nodeid().eq.0) then
         do ipass=1,npass
            lfacmn=minfac(ipass)
            lfacmx=maxfac(ipass)
            write(luout,9989) ipass
            write(luout,*) '-lfacmn- = ',lfacmn
            write(luout,*) '-lfacmx- = ',lfacmx
            write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
         enddo
      endif
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9993) npass,napex,minfac(npass),maxfac(npass),
     1                  radapp,raderr,srfapp,srferr,volapp,volerr
      endif
c
c     ----- at this point each of the faces is assigned to one -----
c           segment. now define centers of segments ...
c
      third =one/three
      lfacmn= minfac(minbem)
      lfacmx= maxfac(minbem)
      do lface=lfacmn,lfacmx
         mface=lface-lfacmn+1
         ijkseg(mface)=mface
         i=ijkfac(1,lface)
         j=ijkfac(2,lface)
         k=ijkfac(3,lface)
         do m=1,3
            xyzseg(m,mface)=(apex(m,i)+apex(m,j)+apex(m,k))*third
         enddo
         dseg=one/dist(xyzseg(1,mface),xyzseg(2,mface),xyzseg(3,mface),
     1                 zero,zero,zero)
         do m=1,3
            xyzseg(m,mface)=xyzseg(m,mface)*dseg
         enddo
      enddo
      nseg=(lfacmx-lfacmn+1)
c
      if(dbug.and.ga_nodeid().eq.0) then
         lfacmn=1
         lfacmx=nseg
         write(luout,*)    'segment to segment mapping = '
         write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
      endif
c
c     ----- now the faces ... -----
c
      if(npass.gt.minbem) then
         lfacmn=minfac(minbem+1)
         lfacmx=maxfac(npass   )
         do lface=lfacmn,lfacmx
            mface=lface-lfacmn+1    
     1                        +(maxfac(minbem)-minfac(minbem)+1)
            ijkseg(mface)=ijkseg(lface)
     1                        -(               minfac(minbem)-1)
            i=ijkfac(1,lface)
            j=ijkfac(2,lface)
            k=ijkfac(3,lface)
            do m=1,3
               xyzseg(m,mface)=(apex(m,i)+apex(m,j)+apex(m,k))*third
            enddo
            dseg=one/dist(xyzseg(1,mface),
     1                    xyzseg(2,mface),
     2                    xyzseg(3,mface),zero,zero,zero)
            do m=1,3
               xyzseg(m,mface)=xyzseg(m,mface)*dseg
            enddo
         enddo
         nfac=(lfacmx-lfacmn+1)
c
c        ----- only keep the faces at granularity maxbem -----
c        ----- discard all other faces -----
c
         do lface=1,maxfac(maxbem-1)-minfac(minbem+1)+1
           ijkseg(nseg+lface) = 0
         enddo
      else
         do iseg=1,nseg
            ifac=iseg+nseg
            ijkseg(ifac)=ijkseg(iseg)
            do m=1,3
               xyzseg(m,ifac)=xyzseg(m,iseg)
            enddo
         enddo
         nfac=nseg
      endif
c
      if(dbug.and.ga_nodeid().eq.0) then
         lfacmn=nseg+1
         lfacmx=nseg+nfac
         write(luout,*)    ' face   to segment mapping = '
         write(luout,9988) (ijkseg(lface),lface=lfacmn,lfacmx)
      endif
c
c     ----- calculate -dsurf dvol- for the -cosmo- theory -----
c
      nfac =nseg
      ndiv =nfac/nseg
      dsurf=srf/dble(nfac)
      dvol =vol/dble(nfac)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9986) nseg,nfac,ndiv,dsurf,dvol
      endif
      if(out.and.ga_nodeid().eq.0) then
         write(luout,9985)
         do i=1,nseg
            done=dist(xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     1                zero,zero,zero)
            write(luout,9984) i,
     1                     xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     2                     ijkseg(i),done
         enddo
      endif
      if(dbug.and.ga_nodeid().eq.0) then
         write(luout,9985)
         do i=nseg+1,nseg+nfac
            done=dist(xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     1                zero,zero,zero)
            write(luout,9984) (i-nseg),
     1                     xyzseg(1,i),xyzseg(2,i),xyzseg(3,i),
     2                     ijkseg(i),done
         enddo
      endif
c
c     ----- calculate -adiag- of the -cosmo- theory -----
c
      avgdia=zero
      avgfac=zero
      do mseg=1,nseg
         sum=zero
         do lseg=1,nseg
            if(lseg.ne.mseg) then
               l1=mseg
               l2=lseg
         sum=sum+rad/dist(xyzseg(1,l2),xyzseg(2,l2),xyzseg(3,l2),
     1                    xyzseg(1,l1),xyzseg(2,l1),xyzseg(3,l1))
            endif
         enddo
         fac=(dble(nseg)-sum)/sqrt(dble(nseg))
         adiag=sqrt(four*pi)*fac
         if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
            write(luout,9983) mseg,adiag,fac,dble(nseg),sum
         endif
         avgdia=avgdia+adiag
         avgfac=avgfac+fac
      enddo
      adiag=avgdia/dble(nseg)
      fac  =avgfac/dble(nseg)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9978)      adiag,fac               
      endif
c
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(' vol, approx., ratio, error = ',2f12.8,2 f8.4)
 9996 format(' pass, napex, nface, mnfac, mxfac = ',i3,4i8)
 9995 format('       napex, nface               = ',3x,2i8)
 9994 format(1x,'sphere from -icosahedron-',/,1x,25(1h-))
 9993 format(' npass = ',i2,' napex = ',i8,
     1       ' minfac = ',i8,' maxfac = ',i8,/,
     2       ' rad = ',f10.6,' error = ',f8.4,/,
     3       ' srf = ',f10.6,' error = ',f8.4,/,
     4       ' vol = ',f10.6,' error = ',f8.4)
 9992 format(' srf, approx., ratio, error = ',2f12.8,2 f8.4)
 9991 format(' rad, approx., ratio, error = ',2f12.8,2 f8.4)
 9990 format(' absolute -ijkseg- , for -ipass- = ',i5)
 9989 format(' relative -ijkseg- , for -ipass- = ',i5)
 9988 format(12i6)
 9987 format(' assigned segment for -lface- = ',i7,
     1       ' is = ',i7,' ( greater than -maxseg- = ',i4,' )')
 9986 format(' nseg,nfac,ndiv=nfac/nseg,dsurf,dvol = ',3i7,2f10.6)
 9985 format('   pt  ','      x     ','      y     ','      z     ',
     1       ' seg ','    norm    ',/,1x,59(1h-))
 9984 format(i7,3f12.8,i5,f12.8)
 9983 format(' mseg,adiag,fac,m,sum = ',i7,4f12.6)
 9982 format(1x,'sphere from -octahedron-',/,1x,24(1h-))
 9981 format(' pass # = ',i5)
 9980 format(' minfac = ',10i5)
 9979 format(' maxfac = ',10i5)
 9978 format('      adiag,fac       = ',   2f12.6)
      end
C>
C> \brief Setup the data for an Octahedron
C>
C> This routine initiates the segments of an octahedron. The output
C> is split over a few arrays. One array `apex` holds the coordinates
C> of the corners of the triangles, another array `ijkfac` lists 
C> for each triangle which points in the `apex` array hold the
C> corresponding corner points, and finally `ijkseg` record the 
C> mapping between faces and segments.
C>
      subroutine hnd_sphoct(apex,napex,ijkfac,ijkseg,nface)
      implicit none
#include "global.fh"
#include "stdio.fh"
c
      logical out
      integer            napex     !< [Output] The number of apeces
      integer            nface     !< [Output] The number of faces
      double precision    xyz(3,6)
      integer             ijk(3,8)
      double precision   apex(3,*) !< [Output] The corners of the 
                                   !< triangles
      integer          ijkfac(3,*) !< [Output] For each triangle which
                                   !< points make up its corners
      integer          ijkseg(  *) !< [Output] For each face the number
                                   !< of the segment it belongs to
      integer iapex, lface
      data xyz / 1.0d+00, 0.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00,
     1          -1.0d+00, 0.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     2           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00/
      data ijk / 5, 1, 2, 5, 2, 3, 5, 3, 4, 5, 4, 1,
     1           6, 1, 2, 6, 2, 3, 6, 3, 4, 6, 4, 1/
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
      if(out) then
         write(luout,9997)
      endif
c
c     ----- set the 6 apeces of an octahedron -----
c
c     1     1.     0.     0.
c     2     0.     1.     0.
c     3    -1.     0.     0.
c     4     0.    -1.     0.
c     5     0.     0.     1.
c     6     0.     0.    -1.
c
      napex=6
      do iapex=1,napex
         apex(1,iapex)=xyz(1,iapex)
         apex(2,iapex)=xyz(2,iapex)
         apex(3,iapex)=xyz(3,iapex)
      enddo
      if(out) then
         write(luout,9999)
         do iapex=1,napex
            write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                           apex(3,iapex)
         enddo
      endif
c
      nface=8
      do lface=1,nface
         ijkfac(1,lface)=ijk(1,lface)
         ijkfac(2,lface)=ijk(2,lface)
         ijkfac(3,lface)=ijk(3,lface)
         ijkseg(  lface)=      lface
      enddo
c
      if(out) then
         write(luout,*) '...... end of -sphoct- ......'
      endif
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(/,1x,'octahedron',/,1x,10(1h-))
      end
c
      subroutine hnd_sphico(apex,napex,ijkfac,ijkseg,nface)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
      logical out
      dimension      c(3,12)
      dimension      s(3,12)
      dimension    ijk(3,20)
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      data c   / 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00,
     1           0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     2           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00,
     3           0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     4           1.0d+00, 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,
     5          -1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00/
      data s   / 0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     1           0.0d+00, 0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00,
     2           1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,
     3           1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00, 0.0d+00,
     4           0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00,
     5           0.0d+00, 1.0d+00, 0.0d+00, 0.0d+00,-1.0d+00, 0.0d+00/
      data ijk / 1, 2, 9, 1, 9, 5, 1, 5, 6, 1, 6,11, 1,11, 2,
     1                    2, 9, 7, 2, 7, 8, 2, 8,11,
     2           3, 4,10, 3,10, 5, 3, 5, 6, 3, 6,12, 3,12, 4,
     3                    4,10, 7, 4, 7, 8, 4, 8,12,
     4           9,10, 7, 9, 5,10,11,12, 8,11, 6,12/
      data one  /1.0d+00/
      data two  /2.0d+00/
      data five /5.0d+00/
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
      if(out) then
         write(luout,9997)
      endif
c
c     ----- set the 12 apeces of an icosahedron -----
c
c     1     0.     cosa   sina
c     2     0.     cosa  -sina
c     3     0.    -cosa   sina
c     4     0.    -cosa  -sina
c     5     sina   0.     cosa
c     6    -sina   0.     cosa
c     7     sina   0.    -cosa
c     8    -sina   0.    -cosa
c     9     cosa   sina   0.
c    10     cosa  -sina   0.
c    11    -cosa   sina   0.
c    12    -cosa  -sina   0.
c
      ang=acos(one/sqrt(five))/two
      cosa=cos(ang)
      sina=sin(ang)
      napex=12
      do iapex=1,napex
         apex(1,iapex)=cosa*c(1,iapex)+sina*s(1,iapex)
         apex(2,iapex)=cosa*c(2,iapex)+sina*s(2,iapex)
         apex(3,iapex)=cosa*c(3,iapex)+sina*s(3,iapex)
      enddo
      if(out) then
         write(luout,9999)
         do iapex=1,napex
            write(luout,9998) iapex,apex(1,iapex),apex(2,iapex),
     1                           apex(3,iapex)
         enddo
      endif
c
      nface=20
      do lface=1,nface
         ijkfac(1,lface)=ijk(1,lface)
         ijkfac(2,lface)=ijk(2,lface)
         ijkfac(3,lface)=ijk(3,lface)
         ijkseg(  lface)=      lface
      enddo
c
      if(out) then
         write(luout,*) '...... end of -sphico- ......'
      endif
      return
 9999 format(/,1x,'  apex',5x,'x',6x,5x,'y',6x,5x,'z',6x,/,1x,42(1h-))
 9998 format(1x,i6,f12.8,f12.8,f12.8)
 9997 format(/,1x,'icosahedron',/,1x,11(1h-))
      end
C>
C> \brief Partition a given triangle into four triangles and project
C> them outward onto the unit sphere
C>
      subroutine hnd_sphapx(apex,mapex,ijkfac,ijkseg,mface,lface,
     1                             ijkold,ijknew,dmin)
      implicit double precision (a-h,o-z)
#include "global.fh"
#include "stdio.fh"
c
      logical out
      logical duplic
      dimension   apex(3,*)
      dimension ijkfac(3,*)
      dimension ijkseg(  *)
      dimension ijkold(3)
      dimension ijknew(3)
      dimension    xyz(3,3)
      dimension      d(3)
      dimension xyzijk(3)
      equivalence (xyz(1,1),xij),(xyz(2,1),yij),(xyz(3,1),zij),
     1            (xyz(1,2),xjk),(xyz(2,2),yjk),(xyz(3,2),zjk),
     2            (xyz(1,3),xki),(xyz(2,3),yki),(xyz(3,3),zki)
      data zero  /0.0d+00/
      data pt5   /0.5d+00/
      data one   /1.0d+00/
      data three /3.0d+00/
      data tol   /1.0d-04/
c
      dist(x1,y1,z1,x2,y2,z2)=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
c
      out=.false.
      out=out.and.ga_nodeid().eq.0
c
c     ----- create mid-point of the 3 edges -----
c
      iold=ijkold(1)
      jold=ijkold(2)
      kold=ijkold(3)
      do m=1,3
         xyz(m,1)=(apex(m,iold)+apex(m,jold))*pt5
         xyz(m,2)=(apex(m,jold)+apex(m,kold))*pt5
         xyz(m,3)=(apex(m,kold)+apex(m,iold))*pt5
      enddo
c
c     ----- project onto sphere -----
c
      d(1)=dist(xij,yij,zij,zero,zero,zero)
      d(2)=dist(xjk,yjk,zjk,zero,zero,zero)
      d(3)=dist(xki,yki,zki,zero,zero,zero)
      d(1)=one/d(1)
      d(2)=one/d(2)
      d(3)=one/d(3)
      do l=1,3
         do m=1,3
            xyz(m,l)=xyz(m,l)*d(l)
         enddo
      enddo
c
c     ----- check for duplicate apeces -----
c
      newapx=0
      do iapx=1,3
         duplic=.false.
         lduplc=0
         do lapex=1,mapex
            dd=dist(xyz(1,  iapx),xyz(2,  iapx),xyz(3,  iapx),
     1              apex(1,lapex),apex(2,lapex),apex(3,lapex))
            if(abs(dd).lt.tol) then
               duplic=.true.
               lduplc=lapex
            endif
         enddo
         if(duplic) then
            ijknew(iapx)=lduplc
            if(out) then
               write(luout,9999) iapx,ijkold,lduplc
            endif
         else
            newapx=newapx+1
            japx=mapex+newapx
            ijknew(iapx)=japx
            do m=1,3
               apex(m,japx)=xyz(m,iapx)
            enddo
            if(out) then
               write(luout,9998) iapx,ijkold,japx,
     1                        apex(1,japx),apex(2,japx),apex(3,japx)
            endif
         endif
      enddo
      mapex=mapex+newapx
c
c     ----- make up new faces and their centers -----
c
      third=one/three
      dmin =one
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(1)
      ijkfac(2,mface)=ijknew(1)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(2)
      ijkfac(2,mface)=ijknew(1)
      ijkfac(3,mface)=ijknew(2)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijkold(3)
      ijkfac(2,mface)=ijknew(2)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      mface=mface+1
      ijkseg(  mface)=lface
      ijkfac(1,mface)=ijknew(1)
      ijkfac(2,mface)=ijknew(2)
      ijkfac(3,mface)=ijknew(3)
      do m=1,3
         xyzijk(m)=(apex(m,iold)+apex(m,jold)+apex(m,kold))*third
      enddo
      dijk=dist(xyzijk(1),xyzijk(2),xyzijk(3),zero,zero,zero)
      dmin=min(dmin,dijk)
c
      if(out) then
         write(luout,9997) dmin,mface
      endif
c
      return
 9999 format(' duplicated apex =',i2,' for face ',3i5,'. same as = ',i5)
 9998 format('    new     apex =',i2,' for face ',3i5,'.  newapx = ',i5,
     1       /,7x,3f12.8)
 9997 format(' --- dmin = ',f12.8,' --- mface = ',i10)
      end
c
C> \brief Calculate the screened COSMO charges
C>
C> ## Introduction ##
C>
C> In the COSMO model [1] the charges to represent the solvation effects
C> are obtained from solving a linear system of equations \f$ Ax=b \f$.
C> In this system \f$ b \f$ is the electrostatic potential at the 
C> point charge positions. The matrix \f$ A \f$ consists of the 
C> electrostatic interaction between two unit charges at the point
C> charge position, i.e.
C> \f{eqnarray*}{
C>   A_{\mu\nu} &=& ||t_\mu - t_\nu||^{-1} \\\\
C>   A_{\mu\mu} &=& 3.8 |S_\mu|^{-1/2} \\\\
C>              &=& a_{\mathrm{diag}} |S_\mu|^{-1/2}
C> \f}
C> where \f$ S_\mu \f$ is the surface area associated with the COSMO
C> charge (see [1] Eqs. 7a and 7b). The quantity \f$a_{\mathrm{diag}}\f$
C> is evaluated using Eq.(B1). This done in two steps in that the
C> this quantity is evaluated for the unit sphere in `hnd_cossph`, any
C> remaining scale factors are applied in the evaluation of \f$A\f$.
C>
C> The original COSMO approach has problems when charge approach each
C> other and the interaction becomes singular. Therefore we have 
C> replaced this with a model where the interactions are smooth [4]
C> which is outlined below.
C>
C> ## Surface charges satisfying Gauss's theorem ##
C>
C> Solving these equations gives the "raw" COSMO charges \f$ x \f$.
C> These charges should sum up to the total charge contained within 
C> the Solvent Accessible Surface. For 2 reasons that will be in
C> practice not be exactly true:
C>
C> - The discretization of the SAS is not perfect
C>
C> - The electron distribution being represented with Gaussian functions
C>   extends beyond the SAS
C>
C> Therefore the raw COSMO charges are corrected by adding corrections
C> based on a Lagrange multiplier technique [4]. The corresponding
C> equations can be derived by starting from an energy expression
C> in terms of the solute charge distribution and the surface charges.
C> This energy expression including the Lagrange term is
C> \f{eqnarray*}{
C>   E(Q,q,\lambda) &=& \frac{1}{2}Q^TCQ + Q^TBq
C>                   +  \frac{1}{2f(\epsilon)}q^TAq
C>                   +  \lambda\left(f(\epsilon)Q_{in}+\sum_i q_i\right)
C> \f}
C> where \f$Q\f$ is the charge distribution of the solute (including
C> both nucleii and electrons), \f$q\f$ is the surface charge
C> distribution, \f$f(\epsilon)\f$ is the dielectric screening constant
C> as discussed below, \f$Q_{in}\f$ is the charge within the cavity
C> or equivalently the solute charge, and finally \f$A\f$, \f$B\f$,
C> and \f$C\f$ are Coulomb matrices.
C>
C> From this equation the surface charges can be derived by minimizing
C> \f$E\f$ wrt \f$q\f$ and \f$\lambda\f$. This yields
C> \f{eqnarray*}{
C>    \frac{\partial E}{\partial q} 
C>       &=& BQ + \frac{1}{f(\epsilon)}Aq + \Lambda = 0\\\\
C>    \frac{\partial E}{\partial \lambda}
C>       &=& f(\epsilon)Q_{in} + \sum_i q_i = 0
C> \f}
C> where \f$\Lambda\f$ is a vector of which each element is
C> \f$\lambda\f$, i.e. \f$\forall_i, \Lambda_i = \lambda\f$. 
C> Next we get
C> \f{eqnarray*}{
C>   q &=& -f(\epsilon)A^{-1}\left(BQ+\Lambda\right) \\\\
C>   \lambda &=& \frac{Q_{in}-\sum_i\left[A^{-1}BQ\right]_i}{
C>                     \sum_{ij}A^{-1}_{ij}}
C> \f}
C> Because \f$E\f$ is variationally optimized wrt \f$Q\f$, \f$q\f$,
C> and \f$\lambda\f$ the gradient expression only involves derivatives
C> of \f$A\f$, \f$B\f$ and \f$C\f$ just like the original COSMO
C> gradients [1].
C>
C> Previously the surface charge correction was implemented by scaling
C> the raw COSMO charges.
C> However, this led to complications with neutral molecules
C> where the correct integrated surface charge is 0. Hence the 
C> correction factor would be 0 as well, eliminating the COSMO charges
C> and hence all solvation effects. This problem had been patched by
C> calculating the COSMO charges for the nucleii and the electrons 
C> separately. This led to a cumbersome and expensive algorithm. So
C> the use of a Lagrange constraint is a solution that is
C> preferred over scaling the charges. The difference between the two
C> approaches should be small provided the discretization is fine
C> enough.
C>
C> In the COSMO model [1] it is realized that dielectric screening
C> scales as
C> \f{eqnarray*}{
C>   f(\epsilon) &= \frac{\epsilon-1}{\epsilon+a}, & 0\le a \le 2
C> \f}
C> Klamt and Sch&uuml;&uuml;rmann suggested to use \f$ a = 1/2 \f$, 
C> essentially based on an argument that the true value should not 
C> exceed 1 (see appendix A). Stevanovich and Truong [3] realized that
C> the screened charges should be such that the Gauss theorem holds,
C> which requires \f$ a = 0 \f$ (see Eq. (5)). Based on this physical
C> motivation the latter value is used by default.
C>
C> The linear system of equations to be solved may be tackled in 2
C> different ways. For small systems a direct solver is appriopriate,
C> whereas for large systems an iterative solver is used. These solvers
C> have different implications. If \f$ N \f$ is the number of COSMO
C> charges then
C>
C> - the direct solver uses \f$ O(N^2) \f$ memory and \f$ O(N^3) \f$
C>   instructions
C>
C> - the iterative solver uses \f$ O(N) \f$ memory and \f$ k*O(N^2) \f$
C>   instructions (\f$ k \f$ is the number of iterations to convergence)
C>
C> The costs are based on the assumptions that for the direct solver
C> the matrix is stored explicitly and the inverse is not stored but
C> recalculated every time. For the iterative solver only a 
C> matrix-vector multiplication is implemented that regenerates the
C> matrix elements every time, and the number of iterations is roughly
C> independent from the values of the matrix. 
C>
C> Based on these assumptions the iterative solver is optimal for large
C> systems both with respect to memory requirements as well as compute
C> requirements. The case for the iterative solver can be improved 
C> further by parallelizing the matrix-vector multiplication which
C> reduces the compute cost per processor to \f$ k*O(N^2)/N_{proc} \f$.
C> Furthermore if we start the iterative solver in each SCF cycle with
C> the solution from the previous iteration rather than \f$ x=0 \f$ then
C> \f$ k \f$ may be reduced as well. In practice \f$ k \f$ does not
C> depend strongly on the initial value of \f$ x \f$, reductions by at
C> most a factor 2 are seen when the SCF is nearly converged.
C>
C> ## Singularity free surface charge self-interaction ##
C>
C> The cavity in continuum solvation models is constructed by creating
C> spherical cavities around all atoms and merging these volumes. The
C> Solvent Accessible Surface (SAS) is created representing the
C> spherical surface around every atom with points and eliminating the 
C> points that fall inside the sphere around a neighboring atom.
C> Klamt et al. [1] suggested representing the surface charge by point
C> charges at the surface discretization points. This leads to
C> singularities in the solvation energy when some of these points come
C> together, typically in the vicinity of the boundary between spheres.
C>
C> The singularities in the surface self-interaction energy need to be
C> addressed to ensure that sensible geometry optimizations are 
C> possible. York and Karplus [4] suggested formally replacing the 
C> surface point charges by Gaussian charge distributions. This leads
C> to an interaction of the form
C> \f{eqnarray*}{
C>   A(r_i,r_j) &=& \frac{\mathrm{erf}(\zeta_{ij}r_{ij})}{r_{ij}}
C> \f}
C> This interaction is relatively easily implemented in force fields
C> where all charges are point charges. Implementing this for the
C> interaction between a point charge and an electron distribution
C> is conceptually harder to do.
C>
C> A consistent implementation would need to use this expression for
C> both the \f$A\f$ and \f$B\f$ matrix. Scalmani and Frisch [5] have
C> sought to do this exactly by explicitly representing the surface
C> charges as Gaussians. This could be done of course using the 
C> charge density fitting integrals but it would require some 
C> engineering as the surface charges would need to be stored in a
C> geometry object. Lange and Herbert [6,7] have followed York and
C> Karplus more closely and applied this approximation only to
C> matrix \f$A\f$ and not to the surface charge-electron interaction.
C> However, even applying this interaction just for \f$A\f$ is 
C> problematic as the expression still contains singularities at
C> \f$r_{ij} = 0\f$ that have to be handled explicitly to avoid floating
C> point exceptions.
C>
C> Considering the use of the potential above in more detail we have
C> \f{eqnarray*}{
C>   A(r_i,r_j) &=& \frac{\mathrm{erf}(\zeta_{ij}r_{ij})}{r_{ij}} \\\\
C>   A(r_i,r_i) &=& \lim_{r_{ij}\to 0}\frac{\mathrm{erf}(\zeta_{ii}r_{ij})}{r_{ij}}(F(r_i))^{-1} \\\\
C>              &=& \zeta_i\sqrt{2/\pi}(F(r_i))^{-1} \\\\
C>   \zeta_{ij} &=& \frac{\zeta_i\zeta_j}{
C>                        \left(\zeta_i^2+\zeta_j^2\right)^{1/2}} \\\\
C>   \zeta_i    &=& \frac{\zeta}{R_I\sqrt{w_i}} \\\\
C>              &=& \frac{\zeta}{\sqrt{|S_i|}} \\\\
C>              &=& \frac{\zeta\sqrt{M}}{R_I\sqrt{2\pi}}
C> \f}
C> where \f$R_I\f$ is the radius of the cavity around atom \f$I\f$,
C> \f$w_i\f$ is the weight of point \f$i\f$ on the unit sphere, hence 
C> \f$R_I\sqrt{w_i}\f$ is equivalent to the surface \f$|S_i|\f$ of the
C> point, \f$\zeta\f$ is a width parameter for the Gaussian distribution
C> that has been optimized to reproduce the Born solvation energy,
C> \f$M\f$ is the number of discretization points on the sphere.
C> In [4] Table 1 it is shown that \f$\zeta\f$ is essentially
C> \f$4.90\f$ for Lebedev grids (for the Boundary Element Mesh we
C> use it is \f$1.00\f$). Finally \f$F_i\f$ is the switching function
C> defined below.
C> From this the limit of two point charges approaching eachother can
C> be established as
C> \f{eqnarray*}{
C>   \lim_{r_j \to r_i} A(r_i,r_j) &=& \zeta_{ij}\frac{2}{\sqrt{\pi}}
C> \f}
C> When two point charges come so close together that \f$r_{ij} < C\f$
C> then this last expression has to be used, otherwise the regular
C> expression for \f$A(r_i,r_j)\f$ should be used.
C>
C> The gradient of this expression is given by
C> \f{eqnarray*}{
C>   \nabla_M A(r_i,r_j)
C>   &=& -\left(\mathrm{erf}(\zeta_{ij}r_{ij})-
C>        \frac{2\zeta_{ij}}{\sqrt{\pi}}
C>        e^{-\zeta_{ij}^2r_{ij}^2}\right)
C>        \frac{\nabla_M r_{ij}}{r_{ij}^2} \\\\
C>   \nabla_M A(r_i,r_i)
C>   &=& -A(r_i,r_i)\sum_B\frac{\partial F}{\partial R_B}\nabla_M r_i
C> \f}
C>
C> ## Continuous switching functions ##
C>
C> In order to obtain a smooth function wrt the nuclear coordinates
C> it is necessary that when the atoms move the surface areas associated
C> with a point charge change smoothly. The approach suggested by 
C> York and Karplus [4] proposes to use multiple radii around an atom. 
C> Each atom has an inner radius \f$R_{in}\f$ and an outer radius
C> \f$R_{out}\f$. The difference between them is the switching radius 
C> \f$R_{sw} = R_{out} - R_{in}\f$. The areas of cavity surface points
C> are multiplied with a weighting factor based on their relative
C> position. In practice the diagonal elements of \f$A\f$ are scaled 
C> by the inverse of the surface areas generating corresponding weights
C> as
C> \f{eqnarray*}{
C>   r_{A,i} - R_B \ge R_{out},&& W_{Ai,B} = 1 \\\\
C>   r_{A,i} - R_B \le R_{in}, && W_{Ai,B} = \infty \\\\
C>   R_{in} < r_{A,i}-R_B < R_{out},
C>      && W_{Ai,B} = 1/f\left(\frac{r_{A,i}-R_B-R_{in}}{R_{SW}}\right)
C> \f}
C> where the function \f$f(r)\f$ is given by
C> \f{eqnarray*}{
C>   f(r) &=& r^3\left(10-15r+6r^2\right) \\\\
C>   \frac{\partial f(r)}{\partial r} &=& 30r^2(r-1)^2
C> \f}
C> The weighting function overall is given by
C> \f{eqnarray*}{
C>    F_{Ai} = \prod_{B\neq A} W_{Ai,B}
C> \f}
C> It can be shown that if \f$F_{Ai}\f$ is close to \f$0\f$ then
C> the corresponding point will not contribute to the energy
C> expression and can be eliminated. This is most easily shown by
C> considering the energy expression and substituting the surface
C> charge expression. If \f$F_{Ai}\f$ goes to 0 then the diagonal
C> element \f$A(r_i,r_i)\f$ approaches infinity. This means that the
C> corresponding row and column in \f$A^{-1}\f$ goes to zero, and
C> the energy expression depends only on \f$A^{-1}\f$ hence the 
C> corresponding point will not contribute to the energy.
C>
C> One interesting observation is that matrix \f$A\f$ is used in
C> linear system of equations \f$Aq = BQ\f$. The condition number of
C> the matrix \f$A\f$ determines the ratio of the relative error in
C> \f$x\f$ and the relative error in \f$BQ\f$. In particular, the 
C> larger the condition number the larger the relative error in \f$x\f$
C> for a given relative error in \f$BQ\f$. The condition number
C> is given by
C> \f{eqnarray*}{
C>    K(A) = \left|\frac{\lambda_{\mathrm{max}}(A)}{
C>                       \lambda_{\mathrm{min}}(A)}\right|
C> \f}
C> where \f$\lambda_{\mathrm{max}}(A)\f$ is the maximum eigenvalue, 
C> and \f$\lambda_{\mathrm{min}}(A)\f$ is the minimal eigenvalue of 
C> \f$A\f$. The approach by York and Karplus for eliminating surface
C> charges is based on raising the condition number to infinity, at
C> which point the vector \f$q\f$ must become very inaccurate. Hence
C> the cutoff for eliminating charges must be chosen carefully to
C> limit the condition number of the remaining part of \f$A\f$.
C>
C> ### References ###
C>
C>   [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>   "COSMO: a new approach to dielectric screening in solvents with
C>    explicit expressions for the screening energy and its gradient",
C>   <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>   <a href="http://dx.doi.org/10.1039/P29930000799">
C>   10.1039/P29930000799</a>.
C>
C>   [2] M.A. Aguilar, F.J. Olivares del Valle, J. Tomasi,
C>   "Nonequilibrium solvation: An ab initio quantummechanical method
C>    in the continuum cavity model approximation",
C>   <i>J. Chem. Phys.</i> (1993) <b>98</b>, pp 7375-7384, DOI:
C>   <a href="http://dx.doi.org/10.1063/1.464728">
C>   10.1063/1.464728</a>.
C>
C>   [3] E.V. Stefanovich, T.N. Truong,
C>   "Optimized atomic radii for quantum dielectric continuum solvation
C>    models", <i>Chem. Phys. Lett.</i> (1995) <b>244</b>, pp 65-74,
C>   DOI:
C>   <a href="http://dx.doi.org/10.1016/0009-2614(95)00898-E">
C>   10.1016/0009-2614(95)00898-E</a>.
C>
C>   [4] D.M. York, M. Karplus,
C>   "A smooth solvation potential based on the conductor-like 
C>    screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>   pp 11060-11079, DOI:
C>   <a href="http://dx.doi.org/10.1021/jp992097l">
C>   10.1021/jp992097l</a>.
C>
C>   [5] G. Scalmani, M.J. Frisch,
C>   "Continuous surface charge polarizable continuum models of
C>    solvation. I. General formalism", <i>J. Chem. Phys.</i> (2010)
C>   <b>132</b>, 114110, DOI:
C>   <a href="http://dx.doi.org/10.1063/1.3359469">
C>   10.1063/1.3359469</a>.
C>
C>   [6] A.W. Lange, J.M. Herbert, 
C>   "Polarizable continuum reaction-field solvation models affording
C>    smooth potential energy surfaces", <i>J. Phys. Chem. Lett.</i>
C>   (2010) <b>1</b>, pp 556-561, DOI:
C>   <a href="http://dx.doi.org/10.1021/jz900282c">
C>   10.1021/jz900282c</a>.
C>
C>   [7] A.W. Lange, J.M. Herbert,
C>   "A smooth, nonsingular, and faithful discretization scheme for
C>    polarizable continuum models: The switching/Gaussian approach",
C>   <i>J. Chem. Phys.</i> (2010) <b>133</b>, 244111, DOI:
C>   <a href="http://dx.doi.org/10.1063/1.3511297">
C>   10.1063/1.3511297</a>.
C>
      subroutine hnd_coschg(g_dens,ndens,rtdb,geom,basis,nat,nefc,
     &                      efcc,efcs,efcz,efciat,ratm,ecos,
     &                      cosmo_file)
      implicit none
#include "errquit.fh"
#include "global.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "nwc_const.fh"
#include "geom.fh"
#include "bq.fh"
#include "stdio.fh"
#include "prop.fh"
#include "util.fh"
#include "inp.fh"
c
      integer ndens         !< [Input] the number of density matrices
      integer g_dens(ndens) !< [Input] the handle for the density
                            !< matrices
      integer rtdb          !< [Input] the RTDB handle
      integer geom          !< [Input] the molecular geometry handle
      integer basis         !< [Input] the molecular basis set handle
      integer nat           !< [Input] the number of atoms
      integer nefc          !< [Input] the number of COSMO charges
c
      double precision efcc(3,nefc) !< [Input] the COSMO charge
                                    !< coordinates
      double precision efcs(nefc)   !< [Input] the COSMO charge
                                    !< surface area
c
      double precision efcz(nefc)   !< [Output] the COSMO charges
      integer          efciat(nefc) !< [Input] the atom associated
                                    !< with each surface charge
      double precision ratm(nat)    !< [Input] the atom radii
      double precision ecos !< [Output] the energy contribution due to
                            !< the COSMO charges
      logical  status
c
      logical  dbug,more,out,direct,noall,all,elec,nucl,iefc_done
      character*16 at_tag
      integer istrlen
      character*255 cosmo_file
      integer fn
      integer lineq ! 0: fast direct solver, 1: slow iterative solver
      integer minbem 
      integer maxbem 
      integer ificos ! 0 use octahedron, 1 use icosahedron tesselation
      integer ifscrn ! do screening/correction?
      integer iat ! counter over atoms
      integer jef ! counter over COSMO charges
      integer l_i10, i10
      integer l_i11, i11
      integer l_i12, i12
      integer l_i20, i20
      integer l_i21, i21
      integer l_i22, i22
      integer l_i30, i30
      integer l_i40, i40
      integer l_i50, i50
      integer l_i60, i60
      integer l_i70, i70
      integer l_i80, i80
      integer l_i90, i90
      integer i,ipt,ief ! counters
      integer i_init ! number of ints in memory requirement vector
      integer init ! memory requirement vector
      integer ierr ! error flag
      integer iep ! memory offset of b from Ax=b
      integer nodcmp ! flag specifying how to handle errors
      integer need ! the amount of memory needed
c
      integer l_epot, l_xyzpt, l_zanpt ! memory handles
      integer k_epot, k_xyzpt, k_zanpt ! memory offsets
c
      double precision adiag
      double precision dsurf
      double precision dvol
      double precision srfmol
      double precision volmol
      double precision dielec    ! dielectric constant (slow component)
      double precision dielecinf ! optical dielectric constant (fast component)
      double precision screen ! dielectric screening factor
      double precision rsolv
      double precision zeta,gammas,swtol
      double precision ptspatm
      common/hnd_cosmod/lineq,minbem,maxbem,ificos,ifscrn
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision charge ! the total QM region charge
      double precision chgfac ! scale factor for COSMO charges
      double precision chgcos ! the total COSMO surface charge
      double precision chgcvg ! the convergence of the COSMO charges
      double precision chgina ! the inv(A) COSMO charge
      double precision corcos ! the COSMO charge correction
      double precision errcos ! the COSMO charge error
      double precision delchg ! charge difference
      double precision aii,aij,bij,chgief,dij,deta,dum,oldief
      double precision atmefc ! atom - COSMO charge interaction
      double precision efcefc ! COSMO charge - COSMO charge interaction
      double precision elcefc ! electron - COSMO charge interaction
      double precision allefc ! total QM - COSMO charge interaction
      double precision zan ! charges
      double precision xi, xj, xn, xp ! X-coordinates
      double precision yi, yj, yn, yp ! Y-coordinates
      double precision zi, zj, zn, zp ! Z-coordinates
      double precision qi, qj ! charges
      double precision rr ! distance
      double precision solnrg ! solvation energy
      double precision dlambda ! lambda (surface charge correction)
      double precision pi
      double precision x(1)
      equivalence (x(1),dbl_mb(1))
      double precision zero, pt5, one, two ! constants
      data zero   /0.0d+00/
      data pt5    /0.5d+00/
      data one    /1.0d+00/
      data two    /2.0d+00/
      double precision zetai, zetaj, zetaij
      double precision bohr
      parameter (bohr=0.529177249d0)
      logical stat
c
      double precision util_erf
      logical util_io_unit
      external util_io_unit, util_erf
c
      more=.false.
      dbug=.false.
      dbug=dbug.or.more
      out =.false.
      out =out.or.dbug
c
      pi = acos(-1.0d0)
c
      dbug=dbug.and.ga_nodeid().eq.0
      more=more.and.ga_nodeid().eq.0
      out =out .and.ga_nodeid().eq.0
c
      iefc_done=.false.
      all=.false.
      ecos=zero
      istrlen = 0
c
c     ----- get electrostatic potential at surface points -----
c
c     --- total field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i10",l_i10,i10))
     &   call errquit("hnd_coschg: malloc i10 failed",913,MA_ERR)
c     --- nuclear field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i11",l_i11,i11))
     &   call errquit("hnd_coschg: malloc i11 failed",914,MA_ERR)
c     --- electron field
      if (.not.ma_push_get(mt_dbl,nefc,"hnd_coschg:i12",l_i12,i12))
     &   call errquit("hnd_coschg: malloc i12 failed",914,MA_ERR)
c
c     ----- calculate electronic contribution at all points -----
c
      call hnd_elfcon(basis,geom,g_dens(ndens),efcc,nefc,x(i12),0)
      do ipt=1,nefc
         x(i12+ipt-1) = -x(i12+ipt-1)
      enddo
c
c     ----- nuclear contribution -----
c
      if (.not.geom_ncent(geom,nat)) call
     &    errquit('hnd_coschg: geom_ncent',911,GEOM_ERR)
      if (.not. ma_push_get(mt_dbl,3*nat,'xyz pnt',l_xyzpt,k_xyzpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'epot pnt',l_epot,k_epot))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      if (.not. ma_push_get(mt_dbl,nat,'zan pnt',l_zanpt,k_zanpt))
     &    call errquit('hnd_coschg: ma failed',911,MA_ERR)
      do iat=1,nat
        if(.not.geom_cent_get(geom,iat,at_tag,dbl_mb(k_xyzpt+3*(iat-1)),
     &     dbl_mb(k_zanpt+iat-1))) call
     &     errquit('hnd_coschg: geom_cent_get',911,GEOM_ERR)
      enddo ! iat
c
      do ipt=1,nefc
         xp = efcc(1,ipt)
         yp = efcc(2,ipt)
         zp = efcc(3,ipt)
         x(i11+ipt-1) = 0.0d0
         do i = 1,nat
            xn  = dbl_mb(k_xyzpt  +3*(i-1)) - xp
            yn  = dbl_mb(k_xyzpt+1+3*(i-1)) - yp
            zn  = dbl_mb(k_xyzpt+2+3*(i-1)) - zp
            zan = dbl_mb(k_zanpt+i-1)
            rr =  sqrt(xn*xn + yn*yn + zn*zn)
            x(i11+ipt-1) = x(i11+ipt-1) + zan/rr
         enddo ! i
         x(i10+ipt-1) = x(i11+ipt-1) + x(i12+ipt-1)
      enddo ! ipt
c
c     ----- get surface charges -----
c
c     ----- set up the memory based on the solver -----
c     lineq = 0: fast direct solver, lineq = 1: slow iterative solver
c
      if (lineq.eq.0) then
c
        stat = .true.
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i20",
     &                                l_i20,i20)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i21",
     &                                l_i21,i21)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i22",
     &                                l_i22,i22)
        stat = stat .and. ma_push_get(mt_dbl,nefc*nefc,"hnd_coschg i30",
     &                                l_i30,i30)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i40",
     &                                l_i40,i40)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i50",
     &                                l_i50,i50)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i80",
     &                                l_i80,i80)
c
c      check memory
c
       if (.not.stat) then
        call errquit("hnd_coschg: out of memory: lineq = 0 ",950,MA_ERR)
       endif
c
      else if (lineq.eq.1) then
c
        stat = .true.
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i20",
     &                                l_i20,i20)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i21",
     &                                l_i21,i21)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i22",
     &                                l_i22,i22)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i30",
     &                                l_i30,i30)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i40",
     &                                l_i40,i40)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i50",
     &                                l_i50,i50)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i60",
     &                                l_i60,i60)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i70",
     &                                l_i70,i70)
        stat = stat .and. ma_push_get(mt_dbl,nefc,"hnd_coschg i80",
     &                                l_i80,i80)
c
c      check memory
c
       if (.not.stat) then
        call errquit("hnd_coschg: out of memory: lineq = 1 ",950,MA_ERR)
       endif
c
      else 
       call errquit("hnd_coschg: unknown solver",911,INPUT_ERR)
      end if ! lineq 
c
      if(out) then
         if(lineq.eq.1) then
            write(luout,*) 'use iterative method for -lineq-'
         else
            write(luout,*) 'use in-memory method for -lineq-'
         endif
      endif
c
      if(lineq.eq.0) then
c
c     ----- in memory -----
c
         if(out) then
            write(luout,*) '-lineq- in memory'
         endif
c
c     ----- calculate q* = A^{-1}BQ -----
c
         iep=i10

         do ief=1,nefc
            x(ief+i20-1)=x(ief+iep-1)
         enddo

         call hnd_cosmata(nat,nefc,efcc,efcs,efciat,ratm,x(i30))
c
         nodcmp=1
         call hnd_linequ(x(i30),nefc,x(i20),nefc,
     1                   x(i40),x(i50),deta,ierr,nodcmp)
c
c     ----- calculate t = A^{-1}1 -----
c
         do ief=1,nefc
            x(ief+i80-1)=one
         enddo

         call hnd_cosmata(nat,nefc,efcc,efcs,efciat,ratm,x(i30))
c
         nodcmp=1
         call hnd_linequ(x(i30),nefc,x(i80),nefc,
     1                   x(i40),x(i50),deta,ierr,nodcmp)
c
      else
c
c     ----- iterative -----
c
         if(out) then
            write(luout,*) 'iterative -lineq-'
         endif
c
c     ----- calculate qraw = q* = A^{-1}BQ -----
c
         if (.not.rtdb_get(rtdb,'cosmo:qraw',mt_dbl,nefc,x(i20))) then
c
c           If no raw COSMO charges were found initialize x=0
c
            do ief=1,nefc
               x(ief+i20-1)=zero
            enddo
         endif
c
         direct=.true.
c
c     ----- solve ... -----
c
         iep=i10

         call hnd_cosequ(nat,x(iep),x(i20),nefc,
     1                   x(i40),x(i50),x(i60),x(i70),
     2                   efcc,efcs,efciat,ratm)
c
         if (.not.rtdb_put(rtdb,'cosmo:qraw',mt_dbl,nefc,x(i20))) then
            call errquit('hnd_coschg: could not store raw COSMO charge',
     1                   nefc,RTDB_ERR)
         endif
c
c     ----- calculate t = A^{-1}1 -----
c
         if (.not.rtdb_get(rtdb,'cosmo:rawt',mt_dbl,nefc,x(i80))) then
c
c           If no raw A^{-1}1 charges were found initialize x=0
c
            do ief=1,nefc
               x(ief+i80-1)=zero
            enddo
         endif
c
         direct=.true.
c
c     ----- solve ... -----
c
         do ief=1,nefc
            x(ief+i21-1)=one
         enddo
         iep=i21

         call hnd_cosequ(nat,x(iep),x(i80),nefc,
     1                   x(i40),x(i50),x(i60),x(i70),
     2                   efcc,efcs,efciat,ratm)
c
         if (.not.rtdb_put(rtdb,'cosmo:rawt',mt_dbl,nefc,x(i80))) then
            call errquit('hnd_coschg: could not store raw A^{-1}1',
     1                   nefc,RTDB_ERR)
         endif
      endif
c
c     ----- correct the COSMO charges ... -----
c
      chgcos=zero
      chgina=zero ! inverse of A
      do ief=1,nefc
         chgief=x(ief+i20-1)
         chgcos=chgcos+chgief
         chgina=chgina+x(i80+ief-1)
      enddo
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))
     $     charge = 0.0d0
      errcos=charge-chgcos
      dlambda=errcos/chgina
      if (ifscrn.ne.1) then
         do ief=1,nefc
            x(ief+i20-1)=x(ief+i20-1)+dlambda*x(i80+ief-1)
         enddo
      endif
c
c     ----- charge screening due to the dielectric medium -----
c
c     ----- set screening factor -----
c
      if(ifscrn.eq.0) then
         chgfac=one
      elseif(ifscrn.eq.1) then
         chgfac=screen
      elseif(ifscrn.eq.2) then
         chgfac=screen
      else
         write(luout,*) 'illegal -ifscrn- ... stop. = ',ifscrn
         call errquit('hnd_coschg: illegal -ifscrn-...',911,UERR)
         chgfac=one
      endif

      if(out) then
         write(luout,*) ' ifscrn                       = ',ifscrn
         write(luout,9984) chgfac
      endif
c
c     ----- apply screening factor -----
c
      chgcos=zero
      do ief=1,nefc
         x(ief+i20-1)=-chgfac*x(ief+i20-1)
         chgcos=chgcos+x(ief+i20-1)
      enddo
c
c     ----- store effective charges in -efcz- ... -----
c           check convergence ...
c
      chgcvg=zero
      do ief=1,nefc
         oldief=efcz(ief)
         chgief=   x(ief+i20-1)
         delchg=abs(chgief-oldief)
         if(delchg.gt.chgcvg) then
            chgcvg=delchg
         endif
         efcz(ief)=chgief
      enddo

      if(all) then
         if(out) then
            write(luout,9987) chgcvg
         endif
      endif

      if(dbug) then
         write(luout,9998)
         do ief=1,nefc
            write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
         enddo
      endif
c
c     ----- calculate energy terms from -cosmo- charges -----
c
      allefc=zero
      atmefc=zero
      efcefc=zero
      do jef=1,nefc
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         qj=efcz(  jef)
         zetaj=zeta*sqrt(ptspatm)/(ratm(efciat(jef))*sqrt(2.0d0*pi))
c
         allefc=allefc+qj*x(jef+i10-1)
c
         do iat=1,nat
            xi=dbl_mb(k_xyzpt  +3*(iat-1))
            yi=dbl_mb(k_xyzpt+1+3*(iat-1))
            zi=dbl_mb(k_xyzpt+2+3*(iat-1))
            qi=dbl_mb(k_zanpt+iat-1)
            dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
            bij=one/dij
            atmefc=atmefc+qi*bij*qj
         enddo
c
         do ief=1,nefc
            zetai=zeta*sqrt(ptspatm)/(ratm(efciat(ief))*sqrt(2.0d0*pi))
            xi=efcc(1,ief)
            yi=efcc(2,ief)
            zi=efcc(3,ief)
            qi=efcz(  ief)
            if(ief.eq.jef) then
               aij=efcs(ief)
            else
               dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
               zetaij=zetai*zetaj/sqrt(zetai**2+zetaj**2)
               if (dij.lt.1.0d-6) then
                 aij=2.0d0*zetaij/sqrt(pi)
               else
                 aij=util_erf(zetaij*dij)/dij
               endif
            endif
            efcefc=efcefc+qi*aij*qj
         enddo
      enddo
      efcefc= efcefc/(two*chgfac)
      elcefc= allefc-atmefc
      solnrg= allefc+efcefc
      ecos  = atmefc+efcefc
      if (out) then
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9995) efcefc
         write(luout,9988) solnrg
         write(luout,9989) allefc,(-two*efcefc)
         write(luout,9994) ecos  
      endif
c
c     ----- other form of the solvation energy ... -----
c
      allefc=zero
      atmefc=zero
      elcefc=zero
      do jef=1,nefc
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         qj=efcz(  jef)
c
         allefc=allefc+qj*x(jef+i10-1)
         atmefc=atmefc+qj*x(jef+i11-1)
         elcefc=elcefc+qj*x(jef+i12-1)
      enddo
      solnrg= pt5* allefc
c     no idea what the following line is supposed to mean
c     ecos  = pt5*(atmefc-elcefc)
      if (out) then
         write(luout,9991) atmefc
         write(luout,9990) elcefc
         write(luout,9989) allefc
         write(luout,9988) solnrg
         write(luout,9994) ecos
      endif

      if(dbug) then
         write(luout,9998)
         do ief=1,nefc
            write(luout,9997) ief,(efcc(i,ief),i=1,3),efcz(ief)
         enddo
      endif
c
      if(out) then
         write(luout,9993)
      endif
c
c     ----- save -cosmo- charges and energy to -rtdb- -----
c
      if (.not. rtdb_put(rtdb,'cosmo:energy',mt_dbl,1,ecos))
     &   call errquit('hnd_coschg: rtdb put failed for ecos',911,
     &       rtdb_err)
      if(.not.rtdb_put(rtdb,'cosmo:efcz',mt_dbl,  nefc,efcz))
     $   call errquit('hnd_coschg: rtdb put failed for efcz',912,
     &       rtdb_err)
c
c     ----- for the time being, save in 'geometry' object -----
c
      if(out) then
         write(luout,*) 'in -hnd_coschg ... for -efc- geom = ',geom
      endif
      status=bq_set(cosmo_bq_efc,nefc,efcz,efcc)
      if (.not.status) then
         call errquit('hnd_coschg: bq_set failed !', 0,
     &       geom_err)
      endif
      status=bq_rtdb_store(rtdb,cosmo_bq_efc)
      if (.not.status) then
         call errquit('hnd_coschg: bq_rtdb_store failed !', 0,
     &       geom_err)
      endif
     
c
c     ----- printing cosmo charges for bq module -----
c
      istrlen = inp_strlen(cosmo_file)
      if (istrlen.le.0) cosmo_file = "cosmo.xyz"
      call util_file_name_resolve(cosmo_file,.false.)
      if(ga_nodeid().eq.0) then
        if(.not.util_io_unit(80,90,fn))
     +     call errquit('cannot get free unit', 0,
     +       0)
c
        open(unit=fn,form="formatted",file=cosmo_file)
        if (dbug) then
          write(*,*) "printing cosmo charges for bq module",
     +     cosmo_file
        end if
c
        write(fn,*) nefc
        write(fn,*) "cosmo charges (= -Bq-charge) (coordinates in Angstr
     +om) format: Bq x y z q"
        do ief=1,nefc
           write(fn,*) 
     +      "Bq",
     +      efcc(1,ief)*bohr,
     +      efcc(2,ief)*bohr,
     +      efcc(3,ief)*bohr,
     +     -efcz(  ief)
        end do
        close(fn)
      end if
c
c     ----- release memory block -----
c
      if(.not.ma_chop_stack(l_i10)) call
     &  errquit('hnd_coschg, ma_pop_stack of init failed',911,MA_ERR)
c
      return
 9999 format(/,10x,15(1h-),
     1       /,10x,'-cosmo- charges',
     2       /,10x,15(1h-))
 9998 format(4x,'iefc',6x,'x',5x,6x,'y',5x,6x,'z',5x,5x,'q',4x,
     1     /,1x,53(1h-))
 9997 format(1x,i7,3f12.6,f10.6)
 9995 format(' -efcefc- energy = ',f20.12)
 9994 format(' -ecos  - energy = ',f20.12)
 9993 format(' ...... end of -coschg- ......')
 9992 format(' cosmo potential at -ief = ',i6,f16.10)
 9991 format(' -atmefc- energy = ',f20.12)
 9990 format(' -elcefc- energy = ',f20.12)
 9989 format(' -allefc- energy = ',f20.12,f20.12)
 9988 format(' -solnrg- energy = ',f20.12)
 9987 format(' -cosmo- charges convergence = ',f10.6)
 9986 format(' -wfntyp- = ',a8)
 9985 format(' -scftyp- = ',a8)
 9984 format(' applied screening factor -chgfac- = ',f10.6)
      end
C>
C> \brief Compute matrix A
C> 
C> Compute matrix `A` as needed for the in-memory COSMO
C> charge fitting.
C>
      subroutine hnd_cosmata(nat,nefc,efcc,efcs,efciat,ratm,a)
      implicit none
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      integer          efciat(nat)  !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision a(nefc,nefc) !< [Output] Matrix `A`
c
      integer jef, ief
      double precision aii, xi, xj, yi, yj, zi, zj, dij, one
      double precision zetai, zetaj, zetaij
      parameter (one = 1.0d0)
c
      double precision dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision adiag,dsurf,dvol,srfmol,volmol,ptspatm
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
c
      double precision util_erf
      external util_erf
c
      double precision pi
      pi = acos(-1.0d0)
c
      do jef=1,nefc
         zetaj=zeta*sqrt(ptspatm)/(ratm(efciat(jef))*sqrt(2.0d0*pi))
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         do ief=1,nefc
            if(ief.eq.jef) then
               aii=efcs(ief)
               a(ief,jef)=aii
            else
               zetai=zeta*sqrt(ptspatm) /
     &               (ratm(efciat(ief))*sqrt(2.0d0*pi))
               zetaij=zetai*zetaj/sqrt(zetai**2+zetaj**2)
               xi=efcc(1,ief)
               yi=efcc(2,ief)
               zi=efcc(3,ief)
               dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
               if (dij.lt.1.0d-6) then
                 a(ief,jef)=2.0d0*zetaij*sqrt(1.0d0/pi)
               else
                 a(ief,jef)=util_erf(zetaij*dij)/dij
               endif
            endif
         enddo
      enddo
c
      end
C>
C> \brief On the fly matrix-vector multiplication
C>
C> This routine multiplies the COSMO matrix with the current guess
C> for the charge vector using a dot-product based algorithm. The matrix
C> is generated on the fly. For performance reasons the routine is
C> replicated data parallel.
C>
      subroutine hnd_cosaxd(nat,x,ax,nefc,efcc,efcs,efciat,ratm)
      implicit none
#include "global.fh"
#include "msgids.fh"
c
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision x(nefc)      !< [Input] Vector `x`
      double precision ax(nefc)     !< [Output] Matrix-vector product
                                    !< `Ax`
c
      double precision zetai, zetaj, zetaij, pi, aij, dij, dum, xj
      double precision zero, one
      parameter (zero=0.0d+00, one=1.0d+00)
c
      double precision dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      double precision adiag,dsurf,dvol,srfmol,volmol,ptspatm
      common/hnd_cospar/dielec,dielecinf,screen,rsolv,zeta,gammas,swtol
      common/hnd_cosdat/adiag,dsurf,dvol,srfmol,volmol,ptspatm
c
      double precision util_erf
      external util_erf
c
c     Introduced a trivial replicated data parallelization of this
c     matrix-vector multiplication
c
      double precision r
      integer i, j
      r (i,j)=sqrt((efcc(1,i)-efcc(1,j))**2+
     1             (efcc(2,i)-efcc(2,j))**2+
     2             (efcc(3,i)-efcc(3,j))**2)
c     d (i  )=adiag/sqrt(efcs(i))
c
      pi=acos(-1.0d0)
      call dfill(nefc,0.0d0,ax,1)
      do i=ga_nodeid()+1,nefc,ga_nnodes()
         zetai=zeta*sqrt(ptspatm) /
     &         (ratm(efciat(i))*sqrt(2.0d0*pi))
         dum=zero
         do j=1,nefc
            if(j.eq.i) then
               aij=efcs(i)
            else
               zetaj=zeta*sqrt(ptspatm) /
     &               (ratm(efciat(j))*sqrt(2.0d0*pi))
               zetaij=zetai*zetaj/sqrt(zetai**2+zetaj**2)
               dij=r(i,j)
               if (dij.lt.1.0d-6) then
                 aij=2.0d0*zetaij*sqrt(1.0d0/pi)
               else
                 aij=util_erf(zetaij*dij)/dij
               endif
            endif
            xj=x(j)
            dum=dum+aij*xj
         enddo
         ax(i)=dum
      enddo
      call ga_dgop(msg_cosaxd,ax,nefc,'+')
c
      return
      end
C>
C> \brief On the fly vector-matrix multiplication
C>
C> This routine multiplies the current guess for the COSMO charges
C> with the matrix. The routine is replicated data parallel.
C> The matrix `A` is symmetric so we simply call the matrix-vector
C> product.
C>
      subroutine hnd_cosxad(nat,x,xa,nefc,efcc,efcs,efciat,ratm)
      implicit none
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision x(nefc)      !< [Input] Vector `x`
      double precision xa(nefc)     !< [Output] Vector-matrix product
                                    !< `xA`
c
      call hnd_cosaxd(nat,x,xa,nefc,efcc,efcs,efciat,ratm)
c
      return
      end
C>
C> \brief Solve a linear system of equations using an iterative
C> procedure
C>
C> This routine solves a linear system of equations \f$A\cdot x = b\f$
C> using an iterative procedure based on [1] page 70.
C>
C> ### References ###
C>
C> [1] W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,
C>     "Numerical Recipes: in Fortran 77", 2nd edition, Volume 1, 
C>     Cambridge University Press, 1992, ISBN: 0-521-43064-X,
C>     <a href="http://apps.nrbook.com/fortran/index.html">nrbook</a>.
C>
C> [2]  A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="http://dx.doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
      subroutine hnd_cosequ(nat,b,x,nefc,g,h,xi,xj,efcc,efcs,efciat,
     &                      ratm)
      implicit none
c
c     ----- solve a * x = b , using an iterative procedure       -----
c
c     ----- numerical recipes (p.70), cambridge university press -----
c          w.h.press, b.p.flannery, s.a.teukolsky, w.t.vetterling
c
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
c
      logical     dbug
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision b(nefc)      !< [Input] The RHS of Ax=b
      double precision x(nefc)      !< [Output] The solution of Ax=b
c
      double precision g(nefc) !< [Scratch]
      double precision h(nefc) !< [Scratch]
      double precision xi(nefc) !< [Scratch]
      double precision xj(nefc) !< [Scratch]
c
      double precision rp, bsq, anum, aden, rsq, gg, dgg, gam
      double precision zero, eps, eps2
      data zero   /0.0d+00/
      data eps    /1.0d-07/
c
      integer i, j, irst, iter
c
      dbug=.false.
      if(dbug) then
         write(luout,*) 'in -cosequ-'
         do i=1,nefc
            write(luout,9999) i,b(i)
         enddo
      endif
c
      eps2=nefc*eps*eps
      irst=0
   10 irst=irst+1
      call hnd_cosaxd(nat,x,xi,nefc,efcc,efcs,efciat,ratm)
      rp=zero
      bsq=zero
      do j=1,nefc
         bsq=bsq+b(j)*b(j)
         xi(j)=xi(j)-b(j)
         rp=rp+xi(j)*xi(j)
      enddo ! j
      call hnd_cosxad(nat,xi,g,nefc,efcc,efcs,efciat,ratm)
      do j=1,nefc
         g(j)=-g(j)
         h(j)= g(j)
      enddo ! j
      do iter=1,10*nefc
         call hnd_cosaxd(nat,h,xi,nefc,efcc,efcs,efciat,ratm)
         anum=zero
         aden=zero
         do j=1,nefc
            anum=anum+g(j)*h(j)
            aden=aden+xi(j)*xi(j)
         enddo ! j
         if(aden.eq.zero) then
            write(luout,*) 'very singular matrix'
            call errquit('hnd_cosequ: singular matrix',911,UERR)
         endif
         anum=anum/aden
         do j=1,nefc
            xi(j)=x(j)
            x(j)=x(j)+anum*h(j)
         enddo ! j
         call hnd_cosaxd(nat,x,xj,nefc,efcc,efcs,efciat,ratm)
         rsq=zero
         do j=1,nefc
            xj(j)=xj(j)-b(j)
            rsq=rsq+xj(j)*xj(j)
         enddo ! j
         if (iter.gt.10*nefc-min(10,nefc/10)) then
           if (ga_nodeid().eq.0) then
             write(luout,'(" hnd_cosequ: it,residue,thresh = ",
     &                     i5,2f12.5)')iter,rsq,bsq*eps2
           endif
         endif
         if(rsq.eq.rp.or.rsq.le.bsq*eps2) return
         if(rsq.gt.rp) then
            do j=1,nefc
               x(j)=xi(j)
            enddo ! j
            if(irst.ge.3) return
            go to 10
         endif
         rp=rsq
         call hnd_cosxad(nat,xj,xi,nefc,efcc,efcs,efciat,ratm)
         gg=zero
         dgg=zero
         do j=1,nefc
            gg=gg+g(j)*g(j)
            dgg=dgg+(xi(j)+g(j))*xi(j)
         enddo ! j
         if(gg.eq.zero) return
         gam=dgg/gg
         do j=1,nefc
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
         enddo ! j
      enddo ! iter
      write(luout,*) 'too many iterations'
      call errquit('hnd_cosequ: too many iters',911,UERR)
      return
 9999 format(i8,f16.8)
      end
C>
C> \brief Direct linear system solver
C>
C> This is direct linear system solver. On return the matrix A has been
C> destroyed.
C>
      subroutine hnd_linequ(a,lda,b,n,ib,t,deta,ierr,nodcmp)
      implicit real*8 (a-h,o-z)
      dimension a(lda,1),b(n),ib(n),t(n)
c
c     ----- solve a * x = b , with x returned in b -----
c
      if(nodcmp.ne.1) go to 20
c
      call hnd_ludcmp(a,lda,n,ib,t,deta,ierr)
      do 10 j=1,n
   10 deta=deta*a(j,j)
c
   20 continue
      call hnd_lubksb(a,lda,n,ib,b)
c
      return
      end
c
      subroutine hnd_lubksb(a,lda,n,ib,b)
      implicit real*8 (a-h,o-z)
c
      dimension a(lda,1),ib(n),b(n)
      data zero /0.0d+00/
c
      ii=0
      do 12 i=1,n
         ll=ib(i)
         sum=b(ll)
         b(ll)=b(i)
         if(ii.ne.0) then
            do 11 j=ii,i-1
               sum=sum-a(i,j)*b(j)
   11       continue
         else if (sum.ne.zero) then
            ii=i
         endif
         b(i)=sum
   12 continue
      do 14 i=n,1,-1
         sum=b(i)
         if(i.lt.n) then
            do 13 j=i+1,n
               sum=sum-a(i,j)*b(j)
   13       continue
         endif
         b(i)=sum/a(i,i)
   14 continue
      return
      end
c
      subroutine hnd_ludcmp(a,lda,n,ib,vv,d,ierr)
      implicit real*8 (a-h,o-z)
c
      dimension a(lda,1),ib(n),vv(n)
      data tiny     /1.0d-20/
      data zero,one /0.0d+00,1.0d+00/
c
      ierr=0
      d=one
      do 12 i=1,n
         aamax=zero
         do 11 j=1,n
            if( abs(a(i,j)).gt.aamax) aamax= abs(a(i,j))
   11    continue
         if(aamax.eq.zero) then
            ierr=1
            return
         endif
         vv(i)=one/aamax
   12 continue
      do 19 j=1,n
         if(j.gt.1) then
            do 14 i=1,j-1
               sum=a(i,j)
               if(i.gt.1) then
                  do 13 k=1,i-1
                     sum=sum-a(i,k)*a(k,j)
   13             continue
                  a(i,j)=sum
               endif
   14       continue
         endif
         aamax=zero
         imax=0
         do 16 i=j,n
            sum=a(i,j)
            if(j.gt.1) then
               do 15 k=1,j-1
                  sum=sum-a(i,k)*a(k,j)
   15          continue
               a(i,j)=sum
            endif
            dum=vv(i)* abs(sum)
            if(dum.ge.aamax) then
               imax=i
               aamax=dum
            endif
   16    continue
         if(j.ne.imax) then
            do 17 k=1,n
               dum=a(imax,k)
               a(imax,k)=a(j,k)
               a(j,k)=dum
   17       continue
            d=-d
            vv(imax)=vv(j)
         endif
         ib(j)=imax
         if(j.ne.n) then
            if(a(j,j).eq.zero) a(j,j)=tiny
            dum=one/a(j,j)
            do 18 i=j+1,n
               a(i,j)=a(i,j)*dum
   18       continue
         endif
   19 continue
      if(a(n,n).eq.zero) a(n,n)=tiny
      return
      end
C>
C> \brief Compute the function \f$f(r)\f$
C>
C> Computes the function \f$f(r)\f$ as discussed with function
C> `hnd_coschg`. This function is Eq.(64) in [1].
C>
C> ### References ###
C>
C> [1] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="http://dx.doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
      double precision function cosff(r)
      implicit none
#include "errquit.fh"
      double precision r !< [Input] penetration fraction
c
      if (r.lt.0.0d0) then
        cosff = 0.0d0
      else if (r.gt.1.0d0) then
        cosff = 1.0d0
      else
        cosff = r**3*(10.0d0-15.0d0*r+6.0d0*r**2)
      endif
c
      return
      end
C>
C> \brief Compute the function \f$\frac{\partial f(r)}{\partial r}\f$
C>
C> Computes the function \f$\frac{\partial f(r)}{\partial r}\f$ as
C> discussed with function `hnd_coschg`. This function is the 
C> derivative of Eq.(64) in [1].
C>
C> ### References ###
C>
C> [1] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="http://dx.doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
      double precision function cosdff(r)
      implicit none
#include "errquit.fh"
      double precision r !< [Input] penetration fraction
c
      if (r.lt.0.0d0) then
        cosdff = 0.0d0
      else if (r.gt.1.0d0) then
        cosdff = 0.0d0
      else
        cosdff = 30.0d0*(r**2)*((1.0d0-r)**2)
      endif
c
      return
      end
