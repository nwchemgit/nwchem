C>
C> \ingroup cosmo
C> @{
C>
C> \file hnd_cosmo_lib.F
C> Library of HONDO routines used to implement COSMO
C>
C> \brief Setup the COSMO cavity surface
C>
      subroutine hnd_cosset(rtdb,nat,c,radius,geom)
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
#include "cosmo_params.fh"
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
      integer geom
c
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
      integer l_i140, i140
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
      if(.not.ma_push_get(mt_dbl,3*nat,"xyzatm",l_i10,i10))
     c     call errquit('hndcosset: not enuf mem',0,MA_ERR)
      if(.not.ma_push_get(mt_dbl,  nat,"ratm",l_i20,i20))
     c     call errquit('hndcosset: not enuf mem',1,MA_ERR)
      if(.not.ma_push_get(mt_int,  nat,"nspa",l_i30,i30))
     c     call errquit('hndcosset: not enuf mem',2,MA_ERR)
      if(.not.ma_push_get(mt_int,  nat,"nppa",l_i40,i40))
     c     call errquit('hndcosset: not enuf mem',3,MA_ERR)
      if(.not.ma_push_get(mt_int,3*mxface,"ijkfac",l_i50,i50))
     c     call errquit('hndcosset: not enuf mem',4,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*mxface,"xyzseg",l_i60,i60))
     c     call errquit('hndcosset: not enuf mem',5,MA_ERR)
      if(.not.ma_push_get(mt_int,  mxface,"ijkseg",l_i70,i70))
     c     call errquit('hndcosset: not enuf mem',6,MA_ERR)
      if(.not.ma_push_get(mt_log,  mxface*nat,"insseg",l_i80,i80))
     c     call errquit('hndcosset: not enuf mem',7,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*mxface*nat,"xyzspa",l_i90,i90))
     c     call errquit('hndcosset: not enuf mem',8,MA_ERR)
      if(.not.ma_push_get(mt_int,  mxface*nat,"ijkspa",l_i100,i100))
     c     call errquit('hndcosset: not enuf mem',9,MA_ERR)
      if(.not.ma_push_get(mt_int,  mxface*nat,"numpps",l_i110,i110))
     c     call errquit('hndcosset: not enuf mem',10,MA_ERR)
      if(.not.ma_push_get(mt_dbl,3*mxapex    ,"apex",l_i120,i120))
     c     call errquit('hndcosset: not enuf mem',11,MA_ERR)
      if(.not.ma_push_get(mt_dbl,  mxface*nat,"xyzff",l_i130,i130))
     c     call errquit('hndcosset: not enuf mem',12,MA_ERR)
      if(.not.ma_push_get(mt_dbl,mxface,"Lebedev w",l_i140,i140))
     c     call errquit('hndcosset: not enuf mem',13,MA_ERR)
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
     2                int_mb(i50),dbl_mb(i60),int_mb(i70),log_mb(i80),
     3                dbl_mb(i90),int_mb(i100),int_mb(i110),
     4                dbl_mb(i120),dbl_mb(i130),dbl_mb(i140),
     5                geom,rtdb)

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
     3                  xyzspa,ijkspa,numpps,apex,xyzff,wleb,
     4                  geom,rtdb)
      implicit none
c
#include "nwc_const.fh"
#include "cosmo_params.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
#include "cosmoP.fh"
#include "mafdecls.fh"
#include "util_params.fh"
#include "util.fh"
      integer rtdb, nat
      integer mxatm
      integer mxfac
      integer mxapx
      integer geom
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
      double precision  wleb(mxfac)
      double precision  energy
c
      logical    oprint_out
      logical    oprint_more
      logical    oprint_debug
c
      integer i, iat, lfac, lseg, ndiv, nfac, nseg
      integer mfac
      integer jat
c
      oprint_out = util_print("cosmo_cossrf_out",print_high).and.
     D     (ga_nodeid().eq.0)
      oprint_debug = util_print("cosmo_cossrf_debug",print_debug).and.
     D     (ga_nodeid().eq.0)
      oprint_more = util_print("cosmo_cossrf_more",print_never).and.
     D     (ga_nodeid().eq.0)
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

      if ((iangleb.lt.1).and.(fibonacci.lt.1)) then
        call hnd_cossph(nseg,nfac,ndiv,
     1                  ijkfac,xyzseg,ijkseg,mxfac,apex,mxapx,
     2                  dsurf,dvol,adiag)
      elseif (fibonacci.lt.1) then
        call cosmo_lebedev(iangleb, nseg, wleb, xyzseg)
        ijkseg(1:nseg) = 0
      else
        nseg = fibonacci
        call cosmo_fibonacci(nseg, xyzseg, dsurf, dvol)
        ijkseg(1:nseg) = 0
      endif

      if (thomson.gt.0) then
        call cosmo_thomson(nseg,xyzseg)
      endif

      ptspatm = dble(nseg)
c
c     ----- debug printing -----
c
      if(oprint_out) then
         write(luout,9999) nseg,nfac,ndiv,dsurf,dvol
         write(luout,9995) adiag
         if(oprint_more) then
            write(luout,9998)
            do lseg=1,nseg
               write(luout,9997) lseg,xyzseg(1,lseg),xyzseg(2,lseg),
     1                             xyzseg(3,lseg),ijkseg(  lseg)
            enddo
         endif
         if(oprint_debug.and. (iangleb.lt.1).and.(fibonacci.lt.1)) then
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
            if (do_cosmo_model.eq.DO_COSMO_KS) then
              ratm(iat)=(radius(iat)+rsolv)/cau2ang
            else
              ratm(iat)=radius(iat)/cau2ang
            endif
         endif
      enddo
c
c     ----- create -solvent accessible surface- of the molecule -----
c

      call hnd_cossas(nat,xyzatm,ratm,mxatm,
     1                nspa,nppa,xyzspa,ijkspa,
     2                nseg,nfac,xyzseg,ijkseg,insseg,
     3                numpps,xyzff,wleb,mxfac,geom,rtdb)
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
C> significantly (maxbem and rsolv are not used anymore, also the
C> elimination of point charges is no longer based on reducing the
C> segment surface area until it vanishes, instead the surface charge
C> of a segment is quenched with a switching function to eliminate the
C> contribution of a surface point).
C>
C> ## Popular demand ##
C>
C> Due to popular demand this routine can do either the original
C> Klamt-Schuurmann approach or the York-Karplus approach. The approach
C> used is dictated by the `do_cosmo_model` variable.
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="https://doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
      subroutine hnd_cossas(nat,xyzatm,ratm,mxatom,
     1                      nspa,nppa,xyzspa,ijkspa,
     2                      nseg,nfac,xyzseg,ijkseg,insseg,
     3                      numpps,xyzff,wleb,mxface,geom,
     4                      rtdb)
      implicit none
#include "cosmo_params.fh"
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "bq.fh"
#include "prop.fh"
cnew
#include "cosmoP.fh"
#include "util_params.fh"
#include "util.fh"
c
      integer rtdb    !< [Input] The RTDB handle
      integer nat     !< [Input] The actual number of atoms
      integer mxface  !< [Input] The maximum number of faces
      integer mxatom  !< [Input] The maximum number of atoms
      integer nseg    !< [Input] The actual number of segments
      integer nfac    !< [Input] The actual number of faces

      integer geom
c
      logical     dbug
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
      double precision   wleb(         mxface)
      double precision zero, one
      data one     /1.0d+00/
      integer l_efcc, k_efcc, l_efcs, k_efcs, l_efcz, k_efcz
      integer l_efclb, k_efclb, k_efciat, l_efciat
      integer l_efczz, k_efczz
      integer l_msrf,k_msrf
      double precision ratm_real,dij,dum,cavdsp,pi,zetai,zetaii
      integer m,mfac,mseg
      integer nefc,iat,jat,npp,i,iseg,ifac,ief,ipp

      integer l_dum1, l_dum2, l_attag, k_dum1, k_dum2, k_attag
      character(len=16) aname

      integer op, iat_new
      double precision vector(3), symfact
c
      double precision cosff
      external         cosff

      logical  sym_atom
      external sym_atom

      integer sym_number_ops, sym_center_map
      external sym_number_ops, sym_center_map

      logical geom_cart_get
      external geom_cart_get
c
c     MN solvation models -->
c
      logical do_cosmo_smd
c
c     <-- MN solvation models
c
      double precision dist, xi, yi, zi, xj, yj, zj, rin, rout, alphai
      parameter (alphai = 0.5d0)
      logical oprint_debug
      logical oprint_molsurf
c
      dist(xi,yi,zi,xj,yj,zj)=sqrt((xj-xi)**2+(yj-yi)**2+(zj-zi)**2)
      rin(iat)=ratm(iat)*(1.0d0-alphai*gammas*sqrt(0.25d0**minbem))
      rout(iat)=ratm(iat)*(1.0d0+(1.0d0-alphai)*gammas*
     &                     sqrt(0.25d0**minbem))
c
      oprint_debug = util_print("cosmo_cossas_debug",print_debug).and.
     D     (ga_nodeid().eq.0)
      oprint_molsurf = util_print("cosmo_mol_surface",print_never).and.
     D     (ga_nodeid().eq.0)
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
          if (do_cosmo_model.eq.DO_COSMO_KS) then
            write(luout,9997) iat,xyzatm(1,iat),xyzatm(2,iat),
     1                                       xyzatm(3,iat),
     2                    (ratm(iat)*cau2ang-rsolv)
          else if (do_cosmo_model.eq.DO_COSMO_YK) then
            write(luout,9997) iat,xyzatm(1,iat),xyzatm(2,iat),
     1                                       xyzatm(3,iat),
     2                    (ratm(iat)*cau2ang)
          endif
        enddo
      endif
c
c     ----- clear arrays ..... -----
c
      do iat=1,nat
         do i=1,mxface
            ijkspa(i,iat)=0
            numpps(i,iat)=0
            xyzff(i,iat)=0d0
         enddo
      enddo
c
c     ----- sift through atomic centers and decide if a face -----
c           belongs to the -sas- or is inside the molecule.
c
      do iat=1,nat
         if(ratm(iat).ne.0d0) then
            if (sym_number_ops(geom).gt.0) then
              if (.not.sym_atom(geom, iat, symfact)) cycle

              ! Identity operator first
              do iseg=1,nseg
                ijkspa(iseg,iat)=ijkseg(iseg)
                xyzff(iseg,iat)=one
                vector(:)=xyzseg(:,iseg)
                xyzspa(:,iseg,iat)=vector(:)*ratm(iat)+xyzatm(:,iat)
              enddo
              do ifac=1,nfac
                ijkspa(ifac+nseg,iat) = ijkseg(ifac+nseg)
                vector(:)=xyzseg(:,ifac+nseg)
                xyzspa(:,ifac+nseg,iat)=vector(:)*ratm(iat)
     $                                 +xyzatm(:,iat)
              enddo

              ! Loop over rest of operators
              do op=1, sym_number_ops(geom)
                iat_new = sym_center_map(geom, iat, op)
                if (iat_new.eq.iat) cycle
                do iseg=1,nseg
                  ijkspa(iseg,iat_new)=ijkseg(iseg)
                  xyzff(iseg,iat_new)=one
                  call sym_apply_op(geom, op, xyzseg(1,iseg),vector)
                  xyzspa(1:3,iseg,iat_new)=vector(:)*ratm(iat_new) +
     $                                     xyzatm(:,iat_new)
                enddo
                do ifac=1,nfac
                  ijkspa(ifac+nseg,iat_new) = ijkseg(ifac+nseg)
                  call sym_apply_op(geom,op,xyzseg(1,ifac+nseg),vector)
                  xyzspa(1:3,ifac+nseg,iat_new)=vector(:)*ratm(iat_new) 
     $                                         + xyzatm(:,iat_new)
                enddo
              enddo
            else
              do iseg=1,nseg
                 ijkspa(iseg,iat)=ijkseg(iseg)
                 xyzff(iseg,iat)=one
                 do m=1,3
                    xyzspa(m,iseg,iat)=xyzseg(m,iseg)*ratm(iat)
     1                                +xyzatm(m,iat)
                 enddo
              enddo
              do ifac=1,nfac
                 ijkspa(ifac+nseg,iat)=ijkseg(ifac+nseg)
                 do m=1,3
                   xyzspa(m,ifac+nseg,iat)=xyzseg(m,ifac+nseg)*ratm(iat)
     1                                    +xyzatm(m,iat)
                 enddo
              enddo
            endif
         endif
      enddo

      do iat=1,nat
         if(ratm(iat).ne.0d0) then
            if(oprint_debug) then
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
               if (do_cosmo_model.eq.DO_COSMO_KS) then
                 if(jat.ne.iat.and.(dij.lt.(ratm(iat)+ratm(jat)))) then
                   do ifac=1,nfac
                     dum=dist(xyzspa(1,ifac+nseg,iat),
     1                        xyzspa(2,ifac+nseg,iat),
     2                        xyzspa(3,ifac+nseg,iat),
     3                        xyzatm(1,jat),
     4                        xyzatm(2,jat),
     5                        xyzatm(3,jat))
                     if(dum.lt.ratm(jat)) then
                        ijkspa(ifac+nseg,iat)=0
                     endif
                   enddo
                 endif
               else if (do_cosmo_model.eq.DO_COSMO_YK) then
                 if((jat.ne.iat).and.(ratm(jat).ne.0d0)
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
               endif
            enddo
            if(dbug.and.ga_nodeid().eq.0.and.
     $         iangleb.lt.1.and.fibonacci.lt.1) then
               write(luout,9996) iat
               write(luout,9995) (ijkspa(ifac+nseg,iat),ifac=1,nfac)
            endif
c
c     ----- check segments belonging to -sas- -----
c
            if (do_cosmo_model.eq.DO_COSMO_KS) then
              do ifac=1,nseg+nfac
                insseg(ifac,iat)=.true.
              enddo
              do ifac=1,nfac
                iseg=ijkspa(ifac+nseg,iat)
                if(iseg.ne.0) then
                   insseg(ifac+nseg,iat)=.false.
                   insseg(     iseg,iat)=.false.
                endif
              enddo
            else if (do_cosmo_model.eq.DO_COSMO_YK) then
              do iseg=1,nseg
                insseg(iseg,iat)=.not.(xyzff(iseg,iat).ge.swtol)
              enddo
            endif
            if(dbug.and.ga_nodeid().eq.0) then
               write(luout,9994) iat
               if (do_cosmo_model.eq.DO_COSMO_KS) then
                 write(luout,9993) (insseg(ifac,iat),ifac=1,nseg+nfac)
               else if (do_cosmo_model.eq.DO_COSMO_YK) then
                 write(luout,9993) (insseg(iseg,iat),iseg=1,nseg)
               endif
            endif
            mseg=0
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) mseg=mseg+1
            enddo
            mfac=0
            if (do_cosmo_model.eq.DO_COSMO_KS) then
              do ifac=1,nfac
                if(.not.insseg(ifac+nseg,iat)) mfac=mfac+1
              enddo
            endif
            nspa(iat)=mseg
            nppa(iat)=mfac
c
c           ----- surface area of segments -----
c
            if (do_cosmo_model.eq.DO_COSMO_KS) then
              do iseg=1,nseg
                numpps(iseg,iat)=0
              enddo
              do ifac=1,nfac
                iseg=ijkspa(ifac+nseg,iat)
                if(iseg.ne.0) numpps(iseg,iat)=numpps(iseg,iat)+1
              enddo
            endif
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
      if(oprint_debug) then
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
         if(ratm(iat).ne.0d0) then
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
      if(.not.ma_push_get(mt_dbl,nefc,'cosmo zeta',l_efczz,k_efczz))
     & call errquit('cosmo_cossas malloc k_efczz failed',911,MA_ERR) 
      if(.not.ma_push_get(mt_byte,nefc*8,'cosmo tags',l_efclb,k_efclb))
     &     call errquit('cosmo_cossas malloc k_tag failed',911,MA_ERR)
      if(oprint_molsurf)  then
         if(.not.ma_push_get(mt_dbl,nefc,'molsurf',l_msrf,k_msrf))
     &        call errquit('cosmo_cossas malloc k_tag failed',1,MA_ERR)
         call dfill(nefc,0.0d0,dbl_mb(k_msrf),1)
      endif

      if(oprint_debug) then
        if(.not.ma_push_get(mt_dbl,nat*3,'coord',l_dum1,k_dum1))
     &   call errquit('cosmo_cossas malloc dum1 failed',911,MA_ERR)
        if(.not.ma_push_get(mt_dbl,nat,'cosmo z',l_dum2,k_dum2)) call
     &    errquit('cosmo_cossas malloc k_dum2 failed',911,MA_ERR)
        if(.not.ma_push_get(mt_byte,nat*16,'tags',l_attag,k_attag)) call
     &    errquit('cosmo_cossas malloc k_attag failed',911,MA_ERR)
c
        if(.not.geom_cart_get(geom,nat,byte_mb(k_attag),
     &     dbl_mb(k_dum1),dbl_mb(k_dum2))) call errquit
     $      (' cosmo_cossas: geom_cart_get failed.',911, GEOM_ERR)
      endif
c
c     ----- save coordinates of surface points -----
c           save segment surfaces
c           save segment to atom mapping
c
      srfmol=0d0
      volmol=0d0
      ief   =0


      do iat=1,nat
         if(ratm(iat).ne.0d0) then
            ratm_real=-1d99
            if (do_cosmo_model.eq.DO_COSMO_KS) then
              ratm_real=ratm(iat)-rsolv/cau2ang
            else if (do_cosmo_model.eq.DO_COSMO_YK) then
              ratm_real=ratm(iat)
            endif 
            
            do iseg=1,nseg
               if(.not.insseg(iseg,iat)) then
                  ief=ief+1

                  if (do_cosmo_model.eq.DO_COSMO_KS) then
                    do i=1,3
                      dbl_mb(k_efcc+3*(ief-1)+i-1)=
     1                     (xyzspa(i,iseg,iat)-xyzatm(i,iat))*
     2                     ratm_real/ratm(iat) + xyzatm(i,iat)
                    enddo
                  elseif (do_cosmo_model.eq.DO_COSMO_YK) then
c                    do i=1,3
c                     dbl_mb(k_efcc+3*(ief-1)+i-1)=xyzspa(i,iseg,iat)
                     dbl_mb(k_efcc+3*(ief-1))  =xyzspa(1,iseg,iat)
                     dbl_mb(k_efcc+3*(ief-1)+1)=xyzspa(2,iseg,iat)
                     dbl_mb(k_efcc+3*(ief-1)+2)=xyzspa(3,iseg,iat)
c                    enddo
                  endif
                  ipp=numpps(iseg,iat)
                  if (do_cosmo_model.eq.DO_COSMO_KS) then
                    dbl_mb(k_efcs+ief-1) = dble(ipp)*dsurf*ratm_real**2
                    srfmol = srfmol + dble(ipp)*dsurf*ratm_real**2
                    volmol = volmol + dble(ipp)*dvol *ratm_real**3
                  else if (do_cosmo_model.eq.DO_COSMO_YK) then
c
c                   --- eval eq.(67) from [2] ---
c
                    if (iangleb.lt.1) then
                      zetai=zeta/sqrt(dsurf)/ratm_real
                      zetaii=zetai*sqrt(2d0/pi)
                      dbl_mb(k_efcs+ief-1) = (zetaii/xyzff(iseg,iat))
c     769
                      if(oprint_molsurf)  then
                         dbl_mb(k_msrf+ief-1)=
     x                        xyzff(iseg,iat)*dsurf*ratm_real**2
                      endif
                      srfmol = srfmol+xyzff(iseg,iat)*dsurf*ratm_real**2
                      volmol = volmol+xyzff(iseg,iat)*dvol *ratm_real**3
                    else
                      zetai=zeta/(ratm_real*sqrt(wleb(iseg)))
                      zetaii=zetai*sqrt(2d0/pi)
                      dbl_mb(k_efcs+ief-1) = zetaii/xyzff(iseg,iat)
c 775
                      if(oprint_molsurf)  then
                         dbl_mb(k_msrf+ief-1)=
     x                        xyzff(iseg,iat)*wleb(iseg)*ratm_real**2
                      endif
                      srfmol = srfmol + xyzff(iseg,iat)*
     $                                  wleb(iseg)*ratm_real**2
                      volmol = volmol + xyzff(iseg,iat)*
     $                                  wleb(iseg)*ratm_real**3/3d0
                    endif
                  endif
                  int_mb(k_efciat+ief-1)=iat
                  dbl_mb(k_efczz+ief-1) =zetai
               endif
            enddo
         endif
      enddo
      srfmol=srfmol*(cau2ang**2)
      volmol=volmol*(cau2ang**3)
c
      if(oprint_debug) then
        write(luout,*)
        write(luout,'(" --- Geometry including SAS --- ")')
        write(luout,*) nat+nefc
        write(luout,*)
        do iat=1,nat
          call util_set_ma_char_arrays(16,
     $                                byte_mb(k_attag+16*(iat-1)),
     $                                aname)
          write(luout,'(A10,3F12.6)') aname(1:2),
     $                                (xyzatm(i,iat)*cau2ang,i=1,3)
        enddo
        do ief=1,nefc
          write(luout,'(A10,3F12.6)') "X",
     1    (dbl_mb(k_efcc+3*(ief-1)+i-1)*cau2ang,i=1,3)
        enddo
        write(luout,'(" --- End of Geometry including SAS --- ")')
        write(luout,*)
        if (.not.ma_chop_stack(l_dum1))
     $    call errquit('hnd_cossas: could not chop stack',0,MA_ERR)   
      endif
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
c
c MN solvation models -->
c
c      if(ga_nodeid().eq.0) then
c         write(luout,9981) cavdsp
c      endif
      if (.not.
     $ rtdb_get(rtdb,'cosmo:do_cosmo_smd',mt_log,1,do_cosmo_smd))
     $ call errquit('hnd_cossas: cannot get do_cosmo_smd from rtdb',
     $ 0,rtdb_err)
      if(ga_nodeid().eq.0) then
       if (.not.do_cosmo_smd) write(luout,9981) cavdsp
      endif
c
c <-- MN solvation models
c
c     ----- print segment surfaces -----
c
      if(oprint_molsurf) then
        if (do_cosmo_model.eq.DO_COSMO_KS) then
           write(luout,9989)
           do ief=1,nefc
              write(luout,9988) ief,dbl_mb(k_efcs+ief-1),
     &                          int_mb(k_efciat+ief-1)
           enddo
        else if (do_cosmo_model.eq.DO_COSMO_YK) then
           write(luout,*) ' Contributions to the molecular surface '
           do ief=1,nefc
              write(luout,*) ief,
     W           dbl_mb(k_msrf+ief-1)
           enddo

        endif
      endif
c
      do ief=1,nefc
         dbl_mb(k_efcz+ief-1)=0d0
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
      if(.not.rtdb_put(rtdb,'cosmo:efczz',mt_dbl,nefc,dbl_mb(k_efczz)))
     &   call errquit('hnd_cossas: rtdb put failed for efczz ',914,
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
      if(.not.ma_chop_stack(l_efcc)) call
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
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
      subroutine hnd_cossph(nseg,nfac,ndiv,
     1                  ijkfac,xyzseg,ijkseg,mxface,apex,mxapex,
     2                  dsurf_in,dvol_in,adiag_in)
      implicit double precision (a-h,o-z)
#include "cosmo_params.fh"
#include "global.fh"
#include "stdio.fh"
cnew
#include "cosmoP.fh"
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
c        do lface=1,maxfac(maxbem-1)-minfac(minbem+1)+1
c          ijkseg(nseg+lface) = 0
c        enddo
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
      if (do_cosmo_model.eq.DO_COSMO_YK) nfac =nseg
      ndiv =nfac/nseg
      dsurf_in=srf/dble(nfac)
      dvol_in =vol/dble(nfac)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9986) nseg,nfac,ndiv,dsurf_in,dvol_in
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
         adiag_in=sqrt(four*pi)*fac
         if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
            write(luout,9983) mseg,adiag_in,fac,dble(nseg),sum
         endif
         avgdia=avgdia+adiag_in
         avgfac=avgfac+fac
      enddo
      adiag_in=avgdia/dble(nseg)
      fac  =avgfac/dble(nseg)
      if(some.or.out.or.dbug.and.ga_nodeid().eq.0) then
         write(luout,9978)      adiag_in,fac               
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
C>
C> \brief Evaluate COSMO related energy terms
C>
C> Based on the COSMO charges a number of energy terms can be evaluated,
C> including:
C>
C> - the nuclear - COSMO charge interaction energy
C>
C> - the electron - COSMO charge interaction energy
C>
C> - the total solute charge density - COSMO charge interaction energy
C>
C> - the COSMO charge - COSMO charge interaction energy
C>
C> These terms are inherent in Ref.[1] Eq.(8) and in Ref.[2] Eq.(39).
C>
C> ### References ###
C>
C> [1] A. Klamt, G. Sch&uuml;&uuml;rmann,
C>     "COSMO: a new approach to dielectric screening in solvents with
C>      explicit expressions for the screening energy and its gradient",
C>     <i>J. Chem. Soc., Perkin Trans. 2</i>, 1993, pp 799-805, DOI:
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
C> [2] D.M. York, M. Karplus,
C>     "A smooth solvation potential based on the conductor-like
C>      screening model", <i>J. Phys. Chem. A</i> (1999) <b>103</b>,
C>     pp 11060-11079, DOI:
C>     <a href="https://doi.org/10.1021/jp992097l">
C>     10.1021/jp992097l</a>.
C>
      subroutine hnd_cos_energy(nat,nefc,chgscr,efcc,efcs,efcz,efczz,
     &           efciat,ratm,catm,zatm,pot,allefc,atmefc,elcefc,efcefc)
      implicit none
c
#include "cosmoP.fh"
#include "cosmo_params.fh"
#include "global.fh"
c
      integer nat  !< [Input] The number of electrons
      integer nefc !< [Input] The number of COSMO charges
      integer efciat(nefc) !< [Input] Mapping of COSMO charges to the
                           !< atom that generated the corresponding 
                           !< solvent accessible surface area segment
c
      double precision chgscr !< [Input] The dielectric screening
                              !< factor \f$f(\epsilon)\f$ (see Ref.[2]
                              !< Eq.(46) with \f$\epsilon_1=1\f$ and
                              !< \f$\epsilon_2=\epsilon\f$)
      double precision efcc(3,nefc) !< [Input] The COSMO charge 
                                    !< coordinates
      double precision efcs(nefc)   !< [Input] The COSMO charge surface
                                    !< area \f$|S_\mu|\f$ in Ref.[1]
                                    !< (see e.g. Eq.(7b)) or Ref.[2]
                                    !< Eq.(67).
      double precision efcz(nefc)   !< [Input] The COSMO charges
      double precision efczz(nefc)
      double precision ratm(nat)    !< [Input] The atomic radii
      double precision catm(3,nat)  !< [Input] The atomic coordinates
      double precision zatm(nat)    !< [Input] The nuclear charges
      double precision pot(nefc)    !< [Input] The molecular potential
                                    !< at the COSMO charge positions
c
      double precision allefc !< [Output] The total solute charge 
                              !< density - COSMO charge interaction
                              !< energy
      double precision atmefc !< [Output] The nuclear - COSMO charge
                              !< interaction energy
      double precision elcefc !< [Output] The electron - COSMO charge
                              !< interaction energy
      double precision efcefc !< [Output] The COSMO charge - COSMO
                              !< charge interaction energy
c
c
      integer ief, jef !< Counters over COSMO charges
      integer iat      !< Counter over atoms
c
      double precision xi,  yi,  zi,  qi
      double precision xj,  yj,  zj,  qj
      double precision aii, aij, bij, dij
c
      double precision zetai, zetaj, zetaij
c
      double precision zero, one, two
      parameter (zero=0.0d0, one=1.0d0, two=2.0d0)
      double precision pi
      double precision derf
c
      pi    =acos(-1.0d0)
      allefc=zero
      atmefc=zero
      efcefc=zero
      do jef=1+ga_nodeid(),nefc,ga_nnodes()
         xj=efcc(1,jef)
         yj=efcc(2,jef)
         zj=efcc(3,jef)
         qj=efcz(  jef)
c
         allefc=allefc+qj*pot(jef)
c
         do iat=1,nat
            xi=catm(1,iat)
            yi=catm(2,iat)
            zi=catm(3,iat)
            qi=zatm(  iat)
            dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
            bij=one/dij
            atmefc=atmefc+qi*bij*qj
         enddo
      enddo
c
      if (do_cosmo_model.eq.DO_COSMO_KS) then
        do jef=1+ga_nodeid(),nefc,ga_nnodes()
           xj=efcc(1,jef)
           yj=efcc(2,jef)
           zj=efcc(3,jef)
           qj=efcz(  jef)
           efcefc=efcefc+qj*adiag*qj/sqrt(efcs(jef))
           do ief=jef+1,nefc
              qi=efcz(  ief)
              xi=efcc(1,ief)
              yi=efcc(2,ief)
              zi=efcc(3,ief)
              dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
              aij=one/dij
              efcefc=efcefc+2*qi*aij*qj
           enddo
        enddo
      else if (do_cosmo_model.eq.DO_COSMO_YK) then
        do jef=1+ga_nodeid(),nefc,ga_nnodes()
           xj=efcc(1,jef)
           yj=efcc(2,jef)
           zj=efcc(3,jef)
           qj=efcz(  jef)
           zetaj=efczz(jef)
           efcefc=efcefc+qj*efcs(jef)*qj
c
           do ief=jef+1,nefc
              zetai=efczz(ief)
              xi=efcc(1,ief)
              yi=efcc(2,ief)
              zi=efcc(3,ief)
              qi=efcz(  ief)
              dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
              zetaij=zetai*zetaj/sqrt(zetai**2+zetaj**2)
              if (dij.lt.1.0d-8) then
                aij=2.0d0*zetaij/sqrt(pi)
              else
                aij=derf(zetaij*dij)/dij
              endif
              efcefc=efcefc+2*qi*aij*qj
           enddo
        enddo
      endif
      call ga_dgop(0,efcefc,1,'+')
      call ga_dgop(0,allefc,1,'+')
      call ga_dgop(0,atmefc,1,'+')

      efcefc= efcefc/(two*chgscr)
      elcefc= allefc-atmefc
c
      end
C>
C> \brief Compute matrix A
C> 
C> Compute matrix `A` as needed for the in-memory COSMO
C> charge fitting.
C>
      subroutine hnd_cosmata(nat,nefc,efcc,efcs,efczz,efciat,ratm,a)
      implicit none
#include "cosmoP.fh"
#include "cosmo_params.fh"
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      double precision efczz(nefc)
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
c
      double precision derf
      double precision factor, factor2

      double precision pi
      pi = acos(-1.0d0)
c
      call dfill(nefc**2,0d0,a,1)
c
      if (do_cosmo_model.eq.DO_COSMO_KS) then
        do jef=1,nefc
          xj=efcc(1,jef)
          yj=efcc(2,jef)
          zj=efcc(3,jef)
          do ief=1,nefc
            if(ief.eq.jef) then
              aii=adiag/sqrt(efcs(ief))
              a(ief,jef)=aii
            else
              xi=efcc(1,ief)
              yi=efcc(2,ief)
              zi=efcc(3,ief)
              dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
              a(ief,jef)=one/dij
            endif
          enddo
        enddo
      else if (do_cosmo_model.eq.DO_COSMO_YK) then
        factor = zeta*dsqrt(ptspatm/(4d0*pi))
        factor2 = 2d0/sqrt(pi)
        do jef=1,nefc
          zetaj=efczz(jef)
          xj=efcc(1,jef)
          yj=efcc(2,jef)
          zj=efcc(3,jef)
          a(jef,jef)=efcs(jef)
          do ief=1,jef-1
            zetai=efczz(ief)
            zetaij=zetai*zetaj/sqrt(zetai**2+zetaj**2)
            xi=efcc(1,ief)
            yi=efcc(2,ief)
            zi=efcc(3,ief)
            dij=sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
            if (dij.lt.1.0d-5) then
              a(ief,jef)=factor2*zetaij*(1d0 - (zetaij*dij)**2/3d0)
            else
              a(ief,jef)=derf(zetaij*dij)/dij
            endif
            a(jef,ief) = a(ief,jef)
          enddo
        enddo
      endif
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
      subroutine hnd_cosaxd(nat,x,ax,nefc,efcc,efcs,efczz,efciat,ratm)
      implicit none
#include "global.fh"
#include "msgids.fh"
#include "cosmoP.fh"
#include "cosmo_params.fh"
c
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      double precision efczz(nefc)
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
c
      double precision derf
c
c     Introduced a trivial replicated data parallelization of this
c     matrix-vector multiplication
c
      double precision r, d, factor, factor2
      integer i, j
      r (i,j)=sqrt((efcc(1,i)-efcc(1,j))**2+
     1             (efcc(2,i)-efcc(2,j))**2+
     2             (efcc(3,i)-efcc(3,j))**2)
      d (i  )=adiag/sqrt(efcs(i))
c
      pi=acos(-1.0d0)
      call dfill(nefc,0.0d0,ax,1)
c
      if (do_cosmo_model.eq.DO_COSMO_KS) then
        do i=ga_nodeid()+1,nefc,ga_nnodes()
          ax(i) = ax(i) + d(i)*x(i)
          do j=i+1,nefc
            aij = one/r(i,j)
            ax(i) = ax(i) + x(j)*aij
            ax(j) = ax(j) + x(i)*aij
          enddo
        enddo
      else if (do_cosmo_model.eq.DO_COSMO_YK) then
        factor = zeta*dsqrt(ptspatm/(4d0*pi))
        factor2 = 2d0/sqrt(pi)
        do i=ga_nodeid()+1,nefc,ga_nnodes()
          zetai = efczz(i)
          ax(i) = ax(i) + efcs(i)*x(i)
          do j=i+1,nefc
            zetaj = efczz(j)
            zetaij= zeta*zetaj/sqrt(zetai**2+zetaj**2)
            dij=r(i,j)
            if (dij.lt.1.0d-8) then
              aij=factor2*zetaij*(1d0 - (zetaij*dij)**2/3d0)
            else
              aij=derf(zetaij*dij)/dij
            endif
            ax(i)=ax(i)+aij*x(j)
            ax(j)=ax(j)+aij*x(i)
          enddo
        enddo
      endif
c
      call ga_dgop(msg_cosaxd,ax(1),nefc,'+')
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
      subroutine hnd_cosxad(nat,x,xa,nefc,efcc,efcs,efczz,efciat,ratm)
      implicit none
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      double precision efczz(nefc)
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision x(nefc)      !< [Input] Vector `x`
      double precision xa(nefc)     !< [Output] Vector-matrix product
                                    !< `xA`
c
      call hnd_cosaxd(nat,x,xa,nefc,efcc,efcs,efczz,efciat,ratm)
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
C>     <a href="https://doi.org/10.1039/P29930000799">
C>     10.1039/P29930000799</a>.
C>
      subroutine hnd_cosequ(nat,b,x,nefc,g,h,xi,xj,efcc,efcs,efczz,
     &                      efciat,ratm)
      implicit none
c
c     ----- solve A * x = b , using an iterative procedure       -----
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
      double precision efczz(nefc)
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
      call hnd_cosaxd(nat,x,xi,nefc,efcc,efcs,efczz,efciat,ratm)
      rp=zero
      bsq=zero
      do j=1,nefc
         bsq=bsq+b(j)*b(j)
         xi(j)=xi(j)-b(j)
         rp=rp+xi(j)*xi(j)
      enddo ! j
      call hnd_cosxad(nat,xi,g,nefc,efcc,efcs,efczz,efciat,ratm)
      do j=1,nefc
         g(j)=-g(j)
         h(j)= g(j)
      enddo ! j
      do iter=1,10*nefc
         call hnd_cosaxd(nat,h,xi,nefc,efcc,efcs,efczz,efciat,ratm)
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
         call hnd_cosaxd(nat,x,xj,nefc,efcc,efcs,efczz,efciat,ratm)
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
         call hnd_cosxad(nat,xj,xi,nefc,efcc,efcs,efczz,efciat,ratm)
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
C> This is direct linear system solver, solving \f$Ax=b\f$ for \f$x\f$.
C> On return the matrix \f$A\f$ has been destroyed, and \f$b\f$ has been
C> overwritten with the answer \f$x\f$.
C>
      subroutine hnd_linequ(a,lda,b,n,ib,t,ierr,nodcmp)
#ifdef USE_OPENMP
      use omp_lib
#endif
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer lda !< [Input] The leading dimension of \f$A\f$
      integer n   !< [Input] The dimension of matrix \f$A\f$
      double precision a(lda,n) !< [In/Output] On input the matrix
      !< \f$A\f$; On output the matrix has been destroyed
      double precision b(n,*) !< [In/Output] On input the right-hand-side 
      !< \f$b\f$; on output the solution \f$x\f$
      double precision t(n)
      integer ib(n)
      integer ierr !< [Output] Error code
      integer nodcmp !< [Input] Flag, if \f$\mathrm{nodcmp}\le 1\f$
      !< then skip the decomposition
      integer j !< Counter
      integer lwork, iwork, kwork
      double precision t0
      integer iMaxThreads

#ifdef USE_OPENMP
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif
c
c     ----- solve a * x = b , with x returned in b -----
c
      ierr = 0

      ! Save diagonal in case dposv fails
      do j=1,n
        t(j) = a(j,j)
      enddo

      call dposv('u', n, 2, a, lda, b, n, ierr)

      if (ierr.ne.0) then

        lwork = 64*n
        if (.not.ma_push_get(mt_dbl,lwork,"dsysv work",iwork,kwork))
     $  call errquit("hnd_linequ: could not allocate work",lwork,MA_ERR)

        ! Restore diagonal
        do j=1,n
          a(j,j) = t(j)
        enddo

        call dsysv('l',n,2,a,lda,ib,b,n,dbl_mb(kwork),lwork,ierr)

        if (.not.ma_chop_stack(iwork))
     $    call errquit("hnd_linequ: could not chop stack",0,MA_ERR)   

      endif
c
#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif
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
C>     <a href="https://doi.org/10.1021/jp992097l">
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
C>     <a href="https://doi.org/10.1021/jp992097l">
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

      subroutine hnd_cg(nat,b,x,nefc,y,r,p,ap,efcc,efcs,efczz,efciat,
     $                  ratm)
#ifdef USE_OPENMP
      use omp_lib
#endif
      implicit none
#include "cosmoP.fh"
#include "cosmo_params.fh"
#include "stdio.fh"
#include "global.fh"
      logical     dbug
      integer          nat  !< [Input] The number of atoms
      integer          nefc !< [Input] The number of surface charges
      double precision efcc(3,nefc) !< [Input] The surface charge coords
      double precision efcs(nefc)   !< [Input] The surface areas
      double precision efczz(nefc)
      integer          efciat(nefc) !< [Input] The atom of a surface 
                                    !< charge
      double precision ratm(nat)    !< [Input] The atom radii
      double precision b(nefc)      !< [Input] The RHS of Ax=b
      double precision x(nefc)      !< [Output] The solution of Ax=b
c
      double precision y(nefc) !< [Scratch]
      double precision ap(nefc) !< [Scratch]
      double precision p(nefc) !< [Scratch]
      double precision r(nefc) !< [Scratch]

      double precision rtol2, rsqnew, alpha, beta
      double precision xnorm, facold, facnew
c
      integer iter, maxiter
c
      data rtol2   / 1d-13 /
      data maxiter / 100 /
      integer iMaxThreads, i

#ifdef USE_OPENMP
      iMaxThreads = omp_get_max_threads()
      call util_blas_set_num_threads(iMaxThreads)
#endif

      iter = 0

      xnorm = dot_product(x,x)
      if (xnorm.eq.0d0) then
        r(:) = b(:)
      else
        call hnd_cosaxd(nat,x,r,nefc,efcc,efcs,efczz,efciat,ratm)
        r(:) = b(:) - r(:)
      endif

      rsqnew = dot_product(r,r)
      if (rsqnew.lt.rtol2) return


      if (do_cosmo_model.eq.DO_COSMO_YK) then
        p(:) = r(:)/efcs(:)
      elseif (do_cosmo_model.eq.DO_COSMO_KS) then
        p(:) = dsqrt(efcs(:))*r(:)/adiag
      endif

      facnew = dot_product(r,p)

      do while(iter.le.nefc)

        iter = iter + 1
        call hnd_cosaxd(nat,p,ap,nefc,efcc,efcs,efczz,efciat,ratm)

        alpha = facnew/dot_product(p,ap)
        call daxpy(nefc,alpha,p,1,x,1)
        call daxpy(nefc,-alpha,ap,1,r,1)

        rsqnew = dot_product(r,r)
        if (rsqnew.lt.rtol2) goto 2000

        if (do_cosmo_model.eq.DO_COSMO_YK) then
          y(:) = r(:)/efcs(:)
        elseif (do_cosmo_model.eq.DO_COSMO_KS) then
          y(:) = dsqrt(efcs(:))*r(:)/adiag
        endif

        facold = facnew
        facnew = dot_product(r,y)
        beta = facnew/facold

        p(:) = y(:) + beta*p(:)

      enddo

      if (ga_nodeid().eq.0) then
        write(luout,1000) dsqrt(rsqnew)
      endif
 1000 format(2x,"Warning! hnd_cg residual: ",G8.3)

 2000 continue

#ifdef USE_OPENMP
      call util_blas_set_num_threads(1)
#endif
      return
      end

C>
C> @}
c $Id$
