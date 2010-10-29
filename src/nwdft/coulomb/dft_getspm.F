      Subroutine dft_getspm(geom,lmax,finest,iga_dens,basis)

C$Id$
      implicit none
#include "errquit.fh"
      integer lmax               ! max value of ang momentum
      integer iga_dens           ! GA handle for DM
      integer basis              ! basis handle
      integer finest             ! finest level for CMM
      integer  geom
C
C     local
C
      double precision centerl(3),q
      double precision total_nuc_charge
      integer natoms,iatom
      integer max_at_bf,max_at_bf2
      integer lcoord,icoord,lcharge,icharge,itags,ltags
      integer ldipole,idipole,lPmat,iPmat
      integer lxcrd,ixcrd,lycrd,iycrd,lzcrd,izcrd
      integer i,ihi,istart,ilo,icent

#include "bas.fh"
#include "geom.fh"

#include "mafdecls.fh"
#include "stdio.fh"


      if(.not. geom_ncent(geom, natoms))
     &  call errquit('dft_scf: geom_ncent failed',0, GEOM_ERR)

      if(.not.Ma_Alloc_Get(MT_Dbl,natoms*3,'coordinates',lcoord,icoord))
     &  call errquit('dft_getspm: cannot allocate coordinates',0,
     &       MA_ERR)
      if(.not.Ma_Push_Get(MT_Dbl,natoms,'charges',lcharge,icharge))
     &  call errquit('dft_getspm: cannot allocate charges',0,
     &       MA_ERR)
      if(.not.Ma_Push_Get(MT_Byte,natoms*16,'center tags',ltags,itags))
     &  call errquit('dft_getspm: cannot allocate center tags',0,
     &       MA_ERR)

      if(.not. geom_cart_get(geom, natoms, Byte_MB(itags),
     &  Dbl_MB(icoord), Dbl_MB(icharge)))
     &  call errquit('dft_getspm: geom_cart_get failed',1, GEOM_ERR)

c**** 
c**** compute center of charge
c**** 

      centerl(1) = 0.d0
      centerl(2) = 0.d0
      centerl(3) = 0.d0

      do icent = 0,natoms-1
        q =dbl_mb(icharge+icent)
        centerl(1) = centerl(1) + q*dbl_mb(icoord+3*icent  )
        centerl(2) = centerl(2) + q*dbl_mb(icoord+3*icent+1)
        centerl(3) = centerl(3) + q*dbl_mb(icoord+3*icent+2)
      enddo 
      if(.not. geom_nuc_charge(geom,total_nuc_charge)) call errquit
     &  ('dft_getspm: geom_nuc_charge failed',20, GEOM_ERR)

      centerl(1) = centerl(1)/total_nuc_charge
      centerl(2) = centerl(2)/total_nuc_charge
      centerl(3) = centerl(3)/total_nuc_charge

      if(.not.ma_pop_stack(ltags))
     &  call errquit('dft_getspm: cannot pop stack',0, MA_ERR)

      if(lmax.gt.0) then
        if(.not.Ma_Push_Get(MT_Dbl,3*natoms,'dipole',ldipole,
     &    idipole))
     &    call errquit('dft_getspm: cannot allocate dipole',0, MA_ERR)
      endif

      max_at_bf = 0
      do iatom = 1, natoms
        if(.not. bas_ce2bfr(basis, iatom, ilo, ihi))
     $    call errquit('quadvxc0: bas_ce2bfr failed', iatom, BASIS_ERR)
        max_at_bf = max(max_at_bf, ihi-ilo+1)
      enddo
      max_at_bf2=max_at_bf*max_at_bf

      if(.not.Ma_Push_Get(MT_Dbl,max_at_bf2,'local DM',lPmat,iPmat))
     &  call errquit('dft_getspm: cannot allocate local DM',0, MA_ERR)

      call dft_genspm(lmax,iga_dens,basis,
     &  natoms,centerl,dbl_MB(iPmat),max_at_bf,
     &  DBL_MB(icharge),DBL_MB(idipole),dbl_mb(icoord))
c
c     move the coordinates form coord(3,nat) to x y z arrays
c
      if(.not.Ma_Push_Get(MT_Dbl,natoms,'x coord',lxcrd,ixcrd))
     &  call errquit('dft_getspm: cannot allocate x coord',0, MA_ERR)
      if(.not.Ma_Push_Get(MT_Dbl,natoms,'y coord',lycrd,iycrd))
     &  call errquit('dft_getspm: cannot allocate y coord',0, MA_ERR)
      if(.not.Ma_Push_Get(MT_Dbl,natoms,'z coord',lzcrd,izcrd))
     &  call errquit('dft_getspm: cannot allocate z coord',0, MA_ERR)
      istart=0
      do i=0,natoms-1
        dbl_mb(ixcrd+i)=dbl_mb(icoord+istart)
        dbl_mb(iycrd+i)=dbl_mb(icoord+istart+1)
        dbl_mb(izcrd+i)=dbl_mb(icoord+istart+2)
        istart=istart+3
      enddo

      if(.not.ma_free_heap(lcoord))
     &  call errquit('dft_getspm: cannot free heap',0, MA_ERR)
c
c     re-center using geometrical center
c
      call dft_geocent(natoms,
     &  dbl_mb(ixcrd),dbl_mb(iycrd),dbl_mb(izcrd))
c
c     call to CMM code
c
C      call cmm_cmm(natoms,dbl_mb(ixcrd),dbl_mb(iycrd),dbl_mb(izcrd),
C     &  dbl_mb(icharge),finest)

      if(.not.ma_pop_stack(lzcrd))
     &  call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lycrd))
     &  call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lxcrd))
     &  call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
      if(.not.ma_pop_stack(lPmat))
     &  call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
      if(lmax.gt.0) then
        if(.not.ma_pop_stack(ldipole))
     &    call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
      endif
      if(.not.ma_pop_stack(lcharge))
     &    call errquit('dft_getspm: cannot pop stack',0, MA_ERR)
c
c       compute geometrical center for CMM
c

      return
      end
      subroutine dft_geocent(natoms,xcrd,ycrd,zcrd)
      implicit none
c-----------------------------------------------------------------------
c     Find the ranges of the molecular system in the x, y, and z directions
c-----------------------------------------------------------------------
#include "stdio.fh"
#include "global.fh"
#include "tcgmsg.fh"
      integer natoms
      double precision xcrd(natoms),ycrd(natoms),zcrd(natoms)
c     
      double precision xmin,xmax
      double precision ymin,ymax
      double precision zmin,zmax
      double precision xorigin,yorigin,zorigin
      double precision xrange,yrange,zrange
      double precision two
      parameter (two=2.d0)
      
      integer i

      xmin = 1.e+30
      ymin = 1.e+30
      zmin = 1.e+30
      xmax = -1.e+30
      ymax = -1.e+30
      zmax = -1.e+30

      do i = 1 , natoms

        xmin = min(xmin,xcrd(i))
        ymin = min(ymin,ycrd(i))
        zmin = min(zmin,zcrd(i))
        xmax = max(xmax,xcrd(i))
        ymax = max(ymax,ycrd(i))
        zmax = max(zmax,zcrd(i))

      end do

      xrange = dabs(xmin-xmax)
      yrange = dabs(ymin-ymax)
      zrange = dabs(zmin-zmax)

c-----------------------------------------------------------------------
c     Compute the geometrical center of molecular system
c-----------------------------------------------------------------------

      xorigin = (xmin + xmax)/TWO
      yorigin = (ymin + ymax)/TWO
      zorigin = (zmin + zmax)/TWO

c-----------------------------------------------------------------------
c     Refer the coordinates to the geometrical center of molecular system
c-----------------------------------------------------------------------

      do i = 1 , natoms

        xcrd(i) =  xcrd(i) - xorigin
        ycrd(i) =  ycrd(i) - yorigin
        zcrd(i) =  zcrd(i) - zorigin

      end do

c-----------------------------------------------------------------------
c     Output for checking
c-----------------------------------------------------------------------
      if(ga_nodeid().eq.0) then
        write(LuOut, '("minima",3(1x,f8.3))') xmin, ymin, zmin
        write(LuOut, '("maxima",3(1x,f8.3)/)') xmax, ymax, zmax

        write(LuOut, '("xwidth",f8.3)')  xrange
        write(LuOut, '("ywidth",f8.3)')  yrange
        write(LuOut, '("zwidth",f8.3/)') zrange 

        write(LuOut, '("origin",3(f8.3,1x)//)') 
     &       xorigin, yorigin, zorigin
      endif
      
      return
      end
