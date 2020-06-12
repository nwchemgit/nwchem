c     $Id$
      subroutine xc_xdm_init(rtdb,iixdm_v,iixdm_ml)
      implicit none

#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cdft.fh"

      integer rtdb  ! input
      integer natoms ! input
      integer iixdm_v, iixdm_ml ! output

      logical savedvml

c     xdm common
      integer nxdm
      integer ixdm_v, ixdm_ml 
      integer lxdm_v, lxdm_ml
      common /xdmd/ nxdm, ixdm_v, lxdm_v, ixdm_ml, lxdm_ml
      data ixdm_v / 0 /
      data ixdm_ml / 0 /

      if (.not. rtdb_get(rtdb, 'dft:xdm', mt_int, 1,
     &   lxdm)) lxdm = 0
      if(lxdm.eq.0) return

c     are the volumes and moments available from a previous iteration?
      if (.not. rtdb_get(rtdb,'dft:xdmsave', mt_log, 1, savedvml)) 
     &     savedvml = .false.

c     get memory for volume, alpha, ml
      nxdm = ntypes

      if (savedvml) then
c     check if v and ml have been allocated
         if(ixdm_v.eq.0.or.ixdm_ml.eq.0)  then
            if(ga_nodeid().eq.0) then
               write(6,*) ' xdm: v and ml not allocated '
               write(6,*) ' xdm: resetting savedvml '
            endif
            savedvml=.false.
            if (.not. rtdb_put(rtdb,'dft:xdmsave', mt_log, 1, savedvml)) 
     c       call errquit('xc_xdm_init: cant rtdb_put',0,0)
         else
         endif
      endif
      if (.not.savedvml) then
         if (.not.MA_alloc_get(mt_dbl, ntypes,
     &        'xdm_v', lxdm_v, ixdm_v)) then
            call errquit('xc_xdm_init: cant alloc xdm_v',0,0)
         endif
         if (.not.MA_alloc_get(mt_dbl, 3*ntypes,
     &        'xdm_ml', lxdm_ml, ixdm_ml)) then
            call errquit('xc_xdm_init: cant alloc xdm_ml',0,0)
         endif
      endif

      iixdm_v = ixdm_v
      iixdm_ml = ixdm_ml

      return
      end 

      subroutine xc_xdm(rtdb,g_dens,g_vxc,n,nexc,exdm,fxdm,v,ml,what)
      implicit none

#include "errquit.fh"
#include "rtdb.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "util_params.fh"

      integer rtdb, g_dens(2), g_vxc(4) ! input
      integer n, nexc ! input
      double precision exdm, fxdm(3,n) ! output
      double precision v(ntypes) ! inout, atomic volumes
      double precision ml(3,ntypes) ! inout, moments
      character*(*) what ! energy/forces flag

      integer i, j, me, nmu
      double precision rho_n, fxx
      double precision dum1, ftemp(3,n)
      logical setparam, savedvml
      integer izz, ixx, jxx

c     coefficents and vdw radii
      double precision c6(ntypes,ntypes), c8(ntypes,ntypes)
      double precision c10(ntypes,ntypes), rc(ntypes,ntypes)
      double precision rvdw(ntypes,ntypes)
      double precision a(ntypes)

c     geometry get
      logical geom_cent_get
      external geom_cent_get
      character*16 tag
      double precision x1(3), x2(3), q, r

c     input parameters
      double precision a1, a2
      logical onlyc, varc

c     xdm common
      integer nxdm
      integer ixdm_v, ixdm_ml, lxdm_v, lxdm_ml
      common /xdmd/ nxdm, ixdm_v, lxdm_v, ixdm_ml, lxdm_ml

c     free volumes and polarizabilities
C     DKH LSDA/UGBS Free Atomic Volumes
      double precision frevol(103)
      DATA (FREVOL(I),I=1,103)/
     +  9.194D0, 4.481D0,  91.957D0, 61.357D0, 49.813D0, 36.728D0,
     +  27.633D0, 23.517D0, 19.322D0, 15.950D0, 109.359D0, 103.064D0,
     +  120.419D0, 104.229D0, 86.782D0, 77.133D0, 66.372D0, 57.336D0,
     +  203.093D0, 212.202D0, 183.101D0, 162.278D0, 143.250D0,
     +  108.209D0, 123.098D0, 105.735D0, 92.944D0, 83.794D0, 75.750D0,
     +  81.177D0, 118.371D0, 116.334D0, 107.474D0, 103.221D0, 95.111D0,
     +  87.605D0, 248.772D0, 273.748D0, 249.211D0, 223.801D0,  
     +  175.809D0, 156.831D0, 160.042D0, 136.654D0, 127.754D0, 
     +  97.024D0, 112.778D0, 121.627D0, 167.906D0, 172.030D0, 
     +  165.500D0, 163.038D0, 153.972D0, 146.069D0, 341.992D0,
     +  385.767D0, 343.377D0, 350.338D0, 334.905D0, 322.164D0,  
     +  310.337D0,  299.537D0,  289.567D0,  216.147D0,  268.910D0,
     +  259.838D0,  251.293D0,  243.174D0,  235.453D0,  228.284D0,
     +  229.617D0,  209.971D0,  197.541D0,  183.236D0,  174.685D0,
     +  164.139D0,  150.441D0,  135.765D0,  125.297D0,  131.258D0,
     +  185.769D0,  195.671D0,  193.036D0,  189.142D0,  185.919D0,
     +  181.089D0,  357.787D0,  407.283D0,  383.053D0,  362.099D0,
     +  346.565D0,  332.462D0,  319.591D0,  308.095D0,  297.358D0,
     +  300.572D0,  275.792D0,  266.317D0,  257.429D0,  209.687D0,
     +  203.250D0,  230.248D0,  236.878D0/
      double precision frepol(103)
C     free atom polariz.: CRC Handbook of Chemistry and Physics, 88th Ed.
      DATA (FREPOL(I),I=1,102)/
     +  0.6668D0,  0.2051D0,  24.3300D0,  5.6000D0,  3.0300D0,
     +  1.7600D0,  1.1000D0,  0.8020D0,  0.5570D0,  0.3956D0,
     +  24.1100D0,  10.6000D0,  6.8000D0,  5.3800D0,  3.6300D0,
     +  2.9000D0,  2.1800D0,  1.6411D0,  43.4000D0,  22.8000D0,
     +  17.8000D0,  14.6000D0,  12.4000D0,  11.6000D0,  9.4000D0,
     +  8.4000D0,  7.5000D0,  6.8000D0,  6.2000D0,  5.7500D0,
     +  8.1200D0,  6.0700D0,  4.3100D0,  3.7700D0,  3.0500D0,
     +  2.4844D0,  47.3000D0,  27.6000D0,  22.7000D0,  17.9000D0,
     +  15.7000D0,  12.8000D0,  11.4000D0,  9.6000D0,  8.6000D0,
     +  4.8000D0,  7.2000D0,  7.3600D0,  10.2000D0,  7.7000D0,
     +  6.6000D0,  5.5000D0,  5.3500D0,  4.0440D0,  59.4200D0,
     +  39.7000D0,  31.1000D0,  29.6000D0,  28.2000D0,  31.4000D0,
     +  30.1000D0,  28.8000D0,  27.7000D0,  23.5000D0,  25.5000D0,
     +  24.5000D0,  23.6000D0,  22.7000D0,  21.8000D0,  21.0000D0,
     +  21.9000D0,  16.2000D0,  13.1000D0,  11.1000D0,  9.7000D0,
     +  8.5000D0,  7.6000D0,  6.5000D0,  5.8000D0,  5.0200D0,
     +  7.6000D0,  6.8000D0,  7.4000D0,  6.8000D0,  6.0000D0,
     +  5.3000D0,  48.6000D0,  38.3000D0,  32.1000D0,  32.1000D0,
     +  25.4000D0,  24.9000D0,  24.8000D0,  24.5000D0,  23.3000D0,
     +  23.0000D0,  22.7000D0,  20.5000D0,  19.7000D0,  23.8000D0,
     +  18.2000D0,  17.5000D0/
c
c     atomic densities
      integer ngau, k, ii
      double precision pi, vfree, rmid, qq, rr, h, rwei, rhofree
      parameter (ngau = 251)
c
      pi=acos(-1d0)
c     
      
c     read input parameters
      setparam = .not.rtdb_get(rtdb,'dft:xdm_a1',mt_dbl,1,a1)
      setparam = setparam .or.
     &     .not.rtdb_get(rtdb,'dft:xdm_a2',mt_dbl,1,a2)
      if (setparam) call xc_xdm_setdefaults(a1,a2)

      if (.not.rtdb_get(rtdb,'dft:xdm_onlyc',mt_log,1,onlyc)) 
     &     onlyc = .false.
      if (.not.rtdb_get(rtdb,'dft:xdm_varc',mt_log,1,varc)) 
     &     varc = .false.

      me = ga_nodeid()

c     redo the volume and moment calculation?
      if (.not. rtdb_get(rtdb,'dft:xdmsave', mt_log, 1, savedvml)) 
     &     savedvml = .false.
      if (varc .and. what.eq.'energy') then
         savedvml = .false.
      endif

c     header and convert a2 to au
      if (me.eq.0 .and. what.eq.'energy' .and..not.savedvml) then
         write (luout,*)
         write (luout,'("+ XDM: a1       = ",F12.5)') a1
         write (luout,'("       a2 (ang) = ",F12.5)') a2
         write (luout,'("       onlyc?   = ",L12)') onlyc
      endif
      a2 = a2 /cau2ang 

c     do the v and ml calculation
      if (.not.savedvml) then
c     initialize common values
         call dfill(ntypes, 0d0, v, 1)
         call dfill(3*ntypes, 0d0, ml, 1)

c     calculate volumes and moments using the dft grid
         rho_n = 0
         exdm = 0

         call grid_quadv0_gen(rtdb,g_dens,g_vxc,nexc,rho_n,exdm,1,6,
     &        dum1,.false.,.false.)
         call ga_dgop(msg_xdm1,v,ntypes,'+')
         call ga_dgop(msg_xdm3,ml,3*ntypes,'+')
         call ga_sync()

         do i = 1, ntypes
c     skip ghost atoms
            izz = znuc_atom_type(i)
            if (izz.eq.0) cycle

c     Only half of it if spin-polarized calculation
            if (ipol .eq. 2) then
               v(i) = v(i) / 2d0
               ml(1,i) = ml(1,i) / 2d0
               ml(2,i) = ml(2,i) / 2d0
               ml(3,i) = ml(3,i) / 2d0
            endif

c     Divide by the multiplicity
            nmu = 0
            do j = 1, n
               if (iatype(j).eq.i) nmu = nmu + 1
            enddo
            v(i) = v(i) / nmu
            ml(1,i) = ml(1,i) / nmu
            ml(2,i) = ml(2,i) / nmu
            ml(3,i) = ml(3,i) / nmu
         enddo

c     Calculate volumes and moments only once, then save for future use.
         if (.not. rtdb_put(rtdb,'dft:xdmsave', mt_log, 1, .true.))
     &        call errquit('xc_xdm: rtdb_put dft:xdmsave failed', 0,
     &        RTDB_ERR)
      endif
c     fill polarizabilities vector 
      do i = 1, ntypes
c        skip ghost atoms
         izz = znuc_atom_type(i)
         if (izz.eq.0) cycle
         a(i) = frepol(izz)/0.148184D0 * min(v(i)/frevol(izz),1d0)
      enddo

      if (me .eq. 0 .and..not.savedvml) then
         write (luout,*)
         write (luout,'("+ XDM: volume and moments")')
         write (luout,'("  All results in atomic units.")')
         write (luout,'("  i V alpha M1 M2 M3")')
         do i = 1, n
            ixx = iatype(i)
            izz = znuc_atom_type(ixx)
            if (izz.eq.0) cycle
            write (luout,'(I3,5(X,F18.10))') i, v(ixx), a(ixx),
     &           ml(1,ixx), ml(2,ixx), ml(3,ixx)
         end do
         write (luout,*)
      endif

c     calculate interaction coefficients
      if (me .eq. 0 .and. what.eq.'energy' .and..not.savedvml) then
         write (luout,'("+ XDM: dispersion coefficients")')
         write (luout,'(" All results in atomic units.")')
         write (luout,'("  i j C6 C8 C10 Rc Rvdw")')
      endif
      do i = 1, ntypes
         izz = znuc_atom_type(i)
         if (izz.eq.0) cycle

         do j = 1, i
            izz = znuc_atom_type(j)
            if (izz.eq.0) cycle
            c6(i,j) = a(i)*a(j)*ml(1,i)*ml(1,j) / 
     &           (ml(1,i)*a(j) + ml(1,j)*a(i))
            c6(j,i) = c6(i,j)
            c8(i,j) = 3d0/2d0 * (a(i)*a(j)*(ml(1,i)*
     &           ml(2,j)+ml(2,i)*ml(1,j))) /
     &           (ml(1,i)*a(j)+ml(1,j)*a(i))
            c8(j,i) = c8(i,j)
            c10(i,j) = 2 * a(i)*a(j) * (ml(1,i)*
     &           ml(3,j) + ml(3,i)*ml(1,j)) /
     &           (ml(1,i)*a(j) + ml(1,j)*a(i)) + 
     &           21d0/5d0 * a(i)*a(j)*
     &           ml(2,i)*ml(2,j) / (a(j)*ml(1,i)+
     &           a(i)*ml(1,j))
            c10(j,i) = c10(i,j)
            rc(i,j) = (dsqrt(c8(i,j)/c6(i,j)) + 
     &           dsqrt(c10(i,j)/c8(i,j)) + 
     &           (c10(i,j)/c6(i,j))**(0.25d0)) / 3
            rc(j,i) = rc(i,j)
            rvdw(i,j) = a1 * rc(i,j) + a2
            rvdw(j,i) = rvdw(i,j)
         enddo
      end do

      if (what.eq.'energy'.and..not.savedvml) then
         do i = 1, n
            do j = 1, i
               if (me .eq. 0) then
                  ixx = iatype(i)
                  izz = znuc_atom_type(ixx)
                  if (izz.eq.0) cycle
                  jxx = iatype(j)
                  izz = znuc_atom_type(jxx)
                  if (izz.eq.0) cycle
                  write (luout,'(I3,X,I3,1p,5(X,E18.10))') i, j, c6(ixx
     &                 ,jxx),c8(ixx,jxx), c10(ixx,jxx), rc(ixx,jxx)
     &                 ,rvdw(ixx,jxx)
               endif
            enddo
         enddo
      endif

      if (onlyc) return

c     Sum the dispersion energy
      exdm = 0d0
      do i = 1, n
         do j = 1,3
            ftemp(j,i) = 0d0
         enddo
      enddo

      do i = 1, n
         ixx = iatype(i)
         izz = znuc_atom_type(ixx)
         if (izz.eq.0) cycle
         if (.not.geom_cent_get(geom,i,tag,x1,q)) 
     &        call errquit('xc_xdm: geom_cent_get failed',geom,GEOM_ERR)

         do j = 1, i-1
            jxx = iatype(j)
            izz = znuc_atom_type(jxx)
            if (izz.eq.0) cycle
            if (.not.geom_cent_get(geom,j,tag,x2,q)) 
     &           call errquit('xc_xdm: geom_cent_get failed',geom
     &           ,GEOM_ERR)

            r = dsqrt((x1(1)-x2(1))**2+(x1(2)-x2(2))**2+
     &           (x1(3)-x2(3))**2)
            if (what.eq.'energy') then
               exdm = exdm - c6(ixx,jxx) / (rvdw(ixx,jxx)**6 + r**6)
     -           - c8(ixx,jxx) / (rvdw(ixx,jxx)**8 + r**8) 
     -           - c10(ixx,jxx) / (rvdw(ixx,jxx)**10 + r**10) 
            else
               fxx = -(
     +    6 * c6(ixx,jxx) * r**4 / (rvdw(ixx,jxx)**6 + r**6)**2+
     +    8 * c8(ixx,jxx) * r**6 / (rvdw(ixx,jxx)**8 + r**8)**2+
     +  10 * c10(ixx,jxx) * r**8 / (rvdw(ixx,jxx)**10 + r**10)**2)

               do k = 1, 3
                  ftemp(k,i) = ftemp(k,i) + (x2(k) - x1(k)) * fxx
                  ftemp(k,j) = ftemp(k,j) - (x2(k) - x1(k)) * fxx
               enddo
            endif
         enddo
      enddo
      if(what.eq.'forces') then
c     bit needed for forces
         do i = 1, n
            do k = 1, 3
               fxdm(k,i) = fxdm(k,i) + ftemp(k,i)
            end do
         end do
      endif

      if (me .eq. 0) then
         if (what.eq.'energy') then
            write (luout,
     &           '("+ XDM dispersion energy (Hy) = ",1p,E22.12)')
     &           exdm
         else
            do i = 1, n
               write (luout,
     &              '("  FXDM(",I2.2,")=",1p,3(E20.12,X),"a.u.")') 
     &              i, ftemp(1:3,i)
            end do
         endif
         write (luout,*)
      endif
      call ga_sync()
      return
      end

      subroutine xc_eval_xdm(rho,delrho,lap,nq,
     &     qxyz,qwght,ttau,
     &     natoms,xyz,v,ml)
      implicit none

#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "cdft.fh"
#include "mafdecls.fh"

      integer nq  ! input, no. of quadrature weights
      double precision rho(nq,ipol*(ipol+1)/2) ! input, electron density
      double precision delrho(nq,3,ipol) ! input, gradient of the electron density
      double precision lap(nq,ipol*(ipol+1)/2) ! input, laplacian
      double precision qxyz(3,nq) ! input, quadrature node positions
      double precision qwght(nq)  ! input, quadrature weights
      double precision ttau(nq,ipol) ! input, ked
      integer natoms ! input, no. of atoms
      double precision xyz(3,natoms) ! input, atomic positions
      double precision v(ntypes) ! inout, atomic volumes
      double precision a(ntypes) ! inout, atomic polarizabilities
      double precision ml(3,ntypes) ! inout, moments

      integer iz, it

c     parameters
      double precision pi
c      
      integer i, j, k, iat, isp
      double precision wei, rhoat, rhoi, r, x, y, z, vsum
      double precision rhot, rhos, grho, taus, laps
      double precision ds, qs, rhs, xroot, xshift, xold
      double precision expx, gx, fx, ffx, db, db2, db3
      double precision r2, r3
      double precision eps,thresh,mtwo3rds
      parameter(eps=1d-40,thresh=1d-12,mtwo3rds=-2d0/3d0)

c     atomic densities
      double precision rhop(nq), b(nq)
      double precision xdm_rhop
      external xdm_rhop
c
      pi=acos(-1d0)

      do i = 1, nq
c        calculate promolecular density on the grid
         rhop(i) = 0d0
         do j = 1, natoms
            x = qxyz(1,i) - xyz(1,j)
            y = qxyz(2,i) - xyz(2,j)
            z = qxyz(3,i) - xyz(3,j)
            r = dsqrt(x*x + y*y + z*z)

            iz = znuc_atom_type(iatype(j))
c           iz(j) == 0 is a ghost atom (bsse)
            if (iz.ge.1 .or. iz.le.94) then
               rhop(i) = rhop(i) + xdm_rhop(iz,r)
            elseif (iz.lt.0 .or. iz.gt.94) then
               call errquit('xc_xdm: wrong atomic number',j,iz)
            endif
         enddo
      end do

      do isp = 1, ipol
         do i = 1, nq
c     calculate dipole moment
            if (ipol .eq. 1) then
               rhot = max(rho(i,1),eps)
               rhos = max(rhot / 2d0,eps/2d0)
               grho = dsqrt(delrho(i,1,1)**2 + delrho(i,2,1)**2 + 
     +              delrho(i,3,1)**2) / 2d0
               taus = ttau(i,1)
               laps = lap(i,1) / 2d0
            else
               rhot = max(rho(i,2),eps)
               if (isp .eq. 1) then
                  rhos = rho(i,2) - rho(i,3)
                  grho = dsqrt((delrho(i,1,1)-delrho(i,1,2))**2 + 
     +                         (delrho(i,2,1)-delrho(i,2,2))**2 + 
     +                         (delrho(i,3,1)-delrho(i,3,2))**2)
                  taus = (ttau(i,1)-ttau(i,2))*2d0
                  laps = (lap(i,1)-lap(i,2))
               else
                  rhos = rho(i,3)
                  grho = dsqrt(delrho(i,1,2)**2 + delrho(i,2,2)**2 + 
     +                 delrho(i,3,2)**2)
                  taus = ttau(i,2)*2d0
                  laps = lap(i,2)
               endif
               rhos = max(rhos,0.5*eps)
            endif

            if(abs(rhos).gt.eps) then
               ds = taus - 0.25d0 * grho**2 / rhos
               qs = 1d0/6d0 * (laps - 2d0 * ds)
               rhs = 2d0/3d0*pi**(2d0/3d0)*(rhos)**(5d0/3d0)/qs
            else
               rhs= 0d0
            endif
            if (rhs > 0d0) then
               xroot = 3d0
               xshift = 1d0
               do while ((xroot*exp(-2d0*xroot/3d0))/(xroot-2d0)<rhs)
                  xshift = xshift * 0.1d0
                  xroot = 2d0 + xshift
               end do
            else
               xroot = 1d0
               xshift = 1d0
               do while ((xroot*exp(-2d0*xroot/3d0))/(xroot-2d0)>rhs)
                  xshift = xshift * 0.1d0
                  xroot = 2d0 - xshift
               end do
            end if
            if(abs(xroot).gt.eps) then
            xold = xroot + 1d0
            do while (abs(xroot - xold) > thresh)
               xold = xroot
               expx = exp(-2d0 * xroot / 3d0)
               gx = (xroot * expx) / (xroot - 2d0)
               fx = gx - rhs
               ffx = gx * (1d0 / xroot - 2d0/3d0 - 1d0 / (xroot - 2d0))
               xroot = xroot - fx / ffx
            end do
            endif
            b(i) = xroot * (exp(-xroot) / (8d0*pi*rhos))**(1d0/3d0)
         enddo      

c     calculate contribution to atomic volumes, alphas and moments
c     skip ghost atoms
         do iat = 1, natoms
            it = iatype(iat)
            iz = znuc_atom_type(it)
            if (iz.eq.0) cycle

            vsum = 0d0
            do i = 1, nq
               if (ipol .eq. 1) then
                  rhot = rho(i,1)
               else
                  rhot = rho(i,2)
               endif

c     atomic density of iat at the grid node
               x = qxyz(1,i) - xyz(1,iat)
               y = qxyz(2,i) - xyz(2,iat)
               z = qxyz(3,i) - xyz(3,iat)
               r = dsqrt(x*x + y*y + z*z)
               rhoi = xdm_rhop(iz,r)
               wei = rhoi / max(rhop(i),eps) * qwght(i)

c     contribution to atomic volume
               vsum = vsum + r**3 * rhot * wei

c     contribution to atomic moments
               db = max(r-b(i),0d0)
               db2 = db * db
               db3 = db2 * db
               r2 = r * r
               r3 = r2 * r
               ml(1,it) = ml(1,it) + 
     +              wei * rhot * (r - db)**2
               ml(2,it) = ml(2,it) + 
     +              wei * rhot * (r2 - db2)**2
               ml(3,it) = ml(3,it) + 
     +              wei * rhot * (r3 - db3)**2
            enddo
            v(it) = v(it) + vsum
         enddo
      enddo
      end

      subroutine xc_xdm_cleanup(rtdb)

#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c     xdm common
      logical savedvml
      integer nxdm, lxdm
      integer ixdm_v, ixdm_ml, lxdm_v, lxdm_ml
      common /xdmd/ nxdm, ixdm_v, lxdm_v, ixdm_ml, lxdm_ml

      if (.not. rtdb_get(rtdb, 'dft:xdm', mt_int, 1,
     &   lxdm)) lxdm = 0
      if(lxdm.eq.0) return


      if (.not. rtdb_get(rtdb,'dft:xdmsave', mt_log, 1, savedvml)) 
     &     savedvml = .false.
      if (.not.savedvml) then
      if (.not.MA_free_heap(lxdm_ml))
     &     call errquit('xc_xdm_cleanup: cannot clean ml',1,MA_ERR)
      if (.not.MA_free_heap(lxdm_v))
     &     call errquit('xc_xdm_cleanup: cannot clean v',1,MA_ERR)
      endif
      
      end

      subroutine xc_xdm_setdefaults(a1,a2)
#include "cdft.fh"   
c X          C        a1          a2
c becke86b   pbe96    0.82        1.16

      double precision a1, a2 ! output

      if (abs(xfac(55)).gt.1d-2 .and. abs(cfac(12)).gt.1d-2) then
         write (*,*) "WARNING: B86b!"
         a1 = 0.82d0
         a2 = 1.16d0
      else
         write(6,*)'WARNING:'
         write(6,*)'untested DF in XDM!!'
         a1 = 0.82d0
         a2 = 1.16d0
      endif

      end 

      function xdm_rhop(iz,r)

      integer iz ! input atomic number
      double precision r ! input distance to atom
      double precision xdm_rhop ! output promolecular density

c LDA atomic densities
c rho(r) = sum_i C_i * exp(-r/zeta_i)
      double precision rc1(94), rc2(94), rc3(94), rc4(94),
     + rc5(94), rc6(94), rc7(94), zeta1(94), zeta2(94),
     + zeta3(94), zeta4(94), zeta5(94), zeta6(94), zeta7(94)
      DATA (rc1(I),I=1,94)/
     + 2.38247581d-01, 2.65077191d+00, 1.35064937d+01, 3.52085244d+01,
     + 7.05975356d+01, 1.25179296d+02, 2.02122544d+02, 3.05256892d+02,
     + 4.38665900d+02, 6.07145508d+02, 8.39844590d+02, 1.10364452d+03,
     + 1.42251250d+03, 1.79552242d+03, 2.23091423d+03, 2.73461374d+03,
     + 3.31278653d+03, 3.97189405d+03, 4.70704295d+03, 5.55728188d+03,
     + 6.51013713d+03, 7.57599544d+03, 8.76508298d+03, 1.00881454d+04,
     + 1.15553877d+04, 1.31780299d+04, 1.49681265d+04, 1.69392910d+04,
     + 1.91364487d+04, 2.14808938d+04, 2.42414395d+04, 2.69764334d+04,
     + 3.00224115d+04, 3.33619034d+04, 3.70058615d+04, 4.09716964d+04,
     + 4.52399350d+04, 1.44691990d+04, 1.75811295d+04, 1.98707444d+04,
     + 2.74722280d+04, 3.04390736d+04,-2.32888935d+00, 3.87816243d+04,
     + 4.37829818d+04, 4.98625265d+04, 5.57415421d+04, 3.56643213d+03,
     + 7.17222714d+04,-1.14919886d+01, 8.97267971d+04, 1.00937307d+05,
     + 1.13620429d+05, 1.27857026d+05, 1.40502758d+05, 1.57775964d+05,
     + 1.78064645d+05, 2.00731613d+05, 2.26118301d+05, 2.54550489d+05,
     + 2.86238347d+05, 3.21681109d+05, 3.61267800d+05, 4.05477009d+05,
     + 4.54870353d+05, 5.09923518d+05, 5.71460462d+05, 6.40210146d+05,
     + 7.16876859d+05, 8.02519272d+05, 8.94824315d+05, 1.00284493d+06,
     + 1.12164004d+06, 1.25799024d+06, 1.41707468d+06, 1.56657370d+06,
     + 1.75110216d+06, 1.95685098d+06, 2.18876560d+06, 2.45050333d+06,
     + 2.73584853d+06, 3.06874608d+06, 3.44077168d+06, 3.85778178d+06,
     + 4.69176541d+06, 5.33909792d+06, 6.13639622d+06, 6.98851575d+06,
     + 7.92638571d+06, 9.01616115d+06, 1.02412946d+07, 1.16535915d+07,
     + 1.32723280d+07, 1.51679680d+07/
      DATA (zeta1(I),I=1,94)/
     + 5.94389750d-01, 3.63728528d-01, 1.72139700d-01, 1.21899504d-01,
     + 1.02003256d-01, 8.49661099d-02, 7.32065807d-02, 6.45649696d-02,
     + 5.78507920d-02, 5.22970367d-02, 4.14514521d-02, 3.73253761d-02,
     + 3.46867228d-02, 3.16764946d-02, 2.92112390d-02, 2.71143558d-02,
     + 2.52900857d-02, 2.36835198d-02, 2.20789393d-02, 2.08540286d-02,
     + 1.97256033d-02, 1.87075897d-02, 1.77753759d-02, 1.69434569d-02,
     + 1.61518812d-02, 1.54331276d-02, 1.47661483d-02, 1.41466847d-02,
     + 1.36293029d-02, 1.30268084d-02, 1.27388256d-02, 1.20747615d-02,
     + 1.15587091d-02, 1.10931067d-02, 1.06586168d-02, 1.02470835d-02,
     + 9.84539854d-03, 6.81246963d-04, 8.60318829d-04, 7.16123706d-04,
     + 1.82083380d-03, 1.49186880d-03, 8.37924126d-07, 1.20208187d-03,
     + 1.09125302d-03, 1.05045931d-03, 9.15785226d-04, 4.45100103d-03,
     + 8.42267282d-04, 4.55784641d-06, 6.49666842d-04, 6.03861716d-04,
     + 5.65447958d-04, 5.33042747d-04, 3.92130359d-04, 3.61775974d-04,
     + 3.55549092d-04, 3.48815329d-04, 3.42580721d-04, 3.36626058d-04,
     + 3.30233653d-04, 3.24004067d-04, 3.17860234d-04, 3.12193286d-04,
     + 3.05756593d-04, 3.00104329d-04, 2.94369895d-04, 2.88669044d-04,
     + 2.83327379d-04, 2.78093140d-04, 3.59124379d-04, 3.01791198d-04,
     + 2.88675155d-04, 2.57064717d-04, 1.32937491d-04, 2.70674086d-04,
     + 2.64755263d-04, 2.60031593d-04, 2.52814234d-04, 2.44193852d-04,
     + 2.43300956d-04, 2.32084917d-04, 2.23676582d-04, 2.16475971d-04,
     + 9.83996587d-05, 8.86389763d-05, 7.29954688d-05, 6.85803806d-05,
     + 6.78189840d-05, 6.58140892d-05, 6.50938107d-05, 6.38443503d-05,
     + 6.26090037d-05, 6.04309983d-05/
      DATA (rc2(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 4.68434853d-02, 2.52134142d-01,
     + 2.62330732d-01, 5.80715241d-01, 1.09645788d+00, 1.86382178d+00,
     + 2.97428667d+00, 4.62608071d+00, 2.21437182d+01, 3.56786858d+01,
     + 4.36941392d+01, 6.48665951d+01, 9.00636335d+01, 1.20264312d+02,
     + 1.56401427d+02, 1.99211797d+02, 2.67125067d+02, 3.21464903d+02,
     + 3.85771242d+02, 4.57510139d+02, 5.37648439d+02, 6.20149697d+02,
     + 7.18923103d+02, 8.22091519d+02, 9.34290727d+02, 1.05514023d+03,
     + 1.14249171d+03, 1.32583496d+03, 1.28645677d+03, 1.60012205d+03,
     + 1.84950772d+03, 2.09753208d+03, 2.35990684d+03, 2.64509680d+03,
     + 2.97636124d+03, 4.49374019d+04, 4.79280798d+04, 5.31722423d+04,
     + 5.00658603d+04, 5.59792506d+04, 6.72129157d+04, 6.73325347d+04,
     + 7.36240653d+04, 7.95039519d+04, 8.76089515d+04, 1.02661441d+05,
     + 1.02056673d+05, 1.20930700d+05, 1.23237713d+05, 1.33810482d+05,
     + 1.44995277d+05, 1.56862063d+05, 1.77417216d+05, 1.90991139d+05,
     + 2.05075073d+05, 2.20185428d+05, 2.36331156d+05, 2.53473962d+05,
     + 2.71763871d+05, 2.91260169d+05, 3.12044059d+05, 3.34190817d+05,
     + 3.57837847d+05, 3.82982425d+05, 4.09826224d+05, 4.38465778d+05,
     + 4.68930154d+05, 5.01396788d+05, 4.85613275d+05, 5.51030851d+05,
     + 5.93454095d+05, 6.54789149d+05, 5.99621174d+05, 7.21896735d+05,
     + 7.70758213d+05, 8.21242604d+05, 8.79338995d+05, 9.44887717d+05,
     + 9.99564233d+05, 1.08189801d+06, 1.16519392d+06, 1.25282932d+06,
     + 1.57499234d+06, 1.79801905d+06, 2.16513590d+06, 2.44913677d+06,
     + 2.69497434d+06, 2.99908623d+06, 3.30131969d+06, 3.65338787d+06,
     + 4.04499005d+06, 4.52885470d+06/
      DATA (zeta2(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.01441639d+00, 7.08517236d-01,
     + 7.30003765d-01, 6.15383618d-01, 5.33576000d-01, 4.72577360d-01,
     + 4.24378943d-01, 3.83754227d-01, 2.38826415d-01, 1.99477996d-01,
     + 1.92829772d-01, 1.67225112d-01, 1.48741905d-01, 1.34281775d-01,
     + 1.22429085d-01, 1.12511023d-01, 9.82839214d-02, 9.20130613d-02,
     + 8.66687270d-02, 8.18541208d-02, 7.75184939d-02, 7.44554398d-02,
     + 7.02827730d-02, 6.71643705d-02, 6.43177023d-02, 6.17266464d-02,
     + 6.13921845d-02, 5.71521336d-02, 6.13478031d-02, 5.43903881d-02,
     + 5.05712825d-02, 4.76791348d-02, 4.52299778d-02, 4.30400573d-02,
     + 4.01557534d-02, 1.12290074d-02, 1.13416879d-02, 1.07424855d-02,
     + 1.19482805d-02, 1.14350394d-02, 1.00142006d-02, 1.06774861d-02,
     + 1.03201321d-02, 1.00700543d-02, 9.64166067d-03, 8.51583692d-03,
     + 9.15978091d-03, 7.98513470d-03, 8.34889598d-03, 8.04927525d-03,
     + 7.77194190d-03, 7.51012133d-03, 6.60822060d-03, 6.22807064d-03,
     + 6.07149684d-03, 5.91742394d-03, 5.77498801d-03, 5.63905896d-03,
     + 5.50038331d-03, 5.36796683d-03, 5.23954025d-03, 5.11839482d-03,
     + 4.99401318d-03, 4.87751356d-03, 4.76426996d-03, 4.65444537d-03,
     + 4.54930660d-03, 4.44815005d-03, 5.12909519d-03, 4.61664077d-03,
     + 4.44800980d-03, 4.05938174d-03, 1.52242450d-03, 4.15624532d-03,
     + 4.06185185d-03, 3.98192150d-03, 3.87525947d-03, 3.75252197d-03,
     + 3.71912608d-03, 3.56493520d-03, 3.44266883d-03, 3.33335862d-03,
     + 1.04842450d-03, 8.93405393d-04, 6.57894932d-04, 6.00392630d-04,
     + 5.94698138d-04, 5.71622306d-04, 5.66927029d-04, 5.54520209d-04,
     + 5.42218623d-04, 5.16230841d-04/
      DATA (rc3(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 6.20076037d-02, 2.76051005d-01,
     + 1.61014857d-01, 3.93461956d-01, 7.64466624d-01, 1.30215679d+00,
     + 2.04180763d+00, 3.00984485d+00, 6.98610402d+00, 1.00834127d+01,
     + 1.23076665d+01, 1.50188788d+01, 1.82395980d+01, 1.97075623d+01,
     + 2.64041322d+01, 3.15174107d+01, 3.73934453d+01, 4.40633918d+01,
     + 4.26128540d+01, 6.00190607d+01, 4.45592617d+01, 7.46129815d+01,
     + 1.02455412d+02, 1.31470485d+02, 1.62925314d+02, 1.97827848d+02,
     + 2.65304373d+02, 1.57117900d+03, 1.38597704d+03, 1.71450052d+03,
     + 1.10431209d+03, 1.23423271d+03, 1.74411089d+03, 1.49332905d+03,
     + 1.64723180d+03, 1.75852626d+03, 2.01523534d+03, 2.69723853d+03,
     + 2.31628288d+03, 3.21014982d+03, 3.16820224d+03, 3.57493869d+03,
     + 4.00748342d+03, 4.48040398d+03, 8.18987220d+03, 1.07287229d+04,
     + 1.14521494d+04, 1.22608234d+04, 1.30255796d+04, 1.37766798d+04,
     + 1.46837652d+04, 1.55931040d+04, 1.65429433d+04, 1.74661109d+04,
     + 1.85712156d+04, 1.96379950d+04, 2.07455462d+04, 2.18891141d+04,
     + 2.30389584d+04, 2.41952252d+04, 1.13020962d+04, 1.75777161d+04,
     + 1.98740166d+04, 2.95731860d+04, 4.09107217d+05, 2.31029894d+04,
     + 2.43065901d+04, 2.51685750d+04, 2.69880676d+04, 2.96363456d+04,
     + 2.90462255d+04, 3.34752965d+04, 3.71679846d+04, 4.07226602d+04,
     + 7.55624974d+05, 8.77788402d+05, 1.10863240d+06, 1.22508590d+06,
     + 1.29264761d+06, 1.38835313d+06, 1.46321597d+06, 1.55554758d+06,
     + 1.65439402d+06, 1.79235741d+06/
      DATA (zeta3(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00835672d+00, 7.69492311d-01,
     + 8.95052085d-01, 7.49562876d-01, 6.50577105d-01, 5.78352607d-01,
     + 5.22541213d-01, 4.77987125d-01, 3.65799244d-01, 3.09336153d-01,
     + 2.97043149d-01, 2.87814545d-01, 2.79727989d-01, 2.91308753d-01,
     + 2.64243561d-01, 2.56546792d-01, 2.48992596d-01, 2.41688705d-01,
     + 2.55650963d-01, 2.27972392d-01, 2.57800251d-01, 2.17760648d-01,
     + 1.95357993d-01, 1.79039271d-01, 1.65938044d-01, 1.54814449d-01,
     + 1.34999942d-01, 6.51041574d-02, 7.21674364d-02, 6.50557950d-02,
     + 8.90037912d-02, 8.60216892d-02, 7.34617487d-02, 8.14402358d-02,
     + 7.89401490d-02, 7.80479431d-02, 7.36580206d-02, 6.49836533d-02,
     + 7.11605517d-02, 6.14493259d-02, 6.10610277d-02, 5.77205998d-02,
     + 5.47391222d-02, 5.19382050d-02, 3.30078991d-02, 2.52817223d-02,
     + 2.45802867d-02, 2.40980268d-02, 2.39041263d-02, 2.36226590d-02,
     + 2.32201481d-02, 2.28792384d-02, 2.25607988d-02, 2.23758534d-02,
     + 2.19870231d-02, 2.17465623d-02, 2.15260243d-02, 2.13292973d-02,
     + 2.11855810d-02, 2.10847684d-02, 4.09382448d-02, 3.03603604d-02,
     + 2.82918344d-02, 2.03824377d-02, 5.62977448d-03, 2.72440875d-02,
     + 2.68637757d-02, 2.67959682d-02, 2.60212616d-02, 2.48170625d-02,
     + 2.57591136d-02, 2.36991792d-02, 2.23967941d-02, 2.13654035d-02,
     + 4.55423037d-03, 4.25297993d-03, 3.72189295d-03, 3.53030502d-03,
     + 3.47317821d-03, 3.36862712d-03, 3.31841195d-03, 3.24626475d-03,
     + 3.17505028d-03, 3.05979289d-03/
      DATA (rc4(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 6.56020424d-02, 2.77840007d-01,
     + 3.64029769d-01, 4.23460753d-01, 4.63193641d-01, 2.89810882d-01,
     + 5.14790980d-01, 5.33655411d-01, 5.50170763d-01, 5.65407164d-01,
     + 2.84096529d-01, 5.93010550d-01, 1.19942620d-01, 3.60612540d-01,
     + 7.40570005d-01, 1.28553223d+00, 2.02073112d+00, 2.97615952d+00,
     + 8.39557269d+00, 6.02606011d+01, 4.87023540d+01, 6.80933520d+01,
     + 1.62265279d+01, 1.72400236d+01, 3.36917884d+01, 1.83723888d+01,
     + 1.98403203d+01, 1.85143190d+01, 2.48910274d+01, 4.01353587d+01,
     + 2.64509614d+01, 4.78655062d+01, 6.01556846d+01, 7.90020599d+01,
     + 1.01162437d+02, 1.28370348d+02, 7.11169654d+02, 1.64398006d+03,
     + 1.80851978d+03, 1.91100792d+03, 1.94620437d+03, 2.02596238d+03,
     + 2.12145362d+03, 2.21085640d+03, 2.29697962d+03, 2.34539613d+03,
     + 2.45263808d+03, 2.51486491d+03, 2.56861396d+03, 2.61136834d+03,
     + 2.63029593d+03, 2.62765250d+03, 1.95848638d+02, 6.52135512d+02,
     + 8.28259407d+02, 2.73088861d+03, 1.26198542d+04, 8.86503684d+02,
     + 9.16101499d+02, 9.11293836d+02, 9.93328324d+02, 1.14814275d+03,
     + 1.00157533d+03, 1.30190843d+03, 1.54552793d+03, 1.77816217d+03,
     + 2.04357126d+04, 2.35854277d+04, 3.42350892d+04, 4.03227664d+04,
     + 4.16645314d+04, 4.53355366d+04, 4.65504326d+04, 4.90202580d+04,
     + 5.16318071d+04, 5.74706407d+04/
      DATA (zeta4(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.09950352d+00, 8.70037206d-01,
     + 8.15822153d-01, 7.79528672d-01, 7.53236219d-01, 7.85765325d-01,
     + 7.14317589d-01, 6.98660493d-01, 6.84606242d-01, 6.71753041d-01,
     + 7.16027394d-01, 6.48884245d-01, 9.09499464d-01, 7.63611851d-01,
     + 6.70634196d-01, 6.03533737d-01, 5.51786971d-01, 5.10051860d-01,
     + 3.79712304d-01, 2.34159182d-01, 2.43754955d-01, 2.22148073d-01,
     + 3.31676536d-01, 3.30775898d-01, 2.71976660d-01, 3.39485392d-01,
     + 3.39067921d-01, 3.57675223d-01, 3.28498848d-01, 2.82667731d-01,
     + 3.28954437d-01, 2.74280814d-01, 2.50819436d-01, 2.27953575d-01,
     + 2.08767308d-01, 1.91777539d-01, 1.05361706d-01, 7.68835000d-02,
     + 7.41638307d-02, 7.29762271d-02, 7.29383574d-02, 7.20265634d-02,
     + 7.10019314d-02, 7.00798829d-02, 6.92361234d-02, 6.88348560d-02,
     + 6.78372320d-02, 6.73390280d-02, 6.69450435d-02, 6.66732797d-02,
     + 6.66395164d-02, 6.68229001d-02, 2.07796306d-01, 1.22939973d-01,
     + 1.09882989d-01, 6.61756802d-02, 3.89737289d-02, 1.16012046d-01,
     + 1.17207641d-01, 1.20249726d-01, 1.17166773d-01, 1.10793849d-01,
     + 1.19240694d-01, 1.06979518d-01, 9.93514855d-02, 9.33089232d-02,
     + 3.36075409d-02, 3.12070010d-02, 2.45901087d-02, 2.22122582d-02,
     + 2.21170983d-02, 2.11671631d-02, 2.11828525d-02, 2.07857633d-02,
     + 2.03904909d-02, 1.91739084d-02/
      DATA (rc5(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 7.41352221d-02, 3.29043651d-01, 4.03366961d-01, 6.16991180d-01,
     + 3.28706457d-01, 3.89407305d-01, 7.60824501d-01, 3.56433510d-01,
     + 3.21154553d-01, 2.12638921d-01, 2.68825717d-01, 6.74790071d-01,
     + 1.16080256d-01, 3.74840652d-01, 7.76154174d-01, 1.36358170d+00,
     + 2.15461358d+00, 3.17719238d+00, 1.79344549d+01, 4.54639997d+01,
     + 5.36279623d+01, 5.66894253d+01, 5.62578366d+01, 6.02089260d+01,
     + 6.47699752d+01, 6.97534168d+01, 7.51227593d+01, 8.08249928d+01,
     + 8.68693210d+01, 9.31064937d+01, 9.95858875d+01, 1.06259016d+02,
     + 1.12891772d+02, 1.19452938d+02, 4.86622657d-01, 3.75386966d+01,
     + 6.13420753d+01, 1.64863797d+02, 2.83371154d+02, 3.37760283d+01,
     + 2.55288121d+01, 1.66310473d+01, 1.76333146d+01, 2.35863127d+01,
     + 1.44663326d+01, 2.67219996d+01, 4.10151210d+01, 5.89740153d+01,
     + 4.56772818d+02, 5.85287609d+02, 1.29543017d+03, 1.78567564d+03,
     + 1.79552725d+03, 2.05105206d+03, 2.02552823d+03, 2.13209240d+03,
     + 2.24557320d+03, 2.69991263d+03/
      DATA (zeta5(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.11162276d+00, 8.92973194d-01, 8.48665299d-01, 7.81782861d-01,
     + 8.20458147d-01, 7.79621642d-01, 7.10637387d-01, 7.49615537d-01,
     + 7.43606586d-01, 7.03198459d-01, 7.34936372d-01, 6.65219859d-01,
     + 9.33059092d-01, 7.92978117d-01, 7.04257390d-01, 6.39255341d-01,
     + 5.88851087d-01, 5.47985182d-01, 3.64082401d-01, 2.77824918d-01,
     + 2.62876282d-01, 2.61606275d-01, 2.68464970d-01, 2.65686159d-01,
     + 2.62414808d-01, 2.58938476d-01, 2.55366523d-01, 2.49843912d-01,
     + 2.48193305d-01, 2.44727449d-01, 2.41346773d-01, 2.38072521d-01,
     + 2.34997746d-01, 2.32105807d-01, 8.03394712d-01, 2.85600858d-01,
     + 2.53758662d-01, 2.08195570d-01, 1.84383408d-01, 2.84437965d-01,
     + 3.08075946d-01, 3.63374585d-01, 3.64788118d-01, 3.33552690d-01,
     + 3.95115818d-01, 3.33602667d-01, 2.93606072d-01, 2.62490124d-01,
     + 1.63387243d-01, 1.50278726d-01, 1.11071589d-01, 9.65712781d-02,
     + 9.73247043d-02, 9.18921858d-02, 9.40757894d-02, 9.28564165d-02,
     + 9.15796780d-02, 8.47922939d-02/
      DATA (rc6(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 8.41198181d-02, 3.30444084d-01,
     + 4.91705562d-01, 4.94765259d-01, 3.69753458d-01, 3.69539339d-01,
     + 3.69853666d-01, 3.70622471d-01, 3.71721718d-01, 4.34972852d-01,
     + 3.75404514d-01, 3.77695171d-01, 3.80389931d-01, 3.83527706d-01,
     + 3.86898212d-01, 3.90604654d-01, 2.97047551d-03, 3.99986047d-01,
     + 7.20996101d-01, 1.06952842d+00, 1.37569576d+00, 1.25089939d+00,
     + 1.23420443d+00, 6.29685108d-01, 5.66212693d-01, 1.04389432d+00,
     + 7.78008312d-02, 2.99789533d-01, 6.88443909d-01, 1.27022912d+00,
     + 2.43994094d+00, 3.56059743d+00, 1.77231239d+01, 4.17467686d+01,
     + 3.83943193d+01, 5.43331103d+01, 3.96156778d+01, 3.92950589d+01,
     + 3.94476243d+01, 5.61524915d+01/
      DATA (zeta6(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.14652099d+00, 9.40419863d-01,
     + 8.71830139d-01, 8.63645308d-01, 8.99781351d-01, 8.92847773d-01,
     + 8.86263534d-01, 8.79891833d-01, 8.73693204d-01, 8.44849538d-01,
     + 8.61492212d-01, 8.55512154d-01, 8.49557770d-01, 8.43601596d-01,
     + 8.37682504d-01, 8.31763174d-01, 1.16006472d+00, 8.07528576d-01,
     + 7.31082249d-01, 6.83425161d-01, 6.52505241d-01, 6.40350000d-01,
     + 6.28116671d-01, 6.50985392d-01, 6.43444641d-01, 6.04768070d-01,
     + 9.74028988d-01, 8.20486767d-01, 7.26711714d-01, 6.60125651d-01,
     + 6.03396512d-01, 5.63782867d-01, 3.80991010d-01, 2.97211769d-01,
     + 3.02765896d-01, 2.71958506d-01, 3.05459861d-01, 3.09950554d-01,
     + 3.13551092d-01, 2.87960597d-01/
      DATA (rc7(I),I=1,94)/
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 0.00000000d+00, 0.00000000d+00,
     + 0.00000000d+00, 0.00000000d+00, 1.01580194d-01, 3.80390160d-01,
     + 4.13217274d-01, 7.39534177d-01, 4.50487346d-01, 4.21044941d-01,
     + 3.82783864d-01, 5.04264000d-01/
      DATA (zeta7(I),I=1,94)/
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.00000000d+00, 1.00000000d+00,
     + 1.00000000d+00, 1.00000000d+00, 1.10753948d+00, 9.23783778d-01,
     + 9.01059856d-01, 8.16067064d-01, 8.69278115d-01, 8.68260354d-01,
     + 8.71564352d-01, 8.34358766d-01/

      xdm_rhop = rc1(iz) * exp(-r/zeta1(iz))
      if (zeta2(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc2(iz) * exp(-r/zeta2(iz))
      if (zeta3(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc3(iz) * exp(-r/zeta3(iz))
      if (zeta4(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc4(iz) * exp(-r/zeta4(iz))
      if (zeta5(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc5(iz) * exp(-r/zeta5(iz))
      if (zeta6(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc6(iz) * exp(-r/zeta6(iz))
      if (zeta7(iz)<1d-12) return
      xdm_rhop = xdm_rhop + rc7(iz) * exp(-r/zeta7(iz))

      end
      integer function xc_xdm_lxdm()
      implicit none

#include "cdft.fh"
      xc_xdm_lxdm=lxdm
      return
      end
