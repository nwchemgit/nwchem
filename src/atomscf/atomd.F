      subroutine atomd(oprin,iwr,znps,ic,isymax,hatom,
     + pcap,qcap,fc,fo,s,u,t,h,dc,dos,dt,dold,ss,
     + c, copn, smin, qmin, transf, cc , nbb, element, tag,spec)
C$Id$
      implicit none
#include "errquit.fh"
#include "inp.fh"
c
      logical oprin,spec
      integer iwr, nbb, isymax
      integer ic(4,*)
      real *8  pcap(*), qcap(*), fc(*), fo(*), s(*), u(*), t(*)
      real *8  h(*), dc(*), dos(*), dt(*), dold(*), ss(*), hatom(*)
      real *8  c(*),copn(*),smin(nbb,*),qmin(nbb,*),transf(*),cc(*)
      real *8 znps
      character*(*) element, tag ! Tag and element of atom
      integer tagl, elel        !  Length of strings minus trailing blanks
c
c.......................................................................
c     atomic r h f code for gtos. uses roothaan double diagonalization.
c.......................................................................
c
c      routine datoms .. tramad are a totally separate unit
c      they use commom cguess to communicate
c      oprin,iwr are new parameters to control printing
c      hatom are the pseudo-corrections in the contracted basis from gam
c            they replace atoms own integrals
c.......................................................................
c
#include "cguess.fh"
c
      integer i, j, l, m
      integer nbcdim, jcount, k, iant, ksum
      integer noff, nbci, ll, mm, ndim, nsqt, ns, maxitr, lm
      integer nstep1, nstep2, knteps, nbc1, nbc2, nbc3, icount
      real *8 znsave, tlarge
c
c     nsht      = total number of shells
c     n1(i)     = nbas(i) * (nbas(i) + 1 ) / 2
c     nbc(i)    = number of cont. orbitals in symmetry i
c     cont(i)   = contraction coeff. assosiated with primitive no. i
c     nstrt(i)  = number for first primitive in cont. no. i
c     nbct      = total number of cont. basis functions.
c.......................................................................
      if (dabs(zn).ge.1.d-8) then
c.......................................................................
         call atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps,spec)
c
c     distribute electrons according to aufbau
c.......................................................................

crz
c     fix for use of pseudopotentials
         znsave = zn
         zn = znps
c.......................................................................
c
c     move basis set information from transfer variables to working
c     variables. ** not necessary in present version
c.......................................................................
         nbct = 0
         ndim = 0
         nsqt = 0
         nbcdim = 0
         jcount = 0
         nsht = 0
         k = 0
         do i = 1 , nsym
            n1(i) = nbas(i)*(nbas(i)+1)/2
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
            ndim = ndim + n1(i)
            nsht = nsht + ncsh(i) + nosh(i)
            nsqt = nsqt + nbas(i)**2
            nbct = nbct + nbc(i)
            do j = 1 , nbas(i)
               jcount = jcount + 1
            enddo
            do j = 1 , n1(i)
               k = k + 1
               dold(k) = 0.0d0
            enddo
         enddo
         ns = 1
         nstrt(ns) = 1
         do l = 1 , nsym
            iant = nbc(l)
            ksum = 0
            do i = 1 , iant
               nstrt(ns+i) = nstrt(ns+i-1) + ic(l,i)
               ksum = ksum + ic(l,i)
            enddo
            if (ksum.ne.nbas(l)) write (iwr,6020) l , ksum , nbas(l)
            ns = ns + iant
         enddo
         do i = 1 , nsqt
            c(i) = 0.d0
            copn(i) = 0.0d0
            cc(i) = 0.d0
         enddo
         do i = 1 , ndim
            dold(i) = 0.0d0
         enddo
cjvl  few extra checks
         do i=1,nsym
            if (nbc(i).lt.ncsh(i)+nosh(i)) then
               oprin = .true.
               ndim = 0
            end if
         enddo
cjvl      
         tagl = inp_strlen(tag)
         elel = inp_strlen(element)
         if (oprin) write (iwr,6030) zn, tag(1:tagl), element(1:elel)
         if (oprin) write (iwr,6040) (nbas(i),i=1,nsym)
         if (oprin) write (iwr,6050) (nbc(i),i=1,nsym)
         if (oprin) write (iwr,6060) (ncsh(i),i=1,nsym)
         if (oprin) write (iwr,6070) (nosh(i),i=1,nsym)
         if (oprin) write (iwr,6080) (nccup(i),i=1,nsym)
cjvl
         if (ndim.eq.0)
     &       call errquit('atomscf:atomd: error in atom scf',0,
     &       UNKNOWN_ERR)
cjvl
         maxitr = 100
c
c     debug output of exponents/coeffs/...
c
         if (oprin) call printbasis
c
c..
c..     calculate 1-electron ints
c..
         call oeigd(fc,s,u,t,h)
c.......................................................................
c
c     copy overlap matrix to ss
c.......................................................................
         do i = 1 , ndim
            ss(i) = s(i)
            fc(i) = s(i)
            fo(i)=0.d0
         enddo
c.......................................................................
c
c     now transform ss to contracted basis, then set up transformation
c     matrix to o.n. contracted basis.
c.......................................................................
         call trafsd(nsym,nbas,ndim,ss,nbc,cont,nstrt,dc)
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)
c...
         nstep1 = 1
         nstep2 = 1
         do i = 1 , nsym
crz      to surpress problems when there is one type of shell missing...
            if (nbc(i).ne.0) then
               call shalfd(fc(nstep1),transf(nstep2),nbc(i))
               call starcd(cc(nstep2),ss(nstep1),nbc(i),ncsh(i),
     +                     nosh(i))
               nstep1 = nstep1 + nbc(i)*(nbc(i)+1)/2
               nstep2 = nstep2 + nbc(i)**2
            end if
         enddo
         nitscf = 0
         nconv = 0
         damp = .30d0
         if (znsave.eq.30.0d0) then
            damp = .9d0
            maxitr = 200
c         print *,' zn is special ',damp,maxitr
         end if
 110     nitscf = nitscf + 1
c.......................................................................
c
c     transform vectors and set up matrices in primitive basis,
c     then transform fock matrices to contracted basis.
c.......................................................................
         call tracd(c,cc,nsqt,nsym, nbas, ncsh, nosh, nbc,nstrt,
     +                 cont)
         call densid(dt,dold,dos,nsym,nosh,ncsh,nccup,c,damp,nconv,nbas,
     +               nitscf,tlarge)
c
c... check for convergence on tlarge (max change of d-matrix)
c
         if (tlarge.le.1.0d-5) nconv = 1
         if (nitscf+20.ge.maxitr) then
            write (iwr,6010) nitscf , energ , tlarge
         end if
c
         call hamild(pcap,qcap,fc,fo,s,u,t,h,dos,dt,c,smin,qmin,nbb)
c...
c...    now add the h-pseudo-contributions (from xpsnld)
c...    the matrices are in the contracted basis here
c...
         do i = 1 , ndim
            fc(i) = fc(i) + hatom(i)  ! hatom is in the contracted basis
            fo(i) = fo(i) + hatom(i)
         enddo
         call trafsd(nsym,nbas,ndim,fc,nbc,cont,nstrt,dc)  ! convert fc to contracted basis
         call trafsd(nsym,nbas,ndim,fo,nbc,cont,nstrt,dc)  ! convert fo to contracted basis
c
c.......................................................................
c
c     do double diagonalization by symmetries:
c         1. transform block to o.n.basis (contracted).
c         2. store o.n. transformation matrix in vector matrix.
c         3. diagonalize.
c         4. order eigenvectors by eigenvalue.
c         5. if necessary, merge open and closed vectors.
c.......................................................................
         nstep1 = 1
         nstep2 = 1
         knteps = 0
         do i = 1 , nsym
            nbc1 = nbc(i)
            if(nbc1.ne.0) then
            nbc2 = nbc1**2
            nbc3 = (nbc2+nbc1)/2
            if (ncsh(i).ne.0) then
               call tramad(fc(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,cc(nstep2),1)
               call jacod(fc(nstep1),cc(nstep2),nbc1,n1(i),nbc2,1,nbc1,
     +                    dc,dt,nbc1)
            end if
            if (nosh(i).ne.0) then
               call tramad(fo(nstep1),transf(nstep2),dc,nbc3,nbc1,dt)
               call dcopy(nbc2,transf(nstep2),1,copn(nstep2),1)
               call jacod(fo(nstep1),copn(nstep2),nbc1,n1(i),nbc2,1,
     +                    nbc1,dc,dt,nbc1)
            end if
            icount = nstep1
            do j = 1 , nbc1
               dc(j) = fc(icount)
               dos(j) = fo(icount)
               icount = icount + 1 + j
            enddo
            call orderd(cc(nstep2),nbc1,nbc1,dt,dc)
            if (nosh(i).gt.0) then
               call orderd(copn(nstep2),nbc1,nbc1,dt,dos)
               call cmergd(cc(nstep2),copn(nstep2),ncsh(i),nbc1,nosh(i))
            end if
            nstep1 = nstep1 + nbc1*(nbc1+1)/2
            nstep2 = nstep2 + nbc2
            if (nconv.gt.0) then
               if (ncsh(i).gt.0) then
                  do j = 1 , ncsh(i)
                     knteps = knteps + 1
                     eps(knteps) = dc(j)
                  enddo
               end if
               if (nosh(i).gt.0) then
                  do j = 1 , nosh(i)
                     knteps = knteps + 1
                     eps(knteps) = dos(ncsh(i)+j)
                  enddo
               end if
            end if
            endif
         enddo
         if (nitscf.ge.maxitr) nconv = 1
         if (nconv.le.0) go to 110
         if (oprin) then
            call outpud(copn,cc,1,iwr)
         else
            call outpud(copn,cc,0,iwr)
         end if
         call densid_z(dt,dold,dos,nsym,nosh,ncsh,nccup,cc,damp,nconv,
     &       nbc,nitscf,tlarge,tag,eps)
*         call densid(dt,dold,dos,nsym,nosh,ncsh,nccup,cc,damp,nconv,nbc,
*     +               nitscf,tlarge)
         lm = 0
         do i = 1 , nsym
            noff = lm
            nbci = nbc(i)
            do l = 1 , nbci
               ll = noff + l*(l+1)/2
CDIR$  NEXTSCALAR
               do m = 1 , l
                  mm = noff + m*(m+1)/2
                  lm = lm + 1
                  dt(lm) = dt(lm)*sqrt(ss(ll)*ss(mm))
                  if (m.ne.l) dt(lm) = dt(lm)/2.0d0
               enddo
            enddo
         enddo
      else
c.......................................................................
c
c     special section for handling the case of floating functions
c     on centers with no charge.
c.......................................................................
         nbcdim = 0
         do i = 1 , isymax
            nbcdim = nbcdim + nbc(i)*(nbc(i)+1)/2
         enddo
         do i = 1 , nbcdim
            dt(i) = 0.0d0
         enddo
         nsym = isymax
         energ = 0.0d0
      end if
c
      return
 6010 format(' it.',i4,'  energy',d19.10,'  div.',d13.5)
 6020 format('-',' wrong contraction in symmetry   ',3i5)
 6030 format(/'      charge =',f10.1,1x,a,' (',a,')'//
     $     '      symmetry species',12x, 's',5x,'p',5x,'d',5x,'f')
 6040 format(6x,'number of basis functions =',4(i2,4x))
 6050 format(6x,'number of cont. functions =',4(i2,4x))
 6060 format(6x,'number of closed shells   =',4(i2,4x))
 6070 format(6x,'number of open shells     =',4(i2,4x))
 6080 format(6x,'open shell occupation     =',4(i2,4x))
      end
