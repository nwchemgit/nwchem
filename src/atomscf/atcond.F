      subroutine atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps,spec) 
* $Id$
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c.......................................................................
c
c     for atom of nuclear charge zn find electron configuration
c     and k(l,l,0) coupling coefficient. simple aufbau is
c     assumed throughout the periodic system. the algorithm used
c     will work for zn less than 119.
c
c     to allow for pseudo-potentials the effective charge znps
c     is used to determine the number of closed shells filled
c     jvl  (daresbury 1988): corrected July 93
c.......................................................................
      double precision zn
      double precision znps
      double precision ajmn(*)
      integer nsym
      integer ncsh(4), nccup(4), nosh(4)
c
c... local
c
      integer maxbas, ind
      integer i, j, k, nz, nzps, idiff, ihw, nlast, ksym
      integer nshell, nelec, isymm
      integer ncore
      double precision t, f
c
      logical spec
c
c...  Begin ECP-data
c
      integer nelecp(26), iecp(26)
      integer isymhw(17)
      integer isymce(17)
      integer isymab(17)
c
c...  Different schemes to fill the core orbitals in the ECP
c...  (1=s, 2=p, 3=d, 4=f)
c
      data isymhw / 1, 1,2, 1,2, 3, 1,2, 3 ,1,2, 4, 3, 1,2, 4, 3/
      data isymce / 1, 1,2, 1,2,3, 1,2,3,4, 1,2,3,4, 1,2,3/
      data isymab / 1, 1, 2, 1, 3, 2, 1, 3, 2, 1, 4, 3, 2, 1, 4, 3, 2/
c
c...  nelecp(i) gives the number of electrons in ECP i.
c...  iecp(i)   gives a filling scheme for ECP i.
c
      data nelecp/
     &    0,2,4,10,12,18,22,28,30,36,40,46,48,54,60,62,68,72,78,80,
     &    86,92,94,100,104,110/
      data iecp/
     &    1,1,1,1 ,1 ,1 ,2 ,1 ,1 ,1 ,2 ,1 ,2 ,3 ,2 ,1 ,1 ,2 ,1 ,2 ,
     &    3 ,1 ,1 ,1  ,2  ,1  /
c
c...  End ECP-Data
c.......................................................................
c
c     initialize .
c.......................................................................
      nsym = 0
      maxbas = 4
      do 20 i = 1 , maxbas
        ncsh(i) = 0
        nosh(i) = 0
        nccup(i) = 0
 20   continue
c.......................................................................
      nz = zn + 0.1d0
      nzps = znps + 0.1d0
c
      if (nz.ge.58.and.nz.le.71.and.spec) then
         ncsh(1)=1
         ncsh(2)=1
         ncsh(3)=0
         ncsh(4)=0
         nsym=2
         return 
      endif
c
      idiff = nz - nzps
      ncore = idiff
      ihw = 0
      if (idiff.ne.0) then
c
c...     determine pseudopotential type (hw or cep or ab)
c
        do i = 1, 26
          if (nelecp(i).eq.idiff) then
            ihw = iecp(i)
            goto 100
          endif
        enddo
        write(6,600)idiff
 600    format(
     &      '*** The number of electrons ',i3,' in the ECP does not',
     &      /,'accord with a known set of completely filled shells !!!')
        call errquit('unrecognised pseudopotential',911, INPUT_ERR)
 100    continue
      endif
      nzps = idiff
c.......................................................................
c
c     fill up shells - all electron first
c.......................................................................
      nlast = 0
      do 50 i = 1 , maxbas
        do 40 j = 1 , 2
          ksym = i
          do 30 k = 1 , i
            nshell = nlast + 4*ksym - 2
            if (nz.lt.nshell) go to 60
            nsym = max0(nsym,ksym)
            ncsh(ksym) = ncsh(ksym) + 1
            ksym = ksym - 1
            nlast = nshell
 30       continue
 40     continue
 50   continue
c.......................................................................
c
c     check if open shell atom. test for la and ac.
c.......................................................................
 60   if (nz.eq.57 .or. nz.eq.89) then
        ncsh(4) = ncsh(4) - ncsh(4)/2
        ksym = 3
      end if
c.......................................................................
c
c ... now consider pseudopotentials ...
c ... decrease ncsh according to nelcor (nzps)
c
      if (nzps.gt.0) then
        nelec = 0
        if(ihw.eq.1) then
          do 80 i = 1,17
            nelec = nelec + 4*isymhw(i) - 2
            if (nelec.le.nzps) then
              isymm = isymhw(i)
              ncsh(isymm) = ncsh(isymm) - 1
            else
              go to 90
            endif
 80       continue
        else if(ihw.eq.2) then
          do 800 i = 1,17
            nelec = nelec + 4*isymce(i) - 2
            if (nelec.le.nzps) then
              isymm = isymce(i)
              ncsh(isymm) = ncsh(isymm) - 1
            else
              go to 90
            endif
 800      continue
        else
          do 810 i = 1,17
            nelec = nelec + 4*isymab(i) - 2
            if (nelec.le.nzps) then
              isymm = isymab(i)
              ncsh(isymm) = ncsh(isymm) - 1
            else
              go to 90
            endif
 810      continue
        endif
      endif
c
c     check for cases where pseudo-potential changes  state
c     (cu  pseudo incl d10 =>s1)
c
 90   if (nz.eq.29.and.nzps.eq.28) then
        nlast = 28
        ksym = 1
      end if
c
      do 70 i = 1 , 24
        ajmn(i) = 0.0d0
 70   continue
      if (nz.ne.nlast) then
        nosh(ksym) = 1
        nsym = max0(nsym,ksym)
        nccup(ksym) = nz - nlast
c.......................................................................
c
c     set k(l,l,0).
c.......................................................................
        ind = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
        t = dble(nccup(ksym))
        f = dble(4*ksym - 4)
        ajmn(ind) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
      end if
      return
      end
*_old_version:      subroutine old_atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps)
*_old_version:      implicit none
*_old_version:#include "errquit.fh"
*_old_version:c.......................................................................
*_old_version:c     
*_old_version:c     for atom of nuclear charge zn find electron configuration
*_old_version:c     and k(l,l,0) coupling coefficient. simple aufbau is
*_old_version:c     assumed throughout the periodic system. the algorithm used
*_old_version:c     will work for zn less than 119.
*_old_version:c     
*_old_version:c     to allow for pseudo-potentials the effective charge znps
*_old_version:c     is used to determine the number of closed shells filled
*_old_version:c     jvl  (daresbury 1988): corrected July 93
*_old_version:c.......................................................................
*_old_version:*
*_old_version:*
*_old_version:*  This is where you want to read user input of the occupancy.
*_old_version:*  results from this routine are ncsh, noch, nccup, ajmn, and nsym
*_old_version:*  ncsh number of closed shells
*_old_version:*  nosh number of open shells
*_old_version:*  nccup occupancy of open shells
*_old_version:*  nsym number of symmetrires (s,p,d,f) that are occupied
*_old_version:*  ajmn K integrals (see below) (may assume only one open shell?)
*_old_version:*  
*_old_version:*
*_old_version:      real *8 zn, znps
*_old_version:      real *8 ajmn(*)
*_old_version:      integer ncsh(*),nccup(*),nosh(*)
*_old_version:      integer nelhw(103),nelcep(103)
*_old_version:      integer isymhw(13)
*_old_version:      integer isymce(13)
*_old_version:      integer nsym
*_old_version:c     
*_old_version:      integer i, j, k, nz, nlast, ksym, ihw
*_old_version:      integer nzps, idiff, nshell, isymm, nelec
*_old_version:      real *8 t, f
*_old_version:c     
*_old_version:      data isymhw / 1,1,2,1,2, 3,1,2, 3,1,2, 4, 3/
*_old_version:      data isymce / 1,1,2,1,2, 3,1,2, 3, 4,1,2, 3/
*_old_version:      data nelhw/
*_old_version:     $     0, 0,
*_old_version:     $     2, 2, 2, 2, 2, 2, 2, 2,
*_old_version:     $     10, 10, 10, 10, 10, 10, 10, 10,
*_old_version:     $     18, 18,
*_old_version:     $     18, 18, 18, 18, 18,
*_old_version:     $     18, 18, 18, 18, 18,
*_old_version:     $     28, 28, 28, 28, 28, 28,
*_old_version:     $     36,36,36,36,36,36,36,36,36,36,36,36,
*_old_version:     $     46,46,46,46,46,46,54,54,54,0,0,0,
*_old_version:     $     0,0,0,0,0,0,0,0,0,0,0,68,
*_old_version:     $     68,68,68,68,68,68,68,68,68,78,78,0,
*_old_version:     $     0,0,0,0,0,0,0,0,0,0,0,0,
*_old_version:     $     0,0,0,0,0,0,0   /
*_old_version:      data nelcep /
*_old_version:     $     0, 0,
*_old_version:     $     2, 2, 2, 2, 2, 2, 2, 2,
*_old_version:     $     10, 10, 10, 10, 10, 10, 10, 10,
*_old_version:     $     18, 18,
*_old_version:     $     10, 10, 10, 10, 10,
*_old_version:     $     10, 10, 10, 10, 10,
*_old_version:     $     10, 28, 28, 28, 28, 28,
*_old_version:     $     36,36,28,28,28,28,28,28,28,28,28,28,
*_old_version:     $     28,46,46,46,46,46,54,54,46,0,0,0,
*_old_version:     $     0,0,0,0,0,0,0,0,0,0,0,60,
*_old_version:     $     60,60,60,60,60,60,60,60,60,78,78,78,
*_old_version:     $     78,78,0,0,0,0,0,0,0,0,0,0,
*_old_version:     $     0,0,0,0,0,0,0   /
*_old_version:c.......................................................................
*_old_version:c     initialize .
*_old_version:c.......................................................................
*_old_version:      nsym = 0
*_old_version:      do i = 1 , 4
*_old_version:         ncsh(i) = 0
*_old_version:         nosh(i) = 0
*_old_version:         nccup(i) = 0
*_old_version:      enddo
*_old_version:c.......................................................................
*_old_version:      nz = zn + 0.1d0
*_old_version:      nzps = znps + 0.1d0
*_old_version:      idiff = nz - nzps
*_old_version:      ihw = 0
*_old_version:      if (idiff.ne.0) then
*_old_version:c     determine pseudopotential type (hw or cep)
*_old_version:         if (nelhw(nz).eq.idiff) then
*_old_version:            ihw = 1
*_old_version:         else if (nelcep(nz).eq.idiff) then
*_old_version:            ihw = 2
*_old_version:         else
*_old_version:            call errquit
*_old_version:     &         ('atomscf:old_atcond:unrecognised pseudopotential',0)
*_old_version:         endif
*_old_version:      endif
*_old_version:      nzps = idiff
*_old_version:c.......................................................................
*_old_version:c     
*_old_version:c     fill up shells - all electron first
*_old_version:c.......................................................................
*_old_version:      nlast = 0
*_old_version:      do i = 1 , 4
*_old_version:         do j = 1 , 2
*_old_version:            ksym = i
*_old_version:            do k = 1 , i
*_old_version:               nshell = nlast + 4*ksym - 2
*_old_version:               if (nz.lt.nshell) go to 60
*_old_version:               nsym = max0(nsym,ksym)
*_old_version:               ncsh(ksym) = ncsh(ksym) + 1
*_old_version:               ksym = ksym - 1
*_old_version:               nlast = nshell
*_old_version:            enddo
*_old_version:         enddo
*_old_version:      enddo
*_old_version:c.......................................................................
*_old_version:c     
*_old_version:c     check if open shell atom. test for la and ac.
*_old_version:c.......................................................................
*_old_version: 60   if (nz.eq.57 .or. nz.eq.89) then
*_old_version:         ncsh(4) = ncsh(4) - ncsh(4)/2
*_old_version:         ksym = 3
*_old_version:      end if
*_old_version:c     
*_old_version:c     ... now consider pseudopotentials ...
*_old_version:c     ... decrease ncsh according to nelcor (nzps)
*_old_version:c     
*_old_version:      if (nzps.gt.0) then
*_old_version:         nelec = 0
*_old_version:         if(ihw.eq.1) then
*_old_version:            do i = 1,13
*_old_version:               nelec = nelec + 4*isymhw(i) - 2
*_old_version:               if (nelec.le.nzps) then
*_old_version:                  isymm = isymhw(i)
*_old_version:                  ncsh(isymm) = ncsh(isymm) - 1
*_old_version:               else
*_old_version:                  go to 90
*_old_version:               endif
*_old_version:            enddo
*_old_version:         else
*_old_version:            do i = 1,13
*_old_version:               nelec = nelec + 4*isymce(i) - 2
*_old_version:               if (nelec.le.nzps) then
*_old_version:                  isymm = isymce(i)
*_old_version:                  ncsh(isymm) = ncsh(isymm) - 1
*_old_version:               else
*_old_version:                  go to 90
*_old_version:               endif
*_old_version:            enddo
*_old_version:         endif
*_old_version:      endif
*_old_version:c     
*_old_version:c     check for cases where pseudo-potential changes  state
*_old_version:c     (cu  pseudo incl d10 => s1)
*_old_version:c     
*_old_version: 90   if (nz.eq.29.and.nzps.eq.28) then
*_old_version:         nlast = 28
*_old_version:         ksym = 1
*_old_version:      end if
*_old_version:c     
*_old_version:      do i = 1 , 24
*_old_version:         ajmn(i) = 0.0d0
*_old_version:      enddo
*_old_version:      if (nz.ne.nlast) then
*_old_version:         nosh(ksym) = 1
*_old_version:         nsym = max(nsym,ksym)
*_old_version:         nccup(ksym) = nz - nlast
*_old_version:c.......................................................................
*_old_version:c     
*_old_version:c     set k(l,l,0).
*_old_version:c.......................................................................
*_old_version:         i = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
*_old_version:         t = dble(nccup(ksym))
*_old_version:         f = 4*ksym - 4
*_old_version:         ajmn(i) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
*_old_version:      end if
*_old_version:      return
*_old_version:      end
