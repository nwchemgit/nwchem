      subroutine atcond(zn,ncsh,nosh,nccup,ajmn,nsym,znps)
C$Id: atcond.F,v 1.4 1997-04-29 19:50:28 d3e129 Exp $
      implicit none
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
*
*
*  This is where you want to read user input of the occupancy.
*  results from this routine are ncsh, noch, nccup, ajmn, and nsym
*  ncsh number of closed shells
*  nosh number of open shells
*  nccup occupancy of open shells
*  nsym number of symmetrires (s,p,d,f) that are occupied
*  ajmn K integrals (see below) (may assume only one open shell?)
*  
*
      real *8 zn, znps
      real *8 ajmn(*)
      integer ncsh(*),nccup(*),nosh(*)
      integer nelhw(103),nelcep(103)
      integer isymhw(13)
      integer isymce(13)
      integer nsym
c     
      integer i, j, k, nz, nlast, ksym, ihw
      integer nzps, idiff, nshell, isymm, nelec
      real *8 t, f
c     
      data isymhw / 1,1,2,1,2, 3,1,2, 3,1,2, 4, 3/
      data isymce / 1,1,2,1,2, 3,1,2, 3, 4,1,2, 3/
      data nelhw/
     $     0, 0,
     $     2, 2, 2, 2, 2, 2, 2, 2,
     $     10, 10, 10, 10, 10, 10, 10, 10,
     $     18, 18,
     $     18, 18, 18, 18, 18,
     $     18, 18, 18, 18, 18,
     $     28, 28, 28, 28, 28, 28,
     $     36,36,36,36,36,36,36,36,36,36,36,36,
     $     46,46,46,46,46,46,54,54,54,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0,68,
     $     68,68,68,68,68,68,68,68,68,78,78,0,
     $     0,0,0,0,0,0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0   /
      data nelcep /
     $     0, 0,
     $     2, 2, 2, 2, 2, 2, 2, 2,
     $     10, 10, 10, 10, 10, 10, 10, 10,
     $     18, 18,
     $     10, 10, 10, 10, 10,
     $     10, 10, 10, 10, 10,
     $     10, 28, 28, 28, 28, 28,
     $     36,36,28,28,28,28,28,28,28,28,28,28,
     $     28,46,46,46,46,46,54,54,46,0,0,0,
     $     0,0,0,0,0,0,0,0,0,0,0,60,
     $     60,60,60,60,60,60,60,60,60,78,78,78,
     $     78,78,0,0,0,0,0,0,0,0,0,0,
     $     0,0,0,0,0,0,0   /
c.......................................................................
c     initialize .
c.......................................................................
      nsym = 0
      do i = 1 , 4
         ncsh(i) = 0
         nosh(i) = 0
         nccup(i) = 0
      enddo
c.......................................................................
      nz = zn + 0.1d0
      nzps = znps + 0.1d0
      idiff = nz - nzps
      ihw = 0
      if (idiff.ne.0) then
c     determine pseudopotential type (hw or cep)
         if (nelhw(nz).eq.idiff) then
            ihw = 1
         else if (nelcep(nz).eq.idiff) then
            ihw = 2
         else
            call errquit
     &         ('atomscf:atcond:unrecognised pseudopotential',0)
         endif
      endif
      nzps = idiff
c.......................................................................
c     
c     fill up shells - all electron first
c.......................................................................
      nlast = 0
      do i = 1 , 4
         do j = 1 , 2
            ksym = i
            do k = 1 , i
               nshell = nlast + 4*ksym - 2
               if (nz.lt.nshell) go to 60
               nsym = max0(nsym,ksym)
               ncsh(ksym) = ncsh(ksym) + 1
               ksym = ksym - 1
               nlast = nshell
            enddo
         enddo
      enddo
c.......................................................................
c     
c     check if open shell atom. test for la and ac.
c.......................................................................
 60   if (nz.eq.57 .or. nz.eq.89) then
         ncsh(4) = ncsh(4) - ncsh(4)/2
         ksym = 3
      end if
c     
c     ... now consider pseudopotentials ...
c     ... decrease ncsh according to nelcor (nzps)
c     
      if (nzps.gt.0) then
         nelec = 0
         if(ihw.eq.1) then
            do i = 1,13
               nelec = nelec + 4*isymhw(i) - 2
               if (nelec.le.nzps) then
                  isymm = isymhw(i)
                  ncsh(isymm) = ncsh(isymm) - 1
               else
                  go to 90
               endif
            enddo
         else
            do i = 1,13
               nelec = nelec + 4*isymce(i) - 2
               if (nelec.le.nzps) then
                  isymm = isymce(i)
                  ncsh(isymm) = ncsh(isymm) - 1
               else
                  go to 90
               endif
            enddo
         endif
      endif
c     
c     check for cases where pseudo-potential changes  state
c     (cu  pseudo incl d10 => s1)
c     
 90   if (nz.eq.29.and.nzps.eq.28) then
         nlast = 28
         ksym = 1
      end if
c     
      do i = 1 , 24
         ajmn(i) = 0.0d0
      enddo
      if (nz.ne.nlast) then
         nosh(ksym) = 1
         nsym = max(nsym,ksym)
         nccup(ksym) = nz - nlast
c.......................................................................
c     
c     set k(l,l,0).
c.......................................................................
         i = ksym*(ksym+1)*(ksym+2)/6 - ksym + 1
         t = dble(nccup(ksym))
         f = 4*ksym - 4
         ajmn(i) = -(f+1.0d0)/(f+2.0d0) + (t-1.0d0)/t
      end if
      return
      end
