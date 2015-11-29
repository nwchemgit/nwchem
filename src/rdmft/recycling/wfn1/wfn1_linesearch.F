C> \ingroup wfn1
C> @{
C>
C> \brief Perform a line search for minimum
C>
C> This routine performs a line search the find a minimum of the 
C> energy. In this implementation a second order approximation is
C> used in that a parabola is fitted to 3 data elements. The minimum
C> of the parabola is taken as a guess for the true minimum. 
C> This general idea is repeated in a cycle until the gradient at the
C> last point is small enough.
C>
C> The algorithm behind this uses a start up phase and an iterative
C> process. In both phases a parabola is fitted to available data and
C> its minimum established. In the start up phase a parabola is
C> fitted to f(0), df(0) and f(x1). In the iterative phases a different
C> appraoch based on known function values is used. I.e. the input
C> data is f(x1), f(x2), and f(x3). 
C>
C> In short the start up phase works as follows:
C> - Given is f(0) and df(0)
C> - x1 = -df(0), calc f(x1)
C> - fit parabola to f(0), df(0), and f(x1), and find its minimum xm
C> - Calculate f(xm)
C> - Add (0,f(0)), (x1,f(x1)), and (xm,f(xm)) to an ordered table
C>
C> The iterative process then proceeds:
C> - Find the point xg with the lowest energy f(xg)
C> - Find two points closest to xg in position
C> - Call the sequences of these three points x1, x2, and x3, along
C>   with associated values f(x1), f(x2), and f(x3)
C> - Fit a parabola to the three data points and establish the 
C>   corresponding minum xm
C> - Calculate f(xm)
C> - Add (xm,f(xm)) to the table
C>
      subroutine wfn1_linesearch(rtdb,geom,uwfn1,ist,nbf,nea,neb,
     +           h1,eri,erix,ov,e0,
     +           lambda,temp,pow,tol,nperma,npermb,perma,permb,
     +           signa,signb,ovla,ovlb,
     +           r_noa,r_nob,r_ofa,r_ofb,r_noa2,r_nob2,r_ofa2,r_ofb2,
     +           d_noa,d_nob,d_ofa,d_ofb,d_noa2,d_nob2,d_ofa2,d_ofb2,
     +           rnoa,rnob,rofa,rofb,m1,m2,m3,step)
      implicit none
      integer rtdb !< [Input] The RTDB handle
      integer geom !< [Input] The geometry handle
      logical uwfn1 !< [Input] If .true. do unrestricted WFN1
      integer ist  !< [Input] The state number
      integer nbf  !< [Input] The number of basis functions
      integer nea  !< [Input] The number of alpha electrons
      integer neb  !< [Input] The number of beta electrons
      double precision h1(nbf,nbf) !< [Input] The 1-electron Hamiltonian
      double precision eri(nbf,nbf,nbf,nbf) !< [Input] The Coulomb
      !< 2-electron repulsion integrals
      double precision erix(nbf,nbf,nbf,nbf) !< [Input] The exchange
      !< 2-electron integrals
      double precision ov(nbf,nbf) !< [Input] The overlap integrals
      double precision e0 !< [Input] Energy of current wavefunction
      double precision lambda(ist) !< [Input] Lambda of current wavefunction
      double precision temp !< [Input] The temperature (Fermi-smearing)
      double precision pow !< [Input] The density matrix power (DMFT)
      double precision tol !< [Input] The convergence tolerance
      integer nperma !< [Input] The number of alpha electron permutations
      integer npermb !< [Input] The number of beta electron permutations
      integer perma(nea,nperma) !< [Input] The alpha electron permutations
      integer permb(neb,npermb) !< [Input] The beta electron permutations
      double precision signa(nperma) !< [Input] The sign of the alpha
                                     !< permutations
      double precision signb(npermb) !< [Input] The sign of the beta
                                     !< permutations
      double precision ovla(nea,nea) !< [Scratch] Overlap of alpha
                                     !< orbitals
      double precision ovlb(neb,neb) !< [Scratch] Overlap of beta
                                     !< orbitals
      double precision r_noa(nbf,nbf,ist) !< [Input] The alpha natural 
      !< orbitals
      double precision r_nob(nbf,nbf,ist) !< [Input] The beta natural 
      !< orbitals
      double precision r_ofa(nbf,nbf,ist) !< [Input] The alpha
      !< occupation functions
      double precision r_ofb(nbf,nbf,ist) !< [Input] The beta occupation
      !< functions
      double precision r_noa2(nbf,nbf,ist) !< [Scratch] The guess alpha
      !< natural orbitals
      double precision r_nob2(nbf,nbf,ist) !< [Scratch] The guess beta
      !< natural orbitals
      double precision r_ofa2(nbf,nbf,ist) !< [Scratch] The guess alpha
      !< occupation functions
      double precision r_ofb2(nbf,nbf,ist) !< [Scratch] The guess beta
      !< occupation functions
      double precision d_noa(nbf,nbf) !< [Input] The alpha natural 
      !< orbitals gradient
      double precision d_nob(nbf,nbf) !< [Input] The beta natural 
      !< orbitals gradient
      double precision d_ofa(nbf,nbf) !< [Input] The alpha occupation 
      !< functions gradient
      double precision d_ofb(nbf,nbf) !< [Input] The beta occupation 
      !< functions gradient
      double precision d_noa2(nbf,nbf) !< [Scratch] The guess alpha
      !< natural orbitals step
      double precision d_nob2(nbf,nbf) !< [Scratch] The guess beta
      !< natural orbitals step
      double precision d_ofa2(nbf,nbf) !< [Scratch] The guess alpha
      !< occupation functions step
      double precision d_ofb2(nbf,nbf) !< [Scratch] The guess beta
      !< occupation functions step
      double precision rnoa(nbf,nbf) !< [Scratch] The alpha
      !< natural orbitals rotation
      double precision rnob(nbf,nbf) !< [Scratch] The beta
      !< natural orbitals rotation
      double precision rofa(nbf,nbf) !< [Scratch] The alpha
      !< occupation functions rotation
      double precision rofb(nbf,nbf) !< [Scratch] The beta
      !< occupation functions rotation
      double precision m1(nbf,nbf) !< [Scratch]
      double precision m2(nbf,nbf) !< [Scratch]
      double precision m3(nbf,nbf) !< [Scratch]
      double precision step !< [Output] The scale factor that defines
      !< the recommended step
c
      integer ii,jj,jst
      double precision dnrm
      double precision en0,on0(ist)
      integer mxpts, npts
      parameter(mxpts = 100)
      double precision xx(mxpts) !< The coordinate
      double precision fx(mxpts) !< The function value f(xx(i))
      double precision trst      !< The trust region
      double precision xm        !< The coordinate of the predicted minimum
      double precision fxm       !< The predicted minimum value
      double precision a         !< Polynomial coefficient
      double precision b         !< Polynomial coefficient
      double precision c         !< Polynomial coefficient
      double precision df        !< Estimate of gradient at a point
      integer ix1, ix2, ix3, indxmin
      double precision dum1,dum2,dum3
      integer wfn1_idmin
      external wfn1_idmin
c
      trst = 2.0d0*acos(-1.0d0)
c
      npts = 1
      xx(npts) = 0.0d0
      fx(npts) = e0
c
      dnrm = 0.0d0
      do jj = 1, nbf
        do ii = 1, nbf
          dnrm = dnrm + d_noa(ii,jj)**2
          dnrm = dnrm + d_nob(ii,jj)**2
          dnrm = dnrm + d_ofa(ii,jj)**2
          dnrm = dnrm + d_ofb(ii,jj)**2
        enddo
      enddo
      dnrm = sqrt(dnrm*1.000d0)
      do jj = 1, nbf
        do ii = 1, nbf
          d_noa2(ii,jj) = -d_noa(ii,jj)
          d_nob2(ii,jj) = -d_nob(ii,jj)
          d_ofa2(ii,jj) = -d_ofa(ii,jj)
          d_ofb2(ii,jj) = -d_ofb(ii,jj)
        enddo
      enddo
c
      call wfn1_exp(nbf,d_noa2,m1,m2,rnoa)
      call wfn1_exp(nbf,d_nob2,m1,m2,rnob)
      call wfn1_exp(nbf,d_ofa2,m1,m2,rofa)
      call wfn1_exp(nbf,d_ofb2,m1,m2,rofb)
      if (.not.uwfn1) then
        call dcopy(nbf*nbf,rnoa,1,rnob,1)
        call dcopy(nbf*nbf,rofa,1,rofb,1)
      endif
c
      call dcopy(nbf*nbf*ist,r_noa,1,r_noa2,1)
      call dcopy(nbf*nbf*ist,r_nob,1,r_nob2,1)
      call dcopy(nbf*nbf*ist,r_ofa,1,r_ofa2,1)
      call dcopy(nbf*nbf*ist,r_ofb,1,r_ofb2,1)
c
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,rnoa,nbf,
     +           0.0d0,r_noa2(1,1,ist),nbf)
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,rnob,nbf,
     +           0.0d0,r_nob2(1,1,ist),nbf)
c
      if (ist.eq.1) then
        call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnoa,nbf,
     +             r_ofa(1,1,ist),nbf,0.0d0,m1,nbf)
      else
        call dcopy(nbf*nbf,r_ofa(1,1,ist),1,m1,1)
      endif
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,rofa,nbf,0.0d0,
     +           r_ofa2(1,1,ist),nbf)
c
      if (ist.eq.1) then
        call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnob,nbf,
     +             r_ofb(1,1,ist),nbf,0.0d0,m1,nbf)
      else
        call dcopy(nbf*nbf,r_ofb(1,1,ist),1,m1,1)
      endif
      call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,rofb,nbf,0.0d0,
     +           r_ofb2(1,1,ist),nbf)
c
      call wfn1_engrad(rtdb,geom,ist,nbf,nea,neb,h1,eri,erix,ov,
     +                 dum1,dum2,
     +                 r_noa2,r_nob2,r_ofa2,r_ofb2,m1(1,1),m1(1,2),
     +                 m2(1,1),m2(1,2),m3(1,1),m3(1,2),temp,pow,
     +                 nperma,npermb,perma,permb,signa,signb,ovla,ovlb,
     +                 en0,dum3,on0,.false.)
      do jst = 1, ist-1
        en0 = en0 + lambda(jst)*on0(jst)
      enddo
      call wfn1_insert(mxpts,npts,xx,fx,dnrm,en0)
      call wfn1_f0df0f1(dnrm,e0,-dnrm,en0,trst,a,b,c,xm,fxm)
      df = 2.0d0*a*dnrm+b
cDEBUG
c     if (a.lt.0.0d0.and.xm.lt.0.0d0) then
c       xm = -xm
c     endif
c     write(*,*)"*** f0df0f1 a:      en0,a,b,c= ",en0,a,b,c
c     write(*,*)"*** f0df0f1 b: dnrm,xm,fxm,df= ",dnrm,xm,fxm,df
cDEBUG
      step = dnrm/dnrm
c
      do while (abs(df).gt.tol.and.xm.ge.1.0d-16.and.npts.lt.mxpts)
        do jj = 1, nbf
          do ii = 1, nbf
            d_noa2(ii,jj) = -xm/dnrm*d_noa(ii,jj)
            d_nob2(ii,jj) = -xm/dnrm*d_nob(ii,jj)
            d_ofa2(ii,jj) = -xm/dnrm*d_ofa(ii,jj)
            d_ofb2(ii,jj) = -xm/dnrm*d_ofb(ii,jj)
          enddo
        enddo
c
        call wfn1_exp(nbf,d_noa2,m1,m2,rnoa)
        call wfn1_exp(nbf,d_nob2,m1,m2,rnob)
        call wfn1_exp(nbf,d_ofa2,m1,m2,rofa)
        call wfn1_exp(nbf,d_ofb2,m1,m2,rofb)
        if (.not.uwfn1) then
          call dcopy(nbf*nbf,rnoa,1,rnob,1)
          call dcopy(nbf*nbf,rofa,1,rofb,1)
        endif
c
        call dcopy(nbf*nbf*ist,r_noa,1,r_noa2,1)
        call dcopy(nbf*nbf*ist,r_nob,1,r_nob2,1)
        call dcopy(nbf*nbf*ist,r_ofa,1,r_ofa2,1)
        call dcopy(nbf*nbf*ist,r_ofb,1,r_ofb2,1)
c
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_noa(1,1,ist),nbf,
     +             rnoa,nbf,0.0d0,r_noa2(1,1,ist),nbf)
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,r_nob(1,1,ist),nbf,
     +             rnob,nbf,0.0d0,r_nob2(1,1,ist),nbf)
c
        if (ist.eq.1) then
          call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnoa,nbf,
     +               r_ofa(1,1,ist),nbf,0.0d0,m1,nbf)
        else
          call dcopy(nbf*nbf,r_ofa(1,1,ist),1,m1,1)
        endif
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,rofa,nbf,0.0d0,
     +             r_ofa2(1,1,ist),nbf)
c
        if (ist.eq.1) then
          call dgemm('t','n',nbf,nbf,nbf,1.0d0,rnob,nbf,
     +               r_ofb(1,1,ist),nbf,0.0d0,m1,nbf)
        else
          call dcopy(nbf*nbf,r_ofb(1,1,ist),1,m1,1)
        endif
        call dgemm('n','n',nbf,nbf,nbf,1.0d0,m1,nbf,rofb,nbf,0.0d0,
     +             r_ofb2(1,1,ist),nbf)
c
        call wfn1_engrad(rtdb,geom,ist,nbf,nea,neb,h1,eri,erix,ov,
     +                   dum1,dum2,
     +                   r_noa2,r_nob2,r_ofa2,r_ofb2,m1(1,1),m1(1,2),
     +                   m2(1,1),m2(1,2),m3(1,1),m3(1,2),temp,pow,
     +                   nperma,npermb,perma,permb,signa,signb,
     +                   ovla,ovlb,en0,dum3,on0,.false.)
        do jst = 1, ist-1
          en0 = en0 + lambda(jst)*on0(jst)
        enddo
        call wfn1_insert(mxpts,npts,xx,fx,xm,en0)
        indxmin = wfn1_idmin(npts,fx,1)
cDEBUG
c       do ii = 1, npts
c         write(*,*)'*** pts: ',ii,xx(ii),fx(ii)
c       enddo
c       write(*,*)'*** indxmin = ',indxmin
cDEBUG
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
cDEBUG
c     write(*,*)"*** closest a: ",ix1,ix2,ix3
c     write(*,*)"*** closest b: ",fx(ix1),fx(ix2),fx(ix3)
cDEBUG
        if (indxmin.eq.1) then
          call wfn1_f0df0f1(xx(ix2),e0,-dnrm,fx(ix2),trst,a,b,c,xm,fxm)
        else
          call wfn1_f1f2f3(xx(ix1),xx(ix2),xx(ix3),
     +                     fx(ix1),fx(ix2),fx(ix3),trst,a,b,c,xm,fxm)
        endif
        df = 2.0d0*a*xx(indxmin)+b
cDEBUG
c     write(*,*)"*** f1f2f3 a:    en0,a,b,c= ",en0,a,b,c
c     write(*,*)"*** f1f2f3 b: xx,xm,fxm,df= ",xx(indxmin),xm,fxm,df
cDEBUG
        step = xx(indxmin)/dnrm
      enddo
cDEBUG
c     write(*,*)'*** step: ',step
c     write(*,*)
cDEBUG
c
      end
C>
C> @}
