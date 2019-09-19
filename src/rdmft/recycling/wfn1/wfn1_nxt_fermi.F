C> \ingroup wfn1_nxt
C> @{
C>
C> \brief Calculate the Fermi electron distribution over orbitals
C>
C> Calculate the occupation numbers for the natural orbitals using
C> Fermi-Dirac smearing. The occupation numbers are given by
C> \f{eqnarray*}{
C>   f_i &=& \frac{1}{1+\exp(\beta[\epsilon_i - \mu])} \\\\
C>   \beta &=& 1/T
C> \f}
C> The chemical potential \f$\mu\f$ is to be optimized so that
C> \f{eqnarray*}{
C>   N_e &=& \sum_i f_i
C> \f}
C> For a given \f$\mu\f$ and \f$T\f$ we can define the error in terms
C> of number of electrons as
C> \f{eqnarray*}{
C>    N_{d}        &=& \sum_i f_i \\\\
C>    d            &=& \left(N_e - N_d\right)^2 \\\\
C>    \mathrm{err} &=& \sqrt{d}
C> \f}
C> The gradient of \f$d\f$ wrt. \f$\mu\f$ is
C> \f{eqnarray*}{
C>   \frac{\partial d}{\partial\mu} &=& -2\left(N_e-N_d\right)
C>   \frac{\partial N_d}{\partial\mu} \\\\
C>   \frac{\partial N_d}{\partial\mu} &=&
C>   \sum_i\frac{\beta\exp\left(\beta[\epsilon_i-\mu]\right)}
C>         {\left(1+\exp\left(\beta[\epsilon_i-\mu]\right)\right)^2}
C> \f}
C> The Hessian of \f$d\f$ wrt. \f$\mu\f$ is
C> \f{eqnarray*}{
C>   \frac{\partial^2 d}{\partial\mu^2}
C>   &=& 2\left(\frac{\partial N_d}{\partial\mu}\right)^2
C>    -  2\left(N_e-N_d\right)\frac{\partial^2 N_d}{\partial\mu^2} \\\\
C>   \frac{\partial^2 N_d}{\partial\mu^2} &=&
C>   \sum_i 2\frac{\left(\beta\exp\left(\beta[\epsilon_i-\mu]\right)\right)^2}
C>                {\left(1+\exp\left(\beta[\epsilon_i-\mu]\right)\right)^3}
C>  -\sum_i \frac{\beta^2\exp\left(\beta[\epsilon_i-\mu]\right)}
C>                {\left(1+\exp\left(\beta[\epsilon_i-\mu]\right)\right)^2}
C> \f}
C> With the equations above we can define a quasi-Newton update for
C> \f$\mu\f$ as
C> \f{eqnarray*}{
C>   \mu_{i+1} &=& \mu_i - \frac{\partial d}{\partial\mu}/\frac{\partial^2 d}{\partial\mu^2}
C> \f}
C> A good starting value for \f$\mu\f$ is 
C> \f{eqnarray*}{
C>    \mu_0 &=& \frac{1}{2}\left(\epsilon_{\mathrm{lumo}}-\epsilon_{\mathrm{homo}}\right)
C> \f}
C>
C>
C> ## References ##
C>
C> [1] R.W. Warren, B.I. Dunlap
C>     "Fractional occupation numbers and density functional
C>     energy gradients with the linear combination of Gaussian-
C>     type orbitals approach"
C>     Chemical Physics Letters 262 (1996) 384-392, DOI:
C>     <a href="https://doi.org/10.1016/0009-2614(96)01107-4">
C>     10.1016/0009-2614(96)01107-4</a>.
C>
      subroutine wfn1_nxt_fermi(e,f,em,nfocc,noccmx,t,ne,nmo,iwr,out)
      implicit none
#include "global.fh"
c
c     Description:
c     ------------
c
c     Calculate the occupation numbers for the orbitals using
c     Fermi-Dirac smearing:
c
c                 1
c        f_i = -------------------------, b = 1/t
c              1 + exp( b * [e_i - mu] )
c
c     mu is to be optimised such that sum_i f_i = Ne. 
c
c     Define Nd = sum_i f_i then the error is
c
c        d = (Ne - Nd) ** 2
c
c        error = sqrt(d)
c
c     The gradient of d with respect to mu is
c
c        du = -2 (Ne - Nd) dNdu
c
c                     b * exp( b * [e_i - mu] )
c        dNdu = sum_i -------------------------------
c                     ( 1 + exp( b * [e_i - mu] ))**2
c
c     The second derivative of d with respect to mu is
c
c        d2u = 2 dNdu ** 2 - 2 (Ne - Nd) d2Ndu
c
c                         ( b * exp( b * [e_i - mu] )) ** 2
c        d2dNdu = sum_i 2 ---------------------------------
c                         ( 1 + exp( b * [e_i - mu] )) ** 3
c
c                         b**2 * exp( b * [e_i - mu] )
c               - sum_i   ---------------------------------
c                         ( 1 + exp( b * [e_i - mu] )) ** 2
c
c     Now we can use a quasi newton update:
c
c         mu_i+1 = mu_i - du / d2dNdu
c
c     Start with mu = (E_lumo + E_homo)/2.
c
c
c     Also calculate the Mermin term to the total energy:
c
c        Em = 1/b * sum_i {(1-f_i)*ln(1-f_i) + f_i*ln(f_i)}
c
c     Adding this term to the total energy will give the correct
c     finite temperature energy. The gradient expressions should not
c     need any modification [1].
c
c     Comments:
c     ---------
c
c     1) The Fermi-Dirac smearing enforces strict aufbau ordering
c        of the orbitals through the occupations. Thus it cannot be
c        used together with options that may break the aufbau ordering
c        such as (large) level shifters, or locking.
c
c     2) The Fermi-Dirac smearing breaks the strict distinction 
c        between occupied and virtual orbitals. This has consequences
c        for the definition of the tester.
c
c     References:
c     -----------
c
c     [1] R.W. Warren, B.I. Dunlap
c         "Fractional occupation numbers and density functional
c          energy gradients with the linear combination of Gaussian-
c          type orbitals approach"
c         Chemical Physics Letters 262 (1996) 384-392.
c
c     Input
c
      integer nmo ! number of MO's
      integer ne  ! number of electrons
      double precision e(nmo) ! orbital energies
      double precision t      ! smearing parameter
      integer iwr ! unit number of standard output
      logical out ! .true. if output requested
c
c     Output
c
      double precision f(nmo)    ! occupation numbers
      double precision Em        ! the Mermin energy term
      integer nfocc  ! the number of fully occupied orbitals
      integer noccmx ! the number of (fully or partially) occupied 
                     ! orbitals
c
c     Local variables
c
cDEBUG
      double precision temp(nmo)
cDEBUG
      logical onotaufbau
      double precision mu
      double precision b,b2
      double precision Nd,dNdu,d2Ndu
      double precision z,z2,zp1,zp12,zp13
      double precision d,du,d2u
      double precision p,tt
      double precision Ehomo, Elumo
      integer i,j,k
      double precision tol
      parameter(tol=1.0d-8)
      integer maxit
      parameter(maxit=1000)
c
      onotaufbau = .false.
      if (ne.le.0) then
         do i = 1, nmo
            f(i) = 0.0d0
         enddo
         return
      endif
c     write(*,*)'t=',t
      tt = max(t,tol)
      b  = 1.0d0/tt
      b2 = b*b
c
c     First we need to calculate the starting mu for
c     which we need the LUMO and HOMO energy levels
c     We need to do this because if level shifters are
c     used the energy levels may be in the "wrong" order.
c
      do i = 1, nmo
         f(i) = 0.0d0
      enddo
      do i = 1, Ne
         k = nmo
         do j = 1, nmo-1
            if (f(j).lt.0.5d0.and.e(j).lt.e(k)) then
               k=j
            endif
         enddo
         f(k) = 1.0d0
      enddo
      Elumo = e(nmo)
      Ehomo = e(1)
      do j = 2, nmo-1
         if (f(j).gt.0.5d0) then
            Ehomo = max(Ehomo,e(j))
         else
            Elumo = min(Elumo,e(j))
         endif
      enddo
      mu = 0.5d0 * (Ehomo + Elumo)
c
c     Now optimise mu to conserve the number of electrons
c
      if (out.and.ga_nodeid().eq.0) then
         write(iwr,"(2(/1x,a4,1x,a5,10x,a12,8x,a6,1x,a9))")
     +        "iter","Nelec","Nelec_approx","Nerror","chem.pot.",
     +        "----","-----","------------","------","---------"
      endif
      j      = 0 
 10   j      = j+1
      Nd     = 0.0d0
      dNdu   = 0.0d0
      d2Ndu  = 0.0d0
      Em     = 0.0d0
      nfocc  = 0
      noccmx = 0
      do i = 1, nmo
         p = b*(e(i)-mu)
         if (p.gt.-dlog(tol)) then
            f(i)   = 0.0d0
         else if (p.lt.dlog(tol)) then
            f(i)   = 1.0d0
            Nd     = Nd+1.0d0
            nfocc  = nfocc+1
c           noccmx = noccmx+1
            noccmx = i
         else
            z      = dexp(p)
            z2     = z*z
            zp1    = 1.0d0+z
            zp12   = zp1*zp1
            zp13   = zp1*zp12
            Nd     = Nd + 1.0d0/zp1
            dNdu   = dNdu + b*z/zp12
            d2Ndu  = d2Ndu + 2.0d0*b2*z2/zp13 - b2*z/zp12
            f(i)   = 1.0d0/zp1
            Em     = Em + tt*((1.0d0-f(i))*dlog(1.0d0-f(i))+
     &                        f(i)*dlog(f(i)))
c           noccmx = noccmx+1
            noccmx = i
         endif
      enddo
cDEBUG
c     Em = 0.0d0
cDEBUG
      d   =  (Ne-Nd)
      du  = -2.0d0*(Ne-Nd)*dNdu
      d2u =  2.0d0*dNdu*dNdu-2.0d0*(Ne-Nd)*d2Ndu
      if (out.and.ga_nodeid().eq.0) then
         write(iwr,"(1x,i4,1x,i5,f22.8,f14.8,f10.4)")j,Ne,Nd,d,mu
      endif
      if (dabs(d)/Ne.gt.tol) then
         if (j.lt.maxit) then
            mu = mu - du/d2u
            goto 10
         endif
      endif
c
c     Check that the orbitals energies are in aufbau
c     ordering. Otherwise the Fermi-Dirac smearing will mess
c     everything up as it insists to occuppy the orbitals 
c     according to aufbau.
c
c     Aufbau ordering is checked against the occupations.
c     Because interchanges of equally occupied orbitals
c     are insignificant and the orbital occupations are 
c     energy dependent this seems the best way to do it.
c
      do i=1,nmo-1
         if (f(i+1).gt.f(i)+tol) then
c           write(iwr,"(2(/1x,a3,a22))")
c    +                "mo","energy",
c    +                "--","------"
c           do j = 1, nmo
c              write(iwr,"(1x,i3,f22.8)")j,e(j)
c           enddo
            if (ga_nodeid().eq.0) then
              write(iwr,"(1x,'WARNING: fermi_smear: aufbau violation',
     +                    ' detected')")
              write(iwr,"(1x,'orbitals ',i5,' and ',i5,' in non-aufbau',
     +                    ' order')")i,i+1
            endif
            onotaufbau = .true.
            call errquit("fermi_smear: aufbau violation detected")
         endif
      enddo
      if ((out.or.j.ge.maxit.or.onotaufbau).and.ga_nodeid().eq.0) then
         write(iwr,"(/1x,'fully occupied orbitals: ',i5)")nfocc
         write(iwr,"( 1x,'all   occupied orbitals: ',i5)")noccmx
         write(iwr,"(2(/1x,a3,2a22))")
     +             "mo","energy","occupation",
     +             "--","------","----------"
         do i = 1, nmo
            write(iwr,"(1x,i3,2f22.8)")i,e(i),f(i)
         enddo
         write(iwr,*)
      endif
      if (j.ge.maxit) then
         write(*,*)"fermi_smear: too many iterations"
         write(iwr,"(/1x,'fully occupied orbitals: ',i5)")nfocc
         write(iwr,"( 1x,'all   occupied orbitals: ',i5)")noccmx
         write(iwr,"(2(/1x,a3,2a22))")
     +             "mo","energy","occupation",
     +             "--","------","----------"
         do i = 1, nmo
            write(iwr,"(1x,i3,2f22.8)")i,e(i),f(i)
         enddo
         write(iwr,*)
c        call errquit("fermi_smear: too many iterations",0,0)
      endif
      end
C>
C> @}
