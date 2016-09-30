!****************************************************************************
!
! Find real roots of quadratic polynomial aa*x**2+bb*x+cc=0
! nroots = number of distinct real roots
! for nroots = 2 or 1, solutions are x = root1, x = root2
! for nroots = 0, solutions are x = root1 +/- i root2

subroutine rquadroots(aa,bb,cc,nroots,root1,root2)
implicit none
double precision :: aa,bb,cc
double precision :: root1,root2
integer :: nroots
double precision :: qq,sqq,tt

  qq=bb**2-4*aa*cc
  if (qq.lt.0.d0) then
    nroots=0
    root1 = -bb/(2*aa)  ! real part of roots
    root2 = sqrt(abs(qq))/(2*aa) ! +/- imaginary part of roots
    return
  elseif (qq.eq.0.d0) then
    nroots=1
    root1=-bb/(2*aa)
    root2=root1
  else
    nroots=2
    sqq=sqrt(qq)
    if (bb.ge.0.d0) then
      tt=-0.5*(bb+sqq)
    else
      tt=-0.5*(bb-sqq)
    endif
    root1=tt/aa
    root2=cc/tt
!  root1=(-bb-sqq)/(2*aa)
!  root2=(-bb+sqq)/(2*aa)
  endif

  return
end subroutine rquadroots

!****************************************************************************
!
! Find roots of cubic polynomial aa*x**3 + bb*x**2 + cc*x + dd
! nroots = number of real roots
! for nroots = 3, solutions are x = root1, x = root2, x = root3
! for nroots = 1, solutions are x = root1, x = root2 +/- i root3

subroutine rcubicroots(aa,bb,cc,dd,nroots,root1,root2,root3)
implicit none
double precision :: aa,bb,cc,dd
double precision :: root1,root2,root3
integer :: nroots
double precision :: bp,cp,dp
double precision :: QQ,RR,theta,adscr,bdscr,QQ3,sqQQ3,sqQQ,twopi
integer :: ii,jj,kk,ll
double precision :: rt(4),tt(4)
common /roottest/ rt

  twopi = 2*acos(-1.d0)
  bp = bb/aa
  cp = cc/aa
  dp = dd/aa
  QQ = (bp*bp - 3*cp)/9.d0
  RR = (2*bp**3-9*bp*cp+27*dp)/54.d0
  QQ3 = QQ**3
  if (RR*RR .lt. QQ3) then                   ! three real roots
    sqQQ3 = sqrt(QQ3)
    sqQQ = sqrt(QQ)
    theta = acos(-RR/sqQQ3)
    nroots = 3
    root1 = 2*sqQQ*cos(theta/3.d0) - bp/3.d0
    root2 = 2*sqQQ*cos((theta+twopi)/3.d0) - bp/3.d0
    root3 = 2*sqQQ*cos((theta-twopi)/3.d0) - bp/3.d0
  else                                       ! one real and two complex roots
    nroots = 1
    if (RR .ge. 0.d0) then
      adscr = - (abs(RR) + sqrt(RR*RR-QQ3))**(1.d0/3.d0)
    else
      adscr = (abs(RR) + sqrt(RR*RR-QQ3))**(1.d0/3.d0)
    endif
    if (adscr .eq. 0.d0) then
      bdscr = 0.d0
    else
      bdscr = QQ/adscr
    endif
    root1 = adscr + bdscr - bp/3.d0          ! real root
    root2 = - 0.5*(adscr + bdscr) - bp/3.d0  ! real part of complex roots
    root3 = 0.5*sqrt(3.d0)*(adscr-bdscr)     ! imaginary part of complex roots
  endif

  return
end subroutine rcubicroots

!****************************************************************************
!
! Find roots of quartic polynomial aa*x**4 + bb*x**3 + cc*x**2 + dd*x + ee
! nroots = number of real roots
! for nroots = 4, solutions are x = root1, x = root2, x = root3, x = root4
! for nroots = 2, solutions are x = root1, x = root2, x = root3 +/- i root4
! for nroots = 0, solutions are x = root1 +/- root2,  x = root3 +/- i root4

subroutine rquarticroots(aa,bb,cc,dd,ee,nroots,root1,root2,root3,root4)
implicit none
double precision :: aa,bb,cc,dd,ee
double precision :: root1,root2,root3,root4
integer :: nroots
double precision :: bp,cp,dp,ep
double precision :: alpha,beta,gam,yy,ytrm(0:2),temp,gg,hh,sgg,shh
integer :: iroots,jroots
double complex :: uu2,uu
! DEBUG vars
integer :: ii,jj,kk,ll
double precision :: PP,QQ,RR,UD
double complex :: vv1,vv2,ww1,ww2,ww3
double precision :: rt(4),ut(4),yt(3)
common /roottest/ rt

  bp = bb/aa
  cp = cc/aa
  dp = dd/aa
  ep = ee/aa
  alpha = -3*bp*bp/8.d0 + cp
  beta = bp**3/8.d0 - bp*cp/2.d0 + dp
  gam = -(3*bp**4)/256.d0 + cp*bp*bp/16.d0 - bp*dp/4.d0 + ep
! now have u^4 + alpha u^2 + beta u + gam = 0
! x = u - bb/(4*aa)
  if (beta.eq.0.d0) then
!   u^4 + alpha u^2 + gam = 0, quadratic in uu^2
    call rquadroots(1.d0,alpha,gam,nroots,root1,root2)
    if (nroots.gt.0) then
      if (root1.lt.root2) then
        temp = root2
        root2 = root1
        root1 = temp
      endif
      nroots = 4
      if (root2.ge.0.d0) then
        root4 = -sqrt(root2) - bb/(4*aa)
        root3 =  sqrt(root2) - bb/(4*aa)
      else
        nroots = 2
        root3 = - bb/(4*aa)
        root4 = sqrt(abs(root2))
      endif
      if (root1.ge.0.d0) then
        root2 = -sqrt(root1) - bb/(4*aa)
        root1 =  sqrt(root1) - bb/(4*aa)
      else
        nroots = nroots - 2
        root1 = - bb/(4*aa)
        root2 = sqrt(abs(root1))
      endif
    else
      uu2 = root1 + (0.d0,1.d0)*root2
      uu = sqrt(uu2)
      nroots = 0
      root1 = dble(uu) - bb/(4*aa)
      root2 = dimag(uu)
      root3 = -dble(uu) - bb/(4*aa)
      root4 = root2
    endif
  else if (gam.eq.0.d0) then
! (u^3 + alpha u + beta) * u = 0
    root1 = - bb/(4*aa)
    call rcubicroots(1.d0,0.d0,alpha,beta,nroots,root2,root3,root4)
    nroots = nroots + 1
    root2 = root2 - bb/(4*aa)
    root3 = root3 - bb/(4*aa)
    if (nroots.eq.4) root4 = root4 - bb/(4*aa)
  else
    ytrm(2) = 2.5*alpha
    ytrm(1) = 2*alpha**2-gam
    ytrm(0) = 0.5*(alpha**3-alpha*gam-beta*beta/4)
    call rcubicroots(1.d0,ytrm(2),ytrm(1),ytrm(0),nroots,root1,root2,root3)
    yy = root1
    if (nroots.gt.1) yy = max(root1,root2,root3)
    gg = sqrt(max(3.d-308,alpha+2*yy))
    hh = alpha + yy - beta/(2*gg)
    call rquadroots(1.d0,gg,hh,iroots,root1,root2)
    gg = -gg
    hh = alpha + yy - beta/(2*gg)
    call rquadroots(1.d0,gg,hh,jroots,root3,root4)
    nroots = 0
    if (iroots.gt.0) then
      nroots = nroots+2
      root1 = root1 - bb/(4*aa)
      root2 = root2 - bb/(4*aa)
    else
      root1 = root1 - bb/(4*aa)
    endif
    if (jroots.gt.0) then
      nroots = nroots+2
      root3 = root3 - bb/(4*aa)
      root4 = root4 - bb/(4*aa)
    else
      root3 = root3 - bb/(4*aa)
    endif
    if (jroots.gt.iroots) then
      temp = root1
      root1 = root3
      root3 = temp
      temp = root2
      root2 = root4
      root4 = temp
    endif
  endif

  return
end subroutine rquarticroots

!****************************************************************************
!
! find sqrt(1+xx)-1
! useful when xx is small to avoid numerical error

double precision function sqrtm1(xx)
double precision :: xx
double precision :: sum1,sum1_old,prod1
integer :: kk

  if (abs(xx).gt.0.03) then
    sqrtm1 = sqrt(1.d0+xx)-1.d0
  else
    sqrtm1 = 0.d0
    sum1 = 0.d0
    prod1 = 1.d0
    kk = 1
    do
      prod1 = prod1*(-1*(2*kk)*(2*kk-1)*xx)/dble(kk*kk*4)
      sum1_old = sum1
      sum1 = sum1 + prod1/dble(1-2*kk)
      if (sum1.eq.sum1_old) exit
      kk = kk+1
    enddo
    sqrtm1 = sum1
  endif
  return

end function sqrtm1
