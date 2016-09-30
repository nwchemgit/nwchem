!**********************************************************************
!   This is Steve White's rewrite of Mike Teter's integration routine.
!   Modified by J. Rehr for complex integration.
!   Modified by L. Campbell for Fortran 90
!   The following is a listing of the arguments in the initial function
!   statement:
!      fn    -- routine requires external function statement in MAIN
!      xmin  -- lower limit
!      xmax  -- upper limit
!      abr   -- absolute tolerable error
!      rlr   -- relative tolerable error
!      nsing -- number of singularities or regions requiring
!                   special attention
!      xsing -- array of locations of singularities or endpoints
!                   of special regions
!      error -- output for routine error messages
!      numcal-- the number of times fn was called
!      maxns -- the maximum number of regions being considered simultaneously

double precision function grater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
implicit none
integer :: nsing,numcal,maxns,mx,nstack,ii,jj,kk,icount
parameter (mx=1500)
double precision :: fn,value,valu,fval(3,mx)
double precision :: xmin,xmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: xleft(mx),dx(3),wt(3)
double precision :: wt9(9),xsing(30)
external fn
logical :: atsing
save dx,wt,wt9
data dx/0.1127016653792583,0.5,0.8872983346207417/
data wt/0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9/0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater) and the region
! is removed from consideration (and a new region is the right-most).
  nstack=nsing+1
  maxns = nstack
  error=0.d0
  grater=0.d0
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
  xleft(1)=xmin
  xleft(nsing+2)=xmax
  if(nsing.gt.0) then
    do jj=1,nsing
      xleft(jj+1)=xsing(jj)
    enddo
  endif
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=xleft(ii+1)-xleft(ii)
    do jj=1,3
      fval(jj,ii)=fn(xleft(ii)+del*dx(jj))
    enddo
  enddo
  numcal = nstack * 3
  do
    if(nstack+3.ge.mx) then
      write(*,*) 'TOO MANY REGIONS'
      stop 0006
    endif
! Divide the rightmost region into three subregions.
    del=xleft(nstack+1)-xleft(nstack)
    xleft(nstack+3)=xleft(nstack+1)
    xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
    xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=0.d0
    valu=0.d0
    do jj=nstack,nstack+2
      del1=xleft(jj+1)-xleft(jj)
!      print*, 'fn call 2'
      fval(1,jj)=fn(xleft(jj)+dx(1)*del1)
      fval(3,jj)=fn(xleft(jj)+dx(3)*del1)
!      print*, 'fn call 2'
      numcal = numcal + 2
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac = del / (xmax - xmin)
    atsing = .false.
    if(frac .le. 1.0e-8) atsing = .true.
    if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.           &
&       (atsing .and. (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
      grater=grater+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo

end function grater


double complex function cgrater(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
implicit none
integer :: nsing,numcal,maxns,mx,nstack,ii,jj,kk,icount
parameter (mx=1500)
double complex :: fn,value,valu,fval(3,mx)
double precision :: xmin,xmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: xleft(mx),dx(3),wt(3)
double precision :: wt9(9),xsing(20)
external fn
logical :: atsing
save dx,wt,wt9
data dx/0.1127016653792583,0.5,0.8872983346207417/
data wt/0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9/0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in cgrater) and the region
! is removed from consideration (and a new region is the right-most).
  nstack=nsing+1
  maxns = nstack
  error=0.d0
  cgrater=(0.d0,0.d0)
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
  xleft(1)=xmin
  xleft(nsing+2)=xmax
  if(nsing.gt.0) then
    do jj=1,nsing
      xleft(jj+1)=xsing(jj)
    enddo
  endif
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=xleft(ii+1)-xleft(ii)
    do jj=1,3
      fval(jj,ii)=fn(xleft(ii)+del*dx(jj))
    enddo
  enddo
  numcal = nstack * 3
  do
    if(nstack+3.ge.mx) then
      write(*,*) 'TOO MANY REGIONS'
      stop 0006
    endif
! Divide the rightmost region into three subregions.
    del=xleft(nstack+1)-xleft(nstack)
    xleft(nstack+3)=xleft(nstack+1)
    xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
    xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=(0.d0,0.d0)
    valu=(0.d0,0.d0)
    do jj=nstack,nstack+2
      del1=xleft(jj+1)-xleft(jj)
!      print*, 'fn call 2'
      fval(1,jj)=fn(xleft(jj)+dx(1)*del1)
      fval(3,jj)=fn(xleft(jj)+dx(3)*del1)
!      print*, 'fn call 2'
      numcal = numcal + 2
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac = del / (xmax - xmin)
    atsing = .false.
    if(frac .le. 1.0e-8) atsing = .true.
    if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.           &
&       (atsing .and. (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
      cgrater=cgrater+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo

end function cgrater


! Copy of grater for two-dimensional integration
double precision function grater2(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
implicit none
integer :: nsing,numcal,maxns,mx,nstack,ii,jj,kk,icount
parameter (mx=1500)
double precision :: fn,value,valu,fval(3,mx)
double precision :: xmin,xmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: xleft(mx),dx(3),wt(3)
double precision :: wt9(9),xsing(20)
external fn
logical :: atsing
save dx,wt,wt9
data dx/0.1127016653792583,0.5,0.8872983346207417/
data wt/0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9/0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater2) and the region
! is removed from consideration (and a new region is the right-most).
  nstack=nsing+1
  maxns = nstack
  error=0.d0
  grater2=0.d0
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
  xleft(1)=xmin
  xleft(nsing+2)=xmax
  if(nsing.gt.0) then
    do jj=1,nsing
      xleft(jj+1)=xsing(jj)
    enddo
  endif
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=xleft(ii+1)-xleft(ii)
    do jj=1,3
      fval(jj,ii)=fn(xleft(ii)+del*dx(jj))
    enddo
  enddo
  numcal = nstack * 3
  do
    if(nstack+3.ge.mx) then
      write(*,*) 'TOO MANY REGIONS'
      stop 0006
    endif
! Divide the rightmost region into three subregions.
    del=xleft(nstack+1)-xleft(nstack)
    xleft(nstack+3)=xleft(nstack+1)
    xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
    xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=0.d0
    valu=0.d0
    do jj=nstack,nstack+2
      del1=xleft(jj+1)-xleft(jj)
!      print*, 'fn call 2'
      fval(1,jj)=fn(xleft(jj)+dx(1)*del1)
      fval(3,jj)=fn(xleft(jj)+dx(3)*del1)
!      print*, 'fn call 2'
      numcal = numcal + 2
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac = del / (xmax - xmin)
    atsing = .false.
    if(frac .le. 1.0e-8) atsing = .true.
    if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.           &
&       (atsing .and. (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
      grater2=grater2+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo

end function grater2


! Copy of grater for three-dimensional integration
double precision function grater3(fn,xmin,xmax,abr,rlr,nsing,xsing,error,numcal,maxns)
implicit none
integer :: nsing,numcal,maxns,mx,nstack,ii,jj,kk,icount
parameter (mx=1500)
double precision :: fn,value,valu,fval(3,mx)
double precision :: xmin,xmax,del,del1,dif,frac
double precision :: abr,rlr,error
double precision :: xleft(mx),dx(3),wt(3)
double precision :: wt9(9),xsing(20)
external fn
logical :: atsing
save dx,wt,wt9
data dx/0.1127016653792583,0.5,0.8872983346207417/
data wt/0.277777777777777778,0.4444444444444444444,0.2777777777777777778/
data wt9/0.0616938806304841571,0.108384229110206161,         &
&           0.0398463603260281088,0.175209035316976464,      &
&           0.229732989232610220,0.175209035316976464,       &
&           0.0398463603260281088,0.108384229110206161,      &
&           0.0616938806304841571  /

! nstack is the number of different intervals into which the
! integration region is currently divided. The number of regions can
! grow if more accuracy is needed by dividing the right-most region
! into three regions. The number of regions shrinks when the integral
! over the right-most region is accurate enough, in which case that
! integral is added to the total (stored in grater3) and the region
! is removed from consideration (and a new region is the right-most).
  nstack=nsing+1
  maxns = nstack
  error=0.d0
  grater3=0.d0
! The array xleft stores the boundary points of the regions.
! The singular points just govern the initial placement of the regions.
  xleft(1)=xmin
  xleft(nsing+2)=xmax
  if(nsing.gt.0) then
    do jj=1,nsing
      xleft(jj+1)=xsing(jj)
    enddo
  endif
! For each region, calculate the function and store at three selected points.
  do ii=1,nstack
    del=xleft(ii+1)-xleft(ii)
    do jj=1,3
      fval(jj,ii)=fn(xleft(ii)+del*dx(jj))
    enddo
  enddo
  numcal = nstack * 3
  do
    if(nstack+3.ge.mx) then
      write(*,*) 'TOO MANY REGIONS'
      stop 0006
    endif
! Divide the rightmost region into three subregions.
    del=xleft(nstack+1)-xleft(nstack)
    xleft(nstack+3)=xleft(nstack+1)
    xleft(nstack+1)=xleft(nstack)+del*dx(1)*2.
    xleft(nstack+2)=xleft(nstack+3)-del*dx(1)*2.
! The three data points already found for the region become the
! middle data points (number 2 in first index of fval) for each region.
    fval(2,nstack+2)=fval(3,nstack)
    fval(2,nstack+1)=fval(2,nstack)
    fval(2,nstack)=fval(1,nstack)
! Now do the integral over the right-most region in two different ways-
! a three point integral (valu) over each of the three subregions
! and a more accurate nine-point integral (value) over whole region.
! valu is used only for the error estimate.
    icount=0
    value=0.d0
    valu=0.d0
    do jj=nstack,nstack+2
      del1=xleft(jj+1)-xleft(jj)
!      print*, 'fn call 2'
      fval(1,jj)=fn(xleft(jj)+dx(1)*del1)
      fval(3,jj)=fn(xleft(jj)+dx(3)*del1)
!      print*, 'fn call 2'
      numcal = numcal + 2
      do kk=1,3
        icount=icount+1
        value=value+wt9(icount)*fval(kk,jj)*del
        valu=valu+fval(kk,jj)*wt(kk)*del1
      enddo
    enddo
    dif=abs(value-valu)
! If the following condition is true, add in this integral to the total,
! and reduce the number of regions under consideration.
    frac = del / (xmax - xmin)
    atsing = .false.
    if(frac .le. 1.0e-8) atsing = .true.
    if(dif .le. abr*frac .or. dif.le.rlr*abs(value) .or.           &
&       (atsing .and. (frac .le. 1.0e-15 .or. dif .le. abr*0.1  ))) then
      grater3=grater3+value
      error=error+abs(dif)
      nstack=nstack-1
! If no more regions, we are done.
      if(nstack.le.0) return
    else
! If the integration is insufficiently accurate, make each of the
! three subregions of the right-most region into regions.
! On next pass the right-most of these is the new current region.
      nstack=nstack+2
      maxns = max(maxns,nstack)
    endif
  enddo

end function grater3


