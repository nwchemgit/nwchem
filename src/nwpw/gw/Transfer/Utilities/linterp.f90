function linterp(data,grid,ndata,coord)
! For a set value of gridpoints GRID on the x-axis, and known function values 
! f(x) stored in DATA, find the linear interpolation of the value f(coord).
implicit none
integer :: ndata
double precision :: linterp,data(ndata),grid(ndata),coord
integer :: iihi,iilo,iimid,ii

if (coord.lt.grid(1).or.coord.gt.grid(ndata)) then
  write(6,*) 'input coordinate out of range in function linterp'
  stop
endif
iihi=ndata
iilo=1
do
  if (iihi-iilo.gt.1) then
    iimid=(iihi+iilo)/2
    if (grid(ndata).ge.grid(1)) then
      if (coord.ge.grid(iimid)) then
        iilo=iimid
      else
        iihi=iimid
      endif
    else 
      if (coord.le.grid(iimid)) then
        iilo=iimid
      else
        iihi=iimid
      endif
    endif
  else
    exit
  endif
enddo
if (coord.eq.grid(1)) then
  ii=iilo
elseif (coord.eq.grid(ndata)) then
  ii=iihi-1
else
  ii=iimid
endif
linterp=data(ii)+(data(ii+1)-data(ii))*(coord-grid(ii))/(grid(ii+1)-grid(ii))

return
end function linterp
