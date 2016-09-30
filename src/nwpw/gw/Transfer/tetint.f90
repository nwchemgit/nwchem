! Tetrahedron integration

! Two sets of routines, one dealing with the contribution from a single energy plane, and one dealing with the contribution between two energy planes.
! The latter set avoids certain artifacts that can occur when you have near but not quite degenerate energies, and are evaluating at an energy point between those two near-degenerate energies.
! The two-energy-plane sets start with "v" (for volumetric integration).
! Both sets use subroutines fval and fpol.


!***************************************************************************

! fval finds the energies at the corners of a cube given the values on a ngkpt(1) x ngkpt(2) x ngkpt(3) grid of energies.

subroutine fval(omega,ikk,ikkp,ngkpt,val)
implicit none
integer :: ikk(3),ikkp(3),ngkpt(3),ii
double precision :: omega(ngkpt(1),ngkpt(2),ngkpt(3))
double precision :: val(8)

do ii=1,3
  ikkp(ii)=ikk(ii)+1
  if (ikkp(ii).gt.ngkpt(ii)) ikkp(ii)=ikkp(ii)-ngkpt(ii)
enddo

val(1)=omega( ikk(1), ikk(2), ikk(3))
val(2)=omega(ikkp(1), ikk(2), ikk(3))
val(3)=omega( ikk(1),ikkp(2), ikk(3))
val(4)=omega(ikkp(1),ikkp(2), ikk(3))
val(5)=omega( ikk(1), ikk(2),ikkp(3))
val(6)=omega(ikkp(1), ikk(2),ikkp(3))
val(7)=omega( ikk(1),ikkp(2),ikkp(3))
val(8)=omega(ikkp(1),ikkp(2),ikkp(3))

return
end subroutine fval

!***************************************************************************

! fpol finds the complex function values at the corners of a cube given the values on a ngkpt(1) x ngkpt(2) x ngkpt(3) grid of energies.
! This is the function which will be integrated over.

subroutine fpol(poli,ngkpt,ikk,ikkp,polv)
implicit none
integer :: ikk(3),ikkp(3),ngkpt(3)
double complex :: poli(ngkpt(1),ngkpt(2),ngkpt(3))
double complex :: polv(8)

polv(1)=poli( ikk(1), ikk(2), ikk(3))
polv(2)=poli(ikkp(1), ikk(2), ikk(3))
polv(3)=poli( ikk(1),ikkp(2), ikk(3))
polv(4)=poli(ikkp(1),ikkp(2), ikk(3))
polv(5)=poli( ikk(1), ikk(2),ikkp(3))
polv(6)=poli(ikkp(1), ikk(2),ikkp(3))
polv(7)=poli( ikk(1),ikkp(2),ikkp(3))
polv(8)=poli(ikkp(1),ikkp(2),ikkp(3))

return
end subroutine fpol

!***************************************************************************

! The single-surface integration over a cube.  The cube is broken into six tetrahedrons, and each is integrated over via subroutine tetint.

subroutine cubeint(val,ww,polv,xint)
implicit none
double precision :: val(8),ww
double complex :: polv(8)
double complex :: xint,xintval,fpyr(4)
double precision :: rr(3,8),rrpyr(3,4),valpyr(4)
integer :: ivndx(4,6)
integer :: itet,iv
data rr(1:3,1) /0.0d0,0.0d0,0.0d0/
data rr(1:3,2) /1.0d0,0.0d0,0.0d0/
data rr(1:3,3) /0.0d0,1.0d0,0.0d0/
data rr(1:3,4) /1.0d0,1.0d0,0.0d0/
data rr(1:3,5) /0.0d0,0.0d0,1.0d0/
data rr(1:3,6) /1.0d0,0.0d0,1.0d0/
data rr(1:3,7) /0.0d0,1.0d0,1.0d0/
data rr(1:3,8) /1.0d0,1.0d0,1.0d0/
data ivndx(1:4,1) /1,2,3,5/
data ivndx(1:4,2) /3,5,6,7/
data ivndx(1:4,3) /2,3,5,6/
data ivndx(1:4,4) /2,3,4,6/
data ivndx(1:4,5) /3,4,6,7/
data ivndx(1:4,6) /4,6,7,8/
save rr,ivndx

xint=0.d0
do itet=1,6
!do itet=1,1
  do iv=1,4
    valpyr(iv)=val(ivndx(iv,itet))  ! energies at tetrahedron corners
    fpyr(iv)=polv(ivndx(iv,itet))   ! function values at tetrahedron corners
    rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))  ! coordinates of tetrahedron corners
  enddo
  call tetint(rrpyr,valpyr,ww,fpyr,xintval)
! debug
!write(27,'(i5,3x,2f19.6)') itet,xintval
  xint=xint+xintval
enddo
return
end subroutine cubeint

!***************************************************************************

subroutine tetint(kvtx,evtx,engy,fvtx,cint)
! integrate sint0=int_tetrahedron d^3k A(k) delta(e(k)-engy)
! where A(k)=fvtx(indxe(1))+fgrad dot k
! see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972)
implicit none
double precision :: kvtx(3,4),evtx(4),engy,xk(3,4)
double complex :: fvtx(4),cint
integer :: indxe(4),ii
double precision :: vol,sint0,sint1(3)
double complex :: fgrad(3)

call indxhpsort(4,4,evtx,indxe)
call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
xk(1:3,1)=(/0.d0,0.d0,0.d0/)
call mksint(kvtx(1:3,indxe),xk,vol,evtx(indxe),engy,sint0,sint1)
cint=fvtx(indxe(1))*sint0+fgrad(1)*sint1(1)+fgrad(2)*sint1(2)+fgrad(3)*sint1(3)

!write(6,*) "Area"
!write(6,'(4x,2e14.6)') sint0
!write(6,*)
!write(6,*) "First moment"
!do ii=1,3
!  write(6,'(i4,2e14.6)') ii,sint1(ii)
!enddo
!write(6,*)
!write(6,*) "Function value"
!write(6,'(4x,2e14.6)') fvtx(indxe(1))
!write(6,*)
!write(6,*) "Function gradient"
!do ii=1,3
!  write(6,'(i4,2e14.6)') ii,fgrad(ii)
!enddo
!write(6,*)
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(1)),evtx(indxe(1))-engy,fvtx(indxe(1))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(2)),evtx(indxe(2))-engy,fvtx(indxe(2))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(3)),evtx(indxe(3))-engy,fvtx(indxe(3))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(4)),evtx(indxe(4))-engy,fvtx(indxe(4))
!write(6,'(4i4)') indxe
!write(6,'(a,f10.6)') 'vol    ',vol
!!write(6,'(a,es14.6)') 'sint0  ',sint0
! debug
!write(27,'(a,f14.6)') 'sint0  ',sint0
!write(27,'(a,3f14.6)') 'sint1 ',sint1
!write(27,'(a,"(",f10.6,",",f10.6,")")') 'fvtx(1) ',fvtx(indxe(1))
!write(27,'(a,3("(",f10.6,",",f10.6,")"))') 'fgrad ',fgrad
!write(6,'(a,f14.6)') 'sint0*fvtx(1)  ',sint0*fvtx(indxe(1))
!write(6,'(a,f14.6)') 'sint1*fgrad  ',dot_product(sint1,fgrad)
!write(6,'(a,"(",f14.6,",",f14.6,")")') 'cint  ',cint
!write(6,'(3f10.6)') sint1/sint0
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,sint1/sint0)
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(2))-kvtx(1:3,indxe(1)))
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(3))-kvtx(1:3,indxe(1)))
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(4))-kvtx(1:3,indxe(1)))

return
end subroutine tetint

!***************************************************************************

! The volumetric two-surface integration over a cube.  The cube is broken into six tetrahedrons, and each is integrated over via subroutine tetint.

subroutine vcubeint(val,ww,dw,polv,xint)
implicit none
double precision :: val(8),ww,dw
double complex :: polv(8)
double complex :: xint,xintval,fpyr(4)
double precision :: rr(3,8),rrpyr(3,4),valpyr(4)
integer :: ivndx(4,6)
integer :: itet,iv
data rr(1:3,1) /0.0d0,0.0d0,0.0d0/
data rr(1:3,2) /1.0d0,0.0d0,0.0d0/
data rr(1:3,3) /0.0d0,1.0d0,0.0d0/
data rr(1:3,4) /1.0d0,1.0d0,0.0d0/
data rr(1:3,5) /0.0d0,0.0d0,1.0d0/
data rr(1:3,6) /1.0d0,0.0d0,1.0d0/
data rr(1:3,7) /0.0d0,1.0d0,1.0d0/
data rr(1:3,8) /1.0d0,1.0d0,1.0d0/
data ivndx(1:4,1) /1,2,3,5/
data ivndx(1:4,2) /3,5,6,7/
data ivndx(1:4,3) /2,3,5,6/
data ivndx(1:4,4) /2,3,4,6/
data ivndx(1:4,5) /3,4,6,7/
data ivndx(1:4,6) /4,6,7,8/
save rr,ivndx

xint=0.d0
do itet=1,6
!do itet=1,1
  do iv=1,4
    valpyr(iv)=val(ivndx(iv,itet))  ! energies at tetrahedron corners
    fpyr(iv)=polv(ivndx(iv,itet))   ! function values at tetrahedron corners
    rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))  ! coordinates of tetrahedron corners
  enddo
  call vtetint(rrpyr,valpyr,ww,dw,fpyr,xintval)
! debug
!write(27,'(i5,3x,2f19.6)') itet,xintval
  xint=xint+xintval
enddo
return
end subroutine vcubeint

!***************************************************************************

subroutine vtetint(kvtx,evtx,engy,de,fvtx,cint)
! volumetric tetrahedron integration, find volume between engy and engy+dw
! extension of tetint to case where need [int dw tetint(kvtx,evtx,w,fvtx,cint)]
! helps avoid artifacts when tetrahedron vertex energies are very closely spaced and engy just happens to lie between them
! replace area of plane of constant energy with volume of truncated pyramid ofallowed energies, center of mass of plane with center of mass of truncated pyramid.
! otherwise, as notes for tetint
implicit none
double precision :: kvtx(3,4),evtx(4),engy,de,xk(3,4)
double complex :: fvtx(4),cint
integer :: indxe(4),ii
double precision :: vol,sint0,sint1(3)
double complex :: fgrad(3)

call indxhpsort(4,4,evtx,indxe)
call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
xk(1:3,1)=(/0.d0,0.d0,0.d0/)
call mkvsint(kvtx(1:3,indxe),xk,vol,evtx(indxe),engy,de,sint0,sint1)
cint=fvtx(indxe(1))*sint0+fgrad(1)*sint1(1)+fgrad(2)*sint1(2)+fgrad(3)*sint1(3)
! debug
!write(27,'(a,f14.6)') 'sint0  ',sint0/de
!write(27,'(a,3f14.6)') 'sint1 ',sint1(1)/de,sint1(2)/de,sint1(3)/de
!write(27,'(a,f14.6)') 'sint0  ',sint0
!write(27,'(a,3f14.6)') 'sint1 ',sint1
!write(27,'(a,"(",f10.6,",",f10.6,")")') 'fvtx(1) ',fvtx(indxe(1))
!write(27,'(a,3("(",f10.6,",",f10.6,")"))') 'fgrad ',fgrad
!write(27,'(a,"(",f10.6,",",f10.6,")")') 'sint0*fvtx(1) ',sint0*fvtx(indxe(1))
!write(27,'(a,3("(",f10.6,",",f10.6,")"))') 'sint1*fgrad ',sint1(1)*fgrad(1),sint1(2)*fgrad(2),sint1(3)*fgrad(3)

return
end subroutine vtetint

!***************************************************************************

subroutine rvtetint(kvtx,evtx,engy,de,fvtx,cint)
! volumetric tetrahedron integration of real-valued function fvtx
! find volume between engy and engy+dw
! extension of tetint to case where need [int dw tetint(kvtx,evtx,w,fvtx,cint)]
! helps avoid artifacts when tetrahedron vertex energies are very closely spaced and engy just happens to lie between them
! replace area of plane of constant energy with volume of truncated pyramid ofallowed energies, center of mass of plane with center of mass of truncated pyramid.
! otherwise, as notes for tetint
implicit none
double precision :: kvtx(3,4),evtx(4),engy,de,xk(3,4)
double precision :: fvtx(4),cint
integer :: indxe(4),ii
double precision :: vol,sint0,sint1(3)
double precision :: fgrad(3)

!do ii=1,4
!write(6,*) kvtx(:,ii)
!!write(6,*) evtx(ii)
!enddo
!write(6,*)

call indxhpsort(4,4,evtx,indxe)
call mkgradr(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
xk(1:3,1)=(/0.d0,0.d0,0.d0/)
call mkvsint(kvtx(1:3,indxe),xk,vol,evtx(indxe),engy,de,sint0,sint1)
cint=fvtx(indxe(1))*sint0+fgrad(1)*sint1(1)+fgrad(2)*sint1(2)+fgrad(3)*sint1(3)
! debug
!write(27,'(a,f14.6)') 'sint0  ',sint0/de
!write(27,'(a,3f14.6)') 'sint1 ',sint1(1)/de,sint1(2)/de,sint1(3)/de
!write(6,'(a,f14.6)') 'sint0  ',sint0
!write(6,'(a,3f14.6)') 'sint1 ',sint1
!write(6,'(a,f10.6)') 'fvtx(1) ',fvtx(indxe(1))
!write(6,'(a,3f10.6)') 'fgrad ',fgrad
!write(27,'(a,"(",f10.6,",",f10.6,")")') 'sint0*fvtx(1) ',sint0*fvtx(indxe(1))
!write(27,'(a,3("(",f10.6,",",f10.6,")"))') 'sint1*fgrad ',sint1(1)*fgrad(1),sint1(2)*fgrad(2),sint1(3)*fgrad(3)

return
end subroutine rvtetint

!***************************************************************************

subroutine vnitetint(kvtx,evtx,engy,de,ftens,bmet,ivtx,cint)
! volumetric tetrahedron integration for non-isotropic function at one vertex
! ftens is tensor response at vertex ivtx
! response near vertex ivtx is khat * ftens * khat for khat the direction vector, and * represents the dot product.
! Dot product is evaluated with metric tensor bmet; A * B = sum_{ii,jj} A(ii) bmet(ii,jj) B(jj), to allow for reduced coordinate description (in cartesian coordinates, pass unit matrix for bmet).
! due to linearity, can treat only this vertex with function values at other vertices set to 0,
! then do normal tetrahedral integration over other vertices with value at ivtx set to 0
implicit none
double precision :: kvtx(3,4),evtx(4),engy,de,xk(3,4),bmet(3,3)
double complex :: ftens(3,3),cint,twvec(3),tuvec(3)
integer :: ivtx
integer :: indxe(4),ii,jj,kk,ll
double precision :: vol,rr(3,3),fside(3),khat(3,3),mag,tvvec(3)
double precision :: engy2, ehi, elo, de1
double precision :: dewl1,dewl2,dewl3,dewl4
double precision :: dev21, dev31, dev41, dev32, dev42, dev43
double precision :: qq(3,6),qp(3,6),dq(3,6)
double precision :: denom1,denom2,term1,term2,term3
double precision :: ll000,ll001,ll002,ll003,ll033,ll102,ll103,ll112,ll122,ll123,ll203,ll330,ll331,ll332,ll333

call indxhpsort(4,4,evtx,indxe)

engy2 = engy + de
xk(1:3,indxe(1))=(/0.d0,0.d0,0.d0/)
do ii=1,3
  do jj=1,3
    xk(jj,indxe(ii+1))=kvtx(jj,indxe(ii+1))-kvtx(jj,indxe(1))
  enddo
enddo
do jj=1,3
  call cross(xk(1:3,indxe(mod(jj,3)+2)),xk(1:3,indxe(mod(jj+1,3)+2)),rr(1:3,indxe(jj)))
enddo
vol=dot_product(xk(1:3,indxe(2)),rr(1:3,indxe(1)))
rr=rr/vol
vol=abs(vol)
do ii=1,3
  if (ii.lt.ivtx) then 
    jj=ii
  else
    jj=ii+1
  endif
  khat(:,ii)=xk(:,jj)-xk(:,ivtx)
  mag=sqrt(dot_product(khat(:,ii),khat(:,ii)))
  khat(:,ii)=khat(:,ii)/mag
enddo
do ii=1,3
  tvvec(:) = bmet(:,1)*khat(1,ii) + bmet(:,2)*khat(2,ii) + bmet(:,3)*khat(3,ii)
  twvec(:) = ftens(:,1)*tvvec(1) + ftens(:,2)*tvvec(2) + ftens(:,3)*tvvec(3)
  tuvec(:) = bmet(:,1)*twvec(1) + bmet(:,2)*twvec(2) + bmet(:,3)*twvec(3)
  fside(ii) = khat(1,ii)*tuvec(1) + khat(2,ii)*tuvec(2) + khat(3,ii)*tuvec(3)
enddo

cint = (0.d0,0.d0)

dev21 = evtx(indxe(2))-evtx(indxe(1))
dev31 = evtx(indxe(3))-evtx(indxe(1))
dev41 = evtx(indxe(4))-evtx(indxe(1))
dev32 = evtx(indxe(3))-evtx(indxe(2))
dev42 = evtx(indxe(4))-evtx(indxe(2))
dev43 = evtx(indxe(4))-evtx(indxe(3))

if (engy2.le.evtx(indxe(1)) .or. engy.ge.evtx(indxe(4))) then
! debug
!write(6,*) "n 0"
  return
else
  if ((engy.lt.evtx(indxe(2))).and.(evtx(indxe(2)).gt.evtx(indxe(1)))) then
! debug
!write(6,*) "p 1"
    if (engy.lt.evtx(indxe(1))) then
      elo = evtx(indxe(1))
      if (engy2.lt.evtx(indxe(2))) then
        ehi = engy2
      else
        ehi = evtx(indxe(2))
      endif
      de1 = ehi - elo
    else
      elo = engy
      if (engy2.gt.evtx(indxe(2))) then
        ehi = evtx(indxe(2))
        de1 = ehi - elo
      else
        de1 = de
      endif
    endif
    denom1=6.d0*dev21*dev31*dev41
    if (ivtx.eq.1) then
      term2=evtx(indxe(1))**2+2.d0*evtx(indxe(1))*evtx(indxe(2))
      term1=2.d0*evtx(indxe(1))+evtx(indxe(2))
      ll001=de1*(elo-evtx(indxe(1)))**2*(evtx(indxe(2))-elo)                  &
&          +de1**2*(elo*term1 - 1.5d0*elo**2 - 0.5d0*term2)                   &
&          +de1**3*(term1/3.d0 - elo)                                         &
&          -de1**4/4.d0
      term2=evtx(indxe(1))**2+2.d0*evtx(indxe(1))*evtx(indxe(3))
      term1=2.d0*evtx(indxe(1))+evtx(indxe(3))
      ll002=de1*(elo-evtx(indxe(1)))**2*(evtx(indxe(3))-elo)                  &
&          +de1**2*(elo*term1 - 1.5d0*elo**2 - 0.5d0*term2)                   &
&          +de1**3*(term1/3.d0 - elo)                                         &
&          -de1**4/4.d0
      term2=evtx(indxe(1))**2+2.d0*evtx(indxe(1))*evtx(indxe(4))
      term1=2.d0*evtx(indxe(1))+evtx(indxe(4))
      ll003=de1*(elo-evtx(indxe(1)))**2*(evtx(indxe(4))-elo)                  &
&          +de1**2*(elo*term1 - 1.5d0*elo**2 - 0.5d0*term2)                   &
&          +de1**3*(term1/3.d0 - elo)                                         &
&          -de1**4/4.d0
      cint=cint+(vol/denom1)*(ll001*fside(1)/dev21+ll002*fside(2)/dev31+ll003*fside(3)/dev41)
    elseif (ivtx.eq.2) then
      term1=elo-evtx(indxe(1))
      ll000=de1*term1**3+1.5d0*de1**2*term1**2+de1**3*term1+de1**4/4.d0
      cint=cint+(vol/denom1)*ll000*fside(1)/dev21
    elseif (ivtx.eq.3) then
      term1=elo-evtx(indxe(1))
      ll000=de1*term1**3+1.5d0*de1**2*term1**2+de1**3*term1+de1**4/4.d0
      cint=cint+(vol/denom1)*ll000*fside(1)/dev31
    else ! ivtx=4
      term1=elo-evtx(indxe(1))
      ll000=de1*term1**3+1.5d0*de1**2*term1**2+de1**3*term1+de1**4/4.d0
      cint=cint+(vol/denom1)*ll000*fside(1)/dev41
    endif
  endif
  if ((engy.lt.evtx(indxe(3))).and. &
&     (engy2.gt.evtx(indxe(2))).and. &
&     (evtx(indxe(3)).gt.evtx(indxe(2)))) then
! debug
!write(6,*) "p 2"
    denom1 = 6*dev41*dev42*dev31
    denom2 = 6*dev32*dev42*dev31
    if (engy.lt.evtx(indxe(2))) then
      elo = evtx(indxe(2))
      if (engy2.lt.evtx(indxe(3))) then
        ehi = engy2
      else
        ehi = evtx(indxe(3))
      endif
      de1 = ehi - elo
    else
      elo = engy
      if (engy2.gt.evtx(indxe(3))) then
        ehi = evtx(indxe(3))
        de1 = ehi - elo
      else
        de1 = de
      endif
    endif
    dewl1 = elo - evtx(indxe(1))
    dewl2 = elo - evtx(indxe(2))
    dewl3 = elo - evtx(indxe(3))
    dewl4 = elo - evtx(indxe(4))
    if (ivtx.eq.1.or.ivtx.eq.2) then
      ll122=de1*dewl2*dewl3*dewl3                                             &
&          +de1**2*(dewl2*dewl3*2.d0+dewl3**2)/2.d0                           &
&          +de1**3*(dewl2+dewl3+dewl3)/3.d0                                   &
&          +de1**4/4.d0
      ll033=de1*dewl1*dewl4*dewl4                                             &
&          +de1**2*(dewl1*dewl4*2.d0+dewl4**2)/2.d0                           &
&          +de1**3*(dewl1+dewl4+dewl4)/3.d0                                   &
&          +de1**4/4.d0
      if (ivtx.eq.1) then
        ll203=de1*dewl1*dewl3*dewl4                                           &
&            +de1**2*(dewl1*dewl3+dewl1*dewl4+dewl3*dewl4)/2.d0               &
&            +de1**3*(dewl1+dewl3+dewl4)/3.d0                                 &
&            +de1**4/4.d0
        cint=cint+vol*((ll203*fside(2)/dev31+ll033*fside(3)/dev41)/denom1+ll122*fside(2)/(denom2*dev31))
      else ! ivtx.eq.2
        ll123=de1*dewl2*dewl3*dewl4                                           &
&            +de1**2*(dewl2*dewl3+dewl2*dewl4+dewl3*dewl4)/2.d0               &
&            +de1**3*(dewl2+dewl3+dewl4)/3.d0                                 &
&            +de1**4/4.d0
        cint=cint+vol*(ll033*fside(3)/(denom1*dev42)+(ll122*fside(2)/dev32+ll123*fside(3)/dev42)/denom2)
      endif
    else ! ivtx.eq.3.or.ivtx.eq.4
      ll003=de1*dewl1*dewl1*dewl4                                             &
&          +de1**2*(dewl1**2+dewl1*dewl4*2.d0)/2.d0                           &
&          +de1**3*(dewl1+dewl1+dewl4)/3.d0                                   &
&          +de1**4/4.d0
      ll112=de1*dewl2*dewl2*dewl3                           &
&        +de1**2*(dewl2**2+dewl2*dewl3*2.d0)/2.d0           &
&        +de1**3*(dewl2+dewl2+dewl3)/3.d0                   &
&        +de1**4/4.d0
      if (ivtx.eq.3) then
        ll102=de1*dewl1*dewl2*dewl3                                           &
&            +de1**2*(dewl1*dewl2+dewl1*dewl3+dewl2*dewl3)/2.d0               &
&            +de1**3*(dewl1+dewl2+dewl3)/3.d0                                 &
&            +de1**4/4.d0
        cint=cint-vol*(ll003*fside(1)/(denom1*dev31)+(ll102*fside(1)/dev31+ll112*fside(2)/dev32)/denom2)
      else ! ivtx.eq.4
        ll103=de1*dewl1*dewl2*dewl4                                           &
&            +de1**2*(dewl1*dewl2+dewl1*dewl4+dewl2*dewl4)/2.d0               &
&            +de1**3*(dewl1+dewl2+dewl4)/3.d0                                 &
&            +de1**4/4.d0
        cint=cint-vol*((ll003*fside(1)/dev41+ll103*fside(2)/dev42)/denom1+ll112*fside(2)/(denom2*dev42))
      endif
    endif
  endif
  if ((engy2.gt.evtx(indxe(3))).and.(evtx(indxe(4)).gt.evtx(indxe(3)))) then
    denom1 = 6*dev41*dev42*dev43
! debug
!write(6,*) "p 3"
    if (engy2.gt.evtx(indxe(4))) then
      ehi = evtx(indxe(4))
      if (engy.lt.evtx(indxe(3))) then
        elo = evtx(indxe(3))
      else
        elo = engy
      endif
      de1 = ehi - elo
    else 
      ehi = engy2
      if (engy.lt.evtx(indxe(3))) then
        elo = evtx(indxe(3))
        de1 = ehi - elo
      else
        elo = engy
        de1 = de
      endif
    endif
    dewl4 = elo - evtx(indxe(4))
    if (ivtx.eq.1.or.ivtx.eq.2.or.ivtx.eq.3) then
      term1=evtx(indxe(4))-elo
      ll333=de1*dewl4**3+1.5d0*de1**2*dewl4**2+de1**3*dewl4+de1**4/4.d0
      if (ivtx.eq.1) then
        cint=cint-(vol/denom1)*ll333*fside(1)/dev41
      elseif (ivtx.eq.2) then
        cint=cint-(vol/denom1)*ll333*fside(2)/dev42
      else ! ivtx.eq.3
        cint=cint-(vol/denom1)*ll333*fside(3)/dev43
      endif
    else ! ivtx.eq.4
      dewl1 = elo - evtx(indxe(1))
      dewl2 = elo - evtx(indxe(2))
      dewl3 = elo - evtx(indxe(3))
      ll330=de1*dewl1*dewl4*dewl4                                 &
&          +de1**2*(dewl1*dewl4*2.d0+dewl4**2)/2.d0               &
&          +de1**3*(dewl1+2.d0*dewl4)/3.d0                        &
&          +de1**4/4.d0
      ll331=de1*dewl2*dewl4*dewl4                                 &
&          +de1**2*(dewl2*dewl4*2.d0+dewl4**2)/2.d0               &
&          +de1**3*(dewl2+2.d0*dewl4)/3.d0                        &
&          +de1**4/4.d0
      ll332=de1*dewl3*dewl4*dewl4                                 &
&          +de1**2*(dewl3*dewl4*2.d0+dewl4**2)/2.d0               &
&          +de1**3*(dewl3+2.d0*dewl4)/3.d0                        &
&          +de1**4/4.d0
      cint=cint+(vol/denom1)*(ll330*fside(1)/dev41+ll331*fside(2)/dev42+ll332*fside(3)/dev43)
    endif
  endif
endif
return
end subroutine vnitetint

!***************************************************************************

subroutine mkvsint(kvtx,xk,vol,evtx,engy,de,sint0,sint1)
! integrate sint0=int_engy^(engy+de) dw int_tetrahedron d^3k delta(e(k)-w)
! and       sint1=int_engy^(engy+de) dw int_tetrahedron d^3k k delta(e(k)-w)
! where e(k)=evtx(1)+grad(e) dot k
! vertices are assumed ordered so evtx(1)<=evtx(2)<=evtx(3)<=evtx(4)
implicit none
double precision :: kvtx(3,4),evtx(4),engy,de,vol,xk(3,4)
double precision :: vol_t
double precision :: sint0,sint1(3)
integer :: indxe(4),ii,jj
double precision :: engy2, ehi, elo, de1
double precision :: dewl1,dewl2,dewl3,dewl4,dewh1,dewh2,dewh3,dewh4
double precision :: dev21, dev31, dev41, dev32, dev42, dev43
double precision :: qq(3,6),qp(3,6),dq(3,6)
double precision :: denom1,denom2,term1,term2
double precision :: vv1,vv2,vv3,vv4,vv5,vv6
! debug
!double precision :: qtest(3),qtest2(3)

vol_t = vol/6.d0

engy2 = engy + de
sint0=0.d0
do ii=1,3
  sint1(ii)=0.d0
enddo

dev21 = evtx(2)-evtx(1)
dev31 = evtx(3)-evtx(1)
dev41 = evtx(4)-evtx(1)
dev32 = evtx(3)-evtx(2)
dev42 = evtx(4)-evtx(2)
dev43 = evtx(4)-evtx(3)

if (engy2.le.evtx(1) .or. engy.ge.evtx(4)) then
! debug
!write(6,*) "n 0"
  return
else if (engy.le.evtx(1) .and. engy2.ge.evtx(4)) then
! debug
!write(6,*) "a 0"
  sint0 = vol_t
  sint1(:) = 0.25*vol_t*(xk(:,2) + xk(:,3) + xk(:,4))
  return
else
  if ((engy.lt.evtx(2)).and.(evtx(2).gt.evtx(1))) then
    denom1 = dev21*dev31*dev41
    if (engy.lt.evtx(1)) then
      elo = evtx(1)
      if (engy2.lt.evtx(2)) then
! debug
!write(6,*) "p 1 : 0 1"
        ehi = engy2
      else
! debug
!write(6,*) "p 1 : 0 0"
        ehi = evtx(2)
      endif
      de1 = ehi - elo
      sint0 = de1**3*vol_t/denom1
    else
      elo = engy
      dewl1 = elo - evtx(1)
      if (engy2.gt.evtx(2)) then
! debug
!write(6,*) "p 1 : 1 0"
        ehi = evtx(2)
        de1 = ehi - elo
      else
! debug
!write(6,*) "p 1 : 1 1"
        ehi = engy2
        de1 = de
      endif
      sint0 = (3*de1*dewl1**2 + 3*de1**2*dewl1 + de1**3)*vol_t/denom1
      term1 = 0.25d0*de1/dev21
      dq(:,1) = term1*xk(:,2)
      term1 = 0.25d0*de1/dev31
      dq(:,2) = term1*xk(:,3)
      term1 = 0.25d0*de1/dev41
      dq(:,3) = term1*xk(:,4)
      vv1 = dewl1**3*vol_t/denom1
      sint1(:) = vv1*(dq(:,1)+dq(:,2)+dq(:,3))
    endif
    dewh1 = ehi - evtx(1)
    term1 = 0.25d0*dewh1/dev21
    qp(:,1) = term1*xk(:,2)
    term1 = 0.25d0*dewh1/dev31
    qp(:,2) = term1*xk(:,3)
    term1 = 0.25d0*dewh1/dev41
    qp(:,3) = term1*xk(:,4)
    sint1(:) = sint1(:) + sint0*(qp(:,1)+qp(:,2)+qp(:,3))
  endif
  if (((engy.gt.evtx(2).and.engy.lt.evtx(3)).or. &
&     (engy2.gt.evtx(2).and.engy2.lt.evtx(3)).or. &
&     (engy.le.evtx(2).and.engy2.ge.evtx(3))).and.(evtx(3).gt.evtx(2))) then
    denom1 = dev41*dev32*dev42
    denom2 = dev31*dev41*dev32
    if (engy.lt.evtx(2)) then
      elo = evtx(2)
      if (engy2.lt.evtx(3)) then
! debug
!write(6,*) "p 2 : 0 1"
        ehi = engy2
      else
! debug
!write(6,*) "p 2 : 0 0"
        ehi = evtx(3)
      endif
      de1 = ehi - elo
    else
      elo = engy
      if (engy2.gt.evtx(3)) then
! debug
!write(6,*) "p 2 : 1 0"
        ehi = evtx(3)
        de1 = ehi - elo
      else
! debug
!write(6,*) "p 2 : 1 1"
        ehi = engy2
        de1 = de
      endif
    endif
    dewl1 = elo - evtx(1)
    dewl2 = elo - evtx(2)
    dewl3 = evtx(3) - elo
    dewl4 = evtx(4) - elo
    dewh1 = ehi - evtx(1)
    dewh2 = ehi - evtx(2)
    dewh3 = evtx(3) - ehi
    dewh4 = evtx(4) - ehi
    vv4 = de1*dewl3*dewl1*vol_t/denom2
    vv5 = abs((de1/dev32)*(dewh1*dewh4/(dev41*dev42) &
&                     +dewh2/dev42-dewh1/dev41)*vol_t)
    vv6 = abs((de1/dev32)*(dewh1*dewl1/(dev41*dev31) &
&                     -(dewh1/dev32)*(dewl2/dev31+dewl3/dev41))*vol_t)
    sint0 = sint0 + vv4 + vv5 + vv6
    term1 = 0.25d0*dewl1/dev31
    qq(:,2) = term1*xk(:,3)
    term1 = 0.25d0*dewl1/dev41
    qq(:,3) = term1*xk(:,4)
    term1 = 0.25d0*dewl2/dev32
    term2 = 0.25d0*dewl3/dev32
    qq(:,4) = term1*xk(:,3) + term2*xk(:,2)
    term1 = 0.25d0*dewh1/dev41
    qp(:,3) = term1*xk(:,4)
    term1 = 0.25d0*dewh2/dev32
    term2 = 0.25d0*dewh3/dev32
    qp(:,4) = term1*xk(:,3) + term2*xk(:,2)
    term1 = 0.25d0*dewh2/dev42
    term2 = 0.25d0*dewh4/dev42
    qp(:,5) = term1*xk(:,4) + term2*xk(:,2)
    sint1(:) = sint1(:) + vv4*(qq(:,3)+qq(:,2)+qq(:,4)+qp(:,3)) &
&                       + vv5*(qq(:,4)+qp(:,5)+qp(:,4)+qp(:,3)) &
&                       + vv6*(qq(:,4)+qq(:,2)+qp(:,3)+qp(:,4))
    if (engy.gt.evtx(2)) then
      vv1 = abs((de1/dev42)*(dewl1*dewl2/(dev41*dev32)-dewl2/dev32)*vol_t)
      vv3 = de1*dewh4*dewl2*vol_t/denom1
      sint0 = sint0 + vv1 + vv3
      term1 = 0.25d0*dewl2/dev42
      term2 = 0.25d0*dewl4/dev42
      qq(:,5) = term1*xk(:,4) + term2*xk(:,2)
      sint1(:) = sint1(:) + vv1*(qq(:,3)+qq(:,4)+qq(:,5)+qp(:,5)) &
&                         + vv3*(qq(:,3)+qq(:,4)+qp(:,5)+qp(:,3)) 
    endif
    if (engy2.lt.evtx(3)) then
      vv2 = de1*dewh3*dewh1*vol_t/denom2
      sint0 = sint0 + vv2
      term1 = 0.25d0*dewh1/dev31
      qp(:,2) = term1*xk(:,3)
      sint1(:) = sint1(:) + vv2*(qq(:,2)+qp(:,2)+qp(:,4)+qp(:,3)) 
    endif
! debug
!qq(:,2)=qq(:,2)*4.d0;
!qq(:,3)=qq(:,3)*4.d0;
!qq(:,4)=qq(:,4)*4.d0;
!qq(:,5)=qq(:,5)*4.d0;
!qp(:,2)=qp(:,2)*4.d0;
!qp(:,3)=qp(:,3)*4.d0;
!qp(:,4)=qp(:,4)*4.d0;
!qp(:,5)=qp(:,5)*4.d0;
!call cross(qp(:,3)-qq(:,4),qq(:,2)-qq(:,4),qtest)
!qtest2=qp(:,4)-qq(:,4)
!write(27,*) dot_product(qtest2,qtest)/6.d0
!write(27,*) vv6
! debug
  endif
  if ((engy2.gt.evtx(3)).and.(evtx(4).gt.evtx(3))) then
    denom1 = dev41*dev42*dev43
    if (engy2.gt.evtx(4)) then
      ehi = evtx(4)
      if (engy.lt.evtx(3)) then
! debug
!write(6,*) "p 3 : 0 1"
        elo = evtx(3)
      else
! debug
!write(6,*) "p 3 : 0 0"
        elo = engy
      endif
      de1 = ehi - elo
      vv1 = de1**3*vol_t/denom1
    else 
      ehi = engy2
      dewh4 = evtx(4) - ehi
      if (engy.lt.evtx(3)) then
! debug
!write(6,*) "p 3 : 1 0"
        elo = evtx(3)
        de1 = ehi - elo
      else
! debug
!write(6,*) "p 3 : 1 1"
        elo = engy
        de1 = de
      endif
      vv1 = (3*de1*dewh4**2  + 3*de1**2*dewh4 + de1**3)*vol_t/denom1
      term1 = 0.25d0*de1/dev41
      dq(:,3) = term1*xk(:,4)
      term1 = 0.25d0*de1/dev42
      dq(:,5) = term1*(xk(:,4)-xk(:,2))
      term1 = 0.25d0*de1/dev43
      dq(:,6) = term1*(xk(:,4)-xk(:,3))
      dewh4 = evtx(4) - ehi
      vv2 = dewh4**3*vol_t/denom1
      sint1(:) = sint1(:) - vv2*(dq(:,3)+dq(:,5)+dq(:,6))
    endif
    sint0 = sint0 + vv1
    dewl1 = elo - evtx(1)
    dewl2 = elo - evtx(2)
    dewl3 = elo - evtx(3)
    dewl4 = evtx(4) - elo
    term1 = 0.25d0*dewl1/dev41
    qq(:,3) = term1*xk(:,4)
    term1 = 0.25d0*dewl2/dev42
    term2 = 0.25d0*dewl4/dev42
    qq(:,5) = term1*xk(:,4) + term2*xk(:,2)
    term1 = 0.25d0*dewl3/dev43
    term2 = 0.25d0*dewl4/dev43
    qq(:,6) = term1*xk(:,4) + term2*xk(:,3)
    sint1(:) = sint1(:) + vv1*(qq(:,3)+qq(:,5)+qq(:,6)+0.25*xk(:,4))
  endif
endif

return
end subroutine mkvsint

!***************************************************************************

subroutine mksint(kvtx,xk,vol,evtx,engy,sint0,sint1)
! integrate sint0=int_tetrahedron d^3k delta(e(k)-engy)
! and       sint1=int_tetrahedron d^3k k delta(e(k)-engy)
! where e(k)=evtx(1)+grad(e) dot k
! vertices are assumed ordered so evtx(1)<=evtx(2)<=evtx(3)<=evtx(4)
! see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972)
implicit none
double precision :: kvtx(3,4),evtx(4),engy,vol,xk(3,4)
double precision :: sint0,sint1(3)
integer :: indxe(4),ii,jj
double precision :: f0b,f1b,f2b,f3b,ss(3),s0(3),s1(3),s2(3),s3(3),xfact,xfact0,xfact1,xfact2,xfact3
double precision :: de43,de21
double precision :: xkt(3,3),avec(3),xdum(3,3),cm(3)
double precision :: bgrad(3)

if (engy.le.evtx(1)) then
  sint0=0.d0
  do ii=1,3
    sint1(ii)=0.d0
  enddo
elseif (engy.le.evtx(2)) then
! debug
!write(27,*) "p 1"
  sint0=vol*(engy-evtx(1))**2 &
&       /(2.d0*(evtx(2)-evtx(1))*(evtx(3)-evtx(1))*(evtx(4)-evtx(1)))
  xfact=(engy-evtx(1))/3.d0
  do ii=1,3
    ss(ii)=0.d0
    do jj=2,4
      ss(ii)=ss(ii)+xfact*xk(ii,jj)/(evtx(jj)-evtx(1))
!      ss(ii)=ss(ii)+xfact*(kvtx(ii,jj)-kvtx(ii,1))/(evtx(jj)-evtx(1))
    enddo
    sint1(ii)=sint0*ss(ii)
  enddo
elseif (engy.lt.evtx(3)) then
! debug
!write(27,*) "p 2"
  de21=evtx(2)-evtx(1)
  de43=evtx(4)-evtx(3)
  if (de21.lt.1.d-5*(evtx(4)-evtx(1)).and.de43.lt.1.d-5*(evtx(4)-evtx(1))) then
    xkt(:,1)=xk(:,3)*(engy-evtx(1))/(evtx(3)-evtx(1))
    xkt(:,2)=xk(:,4)*(engy-evtx(1))/(evtx(4)-evtx(1))
    xkt(:,3)=(xk(:,4)-xk(:,2))*(engy-evtx(2))/(evtx(4)-evtx(2))+xk(:,2)
    call cross(xkt(:,1)-xkt(:,2),xkt(:,3)-xkt(:,2),avec)
    call mkgradr(kvtx,evtx,vol,xdum,bgrad)
    sint0=sqrt(dot_product(avec,avec)/dot_product(bgrad,bgrad))
    cm=(xkt(1:3,1)+xkt(1:3,3))/2.d0
    sint1=cm*sint0
  elseif (de21.gt.de43) then
    f0b=vol*(engy-evtx(1))**2 &
&       /(2.d0*de21*(evtx(3)-evtx(1))*(evtx(4)-evtx(1)))
    f1b=vol*(engy-evtx(2))**2 &
&       /(2.d0*de21*(evtx(3)-evtx(2))*(evtx(4)-evtx(2)))
    sint0=f0b-f1b
    xfact0=(engy-evtx(1))/3.d0
    xfact1=(engy-evtx(2))/3.d0
    do ii=1,3
      s0(ii)=0.d0
      do jj=2,4
        s0(ii)=s0(ii)+xfact0*xk(ii,jj)/(evtx(jj)-evtx(1))
!        s0(ii)=s0(ii)+xfact0*(kvtx(ii,jj)-kvtx(ii,1))/(evtx(jj)-evtx(1))
      enddo
      s1(ii)=xk(ii,2)
!      s1(ii)=kvtx(ii,2)-kvtx(ii,1)
      do jj=1,4
        if (jj.eq.2) cycle
        s1(ii)=s1(ii)+xfact1*(xk(ii,jj)-xk(ii,2))/(evtx(jj)-evtx(2))
!        s1(ii)=s1(ii)+xfact1*(kvtx(ii,jj)-kvtx(ii,2))/(evtx(jj)-evtx(2))
      enddo
      sint1(ii)=s0(ii)*f0b-s1(ii)*f1b
    enddo
  else
    f3b=vol*(engy-evtx(4))**2 &
&         /(2.d0*(evtx(4)-evtx(1))*(evtx(4)-evtx(2))*de43)
    f2b=vol*(engy-evtx(3))**2 &
&         /(2.d0*(evtx(3)-evtx(1))*(evtx(3)-evtx(2))*de43)
    sint0=f3b-f2b
    xfact3=(engy-evtx(4))/3.d0
    xfact2=(engy-evtx(3))/3.d0
    do ii=1,3
      s3(ii)=xk(ii,4)
!      s3(ii)=kvtx(ii,4)-kvtx(ii,1)
      do jj=1,3
        s3(ii)=s3(ii)+xfact3*(xk(ii,jj)-xk(ii,4))/(evtx(jj)-evtx(4))
!        s3(ii)=s3(ii)+xfact3*(kvtx(ii,jj)-kvtx(ii,4))/(evtx(jj)-evtx(4))
      enddo
      s2(ii)=xk(ii,3)
!      s2(ii)=kvtx(ii,3)-kvtx(ii,1)
      do jj=1,4
        if (jj.eq.3) cycle
        s2(ii)=s2(ii)+xfact2*(xk(ii,jj)-xk(ii,3))/(evtx(jj)-evtx(3))
!        s2(ii)=s2(ii)+xfact2*(kvtx(ii,jj)-kvtx(ii,3))/(evtx(jj)-evtx(3))
      enddo
      sint1(ii)=s3(ii)*f3b-s2(ii)*f2b
    enddo
  endif
elseif (engy.le.evtx(4)) then
! debug
!write(27,*) "p 3"
  sint0=vol*(engy-evtx(4))**2 &
&         /(2.d0*(evtx(4)-evtx(1))*(evtx(4)-evtx(2))*(evtx(4)-evtx(3)))
  xfact=(engy-evtx(4))/3.d0
  do ii=1,3
    ss(ii)=xk(ii,4)
!    ss(ii)=kvtx(ii,4)-kvtx(ii,1)
    do jj=1,3
      ss(ii)=ss(ii)+xfact*(xk(ii,jj)-xk(ii,4))/(evtx(jj)-evtx(4))
!      ss(ii)=ss(ii)+xfact*(kvtx(ii,jj)-kvtx(ii,4))/(evtx(jj)-evtx(4))
    enddo
    sint1(ii)=sint0*ss(ii)
  enddo
else
  sint0=0.d0
  do ii=1,3
    sint1(ii)=0.d0
  enddo
endif

return
end subroutine mksint

!***************************************************************************

! Find the gradient of a real valued function given the values at the corners of a tetrahedron

subroutine mkgradr(kvtx,valvtx,vol,xk,grad)
implicit none
double precision :: kvtx(3,4)
double precision :: valvtx(4)
double precision :: grad(3)
double precision :: rr(3,3),xk(3,4),vol
double precision :: valdiff
double precision :: dot
external dot
integer ii,jj

do ii=1,3
  do jj=1,3
    xk(jj,ii)=kvtx(jj,ii+1)-kvtx(jj,1)
  enddo
enddo
do jj=1,3
  call cross(xk(1:3,mod(jj,3)+1),xk(1:3,mod(jj+1,3)+1),rr(1:3,jj))
enddo
vol=dot_product(xk(1:3,1),rr(1:3,1))
rr=rr/vol
vol=abs(vol)
do ii=1,3
  grad(ii)=0.d0
enddo
do jj=1,3
  valdiff=(valvtx(jj+1)-valvtx(1))
  do ii=1,3
    grad(ii)=grad(ii)+valdiff*rr(ii,jj)
  enddo
enddo

return
end subroutine mkgradr

!***************************************************************************

! Find the gradient of a complex valued function given the values at the corners of a tetrahedron

subroutine mkgrad(kvtx,valvtx,vol,xk,grad)
implicit none
double precision :: kvtx(3,4)
double complex :: valvtx(4)
double complex :: grad(3)
double precision :: rr(3,3),xk(3,4),vol
double complex :: valdiff
double precision :: dot
external dot
integer ii,jj

do ii=1,3
  do jj=1,3
    xk(jj,ii)=kvtx(jj,ii+1)-kvtx(jj,1)
  enddo
enddo
do jj=1,3
  call cross(xk(1:3,mod(jj,3)+1),xk(1:3,mod(jj+1,3)+1),rr(1:3,jj))
enddo
vol=dot_product(xk(1:3,1),rr(1:3,1))
rr=rr/vol
vol=abs(vol)
do ii=1,3
  grad(ii)=0.d0
enddo
do jj=1,3
  valdiff=(valvtx(jj+1)-valvtx(1))
  do ii=1,3
    grad(ii)=grad(ii)+valdiff*rr(ii,jj)
  enddo
enddo

return
end subroutine mkgrad

!***************************************************************************

! vector cross product

subroutine cross(v1,v2,v3)
! compute v3 = v1 x v2
implicit none
double precision :: v1(3),v2(3)
double precision :: v3(3)
v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
return
end subroutine cross

!***************************************************************************

! Tetrahedron integration over cube between wlo and whi

subroutine cubeint_d(val,valp,ww,polv,wlo,whi,xint)
implicit none
double precision :: val(8),valp(8),ww,wlo,whi
double complex :: polv(8)
double complex :: xint,xintval1,xintval2,fpyr(4)
double precision :: rr(3,8),rrpyr(3,4),valpyr(4),valpyrp(4)
integer :: ivndx(4,6)
integer :: itet,iv
integer :: full1,full2
data rr(1:3,1) /0.0d0,0.0d0,0.0d0/
data rr(1:3,2) /1.0d0,0.0d0,0.0d0/
data rr(1:3,3) /0.0d0,1.0d0,0.0d0/
data rr(1:3,4) /1.0d0,1.0d0,0.0d0/
data rr(1:3,5) /0.0d0,0.0d0,1.0d0/
data rr(1:3,6) /1.0d0,0.0d0,1.0d0/
data rr(1:3,7) /0.0d0,1.0d0,1.0d0/
data rr(1:3,8) /1.0d0,1.0d0,1.0d0/
data ivndx(1:4,1) /1,2,3,5/
data ivndx(1:4,2) /3,5,6,7/
data ivndx(1:4,3) /2,3,5,6/
data ivndx(1:4,4) /2,3,4,6/
data ivndx(1:4,5) /3,4,6,7/
data ivndx(1:4,6) /4,6,7,8/
save rr,ivndx

xint=(0.d0,0.d0)
do itet=1,6
!do itet=1,1
  do iv=1,4
    valpyr(iv)=val(ivndx(iv,itet))
    valpyrp(iv)=valp(ivndx(iv,itet))
    fpyr(iv)=polv(ivndx(iv,itet))
    rrpyr(1:3,iv)=rr(1:3,ivndx(iv,itet))
  enddo
  call tetint_f(rrpyr,valpyr,valpyrp,ww,fpyr,whi,xintval1,full1)
  call tetint_f(rrpyr,valpyr,valpyrp,ww,fpyr,wlo,xintval2,full2)
!write(6,*) itet,full1,full2,dble(xintval1),dble(xintval2)
!write(6,*) itet,xintval2-xintval1
  if ((full1.ne.full2).or.(full1.eq.0)) xint=xint+xintval1-xintval2
!write(6,*) itet,xint
!write(6,*)
enddo
return
end subroutine cubeint_d

!***************************************************************************

! Tetrahedron integration with a Fermi cut-off

subroutine tetint_f(kvtx,evtx,evtx2,engy,fvtx,wcut,cint,full)
! integrate sint0=int_tetrahedron d^3k A(k) delta(e(k)-engy) theta(e2(k)-wcut)
! where A(k)=fvtx(indxe(1))+fgrad dot k
! see G. Lehmann and M. Taut, Phys. stat. sol. (b) 54 469 (1972)
implicit none
double precision :: kvtx(3,4),kvtxp(3,4),kvtxm(3,4),kvtxn(3,4)
double precision :: xk(3,4),xkp(3,4),xkm(3,4),xkn(3,4),xkcond(3,4)
double precision :: evtx(4),evtx2(4),evtxp(4),evtxm(4),evtxn(4),wmax,wmin,wcond(4)
double precision :: engy,wcut
double complex :: fvtx(4),fvtxp(4),fvtxm(4),fvtxn(4),cint,tint
integer :: indxe(4),indxe2(4),indxep(4),indxem(4),indxen(4),ii,jj
double precision :: vol,volp,sint0,sint1(3),tint0,tint1(3),rr(3)
double complex :: fgrad(3), fgradp(3)
double precision :: bgrad(3)
double precision :: delta
integer :: full
integer :: ncint ! number of intersections of conduction plane with tetrahedron

call indxhpsort(4,4,evtx,indxe)
call indxhpsort(4,4,evtx2,indxe2)
xk(1:3,1)=(/0.d0,0.d0,0.d0/)
xkp(1:3,1)=(/0.d0,0.d0,0.d0/)
!write(6,'(f10.6,3x,4f10.6)') wcut,evtx2(indxe2)
!write(6,'(f10.6,3x,4f10.6)') engy,evtx(indxe)
if (wcut.ge.evtx2(indxe2(4))) then
!write(6,*) 4
  call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
  call mksint(kvtx(1:3,indxe),xk,vol,evtx(indxe),engy,sint0,sint1)
  cint=fvtx(indxe(1))*sint0+fgrad(1)*sint1(1)+fgrad(2)*sint1(2)+fgrad(3)*sint1(3)
  full=1
elseif (wcut.ge.evtx2(indxe2(3))) then
!write(6,*) 3
! values at intersections of cut-off plane with tetrahedron
  call mkgradr(kvtx(1:3,indxe),evtx(indxe),vol,xk(1:3,2:4),bgrad)
!write(6,'(4f10.6)') evtx(indxe)
  xkp(:,2)=(kvtx(:,indxe2(2))-kvtx(:,indxe2(4)))*(wcut-evtx2(indxe2(4)))/(evtx2(indxe2(2))-evtx2(indxe2(4)))
  xkp(:,3)=(kvtx(:,indxe2(3))-kvtx(:,indxe2(4)))*(wcut-evtx2(indxe2(4)))/(evtx2(indxe2(3))-evtx2(indxe2(4)))
  xkp(:,4)=(kvtx(:,indxe2(1))-kvtx(:,indxe2(4)))*(wcut-evtx2(indxe2(4)))/(evtx2(indxe2(1))-evtx2(indxe2(4)))
  do ii=1,4
    kvtxp(:,ii)=xkp(:,ii)+kvtx(:,indxe2(4))
  enddo
  evtxp(1)=evtx(indxe2(4))
  evtxp(2)=evtx(indxe(1))+dot_product(bgrad,kvtxp(:,2)-kvtx(:,indxe(1)))
  evtxp(3)=evtx(indxe(1))+dot_product(bgrad,kvtxp(:,3)-kvtx(:,indxe(1)))
  evtxp(4)=evtx(indxe(1))+dot_product(bgrad,kvtxp(:,4)-kvtx(:,indxe(1)))
  call indxhpsort(4,4,evtxp,indxep)
! find energies at intersections of conduction plane with tetrahedron
! wmin,wmax,wtest: valence energies at intersection points
  ncint=0
  wmax=evtx2(indxe2(1))
  wmin=evtx2(indxe2(4))
  do ii=1,3
  do jj=ii+1,4
    if ((engy.ge.evtx(indxe(ii))).and.(engy.le.evtx(indxe(jj)))) then
      ncint=ncint+1
      delta=(engy-evtx(indxe(ii)))/(evtx(indxe(jj))-evtx(indxe(ii)))
      wcond(ncint)=evtx2(indxe(ii))+(evtx2(indxe(jj))-evtx2(indxe(ii)))*delta
      xkcond(:,ncint)=kvtx(:,indxe(ii))+(kvtx(:,indxe(jj))-kvtx(:,indxe(ii)))*delta
      wmax=max(wmax,wcond(ncint))
      wmin=min(wmin,wcond(ncint))
    endif
  enddo
  enddo
!write(6,*) "Boundary vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,2f10.6)') ii,kvtx(1:3,indxe2(ii)),evtx2(indxe2(ii)),evtx(indxe2(ii))
!enddo
!write(6,*)
!write(6,*) "Valence intersection energy",wcut
!write(6,*) "Valence intersection plane"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,kvtxp(1:3,indxe2(ii)),evtxp(indxe2(ii))
!enddo
!write(6,*)
!write(6,'(a,f10.6)') "Conduction intersection energy",engy
!write(6,*) "Conduction intersection vertices"
!do ii=1,ncint
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,xkcond(:,ii),wcond(ii)
!enddo
!write(6,*)
  if (wmax.le.wcut) then
!write(6,*) "3a"
! cut-off plane above plane of energy integration - full intergal
    call tetint(kvtx,evtx,engy,fvtx,cint)
    full=1
  elseif (wmin.ge.wcut) then
!write(6,*) "3b"
! cut-off plane below plane of energy integration - no contribution
    cint=(0.d0,0.d0)
    full=-1
  else
!write(6,*) "3c"
! cut-off plane intersects plane of energy integration - partial contribution
! integral over entire tetrahedron
    call tetint(kvtx,evtx,engy,fvtx,cint)
! subtract off integral above cut-off plane
    call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
    fvtxp(1)=fvtx(indxe2(4))
    fvtxp(2)=fvtx(indxe(1))+dot_product(fgrad,(kvtxp(:,2)-kvtx(:,indxe(1))))
    fvtxp(3)=fvtx(indxe(1))+dot_product(fgrad,(kvtxp(:,3)-kvtx(:,indxe(1))))
    fvtxp(4)=fvtx(indxe(1))+dot_product(fgrad,(kvtxp(:,4)-kvtx(:,indxe(1))))
    call tetint(kvtxp,evtxp,engy,fvtxp,tint)
    cint=cint-tint
    full=0
  endif
elseif (wcut.gt.evtx2(indxe2(2))) then
!write(6,*) 2
  call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
  call mkgradr(kvtx(1:3,indxe),evtx(indxe),vol,xk(1:3,2:4),bgrad)
! kvtxm, evtxm values for cut-off plane crossing tetrahedron
  kvtxm(:,4)=(kvtx(:,indxe2(4))-kvtx(:,indxe2(1)))*(wcut-evtx2(indxe2(1)))/(evtx2(indxe2(4))-evtx2(indxe2(1)))+kvtx(:,indxe2(1))
  kvtxm(:,2)=(kvtx(:,indxe2(2))-kvtx(:,indxe2(4)))*(wcut-evtx2(indxe2(4)))/(evtx2(indxe2(2))-evtx2(indxe2(4)))+kvtx(:,indxe2(4))
  kvtxm(:,3)=(kvtx(:,indxe2(3))-kvtx(:,indxe2(1)))*(wcut-evtx2(indxe2(1)))/(evtx2(indxe2(3))-evtx2(indxe2(1)))+kvtx(:,indxe2(1))
  kvtxm(:,1)=(kvtx(:,indxe2(3))-kvtx(:,indxe2(2)))*(wcut-evtx2(indxe2(2)))/(evtx2(indxe2(3))-evtx2(indxe2(2)))+kvtx(:,indxe2(2))
  do ii=1,4
    evtxm(ii)=evtx(indxe(1))+dot_product(bgrad,kvtxm(:,ii)-kvtx(:,indxe(1)))
    fvtxm(ii)=fvtx(indxe(1))+dot_product(fgrad,kvtxm(:,ii)-kvtx(:,indxe(1)))
  enddo
  call indxhpsort(4,4,evtxm,indxem)
  ncint=0
  wmax=evtx2(indxe2(1))
  wmin=evtx2(indxe2(4))
  do ii=1,3
  do jj=ii+1,4
    if ((engy.ge.evtx(indxe(ii))).and.(engy.le.evtx(indxe(jj)))) then
      ncint=ncint+1
      delta=(engy-evtx(indxe(ii)))/(evtx(indxe(jj))-evtx(indxe(ii)))
      wcond(ncint)=evtx2(indxe(ii))+(evtx2(indxe(jj))-evtx2(indxe(ii)))*delta
      xkcond(:,ncint)=kvtx(:,indxe(ii))+(kvtx(:,indxe(jj))-kvtx(:,indxe(ii)))*delta
      wmax=max(wmax,wcond(ncint))
      wmin=min(wmin,wcond(ncint))
    endif
  enddo
  enddo
!write(6,*) wcut
!write(6,*) "Function values"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtx(indxe2(ii))
!enddo
!write(6,*)
!write(6,*) "Boundary vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,2f10.6)') ii,kvtx(1:3,indxe2(ii)),evtx2(indxe2(ii)),evtx(indxe2(ii))
!enddo
!write(6,*)
!write(6,*) "Function values at valence intersection"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtxm(ii)
!enddo
!write(6,*)
!write(6,*) "Valence intersection plane"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,2f10.6)') ii,kvtxm(1:3,ii),evtxm(ii)
!enddo
!write(6,*)
!write(6,'(a,f10.6)') "Conduction intersection vertices",engy
!write(6,*) "Conduction intersection vertices"
!do ii=2,4
!  write(6,'(2i4,3f10.6)') 1,ii,( (kvtx(:,indxe(ii))-kvtx(:,indxe(1)))*(engy-evtx(indxe(1)))/(evtx(indxe(ii))-evtx(indxe(1))) &
!& +kvtx(:,indxe(1)))
!enddo
!do ii=3,4
!  write(6,'(2i4,3f10.6)') 2,ii,( (kvtx(:,indxe(ii))-kvtx(:,indxe(2)))*(engy-evtx(indxe(2)))/(evtx(indxe(ii))-evtx(indxe(2))) &
!& +kvtx(:,indxe(2)))
!enddo
!write(6,'(2i4,3f10.6)') 3,4,( (kvtx(:,indxe(3))-kvtx(:,indxe(4)))*(engy-evtx(indxe(4)))/(evtx(indxe(3))-evtx(indxe(4))) &
!& +kvtx(:,indxe(4)))
  if (wmax.le.wcut) then
!write(6,*) "2a"
! cut-off plane above plane of energy integration - full intergal
    call tetint(kvtx,evtx,engy,fvtx,cint)
    full=1
  elseif (wmin.ge.wcut) then
!write(6,*) "2b"
    cint=(0.d0,0.d0)
    full=-1
  else
!write(6,*) "2c"
! break tetrahedron in two along plane 'n', integrate lower half-tetrahedron first
    do ii=1,3
      kvtxn(:,ii)=kvtx(:,indxe2(ii))
      evtxn(ii)=evtx(indxe2(ii))
      fvtxn(ii)=fvtx(indxe2(ii))
    enddo
    kvtxn(:,4)=kvtxm(:,4)
    evtxn(4)=evtxm(4)
    fvtxn(4)=fvtxm(4)
!write(6,*) "INTEGRATION 1"
!write(6,*) "Function values at integration tetrahedron"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtxn(ii)
!enddo
!write(6,*)
!write(6,*) "Energy of integration plane"
!write(6,'(f10.6)') engy
!write(6,*)
!write(6,*) "Integration tetrahedron vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,kvtxn(1:3,ii),evtxn(ii)
!enddo
!write(6,*)
    call tetint(kvtxn,evtxn,engy,fvtxn,cint)
!write(6,*) "integration result",cint
! subtract off part above cut-off energy plane
    kvtxp(:,1)=kvtx(:,indxe2(3))
    kvtxp(:,2)=kvtxm(:,1)
    kvtxp(:,3)=kvtxm(:,3)
    kvtxp(:,4)=kvtxm(:,4)
    evtxp(1)=evtx(indxe2(3))
    evtxp(2)=evtxm(1)
    evtxp(3)=evtxm(3)
    evtxp(4)=evtxm(4)
    fvtxp(1)=fvtx(indxe2(3))
    fvtxp(2)=fvtxm(1)
    fvtxp(3)=fvtxm(3)
    fvtxp(4)=fvtxm(4)
!write(6,*) "INTEGRATION 2"
!write(6,*) "Function values at integration tetrahedron"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtxp(ii)
!enddo
!write(6,*)
!write(6,*) "Energy of integration plane"
!write(6,'(f10.6)') engy
!write(6,*)
!write(6,*) "Integration tetrahedron vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,kvtxp(1:3,ii),evtxp(ii)
!enddo
!write(6,*)
    call tetint(kvtxp,evtxp,engy,fvtxp,tint)
!write(6,*) "integration result",tint
    cint=cint-tint
! add part in upper half-tetrahedron
    kvtxp(:,1)=kvtx(:,indxe2(2))
    kvtxp(:,2)=kvtxm(:,2)
    kvtxp(:,3)=kvtxm(:,1)
    kvtxp(:,4)=kvtxm(:,4)
    evtxp(1)=evtx(indxe2(2))
    evtxp(2)=evtxm(2)
    evtxp(3)=evtxm(1)
    evtxp(4)=evtxm(4)
    fvtxp(1)=fvtx(indxe2(2))
    fvtxp(2)=fvtxm(2)
    fvtxp(3)=fvtxm(1)
    fvtxp(4)=fvtxm(4)
!write(6,*) "INTEGRATION 3"
!write(6,*) "Function values at integration tetrahedron"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtxp(ii)
!enddo
!write(6,*)
!write(6,*) "Energy of integration plane"
!write(6,'(f10.6)') engy
!write(6,*)
!write(6,*) "Integration tetrahedron vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,kvtxp(1:3,ii),evtxp(ii)
!enddo
!write(6,*)
    call tetint(kvtxp,evtxp,engy,fvtxp,tint)
!write(6,*) "integration result",tint
    cint=cint+tint
    full=0
  endif
elseif (wcut.gt.evtx2(indxe2(1))) then
!write(6,*) 1
  call mkgrad(kvtx(1:3,indxe),fvtx(indxe),vol,xk(1:3,2:4),fgrad)
  call mkgradr(kvtx(1:3,indxe),evtx(indxe),vol,xk(1:3,2:4),bgrad)
  xkp(:,2)=(kvtx(:,indxe2(2))-kvtx(:,indxe2(1)))*(wcut-evtx2(indxe2(1)))/(evtx2(indxe2(2))-evtx2(indxe2(1)))
  xkp(:,3)=(kvtx(:,indxe2(3))-kvtx(:,indxe2(1)))*(wcut-evtx2(indxe2(1)))/(evtx2(indxe2(3))-evtx2(indxe2(1)))
  xkp(:,4)=(kvtx(:,indxe2(4))-kvtx(:,indxe2(1)))*(wcut-evtx2(indxe2(1)))/(evtx2(indxe2(4))-evtx2(indxe2(1)))
  do ii=1,4
    kvtxp(:,ii)=xkp(:,ii)+kvtx(:,indxe2(1))
  enddo
  evtxp(1)=evtx(indxe2(1))
  fvtxp(1)=fvtx(indxe2(1))
  do ii=2,4
    evtxp(ii)=evtx(indxe(1))+dot_product(bgrad,kvtxp(:,ii)-kvtx(:,indxe(1)))
    fvtxp(ii)=fvtx(indxe(1))+dot_product(fgrad,kvtxp(:,ii)-kvtx(:,indxe(1)))
  enddo
  call indxhpsort(4,4,evtxp,indxep)
  ncint=0
  wmax=evtx2(indxe2(1))
  wmin=evtx2(indxe2(4))
  do ii=1,3
  do jj=ii+1,4
    if ((engy.ge.evtx(indxe(ii))).and.(engy.le.evtx(indxe(jj)))) then
      ncint=ncint+1
      delta=(engy-evtx(indxe(ii)))/(evtx(indxe(jj))-evtx(indxe(ii)))
      wcond(ncint)=evtx2(indxe(ii))+(evtx2(indxe(jj))-evtx2(indxe(ii)))*delta
      xkcond(:,ncint)=kvtx(:,indxe(ii))+(kvtx(:,indxe(jj))-kvtx(:,indxe(ii)))*delta
      wmax=max(wmax,wcond(ncint))
      wmin=min(wmin,wcond(ncint))
    endif
  enddo
  enddo
!write(6,*) wcut
!write(6,*) "Function values"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtx(indxe2(ii))
!enddo
!write(6,*)
!write(6,*) "Boundary vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,2f10.6)') ii,kvtx(1:3,indxe2(ii)),evtx2(indxe2(ii)),evtx(indxe2(ii))
!enddo
!write(6,*) "Function values at integration plane"
!do ii=1,4
!  write(6,'(i4,2e14.6)') ii,fvtxp(ii)
!enddo
!write(6,*)
!write(6,*) "Energy of integration plane"
!write(6,'(f10.6)') engy
!write(6,*)
!write(6,*) "Valence intersection vertices"
!do ii=1,4
!  write(6,'(i4,3f10.6,3x,f10.6)') ii,kvtxp(1:3,ii),evtxp(ii)
!enddo
!write(6,*)
!write(6,*) "Conduction intersection vertices"
!do ii=2,4
!  write(6,'(2i4,3f10.6)') 1,ii,( (kvtx(:,indxe(ii))-kvtx(:,indxe(1)))*(engy-evtx(indxe(1)))/(evtx(indxe(ii))-evtx(indxe(1))) &
!& +kvtx(:,indxe(1)))
!enddo
!do ii=3,4
!  write(6,'(2i4,3f10.6)') 2,ii,( (kvtx(:,indxe(ii))-kvtx(:,indxe(2)))*(engy-evtx(indxe(2)))/(evtx(indxe(ii))-evtx(indxe(2))) &
!& +kvtx(:,indxe(2)))
!enddo
!write(6,'(2i4,3f10.6)') 3,4,( (kvtx(:,indxe(3))-kvtx(:,indxe(4)))*(engy-evtx(indxe(4)))/(evtx(indxe(3))-evtx(indxe(4))) &
!& +kvtx(:,indxe(4)))
!write(6,*)
  if (wmin.lt.wcut) then
    call tetint(kvtxp,evtxp,engy,fvtxp,cint)
    if (wmax.le.wcut) then
!write(6,*) '1a'
      full=1
    else
!write(6,*) '1b'
      full=0
    endif
  else
!write(6,*) '1c'
    cint=(0.d0,0.d0)
    full=-1
  endif
else
!write(6,*) 0
  cint=(0.d0,0.d0)
  full=-1
endif
!write(6,*) cint

!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(1)),evtx(indxe(1))-engy,fvtx(indxe(1))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(2)),evtx(indxe(2))-engy,fvtx(indxe(2))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(3)),evtx(indxe(3))-engy,fvtx(indxe(3))
!write(6,'(3f10.6,3x,f12.6,2x,"(",f10.6,",",f10.6,")")') kvtx(1:3,indxe(4)),evtx(indxe(4))-engy,fvtx(indxe(4))
!write(6,'(4i4)') indxe
!write(6,'(a,f10.6)') 'vol    ',vol
!write(6,'(a,f10.6)') 'sint0  ',sint0
!write(6,'(a,es14.6)') 'sint0  ',sint0
!write(6,'(a,"(",f10.6,",",f10.6,")")') 'fvtx(1) ',fvtx(indxe(1))
!write(6,'(a,3f10.6)') 'sint1 ',sint1
!write(6,'(a,3("(",f10.6,",",f10.6,")"))') 'fgrad ',fgrad
!write(6,'(a,"(",f10.6,",",f10.6,")")') 'cint  ',cint
!write(6,'(3f10.6)') sint1/sint0
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,sint1/sint0)
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(2))-kvtx(1:3,indxe(1)))
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(3))-kvtx(1:3,indxe(1)))
!write(6,'("(",f10.6,",",f10.6,")")') fvtx(indxe(1))+dot_product(fgrad,kvtx(1:3,indxe(4))-kvtx(1:3,indxe(1)))

return
end subroutine tetint_f

!***************************************************************************

subroutine qsquare(qq,bmet,qq2)
double precision :: qq(3),bmet(3,3),qq2
integer :: ii,jj

  qq2=0.d0
  do ii=1,3
  do jj=1,3
    qq2 = qq2 + qq(ii)*bmet(ii,jj)*qq(jj)
  enddo
  enddo

end subroutine qsquare

!***************************************************************************

subroutine singular_adaptive_tetint(qvtx,evtx,engy,de,ivtx,bmet,abr,rlr,rlr0,cint)
! Integrate q(x)q/q^4 over volume between two constant-energy planes
! over tetrahedron specified by vertices q = qvtx and energies evtx at vertices
! where (x) is the outer product and one vertex has qvtx=(0,0,0).
! ivtx identifies vertex with qvtx=(0,0,0).
! Sub-divides a tetrahedron containing a singular vertex into 8 sub-tetrahedra
! and sums over the tetrahedral integrals of the 7 non-singular sub-tetrahedra.
! Repeats process until converged.
! Vectors in reduced coordinates; for Cartesian coordinates set bmet to unit matrix
integer :: it,jt
double precision :: qvtx(3,4),evtx(4),engy,de,bmet(3,3),abr,rlr,rlr0
double precision :: cint(3,3)
integer :: ngkpt(3)
integer :: ivtx ! the index of the singular vertex
double precision :: qq(3),qq2
double precision :: qqval(3,9),qq2val(9) ! coordinates on sub-divided grid, and their squares
double precision :: eeval(9) ! energy at above coordinates
double precision :: fsum(3,3)
integer :: vindx(3) ! tells which vertices are non-singular
double precision :: rqvtx(3,4),revtx(4)
double precision :: rfvtx(4),rcint(3,3)
integer :: ii,jj,kk,iv,itet
double precision :: step
double precision :: eemin,eemax
logical :: converged
integer :: ivndx(4,7)
data ivndx(1:4,1) /1,4,7,8/
data ivndx(1:4,2) /2,5,7,9/
data ivndx(1:4,3) /3,6,8,9/
data ivndx(1:4,4) /4,5,6,7/
data ivndx(1:4,5) /4,6,7,8/
data ivndx(1:4,6) /5,6,7,9/
data ivndx(1:4,7) /6,7,8,9/

  do ii=1,3
  do jj=1,3
    cint(ii,jj) = (0.d0, 0.d0)
  enddo
  enddo

!write(6,*) "energies: ",engy,engy+de
  jj = 1
  do ii=1,4
    if (ii.ne.ivtx) then
      vindx(jj) = ii
      jj = jj+1
    endif
  enddo
!write(6,*) "ivtx",ivtx
!write(6,*) "vindx",vindx
  do ii=1,3
    qqval(:,ii) = qvtx(:,vindx(ii))
    call qsquare(qqval(:,ii),bmet,qq2val(ii))
    eeval(ii) = evtx(vindx(ii))
  enddo
!write(6,*) "qq2"
!write(6,*) qq2
!write(6,*) "bmet"
!do ii=1,3
!write(6,*) (bmet(ii,jj),jj=1,3)
!enddo
!write(6,*) "evtx"
!do ii=1,4
! write(6,*) evtx(ii)
!enddo
!write(6,*) "qvtx"
!do ii=1,4
! write(6,*) qvtx(:,ii)
!enddo
!read(*,*)
  step = 0.5d0
  do
    eeval(4) = step * evtx(vindx(1)) + (1.d0-step) * evtx(ivtx)
    eeval(5) = step * evtx(vindx(2)) + (1.d0-step) * evtx(ivtx)
    eeval(6) = step * evtx(vindx(3)) + (1.d0-step) * evtx(ivtx)
    eeval(7) = (eeval(1) + eeval(2))/2.d0
    eeval(8) = (eeval(1) + eeval(3))/2.d0
    eeval(9) = (eeval(2) + eeval(3))/2.d0
    qqval(:,4) = step * qvtx(:,vindx(1)) ! value of qq at ivtx = (0,0,0)
    qqval(:,5) = step * qvtx(:,vindx(2))
    qqval(:,6) = step * qvtx(:,vindx(3))
    qqval(:,7) = (qqval(:,1) + qqval(:,2))/2.d0
    qqval(:,8) = (qqval(:,1) + qqval(:,3))/2.d0
    qqval(:,9) = (qqval(:,2) + qqval(:,3))/2.d0
    do ii=4,9
      call qsquare(qqval(:,ii),bmet,qq2val(ii))
    enddo
!write(6,*) "qqval"
!do ii=1,9
!write(6,*) qqval(:,ii)
!enddo
!write(6,*)

    do it=1,3
    do jt=1,3
      fsum(it,jt) = 0.d0
    enddo
    enddo
    do itet=1,7
!write(6,*) "Tetrahedron ",itet
      do iv=1,4
        revtx(iv)=eeval(ivndx(iv,itet))  ! energies at tetrahedron corners
        rqvtx(:,iv)=qqval(:,ivndx(iv,itet))  ! coordinates of tetrahedron corners
!write(6,'(5x,3f10.5)') rqvtx(:,iv)
      enddo
!write(6,*)
      call q_adaptive_tetint(rqvtx,revtx,engy,de,bmet,abr,rlr0,rcint)
      fsum = fsum + rcint
!write(6,*)
!do ii=1,3
!write(6,*) (dble(rcint(ii,jj)),jj=1,3)
!enddo
!read(*,*)
!write(6,*) engy,de
!write(6,*) itet,rcint,fsum
!write(6,*)
    enddo
    cint = cint+fsum
!write(6,*)
!do ii=1,3
!write(6,*) (dble(cint(ii,jj)),jj=1,3)
!enddo
!read(*,*)
    eemin = min(eeval(4),eeval(5),eeval(6),eeval(ivtx))
    eemax = max(eeval(4),eeval(5),eeval(6),eeval(ivtx))
    if (engy.gt.eemax.or.engy+de.lt.eemin) then
       exit  ! energy slice does not exist in singular sub-tetrahedron - nothing more to do
    else 
      converged = .true.
      do it=1,3
      do jt=1,3
        if (abs(fsum(it,jt)).gt.(abs(rlr*cint(it,jt))+abr)) converged = .false.
      enddo
      enddo
      if (converged) exit ! converged
    endif

    eeval(1) = eeval(4)
    eeval(2) = eeval(5)
    eeval(3) = eeval(6)
    qqval(:,1) = qqval(:,4)
    qqval(:,2) = qqval(:,5)
    qqval(:,3) = qqval(:,6)
    step = step/2.d0
  enddo

end subroutine singular_adaptive_tetint

!***************************************************************************

subroutine q_adaptive_tetint(qvtx,evtx,engy,de,bmet,abr,rlr,cint)
! Integrate q(x)q/q^4 over volume between two constant-energy planes
! over tetrahedron specified by vertices q = qvtx and energies evtx at vertices
! where (x) is the outer product and q=(0,0,0) does not exist in or on tetrahedron.
! Sub-divides a tetrahedron containing a singular vertex into 8 sub-tetrahedra
! and sums over the tetrahedral integrals.
! Repeats process until converged.
! Vectors in reduced coordinates; for Cartesian coordinates set bmet to unit matrix
double precision :: qvtx(3,4),evtx(4),engy,de,bmet(3,3),abr,rlr
double precision :: cint(3,3)
integer :: ngkpt(3)
integer :: nstack, istack, mx
parameter (mx=500)
integer :: recursion_depth(mx), recursion_store
double precision :: qq(3),qq2
double precision :: qqsub(3,4),qq2sub(4) ! coordinates on sub-divided tetrahedron, and their squares
double precision :: eesub(4) ! energy at above coordinates
double precision :: qqgrid(3,0:9),qq2grid(0:9) ! coordinates on sub-division grid, and their squares
double precision :: eegrid(0:9) ! energy at above coordinates
double precision :: qqval(3,4,mx),qq2val(4,mx) ! coordinates of tetrahedron stack, and their squares
double precision :: eeval(4,mx) ! energy at above coordinates
double precision :: ffval(3,3,mx) ! integrated value of tetrahedra on stack
double precision :: ffstore(3,3), ffsum(3,3), ffmax
double precision :: fsum, cint_old
double precision :: rfvtx(4),rcint(3,3)
integer :: ii,jj,kk,iv,itet,it,jt
double precision :: step
double precision :: eemin,eemax
logical :: converged
integer :: ivndx(4,8),gridndx(2,4:9)
data ivndx(1:4,1) /1,4,7,8/
data ivndx(1:4,2) /2,5,7,9/
data ivndx(1:4,3) /3,6,8,9/
data ivndx(1:4,4) /4,5,6,7/
data ivndx(1:4,5) /4,6,7,8/
data ivndx(1:4,6) /5,6,7,9/
data ivndx(1:4,7) /6,7,8,9/
data ivndx(1:4,8) /4,5,6,0/
data gridndx(1:2,4) /1,0/
data gridndx(1:2,5) /2,0/
data gridndx(1:2,6) /3,0/
data gridndx(1:2,7) /1,2/
data gridndx(1:2,8) /1,3/
data gridndx(1:2,9) /2,3/

  do ii=1,3
  do jj=1,3
    cint(ii,jj) = 0.d0
  enddo
  enddo

!write(6,*) "energies: ",engy,engy+de
  nstack = 1

  qqval(:,:,1) = qvtx
  eeval(:,1) = evtx
!do iv=1,4
!write(6,*) qvtx(:,iv)
!enddo
  do iv=1,4
    call qsquare(qqval(:,iv,1),bmet,qq2val(iv,1))
  enddo
!write(6,*) "initial tetrahedron"
!do iv=1,4
!write(6,*) qqval(:,iv,1)
!enddo
!do iv=1,4
!write(6,*) qq2val(iv,1),eeval(iv,1)
!enddo
!write(6,*)
!write(6,*) engy,de
!write(6,*)
  do it=1,3
  do jt=it,3
    do iv=1,4
      rfvtx(iv) = qqval(it,iv,1)*qqval(jt,iv,1)/qq2val(iv,1)**2
    enddo
    call rvtetint(qqval(:,:,1),eeval(:,1),engy,de,rfvtx,ffval(it,jt,1))
!write(6,*) rfvtx
!write(6,*) ffval(it,jt,1)
  enddo
  enddo

  recursion_depth(1) = 0
  do 
    if (nstack+8.ge.mx) then
      write(*,*) "q_adaptive_tetint: TOO MANY REGIONS"
      stop
    endif
!write(6,*) "initial tetrahedron"
!do iv=1,4
!write(6,*) qqval(:,iv,nstack)
!enddo
!do iv=1,4
!write(6,*) qq2val(iv,nstack),eeval(iv,nstack)
!enddo
!read(*,*)
! divide last tetrahderon on stack into 8 sub-tetrahedra
    do ii=1,4
      jj = mod(ii,4)
      qqgrid(:,jj) = qqval(:,ii,nstack)
      call qsquare(qqgrid(:,jj),bmet,qq2grid(jj))
      eegrid(jj) = eeval(ii,nstack)
    enddo
    do ii=4,9
      qqgrid(:,ii) = (qqgrid(:,gridndx(1,ii))+qqgrid(:,gridndx(2,ii)))*0.5d0
      call qsquare(qqgrid(:,ii),bmet,qq2grid(ii))
      eegrid(ii) = (eegrid(gridndx(1,ii))+eegrid(gridndx(2,ii)))*0.5d0
    enddo
! store value of original tetrahedron
    ffstore = ffval(:,:,nstack)
    do it=1,3
    do jt=1,3
      ffsum(it,jt) = 0.d0
    enddo
    enddo
    recursion_store = recursion_depth(nstack)
! compute values of sub-tetrahedra
    do itet=1,8
      istack = nstack+itet-1
      recursion_depth(istack) = recursion_store+1
      do iv=1,4
        qqval(:,iv,istack) = qqgrid(:,ivndx(iv,itet))
        qq2val(iv,istack) = qq2grid(ivndx(iv,itet))
        eeval(iv,istack) = eegrid(ivndx(iv,itet))
      enddo
!write(6,*) "tetrahedron ",itet
!do iv=1,4
!write(6,*) qqval(:,iv,istack)
!enddo 
!do iv=1,4
!write(6,*) qq2val(iv,istack),eeval(iv,istack)
!enddo
      do it=1,3
      do jt=it,3
        do iv=1,4
          rfvtx(iv) = qqval(it,iv,istack)*qqval(jt,iv,istack)/qq2val(iv,istack)**2
        enddo
!do iv=1,4
!write(6,*) qqval(:,iv,istack)
!!write(6,*) eeval(iv,istack)
!enddo
!write(6,*)
        call rvtetint(qqval(:,:,istack),eeval(:,istack),engy,de,rfvtx,ffval(it,jt,istack))
!do iv=1,4
!write(6,*) qqval(:,iv,istack)
!!write(6,*) eeval(iv,istack)
!enddo
!write(6,*)
!read(*,*)
        ffsum(it,jt) = ffsum(it,jt) + ffval(it,jt,istack)
!write(6,*) rfvtx
!write(6,*) ffval(it,jt,istack)
      enddo
      enddo
!read(*,*)
    enddo
!write(6,*) nstack, recursion_depth(nstack)
!do it=1,3
!do jt=1,3
!write(6,*) cint(it,jt),ffsum(it,jt),ffstore(it,jt)
!enddo
!enddo
!read(*,*)
! compare sum over sub-tetrahedra to original value to see if converged
    ffmax = 0.d0
    do it=1,3
    do jt=it,3
      if (abs(ffsum(it,jt)).gt.ffmax) ffmax=abs(ffsum(it,jt))
    enddo
    enddo
    converged = .true.
    do it=1,3
    do jt=it,3
      if (abs(ffsum(it,jt)-ffstore(it,jt)).gt.(abs(rlr*ffmax)+abr)) converged = .false.
    enddo
    enddo
    if ((converged.or.recursion_depth(nstack).gt.5).and.recursion_depth(nstack).ge.2) then
! if converged, remove this tetrahedron from the stack and add its value to total integrated value
      nstack = nstack-1
      cint = cint + ffsum
    else 
! otherwise replace tetrahedron with its sub-tetrahedra, expanding the stack
      nstack_old = nstack
      nstack = nstack+7
    endif
    if (nstack.le.0) exit
  enddo

  do it=2,3
  do jt=1,it-1
    cint(it,jt) = cint(jt,it)
  enddo
  enddo
!write(6,*)
!do it=1,3
!do jt=1,3
!write(6,*) cint(it,jt)
!enddo
!enddo
!read(*,*)

end subroutine q_adaptive_tetint
