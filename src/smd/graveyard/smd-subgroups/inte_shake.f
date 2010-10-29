c
c $Id$
c

      SUBROUTINE inte_shake(tstep,ntcons,ekin)

      implicit none

      include 'p_array.inc'
      include 'p_input.inc'
      include 'cm_atom.inc'
      include 'cm_latt.inc'
      include 'cm_cons.inc'

      integer i,j,it,maxit,iatm,iatm1,iatm2,idist
      integer ntcons

      real*8  ncc,nvv,dcc
      real*8  force,rma,rmb
      real*8  ekin,tmpvx,tmpvy,tmpvz
      real*8  tstep,tstepsq,rtstep,tol,mxdiff
      real*8  nrij,orij,nrijsq,dijsq,diffsq,dotprod

      dimension ncc(mxatms,3),nvv(mxatms,3),dcc(mxatms,3)
      dimension nrij(mxcons2,3),orij(mxcons2,3)

      ekin=0.0
      tstepsq=tstep**2
      rtstep=1.0/tstep
      tol=1.0d-8
      maxit=100
      mxdiff = 0.0d0

      do i=1,ntcons

       iatm1=conlist(i,1)
       iatm2=conlist(i,2)

       orij(i,1)=ccc(iatm1,1)-ccc(iatm2,1)
       orij(i,2)=ccc(iatm1,2)-ccc(iatm2,2)
       orij(i,3)=ccc(iatm1,3)-ccc(iatm2,3)

      enddo

      call tool_rebox(ntcons,mxcons2,latt,rlatt,orij)

      do i=1,natms

       iatm=atmtype(i)

       ncc(i,1)=ccc(i,1)
       ncc(i,2)=ccc(i,2)
       ncc(i,3)=ccc(i,3)

       nvv(i,1)=vvv(i,1)+fff(i,1)*tstep/typmass(iatm)
       nvv(i,2)=vvv(i,2)+fff(i,2)*tstep/typmass(iatm)
       nvv(i,3)=vvv(i,3)+fff(i,3)*tstep/typmass(iatm)

       ccc(i,1)=ncc(i,1)+tstep*nvv(i,1)
       ccc(i,2)=ncc(i,2)+tstep*nvv(i,2)
       ccc(i,3)=ncc(i,3)+tstep*nvv(i,3)

      enddo 

      do i=1,maxit

       do j=1,ntcons

        iatm1=conlist(j,1)
        iatm2=conlist(j,2)

        nrij(j,1)=ccc(iatm1,1)-ccc(iatm2,1)
        nrij(j,2)=ccc(iatm1,2)-ccc(iatm2,2)
        nrij(j,3)=ccc(iatm1,3)-ccc(iatm2,3)

       enddo

       call tool_rebox(ntcons,mxcons2,latt,rlatt,nrij) 

       do j=1,natms
        dcc(j,1)=0.0
        dcc(j,2)=0.0
        dcc(j,3)=0.0
       enddo

       do j=1,ntcons

        iatm1=conlist(j,1)
        iatm2=conlist(j,2)
        idist=conlist(j,3)

        nrijsq=nrij(j,1)**2+nrij(j,2)**2+nrij(j,3)**2
        dijsq=consdist(idist)**2
        diffsq=dijsq-nrijsq
        mxdiff=max(mxdiff,abs(diffsq)/consdist(idist))

        dotprod=orij(j,1)*nrij(j,1)
     $         +orij(j,2)*nrij(j,2)
     $         +orij(j,3)*nrij(j,3)

        iatm=atmtype(iatm1)
        rma= tstepsq/typmass(iatm)
        iatm=atmtype(iatm2)
        rmb=-tstepsq/typmass(iatm)
        force=diffsq/(-2.0*(rma-rmb)*dotprod)

        dcc(iatm1,1)=dcc(iatm1,1)-rma*orij(j,1)*force
        dcc(iatm1,2)=dcc(iatm1,2)-rma*orij(j,2)*force
        dcc(iatm1,3)=dcc(iatm1,3)-rma*orij(j,3)*force
        dcc(iatm2,1)=dcc(iatm2,1)-rmb*orij(j,1)*force
        dcc(iatm2,2)=dcc(iatm2,2)-rmb*orij(j,2)*force
        dcc(iatm2,3)=dcc(iatm2,3)-rmb*orij(j,3)*force

       enddo

       do j=1,ntcons

        iatm1=conlist(j,1)
        iatm2=conlist(j,2)

        ccc(iatm1,1)=ccc(iatm1,1)+0.5*dcc(iatm1,1)
        ccc(iatm1,2)=ccc(iatm1,2)+0.5*dcc(iatm1,2)
        ccc(iatm1,3)=ccc(iatm1,3)+0.5*dcc(iatm1,3)
        ccc(iatm2,1)=ccc(iatm2,1)+0.5*dcc(iatm2,1)
        ccc(iatm2,2)=ccc(iatm2,2)+0.5*dcc(iatm2,2)
        ccc(iatm2,3)=ccc(iatm2,3)+0.5*dcc(iatm2,3)

       enddo

       mxdiff=mxdiff*0.5

       if(mxdiff.lt.tol)goto 100

      enddo

100   continue

      do i=1,natms

       iatm=atmtype(i)

       nvv(i,1)=(ccc(i,1)-ncc(i,1))*rtstep
       nvv(i,2)=(ccc(i,2)-ncc(i,2))*rtstep
       nvv(i,3)=(ccc(i,3)-ncc(i,3))*rtstep

       tmpvx=0.5*(nvv(i,1)+vvv(i,1))
       tmpvy=0.5*(nvv(i,2)+vvv(i,2))
       tmpvz=0.5*(nvv(i,3)+vvv(i,3))

       ekin=ekin+typmass(iatm)*(tmpvx**2+tmpvy**2+tmpvz**2)

       fff(i,1)=(nvv(i,1)-vvv(i,1))*typmass(iatm)*rtstep
       fff(i,2)=(nvv(i,2)-vvv(i,2))*typmass(iatm)*rtstep
       fff(i,3)=(nvv(i,3)-vvv(i,3))*typmass(iatm)*rtstep

       vvv(i,1)=nvv(i,1)
       vvv(i,2)=nvv(i,2)
       vvv(i,3)=nvv(i,3)

      enddo

      call tool_rebox(natms,mxatms,latt,rlatt,ccc) 

      ekin=0.5*ekin

      return

      END
