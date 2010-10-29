c
c $Id$
c

      subroutine smd_leapf_shake(natms,
     >                           ntcons,
     >                           tstep,
     >                           ekin,
     >                           mass,
     >                           icon1,
     >                           icon2,
     >                           consdist,
     >                           ncc,
     >                           nvv,
     >                           dcc,
     >                           nrij,
     >                           orij,
     >                           fff,
     >                           vvv,
     >                           ccc)

      implicit none
c
      integer natms
      integer ntcons
      double precision tstep
      double precision ekin
      double precision mass(natms)
      integer          icon1(ntcons)
      integer          icon2(ntcons)
      double precision consdist(ntcons)
      double precision ncc(natms,3),nvv(natms,3),dcc(natms,3)
      double precision nrij(ntcons,3),orij(ntcons,3)
      double precision fff(natms,3)
      double precision vvv(natms,3)
      double precision ccc(natms,3)

      integer i,j,it,maxit,iatm1,iatm2

      double precision   force,rma,rmb
      double precision   tmpvx,tmpvy,tmpvz
      double precision   tstepsq,rtstep,tol,mxdiff
      double precision   nrijsq,dijsq,diffsq,dotprod

      ekin=0.0
      tstepsq=tstep**2
      rtstep=1.0/tstep
      tol=1.0d-8
      maxit=100
      mxdiff = 0.0d0

      do i=1,ntcons

       iatm1=icon1(i)
       iatm2=icon2(i)

       orij(i,1)=ccc(iatm1,1)-ccc(iatm2,1)
       orij(i,2)=ccc(iatm1,2)-ccc(iatm2,2)
       orij(i,3)=ccc(iatm1,3)-ccc(iatm2,3)

      enddo

      call smd_lat_rebox(ntcons,orij)

      do i=1,natms

       ncc(i,1)=ccc(i,1)
       ncc(i,2)=ccc(i,2)
       ncc(i,3)=ccc(i,3)

       nvv(i,1)=vvv(i,1)+fff(i,1)*tstep/mass(i)
       nvv(i,2)=vvv(i,2)+fff(i,2)*tstep/mass(i)
       nvv(i,3)=vvv(i,3)+fff(i,3)*tstep/mass(i)

       ccc(i,1)=ncc(i,1)+tstep*nvv(i,1)
       ccc(i,2)=ncc(i,2)+tstep*nvv(i,2)
       ccc(i,3)=ncc(i,3)+tstep*nvv(i,3)

      enddo 


      do i=1,maxit

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrij(j,1)=ccc(iatm1,1)-ccc(iatm2,1)
        nrij(j,2)=ccc(iatm1,2)-ccc(iatm2,2)
        nrij(j,3)=ccc(iatm1,3)-ccc(iatm2,3)

       enddo

       call smd_lat_rebox(ntcons,nrij) 

       do j=1,natms
        dcc(j,1)=0.0
        dcc(j,2)=0.0
        dcc(j,3)=0.0
       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

        nrijsq=nrij(j,1)**2+nrij(j,2)**2+nrij(j,3)**2
        dijsq=consdist(j)**2
        diffsq=dijsq-nrijsq
        mxdiff=max(mxdiff,abs(diffsq)/consdist(j))

        dotprod=orij(j,1)*nrij(j,1)
     $         +orij(j,2)*nrij(j,2)
     $         +orij(j,3)*nrij(j,3)

        rma= tstepsq/mass(iatm1)
        rmb=-tstepsq/mass(iatm2)
        force=diffsq/(-2.0*(rma-rmb)*dotprod)

        dcc(iatm1,1)=dcc(iatm1,1)-rma*orij(j,1)*force
        dcc(iatm1,2)=dcc(iatm1,2)-rma*orij(j,2)*force
        dcc(iatm1,3)=dcc(iatm1,3)-rma*orij(j,3)*force
        dcc(iatm2,1)=dcc(iatm2,1)-rmb*orij(j,1)*force
        dcc(iatm2,2)=dcc(iatm2,2)-rmb*orij(j,2)*force
        dcc(iatm2,3)=dcc(iatm2,3)-rmb*orij(j,3)*force


       enddo

       do j=1,ntcons

        iatm1=icon1(j)
        iatm2=icon2(j)

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

       nvv(i,1)=(ccc(i,1)-ncc(i,1))*rtstep
       nvv(i,2)=(ccc(i,2)-ncc(i,2))*rtstep
       nvv(i,3)=(ccc(i,3)-ncc(i,3))*rtstep

       tmpvx=0.5*(nvv(i,1)+vvv(i,1))
       tmpvy=0.5*(nvv(i,2)+vvv(i,2))
       tmpvz=0.5*(nvv(i,3)+vvv(i,3))

       ekin=ekin+mass(i)*(tmpvx**2+tmpvy**2+tmpvz**2)

       fff(i,1)=(nvv(i,1)-vvv(i,1))*mass(i)*rtstep
       fff(i,2)=(nvv(i,2)-vvv(i,2))*mass(i)*rtstep
       fff(i,3)=(nvv(i,3)-vvv(i,3))*mass(i)*rtstep

       vvv(i,1)=nvv(i,1)
       vvv(i,2)=nvv(i,2)
       vvv(i,3)=nvv(i,3)

      enddo


      ekin=0.5*ekin
    
      return

      END
