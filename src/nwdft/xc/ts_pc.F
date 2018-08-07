c Copyright 2018 (C) Orbital-free DFT group at University of Florida
c Licensed under the Educational Community License, Version 2.0 
c (the "License"); you may not use this file except in compliance with 
c the License. You may obtain a copy of the License at
c
c    http://www.osedu.org/licenses/ECL-2.0
c
c Unless required by applicable law or agreed to in writing,
c software distributed under the License is distributed on an "AS IS"
c BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express
c or implied. See the License for the specific language governing
c permissions and limitations under the License.
c
c ---------------------------------------------------------------------
c
c Perdew & Constantin kinetic energy functional enhancement factor
c and derivatives
c
c           META GGA
C         utilizes ingredients:
c                              rho   -  density
c                              delrho - gradient of density
c                              laprho - laplacian of density
c
c     Written by:
c     Daniel Mejia-Rodriguez
c     QTP, Department of Physics, University of Florida

      Subroutine ts_pc(tol_rho, rho, delrho, laprho,
     &                 dfsdp, dfsdq, fs, facttwo)
      implicit none
c
      double precision fs !< [Output]
      double precision dfsdp !< [Output]
      double precision dfsdq !< [Output]
c
c     Charge Density
c
      double precision rho !< [Input]
c
c     Charge Density Gradient
c
      double precision delrho(3) !< [Input]
c
c     Charge Density Laplacian
c
      double precision laprho !< [Input]
c
c
      double precision facttwo ! 2 for o.s. 1 for c.s.
c
      double precision tol_rho, pi
      double precision rho83, rho53
      double precision p, g2
      double precision F83, F53, F23
      double precision afact2, rhoval, Pconst
c
      double precision amgga, bmgga
      double precision ge2,ge4,mge4,zpc,fpc,bfpc,q,lapval
      double precision exp8,exp9,factor,ftw,mge4den,opftw

c     functional derivatives below FFFFFFFFFFFF

      double precision dge2dp,dge2dq,dge4dp,dge4dq,dmge4dp,dmge4dq
      double precision dfactordp,dfactordq
      double precision dfpc,dfsdz
      
c     functional derivatives above FFFFFFFFFFFF
      
      parameter (F83=8.d0/3.d0, F53=5.d0/3.d0, F23=2d0/3d0)
      parameter (amgga=1.784720d0, bmgga=0.258304d0)
c
      pi=acos(-1d0)
      Pconst = (3.d0*pi**2)**F23
      afact2=1d0/facttwo
c
      if (rho.ge.tol_rho) then

            rhoval=rho*facttwo
            rho53 = rhoval**F53
            rho83 = rho53*rhoval
            g2 = delrho(1)*delrho(1) +
     &           delrho(2)*delrho(2) +
     &           delrho(3)*delrho(3)
            g2 = g2*facttwo*facttwo
            lapval = laprho*facttwo
c
            p = g2/(4d0*Pconst*rho83)
            q = lapval/(4d0*Pconst*rho53)

            ge2 = 5.d0/27.d0*p + 20.d0/9.d0*q
            dge2dp = 5.d0/27.d0
            dge2dq = 20.d0/9.d0

            ge4 = q*(8.d0/81.d0*q - 1.d0/9.d0*p) + 8.d0/243.d0*p*p
            dge4dp = 16.d0/243.d0*p - q/9.d0
            dge4dq = 16.d0/81.d0*q - p/9.d0

            ftw = F53*p
            opftw = 1d0 + ftw
            factor = ge4/opftw
            dfactordp = 1d0/opftw * (dge4dp - F53*factor)
            dfactordq = dge4dq/opftw

            if (factor.gt.1.d5) then
              mge4 = (1d0/ge4 + ge2/ge4 + 1d0)*opftw
              dmge4dp = (dge2dp + dge4dp - mge4*dfactordp)/factor
              dmge4dq = (dge2dq + dge4dq - mge4*dfactordq)/factor
            else
              mge4den = dsqrt(1.d0 + factor**2)
              mge4 = (1d0 + ge2 + ge4)/mge4den
              dmge4dp = (dge2dp + dge4dp - 
     &                   mge4*factor*dfactordp/mge4den)/mge4den
              dmge4dq = (dge2dq + dge4dq - 
     &                   mge4*factor*dfactordq/mge4den)/mge4den
            endif
c
            zpc = mge4 - ftw
            fpc = 0d0
            dfpc = 0d0
c
            if (zpc.le.0.001d0) then
              fpc = 0d0
              dfpc = 0d0
            elseif (zpc.lt.amgga) then
              exp8 = dexp(-amgga/zpc)
              exp9 = dexp(-amgga/(amgga-zpc))
              bfpc = (exp8*exp9 + exp8)/(exp8 + exp9)
              fpc = bfpc**bmgga
              dfpc = amgga*bmgga*fpc/(exp8 + exp9)
              if (amgga-zpc.gt.0.0013d0) then
                dfpc = dfpc*exp9*(1d0/zpc**2 +
     &                 (1d0-exp8)/(1d0+exp9*(amgga-zpc)**2))
              else
                dfpc = dfpc*exp9/zpc**2
              end if
            else
              fpc=1d0
              dfpc=0d0
            endif

            fs=ftw+zpc*fpc

            dfsdz = fpc + zpc*dfpc
            dfsdp = F53 + dfsdz*(dmge4dp - F53)
            dfsdq = dfsdz*dmge4dq

      else
c
        fs = 0d0
        dfsdp = 0d0
        dfsdq = 0d0
c
      end if
c
      return
      end
