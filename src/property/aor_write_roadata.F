      subroutine aor_write_roadata
     &   (luout, alfa, beta, gp, dipquad, dqpol)
      
      integer luout
      double precision alfa(3,3), beta(3,3), gp(3,3) 
      double precision dipquad(3,6), dqpol(3,3,3)

      write (luout,'(/1x,a)') 'roa begin'
      write (luout,'(/1x,a)') 'alpha:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', alfa(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', alfa(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', alfa(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', alfa(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', alfa(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', alfa(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', alfa(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', alfa(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', alfa(3,3)
      
      write (luout,'(/1x,a)') 'beta:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', beta(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', beta(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', beta(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', beta(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', beta(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', beta(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', beta(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', beta(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', beta(3,3)

      write (luout,'(/1x,a)') 'G-prime = -omega beta:'
      write (luout,'(1x,a,t10,e20.13)') 'XX', gp(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XY', gp(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XZ', gp(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YX', gp(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YY', gp(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YZ', gp(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZX', gp(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZY', gp(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZZ', gp(3,3)
      
c     second index of dipquadre:
c     XX=1, XY=YX=2, XZ=ZX=3, YY=4, YZ=ZY=5, ZZ=6
      write (luout,'(/1x,a)') 'dipole-quadrupole (Cartesian):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dipquad(1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dipquad(1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dipquad(1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dipquad(1,4)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dipquad(1,5)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dipquad(1,6)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dipquad(2,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dipquad(2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dipquad(2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dipquad(2,4)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dipquad(2,5)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dipquad(2,6)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dipquad(3,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dipquad(3,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dipquad(3,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dipquad(3,4)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dipquad(3,5)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dipquad(3,6)

      write (luout,'(/1x,a)') 'dipole-quadrupole (traceless):'
      write (luout,'(1x,a,t10,e20.13)') 'XXX', dqpol(1,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'XXY', dqpol(1,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'XXZ', dqpol(1,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'XYY', dqpol(1,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'XYZ', dqpol(1,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'XZZ', dqpol(1,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'YXX', dqpol(2,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'YXY', dqpol(2,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'YXZ', dqpol(2,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'YYY', dqpol(2,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'YYZ', dqpol(2,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'YZZ', dqpol(2,3,3)

      write (luout,'(1x,a,t10,e20.13)') 'ZXX', dqpol(3,1,1)
      write (luout,'(1x,a,t10,e20.13)') 'ZXY', dqpol(3,1,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZXZ', dqpol(3,1,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZYY', dqpol(3,2,2)
      write (luout,'(1x,a,t10,e20.13)') 'ZYZ', dqpol(3,2,3)
      write (luout,'(1x,a,t10,e20.13)') 'ZZZ', dqpol(3,3,3)

      write (luout,'(/1x,a/)') 'roa end'

      return
      end
c $Id$
