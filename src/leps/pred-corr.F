c    *******************************************************************
c    ** gear 5-value predictor-corrector algorithm for translation.   **
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** gear, numerical initial value problems in ordinary            **
c    ** differential equations (prentice-hall,1971).                  **
c    **                                                               **
c    ** supplied routines:                                            **
c    **                                                               **
c    ** subroutine predic ( dt )                                      **
c    **    predicts the new positions, velocities, etc.               **
c    ** subroutine correc ( dt, m, k )                                **
c    **    corrects the positions, velocities etc. using gear method  **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** integer n                  number of molecules                **
c    ** real    dt                 timestep                           **
c    ** real    rx(n),ry(n),rz(n)  positions                          **
c    ** real    vx(n),vy(n),vz(n)  velocities                         **
c    ** real    ax(n),ay(n),az(n)  accelerations                      **
c    ** real    bx(n),by(n),bz(n)  third derivatives                  **
c    ** real    cx(n),cy(n),cz(n)  fourth derivatives                 **
c    ** real    fx(n),fy(n),fz(n)  forces                             **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** at each timestep, call predic, force, correc in order         **
c    ** followed by accumulation of thermodynamic quantities.         **
c    ** the force routine (not supplied here: see f.17) calculates    **
c    ** potential energy and forces on all atoms.                     **
c    *******************************************************************



        subroutine predic (n, dt )

c    *******************************************************************
c    ** predictor routine                                             **
c    **                                                               **
c    ** in timestep-scaled variables the predictor is the pascal      **
c    ** triangle matrix.  in unscaled variables it is a taylor series **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** predic is called to advance the coordinates, velocities etc.  **
c    ** by one timestep dt, prior to force evaluation.                **
c    *******************************************************************

        integer     n 
        double precision r(3,n), v(3,n), a(3,n), b(3,n), c(3,n)
        double precision dt
        integer     i,j
C ---- dt# are the Taylor expansion series terms dt1=dt, dt2=dt^2/2, dt3=dt^3/6 etc. 
        real        dt1, dt2, dt3, dt4

c    *******************************************************************

       dt1 = dt
       dt2 = dt1 * dt / 2.0
       dt3 = dt2 * dt / 3.0
       dt4 = dt3 * dt / 4.0

       do i = 1,n 
        do j=1,3
         r(j,i) = r(j,i) + dt1*v(j,i) + dt2*a(j,i)  
     1          +  dt3*b(j,i) + dt4*c(j,i)
         v(j,i) = v(j,i) + dt1*a(j,i) + dt2*b(j,i) + dt3*c(j,i)
         a(j,i) = a(j,i) + dt1*b(j,i) + dt2*c(j,i)
         b(j,i) = b(j,i) + dt1*c(j,i)
        enddo
       enddo 

        return
        end



        subroutine correc (n, dt, mass, ek )

c    *******************************************************************
c    ** corrector routine                                             **
c    **                                                               **
c    ** corrects positions, velocities etc. after force evaluation.   **
c    ** in timestep-scaled variables the numerical coefficients are   **
c    ** given by gear (ref above): 19/120, 3/4, 1, 1/2, 1/12.         **
c    ** in unscaled form these must be multiplied by factors          **
c    ** involving the timestep as shown here.                         **
c    **                                                               **
c    ** principal variables:                                          **
c    **                                                               **
c    ** real    m                        atomic mass                  **
c    ** real    ek                        kinetic energy per atom      **
c    ** real    c0,c1,c2,c3              gear coefficients            **
c    **                                                               **
c    ** usage:                                                        **
c    **                                                               **
c    ** it is assumed that intermolecular forces have been calculated **
c    ** and stored in force (gradient of the potential). correc simply
c       applies the corrector   **
c    ** equations based on the differences between predicted and      **
c    ** evaluated accelerations.  it also calculates kinetic energy.  **
c    *******************************************************************

        integer     n
        double precision r(3,n), v(3,n), a(3,n), b(3,n), c(3,n)
        double precision dt, force(3,n)
        integer     i,j
C ---- dt# are the Taylor expansion series terms dt1=dt, dt2=dt^2/2, dt3=dt^3/6 etc.
        double precision        dt1, dt2, dt3, dt4
        double precision       mass(n), ek
        double precision ai(3),corr(3)
        double precision  cr, cv, cb, cc
        double precision  c0, c1, c3, c4
        parameter ( c0 = 19.0 / 120.0, c1 = 3.0 / 4.0,
     1              c3 = 1.0 / 2.0,    c4 = 1.0 / 12.0 )

c    *******************************************************************

        dt1 = dt
        dt2 = dt1 * dt / 2.0
        dt3 = dt2 * dt / 3.0
        dt4 = dt3 * dt / 4.0

        cr = c0 * dt2
        cv = c1 * dt2 / dt1
        cb = c3 * dt2 / dt3
        cc = c4 * dt2 / dt4

        ek = 0.0

        do i = 1, n
         do j=1,3  
           
           ai(j) = force(j,i) / mass(i)
           corr(j) = ai(j) - a(j,i)

           r(j,i) = r(j,i) + cr * corr(j)
           v(j,i) = v(j,i) + cv * corr(j)
           a(j,i) = ai(j)
           b(j,i) = b(j,i) + cb * corr(j)
           c(j,i) = c(j,i) + cc * corr(j)
         
        
           ek = ek + v(j,i) ** 2 
         enddo 
       enddo 

       ek = 0.5 * mass(i) * ek

       return
       end




c $Id$
