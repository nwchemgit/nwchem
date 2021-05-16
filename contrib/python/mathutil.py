from array import *
from math import sqrt

def zerovector(n):
    a = array('d',range(n))
    for i in range(n):
        a[i] = 0.0
    return a

def zeromatrix(n,m):
    a = list(range(n))
    for i in range(n):
        a[i] = zerovector(m)
    return a

def copyvector(x):
    a = array('d',range(len(x)))
    for i in range(len(x)):
        a[i] = x[i]
    return a

def copymatrix(x):
    n = len(x)
    m = len(x[0])
    a = zeromatrix(n,m)
    for i in range(n):
        for j in range(m):
            a[i][j] = x[i][j]
    return a

def transpose(x):
    n = len(x)
    m = len(x[0])
    a = zeromatrix(m,n)
    for i in range(n):
        for j in range(m):
            a[j][i] = x[i][j]
    return a

def dot(a,b):
    sum = 0.0
    for i in range(len(a)):
        sum = sum + a[i]*b[i]
    return sum

def mxm(a,b):
    n = len(a)
    k = len(a[0])
    kk= len(b)
    m = len(b[0])
    if kk != k:
        raise "matrices do not conform for multiplication"
    c = zeromatrix(n,m)
    for i in range(n):
        for j in range(m):
            sum = 0.0
            for l in range(k):
                sum = sum + a[i][l]*b[l][j]
            c[i][j] = sum
    return c

def mxv(a,b):
    n = len(a)
    k = len(a[0])
    kk = len(b)
    if k != kk:
        raise "matrix and vector do not conform for multiplication"
    c = zerovector(n)
    for i in range(n):
        c[i] = dot(a[i],b)
    return c

def printvector(a):
    n = len(a)
    for i in range(n):
        print ("%12.5e "%a[i]),
    print(" ")

def printmatrix(a):
    n = len(a)
    for i in range(n):
        printvector(a[i])

def numderiv(func,x,step,eps):
    '''
    Use central differences to compute the gradient and diagonal
    elements of the Hessian. 
    func(x) = function to be differentiated
    x[] = (array) point at which to differentiate
    step[] = (array) remembers finite difference step between
    .        successive calls.  Set to zero on first call 
    .        or set close to appropriate value 
    eps = expected precision in func
    
    Some care is taken to adjust the step so that the gradient and
    Hessian diagonal are estimated with about 4 digits of precision
    but some noise is unavaoidable due either to the noise in the
    function or cubic/higher terms in the Taylor expansion.
    '''
    
    n = len(x)
    g = zerovector(n)
    h = zerovector(n)
    f0 = func(x)
    for i in range(n):
        if step[i] == 0.0:
            step[i] = max(abs(x[i])*0.01,0.0001)
        xi = x[i]
        while 1:
            x[i] = xi + step[i]
            f1 = func(x)
            if abs(f1-f0) < (1e4*eps):
                #print ' Increasing step ',i,step[i],abs(f1-f0)
                step[i] = step[i]*2.0
            elif abs(f1-f0) > (1e5*eps):
                #print ' Decreasing step ',i,step[i],abs(f1-f0)
                step[i] = step[i]/3.0
            else:
                break
        x[i] = xi - step[i]
        fm1 = func(x)
        x[i] = xi
        g[i] = (f1 - fm1)/(2*step[i])
        h[i] = (f1 + fm1 - 2.0*f0)/(step[i]*step[i])

    return (f0,g,h)


def quadfit(alpha0, f0, alpha1, f1, alpha2, f2):
    '''
    Given 3 points compute the gradient and hessian at point 0
    using a quadratic fit.
    '''
    delta1 = alpha1 - alpha0
    delta2 = alpha2 - alpha0
    d1 = (f1 - f0)/delta1
    d2 = (f2 - f0)/delta2
    h0 = 2.0*(d1 - d2)/(delta1-delta2)
    g0 = d1 - 0.5*h0*delta1

    test1 = f0 + g0*delta1 + 0.5*h0*delta1*delta1
    test2 = f0 + g0*delta2 + 0.5*h0*delta2*delta2
    return (f0, g0, h0)

def takestep(x0, s, alpha):
    x = zerovector(len(x0))
    for j in range(len(x)):
        x[j] = x0[j] + s[j]*alpha
    return x

def quadratic_step(trust, g0, h0):
    if h0 > 0:
        delta2 = -g0/h0
        if abs(delta2) > trust:
            print ("                    Step restriction: %12.5f %12.5f " % (delta2, trust))
            delta2 = abs(trust*delta2)/delta2
    else:
        print ("                    Negative curvature ")
        delta2 = -abs(trust*g0)/g0
    return delta2


def linesearch(func, x0, s, lsgrad, eps):
    # Assume here that some halfway OK preconditioning
    # is being used so we expect a step around unity.
    # Nevertheless, must exercise some caution.

    # First step in small increments until we've either
    # bracketed the minimum or gone downhil with enough
    # energy difference to start fitting

    print (" Line search: step   alpha     grad     hess        value")
    print ("              ---- --------- -------- --------  ----------------")

    trust = 0.2
    alpha0 = 0.0
    f0 = func(x0)
    print ("                   %9.2e %8.1e          %16.8f" % (alpha0 , lsgrad , f0))

    if lsgrad < 0:
        alpha1 = alpha0 + trust
    else:
        alpha1 = alpha0 - trust
    f1 = func(takestep(x0,s,alpha1))
    print ("                   %9.2e                   %16.8f" % \
          (alpha1, f1))
    
    while f1 > f0:
        if trust < 0.00125:
            print (" system is too badly conditioned for initial step")
            return (alpha0,f0)          # Cannot seem to find my way
        trust = trust * 0.5
        if lsgrad < 0:
            alpha1 = alpha0 + trust
        else:
            alpha1 = alpha0 - trust
        f1 = func(takestep(x0,s,alpha1))
        print ("                   %9.2e                   %16.8f" % \
              (alpha1, f1))
        
    g0 = lsgrad
    h0 = (f1-f0-alpha1*g0)/alpha1**2
    if f1 < f0:
        g0 = g0 + h0*(alpha1 - alpha0)
        alpha0, alpha1, f0, f1 = alpha1, alpha0, f1, f0
    
    alpha2 = alpha0 + quadratic_step(trust,g0,h0)

    nbackstep =0

    for iter in range(1,10):
        f2 = func(takestep(x0,s,alpha2))
        #print ' alphas ', alpha0, alpha1, alpha2
        #print ' fs     ', f0, f1, f2
        
        if iter == 1:
            f2prev = f2
            
        print ("                   %9.2e                   %16.8f" % \
              (alpha2, f2))

        # Check for convergence or insufficient precision to proceed further
        if (abs(f0-f1)<(10*eps)) and (abs(f1-f2)<(10*eps)):
            print ("                  "),
            print (" Insufficient precision ... terminating LS")
            break
        
        if (f2-f2prev) > 0:
            # New point is higher than previous worst
            if  nbackstep < 3:
                nbackstep = nbackstep + 1
                print ("                  "),
                print (" Back stepping due to uphill step")
                trust = max(0.01,0.2*abs(alpha2 - alpha0))  # Reduce trust radius
                alpha2 = alpha0 + 0.2*(alpha2 - alpha0)
                continue
        elif (f2-f0) < 0:
            trust = min(4.0,trust*2.0)  # Seem to have narrowed the search

        nbackstep = 0

        f2prev = f2

        # Order in increasing energy
        if f1 < f0:
            alpha0, alpha1, f0, f1 = alpha1, alpha0, f1, f0
        if f2 < f0:
            alpha0, alpha2, f0, f2 = alpha2, alpha0, f2, f0
        if f2 < f1:
            alpha1, alpha2, f1, f2 = alpha2, alpha1, f2, f1

        (f0, g0, h0) = quadfit(alpha0, f0, alpha1, f1, alpha2, f2)
        print ("              %4i %9.2e %8.1e %8.1e %16.8f" % \
              (iter, alpha0, g0, h0, f0))

        if (h0>0.0) and (abs(g0) < 0.03*abs(lsgrad)):
            print ("                  "),
            print (" gradient reduced 30-fold ... terminating LS")
            break
        
        # Determine the next step
        delta = quadratic_step(trust,g0,h0)
        alpha2 = alpha0 + delta
        df = g0*delta + 0.5*h0*delta*delta
        if abs(df) < 10.0*eps:
            print ("                  "),
            print (" projected energy reduction < 10*eps ... terminating LS")
            break
        
        

    return (alpha0, f0)

def jacobi(ainput):
    '''
    Diagonalize a real symmetric matrix using the variable threshold
    cyclic Jacobi method.
    
    (v,e) = jacobi(a)
    
    Input: a[n][m] is a real symmetric matrix

    Returns: (v,e) where v is the list of eigenvectors and e is an
    array of the corresponding eigenvalues in ascending order.
    v[k] is a vector containing the kth eigenvector.  These satisfy

    A*Vt = Vt*e

    or

    V*A = e*V

    or
    
    sum(j)(a[i][j]v[k][j]) = e[k]*v[k][i]
    '''
    a = copymatrix(ainput)
    n = len(a)
    m = len(a[0])
    if n != m:
        raise 'Matrix must be square'
    for i in range(n):
        for j in range(m):
            if a[i][j] != a[j][i]:
                raise ' Matrix must be symmetric'

    tolmin = 1e-14
    tol = 1e-4

    v = zeromatrix(n,n)
    for i in range(n):
        v[i][i] = 1.0

    maxd = 0.0
    for i in range(n):
        maxd = max(abs(a[i][i]),maxd)

    for iter in range(50):
        nrot = 0
        for i in range(n):
            for j in range(i+1,n):
                aii = a[i][i]
                ajj = a[j][j]
                daij = abs(a[i][j])
                if daij > tol*maxd: # Screen small elements
                    nrot = nrot + 1
                    s = aii - ajj
                    ds = abs(s)
                    if daij > (tolmin*ds): # Check for sufficient precision
                        if (tol*daij) > ds:
                            c = s = 1/sqrt(2.)
                        else:
                            t = a[i][j]/s
                            u = 0.25/sqrt(0.25+t*t)
                            c = sqrt(0.5+u)
                            s = 2.*t*u/c
                            
                        for k in range(n):
                            u = a[i][k]
                            t = a[j][k]
                            a[i][k] = s*t + c*u
                            a[j][k] = c*t - s*u

                        for k in range(n):
                            u = a[k][i]
                            t = a[k][j]
                            a[k][i] = s*t + c*u
                            a[k][j]= c*t - s*u

                        for k in range(n):
                            u = v[i][k]
                            t = v[j][k]
                            v[i][k] = s*t + c*u
                            v[j][k] = c*t - s*u

                        a[j][i] = a[i][j] = 0.0
                        maxd = max(maxd,abs(a[i][i]),abs(a[j][j]))
        if nrot == 0 and tol <= tolmin:
            break
        tol = max(tolmin,tol*0.99e-2)

    if nrot != 0:
        raise "Jacobi iteration did not converge in 50 passes"

    # Sort eigenvectors and values into increasing order
    e = zerovector(n)
    for i in range(n):
        e[i] = a[i][i]
        for j in range(i):
            if e[j] > e[i]:
                (e[i],e[j]) = (e[j],e[i])
                (v[i],v[j]) = (v[j],v[i])

    return (v,e)

def hessian_update_bfgs(hp, dx, g, gp):
    '''
    Apply the BFGS update to the approximate Hessian h[][].

    hp[][] = Hessian matrix from previous iteration
    dx[]  = Step from previous iteration
    .       (dx[] = x[] - xp[] where xp[] is the previous point)
    g[]   = gradient at current point
    gp[]  = gradient at previous point

    Returns the updated hessian
    '''

    n = len(hp)
    hdx  = mxv(hp,dx)
    dg = zerovector(n)
    for i in range(n):
        dg[i] = g[i] - gp[i]
    
    dxhdx = dot(dx,hdx)
    dxdx  = dot(dx,dx)
    dxdg  = dot(dx,dg)
    dgdg  = dot(dg,dg)
    h = copymatrix(hp)

    if (dxdx > 0.0) and (dgdg > 0.0) and (abs(dxdg/sqrt(dxdx*dgdg)) > 1.e-4):
        for i in range(n):
            for j in range(n):
                h[i][j] = h[i][j] + dg[i]*dg[j]/dxdg - hdx[i]*hdx[j]/dxhdx
    else:
        print (' BFGS not updating dxdg (%e), dgdg (%e), dxhdx (%f), dxdx(%e)' % (dxdg, dgdg, dxhdx, dxdx))

    return h

def quasinr(func, guess, tol, eps, printvar=None):
    '''
    Unconstrained minimization of a function of n variables
    without analytic derivatives using quasi-Newtwon with BFGS update
    and numerical gradients.
    
    func(x) is a function that takes an array of n values and
    returns the function value
    
    guess[] is an array of n values for the initial guess
    
    tol is the convergence criterion for the maximum value
    of the gradient
    
    eps is the expected precision in the function value
    
    printvar(x) is an optional user function to print the values of
    parameters each macro iteration
    '''
    
    n = len(guess)
    x = copyvector(guess)
    s = zerovector(n)
    g = zerovector(n)
    gp = zerovector(n)
    step = zerovector(n)
    hessian = zeromatrix(n,n)

    alpha = 0.0

    for iter in range(50*n):
        (value,g,h) = numderiv(func, x, step, eps)
        gmax = max(map(abs,g))
        
        print (' ')
        print (' iter    gmax         value ')
        print (' ---- --------- ----------------')
        print ("%4i %9.2e %16.8f" % (iter,gmax,value))

        if (printvar):
            printvar(x)

        if gmax < tol:
            print (' Converged!')
            break

        if iter == 0:
            for i in range(n):
                hessian[i][i] = max(abs(h[i]),1e-4)
        else:
            hessian = hessian_update_bfgs(hessian, s, g, gp)

        (v,e) = jacobi(hessian)
        emax = max(map(abs,e))
        emin = emax*1e-4   # Control noise in small eigenvalues
        print ('\n Eigenvalues of the Hessian:')
        printvector(e)

        # Transform to spectral form, take step, transform back
        gs = mxv(v,g)
        for i in range(n):
            if e[i] < emin:
                print (' Mode %d: small/negative eigenvalue (%f).' % (i, e[i]))
                s[i] = -gs[i]/emin
            else:
                s[i] = -gs[i]/e[i]
        s = mxv(transpose(v),s)

        # Apply overall step restriction ... better LS obviates this
        scale = 1.0
        for i in range(n):
            trust = max(abs(x[i]),abs(x[i]/sqrt(max(1e-4,abs(hessian[i][i])))))
            if abs(s[i]) > trust:
                print (' restricting ', i, trust, abs(x[i]), \
                      abs(x[i]/sqrt(abs(hessian[i][i]))), s[i])
                scale = min(scale,trust/abs(s[i]))
        if scale != 1.0:
            for i in range(n):
                s[i] = s[i]*scale
        
        (alpha,value) = linesearch(func, x, s, dot(s,g), eps)

        if alpha == 0.0:
            print (' Insufficient precision to proceed further')
            break

        for i in range(n):
            s[i] = s[i]*alpha
            x[i] = x[i] + s[i]
            gp[i]= g[i]
    return (value,x)


def cgminold(func, dfunc, guess, tol):
    '''
    Simple conjugate gradient assuming analtyic derivatives.
    '''
    n = len(guess)
    x = copyvector(guess)
    s = zerovector(n)
    g = zerovector(n)
    gp= zerovector(n)
    value = func(x)

    for iter in range(10*n):
        g = dfunc(x)
        gmax = max(map(abs,g))
        
        print (' ')
        print (' iter    gmax         value ')
        print (' ---- --------- ----------------')
        print ("%4i %9.2e %16.8f" % (iter,gmax,value))
        
        if gmax < tol:
            print (' Converged!')
            break
        
        if (iter == 0) or ((iter%20) == 0):
            beta = 0.0
        else:
            beta = (dot(g,g) - dot(g,gp))/(dot(s,g)-dot(s,gp))

        for i in range(n):
            s[i] = -g[i] + beta*s[i]

        (alpha,value) = linesearch(func, x, s, dot(s,g), 1e-12)

        for i in range(n):
            s[i] = s[i]*alpha
            x[i] = x[i] + s[i]
            gp[i]= g[i]

    return (value,x)


def cgmin(func, dfunc, guess, tol, precond=None, reset=None):
    '''
    Conjugate gradient with optional preconditioning and
    use of analytic gradients.
    '''
    n = len(guess)
    x = copyvector(guess)
    s = zerovector(n)
    g = zerovector(n)
    gp= zerovector(n)
    value = func(x)

    if not reset:
        reset = n
    reset = min(reset,n)

    for iter in range(10*n):
        g = dfunc(x)
        gmax = max(map(abs,g))
        
        print (' ')
        print (' iter    gmax         value ')
        print (' ---- --------- ----------------')
        print ("%4i %9.2e %16.8f" % (iter,gmax,value))

        if gmax < tol:
            print (' Converged!')
            break
        
        if precond:
            precondg = precond(g)
        else:
            precondg = g
        
        if (iter % reset) == 0:
            beta = 0.0
        else:
            beta = (dot(precondg,g) - dot(precondg,gp))/(dot(s,g)-dot(s,gp))

        for i in range(n):
            s[i] = -precondg[i] + beta*s[i]

        (alpha,value) = linesearch(func, x, s, dot(s,g),
                                   max(1e-16,abs(value)*1e-12))

        for i in range(n):
            s[i] = s[i]*alpha
            x[i] = x[i] + s[i]
            gp[i]= g[i]
        
    return (value,x)

def cgmin2(func, guess, tol, eps, printvar=None,reset=None):
    '''
    Unconstrained minimization of a function of n variables
    without analytic derivatives using conjugate gradient with
    diagonal preconditioning.
    
    func(x) is a function that takes an array of n values and
    returns the function value
    
    guess[] is an array of n values for the initial guess
    
    tol is the convergence criterion for the maximum value
    of the gradient
    
    eps is the expected precision in the function value
    
    printvar(x) is an optional user function to print the values of
    parameters each iteration
    
    reset is the number of iterations between forced resets of the
    conjugacy.  In principle this could be n but noise in the
    numerical gradients makes a smaller number a better choice.
    '''
    
    n = len(guess)
    x = copyvector(guess)
    s = zerovector(n)
    g = zerovector(n)
    gp = zerovector(n)
    step = zerovector(n)
    precondg = zerovector(n)

    alpha = 0.0

    if not reset:
        reset = n
    reset = min(reset,n)

    for iter in range(50*n):
        (value,g,hh) = numderiv(func, x, step, eps)
        gmax = max(map(abs,g))
        
        print (' ')
        print (' iter    gmax         value ')
        print (' ---- --------- ----------------')
        print ("%4i %9.2e %16.8f" % (iter,gmax,value))

        if (printvar):
            printvar(x)

        if gmax < tol:
            print (' Converged!')
            break

        if (iter % reset) == 0:
            # On the first iteration or if not applying conjugacy
            # we can recompute the diagonal preconditioner
            h = copyvector(hh)
            for i in range(n):
                h[i] = max(abs(h[i]),1e-6)
        
        # Preconditioning with the diagonal of the Hessian
        for i in range(n):
            precondg[i] = g[i] / h[i]
            
        # Should be able to reset every n steps but noisy gradients
        # means that we don't have enough info.
        if (iter % reset) == 0:
            if iter != 0:
                print(" Resetting conjugacy")
            beta = 0.0
        else:
            beta = (dot(precondg,g) - dot(precondg,gp))/(dot(s,g)-dot(s,gp))

        for i in range(n):
            s[i] = -precondg[i] + beta*s[i]

        (alpha,value) = linesearch(func, x, s, dot(s,g), eps)

        if alpha == 0.0:
            # LS failed, probably due to lack of precision.
            if beta != 0.0:
                print ("LS failed - trying preconditioned steepest descent direction")
                for i in range(n):
                    s[i] = -g[i]
                (alpha,value) = linesearch(func, x, s, dot(s,g), eps)
            if alpha == 0.0:
                print (" Insufficient precision to proceed further")
                break

        for i in range(n):
            s[i] = s[i]*alpha
            x[i] = x[i] + s[i]
            gp[i]= g[i]
    return (value,x)
    

if __name__ == '__main__':

    def precond(g):
        # Used to test optional preconditioner for cgmin().
        precondg = copyvector(g)
        for i in range(len(g)):
            precondg[i] = precondg[i]/(i+2.0)
        return precondg

    def df(x):
        d = zerovector(len(x))
        for i in range(len(x)):
            d[i] = x[i]*(i+1)
            for j in range(len(x)):
                d[i] = d[i] + x[j]
        return d

    def f(x):
        sum = 0.0
        for i in range(len(x)):
            for j in range(len(x)):
                sum = sum + 0.5*x[i]*x[j]
        for i in range(len(x)):
            sum = sum + 0.5*x[i]*x[i]*(i+1)
        return sum


    print ('\n\n   TESTING QUASI-NR SOLVER \n\n')
    quasinr(f, [1.,0.5,0.3,-0.4], 1e-4, 1e-10)

    print ('\n\n   TESTING GC WITH NUM. GRAD. AND DIAG. PRECOND.\n\n')
    cgmin2(f, [1.,0.5,0.3,-0.4], 1e-4, 1e-10, reset=20)

    print ('\n\n   TESTING GC WITH ANAL. GRAD. AND WITHOUT OPTIONAL PRECOND.\n\n')
    cgmin(f, df, [1.,0.5,0.3,-0.4], 1e-4)

    print ('\n\n   TESTING GC WITH ANAL. GRAD. AND WITH OPTIONAL PRECOND.\n\n')
    cgmin(f, df, [1.,0.5,0.3,-0.4], 1e-4, precond=precond)

    print ('\n\n   TESTING GC WITH ANAL. GRAD. AND NO PRECOND.\n\n')
    cgminold(f, df, [1.,0.5,0.3,-0.4], 1e-4)

    print ('\n\n   TESTING JACOBI EIGENSOLVER\n\n')
    n = 50
    a = zeromatrix(n,n)
    for i in range(n):
        for j in range(i,n):
            a[j][i] = a[i][j] = (i*j+1.)/(i+j+1.)

    (v,e)= jacobi(a)

    print (' eigenvalues')
    printvector(e)
    #print ' v '
    #printmatrix(v)
    ev = mxm(v,a)
    for i in range(n):
        err = 0.0
        for j in range(n):
            err = max(err,abs(ev[i][j] - e[i]*v[i][j]))
        err = err/(n*max(1.0,abs(e[i])))
        if err > 1e-12:
            print (' Error in eigenvector ', i, err)
