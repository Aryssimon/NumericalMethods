import numpy as np
import math
import matplotlib.pyplot as plt


def bracket(f, x1, h, c):
    x2 = x1 + h
    fx2 = f(x2)
    if fx2 > f(x1):
        h = -h
        x2 = x1 + h
        fx2 = f(x2)
    for i in range(1000):
        h *= c
        x3 = x2 + h
        fx3 = f(x3)
        if fx3 > fx2:
            return x1, x3
        x1 = x2
        x2 = x3
        fx2 = fx3


"""
def func(x):
    return np.sin(3 * x) + np.cos(x)
print(bracket(func, 1.5, 0.001, 1.1))
x = np.arange(0, 5, 0.01)
y = func(x)
plt.plot(x, y, 'r-')
plt.show()
"""


def golden(f, a, b, tol=1.0e-9):
    nIter = int(np.ceil(-2.078087 * np.log(tol / abs(b - a))))
    R = (-1 + np.sqrt(5)) / 2
    C = 1.0 - R
    x1 = R * a + C * b
    x2 = C * a + R * b
    fx1 = f(x1)
    fx2 = f(x2)
    for i in range(nIter):
        if fx1 > fx2:
            a = x1
            x1 = x2
            fx1 = fx2
            x2 = C * a + R * b
            fx2 = f(x2)
            if abs(x1 - x2) < tol:
                return x2, fx2
        else:
            b = x2
            x2 = x1
            fx2 = fx1
            x1 = R * a + C * b
            fx1 = f(x1)
            if abs(x2 - x1) < tol:
                return x1, fx1
    if fx1 < fx2:
        return x1, fx1
    return x2, fx2







def powell(F,x,h=0.1,tol=1.0e-6, sequences=None):

    def f(s): return F(x + s*v)    # F in direction of v

    if(sequences is not None):sequences+=[x.copy()]
    n = len(x)                     # Number of design variables
    df = np.zeros(n)               # Decreases of F stored here
    u = np.identity(n)             # Vectors v stored here by rows
    for j in range(30):            # Allow for 30 cycles:

        xOld = x.copy()            # Save starting point
        fOld = F(xOld)
      # First n line searches record decreases of F
        for i in range(n):
            v = u[i]
            a,b = bracket(f,0.0,h)
            s,fMin = search(f,a,b)
            df[i] = fOld - fMin
            fOld = fMin
            x = x + s*v
      # Last line search in the cycle
        v = x - xOld
        a,b = bracket(f,0.0,h)
        s,fLast = search(f,a,b)
        x = x + s*v
        if(sequences is not None):sequences+=[x.copy()]
      # Check for convergence
        if math.sqrt(np.dot(x-xOld,x-xOld)/n) < tol: return x,j+1
      # Identify biggest decrease & update search directions
        iMax = np.argmax(df)
        for i in range(iMax,n-1):
            u[i] = u[i+1]
        u[n-1] = v
    print("Powell did not converge")



## module downhill
''' x = downhill(F,xStart,side=0.1,tol=1.0e-6)
    Downhill simplex method for minimizing the user-supplied
    scalar function F(x) with respect to the vector x.
    xStart = starting vector x.
    side   = side length of the starting simplex (default is 0.1)
'''

def downhill(F,xStart,side=0.1,tol=1.0e-6, sequences=None):
    n = len(xStart)                 # Number of variables
    x = np.zeros((n+1,n))
    f = np.zeros(n+1)

  # Generate starting simplex
    x[0] = xStart
    for i in range(1,n+1):
        x[i] = xStart
        x[i,i-1] = xStart[i-1] + side
  # Compute values of F at the vertices of the simplex
    for i in range(n+1): f[i] = F(x[i])


  # Main loop
    for k in range(500):
        if(sequences is not None):
            sequences+=[x.copy()]

      # Find highest and lowest vertices
        iLo = np.argmin(f)
        iHi = np.argmax(f)
      # Compute the move vector d
        d = (-(n+1)*x[iHi] + np.sum(x,axis=0))/n
      # Check for convergence
        if math.sqrt(np.dot(d,d)/n) < tol: return x[iLo]

      # Try reflection
        xNew = x[iHi] + 2.0*d
        fNew = F(xNew)
        if fNew <= f[iLo]:        # Accept reflection
            x[iHi] = xNew
            f[iHi] = fNew
          # Try expanding the reflection
            xNew = x[iHi] + d
            fNew = F(xNew)
            if fNew <= f[iLo]:    # Accept expansion
                x[iHi] = xNew
                f[iHi] = fNew
        else:
          # Try reflection again
            if fNew <= f[iHi]:    # Accept reflection
                x[iHi] = xNew
                f[iHi] = fNew
            else:
              # Try contraction
                xNew = x[iHi] + 0.5*d
                fNew = F(xNew)
                if fNew <= f[iHi]: # Accept contraction
                    x[iHi] = xNew
                    f[iHi] = fNew
                else:
                  # Use shrinkage
                    for i in range(len(x)):
                        if i != iLo:
                            x[i] = (x[i] - x[iLo])*0.5
                            f[i] = F(x[i])
    print("Too many iterations in downhill")
    return x[iLo]
