__author__ = 'pja'

from numpy import *
from numpy.linalg import *
from numpy.polynomial import *
from math import isinf
from operator import mod


def characteristic( A , v ):
  return det(v*eye(A.shape[0]) - A)

def companion( coefs ):
  N = len(coefs) - 1
  A = matrix(zeros((N,N)))
  for i in range(0,N):
      A[i,N-1] = -coefs[i]/coefs[N]
  for i in range(0,N-1):
      A[i+1,i] = 1
  return A

def polyDerivative( poly ):
  N = len(poly)-1
  deriv = [0]*N
  for power,value in enumerate(poly[1:]):
      deriv[power] = value*(power+1)
  return deriv


def polyDegree( poly ):
    degree = len(poly)-1
    while degree >= 0:
        if poly[degree]: break
        degree -= 1
    return degree

def polyDivide( numerator , denominator ):
    nn = len(numerator)-1; nd = len(denominator)-1
    while nd >= 0 and denominator[nd] == 0: nd -= 1
    if nd < 0: raise ValueError('Divided by zero')

    r = numerator[:]
    q = [0]*len(numerator)

    for k in range(nn-nd,-1,-1):
        q[k]=r[nd+k]/float(denominator[nd])
        for j in range(k+nd,k-1,-1): r[j] -= q[k]*denominator[j-k]
    # Set it exactly to zero so it knows to skip over these later on
    for j in range(nd,nn+1): r[j] = 0

    return q,r

def polyEvaluate( poly , x ):
    # for infinity only the largest power needs to be evaluated
    # also the normal code below won't work since 0*inf = NaN
    if isinf(x):
        degree = polyDegree(poly)
        if degree == -1: return 0
        if degree == 0: return poly[0]
        if mod(degree,2) == 0: x = float('inf')
        return sign(poly[degree])*x

    # Evaluate using Horner's Method
    ret = 0

    for c in reversed(poly):
        ret = ret*x + c

    return ret


def countChange( sequence ):
    # find first non-zero entry
    isPlus = None
    for c in sequence:
        if c > 0:
            isPlus = True;break
        elif c < 0:
            isPlus = False;break
    if isPlus is None: return 0

    signChanges = 0
    for c in sequence:
        if isPlus:
            if c < 0: isPlus = False; signChanges += 1
        else:
            if c > 0: isPlus = True; signChanges += 1
    return signChanges

def numRealRoots( poly , minVal , maxVal ):
    # construct the sturm sequence recursively
    def s( polyA , polyB  , x ):
        # divide the polynomial and take the remainder
        q,r = polyDivide(polyA,polyB)

        # see if only a constant term remains
        if polyDegree(r) <= 0: return [-r[0]]

        # calculate the negative of the remainder
        polyC = [ -c for c in r ]

        return [polyEvaluate(polyC,x)] + s( polyB, polyC , x )

    polyD = polyDerivative( poly )

    sequenceMin = [polyEvaluate(poly,minVal) , polyEvaluate(polyD,minVal)] + s(poly,polyD,minVal)
    sequenceMax = [polyEvaluate(poly,maxVal) , polyEvaluate(polyD,maxVal)] + s(poly,polyD,maxVal)

    return countChange(sequenceMin) - countChange(sequenceMax)


def findRealRoots( poly , searchW = 20 , tol=1e-10 , maxIterations = 2000 ):

    numRoots = numRealRoots( poly , float('-inf'),float('inf'))
    if not numRoots: return []

    # Find a bound which contains all the real roots
    iter = 0
    while True and iter < maxIterations:
        if numRoots == numRealRoots( poly , -searchW , searchW ):
            break
        searchW *= 2
        iter += 1

    if iter >= maxIterations: raise RuntimeError('Too many iterations bounding all roots')

    # SIMPLIFY CODE LIKE THE JAVA VERSION

    # tighten the bounds, required for the next step to work
    l = -searchW; u = searchW
    while True:
        m = (l+u)/2.0
        if numRealRoots(poly,l,m) == 0:
            l = m
        else: break;
    while True:
        m = (l+u)/2.0
        if numRealRoots(poly,m,u) == 0:
            u = m
        else: break;

    # find upper and lower bounds for individual roots
    bounds = [l,u]
    iter = 0
    while len(bounds) < numRoots+1 and iter < maxIterations:
        m = (l+u)/2.0
        if numRealRoots(poly,l,m) == 1:
            bounds.insert(-1,m)
            l = m; u = bounds[-1]
        else:
            u = m
        iter += 1

    if iter >= maxIterations: raise RuntimeError('Too many iterations finding upper and lower bounds')

     # use bisection to refine the estimate for each root
    roots = [0]*numRoots
    for i in xrange(0,numRoots):
        l = bounds[i]; u = bounds[i+1]
        iter = 0

        # Use a relative threshold to determine convergence
        while (u-l)/abs(float(u)) > tol and iter < maxIterations:
            m = (l+u)/2.0
            if numRealRoots(poly,l,m) == 1:
                u = m
            else:
                l = m
            iter += 1
        if iter >= maxIterations: raise RuntimeError('Too many iterations for a single root')
        roots[i] = (u+l)/2.0


    return roots


def polishRoot( poly , root , tol , maxIterations=500 ):
    polyD = polyDerivative(poly)

    n = abs(root)
    if n == 0:  n = 1

    for i in xrange(0,maxIterations):
        p = polyEvaluate(poly,root)
        d = polyEvaluate(polyD,root)
        if d == 0:
            print 'Derivative vanished'
            return
        delta = p/d
        root -= delta
        if abs(delta)/n < tol:
            break
    return root


coefs = [-1.322309e+02 , 3.713984e+02 , -5.007874e+02 , 3.744386e+02 ,-1.714667e+02  , 4.865014e+01 ,-1.059870e+01  ,  1.642273e+00 ,-2.304341e-01,2.112391e-03,-2.273737e-13 ]

# TODO Is there a way to avoid overflow by evaluating sign changes using horn's method?

#print polyDiv(a,polyDerivative(a))

#print numRealRoots(a,float('-inf'),float('inf'))
#print numRealRoots(a,-3,-2)
#print numRealRoots(a,2,3)
#print numRealRoots(a,-3,3)


print numRealRoots([2,3,4],float('-inf'),float('inf'))
#
#print findRealRoots(coefs,50,1e-10)
#print eigvals(companion(coefs))

print 'Companion'
for n in eigvals(companion(coefs)):
    if n.imag == 0:
        n = n.real
        print 'Root = ',n,'  zero = ',polyEvaluate(coefs,n)

print 'Sturm'
for n in findRealRoots(coefs,50,1e-10):
    n = polishRoot(coefs,n,1e-15)
    print 'Root = ',n,'  zero = ',polyEvaluate(coefs,n)

print 'NumPy' # Used companion matrix technique
for n in polyroots(coefs):
    if n.imag == 0:
        n = n.real
        print 'Root = ',n,'  zero = ',polyEvaluate(coefs,n)

#r = isolateRoot(coefs,100,1e-4)
#print 'Estimated root r = ',r,' value = ',polyEvaluate(coefs,r)
#r = polishRoot(coefs,r,1e-20,100)
#print 'Polished root r = ',r,' value = ',polyEvaluate(coefs,r)
#
#found = eigvals(companion(coefs))[8].real
#print 'Companion root r = ',found,' value = ',polyEvaluate(coefs,found)
#r = polishRoot(coefs,found,1e-20,100)
#print 'Polished root r = ',r,' value = ',polyEvaluate(coefs,r)
