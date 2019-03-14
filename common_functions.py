#!/usr/bin/env python

from math import gamma,e
import sys
#import numpy as np

def eprint(*args,**kwargs):
    print(*args,**kwargs, file=sys.stderr)
        
def vprint(verbose,*args,**kwargs):
    if verbose:
        print(*args,**kwargs, file=sys.stderr)
       
# Binary search 
def pfunc(plocal,r,N):
    q = 1.0-plocal
    
    # Find x10
    x = 0.0
    for j in range(1,11):
        x = 1.0 + (q*(plocal**r)*(x**(r+1.0)))

    # do the equation
    result = (1.0 - plocal*x)
    result = result/((r+1.0 - (r*x))*q)
    try:
        result = result/(x**(N+1))
    except OverflowError:
        # catch OverflowError resulting from large N.
        result = 0.0
    return result

def search_for_p(r,N,iterations=1000, min_plocal=0.0, max_plocal=1.0, tolerance=0.00000001,verbose=False):
    # Binary chop search for Plocal
    iteration = 0
    found = False
    
    #vprint(verbose,"SEARCH FOR P")
    #vprint(verbose,f'min {min_plocal}  max {max_plocal} verbose={verbose} r={r} N={N}')
    while (iteration < iterations):
        candidate = (min_plocal + max_plocal)/2.0 # start in the middle of the range
        result = pfunc(candidate,r,N)
        #print ("iteration =",iteration)
        #if verbose:
        #    vprint(verbose,f'candidate {candidate}  min {min_plocal}  max {max_plocal}')
        iteration += 1
        if iteration > iterations:
            found = False
            break
        elif (result > (0.99-tolerance)) and (result < (0.99+tolerance)):
            found = True
            P_local = candidate
            break
        elif result > 0.99:
            min_plocal = candidate
        else:
            max_plocal = candidate

    if (found == False):
        print ("Warning: P_local not found")

    return P_local

# Binary search provided by Joshua Hill

## Binary search algorithm to find value of p s.t. the expected value
## of the Maurer Universal Statistic equals mu_bar within tolerance
## presumes a decreasing function
#def solve_for_p(mu_bar, n, v, tolerance=1e-09):
#    assert n > 0
#
#    #This is a hackish way of checking to see if the difference is within approximately 4 ULPs
#    absEpsilon = 4.0 * max((np.nextafter(mu_bar, mu_bar+1.0) - mu_bar), (mu_bar - np.nextafter(mu_bar, mu_bar-1.0)))
#    #If we don't have numpy, then this will work for most of the ranges we're concerned with
#    #absEpsilon = 4.0*sys.float_info.epsilon
#
#    if mu_bar > EppM(1.0/float(n), n, v):
#        return False, 0.0
#
#    ldomain = 1.0/float(n)
#    hdomain = 1.0
#
#    lbound = ldomain
#    lvalue = float("inf")
#    hbound = hdomain
#    hvalue = float("-inf")
#
#    #Note that the bounds are in the interval [0, 1], so underflows
#    #are an issue, but overflows are not
#    center = (lbound + hbound) / 2
#    assert (center > ldomain) and (center < hdomain)
#
#    centerVal = EppM(center, n, v)
#
#    for rounds in range(1076):
#        if isclose(mu_bar, centerVal, tolerance, absEpsilon):
#            return True, center
#
#        if mu_bar < centerVal:
#            lbound = center
#            lvalue = centerVal
#        else:
#            hbound = center
#            hvalue = centerVal
##We now verify that ldomain <= lbound < center < hbound <= hdomain
##and that target in [ hvalue, lvalue ]
#        if lbound >= hbound:
#            print ("Bounds have converged after %d rounds and target was not found. Returning largest bound." % rounds)
#            return True, min(max(lbound, hbound), hdomain)
#
#        if (lbound < ldomain) or (lbound > hdomain) or (hbound < ldomain) or (hbound > hdomain):
#            print ("The current search interval is not a subset of the domain after %d rounds and target was not found." % rounds)
#            return False, 0.0
#
#        if (mu_bar > lvalue) or (mu_bar < hvalue):
#            print ("Target is not within the search interval after %d rounds" % rounds)
#            return False, 0.0
#
#        lastCenter = center
#        center = (lbound + hbound) / 2.0
#
#        if (center <= lbound) or (center >= hbound):
#            print ("The next center is outside of the search interval after %d rounds" % rounds)
#            return False, 0.0
#
#        if lastCenter == center:
#            print ("Detected cycle after %d rounds. Returning upper bound." % rounds)
#            return True, hbound
#
#        centerVal = EppM(center, n, v)
#
#        #invariant: if this isn't true, then this isn't loosely monotonic
#        if (centerVal < hvalue) or (centerVal > lvalue):
#            print ("CenterVal is not within the search value interval after %d rounds. Returning upper bound." % rounds)
#            return True, hbound
#
#    #We ran out of rounds for the binary search
#    if isclose(mu_bar, centerVal, tolerance, absEpsilon):
#        return True, p
#    else:
#        print ("Ran out of search rounds. Returning upper bound")
#        return True, min(hbound, hdomain)


# Continued Fraction Computation
# 6.5.31 Handbook of Mathematical Functions, page 263
#    Recursive implementation
def upper_incomplete_gamma(a,x,d=0,iterations=100):
    if d == iterations:
        if ((d % 2) == 1):
            return 1.0 # end iterations
        else:
            m = d/2
            return x + (m-a)
    if d == 0:
        result = ((x**a) * (e**(-x)))/upper_incomplete_gamma(a,x,d=d+1)
        return result
    elif ((d % 2) == 1):
        m = 1.0+((d-1.0)/2.0)
        return x+ ((m-a)/(upper_incomplete_gamma(a,x,d=d+1)))
    else:
        m = d/2
        return 1+(m/(upper_incomplete_gamma(a,x,d=d+1)))

# 6.5.31 Handbook of Mathematical Functions, page 263
#    Recursive implementation
def upper_incomplete_gamma2(a,x,d=0,iterations=100):
    if d == iterations:
        return 1.0 
    if d == 0:
        result = ((x**a) * (e**(-x)))/upper_incomplete_gamma2(a,x,d=d+1)
        return result
    else:
        m = (d*2)-1
        return (m-a)+x+ ((d*(a-d))/(upper_incomplete_gamma2(a,x,d=d+1)))

def lower_incomplete_gamma(a,x,d=0,iterations=100):
    if d == iterations:
        if ((d % 2) == 1):
            return 1.0 # end iterations
        else:
            m = d/2
            return x + (m-a)
    if d == 0:
        result = ((x**a) * (e**(-x)))/lower_incomplete_gamma(a,x,d=d+1)
        return result
    elif ((d % 2) == 1):
        m = d - 1
        n = (d-1.0)/2.0
        return a + m - (((a+n)*x)/lower_incomplete_gamma(a,x,d=d+1))
    else:
        m = d-1
        n = d/2.0
        return a+m+((n*x)/(lower_incomplete_gamma(a,x,d=d+1)))

def lower_incomplete_gamma2(a,x):
    return gamma(a)-upper_incomplete_gamma2(a,x)

def complimentary_incomplete_gamma(a,x):
    return 1.0-upper_incomplete_gamma(a,x)

# Scipy name mappings
def gammainc(a,x):
    return lower_incomplete_gamma(a,x)/gamma(a)

def gammaincc(a,x):
    return upper_incomplete_gamma(a,x)/gamma(a)
   
