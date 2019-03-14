#!/usr/bin/env python

# sp_800_90b_multi_mwc.py
#

from __future__ import print_function
from __future__ import division

import math
import operator as op
from functools import reduce
from common_functions import *

def nCr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom
    
def bad_nCr(n,r):
    result = math.factorial(n)
    result = result / math.factorial(r)
    result = result / math.factorial(n-r)
    return result

def bits_to_int(bits):
    theint = 0
    for c,i in enumerate(range(len(bits))):
        theint = theint  + (bits[i] << c)
    return theint

def int_to_bits(s,l):
    thebits=list()
    for i in range(l):
        thebits.append(s & 0x01)
        s = s >> 1
    return thebits

#def pfunc(plocal,r,N):
#    q = 1.0-plocal
#    
#    # Find x10
#    x = 0.0
#    for j in range(1,11):
#        x = 1.0 + (q*(plocal**r)*(x**(r+1.0)))
#
#    # do the equation
#    result = (1.0 - plocal*x)
#    result = result/((r+1.0 - (r*x))*q)
#    result = result/(x**(N+1))
#
#    return result

def lag_prediction(bits,symbol_length=1,verbose=True, D=128):
    vprint(verbose,"LAG PREDICTION Test")
    bitcount = len(bits)
    L = bitcount//symbol_length

    #vprint(verbose,bits)
    vprint(verbose,"   Symbol Length           ",symbol_length)
    vprint(verbose,"   Number of bits          ",(L * symbol_length))
    vprint(verbose,"   Number of Symbols       ",L)

    # Split bits into integer symbols
    # Prefix with 0 to start index at 1.
    s = [0,]+[ bits_to_int(bits[symbol_length*i:symbol_length*(i+1)]) for i in range(L)]
    #vprint(verbose,symbols) 

    #Steps 1
    #w = ws # Window Sizes
    N = L-1
    vprint(verbose,"   N                       ",N)
    lag = [None for i in range(D+1)] # add to to base from 1.
    correct = [0 for i in range(N+1)]
    
    # Step 2
    scoreboard = [0 for i in range(D+1)]
    frequent = [0,None, None, None, None]
    winner = 1
    #prediction = None

    # Step 3
    for i in range(2,L+1):
        for d in range(1,D+1):
            if (d < i):
                lag[d] = s[i-d]
            else:
                lag[d] = None
        prediction = lag[winner]
        if (prediction == s[i]):
            correct[i-1] = 1
        for d in range(1,D+1):
            if lag[d]==s[i]:
                scoreboard[d]=scoreboard[d]+1
                if scoreboard[d] >= scoreboard[winner]:
                    winner = d
    
    # Step 4
    C = 0
    for i in correct:
        if i==1:
            C += 1
    #print ("correct = ",correct)
    # Step 5
    P_global = C/N
    if P_global == 0:
        P_prime_global = 1.0 -(0.01**(1.0/N))
    else:
        P_prime_global = min(1.0,P_global + (2.576*math.sqrt((P_global*(1.0-P_global)/(N-1.0))))) 

    vprint(verbose,"   P_global                ",P_global)
    vprint(verbose,"   P_prime_global          ",P_prime_global)
    
    # Step 6
    
    # find longest run of ones in correct[]
    runlength = 0
    max_runlength = 0
    for c in correct:
        if c == 1:
            runlength += 1
        else:
            runlength = 0

        if runlength > max_runlength:
            max_runlength = runlength

    r = max_runlength+1
    vprint(verbose,"    C                    ",C)
    vprint(verbose,"    r                    ",r)
    
    #solve_for_p(mu_bar=0.99, n=N, v=r, tolerance=1e-09)
    P_local = search_for_p(r,N,verbose=verbose)
    
    if False:
        # Binary chop search for Plocal
        iterations = 1000
        iteration = 0
        min_plocal = -0.1
        max_plocal = 1.1
        found = False
        while (iteration < 1000):
            candidate = (min_plocal + max_plocal)/2.0
            result = pfunc(candidate,r,N)
            iteration += 1
            if iteration > iterations:
                found = False
                break
            elif (result > (0.99-0.00000001)) and (result < (0.99+0.00000001)):
                found = True
                P_local = candidate
                break
            elif result > 0.99:
                min_plocal = candidate
            else:
                max_plocal = candidate

        if (found == False):
            print ("Warning: P_local not found")

    vprint(verbose,"   P_local                 ",P_local)
    k = 2.0**symbol_length
    min_entropy = -math.log(max(P_prime_global,P_local,1.0/k),2)
    min_entropy_per_bit = min_entropy/symbol_length

    vprint(verbose,"   Min Entropy per symbol  ",min_entropy)
    vprint(verbose,"   Min Entropy per bit     ",min_entropy_per_bit)
    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [2,1,3,2,1,3,1,3,1,2]
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = lag_prediction(bits,symbol_length=2,D=3)
    
    vprint(verbose,"min_entropy = ",min_entropy)
