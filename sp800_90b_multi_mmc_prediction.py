#!/usr/bin/env python

# sp_800_90b_lz78y.py
#

from __future__ import print_function
from __future__ import division

import math
from common_functions import *

precision = 300

def bits_to_int(bits):
    theint = 0
    for i in range(len(bits)):
        theint = (theint << 1) + bits[i]
    return theint
 
def bits_to_int(bits):
    theint = 0
    for i in range(len(bits)):
        theint = theint + (bits[i] << i)
        #theint = (theint << 1) + bits[i]
    return theint

def int_to_bits(s,l):
    thebits=list()
    for i in range(l):
        thebits.append(s & 0x01)
        s = s >> 1
    return thebits

def p_local_func(p,r,N):
    q = (1.0-p)
    
    x = 1.0
    for i in range(1,11):
        x = 1.0 + q*(p**r)*(x**(r+1))
    #vprint(verbose,"     x : ",x)
    result = (1.0-(p*x))/((r+1.0-(r*x))*q)
    result = result / (x**(N+1))
    return result

    
def multi_mmc_prediction(bits,symbol_length=1,verbose=True,D=16):
    vprint(verbose,"MULTI MMC PREDICTION Test")
    bitcount = len(bits)
    L = bitcount//symbol_length
    
    #vprint(verbose,bits)
    vprint(verbose,"    Symbol Length        ",symbol_length)
    vprint(verbose,"    Number of bits       ",(L * symbol_length))
    vprint(verbose,"    Number of Symbols    ",L)

    # Split bits into integer symbols
    #   prepend with 0, so the symbols are indexed from 1
    #vprint(verbose,bits)
    S = [0,] + [ bits_to_int(bits[symbol_length*i:symbol_length*(i+1)]) for i in range(L)]
    #vprint(verbose,S)
    #Step 1
    N = L-2
    subpredict = [None for x in range(D+1)] # add one to start index at one
    entries = [0 for x in range(D+1)]
    maxEntries = 100000
    correct = [0 for x in range(N+1)]
    
    vprint(verbose,"    D                    ",D)
    vprint(verbose,"    L                    ",L)
    vprint(verbose,"    N                    ",N)
    
    #step 2
    M = [dict() for x in range(D+1)]
    
    #step 3
    scoreboard = [0 for x in range(D+1)]
    winner = 1
    
    vprint(verbose,"    STEP 4")
    # step 4
    ys = list()
    for i in range(3,L+1):
        for d in range(1,D+1):
            if d < (i-1):
                x = S[i-d-1:i-1]
                y = S[i-1]
                atuple = (tuple(x),y)
                
                if atuple in M:
                    M[d][atuple] += 1
                else:
                    if entries[d] < maxEntries:
                        M[d][atuple]=1
                        entries[d] += 1
                        ys.append(y)
        for d in range(1,D+1):
            if d<i:
                # find y corresponding to highest M[Si-d,...,Si-2,y]
                ymax = -10
                maxtuple = None
                for atuple in M[d]:
                    if M[d][atuple] > ymax:
                        maxtuple = atuple
                        ymax = M[d][atuple]
                    else:
                        M[d][atuple] == ymax
                        if atuple[1] > maxtuple[1]:
                            maxtuple = atuple
                            ymax = M[d][atuple]
                subpredict[d] = ymax
                allzero = True
                for atuple in M[d]:
                    if M[d][atuple] != 0:
                        allzero = False
                        break
                if allzero:
                    subpredict[d] = None
                
        prediction = subpredict[winner]
        if prediction == S[i]:
            correct[i-2] = 1
        
        #update scoreboard
        for d in range(1,D+1):
            if subpredict[d] == S[i]:
                scoreboard[d] += 1
                if scoreboard[d] >= scoreboard[winner]:
                    winner = d
    vprint(verbose,"    STEP 5")
    # step 5
    C = 0
    for c in correct:
        if c == 1:
            C+=1
    
    vprint(verbose,"    STEP 6")
    # step 6
    p_global = float(C)/float(N)
    if (p_global == 0):
        p_prime_global = 1-(0.001**(1.0/N))
    else:
        p_prime_global = min(1.0,p_global + (2.576 * math.sqrt( (p_global*(1.0-p_global))/(N-1.0))))
    
    vprint(verbose,"    p_global             ", p_global)
    vprint(verbose,"    p_prime_global       ", p_prime_global)
    
    rlen = 0
    currentlen = 0
    for x in correct:
        if (x!=1):
            currentlen = 0
        else:
            currentlen += 1
            if currentlen > rlen:
                rlen = currentlen
    r = 1+rlen  
    vprint(verbose,"    C                    ",C)
    vprint(verbose,"    r                    ",r)
    
    vprint(verbose,"    STEP 7")
    # Step 7
    p_local = search_for_p(r,N,iterations=1000, min_plocal=0.0, max_plocal=1.0, tolerance=0.00000001,verbose=False) 
    
    vprint(verbose,"    p_local              ", p_local)

    vprint(verbose,"    STEP 8")
    # Step 8
    pu = max(p_prime_global,p_local, 1.0/(2**symbol_length))
    min_entropy_per_symbol = -math.log(pu,2)
    min_entropy_per_bit = min_entropy_per_symbol/symbol_length

    vprint(verbose,"    pu                   ",pu)
    vprint(verbose,"    Symbol Min Entropy   ",min_entropy_per_symbol)
    vprint(verbose,"    Min Entropy per bit  ",min_entropy_per_bit)

    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [2, 1, 3, 2, 1, 3, 1, 3, 1]
    
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = multi_mmc_prediction(bits,verbose=True, symbol_length=2,D=3)
    
    print("min_entropy = ",min_entropy)
          

