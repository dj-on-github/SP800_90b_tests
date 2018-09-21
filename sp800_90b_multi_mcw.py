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

def pfunc(plocal,r,N):
    q = 1.0-plocal
    
    # Find x10
    x = 0.0
    for j in range(1,11):
        x = 1.0 + (q*(plocal**r)*(x**(r+1.0)))

    # do the equation
    result = (1.0 - plocal*x)
    result = result/((r+1.0 - (r*x))*q)
    result = result/(x**(N+1))

    return result

def multi_mcw(bits,symbol_length=1, verbose=True, ws = [0,63,255,1023,4095]):
    vprint(verbose,"MULTI MCW Test")
    bitcount = len(bits)
    L = bitcount//symbol_length

    #vprint(verbose,bits)
    vprint(verbose,"   Symbol Length           ",symbol_length)
    vprint(verbose,"   Number of bits          ",(L * symbol_length))
    vprint(verbose,"   Number of Symbols       ",L)

    # Split bits into integer symbols
    symbols = [ bits_to_int(bits[symbol_length*i:symbol_length*(i+1)]) for i in range(L)]
    #vprint(verbose,symbols) 

    #Steps 1
    w = ws # Window Sizes
    N = L-w[1]
    vprint(verbose,"   N                       ",N)
    correct = [0 for i in range(N+1)]
    
    # Step 2
    scoreboard = [0,0,0,0,0]
    frequent = [0,None, None, None, None]
    winner = 1
    prediction = None

    # Step 3
    symbols = [0,]+symbols
    #for i in range(w[1]+1,L+1):
    #vprint(verbose,"   i     frequent                    scoreboard3b    winner  prediction   si      correct[i-w[1]] scoreboard3d")
    for i in range(w[1]+1,L+1):
        for j in [1,2,3,4]:
            if (i > w[j]):
                counts = dict()
                tiebreaker = 1
                #print ("RANGE: ",list(range(i-w[j],i)),"  Bits :",[symbols[x] for x in range(i-w[j],i)])
                for index in range(i-w[j],i):
                    s = symbols[index]
                    if s in counts:
                        (c,t) = counts[s]
                        c += 1
                        t = tiebreaker
                        tiebreaker += 1
                        counts[s] = (c,t)
                    else:
                        t = tiebreaker
                        tiebreaker += 1
                        counts[s] = (1,t)
                #vprint(verbose,"Counts : ",counts)
                # find max frequency
                themax = 0
                for s in counts:
                    (c,t) = counts[s]
                    if c > themax:
                        themax = c
                #vprint(verbose,"MAX COUNT:",themax)
                # use the tiebreaker
                themax_tiebreaker = 0
                for s in counts:
                    (c,t) = counts[s]
                    if c == themax:
                        #vprint(verbose,"IF ",t,">",themax_tiebreaker, "answer=",(t > themax_tiebreaker))
                        if t > themax_tiebreaker:
                            #vprint(verbose,"  T > THEMAX_TIEBREAKER  t:",t,"   tmt:",themax_tiebreaker)
                            themax_tiebreaker = t
                            most_frequent_symbol = s
                            #vprint(verbose,"  NOW T = THEMAX_TIEBREAKER  t:",t,"  tmt:",themax_tiebreaker)

                    #vprint(verbose," TIEBREAKER: s=",s,"  count = ",c,"  t=",t," max_tieb:",themax_tiebreaker," most_freq_s:",most_frequent_symbol)
                # set frequent[j] to the most frequent and recent symbol
                frequent[j] = most_frequent_symbol
            else:
                frequent[j] = None
        
        prediction = frequent[winner]
        #scoreboard3b = scoreboard[:]
        if (prediction == symbols[i]):
            correct[i-w[1]] = 1
        for j in [1,2,3,4]:
            if (frequent[j] == symbols[i]):
                scoreboard[j] += 1
                if scoreboard[j] >= scoreboard[winner]:
                    winner = j
        
        #scoreboard3d = scoreboard[:]
        #vprint(verbose,"  ",str(i).ljust(5),str(frequent[1:]).ljust(27),str(scoreboard3b[1:]).ljust(15),
        #      str(winner).ljust(7),str(prediction).ljust(12),str(symbols[i]).ljust(7),
        #      str(correct[i-w[1]]).ljust(15),str(scoreboard3d[1:]).ljust(12),)
    #vprint(verbose,"   Correct                 ",correct)
    # Step 4
    C = 0
    for i in correct:
        if i==1:
            C += 1
 
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
    
    # Binary chop search for Plocal
    P_local = search_for_p(r,N,verbose=verbose)

    vprint(verbose,"   P_local                 ",P_local)
    k = 2.0**symbol_length
    min_entropy = -math.log(max(P_prime_global,P_local,1.0/k),2)
    min_entropy_per_bit = min_entropy/symbol_length

    vprint(verbose,"   Min Entropy per symbol  ",min_entropy)
    vprint(verbose,"   Min Entropy per bit     ",min_entropy_per_bit)
    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [1,2,1,0,2,1,1,2,2,0,0,0]
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = multi_mcw(bits,symbol_length=2,ws=[0,3,5,7,9])
    
    vprint(verbose,"min_entropy = ",min_entropy)
