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

    
def lz78y(bits,symbol_length=1,verbose=True,B=16):
    vprint(verbose,"LZ78Y Test")
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
    N = L-B-1
    
    vprint(verbose,"    B                    ",B)
    vprint(verbose,"    N                    ",N)
    correct = [0 for x in range(N+1)]
    maxDictionarySize = 65536
    
    # Step 2
    D = dict()
    dictionarySize = 0
    
    # Step 3
    if verbose:
        vprint(verbose,"    ","i".ljust(4),"Add to D".ljust(20),"prev".ljust(14),"Max D[prev]".ljust(16),
                "prediction".ljust(12),"Si".ljust(4),"Correct_i-b-1")
    for i in range(B+2,L+1):
        add_to_d = list()
        prevlist = list()
        maxdlist = list()
        
        for j in range(B,0,-1):
            # 3a
            ss = tuple(S[i-j-1:(i-2)+1])
            
            if (ss not in D) and (dictionarySize < maxDictionarySize):
                D[ss] = dict()
                D[ss][S[i-1]] = 0
                add_to_d.append("D["+str(ss)+"]["+str(S[i-1])+"]")
                dictionarySize += 1
            if (ss in D):
                if S[i-1] not in D[ss]:
                    D[ss][S[i-1]] = 0 
                    add_to_d.append("D["+str(ss)+"]["+str(S[i-1])+"]")
                D[ss][S[i-1]] = D[ss][S[i-1]]+1
            
        # 3b
        prediction = None
        maxcount = None
        for j in range(B,0,-1):
            prev = tuple(S[i-j:(i-1)+1])
            prevlist.append(str(prev))
            if prev in D:
                maxyval = 0
                
                for cy in range(2**symbol_length):
                    if cy in D[prev]:
                        if D[prev][cy] > maxyval:
                            maxyval = D[prev][cy]
                            y = cy
                if (maxcount == None) or (D[prev][y] > maxcount):
                    prediction = y
                    maxcount = D[prev][y]
                    maxdlist.append(maxcount)
        if prediction == S[i]:
            correct[i-B-1] = 1
        
        # print out table line
        #vprint(verbose,add_to_d)
        #vprint(verbose,prevlist)
        #vprint(verbose,maxdlist)
        #if verbose:
        #    for pad in range(20):
        #        add_to_d.append("-")
        #        prevlist.append("-")
        #        maxdlist.append("-")
        #    for line in range(4):
        #        if line == 0:
        #            vprint(verbose,"    ",str(i).ljust(4),add_to_d[line].ljust(20), prevlist[line].ljust(14), str(maxcount).ljust(16),
        #                        str(prediction).ljust(12),str(S[i]).ljust(4), correct[i-B-1])
        #        else:
        #            vprint(verbose,"    "," ".ljust(4),add_to_d[line].ljust(20), prevlist[line].ljust(14), str(maxcount).ljust(16),
        #                        " ".ljust(12)," ".ljust(4), " ")
    # step 4
    C = sum(correct)
    #vprint(verbose,"    correct              ",correct)
    p_global = float(C)/float(N)
    if (p_global == 0):
        p_prime_global = 1-(0.001**(1.0/N))
    else:
        p_prime_global = min(1.0,p_global + (2.576 * math.sqrt( (p_global*(1.0-p_global))/(N-1.0))))
    
    vprint(verbose,"    p_global             ", p_global)
    vprint(verbose,"    p_prime_global       ", p_prime_global)
    
    # Step 5
     #  Find run of longest ones in correct, to find r
    
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
    
    vprint(verbose,"    r                    ", r)
    
    #   iteratively find Plocal 
    p_local = search_for_p(r,N,iterations=1000, min_plocal=0.0, max_plocal=1.0, tolerance=0.00000001,verbose=False)  
                
    vprint(verbose,"    p_local              ", p_local)

    # Step 6
    pu = max(p_prime_global,p_local, 1.0/(2**symbol_length))
    min_entropy_per_symbol = -math.log(pu,2)
    min_entropy_per_bit = min_entropy_per_symbol/symbol_length

    vprint(verbose,"    pu                   ",pu)
    vprint(verbose,"    Symbol Min Entropy   ",min_entropy_per_symbol)
    vprint(verbose,"    Min Entropy per bit  ",min_entropy_per_bit)

    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [2,1,3,2,1,3,1,3,1,2,1,3,2]
    
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = lz78y(bits,symbol_length=2,B=4)
    
    vprint(verbose,"min_entropy = ",min_entropy)
          

