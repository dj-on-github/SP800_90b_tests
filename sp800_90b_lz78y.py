#!/usr/bin/env python

# sp_800_90b_lz78y.py
#

from __future__ import print_function
from __future__ import division

import math
from mpmath import *
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
    #print("     x : ",x)
    result = (1.0-(p*x))/((r+1.0-(r*x))*q)
    result = result / (x**(N+1))
    return result

    
def lz78y(bits,symbol_length=1,B=16, quiet=True):
    print("LZ78Y Test")
    bitcount = len(bits)
    L = bitcount//symbol_length
    
    #print(bits)
    print("    Symbol Length        ",symbol_length)
    print("    Number of bits       ",(L * symbol_length))
    print("    Number of Symbols    ",L)

    # Split bits into integer symbols
    #   prepend with 0, so the symbols are indexed from 1
    #print(bits)
    S = [0,] + [ bits_to_int(bits[symbol_length*i:symbol_length*(i+1)]) for i in range(L)]
    #print(S)
    #Step 1
    N = L-B-1
    
    print("    B                    ",B)
    print("    N                    ",N)
    correct = [0 for x in range(N+1)]
    maxDictionarySize = 65536
    
    # Step 2
    D = dict()
    dictionarySize = 0
    
    # Step 3
    if not quiet:
        print("    ","i".ljust(4),"Add to D".ljust(20),"prev".ljust(14),"Max D[prev]".ljust(16),
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
        #print(add_to_d)
        #print(prevlist)
        #print(maxdlist)
        if not quiet:
            for pad in range(20):
                add_to_d.append("-")
                prevlist.append("-")
                maxdlist.append("-")
            for line in range(4):
                if line == 0:
                    print("    ",str(i).ljust(4),add_to_d[line].ljust(20), prevlist[line].ljust(14), str(maxcount).ljust(16),
                                str(prediction).ljust(12),str(S[i]).ljust(4), correct[i-B-1])
                else:
                    print("    "," ".ljust(4),add_to_d[line].ljust(20), prevlist[line].ljust(14), str(maxcount).ljust(16),
                                " ".ljust(12)," ".ljust(4), " ")
    # step 4
    C = sum(correct)
    #print("    correct              ",correct)
    p_global = float(C)/float(N)
    if (p_global == 0):
        p_prime_global = 1-(0.001**(1.0/N))
    else:
        p_prime_global = min(1.0,p_global + (2.576 * math.sqrt( (p_global*(1.0-p_global))/(N-1.0))))
    
    print("    p_global             ", p_global)
    print("    p_prime_global       ", p_prime_global)
    
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
    
    print("    r                    ", r)
    
    #   iteratively fine Plocal   
    iterations = 1000
    iteration = 0
    last_p_mid = -1.0
    p_min = 0.0
    p_mid = 0.5
    p_max = 1.0

    found = False
    
    while (not(found)):
        candidate = p_local_func(p_mid,r,N)
        if candidate > 0.99:
            p_min = p_mid
            p_mid = (p_min+p_max)/2.0
            #print("   G Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
        elif candidate < 0.99:
            p_max = p_mid
            p_mid = (p_min+p_max)/2.0
            #print("   L Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
        elif (candidate == 0.99) or (p_mid == last_p_mid):
            found = True
            p_local = p_mid
            #print("   M Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
            break

        iteration += 1
        if iteration > iterations:
            found = False
            p_local = p_mid
            break
                
    print("    p_local              ", p_local)

    # Step 6
    pu = max(p_prime_global,p_local, 1.0/(2**symbol_length))
    min_entropy_per_symbol = -math.log(pu,2)
    min_entropy_per_bit = min_entropy_per_symbol/symbol_length

    print("    pu                   ",pu)
    print("    Symbol Min Entropy   ",min_entropy_per_symbol)
    print("    Min Entropy per bit  ",min_entropy_per_bit)

    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [2,1,3,2,1,3,1,3,1,2,1,3,2]
    
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = lz78y(bits,symbol_length=2,B=4)
    
    print("min_entropy = ",min_entropy)
          

