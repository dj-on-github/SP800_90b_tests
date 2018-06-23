#!/usr/bin/env python

# sp_800_90b_lrs.py
#

from __future__ import print_function
from __future__ import division

import math
import operator as op
from functools import reduce

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
    for i in range(len(bits)):
        theint = (theint << 1) + bits[i]
    return theint

def int_to_bits(s,l):
    thebits=list()
    for i in range(l):
        thebits.append(s & 0x01)
        s = s >> 1
    return thebits
 
def lrs(bits,symbol_length=1, threshold=35):
    print("LRS Test")
    bitcount = len(bits)
    L = bitcount//symbol_length

    #print(bits)
    print("   Symbol Length        ",symbol_length)
    print("   Number of bits       ",(L * symbol_length))
    print("   Number of Symbols    ",L)
    print("   t-threshold = ",threshold)

    # Split bits into integer symbols
    symbols = [ bits_to_int(bits[symbol_length*i:symbol_length*(i+1)]) for i in range(L)]
    #print(symbols) 

    #Steps 1
    # Find-
    # The smallest u-tuple length for which the count is less than 35
    tuple_dict = dict()
    max_count = None
    max_tuple = None
    for u in range(1,min(L+1,128)):  # (max_count == None) or (max_count > threshold):
        max_count = 0
        max_tuple = None
        print ("   Testing u=",u,end="") 
        tuple_position_count = 1+L-u
        #print("   Searching through ",tuple_position_count," positions")

        for i in range(tuple_position_count):
            the_tuple = tuple(symbols[i:i+u])
            if the_tuple in tuple_dict:
                tuple_dict[the_tuple] += 1
            else:
                tuple_dict[the_tuple] = 1

            if tuple_dict[the_tuple] > max_count:
                max_count = tuple_dict[the_tuple]
                max_tuple = the_tuple
        print("   max tuple count: ",max_count)
        
        # Breakout condition
        if max_count < threshold:
            found_u = u
            break
    max_count = 0
       
    # Step 2
    tuple_dict = dict()
    last_max=threshold+100
    for v in range(1,min(L+1,128)):
        last_max = max_count
        max_count = 0
        max_tuple = None
        print ("   Testing v=",v,end="") 
        tuple_position_count = 1+L-v
        #print("   Searching through ",tuple_position_count," positions")

        for i in range(tuple_position_count):
            the_tuple = tuple(symbols[i:i+v])
            if the_tuple in tuple_dict:
                tuple_dict[the_tuple] += 1
            else:
                tuple_dict[the_tuple] = 1

            if tuple_dict[the_tuple] > max_count:
                max_count = tuple_dict[the_tuple]
                max_tuple = the_tuple
        print("   max tuple count: ",max_count)
        
        # Breakout condition
        if (last_max > 1) and (max_count==1):
            found_v = v
            break
    v = found_v
    
    # Step 3
    P = [0.0 for x in range(v+1)]
    P_max = [0.0 for x in range(v+1)]
    for W in range(u,v+1):
        C=list()
        C.append(0) # Zeroth element ignored
        ith_unique_W_tuple = list()
        ith_unique_W_tuple_count = list()
        tuple_dict=dict()
        
        for i in range(1+L-W):
            the_tuple = tuple(symbols[i:i+W])
            if the_tuple in tuple_dict:
                tuple_dict[the_tuple] += 1
                #print(ith_unique_W_tuple_count)
                ith_unique_W_tuple_count[-1] += 1
            else:
                tuple_dict[the_tuple] = 1
                ith_unique_W_tuple.append(the_tuple)
                ith_unique_W_tuple_count.append(1)

        C = ith_unique_W_tuple_count
        #print("   C = ",C)
        p_max = [0.0 for x in range(W)]
        for c in C:
            if (2 <= c):
                P[W] += nCr(c,2)
        P[W] = P[W]/(nCr(L-W+1,2))
        P_max[W]=P[W]**(1.0/W)
    #print("pmax[u,v]       ",P_max[u:v+1])
    p_hat = max(P_max[u:v+1])
    print("   p_hat                ",p_hat)
    
    # Step 4
    pu = min(1.0,p_hat + (2.576*math.sqrt((p_hat*(1.0-p_hat)/(L-1.0))))) 
        
    # Step 5
    min_entropy_per_symbol = -math.log(pu,2)
    min_entropy_per_bit = min_entropy_per_symbol/symbol_length

    print("   pu                   ",pu)
    print("   Symbol Min Entropy   ",min_entropy_per_symbol)
    print("   Min Entropy per bit  ",min_entropy_per_bit)

    return (False, None, min_entropy_per_bit)

if __name__ == "__main__":
    bits = list()
    symbols = [2, 2, 0, 1, 0, 2, 0, 1, 2, 1, 2, 0, 1, 2, 1, 0, 0, 1, 0, 0, 0]
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = lrs(bits,symbol_length=2,threshold=3)
    
    print("min_entropy = ",min_entropy)
