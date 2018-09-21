#!/usr/bin/env python

# sp_800_90b_compression.py
#

from __future__ import print_function
from __future__ import division

import math
from common_functions import *

def bits_to_int(bits):
    theint = 0
    for i in range(len(bits)):
        theint = (theint << 1) + bits[i]
    return theint
    
def F(z,t,u):
    if u < t:
        return (z**2.0)*((1.0-z)**(u-1.0))
    if u == t:
        return z*((1.0-z)**(t-1.0))

# The equations in step 7 of 6.3.4 are downright misleading and do not work.
# This function more or less follows what NIST did in their code but it looks
# nothing like the equations in the spec.
def G(z, v, d, L):
    g_sum = 0.0
    st = [math.log(u, 2.0) * ((1.0-z)**(u-1.0)) for u in range((d+1), v+d+1)]
    g_sum = v*z*z * sum([math.log(u, 2.0) * ((1.0-z)**(u-1.0)) for u in range(1,(d+1))])
    g_sum += z*z * sum([(v-t-1) * st[t] for t in range(v-1)])
    g_sum += z * sum(st)
    return g_sum/v
        
def compression(bits,symbol_length=1,verbose=True, d=1000):
    vprint(verbose,"COMPRESSION Test")
    L = len(bits)

    if symbol_length != 1:
        vprint(verbose,"   Warning, Compression test treats data at 1 bit symbols. Setting symbol length to 1")

    
    #vprint(verbose,bits)
    vprint(verbose,"   Symbol Length        1")
    vprint(verbose,"   Number of bits      ",L)

    # step 1
    b = 6
    blocks = L//b
    s_prime = [0,]+[bits_to_int(bits[b*i:b*(i+1)]) for i in range(blocks)]

    vprint(verbose,"   Number of blocks    ",blocks)

    # Step 2
    dict_data = s_prime[1:d+1]
    v = blocks-d
    test_data=s_prime[d+1:]

    vprint(verbose,"   v                   ",v)

    # Step 3
    dictionary = [0 for i in range((2**b)+1)] # Make it 1 bigger and leave the zero element dangling
                                              # so the indexes match the spec which uses 1 based arrays.
    for i in range(1,d+1):
        dictionary[s_prime[i]]=i

    # Step 4
    D = [0,]+[0 for i in range(v)]
    for i in range(d+1,blocks+1):
        #vprint(verbose,"  i = ",i,end="")
        #vprint(verbose,"  s_prime[%d]=" % i,s_prime[i])
        if dictionary[s_prime[i]] != 0:
            #print ("D[i-d] = D[%d - %d] = D[%d]" % (i,d,i-d))
            D[i-d] = i-dictionary[s_prime[i]]
            dictionary[s_prime[i]] = i
        if dictionary[s_prime[i]] == 0:
            dictionary[s_prime[i]] = i
            D[i-d] = i

    # Step 5

    x_sum = 0.0
    for i in range(1,v+1):
        #vprint(verbose,"   D[",i,"] = ",D[i], "log2(D[i])=",math.log(D[i],2))
        x_sum += math.log(D[i],2)
    x_bar = x_sum/v

    vprint(verbose,"   x_bar               ",x_bar)

    c = 0.5907

    s_sum = 0.0
    for i in range(1,v+1):
        s_sum += (math.log(D[i],2)**2)
    s_sum = s_sum/(v-1.0)
    s_sum = s_sum - (x_bar**2)
    sigma_hat = c * math.sqrt(s_sum)
    
    vprint(verbose,"   sigma_hat           ",sigma_hat)

    # Step 6
    
    x_bar_prime = x_bar - ((2.576*sigma_hat)/math.sqrt(v))
    vprint(verbose,"   x_bar_prime         ",x_bar_prime)

    # Step 7

    p_min = 2.0 ** -b  # binary search bounds
    p_max = 1.0
    p_mid = (p_min+p_max)/2.0
   
    vprint(verbose,"   p_min               ",p_min) 
    vprint(verbose,"   p_max               ",p_max) 
    iterations = 1000
    iteration = 0

    found = False
    while (iteration < iterations):
        q = (1.0-p_mid)/((2.0**b)-1.0)
        candidate = G(p_mid,v,d,L) + (((2.0**b)-1.0)*G(q,v,d,L))

        if abs(candidate -x_bar_prime) < 0.00000000001:
            found = True
            break
        elif candidate > x_bar_prime:
            p_min = p_mid
            p_mid = (p_min+p_max)/2.0
        elif candidate < x_bar_prime:
            p_max = p_mid
            p_mid = (p_min+p_max)/2.0

        iteration += 1

    # Step 8
    if found:
        min_entropy = -math.log(p_mid,2)/b
        vprint(verbose,"   min_entropy =",min_entropy)
        return(False,None,min_entropy)
    else:
        min_entropy = 1.0
        vprint(verbose,"   min_entropy = 1.0")
        return(False,None,min_entropy)

if __name__ == "__main__":
    bits = [1,0,0,0,1,1,1,0,
            0,1,0,1,0,1,0,1,
            1,1,0,0,1,1,0,0,
            0,1,1,1,0,0,1,0,
            1,0,1,0,1,1,1,0,
            1,1,1,0,0,0,1,1]

    (iid_assumption,T,min_entropy) = compression(bits,1,d=4)
    
    vprint(verbose,"min_entropy = ",min_entropy)
