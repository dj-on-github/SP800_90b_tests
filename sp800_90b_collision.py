#!/usr/bin/env python

# sp_800_90b_collision.py
#

from __future__ import print_function
from __future__ import division

import math

from gamma_functions import *

def pq_func(p):
    q = 1.0-p

    z = 1.0/q

    fq = upper_incomplete_gamma(3,z) * (z**(-3.0)) * (e**(z))

    result = ((p*(q**-2.0)) * (1.0 + (0.5*((1/p)-(1/q)))) * fq)
    result = result - ((p*(q**-1.0)) * 0.5 * ((1/p)-(1/q)))

    return result

def collision(bits,symbol_length=1):
    print("Collision Test")
    if symbol_length > 1:
        print("Warning: Collision test is only to be run with symbol length of 1")

    t = list()

    # Step 1
    v = 0

    # Step 2
    index = 1
    found = False
    while True:
        j = index
        if bits[index-1]==bits[index]:
            found = True
            j=index+1
            #print("  ",bits[index-1:j])
        elif j > (len(bits)-2):
            found = False
            break
        else:
            found = True
            j = index+2
            #print("  ",bits[index-1:j])

        # Step 3
        if (found == True):
            v = v + 1
            t.append(j - index + 1)
            index = j+1

        # Step 4
        if index > (len(bits)-2):
            break

    #print("   T = ",t)
    # Step 5
    t_sum = 0.0
    sq_sum = 0.0
    t_sum = sum(t)
    x_bar = t_sum/v 
    print("   x_bar         ",x_bar)
 
    for ti in t:
        sq_sum += (ti-x_bar)**2
    sigma_hat = math.sqrt((1.0/(v-1))*sq_sum)
    print("   sigma_hat     ",sigma_hat)
    # Step 6
    x_bar_prime = x_bar - 2.576*(sigma_hat/math.sqrt(v))
    print("   x_bar_prime   ",x_bar_prime)

    # Step 7
    iterations = 1000
    iteration = 0
    last_p_mid = -1.0
    p_min = 0.5
    p_mid = 0.75
    p_max = 1.0

    found = False
    while (not(found)):
        candidate = pq_func(p_mid)
        if candidate > x_bar_prime:
            p_min = p_mid
            p_mid = (p_min+p_max)/2.0
            #print("   G Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
        elif candidate < x_bar_prime:
            p_max = p_mid
            p_mid = (p_min+p_max)/2.0
            #print("   L Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
        elif (candidate == x_bar_prime) or (p_mid == last_p_mid):
            found = True
            p = p_mid
            #print("   M Last =",last_p_mid," Pmid =",p_mid, " Candidate = ",candidate," tgt = ",x_bar_prime)
            break
        last_p_mid = p_mid

        iteration += 1
        if iteration > iterations:
            found = False
            break 
    # step 8

    if found:
        if (p<0.5):
            p = 1.0-p
        print("   p =",p)
        min_entropy = -math.log(p,2.0)
        print("   min_entropy =",min_entropy)
    else:
        print("   p = 0.5")
        min_entropy = 1.0
        print("   min_entropy = 1.0")

    return(False, None, min_entropy)

if __name__ == "__main__":
    bits = [1,0,0,0,1,1,1,0,0,
            1,0,1,0,1,0,1,1,1,
            0,0,1,1,0,0,0,1,1,
            1,0,0,1,0,1,0,1,0,
            1,1,1,0]

    (iid_assumption,T,min_entropy) = collision(bits,1)
    
    print("min_entropy = ",min_entropy)
