#!/usr/bin/env python

# sp_800_90b_ttuple.py
#


from __future__ import print_function
from __future__ import division

import math

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
 
def ttuple(bits,symbol_length=1, threshold=35):
    print("T-TUPLE Test")
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

    #Steps 1 and 2
    # Find-t
    # The t-tuple length for which the count is at least 35
    tuple_dict = dict()
    max_count = None
    max_tuple = None
    Q = [0 for x in range(1024)]  # Large enough to always be big enough
    P = [0 for x in range(1024)]  # Large enough to always be big enough
    P_max_array = [0 for x in range(1024)]  # Large enough to always be big enough
    last_five_maxes = [threshold+100 for i in range(5)]  # Keep track of the last 10. If they were all one,
                                             # end to loop to save compute time.
    for t in range(1,min(L+1,128)):  # (max_count == None) or (max_count > threshold):
        max_count = 0
        max_tuple = None
        print ("   Testing t=",t,end="") 
        tuple_position_count = 1+L-t
        #print("   Searching through ",tuple_position_count," positions")


        for i in range(tuple_position_count):
            the_tuple = tuple(symbols[i:i+t])
            if the_tuple in tuple_dict:
                tuple_dict[the_tuple] += 1
            else:
                tuple_dict[the_tuple] = 1

            if tuple_dict[the_tuple] > max_count:
                max_count = tuple_dict[the_tuple]
                max_tuple = the_tuple
            #print ("   Found ",the_tuple," at location ",i," count = ",tuple_dict[the_tuple])
        Q[t]=max_count
        last_five_maxes = last_five_maxes[1:]
        last_five_maxes.append(max_count)
        print("   max tuple count: ",max_count)
        if (max(last_five_maxes)==1) or (max(last_five_maxes) < (threshold-10)):
            break
        #print("   Q[t] = ",max_count, "  Q[i]=",Q[1:t+1])

    for pos,qt in reversed(list(enumerate(Q[:L+1]))):
        #print("   pos=",pos, "  qt=",qt)
        if qt >= threshold:
            found = True
            t = pos
            break

    if found:
        print("   Found t = ",t)
    else:
        print("   Error, no t found")
        exit()

    # Step 2
    for i in range(1,t+1):
            P[i] = Q[i]/(L-i+1.0)
            P_max_array[i] = P[i]**(1.0/i)
    p_max = max(P_max_array)

    pu = min(1.0,p_max + (2.576*math.sqrt((p_max*(1.0-p_max)/(L-1.0))))) 

    min_entropy_per_symbol = -math.log(pu,2.0)
    min_entropy = min_entropy_per_symbol/symbol_length

    print("   pu                   ",pu)
    print("   Symbol Min Entropy   ",min_entropy_per_symbol)
    print("   Min Entropy per bit  ",min_entropy)

    return (False, None, min_entropy)

if __name__ == "__main__":
    bits = list()
    symbols = [2, 2, 0, 1, 0, 2, 0, 1, 2, 1, 2, 0, 1, 2, 1, 0, 0, 1, 0, 0, 0]
    for s in symbols:
        bits = bits + int_to_bits(s,2)
    (iid_assumption,T,min_entropy) = ttuple(bits,symbol_length=2,threshold=3)
    
    print("min_entropy = ",min_entropy)
