#!/usr/bin/env python

# sp_800_90b_mcv    .py
#


from __future__ import print_function
from __future__ import division

import math

def bits_to_int(bits):
    theint = 0
    for i in range(len(bits)):
        theint = (theint << 1) + bits[i]
    return theint
       
def mcv(bits,symbol_length):
    print("MCV Test")
    bitcount = len(bits)
    L = bitcount//symbol_length

    #print(bits)
    print("  Symbol Length        ",symbol_length)
    print("  Number of bits       ",(L * symbol_length))
    print("  Number of Symbols    ",L)
    # Make Frequency Table
    freq_table = list()
    for i in range(2**symbol_length):
        freq_table.append(0)

    # Build the frequency table
    # Keep track of the most frequent symbol
    biggest = 0
    biggest_symbol = 0
    for i in range(L):
        symbol_bits = bits[i*symbol_length:((i+1)*symbol_length)]
        symbol = bits_to_int(symbol_bits)
        #print (" symbol:",symbol," symbol_bits",symbol_bits)
        #print(symbol_bits,symbol)
        freq_table[symbol] += 1
        if freq_table[symbol] > biggest:
            biggest = freq_table[symbol]
            biggest_symbol = symbol

    print("  Most common symbol   ",biggest_symbol)

    # do the SP800-90b section 6.3.1 sums
    p_hat = biggest/L
    print("  p_hat                ",p_hat)

    pu = p_hat + (2.576*(math.sqrt((p_hat*(1.0-p_hat))/(L-1.0))))
    if pu > 1.0:
        pu = 1.0
    min_entropy_per_symbol = (-math.log(pu,2.0))
    min_entropy = (-math.log(pu,2.0))/symbol_length
    print("  pu                   ",pu)
    print("  Symbol Min Entropy   ",min_entropy_per_symbol)
    print("  Min Entropy per bit  ",min_entropy)

    return (False, None, min_entropy)

if __name__ == "__main__":
    bits = [0,0,0,1,0,1,1,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,0,0,1,0,0,1]
    (iid_assumption,T,min_entropy) = mcv(bits,2)
    
    print("min_entropy = ",min_entropy)
