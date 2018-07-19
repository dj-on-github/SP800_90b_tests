#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import math
import operator as op
from functools import reduce

import random

def new_nCr(n,r):
    a=1
    b=1
    k = r
    if r > (n-r):
        k = n-r
    i = n
    while k >= i:
        a = a*i
        b = b*i
        if (a % b)==0:
            a=a/b
            b = 1
        k -= 1
        i -= 1
    return a/b
        
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

for i in range(100):
    n = random.randrange(100)
    r = random.randrange(100)
    if (r > n):
        x = r
        r = n
        n = x
    a = new_nCr(n,r)
    b = nCr(n,r)
    c = bad_nCr(n,r)
    
    if (a != b) or (b != c) or (a != c):
        print(a,b,c," mismatch")
    else:
        print(a,b,c)