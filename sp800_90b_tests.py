#!/usr/bin/env python3

# sp800_90_tests.py
# 
# This file is part of sp800_90_tests.
# 

from __future__ import print_function

import argparse
import sys

def read_bits_from_file(filename,bigendian,symbol_length=1,symbols=(1024*1024)):
    bitlist = list()
    bitcount = symbol_length*symbols
    bits_so_far = 0
    if filename == None:
        f = sys.stdin
    else:
        f = open(filename, "rb")
    while bits_so_far < bitcount:
        bytes = f.read(16384)
        if bytes:
            for bytech in bytes:
                if sys.version_info > (3,0):
                    byte = bytech
                else:
                    byte = ord(bytech)
                for i in range(8):
                    if bigendian:
                        bit = (byte & 0x80) >> 7
                        byte = byte << 1
                    else:
                        bit = (byte >> i) & 1
                    bitlist.append(bit)    
                    bits_so_far += 1
        else:
            break
    return bitlist

parser = argparse.ArgumentParser(description='Test data to establish an entropy estimate, using NIST SP800-90B algorithms.')
parser.add_argument('filename', type=str, nargs='?', help='Filename of binary file to test')
parser.add_argument('--be', action='store_false',help='Treat data as big endian bits within bytes. Defaults to little endian')
parser.add_argument('-t', '--testname', default=None,help='Select the test to run. Defaults to running all tests. Use --list_tests to see the list')
parser.add_argument('-l', '--symbol_length', type=int, default=1,help='Indicate the length of each symbol in bits')
parser.add_argument('--list_tests', action='store_true',help='Display the list of tests')
parser.add_argument('--test_iid', action='store_true',default=False,help='Run Tests of IID Assumption (section 5)')

args = parser.parse_args()

bigendian = args.be
filename = args.filename

# IID Assumption Testing

#   5.1.1  Excursion Test Statistic
#   5.1.2  Number of Directional Runs
#   5.1.3  Length of Directional Runs
#   5.1.4  Number of Increases and Decreases
#   5.1.5  Number of Runs Based on the Median
#   5.1.6  Length of Runs Based on Median
#   5.1.7  Average Collision Test Statistic
#   5.1.8  Maximum Collision Test Statistic
#   5.1.9  Periodicity Test Statistic
#   5.1.10 Covariance Test Statistic
#   5.1.11 Compression Test Statistic

# IID Assumption Testing: Additional Chi-squared Tests

#   5.2.1  Testing Independence for Non-Brinary Data
#   5.2.2  Testing Goodness-of-fit for Non-Binary Data
#   5.2.3  Testing Independence for Binary Data
#   5.2.4  Testing Goodness-of-fit for Binary Data
#   5.2.5  Length of the ongest Repeated Substring Test

# Entropy Estimation Tests for Non_IID Data

#   6.3.1  The Most Common Value Estimate
#   6.3.2  The Collision Estimate    (binary only)
#   6.3.3  The Markov Estimate       (binary only)
#   6.3.4  The Compression Estimate  (binary only)
#   6.3.5  The t-Tuple Estimate
#   6.3.6  The Longest Repeated Substring (LRS) Estimate
#   6.3.7  The Multi Most Common in Window Prediction Estimate
#   6.3.8  The Lag Prediction Estimate
#   6.3.9  The MultiMMC Prediction Estimate
#   6.3.10 The LZ78Y Prediction Estimate

non_iid_testlist = [
        'mcv',
        'collision',
        'markov',
        'compression',
        'ttuple',
        'lrs',
        'multi_mcv',
        #'lag_prediction',
        #'multi_mmc_prediction',
        'lz78y'
        ]


iid_testlist = [
        'excursion',
        'number_of_directional_runs',
        'length_of_directional_runs',
        'number_of_increases_and_decreases',
        'number_of_runs_based_on_median',
        'length_of_runs_based_in_median',
        'average_collision',
        'max_collision',
        'periodicity',
        'covariance',
        'compression_test_statistic',
        'mcv']   # MCV is the only min entropy estimator for IID (see section 6.1)

if args.list_tests:
    print("Testing the IID Assumption")
    for i,testname in zip(range(len(iid_testlist)),iid_testlist):
        print(str(i+1).ljust(4)+": "+testname)
    
    print("Non-IID Entropy Estimators")
    for i,testname in zip(range(len(non_iid_testlist)),non_iid_testlist):
        print(str(i+1).ljust(4)+": "+testname)
    exit()

if args.test_iid:
    testlist = iid_testlist
else:
    testlist = non_iid_testlist

symbol_length = int(args.symbol_length)
bits = read_bits_from_file(filename,bigendian)    

if args.testname:
    if args.testname in testlist:    
        m = __import__ ("sp800_90b_"+args.testname)
        func = getattr(m,args.testname)
        print("TEST: %s" % args.testname)
        (iid_assumption,T,entropy_estimate) = func(bits,symbol_length)

        if iid_assumption:
            print("IID Assumption : T = ",str(T))
        else:
            print("Min Entropy Estimate : H_inf(X) = ",str(entropy_estimate))
    else:
        print("Test name (%s) not known" % args.ttestname)
        exit()
else:
    results = list()
    me_list=list()
    t_list=list()
    for testname in testlist:
        print("TEST: %s" % testname)
        if (testname=="markov" or testname=="collision") and (symbol_length > 1):
            print("  Skipping test, it only runs on 1 bit symbols")
        else:
            m = __import__ ("sp800_90b_"+testname)
            func = getattr(m,testname)
        
            (iid_assumption,T,min_entropy) = func(bits,symbol_length)

            summary_name = testname

            if T != None:
                print("  T="+str(T))
                summary_t = str(T)
                t_list.append(T)
            else:
                summary_t = None

            if min_entropy != None:
                #print("  min_entropy="+str(min_entropy))
                summary_me = str(min_entropy)
                me_list.append(min_entropy)
            else:
                summary_me = None
        
            results.append((summary_name,summary_t, summary_me))
        
    print
    print("SUMMARY")
    print("-------")
    print("NAME".ljust(40),"T".ljust(18),"MIN ENTROPY") 
    for result in results:
        (summary_name,summary_t, summary_me) = result
        print(summary_name.ljust(40),str(summary_t).ljust(18),str(summary_me))
    min_min_entropy = min(me_list)

    print("Minimum Min Entropy = ",min_min_entropy)
    print("COMPLETE")

