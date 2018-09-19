# SP800_90b_tests

**A Python implementation of the non_iid tests in SP800-90B**

This is very rough right now. I've got as far as coding all the non-iid tests and they match the example vectors. I've opened it up so that the two other people in the world who like doing this sort of this can play with it and maybe be motivated to improve it.

**to run**

python sp800_90b_tests.py mybinarydatafile.bin

If you execute the individual test files, they will run the example vectors from the spec.

**Issues**
There are problems:
1) It's slow. There are optimizations that can be made.
2) The LRS test runs out of memory on windows when the dictionary reaches 11184810 keys (this is on a Win10 machine with 16GB of memory). It did fine on Linux. This is really a shortcoming of the LRS test which for reasonable input data sizes, creates a huge dictionary.[Fixed]
3) Its probably buggy.
4) Actually it is buggy LRS, MMC and so on. Changes are afoot.

**A note on IID vs. Non-IID**
I might get around to implementing the IID tests (these are tests to decide if you claim of having IID data is true). However I doubt it since I am very non-motivated. This is because no entropy source in this universe is IID and so claiming it and testing for it is just stupid. All that IID test suite and the shuffling nonsense needs to be removed from the spec.

**Why does this exist? Why not just use the NIST code?**
Because the NIST 90B software sucks. The command line interface doesn't accept multiple files and doesn't output CSV. The intent here is to mirror the very useful file handling in djent. This implementation has a normal command line interface. Point it as your binary random file. It runs the tests and gives you a summary of the entropy estimates at the end.

**Command Line Options**
There are options and stuff:

```
$ ./sp800_90b_tests.py -h
usage: sp800_90b_tests.py [-h] [--be] [-t TESTNAME] [-l SYMBOL_LENGTH]
                          [--list_tests] [--test_iid]
                          [filename]

Test data to establish an entropy estimate, using NIST SP800-90B algorithms.

positional arguments:
  filename              Filename of binary file to test

optional arguments:
  -h, --help            show this help message and exit
  --be                  Treat data as big endian bits within bytes. Defaults
                        to little endian
  -t TESTNAME, --testname TESTNAME
                        Select the test to run. Defaults to running all tests.
                        Use --list_tests to see the list
  -l SYMBOL_LENGTH, --symbol_length SYMBOL_LENGTH
                        Indicate the length of each symbol in bits
  --list_tests          Display the list of tests
  --test_iid            Run Tests of IID Assumption (section 5)
```
