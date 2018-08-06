# SP800_90b_tests
A Python implementation of the non_iid tests in SP800-90B. 

This is very rough right now. I've got as far as coding all the non-iid tests and they match the example vectors. I've opened it up so that the two other people in the world who like doing this sort of this can play with it and maybe be motivated to improve it.

to run:

python sp800_90b_tests.py mybinarydatafile.bin

If you execute the individual test files, they will run the example vectors from the spec.

There are problems:
1) It's slow. There are optimizations that can be made.
2) The LRS test runs out of memory on windows when the dictionary reaches 11184810 keys (this is on a Win10 machine with 16GB of memory). It did fine on Linux. This is really a shortcoming of the LRS test which for reasonable input data sizes, creates a huge dictionary.
3) Its probably buggy.

A note on IID vs. Non-IID:
I might get around to implementing the non-IID tests (these are tests to decide if you claim of having IID data is true). However I doubt it since I am very non-motivated. This is because no entropy source in this universe is IID and so claiming it and testing for it is just stupid. All that non-IID test suite and the shuffling nonsense needs to be removed from the spec.

Why does this exist? Why not just use the NIST code? 
Because the NIST 90B software sucks. It is buggy and has a terrible text mode menus, questions and answers user interface. This implementation has a normal command line interface. Point it as your random file. It runs the tests and gives you a summary of the entropy estimates at the end.
