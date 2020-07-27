# Settlement violation estimates for Proof-of-Stake blockchains with concurrent honest slot leaders under the longest-chain rule
Compute the relative margin in time and space O(T^3) where T is the length of the characteristic string. 
See Section 6.6 of [the paper](https://eprint.iacr.org/2020/041.pdf) for reference. 

We remark that computing these probabilities on a set of size 2^T would take time/space exponential in T if done naively.

## Author
Saad Quader and Alexander Russell

## Description
The code is in C++. The following files should be built as executables:
* `forkability_table.cpp` 
* `prob_forkable.cpp`

Look at the comments around the `main()` function in each file. They are fairly self-explanatory. 

The file `prob_forkable.cpp`, when built into an executable, is an interactive program which asks for
* `N`, a positive integer, the length of the characteristic string `w = w_1 ... w_N` where the symbol `w_i` is one of `h, H, A`; 
* `eps`, a real between 0 and 1, so that `Pr[w_i = A] = (1 - eps)/2` independently for each `i = 1, ... , N`; and
* `r`, a real between 0 and 1, so that `r = Pr[w_i = h]/Pr[w_i != A]`.
It outputs an upper bound on the probability that `w` is forkable, i.e., slot 1 is not `N`-settled.

The file `forkability_table.cpp`, when built into an executable, is an interactive program which requests an output file name and then reproduces the table in Section 6.6 of [the paper](https://eprint.iacr.org/2020/041.pdf). Look inside the `main()` function to see which values of `N`, `eps`, and `r` are being used.


If building on Windows in debug mode, the header `win_mem_leak.h` helps detect memory leak, if any.
