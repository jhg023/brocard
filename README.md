# brocard

The code in this repository was used in an attempt to find more solutions to [Brocard's problem](https://en.wikipedia.org/wiki/Brocard%27s_problem).

Up until now, only the first `1x10^12` (1 trillion) numbers have been tested in an attempt to find an additional solution. One of our primary goals was to increase the amount of values tested by three orders of magnitude.
 - A non-overclocked [Threadripper 3970x](https://www.amd.com/en/products/cpu/amd-ryzen-threadripper-3970x) (32 cores, 64 threads) was used to run this application over that time period.

# Prerequisites
 - [flint2](https://github.com/wbhart/flint2)
 - [GMP](https://gmplib.org/)
 - [Celero](https://github.com/DigitalInBlue/Celero) (optional, to run benchmarks)
 - [gtest](https://github.com/google/googletest) (optional, to run `mulmod` test)
 
 # Findings
 We tested the first `1x10^15` (1 quadrillion) values over a period of ~5 months (January-May, 2020), but no additional solutions were found.
 
 We used a similar algorithm to the one described in [this paper](http://unsolvedproblems.org/S99.pdf), but added additional optimizations (short-circuiting, Zen2 architecture-friendly code, etc.)
 
 The number which passed the most tests (`49/50`) with the primes that were chosen was `602,723,832,772,967`.
  - The odds of finding another number that passes that many tests is 1 in 500 trillion.
