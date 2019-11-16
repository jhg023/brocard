#include <celero/Celero.h>
#include <cstdint>
#include <flint/flint.h>
#include <flint/ulong_extras.h>
#include <vector>

using celero::DoNotOptimizeAway;
using std::vector;

uint64_t mulmod64( uint64_t a, uint64_t b, uint64_t n ) {
  uint64_t d;
  uint64_t unused; // dummy output, unused, to tell GCC that RAX register is modified by this snippet
  asm( "mulq %3\n\t"
       "divq %4"
       : "=a"( unused ), "=&d"( d )
       : "a"( a ), "rm"( b ), "rm"( n )
       : "cc" );
  return d;
}

uint64_t largestPower( uint64_t n, uint64_t p ) {
  // Initialize result
  uint64_t x = 0;

  // Calculate x = n/p + n/(p^2) + n/(p^3) + ....
  while( n ) {
    n /= p;
    x += n;
  }
  return x;
}

// Utility function to do modular exponentiation.
// It returns (x^y) % p
uint64_t power( uint64_t x, uint64_t y, uint64_t p ) {
  uint64_t res = 1; // Initialize result
  x = x % p;        // Update x if it is more than or
  // equal to p
  while( y > 0 ) {
    // If y is odd, multiply x with result
    if( y & 1 )
      res = ( res * x ) % p;

    // y must be even now
    y = y >> 1; // y = y/2
    x = ( x * x ) % p;
  }
  return res;
}

// Function to find modular inverse of a under modulo p
// using Fermat's method. Assumption: p is prime
uint64_t modInverse( uint64_t a, uint64_t p ) {
  return power( a, p - 2, p );
}

uint64_t naive_mulmod( uint64_t n, uint64_t p ) {
  uint64_t result = 1;
  for( uint64_t i; i <= n; ++i ) {
    result = mulmod64( result, i, p );
  }
  return result;
}

// Returns n! % p
uint64_t naive_sieve( uint64_t n, uint64_t p ) {
  vector<bool> A( n, true );

  // Use Sieve of Eratosthenes to find all primes
  // smaller than n
  for( uint64_t i = 2; i <= static_cast<uint64_t>( sqrt( n ) ); ++i ) {
    if( A[i] ) {
      for( uint64_t j = i * i; j < n; j += i ) {
        A[j] = false;
      }
    }
  }

  uint64_t res = 1;

  // Consider all primes found by Sieve
  for( uint64_t i = 2; i <= n; i++ ) {
    if( A[i] ) {
      // Find the largest power of prime 'i' that divides n
      uint64_t k = largestPower( n, i );

      // Multiply result with (i^k) % p
      res = ( res * power( i, k, p ) ) % p;
    }
  }
  return res;
}

// Returns n! % p using Wilson's Theorem
uint64_t naive_wilson( uint64_t n, uint64_t p ) {
  // Initialize result as (p-1)! which is -1 or (p-1)
  uint64_t res = ( p - 1 );

  // Multiply modulo inverse of all numbers from (n+1)
  // to p
  for( uint64_t i = n + 1; i < p; i++ )
    res = ( res * modInverse( i, p ) ) % p;
  return res;
}

#define I_START 1'000'000'000
#define P_START 1'000'004'981

#define I_END 1'000'000'000
#define P_END 1'000'004'981

#define NUM_RUNS 1
#define NUM_ITERATIONS 3

BASELINE( FacBench, naive_mulmod, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t i = I_START; i <= I_END; ++i ) {
    for( uint64_t p = P_START; p <= P_END; ++p ) {
      DoNotOptimizeAway( naive_mulmod( i, p ) );
    }
  }
}

BENCHMARK( FacBench, flint, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t i = I_START; i <= I_END; ++i ) {
    for( uint64_t p = P_START; p <= P_END; ++p ) {
      DoNotOptimizeAway( n_factorial_fast_mod2_preinv( i, p, n_preinvert_limb( p ) ) );
    }
  }
}

/*
BENCHMARK( FacBench, naive_sieve, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t i = I_START; i <= I_END; ++i ) {
    for( uint64_t p = P_START; p <= P_END; ++p ) {
      DoNotOptimizeAway( naive_sieve( i, p ) );
    }
  }
}
*/

BENCHMARK( FacBench, naive_wilson, NUM_RUNS, NUM_ITERATIONS ) {
  for( uint64_t i = I_START; i <= I_END; ++i ) {
    for( uint64_t p = P_START; p <= P_END; ++p ) {
      DoNotOptimizeAway( naive_wilson( i, p ) );
    }
  }
}

CELERO_MAIN
