#include <cstdint>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <locale.h>
#include <unistd.h>

constexpr uint64_t STARTING_N = 1ULL;
constexpr uint64_t ENDING_N = 1'000'000'000ULL;

// Milestone used for printing progress
constexpr uint64_t MILESTONE = 100'000'000;

// The number of threads to use when computing the initial factorial values.
// This is a memory hog, so a small amount of threads should be used here to
// avoid running out of memory.
constexpr uint FACTORIAL_NUM_THREADS = 8;

// The name of the file to write potential solutions to.
#define SOLUTION_FILE_NAME "brocard_solutions.txt"

// The number of primes to use when testing.
// 
// 30 = 1 in 1 billion of finding a potential solution
// 40 = 1 in 1 trillion of finding a potential solution
// 50 = 1 in 1 quadrillion of finding a potential solution
constexpr uint NUM_PRIMES = 40;

// The amount of sub-ranges that the range (ENDING_N - STARTING_N) should be partitioned into.
constexpr uint NUM_SUB_RANGES = 32; //( ENDING_N - STARTING_N ) / 39'062'500ULL;

// If 'last_n[i] - n >= MULMOD_DIFFERENCE', then a more efficient method will be used
// to catch up 'last_n[i]' instead of repeatedly calling 'mulmod_preinv'.
constexpr uint MULMOD_DIFFERENCE = 2'000'000;

struct range_struct {
  uint tid;
  uint64_t start;
  uint64_t end;
  uint64_t *factorials;
  const uint64_t *bits;
  const uint64_t *norms;
  const uint64_t *pinvs;
  const uint64_t *primes;
  const uint64_t *primes_shifted;
};

static inline uint64_t mulmod_preinv( uint64_t a, uint64_t b, uint64_t n, uint64_t ninv, uint64_t norm ) {
  __uint128_t prod = ( __uint128_t ) a * b;
  
  uint64_t a_hi = prod >> 64;
  uint64_t a_lo = ( uint64_t ) prod;

  a_hi <<= norm;

  const __uint128_t u1 = a_hi + ( a_lo >> ( FLINT_BITS - norm ) );
  const uint64_t u0 = ( a_lo << norm );

  __uint128_t u = u1 << 64 | u0;

  prod = ( ( ninv * u1 ) + u ) >> 64;

  uint64_t r = ( u0 - ( ( prod + 1 ) * n ) ) + n;

  return ( r < n ) ? r >> norm : ( r - n ) >> norm;
}

static inline uint64_t factorial_fast_mod2_preinv( uint64_t n, uint64_t p, uint64_t pinv ) {
  slong i, m;
  nmod_t mod;
  mp_ptr t, u, v;
  uint64_t r, s;

  nmod_init( &mod, p );

  m = n_sqrt( n );

  t = _nmod_vec_init( m + 1 );
  u = _nmod_vec_init( m + 1 );
  v = _nmod_vec_init( m + 1 );

  t[0] = UWORD( 0 );

  for( i = 1; i < m; ++i ) {
    t[i] = n_submod( t[i - 1], UWORD( 1 ), p );
  }

  _nmod_poly_product_roots_nmod_vec( u, t, m, mod );

  for( i = 0; i < m; ++i ) {
    t[i] = n_mod2_preinv( i * m + 1, p, pinv );
  }

  _nmod_poly_evaluate_nmod_vec_fast( v, u, m + 1, t, m, mod );

  r = 1;

  for( i = 0; i < m; ++i ) {
    r = n_mulmod2_preinv( r, v[i], mod.n, mod.ninv );
  }

  for( s = m * m + 1; s <= n; ++s ) {
    r = n_mulmod2_preinv( r, s, mod.n, mod.ninv );
  }

  _nmod_vec_clear( t );
  _nmod_vec_clear( u );
  _nmod_vec_clear( v );

  return r;
}

static inline uint64_t initialize_factorial( uint64_t n, uint64_t prime, uint64_t pinv ) {
  if( n < ( prime >> 1 ) ) {
    return factorial_fast_mod2_preinv( n, prime, pinv );
  }

  uint64_t factorial = factorial_fast_mod2_preinv( prime - n - 1, prime, pinv );

  factorial = n_invmod( factorial, prime );

  if( ( n & 1 ) == 0 ) {
    factorial = -factorial + prime;
  }

  return factorial % prime;
}

// A modified version of the jacobi symbol (a/b).
// Returns 0 if 'a == b' or '(a/b) == 1', and 1 if '(a/b) == -1'.
// Assertions: a <= b, b is an odd prime
static inline int jacobi_modified( uint64_t a, uint64_t b, uint64_t bit ) {
  if ( __builtin_expect( a == b, 0 ) ) {
    return 0;
  }

  b >>= 1;

  uint c = __builtin_ctzll( a );
  bit &= c;

  a >>= c;
  a >>= 1;

  do {
    #pragma unroll( 1 )
    for ( uint i = 0; i < 1; ++i ) {
      int64_t t = a - b;

      /* If b > a, invoke reciprocity */
      bit ^= ( a >= b ? 0 : a & b );

      /* b <-- min (a, b) */
      b = a < b ? a : b;

      /* a <-- |a - b| */
      a = ( ( t < 0 ) ? -t : t );

      c = __builtin_ctzll( a ) + 1;

      bit ^= c & ( b ^ ( b >> 1 ) );

      a >>= c;
    }
  } while( b > 0 );

  return bit & 1;
}

static inline void *brocard( void *arguments ) {
  auto *range = static_cast<struct range_struct *>( arguments );

  const uint tid = range->tid;
  const uint64_t start = range->start;
  const uint64_t end = range->end;
  const uint64_t *bits = range->bits;
  const uint64_t *norms = range->norms;
  const uint64_t *pinvs = range->pinvs;
  const uint64_t *primes = range->primes;
  const uint64_t *primes_shifted = range->primes_shifted;

  uint64_t *factorials = range->factorials;
  uint64_t last_n[NUM_PRIMES] = { 0 };

  for( uint64_t &i: last_n ) {
    i = start - 1;
  }

  uint best_i = 25, i, result;
  uint64_t n, norm, prime, pinv, prime_shifted;

  for( n = start; n <= end; ++n ) {
    for( i = 0; i < NUM_PRIMES; ++i ) {
      norm = norms[i];
      pinv = pinvs[i];
      prime = primes[i];
      prime_shifted = primes_shifted[i];

      if( n - last_n[i] <= MULMOD_DIFFERENCE ) { // Allow for underflow
        for( uint64_t j = last_n[i] + 1; j <= n; ++j ) {
          factorials[i] = mulmod_preinv( factorials[i], j, prime_shifted, pinv, norm );
        }
      } else {
        factorials[i] = initialize_factorial( n, prime, pinv );
      }

      last_n[i] = n;

      result = jacobi_modified( factorials[i] + 1, prime, bits[i] );

      if( result ) {
        break;
      }
    }

    if( __builtin_expect( i == NUM_PRIMES, 0 ) ) {
      printf( "[Sub Range #%'d] Potential Solution: %'llu - primes[0] = %'llu - factorials[0] = %'llu\n", tid, n, primes[0], factorials[0] );
      FILE *fp = fopen( SOLUTION_FILE_NAME, "ae" );
      fprintf( fp, "%'llu\n", n );
      fclose( fp );
    } else if( __builtin_expect( i >= best_i, 0 ) ) {
      best_i = i;
      printf( "[Sub Range #%'d] Progress: %'llu (%.2f%%), Tests Passed: %d\n", tid, n, 100.0 * tid / NUM_SUB_RANGES, best_i );
    } else if( __builtin_expect( n % MILESTONE == 0, 0 ) ) {
      printf( "[Sub Range #%'d] Progress: %'llu (%.2f%%)\n", tid, n, 100.0 * tid / NUM_SUB_RANGES );
    }
  }

  //printf( "[Thread #%d] Completed execution!\n", tid );
  return nullptr;
}

/**
 * Generates and returns the first 'NUM_PRIMES' primes after 'start' as a pointer.
 */
auto generate_primes( uint64_t start ) -> uint64_t * {
  n_primes_t iter;
  n_primes_init( iter );
  n_primes_jump_after( iter, start );

  uint64_t *primes = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );

  for( int i = 0; i < NUM_PRIMES; ++i ) {
    primes[i] = n_primes_next( iter );
  }

  n_primes_clear( iter );

  return primes;
}

auto generate_norms( const uint64_t *primes ) -> uint64_t * {
  auto *norms = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );

  for ( int i = 0; i < NUM_PRIMES; ++i ) {
    norms[i] = __builtin_clzll( primes[i] );
  }

  return norms;
}

auto shift_primes( const uint64_t *primes, const uint64_t *norms ) -> uint64_t * {
  auto *primes_shifted = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );

  for ( int i = 0; i < NUM_PRIMES; ++i ) {
    primes_shifted[i] = primes[i] << norms[i];
  }

  return primes_shifted;
}

auto generate_bits( const uint64_t *primes ) -> uint64_t * {
  auto *bits = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );

  for ( int i = 0; i < NUM_PRIMES; ++i ) {
    bits[i] = ( primes[i] >> 1 ) ^ ( primes[i] >> 2 );
  }

  return bits;
}

auto generate_pinvs( const uint64_t *primes ) -> uint64_t * {
  auto *pinvs = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );

  for( int i = 0; i < NUM_PRIMES; ++i ) {
    pinvs[i] = n_preinvert_limb( primes[i] );
  }

  return pinvs;
}

auto main() -> int {
  setlocale( LC_NUMERIC, "" );

  remove( SOLUTION_FILE_NAME );

  uint64_t partition_size = ( ENDING_N - STARTING_N ) / NUM_SUB_RANGES;
  printf( "Number of sub-ranges: %'d\n", NUM_SUB_RANGES );
  printf( "Size of each sub-range: %'llu\n", partition_size );
  printf( "Beginning to calculate %'d factorials...\n", NUM_SUB_RANGES * NUM_PRIMES );

  struct range_struct ranges[NUM_SUB_RANGES];

  for( uint i = 0; i < NUM_SUB_RANGES; ++i ) {
    if ( i > 0 && i % 10 == 0 ) {
      printf( "Finished calculating factorials for %'d sub-ranges out of %'d\n", i, NUM_SUB_RANGES );
    }

    auto *range = static_cast<struct range_struct *>( malloc( sizeof( struct range_struct ) ) );

    range->tid = i;
    range->start = STARTING_N + ( i * partition_size );
    range->end = ( i == NUM_SUB_RANGES - 1 ) ? ENDING_N : range->start + partition_size - 1;
    range->primes = generate_primes( range->end );
    range->norms = generate_norms( range->primes );
    range->primes_shifted = shift_primes( range->primes, range->norms );
    range->bits = generate_bits( range->primes );
    range->pinvs = generate_pinvs( range->primes );

    auto *factorials = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );
    uint64_t n = range->start - 1;

#pragma omp parallel for num_threads( FACTORIAL_NUM_THREADS ) schedule( dynamic )
    for( uint i = 0; i < NUM_PRIMES; ++i ) {
      factorials[i] = initialize_factorial( n, range->primes[i], range->pinvs[i] );
    }

    range->factorials = factorials;

    ranges[i] = *range;
  }

  printf( "Finished calculating factorials for all %'d sub-ranges!\n", NUM_SUB_RANGES );

#pragma omp parallel for schedule( dynamic )
  for( uint i = 0; i < NUM_SUB_RANGES; ++i ) {
    brocard( &ranges[i] );
  }
}
