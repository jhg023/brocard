#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <unistd.h>
#include <cstdint>

constexpr uint64_t STARTING_N = 1;
constexpr uint64_t ENDING_N = 1'000'000'000;

// Milestone used for printing progress
constexpr uint64_t MILESTONE = 100'000'000;

// The number of threads to use when computing the initial factorial values.
// This is a memory hog, so a small amount of threads should be used here to
// avoid running out of memory.
constexpr int FACTORIAL_NUM_THREADS = 8;

// The name of the file to write potential solutions to.
#define SOLUTION_FILE_NAME "brocard_solutions.txt"

// The number of primes to use when testing.
constexpr int NUM_PRIMES = 40;

// The amount of sub-ranges that the range (ENDING_N - STARTING_N) should be partitioned into.
constexpr int NUM_SUB_RANGES = 32;

// If 'last_n[i] - n >= MULMOD_DIFFERENCE', then a more efficient method will be used
// to catch up 'last_n[i]' instead of repeatedly calling 'mulmod_preinv'.
constexpr int MULMOD_DIFFERENCE = 3'000'000;

struct range_struct {
  int tid;
  uint64_t start;
  uint64_t end;
  uint64_t *factorials;
  const uint64_t *primes;
  const uint64_t *pinvs;
};

static inline uint64_t ll_mod_preinv( uint64_t a_hi, uint64_t a_lo, uint64_t n, uint64_t ninv ) {
  const int norm = __builtin_clzll( n );

  n <<= norm;
  a_hi <<= norm;

  // We don't need r_shift, as 'norm' will never be 0
  //const uint64_t u1 = a_hi + r_shift( a_lo, FLINT_BITS - norm );
  const uint64_t u1 = a_hi + ( a_lo >> ( FLINT_BITS - norm ) );
  const uint64_t u0 = ( a_lo << norm );

  uint64_t q0, q1;

  umul_ppmm( q1, q0, ninv, u1 );
  add_ssaaaa( q1, q0, q1, q0, u1, u0 );

  uint64_t r = ( u0 - ( q1 + 1 ) * n ) + n;

  return ( r < n ) ? r >> norm : ( r - n ) >> norm;
}

static inline uint64_t mulmod_preinv( uint64_t a, uint64_t b, uint64_t n, uint64_t ninv ) {
  uint64_t p1, p2;
  umul_ppmm( p1, p2, a, b );
  return ll_mod_preinv( p1, p2, n, ninv );
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
  for( i = 1; i < m; i++ )
    t[i] = n_submod( t[i - 1], UWORD( 1 ), p );

  _nmod_poly_product_roots_nmod_vec( u, t, m, mod );

  for( i = 0; i < m; i++ )
    t[i] = n_mod2_preinv( i * m + 1, p, pinv );

  _nmod_poly_evaluate_nmod_vec_fast( v, u, m + 1, t, m, mod );

  r = 1;
  for( i = 0; i < m; i++ )
    //r = n_mulmod2_preinv(r, v[i], mod.n, mod.ninv);
    r = mulmod_preinv( r, v[i], mod.n, mod.ninv );

  for( s = m * m + 1; s <= n; s++ )
    //r = n_mulmod2_preinv(r, s, mod.n, mod.ninv);
    r = mulmod_preinv( r, s, mod.n, mod.ninv );

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
static inline int jacobi_modified( uint64_t a, uint64_t b ) {
  if ( __builtin_expect( a == b, 0 ) ) {
    return 0;
  }

  b >>= 1;

  int c = __builtin_ctzll( a );
  int bit = c & ( b ^ ( b >> 1 ) );
  
  a >>= c;
  a >>= 1;

  int64_t t;
  uint64_t bgta;

  do {
    t = a - b;
    bgta = t >> 63;
    bit ^= ( bgta & a & b );
    b += ( bgta & t );
    a = ( t ^ bgta ) - bgta;
    c = __builtin_ctzll( t );
    ++c;
    bit ^= c & ( b ^ ( b >> 1 ) );
    a >>= c;
  } while( b != 0 );

  return bit & 1;
}

static inline void *brocard( void *arguments ) {
  auto *range = static_cast<struct range_struct *>( arguments );

  const int tid = range->tid;
  const uint64_t start = range->start;
  const uint64_t end = range->end;
  uint64_t *factorials = range->factorials;
  const uint64_t *primes = range->primes;
  const uint64_t *pinvs = range->pinvs;

  int current_i;
  uint best_i = 25, i;
  uint64_t n, prime, pinv;

  for( n = start; n <= end; ++n ) {
    current_i = -1;

    for( i = 0; i < NUM_PRIMES; ++i ) {
      prime = primes[i];
      pinv = pinvs[i];

      factorials[i] = mulmod_preinv( factorials[i], n, prime, pinv );

      if( current_i < 0 && jacobi_modified( factorials[i] + 1, prime ) ) {
        current_i = i;
      }
    }

    if( __builtin_expect( current_i == -1, 0 ) ) {
      printf( "[Sub Range #%d] Potential Solution: %llu, primes[0] = %llu, factorials[0] = %llu\n", tid, n, primes[0], factorials[0] );
      FILE *fp = fopen( SOLUTION_FILE_NAME, "ae" );
      fprintf( fp, "%llu\n", n );
      fclose( fp );
    } else if( __builtin_expect( current_i >= best_i, 0 ) ) {
      best_i = current_i;
      printf( "[Sub Range #%d] Progress: %llu (%.2f%%), Tests Passed: %d\n", tid, n, 100.0 * tid / NUM_SUB_RANGES, best_i );
    } else if( __builtin_expect( n % MILESTONE == 0, 0 ) ) {
      printf( "[Sub Range #%d] Progress: %llu (%.2f%%)\n", tid, n, 100.0 * tid / NUM_SUB_RANGES );
    }
  }

  //printf( "[Thread #%d] Completed execution!\n", tid );
  return nullptr;
}

/**
 * Generates and returns the first 'num_primes' primes after 'start' as a pointer.
 */
auto generate_primes( uint64_t start, uint64_t num_primes ) -> uint64_t * {
  n_primes_t iter;
  n_primes_init( iter );
  n_primes_jump_after( iter, start );

  auto *primes = static_cast<uint64_t *>( calloc( num_primes, sizeof( uint64_t ) ) );

  for( uint64_t i = 0; i < num_primes; ++i ) {
    primes[i] = n_primes_next( iter );
  }

  n_primes_clear( iter );

  return primes;
}

auto generate_pinvs( const uint64_t *primes, uint64_t num_primes ) -> uint64_t * {
  auto *pinvs = static_cast<uint64_t *>( calloc( num_primes, sizeof( uint64_t ) ) );

  for( uint64_t i = 0; i < num_primes; ++i ) {
    pinvs[i] = n_preinvert_limb( primes[i] );
  }

  return pinvs;
}

auto main() -> int {
  remove( SOLUTION_FILE_NAME );

  uint64_t partition_size = ( ENDING_N - STARTING_N ) / NUM_SUB_RANGES;
  printf( "Partition Size: %llu\n", partition_size );

  struct range_struct ranges[NUM_SUB_RANGES];

  for( int i = 0; i < NUM_SUB_RANGES; ++i ) {
    auto *range = static_cast<struct range_struct *>( malloc( sizeof( struct range_struct ) ) );

    range->tid = i;
    range->start = STARTING_N + ( i * partition_size );
    range->end = ( i == NUM_SUB_RANGES - 1 ) ? ENDING_N : range->start + partition_size - 1;
    range->primes = generate_primes( range->end, NUM_PRIMES );
    range->pinvs = generate_pinvs( range->primes, NUM_PRIMES );

    auto *factorials = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );
    uint64_t n = range->start - 1;

#pragma omp parallel for num_threads( FACTORIAL_NUM_THREADS ) schedule( dynamic )
    for( int i = 0; i < NUM_PRIMES; ++i ) {
      factorials[i] = initialize_factorial( n, range->primes[i], range->pinvs[i] );
    }

    range->factorials = factorials;

    ranges[i] = *range;
  }

#pragma omp parallel for schedule( dynamic )
  for( int i = 0; i < NUM_SUB_RANGES; ++i ) {
    brocard( &ranges[i] );
  }
}
