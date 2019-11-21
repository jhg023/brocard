#include <cfloat>
#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <unistd.h>

// Start at 1 trillion
//#define STARTING_N 1094393000000
constexpr uint64_t STARTING_N = 1;

// End at 5 trillion
//#define ENDING_N 5000000000000
constexpr uint64_t ENDING_N = 1'000'000'000;

// Milestone used for printing progress (5 billion)
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
constexpr int MULMOD_DIFFERENCE = 2'000'000;

struct range_struct {
  int tid;
  uint64_t start;
  uint64_t end;
  uint64_t *factorials;
  const uint64_t *primes;
  const uint64_t *pinvs;
};

static inline uint64_t ll_mod_preinv( uint64_t a_hi, uint64_t a_lo, uint64_t n, uint64_t ninv ) {
  uint64_t q0, q1, r;

  const int norm = __builtin_clzll( n );

  n <<= norm;
  a_hi <<= norm;

  // We don't need r_shift, as 'norm' will never be 0
  //const uint64_t u1 = a_hi + r_shift( a_lo, FLINT_BITS - norm );
  const uint64_t u1 = a_hi + ( a_lo >> ( FLINT_BITS - norm ) );
  const uint64_t u0 = ( a_lo << norm );

  umul_ppmm( q1, q0, ninv, u1 );
  add_ssaaaa( q1, q0, q1, q0, u1, u0 );

  r = ( u0 - ( q1 + 1 ) * n );

  if( r > q0 ) {
    r += n;
  }

  if( __builtin_expect( r < n, 1 ) ) {
    return r >> norm;
  }

  // I've never seen the above condition be false, so this return statement
  // may be able to be removed.
  return ( r - n ) >> norm;
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
  mp_limb_t r, s;

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

static inline int jacobi_unsigned( uint64_t x, uint64_t y ) {
  uint64_t b = x, a = y, temp;
  int s;

  uint exp = __builtin_ctzll( b );
  b >>= exp;

  bool first = (( exp * ( a * a - 1 ) ) & 8) != 0;
  bool second = (( ( a - 1 ) * ( b - 1 ) ) & 4) != 0;

  if (first != second) {
    s = -1;
  } else {
    s = 1;
  }

  while( b != 1 ) {
    if( ( a >> 2 ) < b ) {
      temp = a - b;
      a = b;

      if( temp < b ) {
        b = temp;
      } else if( temp < ( b << 1 ) ) {
        b = temp - a;
      } else {
        b = temp - ( a << 1 );
      }
    } else {
      temp = a % b;
      a = b;
      b = temp;
    }

    if ( b == 0 ) {
      return 0;
    }

    exp = __builtin_ctzll( b );
    b >>= exp;

    first = (( exp * ( a * a - 1 ) ) & 8) != 0;
    second =  (( ( a - 1 ) * ( b - 1 ) ) & 4) != 0;

    if (first != second) {
      s = -s;
    }
  }

  return s;
}

static inline void *brocard( void *arguments ) {
  auto *range = static_cast<struct range_struct *>( arguments );

  int tid = range->tid;
  uint64_t start = range->start;
  uint64_t end = range->end;
  uint64_t *factorials = range->factorials;
  const uint64_t *primes = range->primes;
  const uint64_t *pinvs = range->pinvs;

  uint64_t last_n[NUM_PRIMES] = { 0 };

  for( uint64_t &i : last_n ) {
    i = start - 1;
  }

  int best_i = 25, i;
  uint64_t n, prime, pinv;

  for( n = start; n <= end; ++n ) {
    for( i = 0; i < NUM_PRIMES; ++i ) {
      prime = primes[i];
      pinv = pinvs[i];

      if( last_n[i] == n - 1 ) {
        factorials[i] = mulmod_preinv( factorials[i], n, prime, pinv );
      } else if( n - last_n[i] <= MULMOD_DIFFERENCE ) { // Allow for underflow
        for( uint64_t j = last_n[i] + 1; j <= n; ++j ) {
          factorials[i] = mulmod_preinv( factorials[i], j, prime, pinv );
        }
      } else {
        factorials[i] = initialize_factorial( n, prime, pinv );
      }

      last_n[i] = n;

      if( jacobi_unsigned( factorials[i] + 1, prime ) == -1 ) {
        break;
      }
    }

    if( __builtin_expect( static_cast<long>( i == NUM_PRIMES ), 0 ) != 0 ) {
      printf( "[Sub Range #%d] Potential Solution: %llu, primes[0] = %llu, factorials[0] = %llu\n", tid, n, primes[0], factorials[0] );
      FILE *fp = fopen( SOLUTION_FILE_NAME, "ae" );
      fprintf( fp, "%llu\n", n );
      fclose( fp );
    } else if( __builtin_expect( static_cast<long>( i >= best_i ), 0 ) != 0 ) {
      best_i = i;
      printf( "[Sub Range #%d] Progress: %llu (%.2f%%), Tests Passed: %d\n", tid, n, 100.0 * tid / NUM_SUB_RANGES, best_i );
    } else if( __builtin_expect( static_cast<long>( n % MILESTONE == 0 ), 0 ) != 0 ) {
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

  uint64_t starting_n = STARTING_N;

  for( int i = 0; i < NUM_SUB_RANGES; ++i ) {
    auto *range = static_cast<struct range_struct *>( malloc( sizeof( struct range_struct ) ) );

    range->tid = i;
    range->start = starting_n;
    range->end = ( i == NUM_SUB_RANGES - 1 ) ? ENDING_N : starting_n + partition_size;
    range->primes = generate_primes( range->end, NUM_PRIMES );
    range->pinvs = generate_pinvs( range->primes, NUM_PRIMES );

    auto *factorials = static_cast<uint64_t *>( calloc( NUM_PRIMES, sizeof( uint64_t ) ) );
    uint64_t n = starting_n - 1;

#pragma omp parallel for num_threads( FACTORIAL_NUM_THREADS ) schedule( dynamic )
    for( int i = 0; i < NUM_PRIMES; ++i ) {
      factorials[i] = initialize_factorial( n, range->primes[i], range->pinvs[i] );
    }

    range->factorials = factorials;

    ranges[i] = *range;
    starting_n += partition_size + 1;
  }

#pragma omp parallel for schedule( dynamic )
  for( int i = 0; i < NUM_SUB_RANGES; ++i ) {
    brocard( &ranges[i] );
  }
}