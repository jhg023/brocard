#include <flint/flint.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_poly.h>
#include <flint/ulong_extras.h>
#include <immintrin.h>
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

// The number of primes to use when testing.
constexpr int NUM_PRIMES = 40;

// The amount of '__m256i' elements in the array of primes.
// If using AVX512, this value should be 'NUM_PRIMES / 8'.
constexpr int NUM_PRIMES_VECTORIZED = NUM_PRIMES / 4;

// The amount of sub-ranges that the range (ENDING_N - STARTING_N) should be partitioned into.
constexpr int NUM_SUB_RANGES = 32;

#define LOW_32_BITS_256 = _mm256_set1_epi64x( 0x00000000FFFFFFFFULL );

// The name of the file to write potential solutions to.
#define SOLUTION_FILE_NAME "brocard_solutions.txt"

struct range_struct {
  int tid;
  uint64_t start;
  uint64_t end;
  __m256i *factorial_vectors;
  const __m256i *prime_vectors;
  const __m256i *pinv_vectors;
};

static inline auto mulmod( __m256i a, uint64_t b, __m256i m ) -> __m256i {
  __m256i d = _mm256_setzero_si256();
  __m256i b_vector = _mm256_set1_epi64x( b );
  __m256i mp2 = _mm256_srli_epi64( m, 1 );

  __m256i high_bit = _mm256_set1_epi64x( 0x8000000000000000ULL );

  for ( int i = 0; i < 64; ++i ) {
    __m256i d_mp2_gt = _mm256_cmpgt_epi64( d, mp2 );
    d = _mm256_slli_epi64( d, 1 );
    __m256i d_sub_m = _mm256_sub_epi64( d, m );
    d = _mm256_blendv_pd( d, d_sub_m, d_mp2_gt );
    
    __m256i a_high_bit = _mm256_and_si256( a, high_bit );
    __m256i d_add = _mm256_add_epi64( d, b_vector );
    d = _mm256_blendv_pd( d, d_add, a_high_bit );

    __m256i gt = _mm256_cmpgt_epi64( d, m );
    __m256i eq = _mm256_cmpeq_epi64( d, m );
    __m256i gt_eq = _mm256_or_si256( gt, eq );
    __m256i d_sub = _mm256_sub_epi64( d, m );
    d = _mm256_blendv_pd( d, d_sub, gt_eq );

    a = _mm256_slli_epi64( a, 1 );
  }
  
  return d;
}

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

static inline __m256i initialize_factorial( uint64_t n, __m256i prime_vector, __m256i pinv_vector ) {
  uint64_t primes[4];
  uint64_t pinvs[4];

  // Extract the primes from the vector and store them in 'primes'.
  _mm256_storeu_si256( ( __m256i* ) primes, prime_vector );

  // Extract the inverted limbs from the vector and store them in 'pinvs'.
  _mm256_storeu_si256( ( __m256i* ) pinvs, pinv_vector );

  uint64_t factorials[4];

  for( int i = 0; i < 4; ++i ) {
    uint64_t prime = primes[i];
    uint64_t pinv = pinvs[i];

    if( n < ( prime >> 1 ) ) {
      factorials[i] = factorial_fast_mod2_preinv( n, prime, pinv );
    } else {
      uint64_t factorial = factorial_fast_mod2_preinv( prime - n - 1, prime, pinv );

      factorial = n_invmod( factorial, prime );

      if( ( n & 1 ) == 0 ) {
        factorial = -factorial + prime;
      }

      factorials[i] = factorial % prime;
    }
  }

  return _mm256_set_epi64x( factorials[3], factorials[2], factorials[1], factorials[0] );
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
  __m256i *factorial_vectors = range->factorial_vectors;
  const __m256i *prime_vectors = range->prime_vectors;
  const __m256i *pinv_vectors = range->pinv_vectors;

  int current_i;
  uint best_i = 25, i;
  uint64_t n;

  for( n = start; n <= end; ++n ) {
    current_i = -1;

    for( i = 0; i < NUM_PRIMES_VECTORIZED; ++i ) {
      factorial_vectors[i] = mulmod( factorial_vectors[i], n, prime_vectors[i] );

      if( current_i < 0 && jacobi_modified( factorial_vectors[i] + 1, prime_vectors[i] ) ) {
        current_i = i;
      }
    }

    if( __builtin_expect( current_i == -1, 0 ) ) {
      printf( "[Sub Range #%d] Potential Solution: %llu\n", tid, n );
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
 * Generates and returns an array of size 'NUM_PRIMES_VECTORIZED', containing vectors.
 * Each element (vector) contains four packed primes of the type 'int64_t'.
 */
auto generate_primes( uint64_t start ) -> __m256i * {
  n_primes_t iter;
  n_primes_init( iter );
  n_primes_jump_after( iter, start );

  auto *prime_vectors = static_cast<__m256i *>( calloc( NUM_PRIMES_VECTORIZED, sizeof( __m256i ) ) );

  for( int i = 0; i < NUM_PRIMES_VECTORIZED; ++i ) {
    uint64_t primes[4];

    for ( int j = 0; j < 4; ++j ) {
      primes[j] = n_primes_next( iter );
    }

    // Pack the primes into a vector.
    prime_vectors[i] = _mm256_set_epi64x( primes[3], primes[2], primes[1], primes[0] );
  }

  n_primes_clear( iter );

  return prime_vectors;
}

auto generate_pinvs( const __m256i *prime_vectors ) -> __m256i * {
  auto *pinv_vectors = static_cast<__m256i *>( calloc( NUM_PRIMES_VECTORIZED, sizeof( __m256i ) ) );

  for( uint64_t i = 0; i < NUM_PRIMES_VECTORIZED; ++i ) {
    __m256i prime_vector = prime_vectors[i];

    uint64_t primes[4];
   
    // Extract the primes from the vector and store them in 'primes'.
    _mm256_storeu_si256( ( __m256i* ) primes, prime_vector );

    // Invert limbs of each prime.
    for( uint64_t &prime: primes ) {
      prime = n_preinvert_limb( prime );
    }

    // Pack the inverted limbs back into a vector.
    pinv_vectors[i] = _mm256_set_epi64x( primes[3], primes[2], primes[1], primes[0] );
  }

  return pinv_vectors;
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
    range->prime_vectors = generate_primes( range->end );
    range->pinv_vectors = generate_pinvs( range->prime_vectors );

    auto *factorial_vectors = static_cast<__m256i *>( calloc( NUM_PRIMES_VECTORIZED, sizeof( __m256i ) ) );
    uint64_t n = range->start - 1;

#pragma omp parallel for num_threads( FACTORIAL_NUM_THREADS ) schedule( dynamic )
    for( int i = 0; i < NUM_PRIMES_VECTORIZED; ++i ) {
      factorial_vectors[i] = initialize_factorial( n, range->prime_vectors[i], range->pinv_vectors[i] );
    }

    range->factorial_vectors = factorial_vectors;

    ranges[i] = *range;
  }

#pragma omp parallel for schedule( dynamic )
  for( int i = 0; i < NUM_SUB_RANGES; ++i ) {
    brocard( &ranges[i] );
  }
}
